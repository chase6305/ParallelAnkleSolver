#include "parallel_ankle/constraints.h"

#include <Eigen/SVD>
#include <algorithm>
#include <cmath>
#include <exception>

namespace parallel_ankle::constraints {
namespace {

bool PositiveFinite(double value) { return std::isfinite(value) && value > 0.0; }

bool ValidScales(const Eigen::VectorXd& scales, Eigen::Index size) {
    return scales.size() == 0 ||
           (scales.size() == size && scales.allFinite() && (scales.array() > 0.0).all());
}

Eigen::VectorXd ScalesOrOnes(const Eigen::VectorXd& scales, Eigen::Index size) {
    return scales.size() == 0 ? Eigen::VectorXd::Ones(size).eval() : scales;
}

bool ValidBounds(const Eigen::VectorXd& bounds, Eigen::Index size) {
    return bounds.size() == 0 || (bounds.size() == size && !bounds.array().isNaN().any());
}

struct Evaluation {
    Status status = Status::Success;
    std::string message;
};

Evaluation EvaluateResidual(const EqualityProblem& problem, const Eigen::VectorXd& x, Eigen::VectorXd& value) {
    value = Eigen::VectorXd::Constant(problem.equations, std::numeric_limits<double>::quiet_NaN());
    try {
        if (!problem.residual(x, value)) {
            return {Status::CallbackFailure, "residual callback returned false"};
        }
    } catch (const std::exception& error) {
        return {Status::CallbackFailure, std::string("residual callback threw: ") + error.what()};
    } catch (...) {
        return {Status::CallbackFailure, "residual callback threw an unknown exception"};
    }
    if (value.size() != problem.equations || !value.allFinite()) {
        return {Status::InvalidEvaluation, "residual callback produced a wrong-sized or nonfinite output"};
    }
    return {};
}

Evaluation EvaluateJacobian(const EqualityProblem& problem, const Eigen::VectorXd& x, Eigen::MatrixXd& value) {
    value = Eigen::MatrixXd::Constant(problem.equations, problem.variables, std::numeric_limits<double>::quiet_NaN());
    try {
        if (!problem.jacobian(x, value)) {
            return {Status::CallbackFailure, "Jacobian callback returned false"};
        }
    } catch (const std::exception& error) {
        return {Status::CallbackFailure, std::string("Jacobian callback threw: ") + error.what()};
    } catch (...) {
        return {Status::CallbackFailure, "Jacobian callback threw an unknown exception"};
    }
    if (value.rows() != problem.equations || value.cols() != problem.variables || !value.allFinite()) {
        return {Status::InvalidEvaluation, "Jacobian callback produced a wrong-sized or nonfinite output"};
    }
    return {};
}

void ScaleJacobian(Eigen::MatrixXd& matrix,
                   const Eigen::VectorXd& row_scales,
                   const Eigen::VectorXd& column_scales) {
    matrix.array().colwise() /= row_scales.array();
    matrix.array().rowwise() *= column_scales.transpose().array();
}

double DampedInverse(double singular_value, double damping) {
    // sigma/(sigma^2+lambda^2), without squaring unscaled large values.
    const double scale = std::max(singular_value, damping);
    const double sigma = singular_value / scale;
    const double lambda = damping / scale;
    return (sigma / (sigma * sigma + lambda * lambda)) / scale;
}

}  // namespace

const char* ToString(Status status) {
    switch (status) {
        case Status::Success: return "success";
        case Status::InvalidInput: return "invalid input";
        case Status::CallbackFailure: return "callback failed";
        case Status::InvalidEvaluation: return "invalid callback output";
        case Status::NumericalFailure: return "nonfinite numerical operation";
        case Status::RankDeficient: return "rank-deficient local model";
        case Status::Inconsistent: return "inconsistent linearized constraints";
        case Status::BoundLimited: return "bounds prevent further progress";
        case Status::NoProgress: return "local solve made no progress";
        case Status::MaxIterations: return "iteration limit reached";
    }
    return "unknown status";
}

SolveResult SolveEqualities(const EqualityProblem& problem, const Eigen::VectorXd& initial, const SolveOptions& options) {
    SolveResult result;
    result.x = initial;
    if (problem.variables <= 0 || problem.equations <= 0 || initial.size() != problem.variables ||
        !initial.allFinite() || !problem.residual || !problem.jacobian ||
        !ValidScales(options.residual_scales, problem.equations) ||
        !ValidScales(options.variable_scales, problem.variables) ||
        !ValidBounds(options.lower_bounds, problem.variables) || !ValidBounds(options.upper_bounds, problem.variables) ||
        !PositiveFinite(options.residual_tolerance) || !PositiveFinite(options.step_tolerance) ||
        !PositiveFinite(options.rank_tolerance) || options.rank_tolerance >= 1.0 || !std::isfinite(options.damping) ||
        options.damping < 0.0 || options.max_iterations < 0 || options.max_backtracking <= 0) {
        result.message = "invalid dimensions, callbacks, initial state, scales, bounds, or solver options";
        return result;
    }
    const auto row_scales = ScalesOrOnes(options.residual_scales, problem.equations);
    const auto column_scales = ScalesOrOnes(options.variable_scales, problem.variables);
    const Eigen::VectorXd lower = options.lower_bounds.size() == 0
                                     ? Eigen::VectorXd::Constant(problem.variables, -std::numeric_limits<double>::infinity())
                                     : options.lower_bounds;
    const Eigen::VectorXd upper = options.upper_bounds.size() == 0
                                     ? Eigen::VectorXd::Constant(problem.variables, std::numeric_limits<double>::infinity())
                                     : options.upper_bounds;
    if ((lower.array() > upper.array()).any() || (initial.array() < lower.array()).any() ||
        (initial.array() > upper.array()).any()) {
        result.message = "initial state must lie inside ordered inclusive bounds";
        return result;
    }

    ++result.residual_evaluations;
    auto evaluation = EvaluateResidual(problem, result.x, result.residual);
    if (evaluation.status != Status::Success) {
        result.status = evaluation.status;
        result.message = evaluation.message;
        return result;
    }
    Eigen::VectorXd residual = result.residual.cwiseQuotient(row_scales);
    result.scaled_residual_norm = residual.stableNorm();
    if (!residual.allFinite() || !std::isfinite(result.scaled_residual_norm)) {
        result.status = Status::NumericalFailure;
        result.message = "residual scaling overflowed";
        return result;
    }

    Eigen::MatrixXd jacobian;
    Eigen::JacobiSVD<Eigen::MatrixXd> svd;
    svd.setThreshold(options.rank_tolerance);
    Eigen::VectorXd filter;
    Eigen::VectorXd step;
    Eigen::VectorXd unprojected;
    Eigen::VectorXd candidate;
    Eigen::VectorXd actual_step;
    Eigen::VectorXd candidate_residual;
    Eigen::VectorXd scaled_candidate_residual;
    bool small_step = false;
    for (int iteration = 0; iteration <= options.max_iterations; ++iteration) {
        result.iterations = iteration;
        // A successful linearization at the previous point says nothing about
        // regularity here if this iterate's callback or scaling fails.
        result.rank = 0;
        result.full_column_rank = false;
        ++result.jacobian_evaluations;
        evaluation = EvaluateJacobian(problem, result.x, jacobian);
        if (evaluation.status != Status::Success) {
            result.status = evaluation.status;
            result.message = evaluation.message;
            return result;
        }
        ScaleJacobian(jacobian, row_scales, column_scales);
        if (!jacobian.allFinite()) {
            result.status = Status::NumericalFailure;
            result.message = "Jacobian scaling overflowed";
            return result;
        }
        svd.compute(jacobian, Eigen::ComputeThinU | Eigen::ComputeThinV);
        if (svd.info() != Eigen::Success || !svd.singularValues().allFinite()) {
            result.status = Status::NumericalFailure;
            result.message = "Jacobian decomposition failed";
            return result;
        }
        result.rank = svd.rank();
        result.full_column_rank = result.rank == problem.variables;
        if (result.scaled_residual_norm <= options.residual_tolerance) {
            result.status = Status::Success;
            return result;
        }
        if (iteration == options.max_iterations) {
            result.status = Status::MaxIterations;
            return result;
        }
        if (small_step) {
            result.status = Status::NoProgress;
            result.message = "accepted step is below step tolerance while residual remains too large";
            return result;
        }
        if (result.rank == 0) {
            result.status = Status::RankDeficient;
            result.message = "zero-rank Jacobian cannot provide a local correction";
            return result;
        }

        filter.setZero(svd.singularValues().size());
        for (Eigen::Index index = 0; index < result.rank; ++index) {
            filter(index) = DampedInverse(svd.singularValues()(index), options.damping);
        }
        step = -svd.matrixV() * filter.asDiagonal() * (svd.matrixU().transpose() * residual);
        if (!step.allFinite()) {
            result.status = Status::NumericalFailure;
            result.message = "local correction overflowed";
            return result;
        }

        bool accepted = false;
        bool any_movement = false;
        bool any_valid_trial = false;
        bool any_finite_trial = false;
        bool bounds_restricted = false;
        double alpha = 1.0;
        Evaluation last_failure{Status::CallbackFailure, "no trial could be evaluated"};
        for (int attempt = 0; attempt < options.max_backtracking; ++attempt) {
            unprojected = result.x + (column_scales.array() * (alpha * step.array())).matrix();
            candidate = unprojected.cwiseMax(lower).cwiseMin(upper);
            bounds_restricted = bounds_restricted || (candidate.array() != unprojected.array()).any();
            actual_step = (candidate - result.x).cwiseQuotient(column_scales);
            const double step_norm = actual_step.stableNorm();
            const bool finite_trial = candidate.allFinite() && actual_step.allFinite() && std::isfinite(step_norm);
            any_finite_trial = any_finite_trial || finite_trial;
            if (!finite_trial) {
                last_failure = {Status::NumericalFailure, "trial coordinates or scaled step overflowed"};
            }
            if (finite_trial && step_norm > 0.0) {
                any_movement = true;
                ++result.residual_evaluations;
                evaluation = EvaluateResidual(problem, candidate, candidate_residual);
                if (evaluation.status == Status::InvalidEvaluation) {
                    result.status = evaluation.status;
                    result.message = evaluation.message;
                    return result;
                }
                if (evaluation.status == Status::Success) {
                    scaled_candidate_residual = candidate_residual.cwiseQuotient(row_scales);
                    const double norm = scaled_candidate_residual.stableNorm();
                    if (scaled_candidate_residual.allFinite() && std::isfinite(norm)) {
                        any_valid_trial = true;
                        if (norm <= options.residual_tolerance || norm < result.scaled_residual_norm) {
                            result.x = candidate;
                            result.residual = candidate_residual;
                            result.scaled_residual_norm = norm;
                            residual.swap(scaled_candidate_residual);
                            result.scaled_step_norm = step_norm;
                            small_step = step_norm <= options.step_tolerance;
                            accepted = true;
                            break;
                        }
                    } else {
                        last_failure = {Status::NumericalFailure, "trial residual scaling overflowed"};
                    }
                } else {
                    last_failure = evaluation;
                }
            }
            ++result.backtracking_steps;
            alpha *= 0.5;
        }
        if (!accepted) {
            if (!any_movement) {
                result.status = !any_finite_trial ? Status::NumericalFailure
                                                 : (bounds_restricted ? Status::BoundLimited : Status::NoProgress);
                result.message = "projection, zero correction, or floating-point resolution prevents a usable step";
            } else if (!any_valid_trial) {
                result.status = last_failure.status;
                result.message = last_failure.message;
            } else {
                result.status = Status::NoProgress;
                result.message = "no feasible trial decreases the residual";
            }
            return result;
        }
    }
    return result;
}

DerivativeResult DifferentiateImplicitly(const Eigen::MatrixXd& a,
                                        const Eigen::MatrixXd& b,
                                        const DerivativeOptions& options) {
    DerivativeResult result;
    if (a.rows() <= 0 || a.cols() <= 0 || b.rows() != a.rows() || b.cols() <= 0 || !a.allFinite() || !b.allFinite() ||
        !ValidScales(options.equation_scales, a.rows()) || !ValidScales(options.variable_scales, a.cols()) ||
        !PositiveFinite(options.rank_tolerance) || options.rank_tolerance >= 1.0 ||
        !PositiveFinite(options.consistency_tolerance)) {
        result.message = "invalid implicit derivative dimensions, data, scales, or tolerances";
        return result;
    }
    result.derivative = Eigen::MatrixXd::Constant(a.cols(), b.cols(), std::numeric_limits<double>::quiet_NaN());
    const auto row_scales = ScalesOrOnes(options.equation_scales, a.rows());
    const auto column_scales = ScalesOrOnes(options.variable_scales, a.cols());
    Eigen::MatrixXd normalized_a = a;
    ScaleJacobian(normalized_a, row_scales, column_scales);
    const Eigen::MatrixXd normalized_b = b.array().colwise() / row_scales.array();
    if (!normalized_a.allFinite() || !normalized_b.allFinite()) {
        result.status = Status::NumericalFailure;
        result.message = "implicit derivative scaling overflowed";
        return result;
    }
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(normalized_a, Eigen::ComputeThinU | Eigen::ComputeThinV);
    svd.setThreshold(options.rank_tolerance);
    if (svd.info() != Eigen::Success || !svd.singularValues().allFinite()) {
        result.status = Status::NumericalFailure;
        result.message = "implicit Jacobian decomposition failed";
        return result;
    }
    result.rank = svd.rank();
    if (result.rank != a.cols()) {
        result.status = Status::RankDeficient;
        result.message = "dependent-coordinate partition must have full column rank";
        return result;
    }
    const Eigen::MatrixXd normalized_derivative = svd.solve(-normalized_b);
    const Eigen::MatrixXd error = normalized_a * normalized_derivative + normalized_b;
    result.scaled_residual_norm = error.stableNorm();
    // Check each independent direction separately. A large consistent column
    // must not hide an inconsistent smaller right-hand side.
    result.relative_residual = 0.0;
    for (Eigen::Index column = 0; column < b.cols(); ++column) {
        const double error_norm = error.col(column).stableNorm();
        const double reference_norm = normalized_b.col(column).stableNorm();
        if (!std::isfinite(error_norm) || !std::isfinite(reference_norm)) {
            result.status = Status::NumericalFailure;
            result.message = "implicit derivative residual norm overflowed";
            return result;
        }
        if (reference_norm == 0.0) {
            if ((b.col(column).array() != 0.0).any()) {
                result.status = Status::NumericalFailure;
                result.message = "implicit right-hand-side scaling erased a nonzero column";
                return result;
            }
            if (error_norm != 0.0) {
                result.status = Status::Inconsistent;
                result.relative_residual = std::numeric_limits<double>::infinity();
                result.message = "zero implicit right-hand side produced a nonzero residual";
                return result;
            }
        } else {
            result.relative_residual = std::max(result.relative_residual, error_norm / reference_norm);
        }
    }
    const Eigen::MatrixXd derivative = column_scales.asDiagonal() * normalized_derivative;
    if (!normalized_derivative.allFinite() || !error.allFinite() || !derivative.allFinite() ||
        !std::isfinite(result.scaled_residual_norm) || !std::isfinite(result.relative_residual)) {
        result.status = Status::NumericalFailure;
        result.message = "implicit derivative solve overflowed";
        return result;
    }
    if (result.relative_residual > options.consistency_tolerance) {
        result.status = Status::Inconsistent;
        result.message = "implicit right-hand sides are not in the dependent Jacobian range";
        return result;
    }

    // Validate the actual public result after conversion back to dependent
    // coordinate units. A finite product can underflow to zero and erase a
    // direction even though its normalized solution was consistent.
    const Eigen::MatrixXd represented_normalized = derivative.array().colwise() / column_scales.array();
    const Eigen::MatrixXd represented_error = normalized_a * represented_normalized + normalized_b;
    if (!represented_normalized.allFinite() || !represented_error.allFinite()) {
        result.status = Status::NumericalFailure;
        result.message = "returned implicit derivative cannot represent the normalized solution";
        return result;
    }
    result.scaled_residual_norm = represented_error.stableNorm();
    result.relative_residual = 0.0;
    for (Eigen::Index column = 0; column < b.cols(); ++column) {
        const double error_norm = represented_error.col(column).stableNorm();
        const double reference_norm = normalized_b.col(column).stableNorm();
        if (!std::isfinite(error_norm) || (reference_norm == 0.0 && error_norm != 0.0)) {
            result.status = Status::NumericalFailure;
            result.message = "returned implicit derivative violates a zero or unrepresentable right-hand side";
            return result;
        }
        if (reference_norm != 0.0) {
            result.relative_residual = std::max(result.relative_residual, error_norm / reference_norm);
        }
    }
    if (!std::isfinite(result.scaled_residual_norm) || !std::isfinite(result.relative_residual) ||
        result.relative_residual > options.consistency_tolerance) {
        result.status = Status::NumericalFailure;
        result.message = "returned implicit derivative loses consistency when rescaled to public coordinates";
        return result;
    }
    result.derivative = derivative;
    result.status = Status::Success;
    return result;
}

}  // namespace parallel_ankle::constraints
