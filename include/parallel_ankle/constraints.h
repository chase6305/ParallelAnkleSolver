#pragma once

#include <Eigen/Core>
#include <functional>
#include <limits>
#include <string>

namespace parallel_ankle::constraints {

enum class Status {
    Success,
    InvalidInput,
    CallbackFailure,
    InvalidEvaluation,
    NumericalFailure,
    RankDeficient,
    Inconsistent,
    BoundLimited,
    NoProgress,
    MaxIterations
};

const char* ToString(Status status);

// A local Euclidean equality problem F(x)=0. Callbacks return false when a
// configuration cannot be evaluated. They must fill every output component.
// A malformed/nonfinite output violates the callback contract.
struct EqualityProblem {
    Eigen::Index variables = 0;
    Eigen::Index equations = 0;
    std::function<bool(const Eigen::VectorXd&, Eigen::VectorXd&)> residual;
    std::function<bool(const Eigen::VectorXd&, Eigen::MatrixXd&)> jacobian;
};

struct SolveOptions {
    // Empty vectors mean unit scales or unbounded variables. All scales must
    // be positive and finite. The normalized system is
    //   r = F / residual_scales,
    //   J = diag(1/residual_scales) * F_x * diag(variable_scales).
    Eigen::VectorXd residual_scales;
    Eigen::VectorXd variable_scales;
    Eigen::VectorXd lower_bounds;
    Eigen::VectorXd upper_bounds;
    double residual_tolerance = 1e-10;
    double step_tolerance = 1e-12;
    double rank_tolerance = 1e-10;
    double damping = 1e-8;
    int max_iterations = 50;
    int max_backtracking = 16;
};

struct SolveResult {
    Status status = Status::InvalidInput;
    Eigen::VectorXd x;
    Eigen::VectorXd residual;
    double scaled_residual_norm = std::numeric_limits<double>::infinity();
    double scaled_step_norm = 0.0;
    Eigen::Index rank = 0;
    bool full_column_rank = false;
    int iterations = 0;
    int residual_evaluations = 0;
    int jacobian_evaluations = 0;
    int backtracking_steps = 0;
    std::string message;

    // Success means a sufficiently small residual, not a globally unique
    // solution. full_column_rank reports regularity of the final local model.
    bool success() const { return status == Status::Success; }
};

// Damped SVD steps with local backtracking and box projection. The initial
// state is the warm start and must lie inside all supplied bounds. This does
// not enumerate assembly modes or guarantee global convergence/uniqueness.
SolveResult SolveEqualities(const EqualityProblem& problem,
                           const Eigen::VectorXd& initial,
                           const SolveOptions& options = {});

struct DerivativeOptions {
    Eigen::VectorXd equation_scales;
    Eigen::VectorXd variable_scales;
    double rank_tolerance = 1e-10;
    // Relative to each nonzero scaled right-hand-side column, with no
    // absolute floor. An exactly zero column requires exactly zero residual.
    double consistency_tolerance = 1e-10;
};

struct DerivativeResult {
    Status status = Status::InvalidInput;
    Eigen::MatrixXd derivative;
    Eigen::Index rank = 0;
    double scaled_residual_norm = std::numeric_limits<double>::infinity();
    // Maximum per-column ||A_scaled*X_scaled+B_scaled|| / ||B_scaled||.
    // Exact zero right-hand-side/residual columns contribute zero.
    double relative_residual = std::numeric_limits<double>::infinity();
    std::string message;

    bool success() const { return status == Status::Success; }
};

// For F(x,y)=0, supply A=F_x and B=F_y at an already closed configuration.
// Solves A*(dx/dy)=-B, including redundant equations. A must have full column
// rank, and every right-hand side must be consistent within tolerance.
// The derivative is NaN on rank/consistency failure; no pseudoinverse result
// is presented as a unique derivative for an underdetermined partition.
DerivativeResult DifferentiateImplicitly(const Eigen::MatrixXd& a,
                                        const Eigen::MatrixXd& b,
                                        const DerivativeOptions& options = {});

}  // namespace parallel_ankle::constraints
