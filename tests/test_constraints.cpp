#include <gtest/gtest.h>

#include <cmath>
#include <limits>
#include <stdexcept>

#include "parallel_ankle/constraints.h"

namespace parallel_ankle::constraints {
namespace {

EqualityProblem AffineProblem(const Eigen::MatrixXd& a, const Eigen::VectorXd& b) {
    EqualityProblem problem;
    problem.variables = a.cols();
    problem.equations = a.rows();
    problem.residual = [a, b](const Eigen::VectorXd& x, Eigen::VectorXd& residual) {
        residual = a * x - b;
        return true;
    };
    problem.jacobian = [a](const Eigen::VectorXd&, Eigen::MatrixXd& jacobian) {
        jacobian = a;
        return true;
    };
    return problem;
}

EqualityProblem ScalarAffine(double target) {
    return AffineProblem(Eigen::MatrixXd::Ones(1, 1), Eigen::VectorXd::Constant(1, target));
}

EqualityProblem CircleAndLine() {
    EqualityProblem problem;
    problem.variables = problem.equations = 2;
    problem.residual = [](const Eigen::VectorXd& x, Eigen::VectorXd& residual) {
        residual.resize(2);
        residual << x.squaredNorm() - 1.0, x(0) - 0.6;
        return true;
    };
    problem.jacobian = [](const Eigen::VectorXd& x, Eigen::MatrixXd& jacobian) {
        jacobian.resize(2, 2);
        jacobian << 2.0 * x(0), 2.0 * x(1), 1.0, 0.0;
        return true;
    };
    return problem;
}

TEST(EqualityConstraints, RecoverCircleIntersectionsFromTheirWarmStarts) {
    for (double side : {-1.0, 1.0}) {
        const auto result = SolveEqualities(CircleAndLine(), Eigen::Vector2d(0.7, side * 0.6));
        ASSERT_TRUE(result.success()) << result.message;
        EXPECT_NEAR(result.x(0), 0.6, 1e-9);
        EXPECT_NEAR(result.x(1), side * 0.8, 1e-9);
        EXPECT_LE(result.scaled_residual_norm, 1e-10);
        EXPECT_EQ(result.rank, 2);
        EXPECT_TRUE(result.full_column_rank);
        EXPECT_GT(result.jacobian_evaluations, 0);
    }
}

TEST(EqualityConstraints, RecoverTwoLocalCrankAndRodAssemblies) {
    // Crank end B=(cos(a),sin(a)); fixed rod end C=(2,1); rod length=2.
    EqualityProblem problem;
    problem.variables = problem.equations = 1;
    problem.residual = [](const Eigen::VectorXd& x, Eigen::VectorXd& residual) {
        residual.resize(1);
        residual(0) = 2.0 - 4.0 * std::cos(x(0)) - 2.0 * std::sin(x(0));
        return true;
    };
    problem.jacobian = [](const Eigen::VectorXd& x, Eigen::MatrixXd& jacobian) {
        jacobian.resize(1, 1);
        jacobian(0, 0) = 4.0 * std::sin(x(0)) - 2.0 * std::cos(x(0));
        return true;
    };
    const auto first = SolveEqualities(problem, Eigen::VectorXd::Zero(1));
    const auto second = SolveEqualities(problem, Eigen::VectorXd::Constant(1, 2.0));
    ASSERT_TRUE(first.success()) << first.message;
    ASSERT_TRUE(second.success()) << second.message;
    EXPECT_NEAR(first.x(0), std::atan2(1.0, 2.0) - std::acos(1.0 / std::sqrt(5.0)), 1e-9);
    EXPECT_NEAR(second.x(0), 0.5 * std::acos(-1.0), 1e-9);
}

TEST(EqualityConstraints, AcceptConsistentRedundantEquations) {
    Eigen::MatrixXd a(2, 1);
    a << 1.0, 2.0;
    const auto result = SolveEqualities(AffineProblem(a, Eigen::Vector2d(2.0, 4.0)), Eigen::VectorXd::Zero(1));
    ASSERT_TRUE(result.success()) << result.message;
    EXPECT_NEAR(result.x(0), 2.0, 1e-10);
    EXPECT_EQ(result.rank, 1);
    EXPECT_TRUE(result.full_column_rank);
}

TEST(EqualityConstraints, DoNotCallAnInconsistentLeastSquaresPointASolution) {
    const Eigen::MatrixXd a = Eigen::MatrixXd::Ones(2, 1);
    const auto result = SolveEqualities(AffineProblem(a, Eigen::Vector2d(1.0, 2.0)), Eigen::VectorXd::Zero(1));
    EXPECT_FALSE(result.success());
    EXPECT_EQ(result.status, Status::NoProgress);
    EXPECT_NEAR(result.x(0), 1.5, 1e-8);
    EXPECT_GT(result.scaled_residual_norm, 0.7);
}

TEST(EqualityConstraints, SolveUnderdeterminedFeasibilityWithoutClaimingRegularity) {
    const Eigen::MatrixXd a = Eigen::MatrixXd::Ones(1, 2);
    const auto result = SolveEqualities(AffineProblem(a, Eigen::VectorXd::Constant(1, 3.0)), Eigen::Vector2d::Zero());
    ASSERT_TRUE(result.success()) << result.message;
    EXPECT_NEAR(result.x.sum(), 3.0, 1e-10);
    EXPECT_EQ(result.rank, 1);
    EXPECT_FALSE(result.full_column_rank);
}

TEST(EqualityConstraints, RespectBoundsAndReportWhenTheyPreventClosure) {
    SolveOptions options;
    options.lower_bounds = Eigen::VectorXd::Constant(1, -1.0);
    options.upper_bounds = Eigen::VectorXd::Constant(1, 1.0);
    const auto limited = SolveEqualities(ScalarAffine(2.0), Eigen::VectorXd::Zero(1), options);
    EXPECT_EQ(limited.status, Status::BoundLimited);
    EXPECT_FALSE(limited.success());
    EXPECT_DOUBLE_EQ(limited.x(0), 1.0);
    EXPECT_NEAR(limited.scaled_residual_norm, 1.0, 1e-12);
    EXPECT_EQ(SolveEqualities(ScalarAffine(0.0), Eigen::VectorXd::Constant(1, 2.0), options).status,
              Status::InvalidInput);
    options.lower_bounds(0) = options.upper_bounds(0) = 0.0;
    EXPECT_TRUE(SolveEqualities(ScalarAffine(0.0), Eigen::VectorXd::Zero(1), options).success());
}

TEST(EqualityConstraints, NormalizeDifferentResidualAndVariableUnits) {
    Eigen::Matrix2d a = Eigen::Matrix2d::Zero();
    a.diagonal() << 1e12, 1e-12;
    const Eigen::Vector2d target(2e-6, -3e6);
    SolveOptions options;
    options.variable_scales = Eigen::Vector2d(1e-6, 1e6);
    options.residual_scales = Eigen::Vector2d(1e6, 1e-6);
    const auto result = SolveEqualities(AffineProblem(a, a * target), Eigen::Vector2d::Zero(), options);
    ASSERT_TRUE(result.success()) << result.message;
    EXPECT_LT((result.x - target).cwiseQuotient(options.variable_scales).norm(), 1e-9);
    EXPECT_EQ(result.rank, 2);
}

TEST(EqualityConstraints, BacktrackThroughCallbackDomainFailures) {
    EqualityProblem problem;
    problem.variables = problem.equations = 1;
    problem.residual = [](const Eigen::VectorXd& x, Eigen::VectorXd& residual) {
        if (x(0) <= 0.0) {
            return false;
        }
        residual.resize(1);
        residual(0) = std::log(x(0)) - 1.0;
        return true;
    };
    problem.jacobian = [](const Eigen::VectorXd& x, Eigen::MatrixXd& jacobian) {
        jacobian.resize(1, 1);
        jacobian(0, 0) = 1.0 / x(0);
        return true;
    };
    const auto result = SolveEqualities(problem, Eigen::VectorXd::Constant(1, 10.0));
    ASSERT_TRUE(result.success()) << result.message;
    EXPECT_NEAR(result.x(0), std::exp(1.0), 1e-8);
    EXPECT_GT(result.backtracking_steps, 0);
}

TEST(EqualityConstraints, DistinguishMalformedAndFailedCallbacks) {
    auto problem = ScalarAffine(2.0);
    problem.residual = [](const Eigen::VectorXd&, Eigen::VectorXd&) { return false; };
    EXPECT_EQ(SolveEqualities(problem, Eigen::VectorXd::Zero(1)).status, Status::CallbackFailure);
    problem.residual = [](const Eigen::VectorXd&, Eigen::VectorXd&) -> bool { throw std::runtime_error("domain"); };
    EXPECT_EQ(SolveEqualities(problem, Eigen::VectorXd::Zero(1)).status, Status::CallbackFailure);
    problem.residual = [](const Eigen::VectorXd&, Eigen::VectorXd& residual) {
        residual = Eigen::VectorXd::Zero(2);
        return true;
    };
    EXPECT_EQ(SolveEqualities(problem, Eigen::VectorXd::Zero(1)).status, Status::InvalidEvaluation);
    problem = ScalarAffine(2.0);
    problem.residual = [](const Eigen::VectorXd&, Eigen::VectorXd& residual) {
        residual = Eigen::VectorXd::Constant(1, std::numeric_limits<double>::infinity());
        return true;
    };
    EXPECT_EQ(SolveEqualities(problem, Eigen::VectorXd::Zero(1)).status, Status::InvalidEvaluation);
    problem = ScalarAffine(2.0);
    problem.jacobian = [](const Eigen::VectorXd&, Eigen::MatrixXd&) { return false; };
    EXPECT_EQ(SolveEqualities(problem, Eigen::VectorXd::Zero(1)).status, Status::CallbackFailure);
    problem.jacobian = [](const Eigen::VectorXd&, Eigen::MatrixXd& jacobian) {
        jacobian = Eigen::MatrixXd::Zero(2, 1);
        return true;
    };
    EXPECT_EQ(SolveEqualities(problem, Eigen::VectorXd::Zero(1)).status, Status::InvalidEvaluation);
    problem.jacobian = [](const Eigen::VectorXd&, Eigen::MatrixXd&) { return true; };
    EXPECT_EQ(SolveEqualities(problem, Eigen::VectorXd::Zero(1)).status, Status::InvalidEvaluation);
}

TEST(EqualityConstraints, RejectPartiallyWrittenOutputsAfterASuccessfulEvaluation) {
    for (bool omit_jacobian_entry : {false, true}) {
        SCOPED_TRACE(omit_jacobian_entry);
        int residual_calls = 0;
        int jacobian_calls = 0;
        EqualityProblem problem;
        problem.variables = problem.equations = 2;
        problem.residual = [&](const Eigen::VectorXd& x, Eigen::VectorXd& residual) {
            residual(0) = x(0) * x(0) - 1.0;
            ++residual_calls;
            if (omit_jacobian_entry || residual_calls <= 2) {
                residual(1) = x(1) - 2.0;
            }
            return true;
        };
        problem.jacobian = [&](const Eigen::VectorXd& x, Eigen::MatrixXd& jacobian) {
            jacobian(0, 0) = 2.0 * x(0);
            jacobian(0, 1) = jacobian(1, 0) = 0.0;
            ++jacobian_calls;
            if (!omit_jacobian_entry || jacobian_calls == 1) {
                jacobian(1, 1) = 1.0;
            }
            return true;
        };
        SolveOptions options;
        options.damping = 0.0;
        const auto result = SolveEqualities(problem, Eigen::Vector2d(0.5, 0.0), options);
        EXPECT_EQ(result.status, Status::InvalidEvaluation);
        EXPECT_EQ(result.iterations, 1);
        EXPECT_EQ(result.x, Eigen::Vector2d(1.25, 2.0));
        EXPECT_EQ(result.residual_evaluations, omit_jacobian_entry ? 2 : 3);
        EXPECT_EQ(result.jacobian_evaluations, 2);
    }
}

TEST(EqualityConstraints, InvalidatePreviousRankWhenANewJacobianEvaluationFails) {
    for (bool nonfinite_output : {false, true}) {
        auto problem = ScalarAffine(1.0);
        problem.jacobian = [nonfinite_output](const Eigen::VectorXd& x, Eigen::MatrixXd& jacobian) {
            if (x(0) > 0.5) {
                if (nonfinite_output) {
                    jacobian = Eigen::MatrixXd::Constant(1, 1, std::numeric_limits<double>::infinity());
                    return true;
                }
                return false;
            }
            jacobian = Eigen::MatrixXd::Ones(1, 1);
            return true;
        };
        const auto result = SolveEqualities(problem, Eigen::VectorXd::Zero(1));
        EXPECT_EQ(result.status, nonfinite_output ? Status::InvalidEvaluation : Status::CallbackFailure);
        EXPECT_EQ(result.iterations, 1);
        EXPECT_GT(result.x(0), 0.5);
        EXPECT_EQ(result.rank, 0);
        EXPECT_FALSE(result.full_column_rank);
    }
}

TEST(EqualityConstraints, InvalidatePreviousRankWhenANewJacobianScalingOverflows) {
    EqualityProblem problem;
    problem.variables = problem.equations = 1;
    problem.residual = [](const Eigen::VectorXd& x, Eigen::VectorXd& residual) {
        residual = Eigen::VectorXd::Constant(1, (x(0) - 1.0) * 1e-308);
        return true;
    };
    problem.jacobian = [](const Eigen::VectorXd& x, Eigen::MatrixXd& jacobian) {
        jacobian = Eigen::MatrixXd::Constant(1, 1, x(0) > 0.5 ? 1e308 : 1e-308);
        return true;
    };
    SolveOptions options;
    options.residual_scales = Eigen::VectorXd::Constant(1, 1e-308);
    const auto result = SolveEqualities(problem, Eigen::VectorXd::Zero(1), options);
    EXPECT_EQ(result.status, Status::NumericalFailure);
    EXPECT_EQ(result.iterations, 1);
    EXPECT_GT(result.x(0), 0.5);
    EXPECT_EQ(result.rank, 0);
    EXPECT_FALSE(result.full_column_rank);
}

TEST(EqualityConstraints, ValidateDimensionsOptionsAndRepresentability) {
    auto problem = ScalarAffine(1.0);
    EXPECT_EQ(SolveEqualities(problem, Eigen::Vector2d::Zero()).status, Status::InvalidInput);
    SolveOptions options;
    options.variable_scales = Eigen::VectorXd::Zero(1);
    EXPECT_EQ(SolveEqualities(problem, Eigen::VectorXd::Zero(1), options).status, Status::InvalidInput);
    options = {};
    options.max_backtracking = 0;
    EXPECT_EQ(SolveEqualities(problem, Eigen::VectorXd::Zero(1), options).status, Status::InvalidInput);
    options = {};
    options.lower_bounds = Eigen::VectorXd::Constant(1, std::numeric_limits<double>::quiet_NaN());
    EXPECT_EQ(SolveEqualities(problem, Eigen::VectorXd::Zero(1), options).status, Status::InvalidInput);
    options = {};
    options.residual_scales = Eigen::VectorXd::Constant(1, 1e-308);
    EXPECT_EQ(SolveEqualities(ScalarAffine(1e308), Eigen::VectorXd::Zero(1), options).status, Status::NumericalFailure);
    problem.variables = 0;
    EXPECT_EQ(SolveEqualities(problem, Eigen::VectorXd{}).status, Status::InvalidInput);
}

TEST(EqualityConstraints, DistinguishIterationLimitAndZeroRank) {
    SolveOptions options;
    options.max_iterations = 0;
    EXPECT_EQ(SolveEqualities(ScalarAffine(1.0), Eigen::VectorXd::Zero(1), options).status, Status::MaxIterations);
    EXPECT_TRUE(SolveEqualities(ScalarAffine(1.0), Eigen::VectorXd::Ones(1), options).success());
    const auto zero_rank = AffineProblem(Eigen::MatrixXd::Zero(1, 1), Eigen::VectorXd::Ones(1));
    EXPECT_EQ(SolveEqualities(zero_rank, Eigen::VectorXd::Zero(1)).status, Status::RankDeficient);
}

TEST(EqualityConstraints, AcceptClosedZeroRankWarmStartWithoutClaimingRegularity) {
    EqualityProblem problem;
    problem.variables = problem.equations = 1;
    problem.residual = [](const Eigen::VectorXd& x, Eigen::VectorXd& residual) {
        residual = Eigen::VectorXd::Constant(1, x(0) * x(0));
        return true;
    };
    problem.jacobian = [](const Eigen::VectorXd& x, Eigen::MatrixXd& jacobian) {
        jacobian = Eigen::MatrixXd::Constant(1, 1, 2.0 * x(0));
        return true;
    };
    SolveOptions options;
    options.max_iterations = 0;
    const auto result = SolveEqualities(problem, Eigen::VectorXd::Zero(1), options);
    ASSERT_TRUE(result.success());
    EXPECT_EQ(result.iterations, 0);
    EXPECT_EQ(result.rank, 0);
    EXPECT_FALSE(result.full_column_rank);
    EXPECT_DOUBLE_EQ(result.scaled_residual_norm, 0.0);
}

TEST(ImplicitConstraints, DifferentiateAnalyticCircleAndRedundantConstraints) {
    const Eigen::MatrixXd a = Eigen::MatrixXd::Constant(1, 1, 1.6);
    const Eigen::MatrixXd b = Eigen::MatrixXd::Constant(1, 1, 1.2);
    const auto circle = DifferentiateImplicitly(a, b);
    ASSERT_TRUE(circle.success()) << circle.message;
    EXPECT_NEAR(circle.derivative(0, 0), -0.75, 1e-12);
    Eigen::MatrixXd redundant_a(2, 1);
    redundant_a << 1.0, 2.0;
    Eigen::MatrixXd redundant_b(2, 1);
    redundant_b << 3.0, 6.0;
    const auto redundant = DifferentiateImplicitly(redundant_a, redundant_b);
    ASSERT_TRUE(redundant.success()) << redundant.message;
    EXPECT_NEAR(redundant.derivative(0, 0), -3.0, 1e-12);
    EXPECT_EQ(redundant.rank, 1);
}

TEST(ImplicitConstraints, RejectNonuniqueAndInconsistentPartitions) {
    Eigen::Matrix2d a;
    a << 1.0, 1.0, 2.0, 2.0;
    const auto rank_deficient = DifferentiateImplicitly(a, Eigen::Vector2d(1.0, 2.0));
    EXPECT_EQ(rank_deficient.status, Status::RankDeficient);
    EXPECT_TRUE(rank_deficient.derivative.array().isNaN().all());
    EXPECT_EQ(DifferentiateImplicitly(Eigen::MatrixXd::Ones(1, 2), Eigen::MatrixXd::Ones(1, 1)).status,
              Status::RankDeficient);
    Eigen::MatrixXd full_rank(2, 1);
    full_rank << 1.0, 2.0;
    const auto inconsistent = DifferentiateImplicitly(full_rank, Eigen::Vector2d(3.0, 5.0));
    EXPECT_EQ(inconsistent.status, Status::Inconsistent);
    EXPECT_GT(inconsistent.relative_residual, 1e-3);
    EXPECT_TRUE(inconsistent.derivative.array().isNaN().all());
}

TEST(ImplicitConstraints, CheckEveryIndependentRightHandSide) {
    Eigen::MatrixXd a(2, 1);
    a << 1.0, 0.0;
    Eigen::Matrix2d b;
    b << 1e12, 0.0, 0.0, 1.0;
    const auto result = DifferentiateImplicitly(a, b);
    EXPECT_EQ(result.status, Status::Inconsistent);
    EXPECT_GE(result.relative_residual, 0.99);
}

TEST(ImplicitConstraints, KeepInconsistencyClassificationForTinyRightHandSides) {
    Eigen::MatrixXd a(2, 1);
    a << 1.0, 0.0;
    Eigen::MatrixXd b(2, 1);
    b << 0.0, 1e-12;
    for (double scale : {1e-15, 1.0, 1e15}) {
        SCOPED_TRACE(scale);
        const auto result = DifferentiateImplicitly(scale * a, scale * b);
        EXPECT_EQ(result.status, Status::Inconsistent);
        EXPECT_NEAR(result.relative_residual, 1.0, 1e-12);
        EXPECT_TRUE(result.derivative.array().isNaN().all());
    }
}

TEST(ImplicitConstraints, HandleZeroAndUnderflowedRightHandSidesExplicitly) {
    Eigen::MatrixXd a(2, 1);
    a << 1.0, 2.0;
    const auto zero = DifferentiateImplicitly(a, Eigen::MatrixXd::Zero(2, 2));
    ASSERT_TRUE(zero.success()) << zero.message;
    EXPECT_TRUE(zero.derivative.isZero(0.0));
    EXPECT_DOUBLE_EQ(zero.relative_residual, 0.0);

    DerivativeOptions options;
    options.equation_scales = Eigen::VectorXd::Constant(1, 1e300);
    const auto erased = DifferentiateImplicitly(Eigen::MatrixXd::Ones(1, 1),
                                               Eigen::MatrixXd::Constant(1, 1, 1e-300), options);
    EXPECT_EQ(erased.status, Status::NumericalFailure);
    EXPECT_TRUE(erased.derivative.array().isNaN().all());
}

TEST(ImplicitConstraints, RejectDerivativeDirectionsErasedByPublicUnitConversion) {
    DerivativeOptions options;
    options.variable_scales = Eigen::VectorXd::Constant(1, 1e-308);
    const auto result = DifferentiateImplicitly(Eigen::MatrixXd::Constant(1, 1, 1e308),
                                                Eigen::MatrixXd::Constant(1, 1, 1e-308), options);
    EXPECT_EQ(result.status, Status::NumericalFailure);
    EXPECT_EQ(result.rank, 1);
    EXPECT_TRUE(result.derivative.array().isNaN().all());
    EXPECT_NEAR(result.relative_residual, 1.0, 1e-12);
}

TEST(ImplicitConstraints, PreserveRepresentableSmallPublicDerivatives) {
    DerivativeOptions options;
    options.variable_scales = Eigen::VectorXd::Constant(1, 1e-308);
    const auto result = DifferentiateImplicitly(Eigen::MatrixXd::Constant(1, 1, 1e308),
                                                Eigen::MatrixXd::Ones(1, 1), options);
    ASSERT_TRUE(result.success()) << result.message;
    EXPECT_NE(result.derivative(0, 0), 0.0);
    EXPECT_LT(std::abs(result.derivative(0, 0) / -1e-308 - 1.0), 1e-14);
    EXPECT_LE(result.relative_residual, options.consistency_tolerance);
}

TEST(ImplicitConstraints, NormalizeRowsAndDependentCoordinateUnits) {
    Eigen::Matrix2d a = Eigen::Matrix2d::Zero();
    a.diagonal() << 1e12, 1e-12;
    const Eigen::Vector2d expected(2.0, 3.0);
    DerivativeOptions options;
    options.equation_scales = Eigen::Vector2d(1e12, 1e-12);
    const auto result = DifferentiateImplicitly(a, -a * expected, options);
    ASSERT_TRUE(result.success()) << result.message;
    EXPECT_LT((result.derivative.col(0) - expected).norm(), 1e-12);
    options.equation_scales(0) = 0.0;
    EXPECT_EQ(DifferentiateImplicitly(a, -a * expected, options).status, Status::InvalidInput);
    EXPECT_EQ(DifferentiateImplicitly(a, Eigen::MatrixXd::Zero(1, 1)).status, Status::InvalidInput);
    a(0, 0) = std::numeric_limits<double>::infinity();
    EXPECT_EQ(DifferentiateImplicitly(a, Eigen::Vector2d::Zero()).status, Status::InvalidInput);
}

TEST(ConstraintsScaling, PreserveRepresentableCombinedRowAndColumnScales) {
    for (int direction : {-1, 1}) {
        SCOPED_TRACE(direction);
        const Eigen::MatrixXd a = Eigen::MatrixXd::Constant(1, 1, std::ldexp(1.0, 900 * direction));
        const Eigen::VectorXd rhs = Eigen::VectorXd::Constant(1, std::ldexp(1.0, 400 * direction));
        const double row_scale = std::ldexp(1.0, -500 * direction);
        const double column_scale = std::ldexp(1.0, -700 * direction);
        const double expected = std::ldexp(1.0, -500 * direction);

        SolveOptions solve_options;
        solve_options.residual_scales = Eigen::VectorXd::Constant(1, row_scale);
        solve_options.variable_scales = Eigen::VectorXd::Constant(1, column_scale);
        solve_options.residual_tolerance = std::ldexp(1.0, -950);
        solve_options.damping = 0.0;
        const auto solved = SolveEqualities(AffineProblem(a, rhs), Eigen::VectorXd::Zero(1), solve_options);
        ASSERT_TRUE(solved.success()) << solved.message;
        EXPECT_EQ(solved.x(0), expected);
        EXPECT_TRUE(solved.full_column_rank);

        DerivativeOptions derivative_options;
        derivative_options.equation_scales = solve_options.residual_scales;
        derivative_options.variable_scales = solve_options.variable_scales;
        const auto differentiated = DifferentiateImplicitly(a, rhs, derivative_options);
        ASSERT_TRUE(differentiated.success()) << differentiated.message;
        EXPECT_EQ(differentiated.derivative(0, 0), -expected);
        EXPECT_EQ(differentiated.relative_residual, 0.0);
    }
}

}  // namespace
}  // namespace parallel_ankle::constraints
