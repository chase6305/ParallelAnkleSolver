#include <gtest/gtest.h>

#include <cmath>
#include <limits>
#include <stdexcept>

#include "fixtures.h"

namespace parallel_ankle::test {
namespace {

constexpr double kInfinity = std::numeric_limits<double>::infinity();
constexpr double kQuietNaN = std::numeric_limits<double>::quiet_NaN();
constexpr double kHalfPi = 1.57079632679489661923;

TEST(NumericalJacobianEdges, RejectUnrepresentablePerturbations) {
    const PMAnkle ankle(true, LegacyParams(true));
    for (double step : {1e-30, std::numeric_limits<double>::denorm_min()}) {
        EXPECT_THROW(ankle.ComputeNumericalJacobian(0.1, 0.1, step), std::invalid_argument);
        // The compatibility overload must retain the same validation.
        EXPECT_THROW(ankle.ComputeJacobian(0.1, 0.1, step), std::invalid_argument);
    }
    // One unchanged coordinate is enough to invalidate the stencil.
    EXPECT_THROW(ankle.ComputeNumericalJacobian(0.0, 0.1, 1e-30), std::invalid_argument);
}

TEST(NumericalJacobianEdges, RejectInvalidStepAndAngles) {
    const PMAnkle ankle(true, LegacyParams(true));
    for (double step : {0.0, -0.0, -1e-6, kInfinity, kQuietNaN}) {
        EXPECT_THROW(ankle.ComputeNumericalJacobian(0.1, 0.05, step), std::invalid_argument);
    }
    for (double bad_angle : {kInfinity, -kInfinity, kQuietNaN}) {
        EXPECT_THROW(ankle.ComputeNumericalJacobian(bad_angle, 0.0), std::invalid_argument);
        EXPECT_THROW(ankle.ComputeNumericalJacobian(0.0, bad_angle), std::invalid_argument);
    }
}

TEST(NumericalJacobianEdges, UseTheActualStepOnEachPerturbedAxis) {
    constexpr double kPeriod = 6.28318530717958647693;
    const PMAnkle ankle(true, LegacyParams(true));
    const Eigen::Vector2d pose(1024.0 * kPeriod + 0.2, -64.0 * kPeriod - 0.1);
    constexpr double eps = 1.7e-12;
    const auto base = ankle.GetAnkleMotorAngles(pose.x(), pose.y());
    ASSERT_TRUE(base.success());

    const auto jacobian = ankle.ComputeNumericalJacobian(pose.x(), pose.y(), eps);
    for (int axis = 0; axis < 2; ++axis) {
        Eigen::Vector2d perturbed = pose;
        perturbed(axis) += eps;
        const double step = perturbed(axis) - pose(axis);
        ASSERT_NE(step, eps);
        const auto shifted = ankle.GetAnkleMotorAngles(perturbed.x(), perturbed.y());
        ASSERT_TRUE(shifted.success());
        const Eigen::Vector2d expected = (shifted.motor_angles() - base.motor_angles()) / step;
        EXPECT_EQ(jacobian.col(axis), expected);
    }
}

TEST(InverseKinematicsEdges, PreserveIndependentActuatorFailures) {
    auto params = LegacyParams(true);
    params.l_rod1 = 1.0;
    const PMAnkle ankle(true, params);
    const auto result = ankle.GetAnkleMotorAngles(0.0, 0.0);
    EXPECT_FALSE(result.success());
    EXPECT_FALSE(result.success_1);
    EXPECT_TRUE(result.success_2);
    EXPECT_TRUE(std::isnan(result.theta_1_rad));
    EXPECT_TRUE(std::isfinite(result.theta_2_rad));
}

TEST(ForwardKinematicsEdges, RejectTargetsOutsideSelectedMotorBranchBeforeIterating) {
    const auto params = LegacyParams(true);
    const PMAnkle ankle(true, params);
    const auto neutral = ankle.GetAnkleMotorAngles(0.0, 0.0);
    ASSERT_TRUE(neutral.success());
    for (int motor = 0; motor < 2; ++motor) {
        for (double sign : {-1.0, 1.0}) {
            Eigen::Vector2d target = neutral.motor_angles();
            const double offset = motor == 0 ? params.alpha_1_rad : params.alpha_2_rad;
            target(motor) = offset + sign * (kHalfPi + 0.01);
            const auto result = ankle.SolveForwardKinematics(target(0), target(1));
            EXPECT_EQ(result.status, ForwardStatus::Unreachable);
            EXPECT_EQ(result.iterations, 0);

            // A loose convergence tolerance cannot make an unsupported target reachable.
            ForwardOptions loose;
            loose.tolerance_rad = 100.0;
            const auto loose_result = ankle.SolveForwardKinematics(target(0), target(1), loose);
            EXPECT_EQ(loose_result.status, ForwardStatus::Unreachable);
            EXPECT_EQ(loose_result.iterations, 0);
        }
    }
    for (double extreme : {std::numeric_limits<double>::max(), std::numeric_limits<double>::lowest()}) {
        const auto result = ankle.SolveForwardKinematics(extreme, neutral.theta_2_rad);
        EXPECT_EQ(result.status, ForwardStatus::Unreachable);
        EXPECT_EQ(result.iterations, 0);
    }
}

TEST(ForwardKinematicsEdges, DistinguishInvalidOptionsFromUnreachableFiniteTargets) {
    const PMAnkle ankle(true, LegacyParams(true));
    const auto neutral = ankle.GetAnkleMotorAngles(0.0, 0.0);
    ASSERT_TRUE(neutral.success());
    for (double bad_tolerance : {0.0, -1.0, kInfinity, kQuietNaN}) {
        ForwardOptions options;
        options.tolerance_rad = bad_tolerance;
        const auto result = ankle.SolveForwardKinematics(neutral.theta_1_rad, neutral.theta_2_rad, options);
        EXPECT_EQ(result.status, ForwardStatus::InvalidInput);
        EXPECT_EQ(result.iterations, 0);
    }
    for (double invalid : {kInfinity, -kInfinity, kQuietNaN}) {
        EXPECT_EQ(ankle.SolveForwardKinematics(invalid, neutral.theta_2_rad).status, ForwardStatus::InvalidInput);
        ForwardOptions options;
        options.initial_guess = {0.0, invalid};
        EXPECT_EQ(ankle.SolveForwardKinematics(neutral.theta_1_rad, neutral.theta_2_rad, options).status,
                  ForwardStatus::InvalidInput);
    }
}

TEST(ForwardKinematicsEdges, ReportSingularUnresolvedPoseWithoutReturningNaNState) {
    auto params = LegacyParams(true);
    params.r_a2 = params.r_a1;
    params.r_c2 = params.r_c1;
    params.l_bar2 = params.l_bar1;
    params.l_rod2 = params.l_rod1;
    params.alpha_2_rad = params.alpha_1_rad;
    params.is_left_motor2 = params.is_left_motor1;
    const PMAnkle ankle(true, params);
    const auto neutral = ankle.GetAnkleMotorAngles(0.0, 0.0);
    ASSERT_TRUE(neutral.success());
    const auto result = ankle.SolveForwardKinematics(neutral.theta_1_rad + 0.1, neutral.theta_2_rad - 0.1);
    EXPECT_EQ(result.status, ForwardStatus::Singular);
    EXPECT_EQ(result.iterations, 0);
    EXPECT_TRUE(result.pose.allFinite());
    EXPECT_TRUE(std::isfinite(result.residual_norm));
    EXPECT_GT(result.residual_norm, 0.0);
}

TEST(InverseKinematicsEdges, RejectEncoderOffsetsThatEraseThePhysicalAngle) {
    auto params = LegacyParams(true);
    params.alpha_1_rad = 1e308;
    const PMAnkle ankle(true, params);
    const auto result = ankle.GetAnkleMotorAngles(0.0, 0.0);
    EXPECT_FALSE(result.success_1);
    EXPECT_TRUE(result.success_2);
    EXPECT_EQ(result.status_1, InverseStatus::NumericalFailure);
    EXPECT_TRUE(std::isnan(result.theta_1_rad));
    const auto forward = ankle.SolveForwardKinematics(params.alpha_1_rad, result.theta_2_rad);
    EXPECT_EQ(forward.status, ForwardStatus::NumericalFailure);
    EXPECT_EQ(forward.iterations, 0);
}

TEST(InverseKinematicsEdges, EncoderRoundingCannotSilentlySelectAnotherClosingRoot) {
    AnkleParams params;
    params.r_c1 = params.r_c2 = {0.0, 1.0, 1.0};
    params.l_bar1 = params.l_bar2 = params.l_rod1 = params.l_rod2 = 1.0;
    params.is_left_motor1 = params.is_left_motor2 = true;
    params.alpha_1_rad = 1e308;
    const auto result = PMAnkle(true, params).GetAnkleMotorAngles(0.0, 0.0);
    // gamma=pi/2 and gamma=0 both close this rod, so closure alone cannot
    // detect that a huge offset rounded the selected root into the other one.
    EXPECT_FALSE(result.success_1);
    EXPECT_EQ(result.status_1, InverseStatus::NumericalFailure);
    EXPECT_TRUE(result.success_2);
    EXPECT_NEAR(result.theta_2_rad, kHalfPi, 1e-12);
}

AnkleParams EncoderEndpointFixture(bool identical_actuators, bool positive_motor_direction, double offset) {
    AnkleParams params;
    params.r_c1 = {0.0, 1.0, 1.0};
    params.r_a2 = identical_actuators ? Eigen::Vector3d(0.0, 0.0, 0.0) : Eigen::Vector3d(-1.0, 0.0, 1.0);
    params.r_c2 = identical_actuators ? params.r_c1 : Eigen::Vector3d(-1.0, 1.0, 0.0);
    params.l_bar1 = params.l_bar2 = params.l_rod1 = params.l_rod2 = 1.0;
    params.is_left_motor1 = positive_motor_direction;
    params.is_left_motor2 = true;
    params.alpha_1_rad = offset;
    return params;
}

TEST(ForwardKinematicsEdges, AcceptSuccessfulIkOutputsAtRepresentableEncoderEndpoints) {
    for (bool identical_actuators : {false, true}) {
        for (bool positive_motor_direction : {false, true}) {
            for (double offset : {0.0, -1000.0, 1000.0}) {
                SCOPED_TRACE(identical_actuators);
                SCOPED_TRACE(positive_motor_direction);
                SCOPED_TRACE(offset);
                const auto params = EncoderEndpointFixture(identical_actuators, positive_motor_direction, offset);
                const PMAnkle ankle(true, params);
                const auto ik = ankle.GetAnkleMotorAngles(0.0, 0.0);
                const double endpoint = offset + (positive_motor_direction ? kHalfPi : -kHalfPi);
                ASSERT_TRUE(ik.success_2);
                if (offset == 0.0) {
                    ASSERT_TRUE(ik.success_1) << ToString(ik.status_1);
                }
                if (!ik.success_1) {
                    // A strict rod-relative representability check may reject a
                    // rounded encoder output, but its endpoint is not out of range.
                    EXPECT_EQ(ik.status_1, InverseStatus::NumericalFailure);
                    EXPECT_EQ(ankle.SolveForwardKinematics(endpoint, ik.theta_2_rad).status,
                              ForwardStatus::NumericalFailure);
                    continue;
                }
                ASSERT_TRUE(ik.success());
                EXPECT_DOUBLE_EQ(ik.theta_1_rad, endpoint);
                if (!identical_actuators) {
                    const auto jacobian = ankle.ComputeJacobian(0.0, 0.0);
                    const double determinant = jacobian(0, 0) * jacobian(1, 1) -
                                               jacobian(0, 1) * jacobian(1, 0);
                    EXPECT_GT(std::abs(determinant), 0.5);
                }
                // A closed initial pose is a valid residual solution even if
                // identical actuator rows do not determine a unique local pose.
                const auto forward = ankle.SolveForwardKinematics(ik.theta_1_rad, ik.theta_2_rad);
                EXPECT_TRUE(forward.success()) << ToString(forward.status);
                EXPECT_EQ(forward.iterations, 0);
                EXPECT_EQ(forward.residual_norm, 0.0);
                EXPECT_TRUE(forward.pose.isZero());
            }
        }
    }
}

TEST(ForwardKinematicsEdges, RejectNextRepresentableValueBeyondEncoderEndpoints) {
    const auto params = EncoderEndpointFixture(false, true, 1000.0);
    const PMAnkle ankle(true, params);
    ForwardOptions options;
    options.tolerance_rad = 100.0;
    const double lower = params.alpha_1_rad - kHalfPi;
    const double upper = params.alpha_1_rad + kHalfPi;
    for (double outside : {std::nextafter(lower, -kInfinity), std::nextafter(upper, kInfinity)}) {
        const auto result = ankle.SolveForwardKinematics(outside, 0.0, options);
        EXPECT_EQ(result.status, ForwardStatus::Unreachable);
        EXPECT_EQ(result.iterations, 0);
    }
}

TEST(ForwardKinematicsEdges, PreserveResidualNormBelowTheSquaredDoubleRange) {
    auto params = LegacyParams(true);
    params.alpha_1_rad = params.alpha_2_rad = 0.0;
    const auto uncalibrated = PMAnkle(true, params).GetAnkleMotorAngles(0.0, 0.0);
    ASSERT_TRUE(uncalibrated.success());
    params.alpha_1_rad = -uncalibrated.theta_1_rad;
    params.alpha_2_rad = -uncalibrated.theta_2_rad;
    const PMAnkle ankle(true, params);
    const auto zero = ankle.GetAnkleMotorAngles(0.0, 0.0);
    ASSERT_TRUE(zero.success());
    ASSERT_EQ(zero.theta_1_rad, 0.0);
    ASSERT_EQ(zero.theta_2_rad, 0.0);
    ForwardOptions options;
    options.max_iterations = 0;
    options.tolerance_rad = 1e-250;
    const Eigen::Vector2d target(1e-200, -1e-200);
    const auto result = ankle.SolveForwardKinematics(target(0), target(1), options);
    EXPECT_EQ(result.status, ForwardStatus::MaxIterations);
    EXPECT_EQ(result.iterations, 0);
    EXPECT_GT(result.residual_norm, options.tolerance_rad);
    EXPECT_NEAR(result.residual_norm / std::hypot(target(0), target(1)), 1.0, 1e-14);
}

}  // namespace
}  // namespace parallel_ankle::test
