#include <gtest/gtest.h>

#include <Eigen/LU>
#include <limits>
#include <stdexcept>

#include "fixtures.h"

namespace parallel_ankle::test {
namespace {

constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

TEST(Parameters, RejectInvalidGeometry) {
    EXPECT_THROW((PMAnkle(true, AnkleParams{})), std::invalid_argument);
    for (double bad_length : {0.0, -1.0, kNaN, std::numeric_limits<double>::infinity()}) {
        auto params = LegacyParams(true);
        params.l_rod2 = bad_length;
        EXPECT_THROW((PMAnkle(true, params)), std::invalid_argument);
    }
    auto params = LegacyParams(true);
    params.r_pitch_center(0) = kNaN;
    EXPECT_THROW(params.Validate(), std::invalid_argument);
    params = LegacyParams(true);
    params.alpha_1_rad = kNaN;
    EXPECT_THROW(params.Validate(), std::invalid_argument);
}

TEST(InverseKinematics, PreserveReferenceAnglesAndMirrorSymmetry) {
    // Reference angles for the fixture geometry.
    const PMAnkle right(false, LegacyParams(false));
    const auto right_ik = right.GetAnkleMotorAngles(0.3102, 0.0427);
    ASSERT_TRUE(right_ik.success());
    EXPECT_NEAR(right_ik.theta_1_rad, -0.6463969635890626, 1e-12);
    EXPECT_NEAR(right_ik.theta_2_rad, 0.5588796110904539, 1e-12);
    const PMAnkle left(true, LegacyParams(true));
    const auto left_ik = left.GetAnkleMotorAngles(0.3102, -0.0427);
    ASSERT_TRUE(left_ik.success());
    EXPECT_TRUE(left_ik.motor_angles().isApprox(-right_ik.motor_angles(), 1e-12));
}

TEST(InverseKinematics, DirectionComesFromMotorConfiguration) {
    auto params = LegacyParams(true);
    const auto original = PMAnkle(true, params).GetAnkleMotorAngles(0.2, 0.1);
    const auto relabeled = PMAnkle(false, params).GetAnkleMotorAngles(0.2, 0.1);
    ASSERT_TRUE(original.success());
    EXPECT_TRUE(original.motor_angles().isApprox(relabeled.motor_angles(), 1e-14));
    params.is_left_motor1 = !params.is_left_motor1;
    params.is_left_motor2 = !params.is_left_motor2;
    const auto reversed = PMAnkle(true, params).GetAnkleMotorAngles(0.2, 0.1);
    ASSERT_TRUE(reversed.success());
    EXPECT_NEAR(reversed.theta_1_rad, 2.0 * params.alpha_1_rad - original.theta_1_rad, 1e-12);
    EXPECT_NEAR(reversed.theta_2_rad, 2.0 * params.alpha_2_rad - original.theta_2_rad, 1e-12);
}

TEST(InverseKinematics, SatisfyPhysicalLinkLengths) {
    for (bool is_left : {true, false}) {
        const auto params = LegacyParams(is_left);
        const PMAnkle ankle(is_left, params);
        for (double pitch : {-0.2, 0.0, 0.3102}) {
            for (double roll : {-0.1, 0.0, 0.1}) {
                const auto ik = ankle.GetAnkleMotorAngles(pitch, roll);
                ASSERT_TRUE(ik.success());
                ExpectRodClosure(params, {pitch, roll}, ik.motor_angles());
            }
        }
    }
}

TEST(InverseKinematics, RejectNonfiniteAndUnreachableInputs) {
    const PMAnkle ankle(true, LegacyParams(true));
    const auto invalid = ankle.GetAnkleMotorAngles(kNaN, 0.0);
    EXPECT_FALSE(invalid.success_1);
    EXPECT_FALSE(invalid.success_2);
    EXPECT_TRUE(std::isnan(invalid.theta_1_rad));
    auto params = LegacyParams(true);
    params.l_rod1 = params.l_rod2 = 1.0;
    const PMAnkle unreachable(true, params);
    const auto ik = unreachable.GetAnkleMotorAngles(0.0, 0.0);
    EXPECT_FALSE(ik.success_1);
    EXPECT_FALSE(ik.success_2);
    EXPECT_THROW(unreachable.ComputeJacobian(0.0, 0.0), std::domain_error);
    EXPECT_EQ(unreachable.SolveForwardKinematics(0.0, 0.0).status, ForwardStatus::Unreachable);
}

TEST(ForwardKinematics, RecoverPosesForBothFeet) {
    for (bool is_left : {true, false}) {
        const PMAnkle ankle(is_left, LegacyParams(is_left));
        for (double pitch : {-0.2, 0.0, 0.3102}) {
            for (double roll : {-0.1, 0.0, 0.1}) {
                const auto ik = ankle.GetAnkleMotorAngles(pitch, roll);
                ASSERT_TRUE(ik.success());
                const auto fk = ankle.SolveForwardKinematics(ik.theta_1_rad, ik.theta_2_rad);
                ASSERT_TRUE(fk.success()) << ToString(fk.status);
                EXPECT_LT((fk.pose - Eigen::Vector2d(pitch, roll)).norm(), 2e-8);
                EXPECT_LE(fk.residual_norm, 1e-10);
            }
        }
    }
}

TEST(ForwardKinematics, ReuseTheAcceptedPoseAfterAnUnreachableNewtonStep) {
    const PMAnkle ankle(true, LegacyParams(true));
    const Eigen::Vector2d expected(0.2, -0.4);
    const auto target = ankle.GetAnkleMotorAngles(expected.x(), expected.y());
    ASSERT_TRUE(target.success());

    ForwardOptions options;
    options.initial_guess = {-0.8, 0.2};
    const auto initial = ankle.GetAnkleMotorAngles(options.initial_guess.x(), options.initial_guess.y());
    ASSERT_TRUE(initial.success());
    const Eigen::Matrix2d jacobian = ankle.ComputeJacobian(options.initial_guess.x(), options.initial_guess.y());
    const Eigen::Vector2d step = jacobian.fullPivLu().solve(initial.motor_angles() - target.motor_angles());
    const Eigen::Vector2d full_step = options.initial_guess - step;
    ASSERT_FALSE(ankle.GetAnkleMotorAngles(full_step.x(), full_step.y()).success());

    const auto result = ankle.SolveForwardKinematics(target.theta_1_rad, target.theta_2_rad, options);
    ASSERT_TRUE(result.success()) << ToString(result.status);
    EXPECT_GT(result.iterations, 1);
    EXPECT_LT((result.pose - expected).norm(), 2e-8);
    EXPECT_LE(result.residual_norm, options.tolerance_rad);
    ExpectRodClosure(ankle.parameters(), result.pose, target.motor_angles());
}

TEST(ForwardKinematics, ReportInvalidInputsAndIterationLimit) {
    const PMAnkle ankle(true, LegacyParams(true));
    EXPECT_EQ(ankle.SolveForwardKinematics(kNaN, 0.0).status, ForwardStatus::InvalidInput);
    EXPECT_THROW(ankle.GetAnklePoseFromMotorAngles(kNaN, 0.0), std::runtime_error);
    ForwardOptions options;
    options.max_iterations = 0;
    const auto ik = ankle.GetAnkleMotorAngles(0.2, 0.1);
    ASSERT_TRUE(ik.success());
    EXPECT_EQ(ankle.SolveForwardKinematics(ik.theta_1_rad, ik.theta_2_rad, options).status, ForwardStatus::MaxIterations);
    options.initial_guess = {0.2, 0.1};
    EXPECT_TRUE(ankle.SolveForwardKinematics(ik.theta_1_rad, ik.theta_2_rad, options).success());
    options.max_iterations = -1;
    EXPECT_EQ(ankle.SolveForwardKinematics(0.0, 0.0, options).status, ForwardStatus::InvalidInput);
    options.max_iterations = 30;
    options.tolerance_rad = kNaN;
    EXPECT_EQ(ankle.SolveForwardKinematics(0.0, 0.0, options).status, ForwardStatus::InvalidInput);
}

TEST(ForwardKinematics, DoNotReturnSuccessForImpossibleMotorAngles) {
    const PMAnkle ankle(true, LegacyParams(true));
    EXPECT_FALSE(ankle.SolveForwardKinematics(10.0, 10.0).success());
    EXPECT_THROW(ankle.GetAnklePoseFromMotorAngles(10.0, 10.0), std::runtime_error);
}

TEST(Jacobian, MatchCentralDifference) {
    for (bool is_left : {true, false}) {
        const PMAnkle ankle(is_left, LegacyParams(is_left));
        const Eigen::Vector2d q(0.2, -0.08);
        Eigen::Matrix2d central;
        for (int column = 0; column < 2; ++column) {
            Eigen::Vector2d plus = q;
            Eigen::Vector2d minus = q;
            plus(column) += 1e-6;
            minus(column) -= 1e-6;
            const auto upper = ankle.GetAnkleMotorAngles(plus(0), plus(1));
            const auto lower = ankle.GetAnkleMotorAngles(minus(0), minus(1));
            ASSERT_TRUE(upper.success() && lower.success());
            central.col(column) = (upper.motor_angles() - lower.motor_angles()) / 2e-6;
        }
        EXPECT_LT((ankle.ComputeJacobian(q(0), q(1)) - central).norm(), 1e-5);
        EXPECT_THROW(ankle.ComputeJacobian(0.0, 0.0, 0.0), std::invalid_argument);
        EXPECT_THROW(ankle.ComputeJacobian(kNaN, 0.0), std::invalid_argument);
    }
}

TEST(Mappings, ConservePowerAndRoundTripVelocitiesAndTorques) {
    for (bool is_left : {true, false}) {
        const PMAnkle ankle(is_left, LegacyParams(is_left));
        const Eigen::Matrix2d j = ankle.ComputeJacobian(0.2, -0.08);
        const Eigen::Vector2d task_velocity(0.7, -0.4);
        const Eigen::Vector2d motor_velocity = j * task_velocity;
        const Eigen::Vector2d motor_torque(10.0, -8.0);
        const auto task_torque = ankle.GetTaskTorquesFromMotorTorques(0.2, -0.08, 10.0, -8.0);
        EXPECT_NEAR(task_torque.dot(task_velocity), motor_torque.dot(motor_velocity), 1e-12);
        EXPECT_TRUE(ankle.GetMotorTorquesFromTaskTorques(0.2, -0.08, task_torque(0), task_torque(1))
                        .isApprox(motor_torque, 1e-12));
        EXPECT_TRUE(ankle.GetTaskVelocitiesFromMotorVelocities(0.2, -0.08, motor_velocity(0), motor_velocity(1))
                        .isApprox(task_velocity, 1e-12));
    }
}

TEST(Mappings, PreserveGainEnergyAndCoupling) {
    const PMAnkle ankle(true, LegacyParams(true));
    const Eigen::Matrix2d j = ankle.ComputeJacobian(0.2, -0.08);
    const Eigen::Matrix2d stiffness = ankle.GetMotorStiffnessMatrixFromTaskStiffness(0.2, -0.08, 100.0, 80.0);
    const Eigen::Matrix2d damping = ankle.GetMotorDampingMatrixFromTaskDamping(0.2, -0.08, 5.0, 4.0);
    const Eigen::Matrix2d task_stiffness = Eigen::Vector2d(100.0, 80.0).asDiagonal();
    const Eigen::Matrix2d task_damping = Eigen::Vector2d(5.0, 4.0).asDiagonal();
    EXPECT_TRUE((j.transpose() * stiffness * j).isApprox(task_stiffness, 1e-12));
    EXPECT_TRUE((j.transpose() * damping * j).isApprox(task_damping, 1e-12));
    EXPECT_TRUE(stiffness.isApprox(stiffness.transpose(), 1e-14));
    EXPECT_GT(stiffness.determinant(), 0.0);
    EXPECT_GT(std::abs(stiffness(0, 1)), 1e-3);
    EXPECT_THROW(ankle.GetMotorStiffnessMatrixFromTaskStiffness(0.0, 0.0, -1.0, 1.0), std::invalid_argument);
    EXPECT_THROW(ankle.GetTaskTorquesFromMotorTorques(0.0, 0.0, kNaN, 1.0), std::invalid_argument);
}

TEST(Mappings, PreserveSymmetricSubnormalGainAfterBothFactors) {
    AnkleParams params;
    params.r_c1 = {-4.0, -16.0, 0.0};
    params.r_a1 = {-4.0, -17.0, 1.0};
    params.r_c2 = {0.0, 1.0, 0.0};
    params.r_a2 = {0.0, 0.0, 1.0};
    params.l_bar1 = params.l_bar2 = params.l_rod1 = params.l_rod2 = 1.0;
    params.is_left_motor1 = params.is_left_motor2 = true;
    const PMAnkle ankle(true, params);
    Eigen::Matrix2d jacobian;
    jacobian << 4.0, -16.0, 0.0, 1.0;
    ASSERT_EQ(ankle.ComputeJacobian(0.0, 0.0), jacobian);

    // The first row of J^-1 is [1/4, 4]. Multiplying the gain by 1/4
    // first would erase it before the factor of 4 restores its magnitude.
    const double gain = std::numeric_limits<double>::denorm_min();
    Eigen::Matrix2d expected;
    expected << 0.0, gain, gain, 16.0 * gain;
    EXPECT_EQ(ankle.GetMotorStiffnessMatrixFromTaskStiffness(0.0, 0.0, gain, 0.0), expected);
    EXPECT_EQ(ankle.GetMotorDampingMatrixFromTaskDamping(0.0, 0.0, gain, 0.0), expected);
    EXPECT_THROW(ankle.GetMotorStiffnessMatrixFromTaskStiffness(0.0, 0.0, std::numeric_limits<double>::max(), 0.0),
                 std::overflow_error);
}

TEST(Mappings, RejectSingularInverseMappings) {
    auto params = LegacyParams(true);
    params.r_a2 = params.r_a1;
    params.r_c2 = params.r_c1;
    params.l_bar2 = params.l_bar1;
    params.l_rod2 = params.l_rod1;
    params.alpha_2_rad = params.alpha_1_rad;
    params.is_left_motor2 = params.is_left_motor1;
    const PMAnkle singular(true, params);
    EXPECT_THROW(singular.GetTaskVelocitiesFromMotorVelocities(0.0, 0.0, 1.0, 0.0), std::domain_error);
    EXPECT_THROW(singular.GetMotorTorquesFromTaskTorques(0.0, 0.0, 1.0, 0.0), std::domain_error);
    EXPECT_THROW(singular.GetMotorDampingMatrixFromTaskDamping(0.0, 0.0, 1.0, 1.0), std::domain_error);
    const auto neutral = singular.GetAnkleMotorAngles(0.0, 0.0);
    ASSERT_TRUE(neutral.success());
    EXPECT_EQ(singular.SolveForwardKinematics(neutral.theta_1_rad + 0.1, neutral.theta_2_rad - 0.1).status,
              ForwardStatus::Singular);
}

}  // namespace
}  // namespace parallel_ankle::test
