#include <gtest/gtest.h>

#include <cmath>
#include <limits>
#include <stdexcept>

#include "fixtures.h"
#include "parallel_ankle/second_order.h"

namespace parallel_ankle::test {
namespace {

AnkleParams ScaleSecondOrderGeometry(AnkleParams params, double scale) {
    params.r_roll_center *= scale;
    params.r_pitch_center *= scale;
    params.r_a1 *= scale;
    params.r_a2 *= scale;
    params.r_c1 *= scale;
    params.r_c2 *= scale;
    params.l_bar1 *= scale;
    params.l_bar2 *= scale;
    params.l_rod1 *= scale;
    params.l_rod2 *= scale;
    return params;
}

AnkleParams IdentityJacobianParams() {
    AnkleParams params;
    params.r_c1 = {-1.0, 0.0, 0.0};
    params.r_a1 = {-1.0, -1.0, 1.0};
    params.r_c2 = {0.0, 1.0, 0.0};
    params.r_a2 = {0.0, 0.0, 1.0};
    params.l_bar1 = params.l_bar2 = params.l_rod1 = params.l_rod2 = 1.0;
    params.is_left_motor1 = params.is_left_motor2 = true;
    return params;
}

std::array<Eigen::Matrix2d, 2> DifferenceOfJacobians(const PMAnkle& model, const Eigen::Vector2d& pose, double step) {
    std::array<Eigen::Matrix2d, 2> result;
    for (int column = 0; column < 2; ++column) {
        Eigen::Vector2d plus = pose;
        Eigen::Vector2d minus = pose;
        plus(column) += step;
        minus(column) -= step;
        const Eigen::Matrix2d difference = (model.ComputeJacobian(plus(0), plus(1)) -
                                             model.ComputeJacobian(minus(0), minus(1))) /
                                            (2.0 * step);
        for (int motor = 0; motor < 2; ++motor) {
            result[motor].col(column) = difference.row(motor).transpose();
        }
    }
    return result;
}

double MotorPotential(const PMAnkle& model,
                      const Eigen::Vector2d& pose,
                      const Eigen::Vector2d& motor_reference,
                      const Eigen::Matrix2d& motor_stiffness) {
    const auto ik = model.GetAnkleMotorAngles(pose(0), pose(1));
    if (!ik.success()) {
        throw std::runtime_error("energy test left the supported branch");
    }
    const Eigen::Vector2d error = ik.motor_angles() - motor_reference;
    return 0.5 * error.dot(motor_stiffness * error);
}

Eigen::Matrix2d DifferenceOfEnergy(const PMAnkle& model,
                                  const Eigen::Vector2d& pose,
                                  const Eigen::Vector2d& motor_reference,
                                  const Eigen::Matrix2d& motor_stiffness,
                                  double step) {
    const double center = MotorPotential(model, pose, motor_reference, motor_stiffness);
    Eigen::Matrix2d result;
    for (int j = 0; j < 2; ++j) {
        const Eigen::Vector2d delta = step * Eigen::Vector2d::Unit(j);
        result(j, j) = (MotorPotential(model, pose + delta, motor_reference, motor_stiffness) - 2.0 * center +
                        MotorPotential(model, pose - delta, motor_reference, motor_stiffness)) /
                       (step * step);
    }
    const Eigen::Vector2d dp(step, 0.0);
    const Eigen::Vector2d dr(0.0, step);
    result(0, 1) = (MotorPotential(model, pose + dp + dr, motor_reference, motor_stiffness) -
                    MotorPotential(model, pose + dp - dr, motor_reference, motor_stiffness) -
                    MotorPotential(model, pose - dp + dr, motor_reference, motor_stiffness) +
                    MotorPotential(model, pose - dp - dr, motor_reference, motor_stiffness)) /
                   (4.0 * step * step);
    result(1, 0) = result(0, 1);
    return result;
}

TEST(SecondOrder, MatchIndependentCentralDifferencesOfJacobian) {
    const Eigen::Vector2d pose(0.2, -0.08);
    for (bool is_left : {true, false}) {
        auto params = LegacyParams(is_left);
        params.r_roll_center = {2.0, -3.0, -4.0};
        params.r_pitch_center = {-1.0, 4.0, 7.0};
        const PMAnkle model(is_left, params);
        const auto hessians = ComputeMotorAngleHessians(model, pose(0), pose(1));
        for (double step : {2e-5, 1e-5}) {
            const auto numerical = DifferenceOfJacobians(model, pose, step);
            for (int motor = 0; motor < 2; ++motor) {
                SCOPED_TRACE(motor);
                SCOPED_TRACE(step);
                EXPECT_LT((hessians[motor] - numerical[motor]).norm(), 1e-7);
                EXPECT_TRUE(hessians[motor].isApprox(hessians[motor].transpose(), 1e-14));
            }
        }
    }
}

TEST(SecondOrder, PreserveDerivativesUnderUniformGeometryScaling) {
    const Eigen::Vector2d pose(0.2, -0.08);
    const auto params = LegacyParams(true);
    const auto reference = ComputeMotorAngleHessians(PMAnkle(true, params), pose(0), pose(1));
    for (double scale : {1e-3, 1.0, 1e3}) {
        const auto hessians = ComputeMotorAngleHessians(PMAnkle(true, ScaleSecondOrderGeometry(params, scale)),
                                                       pose(0), pose(1));
        for (int motor = 0; motor < 2; ++motor) {
            EXPECT_LT((reference[motor] - hessians[motor]).norm(), 1e-11);
        }
    }
}

TEST(SecondOrder, KeepEncoderOffsetsOutOfNearFoldGeometricDerivatives) {
    AnkleParams params;
    params.r_c1 = params.r_c2 = Eigen::Vector3d(0.0, 3.0, -4.0);
    params.l_bar1 = params.l_bar2 = 2.0;
    params.l_rod1 = params.l_rod2 = std::nextafter(3.0, std::numeric_limits<double>::infinity());
    params.is_left_motor1 = params.is_left_motor2 = true;
    const PMAnkle original(true, params);
    ASSERT_TRUE(original.GetAnkleMotorAngles(0.0, 0.0).success());
    const auto original_jacobian = original.ComputeJacobian(0.0, 0.0);
    const auto original_hessians = ComputeMotorAngleHessians(original, 0.0, 0.0);

    params.alpha_1_rad = params.alpha_2_rad = 1e4;
    const PMAnkle shifted(true, params);
    ASSERT_TRUE(shifted.GetAnkleMotorAngles(0.0, 0.0).success());
    const auto shifted_hessians = ComputeMotorAngleHessians(shifted, 0.0, 0.0);
    EXPECT_TRUE(shifted.ComputeJacobian(0.0, 0.0).isApprox(original_jacobian, 1e-14));
    for (int motor = 0; motor < 2; ++motor) {
        EXPECT_TRUE(shifted_hessians[motor].isApprox(original_hessians[motor], 1e-13));
    }

    // Pitch enters this collinear fixture quadratically. A much smaller
    // symmetric pitch stencil remains on the same fold-adjacent branch;
    // normal finite-difference steps would leave its tiny feasible interval.
    constexpr double step = 1e-10;
    const auto plus = shifted.ComputeJacobian(step, 0.0);
    const auto minus = shifted.ComputeJacobian(-step, 0.0);
    for (int motor = 0; motor < 2; ++motor) {
        const double numerical_pp = (plus(motor, 0) - minus(motor, 0)) / (2.0 * step);
        const double analytic_pp = shifted_hessians[motor](0, 0);
        ASSERT_GT(std::abs(analytic_pp), 1e6);
        EXPECT_LT(std::abs(numerical_pp - analytic_pp) / std::abs(analytic_pp), 1e-10);
    }
}

TEST(SecondOrder, ShareTheStableShortRodOpeningWithTheCoreJacobian) {
    AnkleParams params;
    const double crank = std::ldexp(1.0, 40);
    params.r_c1 = params.r_c2 = Eigen::Vector3d(0.0, crank, 0.0);
    params.l_bar1 = params.l_bar2 = crank;
    params.l_rod1 = params.l_rod2 = 1.0;
    params.is_left_motor1 = params.is_left_motor2 = true;
    const PMAnkle model(true, params);
    const auto ik = model.GetAnkleMotorAngles(0.0, 0.0);
    ASSERT_TRUE(ik.success());
    EXPECT_GT(ik.theta_1_rad, 0.0);
    const auto jacobian = model.ComputeJacobian(0.0, 0.0);
    EXPECT_NEAR(jacobian(0, 1), 1.0, 1e-13);
    const auto hessians = ComputeMotorAngleHessians(model, 0.0, 0.0);
    for (const auto& hessian : hessians) {
        EXPECT_TRUE(hessian.allFinite());
        EXPECT_LT(hessian.norm(), 1e-12);
    }
}

TEST(SecondOrder, PreserveDerivativesUnderLargeRepresentableTranslations) {
    AnkleParams params;
    params.r_pitch_center = {0.0, 0.0, 2.0};
    params.r_a1 = {0.0, 0.0, 10.0};
    params.r_a2 = {0.0, 0.0, 8.0};
    params.r_c1 = {-2.0, 4.0, 0.0};
    params.r_c2 = {-2.0, -4.0, 0.0};
    params.l_bar1 = params.l_bar2 = 4.0;
    params.l_rod1 = std::sqrt(104.0);
    params.l_rod2 = std::sqrt(68.0);
    const auto reference = ComputeMotorAngleHessians(PMAnkle(true, params), 0.2, -0.1);
    const Eigen::Vector3d translation(1e16, -1e16, 1e16);
    params.r_roll_center += translation;
    params.r_pitch_center += translation;
    params.r_a1 += translation;
    params.r_a2 += translation;
    params.r_c1 += translation;
    params.r_c2 += translation;
    const auto translated = ComputeMotorAngleHessians(PMAnkle(true, params), 0.2, -0.1);
    for (int motor = 0; motor < 2; ++motor) {
        EXPECT_LT((reference[motor] - translated[motor]).norm(), 1e-11);
    }
}

TEST(LoadedStaticStiffness, ZeroPreloadReducesToExistingUnloadedMapping) {
    const PMAnkle model(true, LegacyParams(true));
    const Eigen::Matrix2d task = Eigen::Vector2d(100.0, 80.0).asDiagonal();
    const auto loaded = MapLoadedStaticTaskStiffness(model, 0.2, -0.08, task, Eigen::Vector2d::Zero());
    const auto unloaded = model.GetMotorStiffnessMatrixFromTaskStiffness(0.2, -0.08, 100.0, 80.0);
    EXPECT_TRUE(loaded.motor_stiffness.isApprox(unloaded, 1e-12));
    EXPECT_TRUE(loaded.geometric_stiffness.isZero(0.0));
    EXPECT_TRUE(loaded.task_positive_semidefinite);
    EXPECT_TRUE(loaded.motor_positive_semidefinite);
}

TEST(LoadedStaticStiffness, PreserveMixedScaleAndSubnormalSymmetricDiagonals) {
    const PMAnkle model(true, IdentityJacobianParams());
    ASSERT_TRUE(model.ComputeJacobian(0.0, 0.0).isIdentity(0.0));
    const double subnormal = std::numeric_limits<double>::denorm_min();
    for (const Eigen::Vector2d& diagonal : {Eigen::Vector2d(1e308, 1e-308), Eigen::Vector2d(1.0, subnormal)}) {
        const Eigen::Matrix2d task = diagonal.asDiagonal();
        const auto result = MapLoadedStaticTaskStiffness(model, 0.0, 0.0, task, Eigen::Vector2d::Zero());
        EXPECT_EQ(result.motor_stiffness(0, 0), diagonal(0));
        EXPECT_EQ(result.motor_stiffness(1, 1), diagonal(1));
        EXPECT_EQ(result.motor_stiffness(0, 1), 0.0);
        EXPECT_EQ(result.motor_stiffness(1, 0), 0.0);
        EXPECT_TRUE(result.geometric_stiffness.isZero(0.0));
        EXPECT_TRUE(result.task_positive_semidefinite);
        EXPECT_TRUE(result.motor_positive_semidefinite);
    }
}

TEST(LoadedStaticStiffness, PreserveSubnormalOffDiagonalAndTransposeInvariantMidpoint) {
    const PMAnkle model(true, IdentityJacobianParams());
    ASSERT_TRUE(model.ComputeJacobian(0.0, 0.0).isIdentity(0.0));
    const double subnormal = std::numeric_limits<double>::denorm_min();
    Eigen::Matrix2d task;
    task << 1.0, subnormal, subnormal, 1.0;
    const auto symmetric = MapLoadedStaticTaskStiffness(model, 0.0, 0.0, task, Eigen::Vector2d::Zero());
    EXPECT_EQ(symmetric.motor_stiffness(0, 1), subnormal);
    EXPECT_EQ(symmetric.motor_stiffness(1, 0), subnormal);

    task(1, 0) = 2.0 * subnormal;
    const Eigen::Matrix2d transposed = task.transpose();
    const auto first = MapLoadedStaticTaskStiffness(model, 0.0, 0.0, task, Eigen::Vector2d::Zero());
    const auto second = MapLoadedStaticTaskStiffness(model, 0.0, 0.0, transposed, Eigen::Vector2d::Zero());
    // The exact midpoint is 1.5*denorm_min, which rounds to 2*denorm_min.
    EXPECT_EQ(first.motor_stiffness(0, 1), 2.0 * subnormal);
    EXPECT_EQ(first.motor_stiffness(1, 0), 2.0 * subnormal);
    EXPECT_EQ(second.motor_stiffness(0, 1), first.motor_stiffness(0, 1));
    EXPECT_EQ(second.motor_stiffness(1, 0), first.motor_stiffness(1, 0));
}

TEST(LoadedStaticStiffness, PreserveRepresentableSubnormalAfterBothCongruenceFactors) {
    const double subnormal = std::numeric_limits<double>::denorm_min();
    Eigen::Matrix2d task;
    task << 1.0, subnormal, subnormal, 1.0;
    for (bool reverse : {false, true}) {
        SCOPED_TRACE(reverse);
        const double pitch_factor = reverse ? 1.0 / 1024.0 : 1024.0;
        const double roll_factor = 1.0 / pitch_factor;
        auto params = IdentityJacobianParams();
        params.r_c1.x() = params.r_a1.x() = -pitch_factor;
        params.r_c2.y() = roll_factor;
        params.r_a2.y() = roll_factor - 1.0;
        const PMAnkle model(true, params);
        const Eigen::Matrix2d jacobian = model.ComputeJacobian(0.0, 0.0);
        ASSERT_EQ(jacobian(0, 0), pitch_factor);
        ASSERT_EQ(jacobian(1, 1), roll_factor);
        ASSERT_EQ(jacobian(0, 1), 0.0);
        ASSERT_EQ(jacobian(1, 0), 0.0);

        const auto result = MapLoadedStaticTaskStiffness(model, 0.0, 0.0, task, Eigen::Vector2d::Zero());
        EXPECT_EQ(result.motor_stiffness(0, 0), roll_factor * roll_factor);
        EXPECT_EQ(result.motor_stiffness(1, 1), pitch_factor * pitch_factor);
        EXPECT_EQ(result.motor_stiffness(0, 1), subnormal);
        EXPECT_EQ(result.motor_stiffness(1, 0), subnormal);
        EXPECT_TRUE(result.task_positive_semidefinite);
        EXPECT_TRUE(result.motor_positive_semidefinite);

        const Eigen::Matrix2d overflowing = 1e308 * Eigen::Matrix2d::Identity();
        EXPECT_THROW(MapLoadedStaticTaskStiffness(model, 0.0, 0.0, overflowing, Eigen::Vector2d::Zero()),
                     std::overflow_error);
    }
}

TEST(LoadedStaticStiffness, AverageLargeOffDiagonalsWithoutOverflow) {
    const PMAnkle model(true, IdentityJacobianParams());
    ASSERT_TRUE(model.ComputeJacobian(0.0, 0.0).isIdentity(0.0));
    const double upper = 1.6e308;
    const double lower = std::nextafter(upper, 0.0);
    Eigen::Matrix2d task;
    task << 1e308, upper, lower, 1e308;
    const Eigen::Matrix2d transposed = task.transpose();
    const auto first = MapLoadedStaticTaskStiffness(model, 0.0, 0.0, task, Eigen::Vector2d::Zero());
    const auto second = MapLoadedStaticTaskStiffness(model, 0.0, 0.0, transposed, Eigen::Vector2d::Zero());
    ASSERT_TRUE(first.motor_stiffness.allFinite());
    EXPECT_EQ(first.motor_stiffness(0, 0), 1e308);
    EXPECT_EQ(first.motor_stiffness(1, 1), 1e308);
    EXPECT_GE(first.motor_stiffness(0, 1), lower);
    EXPECT_LE(first.motor_stiffness(0, 1), upper);
    EXPECT_EQ(first.motor_stiffness(1, 0), first.motor_stiffness(0, 1));
    EXPECT_EQ(second.motor_stiffness(0, 1), first.motor_stiffness(0, 1));
}

TEST(LoadedStaticStiffness, MatchComposedPotentialEnergyCurvatureAndRoundTrip) {
    const Eigen::Vector2d pose(0.15, -0.06);
    Eigen::Matrix2d motor_stiffness;
    motor_stiffness << 7.0, 1.2, 1.2, 4.0;
    for (bool is_left : {true, false}) {
        const PMAnkle model(is_left, LegacyParams(is_left));
        const auto ik = model.GetAnkleMotorAngles(pose(0), pose(1));
        ASSERT_TRUE(ik.success());
        const Eigen::Vector2d reference = ik.motor_angles() + Eigen::Vector2d(0.08, -0.05);
        const Eigen::Vector2d preload = motor_stiffness * (reference - ik.motor_angles());
        const auto hessians = ComputeMotorAngleHessians(model, pose(0), pose(1));
        const Eigen::Matrix2d jacobian = model.ComputeJacobian(pose(0), pose(1));
        const Eigen::Matrix2d geometric = -preload(0) * hessians[0] - preload(1) * hessians[1];
        const Eigen::Matrix2d task = jacobian.transpose() * motor_stiffness * jacobian + geometric;

        for (double step : {2e-4, 1e-4}) {
            const Eigen::Matrix2d energy_curvature = DifferenceOfEnergy(model, pose, reference, motor_stiffness, step);
            EXPECT_LT((energy_curvature - task).norm(), 5e-6) << "step=" << step;
        }
        const auto result = MapLoadedStaticTaskStiffness(model, pose(0), pose(1), task, preload);
        EXPECT_TRUE(result.motor_stiffness.isApprox(motor_stiffness, 1e-11));
        EXPECT_TRUE(result.geometric_stiffness.isApprox(geometric, 1e-12));
        EXPECT_TRUE((jacobian.transpose() * result.motor_stiffness * jacobian + result.geometric_stiffness)
                        .isApprox(task, 1e-11));
        EXPECT_TRUE(result.motor_positive_semidefinite);
    }
}

TEST(LoadedStaticStiffness, ReportIndefinitenessWithoutClampingItAway) {
    const PMAnkle model(true, LegacyParams(true));
    for (double scale : {1.0, 1e-30}) {
        const Eigen::Matrix2d task = (scale * Eigen::Vector2d(-1.0, 2.0)).asDiagonal();
        const auto result = MapLoadedStaticTaskStiffness(model, 0.2, -0.08, task, Eigen::Vector2d::Zero());
        EXPECT_FALSE(result.task_positive_semidefinite);
        EXPECT_FALSE(result.motor_positive_semidefinite);
    }
    const auto zero = MapLoadedStaticTaskStiffness(model, 0.2, -0.08, Eigen::Matrix2d::Zero(), Eigen::Vector2d::Zero());
    EXPECT_TRUE(zero.task_positive_semidefinite);
    EXPECT_TRUE(zero.motor_positive_semidefinite);
}

TEST(SecondOrder, RejectUndefinedDerivativesAndSingularLoadedMappings) {
    AnkleParams params;
    params.r_c1 = params.r_c2 = Eigen::Vector3d(0.0, 3.0, -4.0);
    params.l_bar1 = params.l_bar2 = 2.0;
    params.l_rod1 = params.l_rod2 = 3.0;
    params.is_left_motor1 = params.is_left_motor2 = true;
    const PMAnkle transmission_singular(true, params);
    ASSERT_TRUE(transmission_singular.GetAnkleMotorAngles(0.0, 0.0).success());
    EXPECT_THROW(ComputeMotorAngleHessians(transmission_singular, 0.0, 0.0), std::domain_error);

    params = LegacyParams(true);
    params.r_a2 = params.r_a1;
    params.r_c2 = params.r_c1;
    params.l_bar2 = params.l_bar1;
    params.l_rod2 = params.l_rod1;
    params.alpha_2_rad = params.alpha_1_rad;
    params.is_left_motor2 = params.is_left_motor1;
    const PMAnkle parallel_singular(true, params);
    EXPECT_THROW(MapLoadedStaticTaskStiffness(parallel_singular, 0.2, -0.08, Eigen::Matrix2d::Identity(),
                                             Eigen::Vector2d::Zero()), std::domain_error);
}

TEST(LoadedStaticStiffness, RejectNonfiniteAndAsymmetricInputs) {
    const PMAnkle model(true, LegacyParams(true));
    const double nan = std::numeric_limits<double>::quiet_NaN();
    EXPECT_THROW(ComputeMotorAngleHessians(model, nan, 0.0), std::invalid_argument);
    Eigen::Matrix2d asymmetric = Eigen::Matrix2d::Identity();
    asymmetric(0, 1) = 0.1;
    EXPECT_THROW(MapLoadedStaticTaskStiffness(model, 0.0, 0.0, asymmetric, Eigen::Vector2d::Zero()),
                 std::invalid_argument);
    Eigen::Matrix2d nonfinite = Eigen::Matrix2d::Identity();
    nonfinite(0, 0) = nan;
    EXPECT_THROW(MapLoadedStaticTaskStiffness(model, 0.0, 0.0, nonfinite, Eigen::Vector2d::Zero()),
                 std::invalid_argument);
    EXPECT_THROW(MapLoadedStaticTaskStiffness(model, 0.0, 0.0, Eigen::Matrix2d::Identity(), Eigen::Vector2d(nan, 0.0)),
                 std::invalid_argument);
}

}  // namespace
}  // namespace parallel_ankle::test
