#include <gtest/gtest.h>

#include <cmath>
#include <limits>
#include <stdexcept>

#include "fixtures.h"

namespace parallel_ankle::test {
namespace {

AnkleParams IdenticalActuators(const Eigen::Vector3d& foot_point, double bar, double rod) {
    AnkleParams params;
    params.r_c1 = params.r_c2 = foot_point;
    params.l_bar1 = params.l_bar2 = bar;
    params.l_rod1 = params.l_rod2 = rod;
    params.is_left_motor1 = params.is_left_motor2 = true;
    return params;
}

AnkleParams ScaleParams(AnkleParams params, double factor) {
    params.r_roll_center *= factor;
    params.r_pitch_center *= factor;
    params.r_a1 *= factor;
    params.r_a2 *= factor;
    params.r_c1 *= factor;
    params.r_c2 *= factor;
    params.l_bar1 *= factor;
    params.l_bar2 *= factor;
    params.l_rod1 *= factor;
    params.l_rod2 *= factor;
    return params;
}

Eigen::Matrix2d CentralDifference(const PMAnkle& ankle, const Eigen::Vector2d& pose) {
    constexpr double step = 2e-6;
    Eigen::Matrix2d difference;
    for (int column = 0; column < 2; ++column) {
        auto plus = pose;
        auto minus = pose;
        plus(column) += step;
        minus(column) -= step;
        const auto upper = ankle.GetAnkleMotorAngles(plus(0), plus(1));
        const auto lower = ankle.GetAnkleMotorAngles(minus(0), minus(1));
        if (!upper.success() || !lower.success()) {
            throw std::runtime_error("central-difference test fixture left its valid branch");
        }
        difference.col(column) = (upper.motor_angles() - lower.motor_angles()) / (2.0 * step);
    }
    return difference;
}

TEST(Geometry, RejectAsinFoldThatViolatesRodLength) {
    // The old expression gives sin(gamma)=0.9741657 but selected cos(gamma)<0.
    // Folding through asin produces a 2.12053 rod instead of sqrt(5.4).
    const PMAnkle ankle(true, IdenticalActuators({0.0, 1.0, -1.0}, 1.0, std::sqrt(5.4)));
    const auto solution = ankle.GetAnkleMotorAngles(0.0, 0.0);
    EXPECT_FALSE(solution.success_1);
    EXPECT_FALSE(solution.success_2);
    EXPECT_EQ(solution.status_1, InverseStatus::OutsideBranch);
    EXPECT_EQ(solution.status_2, InverseStatus::OutsideBranch);
    EXPECT_TRUE(std::isnan(solution.theta_1_rad));
    EXPECT_THROW(ankle.ComputeJacobian(0.0, 0.0), std::domain_error);
}

TEST(Geometry, RejectFoldedBranchForV02Geometry) {
    // Analytically constructed from actuator 1 at a 92-degree crank angle.
    // Its asin-folded 88-degree solution shortens the rod by 0.147415 mm.
    AnkleParams params;
    params.r_roll_center = {0.0, 0.0, -15.0};
    params.r_pitch_center = {0.0, 0.0, -5.0};
    params.r_a1 = {-29.61, 0.0, 261.0};
    params.r_c1 = {-29.61, 29.544, 0.0};
    params.l_bar1 = 30.0;
    params.l_rod1 = 266.209;
    params.alpha_1_rad = -0.1745;
    params.is_left_motor1 = true;
    params.r_a2 = {-29.61, 0.0, 181.0};
    params.r_c2 = {-29.61, -29.544, 0.0};
    params.l_bar2 = 30.0;
    params.l_rod2 = 186.209;
    params.alpha_2_rad = 0.1745;
    params.is_left_motor2 = false;
    const PMAnkle ankle(true, params);
    const auto solution = ankle.GetAnkleMotorAngles(0.5930940250083878, 0.5);
    EXPECT_FALSE(solution.success_1);
    EXPECT_EQ(solution.status_1, InverseStatus::OutsideBranch);
    EXPECT_TRUE(std::isnan(solution.theta_1_rad));
}

TEST(Geometry, KeepPhysicalCrankSideWhenFootCrossesMotorPlane) {
    // At zero, B=(0,1,0), C=(0,1,-1), and both unit links close.
    // At roll=-2, abs(C.y-A.y) would silently switch the crank to -Y.
    const PMAnkle ankle(true, IdenticalActuators({0.0, 1.0, -1.0}, 1.0, 1.0));
    ASSERT_TRUE(ankle.GetAnkleMotorAngles(0.0, 0.0).success());
    const auto crossed = ankle.GetAnkleMotorAngles(0.0, -2.0);
    EXPECT_FALSE(crossed.success());
    EXPECT_EQ(crossed.status_1, InverseStatus::OutsideBranch);
}

TEST(Geometry, PreserveAnglesAndJacobianUnderUniformLengthScaling) {
    const Eigen::Vector2d pose(0.2, -0.08);
    for (bool is_left : {true, false}) {
        const auto params = LegacyParams(is_left);
        const PMAnkle reference(is_left, params);
        const auto reference_ik = reference.GetAnkleMotorAngles(pose(0), pose(1));
        ASSERT_TRUE(reference_ik.success());
        const auto reference_jacobian = reference.ComputeJacobian(pose(0), pose(1));
        for (double factor : {1e-100, 1e-6, 1e-3, 1.0, 1e3, 1e6, 1e100}) {
            SCOPED_TRACE(factor);
            const PMAnkle scaled(is_left, ScaleParams(params, factor));
            const auto ik = scaled.GetAnkleMotorAngles(pose(0), pose(1));
            ASSERT_TRUE(ik.success());
            EXPECT_LT((ik.motor_angles() - reference_ik.motor_angles()).norm(), 5e-13);
            EXPECT_LT((scaled.ComputeJacobian(pose(0), pose(1)) - reference_jacobian).norm(), 1e-11);
            const auto fk = scaled.SolveForwardKinematics(ik.theta_1_rad, ik.theta_2_rad);
            ASSERT_TRUE(fk.success()) << ToString(fk.status);
            EXPECT_LT((fk.pose - pose).norm(), 2e-8);
        }
    }
}

TEST(Geometry, AnalyticJacobianMatchesIndependentCentralDifference) {
    const Eigen::Vector2d pose(0.17, -0.08);
    for (bool is_left : {true, false}) {
        auto params = LegacyParams(is_left);
        // Nonzero, noncoincident centers exercise all transform derivatives.
        params.r_roll_center = {2.0, -3.0, -4.0};
        params.r_pitch_center = {-1.0, 4.0, 7.0};
        const PMAnkle ankle(is_left, params);
        const auto ik = ankle.GetAnkleMotorAngles(pose(0), pose(1));
        ASSERT_TRUE(ik.success());
        ExpectRodClosure(params, pose, ik.motor_angles());
        const Eigen::Matrix2d analytic = ankle.ComputeJacobian(pose(0), pose(1));
        const Eigen::Matrix2d central = CentralDifference(ankle, pose);
        EXPECT_LT((analytic - central).norm(), 5e-8);
        EXPECT_LT((analytic - ankle.ComputeNumericalJacobian(pose(0), pose(1))).norm(), 1e-5);
        EXPECT_TRUE(ankle.ComputeJacobian(pose(0), pose(1), 2e-7)
                        .isApprox(ankle.ComputeNumericalJacobian(pose(0), pose(1), 2e-7), 1e-14));
    }
}

TEST(Geometry, RejectTransmissionSingularityWhileKeepingReachableIk) {
    // A, B, C are collinear: r=(0,3,-4), bar=2, rod=3.
    // All normalized closure quantities are exactly binary-representable.
    const PMAnkle ankle(true, IdenticalActuators({0.0, 3.0, -4.0}, 2.0, 3.0));
    const auto ik = ankle.GetAnkleMotorAngles(0.0, 0.0);
    ASSERT_TRUE(ik.success());
    ExpectRodClosure(ankle.parameters(), Eigen::Vector2d::Zero(), ik.motor_angles());
    EXPECT_THROW(ankle.ComputeJacobian(0.0, 0.0), std::domain_error);
    const auto closed_at_initial_guess = ankle.SolveForwardKinematics(ik.theta_1_rad, ik.theta_2_rad);
    EXPECT_TRUE(closed_at_initial_guess.success());
    EXPECT_EQ(closed_at_initial_guess.iterations, 0);
    EXPECT_EQ(ankle.SolveForwardKinematics(ik.theta_1_rad + 0.01, ik.theta_2_rad - 0.01).status,
              ForwardStatus::Singular);
}

TEST(Geometry, RejectNumericalStepsThatCannotPerturbTheInput) {
    const PMAnkle ankle(true, LegacyParams(true));
    EXPECT_THROW(ankle.ComputeNumericalJacobian(1e20, 0.0, 1e-6), std::invalid_argument);
    EXPECT_THROW(ankle.ComputeNumericalJacobian(0.0, 1e20, 1e-6), std::invalid_argument);
}

TEST(Geometry, PreserveRelativeGeometryUnderLargeRepresentableTranslations) {
    AnkleParams params;
    params.r_pitch_center = {0.0, 0.0, 2.0};
    params.r_a1 = {0.0, 0.0, 10.0};
    params.r_a2 = {0.0, 0.0, 8.0};
    params.r_c1 = {-2.0, 4.0, 0.0};
    params.r_c2 = {-2.0, -4.0, 0.0};
    params.l_bar1 = params.l_bar2 = 4.0;
    params.l_rod1 = std::sqrt(104.0);
    params.l_rod2 = std::sqrt(68.0);
    const PMAnkle original(true, params);
    const auto reference = original.GetAnkleMotorAngles(0.2, -0.1);
    ASSERT_TRUE(reference.success());
    const auto jacobian = original.ComputeJacobian(0.2, -0.1);

    // Every coordinate offset here is an exactly representable multiple of 2.
    // A world-space transform would still round away the rotated local motion.
    const Eigen::Vector3d translation(1e16, -1e16, 1e16);
    params.r_roll_center += translation;
    params.r_pitch_center += translation;
    params.r_a1 += translation;
    params.r_a2 += translation;
    params.r_c1 += translation;
    params.r_c2 += translation;
    const PMAnkle translated(true, params);
    const auto result = translated.GetAnkleMotorAngles(0.2, -0.1);
    ASSERT_TRUE(result.success());
    EXPECT_LT((result.motor_angles() - reference.motor_angles()).norm(), 1e-13);
    EXPECT_LT((translated.ComputeJacobian(0.2, -0.1) - jacobian).norm(), 1e-12);
}

TEST(Geometry, PreserveRepresentableShortRodOpeningWithLongCranks) {
    const double crank = std::ldexp(1.0, 40);
    const auto params = IdenticalActuators({0.0, crank, 0.0}, crank, 1.0);
    const PMAnkle model(true, params);
    const auto ik = model.GetAnkleMotorAngles(0.0, 0.0);
    ASSERT_TRUE(ik.success());
    const double expected = 2.0 * std::asin(0.5 / crank);
    for (double angle : {ik.theta_1_rad, ik.theta_2_rad}) {
        EXPECT_GT(angle, 0.0);
        EXPECT_NEAR(angle / expected, 1.0, 1e-13);
        const long double y = static_cast<long double>(crank) * (1.0L - std::cos(static_cast<long double>(angle)));
        const long double z = -static_cast<long double>(crank) * std::sin(static_cast<long double>(angle));
        EXPECT_NEAR(static_cast<double>(std::hypot(y, z)), 1.0, 2e-13);
    }
    const auto jacobian = model.ComputeJacobian(0.0, 0.0);
    EXPECT_TRUE(jacobian.allFinite());
    EXPECT_NEAR(jacobian(0, 1), 1.0, 1e-13);
}

TEST(Geometry, ValidateEncodedShortRodAndOutOfPlaneReachInRodUnits) {
    const double crank = std::ldexp(1.0, 40);
    auto params = IdenticalActuators({0.0, crank, 0.0}, crank, 1.0);
    // The physical opening is half an encoder ULP at this offset. An absolute
    // angle-error gate alone would allow it, but zero angle gives a zero rod.
    params.alpha_1_rad = 1e4;
    const auto encoded = PMAnkle(true, params).GetAnkleMotorAngles(0.0, 0.0);
    EXPECT_EQ(encoded.status_1, InverseStatus::NumericalFailure);
    EXPECT_FALSE(encoded.success_1);
    EXPECT_TRUE(encoded.success_2);
    params.alpha_1_rad = 0.0;
    params.r_c1.x() = std::nextafter(1.0, std::numeric_limits<double>::infinity());
    const auto outside = PMAnkle(true, params).GetAnkleMotorAngles(0.0, 0.0);
    EXPECT_EQ(outside.status_1, InverseStatus::Unreachable);
    EXPECT_TRUE(outside.success_2);
}

TEST(Geometry, PreserveSelectedShortRodRootWhenAnotherRootAlsoCloses) {
    const double crank = std::ldexp(1.0, 42);
    auto params = IdenticalActuators({0.0, crank, 1.0}, crank, 1.0);
    const auto selected = PMAnkle(true, params).GetAnkleMotorAngles(0.0, 0.0);
    ASSERT_TRUE(selected.success());
    const double expected = 2.0 * std::atan2(1.0, crank);
    EXPECT_GT(selected.theta_1_rad, 0.0);
    EXPECT_NEAR(selected.theta_1_rad / expected, 1.0, 1e-13);
    const long double phase = selected.theta_1_rad;
    const long double working = static_cast<long double>(crank) * std::sin(phase) - std::cos(phase);
    EXPECT_NEAR(static_cast<double>(working), 1.0, 2e-13);

    // Gamma=0 closes the same rod but belongs to the other working root.
    // At this offset the selected phase is erased, even though its absolute
    // angular loss is below the generic encoder-angle tolerance.
    params.alpha_1_rad = 1e4;
    const auto encoded = PMAnkle(true, params).GetAnkleMotorAngles(0.0, 0.0);
    EXPECT_FALSE(encoded.success_1);
    EXPECT_EQ(encoded.status_1, InverseStatus::NumericalFailure);
    EXPECT_TRUE(std::isnan(encoded.theta_1_rad));
    EXPECT_TRUE(encoded.success_2);
}

TEST(Geometry, MatchIndependentHighPrecisionReferenceNearBranchBoundary) {
    // 90-digit mpmath geometry + numerical differentiation of circle roots,
    // using the exact binary-double inputs (not the production Jacobian).
    // The default 1e-6 forward stencil has ~0.2% truncation error here.
    AnkleParams params;
    params.r_roll_center = {0.0, 0.0, -15.0};
    params.r_pitch_center = {0.0, 0.0, -5.0};
    params.r_a1 = {-29.61, 0.0, 261.0};
    params.r_c1 = {-29.61, -29.544, 0.0};
    params.l_bar1 = params.l_bar2 = 30.0;
    params.l_rod1 = 266.209;
    params.alpha_1_rad = 0.1745;
    params.is_left_motor1 = false;
    params.r_a2 = {-29.61, 0.0, 181.0};
    params.r_c2 = {-29.61, 29.544, 0.0};
    params.l_rod2 = 186.209;
    params.alpha_2_rad = -0.1745;
    params.is_left_motor2 = true;
    const PMAnkle model(false, params);
    Eigen::Matrix2d reference;
    reference << -35.686537329196420370473416548, 42.712179739524455253255576046,
                  1.28318990617409599985867974675, 0.707263279425711855370688046769;
    const auto actual = model.ComputeJacobian(-1.2, 0.18000000000000005);
    EXPECT_LT((actual - reference).norm() / reference.norm(), 1e-9);
}

}  // namespace
}  // namespace parallel_ankle::test
