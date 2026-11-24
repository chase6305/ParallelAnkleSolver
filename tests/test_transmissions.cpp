#include <gtest/gtest.h>

#include <Eigen/Geometry>
#include <cmath>
#include <limits>
#include <stdexcept>

#include "parallel_ankle/transmissions.h"

namespace parallel_ankle::test {
namespace {

constexpr double kPi = 3.14159265358979323846;
constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();
constexpr double kInf = std::numeric_limits<double>::infinity();

void CheckRod(const RotaryRodParams& params, const Eigen::Vector3d& point, const RotaryRodSolution& solution) {
    using Extended = Eigen::Matrix<long double, 3, 1>;
    const Extended delta = point.cast<long double>() - params.pivot.cast<long double>();
    const long double gamma = solution.crank_phase_rad;
    const Extended crank = params.crank_length *
                           (params.u.cast<long double>() * std::cos(gamma) +
                            params.v.cast<long double>() * std::sin(gamma));
    const Extended rod = delta - crank;
    EXPECT_NEAR(static_cast<double>(std::hypot(rod.x(), rod.y(), rod.z()) / params.rod_length), 1.0, 2e-12);

    // Reconstruct the crank from the motor phase without subtracting a huge
    // absolute offset; this independently checks the calibrated rod geometry.
    const long double alpha = params.motor_offset_rad;
    const long double theta = solution.motor_phase_rad;
    const long double crank_cosine = std::cos(theta) * std::cos(alpha) + std::sin(theta) * std::sin(alpha);
    const long double crank_sine = params.motor_sign *
                                  (std::sin(theta) * std::cos(alpha) - std::cos(theta) * std::sin(alpha));
    const Extended calibrated_rod = delta - params.crank_length *
                                               (params.u.cast<long double>() * crank_cosine +
                                                params.v.cast<long double>() * crank_sine);
    EXPECT_NEAR(static_cast<double>(std::hypot(calibrated_rod.x(), calibrated_rod.y(), calibrated_rod.z()) /
                                    params.rod_length),
                1.0, 2e-12);
}

double PhaseForBranch(const RotaryRodResult& result, int branch) {
    EXPECT_TRUE(result.success());
    for (std::size_t index = 0; index < result.count; ++index) {
        if (result.solutions[index].branch == branch) {
            return result.solutions[index].motor_phase_rad;
        }
    }
    ADD_FAILURE() << "missing rotary-rod branch " << branch;
    return kNaN;
}

TEST(RotaryRod, EnumerateBothPhaseFamiliesAndCalibratedGradients) {
    RotaryRodParams params;
    params.motor_sign = -1;
    params.motor_offset_rad = 0.2;
    const Eigen::Vector3d point = Eigen::Vector3d::UnitX();
    const auto result = SolveRotaryRod(params, point);
    ASSERT_EQ(result.status, TransmissionStatus::Success);
    ASSERT_EQ(result.count, 2U);
    for (std::size_t index = 0; index < result.count; ++index) {
        const auto& solution = result.solutions[index];
        const double phase_sign = index == 0 ? 1.0 : -1.0;
        EXPECT_NEAR(solution.crank_phase_rad, phase_sign * kPi / 3.0, 1e-14);
        EXPECT_NEAR(solution.motor_phase_rad, 0.2 - phase_sign * kPi / 3.0, 1e-14);
        EXPECT_EQ(solution.branch, index == 0 ? -1 : 1);
        ASSERT_TRUE(solution.derivative_valid);
        EXPECT_NEAR(solution.crank_gradient.x(), -phase_sign / std::sqrt(3.0), 1e-13);
        EXPECT_NEAR(solution.crank_gradient.y(), 1.0, 1e-13);
        EXPECT_NEAR(solution.crank_gradient.z(), 0.0, 1e-14);
        EXPECT_TRUE(solution.motor_gradient.isApprox(-solution.crank_gradient, 1e-14));
        CheckRod(params, point, solution);
    }
}

TEST(RotaryRod, RetainFullCircleRootsBeyondPrincipalAsin) {
    RotaryRodParams params;
    params.rod_length = std::sqrt(3.0);
    const auto result = SolveRotaryRod(params, Eigen::Vector3d::UnitX());
    ASSERT_TRUE(result.success());
    ASSERT_EQ(result.count, 2U);
    EXPECT_NEAR(result.solutions[0].crank_phase_rad, 2.0 * kPi / 3.0, 1e-14);
    EXPECT_NEAR(result.solutions[1].crank_phase_rad, -2.0 * kPi / 3.0, 1e-14);
    for (std::size_t index = 0; index < result.count; ++index) {
        CheckRod(params, Eigen::Vector3d::UnitX(), result.solutions[index]);
    }
}

TEST(RotaryRod, PreserveBothBranchesAcrossThePhaseCut) {
    RotaryRodParams params;
    params.rod_length = 0.125;
    const Eigen::Vector3d point = -Eigen::Vector3d::UnitX();
    const double opening = 2.0 * std::asin(params.rod_length / 2.0);
    for (int motor_sign : {-1, 1}) {
        params.motor_sign = motor_sign;
        for (double offset : {-0.25, 0.0, 0.25}) {
            params.motor_offset_rad = offset;
            const auto result = SolveRotaryRod(params, point);
            ASSERT_EQ(result.status, TransmissionStatus::Success);
            ASSERT_EQ(result.count, 2U);
            EXPECT_NEAR(result.solutions[0].crank_phase_rad, -kPi + opening, 1e-14);
            EXPECT_NEAR(result.solutions[1].crank_phase_rad, kPi - opening, 1e-14);
            EXPECT_EQ(result.solutions[0].branch, -1);
            EXPECT_EQ(result.solutions[1].branch, 1);
            for (std::size_t index = 0; index < result.count; ++index) {
                const auto& solution = result.solutions[index];
                const double expected_motor = std::remainder(motor_sign * solution.crank_phase_rad + offset,
                                                             2.0 * kPi);
                EXPECT_NEAR(solution.motor_phase_rad, expected_motor, 1e-14);
                EXPECT_TRUE(solution.derivative_valid);
                CheckRod(params, point, solution);
            }
        }
    }

    const auto tangent = SolveRotaryRod(RotaryRodParams{}, {-2.0, 0.0, 0.0});
    ASSERT_EQ(tangent.status, TransmissionStatus::Tangent);
    ASSERT_EQ(tangent.count, 1U);
    EXPECT_EQ(tangent.solutions[0].crank_phase_rad, -kPi);
    EXPECT_EQ(tangent.solutions[0].motor_phase_rad, -kPi);
}

TEST(RotaryRod, SupportArbitraryFixedAxisAndOutOfPlanePoints) {
    RotaryRodParams params;
    params.pivot = {2.0, -1.0, 3.0};
    params.u = Eigen::Vector3d::UnitY();
    params.v = Eigen::Vector3d::UnitZ();
    const Eigen::Vector3d point = params.pivot + params.u;
    const auto result = SolveRotaryRod(params, point);
    ASSERT_TRUE(result.success());
    ASSERT_EQ(result.count, 2U);
    EXPECT_TRUE(result.solutions[0].crank_gradient.isApprox(Eigen::Vector3d(0.0, -1.0 / std::sqrt(3.0), 1.0), 1e-13));
    CheckRod(params, point, result.solutions[0]);

    params = RotaryRodParams{};
    params.rod_length = std::sqrt(2.0);
    const Eigen::Vector3d elevated(1.0, 0.0, 1.0);
    const auto spatial = SolveRotaryRod(params, elevated);
    ASSERT_TRUE(spatial.success());
    EXPECT_NEAR(spatial.solutions[0].crank_phase_rad, kPi / 3.0, 1e-14);
    EXPECT_NEAR(spatial.solutions[0].crank_gradient.z(), -2.0 / std::sqrt(3.0), 1e-13);
    CheckRod(params, elevated, spatial.solutions[0]);
}

TEST(RotaryRod, PreserveRigidFrameAndConsistentUnitChanges) {
    const RotaryRodParams reference;
    const Eigen::Vector3d point(1.0, 0.2, 0.3);
    const auto initial = SolveRotaryRod(reference, point);
    ASSERT_TRUE(initial.success());
    const Eigen::Matrix3d rotation = Eigen::AngleAxisd(0.7, Eigen::Vector3d(1.0, 2.0, -1.0).normalized()).toRotationMatrix();
    const Eigen::Vector3d translation(2.0, -1.0, 3.0);
    auto moved = reference;
    moved.pivot = translation;
    moved.u = rotation * moved.u;
    moved.v = rotation * moved.v;
    const auto transformed = SolveRotaryRod(moved, rotation * point + translation);
    ASSERT_TRUE(transformed.success());
    ASSERT_EQ(transformed.count, initial.count);
    for (std::size_t index = 0; index < initial.count; ++index) {
        EXPECT_NEAR(transformed.solutions[index].crank_phase_rad, initial.solutions[index].crank_phase_rad, 2e-14);
        EXPECT_TRUE(transformed.solutions[index].crank_gradient.isApprox(rotation * initial.solutions[index].crank_gradient,
                                                                       2e-13));
    }
    for (double scale : {1e-200, 1e200}) {
        auto scaled = reference;
        scaled.crank_length *= scale;
        scaled.rod_length *= scale;
        const auto result = SolveRotaryRod(scaled, point * scale);
        ASSERT_TRUE(result.success());
        ASSERT_EQ(result.count, initial.count);
        for (std::size_t index = 0; index < result.count; ++index) {
            EXPECT_NEAR(result.solutions[index].crank_phase_rad, initial.solutions[index].crank_phase_rad, 2e-14);
            EXPECT_TRUE((result.solutions[index].crank_gradient * scale).isApprox(initial.solutions[index].crank_gradient,
                                                                                 2e-13));
            CheckRod(scaled, point * scale, result.solutions[index]);
        }
    }
}

TEST(RotaryRod, MatchPointDifferencesWithinEachBranch) {
    RotaryRodParams params;
    params.motor_sign = -1;
    params.motor_offset_rad = 0.3;
    const Eigen::Vector3d point(1.0, 0.2, 0.3);
    const auto result = SolveRotaryRod(params, point);
    ASSERT_TRUE(result.success());
    for (std::size_t index = 0; index < result.count; ++index) {
        Eigen::Vector3d numerical;
        for (int coordinate = 0; coordinate < 3; ++coordinate) {
            Eigen::Vector3d plus = point;
            Eigen::Vector3d minus = point;
            plus(coordinate) += 1e-6;
            minus(coordinate) -= 1e-6;
            const double upper = PhaseForBranch(SolveRotaryRod(params, plus), result.solutions[index].branch);
            const double lower = PhaseForBranch(SolveRotaryRod(params, minus), result.solutions[index].branch);
            numerical(coordinate) = std::remainder(upper - lower, 2.0 * kPi) / 2e-6;
        }
        EXPECT_TRUE(numerical.isApprox(result.solutions[index].motor_gradient, 2e-8));
    }
}

TEST(RotaryRod, DistinguishTangencyContinuumAndUnreachableTargets) {
    RotaryRodParams params;
    const auto tangent = SolveRotaryRod(params, {2.0, 0.0, 0.0});
    ASSERT_EQ(tangent.status, TransmissionStatus::Tangent);
    ASSERT_EQ(tangent.count, 1U);
    EXPECT_TRUE(tangent.success());
    EXPECT_NEAR(tangent.solutions[0].crank_phase_rad, 0.0, 1e-14);
    EXPECT_FALSE(tangent.solutions[0].derivative_valid);
    EXPECT_FALSE(tangent.solutions[0].crank_gradient.allFinite());
    EXPECT_EQ(tangent.solutions[0].branch, 0);
    EXPECT_EQ(SolveRotaryRod(params, {2.001, 0.0, 0.0}).status, TransmissionStatus::Unreachable);
    EXPECT_EQ(SolveRotaryRod(params, {0.0, 0.0, 2.0}).status, TransmissionStatus::Unreachable);

    params.crank_length = 3.0;
    params.rod_length = 5.0;
    const auto continuum = SolveRotaryRod(params, {0.0, 0.0, 4.0});
    EXPECT_EQ(continuum.status, TransmissionStatus::Continuum);
    EXPECT_EQ(continuum.count, 0U);
    EXPECT_FALSE(continuum.success());
    params.rod_length = 4.0;
    EXPECT_EQ(SolveRotaryRod(params, {0.0, 0.0, 4.0}).status, TransmissionStatus::Unreachable);
}

TEST(RotaryRod, PreserveSmallOpeningsAndRejectUnrepresentableCalibratedPhases) {
    RotaryRodParams params;
    params.crank_length = std::ldexp(1.0, 40);
    const Eigen::Vector3d target(params.crank_length, 0.0, 0.0);
    const auto result = SolveRotaryRod(params, target);
    ASSERT_TRUE(result.success());
    ASSERT_EQ(result.count, 2U);
    const double expected = 2.0 * std::asin(0.5 / params.crank_length);
    EXPECT_NEAR(result.solutions[0].crank_phase_rad / expected, 1.0, 2e-14);
    EXPECT_NEAR(result.solutions[1].crank_phase_rad / expected, -1.0, 2e-14);
    CheckRod(params, target, result.solutions[0]);

    params.crank_length = std::ldexp(1.0, 60);
    params.motor_offset_rad = 1.0;
    const auto unrepresentable = SolveRotaryRod(params, {params.crank_length, 0.0, 0.0});
    EXPECT_EQ(unrepresentable.status, TransmissionStatus::NumericalFailure);
    EXPECT_EQ(unrepresentable.count, 0U);
}

TEST(RotaryRod, ReduceLargeAbsoluteOffsetsBeforeCombiningPhases) {
    RotaryRodParams params;
    params.motor_offset_rad = 1e300;
    const auto result = SolveRotaryRod(params, Eigen::Vector3d::UnitX());
    ASSERT_TRUE(result.success());
    ASSERT_EQ(result.count, 2U);
    for (std::size_t index = 0; index < result.count; ++index) {
        EXPECT_GE(result.solutions[index].motor_phase_rad, -kPi);
        EXPECT_LT(result.solutions[index].motor_phase_rad, kPi);
        CheckRod(params, Eigen::Vector3d::UnitX(), result.solutions[index]);
    }
}

TEST(RotaryRod, RejectOutsidePlaneReachDespiteLongCrank) {
    RotaryRodParams params;
    params.crank_length = 1e12;
    const auto outside = SolveRotaryRod(params, {1e12, 0.0, 1.0 + 1e-14});
    EXPECT_EQ(outside.status, TransmissionStatus::Unreachable);
    EXPECT_EQ(outside.count, 0U);

    // In this rotated plane the normal projection subtracts large products.
    // Its rounding interval overlaps the rod boundary, so it cannot certify
    // a tangent solution even if the computed projection rounds inward.
    const double base = std::ldexp(1.0, 40);
    params.u = {std::sqrt(0.5), std::sqrt(0.5), 0.0};
    params.v = Eigen::Vector3d::UnitZ();
    params.crank_length = std::sqrt(2.0) * (base + 1.0);
    params.rod_length = std::nextafter(std::sqrt(2.0), 0.0);
    const auto uncertain = SolveRotaryRod(params, {base, base + 2.0, 0.0});
    EXPECT_EQ(uncertain.status, TransmissionStatus::NumericalFailure);
    EXPECT_EQ(uncertain.count, 0U);
}

TEST(RotaryRod, ValidateGeometryAndRejectNonfiniteTargets) {
    for (double invalid : {0.0, -1.0, kInf, kNaN}) {
        RotaryRodParams params;
        params.crank_length = invalid;
        EXPECT_THROW(params.Validate(), std::invalid_argument);
        params = RotaryRodParams{};
        params.rod_length = invalid;
        EXPECT_THROW(params.Validate(), std::invalid_argument);
    }
    RotaryRodParams params;
    params.u = {2.0, 0.0, 0.0};
    EXPECT_THROW(params.Validate(), std::invalid_argument);
    params = RotaryRodParams{};
    params.v = params.u;
    EXPECT_THROW(params.Validate(), std::invalid_argument);
    params = RotaryRodParams{};
    params.motor_sign = 0;
    EXPECT_THROW(params.Validate(), std::invalid_argument);
    params = RotaryRodParams{};
    params.motor_offset_rad = kInf;
    EXPECT_THROW(params.Validate(), std::invalid_argument);
    params = RotaryRodParams{};
    params.pivot.z() = kNaN;
    EXPECT_THROW(params.Validate(), std::invalid_argument);
    for (double invalid : {kInf, kNaN}) {
        EXPECT_EQ(SolveRotaryRod(RotaryRodParams{}, {invalid, 0.0, 0.0}).status, TransmissionStatus::InvalidInput);
    }
}

TEST(PrismaticLength, ReturnLengthGradientAndStrokeFeasibility) {
    PrismaticLengthParams params;
    params.minimum_length = params.maximum_length = 5.0;
    const auto result = SolvePrismaticLength(params, {3.0, 4.0, 0.0});
    ASSERT_TRUE(result.success());
    EXPECT_DOUBLE_EQ(result.length, 5.0);
    ASSERT_TRUE(result.derivative_valid);
    EXPECT_TRUE(result.gradient.isApprox(Eigen::Vector3d(0.6, 0.8, 0.0), 1e-14));
    params.maximum_length = 6.0;
    params.minimum_length = 5.1;
    EXPECT_EQ(SolvePrismaticLength(params, {3.0, 4.0, 0.0}).status, TransmissionStatus::Unreachable);
    const auto zero = SolvePrismaticLength(PrismaticLengthParams{}, Eigen::Vector3d::Zero());
    EXPECT_EQ(zero.status, TransmissionStatus::Singular);
    EXPECT_TRUE(zero.success());
    EXPECT_DOUBLE_EQ(zero.length, 0.0);
    EXPECT_FALSE(zero.derivative_valid);
}

TEST(PrismaticLength, RejectInvalidAndUnrepresentableInputs) {
    PrismaticLengthParams params;
    params.maximum_length = -1.0;
    EXPECT_THROW(params.Validate(), std::invalid_argument);
    params = PrismaticLengthParams{};
    params.minimum_length = kNaN;
    EXPECT_THROW(params.Validate(), std::invalid_argument);
    EXPECT_EQ(SolvePrismaticLength(PrismaticLengthParams{}, {kNaN, 0.0, 0.0}).status, TransmissionStatus::InvalidInput);
    params = PrismaticLengthParams{};
    params.anchor = {-std::numeric_limits<double>::max(), 0.0, 0.0};
    EXPECT_EQ(SolvePrismaticLength(params, {std::numeric_limits<double>::max(), 0.0, 0.0}).status,
              TransmissionStatus::NumericalFailure);
}

}  // namespace
}  // namespace parallel_ankle::test
