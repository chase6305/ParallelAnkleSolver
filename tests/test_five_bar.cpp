#include <gtest/gtest.h>

#include <Eigen/LU>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

#include "parallel_ankle/five_bar.h"

namespace parallel_ankle::test {
namespace {

constexpr double kPi = 3.14159265358979323846;
constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();
constexpr double kInf = std::numeric_limits<double>::infinity();

FiveBarParams AsymmetricGeometry() {
    FiveBarParams params;
    params.pivots = {Eigen::Vector2d(-0.4, 0.1), Eigen::Vector2d(0.5, -0.1)};
    params.proximal = {0.8, 0.7};
    params.distal = {1.0, 0.9};
    return params;
}

double Distance(const Eigen::Vector2d& first, const Eigen::Vector2d& second) {
    return std::hypot(first.x() - second.x(), first.y() - second.y());
}

void CheckRods(const FiveBarParams& params, const Eigen::Vector2d& position, const Eigen::Vector2d& angles) {
    for (int leg = 0; leg < 2; ++leg) {
        const Eigen::Vector2d crank = params.proximal[leg] * Eigen::Vector2d(std::cos(angles(leg)), std::sin(angles(leg)));
        const Eigen::Vector2d elbow = params.pivots[leg] + crank;
        EXPECT_NEAR(Distance(elbow, params.pivots[leg]) / params.proximal[leg], 1.0, 2e-12);
        EXPECT_NEAR(Distance(position, elbow) / params.distal[leg], 1.0, 2e-12);
    }
}

Eigen::Vector2d NearestAngles(const FiveBarInverseResult& result, const Eigen::Vector2d& reference) {
    EXPECT_TRUE(result.success());
    Eigen::Vector2d nearest = Eigen::Vector2d::Constant(kNaN);
    double best = kInf;
    for (std::size_t index = 0; index < result.count; ++index) {
        Eigen::Vector2d difference = result.solutions[index].motor_angles - reference;
        for (int leg = 0; leg < 2; ++leg) {
            difference(leg) = std::remainder(difference(leg), 2.0 * kPi);
        }
        if (difference.squaredNorm() < best) {
            best = difference.squaredNorm();
            nearest = reference + difference;
        }
    }
    return nearest;
}

TEST(FiveBarParameters, RejectNonpositiveAndNonfiniteParameters) {
    for (double invalid : {0.0, -1.0, kNaN, kInf}) {
        auto params = AsymmetricGeometry();
        params.proximal[0] = invalid;
        EXPECT_THROW((PlanarFiveBar(params)), std::invalid_argument);
        params = AsymmetricGeometry();
        params.distal[1] = invalid;
        EXPECT_THROW((PlanarFiveBar(params)), std::invalid_argument);
    }
    auto params = AsymmetricGeometry();
    params.pivots[0].x() = kNaN;
    EXPECT_THROW((PlanarFiveBar(params)), std::invalid_argument);
}

TEST(FiveBarGeometry, EnumerateAllAsymmetricWorkingAndAssemblyModes) {
    const auto params = AsymmetricGeometry();
    const PlanarFiveBar mechanism(params);
    const Eigen::Vector2d target(0.1, 1.0);
    const auto inverse = mechanism.InverseAll(target);
    ASSERT_TRUE(inverse.success());
    ASSERT_EQ(inverse.count, 4U);
    unsigned working_modes = 0;
    for (std::size_t index = 0; index < inverse.count; ++index) {
        const auto& solution = inverse.solutions[index];
        CheckRods(params, target, solution.motor_angles);
        ASSERT_NE(solution.branch.working[0], 0);
        ASSERT_NE(solution.branch.working[1], 0);
        const unsigned mode = (solution.branch.working[0] > 0 ? 2U : 0U) +
                              (solution.branch.working[1] > 0 ? 1U : 0U);
        working_modes |= 1U << mode;

        const auto forward = mechanism.ForwardAll(solution.motor_angles);
        ASSERT_TRUE(forward.success());
        ASSERT_EQ(forward.count, 2U);
        EXPECT_EQ(forward.solutions[0].branch.assembly, 1);
        EXPECT_EQ(forward.solutions[1].branch.assembly, -1);
        bool recovered = false;
        for (std::size_t assembly = 0; assembly < forward.count; ++assembly) {
            const auto& output = forward.solutions[assembly];
            CheckRods(params, output.position, solution.motor_angles);
            recovered = recovered || Distance(output.position, target) < 2e-12;
            const auto back = mechanism.InverseAll(output.position);
            ASSERT_TRUE(back.success());
            EXPECT_LT((NearestAngles(back, solution.motor_angles) - solution.motor_angles).norm(), 2e-12);
        }
        EXPECT_TRUE(recovered);
    }
    EXPECT_EQ(working_modes, 15U);
}

TEST(FiveBarGeometry, RemainConsistentAfterRigidBaseTransform) {
    const auto params = AsymmetricGeometry();
    const PlanarFiveBar original(params);
    const Eigen::Vector2d target(0.1, 1.0);
    const auto inverse = original.InverseAll(target);
    ASSERT_TRUE(inverse.success());
    const double angle = 0.7;
    Eigen::Matrix2d rotation;
    rotation << std::cos(angle), -std::sin(angle), std::sin(angle), std::cos(angle);
    const Eigen::Vector2d translation(2.0, -1.0);
    auto moved = params;
    for (auto& pivot : moved.pivots) {
        pivot = rotation * pivot + translation;
    }
    const PlanarFiveBar transformed(moved);
    const Eigen::Vector2d moved_target = rotation * target + translation;
    for (std::size_t index = 0; index < inverse.count; ++index) {
        const Eigen::Vector2d moved_angles = inverse.solutions[index].motor_angles + Eigen::Vector2d::Constant(angle);
        CheckRods(moved, moved_target, moved_angles);
        const auto before = original.Analyze(target, inverse.solutions[index].motor_angles);
        const auto after = transformed.Analyze(moved_target, moved_angles);
        ASSERT_TRUE(before.success());
        ASSERT_TRUE(after.success());
        EXPECT_TRUE((after.motor_jacobian * rotation).isApprox(before.motor_jacobian, 2e-12));
    }
}

TEST(FiveBarGeometry, SupportCoaxialBasesWithFiniteForwardSolutions) {
    FiveBarParams params;
    params.pivots = {Eigen::Vector2d::Zero(), Eigen::Vector2d::Zero()};
    params.proximal = {1.0, 1.0};
    params.distal = {1.5, 1.5};
    const PlanarFiveBar mechanism(params);
    const Eigen::Vector2d motors(kPi / 3.0, -kPi / 3.0);
    const auto forward = mechanism.ForwardAll(motors);
    ASSERT_TRUE(forward.success());
    ASSERT_EQ(forward.count, 2U);
    for (std::size_t index = 0; index < forward.count; ++index) {
        CheckRods(params, forward.solutions[index].position, motors);
        EXPECT_NEAR(forward.solutions[index].position.y(), 0.0, 1e-14);
    }
    // The positive, symmetric extension agrees with the virtual-leg geometry.
    const double expected_length = std::cos(kPi / 3.0) +
                                   std::sqrt(1.5 * 1.5 - std::sin(kPi / 3.0) * std::sin(kPi / 3.0));
    const double maximum_x = std::max(forward.solutions[0].position.x(), forward.solutions[1].position.x());
    EXPECT_NEAR(maximum_x, expected_length, 1e-14);
}

TEST(FiveBarGeometry, DeduplicateTangenciesAndDiagnoseInfiniteSolutionSets) {
    FiveBarParams params;
    params.pivots = {Eigen::Vector2d(0.0, 0.0), Eigen::Vector2d(2.0, 0.0)};
    const auto tangent = PlanarFiveBar(params).ForwardAll(Eigen::Vector2d::Zero());
    ASSERT_TRUE(tangent.success());
    ASSERT_EQ(tangent.count, 1U);
    EXPECT_NEAR(Distance(tangent.solutions[0].position, {2.0, 0.0}), 0.0, 1e-14);
    EXPECT_EQ(tangent.solutions[0].branch.assembly, 0);

    params.pivots[1].setZero();
    const PlanarFiveBar coaxial(params);
    const auto infinite = coaxial.ForwardAll(Eigen::Vector2d::Zero());
    EXPECT_EQ(infinite.status, FiveBarStatus::Singular);
    EXPECT_EQ(infinite.count, 0U);
    EXPECT_EQ(coaxial.InverseAll(Eigen::Vector2d::Zero()).status, FiveBarStatus::Singular);
    const auto stretched = coaxial.InverseAll({2.0, 0.0});
    ASSERT_TRUE(stretched.success());
    ASSERT_EQ(stretched.count, 1U);
    EXPECT_EQ(stretched.solutions[0].branch.working[0], 0);
    EXPECT_EQ(stretched.solutions[0].branch.working[1], 0);

    params.distal[1] = 2.0;
    EXPECT_EQ(PlanarFiveBar(params).ForwardAll(Eigen::Vector2d::Zero()).status, FiveBarStatus::Unreachable);
    params.proximal = {2.0, 2.0};
    params.distal = {1.0, 1.0};
    const auto internal_tangent = PlanarFiveBar(params).InverseAll({1.0, 0.0});
    ASSERT_TRUE(internal_tangent.success());
    EXPECT_EQ(internal_tangent.count, 1U);
    EXPECT_EQ(PlanarFiveBar(params).InverseAll(Eigen::Vector2d::Zero()).status, FiveBarStatus::Unreachable);
    params = FiveBarParams{};
    params.pivots = {Eigen::Vector2d(0.0, 0.0), Eigen::Vector2d(10.0, 0.0)};
    EXPECT_EQ(PlanarFiveBar(params).InverseAll(Eigen::Vector2d::Zero()).status, FiveBarStatus::Unreachable);
}

TEST(FiveBarGeometry, RejectTrulyUnreachableTargetsWithoutBroadClamping) {
    FiveBarParams params;
    params.pivots = {Eigen::Vector2d::Zero(), Eigen::Vector2d::Zero()};
    const double outside = std::nextafter(2.0, kInf);
    EXPECT_EQ(PlanarFiveBar(params).InverseAll({outside, 0.0}).status, FiveBarStatus::Unreachable);
    params.pivots[1] = {outside, 0.0};
    EXPECT_EQ(PlanarFiveBar(params).ForwardAll(Eigen::Vector2d::Zero()).status, FiveBarStatus::Unreachable);
}

TEST(FiveBarDifferential, MatchIndependentDifferencesOnEveryWorkingMode) {
    const PlanarFiveBar mechanism(AsymmetricGeometry());
    const Eigen::Vector2d target(0.1, 1.0);
    const auto inverse = mechanism.InverseAll(target);
    ASSERT_TRUE(inverse.success());
    for (std::size_t index = 0; index < inverse.count; ++index) {
        const Eigen::Vector2d motors = inverse.solutions[index].motor_angles;
        const auto differential = mechanism.Analyze(target, motors);
        ASSERT_TRUE(differential.success());
        Eigen::Matrix2d numerical;
        for (int coordinate = 0; coordinate < 2; ++coordinate) {
            Eigen::Vector2d plus = target;
            Eigen::Vector2d minus = target;
            plus(coordinate) += 1e-6;
            minus(coordinate) -= 1e-6;
            const auto upper = mechanism.InverseAll(plus);
            const auto lower = mechanism.InverseAll(minus);
            ASSERT_TRUE(upper.success());
            ASSERT_TRUE(lower.success());
            numerical.col(coordinate) = (NearestAngles(upper, motors) - NearestAngles(lower, motors)) / 2e-6;
        }
        EXPECT_LT((numerical - differential.motor_jacobian).norm(), 2e-8);
        EXPECT_TRUE((differential.B * differential.motor_jacobian).isApprox(differential.A, 1e-13));
        const Eigen::Vector2d task_velocity(0.3, -0.2);
        const Eigen::Vector2d motor_torque(2.0, -3.0);
        const Eigen::Vector2d motor_velocity = differential.motor_jacobian * task_velocity;
        const Eigen::Vector2d task_force = differential.motor_jacobian.transpose() * motor_torque;
        EXPECT_NEAR(task_force.dot(task_velocity), motor_torque.dot(motor_velocity), 1e-13);
    }
}

TEST(FiveBarDifferential, SeparateActuatorAndAssemblySingularities) {
    FiveBarParams assembly_params;
    assembly_params.pivots = {Eigen::Vector2d(0.0, -1.0), Eigen::Vector2d(2.0, -1.0)};
    const PlanarFiveBar assembly(assembly_params);
    const Eigen::Vector2d vertical(kPi / 2.0, kPi / 2.0);
    const auto forward = assembly.ForwardAll(vertical);
    ASSERT_TRUE(forward.success());
    ASSERT_EQ(forward.count, 1U);
    const auto assembly_result = assembly.Analyze(forward.solutions[0].position, vertical);
    EXPECT_EQ(assembly_result.status, FiveBarStatus::Singular);
    EXPECT_TRUE(assembly_result.assembly_singular);
    EXPECT_FALSE(assembly_result.actuator_singular);
    EXPECT_TRUE(assembly_result.motor_jacobian.allFinite());

    FiveBarParams actuator_params;
    actuator_params.pivots = {Eigen::Vector2d::Zero(), Eigen::Vector2d(2.0, -1.0)};
    actuator_params.proximal = {1.0, 0.5};
    const PlanarFiveBar actuator(actuator_params);
    const Eigen::Vector2d target(2.0, 0.0);
    const auto inverse = actuator.InverseAll(target);
    ASSERT_TRUE(inverse.success());
    ASSERT_EQ(inverse.count, 2U);
    for (std::size_t index = 0; index < inverse.count; ++index) {
        const auto result = actuator.Analyze(target, inverse.solutions[index].motor_angles);
        EXPECT_EQ(result.status, FiveBarStatus::Singular);
        EXPECT_TRUE(result.actuator_singular);
        EXPECT_FALSE(result.assembly_singular);
        EXPECT_FALSE(result.motor_jacobian.allFinite());
    }
}

TEST(FiveBarEdges, HandleInvalidAndExtremeFiniteInputs) {
    const PlanarFiveBar mechanism(AsymmetricGeometry());
    for (double invalid : {kNaN, kInf, -kInf}) {
        EXPECT_EQ(mechanism.ForwardAll({invalid, 0.0}).status, FiveBarStatus::InvalidInput);
        EXPECT_EQ(mechanism.InverseAll({0.0, invalid}).status, FiveBarStatus::InvalidInput);
        EXPECT_EQ(mechanism.Analyze({0.0, 0.0}, {invalid, 0.0}).status, FiveBarStatus::InvalidInput);
    }
    for (double invalid_tolerance : {0.0, -1.0, 1.0, kInf, kNaN}) {
        EXPECT_EQ(mechanism.Analyze({0.1, 1.0}, {0.0, 0.0}, invalid_tolerance).status, FiveBarStatus::InvalidInput);
    }
    EXPECT_EQ(mechanism.Analyze({100.0, 100.0}, {0.0, 0.0}).status, FiveBarStatus::Unreachable);
    auto params = AsymmetricGeometry();
    params.pivots = {Eigen::Vector2d(-std::numeric_limits<double>::max(), 0.0),
                     Eigen::Vector2d(std::numeric_limits<double>::max(), 0.0)};
    const PlanarFiveBar separated(params);
    EXPECT_EQ(separated.ForwardAll(Eigen::Vector2d::Zero()).status, FiveBarStatus::Unreachable);
    EXPECT_EQ(separated.InverseAll(Eigen::Vector2d::Zero()).status, FiveBarStatus::Unreachable);
}

TEST(FiveBarEdges, PreserveGeometryAtVerySmallAndLargeConsistentScales) {
    const auto original = AsymmetricGeometry();
    const Eigen::Vector2d target(0.1, 1.0);
    const PlanarFiveBar reference(original);
    const auto reference_inverse = reference.InverseAll(target);
    ASSERT_TRUE(reference_inverse.success());
    for (double scale : {1e-200, 1e200}) {
        auto params = original;
        for (int leg = 0; leg < 2; ++leg) {
            params.pivots[leg] *= scale;
            params.proximal[leg] *= scale;
            params.distal[leg] *= scale;
        }
        const PlanarFiveBar scaled(params);
        const auto inverse = scaled.InverseAll(target * scale);
        ASSERT_TRUE(inverse.success());
        ASSERT_EQ(inverse.count, reference_inverse.count);
        for (std::size_t index = 0; index < inverse.count; ++index) {
            EXPECT_LT((inverse.solutions[index].motor_angles - reference_inverse.solutions[index].motor_angles).norm(),
                      2e-12);
            CheckRods(params, target * scale, inverse.solutions[index].motor_angles);
            const auto differential = scaled.Analyze(target * scale, inverse.solutions[index].motor_angles);
            const auto reference_differential = reference.Analyze(target, reference_inverse.solutions[index].motor_angles);
            ASSERT_TRUE(differential.success());
            EXPECT_TRUE((differential.motor_jacobian * scale).isApprox(reference_differential.motor_jacobian, 2e-12));
            const auto forward = scaled.ForwardAll(inverse.solutions[index].motor_angles);
            ASSERT_TRUE(forward.success());
            bool recovered = false;
            for (std::size_t assembly = 0; assembly < forward.count; ++assembly) {
                recovered = recovered || (forward.solutions[assembly].position / scale - target).norm() < 2e-12;
            }
            EXPECT_TRUE(recovered);
        }
    }
}

TEST(FiveBarRepresentability, RejectNarrowedPointsThatBreakShortRods) {
    FiveBarParams params;
    params.pivots = {Eigen::Vector2d(1e16, 0.0), Eigen::Vector2d(1e16 + 2.0, 0.0)};
    const auto result = PlanarFiveBar(params).ForwardAll({kPi / 2.0, kPi / 2.0});
    // The tangent has x=1e16+1, which rounds onto a pivot in a double.
    EXPECT_EQ(result.status, FiveBarStatus::NumericalFailure);
    EXPECT_EQ(result.count, 0U);
}

TEST(FiveBarRepresentability, ResolveTinyOpeningBesideLongEqualSides) {
    const double length = std::ldexp(1.0, 40);
    FiveBarParams params;
    params.pivots = {Eigen::Vector2d::Zero(), Eigen::Vector2d::Zero()};
    params.proximal = {length, length};
    const Eigen::Vector2d target(length, 0.0);
    const auto result = PlanarFiveBar(params).InverseAll(target);
    ASSERT_TRUE(result.success());
    ASSERT_EQ(result.count, 4U);
    const double opening = 2.0 * std::asin(1.0 / (2.0 * length));
    for (std::size_t index = 0; index < result.count; ++index) {
        CheckRods(params, target, result.solutions[index].motor_angles);
        for (int leg = 0; leg < 2; ++leg) {
            EXPECT_NEAR(std::abs(result.solutions[index].motor_angles(leg)) / opening, 1.0, 2e-15);
        }
    }
}

TEST(FiveBarRepresentability, DoNotScaleClosureToleranceByLongCranks) {
    const double length = std::ldexp(1.0, 50);
    FiveBarParams params;
    params.pivots = {Eigen::Vector2d::Zero(), Eigen::Vector2d(10.0, 20.0)};
    params.proximal = {length, length};
    const auto result = PlanarFiveBar(params).Analyze({length, 10.0}, Eigen::Vector2d::Zero());
    // Actual distal distances are 10 and sqrt(200), whereas both rods have length 1.
    EXPECT_EQ(result.status, FiveBarStatus::Unreachable);
    EXPECT_FALSE(result.success());
}

TEST(FiveBarReferences, DoggoPositiveExtensionAndVirtualLegJacobian) {
    // Geometry/formulas: Stanford Doggo's pinned ODrive ec25d02, controller.cpp
    // lines 174-195. This checks planar geometry, not raw encoder conventions.
    constexpr double upper = 0.09;
    constexpr double lower = 0.162;
    FiveBarParams params;
    params.pivots = {Eigen::Vector2d::Zero(), Eigen::Vector2d::Zero()};
    params.proximal = {upper, upper};
    params.distal = {lower, lower};
    const PlanarFiveBar mechanism(params);
    Eigen::Matrix2d coordinate_map;
    coordinate_map << -1.0, -1.0, -1.0, 1.0;
    for (double theta : {-0.45, 0.0, 0.35}) {
        for (double gamma : {-1.0, -0.55, 0.55, 1.0}) {
            const double radicand = lower * lower - upper * upper * std::sin(gamma) * std::sin(gamma);
            const double radius = upper * std::cos(gamma) + std::sqrt(radicand);
            const double radius_derivative = -upper * std::sin(gamma) -
                upper * upper * std::sin(gamma) * std::cos(gamma) / std::sqrt(radicand);
            const Eigen::Vector2d position(radius * std::sin(theta), radius * std::cos(theta));
            const Eigen::Vector2d angles(kPi / 2.0 - (theta + gamma), kPi / 2.0 - (theta - gamma));
            CheckRods(params, position, angles);
            const auto forward = mechanism.ForwardAll(angles);
            ASSERT_TRUE(forward.success());
            ASSERT_EQ(forward.count, 2U);
            bool matched_positive_extension = false;
            for (std::size_t assembly = 0; assembly < forward.count; ++assembly) {
                matched_positive_extension = matched_positive_extension ||
                    (forward.solutions[assembly].position - position).norm() < 1e-13;
            }
            EXPECT_TRUE(matched_positive_extension);
            const auto differential = mechanism.Analyze(position, angles);
            ASSERT_TRUE(differential.success());
            Eigen::Matrix2d reference;
            reference << radius * std::cos(theta), radius_derivative * std::sin(theta),
                         -radius * std::sin(theta), radius_derivative * std::cos(theta);
            EXPECT_TRUE((differential.motor_jacobian * reference).isApprox(coordinate_map, 2e-12));
            const int sign = gamma > 0.0 ? 1 : -1;
            EXPECT_EQ(differential.branch.working[0], sign);
            EXPECT_EQ(differential.branch.working[1], -sign);
            EXPECT_EQ(differential.branch.assembly, -sign);
        }
    }
}

}  // namespace
}  // namespace parallel_ankle::test
