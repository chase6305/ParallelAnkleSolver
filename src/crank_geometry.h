#pragma once

#include <Eigen/Core>
#include <algorithm>
#include <cmath>
#include <limits>

#include "parallel_ankle/parallel_ankle.h"

namespace parallel_ankle::detail {

constexpr double kCrankBranchTolerance = 128.0 * std::numeric_limits<double>::epsilon();
constexpr long double kCrankRodTolerance = 512.0L * std::numeric_limits<double>::epsilon();

struct CrankPhase {
    InverseStatus status = InverseStatus::Unreachable;
    double gamma = std::numeric_limits<double>::quiet_NaN();
};

struct CrankPhaseCheck {
    bool rod_closes = false;
    bool selected_working_root = false;
};

// Use rod-relative physical and working-root residuals, never tolerances
// proportional to the square of a possibly much longer crank/displacement.
// Share the long-double trigonometric evaluation for both checks.
inline CrankPhaseCheck CheckCrankPhase(const Eigen::Vector3d& delta,
                                      double bar,
                                      double rod,
                                      double side,
                                      double gamma) {
    if (!std::isfinite(gamma)) {
        return {};
    }
    const long double phase = gamma;
    const long double sine = std::sin(phase);
    const long double cosine = std::cos(phase);
    const long double x = delta.x();
    const long double y = static_cast<long double>(delta.y()) -
                          static_cast<long double>(bar) * side * cosine;
    const long double z = static_cast<long double>(delta.z()) - static_cast<long double>(bar) * sine;
    const long double actual = std::hypot(x, y, z);
    // The positive square-root branch has rho*sin(opening) >= 0. Closure alone cannot
    // distinguish its neighboring root after encoder rounding. Only rod-scale
    // roundoff is admitted at tangency, even with an enormous crank length.
    const long double tolerance = kCrankRodTolerance * rod;
    const long double first = std::abs(static_cast<long double>(delta.y())) * sine;
    const long double second = static_cast<long double>(delta.z()) * cosine;
    const long double working = first - second;
    // Long-double fma may be implemented in software. Away from the decision
    // boundary, bound its difference from ordinary products/subtraction; only
    // an uncertain predicate needs the original fused evaluation. The factor
    // eight includes rounding in the bound and its endpoint comparisons.
    constexpr long double roundoff = 8.0L * std::numeric_limits<long double>::epsilon();
    const long double error = roundoff * std::abs(first) + roundoff * std::abs(second) +
                              8.0L * std::numeric_limits<long double>::min();
    bool selected = false;
    if (std::isfinite(working) && std::isfinite(error) && working - error >= -tolerance) {
        selected = true;
    } else if (!(std::isfinite(working) && std::isfinite(error) && working + error < -tolerance)) {
        const long double fused = std::fma(std::abs(static_cast<long double>(delta.y())), sine, -second);
        selected = std::isfinite(fused) && fused >= -tolerance;
    }
    return {std::isfinite(actual) && std::abs(actual - rod) <= tolerance,
            selected};
}

inline CrankPhase SolveCrankPhase(const Eigen::Vector3d& delta, double bar, double rod, double side) {
    if (!delta.allFinite()) {
        return {};
    }
    if (std::abs(delta.x()) > rod) {
        return {InverseStatus::Unreachable};
    }
    const double scale = std::max({delta.cwiseAbs().maxCoeff(), bar, rod});
    const Eigen::Vector3d r = delta / scale;
    if (side * r.y() < -kCrankBranchTolerance) {
        return {InverseStatus::OutsideBranch};
    }

    // Use the normalized closure formula for well-scaled geometry.
    // Lost short-rod terms are detected by the physical closure check below.
    const double length = bar / scale;
    const double rod_length = rod / scale;
    const double a = r.y();
    const double b = -r.z();
    const double radial_squared = a * a + b * b;
    if (length != 0.0 && rod_length != 0.0 &&
        radial_squared > kCrankBranchTolerance * kCrankBranchTolerance) {
        const double projected = (rod_length * rod_length - length * length - r.squaredNorm()) / (2.0 * length);
        const double discriminant = radial_squared - projected * projected;
        if (std::isfinite(discriminant) && discriminant >= -kCrankBranchTolerance * radial_squared) {
            const double root = std::sqrt(std::max(0.0, discriminant));
            const double sine = (b * projected + std::abs(a) * root) / radial_squared;
            const double cosine = (-std::abs(a) * projected + b * root) / radial_squared;
            if (std::isfinite(sine) && std::isfinite(cosine) && cosine >= -kCrankBranchTolerance) {
                const double gamma = std::atan2(sine, std::max(0.0, cosine));
                const auto checked = CheckCrankPhase(delta, bar, rod, side, gamma);
                if (checked.rod_closes && checked.selected_working_root) {
                    return {InverseStatus::Success, gamma};
                }
            }
        }
    }

    // The motor crank is a YZ circle. Intersect the rod sphere with that
    // plane, then solve the remaining planar triangle using stable half-angle
    // factors. Unscaled long doubles retain a short nonzero rod/opening when
    // subtracting squared normalized lengths would lose it completely.
    const long double x = std::abs(static_cast<long double>(delta.x()));
    const long double y = std::abs(static_cast<long double>(delta.y()));
    const long double z = delta.z();
    const long double crank = bar;
    const long double target_rod = rod;
    const long double projected_radius = std::hypot(y, z);
    if (projected_radius == 0.0L) {
        return {InverseStatus::Singular};
    }
    if (x > target_rod) {
        return {InverseStatus::Unreachable};
    }
    const long double radius = std::sqrt((target_rod - x) * (target_rod + x));
    const long double difference = std::abs(crank - projected_radius);
    const long double sum = crank + projected_radius;
    if (radius < difference || radius > sum) {
        return {InverseStatus::Unreachable};
    }
    const long double sine_square = (radius - difference) * (radius + difference);
    const long double cosine_square = (sum - radius) * (sum + radius);
    if (!std::isfinite(sine_square) || !std::isfinite(cosine_square) ||
        (sine_square == 0.0L && cosine_square == 0.0L)) {
        return {InverseStatus::NumericalFailure};
    }
    const long double opening = 2.0L * std::atan2(std::sqrt(sine_square), std::sqrt(cosine_square));
    const long double phase = std::atan2(z, y) + opening;
    const long double sine = std::sin(phase);
    const long double cosine = std::cos(phase);
    if (!std::isfinite(sine) || !std::isfinite(cosine)) {
        return {InverseStatus::NumericalFailure};
    }
    if (cosine < -kCrankBranchTolerance) {
        return {InverseStatus::OutsideBranch};
    }
    const double gamma = static_cast<double>(std::atan2(sine, std::max(0.0L, cosine)));
    const auto checked = CheckCrankPhase(delta, bar, rod, side, gamma);
    if (!checked.rod_closes || !checked.selected_working_root) {
        return {InverseStatus::NumericalFailure};
    }
    return {InverseStatus::Success, gamma};
}

}  // namespace parallel_ankle::detail
