#pragma once

#include <Eigen/Core>
#include <cmath>
#include <limits>

namespace parallel_ankle::detail {

inline double FiniteMidpoint(double first, double second) {
    if (first == second) {
        return first;
    }
    constexpr double half_max = std::numeric_limits<double>::max() * 0.5;
    constexpr double twice_min = std::numeric_limits<double>::min() * 2.0;
    const double first_magnitude = std::abs(first);
    const double second_magnitude = std::abs(second);
    if (first_magnitude <= half_max && second_magnitude <= half_max) {
        // Add before halving to preserve subnormal midpoints.
        return (first + second) * 0.5;
    }
    if (first_magnitude < twice_min) {
        return first + second * 0.5;
    }
    if (second_magnitude < twice_min) {
        return first * 0.5 + second;
    }
    return first * 0.5 + second * 0.5;
}

inline Eigen::Matrix2d SymmetricAverage(const Eigen::Matrix2d& matrix) {
    Eigen::Matrix2d symmetric = matrix;
    const double off_diagonal = FiniteMidpoint(matrix(0, 1), matrix(1, 0));
    symmetric(0, 1) = off_diagonal;
    symmetric(1, 0) = off_diagonal;
    return symmetric;
}

inline bool TryTransformSymmetricGains(const Eigen::Matrix2d& inverse,
                                       const Eigen::Matrix2d& task,
                                       Eigen::Matrix2d& motor) {
    // A factor can move an entry below the double range while the second
    // factor brings the final result back into range. Narrow only at the end.
    using ExtendedMatrix = Eigen::Matrix<long double, 2, 2>;
    const ExtendedMatrix extended_inverse = inverse.cast<long double>();
    const ExtendedMatrix extended_task = task.cast<long double>();
    const ExtendedMatrix extended_motor = extended_inverse.transpose() * extended_task * extended_inverse;
    constexpr long double double_max = std::numeric_limits<double>::max();
    Eigen::Matrix2d rounded;
    for (int row = 0; row < 2; ++row) {
        for (int column = 0; column < 2; ++column) {
            const long double value = extended_motor(row, column);
            if (!std::isfinite(value) || std::abs(value) > double_max) {
                return false;
            }
            rounded(row, column) = static_cast<double>(value);
        }
    }
    motor = SymmetricAverage(rounded);
    return true;
}

}  // namespace parallel_ankle::detail
