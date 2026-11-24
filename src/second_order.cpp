#include "parallel_ankle/second_order.h"

#include <Eigen/Eigenvalues>
#include <Eigen/LU>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

#include "ankle_derivatives.h"

namespace parallel_ankle {
namespace {

constexpr double kRankThreshold = 1e-10;
constexpr double kSymmetryTolerance = 1e-12;
constexpr double kSemidefiniteTolerance = 128.0 * std::numeric_limits<double>::epsilon();

double FiniteMidpoint(double first, double second) {
    if (first == second) {
        return first;
    }
    constexpr double half_max = std::numeric_limits<double>::max() * 0.5;
    constexpr double twice_min = std::numeric_limits<double>::min() * 2.0;
    const double first_magnitude = std::abs(first);
    const double second_magnitude = std::abs(second);
    if (first_magnitude <= half_max && second_magnitude <= half_max) {
        // Add before halving to preserve subnormal midpoints, including the
        // ties-to-even midpoint of denorm_min and 2*denorm_min.
        return (first + second) * 0.5;
    }
    if (first_magnitude < twice_min) {
        return first + second * 0.5;
    }
    if (second_magnitude < twice_min) {
        return first * 0.5 + second;
    }
    // Both halves are normal here. This also handles opposite large signs
    // without overflowing either their sum or their difference.
    return first * 0.5 + second * 0.5;
}

Eigen::Matrix2d SymmetricAverage(const Eigen::Matrix2d& matrix) {
    Eigen::Matrix2d symmetric = matrix;
    const double off_diagonal = FiniteMidpoint(matrix(0, 1), matrix(1, 0));
    symmetric(0, 1) = off_diagonal;
    symmetric(1, 0) = off_diagonal;
    return symmetric;
}

Eigen::Matrix2d SymmetricInput(const Eigen::Matrix2d& matrix) {
    if (!matrix.allFinite()) {
        throw std::invalid_argument("task stiffness must be finite");
    }
    const double scale = matrix.cwiseAbs().maxCoeff();
    if (scale == 0.0) {
        return matrix;
    }
    const Eigen::Matrix2d normalized = matrix / scale;
    if ((normalized - normalized.transpose()).cwiseAbs().maxCoeff() > kSymmetryTolerance) {
        throw std::invalid_argument("task stiffness must be symmetric");
    }
    // Normalization is only for the relative symmetry check. Reconstructing
    // from it would erase small representable entries in a mixed-scale input.
    return SymmetricAverage(matrix);
}

bool PositiveSemidefinite(const Eigen::Matrix2d& matrix) {
    const double scale = matrix.cwiseAbs().maxCoeff();
    if (scale == 0.0) {
        return true;
    }
    const Eigen::Matrix2d normalized = matrix / scale;
    const Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> solver(normalized, Eigen::EigenvaluesOnly);
    if (solver.info() != Eigen::Success || !solver.eigenvalues().allFinite()) {
        throw std::overflow_error("stiffness spectral diagnostic failed");
    }
    const double spectral_scale = solver.eigenvalues().cwiseAbs().maxCoeff();
    return solver.eigenvalues().minCoeff() >= -kSemidefiniteTolerance * spectral_scale;
}

}  // namespace

std::array<Eigen::Matrix2d, 2> ComputeMotorAngleHessians(const PMAnkle& model, double pitch_rad, double roll_rad) {
    return detail::EvaluateAnkleDerivatives(model.parameters(), pitch_rad, roll_rad).hessians;
}

LoadedStaticStiffnessResult MapLoadedStaticTaskStiffness(const PMAnkle& model,
                                                       double pitch_rad,
                                                       double roll_rad,
                                                       const Eigen::Matrix2d& task_stiffness,
                                                       const Eigen::Vector2d& motor_preload) {
    const Eigen::Matrix2d task = SymmetricInput(task_stiffness);
    if (!motor_preload.allFinite()) {
        throw std::invalid_argument("motor preload must be finite");
    }
    const auto derivatives = detail::EvaluateAnkleDerivatives(model.parameters(), pitch_rad, roll_rad);
    Eigen::FullPivLU<Eigen::Matrix2d> factor(derivatives.jacobian);
    factor.setThreshold(kRankThreshold);
    if (!factor.isInvertible()) {
        throw std::domain_error("loaded stiffness requires an invertible ankle Jacobian");
    }
    const Eigen::Matrix2d inverse = factor.solve(Eigen::Matrix2d::Identity());
    LoadedStaticStiffnessResult result;
    result.geometric_stiffness = -motor_preload(0) * derivatives.hessians[0] -
                                 motor_preload(1) * derivatives.hessians[1];
    const Eigen::Matrix2d corrected = task - result.geometric_stiffness;
    if (!inverse.allFinite() || !result.geometric_stiffness.allFinite() || !corrected.allFinite()) {
        throw std::overflow_error("loaded stiffness mapping produced nonfinite entries");
    }
    // A small entry can underflow in the first double product and become
    // representable again after the second factor. Keep the congruence in
    // extended precision until both factors have been applied. This leaves
    // the preload correction's existing double arithmetic contract intact.
    using ExtendedMatrix = Eigen::Matrix<long double, 2, 2>;
    const ExtendedMatrix extended_inverse = inverse.cast<long double>();
    const ExtendedMatrix extended_corrected = corrected.cast<long double>();
    const ExtendedMatrix extended_motor = extended_inverse.transpose() * extended_corrected * extended_inverse;
    constexpr long double double_max = std::numeric_limits<double>::max();
    Eigen::Matrix2d motor;
    for (int row = 0; row < 2; ++row) {
        for (int column = 0; column < 2; ++column) {
            const long double value = extended_motor(row, column);
            if (!std::isfinite(value) || std::abs(value) > double_max) {
                throw std::overflow_error("loaded stiffness is not representable as a finite double matrix");
            }
            motor(row, column) = static_cast<double>(value);
        }
    }
    result.motor_stiffness = SymmetricAverage(motor);
    result.task_positive_semidefinite = PositiveSemidefinite(task);
    result.motor_positive_semidefinite = PositiveSemidefinite(result.motor_stiffness);
    return result;
}

}  // namespace parallel_ankle
