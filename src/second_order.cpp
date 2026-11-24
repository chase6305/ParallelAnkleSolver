#include "parallel_ankle/second_order.h"

#include <Eigen/Eigenvalues>
#include <Eigen/LU>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

#include "ankle_derivatives.h"
#include "inverse_congruence.h"

namespace parallel_ankle {
namespace {

constexpr double kRankThreshold = 1e-10;
constexpr double kSymmetryTolerance = 1e-12;
constexpr double kSemidefiniteTolerance = 128.0 * std::numeric_limits<double>::epsilon();

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
    return detail::SymmetricAverage(matrix);
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
    if (!detail::TryTransformSymmetricGains(inverse, corrected, result.motor_stiffness)) {
        throw std::overflow_error("loaded stiffness is not representable as a finite double matrix");
    }
    result.task_positive_semidefinite = PositiveSemidefinite(task);
    result.motor_positive_semidefinite = PositiveSemidefinite(result.motor_stiffness);
    return result;
}

}  // namespace parallel_ankle
