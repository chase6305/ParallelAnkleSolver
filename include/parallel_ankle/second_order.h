#pragma once

#include <Eigen/Core>
#include <array>

#include "parallel_ankle/parallel_ankle.h"

namespace parallel_ankle {

// H[i] = d²(theta_i)/d([pitch,roll])² on the model's current smooth branch.
// This optional calculation does not add work to ordinary IK/FK/Jacobian calls.
// An actuator transmission singularity or nonfinite derivative is rejected;
// conditioning can deteriorate rapidly before a singularity is reached.
std::array<Eigen::Matrix2d, 2> ComputeMotorAngleHessians(const PMAnkle& model,
                                                      double pitch_rad,
                                                      double roll_rad);

struct LoadedStaticStiffnessResult {
    Eigen::Matrix2d motor_stiffness;
    // Actual task-space preload contribution: -sum_i(tau_motor_i * H_i).
    // K_task = J^T * K_motor * J + geometric_stiffness.
    Eigen::Matrix2d geometric_stiffness;
    // Relative-tolerance spectral diagnostics of the symmetric matrices.
    // These flags are not certificates of passivity or controller stability.
    bool task_positive_semidefinite = false;
    bool motor_positive_semidefinite = false;
};

// Static conservative restoring convention: tau_motor = -gradient(U_motor),
// K_motor = Hessian(U_motor), and J=d(theta)/dq. The motor preload must be at
// this same pose, with the same motor signs and units as the model.
//
// K_motor = J^(-T) * (K_task - geometric_stiffness) * J^(-1).
//
// task_stiffness must be finite and symmetric. Symmetric indefinite inputs
// and outputs are allowed and reported through the spectral flags. This is
// not a moving-state PD impedance transfer; no velocity terms are included.
// Throws invalid_argument for malformed inputs, domain_error for unsupported
// branches or singular mappings, and overflow_error for nonfinite arithmetic.
LoadedStaticStiffnessResult MapLoadedStaticTaskStiffness(const PMAnkle& model,
                                                       double pitch_rad,
                                                       double roll_rad,
                                                       const Eigen::Matrix2d& task_stiffness,
                                                       const Eigen::Vector2d& motor_preload);

}  // namespace parallel_ankle
