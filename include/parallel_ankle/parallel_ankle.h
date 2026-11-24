#pragma once

#include <Eigen/Core>
#include <limits>
#include <utility>

#include "parallel_ankle/parameters.h"

namespace parallel_ankle {

enum class InverseStatus { Success, InvalidInput, Unreachable, Singular, OutsideBranch, NumericalFailure };

const char* ToString(InverseStatus status);

struct AnkleIKSolution {
    double theta_1_rad = std::numeric_limits<double>::quiet_NaN();
    bool success_1 = false;
    double theta_2_rad = std::numeric_limits<double>::quiet_NaN();
    bool success_2 = false;
    InverseStatus status_1 = InverseStatus::InvalidInput;
    InverseStatus status_2 = InverseStatus::InvalidInput;

    bool success() const { return success_1 && success_2; }
    Eigen::Vector2d motor_angles() const { return {theta_1_rad, theta_2_rad}; }
};

enum class ForwardStatus { Success, InvalidInput, Unreachable, Singular, NoProgress, MaxIterations, NumericalFailure };

const char* ToString(ForwardStatus status);

struct ForwardOptions {
    Eigen::Vector2d initial_guess = Eigen::Vector2d::Zero();
    double tolerance_rad = 1e-10;
    int max_iterations = 30;
};

struct AnkleFKSolution {
    Eigen::Vector2d pose = Eigen::Vector2d::Constant(std::numeric_limits<double>::quiet_NaN());
    ForwardStatus status = ForwardStatus::InvalidInput;
    int iterations = 0;
    double residual_norm = std::numeric_limits<double>::infinity();

    bool success() const { return status == ForwardStatus::Success; }
};

// Rotary-crank ankle solver with two fixed-length rods.
// q = [pitch, roll] and theta = [motor1, motor2], in radians.
// The is_left flag labels the ankle. Signs come only from is_left_motor1/2.
class PMAnkle {
public:
    explicit PMAnkle(bool is_left, const AnkleParams& params);

    bool is_left() const { return is_left_; }
    const AnkleParams& parameters() const { return params_; }

    AnkleIKSolution GetAnkleMotorAngles(double q_pitch_rad, double q_roll_rad) const;

    AnkleFKSolution SolveForwardKinematics(double theta_1_rad,
                                         double theta_2_rad,
                                         const ForwardOptions& options = {}) const;

    // Convenience API that throws if the forward solve fails.
    Eigen::Vector2d GetAnklePoseFromMotorAngles(double theta_1_rad, double theta_2_rad) const;

    // J = d(theta)/d([pitch, roll]), differentiated from the rod closure.
    // Throws for invalid IK or an actuator transmission singularity.
    Eigen::Matrix2d ComputeJacobian(double q_pitch_rad, double q_roll_rad) const;

    // Independent finite-difference reference; each stencil must be valid.
    Eigen::Matrix2d ComputeNumericalJacobian(double q_pitch_rad, double q_roll_rad, double eps = 1e-6) const;

    // Source-compatible explicit-step numerical overload.
    Eigen::Matrix2d ComputeJacobian(double q_pitch_rad, double q_roll_rad, double eps) const;

    Eigen::Vector2d GetTaskTorquesFromMotorTorques(double q_pitch_rad,
                                                double q_roll_rad,
                                                double tau_motor_1,
                                                double tau_motor_2) const;

    Eigen::Vector2d GetMotorTorquesFromTaskTorques(double q_pitch_rad,
                                                double q_roll_rad,
                                                double tau_task_pitch,
                                                double tau_task_roll) const;

    Eigen::Vector2d GetTaskVelocitiesFromMotorVelocities(double q_pitch_rad,
                                                       double q_roll_rad,
                                                       double vel_motor_1,
                                                       double vel_motor_2) const;

    // Local, unloaded linearization; the returned matrices retain coupling.
    Eigen::Matrix2d GetMotorStiffnessMatrixFromTaskStiffness(double q_pitch_rad,
                                                          double q_roll_rad,
                                                          double k_pitch,
                                                          double k_roll) const;

    Eigen::Matrix2d GetMotorDampingMatrixFromTaskDamping(double q_pitch_rad,
                                                      double q_roll_rad,
                                                      double d_pitch,
                                                      double d_roll) const;

private:
    Eigen::Matrix2d MapTaskGains(double q_pitch_rad, double q_roll_rad, double pitch, double roll) const;

    bool is_left_;
    AnkleParams params_;
};

}  // namespace parallel_ankle
