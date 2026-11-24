#include "parallel_ankle/parallel_ankle.h"

#include <Eigen/Geometry>
#include <Eigen/LU>
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

#include "crank_geometry.h"

namespace parallel_ankle {
namespace {

constexpr double kSingularThreshold = 1e-10;
constexpr double kAngleRepresentationTolerance = 1e-12;

struct FootTransform {
    Eigen::Matrix3d pitch;
    Eigen::Matrix3d roll;
};

FootTransform MakeTransform(double pitch, double roll) {
    return {Eigen::AngleAxisd(pitch, Eigen::Vector3d::UnitY()).toRotationMatrix(),
            Eigen::AngleAxisd(roll, Eigen::Vector3d::UnitX()).toRotationMatrix()};
}

struct ActuatorSolution {
    double theta = std::numeric_limits<double>::quiet_NaN();
    InverseStatus status = InverseStatus::Unreachable;
    Eigen::RowVector2d jacobian = Eigen::RowVector2d::Constant(std::numeric_limits<double>::quiet_NaN());
    bool differentiable = false;
};

ActuatorSolution SolveSingle(const FootTransform& transform,
                            const Eigen::Vector3d& motor_center,
                            const Eigen::Vector3d& foot_point,
                            double bar,
                            double rod,
                            double offset,
                            bool is_left_motor,
                            const AnkleParams& geometry,
                            bool differentiate) {
    ActuatorSolution result;
    const Eigen::Vector3d rolled = transform.roll * (foot_point - geometry.r_roll_center);
    // Keep every sum in relative coordinates. Forming a huge world-space C
    // before subtracting A can erase an otherwise representable local motion.
    const Eigen::Vector3d pitch_relative = transform.pitch *
                                          (rolled + (geometry.r_roll_center - geometry.r_pitch_center));
    const Eigen::Vector3d delta = pitch_relative + (geometry.r_pitch_center - motor_center);
    if (!delta.allFinite()) {
        return result;
    }
    const double initial_side = foot_point.y() - motor_center.y();
    if (initial_side == 0.0) {
        result.status = InverseStatus::Singular;
        return result;
    }
    const double side = std::copysign(1.0, initial_side);
    const auto phase = detail::SolveCrankPhase(delta, bar, rod, side);
    if (phase.status != InverseStatus::Success) {
        result.status = phase.status;
        return result;
    }
    const double gamma = phase.gamma;
    const double direction = is_left_motor ? 1.0 : -1.0;
    result.theta = direction * gamma + offset;
    // A finite but very large encoder offset can erase gamma in theta=offset+gamma.
    // Validate the angle that a caller can actually reconstruct from our output.
    const double encoded_gamma = direction * (result.theta - offset);
    if (!std::isfinite(result.theta) || !std::isfinite(encoded_gamma) ||
        std::abs(encoded_gamma - gamma) > kAngleRepresentationTolerance) {
        ActuatorSolution failure;
        failure.status = InverseStatus::NumericalFailure;
        return failure;
    }
    // The phase solver already checked closure and the working root at gamma.
    if (encoded_gamma != gamma) {
        const auto encoded_phase = detail::CheckCrankPhase(delta, bar, rod, side, encoded_gamma);
        if (!encoded_phase.rod_closes || !encoded_phase.selected_working_root) {
            ActuatorSolution failure;
            failure.status = InverseStatus::NumericalFailure;
            return failure;
        }
    }
    result.status = InverseStatus::Success;
    if (differentiate) {
        const double scale = std::max({delta.cwiseAbs().maxCoeff(), bar, rod});
        const Eigen::Vector3d r = delta / scale;
        const double length = bar / scale;
        const double sine = std::sin(gamma);
        const double cosine = std::cos(gamma);
        const Eigen::Vector3d crank(0.0, side * cosine, sine);
        const Eigen::Vector3d link = r - length * crank;
        const Eigen::Vector3d crank_derivative(0.0, -side * sine, cosine);
        const double transmission = length * r.dot(crank_derivative);
        if (std::abs(transmission) <= detail::kCrankBranchTolerance * length * r.stableNorm()) {
            return result;
        }
        const Eigen::Vector3d c_pitch = Eigen::Vector3d::UnitY().cross(pitch_relative / scale);
        const Eigen::Vector3d c_roll = transform.pitch * Eigen::Vector3d::UnitX().cross(rolled / scale);
        result.jacobian << direction * link.dot(c_pitch) / transmission,
                           direction * link.dot(c_roll) / transmission;
        result.differentiable = result.jacobian.allFinite();
    }
    return result;
}

std::pair<ActuatorSolution, ActuatorSolution> SolveActuators(double pitch,
                                                          double roll,
                                                          const AnkleParams& p,
                                                          bool differentiate = false) {
    if (!std::isfinite(pitch) || !std::isfinite(roll)) {
        ActuatorSolution invalid;
        invalid.status = InverseStatus::InvalidInput;
        return {invalid, invalid};
    }
    const auto transform = MakeTransform(pitch, roll);
    return {SolveSingle(transform, p.r_a1, p.r_c1, p.l_bar1, p.l_rod1, p.alpha_1_rad, p.is_left_motor1, p, differentiate),
            SolveSingle(transform, p.r_a2, p.r_c2, p.l_bar2, p.l_rod2, p.alpha_2_rad, p.is_left_motor2, p, differentiate)};
}

Eigen::FullPivLU<Eigen::Matrix2d> FactorJacobian(const Eigen::Matrix2d& matrix) {
    if (!matrix.allFinite()) {
        throw std::domain_error("ankle Jacobian is not finite");
    }
    Eigen::FullPivLU<Eigen::Matrix2d> factor(matrix);
    factor.setThreshold(kSingularThreshold);
    if (!factor.isInvertible()) {
        throw std::domain_error("ankle Jacobian is singular");
    }
    return factor;
}

Eigen::Vector2d CheckedVector(double first, double second) {
    const Eigen::Vector2d vector(first, second);
    if (!vector.allFinite()) {
        throw std::invalid_argument("torques and velocities must be finite");
    }
    return vector;
}

Eigen::Vector2d SolveLinear(const Eigen::Matrix2d& matrix, const Eigen::Vector2d& vector) {
    const Eigen::Vector2d result = FactorJacobian(matrix).solve(vector);
    if (!result.allFinite()) {
        throw std::overflow_error("ankle mapping produced a nonfinite result");
    }
    return result;
}

}  // namespace

const char* ToString(InverseStatus status) {
    switch (status) {
        case InverseStatus::Success: return "success";
        case InverseStatus::InvalidInput: return "invalid input";
        case InverseStatus::Unreachable: return "rod closure is unreachable";
        case InverseStatus::Singular: return "undefined crank angle";
        case InverseStatus::OutsideBranch: return "outside the supported crank branch";
        case InverseStatus::NumericalFailure: return "motor angle cannot represent the rod closure";
    }
    return "unknown status";
}

const char* ToString(ForwardStatus status) {
    switch (status) {
        case ForwardStatus::Success:
            return "success";
        case ForwardStatus::InvalidInput:
            return "invalid input";
        case ForwardStatus::Unreachable:
            return "unreachable iterate or Jacobian stencil";
        case ForwardStatus::Singular:
            return "singular Jacobian";
        case ForwardStatus::NoProgress:
            return "line search made no progress";
        case ForwardStatus::MaxIterations:
            return "iteration limit reached";
        case ForwardStatus::NumericalFailure:
            return "numerical representation failure";
    }
    return "unknown status";
}

PMAnkle::PMAnkle(bool is_left, const AnkleParams& params) : is_left_(is_left), params_(params) {
    params_.Validate();
}

AnkleIKSolution PMAnkle::GetAnkleMotorAngles(double q_pitch_rad, double q_roll_rad) const {
    const auto motors = SolveActuators(q_pitch_rad, q_roll_rad, params_);
    AnkleIKSolution solution;
    solution.theta_1_rad = motors.first.theta;
    solution.theta_2_rad = motors.second.theta;
    solution.status_1 = motors.first.status;
    solution.status_2 = motors.second.status;
    solution.success_1 = solution.status_1 == InverseStatus::Success;
    solution.success_2 = solution.status_2 == InverseStatus::Success;
    return solution;
}

Eigen::Matrix2d PMAnkle::ComputeJacobian(double q_pitch_rad, double q_roll_rad) const {
    if (!std::isfinite(q_pitch_rad) || !std::isfinite(q_roll_rad)) {
        throw std::invalid_argument("Jacobian angles must be finite");
    }
    const auto motors = SolveActuators(q_pitch_rad, q_roll_rad, params_, true);
    if (motors.first.status != InverseStatus::Success || motors.second.status != InverseStatus::Success) {
        throw std::domain_error("cannot compute ankle Jacobian outside the valid IK branch");
    }
    if (!motors.first.differentiable || !motors.second.differentiable) {
        throw std::domain_error("ankle actuator transmission is singular");
    }
    Eigen::Matrix2d jacobian;
    jacobian.row(0) = motors.first.jacobian;
    jacobian.row(1) = motors.second.jacobian;
    return jacobian;
}

Eigen::Matrix2d PMAnkle::ComputeNumericalJacobian(double q_pitch_rad, double q_roll_rad, double eps) const {
    if (!std::isfinite(eps) || eps <= 0.0 || !std::isfinite(q_pitch_rad) || !std::isfinite(q_roll_rad)) {
        throw std::invalid_argument("Jacobian angles must be finite and its step must be positive and finite");
    }
    const double pitch_step = (q_pitch_rad + eps) - q_pitch_rad;
    const double roll_step = (q_roll_rad + eps) - q_roll_rad;
    if (!std::isfinite(pitch_step) || !std::isfinite(roll_step) || pitch_step == 0.0 || roll_step == 0.0) {
        throw std::invalid_argument("Jacobian step does not produce a finite representable perturbation");
    }
    const auto zero = GetAnkleMotorAngles(q_pitch_rad, q_roll_rad);
    const auto pitch = GetAnkleMotorAngles(q_pitch_rad + eps, q_roll_rad);
    const auto roll = GetAnkleMotorAngles(q_pitch_rad, q_roll_rad + eps);
    if (!zero.success() || !pitch.success() || !roll.success()) {
        throw std::domain_error("cannot compute ankle Jacobian at an unreachable IK stencil");
    }
    Eigen::Matrix2d jacobian;
    jacobian.col(0) = (pitch.motor_angles() - zero.motor_angles()) / pitch_step;
    jacobian.col(1) = (roll.motor_angles() - zero.motor_angles()) / roll_step;
    if (!jacobian.allFinite()) {
        throw std::domain_error("ankle Jacobian is not finite");
    }
    return jacobian;
}

Eigen::Matrix2d PMAnkle::ComputeJacobian(double q_pitch_rad, double q_roll_rad, double eps) const {
    return ComputeNumericalJacobian(q_pitch_rad, q_roll_rad, eps);
}

AnkleFKSolution PMAnkle::SolveForwardKinematics(double theta_1_rad,
                                              double theta_2_rad,
                                              const ForwardOptions& options) const {
    AnkleFKSolution result;
    const Eigen::Vector2d target(theta_1_rad, theta_2_rad);
    if (!target.allFinite() || !options.initial_guess.allFinite() || !std::isfinite(options.tolerance_rad) ||
        options.tolerance_rad <= 0.0 || options.max_iterations < 0) {
        return result;
    }
    result.pose = options.initial_guess;
    const double half_pi = 0.5 * std::acos(-1.0);
    if (theta_1_rad < params_.alpha_1_rad - half_pi || theta_1_rad > params_.alpha_1_rad + half_pi ||
        theta_2_rad < params_.alpha_2_rad - half_pi || theta_2_rad > params_.alpha_2_rad + half_pi) {
        result.status = ForwardStatus::Unreachable;
        return result;
    }
    for (int iteration = 0; iteration <= options.max_iterations; ++iteration) {
        result.iterations = iteration;
        const auto ik = GetAnkleMotorAngles(result.pose(0), result.pose(1));
        if (!ik.success()) {
            if (ik.status_1 == InverseStatus::NumericalFailure || ik.status_2 == InverseStatus::NumericalFailure) {
                result.status = ForwardStatus::NumericalFailure;
            } else if (ik.status_1 == InverseStatus::Singular || ik.status_2 == InverseStatus::Singular) {
                result.status = ForwardStatus::Singular;
            } else {
                result.status = ForwardStatus::Unreachable;
            }
            return result;
        }
        const Eigen::Vector2d residual = ik.motor_angles() - target;
        result.residual_norm = residual.stableNorm();
        if (result.residual_norm <= options.tolerance_rad) {
            result.status = ForwardStatus::Success;
            return result;
        }
        if (iteration == options.max_iterations) {
            result.status = ForwardStatus::MaxIterations;
            return result;
        }

        Eigen::Matrix2d jacobian;
        try {
            jacobian = ComputeJacobian(result.pose(0), result.pose(1));
        } catch (const std::domain_error&) {
            result.status = ForwardStatus::Singular;
            return result;
        }
        Eigen::Vector2d step;
        try {
            step = SolveLinear(jacobian, residual);
        } catch (const std::domain_error&) {
            result.status = ForwardStatus::Singular;
            return result;
        } catch (const std::overflow_error&) {
            result.status = ForwardStatus::NoProgress;
            return result;
        }

        bool accepted = false;
        double scale = 1.0;
        for (int attempt = 0; attempt < 12; ++attempt) {
            const Eigen::Vector2d candidate = result.pose - scale * step;
            const auto candidate_ik = GetAnkleMotorAngles(candidate(0), candidate(1));
            if (candidate_ik.success() && (candidate_ik.motor_angles() - target).stableNorm() < result.residual_norm) {
                result.pose = candidate;
                accepted = true;
                break;
            }
            scale *= 0.5;
        }
        if (!accepted) {
            result.status = ForwardStatus::NoProgress;
            return result;
        }
    }
    return result;
}

Eigen::Vector2d PMAnkle::GetAnklePoseFromMotorAngles(double theta_1_rad, double theta_2_rad) const {
    const auto result = SolveForwardKinematics(theta_1_rad, theta_2_rad);
    if (!result.success()) {
        throw std::runtime_error(std::string("ankle forward kinematics failed: ") + ToString(result.status));
    }
    return result.pose;
}

Eigen::Vector2d PMAnkle::GetTaskTorquesFromMotorTorques(double q_pitch_rad,
                                                     double q_roll_rad,
                                                     double tau_motor_1,
                                                     double tau_motor_2) const {
    const Eigen::Vector2d torque = CheckedVector(tau_motor_1, tau_motor_2);
    const Eigen::Vector2d result = ComputeJacobian(q_pitch_rad, q_roll_rad).transpose() * torque;
    if (!result.allFinite()) {
        throw std::overflow_error("ankle torque mapping produced a nonfinite result");
    }
    return result;
}

Eigen::Vector2d PMAnkle::GetMotorTorquesFromTaskTorques(double q_pitch_rad,
                                                     double q_roll_rad,
                                                     double tau_task_pitch,
                                                     double tau_task_roll) const {
    const Eigen::Vector2d torque = CheckedVector(tau_task_pitch, tau_task_roll);
    return SolveLinear(ComputeJacobian(q_pitch_rad, q_roll_rad).transpose(), torque);
}

Eigen::Vector2d PMAnkle::GetTaskVelocitiesFromMotorVelocities(double q_pitch_rad,
                                                          double q_roll_rad,
                                                          double vel_motor_1,
                                                          double vel_motor_2) const {
    const Eigen::Vector2d velocity = CheckedVector(vel_motor_1, vel_motor_2);
    return SolveLinear(ComputeJacobian(q_pitch_rad, q_roll_rad), velocity);
}

Eigen::Matrix2d PMAnkle::MapTaskGains(double q_pitch_rad, double q_roll_rad, double pitch, double roll) const {
    if (!std::isfinite(pitch) || !std::isfinite(roll) || pitch < 0.0 || roll < 0.0) {
        throw std::invalid_argument("stiffness and damping must be finite and nonnegative");
    }
    const Eigen::Matrix2d inverse = FactorJacobian(ComputeJacobian(q_pitch_rad, q_roll_rad))
                                        .solve(Eigen::Matrix2d::Identity());
    const Eigen::Matrix2d result = inverse.transpose() * Eigen::Vector2d(pitch, roll).asDiagonal() * inverse;
    if (!result.allFinite()) {
        throw std::overflow_error("ankle gain mapping produced a nonfinite result");
    }
    return result;
}

Eigen::Matrix2d PMAnkle::GetMotorStiffnessMatrixFromTaskStiffness(double q_pitch_rad,
                                                               double q_roll_rad,
                                                               double k_pitch,
                                                               double k_roll) const {
    return MapTaskGains(q_pitch_rad, q_roll_rad, k_pitch, k_roll);
}

Eigen::Matrix2d PMAnkle::GetMotorDampingMatrixFromTaskDamping(double q_pitch_rad,
                                                          double q_roll_rad,
                                                          double d_pitch,
                                                          double d_roll) const {
    return MapTaskGains(q_pitch_rad, q_roll_rad, d_pitch, d_roll);
}

}  // namespace parallel_ankle
