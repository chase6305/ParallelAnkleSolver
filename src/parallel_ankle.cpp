#include "parallel_ankle/parallel_ankle.h"

#include <Eigen/Geometry>
#include <Eigen/LU>
#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>
#include <string>
#include <utility>

#include "ankle_derivatives.h"
#include "crank_geometry.h"

namespace parallel_ankle {
namespace {

constexpr double kSingularThreshold = 1e-10;
constexpr double kAngleRepresentationTolerance = 1e-12;
constexpr double kTransmissionTolerance = 128.0 * std::numeric_limits<double>::epsilon();

struct FootTransform {
    Eigen::Matrix3d pitch;
    Eigen::Matrix3d roll;
};

FootTransform MakeTransform(double pitch, double roll) {
    return {Eigen::AngleAxisd(pitch, Eigen::Vector3d::UnitY()).toRotationMatrix(),
            Eigen::AngleAxisd(roll, Eigen::Vector3d::UnitX()).toRotationMatrix()};
}

struct ActuatorPose {
    Eigen::Vector3d rolled;
    Eigen::Vector3d pitch_relative;
    Eigen::Vector3d delta;
    double bar;
    double rod;
    double side;
    double direction;
    double gamma = std::numeric_limits<double>::quiet_NaN();
    double theta = std::numeric_limits<double>::quiet_NaN();
    InverseStatus status = InverseStatus::Unreachable;
};

ActuatorPose SolveActuatorPose(const FootTransform& transform,
                               const Eigen::Vector3d& motor_center,
                               const Eigen::Vector3d& foot_point,
                               double bar,
                               double rod,
                               double offset,
                               bool is_left_motor,
                               const AnkleParams& geometry) {
    const Eigen::Vector3d rolled = transform.roll * (foot_point - geometry.r_roll_center);
    // Keep every sum in relative coordinates. Forming a huge world-space C
    // before subtracting A can erase an otherwise representable local motion.
    const Eigen::Vector3d pitch_relative = transform.pitch *
                                          (rolled + (geometry.r_roll_center - geometry.r_pitch_center));
    const Eigen::Vector3d delta = pitch_relative + (geometry.r_pitch_center - motor_center);
    const double initial_side = foot_point.y() - motor_center.y();
    const double side = std::copysign(1.0, initial_side);
    const double direction = is_left_motor ? 1.0 : -1.0;
    ActuatorPose result{rolled, pitch_relative, delta, bar, rod, side, direction};
    if (!delta.allFinite()) {
        return result;
    }
    if (initial_side == 0.0) {
        result.status = InverseStatus::Singular;
        return result;
    }
    const auto phase = detail::SolveCrankPhase(delta, bar, rod, side);
    if (phase.status != InverseStatus::Success) {
        result.status = phase.status;
        return result;
    }
    const double gamma = phase.gamma;
    const double theta = direction * gamma + offset;
    // A finite but very large encoder offset can erase gamma in theta=offset+gamma.
    // Validate the angle that a caller can actually reconstruct from our output.
    const double encoded_gamma = direction * (theta - offset);
    if (!std::isfinite(theta) || !std::isfinite(encoded_gamma) ||
        std::abs(encoded_gamma - gamma) > kAngleRepresentationTolerance) {
        result.status = InverseStatus::NumericalFailure;
        return result;
    }
    // The phase solver already checked closure and the working root at gamma.
    if (encoded_gamma != gamma) {
        const auto encoded_phase = detail::CheckCrankPhase(delta, bar, rod, side, encoded_gamma);
        if (!encoded_phase.rod_closes || !encoded_phase.selected_working_root) {
            result.status = InverseStatus::NumericalFailure;
            return result;
        }
    }
    result.gamma = gamma;
    result.theta = theta;
    result.status = InverseStatus::Success;
    return result;
}

struct AnklePose {
    Eigen::Matrix3d pitch_rotation;
    std::array<ActuatorPose, 2> actuators;

    bool success() const {
        return actuators[0].status == InverseStatus::Success && actuators[1].status == InverseStatus::Success;
    }

    Eigen::Vector2d motor_angles() const { return {actuators[0].theta, actuators[1].theta}; }
};

AnklePose SolveAnklePose(const FootTransform& transform, const AnkleParams& p) {
    return {transform.pitch,
            {SolveActuatorPose(transform, p.r_a1, p.r_c1, p.l_bar1, p.l_rod1, p.alpha_1_rad, p.is_left_motor1, p),
             SolveActuatorPose(transform, p.r_a2, p.r_c2, p.l_bar2, p.l_rod2, p.alpha_2_rad, p.is_left_motor2, p)}};
}

struct ActuatorDifferential {
    Eigen::Vector3d rolled;
    Eigen::Vector3d r;
    Eigen::Vector3d crank;
    Eigen::Vector3d crank_derivative;
    Eigen::Vector3d link;
    Eigen::Matrix<double, 3, 2> first;
    Eigen::RowVector2d jacobian;
    double length;
    double transmission;
};

bool DifferentiateActuator(const ActuatorPose& actuator,
                          const Eigen::Matrix3d& pitch_rotation,
                          ActuatorDifferential& differential) {
    const double scale = std::max({actuator.delta.cwiseAbs().maxCoeff(), actuator.bar, actuator.rod});
    differential.r = actuator.delta / scale;
    differential.length = actuator.bar / scale;
    const double sine = std::sin(actuator.gamma);
    const double cosine = std::cos(actuator.gamma);
    differential.crank = Eigen::Vector3d(0.0, actuator.side * cosine, sine);
    differential.crank_derivative = Eigen::Vector3d(0.0, -actuator.side * sine, cosine);
    differential.transmission = differential.length * differential.r.dot(differential.crank_derivative);
    const double singular_threshold = kTransmissionTolerance * differential.length * differential.r.stableNorm();
    if (std::abs(differential.transmission) <= singular_threshold) {
        return false;
    }
    differential.link = differential.r - differential.length * differential.crank;
    differential.rolled = actuator.rolled / scale;
    differential.first.col(0) = Eigen::Vector3d::UnitY().cross(actuator.pitch_relative / scale);
    differential.first.col(1) = pitch_rotation * Eigen::Vector3d::UnitX().cross(differential.rolled);
    for (int j = 0; j < 2; ++j) {
        differential.jacobian(j) = actuator.direction * differential.link.dot(differential.first.col(j)) /
                                   differential.transmission;
    }
    return differential.jacobian.allFinite();
}

bool DifferentiateAnkle(const AnklePose& pose, Eigen::Matrix2d& jacobian) {
    for (int i = 0; i < 2; ++i) {
        ActuatorDifferential differential;
        if (!DifferentiateActuator(pose.actuators[i], pose.pitch_rotation, differential)) {
            return false;
        }
        jacobian.row(i) = differential.jacobian;
    }
    return true;
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

Eigen::Matrix2d DifferentiateActuatorTwice(const ActuatorPose& actuator,
                                          const Eigen::Matrix3d& pitch_rotation,
                                          const ActuatorDifferential& differential) {
    const auto& first = differential.first;
    const Eigen::Vector3d pp = Eigen::Vector3d::UnitY().cross(first.col(0));
    const Eigen::Vector3d pr = Eigen::Vector3d::UnitY().cross(first.col(1));
    const Eigen::Vector3d rr = pitch_rotation * Eigen::Vector3d::UnitX().cross(
                                                   Eigen::Vector3d::UnitX().cross(differential.rolled));
    const Eigen::Vector2d gradient = actuator.direction * differential.jacobian.transpose();
    Eigen::Matrix2d hessian;
    for (int j = 0; j < 2; ++j) {
        for (int k = j; k < 2; ++k) {
            const Eigen::Vector3d& second = j == 0 && k == 0 ? pp : (j == 1 && k == 1 ? rr : pr);
            const double numerator = first.col(j).dot(first.col(k)) + differential.link.dot(second) -
                                     differential.length *
                                         (first.col(j).dot(differential.crank_derivative) * gradient(k) +
                                          first.col(k).dot(differential.crank_derivative) * gradient(j)) +
                                     differential.length * differential.r.dot(differential.crank) *
                                         gradient(j) * gradient(k);
            const double value = actuator.direction * numerator / differential.transmission;
            hessian(j, k) = value;
            hessian(k, j) = value;
        }
    }
    if (!hessian.allFinite()) {
        throw std::overflow_error("ankle angle Hessian produced nonfinite entries");
    }
    return hessian;
}

}  // namespace

detail::AnkleDerivatives detail::EvaluateAnkleDerivatives(const AnkleParams& params, double pitch, double roll) {
    if (!std::isfinite(pitch) || !std::isfinite(roll)) {
        throw std::invalid_argument("second-order task angles must be finite");
    }
    const auto pose = SolveAnklePose(MakeTransform(pitch, roll), params);
    if (!pose.success()) {
        throw std::domain_error("cannot differentiate outside the supported ankle IK branch");
    }
    // Both orders differentiate the same selected geometric angle, before
    // the encoder's rounding can perturb the transmission near a fold.
    std::array<ActuatorDifferential, 2> differentials;
    AnkleDerivatives derivatives;
    for (int i = 0; i < 2; ++i) {
        if (!DifferentiateActuator(pose.actuators[i], pose.pitch_rotation, differentials[i])) {
            throw std::domain_error("ankle actuator transmission is singular");
        }
        derivatives.jacobian.row(i) = differentials[i].jacobian;
    }
    for (int i = 0; i < 2; ++i) {
        derivatives.hessians[i] = DifferentiateActuatorTwice(pose.actuators[i], pose.pitch_rotation, differentials[i]);
    }
    return derivatives;
}

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
    if (!std::isfinite(q_pitch_rad) || !std::isfinite(q_roll_rad)) {
        return {};
    }
    const auto pose = SolveAnklePose(MakeTransform(q_pitch_rad, q_roll_rad), params_);
    const auto& [first, second] = pose.actuators;
    AnkleIKSolution solution;
    solution.theta_1_rad = first.theta;
    solution.theta_2_rad = second.theta;
    solution.status_1 = first.status;
    solution.status_2 = second.status;
    solution.success_1 = solution.status_1 == InverseStatus::Success;
    solution.success_2 = solution.status_2 == InverseStatus::Success;
    return solution;
}

Eigen::Matrix2d PMAnkle::ComputeJacobian(double q_pitch_rad, double q_roll_rad) const {
    if (!std::isfinite(q_pitch_rad) || !std::isfinite(q_roll_rad)) {
        throw std::invalid_argument("Jacobian angles must be finite");
    }
    const auto pose = SolveAnklePose(MakeTransform(q_pitch_rad, q_roll_rad), params_);
    if (!pose.success()) {
        throw std::domain_error("cannot compute ankle Jacobian outside the valid IK branch");
    }
    Eigen::Matrix2d jacobian;
    if (!DifferentiateAnkle(pose, jacobian)) {
        throw std::domain_error("ankle actuator transmission is singular");
    }
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
    auto evaluation = SolveAnklePose(MakeTransform(result.pose(0), result.pose(1)), params_);
    if (!evaluation.success()) {
        const auto& [first, second] = evaluation.actuators;
        if (first.status == InverseStatus::NumericalFailure || second.status == InverseStatus::NumericalFailure) {
            result.status = ForwardStatus::NumericalFailure;
        } else if (first.status == InverseStatus::Singular || second.status == InverseStatus::Singular) {
            result.status = ForwardStatus::Singular;
        } else {
            result.status = ForwardStatus::Unreachable;
        }
        return result;
    }
    Eigen::Vector2d residual = evaluation.motor_angles() - target;
    result.residual_norm = residual.stableNorm();
    for (int iteration = 0; iteration <= options.max_iterations; ++iteration) {
        result.iterations = iteration;
        if (result.residual_norm <= options.tolerance_rad) {
            result.status = ForwardStatus::Success;
            return result;
        }
        if (iteration == options.max_iterations) {
            result.status = ForwardStatus::MaxIterations;
            return result;
        }

        Eigen::Matrix2d jacobian;
        if (!DifferentiateAnkle(evaluation, jacobian)) {
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
        for (int attempt = 0; attempt < 12; ++attempt, scale *= 0.5) {
            const Eigen::Vector2d candidate = result.pose - scale * step;
            if (!candidate.allFinite()) {
                continue;
            }
            auto candidate_evaluation = SolveAnklePose(MakeTransform(candidate(0), candidate(1)), params_);
            if (!candidate_evaluation.success()) {
                continue;
            }
            const Eigen::Vector2d candidate_residual = candidate_evaluation.motor_angles() - target;
            const double candidate_norm = candidate_residual.stableNorm();
            if (candidate_norm < result.residual_norm) {
                result.pose = candidate;
                result.residual_norm = candidate_norm;
                residual = candidate_residual;
                evaluation = std::move(candidate_evaluation);
                accepted = true;
                break;
            }
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
