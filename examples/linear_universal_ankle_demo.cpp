#include <Eigen/Geometry>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

#include "parallel_ankle/constraints.h"
#include "parallel_ankle/transmissions.h"

namespace {

using parallel_ankle::constraints::DifferentiateImplicitly;
using parallel_ankle::constraints::EqualityProblem;
using parallel_ankle::constraints::SolveEqualities;
using parallel_ankle::constraints::SolveOptions;
using parallel_ankle::constraints::SolveResult;
using parallel_ankle::constraints::Status;

// Independent SI fixture for two UPS length drives and a common-center U joint.
// These are illustrative dimensions, not dimensions imported from HyRoDyn/RH5.
// q=[roll,pitch], R=Rx(roll)*Ry(pitch); all attachment coordinates are in metres.
struct Fixture {
    Eigen::Vector3d center = Eigen::Vector3d::Zero();
    std::array<Eigen::Vector3d, 2> base{Eigen::Vector3d(0.04, 0.03, 0.20),
                                       Eigen::Vector3d(0.04, -0.03, 0.20)};
    std::array<Eigen::Vector3d, 2> platform{Eigen::Vector3d(0.04, 0.03, 0.0),
                                           Eigen::Vector3d(0.04, -0.03, 0.0)};
    Eigen::Vector2d lower = Eigen::Vector2d::Constant(-0.3);
    Eigen::Vector2d upper = Eigen::Vector2d::Constant(0.3);
    double minimum_length = 0.15;
    double maximum_length = 0.25;
    double reference_length = 0.20;
};

struct Attachments {
    std::array<Eigen::Vector3d, 2> position;
    std::array<Eigen::Matrix<double, 3, 2>, 2> derivative;
};

bool CalculateAttachments(const Fixture& model, const Eigen::VectorXd& q, Attachments& output) {
    if (q.size() != 2 || !q.allFinite() || (q.array() < model.lower.array()).any() ||
        (q.array() > model.upper.array()).any()) {
        return false;
    }
    const Eigen::Matrix3d rx = Eigen::AngleAxisd(q(0), Eigen::Vector3d::UnitX()).toRotationMatrix();
    const Eigen::Matrix3d ry = Eigen::AngleAxisd(q(1), Eigen::Vector3d::UnitY()).toRotationMatrix();
    for (int leg = 0; leg < 2; ++leg) {
        const Eigen::Vector3d pitched = ry * model.platform[leg];
        const Eigen::Vector3d relative = rx * pitched;
        output.position[leg] = model.center + relative;
        output.derivative[leg].col(0) = Eigen::Vector3d::UnitX().cross(relative);
        output.derivative[leg].col(1) = rx * Eigen::Vector3d::UnitY().cross(pitched);
    }
    return true;
}

// Separate scalar Rx*Ry and Euclidean-distance oracle, independent of the
// primitive, analytic derivatives and equality solver used by this example.
Eigen::Vector2d ScalarLengthOracle(const Fixture& model, const Eigen::Vector2d& q) {
    const double sr = std::sin(q(0));
    const double cr = std::cos(q(0));
    const double sp = std::sin(q(1));
    const double cp = std::cos(q(1));
    Eigen::Vector2d lengths;
    for (int leg = 0; leg < 2; ++leg) {
        const auto& p = model.platform[leg];
        const double x = cp * p.x() + sp * p.z();
        const double z = -sp * p.x() + cp * p.z();
        const double dx = model.center.x() + x - model.base[leg].x();
        const double dy = model.center.y() + cr * p.y() - sr * z - model.base[leg].y();
        const double dz = model.center.z() + sr * p.y() + cr * z - model.base[leg].z();
        lengths(leg) = std::hypot(dx, dy, dz);
    }
    return lengths;
}

struct LengthIK {
    bool success = false;
    Eigen::Vector2d lengths = Eigen::Vector2d::Constant(std::numeric_limits<double>::quiet_NaN());
    Eigen::Matrix2d jacobian = Eigen::Matrix2d::Constant(std::numeric_limits<double>::quiet_NaN());
    std::string status;
};

LengthIK Inverse(const Fixture& model, const Eigen::Vector2d& q) {
    LengthIK result;
    Attachments attachments;
    if (!CalculateAttachments(model, q, attachments)) {
        result.status = "invalid task coordinates or outside the selected joint chart";
        return result;
    }
    for (int leg = 0; leg < 2; ++leg) {
        parallel_ankle::PrismaticLengthParams params;
        params.anchor = model.base[leg];
        params.minimum_length = model.minimum_length;
        params.maximum_length = model.maximum_length;
        const auto actuator = parallel_ankle::SolvePrismaticLength(params, attachments.position[leg]);
        if (!actuator.success() || !actuator.derivative_valid) {
            result.status = "leg " + std::to_string(leg + 1) + ": " + parallel_ankle::ToString(actuator.status);
            return result;
        }
        result.lengths(leg) = actuator.length;
        result.jacobian.row(leg) = actuator.gradient.transpose() * attachments.derivative[leg];
    }
    result.success = true;
    result.status = "success";
    return result;
}

bool ValidLengths(const Fixture& model, const Eigen::Vector2d& lengths) {
    return lengths.allFinite() && (lengths.array() >= model.minimum_length).all() &&
           (lengths.array() <= model.maximum_length).all();
}

struct Closure {
    Eigen::Vector2d residual;
    Eigen::Matrix2d q_partial;
    Eigen::Matrix2d length_partial = Eigen::Matrix2d::Zero();
};

bool EvaluateClosure(const Fixture& model, const Eigen::VectorXd& q, const Eigen::Vector2d& lengths, Closure& output) {
    Attachments attachments;
    if (!ValidLengths(model, lengths) || !CalculateAttachments(model, q, attachments)) {
        return false;
    }
    for (int leg = 0; leg < 2; ++leg) {
        const Eigen::Vector3d link = attachments.position[leg] - model.base[leg];
        // Raw squared-distance closure; units are m^2 and m^2/rad.
        output.residual(leg) = 0.5 * (link.squaredNorm() - lengths(leg) * lengths(leg));
        output.q_partial.row(leg) = link.transpose() * attachments.derivative[leg];
        output.length_partial(leg, leg) = -lengths(leg);
    }
    return output.residual.allFinite() && output.q_partial.allFinite();
}

SolveResult Forward(const Fixture& model, const Eigen::Vector2d& lengths, const Eigen::Vector2d& warm_start) {
    if (!ValidLengths(model, lengths)) {
        SolveResult rejected;
        rejected.status = Status::InvalidInput;
        rejected.message = "requested actuator length is outside its stroke";
        return rejected;
    }
    EqualityProblem problem;
    problem.variables = problem.equations = 2;
    problem.residual = [&model, lengths](const Eigen::VectorXd& q, Eigen::VectorXd& residual) {
        Closure closure;
        if (!EvaluateClosure(model, q, lengths, closure)) return false;
        residual = closure.residual;
        return true;
    };
    problem.jacobian = [&model, lengths](const Eigen::VectorXd& q, Eigen::MatrixXd& jacobian) {
        Closure closure;
        if (!EvaluateClosure(model, q, lengths, closure)) return false;
        jacobian = closure.q_partial;
        return true;
    };
    SolveOptions options;
    options.lower_bounds = model.lower;
    options.upper_bounds = model.upper;
    options.variable_scales = Eigen::Vector2d::Constant(0.3);
    options.residual_scales = Eigen::Vector2d::Constant(model.reference_length * model.reference_length);
    options.residual_tolerance = 1e-12;
    options.max_iterations = 40;
    return SolveEqualities(problem, warm_start, options);
}

void Require(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

struct Differentials {
    Eigen::Matrix2d length_from_q;
    Eigen::Matrix2d q_from_length;
};

Differentials Differentiate(const Fixture& model, const Eigen::Vector2d& q, const Eigen::Vector2d& lengths) {
    Closure closure;
    Require(EvaluateClosure(model, q, lengths, closure), "closure evaluation failed before differentiation");
    Require(closure.residual.norm() < 1e-12, "implicit differentiation requires an already closed pose");
    parallel_ankle::constraints::DerivativeOptions options;
    options.equation_scales = Eigen::Vector2d::Constant(model.reference_length * model.reference_length);
    options.variable_scales = Eigen::Vector2d::Constant(model.reference_length);
    const auto lengths_from_q = DifferentiateImplicitly(closure.length_partial, closure.q_partial, options);
    Require(lengths_from_q.success(), "d(length)/d(q): " + std::string(parallel_ankle::constraints::ToString(lengths_from_q.status)));
    options.variable_scales = Eigen::Vector2d::Constant(0.3);
    const auto q_from_lengths = DifferentiateImplicitly(closure.q_partial, closure.length_partial, options);
    Require(q_from_lengths.success(), "d(q)/d(length): " + std::string(parallel_ankle::constraints::ToString(q_from_lengths.status)));
    return {lengths_from_q.derivative, q_from_lengths.derivative};
}

void CheckJacobianAndPower(const Fixture& model, const Eigen::Vector2d& q) {
    const auto ik = Inverse(model, q);
    Require(ik.success, "length IK: " + ik.status);
    const auto oracle = ScalarLengthOracle(model, q);
    Require((ik.lengths - oracle).norm() < 2e-14, "length IK differs from scalar RxRy oracle");
    const auto differential = Differentiate(model, q, oracle);
    Require((differential.length_from_q - ik.jacobian).norm() < 2e-13, "implicit and primitive Jacobians disagree");
    Eigen::Matrix2d numerical;
    for (int coordinate = 0; coordinate < 2; ++coordinate) {
        Eigen::Vector2d plus = q;
        Eigen::Vector2d minus = q;
        plus(coordinate) += 1e-6;
        minus(coordinate) -= 1e-6;
        numerical.col(coordinate) = (ScalarLengthOracle(model, plus) - ScalarLengthOracle(model, minus)) /
                                    (plus(coordinate) - minus(coordinate));
    }
    Require((numerical - differential.length_from_q).norm() < 2e-10, "length Jacobian fails independent finite differences");
    const Eigen::Vector2d q_velocity(0.2, -0.1);
    const Eigen::Vector2d force(100.0, 150.0);
    const Eigen::Vector2d length_velocity = differential.length_from_q * q_velocity;
    const Eigen::Vector2d task_torque = differential.length_from_q.transpose() * force;
    Require(std::abs(force.dot(length_velocity) - task_torque.dot(q_velocity)) < 2e-13, "virtual-work mismatch");
    Require((differential.q_from_length * length_velocity - q_velocity).norm() < 2e-12, "implicit velocity inverse mismatch");
}

}  // namespace

int main() {
    try {
        const Fixture model;
        const auto neutral = Inverse(model, Eigen::Vector2d::Zero());
        Require(neutral.success, neutral.status);
        Eigen::Matrix2d expected;
        expected << -0.03, 0.04, 0.03, 0.04;
        Require((neutral.lengths - Eigen::Vector2d::Constant(0.20)).norm() < 2e-14, "neutral length golden mismatch");
        Require((neutral.jacobian - expected).norm() < 2e-14, "neutral Jacobian golden mismatch");
        const Eigen::Vector2d golden_velocity = neutral.jacobian * Eigen::Vector2d(0.2, -0.1);
        const Eigen::Vector2d golden_torque = neutral.jacobian.transpose() * Eigen::Vector2d(100.0, 150.0);
        Require((golden_velocity - Eigen::Vector2d(-0.010, 0.002)).norm() < 2e-14, "neutral velocity golden mismatch");
        Require((golden_torque - Eigen::Vector2d(1.5, 10.0)).norm() < 2e-13, "neutral torque golden mismatch");

        const std::array<Eigen::Vector2d, 6> trajectory{Eigen::Vector2d(0.0, 0.0), Eigen::Vector2d(0.04, -0.03),
                                                       Eigen::Vector2d(0.08, -0.06), Eigen::Vector2d(0.12, -0.08),
                                                       Eigen::Vector2d(0.10, -0.02), Eigen::Vector2d(0.02, 0.06)};
        Eigen::Vector2d preceding = Eigen::Vector2d::Zero();
        int forward_checks = 0;
        for (const auto& truth : trajectory) {
            CheckJacobianAndPower(model, truth);
            const Eigen::Vector2d lengths = ScalarLengthOracle(model, truth);
            const std::array<Eigen::Vector2d, 3> starts{
                Eigen::Vector2d::Zero(), Eigen::Vector2d(-0.05, 0.04), preceding};
            for (const Eigen::Vector2d& start : starts) {
                const auto solved = Forward(model, lengths, start);
                Require(solved.success(), "local FK: " + std::string(parallel_ankle::constraints::ToString(solved.status)) +
                                              ": " + solved.message);
                Require(solved.full_column_rank, "local FK ended at a rank-deficient pose");
                Require((solved.x - truth).norm() < 2e-9, "local FK recovered another pose in this fixture");
                Require((ScalarLengthOracle(model, solved.x) - lengths).norm() < 2e-11, "FK fails independent length closure");
                preceding = solved.x;
                ++forward_checks;
            }
        }

        Fixture restricted_stroke = model;
        restricted_stroke.maximum_length = 0.19;
        const auto outside_ik_stroke = Inverse(restricted_stroke, Eigen::Vector2d::Zero());
        Require(!outside_ik_stroke.success, "outside-stroke IK was not rejected");
        const auto outside_stroke = Forward(model, {0.30, 0.20}, Eigen::Vector2d::Zero());
        Require(outside_stroke.status == Status::InvalidInput, "outside-stroke request was not rejected");
        const auto outside_chart = Forward(model, neutral.lengths, {0.31, 0.0});
        Require(outside_chart.status == Status::InvalidInput && outside_chart.iterations == 0 &&
                    outside_chart.residual_evaluations == 0,
                "outside-chart initial state was not rejected before evaluation");
        // Reverse triangle inequality bounds |ell1-ell2| by
        // ||R(p1-p2)-(A1-A2)|| = .12*|sin(roll/2)| <= .12*sin(.15).
        // The requested .10 m difference is impossible in this chart, even
        // though both lengths separately satisfy the actuator strokes.
        const auto infeasible = Forward(model, {0.15, 0.25}, Eigen::Vector2d::Zero());
        Require(!infeasible.success() && infeasible.status != Status::InvalidInput,
                "chart-infeasible lengths did not reach geometric failure handling");

        const Eigen::Vector2d sample(0.12, -0.08);
        const auto sample_ik = Inverse(model, sample);
        Require(sample_ik.success, "sample IK: " + sample_ik.status);
        const auto differential = Differentiate(model, sample, sample_ik.lengths);
        std::cout << std::setprecision(10)
                  << "2UPS+U reduced spatial ankle; independent SI fixture, not robot hardware dimensions\n"
                  << "q=[roll,pitch] rad; R=Rx(roll)*Ry(pitch); actuator coordinates are lengths in m\n"
                  << "sample q: " << sample.transpose() << ", lengths: " << sample_ik.lengths.transpose() << '\n'
                  << "d(length)/d(q) [m/rad]:\n" << differential.length_from_q << '\n'
                  << "d(q)/d(length) [rad/m]:\n" << differential.q_from_length << '\n'
                  << "neutral virtual-work golden: task torque [1.5,10] Nm; actuator velocities [-.010,.002] m/s\n"
                  << "IK with restricted stroke: " << outside_ik_stroke.status << '\n'
                  << "outside stroke: " << parallel_ankle::constraints::ToString(outside_stroke.status) << '\n'
                  << "outside warm-start chart: " << parallel_ankle::constraints::ToString(outside_chart.status) << '\n'
                  << "infeasible lengths in selected chart: " << parallel_ankle::constraints::ToString(infeasible.status) << '\n'
                  << "self-checks passed: 6 poses, " << forward_checks
                  << " bounded local FK checks, scalar geometry/FD/virtual-work oracles\n";
    } catch (const std::exception& error) {
        std::cerr << "linear_universal_ankle_demo failed: " << error.what() << '\n';
        return 1;
    }
    return 0;
}
