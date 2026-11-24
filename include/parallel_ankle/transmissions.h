#pragma once

#include <Eigen/Core>
#include <array>
#include <cstddef>
#include <limits>

namespace parallel_ankle {

enum class TransmissionStatus { Success, Tangent, Singular, Continuum, Unreachable, InvalidInput, NumericalFailure };
const char* ToString(TransmissionStatus status);

struct RotaryRodParams {
    // All coordinates and lengths share one unit.
    Eigen::Vector3d pivot = Eigen::Vector3d::Zero();
    Eigen::Vector3d u = Eigen::Vector3d::UnitX();
    Eigen::Vector3d v = Eigen::Vector3d::UnitY();
    double crank_length = 1.0;
    double rod_length = 1.0;
    int motor_sign = 1;
    double motor_offset_rad = 0.0;

    // u/v must already be orthonormal within floating-point roundoff. The
    // solver normalizes that residual roundoff, preserving their orientation.
    // Invalid parameters throw std::invalid_argument.
    void Validate() const;
};

struct RotaryRodSolution {
    // Phases in [-pi, pi); each represents a family with arbitrary 2*pi winding.
    double crank_phase_rad = std::numeric_limits<double>::quiet_NaN();
    double motor_phase_rad = std::numeric_limits<double>::quiet_NaN();
    // Sign of (C-A).(-u*sin(gamma)+v*cos(gamma)); zero at tangency.
    int branch = 0;
    // Locally unwrapped phase gradients with respect to C, in rad per length
    // unit; the canonical scalar phase itself has a jump at +/-pi.
    Eigen::Vector3d crank_gradient = Eigen::Vector3d::Constant(std::numeric_limits<double>::quiet_NaN());
    Eigen::Vector3d motor_gradient = Eigen::Vector3d::Constant(std::numeric_limits<double>::quiet_NaN());
    bool derivative_valid = false;
};

struct RotaryRodResult {
    TransmissionStatus status = TransmissionStatus::InvalidInput;
    std::array<RotaryRodSolution, 2> solutions;
    std::size_t count = 0;
    bool success() const { return status == TransmissionStatus::Success || status == TransmissionStatus::Tangent; }
};

// Solves ||C - [A + L*(u*cos(gamma)+v*sin(gamma))]|| = R.
// Motor phase is motor_sign*gamma + motor_offset, modulo 2*pi. This is a
// single local transmission, not a complete mechanism/platform FK solver.
// Tangency returns one valid phase with invalid derivatives. A continuum has
// count=0: the solver does not choose an arbitrary representative angle.
// Returned crank and motor phases are checked against actual rod closure.
RotaryRodResult SolveRotaryRod(const RotaryRodParams& params, const Eigen::Vector3d& point);

struct PrismaticLengthParams {
    Eigen::Vector3d anchor = Eigen::Vector3d::Zero();
    double minimum_length = 0.0;
    double maximum_length = std::numeric_limits<double>::infinity();
    void Validate() const;
};

struct PrismaticLengthResult {
    TransmissionStatus status = TransmissionStatus::InvalidInput;
    double length = std::numeric_limits<double>::quiet_NaN();
    Eigen::Vector3d gradient = Eigen::Vector3d::Constant(std::numeric_limits<double>::quiet_NaN());
    bool derivative_valid = false;
    // A permitted zero length is geometrically valid but direction-singular.
    bool success() const { return status == TransmissionStatus::Success || status == TransmissionStatus::Singular; }
};

// Positive distance from anchor to point, with stroke bounds. Angular bearing
// limits, screw pitch/encoder calibration and multi-leg closure are external.
PrismaticLengthResult SolvePrismaticLength(const PrismaticLengthParams& params, const Eigen::Vector3d& point);

}  // namespace parallel_ankle
