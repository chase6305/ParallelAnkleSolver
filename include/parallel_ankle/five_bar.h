#pragma once

#include <Eigen/Core>
#include <array>
#include <cstddef>
#include <limits>

namespace parallel_ankle {

enum class FiveBarStatus { Success, InvalidInput, Unreachable, Singular, NumericalFailure };
const char* ToString(FiveBarStatus status);

struct FiveBarParams {
    std::array<Eigen::Vector2d, 2> pivots{Eigen::Vector2d(-0.5, 0.0), Eigen::Vector2d(0.5, 0.0)};
    std::array<double, 2> proximal{1.0, 1.0};
    std::array<double, 2> distal{1.0, 1.0};

    // All coordinates/lengths share a unit. Positive finite lengths are required.
    void Validate() const;
};

struct FiveBarBranch {
    // Sign of cross(E2-E1, P-E1); zero at a tangent/coincident elbow geometry.
    int assembly = 0;
    // Each sign is cross(Ei-Ai, P-Ei); zero means a collinear branch boundary.
    std::array<int, 2> working{0, 0};
};

struct FiveBarForwardSolution {
    Eigen::Vector2d position = Eigen::Vector2d::Constant(std::numeric_limits<double>::quiet_NaN());
    FiveBarBranch branch;
};

struct FiveBarForwardResult {
    FiveBarStatus status = FiveBarStatus::InvalidInput;
    std::array<FiveBarForwardSolution, 2> solutions;
    std::size_t count = 0;
    bool success() const { return status == FiveBarStatus::Success; }
};

struct FiveBarInverseSolution {
    Eigen::Vector2d motor_angles = Eigen::Vector2d::Constant(std::numeric_limits<double>::quiet_NaN());
    FiveBarBranch branch;
};

struct FiveBarInverseResult {
    FiveBarStatus status = FiveBarStatus::InvalidInput;
    std::array<FiveBarInverseSolution, 4> solutions;
    std::size_t count = 0;
    bool success() const { return status == FiveBarStatus::Success; }
};

struct FiveBarDifferential {
    FiveBarStatus status = FiveBarStatus::InvalidInput;
    // Unit rod directions form A. B is diagonal: A*x_dot = B*theta_dot.
    Eigen::Matrix2d A = Eigen::Matrix2d::Constant(std::numeric_limits<double>::quiet_NaN());
    Eigen::Matrix2d B = Eigen::Matrix2d::Constant(std::numeric_limits<double>::quiet_NaN());
    // d(theta)/d([x,y]) = B^-1 A, in rad per configured length unit.
    // It remains finite at assembly-only singularities. If F = J^T*tau is used,
    // F has torque/length units: Nm/mm for millimetres, and N only for metres.
    Eigen::Matrix2d motor_jacobian = Eigen::Matrix2d::Constant(std::numeric_limits<double>::quiet_NaN());
    bool actuator_singular = false;
    bool assembly_singular = false;
    FiveBarBranch branch;
    bool success() const { return status == FiveBarStatus::Success; }
};

// Independent planar five-bar geometry, including unequal legs and coaxial bases.
// Motor angles are absolute CCW rotations from global +X, in radians.
// The two passive distal rods meet at the task point P = [x,y].
// Success outputs satisfy a distal-length-relative closure check. If double
// coordinates or motor angles cannot preserve that closure, the solve returns
// NumericalFailure/count=0; extended intermediates do not remove output limits.
class PlanarFiveBar {
public:
    explicit PlanarFiveBar(const FiveBarParams& params);
    const FiveBarParams& parameters() const { return params_; }

    // Returns all isolated finite solutions. A tangent is one Success solution;
    // coincident circles with infinitely many points return Singular/count=0.
    FiveBarForwardResult ForwardAll(const Eigen::Vector2d& motor_angles) const;

    // Enumerates both working modes per leg; tangent branches are deduplicated.
    // A target at a base pivot with equal proximal/distal lengths is nonunique.
    FiveBarInverseResult InverseAll(const Eigen::Vector2d& position) const;

    // Requires a closed configuration. Threshold is a dimensionless sine bound in (0,1).
    // Singular flags distinguish collinear crank/rod from collinear distal rods.
    FiveBarDifferential Analyze(const Eigen::Vector2d& position,
                                const Eigen::Vector2d& motor_angles,
                                double singular_tolerance = 1e-10) const;

private:
    FiveBarParams params_;
};

}  // namespace parallel_ankle
