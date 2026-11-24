#include "parallel_ankle/transmissions.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace parallel_ankle {
namespace {

using ExtendedVector = Eigen::Matrix<long double, 3, 1>;
constexpr long double kBasisTolerance = 64.0L * std::numeric_limits<double>::epsilon();
constexpr long double kArithmeticTolerance = 64.0L * std::numeric_limits<long double>::epsilon();
constexpr long double kRodTolerance = 512.0L * std::numeric_limits<double>::epsilon();

long double Norm(const ExtendedVector& value) { return std::hypot(value.x(), value.y(), value.z()); }

long double Pi() { return std::acos(-1.0L); }

long double CanonicalPhase(long double phase) {
    return std::atan2(std::sin(phase), std::cos(phase));
}

double NarrowPhase(long double phase) {
    double result = static_cast<double>(CanonicalPhase(phase));
    const double pi = static_cast<double>(Pi());
    if (result >= pi) {
        result -= 2.0 * pi;
    }
    return result;
}

bool NarrowGradient(const ExtendedVector& gradient, Eigen::Vector3d& result) {
    const long double bound = std::numeric_limits<double>::max();
    if (!gradient.allFinite() || gradient.cwiseAbs().maxCoeff() > bound) {
        return false;
    }
    result = gradient.cast<double>();
    return result.allFinite();
}

bool RodCloses(const ExtendedVector& residual, long double rod) {
    const long double actual = Norm(residual);
    return std::isfinite(actual) && std::abs(actual - rod) <= kRodTolerance * rod;
}

RotaryRodResult RotaryFailure(TransmissionStatus status) {
    RotaryRodResult result;
    result.status = status;
    return result;
}

}  // namespace

const char* ToString(TransmissionStatus status) {
    switch (status) {
        case TransmissionStatus::Success:
            return "success";
        case TransmissionStatus::Tangent:
            return "tangent configuration; derivative undefined";
        case TransmissionStatus::Singular:
            return "direction undefined";
        case TransmissionStatus::Continuum:
            return "continuum of crank phases";
        case TransmissionStatus::Unreachable:
            return "unreachable transmission target";
        case TransmissionStatus::InvalidInput:
            return "invalid input";
        case TransmissionStatus::NumericalFailure:
            return "result cannot be represented reliably";
    }
    return "unknown status";
}

void RotaryRodParams::Validate() const {
    if (!pivot.allFinite() || !u.allFinite() || !v.allFinite() || !std::isfinite(crank_length) ||
        !std::isfinite(rod_length) || crank_length <= 0.0 || rod_length <= 0.0 ||
        (motor_sign != 1 && motor_sign != -1) || !std::isfinite(motor_offset_rad)) {
        throw std::invalid_argument("rotary-rod geometry must be finite, lengths positive and motor sign +/-1");
    }
    const ExtendedVector first = u.cast<long double>();
    const ExtendedVector second = v.cast<long double>();
    if (std::abs(Norm(first) - 1.0L) > kBasisTolerance || std::abs(Norm(second) - 1.0L) > kBasisTolerance ||
        std::abs(first.dot(second)) > kBasisTolerance) {
        throw std::invalid_argument("rotary-rod plane basis must be orthonormal");
    }
}

RotaryRodResult SolveRotaryRod(const RotaryRodParams& params, const Eigen::Vector3d& point) {
    params.Validate();
    if (!point.allFinite()) {
        return RotaryFailure(TransmissionStatus::InvalidInput);
    }
    ExtendedVector u = params.u.cast<long double>();
    u /= Norm(u);
    ExtendedVector v = params.v.cast<long double>();
    v -= u * u.dot(v);
    v /= Norm(v);
    const ExtendedVector displacement = point.cast<long double>() - params.pivot.cast<long double>();
    if (!displacement.allFinite()) {
        return RotaryFailure(TransmissionStatus::NumericalFailure);
    }
    const long double scale = std::max({displacement.cwiseAbs().maxCoeff(),
                                        static_cast<long double>(params.crank_length),
                                        static_cast<long double>(params.rod_length)});
    const ExtendedVector delta = displacement / scale;
    const long double crank = params.crank_length / scale;
    const long double rod = params.rod_length / scale;
    const long double a = delta.dot(u);
    const long double b = delta.dot(v);
    const long double projected = std::hypot(a, b);
    const long double constraint = delta.squaredNorm() + crank * crank - rod * rod;
    if (projected == 0.0L) {
        const long double error_bound = kArithmeticTolerance *
                                        (delta.squaredNorm() + crank * crank + rod * rod);
        return RotaryFailure(std::abs(constraint) <= error_bound ? TransmissionStatus::Continuum
                                                               : TransmissionStatus::Unreachable);
    }

    ExtendedVector normal;
    normal << u.y() * v.z() - u.z() * v.y(), u.z() * v.x() - u.x() * v.z(),
        u.x() * v.y() - u.y() * v.x();
    normal /= Norm(normal);
    const long double normal_distance = std::abs(delta.dot(normal));
    const long double normal_products = delta.cwiseProduct(normal).cwiseAbs().sum();
    const long double plane_tolerance = kArithmeticTolerance * (normal_products + rod);
    const long double plane_gap = normal_distance - rod;
    if (plane_gap > plane_tolerance) {
        return RotaryFailure(TransmissionStatus::Unreachable);
    }
    if (plane_gap > 0.0L ||
        (std::abs(plane_gap) <= plane_tolerance && normal_products > 2.0L * normal_distance)) {
        // A positive or cancellation-dominated uncertain gap must not be
        // clamped into a feasible tangent circle.
        return RotaryFailure(TransmissionStatus::NumericalFailure);
    }
    const long double in_plane_radius = std::sqrt((rod - normal_distance) * (rod + normal_distance));
    const long double sum = crank + projected;
    const long double difference = std::abs(crank - projected);
    const long double triangle_tolerance = kArithmeticTolerance * std::max({crank, projected, in_plane_radius});
    if (in_plane_radius + triangle_tolerance < difference || in_plane_radius > sum + triangle_tolerance) {
        return RotaryFailure(TransmissionStatus::Unreachable);
    }

    // Half-angle factors retain a tiny nonzero opening when acos(cosine)
    // would round it to zero for a short rod attached to a long crank.
    const long double sine_square = std::max(0.0L, in_plane_radius - difference) * (in_plane_radius + difference);
    const long double cosine_square = std::max(0.0L, sum - in_plane_radius) * (sum + in_plane_radius);
    if (!std::isfinite(sine_square) || !std::isfinite(cosine_square) ||
        (sine_square == 0.0L && cosine_square == 0.0L)) {
        return RotaryFailure(TransmissionStatus::NumericalFailure);
    }
    const long double opening = 2.0L * std::atan2(std::sqrt(sine_square), std::sqrt(cosine_square));
    const long double center = std::atan2(b, a);
    const bool tangent = sine_square == 0.0L || cosine_square == 0.0L;
    const long double offset = std::abs(static_cast<long double>(params.motor_offset_rad)) <= Pi()
                                   ? static_cast<long double>(params.motor_offset_rad)
                                   : CanonicalPhase(params.motor_offset_rad);
    RotaryRodResult result;
    result.status = tangent ? TransmissionStatus::Tangent : TransmissionStatus::Success;
    result.count = tangent ? 1 : 2;
    for (std::size_t index = 0; index < result.count; ++index) {
        auto& solution = result.solutions[index];
        const long double phase = center + (index == 0 ? opening : -opening);
        solution.crank_phase_rad = NarrowPhase(phase);
        solution.motor_phase_rad = NarrowPhase(params.motor_sign * phase + offset);
        const long double gamma = solution.crank_phase_rad;
        const long double sine = std::sin(gamma);
        const long double cosine = std::cos(gamma);
        const ExtendedVector link = delta - crank * (u * cosine + v * sine);
        if (!RodCloses(link, rod)) {
            return RotaryFailure(TransmissionStatus::NumericalFailure);
        }
        const long double recovered_phase = CanonicalPhase(params.motor_sign *
                                                           (static_cast<long double>(solution.motor_phase_rad) - offset));
        const bool phase_changed = recovered_phase != gamma;
        const long double recovered_sine = phase_changed ? std::sin(recovered_phase) : sine;
        const long double recovered_cosine = phase_changed ? std::cos(recovered_phase) : cosine;
        if (phase_changed && !RodCloses(delta - crank * (u * recovered_cosine + v * recovered_sine), rod)) {
            return RotaryFailure(TransmissionStatus::NumericalFailure);
        }
        if (tangent) {
            continue;
        }
        const ExtendedVector tangent_direction = -u * sine + v * cosine;
        const long double denominator = crank * delta.dot(tangent_direction);
        if (denominator == 0.0L || !std::isfinite(denominator)) {
            return RotaryFailure(TransmissionStatus::NumericalFailure);
        }
        if (phase_changed) {
            const ExtendedVector recovered_tangent = -u * recovered_sine + v * recovered_cosine;
            const long double recovered_denominator = crank * delta.dot(recovered_tangent);
            if (recovered_denominator == 0.0L || (recovered_denominator > 0.0L) != (denominator > 0.0L)) {
                return RotaryFailure(TransmissionStatus::NumericalFailure);
            }
        }
        const ExtendedVector gradient = link / (scale * denominator);
        if (!NarrowGradient(gradient, solution.crank_gradient) ||
            !NarrowGradient(params.motor_sign * gradient, solution.motor_gradient)) {
            return RotaryFailure(TransmissionStatus::NumericalFailure);
        }
        solution.branch = denominator > 0.0L ? 1 : -1;
        solution.derivative_valid = true;
    }
    if (result.count == 2 &&
        (result.solutions[0].crank_phase_rad == result.solutions[1].crank_phase_rad ||
         result.solutions[0].motor_phase_rad == result.solutions[1].motor_phase_rad)) {
        return RotaryFailure(TransmissionStatus::NumericalFailure);
    }
    return result;
}

void PrismaticLengthParams::Validate() const {
    if (!anchor.allFinite() || !std::isfinite(minimum_length) || minimum_length < 0.0 ||
        std::isnan(maximum_length) || maximum_length < minimum_length) {
        throw std::invalid_argument("prismatic anchor and nonnegative minimum must be finite; maximum must bound minimum");
    }
}

PrismaticLengthResult SolvePrismaticLength(const PrismaticLengthParams& params, const Eigen::Vector3d& point) {
    params.Validate();
    PrismaticLengthResult result;
    if (!point.allFinite()) {
        return result;
    }
    const ExtendedVector displacement = point.cast<long double>() - params.anchor.cast<long double>();
    const long double length = Norm(displacement);
    if (!std::isfinite(length) || length > std::numeric_limits<double>::max()) {
        result.status = TransmissionStatus::NumericalFailure;
        return result;
    }
    if (length < params.minimum_length || length > params.maximum_length) {
        result.status = TransmissionStatus::Unreachable;
        return result;
    }
    result.length = static_cast<double>(length);
    if (length == 0.0L) {
        result.status = TransmissionStatus::Singular;
        return result;
    }
    if (!NarrowGradient(displacement / length, result.gradient)) {
        result.status = TransmissionStatus::NumericalFailure;
        return result;
    }
    result.derivative_valid = true;
    result.status = TransmissionStatus::Success;
    return result;
}

}  // namespace parallel_ankle
