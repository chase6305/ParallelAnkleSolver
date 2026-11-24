#include "parallel_ankle/five_bar.h"

#include <algorithm>
#include <cmath>
#include <functional>
#include <stdexcept>

namespace parallel_ankle {
namespace {

// Extended intermediates prevent overflow in squared lengths and pivot differences
// while the public coordinates, angles and matrices remain ordinary doubles.
struct Point {
    long double x;
    long double y;
};

Point operator+(Point a, Point b) { return {a.x + b.x, a.y + b.y}; }
Point operator-(Point a, Point b) { return {a.x - b.x, a.y - b.y}; }
Point operator*(Point a, long double scale) { return {a.x * scale, a.y * scale}; }
Point operator/(Point a, long double scale) { return {a.x / scale, a.y / scale}; }
long double Dot(Point a, Point b) { return a.x * b.x + a.y * b.y; }
long double Cross(Point a, Point b) { return a.x * b.y - a.y * b.x; }
long double Norm(Point a) { return std::hypot(a.x, a.y); }
Point Extended(const Eigen::Vector2d& point) { return {point.x(), point.y()}; }

constexpr long double kRoundoff = 64.0L * std::numeric_limits<long double>::epsilon();
constexpr long double kBranchRoundoff = 64.0L * std::numeric_limits<double>::epsilon();

bool Narrow(Point point, Eigen::Vector2d& result) {
    const long double bound = std::numeric_limits<double>::max();
    if (!std::isfinite(point.x) || !std::isfinite(point.y) || std::abs(point.x) > bound || std::abs(point.y) > bound) {
        return false;
    }
    result = {static_cast<double>(point.x), static_cast<double>(point.y)};
    return result.allFinite();
}

Point Crank(long double angle, long double length) {
    return Point{std::cos(angle), std::sin(angle)} * length;
}

std::array<Point, 2> RodVectors(Point position,
                               const std::array<Point, 2>& pivots,
                               const std::array<Point, 2>& cranks) {
    return {(position - pivots[0]) - cranks[0], (position - pivots[1]) - cranks[1]};
}

int CrossSign(Point first, Point second, long double second_length) {
    const long double scale = Norm(first) * second_length;
    const long double value = Cross(first, second);
    if (std::abs(value) <= kBranchRoundoff * scale) {
        return 0;
    }
    return value > 0.0L ? 1 : -1;
}

FiveBarBranch Branch(Point elbow_delta,
                     const std::array<Point, 2>& cranks,
                     const std::array<Point, 2>& rods) {
    FiveBarBranch branch;
    const long double first_rod_length = Norm(rods[0]);
    branch.assembly = CrossSign(elbow_delta, rods[0], first_rod_length);
    branch.working[0] = CrossSign(cranks[0], rods[0], first_rod_length);
    branch.working[1] = CrossSign(cranks[1], rods[1], Norm(rods[1]));
    return branch;
}

// Comparisons allow only arithmetic roundoff in the extended intermediates.
// A target outside the triangle inequalities by an ordinary double step fails.
bool TriangleReachable(long double distance, long double first, long double second) {
    const long double scale = std::max({distance, first, second});
    const long double tolerance = kRoundoff * scale;
    return distance <= first + second + tolerance && distance + tolerance >= std::abs(first - second);
}

// Four times the triangle area, divided by scale^2. The sorted, factored
// Heron expression preserves a tiny side when the other two sides are large.
bool ScaledAreaFour(long double first, long double second, long double third, long double& area_four) {
    std::array<long double, 3> sides{first, second, third};
    std::sort(sides.begin(), sides.end(), std::greater<long double>());
    const long double a = sides[0];
    const long double b = sides[1];
    const long double c = sides[2];
    std::array<long double, 4> factors{
        (a + (b + c)) / a, (c - (a - b)) / a, (c + (a - b)) / a, (a + (b - c)) / a};
    long double product = 1.0L;
    for (long double factor : factors) {
        if (!std::isfinite(factor) || factor < -kRoundoff) {
            return false;
        }
        product *= std::max(0.0L, factor);
    }
    area_four = std::sqrt(product);
    return std::isfinite(area_four);
}

bool RodClosed(long double actual, long double required_length) {
    // A long crank cannot excuse an error comparable to a short distal rod.
    const long double tolerance = 256.0L * std::numeric_limits<double>::epsilon() * required_length;
    return actual > 0.0L && std::isfinite(actual) && std::abs(actual - required_length) <= tolerance;
}

bool RodsClosed(const std::array<Point, 2>& rods, const std::array<double, 2>& required_lengths) {
    return RodClosed(Norm(rods[0]), required_lengths[0]) && RodClosed(Norm(rods[1]), required_lengths[1]);
}

double CanonicalAngle(long double angle) {
    const long double pi = std::acos(-1.0L);
    long double wrapped = std::remainder(angle, 2.0L * pi);
    if (wrapped >= pi) {
        wrapped -= 2.0L * pi;
    }
    return static_cast<double>(wrapped);
}

struct LegSolution {
    double angle = 0.0;
    Point crank{};
    Point rod{};
};

struct LegInverse {
    FiveBarStatus status = FiveBarStatus::Unreachable;
    std::array<LegSolution, 2> solutions;
    std::size_t count = 0;
};

LegInverse InverseLeg(Point target, long double proximal, long double distal) {
    LegInverse result;
    const long double distance = Norm(target);
    if (distance == 0.0L) {
        result.status = proximal == distal ? FiveBarStatus::Singular : FiveBarStatus::Unreachable;
        return result;
    }
    if (!TriangleReachable(distance, proximal, distal)) {
        return result;
    }
    const long double scale = std::max({distance, proximal, distal});
    const long double p = proximal / scale;
    const long double d = distance / scale;
    const long double r = distal / scale;
    long double area_four;
    if (!ScaledAreaFour(proximal, distance, distal, area_four)) {
        result.status = FiveBarStatus::NumericalFailure;
        return result;
    }
    const long double cosine_numerator = (p - r) * (p + r) + d * d;
    const long double center = std::atan2(target.y, target.x);
    const long double opening = std::atan2(area_four, cosine_numerator);
    result.solutions[0].angle = CanonicalAngle(center + opening);
    result.solutions[1].angle = CanonicalAngle(center - opening);
    result.count = area_four == 0.0L || result.solutions[0].angle == result.solutions[1].angle ? 1 : 2;
    for (std::size_t index = 0; index < result.count; ++index) {
        auto& solution = result.solutions[index];
        solution.crank = Crank(solution.angle, proximal);
        solution.rod = target - solution.crank;
        if (!RodClosed(Norm(solution.rod), distal)) {
            result.status = FiveBarStatus::NumericalFailure;
            result.count = 0;
            return result;
        }
    }
    result.status = FiveBarStatus::Success;
    return result;
}

}  // namespace

const char* ToString(FiveBarStatus status) {
    switch (status) {
        case FiveBarStatus::Success:
            return "success";
        case FiveBarStatus::InvalidInput:
            return "invalid input";
        case FiveBarStatus::Unreachable:
            return "unreachable configuration";
        case FiveBarStatus::Singular:
            return "singular or nonunique configuration";
        case FiveBarStatus::NumericalFailure:
            return "result cannot be represented reliably";
    }
    return "unknown status";
}

void FiveBarParams::Validate() const {
    for (int leg = 0; leg < 2; ++leg) {
        if (!pivots[leg].allFinite() || !std::isfinite(proximal[leg]) || !std::isfinite(distal[leg]) ||
            proximal[leg] <= 0.0 || distal[leg] <= 0.0) {
            throw std::invalid_argument("five-bar pivots must be finite and all lengths positive and finite");
        }
    }
}

PlanarFiveBar::PlanarFiveBar(const FiveBarParams& params) : params_(params) { params_.Validate(); }

FiveBarForwardResult PlanarFiveBar::ForwardAll(const Eigen::Vector2d& motor_angles) const {
    FiveBarForwardResult result;
    if (!motor_angles.allFinite()) {
        return result;
    }
    const std::array<Point, 2> pivots{Extended(params_.pivots[0]), Extended(params_.pivots[1])};
    const std::array<Point, 2> cranks{Crank(motor_angles(0), params_.proximal[0]),
                                      Crank(motor_angles(1), params_.proximal[1])};
    const Point first = pivots[0] + cranks[0];
    const Point delta = (pivots[1] - pivots[0]) + (cranks[1] - cranks[0]);
    const long double distance = Norm(delta);
    const long double radius_first = params_.distal[0];
    const long double radius_second = params_.distal[1];
    if (distance == 0.0L) {
        result.status = radius_first == radius_second ? FiveBarStatus::Singular : FiveBarStatus::Unreachable;
        return result;
    }
    if (!TriangleReachable(distance, radius_first, radius_second)) {
        result.status = FiveBarStatus::Unreachable;
        return result;
    }
    const long double scale = std::max({distance, radius_first, radius_second});
    const long double d = distance / scale;
    const long double first_radius = radius_first / scale;
    const long double second_radius = radius_second / scale;
    const long double along = ((first_radius - second_radius) * (first_radius + second_radius) + d * d) / (2.0L * d);
    long double area_four;
    if (!ScaledAreaFour(radius_first, radius_second, distance, area_four)) {
        result.status = FiveBarStatus::NumericalFailure;
        return result;
    }
    const Point direction = delta / distance;
    const Point center = first + direction * (scale * along);
    const Point offset = Point{-direction.y, direction.x} * (scale * area_four / (2.0L * d));
    Eigen::Vector2d positive;
    Eigen::Vector2d negative;
    if (!Narrow(center + offset, positive) || !Narrow(center - offset, negative)) {
        result.status = FiveBarStatus::NumericalFailure;
        return result;
    }
    const auto positive_rods = RodVectors(Extended(positive), pivots, cranks);
    const auto negative_rods = RodVectors(Extended(negative), pivots, cranks);
    if (!RodsClosed(positive_rods, params_.distal) || !RodsClosed(negative_rods, params_.distal)) {
        result.status = FiveBarStatus::NumericalFailure;
        return result;
    }
    result.solutions[0].position = positive;
    result.solutions[0].branch = Branch(delta, cranks, positive_rods);
    result.count = 1;
    if (!(positive.array() == negative.array()).all()) {
        result.solutions[1].position = negative;
        result.solutions[1].branch = Branch(delta, cranks, negative_rods);
        result.count = 2;
    }
    result.status = FiveBarStatus::Success;
    return result;
}

FiveBarInverseResult PlanarFiveBar::InverseAll(const Eigen::Vector2d& position) const {
    FiveBarInverseResult result;
    if (!position.allFinite()) {
        return result;
    }
    const Point target = Extended(position);
    const Point first_pivot = Extended(params_.pivots[0]);
    const Point second_pivot = Extended(params_.pivots[1]);
    const LegInverse first = InverseLeg(target - first_pivot, params_.proximal[0], params_.distal[0]);
    const LegInverse second = InverseLeg(target - second_pivot, params_.proximal[1], params_.distal[1]);
    if (first.status == FiveBarStatus::Unreachable || second.status == FiveBarStatus::Unreachable) {
        result.status = FiveBarStatus::Unreachable;
        return result;
    }
    if (first.status == FiveBarStatus::NumericalFailure || second.status == FiveBarStatus::NumericalFailure) {
        result.status = FiveBarStatus::NumericalFailure;
        return result;
    }
    if (first.status == FiveBarStatus::Singular || second.status == FiveBarStatus::Singular) {
        result.status = FiveBarStatus::Singular;
        return result;
    }
    const Point pivot_delta = second_pivot - first_pivot;
    for (std::size_t left = 0; left < first.count; ++left) {
        for (std::size_t right = 0; right < second.count; ++right) {
            const auto& left_solution = first.solutions[left];
            const auto& right_solution = second.solutions[right];
            auto& solution = result.solutions[result.count++];
            solution.motor_angles = {left_solution.angle, right_solution.angle};
            const std::array<Point, 2> cranks{left_solution.crank, right_solution.crank};
            const Point elbow_delta = pivot_delta + (cranks[1] - cranks[0]);
            solution.branch = Branch(elbow_delta, cranks, {left_solution.rod, right_solution.rod});
        }
    }
    result.status = FiveBarStatus::Success;
    return result;
}

FiveBarDifferential PlanarFiveBar::Analyze(const Eigen::Vector2d& position,
                                         const Eigen::Vector2d& motor_angles,
                                         double singular_tolerance) const {
    FiveBarDifferential result;
    if (!position.allFinite() || !motor_angles.allFinite() || !std::isfinite(singular_tolerance) ||
        singular_tolerance <= 0.0 || singular_tolerance >= 1.0) {
        return result;
    }
    result.B.setZero();
    const std::array<Point, 2> pivots{Extended(params_.pivots[0]), Extended(params_.pivots[1])};
    const std::array<Point, 2> cranks{Crank(motor_angles(0), params_.proximal[0]),
                                      Crank(motor_angles(1), params_.proximal[1])};
    const auto rods = RodVectors(Extended(position), pivots, cranks);
    std::array<Point, 2> units;
    for (int leg = 0; leg < 2; ++leg) {
        const long double distance = Norm(rods[leg]);
        if (!RodClosed(distance, params_.distal[leg])) {
            result.status = FiveBarStatus::Unreachable;
            return result;
        }
        units[leg] = rods[leg] / distance;
        result.A.row(leg) << static_cast<double>(units[leg].x), static_cast<double>(units[leg].y);
        const long double diagonal = Dot(units[leg], Point{-cranks[leg].y, cranks[leg].x});
        result.B(leg, leg) = static_cast<double>(diagonal);
        result.actuator_singular = result.actuator_singular ||
                                   std::abs(diagonal) <= singular_tolerance * params_.proximal[leg];
    }
    const Point elbow_delta = (pivots[1] - pivots[0]) + (cranks[1] - cranks[0]);
    result.branch = Branch(elbow_delta, cranks, rods);
    result.assembly_singular = std::abs(Cross(units[0], units[1])) <= singular_tolerance;
    if (!result.A.allFinite() || !result.B.allFinite()) {
        result.status = FiveBarStatus::NumericalFailure;
        return result;
    }
    if (!result.actuator_singular) {
        for (int leg = 0; leg < 2; ++leg) {
            result.motor_jacobian.row(leg) = result.A.row(leg) / result.B(leg, leg);
        }
        if (!result.motor_jacobian.allFinite()) {
            result.status = FiveBarStatus::NumericalFailure;
            return result;
        }
    }
    result.status = result.actuator_singular || result.assembly_singular ? FiveBarStatus::Singular
                                                                       : FiveBarStatus::Success;
    return result;
}

}  // namespace parallel_ankle
