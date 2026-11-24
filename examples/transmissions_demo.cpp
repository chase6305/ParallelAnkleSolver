#include <iomanip>
#include <iostream>

#include "parallel_ankle/transmissions.h"

int main() {
    const parallel_ankle::RotaryRodParams params;
    const auto rotary = parallel_ankle::SolveRotaryRod(params, {1.0, 0.0, 0.0});
    std::cout << std::setprecision(10) << "single rotary-crank/rod transmission: "
              << parallel_ankle::ToString(rotary.status) << '\n';
    if (!rotary.success()) {
        return 1;
    }
    for (std::size_t index = 0; index < rotary.count; ++index) {
        const auto& solution = rotary.solutions[index];
        std::cout << "  branch " << solution.branch << ": crank phase " << solution.crank_phase_rad
                  << " rad, motor phase " << solution.motor_phase_rad << " rad\n"
                  << "  motor gradient [rad/length]: " << solution.motor_gradient.transpose() << '\n';
    }
    const auto linear = parallel_ankle::SolvePrismaticLength({}, {3.0, 4.0, 0.0});
    std::cout << "single telescoping transmission: " << parallel_ankle::ToString(linear.status)
              << ", length " << linear.length << ", gradient " << linear.gradient.transpose() << '\n';
    return linear.success() ? 0 : 1;
}
