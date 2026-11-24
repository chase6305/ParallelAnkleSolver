#include <iomanip>
#include <iostream>

#include "parallel_ankle/five_bar.h"

int main() {
    parallel_ankle::FiveBarParams parameters;
    parameters.pivots = {Eigen::Vector2d(-0.4, 0.1), Eigen::Vector2d(0.5, -0.1)};
    parameters.proximal = {0.8, 0.7};
    parameters.distal = {1.0, 0.9};
    const parallel_ankle::PlanarFiveBar mechanism(parameters);
    const Eigen::Vector2d target(0.1, 1.0);
    const auto inverse = mechanism.InverseAll(target);
    std::cout << std::setprecision(10) << "This example chooses metres for all geometry coordinates.\n"
              << "target [x,y] [m]: " << target.transpose() << '\n'
              << "IK: " << parallel_ankle::ToString(inverse.status) << ", " << inverse.count << " working modes\n";
    if (!inverse.success()) {
        return 1;
    }
    for (std::size_t index = 0; index < inverse.count; ++index) {
        const auto& solution = inverse.solutions[index];
        const auto differential = mechanism.Analyze(target, solution.motor_angles);
        if (!differential.success()) {
            std::cerr << parallel_ankle::ToString(differential.status) << '\n';
            return 1;
        }
        std::cout << "mode " << index << ": motors [rad] " << solution.motor_angles.transpose()
                  << ", working [" << solution.branch.working[0] << ',' << solution.branch.working[1]
                  << "], assembly " << solution.branch.assembly << '\n'
                  << "  differential: " << parallel_ankle::ToString(differential.status) << '\n'
                  << "  dtheta/d[x,y] [rad/m]:\n" << differential.motor_jacobian << '\n';
        const auto forward = mechanism.ForwardAll(solution.motor_angles);
        if (!forward.success()) {
            std::cerr << parallel_ankle::ToString(forward.status) << '\n';
            return 1;
        }
        std::cout << "  FK: " << parallel_ankle::ToString(forward.status) << ", " << forward.count << " solutions\n";
        for (std::size_t assembly = 0; assembly < forward.count; ++assembly) {
            std::cout << "    " << forward.solutions[assembly].position.transpose() << '\n';
        }
    }
    return 0;
}
