#include <cmath>
#include <iomanip>
#include <iostream>

#include "parallel_ankle/second_order.h"

int main() {
    try {
        // Small illustrative rotary-rod geometry; units are consistent mm.
        parallel_ankle::AnkleParams params;
        params.r_pitch_center = {0.0, 0.0, 2.0};
        params.r_a1 = {0.0, 0.0, 10.0};
        params.r_a2 = {0.0, 0.0, 8.0};
        params.r_c1 = {-2.0, 4.0, 0.0};
        params.r_c2 = {-2.0, -4.0, 0.0};
        params.l_bar1 = params.l_bar2 = 4.0;
        params.l_rod1 = std::sqrt(104.0);
        params.l_rod2 = std::sqrt(68.0);
        const parallel_ankle::PMAnkle model(true, params);
        const Eigen::Vector2d pose(0.2, -0.1);
        const Eigen::Vector2d motor_preload(2.0, -1.0);
        const Eigen::Matrix2d task_stiffness = Eigen::Vector2d(100.0, 80.0).asDiagonal();
        const auto result = parallel_ankle::MapLoadedStaticTaskStiffness(model, pose(0), pose(1), task_stiffness,
                                                                       motor_preload);
        const Eigen::Matrix2d jacobian = model.ComputeJacobian(pose(0), pose(1));
        const Eigen::Matrix2d recovered = jacobian.transpose() * result.motor_stiffness * jacobian +
                                          result.geometric_stiffness;
        std::cout << std::setprecision(10)
                  << "static restoring-torque convention: tau_motor = -gradient(U_motor)\n"
                  << "task stiffness:\n" << task_stiffness << '\n'
                  << "task geometric contribution (-sum tau_motor_i H_i):\n" << result.geometric_stiffness << '\n'
                  << "motor stiffness:\n" << result.motor_stiffness << '\n'
                  << "round-trip error: " << (recovered - task_stiffness).norm() << '\n'
                  << "task/motor PSD diagnostics: " << std::boolalpha << result.task_positive_semidefinite << "/"
                  << result.motor_positive_semidefinite << '\n';
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
    return 0;
}
