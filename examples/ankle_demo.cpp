#include <iomanip>
#include <iostream>
#include <stdexcept>

#include "parallel_ankle/config.h"
#include "parallel_ankle/parallel_ankle.h"

int main(int argc, char** argv) {
    if (argc < 2 || argc > 3) {
        std::cerr << "usage: ankle_demo <config.yaml> [profile]\n";
        return 2;
    }
    try {
        const auto config = parallel_ankle::LoadAnkleConfig(argv[1], argc == 3 ? argv[2] : "");
        std::cout << std::setprecision(10) << "profile: " << config.profile << '\n';
        for (bool is_left : {true, false}) {
            const parallel_ankle::PMAnkle ankle(is_left, is_left ? config.left : config.right);
            const auto ik = ankle.GetAnkleMotorAngles(0.1, 0.05);
            if (!ik.success()) {
                throw std::runtime_error("example pose is unreachable");
            }
            const auto fk = ankle.SolveForwardKinematics(ik.theta_1_rad, ik.theta_2_rad);
            if (!fk.success()) {
                throw std::runtime_error(parallel_ankle::ToString(fk.status));
            }
            std::cout << (is_left ? "left" : "right") << " ankle\n"
                      << "  motors [rad]: " << ik.motor_angles().transpose() << '\n'
                      << "  recovered [pitch, roll] [rad]: " << fk.pose.transpose() << '\n'
                      << "  FK residual [rad]: " << fk.residual_norm << '\n'
                      << "  Jacobian:\n" << ankle.ComputeJacobian(0.1, 0.05) << '\n';
        }
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
    return 0;
}
