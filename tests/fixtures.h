#pragma once

#include <gtest/gtest.h>

#include <Eigen/Geometry>
#include <cmath>

#include "parallel_ankle/parallel_ankle.h"

namespace parallel_ankle::test {

// Copied from the original test's MakeAnkleParams, including its motor signs.
inline AnkleParams LegacyParams(bool is_left) {
    AnkleParams params;
    params.r_pitch_center = {0.0, 0.0, 10.0};
    params.r_a1 = {-46.0, 0.0, 257.94};
    params.r_a2 = {-46.0, 0.0, 177.94};
    params.r_c1 = {-46.0, is_left ? -59.0 : 59.0, 0.0};
    params.r_c2 = {-46.0, is_left ? 59.0 : -59.0, 0.0};
    params.l_bar1 = params.l_bar2 = 60.0;
    params.l_rod1 = 268.354;
    params.l_rod2 = 188.355;
    params.alpha_1_rad = is_left ? 0.1745 : -0.1745;
    params.alpha_2_rad = -params.alpha_1_rad;
    params.is_left_motor1 = is_left;
    params.is_left_motor2 = !is_left;
    return params;
}

inline void ExpectRodClosure(const AnkleParams& params, const Eigen::Vector2d& pose, const Eigen::Vector2d& motors) {
    const Eigen::Matrix3d pitch = Eigen::AngleAxisd(pose(0), Eigen::Vector3d::UnitY()).toRotationMatrix();
    const Eigen::Matrix3d roll = Eigen::AngleAxisd(pose(1), Eigen::Vector3d::UnitX()).toRotationMatrix();
    for (int index = 0; index < 2; ++index) {
        const Eigen::Vector3d& a = index == 0 ? params.r_a1 : params.r_a2;
        const Eigen::Vector3d& c0 = index == 0 ? params.r_c1 : params.r_c2;
        const double bar = index == 0 ? params.l_bar1 : params.l_bar2;
        const double rod = index == 0 ? params.l_rod1 : params.l_rod2;
        const double offset = index == 0 ? params.alpha_1_rad : params.alpha_2_rad;
        const bool left_motor = index == 0 ? params.is_left_motor1 : params.is_left_motor2;
        const double crank_angle = (motors(index) - offset) * (left_motor ? 1.0 : -1.0);
        const double side = c0.y() >= a.y() ? 1.0 : -1.0;
        const Eigen::Vector3d b = a + Eigen::Vector3d(0.0, side * bar * std::cos(crank_angle), bar * std::sin(crank_angle));
        const Eigen::Vector3d c = pitch * (roll * (c0 - params.r_roll_center) + params.r_roll_center -
                                         params.r_pitch_center) + params.r_pitch_center;
        EXPECT_NEAR((b - a).norm(), bar, 1e-9);
        EXPECT_NEAR((c - b).norm(), rod, 1e-7) << "rod " << index + 1 << ", pose " << pose.transpose();
    }
}

}  // namespace parallel_ankle::test
