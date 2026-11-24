#pragma once

#include <Eigen/Core>

namespace parallel_ankle {

// Coordinates and lengths must share a unit; the example profiles use millimetres.
// Joint order is [pitch, roll]. Motor numbers follow wiring, not physical side.
struct AnkleParams {
    Eigen::Vector3d r_roll_center = Eigen::Vector3d::Zero();
    Eigen::Vector3d r_pitch_center = Eigen::Vector3d::Zero();

    bool is_left_motor1 = true;
    Eigen::Vector3d r_a1 = Eigen::Vector3d::Zero();
    Eigen::Vector3d r_c1 = Eigen::Vector3d::Zero();
    double l_bar1 = 0.0;
    double l_rod1 = 0.0;
    double alpha_1_rad = 0.0;

    bool is_left_motor2 = false;
    Eigen::Vector3d r_a2 = Eigen::Vector3d::Zero();
    Eigen::Vector3d r_c2 = Eigen::Vector3d::Zero();
    double l_bar2 = 0.0;
    double l_rod2 = 0.0;
    double alpha_2_rad = 0.0;

    // Throws std::invalid_argument for nonfinite data or nonpositive lengths.
    void Validate() const;
};

}  // namespace parallel_ankle
