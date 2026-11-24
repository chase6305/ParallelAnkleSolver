#include <Eigen/Core>
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

#include "parallel_ankle/config.h"
#include "parallel_ankle/parallel_ankle.h"

namespace {
using Extended = Eigen::Matrix<long double, 3, 1>;
constexpr double kPitch = -1.2;
constexpr double kRoll = 0.18000000000000005;
constexpr long double kPhysicalTolerance = 512.0L * std::numeric_limits<double>::epsilon();

Extended RotateX(const Extended& point, long double angle) {
    return {point.x(), std::cos(angle) * point.y() - std::sin(angle) * point.z(),
            std::sin(angle) * point.y() + std::cos(angle) * point.z()};
}
Extended RotateY(const Extended& point, long double angle) {
    return {std::cos(angle) * point.x() + std::sin(angle) * point.z(), point.y(),
            -std::sin(angle) * point.x() + std::cos(angle) * point.z()};
}

struct CheckedPose {
    Eigen::Vector2d angles = Eigen::Vector2d::Constant(std::numeric_limits<double>::quiet_NaN());
    bool valid = false;
    double rod_relative_error = std::numeric_limits<double>::quiet_NaN();
};

CheckedPose CheckPose(const parallel_ankle::PMAnkle& ankle, double pitch, double roll) {
    CheckedPose result;
    const auto ik = ankle.GetAnkleMotorAngles(pitch, roll);
    if (!ik.success()) {
        return result;
    }
    result.angles = ik.motor_angles();
    result.rod_relative_error = 0.0;
    const auto& p = ankle.parameters();
    for (int motor = 0; motor < 2; ++motor) {
        const auto& a = motor == 0 ? p.r_a1 : p.r_a2;
        const auto& c = motor == 0 ? p.r_c1 : p.r_c2;
        const double bar = motor == 0 ? p.l_bar1 : p.l_bar2;
        const double rod = motor == 0 ? p.l_rod1 : p.l_rod2;
        const double offset = motor == 0 ? p.alpha_1_rad : p.alpha_2_rad;
        const bool left_motor = motor == 0 ? p.is_left_motor1 : p.is_left_motor2;
        const Extended delta = RotateY(RotateX(c.cast<long double>() - p.r_roll_center.cast<long double>(), roll) +
                                          p.r_roll_center.cast<long double>() - p.r_pitch_center.cast<long double>(),
                                      pitch) +
                               p.r_pitch_center.cast<long double>() - a.cast<long double>();
        const long double gamma = (static_cast<long double>(result.angles(motor)) - offset) *
                                  (left_motor ? 1.0L : -1.0L);
        const long double side = c.y() > a.y() ? 1.0L : -1.0L;
        const Extended crank(0.0L, side * bar * std::cos(gamma), bar * std::sin(gamma));
        const Extended link = delta - crank;
        const long double relative = std::abs(std::hypot(link.x(), link.y(), link.z()) - rod) / rod;
        if (!std::isfinite(relative) || relative > kPhysicalTolerance) {
            return result;
        }
        result.rod_relative_error = std::max(result.rod_relative_error, static_cast<double>(relative));
    }
    result.valid = true;
    return result;
}

void VectorJson(std::ostream& file, const char* name, const Eigen::Vector3d& vector, bool comma = true) {
    file << '"' << name << "\":[" << vector.x() << ',' << vector.y() << ',' << vector.z() << ']';
    if (comma) file << ',';
    file << '\n';
}

void WriteMetadata(const std::string& prefix, const parallel_ankle::AnkleParams& p, const CheckedPose& base) {
    std::ofstream file(prefix + ".meta.json");
    file.exceptions(std::ios::badbit | std::ios::failbit);
    file << std::setprecision(std::numeric_limits<double>::max_digits10)
         << "{\n\"profile\":\"rotary_v02\",\"side\":\"right\",\n"
         << "\"pitch\":" << kPitch << ",\"roll\":" << kRoll << ",\n"
         << "\"theta\":[" << base.angles(0) << ',' << base.angles(1) << "],\n";
    VectorJson(file, "r_roll_center", p.r_roll_center);
    VectorJson(file, "r_pitch_center", p.r_pitch_center);
    VectorJson(file, "r_a1", p.r_a1);
    VectorJson(file, "r_c1", p.r_c1);
    VectorJson(file, "r_a2", p.r_a2);
    VectorJson(file, "r_c2", p.r_c2);
    file << "\"l_bar1\":" << p.l_bar1 << ",\"l_rod1\":" << p.l_rod1 << ",\n"
         << "\"l_bar2\":" << p.l_bar2 << ",\"l_rod2\":" << p.l_rod2 << ",\n"
         << "\"alpha_1_rad\":" << p.alpha_1_rad << ",\"alpha_2_rad\":" << p.alpha_2_rad << ",\n"
         << "\"is_left_motor1\":" << (p.is_left_motor1 ? "true" : "false") << ",\n"
         << "\"is_left_motor2\":" << (p.is_left_motor2 ? "true" : "false") << "\n}\n";
}

void WriteRow(std::ostream& file,
              const char* method,
              double eps,
              bool valid,
              double physical_error,
              const Eigen::Matrix2d& matrix,
              const Eigen::Matrix2d& analytic) {
    const double pp = kPitch + eps;
    const double pm = kPitch - eps;
    const double rp = kRoll + eps;
    const double rm = kRoll - eps;
    const bool central = std::string(method) == "central";
    const double dp = central ? pp - pm : pp - kPitch;
    const double dr = central ? rp - rm : rp - kRoll;
    const double reference_norm = analytic.norm();
    const double relative = valid ? (matrix - analytic).norm() / (reference_norm > 0.0 ? reference_norm : 1.0)
                                  : std::numeric_limits<double>::quiet_NaN();
    file << method << ',' << eps << ',' << kPitch << ',' << kRoll << ',' << pp << ',' << pm << ',' << rp << ',' << rm
         << ',' << dp << ',' << dr << ',' << (valid ? "ok" : "invalid_physical_stencil") << ',' << physical_error;
    for (int row = 0; row < 2; ++row) {
        for (int column = 0; column < 2; ++column) {
            file << ',' << matrix(row, column);
        }
    }
    file << ',' << relative << '\n';
}
}  // namespace

int main(int argc, char** argv) {
    if (argc != 3) {
        std::cerr << "usage: jacobian_convergence_probe <config.yaml> <output-prefix>\n";
        return 2;
    }
    try {
        const auto config = parallel_ankle::LoadAnkleConfig(argv[1], "rotary_v02");
        const parallel_ankle::PMAnkle ankle(false, config.right);
        const CheckedPose base = CheckPose(ankle, kPitch, kRoll);
        if (!base.valid) throw std::runtime_error("base pose does not satisfy the physical rod constraints");
        const Eigen::Matrix2d analytic = ankle.ComputeJacobian(kPitch, kRoll);
        WriteMetadata(argv[2], config.right, base);
        std::ofstream csv(std::string(argv[2]) + ".csv");
        csv.exceptions(std::ios::badbit | std::ios::failbit);
        csv << std::setprecision(std::numeric_limits<double>::max_digits10)
            << "method,eps,pitch,roll,pitch_plus,pitch_minus,roll_plus,roll_minus,actual_pitch_step,actual_roll_step,"
               "status,max_rod_relative_error,J00,J01,J10,J11,relative_error_vs_analytic\n";
        WriteRow(csv, "analytic", 0.0, true, base.rod_relative_error, analytic, analytic);
        for (double eps : {1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9}) {
            const auto pp = CheckPose(ankle, kPitch + eps, kRoll);
            const auto pm = CheckPose(ankle, kPitch - eps, kRoll);
            const auto rp = CheckPose(ankle, kPitch, kRoll + eps);
            const auto rm = CheckPose(ankle, kPitch, kRoll - eps);
            const bool forward_valid = pp.valid && rp.valid;
            Eigen::Matrix2d forward = Eigen::Matrix2d::Constant(std::numeric_limits<double>::quiet_NaN());
            if (forward_valid) forward = ankle.ComputeNumericalJacobian(kPitch, kRoll, eps);
            WriteRow(csv, "forward", eps, forward_valid,
                     std::max({base.rod_relative_error, pp.rod_relative_error, rp.rod_relative_error}), forward, analytic);
            const bool central_valid = pp.valid && pm.valid && rp.valid && rm.valid;
            Eigen::Matrix2d central = Eigen::Matrix2d::Constant(std::numeric_limits<double>::quiet_NaN());
            if (central_valid) {
                central.col(0) = (pp.angles - pm.angles) / ((kPitch + eps) - (kPitch - eps));
                central.col(1) = (rp.angles - rm.angles) / ((kRoll + eps) - (kRoll - eps));
            }
            WriteRow(csv, "central", eps, central_valid,
                     std::max({base.rod_relative_error, pp.rod_relative_error, pm.rod_relative_error,
                               rp.rod_relative_error, rm.rod_relative_error}), central, analytic);
        }
        std::cout << "Wrote " << argv[2] << ".csv and " << argv[2] << ".meta.json\n";
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
    return 0;
}
