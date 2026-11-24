#include <cmath>

#include <parallel_ankle/parallel_ankle.h>
#include <parallel_ankle/constraints.h>
#include <parallel_ankle/five_bar.h>
#include <parallel_ankle/transmissions.h>
#include <parallel_ankle/second_order.h>

int main() {
    parallel_ankle::AnkleParams params;
    params.r_a1 = params.r_a2 = {0.0, 0.0, 2.0};
    params.r_c1 = {-0.5, -1.0, 0.0};
    params.r_c2 = {-0.5, 1.0, 0.0};
    params.l_bar1 = params.l_bar2 = 1.0;
    params.l_rod1 = params.l_rod2 = std::hypot(0.5, 2.0);
    const parallel_ankle::PMAnkle ankle(true, params);
    const auto ik = ankle.GetAnkleMotorAngles(0.0, 0.0);
    if (!ik.success() || !ik.motor_angles().allFinite()) {
        return 1;
    }
    const auto loaded = parallel_ankle::MapLoadedStaticTaskStiffness(
        ankle, 0.0, 0.0, Eigen::Matrix2d::Identity(), Eigen::Vector2d::Zero());
    if (!loaded.motor_stiffness.allFinite() || !loaded.motor_positive_semidefinite) {
        return 5;
    }
    const parallel_ankle::PlanarFiveBar five_bar(parallel_ankle::FiveBarParams{});
    if (!five_bar.InverseAll({0.0, 1.0}).success()) {
        return 2;
    }
    const auto transmission = parallel_ankle::SolveRotaryRod(parallel_ankle::RotaryRodParams{}, {1.0, 0.0, 0.0});
    if (!transmission.success() || transmission.count != 2) {
        return 4;
    }
    parallel_ankle::constraints::EqualityProblem problem;
    problem.variables = problem.equations = 1;
    problem.residual = [](const Eigen::VectorXd& x, Eigen::VectorXd& residual) {
        residual = Eigen::VectorXd::Constant(1, x(0) - 2.0);
        return true;
    };
    problem.jacobian = [](const Eigen::VectorXd&, Eigen::MatrixXd& jacobian) {
        jacobian = Eigen::MatrixXd::Ones(1, 1);
        return true;
    };
    const auto solved = parallel_ankle::constraints::SolveEqualities(problem, Eigen::VectorXd::Zero(1));
    return solved.success() && std::abs(solved.x(0) - 2.0) < 1e-10 ? 0 : 3;
}
