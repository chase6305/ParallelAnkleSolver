#include <iomanip>
#include <iostream>

#include "parallel_ankle/constraints.h"

int main() {
    using namespace parallel_ankle::constraints;
    // Two rods of length sqrt(2) connect a moving point to fixed centers
    // (0,0) and (span,0). The upper and lower circle intersections are two
    // assembly modes; a warm start and bounds select the upper local branch.
    constexpr double span = 2.0;
    EqualityProblem problem;
    problem.variables = problem.equations = 2;
    problem.residual = [](const Eigen::VectorXd& x, Eigen::VectorXd& residual) {
        residual.resize(2);
        residual << x.squaredNorm() - 2.0, (x(0) - span) * (x(0) - span) + x(1) * x(1) - 2.0;
        return true;
    };
    problem.jacobian = [](const Eigen::VectorXd& x, Eigen::MatrixXd& jacobian) {
        jacobian.resize(2, 2);
        jacobian << 2.0 * x(0), 2.0 * x(1), 2.0 * (x(0) - span), 2.0 * x(1);
        return true;
    };
    SolveOptions options;
    options.lower_bounds = Eigen::Vector2d(-std::numeric_limits<double>::infinity(), 0.0);
    const auto solution = SolveEqualities(problem, Eigen::Vector2d(0.8, 0.9), options);
    if (!solution.success()) {
        std::cerr << ToString(solution.status) << ": " << solution.message << '\n';
        return 1;
    }

    Eigen::MatrixXd a;
    problem.jacobian(solution.x, a);
    Eigen::MatrixXd b(2, 1);
    b << 0.0, -2.0 * (solution.x(0) - span);  // dF/dspan at a closed point
    const auto derivative = DifferentiateImplicitly(a, b);
    if (!derivative.success()) {
        std::cerr << ToString(derivative.status) << ": " << derivative.message << '\n';
        return 1;
    }
    std::cout << std::setprecision(12)
              << "local upper assembly [x,y]: " << solution.x.transpose() << '\n'
              << "scaled closure residual: " << solution.scaled_residual_norm << '\n'
              << "Jacobian rank: " << solution.rank << '\n'
              << "d[x,y]/dspan: " << derivative.derivative.col(0).transpose() << '\n';
    return 0;
}
