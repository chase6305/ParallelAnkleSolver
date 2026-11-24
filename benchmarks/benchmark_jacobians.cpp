#include <algorithm>
#include <array>
#include <charconv>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>

#include "parallel_ankle/parallel_ankle.h"

namespace {

constexpr std::size_t kPoseCount = 64;
constexpr std::size_t kRepetitions = 5;
constexpr unsigned kDefaultIterations = 2000;
constexpr unsigned kMaximumIterations = 100000;
constexpr double kDifferenceStep = 1e-6;
constexpr double kMaximumFrobeniusError = 1e-4;

using Poses = std::array<Eigen::Vector2d, kPoseCount>;
using Samples = std::array<double, kRepetitions>;

// Every timed result contributes to an observable value, consumed after timing.
volatile double result_sink = 0.0;

parallel_ankle::AnkleParams LegacyLeftParams() {
    // Same geometry and motor convention as configs/rotary_ankle_reference.yaml.
    parallel_ankle::AnkleParams params;
    params.r_roll_center = {0.0, 0.0, 0.0};
    params.r_pitch_center = {0.0, 0.0, 10.0};
    params.r_a1 = {-46.0, 0.0, 257.94};
    params.r_c1 = {-46.0, -59.0, 0.0};
    params.l_bar1 = 60.0;
    params.l_rod1 = 268.354;
    params.alpha_1_rad = 0.1745;
    params.is_left_motor1 = true;
    params.r_a2 = {-46.0, 0.0, 177.94};
    params.r_c2 = {-46.0, 59.0, 0.0};
    params.l_bar2 = 60.0;
    params.l_rod2 = 188.355;
    params.alpha_2_rad = -0.1745;
    params.is_left_motor2 = false;
    return params;
}

Poses MakePoses() {
    Poses poses;
    for (std::size_t pitch_index = 0; pitch_index < 8; ++pitch_index) {
        for (std::size_t roll_index = 0; roll_index < 8; ++roll_index) {
            poses[8 * pitch_index + roll_index] = {
                -0.2 + 0.5 * static_cast<double>(pitch_index) / 7.0,
                -0.1 + 0.2 * static_cast<double>(roll_index) / 7.0};
        }
    }
    return poses;
}

double ValidatePoses(const parallel_ankle::PMAnkle& ankle, const Poses& poses) {
    double maximum_error = 0.0;
    for (std::size_t index = 0; index < poses.size(); ++index) {
        const auto& pose = poses[index];
        if (!ankle.GetAnkleMotorAngles(pose(0), pose(1)).success()) {
            throw std::runtime_error("invalid IK at benchmark pose " + std::to_string(index));
        }
        const Eigen::Matrix2d analytic = ankle.ComputeJacobian(pose(0), pose(1));
        const Eigen::Matrix2d numerical = ankle.ComputeNumericalJacobian(pose(0), pose(1), kDifferenceStep);
        if (!analytic.allFinite() || !numerical.allFinite()) {
            throw std::runtime_error("nonfinite Jacobian at benchmark pose " + std::to_string(index));
        }
        const double error = (analytic - numerical).norm();
        if (!std::isfinite(error) || error > kMaximumFrobeniusError) {
            throw std::runtime_error("Jacobian comparison exceeded its error bound at benchmark pose " +
                                     std::to_string(index));
        }
        maximum_error = std::max(maximum_error, error);
    }
    return maximum_error;
}

template <bool Analytic>
double Measure(const parallel_ankle::PMAnkle& ankle, const Poses& poses, unsigned iterations) {
    double checksum = 0.0;
    const auto begin = std::chrono::steady_clock::now();
    for (unsigned iteration = 0; iteration < iterations; ++iteration) {
        const auto& pose = poses[iteration % poses.size()];
        // Volatile input reads prevent constant-pose folding under aggressive
        // whole-program optimization. Both methods include this same overhead.
        volatile double pitch = pose(0);
        volatile double roll = pose(1);
        Eigen::Matrix2d jacobian;
        if constexpr (Analytic) {
            jacobian = ankle.ComputeJacobian(pitch, roll);
        } else {
            jacobian = ankle.ComputeNumericalJacobian(pitch, roll, kDifferenceStep);
        }
        checksum += jacobian(0, 0) + 0.5 * jacobian(0, 1) + 0.25 * jacobian(1, 0) + 0.125 * jacobian(1, 1);
    }
    const auto end = std::chrono::steady_clock::now();
    if (!std::isfinite(checksum)) {
        throw std::runtime_error("nonfinite benchmark checksum");
    }
    result_sink = result_sink + checksum;
    return std::chrono::duration<double, std::nano>(end - begin).count() / static_cast<double>(iterations);
}

struct Summary {
    double median;
    double minimum;
    double maximum;
};

Summary Summarize(Samples samples) {
    std::sort(samples.begin(), samples.end());
    return {samples[kRepetitions / 2], samples.front(), samples.back()};
}

unsigned ParseIterations(std::string_view argument) {
    unsigned iterations = 0;
    const auto result = std::from_chars(argument.data(), argument.data() + argument.size(), iterations);
    if (result.ec != std::errc{} || result.ptr != argument.data() + argument.size() || iterations == 0 ||
        iterations > kMaximumIterations) {
        throw std::invalid_argument("iterations must be an integer from 1 to " + std::to_string(kMaximumIterations));
    }
    return iterations;
}

void PrintSummary(const char* name, const Summary& summary) {
    std::cout << name << ": median " << summary.median << " ns/call, min " << summary.minimum << ", max "
              << summary.maximum << '\n';
}

void PrintUsage(const char* executable) {
    std::cout << "usage: " << executable << " [iterations]\n"
              << "default iterations per method/repetition: " << kDefaultIterations << ", maximum: "
              << kMaximumIterations << '\n';
}

}  // namespace

int main(int argc, char** argv) {
    if (argc == 2 && std::string_view(argv[1]) == "--help") {
        PrintUsage(argv[0]);
        return 0;
    }
    if (argc > 2) {
        PrintUsage(argv[0]);
        return 2;
    }
    unsigned iterations = kDefaultIterations;
    try {
        if (argc == 2) {
            iterations = ParseIterations(argv[1]);
        }
    } catch (const std::invalid_argument& error) {
        std::cerr << error.what() << '\n';
        return 2;
    }

    try {
        const parallel_ankle::PMAnkle ankle(true, LegacyLeftParams());
        const Poses poses = MakePoses();
        // Validate all 64 poses and warm both methods before starting timers.
        const double maximum_error = ValidatePoses(ankle, poses);
        Samples analytic_samples{};
        Samples numerical_samples{};
        for (std::size_t repetition = 0; repetition < kRepetitions; ++repetition) {
            if (repetition % 2 == 0) {
                analytic_samples[repetition] = Measure<true>(ankle, poses, iterations);
                numerical_samples[repetition] = Measure<false>(ankle, poses, iterations);
            } else {
                numerical_samples[repetition] = Measure<false>(ankle, poses, iterations);
                analytic_samples[repetition] = Measure<true>(ankle, poses, iterations);
            }
        }

        const Summary analytic = Summarize(analytic_samples);
        const Summary numerical = Summarize(numerical_samples);
        std::cout << std::setprecision(10)
                  << "geometry: rotary reference left (mm, rad), " << poses.size() << " deterministic poses\n"
                  << "pose range: pitch [-0.2, 0.3], roll [-0.1, 0.1] rad\n"
                  << "validated maximum Frobenius error: " << maximum_error << " (bound "
                  << kMaximumFrobeniusError << ", finite-difference step " << kDifferenceStep << ")\n"
                  << "calls per method/repetition: " << iterations << ", repetitions: " << kRepetitions
                  << ", alternating method order\n";
        PrintSummary("analytic", analytic);
        PrintSummary("numerical", numerical);
        if (analytic.median > 0.0) {
            std::cout << "measured numerical/analytic median ratio: " << numerical.median / analytic.median << '\n';
        }
        std::cout << "consumed checksum: " << result_sink << '\n';
    } catch (const std::exception& error) {
        std::cerr << "benchmark failed: " << error.what() << '\n';
        return 1;
    }
    return 0;
}
