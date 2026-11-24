#pragma once

#include <array>

#include <Eigen/Core>

#include "parallel_ankle/parameters.h"

namespace parallel_ankle::detail {

struct AnkleDerivatives {
    Eigen::Matrix2d jacobian;
    std::array<Eigen::Matrix2d, 2> hessians;
};

AnkleDerivatives EvaluateAnkleDerivatives(const AnkleParams& params, double pitch, double roll);

}  // namespace parallel_ankle::detail
