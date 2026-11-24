#include "parallel_ankle/parameters.h"

#include <cmath>
#include <stdexcept>

namespace parallel_ankle {

void AnkleParams::Validate() const {
    if (!r_roll_center.allFinite() || !r_pitch_center.allFinite() || !r_a1.allFinite() || !r_c1.allFinite() ||
        !r_a2.allFinite() || !r_c2.allFinite() || !std::isfinite(alpha_1_rad) || !std::isfinite(alpha_2_rad)) {
        throw std::invalid_argument("ankle coordinates and offsets must be finite");
    }
    for (double length : {l_bar1, l_rod1, l_bar2, l_rod2}) {
        if (!std::isfinite(length) || length <= 0.0) {
            throw std::invalid_argument("ankle bar and rod lengths must be finite and positive");
        }
    }
}

}  // namespace parallel_ankle
