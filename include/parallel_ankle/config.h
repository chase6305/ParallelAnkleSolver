#pragma once

#include <string>

#include "parallel_ankle/parameters.h"

namespace parallel_ankle {

struct AnkleConfig {
    std::string profile;
    AnkleParams left;
    AnkleParams right;
};

// Reads <profile>.ankle_params.{left,right} from a YAML file.
// An empty profile selects <robot_type>_<version>, using the two YAML selection fields.
// Missing/malformed fields, duplicate keys or invalid geometry raise
// std::invalid_argument with the file and parameter path in its message.
AnkleConfig LoadAnkleConfig(const std::string& path, const std::string& profile = "");

}  // namespace parallel_ankle
