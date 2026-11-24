#include "parallel_ankle/config.h"

#include <yaml-cpp/yaml.h>

#include <cmath>
#include <set>
#include <stdexcept>

namespace parallel_ankle {
namespace {

[[noreturn]] void InvalidField(const std::string& field, const YAML::Node& node, const std::string& message) {
    const auto mark = node.Mark();
    const std::string location = mark.is_null() ? "" : " (line " + std::to_string(mark.line + 1) +
                                                            ", column " + std::to_string(mark.column + 1) + ")";
    throw std::invalid_argument(field + ": " + message + location);
}

void RequireMap(const YAML::Node& node, const std::string& field) {
    if (!node.IsMap()) {
        InvalidField(field, node, "expected a mapping");
    }
    std::set<std::string> keys;
    for (const auto& item : node) {
        if (!item.first.IsScalar()) {
            InvalidField(field, item.first, "mapping keys must be strings");
        }
        const auto key = item.first.as<std::string>();
        if (!keys.insert(key).second) {
            InvalidField(field + "." + key, item.first, "duplicate key");
        }
    }
}

YAML::Node RequireField(const YAML::Node& node, const std::string& key, const std::string& field) {
    const auto value = node[key];
    if (!value.IsDefined()) {
        InvalidField(field, node, "missing required field");
    }
    return value;
}

template <typename T>
T ParseScalar(const YAML::Node& value, const std::string& field, const std::string& expected) {
    if (!value.IsScalar()) {
        InvalidField(field, value, "expected " + expected);
    }
    try {
        return value.as<T>();
    } catch (const YAML::Exception&) {
        InvalidField(field, value, "expected " + expected);
    }
}

template <typename T>
T ReadScalar(const YAML::Node& node, const std::string& key, const std::string& path, const std::string& expected) {
    const auto field = path.empty() ? key : path + "." + key;
    return ParseScalar<T>(RequireField(node, key, field), field, expected);
}

double ParseFinite(const YAML::Node& value, const std::string& field) {
    const double number = ParseScalar<double>(value, field, "a finite number");
    if (!std::isfinite(number)) {
        InvalidField(field, value, "expected a finite number");
    }
    return number;
}

double ReadFinite(const YAML::Node& node, const std::string& key, const std::string& path, bool positive = false) {
    const auto field = path + "." + key;
    const auto value = RequireField(node, key, field);
    const double number = ParseFinite(value, field);
    if (positive && number <= 0.0) {
        InvalidField(field, value, "length must be positive");
    }
    return number;
}

Eigen::Vector3d ParseVector(const YAML::Node& node, const std::string& key, const std::string& path) {
    const auto field = path + "." + key;
    const auto value = RequireField(node, key, field);
    if (!value.IsSequence() || value.size() != 3) {
        InvalidField(field, value, "expected exactly three coordinates");
    }
    Eigen::Vector3d vector;
    for (int index = 0; index < 3; ++index) {
        vector(index) = ParseFinite(value[index], field + "[" + std::to_string(index) + "]");
    }
    return vector;
}

AnkleParams ParseParams(const YAML::Node& node, const std::string& path) {
    RequireMap(node, path);
    AnkleParams params;
    params.r_roll_center = ParseVector(node, "r_roll_center", path);
    params.r_pitch_center = ParseVector(node, "r_pitch_center", path);
    params.r_a1 = ParseVector(node, "r_a1", path);
    params.r_c1 = ParseVector(node, "r_c1", path);
    params.l_bar1 = ReadFinite(node, "l_bar1", path, true);
    params.l_rod1 = ReadFinite(node, "l_rod1", path, true);
    params.alpha_1_rad = ReadFinite(node, "alpha_1_rad", path);
    params.is_left_motor1 = ReadScalar<bool>(node, "is_left_motor1", path, "a boolean");
    params.r_a2 = ParseVector(node, "r_a2", path);
    params.r_c2 = ParseVector(node, "r_c2", path);
    params.l_bar2 = ReadFinite(node, "l_bar2", path, true);
    params.l_rod2 = ReadFinite(node, "l_rod2", path, true);
    params.alpha_2_rad = ReadFinite(node, "alpha_2_rad", path);
    params.is_left_motor2 = ReadScalar<bool>(node, "is_left_motor2", path, "a boolean");
    try {
        params.Validate();
    } catch (const std::invalid_argument& error) {
        InvalidField(path, node, error.what());
    }
    return params;
}

}  // namespace

AnkleConfig LoadAnkleConfig(const std::string& path, const std::string& profile) {
    try {
        const YAML::Node root = YAML::LoadFile(path);
        RequireMap(root, "document");
        AnkleConfig config;
        config.profile = profile;
        if (config.profile.empty()) {
            const auto robot = ReadScalar<std::string>(root, "robot_type", "", "a nonempty string");
            const auto version = ReadScalar<std::string>(root, "version", "", "a nonempty string");
            if (robot.empty() || version.empty()) {
                throw std::invalid_argument("robot_type and version must be nonempty strings");
            }
            config.profile = robot + "_" + version;
        }
        const auto selected = RequireField(root, config.profile, config.profile);
        RequireMap(selected, config.profile);
        const auto ankle_path = config.profile + ".ankle_params";
        const auto ankles = RequireField(selected, "ankle_params", ankle_path);
        RequireMap(ankles, ankle_path);
        config.left = ParseParams(RequireField(ankles, "left", ankle_path + ".left"), ankle_path + ".left");
        config.right = ParseParams(RequireField(ankles, "right", ankle_path + ".right"), ankle_path + ".right");
        return config;
    } catch (const YAML::Exception& error) {
        throw std::invalid_argument("failed to load ankle config '" + path + "': " + error.what());
    } catch (const std::invalid_argument& error) {
        throw std::invalid_argument("failed to load ankle config '" + path + "': " + error.what());
    }
}

}  // namespace parallel_ankle
