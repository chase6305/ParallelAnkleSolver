#include <gtest/gtest.h>

#include <chrono>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <stdexcept>
#include <string>

#include "fixtures.h"
#include "parallel_ankle/config.h"

namespace parallel_ankle::test {
namespace {

const std::filesystem::path kConfigs = std::filesystem::path(PARALLEL_ANKLE_SOURCE_DIR) / "configs";

class TempConfig {
public:
    explicit TempConfig(const std::string& text) {
        const auto stamp = std::chrono::steady_clock::now().time_since_epoch().count();
        path_ = std::filesystem::temp_directory_path() / ("parallel-ankle-" + std::to_string(stamp) + ".yaml");
        std::ofstream file(path_);
        file.exceptions(std::ios::badbit | std::ios::failbit);
        file << text;
    }
    ~TempConfig() {
        std::error_code error;
        std::filesystem::remove(path_, error);
    }
    std::string path() const { return path_.string(); }

private:
    std::filesystem::path path_;
};

std::string LegacyYaml() {
    std::ifstream file(kConfigs / "rotary_ankle_reference.yaml");
    return {(std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>()};
}

std::string ReplaceField(std::string text, const std::string& original, const std::string& replacement) {
    const auto position = text.find(original);
    if (position == std::string::npos) {
        throw std::logic_error("test YAML field not found: " + original);
    }
    text.replace(position, original.size(), replacement);
    return text;
}

void ExpectConfigError(const std::string& text, const std::string& profile, const std::string& field) {
    const TempConfig config(text);
    try {
        LoadAnkleConfig(config.path(), profile);
        FAIL() << "expected invalid configuration at " << field;
    } catch (const std::invalid_argument& error) {
        const std::string message = error.what();
        EXPECT_NE(message.find(config.path()), std::string::npos) << message;
        EXPECT_NE(message.find(field), std::string::npos) << message;
    }
}

TEST(Config, SelectCurrentV02ByDefaultAndPreserveV01) {
    const auto config = LoadAnkleConfig((kConfigs / "rotary_ankle.yaml").string());
    EXPECT_EQ(config.profile, "rotary_v02");
    EXPECT_DOUBLE_EQ(config.left.l_bar1, 30.0);
    EXPECT_DOUBLE_EQ(config.right.l_rod2, 186.209);
    EXPECT_DOUBLE_EQ(config.left.r_roll_center.z(), -15.0);
    EXPECT_DOUBLE_EQ(config.left.r_pitch_center.z(), -5.0);
    EXPECT_TRUE(config.left.is_left_motor1);
    EXPECT_FALSE(config.right.is_left_motor1);
    const auto old = LoadAnkleConfig((kConfigs / "rotary_ankle.yaml").string(), "rotary_v01");
    EXPECT_DOUBLE_EQ(old.left.l_bar1, 60.0);
    EXPECT_DOUBLE_EQ(old.left.r_roll_center.z(), 3.0);
    EXPECT_DOUBLE_EQ(old.left.r_pitch_center.z(), 3.0);
    EXPECT_FALSE(old.left.is_left_motor1);
    EXPECT_TRUE(old.right.is_left_motor1);
}

TEST(Config, LegacyYamlMatchesOriginalTestFixture) {
    const auto config = LoadAnkleConfig((kConfigs / "rotary_ankle_reference.yaml").string());
    for (bool is_left : {true, false}) {
        const auto& loaded = is_left ? config.left : config.right;
        const auto expected = LegacyParams(is_left);
        EXPECT_TRUE(loaded.r_a1.isApprox(expected.r_a1));
        EXPECT_TRUE(loaded.r_a2.isApprox(expected.r_a2));
        EXPECT_TRUE(loaded.r_c1.isApprox(expected.r_c1));
        EXPECT_TRUE(loaded.r_c2.isApprox(expected.r_c2));
        EXPECT_TRUE(loaded.r_roll_center.isApprox(expected.r_roll_center));
        EXPECT_TRUE(loaded.r_pitch_center.isApprox(expected.r_pitch_center));
        EXPECT_EQ(loaded.is_left_motor1, expected.is_left_motor1);
        EXPECT_EQ(loaded.is_left_motor2, expected.is_left_motor2);
        EXPECT_DOUBLE_EQ(loaded.l_bar1, expected.l_bar1);
        EXPECT_DOUBLE_EQ(loaded.l_bar2, expected.l_bar2);
        EXPECT_DOUBLE_EQ(loaded.l_rod1, expected.l_rod1);
        EXPECT_DOUBLE_EQ(loaded.l_rod2, expected.l_rod2);
        EXPECT_DOUBLE_EQ(loaded.alpha_1_rad, expected.alpha_1_rad);
        EXPECT_DOUBLE_EQ(loaded.alpha_2_rad, expected.alpha_2_rad);
    }
}

TEST(Config, RejectMissingAndMalformedConfigurations) {
    EXPECT_THROW(LoadAnkleConfig((kConfigs / "does-not-exist.yaml").string()), std::invalid_argument);
    EXPECT_THROW(LoadAnkleConfig((kConfigs / "rotary_ankle.yaml").string(), "unknown_profile"), std::invalid_argument);
    TempConfig empty("{}");
    EXPECT_THROW(LoadAnkleConfig(empty.path()), std::invalid_argument);
    TempConfig missing("test:\n  ankle_params:\n    left: {}\n    right: {}\n");
    EXPECT_THROW(LoadAnkleConfig(missing.path(), "test"), std::invalid_argument);

    std::ifstream file(kConfigs / "rotary_ankle_reference.yaml");
    const std::string original((std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());
    ASSERT_FALSE(original.empty());
    std::string invalid_length = original;
    const auto length_position = invalid_length.find("l_bar1: 60.0");
    ASSERT_NE(length_position, std::string::npos);
    invalid_length.replace(length_position, std::string("l_bar1: 60.0").size(), "l_bar1: -1.0");
    TempConfig negative(invalid_length);
    EXPECT_THROW(LoadAnkleConfig(negative.path()), std::invalid_argument);

    std::string invalid_vector = original;
    const auto vector_position = invalid_vector.find("r_a1: [-46.0, 0.0, 257.94]");
    ASSERT_NE(vector_position, std::string::npos);
    invalid_vector.replace(vector_position, std::string("r_a1: [-46.0, 0.0, 257.94]").size(), "r_a1: [0.0, 0.0]");
    TempConfig short_vector(invalid_vector);
    EXPECT_THROW(LoadAnkleConfig(short_vector.path()), std::invalid_argument);
}

TEST(Config, IdentifyMissingFieldsAndWrongRootTypes) {
    for (const std::string root : {"", "null", "[]", "[rotary, v02]", "42"}) {
        ExpectConfigError(root, "", "document");
    }
    ExpectConfigError("robot_type: rotary\n", "", "version");
    ExpectConfigError("robot_type: rotary\nversion: ''\n", "", "nonempty");
    ExpectConfigError("test:\n  ankle_params:\n    left: {}\n", "test", "test.ankle_params.left.r_roll_center");
}

TEST(Config, ReportFieldPathsForNonfiniteAndInvalidScalars) {
    const auto source = LegacyYaml();
    ASSERT_FALSE(source.empty());
    for (const std::string value : {".nan", ".inf", "-.inf", "not-a-number"}) {
        ExpectConfigError(ReplaceField(source, "l_bar1: 60.0", "l_bar1: " + value), "",
                          "rotary_reference.ankle_params.left.l_bar1");
    }
    ExpectConfigError(ReplaceField(source, "l_bar1: 60.0", "l_bar1: -1.0"), "",
                      "rotary_reference.ankle_params.left.l_bar1");
    ExpectConfigError(ReplaceField(source, "is_left_motor1: true", "is_left_motor1: 2"), "",
                      "rotary_reference.ankle_params.left.is_left_motor1");
    ExpectConfigError(ReplaceField(source, "r_a1: [-46.0, 0.0, 257.94]", "r_a1: [-46.0, .nan, 257.94]"), "",
                      "rotary_reference.ankle_params.left.r_a1[1]");
    ExpectConfigError(ReplaceField(source, "alpha_1_rad: 0.1745", "alpha_1_rad: [0.1]"), "",
                      "rotary_reference.ankle_params.left.alpha_1_rad");
}

TEST(Config, RejectAmbiguousDuplicateKeys) {
    ExpectConfigError("robot_type: rotary\nrobot_type: other\nversion: v02\n", "", "document.robot_type: duplicate key");
    const auto source = LegacyYaml();
    ASSERT_FALSE(source.empty());
    ExpectConfigError(ReplaceField(source, "l_bar1: 60.0", "l_bar1: 60.0\n      l_bar1: 30.0"), "",
                      "rotary_reference.ankle_params.left.l_bar1: duplicate key");
}

TEST(Config, ExplicitProfileDoesNotRequireDefaultSelectionFields) {
    auto source = LegacyYaml();
    ASSERT_FALSE(source.empty());
    source = ReplaceField(source, "robot_type: rotary\nversion: reference\n", "");
    const TempConfig config(source);
    EXPECT_EQ(LoadAnkleConfig(config.path(), "rotary_reference").profile, "rotary_reference");
    EXPECT_THROW(LoadAnkleConfig(config.path()), std::invalid_argument);
}

class ProfileTest : public ::testing::TestWithParam<const char*> {};

TEST_P(ProfileTest, BothFeetCloseRodsAndRecoverPoses) {
    const std::string profile = GetParam();
    const auto path = kConfigs / (profile == "rotary_reference" ? "rotary_ankle_reference.yaml" : "rotary_ankle.yaml");
    const auto config = LoadAnkleConfig(path.string(), profile);
    for (bool is_left : {true, false}) {
        const auto& params = is_left ? config.left : config.right;
        const PMAnkle ankle(is_left, params);
        for (double pitch : {-0.2, 0.0, 0.31}) {
            for (double roll : {-0.1, 0.0, 0.1}) {
                const auto ik = ankle.GetAnkleMotorAngles(pitch, roll);
                ASSERT_TRUE(ik.success()) << profile << ", " << is_left << ", " << pitch << ", " << roll;
                ExpectRodClosure(params, {pitch, roll}, ik.motor_angles());
                const auto fk = ankle.SolveForwardKinematics(ik.theta_1_rad, ik.theta_2_rad);
                ASSERT_TRUE(fk.success()) << profile << ": " << ToString(fk.status);
                EXPECT_LT((fk.pose - Eigen::Vector2d(pitch, roll)).norm(), 2e-8);
                EXPECT_LT(fk.residual_norm, 1e-9);
            }
        }
    }
}

INSTANTIATE_TEST_SUITE_P(RotaryAnkleProfiles, ProfileTest, ::testing::Values("rotary_v01", "rotary_v02", "rotary_reference"));

}  // namespace
}  // namespace parallel_ankle::test
