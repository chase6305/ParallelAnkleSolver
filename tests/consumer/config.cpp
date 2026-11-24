#include <parallel_ankle/config.h>
#include <parallel_ankle/parallel_ankle.h>

int main(int argc, char** argv) {
    if (argc != 2) {
        return 2;
    }
    const auto config = parallel_ankle::LoadAnkleConfig(argv[1]);
    const parallel_ankle::PMAnkle ankle(true, config.left);
    return !config.profile.empty() && ankle.GetAnkleMotorAngles(0.0, 0.0).success() ? 0 : 1;
}
