# Independent Jacobian convergence probe

This optional probe checks analytical and finite-difference Jacobians at a
high-curvature pose near the selected crank branch boundary.

The fixed point is the current `rotary_v02` right foot at the exact C++ doubles
`pitch=-1.2`, `roll=0.18000000000000005`. The C++ collector loads the production
YAML, independently checks all used stencil points with long-double scalar
rotations and rod closure, and writes round-trip binary-double metadata and
analytic/forward/central matrices for steps from `1e-4` through `1e-9` rad.
Central denominators use the actual representable coordinate spans.

The Python tool reconstructs each metadata number as its exact binary-float
rational. With mpmath at 90 decimal digits it independently intersects the
crank circle and rod sphere using geometric `atan2`, selects the smooth family
only from the base motor angle, then applies `mp.diff`. It uses no production
Python bindings, NumPy, SciPy or YAML parser. mpmath must already be available;
the commands do not install dependencies.

For an existing verified static, YAML-enabled producer, run from the repository
root under the common resource lock:

```bash
(
    set -euo pipefail
    source scripts/resource_limits.sh
    acquire_heavy_work_lock
    mkdir -p build/jacobian-probe
    timeout 120s nice -n 10 c++ -std=c++17 -O2 \
        -Iinclude $(pkg-config --cflags eigen3 yaml-cpp) \
        tools/validation/jacobian_convergence_probe.cpp \
        build/libparallel_ankle_config.a build/libparallel_ankle.a \
        $(pkg-config --libs yaml-cpp) \
        -o build/jacobian-probe/collector
    timeout 15s nice -n 10 build/jacobian-probe/collector \
        configs/rotary_ankle.yaml build/jacobian-probe/result
    timeout 30s nice -n 10 python3 tools/validation/jacobian_convergence_mpmath.py \
        build/jacobian-probe/result --dps 90
)
```

Compilation is one bounded native compiler process. The subshell stops on a
lock/resource failure and releases its descriptor on exit. Thread limits and
resource checks come from the sourced guard. Adjust library paths for other producer
directories; shared-library runtime lookup is not configured by this example.
The four output files are `.meta.json`, `.csv`, `.mp.csv`, and the collector
executable. Retain input geometry and source revision when exporting results.

The reference run found relative
Frobenius errors of `2.02e-12` for the analytic matrix, `1.96e-3` for a forward
step of `1e-6`, and `8.14e-8` for a central step of `1e-7`. Reducing the central
step to `1e-9` increased error to `9.20e-7`, showing the expected roundoff limit.
This is an accuracy probe at one difficult point, not a performance benchmark
or evidence about every mechanism branch.
