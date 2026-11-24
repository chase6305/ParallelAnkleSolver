#!/usr/bin/env bash
set -euo pipefail

project_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
build_dir="${1:-${project_dir}/build}"
source "${project_dir}/scripts/resource_limits.sh"
acquire_heavy_work_lock

timeout 60s nice -n 10 cmake -S "${project_dir}" -B "${build_dir}" -DCMAKE_BUILD_TYPE=Release \
    -DPARALLEL_ANKLE_WITH_YAML=ON -DPARALLEL_ANKLE_BUILD_TESTS=ON -DPARALLEL_ANKLE_BUILD_EXAMPLES=ON
timeout 180s nice -n 10 cmake --build "${build_dir}" --parallel 2
timeout 60s nice -n 10 ctest --test-dir "${build_dir}" --parallel 1 --output-on-failure
timeout 10s nice -n 10 "${build_dir}/ankle_demo" "${project_dir}/configs/rotary_ankle.yaml"
timeout 10s nice -n 10 "${build_dir}/ankle_demo" "${project_dir}/configs/rotary_ankle.yaml" rotary_v01
timeout 10s nice -n 10 "${build_dir}/ankle_demo" "${project_dir}/configs/rotary_ankle_reference.yaml"
