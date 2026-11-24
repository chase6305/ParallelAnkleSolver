#!/usr/bin/env bash
set -euo pipefail

project_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
build_dir="${project_dir}/build"
with_benchmark=false
seen_build=false
for argument in "$@"; do
    case "${argument}" in
        --help|-h)
            echo "Usage: $0 [build-directory] [--benchmark]"
            echo "Build, test, and verify installed/embedded use under one shared resource lock."
            exit 0
            ;;
        --benchmark) with_benchmark=true ;;
        -*) echo "Unknown option: ${argument}" >&2; exit 2 ;;
        *)
            if ${seen_build}; then
                echo "Only one build directory may be specified." >&2
                exit 2
            fi
            build_dir="${argument}"
            seen_build=true
            ;;
    esac
done
mkdir -p -- "${build_dir}"
build_dir="$(cd -- "${build_dir}" && pwd)"
source "${project_dir}/scripts/resource_limits.sh"
acquire_heavy_work_lock

"${project_dir}/scripts/check.sh" "${build_dir}"
"${project_dir}/scripts/check_package.sh" "${build_dir}" --core-only
if ${with_benchmark}; then
    acquire_heavy_work_lock
    timeout 60s nice -n 10 cmake -S "${project_dir}" -B "${build_dir}" -DPARALLEL_ANKLE_BUILD_BENCHMARKS=ON
    timeout 120s nice -n 10 cmake --build "${build_dir}" --target benchmark_jacobians --parallel 2
    timeout 15s nice -n 10 "${build_dir}/benchmark_jacobians" 5000 --kinematics
fi
