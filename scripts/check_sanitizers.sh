#!/usr/bin/env bash
set -euo pipefail

if [[ ${1:-} == --help || ${1:-} == -h ]]; then
    echo "Usage: $0 [build-directory]"
    echo "Check C++ examples/tests with AddressSanitizer and UndefinedBehaviorSanitizer (GCC/Clang)."
    exit 0
fi
if (( $# > 1 )); then
    echo "Only one build directory may be specified." >&2
    exit 2
fi
project_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
build_dir="${1:-${project_dir}/build-sanitize}"
source "${project_dir}/scripts/resource_limits.sh"
acquire_heavy_work_lock

export ASAN_OPTIONS=detect_leaks=1:halt_on_error=1
export UBSAN_OPTIONS=halt_on_error=1:print_stacktrace=1
timeout 60s nice -n 10 cmake -S "${project_dir}" -B "${build_dir}" -DCMAKE_BUILD_TYPE=Debug \
    -DCMAKE_CXX_FLAGS_DEBUG="-O1 -g1 -fsanitize=address,undefined -fno-omit-frame-pointer" \
    -DCMAKE_EXE_LINKER_FLAGS="-fsanitize=address,undefined" \
    -DPARALLEL_ANKLE_WITH_YAML=ON -DPARALLEL_ANKLE_BUILD_TESTS=ON \
    -DPARALLEL_ANKLE_BUILD_EXAMPLES=ON -DPARALLEL_ANKLE_BUILD_BENCHMARKS=OFF
timeout 240s nice -n 10 cmake --build "${build_dir}" --parallel 2
timeout 120s nice -n 10 ctest --test-dir "${build_dir}" --parallel 1 --output-on-failure --timeout 20
