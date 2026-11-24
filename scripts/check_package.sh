#!/usr/bin/env bash
set -euo pipefail

project_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
build_dir="${project_dir}/build"
with_core_only=false
build_argument_seen=false

usage() {
    echo "Usage: $0 [existing-build-directory] [--core-only]"
    echo "Run scripts/check.sh first to build the YAML-enabled producer."
    echo "--core-only also checks a core-only install, reusing the embedded build."
}

for argument in "$@"; do
    case "${argument}" in
        --help|-h)
            usage
            exit 0
            ;;
        --core-only)
            with_core_only=true
            ;;
        -*)
            echo "Unknown option: ${argument}" >&2
            usage >&2
            exit 2
            ;;
        *)
            if ${build_argument_seen}; then
                echo "Only one producer build directory may be specified." >&2
                usage >&2
                exit 2
            fi
            build_dir="${argument}"
            build_argument_seen=true
            ;;
    esac
done

if [[ ! -f "${build_dir}/CMakeCache.txt" || ! -f "${build_dir}/cmake_install.cmake" ]]; then
    echo "No configured producer at ${build_dir}; run scripts/check.sh first." >&2
    exit 2
fi
build_dir="$(cd -- "${build_dir}" && pwd)"
python3 - "${project_dir}" "${build_dir}/CMakeCache.txt" <<'PY'
import sys
from pathlib import Path

project, cache_file = map(Path, sys.argv[1:])
cache = {}
for line in cache_file.read_text().splitlines():
    if not line or line.startswith(('#', '//')) or '=' not in line or ':' not in line:
        continue
    key, value = line.split('=', 1)
    cache[key.split(':', 1)[0]] = value
if Path(cache.get('CMAKE_HOME_DIRECTORY', '')).resolve() != project.resolve():
    raise SystemExit('Producer build belongs to another source directory; run scripts/check.sh first.')
if cache.get('PARALLEL_ANKLE_WITH_YAML', '').upper() not in ('ON', 'TRUE', 'YES', '1'):
    raise SystemExit('Producer must include the YAML loader; run scripts/check.sh first.')
PY

source "${project_dir}/scripts/resource_limits.sh"
acquire_heavy_work_lock

prefix="${build_dir}/package-prefix"
checks_dir="${build_dir}/package-checks"
mkdir -p -- "${checks_dir}"

check_consumer() {
    local name="$1"
    local find_mode="$2"
    local install_prefix="$3"
    shift 3
    local consumer_build="${checks_dir}/${name}"
    echo "Checking package consumer: ${name}"
    timeout 60s nice -n 10 cmake -S "${project_dir}/tests/consumer" -B "${consumer_build}" \
        -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH="${install_prefix}" \
        -DCONSUMER_FIND_MODE="${find_mode}" "$@"
    timeout 120s nice -n 10 cmake --build "${consumer_build}" --config Release --parallel 2
    timeout 30s nice -n 10 ctest --test-dir "${consumer_build}" -C Release --parallel 1 --output-on-failure
}

expect_configure_failure() {
    local name="$1"
    local expected_message="$2"
    shift 2
    local log_file="${checks_dir}/${name}.log"
    local configure_status=0
    timeout 60s nice -n 10 cmake -S "${project_dir}/tests/consumer" -B "${checks_dir}/${name}" \
        -DCMAKE_BUILD_TYPE=Release "$@" >"${log_file}" 2>&1 || configure_status=$?
    if (( configure_status == 0 )); then
        echo "Expected package configuration to fail: ${name}" >&2
        return 1
    fi
    if (( configure_status == 124 || configure_status == 137 )); then
        cat -- "${log_file}" >&2
        return "${configure_status}"
    fi
    if ! rg --fixed-strings --quiet -- "${expected_message}" "${log_file}"; then
        echo "Package configure failed for an unexpected reason: ${name}" >&2
        cat -- "${log_file}" >&2
        return 1
    fi
    echo "Confirmed expected package rejection: ${name}"
}

echo "Installing existing producer into ${prefix}"
timeout 30s nice -n 10 cmake --install "${build_dir}" --config Release --prefix "${prefix}"
check_consumer installed-core core "${prefix}" -DCMAKE_DISABLE_FIND_PACKAGE_yaml-cpp=TRUE
check_consumer installed-config config "${prefix}"
check_consumer optional-config optional-config "${prefix}" -DCMAKE_DISABLE_FIND_PACKAGE_yaml-cpp=TRUE
expect_configure_failure unsupported-component "Unsupported ParallelAnkle component: bogus" \
    -DCMAKE_PREFIX_PATH="${prefix}" -DCONSUMER_FIND_MODE=bogus

echo "Checking embedded core without yaml-cpp or GoogleTest"
embedded_build="${checks_dir}/embedded"
timeout 60s nice -n 10 cmake -S "${project_dir}/tests/consumer/embedded" -B "${embedded_build}" \
    -DCMAKE_BUILD_TYPE=Release -DPARALLEL_ANKLE_SOURCE="${project_dir}" \
    -DCMAKE_DISABLE_FIND_PACKAGE_yaml-cpp=TRUE -DCMAKE_DISABLE_FIND_PACKAGE_GTest=TRUE
timeout 120s nice -n 10 cmake --build "${embedded_build}" --config Release --parallel 2
timeout 30s nice -n 10 ctest --test-dir "${embedded_build}" -C Release --parallel 1 --output-on-failure

if ${with_core_only}; then
    # The embedded subproject is already an Eigen-only producer. Install it
    # directly instead of compiling the same library in another build tree.
    core_prefix="${build_dir}/package-core-prefix"
    timeout 30s nice -n 10 cmake --install "${embedded_build}/parallel_ankle" --config Release --prefix "${core_prefix}"
    check_consumer core-only-install core "${core_prefix}" -DCMAKE_DISABLE_FIND_PACKAGE_yaml-cpp=TRUE
    expect_configure_failure missing-config-component "ParallelAnkle was built without the config component" \
        -DCMAKE_PREFIX_PATH="${core_prefix}" -DCONSUMER_FIND_MODE=config
fi

echo "Package checks passed."
