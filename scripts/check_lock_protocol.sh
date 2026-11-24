#!/usr/bin/env bash
# Lightweight FD/open-file-description lock protocol check, not a resource
# headroom check. Uses private temporary locks and stubs the resource reports
# and Python admission gate. Never compiles, processes data or takes the real
# shared heavy-work lock.
set -euo pipefail

if [[ ${1:-} == --help || ${1:-} == -h ]]; then
    echo "Usage: $0"
    echo "Check lock FD/OFD inheritance, wrong-inode rejection and exit release."
    echo "Uses private temporary locks; resource reports and the Python gate are stubbed."
    echo "This is not resource-headroom validation and runs no build or data computation."
    exit 0
fi
if (( $# != 0 )); then
    echo "Unexpected argument; use --help." >&2
    exit 2
fi

project_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
lock_test_dir=$(mktemp -d "${TMPDIR:-/tmp}/parallel-ankle-lock-protocol.XXXXXX")
lock_a="${lock_test_dir}/A"
lock_b="${lock_test_dir}/B"
cleanup() {
    exec 9>&-
    [[ ! -e ${lock_a} ]] || unlink -- "${lock_a}"
    [[ ! -e ${lock_b} ]] || unlink -- "${lock_b}"
    rmdir -- "${lock_test_dir}"
}
trap cleanup EXIT
: >"${lock_a}"
: >"${lock_b}"

# Drop only this process's inherited descriptor copy before using private files.
# Override an inherited lock path before invoking the production guard.
exec 9>&-
export PARALLEL_ANKLE_HEAVY_LOCK="${lock_a}" PARALLEL_ANKLE_LOCK_WAIT_SECONDS=0
free() { :; }
uptime() { :; }
ps() { :; }
python3() { :; }
export -f free uptime ps python3
source "${project_dir}/scripts/resource_limits.sh"

expect_competitor() {
    local expected="$1"
    local actual=0
    # Close the inherited locked descriptor before opening an independent OFD.
    timeout 5s bash -c 'exec 9>&-; flock -w 0 -E 75 "$1" true' _ "${lock_a}" || actual=$?
    if [[ ${actual} != "${expected}" ]]; then
        echo "Lock competitor returned ${actual}; expected ${expected}." >&2
        return 1
    fi
}

exec 9>>"${lock_a}"
flock -w 0 -E 75 9
(
    exec 9>>"${lock_b}"
    [[ /proc/$$/fd/9 -ef ${lock_a} && /proc/${BASHPID}/fd/9 -ef ${lock_b} ]]
    status=0
    acquire_heavy_work_lock || status=$?
    [[ ${status} == 75 && /proc/${BASHPID}/fd/9 -ef ${lock_a} ]]
    echo "PASS: wrong-inode subshell reopens A and returns 75 while parent holds A."
)

timeout 5s bash -c '
    set -euo pipefail
    source "$1"
    acquire_heavy_work_lock
    [[ /proc/${BASHPID}/fd/9 -ef ${PARALLEL_ANKLE_HEAVY_LOCK} ]]
' _ "${project_dir}/scripts/resource_limits.sh"
expect_competitor 75
echo "PASS: nested script reuses the inherited OFD; outer exclusion survives its exit."

exec 9>&-
(
    set -euo pipefail
    acquire_heavy_work_lock
    expect_competitor 75
    echo "PASS: ordinary documentation subshell excludes an independent opener."
)
expect_competitor 0
echo "PASS: documentation subshell exit releases its last lock descriptor."
echo "Lock protocol checks passed; resource headroom was not evaluated."
