#!/usr/bin/env bash
# Sourced by verification tools. The calling shell holds descriptor 9 throughout
# its work; never launch a background job and exit while holding this lock.

acquire_heavy_work_lock() {
    local heavy_lock="${PARALLEL_ANKLE_HEAVY_LOCK:-${HOME}/.cache/codex-heavy-work.lock}"
    local wait_seconds="${PARALLEL_ANKLE_LOCK_WAIT_SECONDS:-0}"
    if [[ ! ${wait_seconds} =~ ^[0-9]+$ ]] || (( ${#wait_seconds} > 4 )) || (( 10#${wait_seconds} > 3600 )); then
        echo "PARALLEL_ANKLE_LOCK_WAIT_SECONDS must be an integer from 0 to 3600." >&2
        return 2
    fi
    mkdir -p -- "$(dirname -- "${heavy_lock}")"
    # Reuse an inherited descriptor only when it names this exact lock inode.
    # BASHPID follows (...) subshells; $$ can still name their parent process.
    # flock still takes/verifies an exclusive kernel lock on that descriptor;
    # an environment flag alone can never bypass locking. This lets one outer
    # validation hold the lock through nested build/package scripts.
    if [[ ! /proc/${BASHPID}/fd/9 -ef "${heavy_lock}" ]]; then
        exec 9>>"${heavy_lock}"
    fi
    local lock_status=0
    flock -w "${wait_seconds}" -E 75 9 || lock_status=$?
    if (( lock_status != 0 )); then
        if (( lock_status == 75 )); then
            echo "Heavy-work lock is busy; validation was deferred." >&2
        fi
        return "${lock_status}"
    fi

    export CMAKE_BUILD_PARALLEL_LEVEL=2 CTEST_PARALLEL_LEVEL=1 MAX_JOBS=2 RAYON_NUM_THREADS=2
    export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
    free -h
    uptime
    ps -eo pid,ppid,pcpu,pmem,rss,comm --sort=-pcpu | head -n 12 || true
    python3 - <<'PY'
import os
from pathlib import Path

memory = {}
for line in Path('/proc/meminfo').read_text().splitlines():
    key, value = line.split(':', 1)
    memory[key] = int(value.split()[0]) * 1024
reserve = max(12 * 1024**3, memory['MemTotal'] * 0.2)
cpus = len(os.sched_getaffinity(0))
load = os.getloadavg()[0]
for kind in ('cpu', 'memory', 'io'):
    pressure = Path('/proc/pressure') / kind
    if pressure.exists():
        print(f'{kind} PSI: {pressure.read_text().strip()}')
if memory['MemAvailable'] < reserve or load >= 0.75 * cpus:
    print('Insufficient resource headroom; validation was deferred.')
    raise SystemExit(75)
PY
}
