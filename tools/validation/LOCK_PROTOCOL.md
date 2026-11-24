# Lightweight lock protocol validation

Run from the repository root:

```bash
./scripts/check_lock_protocol.sh
```

This checks file-descriptor (FD) and open-file-description (OFD) lock behavior,
not machine resource headroom. It creates two private temporary files, overrides
the lock-path environment variable, and stubs `free`, `uptime`, `ps` and the
Python admission gate. The production guard's inode comparison and actual
`flock` calls remain active. It never takes the real shared heavy-work lock,
compiles code or performs numerical/data processing. Temporary files are removed
on exit; flock attempts do not wait and independent child checks have a five
second timeout.

## Regression and expected behavior

In Bash, `$$` still identifies the parent shell inside `(...)`; `BASHPID`
identifies the currently executing shell. If a parent holds file A on FD 9 but a
subshell redirects its own FD 9 to B, examining `/proc/$$/fd/9` observes A and can
mistakenly reuse and lock B. The guard must inspect
`/proc/${BASHPID}/fd/9` before deciding whether to reuse the current descriptor.

An inherited FD shares its OFD with the parent. Reacquiring `flock` through that
inherited descriptor preserves reentrant exclusion; closing a child copy does
not release the parent's remaining copy. A separately opened descriptor is an
independent lock contender.

| Check | Expected result |
| --- | --- |
| Parent holds A; subshell redirects FD 9 to B and calls the guard for A | Guard reopens A and returns 75, because the parent still holds A |
| Nested script inherits A and calls the guard | Success; an independent opener still returns 75 after the script exits |
| Ordinary documentation subshell acquires A | Independent opener returns 75 inside the block |
| Documentation subshell exits with no parent copy | Independent opener returns 0 |

On 2026-11-12, all four checks passed with the `BASHPID` guard. A preceding
temporary comparison reproduced the old `$$` failure: the wrong-inode subshell
returned 0 while its actual FD 9 still named B. The retained check exercises the
current guard and fails if that regression returns. These results establish
the tested lock lifecycle only; they do not validate memory/load thresholds or
certify that other running conversations follow the resource convention.
