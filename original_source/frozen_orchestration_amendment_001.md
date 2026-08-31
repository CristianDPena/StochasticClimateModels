# Frozen Orchestration Amendment 001

Date: 2026-07-24

## Scope

This amendment changes orchestration logging and integrity checks only. It
does not modify any simulation, statistic, frequency, block, window,
likelihood, optimizer, transfer profile, parameter order, seed, or estimate.

## Reason

The first frozen one-worker launch showed that Windows starts the production
R calculation through this process tree:

```text
Rscript.exe launcher
  cmd.exe
    x64/Rscript.exe calculation
```

The freeze-time controller recorded the working set of only the launcher.
That would understate per-run peak process memory.

## Change

`work/run_overarching_validation_batch.ps1` now:

- snapshots the Windows process table, identifies all descendants of each
  tracked launcher, and records the sum of their working sets as
  `peak_process_working_set_mb`;
- recomputes every frozen statistical-source, kernel, profile, and design hash
  before launching a batch and aborts on any mismatch.
- calls `WaitForExit`, validates the one-row canonical result, and records an
  unavailable Windows launcher exit code explicitly instead of misclassifying
  a valid completed result as failed.

The global R-process memory trace is unchanged. It supplies the equivalent
one-worker measurement for the first attempt, which was already running when
this amendment was made.

## Integrity

Freeze-time controller SHA-256:

```text
5d51665c0e0067ffd7d71c69c494bc8bc48368eda5551ef98de9f4a26eb95b90
```

Final amended controller SHA-256:

```text
e55675359a02f641a4d0d57277446fa1d906b89daa6f02f6eef19d5f1374ceed
```

All eight core statistical source hashes and all six transfer-profile hashes
remain unchanged.
