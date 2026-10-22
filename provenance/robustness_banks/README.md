# Recovered robustness-bank provenance

This directory contains the compact, source-level recovery of Banks 01--04 used
in the robustness screening. It is included to make settings that are not
expanded in manuscript Table A.1 directly inspectable.

## Authoritative case settings

Each `cases/*_cases.csv` row contains the complete run configuration, including

- `post_t0` (start of the stored AUC window),
- `save_raw`, `save_snapshots`, and `nsave_target`,
- `newton_tol`, `newton_max_iter`, and `blowup_factor`,
- all collision, mesh, time-step, viscosity, method, and final-time parameters.

Banks 01, 02, and 04 use `post_t0=60`. This is a fixed observation-window
start, not necessarily a post-collision time. In particular, the reference CC
centres coincide at `t=140`. Bank 03 uses

```text
post_t0 = min(0.8*Tmax, max(10, (x02-x01)/abs(c1-c2) + 10))
```

The manuscript therefore calls the screening statistic an **observation-window
AUC**. The stored column name `rho_hw_post_auc` is historical.

## MATLAB production implementation

The claim that the four banks used the MATLAB production implementation is
confirmed by the recovered source and launch files:

- the SLURM launchers in `slurm/` load `matlab/R2024b`;
- they call `run_absK_case.m` with the corresponding `cases.csv` and array task;
- `run_absK_case.m` dispatches to the recovered M1, M2, or M3 MATLAB driver.

The source files needed for that chain are in `matlab/`; the exact Python case
generators are in `case_generators/`.

## Newton and amplitude-stop conventions

All four case tables set

```text
newton_tol      = 1e-10
newton_max_iter = 10
blowup_factor   = 50
```

In the recovered MATLAB drivers, the Newton stopping test is the Euclidean norm
of the change between consecutive iterates. The loop counter starts at one and
must remain below `newton_max_iter`, so the stored number of Newton updates is at
most nine. Reaching the cap increments `newton_fail_count`; it does not by itself
change `status`.

The amplitude stop is not `max(abs(U)) > 50`. It is

```text
max(abs(U_n)) > blowup_factor * max(abs(U_0))
```

so 50 is a multiplicative growth factor relative to the initial maximum.

`data/robustness_bank_newton_cap_audit.csv` records how often a completed run
also has one or more cap records. This is supplied for audit transparency and is
not used to alter the historical completion counts.

## Included files

- `cases/`: complete row-level case definitions for Banks 01--04;
- `summaries/`: complete per-case processed summaries;
- `matlab/`: recovered production chain;
- `case_generators/`: generators that reproduce the case tables after only the
  absolute `bank_dir` path is normalized;
- `slurm/`: recovered launch/finalization scripts;
- `RECOVERY_AUDIT.md`: numerical and provenance audit.

The multi-gigabyte raw `.mat` trajectories are not included.
