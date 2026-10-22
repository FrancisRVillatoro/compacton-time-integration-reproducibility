# Robustness-bank recovery audit — 2026-10-18

Source archive: `absK_M1M3_repro_small_20261018.tar.gz`

## Recovery status

The original production code and complete per-case CSV summaries for Banks 01–04 are present.

| Bank | cases.csv rows | summary rows | run-id match | completed by M1/M2/M3 |
|---|---:|---:|---|---|
| Bank 01 — mu/dt | 1440 | 1440 | exact | 454 / 480 / 457 |
| Bank 02 — dx refinement | 1296 | 1296 | exact | 395 / 432 / 398 |
| Bank 03 — collision severity | 1728 | 1728 | exact | 482 / 576 / 473 |
| Bank 04 — long time | 576 | 576 | exact | 123 / 192 / 122 |

The original `make_cases_bank01.py`–`make_cases_bank04.py` scripts were rerun in an isolated temporary HOME. After normalizing only the absolute `bank_dir` path, every generated row matches the archived `cases.csv` exactly for all four banks.

## Published aggregate values reproduced from summary CSVs

### Bank 01

| Method | completed | median AUC | p95 AUC |
|---|---:|---:|---:|
| M1 midpoint | 454/480 | 1.42330511857141e-4 | 4.00913082541211e-1 |
| M2 Pade R1,2 | 480/480 | 3.45052521076372e-6 | 1.83845685525307e-3 |
| M3 composition | 457/480 | 9.98754703208066e-5 | 3.64835614026380e-1 |

### Bank 03

| Method | completed | median AUC | p95 AUC |
|---|---:|---:|---:|
| M1 midpoint | 482/576 | 2.84981153304166e-3 | 5.67462167046206 |
| M2 Pade R1,2 | 576/576 | 8.40448871085146e-6 | 3.72917882774235e-3 |
| M3 composition | 473/576 | 1.59704496921077e-3 | 3.08236711780496 |

These reproduce the manuscript's quoted aggregate numbers to the displayed precision.

## Observation-window issue

The archived production configuration establishes that Banks 01, 02 and 04 use `post_t0 = 60` for all cases. For the CC reference configuration, `x01=30`, `x02=65`, `c1=1.25`, `c2=1`, so center coincidence occurs at

`t_center = (65-30)/(1.25-1) = 140`.

Therefore `[60,T]` is not a post-collision window for the CC subset. The stored variable name `rho_hw_post_auc` is historical. In the manuscript and release documentation, this statistic should be described as an observation-window AUC (or equivalent neutral terminology) unless a bank-specific collision-relative window is actually used.

Bank 03 does use a collision-relative rule from the original generator:

`post_t0 = min(0.8*Tmax, max(10, (x02-x01)/abs(c1-c2) + 10))`.

## Long-time-bank interpretation

Bank 04 has 192 cases per method pooled across three final times: 64 at T=200, 64 at T=500 and 64 at T=1000.

| Method | T=200 | T=500 | T=1000 |
|---|---:|---:|---:|
| M1 midpoint | 56/64 | 36/64 | 31/64 |
| M2 Pade R1,2 | 64/64 | 64/64 | 64/64 |
| M3 composition | 55/64 | 36/64 | 31/64 |

Thus “192 runs up to T=1000” must not be read as 192 distinct T=1000 runs.

## Production code recovered

The archive contains the original MATLAB production chain:

- `run_absK_case.m`
- `deFrutosAbsDissipative_HPC.m` — M1 implicit midpoint
- `Pade21enT_HPC.m` — M2 Pade fractional-step method
- `kpqSimulacionMP4T_HPC.m` — M3 fourth-order composition
- `kpqCompactonesPAR_IMPAR_HPC.m`
- `absK_make_operators.m`
- `absK_metrics_init.m`
- `absK_metrics_add.m`
- `absK_metrics_trim.m`
- `absK_finalize_summary.m`
- `absK_write_summary_csv.m`
- original case generators and analysis scripts.

## Metric implementation note

The bank code stores `I2 = dx*sum(abs(U)**(p+1))`, i.e. without the factor `1/(p+1)` used in the manuscript's normalized definition. Relative I2 drift is unchanged by this constant factor, so all bank claims based on relative drift are unaffected. Absolute I2 values from these bank CSVs should not be compared directly with the normalized manuscript definition without dividing by `p+1`.

## Release implication

Banks 01–04 can now be added to a new public release without the multi-gigabyte `.mat` trajectories. The complete case tables + complete per-case summaries + original production code + original analysis scripts are sufficient to regenerate all aggregate bank claims and to rerun any selected simulation from scratch with MATLAB.

The remaining release blockers are outside these four banks:

1. generalized-alpha: recovered historical PDE implementation does not match the two-variable generalized-alpha algorithm analyzed in the manuscript and must be resolved/recomputed before release;
2. Jacobian diagnostics: numerical CSVs and reports are recovered, but the original eigenvalue-production driver has not yet been recovered;
3. exact current manuscript v37 should be checked against this audit before submission/release.
