# Reproducibility code and data

**Temporal filtering of dispersive stiffness in compacton collisions: coherent-scattering benchmarks and accuracy-cost tradeoffs**

This repository contains the code, raw and processed data, provenance records,
figure sources, and verification scripts used by the numerical study.

Release **v1.1.0** extends v1.0.0 with four reproducibility blocks:

1. the corrected first-order generalized-alpha campaign, including the complete
   microbank, fine-grid boundary test, two Section 4.1 collisions, validation,
   scheduler provenance, and Figure 3 data;
2. the collision-Jacobian and fixed-state mesh-scaling experiments, with
   reconstructed authoritative drivers and full independent numerical checks;
3. the four robustness-screen banks, including every `cases.csv`, every
   per-run summary, the MATLAB R2024b production chain, SLURM launchers,
   parameter audits, and Newton-cap counts;
4. editable sources and CSV inputs for all manuscript figures.

The original Direct-WENO and Petrov-Galerkin benchmark material from v1.0.0 is
retained unchanged.

## Quick verification

Create an environment with NumPy, SciPy, Numba and Matplotlib (see
`requirements.txt` and `ENVIRONMENT.md`), then run

```bash
python scripts/freeze_check.py
```

A successful run ends with

```text
FREEZE_CHECK_V1_1_0_OK
```

The default command verifies the frozen raw-input hashes, reruns the small
corrected generalized-alpha validation, audits the complete Picasso outputs,
checks the Jacobian spectra and mesh-scaling tables, audits the four MATLAB
robustness banks, and regenerates all manuscript figure panels. To recompute
the v1.0.0 Direct-WENO/Petrov-Galerkin benchmark chain from raw trajectories, run

```bash
python scripts/freeze_check.py --full
```

Both the default and `--full` audits run in disposable temporary copies. They do not overwrite
the frozen CSV/PDF artifacts in the working tree.

For the more expensive full Jacobian recomputation from initial states, run

```bash
python reproducibility/jacobian/verify_jacobian_traceability.py --full
```

This full ARPACK recomputation also runs in a disposable copy because last-bit
eigenvalue variations are platform dependent.

## Directory map

- `code/`: Direct-WENO and Petrov-Galerkin benchmark code retained from v1.0.0.
- `raw/`: raw Direct-WENO and Petrov-Galerkin trajectories used by the final
  coherent-core and invariant-balance comparisons.
- `reproducibility/generalized_alpha/`: corrected generalized-alpha solver,
  production drivers, complete case-level outputs, Picasso logs, and verifier.
- `reproducibility/jacobian/`: Jacobian drivers, raw snapshots, frozen and
  regenerated spectra, mesh-scaling data, reports, and verifier.
- `provenance/robustness_banks/`: complete case tables, complete per-run
  summaries, MATLAB production sources, case generators, and SLURM launchers.
- `data/`: compact tables used by the manuscript and release audits, including
  both Section 4.1 collisions and the Newton-cap audit.
- `figure_sources/`: closed source/data chain for Figures 1-3.
- `scripts/fig*.py`: editable sources for all manuscript figure panels.
- `figures/`: current PDF figure panels.
- `scripts/`: verification, regeneration, and release-audit entry points.
- `audit/`: verification documentation and release inventory.

## Provenance qualifications

The original historical Jacobian eigensolver drivers were not recovered
byte-for-byte. The included drivers are authoritative reconstructions from the
published semidiscrete equations, recorded parameters, frozen outputs, and
reports. A complete recomputation reproduces the archived eigenvalues,
filtering factors, and fitted mesh exponents within the documented numerical
tolerances. This limitation is preserved in
`reproducibility/jacobian/README.md` and its verification report.

The multi-gigabyte MATLAB trajectory files of the robustness screens are not
duplicated. The release includes every case definition, every per-run summary,
the complete production and analysis source chain, and sufficient information
to rerun any selected row.

Superseded generalized-alpha outputs produced by the defective historical
implementation are not included as valid data and are not used by any figure
or manuscript result.

## Final independent benchmark

The quantitative coherent-core reference remains the absolute Direct-WENO
N=3200 fitted two-core model:

```text
A_s = 1.34807979416512
X_s = 53.783521434422525
A_f = 3.9971963548196703
X_f = 62.517459194972844
```

The absolute-flow N=1600-to-N=3200 fitted-core difference is
`4.426107800112e-4`. It is a practical last-refinement indicator, not a
rigorous error bound.

## License and citation

Code, data and figures are released under the MIT License. Citation metadata
are provided in `CITATION.cff` and `.zenodo.json`.
