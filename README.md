# Reproducibility code and data
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23340893.svg)](https://doi.org/10.5281/zenodo.23340893)

**Temporal filtering of dispersive stiffness in compacton collisions: coherent-scattering benchmarks and accuracy–cost tradeoffs**

This repository contains code, raw trajectories, processed data, provenance,
and verification scripts for the numerical results associated with the title
above.

## Quick verification

Create an environment with NumPy, SciPy, Numba and Matplotlib (see
`requirements.txt` and `ENVIRONMENT.md`), then run

```bash
python scripts/freeze_check.py
```

A successful run ends with

```text
FREEZE_CHECK_OK
```

The check verifies raw-file hashes, re-fits the Direct-WENO N=3200 core
states, recomputes the absolute-flow spatial-refinement diagnostics, re-fits
the Petrov–Galerkin states, recomputes the sampled I2 corrections, rebuilds
the final comparison table, and regenerates the two figures that depend on
the revised benchmark.

## Contents

- `code/`: Petrov–Galerkin, Direct-WENO, core-fitting and I2 diagnostic code.
- `raw/direct_weno/`: raw Direct-WENO trajectories used by the signed-flux and
  spatial-refinement checks.
- `raw/pg_table3/`: stored Petrov–Galerkin trajectories and fitted states used
  by the final accuracy/invariant comparison.
- `raw/pg_controls/`: auxiliary I2 controls.
- `data/`: frozen numerical tables and derived diagnostics.
- `figures/`: numerical figure panels and regenerated benchmark figures.
- `scripts/`: verification and regeneration scripts.
- `provenance/`: production logs and source hashes.
- `audit/`: compact verification reports.

## Final independent benchmark

The quantitative coherent-core reference is the absolute Direct-WENO
N=3200 fitted two-core model:

```text
A_s = 1.34807979416512
X_s = 53.783521434422525
A_f = 3.9971963548196703
X_f = 62.517459194972844
```

The absolute-flow N=1600-to-N=3200 fitted-core difference is
`4.426107800112e-4`. It is used as a practical last-refinement indicator,
not as a rigorous error bound.

## Interpretation limits

The Direct-WENO refinement rates come from two finite-grid intervals and do
not constitute a proof of continuum convergence. The sampled I2 correction is
a quadrature balance and not a certified temporal/spatial error decomposition.

## License

Code, data and figures are released under the MIT License. See `LICENSE`.

## Citation

Citation metadata are provided in `CITATION.cff` and `.zenodo.json`.
