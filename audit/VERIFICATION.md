# Verification summary

The public release verifies the raw Direct-WENO trajectories by SHA-256,
re-fits the N=3200 classical and absolute core models, recomputes the absolute
800--1600--3200 spatial-refinement diagnostics, re-fits all stored
Petrov--Galerkin core states used in the final comparison, recomputes the
sampled I2 balance from the stored trajectories, rebuilds the final comparison
table, and regenerates the spatial-convergence and accuracy--cost figures.

Run:

```bash
python scripts/freeze_check.py
```

A successful verification ends with:

```text
FREEZE_CHECK_OK
```

Interpretation limits: the two Direct-WENO refinement intervals are finite-grid
diagnostics rather than a proof of continuum convergence, and the sampled I2
balance is a quadrature diagnostic rather than a certified separation of
spatial and temporal error.
