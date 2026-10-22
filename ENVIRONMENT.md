# Computational environments

## Direct-WENO production runs retained from v1.0.0

- Python 3.9.13
- NumPy 1.24.4
- SciPy 1.9.1
- Numba 0.57.1
- Matplotlib 3.5.2

The exact production trajectories are included under `raw/` and checked by
SHA-256 before derived quantities are regenerated.

## Corrected generalized-alpha and Jacobian campaign on Picasso

- Python 3.14.3
- NumPy 2.3.5
- SciPy 1.17.0
- one BLAS/OpenMP thread per task

The retained validation reports record agreement with the analytical 2x2
amplification matrix to `1.11e-16` and recovery of implicit midpoint at
`rho_inf=1` to `3.47e-18`.

## Robustness banks on Picasso

- MATLAB R2024b
- `-singleCompThread`
- SLURM `cpu_parti`

The production launchers, dispatcher, integrator drivers, case generators and
per-run summary writers are included under `provenance/robustness_banks/`.

## Verification environment

The Python verification scripts target standard recent CPython installations
with the package lower bounds in `requirements.txt`. Set BLAS thread counts to
one for the full Jacobian recomputation when deterministic resource use is
important.
