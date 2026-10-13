# Computational environment

The Direct-WENO production runs on Picasso used:

- Python 3.9.13
- NumPy 1.24.4
- SciPy 1.9.1
- Numba 0.57.1
- Matplotlib 3.5.2

The verification scripts are written for standard recent Python/NumPy/SciPy
installations. Exact production trajectories are included under `raw/` and are
checked by SHA-256 before derived quantities are regenerated.
