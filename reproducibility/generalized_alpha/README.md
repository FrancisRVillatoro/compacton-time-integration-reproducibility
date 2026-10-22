# Corrected generalized-alpha campaign

## Defect and correction

The semidiscrete convention is

```text
A Udot + G(U) = 0.
```

For the first-order generalized-alpha method, `V_n` approximates `Udot_n` and must remain an independent algorithmic state. After eliminating `V_{n+1}`, the Newton residual for `X=U_{n+1}` is

```text
R(X) = a A (X-U_n) + c A V_n
       + G((1-alpha_f)U_n + alpha_f X),
a = alpha_m/(gamma*dt),
c = 1-alpha_m/gamma.
```

The superseded implementation replaced `A V_n` by `-G(U_n)` at every step. That substitution is valid only for an equilibrium-consistent initial state and changes the method after the first step whenever `c != 0`. `ga_core.py` implements the corrected residual in both sparse and cyclic-pentadiagonal Newton solvers.

## Validation

`validate_ga.py` checks:

- one-step agreement with the exact 2x2 amplification matrix;
- equivalence of sparse and cyclic-pentadiagonal Newton solves;
- recovery of implicit midpoint at `rho_inf=1`;
- approximately second-order scalar convergence throughout the family.

The audited values are:

```text
max_matrix_error          = 1.1102230246251565e-16
midpoint_endpoint_error   = 3.4694469519536142e-18
final scalar orders       = 2.00007, 2.00007, 1.99807, 1.98785, 1.93830
                            for rho_inf = 1, 0.9, 0.5, 0.25, 0.
```

## Picasso campaign

The corrected production campaign contains:

- 396/396 successful microbank runs;
- 62/64 successful fine-boundary runs;
- 7/7 successful independent `N=1000, c1=2` benchmark runs.

The two fine-boundary failures are the nondissipative `rho_inf=1` controls at `N=1000`, not failures of the damped family. All 47 retained Slurm standard-error files are empty.

`results/fig3_right_corrected.csv` is the sole numerical source of the new Figure 3. The coarse `N=100` summary `fig3_left_corrected.csv` is retained only as provenance and is not plotted.

## Traceability chain

```text
Picasso group CSVs
    -> microbank_corrected_all.csv / boundary_corrected_all.csv
    -> fig3_right_corrected.csv
    -> figure_sources/data/fig3_controlled_damping_fine_corrected.csv
    -> figure_sources/fig03_controlled_damping_tradeoff.py
    -> manuscript/figures/fig3_controlled_damping_tradeoff_{left,right}.pdf
    -> manuscript/main.tex
```

The verifier checks every arrow, the Picasso SHA-256 manifest, aggregate reconstruction, status counts, key benchmark values, absence of superseded manuscript claims, and a fresh validation run.

## Commands

```bash
python verify_generalized_alpha_traceability.py --full
```

The `--full` option reruns `validate_ga.py`; it does not repeat the full Picasso production campaign. The immutable production outputs and scheduler provenance are audited in place.
