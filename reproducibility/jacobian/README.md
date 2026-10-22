# Jacobian experiment: closed traceability

## Provenance status

The archived numerical outputs and reports were recovered, but the original historical eigensolver drivers were not recovered byte-for-byte. The files in this directory are reconstructed authoritative drivers based on:

- the semidiscrete equations and Jacobian printed in the manuscript;
- the documented collision and mesh-refinement parameters;
- the archived collision and mesh-scaling reports;
- the frozen CSV outputs;
- a complete independent recomputation that reproduces those outputs within tight numerical tolerances.

This limitation is stated explicitly in `generated/JACOBIAN_VERIFICATION.txt` and must be preserved in any public release.

## Single numerical core

`jacobian_core.py` is the only source of the Petrov--Galerkin matrices, compacton profiles, Padé `R_{1,2}` trajectory step, matrix-free Jacobian, ARPACK settings, scalar stability functions, and exact 2x2 generalized-alpha amplification matrix.

The Jacobian is

```text
J(U)v = -A^{-1}[mu*h^2*D - c0*B
                + (B+C) diag(p*abs(U)^(p-1))] v.
```

No Jacobian or Figure 1/2 calculation calls a nonlinear generalized-alpha time-step implementation.

## Experiment A: collision trajectory

`recompute_collision_jacobian.py` starts from the documented compacton--anticompacton initial state and advances it with the Padé fractional-step method:

```text
p=3, (c1,c2)=(2,-1), c0=1,
(x1,x2)=(40,60), L=100,
N=1000, h=0.1, mu=1e-3, dt=0.025.
```

At `t=5,20,30`, it stores the state and computes the 60 largest-modulus Jacobian eigenvalues with a deterministic ARPACK starting vector. The recomputed values are

```text
t=5:  max |Im(dt*lambda)| = 2227.05533701
t=20: max |Im(dt*lambda)| = 2410.92393672
t=30: max |Im(dt*lambda)| = 2232.45161208
```

## Experiment B: fixed-state mesh scaling

`recompute_mesh_scaling.py` evaluates separated centers `(40,60)` and overlapping centers `(48,52)` on

```text
N = 250, 500, 1000, 2000, 4000,
```

using the 20 largest-modulus eigenvalues. Fits on the three finest grids give

```text
separated: |lambda_max| ~ h^-3.026910331362
            Padé amp   ~ h^ 3.026910533289
            GA(0) amp  ~ h^ 1.526861733781

overlap:   |lambda_max| ~ h^-3.056280962410
            Padé amp   ~ h^ 3.056281227752
            GA(0) amp  ~ h^ 1.542480487822
```

## Frozen and regenerated artifacts

- `frozen/`: archived CSV outputs used as immutable comparison targets.
- `generated/`: fresh outputs, configurations, stdout, and verification reports.
- `raw/jacobian_collision_snapshots.npz`: independently regenerated states at all three collision times.
- `reports/`: recovered narrative reports; these are provenance documents, not executable sources.

The verifier pairs complex conjugate eigenvalues by minimum-cost matching, because ARPACK is free to return conjugate pairs in either order. It checks spectra and filter factors to relative tolerance `5e-11`, mesh quantities to `5e-11`, and fitted coefficients to absolute tolerance `1e-10`.

## Commands

Quick comparison against the already regenerated outputs:

```bash
python verify_jacobian_traceability.py
```

Complete recomputation from initial conditions and fixed physical states:

```bash
python verify_jacobian_traceability.py --full
```
