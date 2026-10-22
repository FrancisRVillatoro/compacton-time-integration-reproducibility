# Real-Jacobian collision experiment

## Setup
- Equation: absolute |K|(p,p), p = 3.
- Interaction: compacton–anticompacton (CA).
- c1 = 2, c2 = -1, c0 = 1.
- Initial centers: x1 = 40, x2 = 60, L = 100.
- Petrov–Galerkin grid: N = 1000, Δx = 0.1.
- Artificial viscosity: μ = 0.001.
- M2 time step: Δt = 0.025.
- Snapshots: t = 5 (pre-collision), 20 (center collision), 30 (post-collision).

For each snapshot, the actual ODE Jacobian
J(U) = -A^(-1)[μ Δx^2 D - c0 B + (B+C) diag(p |U|^(p-1))]
was formed as a matrix-free linear operator. ARPACK was used to compute the 60 eigenvalues of largest magnitude. Thus the spectral plots diagnose the high-frequency/stiff portion of the real Jacobian, not the full spectrum.

## Main numerical result

| t | stage | max |u| | max |Im(Δt λ)| | median M1 | median M2 | median M3 | median GA(0.5) | median GA(0) |
|---:|---|---:|---:|---:|---:|---:|---:|---:|
| 5 | pre-collision | 1.73204 | 2227.06 | 0.999999 | 0.0025395 | 0.999999 | 0.52367 | 0.0261054 |
| 20 | collision | 1.80533 | 2410.92 | 0.999999 | 0.0028388 | 0.999999 | 0.525063 | 0.0276622 |
| 30 | post-collision | 1.73411 | 2232.45 | 0.999999 | 0.00255354 | 0.999999 | 0.523735 | 0.026178 |

## Interpretation
The real high-frequency Jacobian spectrum is strongly dispersive-stiff already before collision and becomes still more extreme at the collision. M1 and M3 transmit these modes with amplification essentially one. M2 damps the same modes by roughly 2.5e-3 to 2.8e-3 in the median. Even maximally damped generalized-alpha (rho_infty=0) is roughly one order of magnitude less aggressive on these modes.

This supports the mechanism:
nonzero interaction background -> q(B+C), especially q C ~ Δx^-3 -> large imaginary semi-discrete eigenvalues -> selective temporal filtering.

Caveat: these are the 60 largest-|lambda| eigenvalues only. A full-spectrum/pseudospectral analysis and refinement in Δx remain to be done.
