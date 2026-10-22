# Jacobian mesh-scaling test

## Controlled design
The same physical profiles are sampled on successively refined periodic grids.
This deliberately avoids conflating mesh scaling with changes in the evolved PDE state.

Two profiles are used:
1. `separated`: compacton at x=40 and anticompacton at x=60.
2. `overlap`: collision-like overlapping profiles centered at x=48 and x=52.

Common parameters: p=3, c1=2, c2=-1, c0=1, L=100, mu=1e-3, dt=0.025.
For each mesh, the 20 largest-magnitude eigenvalues of the exact semi-discrete Jacobian are computed; the stiffest mode is then passed through M1, M2, M3, and generalized-alpha amplification factors.

## Fine-grid fitted exponents (h <= 0.1)

| state | |lambda_max| vs h | M2 amplification vs h | GA(rho_inf=0) amplification vs h |
|---|---:|---:|---:|
| separated | -3.026910 | 3.026911 | 1.526862 |
| overlap | -3.056281 | 3.056281 | 1.542480 |

The asymptotic predictions are -3, +3, and +3/2, respectively.

## Interpretation
The numerical exponents are very close to the predicted values:
- the real semi-discrete Jacobian stiffness grows as h^-3;
- M2 damps the stiffest mode as h^3 because S2(z) ~ 2/z;
- maximally damped generalized-alpha only damps as approximately h^(3/2), consistent with rho(A_0(z)) = O(|z|^-1/2);
- M1 and M3 tend to unit amplification on the same modes.

The same scaling occurs for both a separated physical state and a controlled overlap state, showing that the h^-3 mechanism is tied to the nonlinear dispersive C operator rather than to one accidental collision snapshot.

Caveat: the overlap profile is a controlled diagnostic state, not an exact evolved collision solution. The preceding real-Jacobian snapshot experiment provides the complementary evidence that the same large-imaginary spectrum is present in an actual M2 collision trajectory.
