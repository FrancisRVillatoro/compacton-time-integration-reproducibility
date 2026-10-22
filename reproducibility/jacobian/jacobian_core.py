#!/usr/bin/env python3
"""Numerical core for the Jacobian experiments in the compacton manuscript.

The semidiscrete system is

    A U_t + G(U) = 0,
    G(U) = (mu*h**2*D-c0*B)U + (B+C)g(U),
    g(U) = |U|**(p-1) U.

The Jacobian of the ODE right-hand side is represented matrix-free as

    J(U)v = -A^{-1} [mu*h**2*D-c0*B+(B+C)diag(g'(U))] v.

This module is the single source of truth for the Petrov--Galerkin matrices,
the M2 trajectory used in the collision experiment, and the exact 2x2
first-order generalized-alpha amplification matrix.  It never calls a
nonlinear generalized-alpha time-step implementation.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Iterable

import numpy as np
from scipy.sparse import bmat, coo_matrix
from scipy.sparse.linalg import LinearOperator, eigs, splu, spsolve


def periodic_matrix(n: int, offsets: Iterable[int], coeffs: Iterable[float]):
    """Return a periodic sparse matrix with the requested stencil."""
    idx = np.arange(n)
    rows, cols, data = [], [], []
    for off, coeff in zip(offsets, coeffs):
        rows.append(idx)
        cols.append((idx + int(off)) % n)
        data.append(np.full(n, float(coeff)))
    return coo_matrix(
        (np.concatenate(data), (np.concatenate(rows), np.concatenate(cols))),
        shape=(n, n),
    ).tocsc()


class AbsKPetrovGalerkin:
    """Five-point periodic Petrov--Galerkin discretization of absolute K(p,p)."""

    def __init__(
        self,
        n: int,
        length: float = 100.0,
        p: float = 3.0,
        c0: float = 1.0,
        mu: float = 1.0e-3,
    ) -> None:
        self.N = int(n)
        self.L = float(length)
        self.h = self.L / self.N
        self.p = float(p)
        self.c0 = float(c0)
        self.mu = float(mu)
        h = self.h

        self.A = periodic_matrix(
            self.N,
            (-2, -1, 0, 1, 2),
            (1 / 120, 26 / 120, 66 / 120, 26 / 120, 1 / 120),
        )
        self.B = periodic_matrix(
            self.N,
            (-2, -1, 1, 2),
            (-1 / (24 * h), -10 / (24 * h), 10 / (24 * h), 1 / (24 * h)),
        )
        self.C = periodic_matrix(
            self.N,
            (-2, -1, 1, 2),
            (-1 / (2 * h**3), 2 / (2 * h**3), -2 / (2 * h**3), 1 / (2 * h**3)),
        )
        self.D = periodic_matrix(
            self.N,
            (-2, -1, 0, 1, 2),
            (1 / h**4, -4 / h**4, 6 / h**4, -4 / h**4, 1 / h**4),
        )
        self.BC = (self.B + self.C).tocsc()
        self.Lop = (self.mu * h * h * self.D - self.c0 * self.B).tocsc()
        self.A_lu = splu(self.A)

    def g(self, u: np.ndarray) -> np.ndarray:
        return np.abs(u) ** (self.p - 1.0) * u

    def gp(self, u: np.ndarray) -> np.ndarray:
        return self.p * np.abs(u) ** (self.p - 1.0)

    def G(self, u: np.ndarray) -> np.ndarray:
        return self.Lop @ u + self.BC @ self.g(u)

    def F(self, u: np.ndarray) -> np.ndarray:
        return self.A_lu.solve(self.G(u))

    def JG(self, u: np.ndarray):
        # scipy sparse broadcasting with a length-N vector scales columns,
        # hence this is (B+C) diag(g'(u)).
        return (self.Lop + self.BC.multiply(self.gp(u))).tocsc()

    def jacobian_operator(self, u: np.ndarray) -> LinearOperator:
        jg = self.JG(u)

        def matvec(v: np.ndarray) -> np.ndarray:
            return -self.A_lu.solve(jg @ v)

        return LinearOperator((self.N, self.N), matvec=matvec, dtype=float)


def compacton(
    x: np.ndarray,
    speed: float,
    center: float,
    p: float = 3.0,
    length: float = 100.0,
) -> np.ndarray:
    """Periodic compacton/anticompacton profile used in both experiments."""
    ell = p * np.pi / (p - 1.0)
    amplitude_base = 2.0 * abs(speed) * p / (p + 1.0)
    beta = (p - 1.0) / (2.0 * p)
    exponent = 1.0 / (p - 1.0)
    xi = (x - center + length / 2.0) % length - length / 2.0
    out = np.zeros_like(x, dtype=float)
    support = np.abs(xi) <= ell
    out[support] = (
        np.sign(speed)
        * (amplitude_base * np.cos(beta * xi[support]) ** 2) ** exponent
    )
    return out


def damped_newton(
    residual_jacobian: Callable[[np.ndarray], tuple[np.ndarray, object]],
    x0: np.ndarray,
    tol: float = 1.0e-9,
    maxit: int = 40,
) -> tuple[np.ndarray, int, bool, float]:
    """Newton solve with the line-search sequence used in the archived run."""
    x = np.asarray(x0, dtype=float).copy()
    for iteration in range(maxit):
        residual, jacobian = residual_jacobian(x)
        norm_inf = float(np.linalg.norm(residual, np.inf))
        if norm_inf < tol:
            return x, iteration, True, norm_inf
        delta = spsolve(jacobian, -residual)
        best_norm, best_x = np.inf, None
        for scale in (1.0, 0.5, 0.25, 0.125, 0.0625, 0.03125, 0.015625, 0.0078125):
            trial = x + scale * delta
            if not np.all(np.isfinite(trial)):
                continue
            trial_norm = float(np.linalg.norm(residual_jacobian(trial)[0], np.inf))
            if trial_norm < best_norm:
                best_norm, best_x = trial_norm, trial
            if trial_norm < norm_inf:
                break
        if best_x is None:
            return x, iteration, False, norm_inf
        x = best_x
    residual, _ = residual_jacobian(x)
    return x, maxit, False, float(np.linalg.norm(residual, np.inf))


def step_m2(
    model: AbsKPetrovGalerkin,
    u: np.ndarray,
    dt: float,
    u_previous: np.ndarray | None = None,
    tol: float = 1.0e-9,
    maxit: int = 40,
) -> tuple[np.ndarray, int, bool, float]:
    """One Padé R_{1,2} fractional-step (M2) update."""
    f_u = model.F(u)
    w = u - dt * f_u / 3.0
    n = model.N
    x = u - dt * f_u if u_previous is None else u + (u - u_previous)
    v = x + dt * model.F(x) / 4.0
    y0 = np.r_[x, v]

    def residual_jacobian(y: np.ndarray):
        x_stage, v_stage = y[:n], y[n:]
        r1 = model.A @ ((v_stage - x_stage) / dt) - 0.25 * model.G(x_stage)
        r2 = model.A @ ((x_stage - w) / dt) + (2.0 / 3.0) * model.G(v_stage)
        jacobian = bmat(
            [
                [-model.A / dt - 0.25 * model.JG(x_stage), model.A / dt],
                [model.A / dt, (2.0 / 3.0) * model.JG(v_stage)],
            ],
            format="csc",
        )
        return np.r_[r1, r2], jacobian

    y, iteration, ok, norm_inf = damped_newton(
        residual_jacobian, y0, tol=tol, maxit=maxit
    )
    return y[:n], iteration, ok, norm_inf


def deterministic_v0(n: int) -> np.ndarray:
    """Fixed nonzero ARPACK starting vector, independent of global RNG state."""
    j = np.arange(n, dtype=float)
    return (
        np.sin(2.0 * np.pi * j / n)
        + 0.123 * np.cos(6.0 * np.pi * j / n)
        + 0.071 * np.sin(14.0 * np.pi * j / n)
    )


def largest_magnitude_eigenvalues(
    model: AbsKPetrovGalerkin,
    u: np.ndarray,
    k: int,
    tol: float = 1.0e-10,
    maxiter: int = 50_000,
) -> np.ndarray:
    """Compute and sort the k largest-modulus Jacobian eigenvalues."""
    ncv = min(model.N, max(2 * k + 1, 3 * k))
    values = eigs(
        model.jacobian_operator(u),
        k=k,
        which="LM",
        tol=tol,
        maxiter=maxiter,
        ncv=ncv,
        v0=deterministic_v0(model.N),
        return_eigenvectors=False,
    )
    # The archived CSV orders the selected stiff sector from smaller to larger
    # modulus.  Conjugate-pair ordering is not scientifically meaningful.
    return values[np.argsort(np.abs(values))]


def s1(z: complex | np.ndarray):
    return (1.0 + z / 2.0) / (1.0 - z / 2.0)


def s2(z: complex | np.ndarray):
    return (1.0 + z / 3.0) / (1.0 - 2.0 * z / 3.0 + z * z / 6.0)


_BETA1 = (2.0 + 2.0 ** (1.0 / 3.0) + 2.0 ** (-1.0 / 3.0)) / 3.0
_BETA2 = 1.0 - 2.0 * _BETA1


def s3(z: complex | np.ndarray):
    numerator = (
        (1.0 + _BETA1 * z / 2.0)
        * (1.0 + _BETA2 * z / 2.0)
        * (1.0 + _BETA1 * z / 2.0)
    )
    denominator = (
        (1.0 - _BETA1 * z / 2.0)
        * (1.0 - _BETA2 * z / 2.0)
        * (1.0 - _BETA1 * z / 2.0)
    )
    return numerator / denominator


def ga_matrix(z: complex, rho_inf: float) -> np.ndarray:
    """Exact first-order generalized-alpha 2x2 amplification matrix."""
    r = float(rho_inf)
    hden = 3.0 + 2.0 * r - r * r - 2.0 * z
    return np.array(
        [
            [
                (3.0 + 2.0 * r - r * r + 2.0 * r * z) / hden,
                (1.0 - r * r) / hden,
            ],
            [
                2.0 * z * (1.0 + r) ** 2 / hden,
                (1.0 - 2.0 * r - 3.0 * r * r + 2.0 * r * z) / hden,
            ],
        ],
        dtype=complex,
    )


def ga_rho(z: complex, rho_inf: float) -> float:
    return float(np.max(np.abs(np.linalg.eigvals(ga_matrix(z, rho_inf)))))


def amplification_row(value: complex, dt: float) -> dict[str, float]:
    z = dt * value
    return {
        "z_real": float(np.real(z)),
        "z_imag": float(np.imag(z)),
        "abs_z": float(abs(z)),
        "M1_amp": float(abs(s1(z))),
        "M2_amp": float(abs(s2(z))),
        "M3_amp": float(abs(s3(z))),
        "GA_rho_0p5": ga_rho(z, 0.5),
        "GA_rho_0": ga_rho(z, 0.0),
    }


@dataclass(frozen=True)
class CollisionConfig:
    n: int = 1000
    length: float = 100.0
    p: float = 3.0
    c0: float = 1.0
    c1: float = 2.0
    c2: float = -1.0
    x1: float = 40.0
    x2: float = 60.0
    mu: float = 1.0e-3
    dt: float = 0.025
    snapshot_times: tuple[float, ...] = (5.0, 20.0, 30.0)
    newton_tol: float = 1.0e-9
    newton_maxit: int = 40
    eigen_count: int = 60
