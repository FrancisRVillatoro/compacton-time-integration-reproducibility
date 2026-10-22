#!/usr/bin/env python3
"""Correct first-order generalized-alpha implementation for the |K|(p,p) PG system.

System convention
-----------------
    A Udot + G(U) = 0,
    G(U) = (mu*h^2*D-c0*B)U + (B+C)g(U),
    g(U)=|U|^(p-1) U.

Generalized-alpha convention
----------------------------
Let V_n approximate Udot_n.  The method enforces

    A V_{n+alpha_m} + G(U_{n+alpha_f}) = 0,

with
    U_{n+1} = U_n + dt*((1-gamma)V_n + gamma V_{n+1}).

After eliminating V_{n+1}, the nonlinear residual for X=U_{n+1} is

    R(X) = a A (X-U_n) + c A V_n + G((1-alpha_f)U_n+alpha_f X),
    a = alpha_m/(gamma*dt), c = 1-alpha_m/gamma.

The historical recovered implementation replaced A V_n by -G(U_n) at every
step.  That replacement is valid only for a consistent initial state (and is
irrelevant at rho_inf=1, where c=0).  This file keeps V_n explicitly.
"""
from __future__ import annotations

import math
import time
from dataclasses import dataclass
from typing import Optional

import numpy as np
from scipy.linalg import solve_banded
from scipy.sparse import bmat, coo_matrix
from scipy.sparse.linalg import splu, spsolve


def periodic_matrix(N, offsets, coeffs):
    i = np.arange(N)
    rr, cc, dd = [], [], []
    for off, c in zip(offsets, coeffs):
        rr.append(i)
        cc.append((i + off) % N)
        dd.append(np.full(N, c))
    return coo_matrix(
        (np.concatenate(dd), (np.concatenate(rr), np.concatenate(cc))),
        shape=(N, N),
    ).tocsc()


class AbsKPetrovGalerkin:
    def __init__(self, N, L=100.0, p=3.0, c0=1.0, mu=1e-7):
        self.N = int(N)
        self.L = float(L)
        self.h = self.L / self.N
        self.p = float(p)
        self.c0 = float(c0)
        self.mu = float(mu)
        h = self.h
        self.A = periodic_matrix(
            self.N, [-2, -1, 0, 1, 2],
            [1/120, 26/120, 66/120, 26/120, 1/120],
        )
        self.B = periodic_matrix(
            self.N, [-2, -1, 1, 2],
            [-1/(24*h), -10/(24*h), 10/(24*h), 1/(24*h)],
        )
        self.C = periodic_matrix(
            self.N, [-2, -1, 1, 2],
            [-1/(2*h**3), 2/(2*h**3), -2/(2*h**3), 1/(2*h**3)],
        )
        self.D = periodic_matrix(
            self.N, [-2, -1, 0, 1, 2],
            [1/h**4, -4/h**4, 6/h**4, -4/h**4, 1/h**4],
        )
        self.BC = (self.B + self.C).tocsc()
        self.Lop = (self.mu*h*h*self.D - self.c0*self.B).tocsc()
        self.Alu = splu(self.A)

    def g(self, U):
        return np.abs(U)**(self.p - 1.0) * U

    def gp(self, U):
        return self.p * np.abs(U)**(self.p - 1.0)

    def G(self, U):
        return self.Lop @ U + self.BC @ self.g(U)

    def F(self, U):
        return self.Alu.solve(self.G(U))

    def JG(self, U):
        return (self.Lop + self.BC.multiply(self.gp(U))).tocsc()


def ga_parameters(rho_inf):
    r = float(rho_inf)
    if not (0.0 <= r <= 1.0):
        raise ValueError("rho_inf must lie in [0,1]")
    alpha_m = (3.0-r) / (2.0*(1.0+r))
    alpha_f = 1.0 / (1.0+r)
    gamma = 0.5 + alpha_m - alpha_f
    return alpha_m, alpha_f, gamma


def ga_matrix(z, rho_inf):
    """Exact 2x2 amplification matrix for state (u, dt*udot)."""
    r = float(rho_inf)
    H = 3.0 + 2.0*r - r*r - 2.0*z
    return np.array([
        [(3.0 + 2.0*r - r*r + 2.0*r*z)/H, (1.0-r*r)/H],
        [2.0*z*(1.0+r)**2/H, (1.0 - 2.0*r - 3.0*r*r + 2.0*r*z)/H],
    ], dtype=complex)


def damped_newton(residual_jacobian, x0, tol=1e-10, maxit=20):
    x = np.asarray(x0, dtype=float).copy()
    last = np.inf
    for it in range(maxit):
        R, J = residual_jacobian(x)
        nrm = float(np.linalg.norm(R, np.inf))
        last = nrm
        if nrm < tol:
            return x, it, True, nrm
        dx = spsolve(J, -R)
        best = (np.inf, None)
        for s in (1.0, .5, .25, .125, .0625, .03125, .015625, .0078125):
            xt = x + s*dx
            if not np.all(np.isfinite(xt)):
                continue
            rt = float(np.linalg.norm(residual_jacobian(xt)[0], np.inf))
            if rt < best[0]:
                best = (rt, xt)
            if rt < nrm:
                break
        if best[1] is None:
            return x, it, False, nrm
        x = best[1]
    R, _ = residual_jacobian(x)
    return x, maxit, False, float(np.linalg.norm(R, np.inf))


def step_M1(model, U, dt, Uprev=None, tol=1e-10, maxit=20):
    guess = U - dt*model.F(U) if Uprev is None else U + (U-Uprev)
    def RJ(X):
        M = 0.5*(X+U)
        R = model.A @ ((X-U)/dt) + model.G(M)
        J = model.A/dt + 0.5*model.JG(M)
        return R, J.tocsc()
    return damped_newton(RJ, guess, tol=tol, maxit=maxit)


def step_GA(model, U, V, dt, rho_inf, Uprev=None, tol=1e-10, maxit=20):
    """Correct generalized-alpha step; V is retained in the equilibrium residual."""
    am, af, gam = ga_parameters(rho_inf)
    a = am/(gam*dt)
    c = 1.0-am/gam
    AV = model.A @ V
    guess = U + dt*V if Uprev is None else U + (U-Uprev)
    def RJ(X):
        Ua = (1.0-af)*U + af*X
        R = a*(model.A @ (X-U)) + c*AV + model.G(Ua)
        J = a*model.A + af*model.JG(Ua)
        return R, J.tocsc()
    X, it, ok, nrm = damped_newton(RJ, guess, tol=tol, maxit=maxit)
    Vnew = (X-U)/(gam*dt) - (1.0-gam)/gam*V
    return X, Vnew, it, ok, nrm


def _coeffs(model, a, af, q):
    h, mu, c0 = model.h, model.mu, model.c0
    A = {-2:1/120, -1:26/120, 0:66/120, 1:26/120, 2:1/120}
    B = {-2:-1/(24*h), -1:-10/(24*h), 0:0.0, 1:10/(24*h), 2:1/(24*h)}
    C = {-2:-1/(2*h**3), -1:2/(2*h**3), 0:0.0, 1:-2/(2*h**3), 2:1/(2*h**3)}
    D = {-2:1/h**4, -1:-4/h**4, 0:6/h**4, 1:-4/h**4, 2:1/h**4}
    out = {}
    for o in (-2, -1, 0, 1, 2):
        lop = mu*h*h*D[o] - c0*B[o]
        bc = B[o] + C[o]
        out[o] = a*A[o] + af*lop + af*bc*q
    return out


def solve_cyclic_pentadiagonal(model, a, af, q, rhs):
    N = model.N
    coeff = _coeffs(model, a, af, q)
    ab = np.zeros((5, N), dtype=float)
    for o in (-2, -1, 0, 1, 2):
        vals = coeff[o]
        js = np.arange(N)
        ii = js-o
        m = (ii >= 0) & (ii < N)
        ab[2-o, js[m]] = vals[js[m]]
    cols = np.array([0, 1, N-2, N-1], dtype=int)
    cmap = {int(c): k for k, c in enumerate(cols)}
    W = np.zeros((N, 4), dtype=float)
    for o in (-2, -1, 1, 2):
        for i in range(N):
            j0 = i+o
            if j0 < 0 or j0 >= N:
                j = j0 % N
                W[i, cmap[j]] += coeff[o][j]
    Bmat = np.column_stack([rhs, W])
    sol = solve_banded((2, 2), ab, Bmat, check_finite=False)
    y, Y = sol[:, 0], sol[:, 1:]
    M = np.eye(4) + Y[cols, :]
    corr = np.linalg.solve(M, y[cols])
    return y - Y @ corr


def step_GA_fast(model, U, V, dt, rho_inf, Uprev=None, tol=1e-10, maxit=20):
    """Correct GA step using the cyclic-pentadiagonal Newton solve."""
    am, af, gam = ga_parameters(rho_inf)
    a = am/(gam*dt)
    c = 1.0-am/gam
    AV = model.A @ V
    x = U + dt*V if Uprev is None else U + (U-Uprev)
    for it in range(maxit):
        Ua = (1.0-af)*U + af*x
        R = a*(model.A @ (x-U)) + c*AV + model.G(Ua)
        nr = float(np.linalg.norm(R, np.inf))
        if nr < tol:
            Vnew = (x-U)/(gam*dt) - (1.0-gam)/gam*V
            return x, Vnew, it, True, nr
        q = model.gp(Ua)
        dx = solve_cyclic_pentadiagonal(model, a, af, q, -R)
        bestnr, bestx = np.inf, None
        for s in (1.0, .5, .25, .125, .0625, .03125, .015625, .0078125):
            xt = x+s*dx
            if not np.all(np.isfinite(xt)):
                continue
            Uat = (1.0-af)*U + af*xt
            Rt = a*(model.A @ (xt-U)) + c*AV + model.G(Uat)
            nrt = float(np.linalg.norm(Rt, np.inf))
            if nrt < bestnr:
                bestnr, bestx = nrt, xt
            if nrt < nr:
                break
        if bestx is None:
            break
        x = bestx
    Ua = (1.0-af)*U + af*x
    R = a*(model.A @ (x-U)) + c*AV + model.G(Ua)
    nr = float(np.linalg.norm(R, np.inf))
    Vnew = (x-U)/(gam*dt) - (1.0-gam)/gam*V
    return x, Vnew, maxit, False, nr


def step_M2(model, U, dt, Uprev=None, tol=1e-10, maxit=20):
    FU = model.F(U)
    W = U-dt*FU/3.0
    N = model.N
    X = U-dt*FU if Uprev is None else U+(U-Uprev)
    V = X+dt*model.F(X)/4.0
    Y0 = np.r_[X, V]
    def RJ(Y):
        X, V = Y[:N], Y[N:]
        R1 = model.A @ ((V-X)/dt) - 0.25*model.G(X)
        R2 = model.A @ ((X-W)/dt) + (2.0/3.0)*model.G(V)
        J = bmat([
            [-model.A/dt-0.25*model.JG(X), model.A/dt],
            [ model.A/dt, (2.0/3.0)*model.JG(V)],
        ], format="csc")
        return np.r_[R1, R2], J
    Y, it, ok, nrm = damped_newton(RJ, Y0, tol=tol, maxit=maxit)
    return Y[:N], it, ok, nrm


def compacton(x, c, x0, p=3.0, L=100.0):
    ell = p*np.pi/(p-1.0)
    amp = 2.0*abs(c)*p/(p+1.0)
    beta = (p-1.0)/(2.0*p)
    expo = 1.0/(p-1.0)
    xi = (x-x0+L/2.0) % L - L/2.0
    u = np.zeros_like(x, dtype=float)
    m = np.abs(xi) <= ell
    u[m] = np.sign(c)*(amp*np.cos(beta*xi[m])**2)**expo
    return u


def rho_hw(U):
    z = np.fft.fft(U)
    theta = 2*np.pi*np.fft.fftfreq(U.size)
    mask = np.abs(theta) >= np.pi/2.0
    den = float(np.sum(np.abs(z)**2))
    return float(np.sum(np.abs(z[mask])**2)/den) if den else 0.0


def I2(U, h, p=3.0):
    return float(h*np.sum(np.abs(U)**(p+1.0)/(p+1.0)))


@dataclass(frozen=True)
class Case:
    N: int
    c1: float
    mu: float
    dt: float
    rho_inf: Optional[float]
    method: str = "GA"
    L: float = 100.0
    p: float = 3.0
    c0: float = 1.0
    c2: float = -1.0
    x1: float = 44.0
    x2: float = 56.0
    T: Optional[float] = None
    tpost: Optional[float] = None
    blow: float = 50.0
    newton_tol: float = 1e-10
    newton_maxit: int = 20


def run_case(case: Case):
    model = AbsKPetrovGalerkin(case.N, case.L, case.p, case.c0, case.mu)
    x = np.arange(case.N)*(case.L/case.N)
    U = compacton(x, case.c1, case.x1, case.p, case.L) + compacton(x, case.c2, case.x2, case.p, case.L)
    Up = None
    V = -model.F(U) if case.method == "GA" else None

    tc = (case.x2-case.x1)/(case.c1-case.c0)
    T = case.T if case.T is not None else tc+15.0
    tpost = case.tpost if case.tpost is not None else tc+3.0
    nsteps = int(math.ceil(T/case.dt))
    Ttarget = nsteps*case.dt

    I20 = I2(U, model.h, case.p)
    pr = rho_hw(U)
    pt = 0.0
    auc = 0.0
    mr = pr
    mx = float(np.max(np.abs(U)))
    its = []
    status = "ok"
    fail_t = float("nan")
    final_nr = float("nan")
    st = time.perf_counter()

    for n in range(nsteps):
        try:
            if case.method == "GA":
                X, V, it, ok, nr = step_GA_fast(
                    model, U, V, case.dt, case.rho_inf, Up,
                    tol=case.newton_tol, maxit=case.newton_maxit,
                )
            elif case.method == "M1":
                X, it, ok, nr = step_M1(model, U, case.dt, Up, tol=case.newton_tol, maxit=case.newton_maxit)
            elif case.method == "M2":
                X, it, ok, nr = step_M2(model, U, case.dt, Up, tol=case.newton_tol, maxit=case.newton_maxit)
            else:
                raise ValueError(case.method)
        except Exception:
            status = "solver_exception"
            fail_t = n*case.dt
            break
        Up, U = U, X
        its.append(it)
        final_nr = nr
        tt = (n+1)*case.dt
        if not ok:
            status = "newton_failure"
            fail_t = tt
            break
        if not np.all(np.isfinite(U)):
            status = "nonfinite"
            fail_t = tt
            break
        uu = float(np.max(np.abs(U)))
        rr = rho_hw(U)
        mx = max(mx, uu)
        mr = max(mr, rr)
        if uu > case.blow:
            status = "amplitude_blowup"
            fail_t = tt
            break
        if tt > tpost:
            left = max(pt, tpost)
            if tt > left:
                if pt < tpost:
                    ff = (tpost-pt)/(tt-pt)
                    rl = pr + ff*(rr-pr)
                else:
                    rl = pr
                auc += 0.5*(rl+rr)*(tt-left)
        pt, pr = tt, rr

    finite = bool(np.all(np.isfinite(U)))
    return dict(
        method=case.method,
        rho_inf=(float("nan") if case.rho_inf is None else case.rho_inf),
        N=case.N,
        dx=case.L/case.N,
        c1=case.c1,
        c2=case.c2,
        x1=case.x1,
        x2=case.x2,
        mu=case.mu,
        dt=case.dt,
        collision_time=tc,
        T_target=Ttarget,
        tpost=tpost,
        newton_tol=case.newton_tol,
        newton_maxit=case.newton_maxit,
        status=status,
        fail_t=fail_t,
        rho_hw_max=mr,
        rho_hw_final=(rho_hw(U) if finite else float("nan")),
        rho_hw_post_AUC=(auc if status == "ok" else float("nan")),
        I2_drift_rel=((I2(U, model.h, case.p)-I20)/I20 if finite else float("nan")),
        max_abs_u=mx,
        mean_newton=(float(np.mean(its)) if its else float("nan")),
        median_newton=(float(np.median(its)) if its else float("nan")),
        max_newton=(int(np.max(its)) if its else 0),
        final_residual_inf=final_nr,
        steps_done=len(its),
        wall_s=time.perf_counter()-st,
    )
