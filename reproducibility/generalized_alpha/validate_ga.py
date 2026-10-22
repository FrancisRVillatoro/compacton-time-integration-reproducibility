#!/usr/bin/env python3
import csv
import math
import os
from pathlib import Path
import numpy as np
from scipy.sparse import eye

from ga_core import (
    AbsKPetrovGalerkin, ga_matrix, step_GA, step_GA_fast, step_M1,
)

ROOT = Path(__file__).resolve().parent
OUT = Path(os.environ.get('GA_OUT', ROOT / 'results' / 'validation_recheck'))
OUT.mkdir(parents=True, exist_ok=True)

class ScalarDecay:
    A = eye(1, format='csc')
    def G(self, u): return np.asarray(u)
    def F(self, u): return np.asarray(u)
    def JG(self, u): return self.A

rows=[]
for rho in [1.0, 0.9, 0.5, 0.25, 0.0]:
    prev=None
    for dt in [0.2,0.1,0.05,0.025,0.0125]:
        u=np.ones(1); v=-u.copy(); up=None
        n=round(1.0/dt)
        for _ in range(n):
            x,v,it,ok,nr=step_GA(ScalarDecay(),u,v,dt,rho,up,tol=1e-13,maxit=8)
            if not ok: raise RuntimeError((rho,dt,nr))
            up,u=u,x
        err=abs(float(u[0])-math.exp(-1.0))
        order=None if prev is None else math.log(prev/err,2)
        rows.append(dict(rho_inf=rho,dt=dt,error=err,order=order))
        prev=err

with (OUT/'validation_scalar_order.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)

# Require near-second order on last refinement for all rho.
for rho in [1.0,0.9,0.5,0.25,0.0]:
    rr=[r for r in rows if r['rho_inf']==rho]
    if rr[-1]['order'] < 1.90:
        raise SystemExit(f'SCALAR_ORDER_FAIL rho={rho} order={rr[-1]["order"]}')

# One-step agreement with exact 2x2 amplification matrix.
max_matrix_err=0.0
for rho in [1.0,0.9,0.5,0.25,0.0]:
  for dt in [0.2,0.05,0.0125]:
    u=np.array([1.0]); v=np.array([-1.0])
    x,v1,it,ok,nr=step_GA(ScalarDecay(),u,v,dt,rho,None,tol=1e-13,maxit=8)
    got=np.array([x[0],dt*v1[0]],dtype=complex)
    ref=ga_matrix(-dt,rho)@np.array([1.0,-dt],dtype=complex)
    max_matrix_err=max(max_matrix_err,float(np.max(np.abs(got-ref))))
if max_matrix_err>5e-12:
    raise SystemExit(f'MATRIX_FAIL {max_matrix_err:.3e}')

# Fast cyclic Newton and sparse Newton must agree on a nonlinear PG step.
rng=np.random.default_rng(20261019)
for rho in [0.0,0.5,0.9,1.0]:
    m=AbsKPetrovGalerkin(100,L=100,p=3,c0=1,mu=1e-7)
    U=0.03*rng.standard_normal(100)
    V=-m.F(U)
    a=step_GA(m,U,V,0.05,rho,None,tol=1e-11,maxit=20)
    b=step_GA_fast(m,U,V,0.05,rho,None,tol=1e-11,maxit=20)
    if not (a[3] and b[3]): raise SystemExit(f'FAST_DENSE_NEWTON_FAIL rho={rho}')
    diff=max(float(np.max(np.abs(a[0]-b[0]))),float(np.max(np.abs(a[1]-b[1]))))
    if diff>2e-9: raise SystemExit(f'FAST_DENSE_FAIL rho={rho} diff={diff:.3e}')

# rho_inf=1 must reproduce midpoint U trajectory.
m=AbsKPetrovGalerkin(100,L=100,p=3,c0=1,mu=1e-7)
x=np.arange(100)
U0=0.02*np.cos(2*np.pi*x/100)+0.01*np.sin(4*np.pi*x/100)
Ug=U0.copy(); Um=U0.copy(); V=-m.F(U0); Upg=None; Upm=None
for _ in range(20):
    Xg,V,itg,okg,nrg=step_GA_fast(m,Ug,V,0.05,1.0,Upg,tol=1e-11,maxit=20)
    Xm,itm,okm,nrm=step_M1(m,Um,0.05,Upm,tol=1e-11,maxit=20)
    if not (okg and okm): raise SystemExit('MIDPOINT_ENDPOINT_NEWTON_FAIL')
    Upg,Ug=Ug,Xg; Upm,Um=Um,Xm
midpoint_err=float(np.max(np.abs(Ug-Um)))
if midpoint_err>2e-9: raise SystemExit(f'MIDPOINT_ENDPOINT_FAIL {midpoint_err:.3e}')

summary = OUT/'validation_summary.txt'
summary.write_text(
    'GA_CORRECTED_VALIDATION_OK\n'
    f'max_matrix_error={max_matrix_err:.16e}\n'
    f'midpoint_endpoint_error={midpoint_err:.16e}\n'
    + '\n'.join(
        f'rho={rho:g} final_order={next(r["order"] for r in rows[::-1] if r["rho_inf"]==rho):.12f}'
        for rho in [1.0,0.9,0.5,0.25,0.0]
    ) + '\n'
)
print(summary.read_text(),end='')
