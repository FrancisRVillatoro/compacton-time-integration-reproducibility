#!/usr/bin/env python3
"""Sampling control on a FIXED TR-BDF2 numerical trajectory.

Integrate delta using accepted nodes + actual trapezoidal-stage nodes, and
on the quadratic interpolant through U_n, Y_(n+gamma), U_(n+1). The latter
controls quadrature on the specified interpolant only: it is NOT a rigorous
error bound for the exact semidiscrete flow or a replacement trajectory.
"""
from __future__ import annotations
import argparse,csv,json
from pathlib import Path
import numpy as np
from check_I2_defect import trapz


def run(path,out):
    with np.load(path,allow_pickle=False) as z:
        U=z['U']; t=z['times']; V=z['stage_U']; st=z['stage_times']
        meta=json.loads(str(z['metadata'].item()))
    if meta['method']!='TR-BDF2' or meta['mu']!=0:
        raise ValueError('This control is defined for the mu=0 TR-BDF2 trajectory only')
    N=int(meta['N']); h=meta['L']/N; p=meta['p']; c0=meta['c0']
    if len(V)!=len(U)-1: raise ValueError('Missing stage states')
    gamma=2-np.sqrt(2)
    if not np.allclose(st,t[:-1]+gamma*np.diff(t),rtol=0,atol=1e-12):
        raise ValueError('Stage times do not match gamma')
    theta=2*np.pi*np.fft.fftfreq(N)
    ah=(66+52*np.cos(theta)+2*np.cos(2*theta))/120
    bh=1j*(10*np.sin(theta)+np.sin(2*theta))/(12*h)
    multiplier=bh/ah
    I20=float(h*np.sum(np.abs(U[0])**(p+1)/(p+1)))
    def delta(u):
        d1=np.fft.ifft(np.fft.fft(u,axis=-1)*multiplier,axis=-1).real
        return c0*h*np.sum(np.abs(u)**(p-1)*u*d1,axis=-1)/I20
    q=np.empty(2*len(U)-1); times=np.empty_like(q)
    q[::2]=delta(U);q[1::2]=delta(V);times[::2]=t;times[1::2]=st
    stage_trap=trapz(q,times)
    rows=[]
    for k in (2,4,8,16,32,64):
        r=np.linspace(0,1,k+1)
        w0=(r-gamma)*(r-1)/gamma
        wg=r*(r-1)/(gamma*(gamma-1))
        w1=r*(r-gamma)/(1-gamma)
        Q=0.
        for j in range(len(V)):
            u=w0[:,None]*U[j]+wg[:,None]*V[j]+w1[:,None]*U[j+1]
            Q+=trapz(delta(u),r*(t[j+1]-t[j]))
        rows.append(dict(subintervals_per_step=k,quadratic_interpolant_trapezoid=Q))
    out=Path(out);out.parent.mkdir(parents=True,exist_ok=True)
    with out.with_suffix('.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)
    report=dict(accepted_nodes=len(U),stage_nodes=len(V),accepted_only_trapezoid=trapz(q[::2],t),
                accepted_plus_stage_trapezoid=stage_trap,quadratic_controls=rows,
                interpretation='Numerical sampling controls on one recorded trajectory; the quadratic interpolant is explicitly defined, not the exact ODE solution. No rigorous error bound is asserted.')
    out.with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))
    return report

if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('trajectory',type=Path); ap.add_argument('--out-prefix',type=Path,required=True)
    a=ap.parse_args();run(a.trajectory,a.out_prefix)
