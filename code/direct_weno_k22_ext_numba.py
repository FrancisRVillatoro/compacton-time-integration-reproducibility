#!/usr/bin/env python3
"""Time-localization audit for classical K(2,2) Direct-WENO5.

Cases:
  fast:   one c=3 compacton at x0=40 (speed c-c0=2)
  slow:   one c=1 compacton at x0=60 (stationary for c0=1)
  cc:     c1=3 at x1=40 plus c2=1 at x2=60

Saves synchronized diagnostic histories and snapshots at nout equally spaced times.
"""
import argparse, json, time
from pathlib import Path
import numpy as np
from numba import njit

@njit(cache=True)
def _mapw(w,d):
    return w*(d+d*d-3*d*w+w*w)/(d*d+w*(1-2*d))

@njit(cache=True)
def _weights(b0,b1,b2,d0,d1,d2,eps):
    a0=d0/(eps+b0)**2; a1=d1/(eps+b1)**2; a2=d2/(eps+b2)**2
    s=a0+a1+a2
    w0=a0/s; w1=a1/s; w2=a2/s
    g0=_mapw(w0,d0); g1=_mapw(w1,d1); g2=_mapw(w2,d2)
    sg=g0+g1+g2
    return g0/sg,g1/sg,g2/sg

@njit(cache=True)
def _conv_plus(a,j,eps):
    n=a.size
    am2=a[(j-2)%n]; am1=a[(j-1)%n]; a0=a[j]; ap1=a[(j+1)%n]; ap2=a[(j+2)%n]
    q0=(1/3)*am2-(7/6)*am1+(11/6)*a0
    q1=-(1/6)*am1+(5/6)*a0+(1/3)*ap1
    q2=(1/3)*a0+(5/6)*ap1-(1/6)*ap2
    b0=(13/12)*(am2-2*am1+a0)**2+0.25*(am2-4*am1+3*a0)**2
    b1=(13/12)*(am1-2*a0+ap1)**2+0.25*(am1-ap1)**2
    b2=(13/12)*(a0-2*ap1+ap2)**2+0.25*(3*a0-4*ap1+ap2)**2
    w0,w1,w2=_weights(b0,b1,b2,0.1,0.6,0.3,eps)
    return w0*q0+w1*q1+w2*q2

@njit(cache=True)
def _conv_minus(a,j,eps):
    n=a.size
    am1=a[(j-1)%n]; a0=a[j]; ap1=a[(j+1)%n]; ap2=a[(j+2)%n]; ap3=a[(j+3)%n]
    bm2=ap3; bm1=ap2; b0v=ap1; bp1=a0; bp2=am1
    q0=(1/3)*bm2-(7/6)*bm1+(11/6)*b0v
    q1=-(1/6)*bm1+(5/6)*b0v+(1/3)*bp1
    q2=(1/3)*b0v+(5/6)*bp1-(1/6)*bp2
    b0=(13/12)*(bm2-2*bm1+b0v)**2+0.25*(bm2-4*bm1+3*b0v)**2
    b1=(13/12)*(bm1-2*b0v+bp1)**2+0.25*(bm1-bp1)**2
    b2=(13/12)*(b0v-2*bp1+bp2)**2+0.25*(3*b0v-4*bp1+bp2)**2
    w0,w1,w2=_weights(b0,b1,b2,0.1,0.6,0.3,eps)
    return w0*q0+w1*q1+w2*q2

@njit(cache=True)
def _disp_plus(a,j,eps):
    n=a.size
    am2=a[(j-2)%n]; am1=a[(j-1)%n]; a0=a[j]
    ap1=a[(j+1)%n]; ap2=a[(j+2)%n]; ap3=a[(j+3)%n]; ap4=a[(j+4)%n]
    c0=-0.25*am2+1.5*am1-2*a0+0.5*ap1+0.25*ap2
    c1=0.25*am1+0.5*a0-2*ap1+1.5*ap2-0.25*ap3
    c2=1.75*a0-4.5*ap1+4*ap2-1.5*ap3+0.25*ap4
    B0=-0.5*am2+am1-ap1+0.5*ap2
    C0=0.5*am2-2*am1+3*a0-2*ap1+0.5*ap2
    B1=-1.5*am1+5*a0-6*ap1+3*ap2-0.5*ap3
    C1=0.5*am1-2*a0+3*ap1-2*ap2+0.5*ap3
    B2=-2.5*a0+9*ap1-12*ap2+7*ap3-1.5*ap4
    C2=0.5*a0-2*ap1+3*ap2-2*ap3+0.5*ap4
    b0=B0*B0+(13/3)*C0*C0; b1=B1*B1+(13/3)*C1*C1; b2=B2*B2+(13/3)*C2*C2
    w0,w1,w2=_weights(b0,b1,b2,4/15,0.5,7/30,eps)
    return w0*c0+w1*c1+w2*c2

@njit(cache=True)
def _disp_minus(a,j,eps):
    n=a.size
    am3=a[(j-3)%n]; am2=a[(j-2)%n]; am1=a[(j-1)%n]; a0=a[j]
    ap1=a[(j+1)%n]; ap2=a[(j+2)%n]; ap3=a[(j+3)%n]
    bm2=ap3; bm1=ap2; b0v=ap1; bp1=a0; bp2=am1; bp3=am2; bp4=am3
    c0=-0.25*bm2+1.5*bm1-2*b0v+0.5*bp1+0.25*bp2
    c1=0.25*bm1+0.5*b0v-2*bp1+1.5*bp2-0.25*bp3
    c2=1.75*b0v-4.5*bp1+4*bp2-1.5*bp3+0.25*bp4
    B0=-0.5*bm2+bm1-bp1+0.5*bp2
    C0=0.5*bm2-2*bm1+3*b0v-2*bp1+0.5*bp2
    B1=-1.5*bm1+5*b0v-6*bp1+3*bp2-0.5*bp3
    C1=0.5*bm1-2*b0v+3*bp1-2*bp2+0.5*bp3
    B2=-2.5*b0v+9*bp1-12*bp2+7*bp3-1.5*bp4
    C2=0.5*b0v-2*bp1+3*bp2-2*bp3+0.5*bp4
    b0=B0*B0+(13/3)*C0*C0; b1=B1*B1+(13/3)*C1*C1; b2=B2*B2+(13/3)*C2*C2
    w0,w1,w2=_weights(b0,b1,b2,4/15,0.5,7/30,eps)
    return w0*c0+w1*c1+w2*c2

@njit(cache=True)
def rhs(u,dx,c0,mu,eps,ext):
    n=u.size; alpha=0.0; beta=0.0
    for j in range(n):
        dq=2.0*abs(u[j]) if ext==1 else 2.0*u[j]
        a=abs(dq-c0); b=abs(dq)
        if a>alpha: alpha=a
        if b>beta: beta=b
    alpha=max(alpha,1e-14); beta=max(beta,1e-14)
    fp=np.empty(n); fm=np.empty(n); gp=np.empty(n); gm=np.empty(n)
    for j in range(n):
        q=abs(u[j])*u[j] if ext==1 else u[j]*u[j]; F=q-c0*u[j]
        fp[j]=0.5*(F+alpha*u[j]); fm[j]=0.5*(F-alpha*u[j])
        gp[j]=0.5*(q+beta*u[j]); gm[j]=0.5*(q-beta*u[j])
    fh=np.empty(n); gh=np.empty(n)
    for j in range(n):
        fh[j]=_conv_plus(fp,j,eps)+_conv_minus(fm,j,eps)
        gh[j]=_disp_plus(gp,j,eps)+_disp_minus(gm,j,eps)
    out=np.empty(n); invdx=1.0/dx; invdx3=invdx**3; visc=mu/(dx*dx)
    for j in range(n):
        jm1=(j-1)%n
        r=-(fh[j]-fh[jm1])*invdx-(gh[j]-gh[jm1])*invdx3
        if mu!=0.0:
            d4=u[(j-2)%n]-4*u[(j-1)%n]+6*u[j]-4*u[(j+1)%n]+u[(j+2)%n]
            r -= visc*d4
        out[j]=r
    return out

@njit(cache=True)
def ssprk3(u,dx,dt,c0,mu,eps,ext):
    k1=rhs(u,dx,c0,mu,eps,ext); u1=u+dt*k1
    k2=rhs(u1,dx,c0,mu,eps,ext); u2=.75*u+.25*(u1+dt*k2)
    k3=rhs(u2,dx,c0,mu,eps,ext)
    return u/3+(2/3)*(u2+dt*k3)

@njit(cache=True)
def choose_dt(u,dx,c0,mu,cfl,remaining,ext):
    alpha=0.0; beta=0.0
    for j in range(u.size):
        dq=2.0*abs(u[j]) if ext==1 else 2.0*u[j]
        a=abs(dq-c0); b=abs(dq)
        if a>alpha: alpha=a
        if b>beta: beta=b
    alpha=max(alpha,1e-14); beta=max(beta,1e-14)
    dt=min(cfl*dx**(5/3)/alpha,cfl*dx**3/beta,remaining)
    if mu>0: dt=min(dt,.1*dx*dx/mu)
    return dt

def compacton(x,c,x0,L):
    xi=(x-x0+L/2)%L-L/2
    u=np.zeros_like(x); m=np.abs(xi)<=2*np.pi
    u[m]=np.sign(c)*(4.0*abs(c)/3.0)*np.cos(xi[m]/4.0)**2
    return u

def diagnostics(u,u0,dx):
    I10=dx*np.sum(u0); I1=dx*np.sum(u)
    A30=dx*np.sum(np.abs(u0)**3/3.0); A3=dx*np.sum(np.abs(u)**3/3.0)
    uh=np.fft.fft(u); th=2*np.pi*np.fft.fftfreq(u.size)
    rho=np.sum(np.abs(uh[np.abs(th)>=np.pi/2])**2)/max(np.sum(np.abs(uh)**2),np.finfo(float).eps)
    neg=np.minimum(u,0.0)
    return (I1-I10,(A3-A30)/A30,float(rho),float(np.max(u)),float(np.min(u)),
            float(np.linalg.norm(neg)/max(np.linalg.norm(u),np.finfo(float).eps)))

def exact_single(x,t,c,x0,c0,L):
    return compacton(x,c,x0+(c-c0)*t,L)

def rel(u,v):
    return float(np.linalg.norm(u-v)/max(np.linalg.norm(v),np.finfo(float).eps))

def run(args):
    N=args.N; L=args.L; dx=L/N; x=np.arange(N)*dx
    if args.case=='fast':
        u=compacton(x,3.0,40.0,L); exact_c=3.0; exact_x0=40.0
    elif args.case=='slow':
        u=compacton(x,1.0,60.0,L); exact_c=1.0; exact_x0=60.0
    elif args.case=='cc':
        u=compacton(x,3.0,40.0,L)+compacton(x,1.0,60.0,L); exact_c=None; exact_x0=None
    else:
        raise ValueError(args.case)
    u0=u.copy()
    _=ssprk3(u.copy(),dx,min(1e-10,args.T),args.c0,args.mu,args.eps,args.ext)
    times=np.linspace(0.0,args.T,args.nout)
    snaps=np.empty((args.nout,N),dtype=np.float64)
    cols=np.empty((args.nout,8),dtype=np.float64)
    t=0.0; steps=0; dtmin=float('inf'); dtmax=0.0; tic=time.perf_counter()
    for k,target in enumerate(times):
        while t < target-1e-14:
            dt=choose_dt(u,dx,args.c0,args.mu,args.cfl,target-t,args.ext)
            u=ssprk3(u,dx,dt,args.c0,args.mu,args.eps,args.ext)
            t += dt; steps += 1
            dtmin=min(dtmin,dt); dtmax=max(dtmax,dt)
            if not np.all(np.isfinite(u)) or np.max(np.abs(u))>50:
                raise RuntimeError(f'breakdown at t={t}')
        snaps[k]=u
        I1d,I2d,rho,umax,umin,negL2=diagnostics(u,u0,dx)
        if exact_c is not None:
            ue=exact_single(x,target,exact_c,exact_x0,args.c0,L)
            eexact=rel(u,ue)
        else:
            eexact=np.nan
        cols[k]=[I1d,I2d,rho,umax,umin,negL2,eexact,float(np.argmax(u))*dx]
    wall=time.perf_counter()-tic
    meta=dict(model=('absolute_K22' if args.ext==1 else 'classical_K22')+'_localization',ext=args.ext,case=args.case,N=N,dx=dx,L=L,c0=args.c0,mu=args.mu,
              cfl=args.cfl,T=args.T,eps=args.eps,nout=args.nout,steps=steps,dt_min=dtmin,dt_max=dtmax,
              wall_s=wall,status='ok')
    np.savez_compressed(args.out+'.npz',x=x,times=times,snapshots=snaps,u0=u0,
                        diag=cols,diag_names=np.array(['I1_drift','I2_drift_rel','rho_hw','max_u','min_u','negL2','single_exact_relL2','argmax_x']),
                        metadata=json.dumps(meta))
    Path(args.out+'.json').write_text(json.dumps(meta,indent=2))
    print(json.dumps(meta,indent=2))
    print('K22_LOCALIZATION_TASK_OK')

def selftest():
    u=np.ones(80)*0.3
    defect=float(np.max(np.abs(rhs(u,.2,0.0,0.0,1e-6,0))))
    if defect>1e-12: raise RuntimeError(f'constant-state RHS defect {defect}')
    print(f'DIRECT_WENO_K22_LOCALIZATION_SELFTEST_OK defect={defect:.3e}')

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--selftest',action='store_true'); ap.add_argument('--run',action='store_true')
    ap.add_argument('--case',choices=['fast','slow','cc'],default='cc'); ap.add_argument('--N',type=int,default=800)
    ap.add_argument('--L',type=float,default=100.); ap.add_argument('--c0',type=float,default=1.); ap.add_argument('--mu',type=float,default=0.)
    ap.add_argument('--cfl',type=float,default=.075); ap.add_argument('--T',type=float,default=10.); ap.add_argument('--eps',type=float,default=1e-6)
    ap.add_argument('--nout',type=int,default=201); ap.add_argument('--ext',type=int,choices=[0,1],default=0); ap.add_argument('--out',default='k22_localization')
    a=ap.parse_args()
    if a.selftest: selftest()
    if a.run: run(a)
    if not a.selftest and not a.run: ap.error('choose --selftest and/or --run')
if __name__=='__main__': main()
