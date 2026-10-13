#!/usr/bin/env python3
import argparse, json, math, time
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares
import pg_lstable_methods as pm

L=100.0; P=2.0; C0=1.0; MU=0.0; T=10.0
C1=3.0; C2=1.0; X1=40.0; X2=60.0
REF=dict(A_slow=1.353093728, X_slow=53.722409556,
         A_fast=3.998591785, X_fast=62.513939623)
REFUNC=dict(A_slow=0.003877038, X_slow=0.03080959,
            A_fast=0.0006052192, X_fast=0.001659703)
HW=2*np.pi

def compacton(x,p,c,x0,L):
    ell=p*np.pi/(p-1.0)
    alpha=2.0*abs(c)*p/(p+1.0)
    beta=(p-1.0)/(2.0*p)
    gamma=1.0/(p-1.0)
    xi=(x-x0+L/2.0)%L-L/2.0
    u=np.zeros_like(x)
    m=np.abs(xi)<=ell
    u[m]=np.sign(c)*(alpha*np.cos(beta*xi[m])**2)**gamma
    return u

def pdiff(a,b):
    return (a-b+L/2)%L-L/2

def phi(x,X):
    d=pdiff(x,X)
    y=np.zeros_like(x); m=np.abs(d)<=HW
    y[m]=np.cos(d[m]/4.0)**2
    return y

def core_model(x,p):
    A1,X1,A2,X2=p
    return A1*phi(x,X1)+A2*phi(x,X2)

def fit_cores(x,u):
    # Late-time cores are near the independent Direct-WENO locations.
    p0=np.array([REF["A_slow"],REF["X_slow"],REF["A_fast"],REF["X_fast"]])
    lo=np.array([0.7,50.0,3.0,59.0])
    hi=np.array([2.0,57.0,4.8,66.0])
    scale=max(0.02,0.01*np.max(np.abs(u)))
    r=least_squares(lambda p: core_model(x,p)-u,p0,bounds=(lo,hi),
                    loss="soft_l1",f_scale=scale,max_nfev=3000,
                    xtol=1e-12,ftol=1e-12,gtol=1e-12)
    return r.x,core_model(x,r.x)

def rhohw(U):
    z=np.fft.fft(U); th=2*np.pi*np.fft.fftfreq(U.size)
    den=np.sum(np.abs(z)**2)
    return float(np.sum(np.abs(z[np.abs(th)>=np.pi/2])**2)/den)

def i2(U,h):
    return float(h*np.sum(np.abs(U)**(P+1)/(P+1)))

def task_table():
    tasks=[]
    for N in (800,1600):
        for method in ("M2","RadauIIA2","TR-BDF2"):
            for dt in (0.4,0.2,0.1,0.05):
                tasks.append((N,method,dt))
    for method in ("M2","RadauIIA2","TR-BDF2"):
        tasks.append((400,method,0.05))
    return tasks

def run_task(idx,outdir):
    N,method,dt=task_table()[idx]
    h=L/N; x=np.arange(N)*h
    U0=compacton(x,P,C1,X1,L)+compacton(x,P,C2,X2,L)
    U=U0.copy(); Uprev=None
    model=pm.AbsKPetrovGalerkin(N,L,P,C0,MU)
    nsteps=int(round(T/dt))
    if abs(nsteps*dt-T)>1e-12:
        raise RuntimeError("dt must divide T")
    I20=i2(U0,h)
    its=[]; maxres=0.0; status="ok"; fail_t=None
    tic=time.perf_counter()
    for n in range(nsteps):
        if method=="M2":
            X,it,ok,nrm=pm.step_M2(model,U,dt,Uprev)
        elif method=="RadauIIA2":
            X,it,ok,nrm=pm.step_RadauIIA2(model,U,dt,Uprev)
        elif method=="TR-BDF2":
            X,it,ok,nrm=pm.step_TRBDF2(model,U,dt,Uprev)
        else:
            raise ValueError(method)
        its.append(int(it)); maxres=max(maxres,float(nrm))
        if not ok or not np.all(np.isfinite(X)) or np.max(np.abs(X))>50:
            status="newton_failure" if not ok else "breakdown"
            fail_t=(n+1)*dt
            U=X
            break
        Uprev,U=U,X
    wall=time.perf_counter()-tic

    row=dict(task_id=idx,N=N,dx=h,method=method,dt=dt,T=T,p=P,c0=C0,mu=MU,
             c1=C1,c2=C2,x1=X1,x2=X2,status=status,fail_t=fail_t,
             steps_done=len(its),wall_s=wall,mean_newton=float(np.mean(its)) if its else np.nan,
             max_newton=int(np.max(its)) if its else 0,max_newton_residual=maxres)
    if status=="ok":
        pars,core=fit_cores(x,U)
        refpars=np.array([REF["A_slow"],REF["X_slow"],REF["A_fast"],REF["X_fast"]])
        refcore=core_model(x,refpars)
        rem=U-core
        row.update(
            A_slow=float(pars[0]),X_slow=float(pars[1]),
            A_fast=float(pars[2]),X_fast=float(pars[3]),
            err_A_slow=abs(float(pars[0]-REF["A_slow"])),
            err_X_slow=abs(float(pdiff(pars[1],REF["X_slow"]))),
            err_A_fast=abs(float(pars[2]-REF["A_fast"])),
            err_X_fast=abs(float(pdiff(pars[3],REF["X_fast"]))),
            core_model_relL2=float(np.linalg.norm(core-refcore)/np.linalg.norm(refcore)),
            remainder_frac=float(np.linalg.norm(rem)/np.linalg.norm(U)),
            remainder_neg_frac=float(np.linalg.norm(np.minimum(rem,0))/np.linalg.norm(U)),
            remainder_pos_frac=float(np.linalg.norm(np.maximum(rem,0))/np.linalg.norm(U)),
            I2_drift_rel=float((i2(U,h)-I20)/I20),
            rho_hw=rhohw(U),max_abs_u=float(np.max(np.abs(U))),
            ref_unc_A_slow=REFUNC["A_slow"],ref_unc_X_slow=REFUNC["X_slow"],
            ref_unc_A_fast=REFUNC["A_fast"],ref_unc_X_fast=REFUNC["X_fast"],
        )
    outdir=Path(outdir); outdir.mkdir(parents=True,exist_ok=True)
    stem=f"pgcore_{method}_N{N}_dt{str(dt).replace('.','p')}"
    (outdir/(stem+".json")).write_text(json.dumps(row,indent=2,allow_nan=True))
    np.savez_compressed(outdir/(stem+".npz"),x=x,u=U,u0=U0,metadata=json.dumps(row,allow_nan=True))
    print(json.dumps(row,indent=2,allow_nan=True))
    return row

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--task",type=int)
    ap.add_argument("--outdir",default="pg_lstable_core_results")
    ap.add_argument("--list",action="store_true")
    ap.add_argument("--scalar-validation",action="store_true")
    args=ap.parse_args()
    if args.list:
        for i,t in enumerate(task_table()): print(i,t)
        return
    if args.scalar_validation:
        for r in pm.scalar_validation(): print(r)
        return
    if args.task is None: raise SystemExit("--task required")
    run_task(args.task,args.outdir)

if __name__=="__main__":
    main()
