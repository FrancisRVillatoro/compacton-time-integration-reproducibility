#!/usr/bin/env python3
"""Replay the archived production integrators, recording accepted states.

Only the recording wrapper is new; production functions are left untouched.
No production result is overwritten. --mode table3 computes all nine listed
N=1600 cases. Refinement mode is diagnostic and is NOT an original HPC run.
"""
from __future__ import annotations
import argparse, contextlib, hashlib, importlib, json, os, platform, sys, time
from pathlib import Path
import numpy as np
import scipy
from check_I2_defect import audit, diagnostic_series, trapz


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source-dir',type=Path,default=Path(__file__).parent)
    ap.add_argument('--outdir',type=Path,default=Path('I2_audit_results'))
    ap.add_argument('--mode',choices=['fine','table3','refinement'],default='fine')
    ap.add_argument('--dt',type=float,nargs='+',default=[0.0125,0.00625,0.003125])
    ap.add_argument('--keep-stages',action='store_true')
    args=ap.parse_args()
    source=args.source_dir.resolve(); args.outdir.mkdir(parents=True,exist_ok=True)
    sys.path.insert(0,str(source))
    pm=importlib.import_module('pg_lstable_methods')
    base=importlib.import_module('run_pg_lstable_core_task')
    closure=importlib.import_module('run_pg_lstable_closure')
    modules=[pm,base,closure]
    for module in modules:
        if Path(module.__file__).resolve().parent!=source:
            raise RuntimeError('Unexpected source path: '+module.__file__)
    if (base.L,base.P,base.C0,base.MU,base.T,base.C1,base.C2,base.X1,base.X2)!=(100.,2.,1.,0.,10.,3.,1.,40.,60.):
        raise RuntimeError('Source benchmark configuration differs from audited case; refusing to substitute it silently')
    hashes={Path(m.__file__).name:hashlib.sha256(Path(m.__file__).read_bytes()).hexdigest() for m in modules}
    if args.mode=='fine': cases=[('TR-BDF2',0.025,'closure',1)]
    elif args.mode=='table3':
        cases=[('M2',.4,'core',12),('M2',.2,'core',13),('M2',.05,'core',15),
               ('RadauIIA2',.4,'closure',2),('RadauIIA2',.2,'core',17),('RadauIIA2',.05,'core',19),
               ('TR-BDF2',.1,'core',22),('TR-BDF2',.05,'core',23),('TR-BDF2',.025,'closure',1)]
    else: cases=[('TR-BDF2',d,'refinement',0) for d in args.dt]
    summaries=[]
    for method,dt,kind,idx in cases:
        if dt<=0 or abs(round(base.T/dt)*dt-base.T)>1e-12: raise ValueError('Invalid dt')
        name=f"{method}_N1600_dt{str(dt).replace('.','p')}"
        out=args.outdir/name
        out.mkdir(parents=True,exist_ok=True)
        if (out/'trajectory.npz').exists(): raise FileExistsError(str(out)+' already contains a trajectory')
        original_newton=pm.damped_newton
        owner=closure if method=='RadauIIA2' and kind=='closure' else pm
        attr='step_radau_reinforced' if owner is closure else {'M2':'step_M2','RadauIIA2':'step_RadauIIA2','TR-BDF2':'step_TRBDF2'}[method]
        original_step=getattr(owner,attr)
        saved=[]; stage_states=[]; iters=[]
        def recording_step(model,U,dt,Uprev=None):
            if not saved: saved.append(U.copy())
            stage_solutions=[]
            def recording_newton(*a,**kw):
                result=original_newton(*a,**kw)
                if result[2]: stage_solutions.append(result[0].copy())
                return result
            if args.keep_stages and method=='TR-BDF2': pm.damped_newton=recording_newton
            try:
                result=original_step(model,U,dt,Uprev)
            finally:
                pm.damped_newton=original_newton
            X,it,ok,res=result[:4]
            if not ok or not np.all(np.isfinite(X)):
                raise RuntimeError(f'{name}: nonlinear failure at step {len(saved)} residual={res}')
            if args.keep_stages and method=='TR-BDF2':
                if len(stage_solutions)!=2: raise RuntimeError('Unexpected TR-BDF2 stage count')
                stage_states.append(stage_solutions[0])
            saved.append(X.copy()); iters.append([it,res])
            return result
        setattr(owner,attr,recording_step)
        original_table=base.task_table
        try:
            with (out/'production_log.txt').open('w') as log, contextlib.redirect_stdout(log):
                if kind=='closure': closure.run_task(idx,out)
                else:
                    if kind=='refinement': base.task_table=lambda:[(1600,'TR-BDF2',dt)]
                    base.run_task(idx,out)
        finally:
            setattr(owner,attr,original_step); pm.damped_newton=original_newton; base.task_table=original_table
        U=np.asarray(saved); times=np.arange(len(U))*dt
        meta=dict(L=base.L,p=base.P,c0=base.C0,mu=base.MU,T=base.T,N=1600,dt=dt,method=method,
                  x1=base.X1,x2=base.X2,c1=base.C1,c2=base.C2,
                  provenance='replay of supplied source in the executing environment; not original Picasso trajectory',
                  mode=kind,source_sha256=hashes,python=sys.version,numpy=np.__version__,scipy=scipy.__version__,
                  hostname=platform.node(),thread_environment={k:os.environ.get(k) for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS')})
        extra={}
        if stage_states:
            extra=dict(stage_times=np.arange(len(stage_states))*dt+(2-np.sqrt(2))*dt,
                       stage_U=np.asarray(stage_states))
        path=out/'trajectory.npz'
        np.savez_compressed(path,U=U,times=times,x=np.arange(1600)*base.L/1600,metadata=json.dumps(meta),iterations=np.asarray(iters),**extra)
        rep=audit(path,source,out/'I2_defect')
        rep['case']=name
        summaries.append(rep)
        print(name,'raw=',rep['raw_drift_signed_relative'],'Qtrap=',rep['semidiscrete_defect_trapezoid_signed_relative'],
              'corrected=',rep['corrected_drift_trapezoid_signed_relative'],flush=True)
    (args.outdir/'summary.json').write_text(json.dumps(summaries,indent=2,allow_nan=False)+'\n')
    print('PG_I2_REPLAY_COMPLETE',flush=True)

if __name__=='__main__': main()
