#!/usr/bin/env python3
from __future__ import annotations
import csv, json, sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'code'))
import run_pg_lstable_core_task as rt

FITS = ROOT/'data'/'dweno_flux_final_fits.csv'
OUT = ROOT/'data'/'table3_absolute_reference_recomputed.csv'
RAW = ROOT/'raw'/'pg_table3'
CASES = [
 ('M2',.4,'M2_N1600_dt0p4','pgcore_M2_N1600_dt0p4.npz'),
 ('M2',.2,'M2_N1600_dt0p2','pgcore_M2_N1600_dt0p2.npz'),
 ('M2',.05,'M2_N1600_dt0p05','pgcore_M2_N1600_dt0p05.npz'),
 ('RadauIIA2',.4,'RadauIIA2_N1600_dt0p4','closure_RadauIIA2_N1600_dt0p4.npz'),
 ('RadauIIA2',.2,'RadauIIA2_N1600_dt0p2','pgcore_RadauIIA2_N1600_dt0p2.npz'),
 ('RadauIIA2',.05,'RadauIIA2_N1600_dt0p05','pgcore_RadauIIA2_N1600_dt0p05.npz'),
 ('TR-BDF2',.1,'TR-BDF2_N1600_dt0p1','pgcore_TR-BDF2_N1600_dt0p1.npz'),
 ('TR-BDF2',.05,'TR-BDF2_N1600_dt0p05','pgcore_TR-BDF2_N1600_dt0p05.npz'),
 ('TR-BDF2',.025,'TR-BDF2_N1600_dt0p025','closure_TR-BDF2_N1600_dt0p025.npz')]

def refpars():
    with FITS.open(newline='') as f:
        rows=list(csv.DictReader(f))
    r=next(r for r in rows if int(r['ext'])==1 and int(r['N'])==3200)
    return np.array([float(r['A_s']),float(r['X_s']),float(r['A_f']),float(r['X_f'])])

def err(p,r,n=12800):
    x=np.arange(n)*100.0/n
    m=rt.core_model(x,p); mr=rt.core_model(x,r)
    return float(np.linalg.norm(m-mr)/np.linalg.norm(mr))

def main():
    ref=refpars(); rows=[]
    for method,dt,case,fn in CASES:
        p=RAW/case/fn
        with np.load(p,allow_pickle=False) as a:
            pars,_=rt.fit_cores(a['x'],a['u'])
            meta=json.loads(str(a['metadata'].item()))
        stored=np.array([meta['A_slow'],meta['X_slow'],meta['A_fast'],meta['X_fast']])
        delta=float(np.max(np.abs(pars-stored)))
        if delta>5e-8:
            raise RuntimeError(f'fit regression failed for {case}: {delta}')
        rows.append(dict(method=method,dt=dt,N=1600,A_slow=pars[0],X_slow=pars[1],A_fast=pars[2],X_fast=pars[3],E_core=err(pars,ref),max_refit_vs_stored=delta))
    with OUT.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=rows[0].keys()); w.writeheader(); w.writerows(rows)
    expected=list(csv.DictReader((ROOT/'data'/'table3_final.csv').open()))
    ex={(r['method'],float(r['dt'])):float(r['E_core']) for r in expected}
    maxdiff=max(abs(r['E_core']-ex[(r['method'],float(r['dt']))]) for r in rows)
    print(f'max difference from frozen final table: {maxdiff:.3e}')
    if maxdiff>2e-9:
        raise RuntimeError('frozen table mismatch')

if __name__=='__main__':
    main()
