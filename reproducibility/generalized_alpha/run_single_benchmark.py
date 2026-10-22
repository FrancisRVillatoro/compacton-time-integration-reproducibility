#!/usr/bin/env python3
import csv, os
from pathlib import Path
from ga_core import Case, run_case

# Repeats the historical N=1000, c1=2, c2=-1 single-collision protocol:
# x1=30, x2=65, mu=1e-7, dt=.1, T=60, diagnostic window t>=45.
RHO=[0.75,0.5,0.447626940936,0.25,0.0]
outdir=Path(os.environ.get('GA_OUT','results'));outdir.mkdir(parents=True,exist_ok=True)
rows=[]
for method,rho in [('M1',None)]+[('GA',r) for r in RHO]+[('M2',None)]:
    d=run_case(Case(
        N=1000,c1=2.0,c2=-1.0,x1=30.0,x2=65.0,mu=1e-7,dt=0.1,
        rho_inf=rho,method=method,T=60.0,tpost=45.0,newton_tol=1e-10,newton_maxit=(40 if method=='M2' else 20),
    ))
    d['label']=method if method!='GA' else f'GA_{rho:g}'
    rows.append(d)
fn=outdir/'single_benchmark_corrected.csv'
with fn.open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0].keys()));w.writeheader();w.writerows(rows)
print(f'SINGLE_BENCHMARK_OK rows={len(rows)} file={fn}',flush=True)
