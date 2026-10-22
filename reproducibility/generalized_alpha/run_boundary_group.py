#!/usr/bin/env python3
import csv, os, sys
from pathlib import Path
from ga_core import Case, run_case

RHO=[1.0,0.9,0.75,0.5,0.25,0.1,0.0]
NVAL=[500,1000]
MU=[0.0,1e-7,1e-3,1e-1]
COMBOS=[(N,mu) for N in NVAL for mu in MU]

def main():
    g=int(sys.argv[1]);N,mu=COMBOS[g]
    outdir=Path(os.environ.get('GA_OUT','results/boundary_groups'));outdir.mkdir(parents=True,exist_ok=True)
    rows=[]
    for rho in RHO:
        d=run_case(Case(N=N,c1=3.0,mu=mu,dt=0.2,rho_inf=rho,newton_tol=1e-10,newton_maxit=20))
        d.update(group_id=g,label=f'GA_{rho:g}')
        rows.append(d)
    d=run_case(Case(N=N,c1=3.0,mu=mu,dt=0.2,rho_inf=None,method='M2',newton_tol=1e-10,newton_maxit=40))
    d.update(group_id=g,label='M2');rows.append(d)
    fn=outdir/f'group_{g:02d}.csv'
    with fn.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0].keys()));w.writeheader();w.writerows(rows)
    print(f'BOUNDARY_GROUP_OK group={g} N={N} mu={mu:g} rows={len(rows)} file={fn}',flush=True)
if __name__=='__main__':main()
