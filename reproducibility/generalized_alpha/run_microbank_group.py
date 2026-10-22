#!/usr/bin/env python3
import csv, os, sys
from pathlib import Path
from ga_core import Case, run_case

RHO=[1.0,0.9,0.75,0.5,0.447626940936,0.25,0.20,0.10,0.05,0.0]
C1=[1.25,2.0,3.0]
MU=[0.0,1e-7,1e-3,1e-1]
DT=[0.05,0.1,0.2]
COMBOS=[(c1,mu,dt) for c1 in C1 for mu in MU for dt in DT]
FIELDS=None

def main():
    g=int(sys.argv[1])
    c1,mu,dt=COMBOS[g]
    outdir=Path(os.environ.get('GA_OUT','results/microbank_groups'))
    outdir.mkdir(parents=True,exist_ok=True)
    rows=[]
    for rho in RHO:
        d=run_case(Case(N=100,c1=c1,mu=mu,dt=dt,rho_inf=rho,newton_tol=1e-10,newton_maxit=20))
        d.update(group_id=g,label=f'GA_{rho:g}')
        rows.append(d)
    # Re-run the Padé control to make the new bank self-contained.
    d=run_case(Case(N=100,c1=c1,mu=mu,dt=dt,rho_inf=None,method='M2',newton_tol=1e-10,newton_maxit=40))
    d.update(group_id=g,label='M2')
    rows.append(d)
    fn=outdir/f'group_{g:02d}.csv'
    with fn.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0].keys()));w.writeheader();w.writerows(rows)
    print(f'MICROBANK_GROUP_OK group={g} rows={len(rows)} file={fn}',flush=True)
if __name__=='__main__': main()
