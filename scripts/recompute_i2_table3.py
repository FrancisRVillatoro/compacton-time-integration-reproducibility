#!/usr/bin/env python3
from __future__ import annotations
import csv, sys, tempfile
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'code'))
from check_I2_defect import audit

RAW=ROOT/'raw'/'pg_table3'
OUT=ROOT/'data'/'table3_I2_corrections_recomputed.csv'
CASES=[
 ('M2',.4,'M2_N1600_dt0p4'),
 ('M2',.2,'M2_N1600_dt0p2'),
 ('M2',.05,'M2_N1600_dt0p05'),
 ('RadauIIA2',.4,'RadauIIA2_N1600_dt0p4'),
 ('RadauIIA2',.2,'RadauIIA2_N1600_dt0p2'),
 ('RadauIIA2',.05,'RadauIIA2_N1600_dt0p05'),
 ('TR-BDF2',.1,'TR-BDF2_N1600_dt0p1'),
 ('TR-BDF2',.05,'TR-BDF2_N1600_dt0p05'),
 ('TR-BDF2',.025,'TR-BDF2_N1600_dt0p025')]

def main():
    rows=[]
    with tempfile.TemporaryDirectory() as td:
        td=Path(td)
        for method,dt,case in CASES:
            rep=audit(RAW/case/'trajectory.npz', ROOT/'code', td/case)
            rows.append(dict(
                method=method,dt=dt,N=1600,
                raw_I2_drift_signed_relative=rep['raw_drift_signed_relative'],
                semidiscrete_defect_trapezoid_signed_relative=rep['semidiscrete_defect_trapezoid_signed_relative'],
                corrected_I2_drift_signed_relative=rep['corrected_drift_trapezoid_signed_relative'],
                corrected_I2_drift_absolute_relative=rep['corrected_drift_trapezoid_absolute_relative']))
    with OUT.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=rows[0].keys()); w.writeheader(); w.writerows(rows)
    frozen=list(csv.DictReader((ROOT/'data'/'table3_I2_corrections.csv').open()))
    F={(r['method'],float(r['dt'])):r for r in frozen}
    md=0.0
    for r in rows:
        fr=F[(r['method'],float(r['dt']))]
        md=max(md,abs(r['corrected_I2_drift_absolute_relative']-float(fr['corrected_I2_drift_absolute_relative'])))
    print(f'max corrected-I2 difference from frozen table: {md:.3e}')
    if md>5e-12:
        raise RuntimeError('I2 correction mismatch')

if __name__=='__main__':
    main()
