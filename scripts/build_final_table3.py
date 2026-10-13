#!/usr/bin/env python3
from pathlib import Path
import csv
ROOT=Path(__file__).resolve().parents[1]
ref=list(csv.DictReader((ROOT/'data'/'table3_absolute_reference_recomputed.csv').open()))
i2=list(csv.DictReader((ROOT/'data'/'table3_I2_corrections.csv').open()))
cost=list(csv.DictReader((ROOT/'data'/'table3_picasso_costs.csv').open()))
R={(r['method'],float(r['dt'])):r for r in ref}; I={(r['method'],float(r['dt'])):r for r in i2}; C={(r['method'],float(r['dt'])):r for r in cost}
order=[('M2',.4),('M2',.2),('M2',.05),('RadauIIA2',.4),('RadauIIA2',.2),('RadauIIA2',.05),('TR-BDF2',.1),('TR-BDF2',.05),('TR-BDF2',.025)]
rows=[]
for k in order:
 r,q,c=R[k],I[k],C[k]
 rows.append(dict(method=k[0],dt=k[1],N=1600,E_core=float(r['E_core']),raw_I2_drift_signed=float(q['raw_I2_drift_signed_relative']),semidiscrete_defect_signed=float(q['semidiscrete_defect_trapezoid_signed_relative']),corrected_I2_drift_signed=float(q['corrected_I2_drift_signed_relative']),corrected_I2_drift_abs=float(q['corrected_I2_drift_absolute_relative']),wall_s=float(c['wall_s']),mean_newton=float(c['mean_newton']),A_slow=float(r['A_slow']),X_slow=float(r['X_slow']),A_fast=float(r['A_fast']),X_fast=float(r['X_fast'])))
out=ROOT/'data'/'table3_final_recomputed.csv'
with out.open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)
orig=list(csv.DictReader((ROOT/'data'/'table3_final.csv').open()))
maxd=max(abs(float(a['E_core'])-float(b['E_core'])) for a,b in zip(rows,orig))
print(f'max E_core difference from frozen final table: {maxd:.3e}')
