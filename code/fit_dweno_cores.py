#!/usr/bin/env python3
"""Fit the two-core model to Direct-WENO final states with the production fit
(run_pg_lstable_core_task.fit_cores) and report the core parameters, the core
discrepancies between successive grids, and the distance to the converged
Petrov--Galerkin cores.

Usage:  python fit_dweno_cores.py --pg-code DIR  file1.npz file2.npz ...
DIR must contain run_pg_lstable_core_task.py and pg_lstable_methods.py.
"""
import argparse, json, sys
from pathlib import Path
import numpy as np

# converged Petrov--Galerkin cores (M2, dt = 0.05, N = 1600; production fit)
PG = np.array([1.348983439, 53.757385902, 3.997856645, 62.516210715])

def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--pg-code', required=True)
    ap.add_argument('files', nargs='+'); a = ap.parse_args()
    sys.path.insert(0, str(Path(a.pg_code).resolve()))
    import run_pg_lstable_core_task as rt
    rows = []
    for f in a.files:
        z = np.load(f, allow_pickle=False)
        meta = json.loads(str(z['metadata'].item()))
        x = z['x']; u = z['snapshots'][-1]
        pars, _ = rt.fit_cores(x, u)
        xf = np.arange(1600) * 100.0 / 1600          # common evaluation grid
        e_pg = np.linalg.norm(rt.core_model(xf, pars) - rt.core_model(xf, PG)) / np.linalg.norm(rt.core_model(xf, PG))
        rows.append(dict(file=Path(f).name, N=meta['N'], ext=meta.get('ext', 0), A_s=pars[0], X_s=pars[1],
                         A_f=pars[2], X_f=pars[3], E_core_vs_PG=e_pg, min_u=float(u.min())))
    rows.sort(key=lambda r: (r['ext'], r['N']))
    for r in rows:
        print(f"ext={r['ext']} N={r['N']:5d}  A_s={r['A_s']:.6f} X_s={r['X_s']:.6f} A_f={r['A_f']:.6f} "
              f"X_f={r['X_f']:.6f}  E_core(vs PG)={r['E_core_vs_PG']:.3e}  min u={r['min_u']:+.3e}")
    json.dump(rows, open('dweno_core_fits.json', 'w'), indent=2)

if __name__ == '__main__':
    main()
