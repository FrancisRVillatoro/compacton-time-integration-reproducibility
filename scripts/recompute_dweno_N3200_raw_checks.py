#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'code'))
import run_pg_lstable_core_task as rt

FITS = ROOT / 'data' / 'dweno_flux_final_fits.csv'
PAIRS = ROOT / 'data' / 'dweno_flux_final_pairs.csv'
OUT = ROOT / 'data' / 'dweno_N3200_raw_refit.csv'

RAW = {
    0: ROOT / 'raw' / 'direct_weno' / 'k22ext_cc_N3200_ext0.npz',
    1: ROOT / 'raw' / 'direct_weno' / 'k22ext_cc_N3200_ext1.npz',
}
EXPECTED_SHA = {
    0: '050b4a9c64f9c10cd90df60765cfcbe457af3d40b59b7f85befe0ce2db302fc8',
    1: 'a70fbfa816b97f07535fc203ca56a43a5e2b184618e51ee3e4f0ab08b0550560',
}


def sha256(path: Path) -> str:
    import hashlib
    h = hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()


def fit_raw(path: Path, ext: int):
    z = np.load(path, allow_pickle=False)
    meta = json.loads(str(z['metadata'].item()))
    assert int(meta['ext']) == ext
    assert int(meta['N']) == 3200
    assert meta['model'] == ('classical_K22_localization' if ext == 0 else 'absolute_K22_localization')
    assert meta['case'] == 'cc'
    for k, v in [('L', 100.0), ('c0', 1.0), ('mu', 0.0), ('T', 10.0), ('eps', 1e-6), ('cfl', 0.075)]:
        if abs(float(meta[k]) - v) > 1e-14:
            raise RuntimeError(f'metadata mismatch {k}: {meta[k]} != {v}')
    pars, _ = rt.fit_cores(z['x'], z['snapshots'][-1])
    return meta, np.asarray(pars, dtype=float)


def core_error(p, q, n=12800):
    x = np.arange(n) * 100.0 / n
    mp = rt.core_model(x, p)
    mq = rt.core_model(x, q)
    return float(np.linalg.norm(mp - mq) / np.linalg.norm(mq))


def main():
    frozen_fits = list(csv.DictReader(FITS.open()))
    frozen_pairs = list(csv.DictReader(PAIRS.open()))
    rows = []
    pars_by_ext = {}

    for ext, path in RAW.items():
        got = sha256(path)
        if got != EXPECTED_SHA[ext]:
            raise RuntimeError(f'raw hash mismatch ext={ext}: {got}')
        meta, pars = fit_raw(path, ext)
        pars_by_ext[ext] = pars
        frozen = next(r for r in frozen_fits if int(r['ext']) == ext and int(r['N']) == 3200)
        frozen_pars = np.array([float(frozen['A_s']), float(frozen['X_s']), float(frozen['A_f']), float(frozen['X_f'])])
        max_delta = float(np.max(np.abs(pars - frozen_pars)))
        if max_delta > 5e-8:
            raise RuntimeError(f'N=3200 raw refit mismatch ext={ext}: {max_delta}')
        rows.append(dict(
            ext=ext, N=3200, model=meta['model'], sha256=got,
            A_s=pars[0], X_s=pars[1], A_f=pars[2], X_f=pars[3],
            max_refit_vs_frozen=max_delta,
        ))

    e = core_error(pars_by_ext[0], pars_by_ext[1], 12800)
    fp = next(r for r in frozen_pairs if r['kind'] == 'flux_swap_same_grid' and int(r['N_1']) == 3200 and int(r['N_2']) == 3200)
    ef = float(fp['E_core_eval12800'])
    diff = abs(e - ef)
    if diff > 5e-10:
        raise RuntimeError(f'N=3200 flux-swap mismatch: recomputed={e}, frozen={ef}, diff={diff}')

    with OUT.open('w', newline='') as f:
        fieldnames = list(rows[0].keys())
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader(); w.writerows(rows)

    report = ROOT / 'audit' / 'dweno_N3200_raw_check.txt'
    report.write_text(
        'DIRECT-WENO N=3200 RAW CHECK\n'
        + '\n'.join(
            f"ext={r['ext']} model={r['model']} A_s={r['A_s']:.12f} X_s={r['X_s']:.12f} "
            f"A_f={r['A_f']:.12f} X_f={r['X_f']:.12f} max_refit_delta={r['max_refit_vs_frozen']:.3e}"
            for r in rows
        )
        + f"\nclassical-vs-absolute fitted-core difference (N=3200, eval12800)={e:.15e}\n"
        + f"difference from frozen pair table={diff:.3e}\n"
        + 'DWENO_N3200_RAW_CHECK_OK\n'
    )
    print(report.read_text(), end='')


if __name__ == '__main__':
    main()
