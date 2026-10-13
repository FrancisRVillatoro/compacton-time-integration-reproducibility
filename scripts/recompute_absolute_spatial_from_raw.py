#!/usr/bin/env python3
from __future__ import annotations

import csv
from pathlib import Path

from analyze_absolute_fullfield_cauchy import load_npz, compare, write_csv

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / 'raw' / 'direct_weno'
FROZEN = ROOT / 'data' / 'absolute_spatial_convergence.csv'
OUT = ROOT / 'data' / 'absolute_spatial_convergence_recomputed.csv'


def final_row(nc, nf, core):
    coarse = load_npz(RAW / ('cc_N800_ext1.npz' if nc == 800 else 'k22ext_cc_N1600_ext1.npz'))
    fine = load_npz(RAW / ('k22ext_cc_N1600_ext1.npz' if nf == 1600 else 'k22ext_cc_N3200_ext1.npz'))
    hist, _ = compare(coarse, fine)
    write_csv(hist, ROOT / 'data' / f'absolute_fullfield_{nc}_{nf}_history.csv')
    r = hist[-1]
    return dict(
        N_coarse=nc, N_fine=nf,
        raw_relL2=r['raw_relL2'],
        registered_relL2=r['registered_relL2'],
        global_shift=r['global_shift'],
        low16_registered=r['low_m16_registered'],
        low32_registered=r['low_m32_registered'],
        low64_registered=r['low_m64_registered'],
        low128_registered=r['low_m128_registered'],
        positive_rel_to_fullfine=r['positive_registered_rel_to_fullfine'],
        negative_rel_to_fullfine=r['negative_registered_rel_to_fullfine'],
        negative_L2_fraction_coarse=r['coarse_negative_L2_fraction'],
        negative_L2_fraction_fine=r['fine_negative_L2_fraction'],
        min_u_coarse=r['coarse_min_u'],
        min_u_fine=r['fine_min_u'],
        E_core=core,
    )


def main():
    pairs = {(800, 1600): 1.1755677210528387e-3, (1600, 3200): 4.426107800112411e-4}
    rows = [final_row(a, b, pairs[(a, b)]) for a, b in pairs]
    with OUT.open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)

    frozen = list(csv.DictReader(FROZEN.open()))
    keys = [
        'raw_relL2','registered_relL2','global_shift','low16_registered','low32_registered',
        'low64_registered','low128_registered','positive_rel_to_fullfine','negative_rel_to_fullfine',
        'negative_L2_fraction_coarse','negative_L2_fraction_fine','min_u_coarse','min_u_fine','E_core'
    ]
    maxdiff = 0.0
    for a, b in zip(rows, frozen):
        if int(b['N_coarse']) != int(a['N_coarse']) or int(b['N_fine']) != int(a['N_fine']):
            raise RuntimeError('pair ordering mismatch')
        for k in keys:
            maxdiff = max(maxdiff, abs(float(a[k]) - float(b[k])))
    print(f'max absolute difference from frozen absolute spatial table: {maxdiff:.3e}')
    if maxdiff > 2e-9:
        raise RuntimeError('absolute spatial raw recomputation mismatch')
    (ROOT / 'audit' / 'absolute_spatial_raw_check.txt').write_text(
        f'max absolute difference from frozen absolute spatial table: {maxdiff:.3e}\n'
        'ABSOLUTE_SPATIAL_RAW_CHECK_OK\n'
    )
    print('ABSOLUTE_SPATIAL_RAW_CHECK_OK')


if __name__ == '__main__':
    main()
