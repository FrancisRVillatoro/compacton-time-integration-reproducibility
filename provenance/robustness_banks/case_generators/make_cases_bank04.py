#!/usr/bin/env python3
import csv
import os

root = os.path.expanduser('~/fscratch/absK_Lstable_picasso')
bank = 'bank04_long_time'
bank_dir = os.path.join(root, 'banks', bank)
cases_dir = os.path.join(bank_dir, 'cases')
os.makedirs(cases_dir, exist_ok=True)

for sub in ['summary', 'raw', 'snapshots', 'logs']:
    os.makedirs(os.path.join(bank_dir, sub), exist_ok=True)

outfile = os.path.join(cases_dir, 'cases.csv')

p_vals = [2, 3]
mu_vals = [0, 1e-9, 1e-7, 1e-5, 1e-3, 1e-1, 1e0, 1e1]
dt_vals = [0.05, 0.1]
Tmax_vals = [200, 500, 1000]
methods = ['M1', 'M2', 'M3']

scenarios = [
    ('CC',  1.0),
    ('CA', -1.0),
]

fieldnames = [
    'run_id','bank','bank_dir','method','scenario',
    'p','dx','dt','mu','Lmax','Tmax','c0','c1','c2','x01','x02',
    'post_t0','save_raw','save_snapshots','nsave_target',
    'newton_tol','newton_max_iter','blowup_factor'
]

run_id = 0

with open(outfile, 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=fieldnames)
    w.writeheader()

    for p in p_vals:
        for mu in mu_vals:
            for dt in dt_vals:
                for Tmax in Tmax_vals:
                    for method in methods:
                        for scenario, c2 in scenarios:
                            run_id += 1

                            row = {
                                'run_id': run_id,
                                'bank': bank,
                                'bank_dir': bank_dir,
                                'method': method,
                                'scenario': scenario,
                                'p': p,
                                'dx': 0.05,
                                'dt': dt,
                                'mu': mu,
                                'Lmax': 100,
                                'Tmax': Tmax,
                                'c0': 1,
                                'c1': 1.25,
                                'c2': c2,
                                'x01': 30,
                                'x02': 65,
                                'post_t0': 60,
                                'save_raw': 1,
                                'save_snapshots': 0,
                                'nsave_target': 1000,
                                'newton_tol': 1e-10,
                                'newton_max_iter': 10,
                                'blowup_factor': 50,
                            }
                            w.writerow(row)

print(outfile)
print(run_id)
