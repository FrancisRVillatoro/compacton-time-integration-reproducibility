#!/usr/bin/env python3
from pathlib import Path
import hashlib, subprocess, sys, os

ROOT = Path(__file__).resolve().parents[1]


def sha(p):
    h = hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()


expected = {
    'cc_N800_ext0.npz': '29babc260877594c9f041ffb05dd0853e0161facc0998d6e26a7057de4d7d921',
    'cc_N800_ext1.npz': 'c214e367bba3e383b85c80881809f900899a34b7a30437662dc375d472be6a5d',
    'k22ext_cc_N1600_ext1.npz': '95603ba5531ae0fd90fdc5a9d08c1c54bba3f039b2db85691ed2cf27567b2eac',
    'k22ext_cc_N3200_ext0.npz': '050b4a9c64f9c10cd90df60765cfcbe457af3d40b59b7f85befe0ce2db302fc8',
    'k22ext_cc_N3200_ext1.npz': 'a70fbfa816b97f07535fc203ca56a43a5e2b184618e51ee3e4f0ab08b0550560',
}

for name, digest in expected.items():
    p = ROOT / 'raw' / 'direct_weno' / name
    got = sha(p)
    print(name, got, 'OK' if got == digest else 'FAIL')
    if got != digest:
        sys.exit(2)

steps = [
    'recompute_dweno_N3200_raw_checks.py',
    'recompute_absolute_spatial_from_raw.py',
    'recompute_table3_absolute_reference.py',
    'recompute_i2_table3.py',
    'build_final_table3.py',
    'make_fig5_spatial_convergence.py',
    'make_fig7_accuracy_cost.py',
]
env = dict(os.environ)
env['PYTHONDONTWRITEBYTECODE'] = '1'
for script in steps:
    subprocess.check_call([sys.executable, str(ROOT / 'scripts' / script)], env=env)

print('FREEZE_CHECK_OK')
