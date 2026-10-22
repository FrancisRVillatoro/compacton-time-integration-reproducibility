#!/bin/bash
#SBATCH --job-name=absK_b02
#SBATCH --partition=cpu_partition
#SBATCH --constraint=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --output=/mnt/home/users/tic_118_uma/frv/fscratch/absK_Lstable_picasso/banks/bank02_dx_refinement/logs/finalize-%j.out

set -u

PROJECT=/mnt/home/users/tic_118_uma/frv/fscratch/absK_Lstable_picasso
BANK=bank02_dx_refinement
BANKDIR="$PROJECT/banks/$BANK"

cd "$PROJECT"

echo "Finalizing $BANK"
date

NCASES=$(($(wc -l < "$BANKDIR/cases/cases.csv")-1))
NCSV=$(ls "$BANKDIR"/summary/run_*.csv 2>/dev/null | wc -l)
NMAT=$(ls "$BANKDIR"/raw/run_*.mat 2>/dev/null | wc -l)

echo "NCASES=$NCASES"
echo "NCSV=$NCSV"
echo "NMAT=$NMAT"

python3 - <<'PY'
import csv, glob, os, re, math
from collections import Counter, defaultdict

project = os.path.expanduser('~/fscratch/absK_Lstable_picasso')
bank = 'bank02_dx_refinement'
bank_dir = os.path.join(project, 'banks', bank)
summary_dir = os.path.join(bank_dir, 'summary')
cases_file = os.path.join(bank_dir, 'cases', 'cases.csv')

files = sorted(glob.glob(os.path.join(summary_dir, 'run_*.csv')))
out = os.path.join(summary_dir, 'summary_bank02_all.csv')

if not files:
    raise SystemExit('No summary CSV files found.')

with open(files[0], newline='') as f:
    header = next(csv.reader(f))

nrows = 0
statuses = Counter()
wall = []
groups = defaultdict(list)

def num(x):
    try:
        return float(x)
    except:
        return float('nan')

def pct(vals,p):
    vals = sorted(v for v in vals if not math.isnan(v))
    if not vals:
        return float('nan')
    k = (len(vals)-1)*p/100.0
    i = int(k)
    j = min(i+1, len(vals)-1)
    a = k-i
    return vals[i]*(1-a)+vals[j]*a

with open(out, 'w', newline='') as g:
    writer = csv.writer(g)
    writer.writerow(header)

    for fn in files:
        with open(fn, newline='') as f:
            r = csv.DictReader(f)
            if r.fieldnames != header:
                print('WARNING: different header in', fn)
                continue
            for row in r:
                writer.writerow([row[h] for h in header])
                nrows += 1
                statuses[row.get('status','')] += 1
                wt = num(row.get('wall_time','nan'))
                wall.append(wt)
                groups[row.get('method','')].append(wt)

expected = set()
with open(cases_file, newline='') as f:
    for row in csv.DictReader(f):
        expected.add(int(row['run_id']))

found = set()
for fn in files:
    m = re.match(r'run_(\d+)_', os.path.basename(fn))
    if m:
        found.add(int(m.group(1)))

missing = sorted(expected - found)

with open(os.path.join(summary_dir,'status_counts_bank02.txt'), 'w') as f:
    for k,v in sorted(statuses.items()):
        f.write('%s,%d\n' % (k,v))

with open(os.path.join(summary_dir,'missing_run_ids_bank02.txt'), 'w') as f:
    for rid in missing:
        f.write('%d\n' % rid)

timing = []
timing.append(f"N = {nrows}")
if wall:
    timing.append(f"total_core_hours = {sum(wall)/3600:.6g}")
    timing.append(f"mean_wall_per_run_s = {sum(wall)/len(wall):.6g}")
    for p in [50,75,90,95,99,100]:
        timing.append(f"p{p}_wall_s = {pct(wall,p):.6g}")
    for m in sorted(groups):
        vals = groups[m]
        timing.append(
            f"{m}: N={len(vals)}, mean={sum(vals)/len(vals):.6g} s, "
            f"p95={pct(vals,95):.6g} s, max={max(vals):.6g} s, "
            f"core_hours={sum(vals)/3600:.6g}"
        )

with open(os.path.join(summary_dir,'timing_summary_bank02.txt'), 'w') as f:
    f.write('\n'.join(timing) + '\n')

print('merged_file=', out)
print('files_merged=', len(files))
print('rows_written=', nrows)
print('expected=', len(expected))
print('found=', len(found))
print('missing=', len(missing))
print('missing_run_ids=', missing)
print('statuses=', dict(statuses))
print('\n'.join(timing))
PY

ERROR_REPORT="$BANKDIR/summary/error_scan_bank02.txt"

grep -R -n \
    -e "Undefined function" \
    -e "Error using" \
    -e "Out of memory" \
    -e "exceeded" \
    -e "License checkout failed" \
    -e "MATLAB:nomem" \
    -e "FAILED" \
    "$BANKDIR/logs"/*.out > "$ERROR_REPORT" 2>/dev/null || true

echo "Error report lines:"
wc -l "$ERROR_REPORT"
head -n 50 "$ERROR_REPORT" || true

echo "Disk usage:"
du -sh "$BANKDIR" || true
du -sh "$BANKDIR/summary" "$BANKDIR/raw" "$BANKDIR/logs" 2>/dev/null || true

ARCHIVE_DIR="$PROJECT/archives"
mkdir -p "$ARCHIVE_DIR"

STAMP=$(date +%Y%m%d_%H%M%S)
ARCHIVE="$ARCHIVE_DIR/${BANK}_COMPLETE_${STAMP}.tar.gz"

tar -czf "$ARCHIVE" -C "$PROJECT/banks" "$BANK"

echo "ARCHIVE=$ARCHIVE"
ls -lh "$ARCHIVE"

echo "Finalized $BANK"
date
