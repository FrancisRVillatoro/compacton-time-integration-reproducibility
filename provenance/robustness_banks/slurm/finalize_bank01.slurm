#!/bin/bash
#SBATCH --job-name=absK_fin01
#SBATCH --partition=cpu_partition
#SBATCH --constraint=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --output=/mnt/home/users/tic_118_uma/frv/fscratch/absK_Lstable_picasso/banks/bank01_mu_dt/logs/finalize-%j.out

set -u

PROJECT=/mnt/home/users/tic_118_uma/frv/fscratch/absK_Lstable_picasso
BANK=bank01_mu_dt
BANKDIR="$PROJECT/banks/$BANK"

cd "$PROJECT"

echo "Finalizing $BANK"
date

echo "---- Counts ----"
NCASES=$(($(wc -l < "$BANKDIR/cases/cases.csv")-1))
NCSV=$(ls "$BANKDIR"/summary/run_*.csv 2>/dev/null | wc -l)
NMAT=$(ls "$BANKDIR"/raw/run_*.mat 2>/dev/null | wc -l)

echo "NCASES=$NCASES"
echo "NCSV=$NCSV"
echo "NMAT=$NMAT"

echo "---- Merge CSV files ----"

python3 - <<'PY'
import csv
import glob
import os
import sys

project = os.path.expanduser('~/fscratch/absK_Lstable_picasso')
bank = 'bank01_mu_dt'
bank_dir = os.path.join(project, 'banks', bank)
summary_dir = os.path.join(bank_dir, 'summary')
files = sorted(glob.glob(os.path.join(summary_dir, 'run_*.csv')))
out = os.path.join(summary_dir, 'summary_bank01_all.csv')

if not files:
    raise SystemExit('No summary CSV files found.')

with open(files[0], newline='') as f:
    reader = csv.reader(f)
    header = next(reader)

nrows = 0
statuses = {}

with open(out, 'w', newline='') as g:
    writer = csv.writer(g)
    writer.writerow(header)

    try:
        status_idx = header.index('status')
    except ValueError:
        status_idx = None

    for fn in files:
        with open(fn, newline='') as f:
            reader = csv.reader(f)
            h = next(reader)
            if h != header:
                print('WARNING: different header in', fn)
                continue

            for row in reader:
                writer.writerow(row)
                nrows += 1
                if status_idx is not None and len(row) > status_idx:
                    statuses[row[status_idx]] = statuses.get(row[status_idx], 0) + 1

print('merged_file=', out)
print('files_merged=', len(files))
print('rows_written=', nrows)
print('statuses=', statuses)

status_out = os.path.join(summary_dir, 'status_counts_bank01.txt')
with open(status_out, 'w') as f:
    for k in sorted(statuses):
        f.write('%s,%d\n' % (k, statuses[k]))
print('status_counts_file=', status_out)
PY

echo "---- Search logs for MATLAB/Slurm errors ----"
ERROR_REPORT="$BANKDIR/summary/error_scan_bank01.txt"

grep -R -n \
    -e "Undefined function" \
    -e "Error using" \
    -e "Out of memory" \
    -e "exceeded" \
    -e "License checkout failed" \
    -e "MATLAB:nomem" \
    -e "FAILED" \
    "$BANKDIR/logs"/*.out > "$ERROR_REPORT" 2>/dev/null || true

echo "Error report:"
wc -l "$ERROR_REPORT"
head -n 50 "$ERROR_REPORT" || true

echo "---- Disk usage before compression ----"
du -sh "$BANKDIR" || true
du -sh "$BANKDIR/summary" "$BANKDIR/raw" "$BANKDIR/logs" 2>/dev/null || true

echo "---- Create manifest ----"
MANIFEST="$BANKDIR/summary/manifest_bank01.txt"
{
    echo "BANK=$BANK"
    echo "DATE=$(date)"
    echo "PROJECT=$PROJECT"
    echo "NCASES=$NCASES"
    echo "NCSV=$NCSV"
    echo "NMAT=$NMAT"
    echo "SUMMARY=$BANKDIR/summary/summary_bank01_all.csv"
    echo "ERROR_REPORT=$ERROR_REPORT"
    echo "DU_BANK=$(du -sh "$BANKDIR" | awk '{print $1}')"
} > "$MANIFEST"

cat "$MANIFEST"

echo "---- Compress complete bank output ----"
ARCHIVE_DIR="$PROJECT/archives"
mkdir -p "$ARCHIVE_DIR"

STAMP=$(date +%Y%m%d_%H%M%S)
ARCHIVE="$ARCHIVE_DIR/${BANK}_${STAMP}.tar.gz"

tar -czf "$ARCHIVE" \
    -C "$PROJECT/banks" "$BANK"

echo "ARCHIVE=$ARCHIVE"
ls -lh "$ARCHIVE"

echo "---- Final checks ----"
if [ "$NCSV" -ne "$NCASES" ]; then
    echo "WARNING: NCSV != NCASES"
    echo "Some runs may be missing."
fi

if [ -s "$ERROR_REPORT" ]; then
    echo "WARNING: possible errors found in logs."
else
    echo "No obvious MATLAB/Slurm errors found in logs."
fi

echo "Finalized $BANK"
date
