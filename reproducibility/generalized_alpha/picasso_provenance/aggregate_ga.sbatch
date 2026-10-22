#!/usr/bin/env bash
#SBATCH --job-name=GA_agg
#SBATCH --partition=cpu_parti
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:10:00
#SBATCH --output=logs/GA_agg-%j.out
#SBATCH --error=logs/GA_agg-%j.err
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
PYBIN=$(./select_python.sh)
export PYTHONPATH="$PWD/code"
"$PYBIN" code/aggregate_results.py
sha256sum results/*.csv > results/SHA256_RESULTS.txt
printf 'GA_CORRECTED_PIPELINE_OK\n'
