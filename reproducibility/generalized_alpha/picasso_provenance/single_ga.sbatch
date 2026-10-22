#!/usr/bin/env bash
#SBATCH --job-name=GA_single
#SBATCH --partition=cpu_parti
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G
#SBATCH --time=01:00:00
#SBATCH --output=logs/GA_single-%j.out
#SBATCH --error=logs/GA_single-%j.err
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
PYBIN=$(./select_python.sh)
export PYTHONPATH="$PWD/code"
export GA_OUT="$PWD/results"
"$PYBIN" code/run_single_benchmark.py
