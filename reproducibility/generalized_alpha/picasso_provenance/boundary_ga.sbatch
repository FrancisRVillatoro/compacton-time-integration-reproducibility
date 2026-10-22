#!/usr/bin/env bash
#SBATCH --job-name=GA_bound
#SBATCH --partition=cpu_parti
#SBATCH --array=0-7%4
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G
#SBATCH --time=01:00:00
#SBATCH --output=logs/GA_bound-%A_%a.out
#SBATCH --error=logs/GA_bound-%A_%a.err
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
PYBIN=$(./select_python.sh)
export PYTHONPATH="$PWD/code"
export GA_OUT="$PWD/results/boundary_groups"
"$PYBIN" code/run_boundary_group.py "$SLURM_ARRAY_TASK_ID"
