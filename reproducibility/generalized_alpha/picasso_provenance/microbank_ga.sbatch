#!/usr/bin/env bash
#SBATCH --job-name=GA_micro
#SBATCH --partition=cpu_parti
#SBATCH --array=0-35%12
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=logs/GA_micro-%A_%a.out
#SBATCH --error=logs/GA_micro-%A_%a.err
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
PYBIN=$(./select_python.sh)
export PYTHONPATH="$PWD/code"
export GA_OUT="$PWD/results/microbank_groups"
"$PYBIN" code/run_microbank_group.py "$SLURM_ARRAY_TASK_ID"
