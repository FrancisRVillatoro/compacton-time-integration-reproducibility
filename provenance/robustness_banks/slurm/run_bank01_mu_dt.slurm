#!/bin/bash
#SBATCH --job-name=absK_b01
#SBATCH --partition=cpu_partition
#SBATCH --constraint=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=03:00:00
#SBATCH --array=1-1440%50
#SBATCH --output=/mnt/home/users/tic_118_uma/frv/fscratch/absK_Lstable_picasso/banks/bank01_mu_dt/logs/slurm-%A.out

module load matlab/R2024b

PROJECT=/mnt/home/users/tic_118_uma/frv/fscratch/absK_Lstable_picasso
CODEDIR=$PROJECT/code
CASES=$PROJECT/banks/bank01_mu_dt/cases/cases.csv

cd "$PROJECT"

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

if [ -n "${LOCALSCRATCH:-}" ] && [ -d "$LOCALSCRATCH" ]; then
    RUN_BASE="$LOCALSCRATCH/frv_absK_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
else
    RUN_BASE="/tmp/frv_absK_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
fi

mkdir -p "$RUN_BASE/home" "$RUN_BASE/tmp" "$RUN_BASE/matlab_pref" "$RUN_BASE/matlab_temp" "$RUN_BASE/java_tmp"

export REAL_HOME="$HOME"
export HOME="$RUN_BASE/home"
export TMPDIR="$RUN_BASE/tmp"
export TEMP="$RUN_BASE/tmp"
export TMP="$RUN_BASE/tmp"
export MATLAB_PREFDIR="$RUN_BASE/matlab_pref"
export MATLAB_LOG_DIR="$RUN_BASE/matlab_temp"
export MATLAB_SHELL_TMP="$RUN_BASE/matlab_temp"
export MATLABPATH="$CODEDIR"
export _JAVA_OPTIONS="-Djava.io.tmpdir=$RUN_BASE/java_tmp -Duser.home=$RUN_BASE/home"

echo "PROJECT=$PROJECT"
echo "CODEDIR=$CODEDIR"
echo "CASES=$CASES"
echo "RUN_BASE=$RUN_BASE"
echo "SLURM_ARRAY_TASK_ID=$SLURM_ARRAY_TASK_ID"

ls -lh "$CODEDIR/run_absK_case.m" || exit 2
ls -lh "$CASES" || exit 3

matlab -singleCompThread -nodisplay -nosplash -batch "try, addpath('$CODEDIR'); run_absK_case('$CASES', str2double(getenv('SLURM_ARRAY_TASK_ID'))); catch ME, disp(getReport(ME,'extended')); exit(1); end"

STATUS=$?

cd /
rm -rf "$RUN_BASE"

exit $STATUS
