#!/usr/bin/env bash
#SBATCH --job-name=media_bench
#SBATCH --array=1-32%20
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=08:00:00
#SBATCH --output=benchmark_%A_%a.out
# Run line N of params.txt (N = SLURM_ARRAY_TASK_ID) through 05_complete_medium.py.
# Submit from this folder:  sbatch benchmarking.sh
# For a params.txt with a different number of lines:  sbatch --array=1-N%20 benchmarking.sh
set -euo pipefail

cd "$SLURM_SUBMIT_DIR"
# source activate micom    # or activate your environment before running sbatch
args=$(sed -n "${SLURM_ARRAY_TASK_ID}p" params.txt)

echo "task $SLURM_ARRAY_TASK_ID: python ../05_complete_medium.py $args"
python ../05_complete_medium.py $args --threads "$SLURM_CPUS_PER_TASK"
