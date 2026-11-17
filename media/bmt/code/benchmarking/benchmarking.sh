#!/usr/bin/env bash
#SBATCH --job-name=media_bench
#SBATCH --array=1-296%40
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=log_benchmark_%A_%a.out
# Run line N of params.txt (N = SLURM_ARRAY_TASK_ID) through 05_complete_medium.py.
# Submit from this folder:  sbatch benchmarking.sh
# For a params.txt with a different number of lines:  sbatch --array=1-N%20 benchmarking.sh
# Runs whose summary file already exists are skipped, so an interrupted array can simply be resubmitted.
set -euo pipefail

# source activate micom    # or activate your environment before running sbatch
args=$(sed -n "${SLURM_ARRAY_TASK_ID}p" params.txt)

# Skip runs that already finished: 05 writes <out-dir>/<diet>_run_summary.json last.
# To redo them anyway:  sbatch --export=ALL,FORCE=1 benchmarking.sh
skeleton=${args%% *}
diet=$(basename "$skeleton" _skeleton.csv)
outdir=$(sed -n 's/.*--out-dir \([^ ]*\).*/\1/p' <<< "$args")
if [[ -n "$outdir" && -s "$outdir/${diet}_run_summary.json" && -z "${FORCE:-}" ]]; then
    echo "task $SLURM_ARRAY_TASK_ID: $outdir/${diet}_run_summary.json exists, skipping (FORCE=1 to rerun)"
    exit 0
fi

echo "task $SLURM_ARRAY_TASK_ID: python ../05_complete_medium.py $args"
python ../05_complete_medium.py $args --threads 1
