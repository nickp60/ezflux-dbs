#!/usr/bin/env bash
# Benchmark the medium-completion parameters of 05_complete_medium.py as a SLURM job array:
# one array task per condition, each running 05_complete_medium.py with a single thread.
#
# Usage
#   bash 06_benchmarking_media_params.sh [options] ../media_out/*_skeleton.csv
#
# Conditions (every skeleton x every parameter set)
#   * every --weights value x every --caps value
#   * --minimize-components once per --caps value (weights are ignored there)
#
# Options
#   --weights "mass None C"          weights values (default: "mass None C")
#   --caps "0.5 1 4 10"              max_added_import values (default: "0.5 1 4 10")
#   --skip-minimize-components       leave out the --minimize-components conditions
#   --growth 0.05                    growth rate every taxon must reach
#   --model-db PATH                  default ../data/agora201_refseq216_species_1.qza
#   --out-dir DIR                    default ../results/benchmark_media
#   --time HH:MM:SS                  time limit per task (default 24:00:00)
#   --mem 16G                        memory per task (default 16G)
#   --partition NAME                 sbatch partition (default: cluster default)
#   --max-parallel N                 run at most N array tasks at once (default: no limit)
#   --setup "COMMAND"                run in every task first, e.g. "conda activate micom"
#                                    (the submitting shell's environment is inherited anyway)
#   --skip-existing                  tasks whose run summary already exists do nothing
#   --dry-run                        print the conditions and the sbatch commands only
#
# Layout of --out-dir
#   conditions.tsv                   one line per array task: diet, tag, skeleton, extra arguments
#   slurm/<jobid>_<task>.out         task logs
#   <diet>/<tag>/                    medium CSV, check pickle, run summary, wall_s.txt, status.txt
#   benchmark_report.csv / .md       written by 06_benchmark_report.py in a final job that waits
#                                    for the array (also runnable by hand:
#                                    python 06_benchmark_report.py ../results/benchmark_media)
#
# The same file is the submitter and the task: when SLURM_ARRAY_TASK_ID is set it runs one
# condition, otherwise it submits.
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if [ -n "${SLURM_ARRAY_TASK_ID:-}" ]; then
    # ---- task mode: run one condition -------------------------------------------------------
    if [ -n "${BENCH_SETUP:-}" ]; then
        eval "$BENCH_SETUP"
    fi
    IFS=$'\t' read -r diet tag skeleton extra < <(sed -n "${SLURM_ARRAY_TASK_ID}p" "$BENCH_CONDITIONS")
    run_dir="$BENCH_OUT/$diet/$tag"
    mkdir -p "$run_dir"
    if [ "${BENCH_SKIP_EXISTING:-0}" = "1" ] && [ -f "$run_dir/${diet}_run_summary.json" ]; then
        echo "$diet $tag: summary exists, skipped"
        exit 0
    fi
    rm -f "$run_dir/status.txt"
    start=$(date +%s)
    status=ok
    # $extra holds the weights / cap / minimize arguments, so it must split into words
    python "$here/05_complete_medium.py" "$skeleton" --growth "$BENCH_GROWTH" --threads 1 \
        --model-db "$BENCH_MODEL_DB" --out-dir "$run_dir" $extra > "$run_dir/log.txt" 2>&1 || status=failed
    echo $(( $(date +%s) - start )) > "$run_dir/wall_s.txt"
    echo "$status" > "$run_dir/status.txt"
    echo "$diet $tag: $status"
    [ "$status" = ok ]
    exit $?
fi

# ---- submit mode ----------------------------------------------------------------------------
weights="mass None C"
caps="0.5 1 4 10"
minimize=1
growth=0.05
model_db="$here/../data/agora201_refseq216_species_1.qza"
out_dir="$here/../results/benchmark_media"
time_limit="24:00:00"
mem="16G"
partition=""
max_parallel=""
setup=""
skip_existing=0
dry_run=0
skeletons=()

while [ $# -gt 0 ]; do
    case "$1" in
        --weights) weights="$2"; shift 2 ;;
        --caps) caps="$2"; shift 2 ;;
        --skip-minimize-components) minimize=0; shift ;;
        --growth) growth="$2"; shift 2 ;;
        --model-db) model_db="$2"; shift 2 ;;
        --out-dir) out_dir="$2"; shift 2 ;;
        --time) time_limit="$2"; shift 2 ;;
        --mem) mem="$2"; shift 2 ;;
        --partition) partition="$2"; shift 2 ;;
        --max-parallel) max_parallel="$2"; shift 2 ;;
        --setup) setup="$2"; shift 2 ;;
        --skip-existing) skip_existing=1; shift ;;
        --dry-run) dry_run=1; shift ;;
        -h|--help) sed -n '2,37p' "${BASH_SOURCE[0]}"; exit 0 ;;
        -*) echo "unknown option $1" >&2; exit 2 ;;
        *) skeletons+=("$1"); shift ;;
    esac
done
if [ ${#skeletons[@]} -eq 0 ]; then
    echo "no skeleton files given (see --help)" >&2
    exit 2
fi

mkdir -p "$out_dir"
out_dir="$(cd "$out_dir" && pwd)"
conditions="$out_dir/conditions.tsv"
: > "$conditions"
for skeleton in "${skeletons[@]}"; do
    skeleton="$(cd "$(dirname "$skeleton")" && pwd)/$(basename "$skeleton")"
    diet="$(basename "$skeleton" _skeleton.csv)"
    for w in $weights; do
        for cap in $caps; do
            printf '%s\t%s\t%s\t%s\n' "$diet" "weights-${w}_cap-${cap}" "$skeleton" "--weights $w --max-added-import $cap" >> "$conditions"
        done
    done
    if [ "$minimize" = 1 ]; then
        for cap in $caps; do
            printf '%s\t%s\t%s\t%s\n' "$diet" "mincomp_cap-${cap}" "$skeleton" "--minimize-components --max-added-import $cap" >> "$conditions"
        done
    fi
done
n=$(wc -l < "$conditions")
echo "$n conditions in $conditions"

array="1-$n"
if [ -n "$max_parallel" ]; then
    array="$array%$max_parallel"
fi
# the tasks read their settings from the environment (sbatch passes it on; --export would split
# values at commas)
export BENCH_CONDITIONS="$conditions" BENCH_OUT="$out_dir" BENCH_GROWTH="$growth"
export BENCH_MODEL_DB="$model_db" BENCH_SETUP="$setup" BENCH_SKIP_EXISTING="$skip_existing"
sbatch_opts=(--parsable --job-name=media_bench --array="$array" --cpus-per-task=1 --ntasks=1
             --mem="$mem" --time="$time_limit" --output="$out_dir/slurm/%A_%a.out")
if [ -n "$partition" ]; then
    sbatch_opts+=(--partition="$partition")
fi

if [ "$dry_run" = 1 ]; then
    cat "$conditions"
    echo "sbatch ${sbatch_opts[*]} $here/06_benchmarking_media_params.sh"
    echo "sbatch --dependency=afterany:<jobid> --wrap 'python $here/06_benchmark_report.py $out_dir'"
    exit 0
fi

mkdir -p "$out_dir/slurm"
jobid=$(sbatch "${sbatch_opts[@]}" "$here/06_benchmarking_media_params.sh")
echo "submitted array job $jobid"

report_opts=(--parsable --job-name=media_bench_report --dependency="afterany:$jobid" --cpus-per-task=1
             --mem=2G --time=00:15:00 --output="$out_dir/slurm/report_%j.out")
if [ -n "$partition" ]; then
    report_opts+=(--partition="$partition")
fi
report_job=$(sbatch "${report_opts[@]}" --wrap "${setup:+$setup; }python $here/06_benchmark_report.py $out_dir")
echo "submitted report job $report_job (runs when the array has finished)"
