# Benchmarking the medium completion

Runs `../05_complete_medium.py` with many parameter sets as a SLURM job array and collects the
results. Not part of `run_all.sh` (it is slow and optional).

## Files

- `params.txt`: one line per run. Each line is the argument string for `05_complete_medium.py`:
  the skeleton CSV, `--weights` or `--minimize-components`, `--max-added-import`, and
  `--out-dir` (one folder per run). Paths are relative to this folder. The default file has 32
  lines: two diets (`agora_s12_high_fiber_bmt`, `agora_s12_western_bmt`) x 16 conditions
  (`weights` mass, None and C, each at caps 0.5, 1, 4 and 10, plus `--minimize-components` at each cap).
- `benchmarking.sh`: runs line N of `params.txt`, where N is `SLURM_ARRAY_TASK_ID`, through
  `05_complete_medium.py` with `--threads` set to the cores SLURM gave the task. Defaults, set in
  the `#SBATCH` lines at the top: at most 20 tasks at once, 1 core, 8 GB, 8 hours per task.
- `report.py`: collects the run summaries into `benchmark_report.csv` and `.md`.

## Usage

Activate the environment with `micom` first (tasks inherit it), then, from this folder:

```
sbatch benchmarking.sh
```

When all tasks have finished:

```
python report.py
```

Each run writes `<diet>_micom_medium.csv`, `<diet>_medium_check.pkl` and `<diet>_run_summary.json`
to its own folder under `results/benchmark_media/<diet>/<condition>/`, so `media_out/` is
untouched. Task logs are `benchmark_<jobid>_<task>.out` in this folder.

The report has one row per finished run: taxa that can grow, growth-rate statistics from
`check_db_medium`, added compounds, added flux (mmol/gDW/h), the skeleton's flux, added flux as a
percentage of it, and run times in seconds (`complete_s`, `check_s`, `total_s`). A run without a
summary failed or has not finished; read its log. `report.py` prints how many of the lines in
`params.txt` have a summary.

## Changing the runs

- Edit or add lines in `params.txt`, then set the array size to the number of lines:
  `sbatch --array=1-N%20 benchmarking.sh` (the `#SBATCH --array=1-32%20` line only fits the default file).
- Run one line again: `sbatch --array=7 benchmarking.sh`.
- Different limits: pass them to `sbatch`, for example `sbatch --mem=16G --time=24:00:00 benchmarking.sh`
  (command-line options override the `#SBATCH` lines). `--minimize-components` runs may need more time.
- Test one line without SLURM, from this folder:
  `SLURM_ARRAY_TASK_ID=1 SLURM_SUBMIT_DIR=$PWD SLURM_CPUS_PER_TASK=1 bash benchmarking.sh`
