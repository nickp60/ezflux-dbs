# Benchmarking the medium completion

Runs `../05_complete_medium.py` with many parameter sets as a SLURM job array and collects the
results. Not part of `run_all.sh` (it is slow and optional).

## Files

- `params.txt`: one line per run. Each line is the argument string for `05_complete_medium.py`:
  the skeleton CSV, `--weights` or `--minimize-components`, `--max-added-import`, `--growth`,
  `--o2` and `--out-dir` (one folder per run). Paths are relative to this folder. The default file
  has 296 lines, all on `agora_s12_western`. The first 240: growth 0.001, 0.005, 0.01, 0.05 x O2 0, 0.001, 0.1 x
  five methods (`--minimize-components`, `--weights` None, mass, C, N) x `--max-added-import`
  0.5, 1, 5, 10. Run folders are named `growth-<g>_o2-<o>_<method>_cap-<c>`.
  Lines 241-296 repeat part of that grid with `--cofactor-max-flux` 0.1, 0.03, 0.01, 0.003, 0.001,
  0.0001 and 0 (O2 0.001; growth 0.001 and 0.01; `--minimize-components` and `--weights mass`; cap
  1 and 10), folders `..._cap-<c>_cof-<f>`. The matching runs without a limit are in the first 240 lines.
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
to its own folder under `results/benchmark_media/<diet>/<condition>/` (`<condition>` as above), so `media_out/` is
untouched. Task logs are `benchmark_<jobid>_<task>.out` in this folder.

The report has one row per finished run: its parameters (growth, o2, weights, minimize_components,
max_added_import), taxa that can grow, growth-rate statistics from
`check_db_medium`, added compounds, added flux (mmol/gDW/h), the skeleton's flux, added flux as a
percentage of it, and run times in seconds (`complete_s`, `check_s`, `total_s`). A run without a
summary failed or has not finished; read its log. `report.py` prints how many of the lines in
`params.txt` have a summary.

## Changing the runs

- Regenerate the grid (bash; edit the lists to change it):
  ```
  d=agora_s12_western
  for g in 0.001 0.005 0.01 0.05; do for o in 0 0.001 0.1; do
   for m in "mincomp:--minimize-components" "weights-None:--weights None" "weights-mass:--weights mass" "weights-C:--weights C" "weights-N:--weights N"; do
    for c in 0.5 1 5 10; do n=${m%%:*}; a=${m#*:}
     echo "../../media_out/${d}_skeleton.csv $a --max-added-import $c --growth $g --o2 $o --out-dir ../../results/benchmark_media/$d/growth-${g}_o2-${o}_${n}_cap-${c}"
   done; done; done; done > params.txt
  ```
- Oxygen: `--o2` replaces the `o2` flux of the skeleton for that run (03 sets 0.001) and makes oxygen
  strict in `complete_db_medium`, so completion cannot raise it above the value given. Runs without
  `--o2` use the skeleton's value and the old non-strict behaviour, so they are not directly
  comparable with runs that have it.

- Cofactor limit: Table 12 gives 86 of the 88 "Minerals, vitamins, other" rows a flux of 1.0, which
  makes B12 forms, quinones, hemes and similar compounds a large carbon source. `--cofactor-max-flux L`
  sets each carbon-bearing row of that category to min(flux, L) and makes those rows strict, so
  completion cannot raise them. Water and methanol are never changed; carbon-free ions and gases
  (`INORGANIC` in 05) are unchanged unless `--cofactor-include-inorganic` is given. `report.py` has a
  `cofactor_max_flux` column (empty = no limit). Look for the lowest limit at which `n_grow_check`
  stays at `n_taxa`, then look at which taxa drop out below it.
- Edit or add lines in `params.txt`, then set the array size to the number of lines:
  `sbatch --array=1-N%20 benchmarking.sh` (the `#SBATCH --array=1-296%40` line only fits the default file).
- Resubmitting is safe: a run whose `<diet>_run_summary.json` already exists is skipped (the task log says so), so after a failure or timeout `sbatch benchmarking.sh` only does the missing runs. To redo finished runs, use `sbatch --export=ALL,FORCE=1 benchmarking.sh`, or delete their folders. Run one line: `sbatch --array=7 benchmarking.sh` (add `FORCE=1` if it already finished).
- Different limits: pass them to `sbatch`, for example `sbatch --mem=16G --time=24:00:00 benchmarking.sh`
  (command-line options override the `#SBATCH` lines). `--minimize-components` runs may need more time.
- Test one line without SLURM, from this folder:
  `SLURM_ARRAY_TASK_ID=1 SLURM_SUBMIT_DIR=$PWD SLURM_CPUS_PER_TASK=1 bash benchmarking.sh`
