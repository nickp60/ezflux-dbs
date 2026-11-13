#!/usr/bin/env python3
"""
Benchmark the medium-completion parameters of 05_complete_medium.py.

Usage
-----
    python 06_benchmarking_media_params.py ../media_out/agora_s12_high_fiber_bmt_skeleton.csv \
        ../media_out/agora_s12_western_bmt_skeleton.csv --threads 12 -o ../results/benchmark_media

For every skeleton the script runs 05_complete_medium.py once per parameter set:

  * every --weights value x every --max-added-import value
  * with --minimize-components, once per --max-added-import value (weights are ignored there;
    switch off with --skip-minimize-components)

Defaults: --weights mass None C, --max-added-import 0.5 1 4 10.

Every run gets its own folder  <out>/<diet>/<tag>/  (medium CSV, check pickle, run summary), so
nothing in media_out/ is overwritten.  The report is written to

  <out>/benchmark_report.csv    one row per diet and parameter set
  <out>/benchmark_report.md     the same as a table

Columns: status (ok / failed / timeout), how many taxa can grow (completion and check),
growth-rate statistics from check_db_medium, number of added compounds, added flux (mmol/gDW/h,
sum over compounds not in the skeleton), skeleton flux, added flux as % of the skeleton flux, and
run times in seconds (complete_s and check_s are measured inside the run, wall_s includes
start-up).

--timeout limits one run (seconds); a run that exceeds it is reported as "timeout" and the
benchmark carries on.  --skip-existing reuses a run summary that is already there.  --dry-run
prints the commands only.
"""
import argparse
import json
import subprocess
import sys
import time
from pathlib import Path

import pandas as pd

here = Path(__file__).resolve().parent
ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("skeletons", nargs="+", help="media_out/<diet>_skeleton.csv files")
ap.add_argument("--weights", nargs="+", default=["mass", "None", "C"])
ap.add_argument("--max-added-import", nargs="+", type=float, default=[0.5, 1, 4, 10])
ap.add_argument("--skip-minimize-components", action="store_true")
ap.add_argument("--growth", type=float, default=0.05)
ap.add_argument("--threads", type=int, default=12)
ap.add_argument("--model-db", default=str(here / "../data/agora201_refseq216_species_1.qza"))
ap.add_argument("--timeout", type=float, default=None, help="seconds per run")
ap.add_argument("--skip-existing", action="store_true")
ap.add_argument("--dry-run", action="store_true")
ap.add_argument("-o", "--out-dir", default=str(here / "../results/benchmark_media"))
args = ap.parse_args()

out = Path(args.out_dir)

# parameter sets: (tag, extra arguments)
params = []
for w in args.weights:
    for cap in args.max_added_import:
        params.append((f"weights-{w}_cap-{cap:g}", ["--weights", w, "--max-added-import", str(cap)]))
if not args.skip_minimize_components:
    for cap in args.max_added_import:
        params.append((f"mincomp_cap-{cap:g}", ["--minimize-components", "--max-added-import", str(cap)]))

rows = []
for skeleton in args.skeletons:
    diet = Path(skeleton).stem.removesuffix("_skeleton")
    for tag, extra in params:
        run_dir = out / diet / tag
        cmd = [sys.executable, str(here / "05_complete_medium.py"), str(skeleton),
               "--growth", str(args.growth), "--threads", str(args.threads),
               "--model-db", str(args.model_db), "--out-dir", str(run_dir), *extra]
        print(f"[{diet} {tag}]", " ".join(cmd), flush=True)
        row = {"diet": diet, "params": tag, "status": "dry-run" if args.dry_run else "ok"}
        summary_file = run_dir / f"{diet}_run_summary.json"
        if args.dry_run:
            rows.append(row)
            continue
        if args.skip_existing and summary_file.exists():
            print("  summary exists, skipped", flush=True)
            row["wall_s"] = None
        else:
            run_dir.mkdir(parents=True, exist_ok=True)
            t0 = time.time()
            try:
                proc = subprocess.run(cmd, capture_output=True, text=True, timeout=args.timeout)
                row["wall_s"] = time.time() - t0
                (run_dir / "log.txt").write_text(proc.stdout + "\n" + proc.stderr)
                if proc.returncode != 0:
                    row["status"] = "failed"
                    row["error"] = (proc.stderr.strip().splitlines() or ["no message"])[-1]
            except subprocess.TimeoutExpired:
                row["wall_s"] = time.time() - t0
                row["status"] = "timeout"
        if row["status"] == "ok" and summary_file.exists():
            summary = json.loads(summary_file.read_text())
            for k in ("n_taxa", "n_can_grow_completion", "n_grow_check", "growth_rate_min", "growth_rate_q25",
                      "growth_rate_median", "growth_rate_mean", "growth_rate_max", "n_added_compounds",
                      "added_flux", "skeleton_flux", "added_pct_of_skeleton", "complete_s", "check_s"):
                row[k] = summary[k]
        print(f"  {row['status']}", flush=True)
        rows.append(row)

report = pd.DataFrame(rows)
out.mkdir(parents=True, exist_ok=True)
report.to_csv(out / "benchmark_report.csv", index=False)

# markdown table (no extra packages needed)
table = report.copy()
for c in table.columns:
    if table[c].dtype == float:
        table[c] = table[c].map(lambda v: "" if pd.isna(v) else f"{v:.6g}")
table = table.fillna("")
lines = ["| " + " | ".join(table.columns) + " |", "|" + "|".join("---" for _ in table.columns) + "|"]
lines += ["| " + " | ".join(str(v) for v in r) + " |" for r in table.itertuples(index=False)]
(out / "benchmark_report.md").write_text("\n".join(lines) + "\n")

print(f"\n{len(report)} run(s); report in {out / 'benchmark_report.csv'} and benchmark_report.md")
print(report.drop(columns=[c for c in ("error",) if c in report]).to_string(index=False))
