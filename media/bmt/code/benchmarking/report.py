#!/usr/bin/env python3
"""
Collect the run summaries written by 05_complete_medium.py into one table.

Usage (from this folder, after the array jobs have finished)
    python report.py [../../results/benchmark_media]

Reads <results>/<diet>/<condition>/<diet>_run_summary.json and writes benchmark_report.csv and
benchmark_report.md to <results>.  One row per finished condition: taxa that can grow, growth-rate
statistics from check_db_medium, added compounds, added flux (mmol/gDW/h), the skeleton's flux,
added flux as % of it, and run times in seconds (complete_s, check_s, total_s).  A condition
without a summary failed or has not finished; look in its benchmark_<jobid>_<task>.out log.
"""
import json
import sys
from pathlib import Path

import pandas as pd

here = Path(__file__).resolve().parent
results = Path(sys.argv[1]) if len(sys.argv) > 1 else here / "../../results/benchmark_media"

rows = []
for f in sorted(results.glob("*/*/*_run_summary.json")):
    row = json.loads(f.read_text())
    row["condition"] = f.parent.name
    row["total_s"] = row["complete_s"] + row["check_s"]
    rows.append(row)
report = pd.DataFrame(rows).drop(columns=["weights", "minimize_components", "max_added_import", "growth"])
report = report[["diet", "condition"] + [c for c in report.columns if c not in ("diet", "condition")]]
report.to_csv(results / "benchmark_report.csv", index=False)

table = report.copy()
for col in table.columns:
    if table[col].dtype == float:
        table[col] = table[col].map(lambda v: f"{v:.6g}")
lines = ["| " + " | ".join(table.columns) + " |", "|" + "|".join("---" for _ in table.columns) + "|"]
lines += ["| " + " | ".join(str(v) for v in r) + " |" for r in table.itertuples(index=False)]
(results / "benchmark_report.md").write_text("\n".join(lines) + "\n")

expected = len((here / "params.txt").read_text().splitlines())
print(f"{len(report)} of {expected} conditions have a summary; report in {results / 'benchmark_report.csv'}")
print(report.to_string(index=False))
