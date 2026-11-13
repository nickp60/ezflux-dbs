#!/usr/bin/env python3
"""
Collect the results of 06_benchmarking_media_params.sh into a report.

Usage
-----
    python 06_benchmark_report.py ../results/benchmark_media

Reads <out-dir>/conditions.tsv and, for every condition, the files that 05_complete_medium.py and
the array task wrote into <out-dir>/<diet>/<tag>/.  Writes <out-dir>/benchmark_report.csv and
benchmark_report.md, one row per condition:

  status        ok / failed (with the last line of the log in the CSV column `error`) /
                no result (the task was killed, for example by the time limit, or never ran)
  taxa          taxa in the database, how many can grow after completion, how many reach the growth
                rate in check_db_medium
  growth rate   min, 25%, median, mean and max from check_db_medium
  added flux    added compounds, added flux (mmol/gDW/h, sum over compounds not in the skeleton),
                the skeleton's flux, and added flux as % of it
  run times     complete_s and check_s (inside 05_complete_medium.py) and wall_s (whole task)
"""
import json
import sys
from pathlib import Path

import pandas as pd

KEYS = ("n_taxa", "n_can_grow_completion", "n_grow_check", "growth_rate_min", "growth_rate_q25",
        "growth_rate_median", "growth_rate_mean", "growth_rate_max", "n_added_compounds",
        "added_flux", "skeleton_flux", "added_pct_of_skeleton", "complete_s", "check_s")


def main():
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    out = Path(sys.argv[1])
    conditions = pd.read_csv(out / "conditions.tsv", sep="\t", header=None,
                             names=["diet", "tag", "skeleton", "extra"])
    rows = []
    for c in conditions.itertuples():
        run_dir = out / c.diet / c.tag
        row = {"diet": c.diet, "params": c.tag}
        summary_file = run_dir / f"{c.diet}_run_summary.json"
        status_file = run_dir / "status.txt"
        wall_file = run_dir / "wall_s.txt"
        row["wall_s"] = float(wall_file.read_text()) if wall_file.exists() else None
        if summary_file.exists():
            row["status"] = "ok"
            summary = json.loads(summary_file.read_text())
            row.update({k: summary[k] for k in KEYS})
        elif status_file.exists() and status_file.read_text().strip() == "failed":
            row["status"] = "failed"
            log = (run_dir / "log.txt").read_text().strip().splitlines()
            row["error"] = log[-1] if log else "no message"
        else:
            row["status"] = "no result"
        rows.append(row)

    report = pd.DataFrame(rows)
    report.to_csv(out / "benchmark_report.csv", index=False)

    # markdown table (no extra packages needed)
    table = report.copy()
    for col in table.columns:
        if table[col].dtype == float:
            table[col] = table[col].map(lambda v: "" if pd.isna(v) else f"{v:.6g}")
    table = table.fillna("")
    lines = ["| " + " | ".join(table.columns) + " |", "|" + "|".join("---" for _ in table.columns) + "|"]
    lines += ["| " + " | ".join(str(v) for v in r) + " |" for r in table.itertuples(index=False)]
    (out / "benchmark_report.md").write_text("\n".join(lines) + "\n")

    print(f"{len(report)} condition(s); report in {out / 'benchmark_report.csv'} and benchmark_report.md")
    print(report.drop(columns=["error"], errors="ignore").to_string(index=False))


if __name__ == "__main__":
    main()
