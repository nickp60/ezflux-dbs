#!/usr/bin/env python3
"""
Extract AGORA Supplementary Table 12 (Western and High-fiber diets, mmol/gDW/h)
from the PDF supplied with Magnusdottir et al. 2017 (Nat Biotechnol 35:81-89,
nbt.3703, MOESM2) into a tidy CSV.

Usage
-----
    python extract_s12_table.py  AGORA_S12_diets_nbt3703_MOESM2.pdf  agora_s12_table12.csv

Needs the `pdftotext` command line tool (poppler-utils).  The PDF is read with
`pdftotext -layout`; every table row looks like

        <pdf id>   <EX_reaction>   <name ...>   <Western>   <High fiber>

and section titles (Sugars / Fiber / Fat / Protein / Minerals, vitamins, other)
appear on their own lines.  Output columns:

    table_category, pdf_id, reaction, metabolite, name, western, high_fiber

`metabolite` is the harmonised ID used by the other scripts (EX_ prefix and
(e)/[e]/_e suffix stripped, '__' -> '_').  The script checks the number of rows
(164 expected) and warns otherwise.
"""
import re
import subprocess
import sys
from pathlib import Path

import pandas as pd

NUM = r"[0-9.]+(?:[eE][+-]?\d+)?"
ROW = re.compile(rf"^\s{{5,}}(\S+)\s+(EX_\S+)\s+(.*?)\s+({NUM})\s+({NUM})\s*$")
SECTIONS = ("Sugars", "Fiber", "Fat", "Protein", "Minerals, vitamins, other")
EXPECTED_ROWS = 164


def normalise_id(x: str) -> str:
    x = re.sub(r"^EX_", "", x)
    x = re.sub(r"(\(e\)|\[e\]|_e)$", "", x)
    return x.replace("__", "_")


def extract(pdf: Path) -> pd.DataFrame:
    text = subprocess.run(["pdftotext", "-layout", str(pdf), "-"], check=True,
                          capture_output=True, text=True).stdout
    rows, cat = [], None
    for line in text.splitlines():
        m = ROW.match(line)
        if m:
            rows.append((cat, m.group(1), m.group(2), normalise_id(m.group(2)),
                         m.group(3).strip(), float(m.group(4)), float(m.group(5))))
        elif line.strip() in SECTIONS:
            cat = line.strip()
    return pd.DataFrame(rows, columns=["table_category", "pdf_id", "reaction", "metabolite",
                                       "name", "western", "high_fiber"])


def main():
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    df = extract(Path(sys.argv[1]))
    if len(df) != EXPECTED_ROWS:
        print(f"WARNING: parsed {len(df)} rows, expected {EXPECTED_ROWS}", file=sys.stderr)
    if df.metabolite.duplicated().any():
        print("WARNING: duplicated metabolites:", df.metabolite[df.metabolite.duplicated()].tolist(), file=sys.stderr)
    df.to_csv(sys.argv[2], index=False)
    print(f"{len(df)} rows -> {sys.argv[2]}")
    print(df.table_category.value_counts().to_string())


if __name__ == "__main__":
    main()
