#!/usr/bin/env python3
"""
Extract the medium table from MICOM medium .qza files (QIIME 2 `MicomMedium[Global]`
artifacts, zip archives containing <uuid>/data/medium.csv) into plain CSV files, so
nothing downstream needs QIIME 2 or q2-micom.

Usage
-----
    python 02_qza_to_csv.py ../data/inputs  ../data/media_csv

Writes <outdir>/<name>.csv for every <name>.qza (an 8-hex-digit upload prefix such as
"a8a87aae-" is stripped from the name).  The columns are those of the original table.
"""
import re
import sys
import zipfile
from pathlib import Path

if len(sys.argv) != 3:
    sys.exit(__doc__)
src, out = Path(sys.argv[1]), Path(sys.argv[2])
out.mkdir(parents=True, exist_ok=True)
for f in sorted(src.glob("*.qza")):
    with zipfile.ZipFile(f) as z:
        member = [n for n in z.namelist() if n.endswith("data/medium.csv")]
        assert len(member) == 1, f"{f.name}: expected one data/medium.csv, found {member}"
        name = re.sub(r"^[0-9a-f]{8}-", "", f.stem)
        (out / f"{name}.csv").write_bytes(z.read(member[0]))
        print(f"{f.name} -> {name}.csv")
