#!/usr/bin/env python3
"""
Sanity checks for the skeleton media written by 03_make_micom_media.py
(media_out/<diet>_skeleton.csv), to be run BEFORE the slow complete_db_medium step.

Usage
-----
    python 04_check_skeletons.py ../media_out/*_skeleton.csv \
        --media-dir ../data/media_csv --s12-csv ../data/agora_s12_table12.csv

Exit status is 1 if any check fails, so it can be used in run_all.sh or CI.

Checks (E = error, W = warning)
-------------------------------
E  required columns exist; flux is finite and > 0; metabolite IDs are unique
E  reaction == EX_<id>_m and global_id == EX_<id>(e)
E  legacy IDs: a diet ID that is on the LEGACY_IDS list (glc -> glc_D), or that has the same
   name as another ID which the provided media use more often (catches IDs the models do not
   have: the Table 12 reaction column says EX_glc(e) but the AGORA exchange is EX_glc_D(e))
E  Table 12 table itself: the harmonised `metabolite` column agrees with the table's own ID
   column (ignoring letter case)
E  Table 12 diets (agora_s12_*): every table compound is present, and
   skeleton flux / table value is either 1 (not diluted) or the diet's dilution factor
   (0.1, or 0.8 for *_bmt); the diluted set is identical in all Table 12 diets, and the *_bmt
   and base diets differ only by the factor 8 for diluted compounds
W  IDs that none of the provided media contain (may be legitimate, but are worth a look)
"""
import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

LEGACY_IDS = {"glc": "glc_D"}           # legacy ID -> current AGORA ID
DILUTION = {"base": 0.1, "bmt": 0.8}    # must match `dilutions` in the notebook
errors, warnings = [], []


def err(diet, msg):
    errors.append(f"E [{diet}] {msg}")


def warn(diet, msg):
    warnings.append(f"W [{diet}] {msg}")


def norm_name(s):
    return re.sub(r"[^a-z0-9]", "", str(s).lower())


def strip_id(g):
    return re.sub(r"^EX_|\(e\)$|_e$", "", g).replace("__", "_")


ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("skeletons", nargs="+", help="media_out/<diet>_skeleton.csv files")
ap.add_argument("--media-dir", required=True, help="folder with the provided media as CSV (02_qza_to_csv.py output)")
ap.add_argument("--s12-csv", required=True, help="Table 12 CSV from 01_extract_s12_table.py")
args = ap.parse_args()

# --- reference: IDs, names and how many provided media use each ID ------------------------------
usage, names = {}, {}
for f in sorted(Path(args.media_dir).glob("*.csv")):
    m = pd.read_csv(f)
    ids = m["global_id"].map(strip_id)
    for i in set(ids):
        usage[i] = usage.get(i, 0) + 1
    if "name" in m:
        for i, n in zip(ids, m["name"]):
            if isinstance(n, str):
                names.setdefault(norm_name(n), set()).add(i)
if not usage:
    sys.exit(f"No media CSV found in {args.media_dir}")

# --- Table 12 itself ----------------------------------------------------------------------------
s12 = pd.read_csv(args.s12_csv)
bad = s12[s12.pdf_id.str.lower() != s12.metabolite.str.lower()]
for _, r in bad.iterrows():
    err("table12", f"metabolite {r.metabolite!r} (from {r.reaction}) differs from the table's own ID {r.pdf_id!r}")
if s12.metabolite.duplicated().any():
    err("table12", f"duplicated metabolites {s12.metabolite[s12.metabolite.duplicated()].tolist()}")

# --- per skeleton ---------------------------------------------------------------------------------
skel = {}
for f in args.skeletons:
    diet = re.sub(r"_skeleton$", "", Path(f).stem)
    d = pd.read_csv(f)
    skel[diet] = d
    missing = {"metabolite", "flux", "reaction", "global_id"} - set(d.columns)
    if missing:
        err(diet, f"missing columns {sorted(missing)}")
        continue
    if d.metabolite.duplicated().any():
        err(diet, f"duplicated metabolites {d.metabolite[d.metabolite.duplicated()].tolist()}")
    if not np.isfinite(d.flux).all() or (d.flux <= 0).any():
        err(diet, f"non-finite or non-positive flux: {d.metabolite[~(np.isfinite(d.flux) & (d.flux > 0))].tolist()}")
    wrong = d[(d.reaction != "EX_" + d.metabolite + "_m") | (d.global_id != "EX_" + d.metabolite + "(e)")]
    if len(wrong):
        err(diet, f"reaction/global_id do not follow the metabolite ID for {wrong.metabolite.tolist()}")

    for _, r in d.iterrows():
        if r.metabolite in LEGACY_IDS:
            err(diet, f"{r.metabolite} is a legacy ID; the AGORA exchange is {LEGACY_IDS[r.metabolite]}")
            continue
        for sib in names.get(norm_name(r.get("name")), ()):
            if sib != r.metabolite and usage.get(sib, 0) > usage.get(r.metabolite, 0):
                err(diet, f"{r.metabolite} has the same name ({r['name']!r}) as {sib}, which the provided media "
                          f"use more often ({usage[sib]} vs {usage.get(r.metabolite, 0)} media): wrong ID?")
    unknown = sorted(set(d.metabolite) - set(usage))
    if unknown:
        warn(diet, f"{len(unknown)} IDs are in none of the provided media: {unknown}")

# --- Table 12 diets -------------------------------------------------------------------------------
diluted = {}
for diet, d in skel.items():
    if not diet.startswith("agora_s12_") or "metabolite" not in d:
        continue
    col = "high_fiber" if "high_fiber" in diet else "western"
    factor = DILUTION["bmt" if diet.endswith("_bmt") else "base"]
    d = d.set_index("metabolite")
    table = s12.set_index("metabolite")[col]
    absent = sorted(set(table.index) - set(d.index))
    if absent:
        err(diet, f"Table 12 compounds missing from the skeleton: {absent}")
    common = table.index.intersection(d.index)
    ratio = d.flux[common] / table[common]
    is1 = np.isclose(ratio, 1, rtol=1e-6)
    isf = np.isclose(ratio, factor, rtol=1e-6)
    odd = common[~(is1 | isf)]
    if len(odd):
        err(diet, f"flux / Table 12 value is neither 1 nor {factor} for {odd.tolist()}")
    diluted[diet] = frozenset(common[isf & ~is1])

for col in ("western", "high_fiber"):
    base, bmt = f"agora_s12_{col}", f"agora_s12_{col}_bmt"
    if base in skel and bmt in skel and not skel[base].empty:
        b, t = skel[base].set_index("metabolite").flux, skel[bmt].set_index("metabolite").flux
        common = [m for m in s12.metabolite if m in b.index and m in t.index]
        r = (t[common] / b[common]).round(6)
        odd = r[~r.isin([1.0, round(DILUTION["bmt"] / DILUTION["base"], 6)])]
        if len(odd):
            err(bmt, f"differs from {base} by something other than x1 or x{DILUTION['bmt'] / DILUTION['base']:g}: {odd.to_dict()}")
sets = {k: v for k, v in diluted.items()}
if len(set(sets.values())) > 1:
    ref_name, ref = next(iter(sets.items()))
    for k, v in sets.items():
        if v != ref:
            err(k, f"diluted compounds differ from {ref_name}: only here {sorted(v - ref)}, only there {sorted(ref - v)}")

# --- report ---------------------------------------------------------------------------------------
for line in warnings + errors:
    print(line)
print(f"\n{len(skel)} skeleton(s) checked: {len(errors)} error(s), {len(warnings)} warning(s)")
sys.exit(1 if errors else 0)
