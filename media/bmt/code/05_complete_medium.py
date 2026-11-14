#!/usr/bin/env python3
"""
Complete ONE diet skeleton into a MICOM medium and validate it.

Usage
-----
    python 05_complete_medium.py ../media_out/agora_s12_high_fiber_bmt_skeleton.csv \
        --weights mass --max-added-import 10 --threads 12

Input is a <diet>_skeleton.csv from 03_make_micom_media.py (columns metabolite,
flux, reaction, global_id, ...).  The script runs, in this order,

  1. micom.batch.db_media.complete_db_medium on the model database
  2. micom.batch.db_media.check_db_medium on the completed medium

and writes to --out-dir (default ../media_out):

  <diet>_micom_medium.csv     completed medium (diet rows keep their skeleton flux, added
                              rows get the largest import any model used)
  <diet>_medium_check.pkl     the check_db_medium result (pickle)
  <diet>_run_summary.json     parameters, growth statistics, added flux and run times

Parameters
----------
--weights               "mass", an element such as "C", or "None" (every metabolite weighs the
                        same).  Ignored by micom if --minimize-components is given.
--minimize-components   minimise the number of added compounds instead of the added flux
                        (can be slow)
--max-added-import      largest import (mmol/gDW/h) completion may add per compound
--growth                growth rate every taxon must reach (default 0.05)
--o2                    import of oxygen (mmol/gDW/h) for this run, replacing the "o2" row of the
                        skeleton (03_make_micom_media.py sets 0.001).  Default: keep the skeleton's
                        value.  When given, oxygen is also made strict in complete_db_medium, so
                        completion cannot raise it above this value (otherwise it could, up to
                        --max-added-import, and the export would then drop the extra).
--threads               worker processes
--model-db              model database .qza (default ../data/agora201_refseq216_species_1.qza)

All work happens inside main(), which keeps the multiprocessing used by micom (threads > 1)
safe on platforms that start worker processes by importing this file.

Annotations for the added compounds are taken from the provided media in --media-dir and from
Table 12 (--s12-csv), as in 03_make_micom_media.py.  Nothing from QIIME 2 is used.
"""
import argparse
import json
import pickle
import time
from pathlib import Path

import pandas as pd
from micom.batch.db_media import check_db_medium, complete_db_medium

here = Path(__file__).resolve().parent


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("skeleton", help="media_out/<diet>_skeleton.csv")
    ap.add_argument("--weights", default="mass", help='"mass", an element like "C", or "None" (default: mass)')
    ap.add_argument("--minimize-components", action="store_true", help="minimise the number of added compounds")
    ap.add_argument("--max-added-import", type=float, default=100.0)
    ap.add_argument("--growth", type=float, default=0.05)
    ap.add_argument("--o2", type=float, default=None,
                    help='oxygen import for this run (default: keep the skeleton\'s "o2" flux)')
    ap.add_argument("--threads", type=int, default=12)
    ap.add_argument("--model-db", default=str(here / "../data/agora201_refseq216_species_1.qza"))
    ap.add_argument("--media-dir", default=str(here / "../data/media_csv"))
    ap.add_argument("--s12-csv", default=str(here / "../data/agora_s12_table12.csv"))
    ap.add_argument("--out-dir", default=str(here / "../media_out"))
    args = ap.parse_args()

    weights = None if args.weights == "None" else args.weights
    diet = Path(args.skeleton).stem.removesuffix("_skeleton")
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    skel = pd.read_csv(args.skeleton)
    strict = []
    if args.o2 is not None:
        is_o2 = skel.metabolite == "o2"
        if not is_o2.any():
            raise SystemExit(f"--o2 given but {args.skeleton} has no 'o2' row")
        skel.loc[is_o2, "flux"] = args.o2
        strict = skel.loc[is_o2, "global_id"].tolist()  # completion must not raise oxygen

    # annotation table for the added compounds (same as in 03_make_micom_media.py)
    media = [pd.read_csv(f) for f in sorted(Path(args.media_dir).glob("*.csv"))]
    annotations = (pd.concat([m for m in media if "name" in m])
                   .loc[:, ["metabolite", "name", "hmdb", "kegg.compound", "pubchem.compound", "inchi", "chebi"]]
                   .dropna(subset=["name"]).drop_duplicates("metabolite"))
    s12 = pd.read_csv(args.s12_csv)
    annotations = pd.concat([annotations, s12.loc[~s12.metabolite.isin(annotations.metabolite), ["metabolite", "name"]]])

    # 1. complete
    t0 = time.time()
    manifest, imports = complete_db_medium(str(args.model_db), skel, growth=args.growth, threads=args.threads,
                                           max_added_import=args.max_added_import,
                                           minimize_components=args.minimize_components, weights=weights,
                                           strict=strict)
    complete_s = time.time() - t0
    print(diet, manifest.can_grow.value_counts().to_dict())

    filled = imports.max()
    added = filled[~filled.index.isin(skel.reaction)]
    print(f"Added flux is {added.sum():.2f}/{filled.sum():.2f} mmol/h.")
    print(added.sort_values(ascending=False)[0:20])

    added_df = added.reset_index()
    added_df.iloc[:, 0] = added_df.iloc[:, 0].str.replace("EX_|_m$", "", regex=True)
    added_df.columns = ["metabolite", "flux"]
    added_df = pd.concat([skel[["metabolite", "flux"]], added_df])
    medium = pd.merge(added_df, annotations, on="metabolite", how="left")
    medium["reaction"] = "EX_" + medium.metabolite + "_m"
    medium["global_id"] = "EX_" + medium.metabolite + "(e)"

    # 2. validate
    t0 = time.time()
    check = check_db_medium(str(args.model_db), medium=medium, threads=args.threads)
    check_s = time.time() - t0
    print(diet)
    print(check.growth_rate.describe())

    # 3. save (sorted, flux rounded to 6 significant figures so that diffs between runs stay readable)
    medium = medium.sort_values("metabolite").reset_index(drop=True)
    medium["flux"] = medium.flux.map(lambda v: float(f"{v:.6g}"))
    medium.to_csv(out / f"{diet}_micom_medium.csv", index=False)
    with open(out / f"{diet}_medium_check.pkl", "wb") as outf:
        pickle.dump(check, outf)

    gr = check.growth_rate
    summary = {
        "diet": diet,
        "weights": None if args.minimize_components else weights,
        "minimize_components": args.minimize_components,
        "max_added_import": args.max_added_import,
        "growth": args.growth,
        "o2": float(skel.loc[skel.metabolite == "o2", "flux"].iloc[0]) if (skel.metabolite == "o2").any() else None,
        "n_taxa": int(len(manifest)),
        "n_can_grow_completion": int(manifest.can_grow.sum()),
        "n_grow_check": int((gr >= args.growth - 1e-6).sum()),
        "growth_rate_min": float(gr.min()),
        "growth_rate_q25": float(gr.quantile(0.25)),
        "growth_rate_median": float(gr.median()),
        "growth_rate_mean": float(gr.mean()),
        "growth_rate_max": float(gr.max()),
        "n_added_compounds": int((added > 0).sum()),
        "skeleton_flux": float(skel.flux.sum()),
        "added_flux": float(added.sum()),
        "added_pct_of_skeleton": float(100 * added.sum() / skel.flux.sum()),
        "complete_s": complete_s,
        "check_s": check_s,
    }
    with open(out / f"{diet}_run_summary.json", "w") as outf:
        json.dump(summary, outf, indent=2)
    print(f"{diet}: wrote medium, check pickle and summary to {out}")


if __name__ == "__main__":
    main()
