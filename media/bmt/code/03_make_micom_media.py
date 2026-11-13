#!/usr/bin/env python3
# MICOM media from AGORA Table 12 and a VMH Diet Designer diet
#
# Usage (from any folder):
#     python 03_make_micom_media.py
#
# Follows the structure of data/inputs/vmh_high_fiber_reference.ipynb, applied to five diets at
# once (stacked in one long table with a `diet` column, so every step runs once):
#
#   * custom_diet: your Diet Designer export (mmol/day; duplicated reactions summed)
#   * agora_s12_western, agora_s12_high_fiber: AGORA Supplementary Table 12 (mmol/gDW/h, from
#     01_extract_s12_table.py)
#   * agora_s12_western_bmt, agora_s12_high_fiber_bmt: the same two Table 12 diets, but absorbed
#     metabolites are multiplied by 0.8 instead of 0.1, representing a less efficient upper GI
#     tract that lets more nutrients reach the colon
#
# Steps: 1. floor very low values, 2. dilute metabolites absorbed in the small intestine
# (Recon3D exchanges), 3. add mucin cores, primary bile acids, cellulose and O2, 4. save each
# skeleton to media_out/<diet>_skeleton.csv and sanity-check it with 04_check_skeletons.py.
#
# The skeletons are then completed and validated one diet at a time with 05_complete_medium.py
# (complete_db_medium, check_db_medium); benchmarking/ (a SLURM job array, see its README) compares completion
# parameters.
#
# Needs cobra and data/Recon3D.xml.gz here; micom and data/agora201_refseq216_species_1.qza are
# needed by 05_complete_medium.py.
import subprocess
import sys
from pathlib import Path

import pandas as pd
from cobra.io import read_sbml_model

root = Path(__file__).resolve().parent.parent
dilutions = {"custom_diet": 0.2,
             "agora_s12_western": 0.1, "agora_s12_high_fiber": 0.1,
             "agora_s12_western_bmt": 0.8, "agora_s12_high_fiber_bmt": 0.8}

# --- Read the diets -------------------------------------------------------------------------------
# The template reads agora_metabolites.csv for annotations. We rebuild an equivalent table from
# the provided AGORA media, converted to CSV by 02_qza_to_csv.py (data/media_csv/), plus the names
# in Table 12. Metabolites missing from it are kept with an empty name.
media = [pd.read_csv(f) for f in sorted((root / "data/media_csv").glob("*.csv"))]
annotations = (pd.concat([m for m in media if "name" in m])
               .loc[:, ["metabolite", "name", "hmdb", "kegg.compound", "pubchem.compound", "inchi", "chebi"]]
               .dropna(subset=["name"]).drop_duplicates("metabolite"))
s12 = pd.read_csv(root / "data/agora_s12_table12.csv")
annotations = pd.concat([annotations, s12.loc[~s12.metabolite.isin(annotations.metabolite), ["metabolite", "name"]]])
print("annotations:", annotations.shape)

# custom diet: sum duplicated reactions, strip IDs, map Diet Designer names to AGORA names
custom = pd.read_csv(root / "data/inputs/custom_diet_fluxes.csv").iloc[:, 0:2]
custom.columns = ["reaction", "flux"]
custom["metabolite"] = custom.reaction.str.replace("^EX_", "", regex=True).str.replace("\\[e\\]|\\(e\\)", "", regex=True)
custom["metabolite"] = custom.metabolite.replace({"oleic": "ocdcea", "palmtl": "hdcea", "tedca": "ttdca",
                                                   "xol": "chsterol", "occa": "octa", "4hpro_LT": "4hpro"})
custom = custom.groupby("metabolite", as_index=False).flux.sum()

# Table 12: one column per diet
west = s12[["metabolite", "western"]].rename(columns={"western": "flux"})
fiber = s12[["metabolite", "high_fiber"]].rename(columns={"high_fiber": "flux"})

diet = pd.concat({"custom_diet": custom,
                  "agora_s12_western": west, "agora_s12_high_fiber": fiber,
                  "agora_s12_western_bmt": west, "agora_s12_high_fiber_bmt": fiber},
                 names=["diet"]).reset_index(level=0)
diet.loc[diet.flux == 0, "flux"] = 1e-4  # bug in VMH designer where everything <1e-4 gets truncated to 0

# --- Adjust for intestinal absorption -------------------------------------------------------------
# Dilute every metabolite that has an exchange in the Recon3D human model. The factor is 0.2 for
# the custom diet (as in the template), 0.1 for the Table 12 diets (the factor seen in
# western_diet_gut_agora.qza) and 0.8 for the two _bmt variants (edit `dilutions` above to change
# any of them). Note this reproduces the provided Western medium only roughly: Recon3D has no
# exchange for several ions and does have one for water.
recon3 = read_sbml_model(root / "data/Recon3D.xml.gz")
exchanges = pd.Series([r.id for r in recon3.exchanges])
exchanges = exchanges.str.replace("__", "_").str.replace("_e$|EX_", "", regex=True)

diet["dilution"] = 1.0
absorbed = diet.metabolite.isin(exchanges)
diet.loc[absorbed, "dilution"] = diet.loc[absorbed, "diet"].map(dilutions)
diet["flux"] = diet["flux"] * diet["dilution"]
print(diet.groupby("diet").dilution.value_counts())

# --- Adding host supplied components --------------------------------------------------------------
# Mucins, primary bile acids and a minuscule amount of oxygen, as in the template. Cellulose is set
# to 0.1 for the custom diet only (Table 12 already has its own cellulose value).
diet.set_index(["diet", "metabolite"], inplace=True)

for d in dilutions:
    # mucin
    for met in annotations.loc[annotations.metabolite.str.contains("core"), "metabolite"]:
        diet.loc[(d, met), "flux"] = 1

    # primary BAs
    for met in ["gchola", "tchola"]:
        diet.loc[(d, met), "flux"] = 1

    # anaerobic
    diet.loc[(d, "o2"), "flux"] = 0.001

# fiber
diet.loc[("custom_diet", "cellul"), "flux"] = 0.1

diet.reset_index(inplace=True)
diet["reaction"] = "EX_" + diet.metabolite + "(e)"

# Merge with the annotations to make the table more accessible.
skeleton = pd.merge(diet, annotations, on="metabolite", how="left")  # left: our annotation table is incomplete

skeleton["global_id"] = skeleton.reaction
skeleton["reaction"] = "EX_" + skeleton.metabolite + "_m"
print(skeleton.groupby("diet").size())

# --- Save the skeleton media (before completion) --------------------------------------------------
# So they can already be used with 07_compare_micom_media.py.
(root / "media_out").mkdir(exist_ok=True)
for d, skel in skeleton.groupby("diet"):
    skel.to_csv(root / f"media_out/{d}_skeleton.csv", index=False)

# --- Sanity-check the skeletons -------------------------------------------------------------------
# Before the slow completion step, run 04_check_skeletons.py. It fails (and stops this script) on
# malformed or duplicated IDs, legacy IDs that the models do not have (for example `glc` instead
# of `glc_D`), and Table 12 values that do not match the diet's dilution.
skeletons = sorted(str(p.resolve()) for p in (root / "media_out").glob("*_skeleton.csv"))
subprocess.run([sys.executable, str(root / "code/04_check_skeletons.py"), *skeletons,
                "--media-dir", str(root / "data/media_csv"),
                "--s12-csv", str(root / "data/agora_s12_table12.csv")],
               check=True)

# --- Complete the medium --------------------------------------------------------------------------
# Completion and validation run in 05_complete_medium.py, one skeleton at a time, so parameters
# are easy to vary:
#
#     python 05_complete_medium.py ../media_out/agora_s12_high_fiber_bmt_skeleton.csv \
#         --weights mass --max-added-import 10 --threads 12
#
# To compare several `weights` and `max_added_import` settings use
#
#     cd benchmarking; sbatch benchmarking.sh
