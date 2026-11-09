#!/usr/bin/env bash
# Regenerates everything from the files in data/inputs.  Run from this folder.
set -euo pipefail
# Run notebooks/make_micom_media.ipynb first (needs micom, qiime2, cobra, data/agora201_refseq216_species_1.qza);
# the skeleton CSVs it writes to media_out/ are used below.
cd code
python extract_s12_table.py ../data/inputs/AGORA_S12_diets_nbt3703_MOESM2.pdf ../data/agora_s12_table12.csv
python compare_micom_media.py ../data/inputs -o ../results/comparison_all_diets \
    --fluxes-csv ../data/inputs/custom_diet_fluxes.csv:custom_diet --s12-csv ../data/agora_s12_table12.csv \
    --medium-csv ../media_out/agora_s12_western_bmt_skeleton.csv:agora_s12_western_bmt \
    --medium-csv ../media_out/agora_s12_high_fiber_bmt_skeleton.csv:agora_s12_high_fiber_bmt
python compare_s12_vs_micom.py --s12-csv ../data/agora_s12_table12.csv --media-dir ../data/inputs \
    -o ../results/s12_vs_micom_provided
python estimate_energy.py ../data/inputs/custom_diet_fluxes.csv ../results/comparison_all_diets/all_media_long.csv \
    > ../results/energy/energy_estimates.txt
# The notebook (notebooks/make_micom_media.ipynb) needs micom, qiime2, cobra and data/agora201_refseq216_species_1.qza; run it separately in Jupyter.
