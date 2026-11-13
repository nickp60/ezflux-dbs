#!/usr/bin/env bash
# Runs the scripts in order (01 .. 06).  Run from this folder:  bash run_all.sh
#
# Step 03 is a Jupyter notebook that needs micom, cobra, data/Recon3D.xml.gz and
# data/agora201_refseq216_species_1.qza, so it is NOT run here.  Run
# notebooks/03_make_micom_media.ipynb yourself between steps 02 and 04; it writes
# media_out/<diet>_skeleton.csv (before completion) and, when complete,
# media_out/<diet>_micom_medium.csv.  Step 04 uses the completed files if they exist
# and falls back to the skeleton files otherwise.
set -euo pipefail
mkdir -p  results/energy/
cd code

python 01_extract_s12_table.py ../data/inputs/AGORA_S12_diets_nbt3703_MOESM2.pdf ../data/agora_s12_table12.csv
python 02_qza_to_csv.py ../data/media/media/ ../data/media_csv

python 07_check_skeletons.py ../media_out/*_skeleton.csv --media-dir ../data/media_csv --s12-csv ../data/agora_s12_table12.csv

suffix=micom_medium
ls ../media_out/*_${suffix}.csv > /dev/null 2>&1 || suffix=skeleton
echo "Using media_out/*_${suffix}.csv"

python 04_compare_micom_media.py ../data/media_csv ../media_out/*_${suffix}.csv -o ../results/comparison_all_diets
python 05_compare_s12_vs_micom.py --s12-csv ../data/agora_s12_table12.csv --media-dir ../data/media_csv \
    -o ../results/s12_vs_micom_provided
python 06_estimate_energy.py ../data/inputs/custom_diet_fluxes.csv ../results/comparison_all_diets/all_media_long.csv \
    > ../results/energy/energy_estimates.txt
