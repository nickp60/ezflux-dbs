#!/usr/bin/env bash
# Runs the scripts in execution order (01 .. 09).  Run from this folder:  bash run_all.sh
#
#   01 extract Table 12 from the PDF          06 benchmark completion parameters (NOT run here,
#   02 provided .qza media -> CSV                 it is slow; run it yourself, see the README)
#   03 build skeletons (also runs 04)         07 compare the completed media
#   04 check skeletons                        08 Table 12 vs the provided media
#   05 complete + validate each skeleton      09 energy estimate
#
# Step 05 needs micom and data/agora201_refseq216_species_1.qza.  If the .qza is missing it is
# skipped, and step 07 then stops because it requires completed media (--require-completed).
set -euo pipefail
mkdir -p  results/energy/
cd code

python 01_extract_s12_table.py ../data/inputs/AGORA_S12_diets_nbt3703_MOESM2.pdf ../data/agora_s12_table12.csv
python 02_qza_to_csv.py ../data/media/media/ ../data/media_csv
python 03_make_micom_media.py

if [ -f ../data/agora201_refseq216_species_1.qza ]; then
    for skeleton in ../media_out/*_skeleton.csv; do
        python 05_complete_medium.py "$skeleton"
    done
else
    echo "data/agora201_refseq216_species_1.qza not found: skipping step 05"
fi

python 07_compare_micom_media.py ../data/media_csv ../media_out/*_micom_medium.csv --require-completed -o ../results/comparison_all_diets
python 08_compare_s12_vs_micom.py --s12-csv ../data/agora_s12_table12.csv --media-dir ../data/media_csv \
    -o ../results/s12_vs_micom_provided
python 09_estimate_energy.py ../data/inputs/custom_diet_fluxes.csv ../results/comparison_all_diets/all_media_long.csv \
    > ../results/energy/energy_estimates.txt
