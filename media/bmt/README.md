# MICOM diet media: comparison, AGORA Table 12 reconciliation, and build notebook

Analysis for the Peled Lab of the micom-dev flux-modelling diets (`micom-dev/media`), a custom VMH Diet Designer diet, and the two diets in AGORA Supplementary Table 12 (Western, High fiber), plus two variants of the Table 12 diets with a less efficient upper GI tract (`_bmt`).

## 1. Files assessed

| Group | Files (in `data/inputs/`) | Notes |
|---|---|---|
| MICOM media (QIIME 2 `MicomMedium[Global]`) | `western_diet_gut_agora.qza`, `western_diet_gut_carveme.qza`, `vmh_eu_average_agora.qza`, `vmh_high_fiber_agora.qza`, `vmh_high_fat_low_carb_agora.qza`, `vmh_fermented_agora.qza`, `vmh_formula_agora.qza`, `breast_milk_agora.qza`, `himalaya.qza`, `guerrero_mountains.qza`, `baobab_honey_antelope.qza` | A `.qza` is a zip; the table is `<uuid>/data/medium.csv`. Step 02 extracts them to `data/media_csv/*.csv`; nothing after that reads a `.qza` except the model database used by the notebook. |
| Custom diet | `custom_diet_fluxes.csv` | VMH Diet Designer export, mmol/(human·day); 443 rows, 328 repeated reactions, **summed**; 106 non-zero metabolites. |
| AGORA Table 12 | `AGORA_S12_diets_nbt3703_MOESM2.pdf` → `data/agora_s12_table12.csv` | Magnúsdóttir et al. 2017; 164 rows, mmol/gDW/h. |
| Reference recipe | `vmh_high_fiber_reference.ipynb` | How `vmh_high_fiber_agora.qza` was built. |

## 2. Conclusions

**Scale and units are the main trap.** The media files carry no units, and the recipes differ.

- AGORA Table 12 and the `western_diet_*` media are per-biomass bounds, mmol/gDW/h.
- The custom diet is mmol/(human·day). The `vmh_*` recipe reads a pre-made Diet Designer table, applies a 0.2 dilution to Recon3D metabolites (a multiplication of the flux), then adds mucin cores, primary bile acids, cellulose and O2, and completes the medium. No /24 is visible in the notebook.
- An energy check points both ways. The custom diet works out to about 2,070–2,275 kcal/day by my estimate (your 2,353 kcal is close; the gap is mostly D-amino-acid duplicates and polymer-size assumptions). The *undiluted* starch in `vmh_eu_average` is about 1/20 of the custom diet's, as if /24 had been applied upstream, but the non-starch compounds in the same medium are far too large for that. **The scale of the `vmh_*` media is therefore unverified.** An earlier statement that all `vmh_*` media are mmol/person/h was wrong and is withdrawn.
- Consequently the comparison uses unit-free views as the primary comparison (within-diet composition, presence/absence, Spearman rank correlation) and shows absolute fluxes only within a unit group.

**Why the overall flux looks low.** Molar sums mislead: starch exchanges are polymers (11–12 glucose residues), the diets contain large water, ion and methanol terms, and the MICOM recipe deliberately dilutes small-intestine-absorbed compounds (x0.1 or x0.2) because only the unabsorbed remainder reaches the colon. Raw totals are therefore not an intake.

**Comparison of 16 diets** (`results/comparison_all_diets/`): the 11 provided media, plus the five diets built by the notebook (custom diet, Table 12 Western and High fiber, and their two `_bmt` variants). The results currently in the folder use the notebook's **skeleton** media (diluted and with host additions, but not yet completed by `complete_db_medium`, which could not be run here). The previous results are kept in `results/comparison_all_diets.bak/`; they used the raw custom export and the tabulated Table 12 values.

- The Table 12 diets (all four versions) and the two `western_diet_gut_*` media form one cluster (Spearman ρ 0.77–1.00 with each other) and are uncorrelated or anti-correlated with every other medium (ρ about −0.4 to 0.0). They are flat designs: most values are one of a few constants. The `_bmt` variants differ from the 0.1 versions only by the dilution factor, so their ranks are almost identical (ρ 0.98).
- The `vmh_*` media cluster together (`vmh_eu_average`, `vmh_fermented`, `vmh_high_fiber`).
- The custom diet (as built by the notebook) is closest by rank to `vmh_eu_average` (ρ 0.91), then `vmh_high_fiber` (0.78) and `vmh_fermented` (0.72). It is 47% starch by energy, has about 26 g/day ethanol, is protein-poor (about 7% of energy) and carries D-amino-acid duplicates.
- Two exceptions to the clustering: `breast_milk_agora` and `vmh_formula_agora` stand apart (lactose-dominated), and `himalaya`, `guerrero_mountains` and `baobab_honey_antelope` (units unverified) resemble the VMH group in rank.

**How `western_diet_gut_agora` relates to Table 12** (`results/s12_vs_micom_provided/`). Of the 164 Table 12 compounds, 123 are multiplied by 0.1, 29 are kept at x1 (starch, all polysaccharide fibres, cellobiose, trehalose, raffinose, water, methanol), 10 are zero in the medium although non-zero in the table (glucose, arachidonate, palmitoleate, linoleate, α- and γ-linolenate, octanoate, Fe(III) dicitrate, D-glutamate, pyridoxal phosphate), and 2 are missing (adocbl, nitrate). Nine compounds appear that are not in the table (`adpcbl`, `amet`, `cytd`, `glyasp`, `glyleu`, `gsn`, `o2`, `sbt_D`, `succ`), probably from the medium-completion step. The reason for the 10 zeros is not known to me.

## 3. MICOM `vmh_high_fiber` versus AGORA Table 12 High fiber

They are two different designs that share a name.

| | AGORA Table 12 High fiber | `vmh_high_fiber_agora` (micom-dev) |
|---|---|---|
| Construction | The Western diet rescaled by compound class | VMH Diet Designer "high fiber" diet, processed with the recipe above |
| Fibres / polysaccharides | **x22.15** (all 24 fibres, including starch) | barely changed relative to `vmh_eu_average` (median x1.4, 3 shared compounds) |
| Sugars and sugar alcohols | **x0.265** | unchanged (median x1.0, 11 shared) |
| Fats | **x0.5** | x0.37 (median) |
| Amino acids, vitamins, minerals | exactly x1 | higher: amino acids x1.35, vitamins x1.16, minerals/ions x1.57 |
| Baseline | Western (Table 12) | `vmh_eu_average` (a different baseline) |
| Direction agreement on shared compounds | n = 98 compounds shared; signs agree for only 35%; Spearman of log2 fold-changes 0.15 | |

In short: Table 12 shifts energy *from sugars toward fibre* and leaves everything else fixed, whereas the MICOM `vmh_high_fiber` medium is a higher-protein, higher-mineral, polyol-rich diet in which fibre is not much higher than in the EU-average diet. The completion step also adds roughly 11.5 mmol/h of extra components in the high-fibre recipe.

**Recommendation for simulating high fibre from the Western baseline.** Raise the complex-carbohydrate and fibre terms (starch, cellulose, inulin, arabinoxylan, pectin and similar), as Table 12 does, rather than copying the VMH pattern of more amino acids, sorbitol, mannitol and lactose. The latter carries a protein/polyol signal that is mostly absorbed in the small intestine (polyols are only about a third absorbed; Lenhart & Chey 2017), and more colonic protein shifts fermentation toward proteolytic products (Diether & Willing 2019). Fibre depletion is what drives mucus foraging in Western-diet models (Desai et al. 2016), and the effect of added fibre is the substrate itself (Louis 2021; Dobranowski & Stintzi 2021). If you want the VMH-style design, build it from the VMH Diet Designer as in the notebook, not by hand-editing Western values. These citations were collected earlier in the session from memory and search results; please verify details before citing them in a paper.

## 4. Archive layout and usage guide

The scripts are numbered in execution order.

```
README.md
requirements.txt
run_all.sh              runs 01, 02, 04, 05, 06 (step 03 is a notebook, see below)
code/
  01_extract_s12_table.py      PDF -> data/agora_s12_table12.csv
  02_qza_to_csv.py             provided .qza media -> data/media_csv/*.csv
notebooks/
  03_make_micom_media.ipynb    build the five MICOM-format diets (writes CSV only)
code/
  04_compare_micom_media.py    parse + compare all diets (plots + tables)
  05_compare_s12_vs_micom.py   Table 12 vs the provided media
  06_estimate_energy.py        rough kcal/day estimate
data/
  agora_s12_table12.csv        parsed Table 12 (164 rows)
  agora_metabolites_from_qza.csv   annotation table rebuilt by the notebook
  Recon3D.xml.gz               human model, used for the absorbed-metabolite list
  media_csv/                   the 11 provided media as CSV (output of step 02)
  inputs/                      all files assessed (section 1)
media_out/                     written by the notebook: <diet>_skeleton.csv now, <diet>_micom_medium.csv when completed
results/
  comparison_all_diets/        plots 01-06 + CSV tables for all 16 diets
  s12_vs_micom_provided/       Table 12 reconciliation (PNG + CSV)
  energy/energy_estimates.txt
```

Files ending in `.bak` are earlier versions that were moved aside, not deleted (`run_all.sh.bak`, `requirements.txt.bak`, `README.md.bak`, `notebooks/make_micom_media.ipynb.bak`, `results/comparison_all_diets.bak/`). They can be removed once you have checked the new ones.

Run the scripts from `code/` (or run `bash run_all.sh` from the archive root).

**`01_extract_s12_table.py`** — `python 01_extract_s12_table.py <pdf> <out.csv>`. Needs `pdftotext`. Writes `table_category, pdf_id, reaction, metabolite, name, western, high_fiber`, with harmonised IDs. Warns if the row count is not 164.

**`02_qza_to_csv.py`** — `python 02_qza_to_csv.py ../data/inputs ../data/media_csv`. Extracts `data/medium.csv` from every `.qza` (no QIIME 2 needed) and strips upload-hash prefixes from the names.

**`03_make_micom_media.ipynb`** — run with Jupyter from `notebooks/`; needs `micom`, `cobra`, `data/Recon3D.xml.gz` (included) and `data/agora201_refseq216_species_1.qza` (not included; place it in `data/`). It follows the `vmh_high_fiber` template step by step, with five diets stacked in one table so each step runs once: `custom_diet`, `agora_s12_western`, `agora_s12_high_fiber`, and the two `_bmt` variants.

1. Sum duplicated reactions in the custom diet, map Diet Designer names to AGORA names (oleic→ocdcea, etc.), floor zeros to 1e-4.
2. Dilute every metabolite with a Recon3D exchange: 0.2 for the custom diet (as in the template), 0.1 for the Table 12 diets (the factor seen in `western_diet_gut_agora`), and 0.8 for `agora_s12_western_bmt` and `agora_s12_high_fiber_bmt`, the same Table 12 diets with a less efficient upper GI tract so that more nutrients reach the colon. All are set in the `dilutions` dict in the first code cell.
3. Add mucin cores, gchola/tchola and O2 0.001 to all five, and cellulose 0.1 to the custom diet only.
4. Merge annotations, save `media_out/<diet>_skeleton.csv`, then `complete_db_medium` (growth 0.05, max_added_import 100, weights="mass") per diet, `check_db_medium`, and save `media_out/<diet>_micom_medium.csv`.

The notebook writes CSV only; QIIME 2 artifacts are not created. The annotation table is rebuilt from the provided media (`data/media_csv/`) and Table 12, since the template's `agora_metabolites.csv` is not available, so the merge is a left join (unannotated compounds keep an empty name).

**`04_compare_micom_media.py`** — `python 04_compare_micom_media.py ../data/media_csv ../media_out/*_micom_medium.csv -o ../results/comparison_all_diets [--fluxes-csv CSV[:NAME]] [--s12-csv CSV]`. Inputs are medium CSV files and/or folders of them; the diet name is the file name without `_micom_medium` or `_skeleton`. The two options add diets directly from raw sources: `--fluxes-csv` a Diet Designer export (duplicates summed, zero rows dropped, aliases applied) and `--s12-csv` the tabulated Table 12 diets (named `agora_s12_western` and `agora_s12_high_fiber`, so do not use it together with the notebook output of the same names). Outputs: `all_media_long.csv`, `all_media_wide.csv`, `summary_per_diet.csv`, `composition_fraction_{all,organic}.csv`, `spearman_rho.csv`, `shared_metabolite_counts.csv` and plots 01–06 (overview and overlap, composition, absolute fluxes by class and unit group, top-metabolite heatmap, rank correlation, key nutrient classes). Metabolite classes are set by hand in `CLASS_SETS`; extend them if the "Other" count is not 0.

**`05_compare_s12_vs_micom.py`** — `python 05_compare_s12_vs_micom.py --s12-csv ../data/agora_s12_table12.csv --media-dir ../data/media_csv -o <outdir>`. Prints the 0.1 / 1 / 0 rule and the High-fiber/Western ratios by category, and saves the comparison figure and a CSV. It imports `04_compare_micom_media.py`, so keep the two together.

**`06_estimate_energy.py`** — `python 06_estimate_energy.py <custom_fluxes.csv> <all_media_long.csv>`. Atwater-factor kcal estimate for a Diet Designer file; for the MICOM media it prints stored-value energy only as a scale check (see its docstring).

**What was run here.** The notebook cells up to and including the skeleton table and the skeleton CSVs were executed (outputs are saved in the notebook). The `complete_db_medium`, `check_db_medium` and final save cells were not run because `micom` and the model database were not available here, so the comparison in `results/` uses the skeleton CSVs. After running the notebook, run step 04 again (or `run_all.sh`, which prefers the completed files) to use the completed media. Using the Recon3D exchange list reproduces the 0.1 / x1 split of the provided Western medium for only about 80% of Table 12 compounds (for example Recon3D has no exchange for several ions, and has one for water, which the provided medium leaves undiluted); the provided medium is what the MICOM authors made, so treat the Table 12 diets here as a close but not exact rebuild.

## 5. Open items

- Units and scaling of the `vmh_*` media, and of `himalaya`, `guerrero_mountains` and `baobab_honey_antelope`.
- Why 10 Table 12 compounds are zero in `western_diet_gut_agora`, and what the 9 extra compounds are (completion step assumed).
- Whether the custom diet should be divided by 24 (mmol/day to mmol/h) before it is used as a medium. The notebook follows the reference notebook and does not; add `custom["flux"] /= 24` after the groupby if you want it on the same footing as the `vmh_*` starch values.
