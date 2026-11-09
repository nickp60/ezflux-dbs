# MICOM diet media: comparison, AGORA Table 12 reconciliation, and build notebook

Analysis for the Peled Lab of the micom-dev flux-modelling diets (`micom-dev/media`), a custom VMH Diet Designer diet, and the two diets in AGORA Supplementary Table 12 (Western, High fiber).

## 1. Files assessed

| Group | Files (in `data/inputs/`) | Notes |
|---|---|---|
| MICOM media (QIIME 2 `MicomMedium[Global]`) | `western_diet_gut_agora.qza`, `western_diet_gut_carveme.qza`, `vmh_eu_average_agora.qza`, `vmh_high_fiber_agora.qza`, `vmh_high_fat_low_carb_agora.qza`, `vmh_fermented_agora.qza`, `vmh_formula_agora.qza`, `breast_milk_agora.qza`, `himalaya.qza`, `guerrero_mountains.qza`, `baobab_honey_antelope.qza` | A `.qza` is a zip; the table is `<uuid>/data/medium.csv`. |
| Custom diet | `custom_diet_fluxes.csv` | VMH Diet Designer export, mmol/(human·day); 443 rows, 328 repeated reactions, **summed**; 106 non-zero metabolites after merging. |
| AGORA Table 12 | `AGORA_S12_diets_nbt3703_MOESM2.pdf` → `data/agora_s12_table12.csv` | Magnúsdóttir et al. 2017; 164 rows, mmol/gDW/h. |
| Reference recipe | `vmh_high_fiber_reference.ipynb` | How `vmh_high_fiber_agora.qza` was built. |

## 2. Conclusions

**Scale and units are the main trap.** The `.qza` files carry no units, and the recipes differ.

- AGORA Table 12 and the `western_diet_*` media are per-biomass bounds, mmol/gDW/h.
- The custom diet is mmol/(human·day). The `vmh_*` recipe reads a pre-made Diet Designer table, applies a 0.2 dilution to Recon3D metabolites (a multiplication of the flux), then adds mucin cores, primary bile acids, cellulose and O2, and completes the medium. No /24 is visible in the notebook.
- An energy check points both ways. The custom diet works out to about 2,070–2,275 kcal/day by my estimate (your 2,353 kcal is close; the gap is mostly D-amino-acid duplicates and polymer-size assumptions). The *undiluted* starch in `vmh_eu_average` is about 1/20 of the custom diet's, as if /24 had been applied upstream, but the non-starch compounds in the same medium are far too large for that. **The scale of the `vmh_*` media is therefore unverified.** An earlier statement that all `vmh_*` media are mmol/person/h was wrong and is withdrawn.
- Consequently the redone comparison uses unit-free views as the primary comparison (within-diet composition, presence/absence, Spearman rank correlation) and shows absolute fluxes only within a unit group.

**Why the overall flux looks low.** Molar sums mislead: starch exchanges are polymers (11–12 glucose residues), the diets contain large water, ion and methanol terms, and the MICOM recipe deliberately dilutes small-intestine-absorbed compounds (x0.1 or x0.2) because only the unabsorbed remainder reaches the colon. Raw totals are therefore not an intake.

**Comparison of all 14 diets** (`results/comparison_all_diets/`):

- The two Table 12 diets and the two `western_diet_gut_*` media form one tight cluster (Spearman ρ 0.88–0.99 with each other) and are anti-correlated with every other medium (ρ ≈ −0.2 to −0.55). They are flat designs: most values are one of a few constants.
- The `vmh_*` media cluster together (`vmh_eu_average`, `vmh_fermented`, `vmh_high_fiber` ρ 0.94–0.98).
- The custom diet is closest by rank to `vmh_eu_average` (ρ 0.91), then `vmh_high_fiber` (0.82). It is 47% starch by energy, has about 26 g/day ethanol, is protein-poor (about 7% of energy) and carries D-amino-acid duplicates.
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

```
README.md
requirements.txt        python packages (+ optional micom, qiime2, cobra)
run_all.sh              regenerate everything below from data/inputs
code/
  extract_s12_table.py      PDF -> CSV
  compare_micom_media.py    parse + compare all diets
  compare_s12_vs_micom.py   Table 12 vs the provided media
  estimate_energy.py        rough kcal/day estimate
data/
  agora_s12_table12.csv     parsed Table 12 (164 rows)
  agora_metabolites_from_qza.csv   annotation table rebuilt by the notebook
  Recon3D.xml.gz            human model, used for the absorbed-metabolite list
  inputs/                   all files assessed (section 1)
notebooks/
  make_micom_media.ipynb    build MICOM-format media (executed)
media_out/                  written by the notebook (empty until it is run with micom/qiime2)
results/
  comparison_all_diets/     plots 01-06 + CSV tables for all 14 diets
  s12_vs_micom_provided/    Table 12 reconciliation (PNG + CSV)
  energy/energy_estimates.txt
```

All commands are run from `code/` (the scripts import each other).

**`extract_s12_table.py`** — `python extract_s12_table.py <pdf> <out.csv>`. Needs `pdftotext`. Writes `table_category, pdf_id, reaction, metabolite, name, western, high_fiber`, with harmonised IDs. Warns if the row count is not 164.

**`compare_micom_media.py`** — `python compare_micom_media.py <dir-or-qza...> -o <outdir> [--fluxes-csv CSV[:NAME]] [--s12-csv CSV] [--medium-csv CSV[:NAME]]`. `--fluxes-csv` adds a Diet Designer diet (duplicates summed, zero rows dropped, aliases such as oleic→ocdcea applied; can repeat). `--s12-csv` adds `agora_s12_western` and `agora_s12_high_fiber`. `--medium-csv CSV[:NAME]` adds a MICOM-format medium written by the notebook (repeatable); `results/comparison_all_diets/` includes the two `_bmt` diets this way, from their skeleton CSVs (before `complete_db_medium`, so they differ from the plain Table 12 diets by the added mucin cores, bile acids and O2 as well as the 0.8 dilution). Outputs: `all_media_long.csv`, `all_media_wide.csv`, `summary_per_diet.csv`, `composition_fraction_{all,organic}.csv`, `spearman_rho.csv`, `shared_metabolite_counts.csv` and plots 01–06 (overview and overlap, composition, absolute fluxes by class and unit group, top-metabolite heatmap, rank correlation, key nutrient classes). Metabolite classes are set by hand in `CLASS_SETS`; extend them if the "Other" count is not 0. The old `--vmh-biomass-gdw` option was removed because the `vmh_*` units are unverified.

**`compare_s12_vs_micom.py`** — `python compare_s12_vs_micom.py --s12-csv ../data/agora_s12_table12.csv --media-dir ../data/inputs -o <outdir>`. Prints the 0.1 / 1 / 0 rule and the High-fiber/Western ratios by category, and saves the comparison figure and a CSV.

**`estimate_energy.py`** — `python estimate_energy.py <custom_fluxes.csv> <all_media_long.csv>`. Atwater-factor kcal estimate for a Diet Designer file; for the MICOM media it prints stored-value energy only as a scale check (see its docstring).

**`make_micom_media.ipynb`** — run with Jupyter from `notebooks/`; needs `micom`, `qiime2`, `cobra`, `data/Recon3D.xml.gz` (included) and `data/agora201_refseq216_species_1.qza` (not included; place it in `data/`). It follows the `vmh_high_fiber` template step by step, with five diets (custom, Table 12 Western, Table 12 High fiber, and the two "_bmt" variants below) stacked in one table so each step runs once:

1. Sum duplicated reactions in the custom diet, map Diet Designer names to AGORA names (oleic→ocdcea, etc.), floor zeros to 1e-4.
2. Dilute every metabolite with a Recon3D exchange: 0.2 for the custom diet (as in the template), 0.1 for the Table 12 diets (the factor seen in `western_diet_gut_agora.qza`), and 0.8 for `agora_s12_western_bmt` and `agora_s12_high_fiber_bmt`, which are the same Table 12 diets with a less efficient upper GI tract (more nutrients reach the colon). All are set in the `dilutions` dict in the first code cell.
3. Add mucin cores, gchola/tchola and O2 0.001 to all five, and cellulose 0.1 to the custom diet only.
4. Merge annotations, `complete_db_medium` (growth 0.05, max_added_import 100, weights="mass") per diet, `check_db_medium`, then save `media_out/<diet>.qza` and a CSV.

The annotation table is rebuilt from the AGORA-format `.qza` files and Table 12, since the template's `agora_metabolites.csv` is not available, so the merge is a left join (unannotated compounds keep an empty name).

**What was run here.** The cells up to and including the skeleton table were executed (outputs are saved in the notebook). The `complete_db_medium`, `check_db_medium` and save cells were not run because `micom`, `qiime2` and the model database were not available, so `media_out/` holds only the skeleton CSVs (`<diet>_skeleton.csv`, written by a cell before completion) until you run them. After running the notebook, re-run the comparison with `--medium-csv media_out/<diet>_micom_medium.csv:<diet>` to use the completed media. Using the Recon3D exchange list reproduces the 0.1 / x1 split of the provided Western medium for only about 80% of Table 12 compounds (for example Recon3D has no exchange for several ions, and has one for water, which the provided medium leaves undiluted); the provided medium is what the MICOM authors made, so treat the Table 12 diets here as a close but not exact rebuild.

## 5. Open items

- Units and scaling of the `vmh_*` media, and of `himalaya`, `guerrero_mountains` and `baobab_honey_antelope`.
- Why 10 Table 12 compounds are zero in `western_diet_gut_agora`, and what the 9 extra compounds are (completion step assumed).
- Whether the custom diet should be divided by 24 (mmol/day to mmol/h) before it is used as a medium. The notebook follows the reference notebook and does not; add `custom["flux"] /= 24` after the groupby if you want it on the same footing as the `vmh_*` starch values.
