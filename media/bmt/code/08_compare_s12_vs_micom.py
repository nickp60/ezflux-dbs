#!/usr/bin/env python3
"""
Compare AGORA Supplementary Table 12 (Western / High fiber, mmol/gDW/h) with the
MICOM-formatted media: how was western_diet_gut_agora derived from the table, and
how does the table's "High fiber" design differ from vmh_high_fiber_agora?

Usage
-----
    python 08_compare_s12_vs_micom.py --s12-csv ../data/agora_s12_table12.csv \
        --media-dir ../data/media_csv -o ../results/s12_vs_micom_provided

Writes  s12_vs_provided_media.csv  and  s12_vs_provided_media.png, and prints:
  * the ratio  western_diet_gut_agora / Table-12 Western  (0.1, 1 or 0)
  * Table-12 High fiber / Western ratios per category
  * log2 fold changes in the table's pair vs the vmh_high_fiber / vmh_eu_average pair
"""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

import importlib.util
spec = importlib.util.spec_from_file_location("cm", Path(__file__).with_name("07_compare_micom_media.py"))
cm = importlib.util.module_from_spec(spec); spec.loader.exec_module(cm)

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--s12-csv", required=True)
ap.add_argument("--media-dir", required=True)
ap.add_argument("-o", "--outdir", default="s12_vs_micom")
a = ap.parse_args()
out = Path(a.outdir); out.mkdir(parents=True, exist_ok=True)
sns.set_theme(style="whitegrid")

T = pd.read_csv(a.s12_csv).rename(columns={"western": "west_pdf", "high_fiber": "hf_pdf", "table_category": "table_cat"})
T["cls"] = T.metabolite.map(cm.classify)
L = cm.load_all([a.media_dir])
w = L.pivot(index="metabolite", columns="diet", values="flux")
T = T.set_index("metabolite")
for c in ["western_diet_gut_agora", "western_diet_gut_carveme", "vmh_high_fiber_agora", "vmh_eu_average_agora"]:
    T[c] = w[c].reindex(T.index)
print(f"Table rows: {len(T)}; present in western_diet_gut_agora: {T.western_diet_gut_agora.notna().sum()}")
print("In western_diet_gut_agora but not in table:", sorted(set(w.index[w.western_diet_gut_agora.notna()]) - set(T.index)))

T["ratio_w"] = T.western_diet_gut_agora / T.west_pdf
r = T.dropna(subset=["ratio_w"])
is01 = (r.ratio_w - 0.1).abs() < 0.002
is1 = (r.ratio_w - 1).abs() < 0.02
print(f"\nwestern qza / Table 12:  x0.1 -> {is01.sum()},  x1 -> {is1.sum()},  other -> {(~(is01 | is1)).sum()}")
print("x1 group:", r.index[is1].tolist())
T["hf_over_w_pdf"] = T.hf_pdf / T.west_pdf
print("\nTable 12 High fiber / Western by category:")
print(T.groupby("table_cat").hf_over_w_pdf.agg(["min", "median", "max"]).round(3))

T["fc_pdf"] = np.log2(T.hf_pdf / T.west_pdf)
T["fc_vmh"] = np.log2(T.vmh_high_fiber_agora / T.vmh_eu_average_agora)
sh = T.dropna(subset=["fc_vmh"]).copy()
sh = sh[sh.cls != "Water & gases"]
print(f"\nShared with the vmh pair: {len(sh)}; sign agreement {(np.sign(sh.fc_pdf) == np.sign(sh.fc_vmh)).mean():.2f}; "
      f"Spearman {sh[['fc_pdf', 'fc_vmh']].corr(method='spearman').iloc[0, 1]:.2f}")
g = sh.groupby("cls").agg(n=("fc_vmh", "size"), pdf=("fc_pdf", "median"), vmh=("fc_vmh", "median"))
g["S12_x"] = 2 ** g.pdf; g["vmh_x"] = 2 ** g.vmh
print(g.round(2).sort_values("n", ascending=False).to_string())

fig, ax = plt.subplots(1, 2, figsize=(15, 6.3))
d = T.dropna(subset=["ratio_w"]).copy()
d["grp"] = np.where((d.ratio_w - 0.1).abs() < 0.002, "x0.1 (absorbed in small intestine)",
                    np.where((d.ratio_w - 1).abs() < 0.02, "x1 (starch, fibre, water, methanol)", "set to 0"))
pal = {"x0.1 (absorbed in small intestine)": "#4c78a8", "x1 (starch, fibre, water, methanol)": "#f58518", "set to 0": "#e45756"}
for k, c in pal.items():
    s = d[d.grp == k]
    ax[0].scatter(s.west_pdf, s.western_diet_gut_agora.clip(lower=1e-9), s=28, color=c, label=f"{k} (n={len(s)})",
                  alpha=.85, edgecolor="white", lw=.4)
lim = [3e-9, 30]
ax[0].plot(lim, lim, "k--", lw=1, label="y = x"); ax[0].plot(lim, [l * .1 for l in lim], "k:", lw=1, label="y = 0.1 x")
ax[0].set_xscale("log"); ax[0].set_yscale("log"); ax[0].set_xlim(lim); ax[0].set_ylim(1e-9, 30)
ax[0].set_xlabel("Table 12 Western uptake rate (mmol/gDW/h)"); ax[0].set_ylabel("western_diet_gut_agora flux (zeros plotted at 1e-9)")
ax[0].set_title("Western: Table 12 vs the provided medium"); ax[0].legend(fontsize=8, loc="upper left")
order = ["Simple sugars & sugar alcohols", "Complex carbs & fibre", "Fatty acids & lipids", "Amino acids", "Minerals & ions", "Vitamins & cofactors"]
long = pd.concat([sh.assign(src="AGORA Table 12: High fiber / Western", fc=sh.fc_pdf),
                  sh.assign(src="Provided media: vmh_high_fiber / vmh_eu_average", fc=sh.fc_vmh)])
long = long[long.cls.isin(order)]
sns.boxplot(data=long, y="cls", x="fc", hue="src", order=order, ax=ax[1], palette=["#54a24b", "#b279a2"], fliersize=0, width=.7)
sns.stripplot(data=long, y="cls", x="fc", hue="src", order=order, ax=ax[1], dodge=True, palette=["#2d6a2a", "#7a4a76"], size=2.5, alpha=.6, legend=False)
ax[1].axvline(0, color="k", lw=1); ax[1].set_xlabel("log2 fold-change (high fibre vs baseline), shared metabolites"); ax[1].set_ylabel("")
ax[1].set_title('Two different "high fibre" designs'); ax[1].legend(fontsize=8, loc="upper center", bbox_to_anchor=(0.5, -0.13), ncol=1, frameon=False)
fig.tight_layout(); fig.savefig(out / "s12_vs_provided_media.png", dpi=200)
T.reset_index().to_csv(out / "s12_vs_provided_media.csv", index=False)
print("\nWrote", out.resolve())
