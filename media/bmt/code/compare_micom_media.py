#!/usr/bin/env python3
"""
Parse and compare MICOM growth media (QIIME 2 `MicomMedium[Global]` .qza files,
e.g. from https://github.com/micom-dev/media) and plot their nutrient content.

Usage
-----
    python compare_micom_media.py  ../data/inputs  -o ../results/comparison_all_diets \
        --fluxes-csv ../data/inputs/custom_diet_fluxes.csv:custom_diet \
        --s12-csv ../data/agora_s12_table12.csv

Inputs may be .qza files and/or directories of them.  Extra diets:
  --fluxes-csv  a VMH Diet Designer export (reaction, flux in mmol/(human*day));
                duplicated reactions are summed.  May be repeated.
  --s12-csv     the CSV from extract_s12_table.py; adds two diets,
                agora_s12_western and agora_s12_high_fiber (mmol/gDW/h).

  --medium-csv  a MICOM-format medium CSV written by make_micom_media.ipynb
                (media_out/<diet>_skeleton.csv or <diet>_micom_medium.csv); repeatable.

A .qza is just a zip archive; the medium lives in  <uuid>/data/medium.csv.
Files are read straight out of the zip (nothing is extracted to disk).

Units (read this before interpreting absolute numbers)
-----------------------------------------------------
The files do not record units, and the micom-dev recipes differ in how they
scale.  What is known:
  * AGORA Table 12 diets (and western_diet_*, which is Table 12 x0.1 for
    small-intestine-absorbed compounds): mmol/gDW/h.
  * VMH Diet Designer exports (the custom diet): mmol/(human*day), kept native here.
  * vmh_* media: built from VMH Diet Designer output (see
    vmh_high_fiber_reference.ipynb): 0.2 dilution of compounds in Recon3D, plus
    mucin/bile-acid/cellulose/O2 additions and a completion step.  The notebook
    reads a pre-made .tsv, so any /24 conversion would happen upstream and is
    not visible.  Energy check (estimate_energy.py): undiluted starch in
    vmh_eu_average is ~1/20 of the custom diet's, i.e. consistent with /24 of a
    ~2000 kcal diet, but its non-starch compounds are not.  Scale: unverified.
  * himalaya / guerrero_mountains / baobab_honey_antelope: unverified.
Absolute fluxes from different groups are therefore NOT put on one axis.  The
script produces unit-free views (within-diet composition, presence/absence,
Spearman rank correlation) for all diets, and absolute views per unit group.
"""
import argparse
import io
import re
import sys
import zipfile
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

# --------------------------------------------------------------------------- #
# Configuration
# --------------------------------------------------------------------------- #
DIET_DESIGNER_DIETS = set()   # names of diets loaded via --fluxes-csv (filled at load time)


def unit_group(name: str) -> str:
    if name.startswith("vmh_") or name in DIET_DESIGNER_DIETS:
        return "VMH-derived (scale unverified)"
    if name.startswith(("western_diet", "breast_milk", "agora_s12")):
        return "MICOM per-biomass (mmol/gDW/h)"
    return "Other diets (units unverified)"


# Metabolites excluded from "composition" fractions (they swamp everything).
NON_NUTRIENT_CLASSES = {"Water & gases"}

AMINO = {"ala", "arg", "asn", "asp", "cys", "gln", "glu", "gly", "his", "ile", "leu",
         "lys", "met", "phe", "pro", "ser", "thr", "trp", "tyr", "val"}

# glycine has no _L/_D suffix; the rest are amino-acid derivatives / polyamines
AMINO_EXTRA = {"4hpro_LT", "gly", "orn", "ala_B", "4abut", "4hpro", "5oxpro", "cgly", "Lcystin",
               "metsox_S_L", "glyb", "phpyr", "26dap_M", "ptrc", "spmd"}

CLASS_SETS = {
    "Water & gases": {"h2o", "h", "h2", "o2", "co2", "n2", "n2o", "no", "no2", "no3",
                      "nh4", "h2s", "ph2s"},
    "Minerals & ions": {"ca2", "cl", "k", "na1", "mg2", "fe2", "fe3", "fe3dcit", "zn2",
                        "cu2", "mn2", "cobalt2", "mobd", "sel", "so4", "so3", "tsul",
                        "pi", "ppi", "tet", "i", "cu2x"},
    "Vitamins & cofactors": {"10fthf", "5mthf", "thf", "fol", "btn", "ribflv", "rbflvrd",
                             "thm", "thmmp", "ncam", "nac", "pnto_R", "pydx", "pydx5p",
                             "pydam", "pydxn", "cbl1", "adpcbl", "ascb_L", "retinol",
                             "retinal", "phyQ", "mqn7", "mqn8", "q8", "2dmmq8", "caro",
                             "coa", "dpcoa", "nadp", "nmn", "pheme", "sheme", "chor",
                             "amet", "4abz", "4hbz", "pime", "hqn", "gthrd", "gthox",
                             "vitEg", "vitEa", "vitEaten", "vitEd", "vitA", "bvite", "vitK1", "vitd3",
                             "vitd2", "adocbl"},
    "Nucleotides & bases": {"ade", "adn", "gua", "gsn", "dgsn", "xan", "xtsn", "hxan",
                            "ura", "uri", "duri", "thymd", "cytd", "dcyt", "csn", "din",
                            "dad_2", "ins", "amp", "cmp", "gtp", "datp", "dtmp",
                            "orot", "urate", "urea"},
    "Short-chain & organic acids": {"ac", "but", "ppa", "for", "lac_L", "lac_D", "succ",
                                    "fum", "mal_L", "cit", "oxa", "pyr", "akg", "oaa",
                                    "acac", "bhb", "2hb", "2obut", "3mop", "glutar",
                                    "glyclt", "glcn", "2ddglcn", "actn_R", "isocapr"},
    "Fatty acids & lipids": {"hdca", "hdcea", "ocdca", "ocdcea", "ttdca", "ttdcea", "ddca",
                             "arach", "arachd", "lnlc", "lnlnca", "lnlncg", "octa", "dha",
                             "12dgr180", "chsterol", "lanost", "sql", "chol", "glyc",
                             # Diet Designer fatty acids with no AGORA counterpart
                             "hca", "dca", "hpdca", "ptdca", "gndc", "vacc", "dlnlncg", "mrgrl",
                             "hpdca10Z", "dlnlc", "crvnc", "dcsa", "clpnd", "tecsa", "epa", "ercc",
                             "nrvnc", "strdnc",
                             "glyc3p", "stys", "tchola", "gchola", "etha", "ethamp"},
    "Simple sugars & sugar alcohols": {"glc", "glc_D", "fru", "gal", "lcts", "sucr", "malt",
                                       "tre", "melib", "raffin", "cellb", "xyl_D", "arab_D",
                                       "arab_L", "rib_D", "man", "fuc_L", "rmn", "mnl",
                                       "sbt_D", "xylt", "rbt", "abt", "inost", "mi4p_D",
                                       "xylu_L", "glyc_R", "acgal", "acgam", "acmana",
                                       "acnam", "gam", "galct_D", "udpglcur", "bglc",
                                       "homogal", "drib"},
    "Complex carbs & fibre": {"starch1200", "strch1", "strch2", "amylose300", "amylopect900",
                              "inulin", "pect", "cellul", "xylan", "arabinoxyl", "arabinogal",
                              "arabinan101", "levan1000", "dextran40", "dextrin", "pullulan1200",
                              "lichn", "lmn30", "galmannan", "glcmannan", "kestopt", "xyluglc",
                              "rhamnogalurI", "rhamnogalurII", "malttr", "malthx", "mantr",
                              "arabttr", "chtbs"},
    "Host glycans (mucin etc.)": {"dsT_antigen", "cspg_b", "hspg", "glygn2", "glygn4", "glygn5"},
    "Solvents": {"meoh", "etoh", "indole", "12ppd_S", "fald", "gcald", "acald"},
}
CLASS_RULES = [  # regex fallbacks, evaluated after exact sets
    ("Host glycans (mucin etc.)", re.compile(r"^(g?n?core\d|MGlcn)")),
    ("Peptides", re.compile(r"^(" + "|".join(AMINO | {"ala", "gly"}) + r")(" + "|".join(AMINO) + r")")),
]

CLASS_ORDER = [
    "Amino acids", "Peptides", "Simple sugars & sugar alcohols", "Complex carbs & fibre",
    "Host glycans (mucin etc.)", "Fatty acids & lipids", "Short-chain & organic acids",
    "Vitamins & cofactors", "Nucleotides & bases", "Minerals & ions", "Solvents",
    "Water & gases", "Other",
]
CLASS_COLORS = dict(zip(CLASS_ORDER, sns.color_palette("tab20", len(CLASS_ORDER))))


# --------------------------------------------------------------------------- #
# Parsing & harmonisation
# --------------------------------------------------------------------------- #
def normalise_id(global_id: str) -> str:
    """EX_ala_L(e) (AGORA) and EX_ala__L_e (CarveMe/BiGG) -> ala_L."""
    g = re.sub(r"^EX_", "", str(global_id))
    g = re.sub(r"(\(e\)|_e)$", "", g)
    return g.replace("__", "_")


def classify(met: str) -> str:
    base = met.split("_")[0]
    if (met.endswith(("_L", "_D")) and base in AMINO) or met in AMINO_EXTRA:
        return "Amino acids"
    for cls, ids in CLASS_SETS.items():
        if met in ids:
            return cls
    for cls, rx in CLASS_RULES:
        if rx.match(met):
            return cls
    return "Other"


def read_qza(path: Path) -> pd.DataFrame:
    with zipfile.ZipFile(path) as z:
        members = [n for n in z.namelist() if n.endswith("data/medium.csv")]
        if len(members) != 1:
            raise ValueError(f"{path.name}: expected one data/medium.csv, found {members}")
        df = pd.read_csv(io.BytesIO(z.read(members[0])))
    if "global_id" not in df or "flux" not in df:
        raise ValueError(f"{path.name}: missing global_id/flux columns ({list(df.columns)})")
    out = pd.DataFrame({
        "metabolite": df["global_id"].map(normalise_id),
        "flux": df["flux"].astype(float),
        "name": df["name"] if "name" in df else np.nan,
    })
    # CarveMe medium repeats a few rows (pheme, thm, spmd): collapse, don't sum.
    return out.groupby("metabolite", as_index=False).agg(flux=("flux", "max"), name=("name", "first"))


# VMH Diet Designer names -> the AGORA exchange IDs used by the .qza media
DIET_DESIGNER_ALIASES = {
    "oleic": "ocdcea",      # oleate, C18:1
    "palmtl": "hdcea",      # palmitoleate, C16:1
    "tedca": "ttdca",       # myristate, C14:0
    "xol": "chsterol",      # cholesterol
    "occa": "octa",         # octanoate, C8:0
    "4hpro_LT": "4hpro",    # trans-4-hydroxy-L-proline
}


def read_fluxes_csv(path: Path) -> pd.DataFrame:
    """Two-column CSV: reaction, flux 'in mmol/(human*day)' (VMH Diet Designer export).

    Duplicated reactions are SUMMED, [e]/EX_ are stripped, IDs are aliased to AGORA
    names, zero-flux rows are dropped (so they don't count as 'present'), and fluxes
    are kept in their native units (mmol/day; no /24).
    """
    df = pd.read_csv(path)
    rxn, flux = df.columns[:2]
    ids = (df[rxn].astype(str).str.replace(r"^EX_", "", regex=True)
           .str.replace(r"\[e\]$", "", regex=True))
    ids = ids.map(lambda x: DIET_DESIGNER_ALIASES.get(x, x))
    out = pd.DataFrame({"metabolite": ids, "flux": df[flux].astype(float)})
    out = out.groupby("metabolite", as_index=False)["flux"].sum()
    out = out[out["flux"] > 0].copy()
    out["name"] = np.nan
    return out


def read_s12_csv(path: Path) -> dict:
    """CSV from extract_s12_table.py -> {diet name: DataFrame(metabolite, flux, name)}."""
    t = pd.read_csv(path)
    out = {}
    for col, nm in [("western", "agora_s12_western"), ("high_fiber", "agora_s12_high_fiber")]:
        d = t[["metabolite", col, "name"]].rename(columns={col: "flux"})
        out[nm] = d[d["flux"] > 0].reset_index(drop=True)
    return out


def load_all(paths, fluxes_csvs=(), s12_csv=None, medium_csvs=()):
    files = []
    for p in map(Path, paths):
        files += sorted(p.glob("*.qza")) if p.is_dir() else [p]
    if not files:
        sys.exit("No .qza files found.")
    long = []
    for f in files:
        diet = re.sub(r"^[0-9a-f]{8}-", "", f.stem)   # strip upload-hash prefixes
        d = read_qza(f)
        d["diet"] = diet
        d["unit_group"] = unit_group(diet)
        long.append(d)
    for spec in fluxes_csvs:                       # "path.csv" or "path.csv:diet_name"
        path, _, nm = spec.partition(":")
        nm = nm or Path(path).stem.split("-", 1)[-1]
        DIET_DESIGNER_DIETS.add(nm)
        d = read_fluxes_csv(Path(path))
        d["diet"], d["unit_group"] = nm, unit_group(nm)
        long.append(d)
    for spec in medium_csvs:                       # MICOM-format medium CSV, "path.csv[:diet_name]"
        path, _, nm = spec.partition(":")
        nm = nm or Path(path).stem
        m = pd.read_csv(path)
        d = pd.DataFrame({"metabolite": m["metabolite"].map(normalise_id), "flux": m["flux"].astype(float),
                          "name": m["name"] if "name" in m else np.nan})
        d = d.groupby("metabolite", as_index=False).agg(flux=("flux", "max"), name=("name", "first"))
        d["diet"], d["unit_group"] = nm, unit_group(nm)
        long.append(d)
    if s12_csv:
        for nm, d in read_s12_csv(Path(s12_csv)).items():
            d["diet"], d["unit_group"] = nm, unit_group(nm)
            long.append(d)
    long = pd.concat(long, ignore_index=True)
    long["class"] = long["metabolite"].map(classify)
    names = long.dropna(subset=["name"]).drop_duplicates("metabolite").set_index("metabolite")["name"]
    long["name"] = long["metabolite"].map(names).fillna(long["metabolite"])
    return long


# --------------------------------------------------------------------------- #
# Plots
# --------------------------------------------------------------------------- #
def order_diets(long):
    return (long.drop_duplicates("diet").sort_values(["unit_group", "diet"])["diet"].tolist())


def group_palette(long):
    groups = sorted(long["unit_group"].unique())
    return dict(zip(groups, sns.color_palette("Set2", len(groups))))


def plot_overview(long, outdir):
    """Metabolite counts per class + Jaccard similarity of component sets."""
    diets = order_diets(long)
    gp = group_palette(long)
    counts = (long.groupby(["diet", "class"]).size().unstack(fill_value=0)
              .reindex(index=diets, columns=[c for c in CLASS_ORDER if c in long["class"].unique()]))
    pres = long.assign(v=1).pivot_table(index="metabolite", columns="diet", values="v", fill_value=0)[diets]
    inter = pres.T.values @ pres.values
    union = pres.sum().values[:, None] + pres.sum().values[None, :] - inter
    jac = pd.DataFrame(inter / union, index=diets, columns=diets)

    fig, axes = plt.subplots(1, 2, figsize=(17, 6.5), gridspec_kw={"width_ratios": [1.1, 1]})
    counts.plot.barh(stacked=True, ax=axes[0], color=[CLASS_COLORS[c] for c in counts.columns], width=0.8)
    axes[0].invert_yaxis(); axes[0].set_xlabel("Number of medium components")
    axes[0].set_title("Medium size and composition by nutrient class")
    axes[0].legend(fontsize=8, bbox_to_anchor=(1.0, 1), loc="upper left", frameon=False)
    for tick in axes[0].get_yticklabels():
        tick.set_color(gp[unit_group(tick.get_text())])
    sns.heatmap(jac, ax=axes[1], vmin=0, vmax=1, cmap="viridis", annot=True, fmt=".2f",
                annot_kws={"size": 7}, cbar_kws={"label": "Jaccard similarity"}, square=True)
    axes[1].set_title("Overlap of component sets"); axes[1].set_xlabel(""); axes[1].set_ylabel("")
    handles = [plt.Line2D([], [], marker="s", ls="", color=c, label=g) for g, c in gp.items()]
    fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False, title="Label colour = unit group")
    fig.tight_layout(rect=(0, 0.07, 1, 1))
    fig.savefig(outdir / "01_overview_size_and_overlap.png", dpi=200); plt.close(fig)


def plot_composition(long, outdir):
    """Unit-free: share of molar flux contributed by each class (two views)."""
    diets = order_diets(long)
    gp = group_palette(long)
    views = [
        ("All components except water & gases", NON_NUTRIENT_CLASSES, "composition_fraction_all.csv"),
        ("Organic nutrients only (also excluding ions & solvents, e.g. the\nWestern diets' flat methanol=10 placeholder)",
         NON_NUTRIENT_CLASSES | {"Minerals & ions", "Solvents"}, "composition_fraction_organic.csv"),
    ]
    fig, axes = plt.subplots(1, 2, figsize=(18, 6.5), sharey=True)
    for ax, (title, excl, csv) in zip(axes, views):
        nut = long[~long["class"].isin(excl)]
        tot = nut.groupby("diet")["flux"].sum()
        cols = [c for c in CLASS_ORDER if c in nut["class"].unique()]
        frac = (nut.groupby(["diet", "class"])["flux"].sum().div(tot, level="diet")
                .unstack(fill_value=0).reindex(index=diets, columns=cols))
        frac.plot.barh(stacked=True, ax=ax, color=[CLASS_COLORS[c] for c in cols], width=0.8, legend=False)
        ax.set_xlim(0, 1); ax.set_xlabel("Fraction of total molar flux"); ax.set_title(title, fontsize=10)
        for tick in ax.get_yticklabels():
            tick.set_color(gp[unit_group(tick.get_text())])
        frac.to_csv(outdir / csv)
    axes[0].invert_yaxis(); axes[0].set_ylabel("")
    h = [plt.Rectangle((0, 0), 1, 1, color=CLASS_COLORS[c]) for c in CLASS_ORDER if c not in NON_NUTRIENT_CLASSES]
    fig.legend(h, [c for c in CLASS_ORDER if c not in NON_NUTRIENT_CLASSES], loc="center right", fontsize=8, frameon=False)
    fig.suptitle("Nutrient composition of each medium (unit-free, comparable across all diets)")
    fig.tight_layout(rect=(0, 0, 0.87, 0.97)); fig.savefig(outdir / "02_composition_fraction.png", dpi=200); plt.close(fig)


def plot_absolute_by_class(long, outdir, rescaled):
    """Absolute class totals, one panel per unit group (or one panel if rescaled)."""
    nut = long[~long["class"].isin(NON_NUTRIENT_CLASSES)]
    tab = nut.groupby(["unit_group", "diet", "class"])["flux"].sum().reset_index()
    groups = sorted(tab["unit_group"].unique())
    if rescaled:   # vmh already divided by biomass; merge VMH with MICOM groups
        pass
    fig, axes = plt.subplots(len(groups), 1, figsize=(11, 3.2 * len(groups) + 1), squeeze=False)
    for ax, g in zip(axes[:, 0], groups):
        t = tab[tab["unit_group"] == g].pivot(index="diet", columns="class", values="flux").fillna(0)
        t = t[[c for c in CLASS_ORDER if c in t.columns]]
        t.plot.barh(stacked=True, ax=ax, color=[CLASS_COLORS[c] for c in t.columns], width=0.8, legend=False)
        ax.invert_yaxis(); ax.set_ylabel(""); ax.set_title(g, loc="left", fontsize=10)
        ax.set_xlabel("Total flux (native units of this group; water/gases excluded)")
    h, l = [], []
    for c in CLASS_ORDER:
        if c in nut["class"].unique():
            h.append(plt.Rectangle((0, 0), 1, 1, color=CLASS_COLORS[c])); l.append(c)
    fig.legend(h, l, loc="center right", fontsize=8, frameon=False)
    fig.suptitle("Absolute nutrient supply by class (compare only within a panel)")
    fig.tight_layout(rect=(0, 0, 0.82, 0.97))
    fig.savefig(outdir / "03_absolute_by_class_per_unit_group.png", dpi=200); plt.close(fig)


def plot_heatmap(long, outdir, top_n=70):
    """log10 flux for the most variable / highest-ranked metabolites, grouped by class."""
    diets = order_diets(long)
    gp = group_palette(long)
    m = long[long["class"] != "Water & gases"].pivot(index="metabolite", columns="diet", values="flux")[diets]
    rank = m.rank(pct=True)                       # within-diet percentile: unit-free
    score = rank.mean(axis=1, skipna=True) * (m.notna().sum(axis=1) / len(diets))
    keep = score.sort_values(ascending=False).head(top_n).index
    cls = long.drop_duplicates("metabolite").set_index("metabolite")["class"]
    nm = long.drop_duplicates("metabolite").set_index("metabolite")["name"]
    order = sorted(keep, key=lambda x: (CLASS_ORDER.index(cls[x]), -score[x]))
    data = rank.loc[order]
    data.index = [f"{nm[i][:34]} ({i})" if nm[i] != i else i for i in order]
    fig, ax = plt.subplots(figsize=(10, 0.23 * len(order) + 3))
    sns.heatmap(data, cmap="mako_r", vmin=0, vmax=1, ax=ax, linewidths=0.3, linecolor="white",
                cbar_kws={"label": "Within-diet percentile of flux (unit-free)", "shrink": 0.5},
                mask=data.isna())
    ax.set_facecolor("#e6e6e6"); ax.set_xlabel(""); ax.set_ylabel("")
    ax.tick_params(axis="y", labelsize=7); plt.setp(ax.get_xticklabels(), rotation=60, ha="right")
    for tick, met in zip(ax.get_yticklabels(), order):
        tick.set_color(CLASS_COLORS[cls[met]])
    for tick in ax.get_xticklabels():
        tick.set_color(gp[unit_group(tick.get_text())])
    ax.set_title(f"Top {len(order)} metabolites: relative abundance in each medium\n"
                 "(grey = absent from that medium; row label colour = nutrient class)")
    fig.tight_layout(); fig.savefig(outdir / "04_heatmap_top_metabolites.png", dpi=200); plt.close(fig)


def plot_rank_correlation(long, outdir):
    """Spearman correlation of flux across shared metabolites: invariant to unit scaling."""
    diets = order_diets(long)
    m = long[long["class"] != "Water & gases"].pivot(index="metabolite", columns="diet", values="flux")[diets]
    rho = pd.DataFrame(index=diets, columns=diets, dtype=float); n = rho.copy()
    for a in diets:
        for b in diets:
            both = m[[a, b]].dropna()
            rho.loc[a, b] = both.corr(method="spearman").iloc[0, 1] if len(both) > 3 else np.nan
            n.loc[a, b] = len(both)
    fig, ax = plt.subplots(figsize=(8.5, 7))
    sns.heatmap(rho, vmin=-0.2, vmax=1, cmap="rocket_r", annot=True, fmt=".2f", annot_kws={"size": 7},
                square=True, ax=ax, cbar_kws={"label": "Spearman ρ (shared metabolites)"})
    ax.set_title("Do two media rank their shared nutrients the same way?\n(rank-based, so unit differences don't matter)")
    ax.set_xlabel(""); ax.set_ylabel("")
    fig.tight_layout(); fig.savefig(outdir / "05_rank_correlation.png", dpi=200); plt.close(fig)
    rho.to_csv(outdir / "spearman_rho.csv"); n.astype(int).to_csv(outdir / "shared_metabolite_counts.csv")


def plot_key_nutrients(long, outdir):
    """Panels for classes people usually care about: each metabolite as share of the diet's total."""
    diets = order_diets(long)
    gp = group_palette(long)
    nut = long[~long["class"].isin(NON_NUTRIENT_CLASSES)].copy()
    nut["share"] = nut["flux"] / nut.groupby("diet")["flux"].transform("sum")
    panels = ["Amino acids", "Vitamins & cofactors", "Short-chain & organic acids",
              "Fatty acids & lipids", "Simple sugars & sugar alcohols", "Complex carbs & fibre"]
    fig, axes = plt.subplots(2, 3, figsize=(21, 12))
    for ax, cls in zip(axes.ravel(), panels):
        s = nut[nut["class"] == cls]
        if s.empty:
            ax.axis("off"); continue
        mat = s.pivot(index="metabolite", columns="diet", values="share").reindex(columns=diets)
        mat = mat.loc[mat.max(axis=1).sort_values(ascending=False).head(25).index]
        nm = s.drop_duplicates("metabolite").set_index("metabolite")["name"]
        mat.index = [nm[i][:26] if nm[i] != i else i for i in mat.index]
        sns.heatmap(np.log10(mat.where(mat > 0)), ax=ax, cmap="mako_r", mask=mat.isna() | (mat <= 0),
                    cbar_kws={"label": "log10 share of diet", "shrink": 0.6}, linewidths=0.2, linecolor="white")
        ax.set_facecolor("#e6e6e6"); ax.set_title(cls); ax.set_xlabel(""); ax.set_ylabel("")
        ax.tick_params(axis="y", labelsize=7)
        plt.setp(ax.get_xticklabels(), rotation=60, ha="right", fontsize=8)
        for t in ax.get_xticklabels():
            t.set_color(gp[unit_group(t.get_text())])
    fig.suptitle("Key nutrient classes: share of each diet's non-water molar flux (grey = absent)", y=1.0)
    fig.tight_layout(); fig.savefig(outdir / "06_key_nutrient_classes.png", dpi=200); plt.close(fig)


def plot_absolute_scatter(long, outdir, rescaled):
    """Only meaningful on one axis if units were harmonised via --vmh-biomass-gdw."""
    if not rescaled:
        return
    diets = order_diets(long)
    m = long[long["class"] != "Water & gases"].pivot(index="metabolite", columns="diet", values="flux")[diets]
    fig, ax = plt.subplots(figsize=(10, 6))
    sns.boxplot(data=np.log10(m.where(m > 0)), orient="h", ax=ax, fliersize=2)
    ax.set_xlabel("log10 flux (mmol/gDW/h, VMH diets rescaled)"); ax.set_title("Flux distribution per medium")
    fig.tight_layout(); fig.savefig(outdir / "07_flux_distributions_harmonised.png", dpi=200); plt.close(fig)


# --------------------------------------------------------------------------- #
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("inputs", nargs="+", help=".qza files and/or directories containing them")
    ap.add_argument("-o", "--outdir", default="media_comparison")
    ap.add_argument("--fluxes-csv", action="append", default=[], metavar="CSV[:NAME]",
                    help="extra diet in VMH Diet Designer format (reaction, flux in mmol/(human*day)); "
                         "duplicates are summed; kept in mmol/day. May be repeated.")
    ap.add_argument("--medium-csv", action="append", default=[], metavar="CSV[:NAME]",
                    help="a MICOM-format medium CSV (columns metabolite, flux), e.g. from make_micom_media.ipynb. May be repeated.")
    ap.add_argument("--s12-csv", default=None, help="CSV from extract_s12_table.py (adds the two AGORA S12 diets)")
    args = ap.parse_args()

    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)
    sns.set_theme(style="white", context="notebook")
    long = load_all(args.inputs, args.fluxes_csv, args.s12_csv, args.medium_csv)
    rescaled = False

    # tidy + wide tables
    long.sort_values(["diet", "class", "metabolite"]).to_csv(outdir / "all_media_long.csv", index=False)
    wide = long.pivot(index="metabolite", columns="diet", values="flux")
    meta = long.drop_duplicates("metabolite").set_index("metabolite")[["name", "class"]]
    meta.join(wide).sort_values(["class", "name"]).to_csv(outdir / "all_media_wide.csv")
    summ = long.groupby("diet").agg(unit_group=("unit_group", "first"), n_components=("metabolite", "size"),
                                    total_flux=("flux", "sum"), median_flux=("flux", "median"),
                                    max_flux=("flux", "max"))
    summ["water_flux"] = long[long.metabolite == "h2o"].set_index("diet")["flux"]
    summ["unclassified_components"] = long[long["class"] == "Other"].groupby("diet").size()
    summ.fillna(0).to_csv(outdir / "summary_per_diet.csv")
    print(summ.to_string(float_format=lambda x: f"{x:.4g}"))
    unc = sorted(long.loc[long["class"] == "Other", "metabolite"].unique())
    print(f"\n{len(unc)} metabolites fell into 'Other': {unc}")

    plot_overview(long, outdir)
    plot_composition(long, outdir)
    plot_absolute_by_class(long, outdir, rescaled)
    plot_heatmap(long, outdir)
    plot_rank_correlation(long, outdir)
    plot_key_nutrients(long, outdir)
    plot_absolute_scatter(long, outdir, rescaled)
    print(f"\nWrote tables and plots to {outdir.resolve()}")


if __name__ == "__main__":
    main()
