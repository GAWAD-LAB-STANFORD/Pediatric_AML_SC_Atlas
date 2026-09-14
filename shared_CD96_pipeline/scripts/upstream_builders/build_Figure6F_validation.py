#!/usr/bin/env python
# ====================================================================
# build_Figure6F_validation.py  |  Upstream builder (regenerates source_data from raw inputs)
# reproducibility header (auto-generated; edit freely below this block)
# Inputs : (reads upstream object / raw matrix — see body)
# Outputs: (writes to ../figures or ../source_data — see body)
# Run    : python build_Figure6F_validation.py   (or: bash run_all.sh)  |  Deps: see REPRODUCIBILITY_CHECKLIST.md
# Source data for this figure is in __SUBMISSION_PACKAGE/08_Source_Data/.
# ====================================================================

"""Figure 6F -- CD96 across AML genotypes, validated in BOTH cohorts.

LEFT  (single cell, 28-patient atlas): % leukemic cells CD96+ for EVERY
       cytogenetic subtype.  Each bar is annotated with its patient n; subtypes
       with <3 patients are faded/italic (shown for completeness but NOT
       statistically robust -- several "cold" ones are single-patient).
RIGHT (TARGET bulk, independent): raw CD96 expression (log2 CPM) by fusion-
       partner group, restricted to samples with >=80% marrow blasts -- i.e.
       near-pure leukemia, so bulk expression ~ leukemic expression with no
       deconvolution/normalization needed.  Groups with >=3 high-blast samples
       are shown; rarer ones fall into Other.

Both cohorts show the same result: CD96 marks most AML genotypes but is COLD in
KMT2Ar/MLL -- the one robust, cross-cohort-validated coverage gap.  The two
cohorts are quantitatively concordant across the robust (>=3-patient) genotypes
(Spearman / Pearson printed + saved to Figure6F_crosscohort_corr.csv); the only
discordant subtype is single-patient NUP98-NSD1 (cold in that one scRNA patient,
but CD96-high in TARGET NUP98r bulk) -- a documented single-patient artifact.

Raw expression only (no normalization).  Pure data; no fabricated values.
"""
import os
import sys
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
import matplotlib as mpl

HERE = os.path.dirname(os.path.abspath(__file__))
SUPP = os.path.join(HERE, "supplements")
sys.path.insert(0, SUPP)
from _fig7_style import PALETTE, RC, clean_ax, save_fig

MERGED = os.path.join(HERE, "Figure7F_TARGET_CD96_merged.csv")
SCRNA = os.path.join(HERE, "Figure6F_AMLvsNormal.csv")            # CD96 %pos by subtype + Normal BM
NPAT = os.path.join(HERE, "Figure6F_scrna_subtype_npatients.csv")  # patients per scRNA subtype
NORMAL_CSV = os.path.join(HERE, "Figure6F_TARGET_normalBM_CD96.csv")  # TARGET type-14A normal BM
BLAST = "Bone marrow leukemic blast percentage (%)"
BLAST_MIN = 80.0
MIN_PATIENTS = 3
COLD_THRESH = 15.0     # %CD96+ below this = CD96-cold (descriptive, for bar colour)

# scRNA subtype -> display label (ALL subtypes shown; n_patients annotated)
SC_LABEL = {
    "BCR/ABL": "BCR::ABL1", "RUNX1/RUNX1T1": "t(8;21)",
    "t(2;3)(p15;q26.2)": "t(2;3) MECOM", "t(7;14)(q21;q32)": "t(7;14)",
    "PML/RARA": "PML-RARA", "CN": "CN", "CBFB/MYH11": "inv(16)",
    "Tri(8)": "Trisomy 8", "Tri(15)": "Trisomy 15", "MLLr": "KMT2Ar",
    "del7q": "del(7q)", "MYB/GATA1": "MYB-GATA1", "NUP98/NSD1": "NUP98-NSD1",
    "Tri(8)/MLLr": "Tri8 + KMT2Ar",
}
# scRNA subtype -> matching TARGET fusion group (for cross-cohort correlation)
SC_TO_TARGET = {
    "RUNX1/RUNX1T1": "t(8;21)", "CBFB/MYH11": "inv(16)", "MLLr": "KMT2Ar",
    "CN": "CN", "PML/RARA": "PML-RARA", "NUP98/NSD1": "NUP98r",
}


def fam(row):
    fz = str(row.get("Gene Fusion", "")).strip()
    if fz and fz.lower() != "nan":
        if fz.startswith("KMT2A-") or fz == "other MLL":
            return "KMT2Ar"
        if fz == "RUNX1-RUNX1T1":
            return "t(8;21)"
        if fz == "CBFB-MYH11":
            return "inv(16)"
        if fz.startswith("NUP98-"):
            return "NUP98r"
        if fz == "CBFA2T3-GLIS2":
            return "CBFA2T3::GLIS2"
        if fz == "MLLT10-PICALM":
            return "PICALM::MLLT10"
        if fz == "DEK-NUP214":
            return "DEK::NUP214"
        if fz == "NPM1-MLF1":
            return "NPM1::MLF1"
        return "Other"
    return {"Normal": "CN", "t(8;21)": "t(8;21)", "inv(16)": "inv(16)",
            "PML-RARA": "PML-RARA", "MLL": "KMT2Ar"}.get(
        str(row.get("Primary Cytogenetic Code", "")).strip(), "Other")


def main():
    # ---- single-cell: ALL cytogenetic subtypes (n_patients annotated) ----
    cd96 = pd.read_csv(SCRNA, index_col=0).loc["CD96"]
    npat = pd.read_csv(NPAT).set_index("subtype")["n_patients"]
    sc_rows = [(s, SC_LABEL.get(s, s), float(cd96[s]), int(npat.get(s, 0)))
               for s in cd96.index if s != "Normal BM"]
    sc = (pd.DataFrame(sc_rows, columns=["sub", "geno", "pct", "npat"])
          .sort_values("pct").reset_index(drop=True))

    # ---- TARGET: >=BLAST_MIN% blasts, raw CD96 log2 CPM, by fusion group,
    #      plus the same-platform normal-BM control (TARGET type 14A) ----
    m = pd.read_csv(MERGED)
    m = m[(m["CD96_log2CPM"].notna()) & (m[BLAST] >= BLAST_MIN)].copy()
    m["grp"] = m.apply(fam, axis=1)
    vc = m["grp"].value_counts()
    keep = set(vc[vc >= MIN_PATIENTS].index)
    m["grp"] = np.where(m["grp"].isin(keep), m["grp"], "Other")
    gvals = {g: m.loc[m["grp"] == g, "CD96_log2CPM"].to_numpy()
             for g in m["grp"].unique()}
    gmed = {g: float(np.median(v)) for g, v in gvals.items()}
    gsize = {g: int(len(v)) for g, v in gvals.items()}
    order = sorted(gvals, key=lambda g: gmed[g])          # ascending by median
    COLD = {"KMT2Ar"}                                      # the validated CD96-cold gap
    print(f"  TARGET >={int(BLAST_MIN)}% blasts n={len(m)}")
    for g in order:
        print(f"    {g:16s} n={gsize[g]:>4d}  median={gmed[g]:.2f}")

    # ---- cross-cohort correlation: single-cell %CD96+ vs TARGET median CD96 ----
    pairs = [(SC_LABEL.get(s, s), float(cd96[s]), int(npat.get(s, 0)),
              gmed[SC_TO_TARGET[s]])
             for s in SC_TO_TARGET
             if s in cd96.index and SC_TO_TARGET[s] in gmed]
    cor = pd.DataFrame(pairs, columns=["geno", "sc_pct", "npat", "tgt_med"])
    rob = cor[cor["npat"] >= MIN_PATIENTS]
    pr_a, pp_a = stats.pearsonr(cor["sc_pct"], cor["tgt_med"])
    sr_a, sq_a = stats.spearmanr(cor["sc_pct"], cor["tgt_med"])
    pr_r, pp_r = stats.pearsonr(rob["sc_pct"], rob["tgt_med"])
    sr_r, sq_r = stats.spearmanr(rob["sc_pct"], rob["tgt_med"])
    cor.round(3).to_csv(os.path.join(HERE, "Figure6F_crosscohort_corr.csv"),
                        index=False)
    print("\n  cross-cohort match (scRNA %CD96+ vs TARGET median log2CPM):")
    print(cor.round(2).to_string(index=False))
    print(f"  ALL matched  (n={len(cor)}): Pearson r={pr_a:.2f} (p={pp_a:.3f}), "
          f"Spearman rho={sr_a:.2f} (p={sq_a:.3f})")
    print(f"  ROBUST >=3pt (n={len(rob)}): Pearson r={pr_r:.2f} (p={pp_r:.3f}), "
          f"Spearman rho={sr_r:.2f} (p={sq_r:.3f})")

    with mpl.rc_context(RC):
        fig, (axL, axR) = plt.subplots(
            1, 2, figsize=(12.4, 5.0),
            gridspec_kw=dict(width_ratios=[0.62, 1.0], wspace=0.34))

        # LEFT: single-cell %CD96+ across ALL cytogenetic subtypes.
        #   crimson = CD96-hot, grey = CD96-cold (<COLD_THRESH); solid = robust
        #   (>=3 pt), faded/italic = <3 pt (single-patient subtypes shown for
        #   completeness, NOT statistically robust).
        y = np.arange(len(sc))
        bars = axL.barh(y, sc["pct"], height=0.74, edgecolor="white",
                        color=["#9fb3c8" if v < COLD_THRESH else PALETTE["CD96"]
                               for v in sc["pct"]])
        for b, n in zip(bars, sc["npat"]):
            b.set_alpha(1.0 if n >= MIN_PATIENTS else 0.42)
        for yi, v, n in zip(y, sc["pct"], sc["npat"]):
            axL.text(v + 1.2, yi, f"{v:.0f}%  (n={n})", va="center", fontsize=6.8,
                     color="#1c2833" if n >= MIN_PATIENTS else "#95a5a6")
        axL.set_yticks(y)
        axL.set_yticklabels(sc["geno"], fontsize=8)
        for lab, n in zip(axL.get_yticklabels(), sc["npat"]):
            if n < MIN_PATIENTS:
                lab.set_color("#7f8c8d"); lab.set_style("italic")
        axL.set_xlabel("% leukemic cells CD96+")
        axL.set_xlim(0, max(sc["pct"]) * 1.34)
        axL.set_title("single cell", fontsize=9.0, loc="left", fontweight="bold")
        clean_ax(axL, "x")

        # RIGHT: TARGET raw CD96 by fusion group (crimson = CD96; grey = cold gap)
        data = [gvals[g] for g in order]
        bp = axR.boxplot(data, positions=range(len(order)), widths=0.6, vert=True,
                         showfliers=False, patch_artist=True,
                         medianprops=dict(color="#1c2833", lw=1.4))
        for patch, g in zip(bp["boxes"], order):
            patch.set_facecolor("#9fb3c8" if g in COLD else PALETTE["CD96"])
            patch.set_alpha(0.85); patch.set_edgecolor("#3b5870")
        for i, d in enumerate(data):
            axR.scatter(np.full(len(d), i) + np.linspace(-0.16, 0.16, len(d)),
                        d, s=4, color="#34495e", alpha=0.22, linewidths=0, zorder=4)
        axR.set_xticks(range(len(order)))
        axR.set_xticklabels([f"{g}\n(n={gsize[g]})" for g in order],
                            rotation=35, ha="right", fontsize=8)
        for lab, g in zip(axR.get_xticklabels(), order):
            if g in COLD:
                lab.set_color("#7f8c8d")
        axR.set_ylabel("CD96 expression (log2 CPM)")
        axR.set_title(f"TARGET bulk (≥{int(BLAST_MIN)}% blasts, n={len(m)})",
                      fontsize=9.5, loc="left", fontweight="bold")
        # cross-cohort concordance (robust >=3-patient genotypes)
        axR.text(0.015, 0.975,
                 f"scRNA vs TARGET, robust genotypes (n={len(rob)}):\n"
                 f"Spearman ρ = {sr_r:.2f},  Pearson r = {pr_r:.2f}",
                 transform=axR.transAxes, ha="left", va="top", fontsize=7.6,
                 color="#1c2833",
                 bbox=dict(boxstyle="round,pad=0.4", facecolor="#eef3f8",
                           edgecolor="#9fb3c8", lw=0.7))
        clean_ax(axR, "y")

        # (figure title/commentary in caption -- added by authors)
        save_fig(fig, HERE, "Figure6F_validation")
    print("\nDone -> Figure6F_validation.{pdf,png}")


if __name__ == "__main__":
    main()
