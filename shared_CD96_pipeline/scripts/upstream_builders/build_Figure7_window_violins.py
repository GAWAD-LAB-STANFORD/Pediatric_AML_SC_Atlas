#!/usr/bin/env python
# ====================================================================
# build_Figure7_window_violins.py  |  Upstream builder (regenerates source_data from raw inputs)
# reproducibility header (auto-generated; edit freely below this block)
# Inputs : (reads upstream object / raw matrix — see body)
# Outputs: (writes to ../figures or ../source_data — see body)
# Run    : python build_Figure7_window_violins.py   (or: bash run_all.sh)  |  Deps: see REPRODUCIBILITY_CHECKLIST.md
# Source data for this figure is in __SUBMISSION_PACKAGE/08_Source_Data/.
# ====================================================================

"""Per-target bulk-expression violins with a targetability call-out.

Builds TWO matched panels (same target order + y-scale, for direct comparison):
  Figure7_window_violins        -- ALL primary AML (>=80% blasts)
  Figure7_window_violins_nonMLL -- with KMT2A-rearranged (MLL) patients REMOVED

For each of the 9 surface targets:
  * violin = distribution of RAW bulk expression across the cohort (TARGET,
    log2 CPM).  No normal-marrow subtraction.
  * horizontal dashed line = the "likely targetable" expression threshold.
  * number ABOVE each violin = % of that cohort's samples above the line.

Violins are COLOURED by the stem/progenitor-toxicity flag from single cell:
  HSPC or Myeloid_Pro positive > TOX_THR%  -> grey  ("likely HSPC/myeloid-toxic")
  at/below                                 -> green ("HSPC/Myeloid-sparing")
  CD96 kept crimson (the lead).  HSPC/Myeloid % printed under each target.

KMT2Ar is CD96's one validated coverage gap, so the non-MLL panel shows how much
of the REST of pediatric AML CD96 reaches.  Inputs already on disk:
  Figure7D_target_by_PAC.csv  -- scRNA % positive per PAC (0_HSPC, Myeloid_Pro)
  TARGET STAR counts (filtered) + clinical blast% / fusion / cytogenetics
No fabricated values.
"""
import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl

HERE = os.path.dirname(os.path.abspath(__file__))
SUPP = os.path.join(HERE, "supplements")
sys.path.insert(0, SUPP)
from _fig7_style import PALETTE, RC, clean_ax, save_fig

PAC = os.path.join(HERE, "Figure7D_target_by_PAC.csv")
FILT = os.path.join(HERE, "..", "ALSF_AML_plot", "H5AD",
                    "NCBI_TARGET_Star_counts_filtered.csv")
MERGED = os.path.join(HERE, "Figure7F_TARGET_CD96_merged.csv")
BLAST = "Bone marrow leukemic blast percentage (%)"
BLAST_MIN = 80.0
TOX_THR = 5.0           # % HSPC/Myeloid_Pro positive above which = likely toxic
TARGET_THR = 5.0        # log2 CPM "likely targetable" line

ENS = {"CD96": "ENSG00000153283", "CD33": "ENSG00000105383",
       "IL3RA": "ENSG00000185291", "CLEC12A": "ENSG00000172322",
       "FLT3": "ENSG00000122025", "MSLN": "ENSG00000102854",
       "CD70": "ENSG00000125726", "ADGRG1": "ENSG00000205336",
       "IL1RAP": "ENSG00000196083"}
GENES = list(ENS)
GLABEL = {"CD96": "CD96", "CD33": "CD33", "IL3RA": "CD123\n(IL3RA)",
          "CLEC12A": "CLL-1\n(CLEC12A)", "FLT3": "FLT3", "MSLN": "MSLN",
          "CD70": "CD70", "ADGRG1": "GPR56\n(ADGRG1)", "IL1RAP": "IL1RAP"}


def is_mll(row):
    """KMT2A-rearranged: any KMT2A/MLL fusion, or (no fusion &) cyto code MLL."""
    fz = str(row.get("Gene Fusion", "")).strip()
    if fz and fz.lower() != "nan":
        return fz.startswith("KMT2A-") or fz == "other MLL"
    return str(row.get("Primary Cytogenetic Code", "")).strip() == "MLL"


def bulk_matrix(exclude_mll=False):
    """per-sample bulk log2 CPM (high-blast primary AML) x 9 genes."""
    raw = pd.read_csv(FILT).set_index("ensembl_gene_id")
    raw.index = raw.index.astype(str).str.split(".").str[0]
    g9 = raw.loc[[ENS[g] for g in GENES]].T
    g9.columns = GENES
    m = pd.read_csv(MERGED).set_index("sample")
    ids = [s for s in g9.index if s in m.index
           and pd.notna(m.loc[s, BLAST]) and m.loc[s, BLAST] >= BLAST_MIN]
    if exclude_mll:
        ids = [s for s in ids if not is_mll(m.loc[s])]
    lib = m.loc[ids, "lib_size"].astype(float)
    return np.log2(g9.loc[ids].div(lib, axis=0) * 1e6 + 1)


def col(g, toxic):
    if g == "CD96":
        return PALETTE["CD96"]
    return "#95a5a6" if toxic[g] else "#1e8449"


def build(L, order, xsafe, ytop, out_stem, cohort_label):
    n = len(L)
    med = L.median()
    toxic = {g: float(xsafe[g]) > TOX_THR for g in GENES}
    pct_tgt = {g: 100.0 * (L[g] > TARGET_THR).mean() for g in GENES}
    pd.DataFrame({"bulk_med": med, "hspc_myeloid_max": xsafe.reindex(GENES),
                  "toxic": pd.Series(toxic), "pct_above_thr": pd.Series(pct_tgt)}
                 ).round(2).to_csv(os.path.join(HERE, out_stem + ".csv"))
    print(f"\n  [{cohort_label}]  n={n}")
    for g in order:
        print(f"    {g:8s} med={med[g]:5.2f}  targetable={pct_tgt[g]:5.1f}%")

    with mpl.rc_context(RC):
        fig, ax = plt.subplots(figsize=(11.6, 6.2))
        pos = np.arange(len(order))
        data = [L[g].values for g in order]
        parts = ax.violinplot(data, positions=pos, widths=0.82,
                              showextrema=False, showmedians=True)
        for body, g in zip(parts["bodies"], order):
            body.set_facecolor(col(g, toxic)); body.set_edgecolor("#3b5870")
            body.set_alpha(0.85); body.set_linewidth(0.6)
        parts["cmedians"].set_color("#1c2833"); parts["cmedians"].set_linewidth(1.4)

        ax.axhline(TARGET_THR, color="#1c2833", ls="--", lw=1.2, zorder=5)

        ax.set_ylim(-0.5, ytop + 1.6)
        for p, g in zip(pos, order):
            ax.text(p, ytop + 0.8, f"{pct_tgt[g]:.1f}%", ha="center",
                    va="center", fontsize=10, fontweight="bold",
                    color=col(g, toxic))
            ax.text(p, -1.7, f"HSPC/Mye\n{xsafe[g]:.0f}%", ha="center",
                    va="center", fontsize=6.6,
                    color="#c0392b" if toxic[g] else "#1e8449")

        ax.set_xticks(pos)
        xt = ax.set_xticklabels([GLABEL[g] for g in order], fontsize=8.5)
        for lab, g in zip(xt, order):
            if g == "CD96":
                lab.set_color(PALETTE["CD96"]); lab.set_fontweight("bold")
        ax.set_ylabel("bulk expression in AML  (TARGET ≥80% blasts, log2 CPM)")
        ax.set_xlim(-0.7, len(order) - 0.3)
        ax.margins(y=0)
        ax.set_title(f"n = {n}", fontsize=9, loc="left")  # cohort n only; rest in caption
        clean_ax(ax, "y")
        ax.tick_params(axis="x", length=0)
        fig.subplots_adjust(bottom=0.20)
        save_fig(fig, HERE, out_stem)
    return pct_tgt


def main():
    pac = pd.read_csv(PAC, index_col=0)
    xsafe = pac[["0_HSPC", "Myeloid_Pro"]].max(axis=1)          # worst-case
    L_all = bulk_matrix(exclude_mll=False)
    L_non = bulk_matrix(exclude_mll=True)
    order = list(L_all.median().sort_values(ascending=False).index)  # shared
    ytop = max(float(L_all.max().max()), float(L_non.max().max()))   # shared
    print(f"  TARGET_THR(log2CPM)={TARGET_THR}; TOX_THR={TOX_THR}%")
    pa = build(L_all, order, xsafe, ytop, "Figure7_window_violins",
               "all primary AML")
    pn = build(L_non, order, xsafe, ytop, "Figure7_window_violins_nonMLL",
               "KMT2A-rearranged (MLL) excluded")
    print("\n  CD96 targetable:  all={:.1f}%   non-MLL={:.1f}%   (+{:.1f} pts)"
          .format(pa["CD96"], pn["CD96"], pn["CD96"] - pa["CD96"]))
    print("Done -> Figure7_window_violins[_nonMLL].{pdf,png,csv}")


if __name__ == "__main__":
    main()
