#!/usr/bin/env python
# ====================================================================
# build_Figure7E_clinical_headtohead.py  |  Upstream builder (regenerates source_data from raw inputs)
# reproducibility header (auto-generated; edit freely below this block)
# Inputs : (reads upstream object / raw matrix — see body)
# Outputs: (writes to ../figures or ../source_data — see body)
# Run    : python build_Figure7E_clinical_headtohead.py   (or: bash run_all.sh)  |  Deps: see REPRODUCIBILITY_CHECKLIST.md
# Source data for this figure is in __SUBMISSION_PACKAGE/08_Source_Data/.
# ====================================================================

"""Figure 7E -- CD96 vs clinical-stage AML surface targets vs other discovery
candidates, on a single safety-vs-coverage scatter.

Each gene becomes one point:
  X = max( %HSPC ,  %Myeloid_Pro )  expressing in scRNA   (safety leak;
      lower = better)
  Y = % of leukemic cells (AML_1-27) expressing in scRNA   (AML coverage;
      higher = better)
  size  = % TARGET patients with above-baseline bulk expression
          (population reach)
  color = group (CD96 / clinical-stage AML target / other discovery candidate)

The shaded green band on the left is the safety gate (HSPC<3%, Myeloid<=5%).
The ideal target lives in the upper-LEFT.

Reads the cached ungated per-gene table built by
supplements/build_CD96_stress_test.py.  If the cache is missing, run that
script first.
"""
import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from adjustText import adjust_text

HERE = os.path.dirname(os.path.abspath(__file__))
SUPP = os.path.join(HERE, "supplements")
sys.path.insert(0, SUPP)
from _fig7_style import PALETTE, RC, clean_ax, save_fig

CACHE = os.path.join(SUPP, "_full_ungated_per_gene_table.csv")

# ---- groups ----
CLINICAL_TARGETS = {
    "CD33":     "gemtuzumab ozogamicin (FDA, AML)",
    "IL3RA":    "tagraxofusp (CD123), AML trials",
    "CLEC12A":  "CLL-1 / CD371, multiple CAR-T trials",
    "FLT3":     "midostaurin/gilteritinib + CAR-T",
    "CD70":     "cusatuzumab (ARGX-110), AML trials",
    "ADGRG1":   "GPR56 (LSC marker, preclinical)",
    "IL1RAP":   "nidanilimab / multiple CAR-T",
    "MSLN":     "mesothelin (preclinical AML)",
}
# discovery candidates = the SAME top-15 composite genes shown in panel A
# (Fig 6A), minus CD96 and the clinical targets -> 6A and 6C stay consistent.
_RANKED = os.path.join(HERE, "Figure7A_composite_ranked.csv")
DISCOVERY = [g for g in pd.read_csv(_RANKED).head(15)["gene"].tolist()
             if g != "CD96"]   # (none of the clinical targets rank in the top 15)

# safety gate values (the only "knobs" -- kept transparent)
GATE_HSPC = 3.0
GATE_MYE = 5.0


def _grp(gene):
    if gene == "CD96":
        return "CD96 (this paper)"
    if gene in CLINICAL_TARGETS:
        return "Clinical-stage AML target"
    return "Other discovery candidate"


def main():
    if not os.path.exists(CACHE):
        sys.exit(f"missing cache: {CACHE}\n"
                 "  run  python Figure7/supplements/build_CD96_stress_test.py  first")
    full = pd.read_csv(CACHE)

    keep = ["CD96"] + list(CLINICAL_TARGETS) + DISCOVERY
    sub = full[full["gene"].isin(keep)].copy()
    missing = set(keep) - set(sub["gene"])
    if missing:
        print(f"  WARNING genes not in cached table: {missing}")
    sub["leak"] = sub[["scRNA_HSPC_pct", "scRNA_Myeloid_pct"]].max(axis=1)
    sub["group"] = sub["gene"].apply(_grp)
    sub = sub.sort_values("scRNA_AML_pct", ascending=False).reset_index(drop=True)

    # save the source data for the panel
    out_cols = ["gene", "group", "scRNA_AML_pct", "scRNA_HSPC_pct",
                "scRNA_Myeloid_pct", "leak", "TARGET_pct_above_Tbaseline",
                "scRNA_PPAC_pct", "scRNA_n_pts_above20pct"]
    sub[out_cols].to_csv(os.path.join(HERE, "Figure7E_clinical_headtohead.csv"),
                         index=False)
    print(sub[out_cols].to_string(index=False))

    # cross-panel scheme: crimson = CD96 (hero), navy = clinical comparators,
    # muted grey = other discovery candidates  (matches panels C / D / F)
    GCOL = {
        "CD96 (this paper)":         PALETTE["CD96"],
        "Clinical-stage AML target": "#2471a3",
        "Other discovery candidate": "#9fb3c8",
    }

    with mpl.rc_context(RC):
        fig, ax = plt.subplots(figsize=(8.8, 6.4))

        # safety-gate shading: HSPC < 3% region (green)
        ax.axvspan(-2, GATE_HSPC, color="#eafaf1", alpha=0.55, zorder=0,
                    label=f"<{GATE_HSPC}% HSPC  +  <={GATE_MYE}% Myeloid_Pro")
        # myeloid threshold marker (no text annotation -- legend covers it)
        ax.axvline(GATE_MYE, color="#bdc3c7", ls=":", lw=0.9, zorder=0)

        # point sizes: TARGET breadth scaled to area  (smaller than before)
        sizes = 32 + 6 * sub["TARGET_pct_above_Tbaseline"].values

        for grp in ["Other discovery candidate", "Clinical-stage AML target",
                    "CD96 (this paper)"]:
            sg = sub[sub["group"] == grp]
            ax.scatter(sg["leak"], sg["scRNA_AML_pct"],
                        s=sizes[sg.index], c=GCOL[grp],
                        edgecolor="black",
                        linewidths=1.2 if grp == "CD96 (this paper)" else 0.6,
                        label=grp, alpha=0.95, zorder=3)

        # label every gene including CD96; let adjustText resolve overlaps.
        texts = []
        for _, r in sub.iterrows():
            if r["gene"] == "CD96":
                color = PALETTE["CD96"]; fw = "bold"
            elif r["gene"] in CLINICAL_TARGETS:
                color = "#1a5276"; fw = "bold"
            else:
                color = "#5d6d7e"; fw = "normal"
            texts.append(
                ax.text(r["leak"], r["scRNA_AML_pct"], r["gene"],
                         fontsize=8.5, color=color, fontweight=fw, zorder=6,
                         ha="center", va="center",
                         bbox=dict(boxstyle="round,pad=0.18",
                                    facecolor="white", edgecolor="none",
                                    alpha=0.55))
            )
        # push labels well clear of the markers; thin leader lines connect them
        adjust_text(
            texts, ax=ax,
            expand=(1.9, 2.3),
            force_text=(1.6, 2.0),
            force_static=(1.4, 1.8),
            force_explode=(0.6, 0.9),
            arrowprops=dict(arrowstyle="-", color="#888", lw=0.5, alpha=0.8),
            min_arrow_len=10.0,
        )

        ax.set_xlabel("% normal HSPCs or myeloid progenitors expressing target  "
                      "(scRNA;  max of the two)")
        ax.set_ylabel("% leukemic cells (AML_1-27) expressing target  (scRNA)")
        ax.set_title("")   # title/commentary in caption (authors)
        # size legend (TARGET breadth)
        size_handles = []
        for v in [10, 25, 40]:
            size_handles.append(
                plt.scatter([], [], s=40 + 9 * v, c="white",
                             edgecolor="black", linewidths=0.7,
                             label=f"{v}% TARGET"))
        legend1 = ax.legend(loc="upper right", fontsize=8, frameon=False,
                             title="group", title_fontsize=8)
        ax.add_artist(legend1)
        ax.legend(handles=size_handles, loc="lower right", fontsize=7,
                   frameon=False, title="TARGET reach", title_fontsize=8,
                   labelspacing=1.4, borderpad=1.0)

        ax.set_xlim(-1.0, max(sub["leak"].max() + 3, 56))
        ax.set_ylim(0, max(sub["scRNA_AML_pct"].max() + 12, 55))
        clean_ax(ax, "both")
        save_fig(fig, HERE, "Figure7E_clinical_headtohead")

    print("\nDone -> Figure7E_clinical_headtohead.{pdf,png,csv}")


if __name__ == "__main__":
    main()
