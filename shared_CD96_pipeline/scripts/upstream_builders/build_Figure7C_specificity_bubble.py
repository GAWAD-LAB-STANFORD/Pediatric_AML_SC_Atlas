#!/usr/bin/env python
# ====================================================================
# build_Figure7C_specificity_bubble.py  |  Upstream builder (regenerates source_data from raw inputs)
# reproducibility header (auto-generated; edit freely below this block)
# Inputs : (reads upstream object / raw matrix — see body)
# Outputs: (writes to ../figures or ../source_data — see body)
# Run    : python build_Figure7C_specificity_bubble.py   (or: bash run_all.sh)  |  Deps: see REPRODUCIBILITY_CHECKLIST.md
# Source data for this figure is in __SUBMISSION_PACKAGE/08_Source_Data/.
# ====================================================================

"""Figure 6C -- AML surface-target specificity across marrow cell-type
compartments (scRNA).  Nine clinical/clinical-stage AML immunotherapy targets
x eight compartments; colour = % cells with >=1 raw read.

The point of the panel is the SAFETY CONTRAST: CD96 (top row) is high on the
leukemic compartment and on T/NK, but LOW on HSPC + myeloid progenitors (the
two shaded "safety" columns) -- whereas CD33, CLL-1, FLT3 and IL1RAP light up
exactly those progenitor compartments.  That HSPC/myeloid sparing is the core
reason CD96 has a wider therapeutic window than the targets already in trials.

"% positive" = fraction of cells with CD96/... raw count > 0 (standard scRNA
convention; same as the rest of Figure 6).  Pure data from the 28-patient
atlas; no fabricated values.
"""
import os
import sys
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
import matplotlib.pyplot as plt
import matplotlib as mpl

HERE = os.path.dirname(os.path.abspath(__file__))
SUPP = os.path.join(HERE, "supplements")
sys.path.insert(0, SUPP)
from _fig7_style import PALETTE, RC, save_fig, SEQ_CMAP, CD96_CMAP

H5AD = os.path.join(HERE, "..", "ALSF_AML_plot", "H5AD",
                    "ALSF_AML_Combo_3500_with_PAC_new.h5ad")

GENES = ["CD96", "CD33", "IL3RA", "CLEC12A", "FLT3", "MSLN",
         "CD70", "ADGRG1", "IL1RAP"]
GLABEL = {
    "CD96": "CD96", "CD33": "CD33", "IL3RA": "CD123 (IL3RA)",
    "CLEC12A": "CLL-1 (CLEC12A)", "FLT3": "FLT3", "MSLN": "MSLN",
    "CD70": "CD70", "ADGRG1": "GPR56 (ADGRG1)", "IL1RAP": "IL1RAP",
}

CELLTYPE_GROUPING = {
    "AML": "AML (leukemic)", "AML-MKI67": "AML (leukemic)", "AML-PCNA": "AML (leukemic)",
    "0_HSPC": "HSPC",
    "Myeloid_Pro": "Myeloid_Pro",
    "CTL": "T cells", "Naïve_CD4T": "T cells", "Naïve_CD8T": "T cells",
    "Activated_CD4T": "T cells", "AML-CTL": "T cells", "AML-CD4T": "T cells",
    "AML-Naïve_CD8T": "T cells",
    "NK": "NK", "AML-NK": "NK",
    "CD34+ProB": "B-lineage", "ProB": "B-lineage", "PreB": "B-lineage",
    "CD20+B": "B-lineage", "PlasmaB": "B-lineage", "AML-B": "B-lineage",
    "CD14_Monocyte": "Mature myeloid", "CD16_Monocyte": "Mature myeloid",
    "Macrophage": "Mature myeloid", "mDC": "Mature myeloid", "pDC": "Mature myeloid",
    "AML-CD14": "Mature myeloid", "AML-CD1C": "Mature myeloid",
    "Erythrocytes": "Erythroid", "AML-Ery": "Erythroid",
}
GROUP_ORDER = ["AML (leukemic)", "HSPC", "Myeloid_Pro",
               "Mature myeloid", "B-lineage", "T cells", "NK", "Erythroid"]
SAFETY_COLS = {"HSPC", "Myeloid_Pro"}
VMAX = 70.0


def _decode(arr):
    return [x.decode() if isinstance(x, bytes) else x for x in arr]


def main():
    print("  loading h5ad ...")
    with h5py.File(H5AD, "r") as f:
        ctype = np.array(_decode(f["obs"]["Cell Type"]["categories"][:])
                         )[f["obs"]["Cell Type"]["codes"][:].astype(int)]
        rx = f["raw"]["X"]
        X = sp.csr_matrix((rx["data"][:], rx["indices"][:], rx["indptr"][:]),
                          shape=tuple(rx.attrs["shape"])).tocsc()
        rk = f["raw"]["var"].attrs.get("_index", b"_index")
        rk = rk.decode() if isinstance(rk, bytes) else rk
        symbols = _decode(f["raw"]["var"][rk][:])
        pos = {}
        for g in GENES:
            if g not in symbols:
                raise KeyError(f"{g} not in h5ad raw var")
            raw = np.asarray(X[:, symbols.index(g)].todense(), np.float32).ravel()
            pos[g] = raw > 0

    group = np.array([CELLTYPE_GROUPING.get(c, "Other") for c in ctype])
    comps = [g for g in GROUP_ORDER if (group == g).any()]
    ncix = {c: int((group == c).sum()) for c in comps}

    M = np.zeros((len(GENES), len(comps)))
    for i, g in enumerate(GENES):
        for j, c in enumerate(comps):
            m = group == c
            M[i, j] = 100.0 * pos[g][m].mean()

    df = pd.DataFrame(M, index=GENES, columns=comps).round(3)
    df.to_csv(os.path.join(HERE, "Figure7C_target_by_compartment.csv"))
    print("\n  --- % positive (target x compartment) ---")
    print(df.round(1).to_string())

    cmap = SEQ_CMAP
    norm = mpl.colors.Normalize(vmin=0, vmax=VMAX)
    ng, nc = len(GENES), len(comps)

    def s_of(v):              # dot area scales with % positive (capped so dots don't overlap)
        return v * 6.0 + 2.0

    with mpl.rc_context(RC):
        fig, ax = plt.subplots(figsize=(10.8, 6.6))
        # safety columns + CD96 row shaded behind the dots
        for j, c in enumerate(comps):
            if c in SAFETY_COLS:
                ax.axvspan(j - 0.5, j + 0.5, color="#eaf5ee", zorder=0)
        ax.axhspan(-0.5, 0.5, color="#fbecea", zorder=0)   # CD96 row
        ax.set_axisbelow(True)
        ax.grid(color="#e9eef4", lw=0.6, zorder=1)

        # other targets -> navy magnitude colormap;  CD96 row (i==0) -> red
        xs, ys, ss, cc = [], [], [], []
        rx, ry, rs, rc = [], [], [], []
        for i in range(ng):
            for j in range(nc):
                if i == 0:
                    rx.append(j); ry.append(i); rs.append(s_of(M[i, j])); rc.append(M[i, j])
                else:
                    xs.append(j); ys.append(i); ss.append(s_of(M[i, j])); cc.append(M[i, j])
        ax.scatter(xs, ys, s=ss, c=cc, cmap=cmap, norm=norm,
                   edgecolor="#3b5870", linewidth=0.4, zorder=3)
        ax.scatter(rx, ry, s=rs, c=rc, cmap=CD96_CMAP, norm=norm,
                   edgecolor="#7d2018", linewidth=0.5, zorder=4)
        ax.set_xlim(-0.6, nc - 0.4); ax.set_ylim(-0.6, ng - 0.4)
        ax.invert_yaxis()      # CD96 (row 0) at the top

        ax.set_xticks(range(nc))
        xt = ax.set_xticklabels([f"{c}\n(n={ncix[c]:,})" for c in comps],
                                rotation=30, ha="right", fontsize=8.5)
        for lab, c in zip(xt, comps):
            if c in SAFETY_COLS:
                lab.set_color(PALETTE["favorable"]); lab.set_fontweight("bold")
            elif c == "AML (leukemic)":
                lab.set_color(PALETTE["CD96"]); lab.set_fontweight("bold")
        ax.set_yticks(range(ng))
        yt = ax.set_yticklabels([GLABEL[g] for g in GENES], fontsize=9.5)
        yt[0].set_color(PALETTE["CD96"]); yt[0].set_fontweight("bold")
        ax.tick_params(length=0)
        for s in ax.spines.values():
            s.set_visible(False)

        # single combined legend: dot grows AND darkens with % positive
        for v in (20, 50, 80):
            ax.scatter([], [], s=s_of(v), color=cmap(norm(v)),
                       edgecolor="#3b5870", linewidth=0.4, label=f"{v}%")
        ax.legend(title="% cells positive\n(≥1 raw read)", loc="center left",
                  bbox_to_anchor=(1.01, 0.5), labelspacing=2.0, frameon=False,
                  fontsize=8.5, title_fontsize=8.5, borderpad=1.2, handletextpad=1.2)

        ax.set_title("", pad=10)   # title/commentary in caption (authors)
        fig.tight_layout()
        save_fig(fig, HERE, "Figure7C_specificity_bubble")
    print("\nDone -> Figure7C_specificity_bubble.{pdf,png} + CSV")


if __name__ == "__main__":
    main()
