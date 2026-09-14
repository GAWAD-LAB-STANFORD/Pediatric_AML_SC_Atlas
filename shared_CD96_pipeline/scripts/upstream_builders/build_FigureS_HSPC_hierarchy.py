#!/usr/bin/env python
# ====================================================================
# build_FigureS_HSPC_hierarchy.py  |  Upstream builder (regenerates source_data from raw inputs)
# reproducibility header (auto-generated; edit freely below this block)
# Inputs : (reads upstream object / raw matrix — see body)
# Outputs: (writes to ../figures or ../source_data — see body)
# Run    : python build_FigureS_HSPC_hierarchy.py   (or: bash run_all.sh)  |  Deps: see REPRODUCIBILITY_CHECKLIST.md
# Source data for this figure is in __SUBMISSION_PACKAGE/08_Source_Data/.
# ====================================================================

"""CD96 vs comparators across normal haematopoietic differentiation, shown as a
DOT PLOT with the two lineage arms as SEPARATE column groups (no connecting
lines -- MEP and CLP are on different branches, so no line should join them):

  MYELOID :  HSC  MPP  CMP  GMP  MEP        (GMP & MEP both branch from CMP)
  LYMPHOID:  HSC  MPP  CLP

Dot size + colour = % cells positive (>=1 raw read).  CD96 (top row, red) stays
near-zero across the whole MYELOID stem/progenitor arm -- the marrow-regenerating
cells -- rising only at the lymphoid CLP; CD33 / CLL-1 / FLT3 / GPR56 light up the
myeloid progenitors (the basis for their myelosuppression).

Source: CZ CELLxGENE Census global cell-type table.  No fabricated values.
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
from _fig7_style import PALETTE, RC, save_fig, SEQ_CMAP, CD96_CMAP, TARGET_COLORS

GLOBAL = os.path.join(HERE, "Figure8_census_normal_celltype_global.csv")
GENES = ["CD96", "CD33", "IL3RA", "CLEC12A", "FLT3", "MSLN", "CD70", "ADGRG1",
         "IL1RAP"]
GLABEL = {"CD96": "CD96", "CD33": "CD33", "IL3RA": "CD123 (IL3RA)",
          "CLEC12A": "CLL-1 (CLEC12A)", "FLT3": "FLT3", "MSLN": "MSLN",
          "CD70": "CD70", "ADGRG1": "GPR56 (ADGRG1)", "IL1RAP": "IL1RAP"}
# (census cell type, short label, x position) -- gap between MEP(4) and CLP(6)
STAGES = [("hematopoietic stem cell", "HSC", 0, "M"),
          ("hematopoietic multipotent progenitor cell", "MPP", 1, "M"),
          ("common myeloid progenitor", "CMP", 2, "M"),
          ("granulocyte monocyte progenitor cell", "GMP", 3, "M"),
          ("megakaryocyte-erythroid progenitor cell", "MEP", 4, "M"),
          ("common lymphoid progenitor", "CLP", 6, "L")]
VMAX = 60.0


def s_of(v):
    return v * 5.0 + 4.0


def main():
    g = pd.read_csv(GLOBAL).set_index("cell_type")
    miss = [c for c, _, _, _ in STAGES if c not in g.index]
    if miss:
        raise SystemExit(f"missing cell types: {miss}")
    xs = [x for _, _, x, _ in STAGES]
    M = np.array([[float(g.loc[c, f"{gn}_pct"]) for c, _, _, _ in STAGES]
                  for gn in GENES])
    norm = mpl.colors.Normalize(vmin=0, vmax=VMAX)

    with mpl.rc_context(RC):
        fig, ax = plt.subplots(figsize=(7.6, 4.6))
        ax.axhspan(-0.5, 0.5, color="#fbecea", zorder=0)          # CD96 row
        ax.axvline(5.0, color="#bbbbbb", lw=1.0, ls="--", zorder=1)  # arm split
        ax.set_axisbelow(True); ax.grid(color="#eef2f6", lw=0.5, zorder=1)
        xsB, ysB, ssB, ccB, xsR, ysR, ssR, ccR = ([] for _ in range(8))
        for i in range(len(GENES)):
            for j, x in enumerate(xs):
                v = M[i, j]
                if i == 0:
                    xsR.append(x); ysR.append(i); ssR.append(s_of(v)); ccR.append(v)
                else:
                    xsB.append(x); ysB.append(i); ssB.append(s_of(v)); ccB.append(v)
        ax.scatter(xsB, ysB, s=ssB, c=ccB, cmap=SEQ_CMAP, norm=norm,
                   edgecolor="#3b5870", linewidth=0.3, zorder=3)
        ax.scatter(xsR, ysR, s=ssR, c=ccR, cmap=CD96_CMAP, norm=norm,
                   edgecolor="#7d2018", linewidth=0.45, zorder=4)
        ax.set_xlim(-0.6, 6.6); ax.set_ylim(-0.6, len(GENES) - 0.4)
        ax.invert_yaxis()
        ax.set_xticks(xs)
        ax.set_xticklabels([f"{lab}\n({int(g.loc[c,'n_cells'])//1000}k)"
                            for c, lab, _, _ in STAGES], fontsize=8)
        ax.set_yticks(range(len(GENES)))
        yt = ax.set_yticklabels([GLABEL[x] for x in GENES], fontsize=9)
        yt[0].set_color(PALETTE["CD96"]); yt[0].set_fontweight("bold")
        ax.tick_params(length=0)
        for s in ax.spines.values():
            s.set_visible(False)
        # arm group labels
        ax.text(2.0, -0.95, "MYELOID arm", ha="center", fontsize=9,
                fontweight="bold", color="#1c2833")
        ax.text(6.0, -0.95, "LYMPHOID", ha="center", fontsize=9,
                fontweight="bold", color="#1c2833")
        # size legend
        for v in (10, 30, 50):
            ax.scatter([], [], s=s_of(v), color=SEQ_CMAP(norm(v)),
                       edgecolor="#3b5870", linewidth=0.4, label=f"{v}%")
        ax.legend(title="% positive", loc="center left", bbox_to_anchor=(1.01, 0.5),
                  frameon=False, fontsize=8, title_fontsize=8, labelspacing=1.3)
        ax.set_title("")   # title/commentary in caption (authors)
        fig.subplots_adjust(left=0.18, right=0.86, top=0.84, bottom=0.16)
        save_fig(fig, HERE, "FigureS_HSPC_hierarchy")
    print("  CD96 across stages:", [round(M[0, j], 1) for j in range(len(xs))])
    print("Done -> FigureS_HSPC_hierarchy.{pdf,png}")


if __name__ == "__main__":
    main()
