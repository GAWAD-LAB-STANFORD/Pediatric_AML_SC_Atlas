#!/usr/bin/env python
# ====================================================================
# build_Figure7_assembled.py  |  Upstream builder (regenerates source_data from raw inputs)
# reproducibility header (auto-generated; edit freely below this block)
# Inputs : (reads upstream object / raw matrix — see body)
# Outputs: (writes to ../figures or ../source_data — see body)
# Run    : python build_Figure7_assembled.py   (or: bash run_all.sh)  |  Deps: see REPRODUCIBILITY_CHECKLIST.md
# Source data for this figure is in __SUBMISSION_PACKAGE/08_Source_Data/.
# ====================================================================

"""Assemble FIGURE 7 — CD96 spares the normal blood-forming system and has the
widest therapeutic window (4 panels).

  A  differentiation dot plot   CD96 across HSC->myeloid/lymphoid progenitors
  B  cell-type specificity      9 targets x marrow cell-type compartments
  C  combination window         CD96 + partner: coverage vs stem/progenitor
                                toxicity (within- and across-patient)
  D  combination bars           per-partner coverage, stacked from the CD96 baseline

Pure layout step: tiles pre-rendered panel PNGs.  Re-run the panel builders, then
this.  No recomputation.
"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

HERE = os.path.dirname(os.path.abspath(__file__))
GRID = [
    ["FigureS_HSPC_hierarchy", "Figure7C_specificity_bubble"],     # A diff, B specificity
    ["Figure_CD96_combination_window"],                            # C (wide)
    ["FigureS_CD96_combinations"],                                 # D (wide)
]
USABLE_W, SIDE, COL_GAP, ROW_GAP, TOP, BOT, DPI = 15.0, 0.18, 0.30, 0.45, 0.60, 0.16, 300


def main():
    imgs, asp = {}, {}
    for row in GRID:
        for s in row:
            p = os.path.join(HERE, s + ".png")
            if not os.path.exists(p):
                raise FileNotFoundError(p)
            im = mpimg.imread(p); imgs[s] = im; asp[s] = im.shape[1] / im.shape[0]
    row_h = [(USABLE_W - COL_GAP * (len(r) - 1)) / sum(asp[s] for s in r)
             for r in GRID]
    fig_w = USABLE_W + 2 * SIDE
    fig_h = TOP + BOT + sum(row_h) + ROW_GAP * (len(GRID) - 1)
    fig = plt.figure(figsize=(fig_w, fig_h), facecolor="white")
    letters = "ABCDEFGH"
    y_top, pi = fig_h - TOP, 0
    for ri, row in enumerate(GRID):
        h = row_h[ri]; y_bot = y_top - h; x = SIDE
        for s in row:
            w = h * asp[s]
            ax = fig.add_axes([x / fig_w, y_bot / fig_h, w / fig_w, h / fig_h])
            ax.imshow(imgs[s], aspect="auto", interpolation="lanczos")
            ax.set_xticks([]); ax.set_yticks([])
            for sp in ax.spines.values():
                sp.set_edgecolor("#d0d0d0"); sp.set_linewidth(0.7)
            fig.text(x / fig_w - 0.004, (y_bot + h) / fig_h + 0.003, letters[pi],
                     ha="left", va="bottom", fontsize=20, fontweight="bold",
                     color="#1c2833")
            pi += 1; x += w + COL_GAP
        y_top = y_bot - ROW_GAP
    # (no figure title in the image -- title/caption added by the authors)
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(HERE, f"Figure7_assembled.{ext}"), dpi=DPI,
                    facecolor="white", pad_inches=0)
    plt.close(fig)
    print("Done -> Figure7_assembled.{png,pdf}")


if __name__ == "__main__":
    main()
