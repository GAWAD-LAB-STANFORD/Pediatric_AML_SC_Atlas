#!/usr/bin/env python
# ====================================================================
# build_supplements.py  |  Upstream builder (regenerates source_data from raw inputs)
# reproducibility header (auto-generated; edit freely below this block)
# Inputs : (reads upstream object / raw matrix — see body)
# Outputs: (writes to ../figures or ../source_data — see body)
# Run    : python build_supplements.py   (or: bash run_all.sh)  |  Deps: see REPRODUCIBILITY_CHECKLIST.md
# Source data for this figure is in __SUBMISSION_PACKAGE/08_Source_Data/.
# ====================================================================

"""Assemble the supporting panels into a few SUPPLEMENTARY figures for the CD96
manuscript (the main story lives in Figure 6 + Figure 7).

  Suppl. Fig S6  -- CD96 specificity / clinical comparison (supports Fig 6):
                    A cell-type specificity bubble   B clinical head-to-head
  Suppl. Fig S8  -- CD96 combination strategy (supports Fig 6/7):
                    A coverage-vs-toxicity window     B per-partner bars

Standalone supplements (one panel each, not re-assembled here):
  Suppl. Fig S7  = FigureS_HSPC_hierarchy.png  (myeloid/lymphoid differentiation)
  Suppl. Fig S9  = FigureS_CD96_correlations.png  (CD96 co-expression)

Pure layout step: tiles existing high-res panel PNGs.  Re-run a panel builder
then re-run this.  No recomputation.
"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

HERE = os.path.dirname(os.path.abspath(__file__))
USABLE_W, SIDE, ROW_GAP, TOP, BOT, DPI = 11.0, 0.16, 0.42, 0.55, 0.14, 300

SUPPS = {
    "FigureS6_specificity": dict(
        title="Suppl. Figure S6   CD96 cell-type specificity and clinical "
              "head-to-head (supports Figure 6)",
        rows=[["Figure7C_specificity_bubble"], ["Figure7E_clinical_headtohead"]]),
    "FigureS7_combinations": dict(
        title="Suppl. Figure S7   CD96 combination strategy: coverage vs "
              "stem/progenitor toxicity (supports Figures 6–8)",
        rows=[["Figure_CD96_combination_window"], ["FigureS_CD96_combinations"]]),
}


def assemble(stem, title, rows):
    imgs, asp = {}, {}
    for row in rows:
        for s in row:
            p = os.path.join(HERE, s + ".png")
            if not os.path.exists(p):
                raise FileNotFoundError(p)
            im = mpimg.imread(p); imgs[s] = im
            asp[s] = im.shape[1] / im.shape[0]
    row_h = [(USABLE_W) / sum(asp[s] for s in row) for row in rows]   # 1 col/row
    fig_w = USABLE_W + 2 * SIDE
    fig_h = TOP + BOT + sum(row_h) + ROW_GAP * (len(rows) - 1)
    fig = plt.figure(figsize=(fig_w, fig_h), facecolor="white")
    letters = "ABCDEFGH"
    y_top, pi = fig_h - TOP, 0
    for ri, row in enumerate(rows):
        h = row_h[ri]; y_bot = y_top - h; x = SIDE
        for s in row:
            w = h * asp[s]
            ax = fig.add_axes([x / fig_w, y_bot / fig_h, w / fig_w, h / fig_h])
            ax.imshow(imgs[s], aspect="auto", interpolation="lanczos")
            ax.set_xticks([]); ax.set_yticks([])
            for sp in ax.spines.values():
                sp.set_edgecolor("#d0d0d0"); sp.set_linewidth(0.7)
            fig.text(x / fig_w - 0.005, (y_bot + h) / fig_h + 0.004, letters[pi],
                     ha="left", va="bottom", fontsize=19, fontweight="bold",
                     color="#1c2833")
            pi += 1; x += w
        y_top = y_bot - ROW_GAP
    fig.text(0.014, 1 - 0.28 / fig_h, title, ha="left", va="center",
             fontsize=12, fontweight="bold", color="#1c2833")
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(HERE, f"{stem}.{ext}"), dpi=DPI,
                    facecolor="white", pad_inches=0)
    plt.close(fig)
    print(f"  wrote {stem}.{{png,pdf}}")


def main():
    for stem, spec in SUPPS.items():
        assemble(stem, spec["title"], spec["rows"])
    print("Done -> supplementary figures")


if __name__ == "__main__":
    main()
