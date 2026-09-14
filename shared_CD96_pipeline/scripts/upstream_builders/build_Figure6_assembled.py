#!/usr/bin/env python
# ====================================================================
# build_Figure6_assembled.py  |  Upstream builder (regenerates source_data from raw inputs)
# reproducibility header (auto-generated; edit freely below this block)
# Inputs : (reads upstream object / raw matrix — see body)
# Outputs: (writes to ../figures or ../source_data — see body)
# Run    : python build_Figure6_assembled.py   (or: bash run_all.sh)  |  Deps: see REPRODUCIBILITY_CHECKLIST.md
# Source data for this figure is in __SUBMISSION_PACKAGE/08_Source_Data/.
# ====================================================================

"""Assemble the six discovery panels into FIGURE 6 -- "how CD96 was nominated
as the lead surface immunotherapy target in pediatric AML".

(Formerly Figure 7; renumbered when CD96 became the paper's central story and
the old Figure 6 was retired.  The companion FIGURE 7 is the body-wide
normal-tissue toxicity atlas, built by build_Figure7_toxicity_atlas.py.)

Pure layout step: reads the high-res panel PNGs and tiles them into a 3x2
montage.  It recomputes nothing, so the composite matches the panels exactly.
Re-run a panel builder then re-run this to refresh.  Panel letters (6A..6F)
are baked into each panel's own title.

Layout (row-major):
    A  discovery funnel        |  B  CD96 scRNA atlas (UMAP)
    C  specificity bubble       |  D  cytogenetic + PAC heatmap
    E  clinical head-to-head    |  F  TARGET bulk + blast% validation
       (= therapeutic window /     (cross-cohort coverage validation)
        selectivity vs clinical
        AML targets incl. MSLN)
"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

HERE = os.path.dirname(os.path.abspath(__file__))

# row-major panel grid -- panel PNG file stems (intermediate artifacts keep
# their Figure7* filenames; their titles now read 6A..6F)
# A funnel, B UMAP, C head-to-head (competitive overview), D specificity,
# E coverage, F validation  -- E-content (head-to-head) precedes the two
# dot-plot drill-downs (D specificity, E coverage), which sit on a diagonal so
# the two dot plots aren't adjacent.
# Final 6-panel main figure:
#   A funnel (nomination)        | B CD96 UMAP (marks the leukemia)
#   C specificity bubble          | D cyto/PAC coverage
#   E cross-cohort validation     | F window-violins (all AML)
# (non-MLL window violins -> supplement; head-to-head -> Figure 7)
GRID = [
    ["Figure7A_discovery_funnel",        "Figure7B_CD96_UMAP"],
    ["Figure7E_clinical_headtohead",     "Figure7D_cyto_PAC_heatmap"],
    ["Figure6F_validation",              "Figure7_window_violins"],
]

USABLE_W = 15.75
SIDE     = 0.18
COL_GAP  = 0.22
ROW_GAP  = 0.42
TOP      = 0.62
BOT      = 0.16
DPI      = 300


def main():
    imgs, asp = {}, {}
    for row in GRID:
        for stem in row:
            path = os.path.join(HERE, stem + ".png")
            if not os.path.exists(path):
                raise FileNotFoundError(
                    f"missing panel PNG: {path}\n"
                    f"  build it first with build_{stem.split('_')[0]}*.py")
            im = mpimg.imread(path)
            h, w = im.shape[:2]
            imgs[stem] = im
            asp[stem] = w / h
            print(f"  {stem:38s} {w}x{h}  (aspect {w/h:.3f})")

    # 2-col rows fill the usable width; a singleton (odd) row is capped to the
    # mean 2-col row height and centred so it doesn't dominate.
    multi_h = [(USABLE_W - COL_GAP * (len(r) - 1)) / sum(asp[s] for s in r)
               for r in GRID if len(r) > 1]
    cap = (sum(multi_h) / len(multi_h)) if multi_h else USABLE_W
    row_h = []
    for row in GRID:
        if len(row) == 1:
            row_h.append(min(USABLE_W / asp[row[0]], cap))
        else:
            gaps = COL_GAP * (len(row) - 1)
            row_h.append((USABLE_W - gaps) / sum(asp[s] for s in row))

    fig_w = USABLE_W + 2 * SIDE
    fig_h = TOP + BOT + sum(row_h) + ROW_GAP * (len(GRID) - 1)
    print(f"\n  composite figure: {fig_w:.2f} x {fig_h:.2f} in  "
          f"-> {int(fig_w*DPI)}x{int(fig_h*DPI)} px @ {DPI} dpi")

    fig = plt.figure(figsize=(fig_w, fig_h), facecolor="white")
    letters = "ABCDEFGHIJ"
    pi = 0
    y_top = fig_h - TOP
    for ri, row in enumerate(GRID):
        h = row_h[ri]
        y_bot = y_top - h
        # centre a singleton row; left-align multi-panel rows
        if len(row) == 1:
            x = SIDE + (USABLE_W - h * asp[row[0]]) / 2.0
        else:
            x = SIDE
        for stem in row:
            w = h * asp[stem]
            ax = fig.add_axes([x / fig_w, y_bot / fig_h, w / fig_w, h / fig_h])
            ax.imshow(imgs[stem], aspect="auto", interpolation="lanczos")
            ax.set_xticks([]); ax.set_yticks([])
            for s in ax.spines.values():
                s.set_edgecolor("#d0d0d0"); s.set_linewidth(0.7)
            # bold panel letter just above the tile's top-left corner
            fig.text(x / fig_w - 0.004, (y_bot + h) / fig_h + 0.004,
                     letters[pi], ha="left", va="bottom",
                     fontsize=21, fontweight="bold", color="#1c2833")
            pi += 1
            x += w + COL_GAP
        y_top = y_bot - ROW_GAP

    # (no figure title in the image -- title/caption added by the authors)

    for ext in ("png", "pdf"):
        out = os.path.join(HERE, f"Figure6_assembled.{ext}")
        fig.savefig(out, dpi=DPI, facecolor="white", bbox_inches=None, pad_inches=0)
        print(f"  wrote {out}")
    plt.close(fig)
    print("\nDone -> Figure6_assembled.{png,pdf}")


if __name__ == "__main__":
    main()
