#!/usr/bin/env python
# ====================================================================
# build_Figure8_assembled.py  |  Upstream builder (regenerates source_data from raw inputs)
# reproducibility header (auto-generated; edit freely below this block)
# Inputs : (reads upstream object / raw matrix — see body)
# Outputs: (writes to ../figures or ../source_data — see body)
# Run    : python build_Figure8_assembled.py   (or: bash run_all.sh)  |  Deps: see REPRODUCIBILITY_CHECKLIST.md
# Source data for this figure is in __SUBMISSION_PACKAGE/08_Source_Data/.
# ====================================================================

"""Assemble FIGURE 8 — body-wide normal-tissue toxicity (4 panels).

Two-column layout:
  LEFT  (full height)  A  detailed organ-grouped cell-type toxicity map
                          (bone-marrow stem/progenitors + non-blood tissues)
  RIGHT (stacked)      B  organ atlas (9 targets x organs; blood removed)
                       C  off-tumor compartment fingerprint (CD96 = T/NK only)
                       D  CD96's only normal liabilities (urothelium, adipocytes)

Pure layout step: tiles pre-rendered panel PNGs (build_Figure7_organ_toxicity.py
for A, build_Figure8_panels.py for B/C/D).  No recomputation.
"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

HERE = os.path.dirname(os.path.abspath(__file__))
DPI = 300
A = "Figure7_organ_toxicity"                 # left (detailed map)
RIGHT = ["Figure8_panelB_organatlas", "Figure8_panelC_fingerprint",
         "Figure8_panelD_cd96liab"]


def load(stem):
    im = mpimg.imread(os.path.join(HERE, stem + ".png"))
    return im, im.shape[1] / im.shape[0]


def main():
    imA, aspA = load(A)
    rights = [load(s) for s in RIGHT]

    # geometry (inches): left column width LW, right column width RW
    FIGH = 15.0
    SIDE, GAP, TOP = 0.16, 0.30, 0.30
    avail_h = FIGH - TOP - 0.16
    LW = avail_h * aspA                       # left panel width from its aspect
    RW = 6.6                                   # right column width
    figw = SIDE + LW + GAP + RW + SIDE
    fig = plt.figure(figsize=(figw, FIGH), facecolor="white")

    # left: panel A spanning full height
    xA = SIDE
    axA = fig.add_axes([xA / figw, 0.16 / FIGH, LW / figw, avail_h / FIGH])
    axA.imshow(imA, aspect="auto", interpolation="lanczos")
    axA.set_xticks([]); axA.set_yticks([])
    for s in axA.spines.values():
        s.set_edgecolor("#d0d0d0"); s.set_linewidth(0.7)
    fig.text(xA / figw - 0.004, (0.16 + avail_h) / FIGH + 0.002, "A",
             ha="left", va="bottom", fontsize=20, fontweight="bold",
             color="#1c2833")

    # right: B/C/D stacked, each at natural aspect, top-aligned with small gaps
    rx = SIDE + LW + GAP
    rgap = 0.30
    # allocate heights proportional to 1/aspect so each fills RW
    hs = [RW / asp for _, asp in rights]
    total = sum(hs) + rgap * (len(rights) - 1)
    scale = min(1.0, avail_h / total)
    hs = [h * scale for h in hs]
    y = FIGH - TOP
    letters = "BCD"
    for (im, asp), h, L in zip(rights, hs, letters):
        w = h * asp
        y -= h
        ax = fig.add_axes([rx / figw, y / FIGH, w / figw, h / FIGH])
        ax.imshow(im, aspect="auto", interpolation="lanczos")
        ax.set_xticks([]); ax.set_yticks([])
        for s in ax.spines.values():
            s.set_edgecolor("#d0d0d0"); s.set_linewidth(0.7)
        fig.text(rx / figw - 0.004, (y + h) / FIGH + 0.002, L, ha="left",
                 va="bottom", fontsize=18, fontweight="bold", color="#1c2833")
        y -= rgap
    # (no figure title in the image -- title/caption added by the authors)
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(HERE, f"Figure8_assembled.{ext}"), dpi=DPI,
                    facecolor="white", pad_inches=0)
    plt.close(fig)
    print("Done -> Figure8_assembled.{png,pdf}")


if __name__ == "__main__":
    main()
