#!/usr/bin/env python
"""Composite Figure 4 (corrected scheme): T/NK landscape re-derived from the count matrix.
A UMAP of the 12 T/NK subsets; B per-sample composition (AML vs healthy BM); C defining-marker dot plot.
Panel D (finer-subset TARGET survival) is added after the CIBERSORTx run. Letters only."""
import os
from PIL import Image, ImageDraw, ImageFont
P   = "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/04_Main_Figures/Figure_4_panels"
OUT = "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/04_Main_Figures"
ROWS = [
  [("Figure_4A_tnk_umap.png","A"), ("Figure_4D_tnk_aml_vs_hbm.png","B"), ("Figure_4G_tnk_prognosis_dot.png","C")],  # same UMAP (subsets / AML-vs-HBM) + prognosis
  [("Figure_4B_tnk_composition.png","D")],          # per-sample composition
  [("Figure_4C_tnk_markers.png","E"), ("Figure_4F_tnk_clonality_perpatient.png","F")],  # markers + TCR Shannon clonality
  [("Figure_4_cd8_functional.png","G")],   # CD8 functional state: GZMB-CD8 = spent terminal effector
]  # GSEA heatmap + un-integrated AML-vs-HBM overlay + CD8-vs-clonality scatter moved to supplementary figures
W, GAP, HEAD, TOP, MARGIN = 2200, 44, 74, 24, 28
MAXH = 0.62 * W
def font(sz):
    for p in ["/System/Library/Fonts/Supplemental/Arial Bold.ttf", "/System/Library/Fonts/Helvetica.ttc"]:
        try: return ImageFont.truetype(p, sz)
        except Exception: pass
    return ImageFont.load_default()
fL = font(58)
rowlays = []
for row in ROWS:
    ims = [Image.open(os.path.join(P, f)).convert("RGB") for f, _ in row]
    ratios = [im.width/im.height for im in ims]; n = len(row); avail = W - GAP*(n-1)
    Hrow = min(avail/sum(ratios), MAXH); widths = [r*Hrow for r in ratios]
    tot = sum(widths) + GAP*(n-1); x = MARGIN + (W - tot)/2; placed = []
    for (f, lab), im, w in zip(row, ims, widths):
        placed.append((im, lab, int(round(x)), int(round(w)))); x += w + GAP
    rowlays.append((int(round(Hrow)), placed))
totalH = TOP + sum(HEAD + h for h, _ in rowlays) + GAP*len(rowlays)
canvas = Image.new("RGB", (W + 2*MARGIN, int(totalH)), "white"); d = ImageDraw.Draw(canvas); y = TOP
for Hrow, placed in rowlays:
    for im, lab, x, w in placed:
        if lab: d.text((x, y), lab, font=fL, fill="black")
        canvas.paste(im.resize((w, Hrow), Image.LANCZOS), (x, y + HEAD))
    y += HEAD + Hrow + GAP
canvas.save(os.path.join(OUT, "Figure_4__T_NK_states.png"))
print("wrote Figure_4__T_NK_states  %dx%d" % (canvas.width, canvas.height))
