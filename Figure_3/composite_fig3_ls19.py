#!/usr/bin/env python
"""Composite Figure 3 (LS19 version): a least-differentiated, interferon/pDC-primed progenitor
leukemic state with poor outcome and high CD96 (targetable, sparing normal HSPC/pDC).
A selection volcano (LS19 highlighted); B LS19 with normal pDC/HSPC on the atlas UMAP;
C marker programs (progenitor / interferon / pDC / differentiation / surface); D whole-cohort EFS KM;
E LS19 fraction across cytogenetic subtypes. Letters only; descriptions in the figure legend."""
import os
from PIL import Image, ImageDraw, ImageFont
P   = "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/04_Main_Figures/Figure_3_panels"
P2  = "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/04_Main_Figures/Figure_2_panels"
OUT = "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/04_Main_Figures"
def find(name):
    for d in (P, P2):
        p = os.path.join(d, name)
        if os.path.exists(p): return p
    raise FileNotFoundError(name)
ROWS = [
  [("Figure_3_LS19_volcano.png","A"), ("Figure_3_LS19_umap.png","B"), ("Figure_3_LS19_ddrtree.png","C")],
  [("Figure_3_LS19_program.png","D")],
  [("Figure_3_LS19_surface.png","E"), ("Figure_3_LS19_pathbar.png","F")],
  [("KM_LS_19_wholecohort_EFS.png","G"), ("Figure_2_LS19_by_cyto.png","H")],
]
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
    ims = [Image.open(find(f)).convert("RGB") for f, _ in row]
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
canvas.save(os.path.join(OUT, "Figure_3__LS19_pDC_progenitor.png"))
print("wrote Figure_3__LS19_pDC_progenitor  %dx%d  (aspect %.2f)" % (canvas.width, canvas.height, canvas.height/canvas.width))
