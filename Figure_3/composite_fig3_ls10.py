#!/usr/bin/env python
"""Composite Figure 3 — the LS10 undifferentiated poor-proliferation state.
A UMAP (LS10 vs NS-proliferation comparators vs normal HSC); B pathway contrast
(LS10 = DNA replication/chromosome segregation, comparators = myeloid differentiation);
C whole-cohort EFS KM; D LS10 fraction across cytogenetic subtypes (adverse-enriched);
E candidate HSC-sparing surface markers. Letters only."""
import os
from PIL import Image, ImageDraw, ImageFont
P   = "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/04_Main_Figures/Figure_2_panels"
OUT = "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/04_Main_Figures"
ROWS = [
  [("Figure_3_state_selection_volcano.png","A"), ("UMAP_LS10_vs_prolifNS_HSC.png","B")],
  [("Figure_3_LS10_gsea.png","C"), ("KM_LS_10_wholecohort_EFS.png","D")],
  [("Figure_2_LS10_by_cyto.png","E"), ("Figure_2_LS10_surface.png","F")],
  [("Figure_3_quiescence_program.png","G")],   # leukemic cells lack HSC quiescence (HLF/AVP), retain HOXA9/MEIS1
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
    ims = [Image.open(os.path.join(P, f)).convert("RGB") for f, _ in row]
    ratios = [im.width/im.height for im in ims]
    n = len(row); avail = W - GAP*(n-1)
    Hrow = min(avail/sum(ratios), MAXH)
    widths = [r*Hrow for r in ratios]
    tot = sum(widths) + GAP*(n-1); x = MARGIN + (W - tot)/2
    placed = []
    for (f, lab), im, w in zip(row, ims, widths):
        placed.append((im, lab, int(round(x)), int(round(w)))); x += w + GAP
    rowlays.append((int(round(Hrow)), placed))
totalH = TOP + sum(HEAD + h for h, _ in rowlays) + GAP*len(rowlays)
canvas = Image.new("RGB", (W + 2*MARGIN, int(totalH)), "white")
d = ImageDraw.Draw(canvas); y = TOP
for Hrow, placed in rowlays:
    for im, lab, x, w in placed:
        if lab: d.text((x, y), lab, font=fL, fill="black")
        canvas.paste(im.resize((w, Hrow), Image.LANCZOS), (x, y + HEAD))
    y += HEAD + Hrow + GAP
canvas.save(os.path.join(OUT, "Figure_3__LS10_proliferation.png"))
canvas.save(os.path.join(OUT, "Figure_3__LS10_proliferation.pdf"), "PDF", resolution=200)
print("wrote Figure_3__LS10_proliferation  %dx%d px  (aspect %.2f)" % (canvas.width, canvas.height, canvas.height/canvas.width))
