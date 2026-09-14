#!/usr/bin/env python
"""Supplementary figure: pediatric AML lacks a primitive quiescent-LSC compartment.
A strict-LSC gating funnel (AML vs healthy BM). B those 26 candidates on the atlas UMAP
(leukemic-localized, not the HSC cluster). C the candidates carry GMP/granule
differentiation markers (aberrant GMP-like, not primitive LSC). Letters only."""
import os
from PIL import Image, ImageDraw, ImageFont
P   = "../figures/Figure_2_panels"
OUT = "../figures"
os.makedirs(OUT, exist_ok=True)
ROWS = [
  [("Figure_3_strict_LSC_funnel.png","A"), ("UMAP_strict_LSC.png","B")],
  [("Figure_3_strictLSC_markers.png","C")],
]
W, GAP, HEAD, TOP, MARGIN = 2200, 44, 70, 22, 28
MAXH = 0.62 * W
def font(sz):
    for p in ["/System/Library/Fonts/Supplemental/Arial Bold.ttf", "/System/Library/Fonts/Helvetica.ttc"]:
        try: return ImageFont.truetype(p, sz)
        except Exception: pass
    return ImageFont.load_default()
fL = font(56)
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
        d.text((x, y), lab, font=fL, fill="black")
        canvas.paste(im.resize((w, Hrow), Image.LANCZOS), (x, y + HEAD))
    y += HEAD + Hrow + GAP
canvas.save(os.path.join(OUT, "Figure_S_LSC_scarcity.png"))
canvas.save(os.path.join(OUT, "Figure_S_LSC_scarcity.pdf"), "PDF", resolution=200)
print("wrote Figure_S_LSC_scarcity  %dx%d px  (aspect %.2f)" % (canvas.width, canvas.height, canvas.height/canvas.width))
