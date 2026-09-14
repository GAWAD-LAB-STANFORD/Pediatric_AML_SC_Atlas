#!/usr/bin/env python
"""Composite Figure 2 as a roughly-square multi-column grid. Rows fill the width
(panels scaled to a shared row height, aspect preserved -> no distortion, no
cropped legends). Whole-cohort panels first (A-E), then CBF (F-H). Letters only."""
import os
from PIL import Image, ImageDraw, ImageFont
P   = "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/04_Main_Figures/Figure_2_new_panels"
OUT = "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/04_Main_Figures"

FILES = {"A":"Figure_2A_combined.png", "B":"Figure_2B_prognosis_umap.png",
         "C":"Figure_2C_states_subtypes.png", "D":"Figure_2CD_stemness_by_prognosis.png",
         "E":"Figure_2G_prognostic_GOBP.png", "F":"Figure_2C_CBF_limited.png",
         "G":"Figure_2_LS12_CBF_KM.png", "H":"Figure_2G_CBF_GOBP.png"}
ROWS = [["A"], ["B","C"], ["D","E"], ["F","G","H"]]     # whole cohort A-E, then CBF F-H

W, GAP, HEAD, TOP, MARGIN = 2500, 46, 76, 26, 30
MAXH = 0.55 * W
def font(sz):
    for p in ["/System/Library/Fonts/Supplemental/Arial Bold.ttf", "/System/Library/Fonts/Helvetica.ttc"]:
        try: return ImageFont.truetype(p, sz)
        except Exception: pass
    return ImageFont.load_default()
fL = font(62)
imgs = {k: Image.open(os.path.join(P, v)).convert("RGB") for k, v in FILES.items()}

rowlays = []
for row in ROWS:
    ratios = [imgs[k].width / imgs[k].height for k in row]
    n = len(row); avail = W - GAP*(n-1)
    Hrow = min(avail / sum(ratios), MAXH)
    widths = [r*Hrow for r in ratios]
    tot = sum(widths) + GAP*(n-1)
    x = MARGIN + (W - tot)/2
    placed = []
    for k, w in zip(row, widths):
        placed.append((k, int(round(x)), int(round(w)))); x += w + GAP
    rowlays.append((int(round(Hrow)), placed))

totalH = TOP + sum(HEAD + h for h, _ in rowlays) + GAP*len(rowlays)
canvas = Image.new("RGB", (W + 2*MARGIN, int(totalH)), "white")
d = ImageDraw.Draw(canvas)
y = TOP
for Hrow, placed in rowlays:
    for k, x, w in placed:
        d.text((x, y), k, font=fL, fill="black")
        canvas.paste(imgs[k].resize((w, Hrow), Image.LANCZOS), (x, y + HEAD))
    y += HEAD + Hrow + GAP

canvas.save(os.path.join(OUT, "Figure_2__leukemic_states.png"))
canvas.save(os.path.join(OUT, "Figure_2__leukemic_states.pdf"), "PDF", resolution=200)
print("wrote Figure_2__leukemic_states  %dx%d px  (aspect %.2f)" % (canvas.width, canvas.height, canvas.height/canvas.width))
