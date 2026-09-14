#!/usr/bin/env python
"""Figure S9 (regulon detail) compositor.

S9 has no in-tree renderer: panels A (per-population regulon rank plots) and B (population x
regulon activity heatmap) come from an external pySCENIC/AUCell run, shipped as a static PDF
(Figure_S9__regulon_detail_AB_source.pdf). This script:
  1. renders that A/B source,
  2. patches the stale baked-in "Figure S5" title -> "Figure S9" (the file predates the S5->S9
     supplementary renumber and the number is baked into the PDF, not added at layout),
  3. appends the full 57-state x 89-regulon activity heatmap (Figure_2_regulon_states.png, the
     panel moved out of Figure 2C when 2C became the HOX-only view) as panel C,
and writes the deployed Figure_S9__regulon_detail.{png,pdf}. Rerun after regenerating either input.
"""
import os, io
import fitz  # pymupdf
from PIL import Image, ImageDraw, ImageFont

BASE = "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster"
SUP  = os.path.join(BASE, "05_Supplementary_Figures")
AB   = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "Figure_S9__regulon_detail_AB_source.pdf")  # external A/B (pySCENIC), kept in the code tree (not in deployed figures)
GEN  = os.path.join(BASE, "04_Main_Figures/Figure_2_panels/Figure_2_regulon_states.png")  # panel C
OUT  = os.path.join(SUP, "Figure_S9__regulon_detail")

NAVY = (30, 68, 110)        # sampled from the source title
def font(sz, bold=True):
    for p in (["/System/Library/Fonts/Supplemental/Arial Bold.ttf"] if bold else []) + \
             ["/System/Library/Fonts/Supplemental/Arial.ttf", "/System/Library/Fonts/Helvetica.ttc"]:
        try: return ImageFont.truetype(p, sz)
        except Exception: pass
    return ImageFont.load_default()

# 1. render A/B source
ab = Image.open(io.BytesIO(fitz.open(AB)[0].get_pixmap(dpi=200).tobytes("png"))).convert("RGB")
# 2. patch stale "Figure S5" -> "Figure S9" (title bbox ~ (25,29)-(222,71) at dpi=200)
dr = ImageDraw.Draw(ab)
dr.rectangle((14, 18, 275, 82), fill="white")
dr.text((25, 24), "Figure S9", font=font(50), fill=NAVY)

# 3. load panel C (full state-level regulon heatmap) and match widths
g = Image.open(GEN).convert("RGB")
W = max(ab.width, g.width)
sw = lambda im: im.resize((W, round(im.height * W / im.width)), Image.LANCZOS)
ab, g = sw(ab), sw(g)

# 4. stack: A/B on top, labelled panel C below
HEAD, GAP = 74, 26
canvas = Image.new("RGB", (W, ab.height + GAP + HEAD + g.height), "white")
canvas.paste(ab, (0, 0))
ImageDraw.Draw(canvas).text((12, ab.height + GAP + 4), "C", font=font(58), fill="black")
canvas.paste(g, (0, ab.height + GAP + HEAD))

canvas.save(OUT + ".png")
canvas.save(OUT + ".pdf", "PDF", resolution=200)
print("wrote Figure_S9__regulon_detail.{png,pdf}  %dx%d  (A/B source + panel C state heatmap)" % canvas.size)
