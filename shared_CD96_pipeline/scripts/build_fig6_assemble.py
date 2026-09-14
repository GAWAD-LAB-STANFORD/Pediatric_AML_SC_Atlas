#!/usr/bin/env python
# ====================================================================
# build_fig6_assemble.py  |  Figure 6 (builder)
# reproducibility header (auto-generated; edit freely below this block)
# Inputs : (reads upstream object / raw matrix — see body)
# Outputs: (writes to ../figures or ../source_data — see body)
# Run    : python build_fig6_assemble.py   (or: bash run_all.sh)  |  Deps: see REPRODUCIBILITY_CHECKLIST.md
# Source data for this figure is in __SUBMISSION_PACKAGE/08_Source_Data/.
# ====================================================================

"""Assemble Figure 6 = A/B/C toxicity dot plots (top) over D radar profiles (bottom).
Stacks the two VECTOR panel PDFs (.panel_fig6ABC.pdf from fig6_toxicity.R, .panel_fig6D.pdf
from fig6_DE.R) onto one page with fitz -> vector Figure6.pdf + 200-dpi Figure6.png.
Replaces fig6_assemble.R (the readRDS+patchwork round-trip cannot re-wrap nested facets)."""
import os
import fitz

# self-contained: assemble from this pipeline's own figures/ dir
FIG = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "figures"))
abc = fitz.open(os.path.join(FIG, ".panel_fig6ABC.pdf"))
rad = fitz.open(os.path.join(FIG, ".panel_fig6D.pdf"))

PAGE_W = 16 * 72                                  # 16 in at 72 pt/in
GAP = 12
aw, ah = abc[0].rect.width, abc[0].rect.height
rw, rh = rad[0].rect.width, rad[0].rect.height
abc_h = PAGE_W * ah / aw                          # A/B/C spans the full width
rad_w = 0.80 * PAGE_W                             # radar centred, a touch narrower
rad_h = rad_w * rh / rw
page_h = GAP + abc_h + GAP + rad_h + GAP

out = fitz.open()
pg = out.new_page(width=PAGE_W, height=page_h)
pg.show_pdf_page(fitz.Rect(0, GAP, PAGE_W, GAP + abc_h), abc, 0)
dx = (PAGE_W - rad_w) / 2
pg.show_pdf_page(fitz.Rect(dx, GAP + abc_h + GAP, dx + rad_w, GAP + abc_h + GAP + rad_h), rad, 0)

out.save(os.path.join(FIG, "Figure6.pdf"), deflate=True, garbage=3)
out[0].get_pixmap(dpi=200).save(os.path.join(FIG, "Figure6.png"))
print(f"wrote Figure6.pdf + Figure6.png  ({PAGE_W/72:.1f} x {page_h/72:.1f} in)")
