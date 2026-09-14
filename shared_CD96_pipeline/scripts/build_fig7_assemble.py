#!/usr/bin/env python
# ====================================================================
# build_fig7_assemble.py  |  Figure 7 (builder)
# reproducibility header (auto-generated; edit freely below this block)
# Inputs : (reads upstream object / raw matrix — see body)
# Outputs: (writes to ../figures or ../source_data — see body)
# Run    : python build_fig7_assemble.py   (or: bash run_all.sh)  |  Deps: see REPRODUCIBILITY_CHECKLIST.md
# Source data for this figure is in __SUBMISSION_PACKAGE/08_Source_Data/.
# ====================================================================

"""Assemble Figure 7 (combination strategy):
  A  target<->cell-type interaction circos (CD96 / ITGAX / TNFRSF4)       -- full width, top
  B  CD96 + partner additive coverage vs marrow toxicity (within/across)  -- middle-left
  C  two strategies to extend CD96 (maximize efficacy / limit toxicity)   -- middle-right
  D  punch line — every 2-target combination landscape (balanced score)   -- full width
  E  winners quantified — ranked bars of the frontier combos (numbers)    -- full width, bottom
Vector panel PDFs stacked with fitz -> Figure7.pdf + 200-dpi Figure7.png. A "Figure 7" title
sits at the very top; panel letters A-E are stamped in clean strips ABOVE each panel/row so they
never overlap the figure content."""
import os
import fitz

# self-contained: assemble from this pipeline's own figures/ dir
FIG = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "figures"))
A = fitz.open(os.path.join(FIG, "cellchat", "fig7A_circos.pdf"))         # 16.5 x 6.2  circos
B = fitz.open(os.path.join(FIG, "Figure7_CD96_additive_window.pdf"))     # 13.5 x 7    additive window
C = fitz.open(os.path.join(FIG, "Figure7D_best_combo.pdf"))              # 9.5 x 7.8   best-combo paths
D = fitz.open(os.path.join(FIG, "Figure7D_punchline.pdf"))              # 15.5 x 6    punch line (all-pairs landscape)
E = fitz.open(os.path.join(FIG, "Figure7E_winners.pdf"))               # 16 x 4.9    winners quantified (ranked bars)
NAVY = (0.063, 0.227, 0.400)

def ar(doc):                                       # aspect = width / height
    r = doc[0].rect; return r.width / r.height

PAGE_W = 16 * 72
TH = 30          # figure-title band ("Figure 7")
LH = 18          # panel-letter strip above each panel/row
G = 12           # gap between rows / panels
PAD = 7

a_h = PAGE_W / ar(A)                               # circos spans the full width
hb = (PAGE_W - G) / (ar(B) + ar(C))                # common height of the B|C row
b_w, c_w = ar(B) * hb, ar(C) * hb
d_h = PAGE_W / ar(D)                               # punch-line landscape spans the full width
e_h = PAGE_W / ar(E)                               # winners-quantified bars span the full width
page_h = TH + LH + a_h + G + LH + hb + G + LH + d_h + G + LH + e_h + G

out = fitz.open()
pg = out.new_page(width=PAGE_W, height=page_h)
pg.insert_text(fitz.Point(10, 22), "Figure 7", fontsize=17, fontname="hebo", color=NAVY)

# A — circos (full width), letter in the strip above it
yA = TH + LH
pg.show_pdf_page(fitz.Rect(0, yA, PAGE_W, yA + a_h), A, 0)
pg.insert_text(fitz.Point(PAD, TH + LH - 3), "A", fontsize=15, fontname="hebo", color=NAVY)

# B (additive window) + C (best-combo), middle row, letters in the strip above the row
yL = yA + a_h + G
yrow = yL + LH
pg.show_pdf_page(fitz.Rect(0, yrow, b_w, yrow + hb), B, 0)
pg.show_pdf_page(fitz.Rect(b_w + G, yrow, b_w + G + c_w, yrow + hb), C, 0)
pg.insert_text(fitz.Point(PAD, yL + LH - 3), "B", fontsize=15, fontname="hebo", color=NAVY)
pg.insert_text(fitz.Point(b_w + G + PAD, yL + LH - 3), "C", fontsize=15, fontname="hebo", color=NAVY)

# D — punch line landscape (full width), the consolidated every-2-target view
yLD = yrow + hb + G
yD = yLD + LH
pg.show_pdf_page(fitz.Rect(0, yD, PAGE_W, yD + d_h), D, 0)
pg.insert_text(fitz.Point(PAD, yLD + LH - 3), "D", fontsize=15, fontname="hebo", color=NAVY)

# E — winners quantified (full width), ranked bars of the frontier combos with numbers
yLE = yD + d_h + G
yE = yLE + LH
pg.show_pdf_page(fitz.Rect(0, yE, PAGE_W, yE + e_h), E, 0)
pg.insert_text(fitz.Point(PAD, yLE + LH - 3), "E", fontsize=15, fontname="hebo", color=NAVY)

out.save(os.path.join(FIG, "Figure7.pdf"), deflate=True, garbage=3)
out[0].get_pixmap(dpi=200).save(os.path.join(FIG, "Figure7.png"))
print(f"wrote Figure7.pdf + Figure7.png  ({PAGE_W/72:.1f} x {page_h/72:.1f} in)  [A circos · B additive · C best-combo · D landscape · E winners]")
