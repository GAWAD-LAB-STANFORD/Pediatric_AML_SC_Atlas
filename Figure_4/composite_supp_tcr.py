#!/usr/bin/env python
"""Supplementary Figure S13 - T/NK TCR repertoire (supports the Fig 4F clonal-focusing finding).
Three vector panels stacked with fitz -> Figure_S13__TNK_TCR_repertoire.pdf (+ 200-dpi .png):
  A  clonotype sharing between T subsets (Jaccard; TCR as a lineage barcode)
  B  exhaustion score, expanded vs unexpanded clones (honest null)
  C  T-cell clonality by cytogenetic subtype (honest null; per-subtype n shown)
Panels are the .pdf outputs of panel_tnk_tcr_extra.R (run that first). No baked figure title/subtitle
(number lives in the filename + legend); panel letters A-C are stamped in clean strips above each panel.
The TRBV usage panel was dropped (HBM normal reference n=2 only)."""
import os, fitz

P   = "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/04_Main_Figures/Figure_4_panels"
OUT = "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/05_Supplementary_Figures"
os.makedirs(OUT, exist_ok=True)
A = fitz.open(os.path.join(P, "Figure_4_tcr_sharing.pdf"))        # 5.6 x 4.8  clonotype sharing
B = fitz.open(os.path.join(P, "Figure_4_tcr_exhaustion.pdf"))     # 4.8 x 4.2  exhaustion (null)
C = fitz.open(os.path.join(P, "Figure_4_tcr_clonality_cyto.pdf")) # 6.2 x 4.2  clonality by cytogenetics (null)
BLACK = (0, 0, 0)

def ar(doc):                                   # aspect = width / height
    r = doc[0].rect; return r.width / r.height

PAGE_W = 11 * 72          # 11-inch-wide supplement
LH = 18                   # panel-letter strip above each panel/row
G  = 14                   # gap between panels / rows
PAD = 6
TOP = 6

# top row: A | B, common height
hb = (PAGE_W - G) / (ar(A) + ar(B))
aw, bw = ar(A) * hb, ar(B) * hb
# bottom row: C, centred, height kept close to the top row so the figure stays balanced
hC = min(PAGE_W / ar(C), hb * 1.02)
cw = ar(C) * hC
cx = (PAGE_W - cw) / 2

page_h = TOP + LH + hb + G + LH + hC + G
out = fitz.open()
pg = out.new_page(width=PAGE_W, height=page_h)

# A + B (top row); letters in the strip above the row
yrow = TOP + LH
pg.show_pdf_page(fitz.Rect(0, yrow, aw, yrow + hb), A, 0)
pg.show_pdf_page(fitz.Rect(aw + G, yrow, aw + G + bw, yrow + hb), B, 0)
pg.insert_text(fitz.Point(PAD, TOP + LH - 4), "A", fontsize=15, fontname="hebo", color=BLACK)
pg.insert_text(fitz.Point(aw + G + PAD, TOP + LH - 4), "B", fontsize=15, fontname="hebo", color=BLACK)

# C (bottom row, centred); letter in the strip above it
yLC = yrow + hb + G
yC = yLC + LH
pg.show_pdf_page(fitz.Rect(cx, yC, cx + cw, yC + hC), C, 0)
pg.insert_text(fitz.Point(PAD, yLC + LH - 4), "C", fontsize=15, fontname="hebo", color=BLACK)

pdf = os.path.join(OUT, "Figure_S13__TNK_TCR_repertoire.pdf")
out.save(pdf, deflate=True, garbage=3)
out[0].get_pixmap(dpi=200).save(os.path.join(OUT, "Figure_S13__TNK_TCR_repertoire.png"))
print("wrote Figure_S13__TNK_TCR_repertoire.pdf + .png  (%.1f x %.1f in)" % (PAGE_W/72, page_h/72))
