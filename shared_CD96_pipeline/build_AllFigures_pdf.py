#!/usr/bin/env python
"""Build the clean combined figure PDF: figures only (no cover / TOC / divider pages),
with each page renumbered in sequence. Legacy manuscript pages have their stale top-left
"Figure.X" label covered with a white box and the new label stamped over it; the CD96
pages (no baked label) get the new label stamped in the top-left corner.
Requires PyMuPDF (fitz)."""
import os, fitz

FIN = "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/AML_2026_Final"
MS  = ("/Users/chuckgawad/Library/CloudStorage/Box-Box/01_Gawad Lab Documents/01_Manuscripts/"
       "AML_Submission_2025/00_Blood?_2025/AML single cell RNA-seq_Figures_2026 copy.pdf")
def f(p): return os.path.join(FIN, p)

ms = fitz.open(MS)
fig5, fig6, fig7 = fitz.open(f("figures/Figure5.pdf")), fitz.open(f("figures/Figure6.pdf")), fitz.open(f("figures/Figure7.pdf"))
figS12 = fitz.open(f("figures/supplementary/FigureS6_nonMLL_violins.pdf"))
figS13 = fitz.open(f("figures/supplementary/FigureS7_correlations.pdf"))
figS14 = fitz.open(f("figures/FigureS8.pdf"))
figS15 = fitz.open(f("figures/supplementary/FigureS15_PAC_window.pdf"))
figS16 = fitz.open(f("figures/supplementary/FigureS16_cyto_window.pdf"))
figS17 = fitz.open(f("figures/supplementary/FigureS17_ppac_specific.pdf"))
figS18 = fitz.open(f("figures/supplementary/FigureS18_ppac_headtohead.pdf"))
figS19 = fitz.open(f("figures/supplementary/FigureS19_ppac_additive.pdf"))

# (source_doc, 0-based page, new label, kind)  kind: 'legacy' (cover old label) | 'cd96'
order = [
    (ms, 0,  "Figure 1",  "legacy"),
    (ms, 3,  "Figure 2",  "legacy"),
    (ms, 6,  "Figure 3",  "legacy"),
    (ms, 10, "Figure 4",  "legacy"),
    (fig5, 0, "Figure 5", "cd96"),
    (fig6, 0, "Figure 6", "cd96"),
    (fig7, 0, "Figure 7", "cd96"),
    (ms, 1,  "Figure S1",  "legacy"),
    (ms, 2,  "Figure S2",  "legacy"),
    (ms, 4,  "Figure S3",  "legacy"),
    (ms, 5,  "Figure S4",  "legacy"),
    (ms, 7,  "Figure S5",  "legacy"),
    (ms, 8,  "Figure S6",  "legacy"),   # was main Figure 4 (hdWGCNA)
    (ms, 9,  "Figure S7",  "legacy"),
    (ms, 11, "Figure S8",  "legacy"),
    (ms, 12, "Figure S9",  "legacy"),
    (ms, 13, "Figure S10", "legacy"),
    (ms, 14, "Figure S11", "legacy"),
    (figS12, 0, "Figure S12", "cd96"),
    (figS13, 0, "Figure S13", "cd96"),
    (figS14, 0, "Figure S14", "cd96"),
    (figS15, 0, "Figure S15", "cd96"),
    (figS16, 0, "Figure S16", "cd96"),
    (figS17, 0, "Figure S17", "cd96"),
    (figS18, 0, "Figure S18", "cd96"),
    (figS19, 0, "Figure S19", "cd96"),
    (ms, 15, "Table S1",   "legacy"),
]

NAVY = (0.063, 0.227, 0.400)
out = fitz.open()
for src, pno, label, kind in order:
    if kind == "legacy":
        # cover ONLY the stale top-left label (detected bbox), then stamp the new one;
        # never reaches the panel 'A' tag, which sits lower.
        out.insert_pdf(src, from_page=pno, to_page=pno)
        pg = out[-1]; w, h = pg.rect.width, pg.rect.height
        labw = [wd for wd in pg.get_text("words") if wd[1] < 0.045 * h and wd[0] < 0.42 * w]
        ybot = max((wd[3] for wd in labw), default=0.032 * h)
        pg.draw_rect(fitz.Rect(0, 0, 0.42 * w, ybot + 4), color=(1, 1, 1), fill=(1, 1, 1))
        pg.insert_text(fitz.Point(0.013 * w, ybot - 2), label, fontsize=16, fontname="hebo", color=NAVY)
    else:
        # CD96 pages have no margin; add a clean white header strip ABOVE the figure
        # (covers no content) and stamp the label there.
        sp = src[pno]; w, h = sp.rect.width, sp.rect.height
        hdr = 0.030 * h
        npg = out.new_page(width=w, height=h + hdr)
        npg.show_pdf_page(fitz.Rect(0, hdr, w, h + hdr), src, pno)
        npg.insert_text(fitz.Point(0.013 * w, hdr * 0.70), label,
                        fontsize=min(22, hdr * 0.68), fontname="hebo", color=NAVY)

OUT = f("AML_2026_AllFigures.pdf")
out.save(OUT, deflate=True, garbage=3)
print("WROTE", OUT, "|", out.page_count, "pages")
