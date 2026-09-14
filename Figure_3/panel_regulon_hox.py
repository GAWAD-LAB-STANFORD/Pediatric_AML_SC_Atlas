#!/usr/bin/env python
"""HOX regulon panel (Figure 2C): homeobox/HOX regulon activity (z-scored AUCell) across the 49
leukemic states + 8 normal compartments. HOX regulons on ROWS (labels coloured by homeobox family);
states on COLUMNS (labels coloured by EFS prognosis). HOXA-family self-renewal fires in the stem/
early-progenitor states, not uniformly across the poor-prognosis states. Data = compute_regulon_landscape.py."""
import os, sys, re, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import seaborn as sns
import matplotlib.patches as mpatches
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import config as C
D3 = C.DATA; D2 = D3.replace("Figure_3", "Figure_2")
OUT = os.path.join(C.BASE, "04_Main_Figures/Figure_2_panels")
z = pd.read_csv(os.path.join(D3, "regulon_landscape_z.csv"), index_col=0)
pr = pd.read_csv(os.path.join(D2, "LS_wholecohort_prognosis.csv")).set_index("LS")
ann = pd.read_csv(os.path.join(D3, "LS_annotation.csv")).set_index("LS")
def collabel(c):
    if c.startswith("NORM_"): return c.replace("NORM_", "norm ")
    p = str(ann.lineage_program.get(c, "")); p = "" if p == "unresolved" else p
    return (c.replace("_", "") + " " + p).strip()
NORMCOL = "#7b3294"
PROG = {"favorable": "#0A9396", "poor": "#AE2012", "ns": "#9aa0a6"}          # manuscript palette
def progcol(c): return NORMCOL if c.startswith("NORM_") else PROG[pr.prognosis.get(c, "ns")]
FAM = [("HOXA", "#D55E00", "HOXA (self-renewal)"), ("HOXB", "#CC79A7", "HOXB"),
       ("HOXC", "#E69F00", "HOXC/D"), ("HOXD", "#E69F00", None),
       ("PBX", "#0072B2", "PBX/MEIS"), ("MEIS", "#0072B2", None)]
OTHERCOL = "#009E73"
def famcol(reg):
    for pre, col, _ in FAM:
        if reg.startswith(pre): return col
    return OTHERCOL
hox = [r for r in z.columns if re.match(r"^(HOXA|HOXB|HOXC|HOXD|MEIS|PBX|CDX|HMX|DLX|PITX|SHOX|LHX|NKX)\d", r)]
H = z[hox].T; H.columns = [collabel(c) for c in z.index]                    # rows = HOX regulons, cols = states
cc = pd.DataFrame({"EFS prognosis": [progcol(c) for c in z.index]}, index=H.columns)
g = sns.clustermap(H, cmap=sns.color_palette("vlag", as_cmap=True), center=0, vmin=-2, vmax=2,
    figsize=(15, 8.5), col_colors=cc, xticklabels=True, yticklabels=True,
    cbar_kws={"label": "regulon activity (z)"}, dendrogram_ratio=(0.06, 0.12), colors_ratio=0.02)
g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), fontsize=7, rotation=90)
g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), fontsize=8)
labcol = {collabel(c): progcol(c) for c in z.index}
for t in g.ax_heatmap.get_xticklabels():                                    # state (col) labels by prognosis
    t.set_color(labcol.get(t.get_text(), "#000000"))
    if t.get_text().startswith("norm "): t.set_fontweight("bold")
for t in g.ax_heatmap.get_yticklabels():                                    # HOX regulon (row) labels by family
    t.set_color(famcol(t.get_text())); t.set_fontweight("bold")
fam_handles = [mpatches.Patch(color=col, label=lab) for pre, col, lab in FAM if lab] + \
              [mpatches.Patch(color=OTHERCOL, label="other (lineage)")]
prog_handles = [mpatches.Patch(color=PROG["favorable"], label="favorable EFS"),
                mpatches.Patch(color=PROG["poor"], label="poor EFS"), mpatches.Patch(color=NORMCOL, label="normal")]
leg1 = g.ax_heatmap.legend(handles=fam_handles, title="HOX family (rows)",
    loc="upper left", bbox_to_anchor=(1.15, 1.0), fontsize=8, title_fontsize=9, frameon=True, framealpha=0.9)
g.ax_heatmap.add_artist(leg1)
leg2 = g.ax_heatmap.legend(handles=prog_handles, title="EFS (states)",
    loc="upper left", bbox_to_anchor=(1.15, 0.5), fontsize=8, title_fontsize=9, frameon=True, framealpha=0.9)
g.fig.suptitle("Homeobox/HOX regulon activity: HOXA-family self-renewal fires in stem/progenitor states, "
               "not uniformly across the poor-prognosis states", fontsize=10, y=1.02)
for e in ("png", "pdf"):
    g.savefig(os.path.join(OUT, f"Figure_regulon_HOX_v3.{e}"), dpi=170,
              bbox_inches="tight", bbox_extra_artists=(leg1, leg2))
print(f"wrote Figure_regulon_HOX_v3 (wide: {len(hox)} HOX regulons x {H.shape[1]} cols)")
