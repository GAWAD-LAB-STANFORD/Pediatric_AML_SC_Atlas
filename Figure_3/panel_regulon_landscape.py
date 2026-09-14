#!/usr/bin/env python
"""[SUPERSEDED 2026-08 — output archived to 04_Main_Figures/Figure_2_panels/_superseded/. This single
full-landscape panel was split into the two deployed panels: the HOX/homeobox subset -> Figure 2C
(panel_regulon_hox.py) and the full state x regulon heatmap -> Figure S9C (panelA2_regulon_states.py,
composited by Figure_S9/scripts/composite_S9.py). Kept for provenance; not part of any deployed figure.]

Regulon landscape panel: SCENIC regulon activity (z-scored AUCell, top-2 regulons per
leukemic state, unioned) across the 49 leukemic states + 8 normal compartments. Column strip:
type (normal vs leukemic) + whole-cohort EFS prognosis. Data = compute_regulon_landscape.py.
No fabricated values."""
import os, sys, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import seaborn as sns
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
NORMCOL = "#7b3294"   # normal-compartment label/annotation color
def progcol(c): return NORMCOL if c.startswith("NORM_") else {"favorable": "#00798c", "poor": "#d1495b", "ns": "#dddddd"}[pr.prognosis.get(c, "ns")]
# vlag: muted soft blue-white-red diverging (seaborn), neutral white center — Nature/Cell single-cell convention
VLAG = sns.color_palette("vlag", as_cmap=True)
top = set()
for s in z.index:
    if s.startswith("LS_"): top |= set(z.loc[s].sort_values(ascending=False).head(2).index)
top = sorted(top)
M = z[top].T; M.columns = [collabel(c) for c in z.index]
cc = pd.DataFrame({"group": [progcol(c) for c in z.index]}, index=M.columns)
g = sns.clustermap(M, cmap=VLAG, center=0, vmin=-2, vmax=2, figsize=(17, 13), col_colors=cc,
    xticklabels=True, yticklabels=True, cbar_kws={"label": "regulon activity (z)"},
    dendrogram_ratio=(0.07, 0.09), colors_ratio=0.012)
g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), fontsize=7, rotation=90)
g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), fontsize=6)
for t in g.ax_heatmap.get_xticklabels():
    if t.get_text().startswith("norm "): t.set_color(NORMCOL); t.set_fontweight("bold")
g.fig.suptitle("SCENIC regulon landscape: 49 leukemic states + 8 normal compartments (purple=normal; teal/red=fav/poor EFS)", fontsize=11, y=1.0)
for e in ("png", "pdf"): g.savefig(os.path.join(OUT, f"Figure_regulon_landscape_v3.{e}"), dpi=160, bbox_inches="tight")
print(f"wrote Figure_regulon_landscape_v3 ({M.shape[0]} regulons x {M.shape[1]} cols)")
