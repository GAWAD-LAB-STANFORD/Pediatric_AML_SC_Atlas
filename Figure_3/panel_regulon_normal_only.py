#!/usr/bin/env python
"""[SUPERSEDED 2026-08 — diagnostic; output archived to 04_Main_Figures/Figure_2_panels/_superseded/.
The normal-compartment regulon validation it provided is covered by Figure S9A/B. Kept for provenance;
not part of any deployed figure.]

Regulon reference panel: SCENIC regulon activity across the 8 normal hematopoietic
compartments only (z-scored within normal). Each compartment recovers its canonical master
regulators (B: PAX5/EBF1; T/NK: TCF7/LEF1; HSPC: HOXA9/ERG; DC: TCF4/BATF3; erythroid: ESRRG)
-- a validation of the AUCell aggregation. Data = compute_regulon_landscape.py. No fabrication."""
import os, sys, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import seaborn as sns
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import config as C
D3 = C.DATA; OUT = os.path.join(C.BASE, "05_Supplementary_Figures")   # normal regulon reference -> supplement
os.makedirs(OUT, exist_ok=True)
mean = pd.read_csv(os.path.join(D3, "regulon_landscape_mean.csv"), index_col=0)
norm = [c for c in mean.index if c.startswith("NORM_")]
mn = mean.loc[norm]; z = (mn - mn.mean(0)) / (mn.std(0) + 1e-9)   # z within normal compartments
top = set()
for c in z.index: top |= set(z.loc[c].sort_values(ascending=False).head(6).index)
M = z[sorted(top)].T; M.columns = [c.replace("NORM_", "") for c in z.index]
VLAG = sns.color_palette("vlag", as_cmap=True)
g = sns.clustermap(M, cmap=VLAG, center=0, vmin=-2, vmax=2, figsize=(8, 13), xticklabels=True, yticklabels=True,
    cbar_kws={"label": "regulon activity (z, within normal)"}, dendrogram_ratio=(0.12, 0.06))
g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), fontsize=9, rotation=45, ha="right")
g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), fontsize=6)
g.fig.suptitle("SCENIC regulon activity across normal hematopoietic compartments (reference)", fontsize=11, y=1.0)
for e in ("png", "pdf"): g.savefig(os.path.join(OUT, f"Figure_S_regulon_normal.{e}"), dpi=170, bbox_inches="tight")
print(f"wrote Figure_S_regulon_normal ({M.shape[0]} regulons x {M.shape[1]} normal compartments)")
