#!/usr/bin/env python
"""Figure 2A (v3, corrected 49 states): full-atlas UMAP, two panels -- CD34 expression
and the 49 leukemic states (non-leukemic cells greyed). Leukemic cells + state labels
come from the corrected reference_assignments (Cell Type in {AML,AML-MKI67,AML-PCNA,
AML-CD1C}); non-stable-cluster leukemic cells are greyed. No fabricated values."""
import os, sys, numpy as np, pandas as pd, h5py, anndata as ad, scanpy as sc, matplotlib, matplotlib.pyplot as plt
from scipy.sparse import csr_matrix
from matplotlib.colors import LinearSegmentedColormap
sys.path.insert(0, "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/07_Code/Figure_3")
import config as C
HEATMAP0 = LinearSegmentedColormap.from_list("heatmap0",
  ["#001219","#005F73","#0A9396","#94D2BD","#E9D8A6","#EE9B00","#CA6702","#AE2012","#9B2226"])
OUT = "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/04_Main_Figures/Figure_2_panels"
os.makedirs(OUT, exist_ok=True)

# --- per-cell state labels (corrected) ---
ra = pd.read_csv(C.REFA); ra["barcode"] = ra["barcode"].astype(str).str.replace(r"^b'|'$", "", regex=True)
raws, lslabs = C.stable_raws(); raw2ls = {rw: lslabs[i] for i, rw in enumerate(raws)}
ra["LS"] = ra["res4.50"].map(raw2ls)
bc2ls = dict(zip(ra["barcode"], ra["LS"]))
NSTATE = len(lslabs)
_qual = []                                                       # QUALITATIVE categorical palette for state IDENTITY:
for _n in ("tab20","tab20b","tab20c"):                           # was a sequential HEATMAP0 ramp sampled by index, which
    _qual += [matplotlib.colors.to_hex(c) for c in matplotlib.colormaps[_n].colors]  # collided with the CD34 expression colorbar.
pal = _qual[:NSTATE]                                             # 49 distinct hues; LS labels carry identity, colour is not a scale.

# --- atlas UMAP + CD34 (CPM, log1p) ---
f = h5py.File(C.H5, "r")
obs_names = [x.decode() if isinstance(x, bytes) else x for x in f["obs"]["_index"][:]]
umap = f["obsm/X_umap"][:]; N = len(obs_names)
genes = np.array([g.decode() if isinstance(g, bytes) else g for g in f["var"]["_index"][:]])
cd34 = int(np.where(genes == "CD34")[0][0]); G = len(genes)
RC = f["layers"]["Raw_Counts"]; dd, ii, ip = RC["data"], RC["indices"], RC["indptr"][:]
cvec = np.zeros(N, dtype=np.float32)
for a in range(0, N, 8000):
    b = min(a+8000, N); s, e = int(ip[a]), int(ip[b])
    X = csr_matrix((dd[s:e].astype(np.float32), ii[s:e], (ip[a:b+1]-ip[a]).astype(np.int64)), shape=(b-a, G))
    rs = np.asarray(X.sum(1)).ravel(); rs[rs == 0] = 1
    cvec[a:b] = np.asarray(X[:, cd34].todense()).ravel()*(1e6/rs)
f.close()

state = np.array([bc2ls.get(b) for b in obs_names], dtype=object)   # LS label or None (non-leukemic / non-stable)
state = np.array([s.replace("_", "") if isinstance(s, str) else s for s in state], dtype=object)  # no underscores in labels
disp = [s.replace("_", "") for s in lslabs]
adata = ad.AnnData(X=np.zeros((N, 1), dtype=np.float32),
    obs=pd.DataFrame({"Leukemic states": pd.Categorical(state, categories=disp),
                      "CD34": np.log1p(cvec)}), var=pd.DataFrame(index=["_dummy"]))
adata.obsm["X_umap"] = umap; adata.uns["Leukemic states_colors"] = pal
sc.settings.figdir = OUT; sc.set_figure_params(dpi=150, dpi_save=300, frameon=False)
sc.pl.umap(adata, color="CD34", title="", frameon=False, show=False, cmap=HEATMAP0)
_fig = plt.gcf()                                             # label the colorbar so the 0-8 scale reads as CD34
if len(_fig.axes) > 1: _fig.axes[-1].set_ylabel("CD34 log(CPM+1)", fontsize=8, rotation=90)
plt.savefig(os.path.join(OUT, "Figure_2A_CD34.pdf"), bbox_inches="tight"); plt.savefig(os.path.join(OUT, "Figure_2A_CD34.png"), bbox_inches="tight", dpi=200); plt.close()
# colour by state WITHOUT scanpy's on-data labels (49 collide); place repelled centroid labels instead
from adjustText import adjust_text
import matplotlib.patheffects as pe
sc.pl.umap(adata, color="Leukemic states", title="", frameon=False, show=False, legend_loc=None)
axc = plt.gca(); texts = []
for s in disp:
    m = state == s
    if m.sum() == 0: continue
    cx, cy = float(np.median(umap[m, 0])), float(np.median(umap[m, 1]))
    texts.append(axc.text(cx, cy, s, fontsize=5, fontweight="bold", ha="center", va="center",
                          path_effects=[pe.withStroke(linewidth=1.6, foreground="white")]))
adjust_text(texts, ax=axc, expand=(1.15, 1.3), force_text=(0.4, 0.6),
            arrowprops=dict(arrowstyle="-", color="#9aa0a6", lw=0.3))
plt.savefig(os.path.join(OUT, "Figure_2A_states.pdf"), bbox_inches="tight"); plt.savefig(os.path.join(OUT, "Figure_2A_states.png"), bbox_inches="tight", dpi=200); plt.close()
print(f"wrote Figure_2A_CD34 + Figure_2A_states ({NSTATE} states, {int((state!=None).sum())} leukemic cells labelled) to {OUT}")
