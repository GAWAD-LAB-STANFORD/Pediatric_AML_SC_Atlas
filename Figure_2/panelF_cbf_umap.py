#!/usr/bin/env python
"""CBF section, opening panel: where core-binding-factor AML maps on the atlas UMAP.
Leukemic cells (corrected 70,108 set) colored by cytogenetic origin: t(8;21) =
RUNX1/RUNX1T1, inv(16) = CBFB/MYH11; other leukemic greyed; non-leukemic faint.
Cytogenetic label is per-cell from the discovery sample. No fabricated values."""
import os, sys, numpy as np, pandas as pd, h5py, matplotlib.pyplot as plt
sys.path.insert(0, "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/07_Code/Figure_3")
import config as C
OUT = "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/04_Main_Figures/Figure_2_panels"

ra = pd.read_csv(C.REFA); ra["barcode"] = ra["barcode"].astype(str).str.replace(r"^b'|'$", "", regex=True)
cbfmap = {"RUNX1/RUNX1T1": "t(8;21)", "CBFB/MYH11": "inv(16)"}
ra["cat"] = ra["Cytogenetic"].map(cbfmap).fillna("other leukemic")
bc2cat = dict(zip(ra["barcode"], ra["cat"]))

f = h5py.File(C.H5, "r")
obs_names = [x.decode() if isinstance(x, bytes) else x for x in f["obs"]["_index"][:]]
umap = f["obsm/X_umap"][:]; f.close()
cat = np.array([bc2cat.get(b, "non-leukemic") for b in obs_names], dtype=object)
print("cells:", {k: int((cat==k).sum()) for k in ["t(8;21)","inv(16)","other leukemic","non-leukemic"]})

order = ["non-leukemic", "other leukemic", "t(8;21)", "inv(16)"]   # CBF on top
colmap = {"non-leukemic": "#eeeeee", "other leukemic": "#c9c9c9", "t(8;21)": "#0A9396", "inv(16)": "#CA6702"}
fig, axx = plt.subplots(figsize=(6.2, 5.6))
for k in order:
    m = cat == k
    axx.scatter(umap[m,0], umap[m,1], s=(1.2 if k in ("non-leukemic","other leukemic") else 3.2),
                c=colmap[k], linewidths=0, rasterized=True, label=k,
                alpha=(0.35 if k=="non-leukemic" else 0.9))
axx.set_xticks([]); axx.set_yticks([]); axx.set_frame_on(False)
axx.legend(loc="lower left", frameon=False, markerscale=3, fontsize=9, handletextpad=0.2)
plt.tight_layout()
for ext in ("png","pdf"):
    fig.savefig(os.path.join(OUT, f"Figure_2F_cbf_umap.{ext}"), dpi=200, bbox_inches="tight")
plt.close()
print("wrote Figure_2F_cbf_umap ->", OUT)
