#!/usr/bin/env python
"""Figure 4 rebuild, step 3: (a) per-subset marker mean + % expressing for the defining dot plot, and
(b) a CIBERSORTx single-cell reference matrix (raw counts, columns labelled by T/NK subset) for the user
to run CIBERSORTx on TARGET bulk -> finer-subset (GZMB/GZMK-CD8, NK, etc.) deconvolution + survival.
Reads the h5ad T/NK cells matched to the annotated subsets. No fabricated values."""
import os, numpy as np, pandas as pd, scanpy as sc
H5 = "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/ALSF_AML_plot/H5AD/ALSF_AML_Combo_3500_with_PAC_new.h5ad"
D  = "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/08_Source_Data/Figure_4/v2"
MARKERS = ["CD3D","CD4","CD8A","CD8B","CCR7","SELL","TCF7","IL7R","FOXP3","CTLA4",
           "GZMK","GZMB","GNLY","PRF1","NKG7","NCAM1","FCGR3A","KLRD1","SLC4A10","KLRB1","MKI67"]
SUBSET_ORDER = ["Naïve CD4 T","Memory CD4 T","Treg","Naïve CD8 T","Memory CD8 T","GZMK CD8 T",
                "GZMB CD8 T","GZMB DNT","MAIT","GZMK NK","GZMB NK","Proliferating T"]
pc = pd.read_csv(os.path.join(D, "tnk_subset_percell.csv")).set_index("barcode")
ad = sc.read_h5ad(H5)
ad = ad[ad.obs_names.isin(pc.index)].copy()
sub = pc["subset"].reindex(ad.obs_names)
# log-norm from raw counts
ad.X = ad.layers["Raw_Counts"].copy(); sc.pp.normalize_total(ad, target_sum=1e4); sc.pp.log1p(ad)
present = [g for g in MARKERS if g in ad.var_names]
E = pd.DataFrame(ad[:, present].X.toarray(), columns=present, index=ad.obs_names); E["subset"] = sub.values
rows = []
for s in SUBSET_ORDER:
    d = E[E.subset == s]
    if len(d) == 0: continue
    for g in present:
        rows.append([s, g, d[g].mean(), (d[g] > 0).mean()])
pd.DataFrame(rows, columns=["subset","gene","mean","pct"]).to_csv(os.path.join(D, "tnk_subset_dotplot.csv"), index=False)
print("wrote tnk_subset_dotplot.csv (%d subsets x %d markers)" % (sub.nunique(), len(present)))

# ---- CIBERSORTx single-cell reference (raw counts, columns = subset labels), subsampled ----
rng = np.random.RandomState(0)
keep = []
for s in SUBSET_ORDER:
    idx = np.where(sub.values == s)[0]
    if len(idx) == 0: continue
    keep.append(idx if len(idx) <= 120 else rng.choice(idx, 120, replace=False))
keep = np.concatenate(keep)
ref = ad[keep].copy()
counts = ref.layers["Raw_Counts"]
mat = pd.DataFrame(counts.toarray().T.astype(int), index=ref.var_names,
                   columns=list(sub.values[keep]))          # genes x cells, header = subset label
mat.index.name = "GeneSymbol"
# drop all-zero genes to shrink the file
mat = mat[mat.sum(axis=1) > 0]
ref_path = os.path.join(D, "tnk_cibersortx_reference.txt")
mat.to_csv(ref_path, sep="\t")
print("wrote CIBERSORTx reference:", ref_path, "| genes x cells =", mat.shape,
      "| cells/subset:", pd.Series(mat.columns).value_counts().to_dict())
