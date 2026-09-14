#!/usr/bin/env python
"""Figure 5A (corrected scheme): surfaceome discovery funnel + composite ranking, recomputed from the
main atlas on the corrected leukemic definition (Cell_Type in {AML, AML-MKI67, AML-PCNA, AML-CD1C},
AML samples). The retired PAC/PPAC axis is removed: the composite is now the mean of just
  C1 = % leukemic cells positive (efficacy)   and   C2 = % of the 28 AML patients with >=20% leukemic positivity (breadth).
Gate (official package): normal HSPC positivity ≤10% AND Myeloid_Pro ≤10%, then AML+ (>20% in >=1 patient),
T-cell-receptor genes dropped. Detection = raw count > 0. No fabricated values.
Writes fig5A_funnel.csv (funnel stage counts) and fig5A_top15.csv (full ranking; fig5.R takes the top 15)
to both the pipeline source_data/ and the 08_Source_Data mirror."""
import os, numpy as np, pandas as pd, scipy.sparse as sp, h5py

H5 = "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/ALSF_AML_plot/H5AD/ALSF_AML_Combo_3500_with_PAC_new.h5ad"
SURF = "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/ref/surfaceome_genes.txt"
DEST = ["/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/07_Code/shared_CD96_pipeline/source_data",
        "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/08_Source_Data/Figure_5"]
LEUK = {"AML", "AML-MKI67", "AML-PCNA", "AML-CD1C"}
GATE_HSPC, GATE_MYE, AML_MIN = 10.0, 10.0, 0.20
TCR = {"CD3D", "CD3E", "CD3G", "CD247", "TRAC", "TRBC1", "TRBC2", "TRGC1", "TRGC2", "TRDC"}

dec = lambda a: np.array([x.decode() if isinstance(x, bytes) else x for x in a])
def cat(f, k):
    o = f["obs"][k]; return dec(o["categories"][:])[o["codes"][:].astype(int)]

f = h5py.File(H5, "r")
ct = cat(f, "Cell Type"); st = cat(f, "SampleType"); sid = cat(f, "SampleID")
rk = f["raw"]["var"].attrs.get("_index", b"_index"); rk = rk.decode() if isinstance(rk, bytes) else rk
sym = list(dec(f["raw"]["var"][rk][:]))
rx = f["raw"]["X"]
X = sp.csr_matrix((rx["data"][:], rx["indices"][:], rx["indptr"][:]), shape=tuple(rx.attrs["shape"]))
f.close()

surf = [g.strip() for g in open(SURF) if g.strip()]
present = [g for g in surf if g in sym]
idx = [sym.index(g) for g in present]
posS = (X[:, idx] > 0).tocsr()                       # 96627 x len(present), sparse boolean
present = np.array(present)

leuk = np.isin(ct, list(LEUK)) & (st == "AML")
hspc = ct == "0_HSPC"; mye = ct == "Myeloid_Pro"
pct = lambda m: 100 * np.asarray(posS[m].mean(0)).ravel()
hspc_pct, mye_pct, aml_pct = pct(hspc), pct(mye), pct(leuk)
patients = sorted(set(sid[leuk]))
npts = np.zeros(len(present), int)
for s in patients:
    m = leuk & (sid == s)
    npts += (np.asarray(posS[m].mean(0)).ravel() >= AML_MIN)

df = pd.DataFrame(dict(gene=present,
    scRNA_HSPC_pct=hspc_pct.round(3), scRNA_Myeloid_pct=mye_pct.round(3),
    scRNA_AML_pct=aml_pct.round(2), scRNA_n_pts_above20pct=npts))
df["C1_pct_AML_cells"] = df.scRNA_AML_pct
df["C2_pct_AML_patients"] = df.scRNA_n_pts_above20pct / len(patients) * 100
df["composite"] = (df.C1_pct_AML_cells + df.C2_pct_AML_patients) / 2      # 2-axis (no PPAC)
df["is_tcr"] = df.gene.isin(TCR)

n_detect = len(df)
n_hspc = int((df.scRNA_HSPC_pct <= GATE_HSPC).sum())
n_mye = int(((df.scRNA_HSPC_pct <= GATE_HSPC) & (df.scRNA_Myeloid_pct <= GATE_MYE)).sum())
gated = df[(df.scRNA_HSPC_pct <= GATE_HSPC) & (df.scRNA_Myeloid_pct <= GATE_MYE) & (~df.is_tcr)].copy()
n_aml = int((gated.scRNA_n_pts_above20pct >= 1).sum())
ranked = gated.sort_values("composite", ascending=False).reset_index(drop=True)
ranked.insert(0, "rank", np.arange(1, len(ranked) + 1))

funnel = pd.DataFrame({"stage": ["surfaceome (in silico)", "HSPC-sparing (≤10%)", "+ Myeloid-sparing (≤10%)",
                                 "+ AML+ (>20%, ≥1 pt)", "top 15 composite"],
                       "n": [n_detect, n_hspc, n_mye, n_aml, 15]})
for d in DEST:
    funnel.to_csv(os.path.join(d, "fig5A_funnel.csv"), index=False)
    ranked.drop(columns=["is_tcr"]).to_csv(os.path.join(d, "fig5A_top15.csv"), index=False)

print("Fig 5A (corrected, 2-axis composite):"); print(funnel.to_string(index=False))
print(f"\n  patients: {len(patients)}   gated (TCR dropped): {len(gated)}   AML+(>20%,>=1pt): {n_aml}")
print("\n  top 15 by composite:")
print(ranked.head(15)[["rank", "gene", "composite", "scRNA_HSPC_pct", "scRNA_Myeloid_pct",
                       "scRNA_AML_pct", "scRNA_n_pts_above20pct"]].to_string(index=False))
cd = ranked[ranked.gene == "CD96"]
print("\n  CD96: rank %d, composite %.2f, leukemic %.1f%%, HSPC %.2f%%, Myeloid %.2f%%, patients>20%% = %d/%d"
      % (cd["rank"].iloc[0], cd.composite.iloc[0], cd.scRNA_AML_pct.iloc[0], cd.scRNA_HSPC_pct.iloc[0],
         cd.scRNA_Myeloid_pct.iloc[0], cd.scRNA_n_pts_above20pct.iloc[0], len(patients)))
