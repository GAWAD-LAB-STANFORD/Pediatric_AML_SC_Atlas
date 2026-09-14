#!/usr/bin/env python
# ====================================================================
# build_fig5CD_leads.py  |  Figure 5 (builder)
# reproducibility header (auto-generated; edit freely below this block)
# Inputs : (reads upstream object / raw matrix — see body)
# Outputs: (writes to ../figures or ../source_data — see body)
# Run    : python build_fig5CD_leads.py   (or: bash run_all.sh)  |  Deps: see REPRODUCIBILITY_CHECKLIST.md
# Source data for this figure is in __SUBMISSION_PACKAGE/08_Source_Data/.
# ====================================================================

"""ALTERNATE Figure 5 panels C + D, for the 5 LEAD targets (CD96, CD9, SUCNR1, AMN, IL1RAP).

C (head-to-head): every gated surface gene as a point, x = normal HSPC/Myeloid leak (max),
   y = % leukemic cells positive, size = % patients with >=20% leukemic positivity. The 5
   leads are coloured/labelled; the rest are grey context. -> fig5C_headtohead.csv
D (coverage): % positive + mean expression (CP10k log1p) of the 5 leads across cytogenetic
   subtypes and across PAC clusters. -> fig5D_cyto.csv, fig5D_pac.csv
All from the bundled count matrix.
"""
import os
import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
CM = os.path.join(ROOT, "count_matrix", "ALSF_AML_raw_counts.h5ad")
SURF = os.path.join(ROOT, "..", "ref", "surfaceome_genes.txt")
SD = os.path.join(ROOT, "source_data")

GATE = 10.0
N_PT = 28
LEADS = ["CD96", "CD9", "SUCNR1", "AMN", "IL1RAP"]
LEUK = [f"AML_{i}" for i in range(1, 28)]
PPACS = ["PPAC_1", "PPAC_2", "PPAC_3", "PPAC_4", "PPAC_5"]
CORE = ["AML", "AML-MKI67", "AML-PCNA"]
TCR = {"CD3D", "CD3E", "CD3G", "CD247", "TRAC", "TRBC1", "TRBC2", "TRGC1", "TRGC2", "TRDC"}
PACORD = ["0_HSPC", "Myeloid_Pro", "FPAC_1", "FPAC_2", "PPAC_1", "PPAC_2", "PPAC_3", "PPAC_4", "PPAC_5"]
CYTO_DROP = {"0_HealthyBM1", "0_HealthyBM2"}

A = ad.read_h5ad(CM)
pac = A.obs["Prognosis-Associated Clusters"].astype(str).values
sub = A.obs["AML Sub-Clusters"].astype(str).values
ct = A.obs["Cell Type"].astype(str).values
cyto = A.obs["Cytogenetic"].astype(str).values
samp = A.obs["Sample"].astype(str).values
pac_anno = A.obs["PAC_anno"].astype(str).values if "PAC_anno" in A.obs else pac
leuk = np.isin(sub, LEUK)
core = np.isin(ct, CORE)
hspc = pac == "0_HSPC"; mye = pac == "Myeloid_Pro"
tot = np.asarray(A.X.sum(1)).ravel().astype(float)          # total counts per cell (CP10k)
tot[tot == 0] = 1.0

# ---------- C: head-to-head over all gated surface genes ----------
surf = [g.strip() for g in open(SURF) if g.strip()]
present = [g for g in surf if g in A.var_names]
Xs = A[:, present].X
poss = (np.asarray(Xs.todense()) > 0) if sp.issparse(Xs) else (np.asarray(Xs) > 0)
present = np.array(present)
hspc_pct = 100 * poss[hspc].mean(0); mye_pct = 100 * poss[mye].mean(0)
aml_pct = 100 * poss[leuk].mean(0); ppac_pct = 100 * poss[np.isin(pac, PPACS)].mean(0)
samples = sorted({s for s in np.unique(samp) if "Healthy" not in s})
npts = np.zeros(len(present), int)
for s in samples:
    m = leuk & (samp == s)
    if m.sum():
        npts += (poss[m].mean(0) >= 0.20)
is_tcr = np.array([g in TCR for g in present])
gate = (hspc_pct <= GATE) & (mye_pct <= GATE) & (~is_tcr)
comp = (aml_pct + npts / N_PT * 100 + ppac_pct) / 3
C = pd.DataFrame(dict(gene=present, leak=np.maximum(hspc_pct, mye_pct).round(2),
                      scRNA_AML_pct=aml_pct.round(2), scRNA_HSPC_pct=hspc_pct.round(2),
                      scRNA_Myeloid_pct=mye_pct.round(2),
                      n_pts_pct=(npts / N_PT * 100).round(1), composite=comp.round(2)))[gate]
C["is_lead"] = C.gene.isin(LEADS)
# keep the 5 leads + top-35 others by composite, for a clean head-to-head backdrop
keep = pd.concat([C[C.is_lead], C[~C.is_lead].sort_values("composite", ascending=False).head(35)])
keep = keep.drop_duplicates("gene").sort_values("composite", ascending=False)
keep.to_csv(os.path.join(SD, "fig5C_headtohead.csv"), index=False)
print(f"C: {len(keep)} genes ({keep.is_lead.sum()} leads + context) -> fig5C_headtohead.csv")

# ---------- D: per-group coverage for the 5 leads ----------
Xl = A[:, LEADS].X
Ld = np.asarray(Xl.todense()) if sp.issparse(Xl) else np.asarray(Xl)
Lpos = Ld > 0
Lnorm = np.log1p(Ld / tot[:, None] * 1e4)                   # CP10k log1p

def coverage(groups, mask_of, cell_subset=None):
    rows = []
    for g in groups:
        m = mask_of(g)
        if cell_subset is not None:
            m = m & cell_subset
        if m.sum() < 20:
            continue
        for i, gene in enumerate(LEADS):
            rows.append(dict(gene=gene, group=g, pct=round(100 * Lpos[m][:, i].mean(), 2),
                             label=gene, mean_expr=round(float(Lnorm[m][:, i].mean()), 4)))
    return pd.DataFrame(rows)

cyto_groups = [c for c in pd.unique(cyto) if c not in CYTO_DROP]
Dc = coverage(cyto_groups, lambda g: cyto == g, cell_subset=core)   # leukemic cells per cytogenetics
Dc.to_csv(os.path.join(SD, "fig5D_cyto.csv"), index=False)
Dp = coverage(PACORD, lambda g: pac == g)                          # all cells per PAC cluster
Dp.to_csv(os.path.join(SD, "fig5D_pac.csv"), index=False)
print(f"D: cyto {Dc.group.nunique()} subtypes, pac {Dp.group.nunique()} clusters -> fig5D_cyto.csv, fig5D_pac.csv")
print("\n  leads head-to-head (leak / AML% / %pts):")
print(keep[keep.is_lead][["gene", "leak", "scRNA_AML_pct", "n_pts_pct"]].to_string(index=False))
print("\n  CD96 sanity — cyto coverage (should match official ~63% BCR/ABL):")
print(Dc[(Dc.gene == "CD96")][["group", "pct"]].sort_values("pct", ascending=False).head(4).to_string(index=False))
