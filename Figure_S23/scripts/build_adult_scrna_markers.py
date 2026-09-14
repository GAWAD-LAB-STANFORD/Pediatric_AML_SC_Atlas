#!/usr/bin/env python
# ====================================================================
# build_adult_scrna_markers.py  |  Upstream builder (regenerates source_data from raw inputs)
# reproducibility header (auto-generated; edit freely below this block)
# Inputs : (reads upstream object / raw matrix — see body)
# Outputs: (writes to ../figures or ../source_data — see body)
# Run    : python build_adult_scrna_markers.py   (or: bash run_all.sh)  |  Deps: see REPRODUCIBILITY_CHECKLIST.md
# Source data for this figure is in __SUBMISSION_PACKAGE/08_Source_Data/.
# ====================================================================

"""Adult AML single-cell: are the 12 panel-5E markers expressed in adult AML *blasts*, at what
% of cells, and across which cytogenetic subtypes?  Source = AML scAtlas (748k cells, 20 studies;
external_data/adult_scrna). Restricted to ADULT + myeloid-leukemia + MALIGNANT cells (van Galen/
Zeng hierarchy: LSC + *-like). Normalised the PEDIATRIC way (CP10k -> log1p) for mean expression;
% positive = raw counts > 0 (same convention as the pediatric atlas).
-> source_data/adult_scrna_{overall,subtype}.csv
"""
import os
import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
SD = os.path.join(ROOT, "source_data")
H5 = os.path.join(ROOT, "external_data", "adult_scrna", "AML_scAtlas_748k.h5ad")

# symbol -> (display, group)  -- the 12 panel-5E genes
GENES = {"CD96": ("CD96", "lead"), "CD9": ("CD9", "lead"), "SUCNR1": ("SUCNR1", "lead"),
         "IL1RAP": ("IL1RAP", "lead"), "FLT3": ("FLT3", "control"), "CD33": ("CD33", "clinical"),
         "IL3RA": ("CD123", "clinical"), "CLEC12A": ("CLL-1", "clinical"), "CD7": ("CD7", "ppac"),
         "TNFRSF4": ("TNFRSF4", "ppac"), "ABCA7": ("ABCA7", "ppac"), "ITGAX": ("ITGAX", "ppac")}
MALIG = {"LSC", "HSC-like", "Prog-like", "GMP-like", "cDC-like", "ProMono-like", "Mono-like"}


def eln_group(x):
    """ELN molecular-genetic classification -> consolidated subtype (richer than fusions)."""
    return {"NPM1": "NPM1", "NPM1 with FLT3-ITD": "NPM1/FLT3-ITD", "FLT3-ITD": "FLT3-ITD",
            "RUNX1-RUNX1T1": "t(8;21)", "CBFB-MYH11": "inv(16)/CBF", "inv(16)": "inv(16)/CBF",
            "TP53": "TP53", "MDS Related Genes": "MDS-related", "CK": "Complex karyotype",
            "-7": "−7/del(7q)", "MLL-AF9": "KMT2Ar", "MLLr-Other": "KMT2Ar",
            "CEBPA": "CEBPA"}.get(str(x), "Other")


A = ad.read_h5ad(H5, backed="r")
o = A.obs
mask = ((o["Main Age Group"] == "adult") & (o["disease"] == "myeloid leukemia")
        & (o["HSPC Cell Type"].isin(MALIG))).to_numpy()
idx = np.where(mask)[0]
print(f"adult AML malignant cells: {len(idx):,}  donors: {o.loc[mask, 'donor_id'].nunique()}")

Asub = A[idx].to_memory()
md = Asub.obs.copy()
md["subtype"] = md["ELN Classification"].map(eln_group)

# locate the 12 genes in raw.var (try var_names then feature_name)
rv = Asub.raw.var
names = rv["feature_name"].astype(str).values if "feature_name" in rv.columns else rv.index.astype(str).values
pos = {}
for g in GENES:
    hit = np.where(names == g)[0]
    if len(hit):
        pos[g] = int(hit[0])
    else:
        print(f"  WARNING gene absent from raw.var: {g}")
syms = list(pos)
raw = Asub.raw.X
raw = raw.tocsr() if sp.issparse(raw) else sp.csr_matrix(raw)
libsize = np.asarray(raw.sum(1)).ravel().astype(float); libsize[libsize == 0] = 1.0
G = raw[:, [pos[s] for s in syms]].toarray().astype(float)   # cells x 12 raw counts
POSm = G > 0
EXPR = np.log1p(G / libsize[:, None] * 1e4)                   # CP10k log1p (pediatric convention)


def summarise(grp_mask):
    n = int(grp_mask.sum())
    row = {"n_cells": n}
    for j, s in enumerate(syms):
        row[s + "_pct"] = round(100.0 * POSm[grp_mask, j].mean(), 2) if n else np.nan
        row[s + "_mean"] = round(float(EXPR[grp_mask, j].mean()), 4) if n else np.nan
    return row


def to_long(df_wide, keycol):
    rows = []
    for _, r in df_wide.iterrows():
        for s in syms:
            rows.append(dict(gene=s, label=GENES[s][0], group=GENES[s][1], **{keycol: r[keycol]},
                             pct=r[s + "_pct"], mean=r[s + "_mean"], n_cells=int(r["n_cells"])))
    return pd.DataFrame(rows)


# ---- overall ----
ov = pd.DataFrame([{**{"scope": "adult AML blasts"}, **summarise(np.ones(len(idx), bool))}])
pd.DataFrame([dict(gene=s, label=GENES[s][0], group=GENES[s][1],
                   pct=ov.iloc[0][s + "_pct"], mean=ov.iloc[0][s + "_mean"],
                   n_cells=int(ov.iloc[0]["n_cells"])) for s in syms]).to_csv(
    os.path.join(SD, "adult_scrna_overall.csv"), index=False)

# ---- by ELN molecular-genetic subtype ----
crows = []
for c, g in md.groupby("subtype", observed=True):
    m = (md["subtype"] == c).to_numpy()
    crows.append({"subtype": c, "n_donors": g["donor_id"].nunique(), **summarise(m)})
cdf = pd.DataFrame(crows)
to_long(cdf, "subtype").merge(cdf[["subtype", "n_donors"]], on="subtype").to_csv(
    os.path.join(SD, "adult_scrna_subtype.csv"), index=False)

print("\n=== OVERALL: % of adult AML blasts positive (raw>0) ===")
ovl = pd.read_csv(os.path.join(SD, "adult_scrna_overall.csv")).sort_values("pct", ascending=False)
print(ovl[["label", "group", "pct", "mean"]].to_string(index=False))
print("\n=== ELN molecular subtypes (cells / donors) ===")
print(cdf[["subtype", "n_cells", "n_donors"]].sort_values("n_donors", ascending=False).to_string(index=False))
print("\nwrote adult_scrna_overall.csv, adult_scrna_subtype.csv")
