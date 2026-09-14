#!/usr/bin/env python
"""Figure 5 panel D (corrected scheme), BOTH sub-panels from the single-cell atlas (counts layer) so
they share one modality, one colour scale, one size scale (TARGET bulk validation is shown separately
in panels E/F). Corrected leukemic set = Cell_Type in config.LEUK_CELLTYPES AND SampleType=="AML".
  * fig5D_cyto.csv  — antigen coverage across cytogenetic subtypes (leukemic cells grouped by the atlas
                      'Cytogenetic' obs); n_pt = distinct patients per subtype (single-patient ones flagged).
  * fig5D_ls.csv    — antigen coverage across ALL 49 leukemic states (favorable->poor by whole-cohort
                      Cox HR; prognosis = sig favorable/poor at FDR<0.05 else n.s.) + normal HSPC/Myeloid_Pro.
Both: pct = % cells positive (counts>0), mean_expr = mean log-norm (counts layer, as Figure 3).
Gene symbols mapped CD123->IL3RA, CLL-1->CLEC12A. No fabricated values."""
import os, sys, numpy as np, pandas as pd, h5py
from scipy.sparse import csr_matrix
sys.path.insert(0, "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/07_Code/Figure_3")
import config as C

DEST = ["/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/07_Code/shared_CD96_pipeline/source_data",
        "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/08_Source_Data/Figure_5"]
PROG = "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/08_Source_Data/Figure_2/v3/LS_wholecohort_prognosis.csv"
GSET = ["CD96","CD9","SUCNR1","AMN","IL1RAP","LMBR1L","NRXN2","FURIN","UMODL1","CD7",
        "CD33","CD123","CLL-1","CD38","CD70","FLT3","TNFRSF4","ABCA7","ITGAX"]
SYM = {"CD123": "IL3RA", "CLL-1": "CLEC12A"}          # figure label -> gene symbol
MIN_CELLS = 30                                         # skip cytogenetic subtypes with too few leukemic cells

dec = lambda a: np.array([x.decode() if isinstance(x, bytes) else x for x in a])
def cat(f, k):
    o = f["obs"][k]; return dec(o["categories"][:])[o["codes"][:].astype(int)]

f = h5py.File(C.H5, "r")
obs = dec(f["obs"]["_index"][:]); var = dec(f["var"]["_index"][:])
ct = cat(f, "Cell Type"); st = cat(f, "SampleType"); sid = cat(f, "SampleID"); cyto = cat(f, "Cytogenetic")
gpos = {g: i for i, g in enumerate(var)}
present = [(lab, SYM.get(lab, lab)) for lab in GSET if SYM.get(lab, lab) in gpos]
gi = [gpos[sym] for _, sym in present]
grp = f["layers"]["counts"]
X = csr_matrix((grp["data"][:], grp["indices"][:], grp["indptr"][:]), shape=tuple(grp.attrs["shape"]))
f.close()
Xs = X[:, gi].tocsc()
leuk = np.isin(ct, C.LEUK_CELLTYPES) & (st == "AML")           # corrected leukemic set

def summarize(mask, group_name, extra=None):
    Xg = Xs[mask]
    pct = np.asarray((Xg > 0).mean(0)).ravel(); mean = np.asarray(Xg.mean(0)).ravel()
    out = []
    for j, (lab, _) in enumerate(present):
        row = [lab, group_name, round(100 * pct[j], 2), round(float(mean[j]), 4), lab]
        if extra: row += extra
        out.append(row)
    return out

# ---------------- LEFT: cytogenetic subtypes (single-cell, leukemic cells) ----------------
MPAL = {"BCR/ABL", "t(2;3)(p15;q26.2)", "t(7;14)(q21;q32)"}
cyto_groups = [g for g in pd.unique(cyto[leuk]) if isinstance(g, str) and (leuk & (cyto == g)).sum() >= MIN_CELLS]
rows = []
for g in cyto_groups:
    m = leuk & (cyto == g)
    npt = len(set(sid[m]))
    rows += summarize(m, g, extra=[npt, int(g in MPAL)])
cyto_df = pd.DataFrame(rows, columns=["gene", "group", "pct", "mean_expr", "label", "n_pt", "mpal"])

# ---------------- RIGHT: ALL 49 leukemic states (single-cell) + normal baseline ----------------
pg = pd.read_csv(PROG).sort_values("HR")
pg["sigdir"] = np.where(pg["fdr"] < 0.05, pg["dir"], "n.s.")
LS_ORDER = pg["LS"].tolist()
ls_dir = dict(zip(pg["LS"], pg["sigdir"])); ls_hr = dict(zip(pg["LS"], pg["HR"]))
ra = pd.read_csv(C.REFA); ra["barcode"] = ra["barcode"].astype(str).str.replace(r"^b'|'$", "", regex=True)
raws, ls = C.stable_raws(); raw2ls = {rw: ls[i] for i, rw in enumerate(raws)}
bc2ls = dict(zip(ra.barcode, ra["res4.50"].map(raw2ls)))
lsc = pd.Series([bc2ls.get(b) for b in obs], index=obs)
grpcol = pd.Series(index=obs, dtype=object)
grpcol[ct == "0_HSPC"] = "HSPC"; grpcol[ct == "Myeloid_Pro"] = "Myeloid_Pro"
for L in LS_ORDER:
    grpcol[lsc[lsc == L].index] = L
gc = grpcol.values
rows = []
for L in ["HSPC", "Myeloid_Pro"] + LS_ORDER:
    m = (gc == L)
    if m.sum() == 0: continue
    prog = "normal" if L in ("HSPC", "Myeloid_Pro") else ls_dir.get(L, "n.s.")
    hr = "" if L in ("HSPC", "Myeloid_Pro") else round(float(ls_hr[L]), 3)
    Xg = Xs[m]; pct = np.asarray((Xg > 0).mean(0)).ravel(); mean = np.asarray(Xg.mean(0)).ravel()
    for j, (lab, _) in enumerate(present):
        rows.append([lab, L, round(100 * pct[j], 2), round(float(mean[j]), 4), lab, prog, hr])
ls_df = pd.DataFrame(rows, columns=["gene", "group", "pct", "mean_expr", "label", "prognosis", "hr"])

for d in DEST:
    cyto_df.to_csv(os.path.join(d, "fig5D_cyto.csv"), index=False)
    ls_df.to_csv(os.path.join(d, "fig5D_ls.csv"), index=False)
nsig = int((pg["fdr"] < 0.05).sum())
print("cytogenetic subtypes (single-cell, >=%d leukemic cells): %d" % (MIN_CELLS, len(cyto_groups)))
print("  n patients per subtype:", {g: int(cyto_df[cyto_df.group == g].n_pt.iloc[0]) for g in cyto_groups})
print("ALL %d leukemic states written; %d whole-cohort EFS-significant (FDR<0.05)" % (len(LS_ORDER), nsig))
print("\nCD96 by cytogenetic subtype (%%pos / mean):")
print(cyto_df[cyto_df.gene == "CD96"][["group", "n_pt", "pct", "mean_expr"]].to_string(index=False))
