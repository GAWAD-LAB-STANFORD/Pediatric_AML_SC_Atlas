#!/usr/bin/env python
# ====================================================================
# build_ls_headtohead.py  |  Upstream builder (regenerates source_data)
# Per-poor-prognosis-STATE (LS) efficacy-vs-toxicity (Figure-5C construction).
# Replaces the retired per-PPAC version (build_ppac_headtohead.py) on the 49-state scheme.
# ====================================================================
"""Per-poor-prognosis-STATE head-to-head (Figure-5C construction), six panels.

For each group G (the 5 multi-patient poor LS + "All poor LS" pooled) plot every
marker as one point in the same toxicity-vs-efficacy space as Figure 5C, with
efficacy made state-specific:
    x = leak     = % of normal HSPC / Myeloid_Pro positive (max of the two)   [toxicity]
    y = cov_pct  = % of G's cells positive                                    [efficacy, state-specific]
  size = rec_pct = % of contributing patients (>=20 cells in G) with >=20% of
                   their G cells positive                                     [recurrence]

Markers per panel = the TOP-5 state-specific markers for G (from funnel_ls_specific.csv,
the S15 ranking) UNION the 12 panel-E antigens (Figure 7B/C). The 12 panel antigens are
shown regardless of the toxicity gate (so the toxic clinical comparators FLT3/CD33/
CD123/CLL-1 appear in the toxic zone). Uses the corrected 49-state assignment (the "49
file"). No fabricated values.

Output -> source_data/figS_ls_headtohead.csv  (long: group x marker).
"""
import os, sys
import numpy as np, pandas as pd
import h5py
from scipy.sparse import csr_matrix

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/07_Code/Figure_3")
import config as C
RANK = os.path.join(HERE, "..", "source_data", "funnel_ls_specific.csv")
COMP = os.path.join(HERE, "..", "source_data", "funnel_ls_specific_composition.csv")
OUT  = os.path.join(HERE, "..", "source_data", "figS_ls_headtohead.csv")

REC_THRESH = 0.20
MIN_CELLS_PER_PT = 20
TOP_K = 5
# 12 panel-E antigens (Figure 7B/C): leads + FLT3 control + clinical + additional candidates
CLIN = {"CD96":"CD96","CD9":"CD9","SUCNR1":"SUCNR1","IL1RAP":"IL1RAP","FLT3":"FLT3",
        "CD33":"CD33","IL3RA":"CD123","CLEC12A":"CLL-1","CD7":"CD7","TNFRSF4":"TNFRSF4",
        "ABCA7":"ABCA7","ITGAX":"ITGAX"}

rank = pd.read_csv(RANK); comp = pd.read_csv(COMP).set_index("group")
GROUPS = comp.index.tolist()   # "All poor LS" + the 5 multi-patient poor LS (builder order)
top_by_group = {g: list(sub.sort_values("cov_pct", ascending=False).head(TOP_K)["gene"])
                for g, sub in rank.groupby("group")}
need = set(CLIN)
for g in GROUPS: need |= set(top_by_group.get(g, []))
need = sorted(need)

# ---- atlas: barcodes, Cell Type (gate), SampleID, raw counts for `need` genes ----
dec = lambda a: np.array([x.decode() if isinstance(x, bytes) else x for x in a])
f = h5py.File(C.H5, "r")
obs = dec(f["obs"]["_index"][:])
def catcol(n):
    g = f["obs"][n]
    return pd.Series(dec(g["categories"][:])[g["codes"][:]], index=obs) if isinstance(g, h5py.Group) \
           else pd.Series(dec(g[:]), index=obs)
ct = catcol("Cell Type"); samp = catcol("SampleID")
rv = dec(f["raw"]["var"]["_index"][:]); rpos = {g: i for i, g in enumerate(rv)}
present = [g for g in need if g in rpos]
missing = [g for g in need if g not in rpos]
if missing: print("WARNING not in matrix:", missing)
gi = [rpos[g] for g in present]
Xr = f["raw"]["X"]; shape = tuple(Xr.attrs["shape"]) if "shape" in Xr.attrs else (len(obs), len(rv))
X = csr_matrix((Xr["data"][:], Xr["indices"][:], Xr["indptr"][:]), shape=shape)
pos = (np.asarray(X[:, gi].todense()) > 0)
f.close()
present = np.array(present); col = {g: i for i, g in enumerate(present)}

# LS assignment (the "49 file")
ra = pd.read_csv(C.REFA); ra["barcode"] = ra["barcode"].astype(str).str.replace(r"^b'|'$", "", regex=True)
raws, ls = C.stable_raws(); raw2ls = {rw: ls[i] for i, rw in enumerate(raws)}
bc2ls = dict(zip(ra.barcode, ra["res4.50"].map(raw2ls)))
lscol = pd.Series([bc2ls.get(b) for b in obs], index=obs)

# poor-LS membership for the pooled group
prog = pd.read_csv("/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/08_Source_Data/Figure_2/v3/LS_wholecohort_prognosis.csv")
poor_all = prog[prog.prognosis == "poor"].LS.tolist()

hspc = (ct == "0_HSPC").values; mye = (ct == "Myeloid_Pro").values
hspc_pct = 100 * pos[hspc].mean(0); mye_pct = 100 * pos[mye].mean(0)

rows = []
for gname in GROUPS:
    gmask = lscol.isin(poor_all).values if gname == "All poor LS" else (lscol == gname).values
    n_cells = int(gmask.sum())
    contrib = [s for s in samp[gmask].value_counts().index
               if (gmask & (samp == s).values).sum() >= MIN_CELLS_PER_PT]
    spflag = bool(comp.loc[gname, "single_patient_flag"]) if gname in comp.index else False
    markers = list(dict.fromkeys(top_by_group.get(gname, []) + list(CLIN)))
    cov_all = 100 * pos[gmask].mean(0)
    for g in markers:
        if g not in col: continue
        i = col[g]; cov = float(cov_all[i])
        rec = 100 * np.mean([pos[gmask & (samp == s).values][:, i].mean() >= REC_THRESH
                             for s in contrib]) if contrib else 0.0
        rows.append(dict(group=gname, gene=g, label=CLIN.get(g, g),
                         cls="clinical" if g in CLIN else "discovered",
                         leak=round(float(max(hspc_pct[i], mye_pct[i])), 2),
                         cov_pct=round(cov, 2), rec_pct=round(float(rec), 2),
                         HSPC=round(float(hspc_pct[i]), 2), Myeloid=round(float(mye_pct[i]), 2),
                         in_top5=g in top_by_group.get(gname, []),
                         n_cells=n_cells, single_patient_flag=spflag))

df = pd.DataFrame(rows); df.to_csv(OUT, index=False)
print(f"wrote {OUT}  ({len(df)} rows, {df.group.nunique()} panels)")
pd.set_option("display.width", 170)
for gname in GROUPS:
    sub = df[df.group == gname].sort_values(["cls", "cov_pct"], ascending=[True, False])
    print(f"\n=== {gname} (n={sub.n_cells.iloc[0]}) ===")
    print(sub[["label","cls","in_top5","leak","cov_pct","rec_pct"]].to_string(index=False))
