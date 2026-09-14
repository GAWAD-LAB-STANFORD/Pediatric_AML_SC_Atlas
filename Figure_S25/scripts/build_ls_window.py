#!/usr/bin/env python
# ====================================================================
# build_ls_window.py  |  Upstream builder — additive therapeutic window per poor-prognosis STATE (LS)
# Replaces the retired per-PAC version (build_pac_window.py) on the corrected 49-state scheme.
# ====================================================================
"""Additive therapeutic window per poor-prognosis leukemic state (LS), computed exactly like
the main Fig-7 draw_window: each point is CD96 alone or CD96 OR one partner. For each point:
  toxicity = max over the safety compartments (0_HSPC, Myeloid_Pro) of % cells positive for CD96 OR partner
  efficacy = % of the LS group's cells positive for CD96 OR partner
The R side sorts by toxicity and draws the Pareto frontier (cummax efficacy). Uses the corrected
49-state assignment (the "49 file"): reference_assignments res4.50 -> config.stable_raws.
Groups: "All poor LS" (all 8 poor pooled) + the 5 multi-patient poor states (LS_10/19/20/21/27);
the 3 single-patient poor states are pooled only. No fabricated values.
Output -> source_data/figS_ls_window.csv
"""
import os, sys
import numpy as np, pandas as pd, h5py
from scipy.sparse import csr_matrix
sys.path.insert(0, "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/07_Code/Figure_3")
import config as C
PROG = "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/08_Source_Data/Figure_2/v3/LS_wholecohort_prognosis.csv"
OUT  = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "source_data", "figS_ls_window.csv")
GENES = ["CD96","CD9","SUCNR1","IL1RAP","FLT3","CD33","IL3RA","CLEC12A","CD7","TNFRSF4","ABCA7","ITGAX"]
GLABEL = {"CD96":"CD96","CD9":"CD9","SUCNR1":"SUCNR1","IL1RAP":"IL1RAP","FLT3":"FLT3","CD33":"CD33",
          "IL3RA":"CD123","CLEC12A":"CLL-1","CD7":"CD7","TNFRSF4":"TNFRSF4","ABCA7":"ABCA7","ITGAX":"ITGAX"}
SINGLE_PT_FRAC = 0.77

prog = pd.read_csv(PROG); poor = prog[prog.prognosis == "poor"].copy()
poor["single_pt"] = poor.top_pt_frac >= SINGLE_PT_FRAC
poor_all = poor.LS.tolist()
poor_multi = poor[~poor.single_pt].sort_values("n_cells", ascending=False).LS.tolist()

dec = lambda a: np.array([x.decode() if isinstance(x, bytes) else x for x in a])
f = h5py.File(C.H5, "r"); obs = dec(f["obs"]["_index"][:])
g = f["obs"]["Cell Type"]; ct = dec(g["categories"][:])[g["codes"][:]]
rv = dec(f["raw"]["var"]["_index"][:]); rp = {x: i for i, x in enumerate(rv)}
gi = [rp[x] for x in GENES]
Xr = f["raw"]["X"]; shape = tuple(Xr.attrs["shape"]) if "shape" in Xr.attrs else (len(obs), len(rv))
X = csr_matrix((Xr["data"][:], Xr["indices"][:], Xr["indptr"][:]), shape=shape)
pos = (np.asarray(X[:, gi].todense()) > 0)
f.close()

ra = pd.read_csv(C.REFA); ra["barcode"] = ra["barcode"].astype(str).str.replace(r"^b'|'$", "", regex=True)
raws, ls = C.stable_raws(); raw2ls = {rw: ls[i] for i, rw in enumerate(raws)}
bc2ls = dict(zip(ra.barcode, ra["res4.50"].map(raw2ls)))
lscol = pd.Series([bc2ls.get(b) for b in obs], index=obs)
cd96 = GENES.index("CD96")
hs = pos[ct == "0_HSPC"]; my = pos[ct == "Myeloid_Pro"]

def tox(cols):
    uh = np.zeros(hs.shape[0], bool); um = np.zeros(my.shape[0], bool)
    for k in cols: uh |= hs[:, k]; um |= my[:, k]
    return 100 * max(uh.mean(), um.mean())

rows = []
def win(group, mask):
    n = int(mask.sum())
    if n == 0: return
    P = pos[mask]
    rows.append(dict(group=group, target="CD96", label="CD96", add="CD96 alone",
                     toxicity=round(tox([cd96]), 2), efficacy=round(100*P[:, cd96].mean(), 2), n_cells=n))
    for k in range(len(GENES)):
        if k == cd96: continue
        eff = 100*(P[:, cd96] | P[:, k]).mean()
        rows.append(dict(group=group, target=GENES[k], label=GLABEL[GENES[k]], add="+"+GLABEL[GENES[k]],
                         toxicity=round(tox([cd96, k]), 2), efficacy=round(eff, 2), n_cells=n))

win("All poor LS", lscol.isin(poor_all).values)
for g_ in poor_multi:
    win(g_, (lscol == g_).values)
pd.DataFrame(rows).to_csv(OUT, index=False)
print(f"wrote {OUT}  (groups: All poor LS + {poor_multi})")
d = pd.DataFrame(rows)
print("\nAll poor LS — CD96 + each partner (by toxicity):")
print(d[d.group == "All poor LS"][["add","toxicity","efficacy"]].sort_values("toxicity").to_string(index=False))
