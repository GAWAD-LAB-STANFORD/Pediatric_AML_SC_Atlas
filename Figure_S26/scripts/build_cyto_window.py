#!/usr/bin/env python
# ====================================================================
# build_cyto_window.py  |  Upstream builder — additive therapeutic window by CYTOGENETIC subtype
# Reworked on the corrected scheme: reads the main atlas (raw_counts.h5ad is GEO-only) and uses the
# corrected 4-type leukemic definition (AML/AML-MKI67/AML-PCNA/AML-CD1C).
# ====================================================================
"""Additive therapeutic window by CYTOGENETIC SUBTYPE (analogue of the poor-LS window, grouped by
cytogenetics). Each point = CD96 alone or CD96 OR one partner. For each point:
  toxicity = max over the normal safety compartments (0_HSPC, Myeloid_Pro) of the % covered;
  efficacy = % of that cytogenetic subtype's leukemic (4-type core-AML) cells covered.
Core-AML = Cell Type in {AML, AML-MKI67, AML-PCNA, AML-CD1C}; subtypes with >=500 core-AML cells kept;
'All AML' pools all core-AML. No fabricated values. Output -> source_data/figS_cyto_window.csv."""
import os, sys
import numpy as np, pandas as pd, h5py
from scipy.sparse import csr_matrix
sys.path.insert(0, "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/07_Code/Figure_3")
import config as C
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "source_data", "figS_cyto_window.csv")
GENES = ["CD96","CD9","SUCNR1","IL1RAP","FLT3","CD33","IL3RA","CLEC12A","CD7","TNFRSF4","ABCA7","ITGAX"]
GLABEL = {"CD96":"CD96","CD9":"CD9","SUCNR1":"SUCNR1","IL1RAP":"IL1RAP","FLT3":"FLT3","CD33":"CD33",
          "IL3RA":"CD123","CLEC12A":"CLL-1","CD7":"CD7","TNFRSF4":"TNFRSF4","ABCA7":"ABCA7","ITGAX":"ITGAX"}
AMLCORE = {"AML","AML-MKI67","AML-PCNA","AML-CD1C"}   # corrected 4-type
MIN = 500

dec = lambda a: np.array([x.decode() if isinstance(x, bytes) else x for x in a])
f = h5py.File(C.H5, "r"); obs = dec(f["obs"]["_index"][:])
def cat(n):
    g = f["obs"][n]
    return dec(g["categories"][:])[g["codes"][:]] if isinstance(g, h5py.Group) else dec(g[:])
ctype = cat("Cell Type"); cyto = cat("Cytogenetic")
rv = dec(f["raw"]["var"]["_index"][:]); rp = {x: i for i, x in enumerate(rv)}
gi = [rp[x] for x in GENES]
Xr = f["raw"]["X"]; shape = tuple(Xr.attrs["shape"]) if "shape" in Xr.attrs else (len(obs), len(rv))
X = csr_matrix((Xr["data"][:], Xr["indices"][:], Xr["indptr"][:]), shape=shape)
pos = (np.asarray(X[:, gi].todense()) > 0)
f.close()
cd96 = GENES.index("CD96")
amlcore = np.isin(ctype, list(AMLCORE))
hs = pos[ctype == "0_HSPC"]; my = pos[ctype == "Myeloid_Pro"]

def tox(cols):
    uh = np.zeros(hs.shape[0], bool); um = np.zeros(my.shape[0], bool)
    for k in cols: uh |= hs[:, k]; um |= my[:, k]
    return 100 * max(uh.mean(), um.mean())

rows = []
def win(group, mask, gtype):
    n = int(mask.sum())
    if n == 0: return
    P = pos[mask]
    rows.append(dict(group=group, group_type=gtype, target="CD96", label="CD96", add="CD96 alone",
                     toxicity=round(tox([cd96]), 2), efficacy=round(100*P[:, cd96].mean(), 2), n_cells=n))
    for k in range(len(GENES)):
        if k == cd96: continue
        rows.append(dict(group=group, group_type=gtype, target=GENES[k], label=GLABEL[GENES[k]],
                         add="+"+GLABEL[GENES[k]], toxicity=round(tox([cd96, k]), 2),
                         efficacy=round(100*(P[:, cd96] | P[:, k]).mean(), 2), n_cells=n))

win("All AML", amlcore, "pooled")
for cv in pd.unique(cyto[amlcore]):
    m = amlcore & (cyto == cv)
    if m.sum() >= MIN and "Healthy" not in str(cv):
        win(str(cv), m, "cyto")
pd.DataFrame(rows).to_csv(OUT, index=False)
d = pd.DataFrame(rows)
print(f"wrote {OUT}  ({d[d.group_type=='cyto'].group.nunique()} cyto subtypes >= {MIN} core-AML cells)")
print(d[d.group == "All AML"][["add","toxicity","efficacy"]].sort_values("toxicity").to_string(index=False))
