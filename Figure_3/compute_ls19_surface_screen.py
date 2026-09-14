#!/usr/bin/env python
"""Data-driven Figure 3E: scan the surfaceome for candidate LS19 immunotherapy targets = surface genes
expressed in the HIGHEST fraction of LS19 cells while nearly ABSENT on normal HSPC (the therapeutic
window that spares normal stem/progenitors). Computes % expressing + mean log-norm for every surfaceome
gene across LS19 / normal HSPC / normal pDC / other leukemic, filters to HSPC-sparing genes, ranks by
LS19 % expressing, and writes the selected genes long-format for the dot plot. CD96 is retained and
flagged for reference. CSR-direct, no fabricated values."""
import os, sys, numpy as np, pandas as pd, h5py
from scipy.sparse import csr_matrix
sys.path.insert(0, "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/07_Code/Figure_3"); import config as C
D3 = C.DATA
SURF = "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/ref/surfaceome_genes.txt"
HSPC_MAX = 0.15     # "almost absent in HSC": <=15% of normal HSPC cells expressing
LS19_MIN = 0.25     # require reasonable prevalence in LS19
TOPN     = 14
surf = [l.strip() for l in open(SURF) if l.strip() and not l.startswith("#")]
f = h5py.File(C.H5, "r"); dec = lambda a: np.array([x.decode() if isinstance(x, bytes) else x for x in a])
obs = dec(f["obs"]["_index"][:]); var = dec(f["var"]["_index"][:])
ct = pd.Series(dec(f["obs"]["Cell_Type"]["categories"][:])[f["obs"]["Cell_Type"]["codes"][:]], index=obs)
gpos = {g: i for i, g in enumerate(var)}
present = [g for g in dict.fromkeys(surf) if g in gpos]; gi = [gpos[g] for g in present]
grp = f["layers"]["counts"]; X = csr_matrix((grp["data"][:], grp["indices"][:], grp["indptr"][:]), shape=tuple(grp.attrs["shape"]))
f.close()
ra = pd.read_csv(C.REFA); ra["barcode"] = ra["barcode"].astype(str).str.replace(r"^b'|'$", "", regex=True)
raws, ls = C.stable_raws(); raw2ls = {rw: ls[i] for i, rw in enumerate(raws)}
bc2ls = dict(zip(ra.barcode, ra["res4.50"].map(raw2ls)))
def group(b):
    v = bc2ls.get(b)
    if v == "LS_19": return "LS19 (poor)"
    if isinstance(v, str) and v.startswith("LS_"): return "other leukemic"
    if ct.get(b) == "pDC": return "normal pDC"
    if ct.get(b) == "0_HSPC": return "normal HSPC"
    return None
g = np.array([group(b) for b in obs])
order = ["normal HSPC", "LS19 (poor)", "normal pDC", "other leukemic"]
Xs = X[:, gi].tocsc()                                   # cells x surface-genes
stat = {}
for grp_ in order:
    m = g == grp_
    Xg = Xs[m]
    stat[(grp_, "pct")]  = np.asarray((Xg > 0).mean(0)).ravel()
    stat[(grp_, "mean")] = np.asarray(Xg.mean(0)).ravel()
tab = pd.DataFrame({"gene": present})
for grp_ in order:
    tab[grp_ + "_pct"]  = stat[(grp_, "pct")]
    tab[grp_ + "_mean"] = stat[(grp_, "mean")]
tab["window"] = tab["LS19 (poor)_pct"] - tab["normal HSPC_pct"]
# HSPC-sparing, LS19-prevalent candidates, ranked by LS19 % expressing
cand = tab[(tab["normal HSPC_pct"] <= HSPC_MAX) & (tab["LS19 (poor)_pct"] >= LS19_MIN)] \
         .sort_values("LS19 (poor)_pct", ascending=False)
sel = list(cand["gene"].head(TOPN))
if "CD96" not in sel and "CD96" in set(tab["gene"]): sel.append("CD96")   # keep the lead target for reference
# long format for the dot plot
rows = []
for gene in sel:
    r = tab[tab.gene == gene].iloc[0]
    for grp_ in order:
        rows.append([gene, grp_, r[grp_ + "_mean"], r[grp_ + "_pct"], r["window"]])
out = pd.DataFrame(rows, columns=["gene", "group", "mean", "pct", "window"])
out.to_csv(os.path.join(D3, "LS19_surface_screen.csv"), index=False)
print("surfaceome genes scanned:", len(present), "| candidates (HSPC<=%.0f%%, LS19>=%.0f%%): %d"
      % (HSPC_MAX*100, LS19_MIN*100, len(cand)))
print("\nselected (LS19%% expr / HSPC%% expr / window):")
for gene in sel:
    r = tab[tab.gene == gene].iloc[0]
    flag = "  <-- CD96" if gene == "CD96" else ""
    print(f"  {gene:10} LS19={r['LS19 (poor)_pct']*100:5.1f}%  HSPC={r['normal HSPC_pct']*100:4.1f}%  window={r['window']*100:5.1f}%{flag}")
print("\nCD96 rank by LS19%% among HSPC-sparing candidates:",
      (cand['gene'].tolist().index("CD96")+1) if "CD96" in cand['gene'].tolist() else "not in candidate set")
