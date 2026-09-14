#!/usr/bin/env python
# ====================================================================
# build_bulk_S17.py  |  Upstream builder — TARGET-bulk half of Figure S17
# Regenerates the BULK CD96-vs-comparator co-expression on the CURRENT 12-antigen
# panel (the metacell/single-cell half is build_metacells_S7.py). The previously
# deployed bulk used a stale comparator set (CD70/GPR56/MSLN) and was missing
# CD9/SUCNR1/TNFRSF4/ABCA7/ITGAX/CD7 -> half-empty bulk row + spurious NA column.
# ====================================================================
"""TARGET bulk: CD96 (x) vs each of the 11 panel comparators (y), log2 CPM, by
cytogenetic subtype. Pure data (TARGET STAR counts + clinical metadata); no
fabricated values. Writes:
  source_data/figS17_points_bulk.csv   (marker, cd96, val, subtype)
  source_data/figS17_corr_r.csv        (bulk rows replaced; single_cell rows kept)
"""
import os
import numpy as np, pandas as pd
from scipy import stats

RD = "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/07_Code/shared_CD96_pipeline/source_data"
FILT = "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/ALSF_AML_plot/H5AD/NCBI_TARGET_Star_counts_filtered.csv"
MERGED = "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/Figure7/Figure7F_TARGET_CD96_merged.csv"

# current 12-antigen panel (matches build_metacells_S7.py)
ENS = {"CD96":"ENSG00000153283","CD9":"ENSG00000010278","SUCNR1":"ENSG00000198829",
       "IL1RAP":"ENSG00000196083","FLT3":"ENSG00000122025","CD33":"ENSG00000105383",
       "IL3RA":"ENSG00000185291","CLEC12A":"ENSG00000172322","CD7":"ENSG00000173762",
       "TNFRSF4":"ENSG00000186827","ABCA7":"ENSG00000064687","ITGAX":"ENSG00000140678"}
GENES = list(ENS); OTHERS = [g for g in GENES if g != "CD96"]
GLABEL = {"CD96":"CD96","CD9":"CD9","SUCNR1":"SUCNR1","IL1RAP":"IL1RAP","FLT3":"FLT3",
          "CD33":"CD33","IL3RA":"CD123","CLEC12A":"CLL-1","CD7":"CD7","TNFRSF4":"TNFRSF4",
          "ABCA7":"ABCA7","ITGAX":"ITGAX"}

def bulk_sub(row):
    fz = str(row.get("Gene Fusion", "")).strip()
    if fz and fz.lower() != "nan":
        if fz.startswith("KMT2A-") or fz == "other MLL": return "KMT2Ar"
        if fz == "RUNX1-RUNX1T1": return "t(8;21)"
        if fz == "CBFB-MYH11": return "inv(16)"
        if fz.startswith("NUP98-"): return "NUP98r"
        return "Other"
    return {"Normal":"CN","t(8;21)":"t(8;21)","inv(16)":"inv(16)","PML-RARA":"PML-RARA",
            "MLL":"KMT2Ar"}.get(str(row.get("Primary Cytogenetic Code","")).strip(), "Other")

def main():
    raw = pd.read_csv(FILT).set_index("ensembl_gene_id")
    raw.index = raw.index.astype(str).str.split(".").str[0]
    g = raw.loc[[ENS[x] for x in GENES]].T
    g.columns = GENES
    m = pd.read_csv(MERGED).set_index("sample")
    ids = [s for s in g.index if s in m.index and pd.notna(m.loc[s, "lib_size"])]
    lib = m.loc[ids, "lib_size"].astype(float)
    d = np.log2(g.loc[ids].div(lib, axis=0) * 1e6 + 1)
    d["subtype"] = [bulk_sub(m.loc[s]) for s in ids]
    print(f"TARGET bulk samples: {len(d)}")

    rrows, pts = [], []
    for x in OTHERS:
        ok = np.isfinite(d["CD96"]) & np.isfinite(d[x])
        r = stats.pearsonr(d["CD96"][ok], d[x][ok])[0]
        rho = stats.spearmanr(d["CD96"][ok], d[x][ok])[0]
        rrows.append(dict(modality="bulk", marker=x, label=GLABEL[x],
                          pearson=round(r, 4), spearman=round(rho, 4), n=int(ok.sum())))
        pts.append(pd.DataFrame({"marker": GLABEL[x], "cd96": d["CD96"].values,
                                 "val": d[x].values, "subtype": d["subtype"].values}))
    pd.concat(pts, ignore_index=True).to_csv(os.path.join(RD, "figS17_points_bulk.csv"), index=False)

    # recompute the single-cell (metacell) r from the CURRENT metacell points too,
    # because the deployed figS17_corr_r.csv single_cell rows were still the stale
    # 8-comparator set (CD70/GPR56/MSLN) while the points were already current.
    scp = pd.read_csv(os.path.join(RD, "figS17_points_singlecell.csv"))
    label2gene = {v: k for k, v in GLABEL.items()}
    scrows = []
    for lab, sub in scp.groupby("marker"):
        ok = np.isfinite(sub["cd96"]) & np.isfinite(sub["val"])
        scrows.append(dict(modality="single_cell", marker=label2gene.get(lab, lab), label=lab,
                           pearson=round(stats.pearsonr(sub["cd96"][ok], sub["val"][ok])[0], 4),
                           spearman=round(stats.spearmanr(sub["cd96"][ok], sub["val"][ok])[0], 4),
                           n=int(ok.sum())))
    corr = os.path.join(RD, "figS17_corr_r.csv")
    pd.concat([pd.DataFrame(scrows), pd.DataFrame(rrows)], ignore_index=True).to_csv(corr, index=False)
    print("\nbulk r (current panel):")
    print(pd.DataFrame(rrows)[["label", "pearson", "spearman"]].to_string(index=False))
    allr = [abs(x["pearson"]) for x in rrows] + [abs(x["pearson"]) for x in scrows]
    print(f"\nmax |pearson| overall (sc {len(scrows)} + bulk {len(rrows)}): {max(allr):.2f}")

if __name__ == "__main__":
    main()
