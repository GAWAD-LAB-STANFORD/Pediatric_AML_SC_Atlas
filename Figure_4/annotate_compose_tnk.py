#!/usr/bin/env python
"""Figure 4 rebuild, step 2: annotate the re-derived T/NK subclusters by marker expression and compute
per-sample composition (AML vs healthy BM). Annotation from tnk_cluster_marker_means.csv (see comments).
The LYZ+ low-quality cluster (15) is dropped. Writes per-cell subset labels, per-sample subset fractions,
and a UMAP table for the panels. No fabricated values."""
import os, numpy as np, pandas as pd
D = "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/08_Source_Data/Figure_4/v2"
# cluster -> subset (from tnk_cluster_marker_means.csv):
# 0 CD3-lo GZMB/GNLY/FCGR3A+ -> GZMB NK ; 1 CD3+CD8-lo GNLY+ -> GZMB DNT ; 2 CD4 IL7R+ CCR7-lo -> Memory CD4 T
# 3 CD4 CCR7/SELL/TCF7+ -> Naive CD4 T ; 4 FOXP3/CTLA4+ -> Treg ; 5 CD8 CCR7/TCF7+ -> Naive CD8 T
# 6 CD8 GZMB/GNLY/PRF1+ -> GZMB CD8 T ; 7 CD3-lo cytotoxic -> GZMB NK ; 8 CD8 CCR7/TCF7+ -> Naive CD8 T
# 9,10 CD4 CCR7/TCF7+ -> Naive CD4 T ; 11 CD8 GZMK+GZMB-lo -> GZMK CD8 T ; 12 SLC4A10/KLRB1/TRAV1-2+ -> MAIT
# 13 CD3-lo GZMK/NCAM1/SELL+ -> GZMK NK ; 14 MKI67+ -> Proliferating T ; 15 LYZ+ -> Low-quality (drop) ; 16 CD8 IL7R+ -> Memory CD8 T
CMAP = {0:"GZMB NK",1:"GZMB DNT",2:"Memory CD4 T",3:"Naïve CD4 T",4:"Treg",5:"Naïve CD8 T",
        6:"GZMB CD8 T",7:"GZMB NK",8:"Naïve CD8 T",9:"Naïve CD4 T",10:"Naïve CD4 T",11:"GZMK CD8 T",
        12:"MAIT",13:"GZMK NK",14:"Proliferating T",15:"Low-quality",16:"Memory CD8 T"}
pc = pd.read_csv(os.path.join(D, "tnk_percell.csv"))
pc["subset"] = pc["tnk_leiden"].map(CMAP)
pc = pc[pc["subset"] != "Low-quality"].copy()
pc["group"] = np.where(pc["SampleID"].astype(str).str.startswith("0_HealthyBM"), "HBM", "AML")
pc.to_csv(os.path.join(D, "tnk_subset_percell.csv"), index=False)
# per-sample composition (fraction of each subset among that sample's T/NK cells)
comp = (pc.groupby(["SampleID","subset"]).size().rename("n").reset_index())
tot = comp.groupby("SampleID")["n"].transform("sum")
comp["frac"] = comp["n"] / tot
comp["group"] = np.where(comp["SampleID"].astype(str).str.startswith("0_HealthyBM"), "HBM", "AML")
comp.to_csv(os.path.join(D, "tnk_composition_bySample.csv"), index=False)
# cohort summary
n_samp = pc["SampleID"].nunique(); n_aml = pc.loc[pc.group=="AML","SampleID"].nunique()
print("T/NK cells (excl low-quality):", len(pc), "| samples:", n_samp, "(AML", n_aml, "+ HBM", n_samp-n_aml, ")")
print("\nsubset totals:"); print(pc["subset"].value_counts().to_string())
print("\nmedian subset fraction AML vs HBM:")
piv = comp.pivot_table(index="subset", columns="group", values="frac", aggfunc="median").fillna(0)
print(piv.round(3).to_string())
