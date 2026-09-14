"""Figure 3 (CBF) shared config for the Python builders. Set FIG3_DATA to swap data."""
import os
BASE = "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster"
H5   = "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/ALSF_AML_plot/H5AD/ALSF_AML_Combo_3500_with_PAC_new.h5ad"
DATA = os.environ.get("FIG3_DATA", os.path.join(BASE, "08_Source_Data/Figure_3/v3"))
OUT  = os.environ.get("FIG3_OUT",  os.path.join(BASE, "04_Main_Figures/Figure_3_panels"))
MEANCPM = os.path.join(DATA, "meancpm.npy"); GENES = os.path.join(DATA, "genes.npy")
REFA    = os.path.join(DATA, "reference_assignments.csv")
JPC     = os.path.join(DATA, "jaccard_percluster.csv")
ANNOT   = os.path.join(DATA, "LS_annotation.csv")
CBF_ORDER = os.path.join(DATA, "LS_cbf_order.csv")   # written by 01_cbf_prognosis.R
LEUK_CELLTYPES = ["AML","AML-MKI67","AML-PCNA","AML-CD1C"]
RES = 4.5; JACC = 0.35
os.makedirs(OUT, exist_ok=True)

def stable_raws():
    """res-4.5 stable clusters (Jaccard>=0.35) sorted by n_cells -> LS_1..LS_K order."""
    import pandas as pd
    jp = pd.read_csv(JPC); jp = jp[jp.res == RES]
    stab = jp[jp.mean_jaccard >= JACC].sort_values("n_cells", ascending=False)
    raws = stab.cluster.astype(int).tolist()
    return raws, [f"LS_{i+1}" for i in range(len(raws))]
