#!/usr/bin/env python
# ====================================================================
# query_census_bladder_eye_CD96.py  |  CD96 figure pipeline component
# reproducibility header (auto-generated; edit freely below this block)
# Inputs : (reads upstream object / raw matrix — see body)
# Outputs: (writes to ../figures or ../source_data — see body)
# Run    : python query_census_bladder_eye_CD96.py   (or: bash run_all.sh)  |  Deps: see REPRODUCIBILITY_CHECKLIST.md
# Source data for this figure is in __SUBMISSION_PACKAGE/08_Source_Data/.
# ====================================================================

"""Per-cell CD96 expression in the two organs where CD96 has appreciable normal
expression (urinary bladder, eye), for the major cell types -> violin source.

Light pull: enumerate obs joinids per (tissue, cell_type), deterministically
downsample at the ID level (<= N_MAX per cell type), then read only that small
CD96 slice of the library-NORMALIZED layer (zeros implied) and convert to
log1p(CP10k).  No streaming of millions of X rows.  No fabricated values.
RUN: ../.census_venv/bin/python -u query_census_bladder_eye_CD96.py
"""
import os, numpy as np, pandas as pd, cellxgene_census

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "Figure8_bladder_eye_CD96_percell.csv")
CENSUS_VERSION = "stable"
GENE = "CD96"
N_MAX = 2000   # cells plotted per cell type (deterministic stride; stats on all)
SEL = {
    "urinary bladder": ["umbrella cell of urothelium", "intermediate cell of urothelium",
                        "basal cell of urothelium", "fibroblast", "smooth muscle cell",
                        "endothelial cell", "T cell"],
    "eye": ["retinal pigment epithelial cell", "non-pigmented ciliary epithelial cell",
            "retinal rod cell", "Mueller cell", "GABAergic amacrine cell",
            "fibroblast", "melanocyte"],
}

with cellxgene_census.open_soma(census_version=CENSUS_VERSION) as census:
    print("  resolving CD96 feature_id ...", flush=True)
    var = (census["census_data"]["homo_sapiens"].ms["RNA"].var
           .read(column_names=["soma_joinid", "feature_id", "feature_name", "n_measured_obs"])
           .concat().to_pandas())
    hits = var[var.feature_name == GENE].sort_values("n_measured_obs", ascending=False)
    vj = int(hits.iloc[0].soma_joinid)
    print(f"  CD96 -> var soma_joinid {vj}", flush=True)
    exp = census["census_data"]["homo_sapiens"]
    Xn = exp.ms["RNA"].X["normalized"]
    rows = []
    for tis, cts in SEL.items():
        print(f"  reading obs for {tis} ...", flush=True)
        ids_by_ct = {ct: [] for ct in cts}
        for batch in exp.obs.read(
                value_filter=(f"tissue_general == '{tis}' and is_primary_data == True "
                              f"and disease == 'normal'"),
                column_names=["soma_joinid", "cell_type", "assay"]):
            b = batch.to_pandas()
            b = b[b["assay"].astype(str).str.contains("10x", case=False, na=False)]
            for ct in cts:
                v = b.loc[b["cell_type"] == ct, "soma_joinid"].to_numpy()
                if len(v):
                    ids_by_ct[ct].append(v)
        for ct in cts:
            ids = np.concatenate(ids_by_ct[ct]) if ids_by_ct[ct] else np.array([], dtype=np.int64)
            n = len(ids)
            if n == 0:
                print(f"    MISSING {tis} / {ct}", flush=True); continue
            ids = np.sort(ids)
            if n > N_MAX:
                step = n // N_MAX
                ids_ds = ids[::step][:N_MAX]
            else:
                ids_ds = ids
            vals = {}
            for tbl in Xn.read(coords=(ids_ds.tolist(), [vj])).tables():
                d0 = tbl.column("soma_dim_0").to_numpy()
                dv = tbl.column("soma_data").to_numpy()
                for jid, val in zip(d0, dv):
                    vals[int(jid)] = float(val)
            for jid in ids_ds:
                v = vals.get(int(jid), 0.0)
                rows.append((tis, ct, round(float(np.log1p(v * 1e4)), 5)))
            pos = 100.0 * np.mean([vals.get(int(j), 0.0) > 0 for j in ids_ds])
            print(f"    {tis} / {ct}: n_total={n} plotted={len(ids_ds)} %pos~{pos:.1f}", flush=True)
    df = pd.DataFrame(rows, columns=["organ", "cell_type", "cd96"])
    df.to_csv(OUT, index=False)
    print(f"  wrote {OUT} ({len(df)} cells)", flush=True)
