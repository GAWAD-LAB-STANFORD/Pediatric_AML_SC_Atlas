#!/usr/bin/env python
# ====================================================================
# query_census_leads.py  |  CD96 figure pipeline component
# reproducibility header (auto-generated; edit freely below this block)
# Inputs : (reads upstream object / raw matrix — see body)
# Outputs: (writes to ../figures or ../source_data — see body)
# Run    : python query_census_leads.py   (or: bash run_all.sh)  |  Deps: see REPRODUCIBILITY_CHECKLIST.md
# Source data for this figure is in __SUBMISSION_PACKAGE/08_Source_Data/.
# ====================================================================

"""Body-wide NORMAL-tissue expression of the Figure-8 targets from the WHOLE
CZ CELLxGENE Census (all human normal datasets) -- NOT Tabula Sapiens alone.
RUN WITH THE ISOLATED CENSUS VENV, UNBUFFERED:
    ../.census_venv/bin/python -u query_census_normal_bodywide.py            # full sweep
    ../.census_venv/bin/python -u query_census_normal_bodywide.py spleen     # smoke test 1 tissue

Why the whole Census, not Tabula Sapiens only:
  The CELLxGENE "Gene Expression" web tool aggregates EVERY dataset in the
  Census, which is what a user sees on the website (e.g. IL3RA expressed in
  ~71% of the 93,489 plasmacytoid dendritic cells census-wide).  Restricting
  to the single Tabula Sapiens dataset is wrong here for two reasons:
    (1) IL3RA raw counts are corrupted/absent specifically in the TS Census
        build (0 reads across 82,785 TS blood cells, including pDCs that are
        canonically CD123-high) -- a TS-ingestion artifact, not biology;
    (2) the body-wide reference we want IS the cross-dataset aggregate the
        website shows.
  Pulling the whole Census fixes IL3RA and matches the website.

Method (no fabricated data; FULLY STREAMING -- never materializes a whole
tissue, so brain/blood at tens of millions of cells stay within memory):
  * resolve the 5 SYMBOLS (CD96, CD33, IL3RA, CLEC12A, FLT3) -> canonical
    Ensembl feature_ids (query BY feature_id; a symbol can map to >1 feature,
    e.g. pseudoautosomal IL3RA)
  * for each tissue_general, stream the obs table (is_primary_data == True,
    disease == 'normal') in Arrow BATCHES; for each batch, read ONLY that
    batch's raw-count X slice for the 5 genes by obs/var coordinates and
    accumulate per-cell-type positive counts.  Peak memory is one obs batch
    plus its X slice -- independent of tissue size.
  * restrict to 10x chemistries (detection rate is not comparable across 10x
    vs Smart-seq2; "% positive" needs one chemistry family)
  * "% positive" = fraction of cells with >=1 raw read (same convention used
    throughout Figure 7/8).  Every stored nonzero in the raw layer is a
    detection.
  * write three CSVs:  by tissue, by (tissue, cell_type), and a global
    per-cell-type roll-up (the website-style dotplot, for QC + callouts).
    by-tissue / by-celltype are flushed incrementally so a partial run still
    leaves usable output.

Caveat carried into the figure:  this is a DIFFERENT platform/processing than
the AML atlas, so it is presented as a standalone normal-tissue reference
(relative across organs / cell types), never numerically compared to the AML
coverage.
"""
import os
import sys
import csv
import time
import numpy as np
import pandas as pd
import cellxgene_census
import tiledbsoma as soma

HERE = os.path.dirname(os.path.abspath(__file__))
# ALTERNATE: pull ONLY the 3 NEW lead targets body-wide (CD96 + IL1RAP already in the
# 9-gene census_raw tables); merged with those downstream for the 5-lead Figure 7.
GENES = ["CD9", "SUCNR1", "AMN"]
CENSUS_VERSION = "stable"
MIN_GROUP = 20  # skip (tissue, cell_type) groups smaller than this (unstable %)
OUT_TISSUE = os.path.join(HERE, "Figure8_census_LEADS_by_tissue.csv")
OUT_CT = os.path.join(HERE, "Figure8_census_LEADS_by_celltype.csv")
OUT_CT_GLOBAL = os.path.join(HERE, "Figure8_census_LEADS_celltype_global.csv")


def resolve_canonical_feature_ids(census, genes):
    """symbol -> canonical Census feature_id (the copy measured in most cells)."""
    var = (census["census_data"]["homo_sapiens"].ms["RNA"].var
           .read(column_names=["feature_id", "feature_name", "n_measured_obs"])
           .concat().to_pandas())
    chosen = {}
    for g in genes:
        hits = var[var["feature_name"] == g]
        if hits.empty:
            print(f"  WARNING symbol absent from Census var: {g}", flush=True)
            continue
        best = hits.sort_values("n_measured_obs", ascending=False).iloc[0]
        chosen[g] = best["feature_id"]
        extra = "" if len(hits) == 1 else \
            f"  (chosen over {list(hits.loc[hits['feature_id'] != best['feature_id'], 'feature_id'])})"
        print(f"  '{g}' -> {best['feature_id']} "
              f"(n_measured_obs={int(best['n_measured_obs']):,}){extra}", flush=True)
    return chosen


def list_labels(census, category):
    # summary_cell_counts schema: organism is lower_snake ('homo_sapiens'),
    # category in {'all','assay','cell_type','disease','tissue','tissue_general',...}
    scc = census["census_info"]["summary_cell_counts"].read().concat().to_pandas()
    sub = scc[(scc["organism"] == "homo_sapiens") & (scc["category"] == category)]
    return sub[["label", "total_cell_count"]].copy()


def stream_tissue(exp, tissue, genes_present, varjoin_to_gidx, var_joinids):
    """Stream one tissue's obs in batches; for each batch read only its X slice
    (raw counts, 5 genes) by coordinate and accumulate per-cell-type positives.
    Peak memory ~ one obs batch.  10x cells only.

    Returns (ct_names, totals[ct], pos[ct, gene_present]) or None if no cells.
    """
    obs_filter = (f"tissue_general == '{tissue}' and "
                  f"is_primary_data == True and disease == 'normal'")
    # read the library-NORMALIZED layer: gives both detection (nonzero == raw>0,
    # so % positive is unchanged) AND the value needed for mean expression.
    Xnorm = exp.ms["RNA"].X["normalized"]
    ng = len(genes_present)

    name2code = {}
    tot_vec = np.zeros((0,), dtype=np.int64)
    pos_mat = np.zeros((0, ng), dtype=np.int64)
    sum_mat = np.zeros((0, ng), dtype=np.float64)   # sum of log1p(normalized) -> mean

    def grow(k):
        nonlocal tot_vec, pos_mat, sum_mat
        if k > tot_vec.shape[0]:
            add = k - tot_vec.shape[0]
            tot_vec = np.concatenate([tot_vec, np.zeros(add, dtype=np.int64)])
            pos_mat = np.vstack([pos_mat, np.zeros((add, ng), dtype=np.int64)])
            sum_mat = np.vstack([sum_mat, np.zeros((add, ng), dtype=np.float64)])

    # Sub-chunk every X read to a bounded number of cell coordinates.  An obs
    # Arrow batch for a giant tissue (brain/blood) can be millions of rows; a
    # single coordinate X read over millions of points spikes memory and is
    # slow, so we cap each X read at CELL_CHUNK cells.  Peak memory stays flat
    # regardless of tissue size.
    CELL_CHUNK = 250_000
    n_kept_total = 0

    def handle(joinids, names):
        nonlocal tot_vec, pos_mat, sum_mat, n_kept_total
        n_kept_total += len(joinids)
        # vectorized cell_type name -> stable global code; the only Python loop
        # is over the (few hundred) UNIQUE names in this chunk, not every cell.
        uniq, inv = np.unique(names, return_inverse=True)
        uniq_codes = np.empty(len(uniq), dtype=np.int64)
        for i, nm in enumerate(uniq):
            c = name2code.get(nm)
            if c is None:
                c = len(name2code)
                name2code[nm] = c
            uniq_codes[i] = c
        codes = uniq_codes[inv]
        grow(len(name2code))
        np.add.at(tot_vec, codes, 1)

        code_by_join = pd.Series(codes, index=joinids)
        for tbl in Xnorm.read(coords=(joinids, var_joinids)).tables():
            d0 = tbl.column("soma_dim_0").to_numpy()
            d1 = tbl.column("soma_dim_1").to_numpy()
            val = tbl.column("soma_data").to_numpy()
            cpe = code_by_join.reindex(d0).to_numpy()  # ct code per nonzero
            gidx = np.full(d1.shape, -1, dtype=np.int64)
            for vj, gi in varjoin_to_gidx.items():
                gidx[d1 == vj] = gi
            m = gidx >= 0
            if m.any():
                rr = cpe[m].astype(np.int64); cc = gidx[m]
                np.add.at(pos_mat, (rr, cc), 1)                       # detections
                # Census 'normalized' = fraction of library (sums to 1/cell); scale
                # to CP10k then log1p -> standard log-normalized expression (0-~6).
                np.add.at(sum_mat, (rr, cc), np.log1p(val[m] * 1e4))  # sum log1p(CP10k)

    for batch in exp.obs.read(
            value_filter=obs_filter,
            column_names=["soma_joinid", "assay", "cell_type"]):
        df = batch.to_pandas()
        is10x = df["assay"].astype(str).str.contains("10x", case=False)
        df = df.loc[is10x.values]
        if df.empty:
            continue
        joinids_all = df["soma_joinid"].to_numpy(dtype=np.int64)
        names_all = df["cell_type"].astype(str).to_numpy()
        for s in range(0, len(joinids_all), CELL_CHUNK):
            handle(joinids_all[s:s + CELL_CHUNK], names_all[s:s + CELL_CHUNK])

    if n_kept_total == 0:
        return None
    inv = [None] * len(name2code)
    for nm, c in name2code.items():
        inv[c] = nm
    return inv, tot_vec, pos_mat, sum_mat


def build_var_map(census, feat_for):
    """genes_present (ordered), {var_soma_joinid: gene_idx}, [var_soma_joinids]."""
    genes_present = [g for g in GENES if g in feat_for]
    fids = [feat_for[g] for g in genes_present]
    fid_list = "[" + ", ".join(f"'{f}'" for f in fids) + "]"
    var = (census["census_data"]["homo_sapiens"].ms["RNA"].var
           .read(value_filter=f"feature_id in {fid_list}",
                 column_names=["soma_joinid", "feature_id"])
           .concat().to_pandas())
    fid_to_gidx = {fids[i]: i for i in range(len(fids))}
    varjoin_to_gidx = {int(r.soma_joinid): fid_to_gidx[r.feature_id]
                       for r in var.itertuples() if r.feature_id in fid_to_gidx}
    return genes_present, varjoin_to_gidx, list(varjoin_to_gidx.keys())


def main():
    only_tissue = sys.argv[1] if len(sys.argv) > 1 else None
    print(f"  opening CZ CELLxGENE Census (version='{CENSUS_VERSION}') ...",
          flush=True)
    with cellxgene_census.open_soma(census_version=CENSUS_VERSION) as census:
        print("\n  resolving gene symbols -> canonical feature_ids ...", flush=True)
        feat_for = resolve_canonical_feature_ids(census, GENES)
        genes_present, varjoin_to_gidx, var_joinids = build_var_map(census, feat_for)
        gi_of = {g: i for i, g in enumerate(GENES)}

        exp = census["census_data"]["homo_sapiens"]

        if only_tissue:
            print(f"\n  SMOKE TEST -- single tissue: {only_tissue}", flush=True)
            res = stream_tissue(exp, only_tissue, genes_present,
                                varjoin_to_gidx, var_joinids)
            if res is None:
                print("  (no normal-primary-10x cells)", flush=True)
                return
            ct_names, totals, pos, ssum = res
            rows = []
            for c in range(len(ct_names)):
                tot_c = int(totals[c])
                rc = {"cell_type": ct_names[c], "n_cells": tot_c}
                for g in genes_present:
                    j = genes_present.index(g)
                    rc[f"{g}_pct"] = round(100.0 * pos[c, j] / tot_c, 2) if tot_c else np.nan
                    rc[f"{g}_mean"] = round(ssum[c, j] / tot_c, 4) if tot_c else np.nan
                rows.append(rc)
            sdf = pd.DataFrame(rows).sort_values("n_cells", ascending=False)
            print(f"\n  {only_tissue}: total kept 10x cells = {int(totals.sum()):,}, "
                  f"{len(ct_names)} cell types", flush=True)
            pdc = sdf[sdf["cell_type"].str.contains("plasmacytoid", case=False, na=False)]
            if not pdc.empty:
                print("  pDC QC:", flush=True)
                print(pdc.to_string(index=False), flush=True)
            print("\n  top cell types by n:", flush=True)
            print(sdf.head(15).to_string(index=False), flush=True)
            return

        # smallest-first: the 69 smaller tissues finish and flush to CSV before
        # the two giants (brain/blood); resilient to a late failure.  Final
        # aggregates are order-independent.
        tis_tbl = list_labels(census, "tissue_general").sort_values(
            "total_cell_count", ascending=True)
        tissues = [t for t in tis_tbl["label"].tolist() if t and t != "na"]
        print(f"\n  {len(tissues)} tissue_general categories to sweep "
              f"(smallest first; brain/blood last)", flush=True)

        ft = open(OUT_TISSUE, "w", newline="")
        wt = csv.writer(ft)
        wt.writerow(["tissue", "n_cells"] + [f"{g}_pct" for g in GENES]
                    + [f"{g}_mean" for g in GENES])
        fc = open(OUT_CT, "w", newline="")
        wc = csv.writer(fc)
        # store EXACT positive counts (for lossless re-aggregation, e.g. the
        # blood-cell-excluded organ atlas) alongside % positive and mean expr
        wc.writerow(["tissue", "cell_type", "n_cells"]
                    + [f"{g}_pos" for g in GENES]
                    + [f"{g}_pct" for g in GENES]
                    + [f"{g}_mean" for g in GENES])

        glob = {}  # cell_type -> {"tot": int, "pos": int[GENES], "sum": float[GENES]}

        skipped = []
        for ti, tis in enumerate(tissues, 1):
            res = "ERR"
            for attempt in range(1, 6):   # retry transient S3/network errors
                try:
                    res = stream_tissue(exp, tis, genes_present,
                                        varjoin_to_gidx, var_joinids)
                    break
                except Exception as e:
                    print(f"  [{ti}/{len(tissues)}] {tis}: attempt {attempt}/5 "
                          f"failed ({type(e).__name__}) -- retry in {15*attempt}s",
                          flush=True)
                    time.sleep(15 * attempt)
            if res == "ERR":
                print(f"  [{ti}/{len(tissues)}] {tis}: FAILED after 5 retries "
                      f"-- SKIPPED", flush=True)
                skipped.append(tis)
                continue
            if res is None:
                print(f"  [{ti}/{len(tissues)}] {tis}: 0 normal-primary-10x cells -- skip",
                      flush=True)
                continue
            ct_names, totals, pos, ssum = res
            n_kept = int(totals.sum())

            tissue_pos = pos.sum(axis=0)
            tissue_sum = ssum.sum(axis=0)
            trow = [tis, n_kept]
            for g in GENES:
                if g in genes_present:
                    j = genes_present.index(g)
                    trow.append(round(100.0 * tissue_pos[j] / n_kept, 3) if n_kept else np.nan)
                else:
                    trow.append("")
            for g in GENES:                                   # n-weighted tissue mean
                if g in genes_present:
                    j = genes_present.index(g)
                    trow.append(round(tissue_sum[j] / n_kept, 5) if n_kept else np.nan)
                else:
                    trow.append("")
            wt.writerow(trow)
            ft.flush()

            n_ct_written = 0
            for c in range(len(ct_names)):
                tot_c = int(totals[c])
                cname = ct_names[c]
                gd = glob.setdefault(
                    cname, {"tot": 0,
                            "pos": np.zeros(len(GENES), dtype=np.int64),
                            "sum": np.zeros(len(GENES), dtype=np.float64)})
                gd["tot"] += tot_c
                for g in genes_present:
                    j = genes_present.index(g)
                    gd["pos"][gi_of[g]] += int(pos[c, j])
                    gd["sum"][gi_of[g]] += float(ssum[c, j])
                if tot_c < MIN_GROUP:
                    continue
                crow = [tis, cname, tot_c]
                for g in GENES:  # exact integer positive counts
                    crow.append(int(pos[c, genes_present.index(g)])
                                if g in genes_present else "")
                for g in GENES:  # % positive
                    crow.append(round(100.0 * pos[c, genes_present.index(g)] / tot_c, 3)
                                if g in genes_present else "")
                for g in GENES:  # mean log1p(normalized) over all cells in group
                    crow.append(round(ssum[c, genes_present.index(g)] / tot_c, 5)
                                if g in genes_present else "")
                wc.writerow(crow)
                n_ct_written += 1
            fc.flush()

            qc_gene = GENES[0]   # progress QC on the first lead gene (IL3RA not in this run)
            il_t = ""
            if qc_gene in genes_present:
                j = genes_present.index(qc_gene)
                il_t = round(100.0 * tissue_pos[j] / n_kept, 2) if n_kept else "NA"
            pdc_str = ""
            for cname, gd in glob.items():
                if "plasmacytoid" in cname.lower() and gd["tot"]:
                    p = 100.0 * gd["pos"][gi_of[qc_gene]] / gd["tot"]
                    pdc_str = f"  [pDC running: {qc_gene}={p:.2f}% of {gd['tot']:,}]"
            print(f"  [{ti}/{len(tissues)}] {tis:22s} kept={n_kept:>9,}  "
                  f"{qc_gene}={il_t}%  ct_rows={n_ct_written}{pdc_str}", flush=True)

        ft.close()
        fc.close()

    grows = []
    for cname, gd in glob.items():
        tot = gd["tot"]
        rc = {"cell_type": cname, "n_cells": int(tot)}
        for g in GENES:
            rc[f"{g}_pct"] = round(100.0 * gd["pos"][gi_of[g]] / tot, 3) if tot else np.nan
        for g in GENES:
            rc[f"{g}_mean"] = round(gd["sum"][gi_of[g]] / tot, 5) if tot else np.nan
        grows.append(rc)
    gdf = (pd.DataFrame(grows).sort_values("n_cells", ascending=False)
           .reset_index(drop=True))
    gdf.to_csv(OUT_CT_GLOBAL, index=False)
    print(f"\n  wrote {OUT_TISSUE}")
    print(f"  wrote {OUT_CT}")
    print(f"  wrote {OUT_CT_GLOBAL}  ({len(gdf)} cell types, census-wide)")

    pdc = gdf[gdf["cell_type"].str.contains("plasmacytoid", case=False, na=False)]
    if not pdc.empty:
        print("\n  QC -- plasmacytoid dendritic cell (census-wide):")
        print(pdc.to_string(index=False))

    if skipped:
        print(f"\n  *** WARNING: {len(skipped)} tissue(s) SKIPPED after retries "
              f"(network/errors): {skipped}", flush=True)
        print("  *** RESULTS ARE INCOMPLETE -- re-run before using ***", flush=True)
    else:
        print("\n  all tissues processed (no skips).", flush=True)

    print("\nDone -> Figure8_census_normal_by_tissue.csv, "
          "Figure8_census_normal_by_celltype.csv, "
          "Figure8_census_normal_celltype_global.csv", flush=True)


if __name__ == "__main__":
    main()
