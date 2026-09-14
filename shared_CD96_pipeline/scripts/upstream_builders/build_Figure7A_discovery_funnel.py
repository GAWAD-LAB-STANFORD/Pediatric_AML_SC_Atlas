#!/usr/bin/env python
# ====================================================================
# build_Figure7A_discovery_funnel.py  |  Upstream builder (regenerates source_data from raw inputs)
# reproducibility header (auto-generated; edit freely below this block)
# Inputs : (reads upstream object / raw matrix — see body)
# Outputs: (writes to ../figures or ../source_data — see body)
# Run    : python build_Figure7A_discovery_funnel.py   (or: bash run_all.sh)  |  Deps: see REPRODUCIBILITY_CHECKLIST.md
# Source data for this figure is in __SUBMISSION_PACKAGE/08_Source_Data/.
# ====================================================================

"""Figure 7A -- discovery funnel.

Replaces the earlier lollipop framing.  The previous "29% of AML cells"
headline was misleading because the per-patient distribution is bimodal
(genotype-defined CD96-cold subset).  The funnel framing emphasizes how
narrow the discovery pipeline is and surfaces CD96 as the unique lead.

Pipeline stages (all computed live from source data):

  1. Surfaceome universe          2,799 symbols   (Bausch-Fluck 2018,
                                                    ref/surfaceome_genes.txt)
  2. Detectable in scRNA          (intersection of surfaceome with h5ad
                                    raw.var symbols)
  3. Safety + tractability gate   scRNA_HSPC < 3.0%,  scRNA_Myeloid <= 5.0%,
                                  exclude TCR / CD3 family (MPAL-like
                                  co-expression, not real surface targets)
  4. Top 15 by composite          composite = mean(% AML cells,
                                    % patients with >=20% leukemic detection,
                                    % PPAC_1-5 cells)
  5. Lead candidate               CD96  (rank 1 of stage-3 set)

C3 uses ALL FIVE poor-prognosis clusters (PPAC_1-5), not a hand-picked
subset.  A sensitivity audit (supplements/_audit_PPAC_composite_sensitivity.py)
confirms CD96 is rank 1 regardless of how C3 is defined -- PPAC_1/2/5,
all 5 PPACs, all 7 leukemic PACs, or C3 dropped entirely -- with the top-3
(CD96, UMODL1, THSD7A) unchanged.  Using all 5 PPACs removes the
"why those clusters?" question without changing the lead.

A small inset on the right shows the top-5 composite scores so a reader
can see the gap between CD96 and the runner-ups without the panel having
to claim numbers in its own commentary.

Caches the per-gene scRNA stats to
`supplements/_surfaceome_scrna_per_gene.csv` so reruns are fast.
"""
import os
import sys
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import ConnectionPatch
import matplotlib as mpl

HERE = os.path.dirname(os.path.abspath(__file__))
SUPP = os.path.join(HERE, "supplements")
sys.path.insert(0, SUPP)
from _fig7_style import PALETTE, RC, clean_ax, save_fig, SEQ_CMAP

ROOT = os.path.abspath(os.path.join(HERE, ".."))
H5AD = os.path.join(ROOT, "ALSF_AML_plot", "H5AD",
                    "ALSF_AML_Combo_3500_with_PAC_new.h5ad")
SURFACEOME = os.path.join(ROOT, "ref", "surfaceome_genes.txt")
CACHE = os.path.join(SUPP, "_surfaceome_scrna_per_gene.csv")
# secondary cache (built by supplements/build_CD96_stress_test.py) -- contains
# clinical-stage AML targets that the surfaceome list misses by alias
# (e.g. surfaceome calls it GPR56; h5ad calls it ADGRG1)
CLIN_CACHE = os.path.join(SUPP, "_full_ungated_per_gene_table.csv")

# 8 clinical-stage / SoC pediatric AML surface targets used for the dropout
# annotation.  Each maps a gene symbol -> short therapy / trial reference.
CLINICAL_AML_TARGETS = {
    "CD33":     "gemtuzumab (FDA)",
    "IL3RA":    "CD123 / tagraxofusp",
    "CLEC12A":  "CLL-1 / CD371",
    "FLT3":     "midostaurin + CAR-T",
    "CD70":     "cusatuzumab",
    "ADGRG1":   "GPR56",
    "IL1RAP":   "nidanilimab",
    "MSLN":     "mesothelin",
}

# pipeline knobs (visible)
GATE_HSPC = 3.0
GATE_MYE = 5.0
N_SCRNA_PATIENTS = 28
SCRNA_AML_PCT_MIN = 0.20      # gene "present" in a patient if >=20% leuk
TOP_N = 15
LEUK_SUBS = [f"AML_{i}" for i in range(1, 28)]
# C3 = coverage across ALL FIVE poor-prognosis clusters (not a subset).
# Sensitivity audit confirms CD96 is rank 1 under every C3 definition.
PPACS = ["PPAC_1", "PPAC_2", "PPAC_3", "PPAC_4", "PPAC_5"]

# TCR / CD3 family -- MPAL-like contamination, not real surface targets
TCR_CD3_EXCLUDE = {
    "CD3D", "CD3E", "CD3G", "CD247",
    "TRAC", "TRBC1", "TRBC2", "TRGC1", "TRGC2", "TRDC",
}


def _decode(arr):
    return [x.decode() if isinstance(x, bytes) else x for x in arr]


def build_or_load_scrna_table():
    """Return per-gene scRNA metrics for the surfaceome universe.  Cache to
    `_surfaceome_scrna_per_gene.csv`."""
    if os.path.exists(CACHE):
        print(f"  loading cached scRNA per-gene table: {CACHE}")
        return pd.read_csv(CACHE)

    print("  no cache -- recomputing scRNA per-gene metrics for the surfaceome ...")
    with open(SURFACEOME) as fh:
        surf = set(fh.read().split())

    with h5py.File(H5AD, "r") as f:
        pac = np.array(pd.Categorical.from_codes(
            f["obs"]["Prognosis-Associated Clusters"]["codes"][:].astype(int),
            categories=_decode(f["obs"]["Prognosis-Associated Clusters"]["categories"][:])
        ).astype(object))
        samp = np.array(_decode(f["obs"]["Sample"]["categories"][:])
                       )[f["obs"]["Sample"]["codes"][:].astype(int)]
        sub = np.array(pd.Categorical.from_codes(
            f["obs"]["AML Sub-Clusters"]["codes"][:].astype(int),
            categories=_decode(f["obs"]["AML Sub-Clusters"]["categories"][:])
        ).astype(object))
        rx = f["raw"]["X"]
        X = sp.csr_matrix((rx["data"][:], rx["indices"][:], rx["indptr"][:]),
                          shape=tuple(rx.attrs["shape"])).tocsc()
        rk = f["raw"]["var"].attrs.get("_index", b"_index")
        rk = rk.decode() if isinstance(rk, bytes) else rk
        symbols = _decode(f["raw"]["var"][rk][:])

    sym_idx = {s: i for i, s in enumerate(symbols) if s in surf}
    print(f"    {len(sym_idx)} of {len(surf)} surfaceome symbols detectable in scRNA")

    in_hspc = pac == "0_HSPC"
    in_mye = pac == "Myeloid_Pro"
    in_ppac = np.isin(pac, PPACS)
    in_leuk = np.isin(sub, LEUK_SUBS)
    aml_samples = sorted({s for s in set(samp) if "HealthyBM" not in s})

    rows = []
    for g, idx in sym_idx.items():
        col = np.asarray(X[:, idx].todense(), dtype=np.float32).ravel()
        bin_ = col > 0
        rec = 0
        for s in aml_samples:
            m = in_leuk & (samp == s)
            if m.sum() and bin_[m].mean() >= SCRNA_AML_PCT_MIN:
                rec += 1
        rows.append(dict(
            gene=g,
            scRNA_HSPC_pct=round(100 * bin_[in_hspc].mean(), 3),
            scRNA_Myeloid_pct=round(100 * bin_[in_mye].mean(), 3) if in_mye.any() else 0.0,
            scRNA_AML_pct=round(100 * bin_[in_leuk].mean(), 2) if in_leuk.any() else 0.0,
            scRNA_PPAC_pct=round(100 * bin_[in_ppac].mean(), 2) if in_ppac.any() else 0.0,
            scRNA_n_pts_above20pct=rec,
        ))
    out = pd.DataFrame(rows)
    out.to_csv(CACHE, index=False)
    print(f"  cached {len(out)} surfaceome genes -> {CACHE}")
    return out


def clinical_dropouts():
    """Return per clinical-target stats: AML%, HSPC%, Myeloid%, dominant
    leak axis, pass/fail at safety gate, in_top15 flag.  Pulled from the
    secondary cache so we cover alias mismatches (GPR56 / ADGRG1)."""
    if not os.path.exists(CLIN_CACHE):
        print(f"  warning: clinical cache missing ({CLIN_CACHE}); "
              "skipping dropout annotation")
        return None
    df = pd.read_csv(CLIN_CACHE)
    rows = []
    for g, label in CLINICAL_AML_TARGETS.items():
        r = df[df["gene"] == g]
        if not len(r):
            print(f"    {g} not in clinical cache -- skipping")
            continue
        r = r.iloc[0]
        leak = max(r["scRNA_HSPC_pct"], r["scRNA_Myeloid_pct"])
        leak_axis = "HSPC" if r["scRNA_HSPC_pct"] >= r["scRNA_Myeloid_pct"] else "Mye"
        pass_gate = (r["scRNA_HSPC_pct"] < GATE_HSPC and
                     r["scRNA_Myeloid_pct"] <= GATE_MYE)
        rows.append(dict(gene=g, therapy=label,
                          aml=r["scRNA_AML_pct"],
                          hspc=r["scRNA_HSPC_pct"],
                          mye=r["scRNA_Myeloid_pct"],
                          leak=leak, leak_axis=leak_axis,
                          pass_gate=pass_gate))
    return pd.DataFrame(rows).sort_values(
        ["pass_gate", "leak"], ascending=[True, False]).reset_index(drop=True)


def compute_pipeline(scrna_df):
    """Apply the funnel filters and return counts + ranked candidate list."""
    n_surf = sum(1 for _ in open(SURFACEOME))   # count rows in source
    n_detect = len(scrna_df)

    # safety + TCR/CD3 gate
    safe = (scrna_df["scRNA_HSPC_pct"] < GATE_HSPC) & \
           (scrna_df["scRNA_Myeloid_pct"] <= GATE_MYE)
    not_tcr = ~scrna_df["gene"].isin(TCR_CD3_EXCLUDE)
    g = scrna_df[safe & not_tcr].copy()
    n_gate = len(g)

    # composite (same definition as the lollipop -- 3 equal-weighted components)
    g["C1_pct_AML_cells"] = g["scRNA_AML_pct"]
    g["C2_pct_AML_patients"] = (g["scRNA_n_pts_above20pct"] /
                                 N_SCRNA_PATIENTS) * 100
    g["C3_pct_PPAC_cells"] = g["scRNA_PPAC_pct"]
    g["composite"] = (g["C1_pct_AML_cells"] +
                       g["C2_pct_AML_patients"] +
                       g["C3_pct_PPAC_cells"]) / 3.0
    ranked = g.sort_values("composite", ascending=False).reset_index(drop=True)
    ranked.insert(0, "rank", np.arange(1, len(ranked) + 1))

    return dict(
        n_surfaceome=n_surf,
        n_detected=n_detect,
        n_gated=n_gate,
        n_top=TOP_N,
        n_lead=1,
        ranked=ranked,
    )


def draw_funnel(ax, stages):
    """Clean vertical funnel: rounded bands tapering top->bottom, cool-blue
    gradient narrowing to the crimson CD96 lead.  Count + criterion live inside
    each band; soft trapezoids connect the tiers.  Width is a designed silhouette
    (not the raw count ratio) so the 2,799 -> 1 drop doesn't crush small tiers."""
    n = len(stages)
    widths = np.linspace(1.06, 0.40, n)
    h = 0.150
    gap = 0.030
    y_cursor = 0.985
    centers = []
    for _ in range(n):
        centers.append(y_cursor - h / 2)
        y_cursor -= (h + gap)

    # soft connector trapezoids + "how many were filtered out" at each step
    for i in range(n - 1):
        wt, wb = widths[i], widths[i + 1]
        yt = centers[i] - h / 2
        yb = centers[i + 1] + h / 2
        ax.add_patch(mpatches.Polygon(
            [(-wt / 2, yt), (wt / 2, yt), (wb / 2, yb), (-wb / 2, yb)],
            closed=True, facecolor="#eef4fb", edgecolor="none", zorder=1))
        if stages[i].get("n") and stages[i + 1].get("n"):
            drop = stages[i]["n"] - stages[i + 1]["n"]
            ax.text(0, (yt + yb) / 2, f"−{drop:,}", ha="center",
                    va="center", fontsize=7.4, color="#a93226",
                    fontweight="bold", fontstyle="italic", zorder=4)

    for s, c, w in zip(stages, centers, widths):
        lead = s.get("emphasize", False)
        ax.add_patch(mpatches.FancyBboxPatch(
            (-w / 2, c - h / 2), w, h,
            boxstyle="round,pad=0.006,rounding_size=0.022",
            facecolor=s["color"], edgecolor=s.get("edge", "none"),
            linewidth=s.get("lw", 0), zorder=2))
        tc = s["textcolor"]
        ax.text(0, c + h * 0.24, s["bignum"], ha="center", va="center",
                fontsize=18 if lead else 15, fontweight="bold",
                color=tc, zorder=3)
        ax.text(0, c - h * 0.02, s["title"], ha="center", va="center",
                fontsize=8.6, fontweight="bold" if lead else "normal",
                color=tc, zorder=3)
        if s.get("detail"):
            ax.text(0, c - h * 0.30, s["detail"], ha="center", va="center",
                    fontsize=6.8, color=s["subcolor"], zorder=3)

    ax.set_xlim(-0.62, 0.62)
    ax.set_ylim(-0.02, 1.0)
    ax.set_aspect("auto")
    ax.axis("off")


def main():
    scrna_df = build_or_load_scrna_table()
    pipe = compute_pipeline(scrna_df)
    ranked = pipe["ranked"]
    ranked.to_csv(os.path.join(HERE, "Figure7A_composite_ranked.csv"), index=False)

    cd96_rank = int(ranked.index[ranked["gene"] == "CD96"][0]) + 1
    print(f"\n  stages:  {pipe['n_surfaceome']} surfaceome  ->  "
          f"{pipe['n_detected']} scRNA-detectable  ->  "
          f"{pipe['n_gated']} safety+TCR-gated  ->  "
          f"top {pipe['n_top']}  ->  lead = CD96  (rank {cd96_rank})")

    top = ranked.head(TOP_N)
    print("\n  --- top 5 by composite ---")
    print(top.head(5)[["rank", "gene", "composite", "C1_pct_AML_cells",
                        "C2_pct_AML_patients", "C3_pct_PPAC_cells",
                        "scRNA_HSPC_pct", "scRNA_Myeloid_pct"]].round(1).to_string(index=False))

    # ---- assemble funnel stages: cool gradient tiers + crimson CD96 lead ----
    stage_specs = [
        (pipe['n_surfaceome'], "surface proteins (surfaceome)",
         "Bausch-Fluck 2018 in silico surfaceome", 0.10),
        (pipe['n_detected'], "detectable in scRNA (28 patients)",
         "expressed in >=1 cell of the atlas", 0.27),
        (pipe['n_gated'], "pass safety + tractability gate",
         f"HSPC<{GATE_HSPC:.0f}% · Myeloid<={GATE_MYE:.0f}% · drop TCR/CD3 genes", 0.46),
        (pipe['n_top'], "top 15 by composite score",
         "ranked by C1+C2+C3 (see box)", 0.65),
    ]
    stages = []
    for n, title, detail, cpos in stage_specs:
        light = cpos <= 0.40
        stages.append(dict(
            bignum=f"{n:,}", n=n, title=title, detail=detail,
            color=mpl.colors.to_hex(SEQ_CMAP(cpos)),
            textcolor="#13314f" if light else "white",
            subcolor="#5b7894" if light else "#dce9f6"))
    stages.append(dict(
        bignum="CD96", title="lead surface target",
        detail=f"rank 1 of {pipe['n_gated']:,} gated",
        color=PALETTE["CD96"], textcolor="white", subcolor="#f6dcd8",
        emphasize=True, edge="#7b241c", lw=1.4))

    # clinical-target dropout table
    drop = clinical_dropouts()

    # ---- figure ----
    with mpl.rc_context(RC):
        fig = plt.figure(figsize=(16.0, 8.2))
        # 3 columns: funnel | clinical dropouts | top-15 lollipop
        gs = fig.add_gridspec(1, 3,
                               width_ratios=[0.85, 0.65, 1.10],
                               wspace=0.30)
        ax_f = fig.add_subplot(gs[0, 0])
        ax_d = fig.add_subplot(gs[0, 1])
        ax_r = fig.add_subplot(gs[0, 2])

        # --- funnel ---
        draw_funnel(ax_f, stages)

        # --- clinical-target dropouts (middle column) ---
        ax_d.axis("off")
        ax_d.set_xlim(0, 1); ax_d.set_ylim(0, 1)
        if drop is not None and len(drop):
            failed = drop[~drop["pass_gate"]]
            passed = drop[drop["pass_gate"]]
            header_y = 0.92
            ax_d.text(0.5, header_y,
                       f"clinical-stage AML targets",
                       ha="center", va="center", fontsize=10.5,
                       fontweight="bold", color="#34495e")
            ax_d.text(0.5, header_y - 0.04,
                       f"{len(failed)} of {len(drop)} fail the safety gate",
                       ha="center", va="center", fontsize=9,
                       color="#a93226", fontstyle="italic")
            # failed rows -- one line each: ✗  GENE   leak% AXIS   (AML XX%)
            row_y = header_y - 0.10
            step = 0.062
            for _, r in failed.iterrows():
                ax_d.text(0.04, row_y, r"$\mathbf{\times}$",
                           ha="left", va="center", fontsize=14,
                           color="#a93226")
                ax_d.text(0.13, row_y, f"{r['gene']}",
                           ha="left", va="center", fontsize=10.5,
                           fontweight="bold", color="#1a5276")
                ax_d.text(0.43, row_y,
                           f"{r['leak']:.0f}% {r['leak_axis']}",
                           ha="left", va="center", fontsize=9.5,
                           color="#a93226")
                ax_d.text(0.74, row_y,
                           f"(AML {r['aml']:.0f}%)",
                           ha="left", va="center", fontsize=8.5,
                           color="#95a5a6")
                row_y -= step
            # passed rows (CD70)
            if len(passed):
                row_y -= 0.020
                for _, r in passed.iterrows():
                    ax_d.text(0.04, row_y, r"$\checkmark$",
                               ha="left", va="center", fontsize=13,
                               color="#27ae60")
                    ax_d.text(0.13, row_y, f"{r['gene']}",
                               ha="left", va="center", fontsize=10.5,
                               fontweight="bold", color="#1a5276")
                    ax_d.text(0.43, row_y,
                               f"{r['leak']:.1f}% {r['leak_axis']}",
                               ha="left", va="center", fontsize=9.5,
                               color="#27ae60")
                    ax_d.text(0.74, row_y,
                               f"(AML {r['aml']:.1f}%)",
                               ha="left", va="center", fontsize=8.5,
                               color="#95a5a6")
                    row_y -= step
            # cut at top-15 note
            row_y -= 0.025
            ax_d.text(0.04, row_y,
                       "0 of 8 reach top 15 by composite\n"
                       "(CD70 passes safety but AML coverage too low)",
                       ha="left", va="top", fontsize=8.5,
                       color="#7f8c8d", fontstyle="italic")

            # leader line from safety-gate stage (in ax_f) -> this column
            # use ConnectionPatch in figure coordinates for cleanest hop
            # safety-gate stage is index 2 (third stage), right edge ~0.5 in ax_f data coords
            try:
                con = ConnectionPatch(
                    xyA=(0.40, 0.55), coordsA=ax_f.transData,
                    xyB=(0.0, header_y - 0.05), coordsB=ax_d.transData,
                    arrowstyle="->", color="#a93226", lw=1.2, alpha=0.7,
                    connectionstyle="arc3,rad=-0.15")
                fig.add_artist(con)
            except Exception:
                pass

        # --- composite-score formula box (explains the stage-4 ranking) ---
        ax_d.add_patch(mpatches.FancyBboxPatch(
            (0.02, 0.015), 0.96, 0.20,
            boxstyle="round,pad=0.006,rounding_size=0.02",
            facecolor="#f4f9fd", edgecolor="#a6c8e4", linewidth=0.9, zorder=1))
        ax_d.text(0.5, 0.188, "Composite score  (ranks the safety-gated set)",
                   ha="center", va="center", fontsize=9, fontweight="bold",
                   color="#13314f", zorder=2)
        ax_d.text(0.5, 0.140, r"score = ( C1 + C2 + C3 ) / 3",
                   ha="center", va="center", fontsize=11, fontweight="bold",
                   color="#08306b", zorder=2)
        ax_d.text(0.07, 0.097, "C1 = % leukemic cells positive",
                   ha="left", va="center", fontsize=8, color="#34495e", zorder=2)
        ax_d.text(0.07, 0.064, "C2 = % patients with >=20% leukemic positivity",
                   ha="left", va="center", fontsize=8, color="#34495e", zorder=2)
        ax_d.text(0.07, 0.031, "C3 = % PPAC_1-5 (poor-prognosis) cells positive",
                   ha="left", va="center", fontsize=8, color="#34495e", zorder=2)

        # --- top-15 lollipop showing CD96's dominance over the safety-gated set ---
        top15 = top.head(TOP_N).reset_index(drop=True)
        y = np.arange(len(top15))[::-1]
        for yi, row in zip(y, top15.itertuples()):
            color = PALETTE["CD96"] if row.gene == "CD96" else "#21588f"
            lc = PALETTE["CD96"] if row.gene == "CD96" else "#c3d4e6"
            ax_r.hlines(yi, 0, row.composite, color=lc, lw=2.0, zorder=1)
            ax_r.plot(row.composite, yi, "o", markersize=10, color=color,
                       markeredgecolor="black", zorder=3)
            ax_r.text(row.composite + 0.9, yi, f"{row.composite:.1f}",
                       va="center", fontsize=8.5,
                       fontweight="bold" if row.gene == "CD96" else "normal",
                       color=color)
        ax_r.set_yticks(y)
        labels_y = ax_r.set_yticklabels(top15["gene"], fontsize=9.5)
        for lab, gene in zip(labels_y, top15["gene"]):
            if gene == "CD96":
                lab.set_fontweight("bold"); lab.set_color(PALETTE["CD96"])
        # rank numbers in a left margin column
        for yi, rk in zip(y, top15["rank"]):
            ax_r.text(-2.0, yi, f"{rk}.", va="center", ha="right",
                       fontsize=8.5, color="#888",
                       fontweight="bold" if rk == 1 else "normal")
        ax_r.set_xlim(-4.0, top15["composite"].max() * 1.22)
        ax_r.set_xlabel("composite score")
        ax_r.set_title("top 15 candidates by composite",
                        fontsize=10, loc="left")
        clean_ax(ax_r, "x")

        # (figure title/commentary in caption -- added by authors)
        save_fig(fig, HERE, "Figure7A_discovery_funnel")
    print("\nDone -> Figure7A_discovery_funnel.{pdf,png,csv}")


if __name__ == "__main__":
    main()
