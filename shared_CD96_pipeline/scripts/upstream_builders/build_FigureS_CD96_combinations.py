#!/usr/bin/env python
# ====================================================================
# build_FigureS_CD96_combinations.py  |  Upstream builder (regenerates source_data from raw inputs)
# reproducibility header (auto-generated; edit freely below this block)
# Inputs : (reads upstream object / raw matrix — see body)
# Outputs: (writes to ../figures or ../source_data — see body)
# Run    : python build_FigureS_CD96_combinations.py   (or: bash run_all.sh)  |  Deps: see REPRODUCIBILITY_CHECKLIST.md
# Source data for this figure is in __SUBMISSION_PACKAGE/08_Source_Data/.
# ====================================================================

"""Supplement: best 2nd target to PAIR with CD96 — two complementary questions,
both toxicity-aware.

(1) WITHIN PATIENT (single cell): of the leukemic cells a CD96 therapy would
    MISS (CD96 raw==0), what fraction does each candidate target capture
    (X raw>0)?  Computed per AML patient, then averaged -> "% of CD96-negative
    blasts rescued".  Answers: which partner mops up escapees inside each tumor.

(2) ACROSS PATIENTS (bulk): a patient is targetable by a gene if bulk log2 CPM
    > TARGET_THR (=5, same line as the violins; TARGET >=80% blasts).  For each
    partner X, % of patients covered by CD96 OR X (union), and the increment
    over CD96 alone.  Answers: which pair covers the most patients.  A greedy
    panel (CD96 first) and a SAFE-only greedy (HSPC/Myeloid-sparing partners)
    are also reported.

TOXICITY (x-axis): % of normal HSPC / Myeloid_Pro cells hit by the COMBINATION
you would actually give -- CD96 OR partner -- taken as the max over the two
safety clusters.  Because it is a union, CD96 alone (CD96 OR CD96) is the floor:
adding any partner can only move a point to the right, never left.  sparing =
combination <= TOX_THR (5%) -> green, else grey.  CD96 = crimson.
Pure data (h5ad + TARGET STAR counts); no fabricated values.
"""
import os
import sys
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
import matplotlib.pyplot as plt
import matplotlib as mpl

HERE = os.path.dirname(os.path.abspath(__file__))
SUPP = os.path.join(HERE, "supplements")
sys.path.insert(0, SUPP)
from _fig7_style import PALETTE, RC, clean_ax, save_fig

H5AD = os.path.join(HERE, "..", "ALSF_AML_plot", "H5AD",
                    "ALSF_AML_Combo_3500_with_PAC_new.h5ad")
FILT = os.path.join(HERE, "..", "ALSF_AML_plot", "H5AD",
                    "NCBI_TARGET_Star_counts_filtered.csv")
MERGED = os.path.join(HERE, "Figure7F_TARGET_CD96_merged.csv")
PAC = os.path.join(HERE, "Figure7D_target_by_PAC.csv")
BLAST = "Bone marrow leukemic blast percentage (%)"
BLAST_MIN = 80.0
TARGET_THR = 5.0
TOX_THR = 5.0
MIN_NEG = 30          # min CD96-negative leukemic cells/patient to score rescue
CD96POS_MIN = 10.0    # within-patient analysis uses only CD96-positive patients
                      # (>= this % of blasts CD96+); CD96-cold patients aren't
                      # CD96 candidates and would drag the within-tumor % down

ENS = {"CD96": "ENSG00000153283", "CD33": "ENSG00000105383",
       "IL3RA": "ENSG00000185291", "CLEC12A": "ENSG00000172322",
       "FLT3": "ENSG00000122025", "MSLN": "ENSG00000102854",
       "CD70": "ENSG00000125726", "ADGRG1": "ENSG00000205336",
       "IL1RAP": "ENSG00000196083"}
GENES = list(ENS)
OTHERS = [g for g in GENES if g != "CD96"]
GLABEL = {"CD33": "CD33", "IL3RA": "CD123", "CLEC12A": "CLL-1", "FLT3": "FLT3",
          "MSLN": "MSLN", "CD70": "CD70", "ADGRG1": "GPR56", "IL1RAP": "IL1RAP"}
LEUK_SUBS = [f"AML_{i}" for i in range(1, 28)]


def _dec(a):
    return [x.decode() if isinstance(x, bytes) else x for x in a]


def _cat(f, k):
    return np.array(pd.Categorical.from_codes(
        f["obs"][k]["codes"][:].astype(int),
        categories=_dec(f["obs"][k]["categories"][:])).astype(object))


SAFETY_PACS = ["0_HSPC", "Myeloid_Pro"]   # normal stem/progenitor safety clusters


def within_patient():
    with h5py.File(H5AD, "r") as f:
        sub = _cat(f, "AML Sub-Clusters")
        samp = _cat(f, "Sample")
        pac = _cat(f, "Prognosis-Associated Clusters")
        rx = f["raw"]["X"]
        X = sp.csr_matrix((rx["data"][:], rx["indices"][:], rx["indptr"][:]),
                          shape=tuple(rx.attrs["shape"])).tocsc()
        rk = f["raw"]["var"].attrs.get("_index", b"_index")
        rk = rk.decode() if isinstance(rk, bytes) else rk
        symbols = _dec(f["raw"]["var"][rk][:])
        pos = {g: np.asarray(X[:, symbols.index(g)].todense(),
                             np.float32).ravel() > 0 for g in GENES}
    leuk = np.isin(sub, LEUK_SUBS)

    # toxicity = % of normal HSPC/Myeloid_Pro cells hit by CD96 alone, and by the
    # COMBINATION (CD96 OR partner); max over the two safety clusters (house
    # metric).  The union guarantees CD96 alone is the floor.
    cd96_tox = max(100.0 * pos["CD96"][pac == c].mean() for c in SAFETY_PACS)
    union_tox = {g: max(100.0 * (pos["CD96"][pac == c] | pos[g][pac == c]).mean()
                        for c in SAFETY_PACS) for g in OTHERS}
    res = {g: [] for g in OTHERS}        # rescue of CD96-neg cells (XY size)
    comb = {g: [] for g in OTHERS}       # per-patient % blasts CD96+ OR X+
    cd96cov = []
    n_excl = 0
    for s in sorted(set(samp[leuk])):
        if str(s).startswith("0_Healthy"):
            continue
        m = leuk & (samp == s)
        c96 = pos["CD96"][m]
        c96frac = 100.0 * c96.mean()
        if c96frac < CD96POS_MIN:        # CD96-negative tumor: not a CD96 candidate
            n_excl += 1
            continue
        cd96cov.append(c96frac)
        for g in OTHERS:
            comb[g].append(100.0 * (c96 | pos[g][m]).mean())
        neg = m & (~pos["CD96"])
        nneg = int(neg.sum())
        if nneg >= MIN_NEG:
            for g in OTHERS:
                res[g].append(100.0 * (neg & pos[g]).sum() / nneg)
    cd96_base = float(np.mean(cd96cov))
    rescue = pd.Series({g: float(np.mean(res[g])) for g in OTHERS})
    combined = pd.Series({g: float(np.mean(comb[g])) for g in OTHERS})
    npat = len(cd96cov)
    print(f"  within-patient: {npat} CD96-positive patients (≥{CD96POS_MIN:.0f}% "
          f"blasts CD96+); excluded {n_excl} CD96-negative; "
          f"CD96 baseline={cd96_base:.0f}%")
    return cd96_base, rescue, combined, comb, npat, cd96_tox, union_tox


def across_patient():
    raw = pd.read_csv(FILT).set_index("ensembl_gene_id")
    raw.index = raw.index.astype(str).str.split(".").str[0]
    g9 = raw.loc[[ENS[g] for g in GENES]].T
    g9.columns = GENES
    m = pd.read_csv(MERGED).set_index("sample")
    ids = [s for s in g9.index if s in m.index
           and pd.notna(m.loc[s, BLAST]) and m.loc[s, BLAST] >= BLAST_MIN]
    lib = m.loc[ids, "lib_size"].astype(float)
    L = np.log2(g9.loc[ids].div(lib, axis=0) * 1e6 + 1)
    tgt = (L > TARGET_THR)
    base = float(tgt["CD96"].mean() * 100)
    union = pd.Series({g: float((tgt["CD96"] | tgt[g]).mean() * 100) for g in OTHERS})
    return base, union, tgt, len(ids)


def greedy(tgt, partners):
    covered = tgt["CD96"].to_numpy().copy()
    order, cum = ["CD96"], [100.0 * covered.mean()]
    rem = list(partners)
    while rem:
        best = max(rem, key=lambda g: (tgt[g].to_numpy() & ~covered).sum())
        covered = covered | tgt[best].to_numpy()
        order.append(best); cum.append(100.0 * covered.mean()); rem.remove(best)
    return order, cum


def main():
    cd96_base, rescue, combined, comb_per, npat, cd96_tox, union_tox = within_patient()
    tox = pd.Series(union_tox)        # x-axis toxicity = combination CD96 OR partner
    sparing = {g: float(tox[g]) <= TOX_THR for g in OTHERS}

    base, union, tgt, nbulk = across_patient()
    print(f"  across-patient: n={nbulk}; CD96 alone covers {base:.1f}%; "
          f"within-patient CD96 baseline {cd96_base:.1f}%")

    summ = pd.DataFrame({
        "within_CD96orX_pct": combined.round(2),
        "within_rescue_pct": rescue.round(2),
        "union_cov_pct": union.round(2),
        "increment_over_CD96_pct": (union - base).round(2),
        "HSPC_Myeloid_pct": tox.reindex(OTHERS).round(2),
        "sparing": pd.Series(sparing)})
    summ = summ.sort_values("increment_over_CD96_pct", ascending=False)
    summ.loc["CD96"] = [cd96_base, np.nan, base, 0.0, cd96_tox, True]  # CD96 ref
    summ.to_csv(os.path.join(HERE, "FigureS_CD96_combinations.csv"))
    print("\n", summ.to_string())

    order_all, cum_all = greedy(tgt, OTHERS)
    safe = [g for g in OTHERS if sparing[g]]
    order_safe, cum_safe = greedy(tgt, safe)
    print("\n  greedy (all partners):", list(zip(order_all, [round(c, 1) for c in cum_all])))
    print("  greedy (SAFE only):   ", list(zip(order_safe, [round(c, 1) for c in cum_safe])))

    def col(g):
        return "#1e8449" if sparing[g] else "#95a5a6"

    with mpl.rc_context(RC):
        fig, (axA, axB) = plt.subplots(1, 2, figsize=(13.0, 5.6),
                                       gridspec_kw=dict(wspace=0.42))
        # ---- A: within-patient blasts covered, CD96 baseline + partner ----
        oa = combined.sort_values().index.tolist()
        y = np.arange(len(oa))
        axA.barh(y, cd96_base, color=PALETTE["CD96"], edgecolor="white",
                 height=0.7, zorder=3)                        # CD96 baseline
        axA.barh(y, (combined[oa] - cd96_base).values, left=cd96_base,
                 color=[col(g) for g in oa], edgecolor="white", height=0.7,
                 zorder=3)                                    # partner adds
        for g, yi in zip(oa, y):                              # per-patient points
            v = np.asarray(comb_per[g])
            axA.scatter(v, yi + np.linspace(-0.16, 0.16, len(v)), s=6,
                        color="#34495e", alpha=0.22, zorder=4, linewidths=0)
            axA.text(combined[g] + 1, yi, f"{combined[g]:.0f}%", va="center",
                     fontsize=8)
        axA.set_yticks(y)
        axA.set_yticklabels([f"+ {GLABEL[g]}  ({tox[g]:.0f}%)" for g in oa],
                            fontsize=8.5)
        for lab, g in zip(axA.get_yticklabels(), oa):
            lab.set_color(col(g))
        axA.set_xlabel("% leukemic cells covered by CD96 OR partner\n"
                       f"(mean across {npat} CD96-positive patients; single cell)")
        axA.set_title("within patients", fontsize=10, loc="left", fontweight="bold")
        clean_ax(axA, "x")

        # ---- B: across-patient coverage, CD96 baseline + partner ----
        ob = union.sort_values().index.tolist()
        y2 = np.arange(len(ob))
        axB.barh(y2, base, color=PALETTE["CD96"], edgecolor="white",
                 height=0.7, zorder=3)                        # CD96 baseline
        axB.barh(y2, (union[ob] - base).values, left=base,
                 color=[col(g) for g in ob], edgecolor="white", height=0.7,
                 zorder=3)                                    # partner adds
        for g, yi in zip(ob, y2):
            axB.text(union[g] + 0.6, yi, f"{union[g]:.0f}%  (+{union[g]-base:.0f})",
                     va="center", fontsize=8)
        axB.set_yticks(y2)
        axB.set_yticklabels([f"+ {GLABEL[g]}  ({tox[g]:.0f}%)" for g in ob],
                            fontsize=8.5)
        for lab, g in zip(axB.get_yticklabels(), ob):
            lab.set_color(col(g))
        axB.set_xlim(0, 100)
        axB.set_xlabel(f"% patients targetable by CD96 OR partner\n"
                       f"(TARGET ≥{int(BLAST_MIN)}% blasts, n={nbulk}; log2CPM>{int(TARGET_THR)})")
        axB.set_title("across patients", fontsize=10, loc="left", fontweight="bold")
        clean_ax(axB, "x")

        from matplotlib.patches import Patch
        fig.legend(handles=[Patch(color=PALETTE["CD96"], label="CD96 baseline"),
                            Patch(color="#1e8449", label=f"sparing partner (≤{int(TOX_THR)}%)"),
                            Patch(color="#95a5a6", label="toxic partner (>5%)")],
                   loc="lower center", ncol=3, frameon=False, fontsize=8.5,
                   bbox_to_anchor=(0.5, -0.02))
        # (figure title/commentary in caption -- added by authors)
        fig.subplots_adjust(bottom=0.20)
        save_fig(fig, HERE, "FigureS_CD96_combinations")

    # ---- MAIN-figure XY: coverage (Y) vs HSPC/Myeloid toxicity (X), drawn
    #      BOTH within patients (cells covered) and across patients (patients
    #      covered).  Same axes structure; CD96-alone = crimson star baseline. ----
    xmax = max(float(tox[OTHERS].max()), 10) * 1.12

    def draw_window(ax, yser, ybase, ylabel, title, ylim):
        ax.axvspan(0, TOX_THR, color="#eafaf1", zorder=0)
        ax.axvspan(TOX_THR, xmax, color="#fdecea", zorder=0)
        ax.axvline(TOX_THR, color="#1e8449", ls="--", lw=1.0, zorder=1)
        pts = sorted([(cd96_tox, ybase)] + [(float(tox[g]), float(yser[g]))
                                            for g in OTHERS])
        front, ym = [], -1.0                              # Pareto frontier
        for x, y in pts:
            if y > ym + 1e-9:
                front.append((x, y)); ym = y
        ax.plot([p[0] for p in front], [p[1] for p in front], color="#7f8c8d",
                ls=":", lw=1.2, zorder=2)
        ax.scatter([cd96_tox], [ybase], marker="o", s=150, c=PALETTE["CD96"],
                   edgecolor="black", linewidth=1.2, zorder=6)
        ax.annotate("CD96 alone", (cd96_tox, ybase), textcoords="offset points",
                    xytext=(11, -3), fontsize=8.5, fontweight="bold",
                    color=PALETTE["CD96"])
        for g in OTHERS:
            c = "#1e8449" if sparing[g] else "#95a5a6"
            ax.scatter([tox[g]], [yser[g]], s=110, c=c, edgecolor="black",
                       linewidth=0.8, alpha=0.92, zorder=5)
            ax.annotate(f"+{GLABEL[g]}", (tox[g], yser[g]),
                        textcoords="offset points", xytext=(7, 5), fontsize=8.2,
                        color=c, fontweight="bold")
        ax.set_xlim(-1.5, xmax); ax.set_ylim(*ylim)
        ax.set_xlabel("% normal HSPC / Myeloid_Pro positive — CD96 OR partner  →  toxicity")
        ax.set_ylabel(ylabel)
        ax.set_title(title, fontsize=10, loc="left", fontweight="bold")
        clean_ax(ax, "both")

    with mpl.rc_context(RC):
        fig2, (axW, axA) = plt.subplots(1, 2, figsize=(13.6, 6.2),
                                        gridspec_kw=dict(wspace=0.30))
        draw_window(axW, combined, cd96_base,
                    "% leukemic cells covered by CD96 + partner\n"
                    f"(mean across {npat} CD96-positive patients; single cell)",
                    "within patients", (20, 80))
        draw_window(axA, union, base,
                    "% patients targetable by CD96 + partner\n"
                    f"(TARGET ≥{int(BLAST_MIN)}% blasts)",
                    "across patients", (45, 100))
        # (figure title/commentary in caption -- added by authors)
        save_fig(fig2, HERE, "Figure_CD96_combination_window")
    print("\nDone -> FigureS_CD96_combinations.{pdf,png,csv} + "
          "Figure_CD96_combination_window.{pdf,png}")


if __name__ == "__main__":
    main()
