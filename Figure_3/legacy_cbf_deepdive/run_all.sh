#!/usr/bin/env bash
# Figure 3 (CBF deep-dive) — full pipeline. Re-runs on whatever data config.R/config.py
# point at (default: 08_Source_Data/Figure_3/v3). To swap datasets: export FIG3_DATA=<dir>.
# The dir must contain: LS_survival.csv (per-patient LS_* fractions + subtype2, from the
# CIBERSORTx stitch), LS_annotation.csv, meancpm.npy, genes.npy,
# reference_assignments.csv, jaccard_percluster.csv.
set -e
HERE="$(cd "$(dirname "$0")" && pwd)"
PY="/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/.venv/bin/python"
RS="Rscript"
echo "[1/8] within-CBF EFS Cox + focus state"; $RS "$HERE/01_cbf_prognosis.R"
echo "[2/8] panel A (CBF states)";            $RS "$HERE/02_panelA_cbf_states.R"
echo "[3/8] panel B (focus KM)";              $RS "$HERE/03_panelB_focus_km.R"
echo "[4/8] per-state Wilcoxon DE";           $PY "$HERE/04_compute_state_DE.py"
echo "[5/8] panel E (CBF GO:BP)";             $RS "$HERE/05_panelE_cbf_go.R"
echo "[6/8] build panel C/D matrices";        $PY "$HERE/06_build_panelCD.py"
echo "[7/8] panels C (targetable) + D (HOX)"; $RS "$HERE/07_panelsCD.R"
echo "[8/8] composite";                       $PY "$HERE/08_composite.py"
echo "DONE -> ../../04_Main_Figures/Figure_3__CBF_states.{png,pdf}"
