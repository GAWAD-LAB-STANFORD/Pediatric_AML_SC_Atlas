#!/usr/bin/env bash
# Master runner for the AML_2026 (alternate) figure package.
#   bash run_all.sh
# Renders the reproducible CD96 figures 5-7 + CD96 supplements from bundled source_data/,
# assembles the multi-panel figures (fitz), then rebuilds the combined PDF. Legacy figures 1-4
# (+ their supplements S1-S11, Table S1) are frozen pages pulled from the manuscript PDF by
# build_AllFigures_pdf.py. Override interpreters with RSCRIPT=... PY=...
set -uo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
RSCRIPT="${RSCRIPT:-Rscript}"
PY="${PY:-python3}"
cd "$ROOT"
mkdir -p figures/supplementary methods figures/legacy
cd "$ROOT/scripts"

# --- (optional) rebuild source_data from raw inputs (count matrix + Census + cBioPortal) ----
# Only needed if the underlying data changed; the bundled source_data/ already contains these.
#   $PY build_fig5E_target_leads.py ; $PY build_adult_bulk.py          # Fig 5E/F bulk cohorts
#   $PY build_fig6_toxicity.py ; $PY build_fig6_DE.py                  # Fig 6 (Census + matrix)
#   $PY build_fig7_combo.py ; $PY build_fig7_allpairs.py               # Fig 7B/C + S26 all-pairs
#   $PY cellchat/build_cellchat_3path_input.py                        # Fig 7A circos input
#   $PY build_adult_scrna_markers.py ; $PY build_adult_scrna_umap.py  # Fig S25
#   $PY build_funnel_ls_specific.py                                   # S15 poor-state surface targets (49-state)
#   $PY build_ls_headtohead.py                                        # S16 per-state efficacy-vs-toxicity (needs S15 output)
#   $PY build_metacells_S7.py ; $PY build_bulk_S17.py                 # S17 CD96-comparator co-expression (metacells + TARGET bulk, 12-antigen panel)
#   $PY build_cd96_axis.py                                            # S24 E/F CD96-nectin checkpoint axis (atlas + T/NK object)
#   $PY build_ls_window.py                                            # S25 additive window per poor-LS (merged former S25+S27; 49-state)
#   $PY build_cyto_window.py                                          # S26 additive window by cytogenetic subtype (atlas; 4-type core-AML)
# (Full list + external-data provenance in REPRODUCIBILITY_CHECKLIST.md / SOURCE_DATA.md.)

run() { echo "   -> $1"; mkdir -p "$ROOT/logs"; "$RSCRIPT" "$1" >"$ROOT/logs/$(basename "$1").log" 2>&1 && echo "      ok" || echo "      FAILED: $1  (see logs/$(basename "$1").log)"; }

echo "== [1/4] Figure 5  (target discovery) =="
run fig5.R

echo "== [2/4] Figure 6  (toxicity disqualification: organ + haematopoiesis + deadly + radar) =="
run fig6_toxicity.R          # panels A/B/C  -> .panel_fig6ABC.pdf
run fig6_DE.R                # panel D radar -> .panel_fig6D.pdf
$PY build_fig6_assemble.py   # stack -> Figure6.pdf + Figure6.png

echo "== [3/4] Figure 7  (combination strategy: circos + additive window + best-combo) =="
run cellchat/fig7_circos_3panel.R   # panel A -> fig7A_circos.pdf
run fig7_combo.R                    # panel B -> Figure7_CD96_additive_window.pdf
run fig7d_bestcombo.R               # panel C -> Figure7D_best_combo.pdf
run fig7_allpairs.R                 # panel D punch line -> Figure7D_punchline.pdf (+ Fig S26)
$PY build_fig7_assemble.py          # stack -> Figure7.pdf + Figure7.png
# fig7_coverage.R (cumulative coverage curve) is retained as an analysis asset but no longer a Fig-7 panel

echo "== [4/4] CD96 supplements =="
# figS_ls_specific.R = deployed S15 (poor-prognosis-STATE surface targets; replaces the retired
#   per-PPAC figS_ppac_specific.R). Its data comes from build_funnel_ls_specific.py (see optional block).
for s in figS.R figS8.R figS_ls_window.R figS_cyto_window.R figS_ls_specific.R \
         figS_ls_headtohead.R figS_leads_celltype.R figS_aml_trials.R \
         figS_venn.R figS_cyto_violin.R fig_adult_scrna.R; do run "$s"; done   # (S28 built in step 3)

echo "== captions =="
$PY make_captions_docx.py    # -> methods/Figure_captions.docx + .txt

echo "== combined figure PDF =="
cd "$ROOT"; [ -f build_AllFigures_pdf.py ] && $PY build_AllFigures_pdf.py || echo "   (skip: build_AllFigures_pdf.py not bundled)"

echo "DONE -> figures/ updated; AML_2026_AllFigures.pdf rebuilt."
