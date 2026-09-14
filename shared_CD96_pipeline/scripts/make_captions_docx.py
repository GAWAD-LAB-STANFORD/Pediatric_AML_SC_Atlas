#!/usr/bin/env python
# ====================================================================
# make_captions_docx.py  |  Upstream builder (regenerates source_data from raw inputs)
# reproducibility header (auto-generated; edit freely below this block)
# Inputs : (reads upstream object / raw matrix — see body)
# Outputs: (writes to ../figures or ../source_data — see body)
# Run    : python make_captions_docx.py   (or: bash run_all.sh)  |  Deps: see REPRODUCIBILITY_CHECKLIST.md
# Source data for this figure is in __SUBMISSION_PACKAGE/08_Source_Data/.
# ====================================================================

"""Generate methods/Figure_captions.docx (+ .txt) — one file, all figures, publication
prose, with the statistical test stated for every panel (n / test / correction / threshold).
"""
import os, re
from docx import Document
from docx.shared import Pt, RGBColor

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "methods")
NAVY = RGBColor(0x10, 0x3a, 0x66)

# (label, title, body, stats)
CAPS = [
("Figure 1", "Single-cell atlas of pediatric acute myeloid leukemia",
 "Single-cell RNA-seq of 28 pediatric AML and 2 healthy bone-marrow (BM) samples "
 "(96,627 cells after QC), integrated and annotated, with the NCI TARGET bulk-AML cohort "
 "used for downstream deconvolution and validation. (A) Study design: scRNA-seq of 28 AML "
 "+ 2 healthy BM and NCI TARGET bulk AML (n=1,046), bulk deconvolution to project "
 "single-cell clusters, and filtering to define Prognosis-Associated Clusters (PAC). "
 "(B) UMAP coloured by lineage. (C) UMAP coloured by cell type (29 annotated types). "
 "(D) UMAPs coloured by FAB class and by cytogenetic group. (E) Per-sample hierarchical-"
 "clustering dendrogram with the karyotype and short-variant presence/absence matrix and "
 "side annotations for FAB, survival, relapse and remission, plus blast %, MRD % and age.",
 "All panels descriptive (low-dimensional embeddings and Ward.D2 hierarchical clustering of "
 "the sample-by-principal-component correlation matrix); no inferential statistical test."),

("Figure 2", "Leukemic sub-clusters and Prognosis-Associated Clusters (PAC)",
 "(A) UMAP of CD34 expression and of the 27 leukemic sub-clusters (AML_1–AML_27). "
 "(B) UMAP of the Prognosis-Associated Clusters (0_HSPC, Myeloid_Pro, FPAC_1–2, PPAC_1–5). "
 "(C) PAC composition across samples and cytogenetic groups (dot size = fraction of cells, "
 "colour = cell count). (D) LSC (Eppert) and (E) HSC gene-set scores across PAC clusters "
 "(violin with box). (F) GO:BP terms enriched among genes up- and down-regulated in each PPAC.",
 "Differential genes by Wilcoxon rank-sum test (presto). GO:BP over-representation "
 "(clusterProfiler enricher) by one-sided hypergeometric test with Benjamini–Hochberg FDR "
 "(q<0.05; DE input thresholds padj<0.001 and |log2FC|>1). (C–E) descriptive "
 "(composition and score distributions)."),

("Figure 3", "Transcription-factor regulon landscape of leukemic cells",
 "Regulon activity inferred with SCENIC/pySCENIC (GRNBoost2 → cisTarget → AUCell). "
 "(A) Heatmap of Z-scored regulon activity across samples with prognosis/relapse annotation "
 "and dendrogram. (B) Regulon activity across cell types / PAC clusters. (C) Homeobox/HOX "
 "regulon module (e.g. HMX3, PITX1, SHOX2, HOXB6, LHX9, HOXA13/10/9/7/6, ZEB1) across PACs.",
 "Regulons and per-cell activity inferred by pySCENIC; values shown as Z-scored AUCell. "
 "Descriptive (no inferential test); ordering by hierarchical clustering."),

("Figure 4", "T- and NK-cell states and association with outcome",
 "(A) UMAP of marrow T/NK cells coloured by sample and (B) by cell type (naïve/activated "
 "CD4, CD8 subsets, GZMK/GZMB CD8 and NK, MAIT, DNT). (C) Exhaustion-marker dot plot "
 "(HAVCR2, PDCD1, LAG3) across T/NK subtypes. (D) TCR repertoire: fraction of cells per "
 "TRAV (top) and TRBV (bottom) V-gene call in each sample, with prognosis/relapse "
 "annotation. (F, G) Overall-survival Kaplan–Meier curves for TARGET patients stratified by "
 "deconvolved GZMB-CD8 (F) and GZMK-CD8 (G) abundance (high vs low).",
 "Overall survival by two-sided log-rank test (GZMK-CD8 high vs low p=0.0018; GZMB-CD8 "
 "p=0.0068; risk tables shown). (A–D) descriptive (embeddings, % positive / mean expression, "
 "and clonotype V-gene frequencies)."),

("Figure 5", "CD96 leads a panel of candidate surface immunotherapy targets in pediatric AML",
 "Twelve candidate surface antigens are compared throughout: four de novo leads (CD96, CD9, "
 "SUCNR1, IL1RAP), the pan-AML control FLT3, three clinical-stage antigens (CD33, CD123/IL3RA, "
 "CLL-1/CLEC12A) and four PPAC-specific markers (CD7, TNFRSF4, ABCA7, ITGAX). (A) In-silico "
 "surfaceome discovery funnel (2,129 surface genes → HSPC-sparing → myeloid-sparing → "
 "AML-expressed → top 15) and composite ranking (mean of % leukemic cells positive, % patients "
 ">20 % positive, and % PPAC cells positive); the four leads are highlighted and CD96 ranks "
 "first. (B) scRNA UMAP of the cell-type compartments and the log-normalized expression of the "
 "four lead targets. (C) Head-to-head of all twelve antigens: leukemic-cell coverage (y, % "
 "positive) vs normal HSPC/Myeloid_Pro expression (x, % positive), point size = % of patients "
 "with ≥20 % positive, each antigen in its group colour. (D) Target coverage (% positive, "
 "colour = mean expression) across cytogenetic subtypes (biological order; the three MPAL "
 "samples flagged in teal at right) and PAC clusters. (E) Pediatric NCI TARGET bulk expression "
 "(n=475, ≥80 % blasts) and (F) adult BeatAML bulk expression (n=104, ≥80 % blasts) for all "
 "twelve antigens as violins, with the % of patients above the targetable line annotated above "
 "each.",
 "Descriptive throughout: (A) composite ranking (no inferential test); (B–D) single-cell % "
 "positive and mean expression; (E, F) proportion of patients above a fixed bulk targetable "
 "threshold (log2 CPM scale), distributions shown as violins with median. Pediatric TARGET and "
 "adult BeatAML are independent bulk cohorts (concordant cross-cohort directionality)."),

("Figure 6", "Disqualifying candidate targets by normal-tissue toxicity",
 "Body-wide normal expression of the twelve candidate antigens, used to disqualify targets with "
 "unacceptable on-target/off-tumour liability. (A) Major-organ atlas (parenchymal cell types, "
 "CZ CELLxGENE Census): dot size = % of normal cells positive, colour = per-gene-scaled mean "
 "expression, antigens clustered by their own expression dendrogram. (B) Normal haematopoietic "
 "differentiation atlas (HSC → progenitors → mature myeloid and lymphoid arms). (C) Each "
 "antigen's top severe-toxicity cell types, grouped by organ (heart, brain, eye, lung, liver, "
 "pancreas, kidney, bladder, gut, skin, vasculature). (D) Radar profiles for all twelve "
 "antigens across five axes — AML coverage, patient breadth, marrow sparing, mature-immune "
 "sparing and vital-organ sparing — each axis min–max scaled across the panel (single-cell "
 "efficacy is dropout-limited, so the absolute axis maximum is printed); antigens are ranked by "
 "polygon area (overall-profile score). CD9, IL1RAP and ABCA7 are disqualified on vital-organ "
 "toxicity; CD123 is retained as a clinically validated antigen despite endothelial expression, "
 "leaving CD96, SUCNR1, CD33, CLL-1, CD7, TNFRSF4 and ITGAX to advance.",
 "Descriptive throughout (% positive and per-gene-scaled mean expression; the Census is a "
 "cross-platform relative reference and is not statistically compared to the AML cohort). Panel "
 "D axes are min–max scaled across the twelve antigens and the polygon-area score is geometric; "
 "no inferential statistical test."),

("Figure 7", "Choosing the best second target to extend CD96",
 "How to extend CD96 with a second target. (A) Target↔cell-type interaction circos for the "
 "three retained targets that are curated ligand–receptor partners (CD96, ITGAX, OX40/TNFRSF4) "
 "across eight marrow compartments; edge weight = co-expression interaction score (% "
 "ligand-positive sender × % receptor-positive receiver ÷ 100, receptor side = maximum over the "
 "axis ligands), the self-loop is the AML→AML autocrine signal, and all three panels share one "
 "weight scale. CD96 (the nectin checkpoint CD96→PVR/NECTIN1/NECTIN2[CD112], AML as sender) "
 "dominates with the strongest autocrine signal (5.9); ITGAX (ICAM1/FCER2A/C3→ITGAX+ITGB2) and "
 "OX40 (OX40L/TNFSF4→TNFRSF4) are weaker and AML-receiving. (B) CD96 + one-partner additive "
 "window: leukemic-cell coverage (y) of CD96 OR partner vs the combination's normal marrow "
 "toxicity (x = % of normal HSPC/Myeloid_Pro positive for CD96 OR partner), within patients "
 "(depth; mean across the 19 CD96-targetable patients) and between patients (breadth; % of "
 "patients with ≥20 % of blasts covered); CD96 alone is the floor and the green zone is "
 "marrow-sparing (<10 %); a red ring marks the antigens organ-toxic in Figure 6 (CD9, CD123, "
 "IL1RAP, ABCA7). (C) Two strategies to extend CD96, integrating the within- and "
 "across-patient coverage of (B) with toxicity in one decision view: within-patient depth (y) vs "
 "across-patient breadth (x) for each CD96 + partner, point fill = marrow toxicity, red ring = "
 "organ-toxic in Figure 6 (CD9, CD123, IL1RAP, ABCA7). Each line is one 2-target combination from "
 "CD96: two green lines mark the limit-toxicity options (CD96 + SUCNR1 for depth, CD96 + ITGAX for "
 "breadth — marrow-sparing and organ-clean) and two red lines the maximize-efficacy options "
 "(CD96 + CD33 or CD96 + CLL-1 — higher coverage, accepting marrow toxicity). Note ITGAX's high combination breadth here "
 "reflects complementarity with CD96 (it covers CD96-cold patients), not high single-agent "
 "breadth (Figure 6D). (D) Punch line — every 2-target combination, not anchored on CD96: all 66 "
 "pairwise combinations of the twelve antigens scored by a consolidated balanced coverage score "
 "(geometric mean of normalized within-patient depth and between-patient breadth) vs marrow "
 "toxicity. The winning combinations (Pareto frontier) are labelled with both partners; the "
 "low-toxicity (<10%) winners carry their depth/breadth/toxicity; a dark-red ring marks any combo "
 "with a toxic member — organ-toxic in Figure 6 (CD9, CD123, IL1RAP, ABCA7) or the marrow-toxic "
 "control FLT3 (so every FLT3 pair is ringed). (E) The same winning combinations quantified and "
 "SEPARATED into two toxicity tiers, each ranked by balanced coverage score: (top) low-toxicity "
 "options that are both organ-clean (no organ-toxic/FLT3 member) AND marrow-sparing (<10% "
 "HSPC/Myeloid_Pro positive) — the deployable combinations; (bottom) higher-coverage combinations "
 "that carry a toxic member (organ-toxic or FLT3) or are marrow-suppressive (≥10%) — the coverage "
 "ceiling shown for reference. Bar fill = marrow toxicity %, bar-tip = score · within-patient "
 "depth / between-patient breadth / marrow toxicity %, printed dark red when the pair carries an "
 "organ-toxic/FLT3 member and dark green when clean (so CD33+CLL-1 has a green tip yet sits in the "
 "higher-coverage tier because it is marrow-suppressive). Five of the six low-toxicity options are "
 "CD96 combinations (best CD96+ITGAX, balanced 71: depth 37%, breadth 82%, toxicity 7%), whereas "
 "five of the six higher-coverage combinations contain FLT3 — i.e. raising coverage past the "
 "low-toxicity tier requires accepting organ/FLT3 toxicity, and CD96 is the anchor for "
 "low-toxicity 2-target therapy. (Full four-view version in Figure S26.)",
 "Descriptive throughout. (A) co-expression interaction score (a descriptive single-cell "
 "connectivity metric — % ligand-positive × % receptor-positive — not a CellChat permutation "
 "probability, which is ≈0 for these sparse pairs); (B, C) % positive / union coverage with the "
 "per-axis marrow toxicity; the toxicity-free and efficacy frontiers are geometric. No inferential test."),
]

SUPP = [
("Figure S1", "Healthy bone-marrow reference",
 "Sample and cell-type embeddings of the 2 healthy BM controls, a canonical-marker dot plot, "
 "and force-directed embeddings of lineage markers (GATA1, PAX5, AZU1, CD1C, CD14, HES4) with "
 "pseudotime.", "Descriptive (embeddings, marker expression)."),
("Figure S2", "AML atlas covariates",
 "UMAPs coloured by sample, age, prognosis, relapse and remission; cell-type marker dot plot; "
 "pseudotime/entropy; and a marker-gene UMAP grid.", "Descriptive."),
("Figure S3", "Bulk deconvolution performance",
 "Per-sample observed vs expected cell-type fractions for AML pseudobulk deconvolution.",
 "Performance reported as per-sample Pearson correlation r and root-mean-square error (RMSE); "
 "no inferential test."),
("Figure S4", "Survival association of leukemic sub-clusters",
 "Hazard-ratio forest plots for the leukemic sub-clusters in two TARGET cohorts.",
 "Cox proportional-hazards regression; per-term hazard ratio with 95 % CI and two-sided Wald "
 "p value; global log-rank p reported per model."),
("Figure S5", "Transcription-factor regulon detail",
 "Per-cell-type regulon rank plots (SCENIC) and a regulon-activity heatmap.",
 "Descriptive (Z-scored AUCell)."),
("Figure S6", "hdWGCNA co-expression modules (formerly main Figure 4)",
 "Module force-directed graph and module dot plot (proliferation / myeloid-progenitor / HSC / "
 "myeloid-differentiation sets), a module fold-change volcano, GO:BP for representative modules "
 "(e.g. Module 9, Module 15), and hub-gene networks.",
 "Differential module-eigengene association (hdWGCNA) for the volcano (FDR-adjusted p and "
 "log2 fold-change); GO:BP over-representation (clusterProfiler, hypergeometric, BH-FDR) for the "
 "GO panels; module graphs descriptive."),
("Figure S7", "hdWGCNA module detail",
 "Per-module kME rank plots (Modules 1–18) and module FDR/log2 fold-change by cell type.",
 "Differential module-eigengene test (hdWGCNA; FDR-adjusted)."),
("Figure S8", "Cellular composition and cell cycle",
 "Per-sample cell-type composition and cell-cycle-phase bars; G1/G2M/S fractions by outcome "
 "group and by PAC class; and T-cell and B-lineage fractions, AML vs healthy BM.",
 "Between-group comparisons by two-sided Welch t-test (ggsignif; significance annotated; "
 "NS / *p<0.05 / **p<0.01 / ***p<0.001)."),
("Figure S9", "T/NK detail and survival",
 "Marker dot plot, CITE-seq surface-protein UMAPs, per-sample composition, and survival hazard "
 "ratios by deconvolved T/NK subset.",
 "Cox proportional-hazards (hazard ratio, 95 % CI, two-sided p) for the survival table; other "
 "panels descriptive."),
("Figure S10", "Surfaceome classification screen",
 "Cell-surface-feature performance vs AML classification, CITE-seq surface UMAPs, and candidate "
 "surface-gene expression maps.",
 "Surface markers ranked by Matthews correlation coefficient (MCC) and F1 classification "
 "metrics; no inferential test."),
("Figure S11", "CITE-seq normalization and gating",
 "Antibody-derived-tag normalization (none vs CLR vs dsb) and per-marker CD45 gating scatters.",
 "Descriptive (normalization comparison and gating thresholds)."),
# --- CD96 supplements in CALL-OUT order: S12-S21 support Figure 5, S22-S25 support Figure 7 ---
("Figure S12", "AML immunotherapy clinical-trial landscape",
 "Landscape of cell-surface antigens under clinical investigation in AML, positioning the "
 "discovered lead antigens against those already being pursued therapeutically (clinical-stage "
 "CAR-T and antibody targets).",
 "Descriptive (trial/target annotation); no inferential test."),
("Figure S13", "Discovered antigens versus antigens already in clinical trials",
 "Venn overlap between the composite-ranked discovery hits (Figure 5A) and the surface antigens "
 "currently in AML clinical trials, separating novel leads from already-pursued antigens.",
 "Descriptive (set overlap); no inferential test."),
("Figure S14", "PPAC-specific surface-target discovery",
 "Surface markers ranked within each individual poor-prognosis cluster (PPAC) by the percent of "
 "that cluster's cells that are positive (no composite). Candidates first pass a toxicity funnel "
 "with the tolerance raised to <=10% on both normal compartments (detectable -> HSPC<=10% AND "
 "Myeloid<=10% -> TCR/CD3 dropped; 1,706 eligible genes) — more permissive than the main Figure-5A "
 "funnel, by design. (A) all PPACs pooled; (B-F) PPAC_1-5. Each panel is a lollipop of the top "
 "markers by % positive; CD96 in crimson, ranked naturally (not forced first).",
 "Descriptive (single % metric; no inferential test). CD96 is the top marker in PPAC_1 (33%), "
 "PPAC_2 (60%) and PPAC_4 (35%) and ranks second in the pooled group (behind CD9). In PPAC_3 CD96 "
 "covers only ~8.5% (ABCA7, ITGAX, LAMP5 lead) — the coverage gap. PPAC_5 is a single-patient "
 "cluster (94.5% one patient), flagged not generalizable."),
("Figure S15", "Per-PPAC efficacy versus toxicity (Figure-5C construction), six panels",
 "The Figure-5C head-to-head — each surface marker as one point: x = % of normal HSPC/Myeloid_Pro "
 "positive (toxicity), y = % of that PPAC's cells positive (efficacy) — applied within each PPAC. "
 "(A) all PPACs pooled; (B-F) PPAC_1-5. Each panel plots the top-5 PPAC-specific markers (from "
 "Figure S14) plus the 12 panel-E antigens (the genes in Figure 7B/C), coloured by group (leads "
 "identity colours, FLT3 control gold, clinical orange, PPAC-specific purple); discovered markers "
 "are grey, CD96 crimson. Green band = sparing zone (leak <=10%).",
 "Descriptive (single-point % metrics; no inferential test). In PPAC_1, 2, 4 and the pooled group "
 "the high-efficacy markers inside the sparing zone are CD96 and CD9 (plus SUCNR1, CD7, TNFRSF4); "
 "clinical comparators (FLT3, CD33, CD123, CLL-1) are pushed to the toxic right. PPAC_3 is the "
 "coverage gap (ABCA7/ITGAX/LAMP5); PPAC_5 flagged not generalizable."),
("Figure S16", "Lead-target cell-type specificity",
 "Expression of the lead antigens across the normal marrow cell types — including the 0_HSPC and "
 "Myeloid_Pro safety compartments — i.e. the fuller cell-type safety view behind the Figure-5C "
 "head-to-head toxicity axis.",
 "Descriptive (% positive and mean expression per cell type); no inferential test."),
("Figure S17", "CD96-comparator co-expression in metacells and bulk",
 "CD96 versus each comparator in AML metacells (n=945; ~75-cell within-patient k-means pools, "
 "dropout-robust) and in TARGET bulk (n=1,330), coloured by cytogenetic subtype.",
 "Pearson r and Spearman ρ per target pair; metacell |r| ≤ 0.45 (highest FLT3, r=0.45, r²≈0.20; "
 "all others |r| ≤ 0.37) and bulk |r| ≤ 0.39 — i.e. complementary, not redundant."),
("Figure S18", "Therapeutic-window violins, KMT2Ar/MLL-excluded",
 "The twelve candidate antigens' TARGET bulk expression with KMT2A-rearranged/MLL cases excluded "
 "(n=340), as log2-CPM violins with the % of patients above the targetable line (log2 CPM=5) "
 "annotated above each — a robustness check that the targetability is not driven by KMT2Ar.",
 "Proportion of patients above a fixed threshold; descriptive (median shown)."),
("Figure S19", "Candidate antigens by cytogenetic subtype — pediatric TARGET",
 "The twelve candidate antigens' bulk expression (pediatric NCI TARGET, log2 CPM) faceted by "
 "cytogenetic subtype, exposing the subtype-specific coverage behind the pooled Figure-5E violins.",
 "Descriptive (proportion above the targetable line per subtype); no inferential test."),
("Figure S20", "Candidate antigens by cytogenetic subtype — adult BeatAML",
 "As Figure S19 for the adult BeatAML bulk cohort (log2 CPM), faceted by cytogenetic / ELN "
 "subtype, for cross-cohort comparison with the pediatric panel (Figure 5F / S19).",
 "Descriptive; no inferential test."),
("Figure S21", "Adult AML single-cell target validation",
 "Whether the twelve candidate antigens are expressed in ADULT AML blasts at single-cell level, "
 "at what % of cells, and across which molecular subtypes — processed through the pediatric "
 "Figure-1 pipeline on the 748k-cell AML scAtlas (adult, malignant blasts only). (A) UMAP, "
 "(B) overall positivity, (C) positivity by ELN molecular-genetic subtype.",
 "Descriptive (% positive / mean expression; restricted to malignant blasts)."),
("Figure S22", "CD96 cell-cell communication (CellChat)",
 "Full sender→receiver CD96 network (including T-cell→AML and NK→AML), the autocrine substrate "
 "(CD96 / PVR / NECTIN1 co-expression on leukemic cells), per-patient heterogeneity, and "
 "signalling strength vs % CD96+ blasts (targetable vs CD96-cold AML vs healthy BM).",
 "Per-sample CellChat one-sided permutation test (p≤0.05); reproducible non-autocrine edges: "
 "T-cell→AML 9/19 and NK→AML 5/14 targetable patients. The blast-positivity panel is "
 "descriptive (medians: 7.5×10⁻³ targetable vs 1.1×10⁻⁵ cold AML vs 4.0×10⁻⁶ healthy)."),
("Figure S23", "Additive therapeutic window per PAC",
 "The Figure-7C additive window applied per PAC: each point is CD96 alone or CD96 OR one partner; "
 "x = normal HSPC/Myeloid_Pro toxicity, y = % of PAC cells covered; points sorted by toxicity with "
 "the Pareto frontier (dotted), green = sparing / red = toxic. (A) all PPACs pooled (n=13,746); "
 "(B) PPACs 1, 2, 5 pooled (n=9,388); (C) the seven leukemic PAC groups, with a target colour key.",
 "Descriptive (union % of cells positive). Note: efficacy here is over all patients' leukemic "
 "cells, whereas the main additive window (Figure 7C) restricts to the 19 CD96-targetable "
 "patients, so CD96-alone coverage is lower here."),
("Figure S24", "Additive therapeutic window by cytogenetic subtype",
 "As Figure S23 but grouped by cytogenetics: efficacy = % of that subtype's leukemic (core-AML) "
 "cells covered by CD96 OR one partner. (A) all AML pooled (n=68,711 core-AML cells); (B) the 14 "
 "cytogenetic subtypes with ≥500 core-AML cells, ordered by CD96 coverage, with a colour key.",
 "Descriptive (union % positive). CD96 covers the most cells in BCR/ABL, RUNX1/RUNX1T1, PML/RARA, "
 "t(2;3) and t(7;14), and is CD96-cold in MLLr, del7q, MYB/GATA1, NUP98/NSD1 and Tri(8)/MLLr. "
 "Several subtypes are single-patient and should be read with caution. As in Figure S23, efficacy "
 "is over all patients' leukemic cells (not identical to the main window, Figure 7C)."),
("Figure S25", "Additive CD96 combinations per PPAC (Figure-7C construction), six panels",
 "The Figure-7C additive therapeutic window applied within each PPAC: CD96 alone is the floor and "
 "each point adds one partner as the union 'CD96 OR partner' (x = max % normal HSPC/Myeloid_Pro "
 "covered; y = % of that PPAC's cells covered). Pareto frontier dotted; green = sparing (leak "
 "<=10%). (A) all PPACs pooled; (B-F) PPAC_1-5. Partners = the 12 panel-E antigens (the genes in "
 "Figure 7B/C), coloured by group (leads identity colours, FLT3 control gold, clinical orange, "
 "PPAC-specific purple).",
 "Descriptive (union % covered; no inferential test). Within each PPAC the best marrow-sparing "
 "CD96 combinations use the low-toxicity partners (e.g. +CD9, +ITGAX, +SUCNR1), while the clinical "
 "antigens (+CD33, +CLL-1) and FLT3 add more coverage only by leaving the sparing zone. PPAC_5 is "
 "a single-patient cluster, flagged not generalizable."),
("Figure S26", "Every two-target combination landscape (generalizes Figure 7B)",
 "All 66 pairwise combinations of the twelve panel-E antigens, scored as the union 'X OR Y' per "
 "cell, with marrow toxicity on x (% of HSPC + Myeloid_Pro positive). (A) within patients (depth; "
 "% of leukemic blasts covered, mean across all 28 AML patients); (B) between patients (breadth; "
 "% of patients with ≥20% of blasts covered); (C) the consolidated balanced coverage score = "
 "100 × √(depth_norm × breadth_norm), each axis normalised to its maximum over the pairs, so a "
 "combination scores high only when strong on BOTH depth and breadth — the single best-combo "
 "metric. (D) The winning combinations quantified as ranked bars, SEPARATED into two toxicity "
 "tiers: (top) low-toxicity options — organ-clean AND marrow-sparing (<10%); (bottom) "
 "higher-coverage combinations carrying a toxic member (organ-toxic/FLT3) or marrow-suppressive "
 "(≥10%). Bar fill = marrow toxicity; bar-tip = score · depth / breadth / toxicity %, dark red "
 "when the pair carries an organ-toxic/FLT3 member and dark green when clean. In the scatter "
 "panels the winners (the Pareto frontier = best efficacy reachable at each toxicity) are filled "
 "by marrow toxicity and labelled with both partners; low-toxicity (<10%) winners also carry their "
 "(depth / breadth / toxicity %). Other pairs are grey; single targets are open diamonds (the "
 "floor); the green zone is marrow-sparing (<10%); a dark-red ring (panels A–C) marks winners with "
 "a toxic member — organ-toxic in Figure 6 (CD9, CD123, IL1RAP, ABCA7) or the marrow-toxic control "
 "FLT3 (so every FLT3 pair is flagged).",
 "Descriptive (union % positive / patient coverage; the balanced score is geometric; no "
 "inferential test). In the marrow-sparing zone the best combinations all contain CD96; the best "
 "fully-clean (marrow- AND organ-sparing) 2-target combination is CD96+ITGAX (balanced 71; depth "
 "37%, breadth 82%, toxicity 7%), then CD96+SUCNR1 / CD7 / TNFRSF4. Raising the toxicity ceiling, "
 "the FLT3- and CD33/CLL-1-containing pairs win on raw coverage (FLT3+CD33 balanced 98) but leave "
 "the sparing zone — i.e. CD96 is the anchor for low-toxicity 2-target therapy."),
("Table S1", "Patient and sample characteristics",
 "Per-sample clinical and molecular table (28 AML + 2 healthy BM): diagnosis, FAB, cytogenetics, "
 "short variants, MRD %, blast %, age, disease status, sex, CNS status, relapse and survival.",
 "Descriptive."),
]

PANEL = re.compile(r"(\([A-G](?:\s*[,–-]\s*[A-G])*\))")

def add_caption(doc, label, title, body, stats):
    p = doc.add_paragraph(); p.paragraph_format.space_after = Pt(10)
    r = p.add_run(f"{label}. "); r.bold = True; r.font.color.rgb = NAVY
    r = p.add_run(f"{title}. "); r.bold = True
    for seg in PANEL.split(body):
        run = p.add_run(seg)
        if PANEL.fullmatch(seg):
            run.bold = True
    r = p.add_run("  Statistical analysis: "); r.italic = True; r.bold = True
    p.add_run(stats).italic = True

doc = Document()
doc.styles["Normal"].font.name = "Calibri"; doc.styles["Normal"].font.size = Pt(10.5)
h = doc.add_paragraph(); rr = h.add_run("Figure legends"); rr.bold = True; rr.font.size = Pt(15)
sub = doc.add_paragraph(
    "Pediatric AML scRNA-seq / CD96 manuscript (2026). Single-cell expression is "
    "log-normalized (log1p of counts-per-10,000); bulk expression is log2(CPM+1). Dot plots: "
    "dot size = percent of cells positive (≥1 count), colour = mean expression. The statistical "
    "test is stated at the end of each legend.")
sub.runs[0].italic = True
for c in CAPS: add_caption(doc, *c)
hh = doc.add_paragraph(); r = hh.add_run("Supplementary figures and tables"); r.bold = True; r.font.size = Pt(13)
for c in SUPP: add_caption(doc, *c)

docx_path = os.path.join(OUT, "Figure_captions.docx"); doc.save(docx_path)
# plain-text mirror
with open(os.path.join(OUT, "Figure_captions.txt"), "w") as fh:
    fh.write("FIGURE LEGENDS — pediatric AML scRNA-seq / CD96 manuscript (2026)\n\n")
    for label, title, body, stats in CAPS + SUPP:
        fh.write(f"{label}. {title}. {body}  Statistical analysis: {stats}\n\n")
print("wrote", docx_path, "and Figure_captions.txt")
