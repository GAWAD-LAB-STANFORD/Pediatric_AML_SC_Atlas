# What's needed for total reproducibility

There are **two levels**. The package already meets level 1 for the CD96 figures and about
half the legacy figures; level 2 (raw reads → everything) has a real gap.

- **Level 1 — Figure reproducibility:** regenerate every figure + derived table from the
  processed count matrix + source-data tables.
- **Level 2 — End-to-end reproducibility:** raw FASTQ → counts → QC/integration/clustering/
  annotation/PAC → all downstream. Requires the upstream pipeline + environments + references.

Status today: ✅ CD96 (5–7, S12–S14) and Figs 1,2,S1,S2,S6,S7,S8 reproduce now;
the rest is gated on the items below.

---

## 1. Input files still to obtain  (finite; see MISSING_AND_EXTRAS.md)
- [ ] **`AML-sample-infor.csv`** (clinical = Table S1) — Box/Yakun. *Unblocks Fig 1E, Fig 4, Table S1.*
- [ ] **SCENIC** `reg-py_3.8.csv` (+ adjacency, TF list) — Oak (`retrieve_from_oak.sh`). *Fig 3, S5.*
- [ ] **`NCI_TARGET_gene_counts_filtered_anno.csv`** — Oak. *Surfaceome/S10.*
- [ ] **CIBERSORTx_T_NK `*Results.csv`** — Yakun or rerun CIBERSORTx. *Fig 4 survival, S3, S9.*
- [ ] **TARGET clinical** (3 `.xlsx`) — public (GDC). *Fig 4/S4/S9 survival.*
- [ ] **CITE-seq ADT object** (`…CITEseq…h5mu`) — Oak. *S9, S11.*
- [ ] surfaceome **expr matrices** (`ALSF_RNA-seq_expr_*`) — or regenerate from counts. *S10.*
- [ ] *(out of scope: B-ALL/ETV6-RUNX1 inputs — only the cross-cohort CCI script.)*

## 2. Upstream primary-processing pipeline  ← the main Level-2 gap
The processed object (`…_with_PAC_new.h5ad`) is currently an **input**, not regenerated. To be
end-to-end reproducible, commit the code + parameters that built it:
- [ ] **Cell Ranger** command + version (FASTQ → `matrix_seurat/` counts), GRCh38 build.
- [ ] **scanpy processing notebook**: exact QC cutoffs, normalization, HVG, batch integration,
      Leiden resolution, UMAP params, Scrublet, cell-cycle, **cell-type annotation**, leukemic
      sub-clustering, and the **PAC definition rule**. *(These define `obs` labels used everywhere;
      the committed Seurat scripts re-import them rather than create them.)*
- [ ] **pySCENIC** nextflow config / CLI parameters (to regenerate `reg-py` from counts).

## 3. Script hygiene + a one-command runner
- [ ] **Re-path** the legacy scripts (remove `/Users/yakun/…`, `setwd("~/Downloads…")`) to
      package-relative paths (`count_matrix/`, `source_data/legacy/`).
- [ ] **Fix the known code errors:** `MCCF.R` line-29 syntax stub; `NCI_PAML_Gene_prognosis.R`
      `CSF` used before defined; `Figure_interaction_Chord_plot_B_ALL.R` `options(scipen=n)`;
      `ALSF_AML_Dendrogram_plotting.R` undefined-variable scratch file.
- [ ] **Split** the monolithic `scripts/legacy/server/*.R` into per-figure scripts (optional).
- [ ] Add a **master runner** (`run_all.sh`/Makefile) that executes scripts in dependency order
      and a **from-scratch check** (re-render, diff vs committed).

## 4. Environment + reference pinning  (partly done)
- [x] **Python** env snapshots: `env/requirements_figures_python.txt`, `env/requirements_census_python.txt`.
- [x] **Local venv note (which env runs what):** scanpy-based analysis scripts (`compute_*.py`, `rank_genes_groups`, pseudobulk) run under the project env **`ALSF_AML_2026_updated/.venv`** (scanpy 1.12.1); the lightweight figure/composite scripts (PIL/pandas/pymupdf/seaborn) run under **`07_Code/shared_CD96_pipeline/.venv`** (no scanpy). Scripts use `#!/usr/bin/env python`, so activate the matching env before running. The body-wide Census screen uses the isolated **`ALSF_AML_2026_updated/.census_venv`** (cellxgene-census 1.17.0 / tiledbsoma 2.3.0).
- [x] **R** package versions: `env/R_packages.txt`.
- [ ] Capture the **Sherlock conda envs** (`pyscenic_3.8`, the scanpy env) as `environment.yml`,
      and the pyscenic Singularity image hash (`aertslab-pyscenic-0.12.0.sif`).
- [ ] Pin **reference data + versions:** pyscenic hg38 cisTarget feather + motif `.tbl`; MSigDB
      release (C5:BP, GOCC_CELL_SURFACE v2023.1.Hs); CellChatDB (CellChat 2.2.0.9001);
      GRCh38 build; CSPA version.

## 5. Methods details to confirm (documentation, not code)  — see Methods.md [confirm] tags
- [ ] IRB/consent; library chemistry, sequencer, read depth; Cell Ranger version; exact scanpy QC
      cutoffs; PAC selection rule; legacy R `sessionInfo()`.

## 6. (Optional) statistical refinements — see Statistical_tests_and_packages.md §E
Principled survival cut-points + continuous Cox + covariate adjustment + FDR; GO expressed-gene
universe + q<0.05; Wilcoxon + FDR for composition; pseudobulk for patient-group DE. These change
figure values, so apply during the legacy rerun.

---

### What I can do now without the handoff
- Re-path + error-fix the legacy scripts whose inputs are already present (Figs 1-dendrogram[needs
  clinical], 2, S6, S7, S8) and add the master runner + CD96 re-render verification.
- Generate `environment.yml`/`renv.lock` stubs from the snapshots above.
### What only the lab/Yakun can provide
- The §1 files, the §2 upstream pipeline + Cell Ranger/scanpy params, and the §5 wet-lab details.

---

## Progress — do-now items completed (no handoff needed)
- [x] **Environment captured:** `env/requirements_figures_python.txt`, `env/requirements_census_python.txt`,
  `env/R_packages.txt`, `env/environment.yml`, `env/install_R_packages.R`, `env/README.md`.
- [x] **Master runner:** `run_all.sh` — renders CD96 figs 5–7 + S12–S14 from `source_data/` and
  rebuilds the combined PDF in one command (verified end-to-end: 5/5 scripts + 22-page PDF).
- [x] **Fig 2 GO re-pathed AND verified** — `scripts/legacy/Figure2/Figure2_clusterprofiler.R` now
  reads bundled `source_data/legacy/`, is msigdbr-version-robust, and **reproduces** all 18 GOBP
  dotplots + Module 15 → `figures/legacy/Figure2/`.
- [x] **Code errors fixed:** B-ALL `options(scipen=n)`→100; NCI premature `CSF` use commented;
  `MCCF.R` line-29 stub made to parse (+ NOTE); scratch `ALSF_AML_Dendrogram_plotting.R` flagged
  (canonical = `Figure_1_Dendrogram_plotting.R`). Both syntax-fixed files now parse.

Remaining (need handoff/lab): re-path the other legacy scripts once their inputs arrive (§1–2),
capture the Sherlock conda envs + reference-DB versions (§4), and the wet-lab/upstream details (§5).
