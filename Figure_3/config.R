# Figure 3 (CBF deep-dive) — shared config. To swap datasets, set FIG3_DATA (a dir
# holding LS_survival.csv, LS_annotation.csv, meancpm.npy, genes.npy,
# reference_assignments.csv, jaccard_percluster.csv) or edit BASE below.
suppressMessages({library(ggplot2); library(scales)})
BASE <- "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster"
C3   <- file.path(BASE, "07_Code/Figure_3")
F3 <- new.env()
F3$DATA   <- Sys.getenv("FIG3_DATA", file.path(BASE, "08_Source_Data/Figure_3/v3"))
F3$OUT    <- Sys.getenv("FIG3_OUT",  file.path(BASE, "04_Main_Figures/Figure_3_panels"))
F3$SURV   <- file.path(F3$DATA, "LS_survival.csv")     # per-patient LS_* fractions + subtype2 (+os/efs); from CIBERSORTx stitch
F3$EFSREF <- file.path(BASE, "08_Source_Data/Figure_2/Figure2E_TARGET_survival_input.csv")  # sample,efs_time,efs_event
F3$ANNOT  <- file.path(F3$DATA, "LS_annotation.csv")   # LS,n_cells,lineage_program,stem_z,top_cyto,...
F3$MEANCPM<- file.path(F3$DATA, "meancpm.npy"); F3$GENES <- file.path(F3$DATA, "genes.npy")
F3$REFA   <- file.path(F3$DATA, "reference_assignments.csv")
F3$JPC    <- file.path(F3$DATA, "jaccard_percluster.csv")
F3$UPGENES<- file.path(F3$DATA, "LS_upgenes_wilcoxon.csv")   # per-state Wilcoxon up-genes (from 04_compute_DE.py)
F3$CBF    <- c("t(8;21)","inv(16)")   # CBF cytogenetic subtypes
F3$PCUT   <- 0.10                       # within-CBF EFS Cox p cutoff for "prognostic" states
dir.create(F3$OUT, showWarnings = FALSE, recursive = TRUE)
source(file.path(C3, "manuscript_palette.R"))
