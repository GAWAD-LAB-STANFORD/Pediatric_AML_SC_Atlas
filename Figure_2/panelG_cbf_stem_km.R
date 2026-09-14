#!/usr/bin/env Rscript
# Figure 2H (v3): within CBF, the HSC/MPP PROGRAM (LS25 up-genes, bulk-scored on TARGET,
# median split) stratifies event-free survival -> HSC/MPP-like CBF relapses. Rigorous,
# program-level version of the signal (the deconvolved normal-HSC fraction is too
# zero-inflated to test). Score = mean z of the HSC/MPP program genes (LS25) on bulk CPM.
suppressMessages({library(survival); library(survminer)})
CODE <- "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/07_Code/Figure_2"
source(file.path(CODE, "manuscript_palette.R"))
D2  <- "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/08_Source_Data/Figure_2/v3"
OUT <- "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/04_Main_Figures/Figure_2_panels"

sc <- read.csv(file.path(D2, "TARGET_bulk_program_scores.csv"))
sv <- read.csv(file.path(D2, "..", "..", "Figure_3", "v3", "LS_survival.csv"), check.names=FALSE)
d <- merge(sv[, c("sample","subtype2","efs_time","efs_event")], sc, by="sample")
cbf <- d[d$subtype2 %in% c("t(8;21)","inv(16)") & !is.na(d$efs_time), ]
# principled median split (pre-specified, not optimized); tertiles underpower the graded effect
cbf$grp <- factor(ifelse(cbf$HSC_program > median(cbf$HSC_program), "High", "Low"), levels=c("Low","High"))
sf <- survfit(Surv(efs_time, efs_event) ~ grp, data = cbf)
sd <- survdiff(Surv(efs_time, efs_event) ~ grp, data = cbf); lr <- 1-pchisq(sd$chisq, length(sd$n)-1)
cx <- summary(coxph(Surv(efs_time,efs_event)~scale(HSC_program), data=cbf))
p <- ggsurvplot(sf, data=cbf, palette=c(unname(PROG["favorable"]), unname(PROG["poor"])), conf.int=FALSE, pval=TRUE,
   pval.coord=c(0.2,0.05), risk.table=TRUE, legend.title="HSC/MPP program", legend.labs=c("Low","High"),
   xlab="Years (event-free survival)", ylab="EFS probability", title="", risk.table.height=0.27, break.time.by=2)
png(file.path(OUT,"Figure_2G_cbf_stem_km.png"), width=6.2, height=5.6, units="in", res=200); print(p); dev.off()
pdf(file.path(OUT,"Figure_2G_cbf_stem_km.pdf"), width=6.2, height=5.6); print(p); dev.off()
cat(sprintf("within-CBF HSC/MPP-program median-split KM: n=%d, log-rank p=%.4f; Cox HR=%.3f p=%.4f; group n: %s\n",
    nrow(cbf), lr, cx$conf.int[1], cx$coefficients[5], paste(table(cbf$grp), collapse="/")))
