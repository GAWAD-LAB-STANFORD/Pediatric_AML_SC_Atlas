#!/usr/bin/env Rscript
# Figure 2G: LS12 (stem-like) tertiles within CBF AML, event-free survival.
# Tertiles of the deconvolved LS12 fraction among CBF patients (t(8;21)+inv(16));
# Kaplan-Meier of EFS with two-sided log-rank.
suppressMessages({library(survival); library(survminer)})
SL  <- Sys.getenv("FIG2_DATA", "/private/tmp/claude-503/-Users-chuckgawad-Desktop-ALSF-AML-2026-updated/f5e93d65-31d5-411e-b754-6f7c504928b7/scratchpad")
PKG <- "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated"
OUT <- file.path(PKG, "__SUBMISSION_PACKAGE_v2_recluster/04_Main_Figures/Figure_2_new_panels")
d  <- read.csv(file.path(SL, "target_ls/LS48_survival.csv"))
ef <- read.csv(file.path(PKG, "__SUBMISSION_PACKAGE/08_Source_Data/Figure_2/Figure2E_TARGET_survival_input.csv"))[, c("sample","efs_time","efs_event")]
d  <- merge(d, ef, by = "sample")
cbf <- d[d$subtype2 %in% c("t(8;21)","inv(16)") & !is.na(d$efs_time), ]
cbf$LS12 <- cut(cbf$LS_12, quantile(cbf$LS_12, c(0, 1/3, 2/3, 1)), include.lowest = TRUE,
                labels = c("Low","Mid","High"))
sf <- survfit(Surv(efs_time, efs_event) ~ LS12, data = cbf)
p <- ggsurvplot(sf, data = cbf, palette = c("#0A9396","#EE9B00","#AE2012"), conf.int = FALSE,
   pval = TRUE, pval.coord = c(0.2, 0.05), risk.table = TRUE, legend.title = "LS12 tertile",
   legend.labs = c("Low","Mid","High"), xlab = "Years (event-free survival)",
   ylab = "EFS probability", title = "", risk.table.height = 0.26, break.time.by = 2)
png(file.path(OUT, "Figure_2_LS12_CBF_KM.png"), width = 6.5, height = 5.6, units = "in", res = 200); print(p); dev.off()
pdf(file.path(OUT, "Figure_2_LS12_CBF_KM.pdf"), width = 6.5, height = 5.6); print(p); dev.off()
cat("wrote Figure_2G LS12 CBF KM (EFS)\n")
