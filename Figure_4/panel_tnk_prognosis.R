#!/usr/bin/env Rscript
# Figure 4 prognosis panels: E) forest of per-SD EFS hazard ratios for the 12 deconvolved T/NK subsets
# (CIBERSORTx, TARGET n=1787), coloured favorable/poor at BH-FDR<0.10; F) EFS KM by GZMK-CD8 tertile
# (the strongest favorable subset) with GZMB-CD8 the opposite direction. Data = tnk_subset_survival.csv +
# the CIBERSORTx jobs merged with TARGET outcome. No fabricated values.
suppressMessages({library(ggplot2); library(dplyr); library(survival); library(survminer)})
source("/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/07_Code/Figure_3/manuscript_palette.R")
V2 <- "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/08_Source_Data/Figure_4/v2"
OUT <- "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/04_Main_Figures/Figure_4_panels"
SURV <- "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/08_Source_Data/Figure_3/v3/LS_survival.csv"

## ---- E: forest (EFS) ----
sv <- read.csv(file.path(V2, "tnk_subset_survival.csv")) %>% filter(outcome == "EFS")
sv$dir <- ifelse(sv$FDR >= 0.10, "n.s.", ifelse(sv$HR > 1, "poor", "favorable"))
sv$dir <- factor(sv$dir, levels = c("favorable","poor","n.s."))
sv <- sv %>% arrange(HR); sv$subset <- factor(sv$subset, levels = sv$subset)
cols <- c(favorable = unname(PROG["favorable"]), poor = unname(PROG["poor"]), `n.s.` = "#9aa0a6")
gE <- ggplot(sv, aes(HR, subset, colour = dir)) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey60") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.25, linewidth = 0.5) +
  geom_point(size = 2.6) +
  scale_colour_manual(values = cols, name = NULL) +
  scale_x_continuous(trans = "log2", breaks = c(0.7,0.85,1,1.2,1.4)) +
  labs(x = "EFS hazard ratio (per SD, TARGET n=1787)", y = NULL) +
  theme_classic(base_size = 11) + theme(legend.position = "top")
ggsave(file.path(OUT, "Figure_4E_tnk_forest.pdf"), gE, width = 5.4, height = 4.4)
ggsave(file.path(OUT, "Figure_4E_tnk_forest.png"), gE, width = 5.4, height = 4.4, dpi = 200)

## ---- F: EFS KM by GZMK-CD8 tertile ----
cib <- bind_rows(lapply(sprintf(file.path(V2,"CIBERSORTx_Job%d_Results.csv"),29:32), read.csv, check.names=FALSE)) %>%
       distinct(Mixture, .keep_all = TRUE)
surv <- read.csv(SURV, check.names = FALSE)
d <- merge(cib, surv[, c("sample","efs_time","efs_event")], by.x="Mixture", by.y="sample")
d <- d[!is.na(d$efs_time) & !is.na(d$efs_event), ]
g <- "GZMK CD8 T"
d$grp <- cut(d[[g]], quantile(d[[g]], c(0,1/3,2/3,1)), include.lowest=TRUE, labels=c("low","mid","high"))
fit <- survfit(Surv(efs_time, efs_event) ~ grp, data = d)
kp <- ggsurvplot(fit, data = d, palette = c(unname(PROG["poor"]), "#9aa0a6", unname(PROG["favorable"])),
                 conf.int = FALSE, pval = TRUE, pval.coord = c(0.2, 0.05), risk.table = TRUE,
                 risk.table.height = 0.28, legend.title = "GZMK CD8 T", legend.labs = c("low","mid","high"),
                 xlab = "Years (event-free survival)", ylab = "EFS probability",
                 ggtheme = theme_classic(base_size = 11))
pdf(file.path(OUT, "Figure_4F_tnk_km_GZMKCD8.pdf"), width = 5.0, height = 5.2); print(kp); dev.off()
png(file.path(OUT, "Figure_4F_tnk_km_GZMKCD8.png"), width = 1000, height = 1040, res = 200); print(kp); dev.off()
cat("wrote Figure_4E_tnk_forest + Figure_4F_tnk_km_GZMKCD8\n")
