#!/usr/bin/env Rscript
# S11 (replaces the stale AML_1..28 subcluster forest): whole-cohort EFS hazard ratios for all 49
# leukemic states from univariate Cox on the z-scored deconvolved state fraction (per SD), ordered
# by HR and coloured by prognosis (BH-FDR<0.10: favorable HR<1 / poor HR>1 / n.s.). Data =
# LS_wholecohort_cox_efs.csv + LS_annotation.csv (lineage labels). No fabricated values.
suppressMessages({library(ggplot2); library(dplyr)})
D2  <- "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/08_Source_Data/Figure_2/v3"
D3  <- "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/08_Source_Data/Figure_3/v3"
OUT <- "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/05_Supplementary_Figures"
cox <- read.csv(file.path(D2,"LS_wholecohort_cox_efs.csv"))
ann <- read.csv(file.path(D3,"LS_annotation.csv"))[,c("LS","lineage_program")]
d <- merge(cox, ann, by="LS", all.x=TRUE)
d$prognosis <- with(d, ifelse(fdr<0.10 & HR>1, "poor", ifelse(fdr<0.10 & HR<1, "favorable", "n.s.")))
d$prognosis <- factor(d$prognosis, levels=c("favorable","poor","n.s."))
d$lineage_program[is.na(d$lineage_program) | d$lineage_program=="unresolved"] <- ""
d$label <- trimws(paste0(sub("_","",d$LS), "  ", d$lineage_program))
d <- d[order(d$HR),]; d$label <- factor(d$label, levels=d$label)
pcol <- c(favorable="#00798c", poor="#d1495b", `n.s.`="#9aa0a6")
p <- ggplot(d, aes(HR, label, colour=prognosis)) +
  geom_vline(xintercept=1, linetype="dashed", colour="grey55") +
  geom_errorbarh(aes(xmin=ci_lo, xmax=ci_hi), height=0, linewidth=0.5) +
  geom_point(size=2) +
  scale_colour_manual(values=pcol, name="EFS prognosis\n(FDR<0.10)") +
  scale_x_log10(breaks=c(0.6,0.8,1.0,1.25,1.5)) +
  labs(x="whole-cohort EFS hazard ratio per SD of deconvolved fraction (log scale)", y=NULL,
       title="EFS hazard ratio for all 49 leukemic states (whole TARGET cohort)") +
  theme_minimal(base_size=10) +
  theme(panel.grid.minor=element_blank(), axis.text.y=element_text(size=6.5),
        plot.title=element_text(size=11, face="bold"), legend.position="right")
ggsave(file.path(OUT,"Figure_S11__state_EFS_forest.pdf"), p, width=7.5, height=9)
ggsave(file.path(OUT,"Figure_S11__state_EFS_forest.png"), p, width=7.5, height=9, dpi=200)
cat(sprintf("wrote S11 49-state EFS forest (favorable=%d, poor=%d, n.s.=%d)\n",
    sum(d$prognosis=="favorable"), sum(d$prognosis=="poor"), sum(d$prognosis=="n.s.")))
