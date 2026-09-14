#!/usr/bin/env Rscript
# Render the whole-cohort EFS Kaplan-Meier for a single leukemic state, using the SAME
# zero-aware split rule as the KM table (tertiles when well-detected; Absent/Low/High when
# >1/3 zeros), and annotate BOTH the log-rank p and the binning-free continuous Cox HR/p.
# Usage: Rscript km_state.R LS_45
suppressMessages({library(survival); library(survminer)})
CODE <- "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/07_Code/Figure_3"
source(file.path(CODE, "manuscript_palette.R"))
D3 <- "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/08_Source_Data/Figure_3/v3"
OUT <- "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/04_Main_Figures/Figure_2_panels"
args <- commandArgs(trailingOnly=TRUE); ST <- if (length(args)>=1) args[1] else "LS_45"

d <- read.csv(file.path(D3,"LS_survival.csv"), check.names=FALSE)
d <- d[!is.na(d$efs_time), ]; x <- d[[ST]]
cx <- summary(coxph(Surv(efs_time,efs_event)~scale(x), data=d))
if (mean(x==0) >= 1/3) {                       # zero-inflated
  med <- median(x[x>0]); d$grp <- factor(ifelse(x==0,"Absent",ifelse(x<=med,"Low","High")), levels=c("Absent","Low","High"))
  pal <- c("#c9c9c9", unname(PROG["favorable"]), unname(PROG["poor"])); split <- "Absent/Low/High"
} else {                                        # tertiles
  d$grp <- cut(x, quantile(x,c(0,1/3,2/3,1)), include.lowest=TRUE, labels=c("Low","Mid","High"))
  pal <- c("#0A9396","#EE9B00","#AE2012"); split <- "tertile"
}
sf <- survfit(Surv(efs_time,efs_event)~grp, data=d)
sd <- survdiff(Surv(efs_time,efs_event)~grp, data=d); lr <- 1-pchisq(sd$chisq, length(sd$n)-1)
sub <- sprintf("continuous Cox: HR=%.3f per SD, p=%.1e  (log-rank shown on plot)", cx$conf.int[1], cx$coefficients[5])
p <- ggsurvplot(sf, data=d, palette=pal, conf.int=FALSE, pval=TRUE, pval.coord=c(0.2,0.05),
   risk.table=TRUE, legend.title=sprintf("%s (%s)", gsub("_","",ST), split), legend.labs=levels(d$grp),
   xlab="Years (event-free survival)", ylab="EFS probability", title="", risk.table.height=0.27, break.time.by=2)
# continuous-Cox stats reported in the figure legend, not on the panel (log-rank p shown on plot)
fn <- file.path(OUT, sprintf("KM_%s_wholecohort_EFS", ST))
png(paste0(fn,".png"), width=6.2, height=5.8, units="in", res=200); print(p); dev.off()
pdf(paste0(fn,".pdf"), width=6.2, height=5.8); print(p); dev.off()
cat(sprintf("%s whole-cohort EFS: %s split; log-rank p=%.4f; continuous Cox HR=%.3f p=%.2e; n/grp: %s\n",
    ST, split, lr, cx$conf.int[1], cx$coefficients[5], paste(table(d$grp), collapse="/")))
