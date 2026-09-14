#!/usr/bin/env Rscript
# Figure 3 panel B: tertile Kaplan-Meier (event-free survival) within CBF for the
# FOCUS state (top within-CBF poor predictor, chosen dynamically in step 1).
suppressMessages({library(survival); library(survminer)})
source(file.path(dirname(sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE))),"config.R"))
focus <- readLines(file.path(F3$DATA,"LS_focus_state.txt"))[1]
d <- read.csv(F3$SURV, check.names=FALSE)
if (!all(c("efs_time","efs_event") %in% names(d))) {
  ef <- read.csv(F3$EFSREF)[,c("sample","efs_time","efs_event")]; d <- merge(d,ef,by="sample") }
cbf <- d[d$subtype2 %in% F3$CBF & !is.na(d$efs_time), ]
x <- cbf[[focus]]
if (mean(x == 0) >= 1/3) {                 # zero-inflated: Absent / Low / High (median of non-zero)
  med <- median(x[x > 0])
  cbf$grp <- factor(ifelse(x == 0, "Absent", ifelse(x <= med, "Low", "High")), levels = c("Absent","Low","High"))
  glabs <- paste0("Absent"); glabs <- c("Absent","Low","High"); split_desc <- "abundance"
} else {                                    # tertiles
  cbf$grp <- cut(x, quantile(x, c(0,1/3,2/3,1)), include.lowest = TRUE, labels = c("Low","Mid","High"))
  glabs <- c("Low","Mid","High"); split_desc <- "tertile"
}
sf <- survfit(Surv(efs_time, efs_event) ~ grp, data = cbf)
lab <- gsub("_","",focus)
p <- ggsurvplot(sf, data=cbf, palette=c("#0A9396","#EE9B00","#AE2012"), conf.int=FALSE, pval=TRUE,
   pval.coord=c(0.2,0.05), risk.table=TRUE, legend.title=paste(lab,split_desc), legend.labs=glabs,
   xlab="Years (event-free survival)", ylab="EFS probability", title="", risk.table.height=0.26, break.time.by=2)
png(file.path(F3$OUT,"panelB_focus_km.png"), width=6.5, height=5.6, units="in", res=200); print(p); dev.off()
pdf(file.path(F3$OUT,"panelB_focus_km.pdf"), width=6.5, height=5.6); print(p); dev.off()
sd <- survdiff(Surv(efs_time,efs_event)~grp, data=cbf)
cat(sprintf("wrote panelB_focus_km for %s (%s split); log-rank p=%.4f; group n: %s\n",
    focus, split_desc, 1-pchisq(sd$chisq, length(sd$n)-1), paste(table(cbf$grp), collapse="/")))
