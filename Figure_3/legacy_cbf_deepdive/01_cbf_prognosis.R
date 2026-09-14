#!/usr/bin/env Rscript
# Figure 3 step 1: within-CBF event-free-survival Cox per leukemic state.
# Picks the CBF-prognostic states (p < PCUT) and the top poor state ("focus state",
# used for the KM and the targetable panel) DYNAMICALLY from the data.
suppressMessages(library(survival))
args <- commandArgs(trailingOnly = TRUE)
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "config.R"))
d <- read.csv(F3$SURV, check.names = FALSE)
if (!all(c("efs_time","efs_event") %in% names(d))) {
  ef <- read.csv(F3$EFSREF)[, c("sample","efs_time","efs_event")]
  d <- merge(d, ef, by = "sample")
}
cbf <- d[d$subtype2 %in% F3$CBF & !is.na(d$efs_time), ]
lscols <- grep("^LS_[0-9]+$", names(cbf), value = TRUE)
cat(sprintf("CBF patients=%d, EFS events=%d, states=%d\n", nrow(cbf), sum(cbf$efs_event), length(lscols)))
res <- data.frame()
for (x in lscols) { fr <- cbf[[x]]; if (sd(fr) == 0 || sum(fr > 0) < 10) next
  s <- summary(coxph(Surv(efs_time, efs_event) ~ scale(fr)[,1], data = cbf))
  res <- rbind(res, data.frame(LS = x, HR = s$conf.int[1], p = s$coefficients[5])) }
res$dir <- ifelse(res$HR > 1, "poor", "favorable")
write.csv(res[order(res$p), ], file.path(F3$DATA, "LS_cbf_cox_efs.csv"), row.names = FALSE)
pr <- res[res$p < F3$PCUT, ]
ann <- read.csv(F3$ANNOT); ann <- ann[, c("LS","lineage_program")]
pr <- merge(pr, ann, by = "LS"); pr$group <- ifelse(pr$dir == "favorable","Favorable","Poor")
pr$LSc <- gsub("_","",pr$LS)
pr$lab <- ifelse(pr$lineage_program %in% c("unresolved", NA), pr$LSc, paste0(pr$LSc,"  ",pr$lineage_program))
pr <- pr[order(pr$dir == "poor", pr$p), ]; pr$lab <- make.unique(pr$lab)
write.csv(pr[, c("LS","lab","group","dir","HR","p")], file.path(F3$DATA, "LS_cbf_order.csv"), row.names = FALSE)
focus <- pr$LS[pr$dir == "poor"][which.min(pr$p[pr$dir == "poor"])]
writeLines(focus, file.path(F3$DATA, "LS_focus_state.txt"))
cat(sprintf("CBF-prognostic states: %d (%d favorable, %d poor); FOCUS (top poor) = %s\n",
    nrow(pr), sum(pr$dir=="favorable"), sum(pr$dir=="poor"), focus))
