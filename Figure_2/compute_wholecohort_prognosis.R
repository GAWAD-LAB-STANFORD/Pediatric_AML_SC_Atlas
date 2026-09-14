#!/usr/bin/env Rscript
# Whole-cohort EFS prognosis for the 49 leukemic states: per-state univariate Cox on
# event-free survival (z-scored deconvolved fraction), BH-FDR. Defines the favorable/poor/
# ns groups used by Figure 2 panels B–E. Merges state annotation. No fabricated values.
suppressMessages(library(survival))
D3 <- "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/08_Source_Data/Figure_3/v3"
D2 <- "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/08_Source_Data/Figure_2/v3"
d  <- read.csv(file.path(D3, "LS_survival.csv"), check.names = FALSE)
ls <- grep("^LS_", names(d), value = TRUE)
dd <- d[!is.na(d$efs_time), ]
r <- data.frame()
for (x in ls) {
  fr <- dd[[x]]; if (sd(fr, na.rm = TRUE) == 0) next
  s <- summary(coxph(Surv(dd$efs_time, dd$efs_event) ~ scale(fr)[,1]))
  r <- rbind(r, data.frame(LS = x, HR = round(s$conf.int[1], 3),
             ci_lo = round(s$conf.int[3], 3), ci_hi = round(s$conf.int[4], 3), p = s$coefficients[5]))
}
r$fdr <- p.adjust(r$p, "BH"); r$dir <- ifelse(r$HR > 1, "poor", "favorable"); r <- r[order(r$p), ]
write.csv(r, file.path(D2, "LS_wholecohort_cox_efs.csv"), row.names = FALSE)
ann <- read.csv(file.path(D3, "LS_annotation.csv"))
m <- merge(r, ann, by = "LS", all.x = TRUE)
m$prognosis <- ifelse(m$fdr < 0.10, m$dir, "ns"); m <- m[order(m$p), ]
write.csv(m, file.path(D2, "LS_wholecohort_prognosis.csv"), row.names = FALSE)
cat(sprintf("n=%d, events=%d; FDR<0.10: %d states (%d favorable, %d poor)\n",
    nrow(dd), sum(dd$efs_event), sum(m$fdr < 0.10),
    sum(m$prognosis == "favorable"), sum(m$prognosis == "poor")))
