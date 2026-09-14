#!/usr/bin/env Rscript
# Figure 4 prognosis: deconvolved T/NK-subset fractions (CIBERSORTx, TARGET bulk, 4 batch jobs) tested
# against event-free and overall survival (per-SD Cox), BH-FDR across the 12 subsets. Merges the CIBERSORTx
# results with the same TARGET outcome table used for the leukemic states. No fabricated values.
suppressMessages({library(survival); library(dplyr)})
V2 <- "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/08_Source_Data/Figure_4/v2"
SURV <- "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/08_Source_Data/Figure_3/v3/LS_survival.csv"
jobs <- lapply(sprintf(file.path(V2, "CIBERSORTx_Job%d_Results.csv"), 29:32),
               read.csv, check.names = FALSE)
cib <- bind_rows(jobs) %>% distinct(Mixture, .keep_all = TRUE)
subsets <- setdiff(colnames(cib), c("Mixture","P-value","Correlation","RMSE"))
surv <- read.csv(SURV, check.names = FALSE)
d <- merge(cib, surv[, c("sample","efs_time","efs_event","os_time","os_event","subtype2")],
           by.x = "Mixture", by.y = "sample")
cat("merged samples:", nrow(d), "| subsets:", length(subsets), "\n\n")
res <- list()
for (oc in c("efs","os")) {
  s <- d[!is.na(d[[paste0(oc,"_time")]]) & !is.na(d[[paste0(oc,"_event")]]), ]
  rows <- lapply(subsets, function(g) {
    s$z <- as.numeric(scale(s[[g]])); s$time <- s[[paste0(oc,"_time")]]; s$event <- s[[paste0(oc,"_event")]]
    co <- summary(coxph(Surv(time, event) ~ z, data = s))$coef
    data.frame(subset = g, outcome = toupper(oc), n = nrow(s), events = sum(s$event),
               HR = exp(co[1,1]), lo = exp(co[1,1]-1.96*co[1,3]), hi = exp(co[1,1]+1.96*co[1,3]), p = co[1,5])
  })
  r <- bind_rows(rows); r$FDR <- p.adjust(r$p, "BH"); res[[oc]] <- r
}
out <- bind_rows(res)
write.csv(out, file.path(V2, "tnk_subset_survival.csv"), row.names = FALSE)
for (oc in c("EFS","OS")) {
  cat("===", oc, "(per-SD Cox, BH-FDR across 12 subsets) ===\n")
  print(out %>% filter(outcome==oc) %>% arrange(p) %>%
        mutate(across(c(HR,lo,hi), ~round(.,2)), p=signif(p,2), FDR=signif(FDR,2)) %>%
        select(subset,HR,lo,hi,p,FDR), row.names = FALSE)
  cat("\n")
}
