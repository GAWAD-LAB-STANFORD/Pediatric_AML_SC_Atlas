#!/usr/bin/env Rscript
# KM (log-rank) + continuous Cox for EVERY state -- 49 leukemic (LS_*) + 7 normal
# compartments (NORM_*) -- in the whole cohort and within CBF, on EFS and OS.
# ONE consistent split rule for all states (no per-state cutoff shopping):
#   >=1/3 of patients zero -> Absent / Low / High (Low|High = median of non-zero);
#   otherwise tertiles (Low/Mid/High). Continuous Cox (per SD, no cutoff) is the primary,
#   binning-free statistic; the KM log-rank is the requested survival-curve test.
# FDR (BH) computed within the 49 leukemic states, separately from the 7 normal, per context.
suppressMessages(library(survival))
D3 <- "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/08_Source_Data/Figure_3/v3"
D2 <- "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/08_Source_Data/Figure_2/v3"
sv  <- read.csv(file.path(D3,"LS_survival.csv"), check.names=FALSE)
ann <- read.csv(file.path(D3,"LS_annotation.csv"))
states <- c(grep("^LS_",names(sv),value=TRUE), grep("^NORM_",names(sv),value=TRUE))

split_km <- function(x){
  if (mean(x==0) >= 1/3){
    nz <- x[x>0]; if (length(nz) < 6) return(list(g=NULL, type="too_sparse"))
    med <- median(nz)
    g <- factor(ifelse(x==0,"Absent", ifelse(x<=med,"Low","High")), levels=c("Absent","Low","High"))
    return(list(g=g, type="Absent/Low/High"))
  }
  q <- quantile(x, c(0,1/3,2/3,1))
  if (length(unique(q)) < 4){
    g <- factor(ifelse(x<=median(x),"Low","High"), levels=c("Low","High")); return(list(g=g, type="median"))
  }
  list(g=cut(x,q,include.lowest=TRUE,labels=c("Low","Mid","High")), type="tertile")
}
logrank_p <- function(g,t,e,dat){ if (is.null(g) || nlevels(droplevels(g))<2) return(NA)
  sd <- tryCatch(survdiff(Surv(dat[[t]],dat[[e]])~g), error=function(er) NULL)
  if (is.null(sd)) return(NA); 1-pchisq(sd$chisq, length(sd$n)-1) }
coxstat <- function(x,t,e,dat){ if (sd(x)==0) return(c(NA,NA))
  s <- tryCatch(summary(coxph(Surv(dat[[t]],dat[[e]])~scale(x), data=dat)), error=function(er) NULL)
  if (is.null(s)) return(c(NA,NA)); c(s$conf.int[1], s$coefficients[5]) }

rowfor <- function(s, dat, t, e){
  d <- dat[!is.na(dat[[t]]) & !is.na(dat[[s]]), ]; x <- d[[s]]
  sp <- split_km(x); cx <- coxstat(x,t,e,d)
  data.frame(n=nrow(d), det=round(mean(x>0)*100,1),
             coxHR=round(cx[1],3), coxP=cx[2],
             KMsplit=sp$type, KM_logrankP=round(logrank_p(sp$g,t,e,d),4))
}
cbf <- function(x) x[x$subtype2 %in% c("t(8;21)","inv(16)"), ]

out <- data.frame()
for (s in states){
  cls <- ifelse(grepl("^LS_",s),"leukemic","normal")
  prog <- if (cls=="leukemic") ann$lineage_program[match(s,ann$LS)] else sub("NORM_","",s)
  npt  <- if (cls=="leukemic") ann$n_patients[match(s,ann$LS)] else NA
  we <- rowfor(s, sv, "efs_time","efs_event"); wo <- rowfor(s, sv, "os_time","os_event")
  ce <- rowfor(s, cbf(sv), "efs_time","efs_event"); co <- rowfor(s, cbf(sv), "os_time","os_event")
  out <- rbind(out, data.frame(state=s, class=cls, program=prog, scRNA_pts=npt,
    WC_EFS_n=we$n, WC_EFS_det=we$det, WC_EFS_coxHR=we$coxHR, WC_EFS_coxP=we$coxP,
    WC_EFS_KMsplit=we$KMsplit, WC_EFS_logrankP=we$KM_logrankP,
    WC_OS_coxHR=wo$coxHR, WC_OS_coxP=wo$coxP, WC_OS_logrankP=wo$KM_logrankP,
    CBF_EFS_n=ce$n, CBF_EFS_coxHR=ce$coxHR, CBF_EFS_coxP=ce$coxP, CBF_EFS_logrankP=ce$KM_logrankP,
    CBF_OS_coxHR=co$coxHR, CBF_OS_coxP=co$coxP, CBF_OS_logrankP=co$KM_logrankP))
}
# FDR within leukemic vs normal, per context (on the continuous Cox p, the primary)
for (grp in c("leukemic","normal")){ m <- out$class==grp
  out$WC_EFS_coxFDR[m]  <- p.adjust(out$WC_EFS_coxP[m],"BH")
  out$CBF_EFS_coxFDR[m] <- p.adjust(out$CBF_EFS_coxP[m],"BH")
  out$WC_EFS_logrankFDR[m]  <- p.adjust(out$WC_EFS_logrankP[m],"BH")
  out$CBF_EFS_logrankFDR[m] <- p.adjust(out$CBF_EFS_logrankP[m],"BH")
}
out$dir_WC_EFS <- ifelse(out$WC_EFS_coxHR>1,"poor","favorable")
num <- sapply(out, is.numeric); out[num] <- lapply(out[num], function(z) signif(z,3))
out <- out[order(out$class, out$WC_EFS_coxP), ]
write.csv(out, file.path(D2,"KM_all_states_EFS_OS.csv"), row.names=FALSE)
cat("wrote KM_all_states_EFS_OS.csv (", nrow(out), " states x whole-cohort+CBF, EFS+OS)\n\n", sep="")
cat("=== WHOLE-COHORT EFS: states significant by KM log-rank (FDR<0.10) ===\n")
w <- out[out$WC_EFS_logrankFDR<0.10 & !is.na(out$WC_EFS_logrankFDR), c("state","class","program","dir_WC_EFS","WC_EFS_coxHR","WC_EFS_coxP","WC_EFS_KMsplit","WC_EFS_logrankP","WC_EFS_logrankFDR")]
print(w, row.names=FALSE)
cat("\n=== WITHIN-CBF EFS: states significant by KM log-rank (raw p<0.05) ===\n")
c2 <- out[out$CBF_EFS_logrankP<0.05 & !is.na(out$CBF_EFS_logrankP), c("state","class","program","CBF_EFS_coxHR","CBF_EFS_coxP","CBF_EFS_logrankP","CBF_EFS_logrankFDR")]
print(c2, row.names=FALSE)
