#!/usr/bin/env Rscript
# Figure 4: EFS Kaplan-Meier for each prognostic T/NK subset, split high vs low deconvolved fraction
# (median split), so the DIRECTION is explicit -- you see whether the high (red) or low (teal) group has
# worse survival, i.e. whether MORE or LESS of that subset is bad. Only the FDR<0.10 (EFS) subsets are
# shown. Data = CIBERSORTx Jobs 29-32 + TARGET outcome. No fabricated values.
suppressMessages({library(survival); library(survminer); library(dplyr); library(ggplot2)})
source("/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/07_Code/Figure_3/manuscript_palette.R")
V2  <- "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/08_Source_Data/Figure_4/v2"
OUT <- "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/04_Main_Figures/Figure_4_panels"
SURV<- "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/08_Source_Data/Figure_3/v3/LS_survival.csv"
cib <- bind_rows(lapply(sprintf(file.path(V2,"CIBERSORTx_Job%d_Results.csv"),29:32), read.csv, check.names=FALSE)) %>%
       distinct(Mixture, .keep_all=TRUE)
surv<- read.csv(SURV, check.names=FALSE)
d <- merge(cib, surv[,c("sample","efs_time","efs_event")], by.x="Mixture", by.y="sample")
d <- d[!is.na(d$efs_time)&!is.na(d$efs_event),]
sv <- read.csv(file.path(V2,"tnk_subset_survival.csv"))
sig <- sv %>% filter(outcome=="EFS", FDR<0.10) %>% arrange(HR) %>% pull(subset)   # favorable -> poor
km <- list(); ann <- list()
for (s in sig) {
  x <- d[[s]]; med <- median(x)
  d$grp <- if (med > 0) ifelse(x > med, "high","low") else ifelse(x > 0, "high","low")
  fit <- survfit(Surv(efs_time, efs_event) ~ grp, data=d)
  ss <- surv_summary(fit, data=d); ss$strata <- gsub("grp=", "", as.character(ss$strata))
  ss$subset <- s; km[[s]] <- ss
  p <- 1 - pchisq(survdiff(Surv(efs_time, efs_event) ~ grp, data=d)$chisq, 1)
  hr <- sv$HR[sv$outcome=="EFS" & sv$subset==s]
  ann[[s]] <- data.frame(subset=s, lab=sprintf("%s\n%s → worse   p=%s",
              s, ifelse(hr>1,"more","less"), format.pval(p, digits=1, eps=1e-4)))
}
kmdf <- bind_rows(km); anndf <- bind_rows(ann)
kmdf$subset <- factor(kmdf$subset, levels=sig); anndf$subset <- factor(anndf$subset, levels=sig)
g <- ggplot(kmdf, aes(time, surv, colour=strata)) +
  geom_step(linewidth=0.6) +
  facet_wrap(~subset, ncol=4) +
  scale_colour_manual(values=c(high=unname(PROG["poor"]), low=unname(PROG["favorable"])),
                      labels=c("high fraction","low fraction"), name=NULL) +
  geom_text(data=anndf, aes(x=0, y=0.06, label=lab), inherit.aes=FALSE, hjust=0, size=2.7, lineheight=0.9) +
  labs(x="Years (event-free survival)", y="EFS probability") +
  coord_cartesian(ylim=c(0,1)) +
  theme_classic(base_size=10) + theme(legend.position="top", strip.text=element_text(size=8, face="bold"))
ggsave(file.path(OUT,"Figure_4F_tnk_km_grid.pdf"), g, width=8.6, height=4.6)
ggsave(file.path(OUT,"Figure_4F_tnk_km_grid.png"), g, width=8.6, height=4.6, dpi=200)
cat("wrote Figure_4F_tnk_km_grid for", length(sig), "significant subsets:", paste(sig, collapse=", "), "\n")
