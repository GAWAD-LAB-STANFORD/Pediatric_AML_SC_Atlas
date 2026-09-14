#!/usr/bin/env Rscript
# Figure 4 TCR panel: clonal expansion per T/NK subset (fraction of TCR+ cells in clones of size >=2,
# clonotype = productive CDR3 kept patient-private), coloured by the subset's whole-cohort EFS prognosis
# direction. Clonal expansion concentrates in the cytotoxic CD8 subsets, and the adverse GZMB-CD8 is the
# most clonally expanded. Only subsets with >=20 TCR+ cells shown. Data = tcr_by_subset.csv +
# tnk_subset_survival.csv. No fabricated values.
suppressMessages({library(ggplot2); library(dplyr)})
source("/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/07_Code/Figure_3/manuscript_palette.R")
V2  <- "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/08_Source_Data/Figure_4/v2"
OUT <- "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/04_Main_Figures/Figure_4_panels"
tcr <- read.csv(file.path(V2, "tcr_by_subset.csv"))
sv  <- read.csv(file.path(V2, "tnk_subset_survival.csv")) %>% filter(outcome=="EFS") %>%
       mutate(dir = ifelse(FDR<0.10, ifelse(HR>1,"poor","favorable"), "n.s.")) %>% select(subset, dir)
d <- merge(tcr, sv, by="subset", all.x=TRUE); d$dir[is.na(d$dir)] <- "n.s."
d <- d %>% filter(n_tcr >= 20) %>% arrange(frac_expanded)
d$subset <- factor(d$subset, levels = d$subset)
d$dir <- factor(d$dir, levels = c("favorable","poor","n.s."))
cols <- c(favorable=unname(PROG["favorable"]), poor=unname(PROG["poor"]), `n.s.`="#9aa0a6")
g <- ggplot(d, aes(frac_expanded*100, subset, fill=dir)) +
  geom_col(width=0.7, colour="grey30", linewidth=0.2) +
  geom_text(aes(label=paste0("n=",n_tcr)), hjust=-0.15, size=2.9, colour="grey30") +
  scale_fill_manual(values=cols, name="EFS prognosis") +
  scale_x_continuous(expand=expansion(mult=c(0,0.18))) +
  labs(x="% of TCR+ cells in expanded clones (size >= 2)", y=NULL) +
  theme_classic(base_size=11) +
  theme(legend.position=c(0.98,0.05), legend.justification=c(1,0),
        legend.background=element_rect(fill="white", colour=NA))
ggsave(file.path(OUT, "Figure_4F_tnk_tcr_expansion.pdf"), g, width=5.6, height=4.4)
ggsave(file.path(OUT, "Figure_4F_tnk_tcr_expansion.png"), g, width=5.6, height=4.4, dpi=200)
cat("wrote Figure_4F_tnk_tcr_expansion\n")
