#!/usr/bin/env Rscript
# Figure 4F: clonal expansion measured PER PATIENT (the correct unit; clonality is a per-repertoire
# property). Because per-patient per-subset TCR counts are small (few tens of cells), we use the fraction
# of a patient's TCR+ cells in expanded clones (clone size >= 2) rather than Shannon entropy, which is
# unreliable at low n. Each point is one patient (>=15 TCR+ cells in that subset); AML circles coloured by
# EFS prognosis, healthy-BM donors as purple diamonds. GZMB-CD8 is the most clonally expanded.
# NOTE (reported in text): overall per-patient CD8/T clonality does NOT differ AML vs HBM - the expansion
# is a property of the effector-CD8 subset, present in both. Data = tcr_expanded_perpatient.csv. No fabrication.
suppressMessages({library(ggplot2); library(dplyr); library(scales)})
source("/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/07_Code/Figure_3/manuscript_palette.R")
V2 <-"/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/08_Source_Data/Figure_4/v2"
OUT<-"/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/04_Main_Figures/Figure_4_panels"
allpp <- read.csv(file.path(V2,"tcr_expanded_perpatient.csv"))
pp <- allpp %>% filter(group=="AML")
sv <- read.csv(file.path(V2,"tnk_subset_survival.csv")) %>% filter(outcome=="EFS") %>%
      mutate(dir=ifelse(FDR<0.10, ifelse(HR>1,"poor","favorable"),"n.s.")) %>% select(subset,dir)
np <- pp %>% count(subset) %>% filter(n>=3)
hbm <- allpp %>% filter(group=="HBM", subset %in% np$subset)
pp <- pp %>% filter(subset %in% np$subset) %>% left_join(sv,by="subset")
pp$dir[is.na(pp$dir)] <- "n.s."; pp$dir <- factor(pp$dir, levels=c("favorable","poor","n.s."))
ord <- pp %>% group_by(subset) %>% summarise(m=median(frac_expanded)) %>% arrange(m) %>% pull(subset)
pp$subset <- factor(pp$subset, levels=ord); hbm$subset <- factor(hbm$subset, levels=ord)
lab <- pp %>% count(subset)
PC <- c(favorable=unname(PROG["favorable"]),poor=unname(PROG["poor"]),`n.s.`="#9aa0a6")
g <- ggplot(pp, aes(frac_expanded*100, subset)) +
  geom_boxplot(aes(fill=dir), outlier.shape=NA, width=.55, alpha=.35, linewidth=.4, colour="grey40") +
  geom_jitter(aes(colour=dir), height=.12, size=1.9, alpha=.85) +
  geom_point(data=hbm, aes(frac_expanded*100, subset), shape=23, fill="#7b3294", colour="grey20", size=2.8, stroke=.4) +
  geom_text(data=lab, aes(x=101, y=subset, label=paste0("n=",n)), hjust=0, size=2.7, colour="grey35") +
  scale_fill_manual(values=PC, name="EFS prognosis") + scale_colour_manual(values=PC, guide="none") +
  scale_x_continuous(limits=c(0,112), breaks=c(0,25,50,75,100), expand=c(0,0)) +
  labs(x="% of TCR+ cells in expanded clones (per patient)", y=NULL,
       caption="AML patients = circles (EFS prognosis); healthy-BM donors = purple diamonds") +
  theme_classic(base_size=11) + theme(legend.position=c(.98,.05), legend.justification=c(1,0),
        plot.caption=element_text(size=8, hjust=0))
ggsave(file.path(OUT,"Figure_4F_tnk_clonality_perpatient.png"), g, width=5.8, height=4.2, dpi=200)
cat("wrote Figure_4F (per-patient % expanded); subsets:", paste(ord,collapse=", "), "\n")
