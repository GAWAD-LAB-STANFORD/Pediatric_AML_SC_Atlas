#!/usr/bin/env Rscript
# Figure 3 (LS19) Monocle3 panel: Monocle3 pseudotime (rooted at normal HSPC) of LS19 vs normal
# HSPC (progenitor ref), normal DC/pDC (mature ref), and other leukemic states. LS19 sits early,
# with the progenitor compartment and well before the mature DC/pDC end - supporting that it is a
# least-differentiated state. Data = monocle_pseudotime.csv. No fabricated values.
suppressMessages({library(ggplot2); library(dplyr)})
D3  <- "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/08_Source_Data/Figure_3/v3"
OUT <- "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/04_Main_Figures/Figure_3_panels"
d <- read.csv(file.path(D3,"monocle_pseudotime.csv")) %>% filter(is.finite(monocle_pt))
d$cat <- ifelse(d$group=="LS_19","LS19 (poor)",
         ifelse(d$group=="NORM_HSPC","normal HSPC",
         ifelse(d$group=="NORM_DC","normal DC/pDC",
         ifelse(grepl("^LS_",d$group),"other leukemic", NA))))
d <- d[!is.na(d$cat),]
ordv <- d %>% group_by(cat) %>% summarise(m=median(monocle_pt)) %>% arrange(m) %>% pull(cat)
d$cat <- factor(d$cat, levels=ordv)
pal <- c("normal HSPC"="#E69F00","LS19 (poor)"="#d1495b","other leukemic"="#9aa0a6","normal DC/pDC"="#0072B2")
p <- ggplot(d, aes(cat, monocle_pt, fill=cat)) +
  geom_boxplot(width=0.6, outlier.size=0.5, alpha=0.85) +
  scale_fill_manual(values=pal, guide="none") +
  labs(x=NULL, y="Monocle3 pseudotime (rooted at HSPC)") +
  theme_minimal(base_size=11) +
  theme(panel.grid.major.x=element_blank(), axis.text.x=element_text(angle=25, hjust=1, size=9))
ggsave(file.path(OUT,"Figure_3_LS19_monocle.pdf"), p, width=4.6, height=4.0)
ggsave(file.path(OUT,"Figure_3_LS19_monocle.png"), p, width=4.6, height=4.0, dpi=200)
cat("median monocle pt:", d %>% group_by(cat) %>% summarise(m=round(median(monocle_pt),1)) %>% pull(m) %>% paste(collapse=" / "), "\n")
