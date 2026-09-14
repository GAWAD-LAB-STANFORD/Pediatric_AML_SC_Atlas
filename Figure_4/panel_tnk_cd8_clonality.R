#!/usr/bin/env Rscript
# Figure 4 (T/NK): per-patient cytotoxic-CD8 fraction vs TCR repertoire clonality. Patients whose T/NK
# compartment is skewed toward expanded effector CD8 have more clonally focused repertoires (higher
# clonality = lower Shannon diversity); naive-heavy patients are more diverse. Ties composition to
# clonality at the patient level. Data = tcr_cd8_vs_clonality.csv. No fabricated values.
suppressMessages({library(ggplot2)})
source("/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/07_Code/Figure_3/manuscript_palette.R")
D  <-"/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/08_Source_Data/Figure_4/v2"
OUT<-"/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/04_Main_Figures/Figure_4_panels"
d <- read.csv(D |> file.path("tcr_cd8_vs_clonality.csv"), check.names=FALSE)
d$cyto <- d[["cytotoxic_CD8"]]*100
rs <- cor.test(d$cyto, d$clonality, method="spearman")
g <- ggplot(d, aes(cyto, clonality)) +
  geom_smooth(method="lm", se=TRUE, colour="grey40", fill="grey85", linewidth=.6) +
  geom_point(size=2.6, shape=21, fill=unname(PROG["poor"]), colour="grey20", stroke=.3) +
  annotate("text", x=min(d$cyto), y=max(d$clonality),
           label=sprintf("Spearman r = %.2f, p = %.3f  (n=%d patients)", rs$estimate, rs$p.value, nrow(d)),
           hjust=0, vjust=1, size=3.4) +
  labs(x="cytotoxic CD8 (GZMK + GZMB) % of T/NK, per patient",
       y="TCR clonality (1 - normalized Shannon)") +
  theme_classic(base_size=11)
ggsave(file.path(OUT,"Figure_4_tcr_cd8_vs_clonality.png"), g, width=5.2, height=4.4, dpi=200)
cat("wrote Figure_4_tcr_cd8_vs_clonality | Spearman r=%.2f p=%.3f\n" |> sprintf(rs$estimate, rs$p.value))
