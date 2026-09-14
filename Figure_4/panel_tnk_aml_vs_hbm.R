#!/usr/bin/env Rscript
# Figure 4B: the SAME Harmony batch-corrected T/NK UMAP as panel A, coloured by origin (AML vs healthy BM).
# After batch correction AML and healthy T/NK intermix across the same states, so AML-associated changes
# are shifts in subset proportion (panel D) rather than new AML-specific T/NK states. The un-integrated
# "do they co-cluster without forcing" view is a supplementary figure. Data = tnk_subset_percell.csv.
suppressMessages({library(ggplot2)})
source("/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/07_Code/Figure_3/manuscript_palette.R")
D   <- "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/08_Source_Data/Figure_4/v2"
OUT <- "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/04_Main_Figures/Figure_4_panels"
d <- read.csv(file.path(D, "tnk_subset_percell.csv"))   # integrated (Harmony) coords, same as panel A
d$group <- factor(d$group, levels = c("HBM", "AML"))
d <- d[sample(nrow(d)), ]                        # shuffle so neither colour is fully on top
lab <- c(AML = sprintf("AML (n=%d)", sum(d$group=="AML")), HBM = sprintf("healthy BM (n=%d)", sum(d$group=="HBM")))
g <- ggplot(d, aes(UMAP1, UMAP2, colour = group)) +
  geom_point(size = 0.35, alpha = 0.55) +
  scale_colour_manual(values = c(HBM = unname(PROG["favorable"]), AML = unname(PROG["poor"])),
                      labels = lab, name = NULL) +
  guides(colour = guide_legend(override.aes = list(size = 2.5, alpha = 1))) +
  labs(x = "UMAP1", y = "UMAP2") +
  theme_classic(base_size = 11) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        legend.position = c(0.02, 0.02), legend.justification = c(0, 0))
ggsave(file.path(OUT, "Figure_4D_tnk_aml_vs_hbm.pdf"), g, width = 5.2, height = 4.8)
ggsave(file.path(OUT, "Figure_4D_tnk_aml_vs_hbm.png"), g, width = 5.2, height = 4.8, dpi = 200)
cat("wrote Figure_4D_tnk_aml_vs_hbm\n")
