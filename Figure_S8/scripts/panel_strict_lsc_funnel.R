#!/usr/bin/env Rscript
# Strict true-LSC funnel: sequential gating of AML-sample vs healthy-BM cells toward the
# literature LSC definition (CD34+CD38-, HLF/AVP+ quiescent, non-cycling, and carrying an
# LSC-aberrant marker TIM3/CD96/CLEC12A that normal HSC lack). Counts are the observed
# gate results (see compute log); log-x. Shows true LSC are vanishingly rare (26 AML vs 2 BM).
suppressMessages({library(ggplot2); library(scales)})
CODE <- "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/07_Code/Figure_3"
source(file.path(CODE, "manuscript_palette.R"))
OUT <- "../figures/Figure_2_panels"
steps <- c("CD34+","+ CD38-","+ quiescent (HLF/AVP+)","+ non-cycling","+ LSC marker (TIM3/CD96/CLEC12A)")
d <- data.frame(
  step  = factor(rep(steps, 2), levels = rev(steps)),
  group = rep(c("AML sample","healthy BM"), each = 5),
  count = c(28304, 13194, 95, 78, 26,   812, 518, 160, 102, 2))
g <- ggplot(d, aes(count, step, fill = group)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.65) +
  geom_text(aes(label = count), position = position_dodge(width = 0.7), hjust = -0.15, size = 3) +
  scale_x_log10(expand = expansion(mult = c(0, 0.18)), labels = comma) +
  scale_fill_manual(values = c("AML sample" = unname(PROG["poor"]), "healthy BM" = "#0A9396"), name = NULL) +
  labs(x = "cells (log scale)", y = NULL,
       caption = "final gate = strict true-LSC: 26 AML-sample cells (2 patients, both FLT3-ITD) vs 2 healthy-BM.") +
  theme_classic(base_size = 10) +
  theme(legend.position = "top", plot.caption = element_text(size = 7.6, hjust = 0),
        axis.text.y = element_text(size = 9))
ggsave(file.path(OUT, "Figure_3_strict_LSC_funnel.pdf"), g, width = 6.0, height = 3.4)
ggsave(file.path(OUT, "Figure_3_strict_LSC_funnel.png"), g, width = 6.0, height = 3.4, dpi = 200)
cat("wrote Figure_3_strict_LSC_funnel\n")
