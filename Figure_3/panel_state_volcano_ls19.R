#!/usr/bin/env Rscript
# Selection volcano: every leukemic state tested for whole-cohort EFS (per-SD Cox).
# x = risk (log2 hazard ratio; right = poor, left = favorable), y = -log10 Cox p.
# Dashed line = the p that corresponds to BH-FDR<0.10 (the selection cutoff). LS10 is the
# top poor state; all FDR<0.10 states labelled. Shows LS10 was selected by testing all LS.
suppressMessages({library(ggplot2); library(ggrepel)})
CODE <- "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/07_Code/Figure_3"
source(file.path(CODE, "manuscript_palette.R"))
D2 <- "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/08_Source_Data/Figure_2/v3"
OUT <- "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/04_Main_Figures/Figure_2_panels"
t <- read.csv(file.path(D2, "KM_all_states_table.csv"))
d <- t[t$class == "leukemic", ]
d$log2HR <- log2(d$WC_coxHR); d$neglogp <- -log10(d$WC_coxP)
d$sig <- ifelse(d$WC_coxFDR < 0.10, ifelse(d$WC_coxHR > 1, "poor", "favorable"), "n.s.")
d$sig <- factor(d$sig, levels = c("favorable","poor","n.s."))
d$LSc <- gsub("_","",d$state)
fdr_line <- -log10(max(d$WC_coxP[d$WC_coxFDR < 0.10]))     # p that maps to FDR = 0.10
cols <- c(favorable = unname(PROG["favorable"]), poor = unname(PROG["poor"]), `n.s.` = GREY)
lab <- d[d$sig != "n.s." | d$state == "LS_19", ]
g <- ggplot(d, aes(log2HR, neglogp)) +
  geom_hline(yintercept = fdr_line, linetype = 2, colour = "grey45") +
  geom_vline(xintercept = 0, linetype = 3, colour = "grey60") +
  geom_point(aes(colour = sig, size = state == "LS_19")) +
  geom_point(data = d[d$state=="LS_19",], shape = 21, size = 5, stroke = 1.1, colour = "black", fill = NA) +
  ggrepel::geom_text_repel(data = lab, aes(label = LSc, colour = sig),
      size = 3, fontface = ifelse(lab$state=="LS_19","bold","plain"), max.overlaps = 30, seed = 1, show.legend = FALSE) +
  scale_colour_manual(values = cols, name = NULL) +
  scale_size_manual(values = c(`FALSE` = 2, `TRUE` = 4), guide = "none") +
  annotate("text", x = min(d$log2HR), y = fdr_line + 0.3, label = "FDR = 0.10", hjust = 0, size = 3, colour = "grey35") +
  annotate("text", x = max(d$log2HR)*0.9, y = 0.2, label = "poor →", hjust = 1, size = 3.2, colour = unname(PROG["poor"])) +
  annotate("text", x = min(d$log2HR)*0.9, y = 0.2, label = "← favorable", hjust = 0, size = 3.2, colour = unname(PROG["favorable"])) +
  labs(x = expression(risk:~log[2]~hazard~ratio~(per~SD)), y = expression(-log[10]~Cox~p)) +
  theme_classic(base_size = 11) + theme(legend.position = "top")
ggsave(file.path(OUT, "Figure_3_LS19_volcano.pdf"), g, width = 6.4, height = 5.2)
ggsave(file.path(OUT, "Figure_3_LS19_volcano.png"), g, width = 6.4, height = 5.2, dpi = 200)
cat(sprintf("wrote volcano; FDR=0.10 line at -log10 p = %.2f; LS10 at log2HR=%.2f, -log10p=%.1f\n",
    fdr_line, log2(d$WC_coxHR[d$state=="LS_19"]), -log10(d$WC_coxP[d$state=="LS_19"])))
