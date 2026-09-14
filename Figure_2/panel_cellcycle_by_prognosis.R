#!/usr/bin/env Rscript
# S9 Panel D (regenerated for the corrected scheme): per-state cell-cycle phase fractions grouped
# by whole-cohort EFS prognosis (favorable / poor / n.s.), replacing the stale FPAC/PPAC panel.
# Unit = leukemic state (each state contributes one phase fraction); groups compared by two-sided
# Mann-Whitney. Poor-prognosis states trend more proliferative but not significantly at the state
# level, consistent with Panel C (cell-cycle vs outcome, n.s.) and with proliferation being driven
# by SPECIFIC states (LS_10, LS_20), not group-level cycling. Data = compute above. No fabrication.
suppressMessages({library(ggplot2); library(tidyr); library(dplyr); library(ggsignif)})
D2  <- "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/08_Source_Data/Figure_2/v3"
OUT <- "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/05_Supplementary_Figures"
d <- read.csv(file.path(D2,"state_cellcycle_by_prognosis.csv"))
pcol <- c(favorable="#00798c", poor="#d1495b", ns="#9aa0a6")
d$prognosis <- factor(d$prognosis, levels=c("favorable","poor","ns"), labels=c("favorable","poor","n.s."))
long <- pivot_longer(d, cols=c("G1","S","G2M"), names_to="phase", values_to="frac")
long$phase <- factor(long$phase, levels=c("G1","S","G2M"))
names(pcol) <- c("favorable","poor","n.s.")
p <- ggplot(long, aes(prognosis, frac, colour=prognosis, fill=prognosis)) +
  geom_boxplot(width=0.55, alpha=0.25, outlier.shape=NA) +
  geom_jitter(width=0.12, size=1.1, alpha=0.8) +
  facet_wrap(~phase, nrow=1) +
  scale_colour_manual(values=pcol, guide="none") + scale_fill_manual(values=pcol, guide="none") +
  geom_signif(comparisons=list(c("favorable","poor"), c("poor","n.s.")),
              map_signif_level=FALSE, test="wilcox.test", colour="grey30",
              textsize=2.8, step_increase=0.10, tip_length=0.01) +
  labs(x=NULL, y="per-state fraction of cells in phase",
       title="Cell-cycle phase by whole-cohort EFS prognosis (per leukemic state)") +
  theme_minimal(base_size=11) +
  theme(panel.grid.minor=element_blank(), strip.text=element_text(face="bold"),
        plot.title=element_text(size=11, face="bold"))
ggsave(file.path(OUT,"Figure_S9_panelD_cellcycle_by_prognosis.pdf"), p, width=7.5, height=3.8)
ggsave(file.path(OUT,"Figure_S9_panelD_cellcycle_by_prognosis.png"), p, width=7.5, height=3.8, dpi=200)
cat("wrote regenerated S9 Panel D (cell cycle by prognosis)\n")
