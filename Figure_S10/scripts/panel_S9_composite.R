#!/usr/bin/env Rscript
# S9 (new number S6): cellular composition and cell cycle, rebuilt for the corrected scheme.
# A per-sample lineage composition; B per-sample cell-cycle phase; C cell-cycle fraction by
# outcome (n.s.); D cycling fraction by EFS prognosis of the leukemic state (poor trend, n.s.);
# E normal-lineage depletion in AML vs healthy BM. Data = compute_S9_composition.py +
# state_cellcycle_by_prognosis.csv. No fabricated values.
suppressMessages({library(ggplot2); library(dplyr); library(tidyr); library(ggsignif); library(patchwork)})
D2  <- "../source_data"
OUT <- "../figures"
lin  <- read.csv(file.path(D2,"s9_lineage_by_sample.csv"))
ph   <- read.csv(file.path(D2,"s9_phase_by_sample.csv"))
norm <- read.csv(file.path(D2,"s9_normal_lineage_by_sample.csv"))
cc   <- read.csv(file.path(D2,"state_cellcycle_by_prognosis.csv"))
linlev <- c("AML","AML-B","AML-Ery","AML-Mono","AML-NK","AML-T","0_HSPC","Myeloid_Pro",
            "Erythrocytes","Monocyte","NK","T","B","PlasmaB")
lincol <- setNames(c("#b2182b","#d6604d","#f4a582","#e08214","#c51b7d","#8073ac",
                     "#2166ac","#4393c3","#0fcfc0","#92c5de","#66c2a5","#3288bd","#5aae61","#1b7837"),
                   linlev)
lin$lineage <- factor(lin$lineage, levels=linlev)
samp_ord <- sort(unique(lin$SampleID)); lin$SampleID <- factor(lin$SampleID, levels=samp_ord)
pA <- ggplot(lin, aes(SampleID, frac, fill=lineage)) + geom_col(width=0.9) +
  scale_fill_manual(values=lincol, name="lineage") + labs(x=NULL, y="fraction", title="A  Per-sample lineage composition") +
  theme_minimal(base_size=9) + theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1,size=6),
    plot.title=element_text(size=11,face="bold"), legend.key.size=unit(0.32,"cm"), legend.text=element_text(size=6.5))

phl <- ph %>% select(SampleID,G1,S,G2M) %>% pivot_longer(-SampleID, names_to="phase", values_to="frac")
phl$phase <- factor(phl$phase, levels=c("G1","S","G2M")); phl$SampleID <- factor(phl$SampleID, levels=samp_ord)
phcol <- c(G1="#a6cee3", S="#fdbf6f", G2M="#fb9a99")
pB <- ggplot(phl, aes(SampleID, frac, fill=phase)) + geom_col(width=0.9) +
  scale_fill_manual(values=phcol, name="phase") + labs(x=NULL, y="fraction", title="B  Per-sample cell-cycle phase") +
  theme_minimal(base_size=9) + theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1,size=6),
    plot.title=element_text(size=11,face="bold"))

phc <- ph %>% select(SampleID,G1,S,G2M,outcome) %>% pivot_longer(c(G1,S,G2M), names_to="phase", values_to="frac")
phc$phase <- factor(phc$phase, levels=c("G1","S","G2M"))
phc$outcome <- factor(phc$outcome, levels=c("HBM","Alive","Deceased"))
ocol <- c(HBM="#2a9d8f", Alive="#4a6fe3", Deceased="#d1495b")
pC <- ggplot(phc, aes(outcome, frac, fill=outcome)) + geom_boxplot(width=0.6, outlier.size=0.6) +
  facet_wrap(~phase, nrow=1) + scale_fill_manual(values=ocol, guide="none") +
  geom_signif(comparisons=list(c("Alive","Deceased")), test="wilcox.test", textsize=2.6, tip_length=0.01, colour="grey30") +
  labs(x=NULL, y="per-sample fraction", title="C  Cell cycle by outcome") +
  theme_minimal(base_size=9) + theme(axis.text.x=element_text(angle=30,hjust=1,size=7),
    strip.text=element_text(face="bold"), plot.title=element_text(size=11,face="bold"))

cc$prognosis <- factor(cc$prognosis, levels=c("favorable","poor","ns"), labels=c("favorable","poor","n.s."))
pcol <- c(favorable="#00798c", poor="#d1495b", `n.s.`="#9aa0a6")
pD <- ggplot(cc, aes(prognosis, cycling, colour=prognosis, fill=prognosis)) +
  geom_boxplot(width=0.55, alpha=0.25, outlier.shape=NA) + geom_jitter(width=0.12, size=1, alpha=0.8) +
  scale_colour_manual(values=pcol, guide="none") + scale_fill_manual(values=pcol, guide="none") +
  geom_signif(comparisons=list(c("favorable","poor")), test="wilcox.test", textsize=2.6, tip_length=0.01, colour="grey30") +
  labs(x=NULL, y="cycling fraction (S+G2M) per state", title="D  Cycling by state prognosis") +
  theme_minimal(base_size=9) + theme(axis.text.x=element_text(angle=30,hjust=1,size=7), plot.title=element_text(size=11,face="bold"))

norm$sample_type <- factor(norm$sample_type, levels=c("HBM","AML"))
norm$lineage <- factor(norm$lineage, levels=c("0_HSPC","Myeloid_Pro","Erythrocytes","Monocyte","NK","T","B","PlasmaB"))
pE <- ggplot(norm, aes(sample_type, frac, fill=sample_type)) + geom_boxplot(width=0.6, outlier.size=0.5) +
  facet_wrap(~lineage, nrow=2, scales="free_y") + scale_fill_manual(values=c(HBM="#2a9d8f", AML="#b2182b"), name="sample") +
  geom_signif(comparisons=list(c("HBM","AML")), test="wilcox.test", textsize=2.4, tip_length=0.01, colour="grey30") +
  labs(x=NULL, y="per-sample fraction", title="E  Normal-lineage abundance, AML vs healthy BM") +
  theme_minimal(base_size=9) + theme(axis.text.x=element_text(size=7), strip.text=element_text(size=7,face="bold"),
    plot.title=element_text(size=11,face="bold"))

fig <- pA / pB / (pC | pD) / pE + plot_layout(heights=c(1,1,1.1,1.4))
ggsave(file.path(OUT,"Figure_S10__composition_cell_cycle.pdf"), fig, width=11, height=15)
ggsave(file.path(OUT,"Figure_S10__composition_cell_cycle.png"), fig, width=11, height=15, dpi=150)
cat("wrote rebuilt S9->S6 composition+cell-cycle composite\n")
