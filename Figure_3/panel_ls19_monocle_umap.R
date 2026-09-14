#!/usr/bin/env Rscript
# Figure 3C (LS19): Monocle3 UMAP (batch-aligned by patient, single HSPC-rooted trajectory) coloured
# by Monocle pseudotime on the shared manuscript SEQ gradient (the dot-plot colour scheme). LS19 cells
# are outlined to highlight WHERE they sit: in the early / least-differentiated region of the trajectory,
# near normal HSPC and far from mature DC/pDC. Data = monocle_umap_coords.csv (from monocle3_trajectory.R).
# No fabricated values.
suppressMessages({library(ggplot2); library(dplyr); library(ggrepel)})
CODE <- "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/07_Code/Figure_3"
source(file.path(CODE, "manuscript_palette.R"))
D3  <- "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/08_Source_Data/Figure_3/v3"
OUT <- "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/04_Main_Figures/Figure_3_panels"
d <- read.csv(file.path(D3, "monocle_umap_coords.csv"))
d <- d[is.finite(d$monocle_pt), ]
d$is19 <- d$group == "LS_19"
anch <- d %>% filter(group %in% c("NORM_HSPC","NORM_DC")) %>% group_by(group) %>%
        summarise(UMAP1=median(UMAP1), UMAP2=median(UMAP2), .groups="drop")
anch$lab <- c(NORM_HSPC="normal HSPC", NORM_DC="normal DC/pDC")[anch$group]
c19 <- d %>% filter(is19) %>% summarise(UMAP1=median(UMAP1), UMAP2=median(UMAP2))
p <- ggplot() +
  geom_point(data=d[!d$is19,], aes(UMAP1, UMAP2, colour=monocle_pt), size=0.35, alpha=0.55) +
  geom_point(data=d[d$is19,],  aes(UMAP1, UMAP2, fill=monocle_pt), shape=21,
             colour="#d1495b", stroke=0.45, size=1.5) +
  scale_colour_seq(name="Monocle\npseudotime") + scale_fill_seq(guide="none") +
  geom_text_repel(data=anch, aes(UMAP1, UMAP2, label=lab), size=3.1, fontface="italic",
                  colour="grey15", bg.color="white", bg.r=0.2, seed=1, box.padding=0.7,
                  min.segment.length=0, segment.colour="grey40") +
  geom_text_repel(data=c19, aes(UMAP1, UMAP2, label="LS19"), size=4, fontface="bold",
                  colour="#9B2226", bg.color="white", bg.r=0.25, seed=1,
                  nudge_x=6, nudge_y=-4, min.segment.length=0, segment.colour="#9B2226") +
  labs(x="Monocle UMAP1", y="Monocle UMAP2") +
  theme_classic(base_size=11) +
  theme(axis.text=element_blank(), axis.ticks=element_blank(),
        legend.position="right", legend.key.height=unit(0.7,"cm"))
ggsave(file.path(OUT, "Figure_3_LS19_monocle_umap.pdf"), p, width=5.8, height=4.6)
ggsave(file.path(OUT, "Figure_3_LS19_monocle_umap.png"), p, width=5.8, height=4.6, dpi=200)
cat("wrote Figure_3_LS19_monocle_umap | LS19 median pt =",
    round(median(d$monocle_pt[d$is19]),2), "| overall range", round(range(d$monocle_pt),2), "\n")
