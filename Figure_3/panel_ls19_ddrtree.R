#!/usr/bin/env Rscript
# Figure 3C (LS19): Monocle2 DDRTree trajectory (the classic branched-"line" reverse-graph embedding).
# Cells sit on the tree skeleton, coloured by DDRTree pseudotime on the shared SEQ gradient (dot-plot
# scheme). Direction arrows run low->high pseudotime along the main branches; LS19 cells are highlighted
# (small red rings) and the normal-HSPC stem-cell root and mature DC/pDC tip are labelled. Shows LS19
# occupies the early / least-differentiated part of the tree near normal HSPC. Data = monocle2_ddrtree.R.
# No fabricated values.
suppressMessages({library(ggplot2); library(dplyr); library(ggrepel)})
CODE <- "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/07_Code/Figure_3"
source(file.path(CODE, "manuscript_palette.R"))
D3  <- "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/08_Source_Data/Figure_3/v3"
OUT <- "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/04_Main_Figures/Figure_3_panels"
cells <- read.csv(file.path(D3, "ddrtree_cells.csv"))
nodes <- read.csv(file.path(D3, "ddrtree_nodes.csv"))
edges <- read.csv(file.path(D3, "ddrtree_edges.csv"))
cells$is19 <- cells$group == "LS_19"
# per-node pseudotime = median of nearest cells (for orienting the direction arrows)
nb <- as.matrix(nodes[, c("K1","K2")]); cb <- as.matrix(cells[, c("DDR1","DDR2")])
nn <- apply(cb, 1, function(p) which.min((nb[,1]-p[1])^2 + (nb[,2]-p[2])^2))
cells$nn <- nodes$node[nn]
nodes$pt <- as.numeric(tapply(cells$Pseudotime, cells$nn, median)[nodes$node])
# edge geometry, oriented low -> high pseudotime
e <- merge(edges, setNames(nodes[,c("node","K1","K2","pt")], c("from","fx","fy","fpt")), by="from")
e <- merge(e,     setNames(nodes[,c("node","K1","K2","pt")], c("to","tx","ty","tpt")),   by="to")
e <- e[is.finite(e$fpt) & is.finite(e$tpt), ]
lo_from <- e$fpt <= e$tpt
e$sx <- ifelse(lo_from, e$fx, e$tx); e$sy <- ifelse(lo_from, e$fy, e$ty)
e$ex <- ifelse(lo_from, e$tx, e$fx); e$ey <- ifelse(lo_from, e$ty, e$fy)
e$dpt <- abs(e$tpt - e$fpt)
# downsample arrows to the strongest-gradient edges so the flow is legible, not cluttered
set.seed(1); arr <- e[order(-e$dpt), ][seq_len(min(45, nrow(e))), ]
# label each branch by the normal compartment that defines it, placed toward the branch tip (median of
# that compartment's most-differentiated cells). The degree-1 leaf nodes are small leukemic sub-tips, so
# labelling by dominant normal compartment (the differentiation reference points) is what orients the tree.
pretty <- c(NORM_HSPC="HSPC (root)", NORM_MyeloidPro="myeloid prog", NORM_Monocyte="monocyte",
            NORM_Macrophage="macrophage", NORM_DC="DC / pDC", NORM_Erythroid="erythroid",
            NORM_B="B lineage", NORM_T_NK="T / NK")
lf <- cells %>% filter(group %in% names(pretty)) %>% group_by(group) %>%
  summarise(K1 = median(DDR1[Pseudotime >= quantile(Pseudotime, 0.6)]),
            K2 = median(DDR2[Pseudotime >= quantile(Pseudotime, 0.6)]),
            pt = median(Pseudotime), nc = n(), .groups = "drop") %>%
  filter(nc >= 20)
lf$lab <- unname(pretty[lf$group])
cx <- mean(nodes$K1); cy <- mean(nodes$K2); rx <- diff(range(cells$DDR1)); ry <- diff(range(cells$DDR2))
ang <- atan2(lf$K2 - cy, lf$K1 - cx); lf$nx <- cos(ang)*rx*0.10; lf$ny <- sin(ang)*ry*0.10
# anchor the LS19 label at the centroid of the root-proximal MAJORITY of LS19 cells (the main-trunk cells
# near the HSPC root; DDR2 > -1 excludes the minority that trail toward the DC/pDC tip). LS19 is the
# least-differentiated prognostic state, so the leader must point to where most LS19 cells actually sit --
# the early, HSPC-proximal part of the tree -- not the differentiated branch. ~55% of LS19 cells are here.
c19d <- cells[cells$is19 & cells$DDR2 > -1, ]
c19 <- data.frame(DDR1 = median(c19d$DDR1), DDR2 = median(c19d$DDR2))
cat("branch (compartment) labels:\n"); print(lf[, c("lab","pt","nc")], row.names=FALSE)
p <- ggplot() +
  geom_segment(data=e, aes(sx, sy, xend=ex, yend=ey), colour="grey70", linewidth=0.25) +
  geom_segment(data=arr, aes(sx, sy, xend=ex, yend=ey), colour="grey35", linewidth=0.3,
               arrow=arrow(length=unit(0.05,"inches"), type="closed")) +
  geom_point(data=cells[!cells$is19,], aes(DDR1, DDR2, colour=Pseudotime), size=0.45, alpha=0.7) +
  geom_point(data=cells[cells$is19,],  aes(DDR1, DDR2, fill=Pseudotime), shape=21,
             colour="#AE2012", stroke=0.4, size=0.9) +
  scale_colour_seq(name="DDRTree\npseudotime") + scale_fill_seq(guide="none") +
  geom_text_repel(data=lf, aes(K1, K2, label=lab), size=2.9, fontface="italic", colour="grey15",
                  bg.color="white", bg.r=0.18, nudge_x=lf$nx, nudge_y=lf$ny,
                  segment.colour="grey55", min.segment.length=0, seed=1, box.padding=0.4, max.overlaps=Inf) +
  geom_text_repel(data=c19, aes(DDR1, DDR2, label="LS19"), size=4, fontface="bold", colour="#AE2012",
                  bg.color="white", bg.r=0.25, nudge_x=-rx*0.05, nudge_y=ry*0.24,
                  min.segment.length=0, segment.colour="#AE2012", seed=1) +
  labs(x="DDRTree component 1", y="DDRTree component 2") +
  theme_classic(base_size=11) +
  theme(axis.text=element_blank(), axis.ticks=element_blank(),
        legend.position="right", legend.key.height=unit(0.7,"cm"))
ggsave(file.path(OUT, "Figure_3_LS19_ddrtree.pdf"), p, width=5.8, height=4.6)
ggsave(file.path(OUT, "Figure_3_LS19_ddrtree.png"), p, width=5.8, height=4.6, dpi=200)
cat("wrote Figure_3_LS19_ddrtree | nodes", nrow(nodes), "edges", nrow(e),
    "| LS19 median pt", round(median(cells$Pseudotime[cells$is19]),2), "\n")
