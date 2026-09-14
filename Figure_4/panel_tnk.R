#!/usr/bin/env Rscript
# Figure 4 rebuild (corrected scheme): T/NK landscape re-derived from the count matrix.
# A UMAP of the 12 T/NK subsets; B per-sample composition (AML vs healthy BM, HBM n=2 -> descriptive,
# no formal test); C defining-marker dot plot. Shared manuscript palette: categorical subsets sampled
# from the SEQ gradient (lineage-ordered), dot-plot fill = scale_fill_seq. No fabricated values.
suppressMessages({library(ggplot2); library(dplyr); library(scales); library(ggrepel)})
CODE <- "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/07_Code/Figure_4"
source("/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/07_Code/Figure_3/manuscript_palette.R")
D   <- "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/08_Source_Data/Figure_4/v2"
OUT <- "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/04_Main_Figures/Figure_4_panels"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
ORD <- c("Naïve CD4 T","Memory CD4 T","Treg","Naïve CD8 T","Memory CD8 T","GZMK CD8 T",
         "GZMB CD8 T","GZMB DNT","MAIT","GZMK NK","GZMB NK","Proliferating T")
# 12 distinct colours for the T/NK subsets (lineage-grouped, colourblind-aware; a qualitative palette is
# used here because 12 nominal subsets are not separable on a single sequential gradient)
pal <- c("Naïve CD4 T"="#4E79A7","Memory CD4 T"="#A0CBE8","Treg"="#9C755F",
         "Naïve CD8 T"="#59A14F","Memory CD8 T"="#8CD17D","GZMK CD8 T"="#B6992D",
         "GZMB CD8 T"="#F28E2B","GZMB DNT"="#FFBE7D","MAIT"="#499894",
         "GZMK NK"="#E15759","GZMB NK"="#B07AA1","Proliferating T"="#79706E")[ORD]

## ---- A: T/NK UMAP (Harmony batch-corrected embedding, shared with panel B) ----
pc <- read.csv(file.path(D, "tnk_subset_percell.csv"))   # integrated (Harmony/SampleID) UMAP coords
pc$subset <- factor(pc$subset, levels = ORD)
cent <- pc %>% group_by(subset) %>% summarise(UMAP1 = median(UMAP1), UMAP2 = median(UMAP2), .groups = "drop")
gA <- ggplot(pc, aes(UMAP1, UMAP2, colour = subset)) +
  geom_point(size = 0.35, alpha = 0.7) +
  scale_colour_manual(values = pal, name = NULL) +
  geom_text_repel(data = cent, aes(label = subset), size = 3, colour = "grey15",
                  bg.color = "white", bg.r = 0.15, seed = 1, max.overlaps = Inf, show.legend = FALSE) +
  guides(colour = guide_legend(override.aes = list(size = 2.5))) +
  labs(x = "UMAP1", y = "UMAP2") +
  theme_classic(base_size = 11) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), legend.position = "none")
ggsave(file.path(OUT, "Figure_4A_tnk_umap.pdf"), gA, width = 5.4, height = 4.8)
ggsave(file.path(OUT, "Figure_4A_tnk_umap.png"), gA, width = 5.4, height = 4.8, dpi = 200)

## ---- B: per-sample composition (AML vs HBM) ----
comp <- read.csv(file.path(D, "tnk_composition_bySample.csv"))
comp$subset <- factor(comp$subset, levels = ORD)
ord_s <- comp %>% filter(subset == "GZMB NK") %>% arrange(frac) %>% pull(SampleID)
ord_s <- unique(c(ord_s, setdiff(comp$SampleID, ord_s)))
comp$SampleID <- factor(comp$SampleID, levels = ord_s)
gB <- ggplot(comp, aes(SampleID, frac, fill = subset)) +
  geom_col(width = 0.9) +
  facet_grid(~group, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = pal, name = NULL) +
  scale_y_continuous(labels = percent, expand = c(0, 0)) +
  labs(x = NULL, y = "T/NK composition") +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5),
        panel.grid = element_blank(), legend.key.size = unit(0.32, "cm"),
        strip.text = element_text(face = "bold"))
ggsave(file.path(OUT, "Figure_4B_tnk_composition.pdf"), gB, width = 8.2, height = 4.0)
ggsave(file.path(OUT, "Figure_4B_tnk_composition.png"), gB, width = 8.2, height = 4.0, dpi = 200)

## ---- C: defining-marker dot plot ----
dp <- read.csv(file.path(D, "tnk_subset_dotplot.csv"))
dp$subset <- factor(dp$subset, levels = rev(ORD))
dp$gene <- factor(dp$gene, levels = unique(dp$gene))
gC <- ggplot(dp, aes(gene, subset, size = pct, fill = mean)) +
  geom_point(shape = 21, colour = "grey30", stroke = 0.3) +
  scale_size(range = c(0.2, 6), name = "% expressing", labels = percent) +
  scale_fill_seq(name = "mean\nlog-norm") +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8, face = "italic"),
        panel.grid.major = element_line(colour = "grey92"), legend.key.size = unit(0.35, "cm"))
ggsave(file.path(OUT, "Figure_4C_tnk_markers.pdf"), gC, width = 8.4, height = 4.2)
ggsave(file.path(OUT, "Figure_4C_tnk_markers.png"), gC, width = 8.4, height = 4.2, dpi = 200)
cat("wrote Figure_4 panels A/B/C to", OUT, "\n")
