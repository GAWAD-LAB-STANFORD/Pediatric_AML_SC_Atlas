#!/usr/bin/env Rscript
# Figure 2C (v2): 14 prognostic states x {cytogenetic | molecular} subtypes, sharing
# the y-axis + one stemness strip. Harmonized ltc palette (heatmap0).
suppressMessages({library(ggplot2); library(patchwork); library(scales)})
SL <- "/private/tmp/claude-503/-Users-chuckgawad-Desktop-ALSF-AML-2026-updated/f5e93d65-31d5-411e-b754-6f7c504928b7/scratchpad"
source(file.path(SL, "manuscript_palette.R"))
OUT <- "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/04_Main_Figures/Figure_2_new_panels"
pr  <- read.csv(file.path(SL, "target_ls/LS48_prognosis.csv")); pr <- pr[pr$fdr < 0.10, ]
ann <- read.csv(file.path(SL, "ls48/LS48_annotation.csv"))[, c("LS", "lineage_program")]
stem<- read.csv(file.path(SL, "ls48/LS48_stem_scores.csv"))[, c("LS", "stem_z")]
pr <- merge(merge(pr, ann, by = "LS"), stem, by = "LS"); pr$LSc <- gsub("_", "", pr$LS)
pr$prog <- ifelse(pr$lineage_program == "unresolved", "", pr$lineage_program)
pr$lab  <- ifelse(pr$prog == "", pr$LSc, paste0(pr$LSc, "  ", pr$prog))
pr <- pr[order(pr$dir == "poor", pr$p), ]; pr$lab <- make.unique(pr$lab)
lv <- pr$lab; nfav <- sum(pr$dir == "favorable"); yline <- nrow(pr) - nfav + 0.5

cyto_lv <- c("PML-RARA","t(8;21)","inv(16)","Normal","Other","KMT2A(MLL)","NUP98-r","CBFA2T3-GLIS2","t(6;9)","mono7/del7q","complex(>=3)")
mol_lv  <- c("FLT3-ITD","NPM1-mut","CEBPA-mut","WT1-mut","Triple-neg")
mkdat <- function(csv, xlv){
  d <- read.csv(file.path(SL, csv)); d <- d[d$LS %in% pr$LS, ]; d <- merge(d, pr[, c("LS","lab")], by="LS")
  d$lab <- factor(d$lab, levels = rev(lv)); d$subtype <- factor(d$subtype, levels = xlv); d }
dc <- mkdat("ls48/LS48_cyto_enrichment.csv", cyto_lv)
dm <- mkdat("ls48/LS48_mol_enrichment.csv",  mol_lv)

dotlayer <- function(d, ylabs){
  ggplot(d, aes(subtype, lab, size = mean_pct, fill = neglog10fdr)) +
    geom_hline(yintercept = yline, colour = "grey40", linewidth = 0.5) +
    geom_point(shape = 21, colour = "grey30", stroke = 0.3) +
    scale_fill_seq(limits = c(0,10), oob = squish, name = expression(-log[10]~FDR)) +
    scale_size(range = c(0.5,9), limits = c(0,27), name = "mean % (deconv.)") +
    labs(x = NULL, y = NULL) + theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
          axis.text.y = if (ylabs) element_text(size = 9) else element_blank(),
          panel.grid.major = element_line(colour = "grey93")) }
pc <- dotlayer(dc, TRUE) + ggtitle("Cytogenetic")
blk <- data.frame(y = c(nrow(pr)-nfav/2+0.5, (nrow(pr)-nfav)/2+0.5), lab = c("Favorable","Poor"))
pblk <- ggplot(blk, aes(1, y, label = lab, colour = lab)) +
  geom_text(angle = 90, fontface = "bold", size = 4.6) +
  scale_colour_manual(values = c(Favorable = unname(PROG["favorable"]), Poor = unname(PROG["poor"])), guide = "none") +
  scale_y_continuous(limits = c(0.5, nrow(pr)+0.5), expand = c(0,0)) +
  scale_x_continuous(limits = c(0.5, 1.5)) + theme_void()
pm <- dotlayer(dm, FALSE) + ggtitle("Molecular")
prs <- pr; prs$lab <- factor(prs$lab, levels = rev(lv))
ps <- ggplot(prs, aes("HSC/LSC", lab, fill = stem_z)) +
  geom_hline(yintercept = yline, colour = "grey40", linewidth = 0.5) +
  geom_tile(colour = "white", linewidth = 0.4) +
  scale_fill_div(midpoint=0, limits=c(-1,2), oob=squish, name="stemness\n(z)") +
  labs(x = NULL, y = NULL) + theme_minimal(base_size = 10) +
  theme(axis.text.y = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        panel.grid = element_blank())
fig <- pblk + pc + pm + ps + plot_layout(widths = c(0.05, 1, 0.5, 0.09), guides = "collect") +
  plot_annotation(theme = theme(plot.margin = margin(12,8,6,8)))
ggsave(file.path(OUT, "Figure_2C_states_subtypes.pdf"), fig, width = 13, height = 6.0)
ggsave(file.path(OUT, "Figure_2C_states_subtypes.png"), fig, width = 13, height = 6.0, dpi = 200)
cat("wrote combined Figure 2C\n")
