# ====================================================================
# fig6_DE.R  |  Figure 6D (12-target radar profiles)
# reproducibility header (auto-generated; edit freely below this block)
# Inputs : fig6D_radar.csv
# Outputs: .panel_fig6D.pdf, _review_fig6D.png
# Run    : Rscript fig6_DE.R   (or: bash run_all.sh)  |  Deps: see REPRODUCIBILITY_CHECKLIST.md
# Source data for this figure is in __SUBMISSION_PACKAGE/08_Source_Data/.
# ====================================================================

# Figure 6 panel D — radar/spider profiles for ALL 12 targets (the panel A/B/C set).
# Each axis is scaled to its own max ACROSS the 12 targets (so the strongest marker reaches the
# edge); the absolute max is printed on the axis. Values are correct single-cell metrics — % of
# leukemic cells positive is dropout-limited, hence the per-axis scaling. (Panel E -> Figure 7.)
source("_style.R")
suppressPackageStartupMessages(library(patchwork))

# corrected scheme: CD7 (a clinical CAR-T antigen) is coloured with the clinical comparators; the
# remaining former-PPAC markers become neutral-slate "additional candidates".
gcol <- function(gene, group) dplyr::case_when(
  group == "lead" ~ as.character(LEAD_COLORS[gene]),
  group == "control" ~ "#b7950b",
  group == "clinical" | gene == "CD7" ~ "#e8590c", TRUE ~ "#7f8c8d")

rad <- sd("fig6D_radar.csv")
AX  <- c("AML_coverage", "Patient_breadth", "Marrow_sparing", "Immune_sparing", "Organ_sparing")
mx  <- sapply(AX, function(a) max(rad[[a]], na.rm = TRUE))
AXL <- c(sprintf("AML cov.\n(≤%.0f%%)", mx["AML_coverage"]),
         sprintf("patient\nbreadth (≤%.0f%%)", mx["Patient_breadth"]),
         sprintf("marrow\nsparing (≤%.0f)", mx["Marrow_sparing"]),
         sprintf("immune\nsparing (≤%.0f)", mx["Immune_sparing"]),
         sprintf("organ\nsparing (≤%.0f)", mx["Organ_sparing"]))
K <- length(AX); ang <- 2 * pi * (0:(K - 1)) / K
xy <- function(v, a) list(x = v * sin(a), y = v * cos(a))
rad$col <- gcol(rad$gene, rad$group)
# per-axis scaled (0..100 of that axis's max)
norm <- as.matrix(sweep(rad[, AX], 2, mx, "/") * 100)
# radar polygon AREA as % of the maximal pentagon (overall-profile score; the sin term cancels)
rad$area <- round(rowSums(norm * norm[, c(2:K, 1)]) / (K * 100), 1)
ord <- order(-rad$area); rad <- rad[ord, ]; norm <- norm[ord, , drop = FALSE]   # best area first
rad$label <- factor(rad$label, levels = as.character(rad$label))

# closed polygon per target
poly <- do.call(rbind, lapply(seq_len(nrow(rad)), function(i) {
  v <- norm[i, ]; p <- xy(c(v, v[1]), c(ang, ang[1]))
  data.frame(label = rad$label[i], col = rad$col[i], x = p$x, y = p$y)
}))
rings <- do.call(rbind, lapply(c(25, 50, 75, 100), function(r) {
  a <- seq(0, 2 * pi, length.out = 90); data.frame(r = r, x = r * sin(a), y = r * cos(a))
}))
spokes <- data.frame(x = 0, y = 0, xe = 100 * sin(ang), ye = 100 * cos(ang))
axtext <- data.frame(x = 128 * sin(ang), y = 128 * cos(ang), lab = AXL)
sl <- setNames(sprintf("<b style='color:%s'>%s</b> · area %.0f", rad$col, rad$label, rad$area),
               as.character(rad$label))

pD <- ggplot() +
  geom_path(data = rings, aes(x, y, group = r), colour = "grey86", linewidth = 0.3) +
  geom_segment(data = spokes, aes(x, y, xend = xe, yend = ye), colour = "grey86", linewidth = 0.3) +
  geom_text(data = axtext, aes(x, y, label = lab), size = 1.85, colour = "grey45", lineheight = 0.78) +
  geom_polygon(data = poly, aes(x, y, fill = col, colour = col), alpha = 0.33, linewidth = 0.6) +
  geom_point(data = poly, aes(x, y, colour = col), size = 0.7) +
  scale_fill_identity() + scale_colour_identity() +
  facet_wrap(~label, ncol = 4, labeller = as_labeller(sl)) +
  coord_equal(clip = "off") + xlim(-145, 145) + ylim(-140, 140) +
  labs(title = "D") +                                # commentary -> caption
  theme_void(base_size = 9) +
  theme(plot.title = element_text(face = "bold", size = 11),
        strip.text = ggtext::element_markdown(face = "bold", size = 8.3),
        panel.spacing = grid::unit(0.2, "lines"), plot.margin = margin(6, 10, 6, 10))
ggsave("../figures/_review_fig6D.png", pD, width = 13, height = 10, dpi = 200, bg = "white",
       device = ragg::agg_png)
ggsave("../figures/.panel_fig6D.pdf", pD, width = 13, height = 10, bg = "white",
       device = .PDF_DEVICE)                        # vector panel for the fitz assembler
cat("wrote _review_fig6D.png + .panel_fig6D.pdf  | radar area ranking:\n")
print(rad[, c("label", "area")], row.names = FALSE)
