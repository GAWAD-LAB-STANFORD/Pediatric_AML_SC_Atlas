# ====================================================================
# fig_adult_scrna.R  |  Figure S21 (adult AML single-cell target validation)
# reproducibility header (auto-generated; edit freely below this block)
# Inputs : adult_scrna_overall.csv, adult_scrna_subtype.csv
# Outputs: (writes to ../figures or ../source_data — see body)
# Run    : Rscript fig_adult_scrna.R   (or: bash run_all.sh)  |  Deps: see REPRODUCIBILITY_CHECKLIST.md
# Source data for this figure is in __SUBMISSION_PACKAGE/08_Source_Data/.
# ====================================================================

# SUPPLEMENT — adult AML single-cell (AML scAtlas, malignant blasts), the 12 panel-5E markers.
#   A  scVI UMAP — compartment + all 12 markers (embedded scanpy render)
#   B  overall % of adult AML blasts positive (single cell)
#   C  by ELN molecular-genetic subtype (NPM1, t(8;21), TP53, complex, −7, MDS-related, …)
source("_style.R")
suppressPackageStartupMessages(library(png))

gcol <- function(gene, group) dplyr::case_when(
  group == "lead" ~ as.character(LEAD_COLORS[gene]),
  group == "control" ~ "#b7950b", group == "clinical" ~ "#e8590c", TRUE ~ PPAC_COLOR)

ov <- sd("adult_scrna_overall.csv"); ov$col <- gcol(ov$gene, ov$group)
GORD <- ov |> dplyr::arrange(dplyr::desc(pct)) |> dplyr::pull(label)
labcol <- setNames(ov$col, ov$label)
gmd <- function(v) sprintf("<b style='color:%s'>%s</b>", labcol[v], v)

## A — embedded scVI UMAP (compartment + 12 genes)
pA <- patchwork::wrap_elements(full = grid::rasterGrob(
  png::readPNG("../assets/adult_scrna_umap.png"), interpolate = TRUE))

## B — overall % positive
ovf <- ov; ovf$label <- factor(ovf$label, levels = rev(GORD))
pB <- ggplot(ovf, aes(pct, label)) +
  geom_segment(aes(x = 0, xend = pct, yend = label, colour = col), linewidth = 0.8) +
  geom_point(aes(colour = col), size = 3) +
  geom_text(aes(label = sprintf("%.0f%%", pct)), hjust = -0.3, size = 2.7) +
  scale_colour_identity() + scale_x_continuous(expand = expansion(c(0, 0.16))) +
  scale_y_discrete(labels = function(v) gmd(v)) +
  labs(x = "% of adult AML blasts positive (single cell)", y = NULL,
       title = "B  Marker positivity in adult AML blasts",
       subtitle = NULL) +
  theme_pub(9) + theme(axis.text.y = ggtext::element_markdown(),
                       plot.subtitle = element_text(size = 7.5, colour = "grey35"))

## C — by ELN molecular-genetic subtype (>= 3 donors), ordered favourable -> adverse
sub <- sd("adult_scrna_subtype.csv")
SUBORD <- c("t(8;21)", "inv(16)/CBF", "NPM1", "NPM1/FLT3-ITD", "FLT3-ITD",
            "MDS-related", "−7/del(7q)", "Complex karyotype", "TP53", "Other")
keep <- sub |> dplyr::distinct(subtype, n_donors) |> dplyr::filter(n_donors >= 3)
SUBORD <- SUBORD[SUBORD %in% keep$subtype]
slab <- setNames(sprintf("%s\n(n=%d)", gsub("−", "-", keep$subtype), keep$n_donors), keep$subtype)
sd2 <- sub |> dplyr::filter(subtype %in% SUBORD)
sd2$gy <- factor(sd2$label, levels = rev(GORD))
sd2$subtype <- factor(sd2$subtype, levels = SUBORD)
pC <- ggplot(sd2, aes(subtype, gy)) +
  geom_point(aes(size = pct, fill = mean), shape = 21, colour = "#3b5870", stroke = 0.2) +
  scale_fill_gradientn(colours = SEQ, name = "mean expr.\n(CP10k log1p)") +
  scale_size_area(max_size = 7, limits = c(0, 100), breaks = c(1, 10, 25, 50, 75), name = "% positive") +
  scale_y_discrete(limits = rev(GORD), labels = function(v) gmd(v)) +
  scale_x_discrete(labels = function(v) slab[v]) +
  labs(x = NULL, y = NULL, title = "C  By ELN molecular-genetic subtype",
       subtitle = NULL) +
  theme_pub(9) + theme(axis.text.y = ggtext::element_markdown(),
                       axis.text.x = element_text(angle = 35, hjust = 1),
                       plot.subtitle = element_text(size = 7.5, colour = "grey35"))

p <- pA / (pB | pC) + plot_layout(heights = c(1.35, 1)) +
  plot_annotation(title = "Adult AML single-cell — 12 targets in malignant blasts (AML scAtlas, adult subset)",
                  theme = theme(plot.title = element_text(face = "bold", size = 12)))
savefig(p, "supplementary/FigureS21_adult_scrna", 15, 13)
cat("wrote supplementary/FigureS21_adult_scrna  (ELN subtypes:", paste(SUBORD, collapse = " "), ")\n")
