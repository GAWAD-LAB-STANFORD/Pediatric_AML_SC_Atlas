# ====================================================================
# _style.R  |  Shared style + helpers (colours, sd() data loader) for all CD96 R figures
# reproducibility header (auto-generated; edit freely below this block)
# Inputs : (reads upstream object / raw matrix — see body)
# Outputs: (writes to ../figures or ../source_data — see body)
# Run    : Rscript _style.R   (or: bash run_all.sh)  |  Deps: see REPRODUCIBILITY_CHECKLIST.md
# Source data for this figure is in __SUBMISSION_PACKAGE/08_Source_Data/.
# ====================================================================

# Shared style + helpers for the CD96 R/ggplot2 figure package.
# Run scripts from this scripts/ directory; data read from ../source_data,
# figures written to ../figures.
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(readr)
  library(forcats); library(scales); library(ggrepel); library(patchwork)
  library(ggnewscale); library(ggtext)
})

SRC <- "../source_data"; FIG <- "../figures"
sd <- function(f) readr::read_csv(file.path(SRC, f), show_col_types = FALSE)

# CD96 = crimson everywhere; navy = comparators/magnitude; green = safety; grey = other
PAL <- list(CD96 = "#c0392b", navy = "#2471a3", green = "#1e8449",
            grey = "#95a5a6", sparing = "#1e8449", toxic = "#95a5a6")
# 9-target categorical colours -- one fixed colour per gene, used consistently in
# every figure + supplement (dot-plot gene labels, Fig 7C/D points & bars, violins).
# Palette = Okabe-Ito + Paul Tol (the Nature-recommended colourblind-safe families),
# CD96 = crimson hero; the Okabe yellow was swapped for purple (yellow is invisible
# on white).  Scored best of the candidates for min pairwise distance under
# deuteranopia/protanopia.  Keyed by GLABEL.
TARGET_COLORS <- c("CD96" = "#c0392b", "CD33" = "#E69F00", "CD123" = "#56B4E9",
                   "CLL-1" = "#009E73", "FLT3" = "#7B3294", "MSLN" = "#0072B2",
                   "CD70" = "#CC79A7", "GPR56" = "#999933", "IL1RAP" = "#882255")
GENES <- names(TARGET_COLORS)                     # display order
# --- ALTERNATE figure set: the top-5 LEAD targets from the relaxed-gate (10%/10%)
# Figure-5A composite (CD96 + the 4 additional candidates). CD96 stays crimson; the
# 4 additional get distinct Okabe-Ito (colourblind-safe) colours. Used as the hero
# target set across the alternate main + supplemental figures.
LEAD_TARGETS <- c("CD96", "CD9", "SUCNR1", "IL1RAP")   # de-novo candidate panel (CD96 top-ranked; CD9/IL1RAP later disqualified on organ toxicity)
# PPAC-specific markers (from Figure-5 supplement S17): surface genes that mark individual PPAC
# subsets rather than pan-AML leukemia. Shown as a distinct group (purple) in panels C and D.
# NB CD7 is also a clinical CAR-T antigen; here it is coloured with this PPAC-specific group.
# corrected scheme: these are now "additional candidate" antigens (CD7 is folded into the clinical
# comparators elsewhere; the group is coloured neutral slate rather than the retired PPAC purple).
PPAC_SPECIFIC <- c("CD7", "TNFRSF4", "ABCA7", "ITGAX")
PPAC_COLOR    <- "#7f8c8d"
# MPAL (mixed-phenotype acute leukemia) samples in the scRNA cohort map 1:1 to these three
# cytogenetic subtypes (AML4363 = BCR/ABL, AML3121 = t(2;3), AML882 = t(7;14)); they are NOT
# pure AML, so they are flagged in a distinct teal wherever cytogenetic subtypes appear.
MPAL_CYTO  <- c("BCR/ABL", "t(2;3)(p15;q26.2)", "t(7;14)(q21;q32)")
MPAL_COLOR <- "#0E7C7B"
LEAD_COLORS  <- c("CD96" = "#c0392b", "CD9" = "#0072B2", "SUCNR1" = "#009E73",
                  "AMN" = "#E69F00", "IL1RAP" = "#56B4E9")   # IL1RAP = muted purple (was pink)
# The 12 panel-E antigens (= the genes in Figure 7B/C) with the shared group colour scheme
# (matches Fig 5E): leads = identity colours, FLT3 control = gold, clinical = orange,
# PPAC-specific = purple. Keyed by display label, CD96 first. Used by the window supplements
# (S23 per-PAC, S24 per-cytogenetic, S25 per-PPAC) so they match the main combination figures.
TARGET12_COLORS <- c(CD96 = "#c0392b", CD9 = "#0072B2", SUCNR1 = "#009E73", IL1RAP = "#56B4E9",
                     FLT3 = "#b7950b", CD33 = "#e8590c", CD123 = "#e8590c", "CLL-1" = "#e8590c",
                     CD7 = "#e8590c", TNFRSF4 = "#7f8c8d", ABCA7 = "#7f8c8d", ITGAX = "#7f8c8d")
# CD9 is broadly expressed on normal B-lineage progenitors + HSPC (see safety panel);
# that caveat rides in the HSPC/Myeloid toxicity row, NOT a label flag (shown plain per
# user preference -- no dagger).
LEAD_FLAG <- character(0)
lead_lab  <- function(v) ifelse(as.character(v) %in% LEAD_FLAG,
                                paste0(as.character(v), "†"), as.character(v))
lead_md   <- function(v) sprintf("<b style='color:%s'>%s</b>",
                                 LEAD_COLORS[as.character(v)], lead_lab(v))
# dot plot keyed to the 5 LEAD targets (CD96 top). size = %positive; colour = mean_expr
# if a `mean_expr` column is present, else %positive. df long: label, <xcol>, pct[, mean_expr].
dotplot_leads <- function(df, xcol, xlevels = NULL, maxsize = 6, xlab = "") {
  df$gene_y <- factor(df$label, levels = rev(LEAD_TARGETS))
  if (!is.null(xlevels)) df[[xcol]] <- factor(df[[xcol]], levels = xlevels)
  ytop <- length(LEAD_TARGETS)
  has_expr <- "mean_expr" %in% names(df)
  df$fillv <- if (has_expr) df$mean_expr else df$pct
  ggplot(df, aes(.data[[xcol]], gene_y)) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = ytop - 0.5, ymax = ytop + 0.5,
             fill = "#fbecea") +
    geom_point(aes(size = pct, fill = fillv), shape = 21, colour = "#3b5870", stroke = 0.2) +
    scale_fill_gradientn(colours = SEQ, oob = scales::squish,
                         name = if (has_expr) "mean expr." else "% positive",
                         guide = if (has_expr) "colourbar" else "none") +
    scale_size_area(max_size = maxsize, limits = c(0, 100), name = "% positive",
                    breaks = c(10, 25, 50, 75)) +
    labs(x = xlab, y = NULL) + theme_pub() +
    theme(axis.text.y = ggtext::element_markdown(),
          axis.text.x = element_text(angle = 35, hjust = 1)) +
    scale_y_discrete(limits = rev(LEAD_TARGETS), labels = function(v) lead_md(v))
}
# Sequential fill for dot plots / heatmaps = the manuscript heatmap0 ramp (matches Figure 3's
# scale_fill_seq exactly, so every dot-plot mean-expression fill reads the same across the paper).
SEQ    <- c("#001219","#005F73","#0A9396","#94D2BD","#E9D8A6","#EE9B00","#CA6702","#AE2012","#9B2226")  # heatmap0
CD96SEQ <- c("#fdf4f2","#f7d2cb","#ed9e90","#df6450","#c0392b","#92271c","#5f160f")
PROG   <- c(favorable = "#0A9396", poor = "#AE2012")   # SEQ endpoints (teal / red) — matches Figure 2/3
# NB: a named base_family (e.g. "Helvetica") breaks gridtext/element_markdown
# rendering on this device, so we use the device default sans family.
BASEFONT <- ""

# NB: base on theme_grey (NOT theme_minimal) -- theme_minimal() breaks
# ggtext::element_markdown rendering on this stack.  Styled to look clean.
theme_pub <- function(base = 10) {
  theme_grey(base_size = base, base_family = BASEFONT) +
    theme(panel.background = element_blank(),
          panel.border = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_line(color = "#eef2f6", linewidth = 0.3),
          axis.ticks = element_line(color = "#cfd6dd", linewidth = 0.3),
          plot.title = element_text(face = "bold", size = base),
          plot.title.position = "plot",
          legend.key = element_blank(),
          legend.key.size = unit(0.35, "cm"),
          axis.title = element_text(size = base - 1))
}

# colour the CD96 tick label crimson+bold via ggtext markdown
md_lab <- function(x) ifelse(x == "CD96",
  "<b style='color:#c0392b'>CD96</b>", x)
# colour every gene tick label by its target colour (bold)
gene_md <- function(v) sprintf("<b style='color:%s'>%s</b>", TARGET_COLORS[v], v)

# Pick the best AVAILABLE vector-PDF device that embeds fonts. The base pdf()
# device substitutes fonts on other machines (shifted text); cairo_pdf embeds
# but needs XQuartz; macOS quartz(type="pdf") embeds natively with no X11.
# capabilities("cairo") can report TRUE while the DLL fails to load, so probe.
.probe_dev <- function(opener) tryCatch({
  tmp <- tempfile(fileext = ".pdf"); opener(tmp); grDevices::dev.off()
  ok <- file.exists(tmp) && file.info(tmp)$size > 1000; unlink(tmp); ok
}, error = function(e) FALSE, warning = function(w) FALSE)
.PDF_DEVICE <- local({
  if (.probe_dev(function(f) grDevices::cairo_pdf(f, width = 4, height = 4)))
    return(grDevices::cairo_pdf)
  if (.probe_dev(function(f) grDevices::quartz(file = f, type = "pdf", width = 4, height = 4)))
    return(function(filename, width, height, ...)
      grDevices::quartz(file = filename, type = "pdf", width = width, height = height))
  grDevices::pdf  # last resort: fonts NOT embedded
})

savefig <- function(p, name, w, h) {
  ggplot2::ggsave(file.path(FIG, paste0(name, ".png")), p, width = w, height = h, create.dir = TRUE,
                  dpi = 300, bg = "white", device = ragg::agg_png)
  ggplot2::ggsave(file.path(FIG, paste0(name, ".pdf")), p, width = w, height = h, create.dir = TRUE,
                  bg = "white", device = .PDF_DEVICE)
  cat("  wrote ", name, " (", w, "x", h, ")\n", sep = "")
}

# Reusable dot plot: genes on Y (CD96 top), categories on X.  Standard single-cell
# encoding: dot SIZE = % positive; dot COLOUR = mean expression where a `mean_expr`
# column is present, else % positive (graceful).  Gene labels carry the per-gene
# identity colour.  CD96 stays flagged by the pink row band + crimson label.
# `df` long with columns: label (gene display), <xcol>, pct, [mean_expr].
dotplot <- function(df, xcol, xlevels = NULL, vmax = 70, maxsize = 7,
                    xlab = "", ang = 35) {
  df$gene_y <- factor(df$label, levels = rev(GENES))   # 'label' is reserved in aes()
  if (!is.null(xlevels)) df[[xcol]] <- factor(df[[xcol]], levels = xlevels)
  ytop <- length(GENES)
  has_expr <- "mean_expr" %in% names(df)
  df$fillv <- if (has_expr) df$mean_expr else df$pct
  ggplot(df, aes(x = .data[[xcol]], y = gene_y)) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = ytop - 0.5, ymax = ytop + 0.5,
             fill = "#fbecea") +
    geom_point(aes(size = pct, fill = fillv), shape = 21, colour = "#3b5870",
               stroke = 0.2) +
    scale_fill_gradientn(colours = SEQ, oob = scales::squish,
                         name = if (has_expr) "mean expr." else "% positive",
                         guide = if (has_expr) "colourbar" else "none") +
    scale_size_area(max_size = maxsize, limits = c(0, 100), name = "% positive",
                    breaks = c(10, 25, 50, 75)) +
    labs(x = xlab, y = NULL) +
    theme_pub() +
    theme(axis.text.y = ggtext::element_markdown(),
          axis.text.x = element_text(angle = ang, hjust = 1)) +
    scale_y_discrete(limits = rev(GENES), labels = function(v) gene_md(v))
}

# Dot plot with genes on X (CD96 first, all gene labels target-coloured),
# categories on Y.  Dot SIZE = % positive; dot COLOUR = mean expression where a
# `mean_expr` column is present, else % positive.  `df` long: label, <ycol>, pct.
dotplot_genesx <- function(df, ycol, ylevels = NULL, vmax = 70, maxsize = 6,
                           ylab = NULL) {
  df$gx <- factor(df$label, levels = GENES)
  if (!is.null(ylevels)) df[[ycol]] <- factor(df[[ycol]], levels = ylevels)
  has_expr <- "mean_expr" %in% names(df)
  df$fillv <- if (has_expr) df$mean_expr else df$pct
  ggplot(df, aes(gx, .data[[ycol]])) +
    annotate("rect", xmin = 0.5, xmax = 1.5, ymin = -Inf, ymax = Inf,
             fill = "#fbecea") +
    geom_point(aes(size = pct, fill = fillv), shape = 21, colour = "#3b5870",
               stroke = 0.2) +
    scale_fill_gradientn(colours = SEQ, oob = scales::squish,
                         name = if (has_expr) "mean expr." else "% positive",
                         guide = if (has_expr) "colourbar" else "none") +
    scale_size_area(max_size = maxsize, limits = c(0, 100), name = "% positive",
                    breaks = c(10, 25, 50, 75)) +
    labs(x = NULL, y = ylab) + theme_pub() +
    theme(axis.text.x = ggtext::element_markdown(angle = 35, hjust = 1)) +
    scale_x_discrete(limits = GENES, labels = function(v) gene_md(v))
}

# Shared bulk-expression violin (used by Fig 5E/F and the non-MLL supplement). Genes ordered
# by `ord`; coloured by group (lead = identity, control = gold, clinical = orange, ppac = purple).
# `vio` cols: gene, log2cpm, label, group.  `meta` cols: gene, label, group, median,
# pct_targetable, hspc_mye.  Top: median + %targetable(>5); bottom: HSPC/Myeloid toxicity.
violin_panel <- function(vio, meta, ord, title, ylab = "bulk expression in AML (log2 CPM)") {
  meta <- meta |>
    dplyr::mutate(col = dplyr::case_when(group == "lead" ~ as.character(LEAD_COLORS[gene]),
                                         group == "control" ~ "#b7950b",
                                         group == "ppac" ~ PPAC_COLOR,
                                         TRUE ~ "#e8590c"))
  labcol <- setNames(meta$col, as.character(meta$label))
  vio <- vio |> dplyr::left_join(meta[, c("label", "col")], by = "label") |>
    dplyr::mutate(label = factor(label, levels = ord))
  meta <- meta |> dplyr::mutate(label = factor(label, levels = ord))
  ytop <- max(vio$log2cpm, na.rm = TRUE)
  ggplot(vio, aes(label, log2cpm)) +
    geom_violin(aes(fill = col), scale = "width", linewidth = 0.3, colour = "#3b5870") +
    stat_summary(fun = median, geom = "crossbar", width = 0.5, linewidth = 0.3) +
    geom_hline(yintercept = 5, linetype = "dashed", linewidth = 0.4) +
    scale_fill_identity() +
    geom_text(data = meta, aes(label, ytop + 1.85, label = sprintf("med %.1f", median)),
              size = 2.4, colour = "grey30") +
    geom_text(data = meta, aes(label, ytop + 0.8, label = sprintf("%.0f%%", pct_targetable),
                               colour = col), fontface = "bold", size = 3) +
    scale_colour_identity() +
    geom_text(data = meta, aes(label, -1.6, label = sprintf("HSPC/Mye\n%.0f%%", hspc_mye)),
              size = 2, colour = ifelse(meta$hspc_mye > 10, "#c0392b", "#1e8449")) +
    labs(x = NULL, y = ylab, title = title) +
    coord_cartesian(ylim = c(-2.2, ytop + 2.4), clip = "off") +
    theme_pub(9) + theme(axis.text.x = ggtext::element_markdown(),
                         plot.margin = margin(6, 6, 14, 6)) +
    scale_x_discrete(labels = function(v) sprintf("<b style='color:%s'>%s</b>", labcol[v], v))
}
