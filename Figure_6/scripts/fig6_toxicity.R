# ====================================================================
# fig6_toxicity.R  |  Figure 6A-C (normal-tissue toxicity dot plots)
# reproducibility header (auto-generated; edit freely below this block)
# Inputs : fig6tox_celltypes.csv, fig6tox_hema.csv, fig6tox_organ.csv
# Outputs: .panel_fig6ABC.pdf, _review_fig6ABC.png
# Run    : Rscript fig6_toxicity.R   (or: bash run_all.sh)  |  Deps: see REPRODUCIBILITY_CHECKLIST.md
# Source data for this figure is in __SUBMISSION_PACKAGE/08_Source_Data/.
# ====================================================================

# NEW Figure 6 (toxicity disqualification) — 3 stacked dot plots, each independently
# clustered by marker (own dendrogram):
#   A  major-organ atlas (parenchymal)
#   B  hematopoiesis atlas
#   C  each marker's top-3 SEVERE-toxicity cell types (12 x 3 -> ~36 specific cell types)
#   dot SIZE = % of normal cells positive ; dot COLOUR = mean expr scaled per gene.
source("_style.R")
suppressPackageStartupMessages({library(patchwork); library(ggdendro)})

GENES12 <- c("CD96", "CD9", "SUCNR1", "IL1RAP", "FLT3", "CD33", "CD123", "CLL-1",
             "CD7", "TNFRSF4", "ABCA7", "ITGAX")
GCOL <- c(CD96 = "#c0392b", CD9 = "#0072B2", SUCNR1 = "#009E73", IL1RAP = "#CC79A7",
          FLT3 = "#b7950b", CD33 = "#e8590c", CD123 = "#e8590c", "CLL-1" = "#e8590c",
          CD7 = "#e8590c", TNFRSF4 = "#7f8c8d", ABCA7 = "#7f8c8d", ITGAX = "#7f8c8d")
gmd <- function(v) sprintf("<b style='color:%s'>%s</b>", GCOL[v], v)

# dot plot: genes on Y (order = gorder, bottom->top), categories on X.
tox_dot <- function(df, xcol, xlevels, title, gorder, xlab = "", sep_after = NULL, ylabs = TRUE) {
  df <- df |> dplyr::filter(.data[[xcol]] %in% xlevels, label %in% gorder)
  df <- df |> dplyr::group_by(label) |>
    dplyr::mutate(fillv = ifelse(max(mean, na.rm = TRUE) > 0, mean / max(mean, na.rm = TRUE), 0)) |>
    dplyr::ungroup()
  df[[xcol]] <- factor(df[[xcol]], levels = xlevels)
  df$gy <- factor(df$label, levels = gorder)
  p <- ggplot(df, aes(.data[[xcol]], gy))
  if (!is.null(sep_after))
    p <- p + geom_vline(xintercept = sep_after + 0.5, colour = "grey70", linewidth = 0.4, linetype = "dashed")
  p <- p + geom_point(aes(size = pct, fill = fillv), shape = 21, colour = "#3b5870", stroke = 0.2) +
    scale_fill_gradientn(colours = SEQ, limits = c(0, 1), oob = scales::squish, name = "expr.\n(per-gene)") +
    scale_size_area(max_size = 10, limits = c(0, 50), breaks = c(5, 10, 25, 50),
                    oob = scales::squish, name = "% normal\ncells positive") +
    scale_x_discrete(limits = xlevels) +
    scale_y_discrete(limits = gorder, labels = function(v) gmd(v)) +
    labs(x = xlab, y = NULL, title = title) + theme_pub(9) +
    theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 6.6))
  if (ylabs) p + theme(axis.text.y = ggtext::element_markdown())
  else p + theme(axis.text.y = element_blank())
}

# cluster the 12 markers on a (gene x category) matrix; return order + aligned dendrogram
matof <- function(df, keycol) {
  d <- df |> dplyr::transmute(label, k = .data[[keycol]], pct)
  w <- tidyr::pivot_wider(d, names_from = k, values_from = pct, values_fn = max)
  m <- as.matrix(w[, -1]); rownames(m) <- w$label; m
}
cluster_genes <- function(mat) {
  mat <- mat[GENES12[GENES12 %in% rownames(mat)], , drop = FALSE]
  mat[is.na(mat)] <- 0
  hc <- hclust(dist(mat), "ward.D2")
  ord <- rownames(mat)[hc$order]
  dseg <- ggdendro::dendro_data(as.dendrogram(hc))$segments
  dend <- ggplot(dseg) +
    geom_segment(aes(x = -y, y = x, xend = -yend, yend = xend), linewidth = 0.3, colour = "#95a5a6") +
    scale_y_continuous(limits = c(0.4, nrow(mat) + 0.6), expand = c(0, 0)) +
    scale_x_continuous(expand = expansion(mult = c(0.06, 0))) + theme_void()
  list(order = ord, dend = dend)
}

# ---------- A: major organs (parenchymal), clustered by marker ----------
org <- sd("fig6tox_organ.csv")
MAJOR <- c("heart" = "Heart", "brain" = "Brain", "lung" = "Lung", "liver" = "Liver",
           "kidney" = "Kidney", "pancreas" = "Pancreas", "stomach" = "Stomach",
           "small intestine" = "Small intestine", "colon" = "Colon", "esophagus" = "Esophagus",
           "endocrine gland" = "Endocrine", "skin of body" = "Skin", "adipose tissue" = "Adipose",
           "musculature" = "Skeletal muscle", "bladder organ" = "Bladder", "prostate gland" = "Prostate",
           "uterus" = "Uterus", "breast" = "Breast", "eye" = "Eye")
org <- org |> dplyr::filter(organ %in% names(MAJOR)) |> dplyr::mutate(organ = unname(MAJOR[organ]))
aord <- unname(MAJOR)[unname(MAJOR) %in% unique(org$organ)]
cA <- cluster_genes(matof(org, "organ"))
pA <- tox_dot(org, "organ", aord, "A", cA$order)

# ---------- B: hematopoiesis, clustered by marker ----------
hem <- sd("fig6tox_hema.csv")
hord <- unique(hem$cell)
cB <- cluster_genes(matof(hem, "cell"))
pB <- tox_dot(hem, "cell", hord, "B", cB$order)

# ---------- C: each marker's top-3 SEVERE-toxicity cell types (specific, ~36) ----------
ctp <- sd("fig6tox_celltypes.csv")
SEV_KW <- c("cardiac muscle", "cardiomyocyte", "neuron", "retina", "cone cell", "rod cell",
            "photoreceptor", "lens fiber", "corneal", "pneumocyte", "alveolar epithelial",
            "alveolar type", "hepatocyte", "pancreatic", "islet", "acinar", "kidney", "renal",
            "nephron", "podocyte", "urothelial", "umbrella cell", "enterocyte", "intestinal epithelial",
            "epithelial cell of small intestine", "epithelial cell of large intestine", "m cell of gut",
            "enteroendocrine", "keratinocyte", "spinous", "epidermis", "endothelial cell")
BLOCK <- c("trophoblast", "placent", "embryo", "fetal", "yolk", "syncytio", "cytotroph", "tip cell",
           "eurydendroid", "be cell", "blastomere", "morula", "primordial", "unknown", "pluripotent",
           "stem cell", "progenitor")
sev <- ctp |> dplyr::filter(n >= 1000,
                            grepl(paste(SEV_KW, collapse = "|"), tolower(cell_type)),
                            !grepl(paste(BLOCK, collapse = "|"), tolower(cell_type)))
top3 <- sev |> dplyr::group_by(label) |> dplyr::slice_max(pct, n = 3, with_ties = FALSE) |> dplyr::ungroup()
shorten <- function(x) {
  x <- gsub("epithelial cell", "epith.", x)
  x <- gsub("GABAergic cortical interneuron", "GABAergic neuron", x)
  x <- gsub(" glutamatergic neuron.*", " glut. neuron", x)
  x <- gsub(" of the primary motor cortex| of Henle", "", x)
  x <- gsub(" cell$", "", x); trimws(substr(x, 1, 30))
}
cdat <- sev |> dplyr::filter(cell_type %in% unique(top3$cell_type)) |>
  dplyr::mutate(cell = shorten(cell_type))
# assign each cell type to an organ system; group the columns by organ (faceted)
ORG_ORDER <- c("Heart", "Brain", "Eye", "Lung", "Liver", "Pancreas", "Kidney",
               "Bladder", "Gut", "Skin", "Vasculature", "Other")
organ_of <- function(ct) dplyr::case_when(
  grepl("cardiac|cardiomyocyte|myocyte", ct) ~ "Heart",
  grepl("neuron|gabaergic|glutamatergic", ct) ~ "Brain",
  grepl("retina|cone|rod|photoreceptor|lens|cornea", ct) ~ "Eye",
  grepl("pneumocyte|alveolar|pulmonary", ct) ~ "Lung",
  grepl("hepatocyte|hepatoblast", ct) ~ "Liver",
  grepl("pancreat|islet|acinar", ct) ~ "Pancreas",
  grepl("kidney|renal|nephron|podocyte|henle", ct) ~ "Kidney",
  grepl("urothelial|umbrella", ct) ~ "Bladder",
  grepl("enterocyte|intestin|m cell of gut|enteroendocrine|paneth", ct) ~ "Gut",
  grepl("keratinocyte|spinous|epidermis", ct) ~ "Skin",
  grepl("endothel", ct) ~ "Vasculature", TRUE ~ "Other")
cdat$organ <- factor(organ_of(tolower(cdat$cell_type)), levels = ORG_ORDER)
cell_order <- cdat |> dplyr::group_by(cell) |>
  dplyr::summarise(org = dplyr::first(organ), mx = max(pct), .groups = "drop") |>
  dplyr::arrange(org, dplyr::desc(mx)) |> dplyr::pull(cell)
cC <- cluster_genes(matof(cdat, "cell"))
cdf <- cdat |> dplyr::group_by(label) |>
  dplyr::mutate(fillv = ifelse(max(mean, na.rm = TRUE) > 0, mean / max(mean, na.rm = TRUE), 0)) |>
  dplyr::ungroup() |>
  dplyr::mutate(cell = factor(cell, levels = cell_order), gy = factor(label, levels = cC$order))
pC <- ggplot(cdf, aes(cell, gy)) +
  geom_point(aes(size = pct, fill = fillv), shape = 21, colour = "#3b5870", stroke = 0.2) +
  scale_fill_gradientn(colours = SEQ, limits = c(0, 1), oob = scales::squish, name = "expr.\n(per-gene)") +
  scale_size_area(max_size = 10, limits = c(0, 50), breaks = c(5, 10, 25, 50), oob = scales::squish,
                  name = "% normal\ncells positive") +
  scale_y_discrete(limits = cC$order, labels = function(v) gmd(v)) +
  facet_grid(. ~ organ, scales = "free_x", space = "free_x") +
  labs(x = NULL, y = NULL, title = "C") +
  theme_pub(9) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 6.4),
        axis.text.y = ggtext::element_markdown(),
        strip.text = element_text(face = "bold", size = 7.2),
        strip.background = element_rect(fill = "#eef2f6", colour = NA),
        panel.spacing = grid::unit(0.12, "lines"))
cat("C cell types:", length(cell_order), "| organs:", paste(levels(droplevels(cdat$organ)), collapse = " "), "\n")

DW <- 0.055   # dendrogram width fraction
p <- (cA$dend + pA + plot_layout(widths = c(DW, 1))) /
     (cB$dend + pB + plot_layout(widths = c(DW, 1))) /
     (cC$dend + pC + plot_layout(widths = c(DW, 1))) +
  plot_layout(heights = c(1, 1, 1))                  # no combined title/subtitle; commentary -> caption
ggsave("../figures/_review_fig6ABC.png", p, width = 16, height = 13, dpi = 200, bg = "white",
       device = ragg::agg_png)
ggsave("../figures/.panel_fig6ABC.pdf", p, width = 16, height = 13, bg = "white",
       device = .PDF_DEVICE)                       # vector panel for the fitz assembler
cat("wrote _review_fig6ABC.png + .panel_fig6ABC.pdf\n")
