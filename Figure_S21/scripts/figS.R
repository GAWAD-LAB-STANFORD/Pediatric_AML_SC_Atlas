# ====================================================================
# figS.R  |  Supplementary figures (correlations S17, non-MLL window S18)
# reproducibility header (auto-generated; edit freely below this block)
# Inputs : fig5F_meta_all.csv, fig5F_violin.csv, figS17_corr_r.csv, figS17_points_bulk.csv, figS17_points_singlecell.csv, figS_cyto_pediatric_violin.csv
# Outputs: (writes to ../figures or ../source_data — see body)
# Run    : Rscript figS.R   (or: bash run_all.sh)  |  Deps: see REPRODUCIBILITY_CHECKLIST.md
# Source data for this figure is in __SUBMISSION_PACKAGE/08_Source_Data/.
# ====================================================================

# Supplementary figures (R/ggplot2):
#   S6  therapeutic-window violins, KMT2Ar/MLL excluded (all 12 candidate antigens)
#   S7  CD96 co-expression with the comparators (single cell + TARGET bulk)
source("_style.R")

## ---- S6: non-MLL window violins (KMT2Ar excluded) -------------------------------------------
# Self-contained: the KMT2Ar sample set comes from the cytogenetic-faceted supplement input
# (figS_cyto_pediatric_violin.csv, which carries a `cyto` column), and the non-MLL meta is
# recomputed here so no orphan file is needed. hspc_mye toxicity is joined from fig6F_meta_all.
cyto <- sd("figS_cyto_pediatric_violin.csv")
mll_samples <- unique(cyto$sample[cyto$cyto == "KMT2Ar"])
tox <- sd("fig5F_meta_all.csv") |> dplyr::select(label, hspc_mye)
vio <- sd("fig5F_violin.csv") |> dplyr::filter(!sample %in% mll_samples)
meta <- vio |> dplyr::group_by(gene, label, group) |>
  dplyr::summarise(median = round(median(log2cpm), 4),
                   pct_targetable = round(100 * mean(log2cpm > 5), 2),
                   n = dplyr::n(), .groups = "drop") |>
  dplyr::left_join(tox, by = "label") |>
  dplyr::arrange(dplyr::desc(median))
ord <- meta$label
s6 <- violin_panel(vio, meta, ord,
                   sprintf("KMT2Ar/MLL excluded (%d of 475 \u2265 80%%-blast TARGET samples)", dplyr::n_distinct(vio$sample)))
savefig(s6, "supplementary/FigureS18_nonMLL_violins", 12, 5.2)

## ---- S7: CD96 co-expression scatter (single cell + bulk) ----
SUBCOL <- c("CN" = "#0072B2", "t(8;21)" = "#009E73", "inv(16)" = "#E69F00",
            "KMT2Ar" = "#c0392b", "PML-RARA" = "#CC79A7", "NUP98r" = "#56B4E9",
            "Other" = "#b9c2cb")
collapse_sub <- function(s) dplyr::case_when(
  s %in% c("CN", "t(8;21)", "inv(16)", "KMT2Ar", "PML-RARA") ~ s,
  grepl("NUP98", s) ~ "NUP98r", TRUE ~ "Other")
sc <- sd("figS17_points_singlecell.csv") |> dplyr::mutate(modality = "AML Metacells")
bk <- sd("figS17_points_bulk.csv") |> dplyr::mutate(modality = "TARGET Bulk")
pts <- dplyr::bind_rows(sc, bk) |>
  dplyr::mutate(subtype = factor(collapse_sub(subtype), levels = names(SUBCOL)),
                marker = factor(marker, levels = names(TARGET12_COLORS)[-1]),
                modality = factor(modality, levels = c("AML Metacells", "TARGET Bulk")))
rr <- sd("figS17_corr_r.csv") |>
  dplyr::mutate(modality = ifelse(modality == "single_cell", "AML Metacells", "TARGET Bulk"),
                modality = factor(modality, levels = c("AML Metacells", "TARGET Bulk")),
                marker = factor(label, levels = names(TARGET12_COLORS)[-1]))
s7 <- ggplot(pts, aes(cd96, val)) +
  # metacells (~75-cell pools) + bulk samples -> few thousand points; plot direct
  geom_point(aes(colour = subtype), size = 0.85, alpha = 0.6) +
  scale_colour_manual(values = SUBCOL, name = "cytogenetic subtype",
                      guide = guide_legend(override.aes = list(size = 2.5, alpha = 1))) +
  geom_text(data = rr, aes(x = -Inf, y = Inf, label = sprintf("r=%+.2f", pearson)),
            hjust = -0.12, vjust = 1.4, size = 2.7, fontface = "bold",
            colour = "#1c2833", inherit.aes = FALSE) +
  # independent x per facet so the metacell row (CD96 ~0-2) fills its width
  # instead of sharing the bulk row's wider CD96 (log2 CPM) axis
  ggh4x::facet_grid2(modality ~ marker, scales = "free", independent = "x") +
  labs(x = "CD96", y = "target") + theme_pub(9) +
  theme(legend.position = "bottom", strip.text = element_text(size = 8),
        panel.spacing = unit(3, "pt"))
savefig(s7, "supplementary/FigureS17_correlations", 15, 5.4)
cat("done supplements\n")
