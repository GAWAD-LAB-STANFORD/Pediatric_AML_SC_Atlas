# ====================================================================
# figS_cyto_violin.R  |  Figures S19/S20 (candidate antigens by cytogenetic subtype)
# reproducibility header (auto-generated; edit freely below this block)
# Inputs : fig5F_meta_all.csv, figS_cyto_adult_violin.csv, figS_cyto_pediatric_violin.csv
# Outputs: (writes to ../figures or ../source_data — see body)
# Run    : Rscript figS_cyto_violin.R   (or: bash run_all.sh)  |  Deps: see REPRODUCIBILITY_CHECKLIST.md
# Source data for this figure is in __SUBMISSION_PACKAGE/08_Source_Data/.
# ====================================================================

# Supplements to Figure-5 panels E/F: the 8 targets violined, FACETED by cytogenetic subtype.
#   pediatric = TARGET bulk (>=80% blasts);  adult = BeatAML bulk (>=80% blasts).
# Same 8 genes, same log2(CPM+1) scale and gene order as the main violins; one facet per
# subtype (strip shows n). Violin fill by group (lead identity / FLT3 gold ctrl / clinical
# orange); jittered points show the actual sample n; crossbar = median; dashed = targetable (5).
source("_style.R")

# shared gene x-order = pediatric panel-E order (descending TARGET median) -> twins line up
GENE_ORD <- sd("fig5F_meta_all.csv") |>
  dplyr::arrange(dplyr::desc(median)) |> dplyr::pull(label)

grp_col <- function(meta_like) dplyr::case_when(
  meta_like$group == "lead"    ~ as.character(LEAD_COLORS[meta_like$gene]),
  meta_like$group == "control" ~ "#b7950b",
  meta_like$group == "ppac"    ~ PPAC_COLOR,
  TRUE                          ~ "#e8590c")

VIOLIN_MIN <- 10   # draw a kernel-density violin only when a facet has >= this many samples;
                   # smaller facets show raw points + median crossbar only (no misleading density)

facet_violin <- function(vio, facet_levels, title, subtitle, ncol = 4) {
  # per-label colour (group-based), and n per subtype for the strip labels
  lab_meta <- vio |> dplyr::distinct(gene, label, group)
  lab_meta$col <- grp_col(lab_meta)
  labcol <- setNames(lab_meta$col, lab_meta$label)
  ncyto <- vio |> dplyr::group_by(cyto) |>
    dplyr::summarise(n = dplyr::n_distinct(sample), .groups = "drop")
  strip <- setNames(sprintf("%s  (n=%d)", ncyto$cyto, ncyto$n), ncyto$cyto)
  big <- ncyto$cyto[ncyto$n >= VIOLIN_MIN]
  vio <- vio |>
    dplyr::left_join(lab_meta[, c("label", "col")], by = "label") |>
    dplyr::mutate(label = factor(label, levels = GENE_ORD),
                  cyto  = factor(cyto, levels = facet_levels))
  ggplot(vio, aes(label, log2cpm)) +
    geom_hline(yintercept = 5, linetype = "dashed", linewidth = 0.3, colour = "grey55") +
    geom_violin(data = dplyr::filter(vio, cyto %in% big),
                aes(fill = col), scale = "width", linewidth = 0.25, colour = "#3b5870") +
    geom_jitter(width = 0.13, height = 0, size = 0.45, alpha = 0.4, colour = "#2c3e50") +
    stat_summary(fun = median, geom = "crossbar", width = 0.55, linewidth = 0.28, colour = "black") +
    scale_fill_identity() +
    facet_wrap(~cyto, ncol = ncol, labeller = as_labeller(strip)) +
    labs(x = NULL, y = "bulk expression in AML (log2 CPM)", title = title, subtitle = subtitle) +
    theme_pub(9) +
    theme(axis.text.x = ggtext::element_markdown(angle = 40, hjust = 1, size = 7),
          strip.text = element_text(face = "bold", size = 8.5),
          strip.background = element_rect(fill = "#eef2f6", colour = NA),
          panel.spacing = grid::unit(0.5, "lines"),
          plot.subtitle = element_text(size = 8, colour = "grey30")) +
    scale_x_discrete(labels = function(v) sprintf("<b style='color:%s'>%s</b>", labcol[v], v))
}

## ---- pediatric (TARGET) ----
ped <- sd("figS_cyto_pediatric_violin.csv")
ped_levels <- c("t(8;21)", "inv(16)", "PML-RARA", "KMT2Ar", "NUP98r", "CN", "Other")
ped_levels <- ped_levels[ped_levels %in% unique(ped$cyto)]
pp <- facet_violin(
  ped, ped_levels,
  "Pediatric AML (TARGET bulk, >=80% blasts): target expression by cytogenetic subtype",
  NULL,
  ncol = 4)
savefig(pp, "supplementary/FigureS19_cyto_violins_pediatric", 15, 8)
cat("wrote supplementary/FigureS19_cyto_violins_pediatric\n")

## ---- adult (BeatAML) ----
adt <- sd("figS_cyto_adult_violin.csv")
adt_levels <- c("CBF [t(8;21)/inv(16)]", "PML-RARA", "KMT2Ar", "Normal karyotype", "Other")
adt_levels <- adt_levels[adt_levels %in% unique(adt$cyto)]
pa <- facet_violin(
  adt, adt_levels,
  "Adult AML (BeatAML bulk, >=80% blasts): target expression by cytogenetic subtype",
  NULL,
  ncol = 3)
savefig(pa, "supplementary/FigureS20_cyto_violins_adult", 13, 8)
cat("wrote supplementary/FigureS20_cyto_violins_adult\n")
