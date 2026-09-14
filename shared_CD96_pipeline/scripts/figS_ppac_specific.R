# ====================================================================
# figS_ppac_specific.R  |  Figure S14 (PPAC-specific surface-target discovery)
# reproducibility header (auto-generated; edit freely below this block)
# Inputs : funnel_ppac_specific.csv, funnel_ppac_specific_composition.csv
# Outputs: (writes to ../figures or ../source_data — see body)
# Run    : Rscript figS_ppac_specific.R   (or: bash run_all.sh)  |  Deps: see REPRODUCIBILITY_CHECKLIST.md
# Source data for this figure is in __SUBMISSION_PACKAGE/08_Source_Data/.
# ====================================================================

# Supplementary figure — PPAC-SPECIFIC discovery (S17).
# Markers are ranked WITHIN each poor-prognosis cluster (PPAC) by the simplest
# PPAC-specific efficacy metric: the % of that cluster's cells that are positive
# (no composite). Candidate pool first passes a toxicity funnel with tolerance
# raised to <=10% on BOTH normal compartments (HSPC<=10%, Myeloid<=10%, TCR/CD3
# dropped). One lollipop per group; CD96 in crimson (not forced to the top).
# PPAC_5 is a single-patient cluster -> flagged; its hits are not generalizable.
# Source: funnel_ppac_specific.csv (+ _composition.csv).  Builder: build_funnel_ppac_specific.py
source("_style.R")

d    <- sd("funnel_ppac_specific.csv")
comp <- sd("funnel_ppac_specific_composition.csv")
TOPN <- 10

# CD96 crimson + bold on the y axis; everything else dark slate
lab_md <- function(x) ifelse(x == "CD96", "<b style='color:#c0392b'>CD96</b>",
                             sprintf("<span style='color:#34495e'>%s</span>", x))

panel <- function(g, title) {
  ci   <- comp[comp$group == g, ]
  flag <- isTRUE(ci$single_patient_flag[1])
  dd   <- d |> dplyr::filter(group == g) |>
    dplyr::arrange(dplyr::desc(cov_pct)) |> head(TOPN)
  dd$gene_f <- factor(dd$label, levels = rev(dd$label))
  dd$is96   <- dd$label == "CD96"
  cd <- dd[dd$is96, ]
  sub <- if (nrow(cd))
    sprintf("CD96: %.0f%% of cells positive · %d/%d pts",
            cd$cov_pct, ci$n_patients_contrib[1], ci$n_patients_any[1])
  else sprintf("%d patients contribute (CD96 below top %d)", ci$n_patients_contrib[1], TOPN)

  p <- ggplot(dd, aes(cov_pct, gene_f))
  if (flag)  # red wash + artifact stamp behind single-patient PPAC
    p <- p + annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,
                      fill = "#fdecea") +
      annotate("text", x = max(dd$cov_pct) * 0.98, y = 1.4,
               label = sprintf("single patient\n(%.0f%% of cluster)\nnot generalizable",
                               ci$dominant_patient_share_pct[1]),
               hjust = 1, vjust = 0, size = 2.5, fontface = "italic", colour = "#a93226")
  p +
    geom_segment(aes(x = 0, xend = cov_pct, y = gene_f, yend = gene_f, colour = is96),
                 linewidth = 1.1) +
    geom_point(aes(colour = is96, size = is96)) +
    geom_text(aes(label = sprintf("%.1f", cov_pct), colour = is96),
              hjust = -0.35, size = 2.7, fontface = ifelse(dd$is96, "bold", "plain")) +
    scale_colour_manual(values = c(`TRUE` = PAL$CD96, `FALSE` = PAL$navy), guide = "none") +
    scale_size_manual(values = c(`TRUE` = 3.6, `FALSE` = 2.4), guide = "none") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.18))) +
    scale_y_discrete(labels = function(v) lab_md(v)) +
    labs(x = "% of PPAC cells positive", y = NULL, title = title, subtitle = sub) +
    theme_pub(8.5) +
    theme(axis.text.y = ggtext::element_markdown(size = 8),
          plot.subtitle = element_text(size = 6.6, colour = "#7f8c8d"))
}

pg  <- panel("All PPAC", "All PPACs (pooled)")
p1  <- panel("PPAC_1", "PPAC_1")
p2  <- panel("PPAC_2", "PPAC_2")
p3  <- panel("PPAC_3", "PPAC_3")
p4  <- panel("PPAC_4", "PPAC_4")
p5  <- panel("PPAC_5", "PPAC_5")

fig <- (pg | p1) / (p2 | p3) / (p4 | p5) +
  patchwork::plot_annotation(
    title = "PPAC-specific surface-target discovery",
    subtitle = paste0("Markers ranked within each poor-prognosis cluster by % of cells in the ",
                      "cluster positive, behind a toxicity funnel (HSPC<=10%, Myeloid<=10%, TCR/CD3 ",
                      "dropped; 1,706 eligible genes). CD96 in crimson (ranked naturally, not forced first)."),
    theme = theme(plot.title = element_text(face = "bold", size = 12),
                  plot.subtitle = ggtext::element_textbox_simple(size = 8, colour = "#34495e",
                                                                 margin = margin(b = 6))))

savefig(fig, "supplementary/FigureS14_ppac_specific", 11, 12)
cat("Done -> figures/supplementary/FigureS14_ppac_specific.{pdf,png}\n")
