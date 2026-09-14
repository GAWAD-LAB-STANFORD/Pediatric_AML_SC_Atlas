# ====================================================================
# figS_ls_specific.R  |  Poor-prognosis-STATE (LS) surface-target discovery
# Replaces the retired per-PPAC version (figS_ppac_specific.R) on the corrected
# 49-state scheme. Markers ranked WITHIN each poor-prognosis leukemic state by
# % of that state's cells positive, behind the same toxicity funnel
# (HSPC<=10%, Myeloid<=10%, TCR/CD3 dropped). CD96 in crimson (ranked naturally).
# All panels are multi-patient (generalizable); the 3 single-patient poor LS
# (LS_23/LS_38/LS_41) are pooled into "All poor LS" but given no panel.
# Source: funnel_ls_specific.csv (+ _composition.csv). Builder: build_funnel_ls_specific.py
# ====================================================================
source("_style.R")

d    <- sd("funnel_ls_specific.csv")
comp <- sd("funnel_ls_specific_composition.csv")
TOPN <- 10

# lineage program per featured LS (for panel titles); order = by cell number
PANELS <- c("All poor LS" = "All poor-prognosis LS (pooled)",
            "LS_10" = "LS_10 (proliferative)", "LS_19" = "LS_19 (pDC)",
            "LS_20" = "LS_20 (proliferative)", "LS_21" = "LS_21 (macrophage)",
            "LS_27" = "LS_27 (GMP)")

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
  if (flag)
    p <- p + annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "#fdecea") +
      annotate("text", x = max(dd$cov_pct) * 0.98, y = 1.4,
               label = sprintf("single patient\n(%.0f%% of state)\nnot generalizable",
                               ci$dominant_patient_share_pct[1]),
               hjust = 1, vjust = 0, size = 2.5, fontface = "italic", colour = "#a93226")
  p +
    geom_segment(aes(x = 0, xend = cov_pct, y = gene_f, yend = gene_f, colour = is96), linewidth = 1.1) +
    geom_point(aes(colour = is96, size = is96)) +
    geom_text(aes(label = sprintf("%.1f", cov_pct), colour = is96),
              hjust = -0.35, size = 2.7, fontface = ifelse(dd$is96, "bold", "plain")) +
    scale_colour_manual(values = c(`TRUE` = PAL$CD96, `FALSE` = PAL$navy), guide = "none") +
    scale_size_manual(values = c(`TRUE` = 3.6, `FALSE` = 2.4), guide = "none") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.18))) +
    scale_y_discrete(labels = function(v) lab_md(v)) +
    labs(x = "% of state cells positive", y = NULL, title = title, subtitle = sub) +
    theme_pub(8.5) +
    theme(axis.text.y = ggtext::element_markdown(size = 8),
          plot.subtitle = element_text(size = 6.6, colour = "#7f8c8d"))
}

plist <- Map(panel, names(PANELS), unname(PANELS))
fig <- (plist[[1]] | plist[[2]]) / (plist[[3]] | plist[[4]]) / (plist[[5]] | plist[[6]]) +
  patchwork::plot_annotation(
    title = "Poor-prognosis-state surface-target discovery",
    subtitle = NULL,
    theme = theme(plot.title = element_text(face = "bold", size = 12),
                  plot.subtitle = ggtext::element_textbox_simple(size = 8, colour = "#34495e",
                                                                 margin = margin(b = 6))))

savefig(fig, "supplementary/FigureS15_ls_specific", 11, 12)
cat("Done -> figures/supplementary/FigureS15_ls_specific.{pdf,png}\n")
