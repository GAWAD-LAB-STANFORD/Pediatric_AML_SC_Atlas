#!/usr/bin/env Rscript
# Figure 2D (v3, corrected 49 states): state-level stemness by whole-cohort EFS
# prognosis group. Three groups (Favorable / No assoc. / Poor); each dot = one leukemic
# state (per-state median of the precomputed per-cell score); one panel per score
# (LSC Eppert, HSC). Kruskal-Wallis + pairwise Wilcoxon. No fabricated values.
suppressMessages({library(ggplot2); library(ggpubr); library(patchwork)})
CODE <- "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/07_Code/Figure_2"
source(file.path(CODE, "manuscript_palette.R"))
D2  <- "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/08_Source_Data/Figure_2/v3"
OUT <- "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/04_Main_Figures/Figure_2_panels"

pc  <- read.csv(file.path(D2, "LS_stem_percell.csv"))
agg <- aggregate(cbind(LSC_Eppert, HSC) ~ LS, data = pc, median)     # per-state median
pr  <- read.csv(file.path(D2, "LS_wholecohort_prognosis.csv"))
pr$group <- ifelse(pr$prognosis == "favorable", "Favorable",
             ifelse(pr$prognosis == "poor", "Poor", "No assoc."))
d <- merge(agg, pr[, c("LS","group")], by = "LS", all.x = TRUE)
d$group[is.na(d$group)] <- "No assoc."
d$group <- factor(d$group, levels = c("Favorable","No assoc.","Poor"))
cat("states per group:\n"); print(table(d$group))
cols <- c(Favorable = unname(PROG["favorable"]), `No assoc.` = GREY, Poor = unname(PROG["poor"]))
cmp  <- list(c("Favorable","Poor"), c("Favorable","No assoc."), c("No assoc.","Poor"))

panel <- function(yvar, ylab){
  rng <- range(d[[yvar]]); pad <- diff(rng)
  ggplot(d, aes(group, .data[[yvar]])) +
    geom_boxplot(aes(fill = group), outlier.shape = NA, width = 0.62, alpha = 0.45, linewidth = 0.4) +
    geom_jitter(aes(colour = group), width = 0.13, height = 0, size = 2, alpha = 0.9) +
    scale_fill_manual(values = cols) + scale_colour_manual(values = cols) +
    stat_compare_means(comparisons = cmp, method = "wilcox.test",
                       tip.length = 0.01, size = 3, step.increase = 0.11) +
    stat_compare_means(method = "kruskal.test", size = 3, label.x = 0.7,
                       label.y = rng[2] + pad*0.62) +
    labs(x = NULL, y = ylab) + theme_classic(base_size = 11) +
    theme(legend.position = "none", axis.text.x = element_text(size = 10)) +
    coord_cartesian(ylim = c(rng[1] - pad*0.05, rng[2] + pad*0.9)) }

pL <- panel("LSC_Eppert", "per-state median LSC (Eppert) score")
pH <- panel("HSC",        "per-state median HSC score")
fig <- pL + pH + plot_layout(guides = "collect")
ggsave(file.path(OUT, "Figure_2D_stemness_by_prognosis.pdf"), fig, width = 8, height = 5)
ggsave(file.path(OUT, "Figure_2D_stemness_by_prognosis.png"), fig, width = 8, height = 5, dpi = 200)
for (v in c("LSC_Eppert","HSC")) {
  kw <- kruskal.test(d[[v]], d$group)$p.value
  cat(sprintf("\n%s: Kruskal-Wallis p=%.3g\n", v, kw))
  for (cc in cmp) cat(sprintf("  %s vs %s: Wilcoxon p=%.3g\n", cc[1], cc[2],
      wilcox.test(d[[v]][d$group==cc[1]], d[[v]][d$group==cc[2]])$p.value)) }
cat("\nwrote Figure_2D_stemness_by_prognosis\n")
