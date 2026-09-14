# ====================================================================
# fig7_coverage.R  |  Analysis asset (cumulative coverage; not a final Fig-7 panel)
# reproducibility header (auto-generated; edit freely below this block)
# Inputs : fig5E_coverage.csv
# Outputs: (writes to ../figures or ../source_data — see body)
# Run    : Rscript fig7_coverage.R   (or: bash run_all.sh)  |  Deps: see REPRODUCIBILITY_CHECKLIST.md
# Source data for this figure is in __SUBMISSION_PACKAGE/08_Source_Data/.
# ====================================================================

# Figure 7 panel B — cumulative patient coverage as targets are added to CD96.
# CD96-anchored greedy: start from CD96, then add the partner covering the most still-uncovered
# patients. Patient "covered" = >=20% of blasts positive for >=1 panel target (n = 28). The curve
# climbs CD96 61% -> +CD33 86% -> +SUCNR1 89% -> +ITGAX 93%, then plateaus (CLL-1/CD7/TNFRSF4 +0).
source("_style.R")
cov <- as.data.frame(sd("fig5E_coverage.csv"))
cov$group[cov$label == "CD7"] <- "clinical"        # corrected scheme: CD7 = clinical CAR-T antigen
GCOL <- c(lead = "#2980b9", clinical = "#e8590c", control = "#b7950b", ppac = "#7f8c8d")
GLAB <- c(lead = "lead candidate", clinical = "clinical antigen", ppac = "additional candidate")
cov$xlab <- ifelse(cov$step == 1, "CD96", paste0("+ ", cov$label))
cov$xlab <- factor(cov$xlab, levels = cov$xlab[order(cov$step)])
cov$group <- factor(cov$group, levels = names(GCOL))
cov$plateau <- cov$added == 0                      # targets that add no new patients

p <- ggplot(cov, aes(step, cum_pct)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 100, ymax = Inf, fill = "white") +
  geom_hline(yintercept = 100, linetype = "dotted", colour = "grey70", linewidth = 0.3) +
  geom_area(fill = "#dcebf7", alpha = 0.7) +
  geom_line(colour = "#2471a3", linewidth = 0.8) +
  geom_point(aes(fill = group, alpha = plateau), shape = 21, size = 5, stroke = 0.5, colour = "white") +
  geom_text(aes(label = sprintf("%.0f%%", cum_pct)), vjust = -1.25, size = 3.2, fontface = "bold") +
  geom_text(data = subset(cov, !plateau & step > 1), aes(label = sprintf("+%d pts", added)),
            vjust = 2.4, size = 2.6, colour = "grey40") +
  geom_text(data = subset(cov, plateau), aes(y = cum_pct, label = "no further\ngain"),
            vjust = 1.9, size = 2.4, colour = "grey55", lineheight = 0.85) +
  scale_alpha_manual(values = c(`FALSE` = 1, `TRUE` = 0.35), guide = "none") +
  scale_fill_manual(values = GCOL, labels = GLAB, name = NULL, drop = FALSE) +
  scale_x_continuous(breaks = cov$step, labels = levels(cov$xlab), expand = expansion(c(0.05, 0.08))) +
  scale_y_continuous(limits = c(0, 108), breaks = seq(0, 100, 20)) +
  labs(title = "Cumulative patient coverage, building from CD96",
       subtitle = "greedy best-add order · patient covered = ≥20% of blasts positive for ≥1 panel target · n = 28 AML patients",
       x = NULL, y = "% patients covered (cumulative)") +
  theme_pub(11) +
  theme(legend.position = c(0.99, 0.04), legend.justification = c(1, 0),
        legend.background = element_rect(fill = "white", colour = "grey85", linewidth = 0.3),
        legend.key.size = unit(0.32, "cm"),
        axis.text.x = element_text(face = "bold"),
        plot.title = element_text(face = "bold", size = 12))

savefig(p, "fig7_coverage", 8, 6.4)
cat("wrote fig7_coverage\n"); print(cov[, c("step", "xlab", "added", "cum_pct")], row.names = FALSE)
