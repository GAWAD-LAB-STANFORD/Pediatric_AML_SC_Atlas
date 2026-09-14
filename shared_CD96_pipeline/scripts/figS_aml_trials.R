# ====================================================================
# figS_aml_trials.R  |  Figure S12 (AML immunotherapy clinical-trial landscape)
# reproducibility header (auto-generated; edit freely below this block)
# Inputs : fig_aml_trials.csv
# Outputs: (writes to ../figures or ../source_data — see body)
# Run    : Rscript figS_aml_trials.R   (or: bash run_all.sh)  |  Deps: see REPRODUCIBILITY_CHECKLIST.md
# Source data for this figure is in __SUBMISSION_PACKAGE/08_Source_Data/.
# ====================================================================

# Supplementary figure — AML surface-antigen immunotherapy clinical landscape (ClinicalTrials.gov).
# Number of interventional AML trials naming each surface antigen as the immunotherapy target
# (CAR-T / bispecific / ADC / mAb). The 4 lead targets are highlighted to show where
# they sit: CD96, CD9 and SUCNR1 are completely untapped; only IL1RAP has AML trials (2).
# Source: source_data/fig_aml_trials.csv (queried 2026-06-08). Methodology + caveats in caption.
source("_style.R")

d <- sd("fig_aml_trials.csv")
d <- d |> dplyr::arrange(n_trials, dplyr::desc(label))
d$label <- factor(d$label, levels = d$label)        # ascending for horizontal bars
# colour: leads get their identity colour (CD9 daggered); established targets navy
d$gene_key <- sub(" .*", "", d$antigen)
d$fill <- ifelse(d$is_lead, as.character(LEAD_COLORS[d$gene_key]), PAL$navy)
d$fill[is.na(d$fill)] <- PAL$navy
ax_lab <- function(v) {
  g <- d$gene_key[match(v, d$label)]
  ifelse(g %in% LEAD_TARGETS,
         sprintf("<b style='color:%s'>%s</b>", LEAD_COLORS[g],
                 v),
         as.character(v))
}

p <- ggplot(d, aes(n_trials, label)) +
  geom_col(aes(fill = fill), width = 0.74) +
  scale_fill_identity() +
  geom_text(aes(label = ifelse(n_trials == 0, "0 — none", as.character(n_trials)),
                colour = fill),
            hjust = -0.12, size = 2.9, fontface = "bold") +
  scale_colour_identity() +
  scale_x_continuous(expand = expansion(c(0, 0.13))) +
  scale_y_discrete(labels = function(v) ax_lab(v)) +
  labs(x = "number of interventional AML trials (ClinicalTrials.gov)", y = NULL,
       title = "AML surface-antigen immunotherapy clinical landscape",
       subtitle = NULL) +
  theme_pub(10) +
  theme(axis.text.y = ggtext::element_markdown(size = 9.5),
        plot.subtitle = ggtext::element_textbox_simple(size = 8, colour = "#34495e",
                                                       margin = margin(b = 8)),
        panel.grid.major.y = element_blank())

savefig(p, "supplementary/FigureS12_aml_trial_landscape", 8.5, 6.5)
cat("Done -> figures/supplementary/FigureS12_aml_trial_landscape.{pdf,png}\n")
