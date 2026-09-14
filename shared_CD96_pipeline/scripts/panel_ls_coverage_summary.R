#!/usr/bin/env Rscript
# Figure 5 summary of the panel-D LS dot plot: for each surface antigen, how many leukemic states (LS)
# it could target, split by prognosis. LEFT = boxplot of % cells positive across the favorable vs poor
# prognostic LS (each point = one LS); RIGHT = count of LS targetable (>= TARGET_THRESH % positive),
# total / favorable / poor. Ordered by poor-LS coverage (the aggressive states we most need to hit).
# Single-cell throughout; TARGET-bulk validation is in panels E/F. Data = fig5D_ls.csv. No fabrication.
source("_style.R")
suppressPackageStartupMessages({library(dplyr); library(tidyr)})
TARGET_THRESH <- 25                                   # % cells positive to call an LS "targetable"
d <- as.data.frame(sd("fig5D_ls.csv")) |> filter(grepl("^LS_", group))    # 49 LS (drop normal baseline)
d$prognosis <- factor(d$prognosis, levels = c("favorable", "poor", "n.s."))
NFAV <- d |> filter(prognosis == "favorable") |> distinct(group) |> nrow()
NPOOR <- d |> filter(prognosis == "poor") |> distinct(group) |> nrow()
NTOT <- d |> distinct(group) |> nrow()

# per-antigen counts of targetable LS (total / favorable / poor)
cnt <- d |> group_by(label) |>
  summarise(total = sum(pct >= TARGET_THRESH),
            favorable = sum(pct >= TARGET_THRESH & prognosis == "favorable"),
            poor = sum(pct >= TARGET_THRESH & prognosis == "poor"),
            poormed = median(pct[prognosis == "poor"]), .groups = "drop") |>
  arrange(desc(poor), desc(poormed))
ord <- cnt$label
d$label <- factor(d$label, levels = ord)

# LEFT: boxplot of % positive per antigen, across ALL leukemic states (total) and the favorable / poor
# prognostic subsets. Colours = manuscript scheme (total grey, favorable teal, poor red = PROG).
GREY <- "#808080"
bx <- bind_rows(d |> mutate(cls = "total"),
                d |> filter(prognosis %in% c("favorable", "poor")) |> mutate(cls = as.character(prognosis)))
bx$cls <- factor(bx$cls, levels = c("total", "favorable", "poor"))
bx$label <- factor(bx$label, levels = ord)
PC <- c(total = GREY, favorable = unname(PROG["favorable"]), poor = unname(PROG["poor"]))
pL <- ggplot(bx, aes(label, pct, fill = cls)) +
  geom_hline(yintercept = TARGET_THRESH, linetype = "dashed", colour = "grey55", linewidth = 0.4) +
  geom_boxplot(outlier.shape = NA, width = 0.72, alpha = 0.35, linewidth = 0.3,
               position = position_dodge(0.8), colour = "grey40") +
  geom_point(aes(colour = cls), position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8),
             size = 0.8, alpha = 0.65) +
  scale_fill_manual(values = PC, name = "leukemic states") +
  scale_colour_manual(values = PC, guide = "none") +
  labs(x = NULL, y = "% cells positive (per leukemic state)",
       title = "Leukemic-state coverage per surface antigen (single-cell)",
       subtitle = sprintf("each point = one leukemic state (all %d; %d favorable, %d poor); dashed line = %d%% targetable threshold",
                          NTOT, NFAV, NPOOR, TARGET_THRESH)) +
  theme_pub(10) + theme(axis.text.x = ggtext::element_markdown(angle = 40, hjust = 1),
                        legend.position = c(0.99, 0.99), legend.justification = c(1, 1))

# RIGHT: number of LS targetable per antigen (total / favorable / poor)
cl <- cnt |> select(label, total, favorable, poor) |>
  pivot_longer(-label, names_to = "class", values_to = "n")
cl$class <- factor(cl$class, levels = c("total", "favorable", "poor"),
                   labels = c(sprintf("total (/%d)", NTOT), sprintf("favorable (/%d)", NFAV), sprintf("poor (/%d)", NPOOR)))
cl$label <- factor(cl$label, levels = rev(ord))
CC <- setNames(c(GREY, unname(PROG["favorable"]), unname(PROG["poor"])), levels(cl$class))
pR <- ggplot(cl, aes(n, label, colour = class)) +
  geom_segment(aes(x = 0, xend = n, yend = label), position = position_dodge(0.6), linewidth = 0.4, alpha = 0.5) +
  geom_point(position = position_dodge(0.6), size = 2) +
  scale_colour_manual(values = CC, name = NULL) +
  labs(x = sprintf("# leukemic states targetable (>= %d%% positive)", TARGET_THRESH), y = NULL,
       title = "Number of targetable leukemic states") +
  theme_pub(10) + theme(axis.text.y = ggtext::element_markdown(), legend.position = "bottom")

fig <- pL / pR + plot_layout(heights = c(1, 1.15)) + plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 14))
savefig(fig, "Figure5_LS_coverage_summary", 11, 12)
cat("\nTargetable LS per antigen (>=", TARGET_THRESH, "% positive), sorted by poor:\n")
print(cnt |> select(label, total, favorable, poor), n = 100)
