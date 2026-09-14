#!/usr/bin/env Rscript
# Figure 4 prognosis summary as a dot plot: 12 T/NK subsets x {EFS, OS}. Dot colour = log2 hazard ratio
# per SD (diverging: blue favorable / red poor), size = -log10 BH-FDR (significance). Deconvolved TARGET
# fractions, n=1314. Data = tnk_subset_survival.csv. No fabricated values.
suppressMessages({library(ggplot2); library(dplyr); library(scales)})
source("/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/07_Code/Figure_3/manuscript_palette.R")
V2  <- "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/08_Source_Data/Figure_4/v2"
OUT <- "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/04_Main_Figures/Figure_4_panels"
d <- read.csv(file.path(V2, "tnk_subset_survival.csv"))
d$log2HR <- log2(d$HR); d$nlFDR <- -log10(d$FDR)
efs <- d %>% filter(outcome=="EFS")
ord <- efs %>% arrange(HR) %>% pull(subset)   # favorable (low HR) at top
d$subset  <- factor(d$subset, levels = rev(ord))
d$outcome <- factor(d$outcome, levels = c("EFS","OS"))
# direction group (from EFS) -> one label per group to the LEFT of the subset names (via facet strip)
grp <- efs %>% mutate(dir = ifelse(FDR<0.10, ifelse(HR>1, "more → worse", "more → better"), "n.s.")) %>%
       select(subset, dir)
d <- merge(d, grp, by = "subset")
d$dir <- factor(d$dir, levels = c("more → better", "n.s.", "more → worse"))
g <- ggplot(d, aes(outcome, subset)) +
  geom_point(aes(size = nlFDR, fill = log2HR), shape = 21, colour = "grey30", stroke = 0.3) +
  geom_point(data = subset(d, FDR < 0.10), aes(size = nlFDR), shape = 21, colour = "black", stroke = 0.9, fill = NA) +
  facet_grid(dir ~ ., scales = "free_y", space = "free", switch = "y") +
  scale_fill_div(midpoint = 0, name = expression(log[2]~HR/SD)) +
  scale_size(range = c(1, 10), name = expression(-log[10]~FDR)) +
  labs(x = NULL, y = NULL, caption = "size = -log10 FDR;  black ring = FDR < 0.10") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.major = element_line(colour = "grey92"),
        axis.text.y = element_text(size = 10), axis.text.x = element_text(face = "bold"),
        strip.placement = "outside", strip.background = element_blank(),
        strip.text.y.left = element_text(angle = 90, face = "bold", size = 9),
        panel.spacing = unit(0.25, "lines"), plot.caption = element_text(size = 8, hjust = 0),
        plot.margin = margin(12, 6, 6, 6))
ggsave(file.path(OUT, "Figure_4G_tnk_prognosis_dot.png"), g, width = 5.2, height = 5.2, dpi = 200)
cat("wrote Figure_4G_tnk_prognosis_dot\n")
