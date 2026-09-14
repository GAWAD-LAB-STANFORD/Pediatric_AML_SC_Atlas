#!/usr/bin/env Rscript
# Panel 3G: program map across normal HSPC, leukemic HSC/MPP states (LS5/LS38), LS10, and
# other leukemic. Four modules -- HSC quiescence (HLF/AVP/CRHBP), leukemic stemness
# (HOXA9/MEIS1), proliferation (MKI67/TOP2A/CDK1), myeloid differentiation (MPO/LYZ/CD14).
# Dot size = % expressing, fill = mean log-CPM. Data = compute_program_modules.py. No fabrication.
suppressMessages({library(ggplot2); library(scales)})
CODE <- "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/07_Code/Figure_3"
source(file.path(CODE, "manuscript_palette.R"))
D3 <- "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/08_Source_Data/Figure_3/v3"
OUT <- "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/04_Main_Figures/Figure_2_panels"
d <- read.csv(file.path(D3, "LS_program_modules.csv"))
d$group  <- factor(d$group, levels = rev(c("normal HSPC","LS5/LS38 (leuk. HSC/MPP)","LS10 (proliferative)","other leukemic")))
d$module <- factor(d$module, levels = c("HSC quiescence","leukemic stemness","proliferation","differentiation"))
d$gene   <- factor(d$gene, levels = c("HLF","AVP","CRHBP","HOXA9","MEIS1","MKI67","TOP2A","CDK1","MPO","LYZ","CD14"))
g <- ggplot(d, aes(gene, group, size = pct, fill = mean_expr)) +
  geom_point(shape = 21, colour = "grey30", stroke = 0.3) +
  facet_grid(~module, scales = "free_x", space = "free_x") +
  scale_size(range = c(0.3, 9), limits = c(0, 90), name = "% expressing") +
  scale_fill_seq(name = "mean\nlog-CPM") +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(face = "bold", size = 9), axis.text.y = element_text(size = 9.5),
        strip.text = element_text(face = "bold", size = 9), panel.grid.major = element_line(colour = "grey92"),
        plot.caption = element_text(size = 7.4, hjust = 0))
ggsave(file.path(OUT, "Figure_3_quiescence_program.pdf"), g, width = 8.4, height = 3.3)
ggsave(file.path(OUT, "Figure_3_quiescence_program.png"), g, width = 8.4, height = 3.3, dpi = 200)
cat("wrote Figure_3_quiescence_program (4 modules)\n")
