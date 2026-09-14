#!/usr/bin/env Rscript
# Monocle3 pseudotime ordering of the 49 leukemic states + 8 normal compartments (reference anchors),
# median + IQR, ordered by median Monocle pseudotime; leukemic coloured by lineage program, normal in
# bold purple. Independent cross-check of the diffusion-pseudotime ordering. Data = monocle3_trajectory.R.
suppressMessages({library(ggplot2); library(dplyr)})
D3  <- "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/08_Source_Data/Figure_3/v3"
OUT <- "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/04_Main_Figures/Figure_3_panels"
d   <- read.csv(file.path(D3,"monocle_pseudotime.csv")) %>% filter(is.finite(monocle_pt))
ann <- read.csv(file.path(D3,"LS_annotation.csv"))[,c("LS","lineage_program")]
q <- d %>% group_by(group) %>% summarise(med=median(monocle_pt), lo=quantile(monocle_pt,.25), hi=quantile(monocle_pt,.75), .groups="drop")
q$normal <- grepl("^NORM_", q$group); q <- left_join(q, ann, by=c("group"="LS"))
q$lineage_program[q$normal] <- "normal"
q$lineage_program[is.na(q$lineage_program)|q$lineage_program==""] <- "unresolved"
q$label <- ifelse(q$normal, paste0("norm ", sub("NORM_","",q$group)),
                  paste0(sub("_","",q$group),"  ",ifelse(q$lineage_program=="unresolved","",q$lineage_program)))
q <- q %>% arrange(med); q$label <- factor(q$label, levels=q$label)
pal <- c("HSC/MPP"="#2166ac","GMP"="#4393c3","Monocyte"="#92c5de","Macrophage"="#0fcfc0","cDC"="#66c2a5",
         "pDC"="#3288bd","Megakaryocyte"="#8073ac","Erythroid"="#d6604d","Neutrophil"="#e08214",
         "B/plasma"="#5aae61","T/NK"="#c51b7d","proliferative"="#b2182b","unresolved"="#9aa0a6","normal"="#7b3294")
p <- ggplot(q, aes(med, label, colour=lineage_program)) +
  geom_segment(aes(x=lo,xend=hi,y=label,yend=label), linewidth=0.5) +
  geom_point(aes(shape=normal, size=normal)) +
  scale_colour_manual(values=pal, name="lineage program") +
  scale_shape_manual(values=c(`FALSE`=16,`TRUE`=18), guide="none") +
  scale_size_manual(values=c(`FALSE`=1.9,`TRUE`=3.4), guide="none") +
  labs(x="Monocle3 pseudotime (rooted at normal HSPC)", y=NULL) +
  theme_minimal(base_size=9) +
  theme(panel.grid.major.y=element_line(colour="grey93"), panel.grid.minor=element_blank(),
        axis.text.y=element_text(size=6.4), legend.position="right", legend.key.size=unit(0.32,"cm"))
ggsave(file.path(OUT,"Figure_3_monocle_ordering.pdf"), p, width=8, height=11)
ggsave(file.path(OUT,"Figure_3_monocle_ordering.png"), p, width=8, height=11, dpi=170)
cat("wrote Monocle3 ordering (", nrow(q), "groups )\n")
