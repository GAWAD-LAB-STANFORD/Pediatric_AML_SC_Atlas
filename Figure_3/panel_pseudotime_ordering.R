#!/usr/bin/env Rscript
# Pseudotime ordering: diffusion pseudotime (dpt, rooted at normal HSPC) of the 49 leukemic states
# and the 8 normal hematopoietic compartments (reference anchors), each shown as median + IQR and
# ordered by median dpt. Leukemic states coloured by lineage program; normal compartments in bold
# purple as reference. Tests where leukemic states sit relative to the normal differentiation axis.
# Data = compute_pseudotime_ordering.py + LS_annotation.csv. No fabricated values.
suppressMessages({library(ggplot2); library(dplyr)})
D3  <- "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/08_Source_Data/Figure_3/v3"
OUT <- "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/04_Main_Figures/Figure_3_panels"
dir.create(OUT, showWarnings=FALSE, recursive=TRUE)
d   <- read.csv(file.path(D3,"pseudotime_by_group.csv"), row.names=1)
ann <- read.csv(file.path(D3,"LS_annotation.csv"))[,c("LS","lineage_program")]
q <- d %>% group_by(group) %>% summarise(med=median(dpt), lo=quantile(dpt,.25), hi=quantile(dpt,.75), n=n(), .groups="drop")
q$normal <- grepl("^NORM_", q$group)
q <- left_join(q, ann, by=c("group"="LS"))
q$lineage_program[q$normal] <- "normal"
q$lineage_program[is.na(q$lineage_program) | q$lineage_program==""] <- "unresolved"
q$label <- ifelse(q$normal, paste0("norm ", sub("NORM_","",q$group)),
                  paste0(sub("_","",q$group), "  ", ifelse(q$lineage_program=="unresolved","",q$lineage_program)))
q <- q %>% arrange(med); q$label <- factor(q$label, levels=q$label)
pal <- c("HSC/MPP"="#2166ac","GMP"="#4393c3","Monocyte"="#92c5de","Macrophage"="#0fcfc0","cDC"="#66c2a5",
         "pDC"="#3288bd","Megakaryocyte"="#8073ac","Erythroid"="#d6604d","Neutrophil"="#e08214",
         "B/plasma"="#5aae61","T/NK"="#c51b7d","proliferative"="#b2182b","unresolved"="#9aa0a6","normal"="#7b3294")
p <- ggplot(q, aes(med, label, colour=lineage_program)) +
  geom_segment(aes(x=lo, xend=hi, y=label, yend=label), linewidth=0.5) +
  geom_point(aes(shape=normal, size=normal)) +
  scale_colour_manual(values=pal, name="lineage program") +
  scale_shape_manual(values=c(`FALSE`=16,`TRUE`=18), guide="none") +
  scale_size_manual(values=c(`FALSE`=1.9,`TRUE`=3.4), guide="none") +
  labs(x="diffusion pseudotime (dpt, rooted at normal HSPC)", y=NULL) +
  theme_minimal(base_size=9) +
  theme(panel.grid.major.y=element_line(colour="grey93"), panel.grid.minor=element_blank(),
        axis.text.y=element_text(size=6.4), legend.position="right", legend.key.size=unit(0.32,"cm"))
ggsave(file.path(OUT,"Figure_3_pseudotime_ordering.pdf"), p, width=8, height=11)
ggsave(file.path(OUT,"Figure_3_pseudotime_ordering.png"), p, width=8, height=11, dpi=170)
cat("wrote pseudotime ordering (", nrow(q), "groups )\n")
