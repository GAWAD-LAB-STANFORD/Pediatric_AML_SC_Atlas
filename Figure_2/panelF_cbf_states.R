#!/usr/bin/env Rscript
# Figure 2F (v3): within-CBF EFS Cox for the leukemic states -- forest of HR (per SD of
# deconvolved fraction) for states passing FDR<0.10, coloured by direction, annotated
# with differentiation program + stemness. Companion to the CBF stem-program KM (2G).
# No fabricated values: HR/CI/p from cox on the corrected CIBERSORTx fractions.
suppressMessages({library(ggplot2)})
CODE <- "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/07_Code/Figure_2"
source(file.path(CODE, "manuscript_palette.R"))
D2  <- "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/08_Source_Data/Figure_2/v3"
D3  <- "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/08_Source_Data/Figure_3/v3"
OUT <- "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/04_Main_Figures/Figure_2_panels"

# recompute CI alongside HR (LS_cbf_cox_efs.csv lacks CI); use corrected fractions
suppressMessages(library(survival))
sv <- read.csv(file.path(D3,"LS_survival.csv"), check.names=FALSE)
cbf <- sv[sv$subtype2 %in% c("t(8;21)","inv(16)") & !is.na(sv$efs_time), ]
ann <- read.csv(file.path(D3,"LS_annotation.csv"))
lsN <- grep("^LS_", names(sv), value=TRUE)
r <- data.frame()
for (x in lsN){ fr <- cbf[[x]]; if (sd(fr)==0 || sum(fr>0)<10) next
  s <- summary(coxph(Surv(efs_time,efs_event)~scale(fr), data=cbf))
  r <- rbind(r, data.frame(LS=x, HR=s$conf.int[1], lo=s$conf.int[3], hi=s$conf.int[4], p=s$coefficients[5])) }
r$fdr <- p.adjust(r$p,"BH"); r$dir <- ifelse(r$HR>1,"poor","favorable")
r <- merge(r, ann[,c("LS","lineage_program","stem_z","n_patients")], by="LS")
sig <- r[r$fdr < 0.10, ]
sig$star <- ifelse(sig$n_patients<=1,"*","")
sig$progn <- ifelse(sig$lineage_program=="unresolved","", sig$lineage_program)
sig$lab <- ifelse(sig$progn=="", paste0(gsub("_","",sig$LS),sig$star),
                  paste0(gsub("_","",sig$LS),sig$star,"  ",sig$progn))
sig <- sig[order(sig$HR), ]; sig$lab <- factor(sig$lab, levels=sig$lab)

g <- ggplot(sig, aes(HR, lab, colour=dir)) +
  geom_vline(xintercept=1, linetype=2, colour="grey55") +
  geom_errorbarh(aes(xmin=lo, xmax=hi), height=0.22, linewidth=0.7) +
  geom_point(size=3.4) +
  geom_text(aes(label=sprintf("p=%.3f", p)), hjust=-0.15, vjust=-0.9, size=3, colour="grey25") +
  scale_colour_manual(values=c(favorable=unname(PROG["favorable"]), poor=unname(PROG["poor"])), guide="none") +
  scale_x_continuous(trans="log2", breaks=c(0.6,0.8,1.0,1.25,1.5)) +
  labs(x="within-CBF EFS hazard ratio (per SD)", y=NULL) +  # "* = single-patient state" defined in the figure legend
  theme_classic(base_size=11) +
  theme(axis.text.y=element_text(size=10), plot.caption=element_text(size=8, hjust=0))
ggsave(file.path(OUT,"Figure_2F_cbf_states.pdf"), g, width=5.6, height=3.6)
ggsave(file.path(OUT,"Figure_2F_cbf_states.png"), g, width=5.6, height=3.6, dpi=200)
cat(sprintf("wrote Figure_2F_cbf_states: %d states FDR<0.10 within CBF (n=%d)\n", nrow(sig), nrow(cbf)))
print(sig[,c("LS","HR","lo","hi","p","fdr","dir","lineage_program")], row.names=FALSE, digits=3)
