#!/usr/bin/env Rscript
# S4: bulk-deconvolution validation for the corrected 49-state scheme (deployed CIBERSORTx
# deconvolution of TARGET bulk, Jobs 24-27). (A) per-sample goodness-of-fit across all 1,787
# deconvolved samples (CIBERSORTx Pearson correlation and RMSE; P<0.05 in 99.8%). (B) per-49-LS
# detection reliability (fraction of samples with deconvolved fraction >1% and >5%), states
# ordered by detectability and coloured by whole-cohort EFS prognosis; states left of the dashed
# line are poorly resolved in bulk. Data = compute_deconv_qc.py. Replaces the legacy figure.
suppressMessages({library(ggplot2); library(patchwork)})
D2  <- "../source_data"
OUT <- "../figures"
qc  <- read.csv(file.path(D2,"deconv_fit_qc.csv"))
det <- read.csv(file.path(D2,"LS_detection_49.csv"))
pr  <- read.csv(file.path(D2,"LS_wholecohort_prognosis.csv"))[,c("LS","prognosis")]
pcol <- c(favorable="#00798c", poor="#d1495b", ns="#9aa0a6")
fracP <- mean(qc$P.value < 0.05, na.rm=TRUE)

pA <- ggplot(qc, aes(Correlation)) +
  geom_histogram(bins=40, fill="#00798c", colour="white", linewidth=0.15) +
  geom_vline(xintercept=median(qc$Correlation,na.rm=TRUE), linetype="dashed", colour="grey30") +
  annotate("text", x=median(qc$Correlation,na.rm=TRUE), y=Inf, vjust=1.5, hjust=-0.05,
           label=sprintf("median %.2f", median(qc$Correlation,na.rm=TRUE)), size=3.2, colour="grey30") +
  annotate("text", x=Inf, y=Inf, hjust=1.05, vjust=2.6,
           label=sprintf("n=%d\nP<0.05 in %.1f%%", nrow(qc), 100*fracP), size=3, colour="grey30") +
  labs(x="per-sample fit correlation (CIBERSORTx)", y="TARGET samples",
       title="A  Deconvolution goodness-of-fit") +
  theme_minimal(base_size=11) + theme(panel.grid.minor=element_blank(), plot.title=element_text(size=11,face="bold"))
pB <- ggplot(qc, aes(RMSE)) +
  geom_histogram(bins=40, fill="#8d99ae", colour="white", linewidth=0.15) +
  geom_vline(xintercept=median(qc$RMSE,na.rm=TRUE), linetype="dashed", colour="grey30") +
  annotate("text", x=median(qc$RMSE,na.rm=TRUE), y=Inf, vjust=1.5, hjust=-0.05,
           label=sprintf("median %.2f", median(qc$RMSE,na.rm=TRUE)), size=3.2, colour="grey30") +
  labs(x="per-sample RMSE (CIBERSORTx)", y="TARGET samples", title="B  Deconvolution residual") +
  theme_minimal(base_size=11) + theme(panel.grid.minor=element_blank(), plot.title=element_text(size=11,face="bold"))

det <- merge(det, pr, by="LS", all.x=TRUE); det$prognosis[is.na(det$prognosis)] <- "ns"
det$prognosis <- factor(det$prognosis, levels=names(pcol))
det <- det[order(det$det_gt1pct),]; det$LS <- factor(det$LS, levels=det$LS)
pC <- ggplot(det, aes(det_gt1pct, LS)) +
  geom_vline(xintercept=0.20, linetype="dashed", colour="grey60") +
  geom_segment(aes(x=0, xend=det_gt1pct, y=LS, yend=LS), colour="grey80", linewidth=0.4) +
  geom_point(aes(colour=prognosis), size=2.4) +
  geom_point(aes(x=det_gt5pct), colour="grey40", size=1, shape=1) +
  scale_colour_manual(values=pcol, name="EFS prognosis") +
  scale_x_continuous(limits=c(0,NA)) +
  labs(x="fraction of TARGET samples detected (filled >1%, open >5%)", y=NULL,
       title="C  Per-state detection reliability in bulk (49 leukemic states)") +
  theme_minimal(base_size=10) +
  theme(panel.grid.minor=element_blank(), axis.text.y=element_text(size=6),
        plot.title=element_text(size=11,face="bold"), legend.position="right")

fig <- (pA | pB) / pC + plot_layout(heights=c(1,2.1))
ggsave(file.path(OUT,"Figure_S4__bulk_deconvolution.pdf"), fig, width=10, height=11)
ggsave(file.path(OUT,"Figure_S4__bulk_deconvolution.png"), fig, width=10, height=11, dpi=200)
cat("wrote S4 deconvolution validation (49 states)\n")
