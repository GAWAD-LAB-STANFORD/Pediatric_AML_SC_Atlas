#!/usr/bin/env Rscript
# CBF-RESTRICTED dot plot: enrichment computed among CBF patients ONLY.
# Panel = CBF cytogenetic (t(8;21) vs inv(16)); + stemness strip + within-CBF
# EFS Cox (HR, p) text column. Rows = states associated with within-CBF EFS
# (nominal Cox p<0.10), favorable/poor blocks, ordered by Cox p.
suppressMessages({library(ggplot2); library(patchwork); library(scales); library(survival)})
SL <- "/private/tmp/claude-503/-Users-chuckgawad-Desktop-ALSF-AML-2026-updated/f5e93d65-31d5-411e-b754-6f7c504928b7/scratchpad"
source(file.path(SL, "manuscript_palette.R"))
OUT <- "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/04_Main_Figures/Figure_2_new_panels"
PKG <- "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated"

d  <- read.csv(file.path(SL, "target_ls/LS48_survival.csv"))
ef <- read.csv(file.path(PKG,"__SUBMISSION_PACKAGE/08_Source_Data/Figure_2/Figure2E_TARGET_survival_input.csv"))[,c("sample","efs_time","efs_event")]
d  <- merge(d, ef, by="sample")
cbf <- d[d$subtype2 %in% c("t(8;21)","inv(16)") & !is.na(d$efs_time), ]

# --- within-CBF EFS Cox to select/annotate rows ---
res <- data.frame()
for (x in paste0("LS_",1:48)) { fr <- cbf[[x]]; if (sd(fr)==0 || sum(fr>0)<10) next
  s <- summary(coxph(Surv(efs_time,efs_event)~scale(fr)[,1], data=cbf))
  res <- rbind(res, data.frame(LS=x, HR=s$conf.int[1], p=s$coefficients[5])) }
res$dir <- ifelse(res$HR>1,"poor","favorable")
pr <- res[res$p < 0.10, ]
ann <- read.csv(file.path(SL, "ls48/LS48_annotation.csv"))[, c("LS", "lineage_program")]
stem<- read.csv(file.path(SL, "ls48/LS48_stem_scores.csv"))[, c("LS", "stem_z")]
pr <- merge(merge(pr, ann, by = "LS"), stem, by = "LS"); pr$LSc <- gsub("_", "", pr$LS)
pr$prog <- ifelse(pr$lineage_program == "unresolved", "", pr$lineage_program)
pr$lab  <- ifelse(pr$prog == "", pr$LSc, paste0(pr$LSc, "  ", pr$prog))
pr <- pr[order(pr$dir == "poor", pr$p), ]; pr$lab <- make.unique(pr$lab)
lv <- pr$lab; nfav <- sum(pr$dir == "favorable"); yline <- nrow(pr) - nfav + 0.5

# --- one-sided MWU cytogenetic enrichment within CBF (state higher in this subtype) ---
enrich <- function(dat, grpcol, levs){
  out <- data.frame()
  for (i in seq_len(nrow(pr))) { x <- pr$LS[i]; fr <- dat[[x]]; g <- dat[[grpcol]]
    for (L in levs) { a <- fr[g==L]; b <- fr[g!=L]
      p <- tryCatch(wilcox.test(a,b,alternative="greater")$p.value, error=function(e) 1)
      out <- rbind(out, data.frame(LS=x, lab=pr$lab[i], col=L, mean_pct=mean(a)*100, p=p)) } }
  out$fdr <- p.adjust(out$p,"BH"); out$neglog10fdr <- -log10(out$fdr)
  out$lab <- factor(out$lab, levels=rev(lv)); out }
dcy <- enrich(cbf, "subtype2", c("t(8;21)","inv(16)"))
dcy$col <- factor(dcy$col, levels=c("t(8;21)","inv(16)"))
szmax <- ceiling(max(dcy$mean_pct))

pc <- ggplot(dcy, aes(col, lab, size = mean_pct, fill = neglog10fdr)) +
  geom_hline(yintercept = yline, colour = "grey40", linewidth = 0.5) +
  geom_point(shape = 21, colour = "grey30", stroke = 0.3) +
  scale_fill_seq(limits = c(0,2.5), oob = squish, name = expression(-log[10]~FDR)) +
  scale_size(range = c(1,11), limits = c(0,szmax), name = "mean % (deconv.)") +
  labs(x = NULL, y = NULL) + ggtitle("CBF cytogenetic") + theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
        axis.text.y = element_text(size = 9), panel.grid.major = element_line(colour = "grey93"))

blk <- data.frame(y = c(nrow(pr)-nfav/2+0.5, (nrow(pr)-nfav)/2+0.5), lab = c("Favorable","Poor"))
pblk <- ggplot(blk, aes(1, y, label = lab, colour = lab)) +
  geom_text(angle = 90, fontface = "bold", size = 4.0) +
  scale_colour_manual(values = c(Favorable = unname(PROG["favorable"]), Poor = unname(PROG["poor"])), guide = "none") +
  scale_y_continuous(limits = c(0.5, nrow(pr)+0.5), expand = expansion(mult = 0.06)) +
  scale_x_continuous(limits = c(0.5, 1.5)) + coord_cartesian(clip = "off") + theme_void()

prs <- pr; prs$lab <- factor(prs$lab, levels = rev(lv))
ps <- ggplot(prs, aes("HSC/LSC", lab, fill = stem_z)) +
  geom_hline(yintercept = yline, colour = "grey40", linewidth = 0.5) +
  geom_tile(colour = "white", linewidth = 0.4) +
  scale_fill_div(midpoint=0, limits=c(-1,2), oob=squish, name="stemness\n(z)") +
  labs(x = NULL, y = NULL) + theme_minimal(base_size = 10) +
  theme(axis.text.y = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        panel.grid = element_blank())

# within-CBF EFS Cox as a text column (HR, p) — LS12 = rank 1
prs$txt <- sprintf("%.2f  (p=%s)", prs$HR, formatC(prs$p, format="fg", digits=2))
prs$best <- prs$p == min(prs$p)
ptxt <- ggplot(prs, aes(1, lab, label = txt, fontface = ifelse(best,"bold","plain"))) +
  geom_hline(yintercept = yline, colour = "grey40", linewidth = 0.5) +
  geom_text(size = 3.1, colour = "grey15") +
  scale_x_continuous(limits = c(0.4, 1.6)) + ggtitle("EFS Cox  HR (p)") +
  coord_cartesian(clip = "off") + theme_void() +
  theme(plot.title = element_text(size = 9, hjust = 0.5))

fig <- pblk + pc + ps + ptxt + plot_layout(widths = c(0.05, 0.5, 0.12, 0.34), guides = "collect") +
  plot_annotation(theme = theme(plot.margin = margin(10,8,6,8)))
ggsave(file.path(OUT, "Figure_2C_CBF_limited.pdf"), fig, width = 9.5, height = 5.3)
ggsave(file.path(OUT, "Figure_2C_CBF_limited.png"), fig, width = 9.5, height = 5.3, dpi = 200)
cat(sprintf("wrote CBF-limited dot plot (cyto + stemness + EFS Cox text): %d states; size max %.0f%%\n", nrow(pr), szmax))
