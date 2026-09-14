# ====================================================================
# figS8.R  |  Figure S22 (CD96 CellChat cell-cell communication)
# reproducibility header (auto-generated; edit freely below this block)
# Inputs : cd96_autocrine_coexpr.csv, cd96_persample_long.csv, cd96_sample_info.csv, figS22EF_cellchat.csv
# Outputs: (writes to ../figures or ../source_data — see body)
# Run    : Rscript figS8.R   (or: bash run_all.sh)  |  Deps: see REPRODUCIBILITY_CHECKLIST.md
# Source data for this figure is in __SUBMISSION_PACKAGE/08_Source_Data/.
# ====================================================================

# Figure S8 (R/ggplot2) — CD96 cell–cell communication, supplement to Fig 7E/7F.
# A  full sender→receiver network (all senders incl. the T/NK→AML checkpoint axis)
# B  autocrine substrate: CD96 + PVR/NECTIN1 co-expression on the leukemic cells
# C  per-patient CD96-on-AML signalling (heterogeneity behind the aggregate)
# D  CD96 autocrine signalling vs blast CD96 positivity (targetable vs cold vs healthy)
# All from cellchat_CD96/cellchat_CD96.R (per-sample CellChat, CD96 pathway).
source("_style.R")
rord <- c("AML (leukemic)", "Erythroid", "Myeloid_Pro", "Mature myeloid",
          "HSPC", "B-lineage", "T cells", "NK")

## ---- A: full sender → receiver network (19 CD96-targetable samples) ----
agg <- sd("figS22EF_cellchat.csv") |>
  dplyr::mutate(lr = sub("CD96_", "CD96–", interaction),
                source = factor(source, levels = rord),
                target = factor(target, levels = rev(rord)))
pA <- ggplot(agg, aes(source, target, size = n_sig, fill = mean_prob)) +
  geom_point(shape = 21, colour = "grey30", stroke = 0.3) +
  facet_wrap(~ lr) +
  scale_fill_gradientn(colours = SEQ, name = "mean\ncomm prob", trans = "sqrt") +
  scale_size_area(max_size = 7, name = "# patients\nsignificant", breaks = c(1, 4, 8, 12)) +
  labs(x = "sender (CD96+)", y = "receiver (PVR / NECTIN1+)") +
  theme_pub(8.5) + ggtitle("All CD96 Senders → Receivers (19 CD96+ Samples)") +
  theme(axis.text.x = element_text(angle = 40, hjust = 1),
        panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
        strip.text = element_text(face = "bold"))

## ---- B: autocrine substrate (co-expression on leukemic cells) ----
acl <- sd("cd96_autocrine_coexpr.csv") |>
  tidyr::pivot_longer(c(pct_CD96, pct_PVR, pct_NECTIN1, pct_CD96_PVR, pct_CD96_NECTIN1),
                      names_to = "cat", values_to = "pct") |>
  dplyr::mutate(cat = dplyr::recode(cat,
      pct_CD96 = "CD96\n(ligand)", pct_PVR = "PVR", pct_NECTIN1 = "NECTIN1",
      pct_CD96_PVR = "CD96 &\nPVR", pct_CD96_NECTIN1 = "CD96 &\nNECTIN1"),
    cat = factor(cat, levels = c("CD96\n(ligand)", "PVR", "NECTIN1",
                                 "CD96 &\nPVR", "CD96 &\nNECTIN1")))
pB <- ggplot(acl, aes(cat, pct)) +
  geom_boxplot(outlier.shape = NA, fill = "grey93", width = 0.6, linewidth = 0.3) +
  geom_jitter(width = 0.12, height = 0, size = 1.4, alpha = 0.55, colour = PAL$CD96) +
  labs(x = NULL, y = "% of leukemic cells positive") +
  theme_pub(9) + ggtitle("Autocrine Substrate On Leukemic Cells (19 CD96+ Samples)")

## ---- C: per-patient CD96-on-AML signalling (heterogeneity) ----
ps <- sd("cd96_persample_long.csv")
aml <- ps |> dplyr::filter(targetable, source == "AML (leukemic)") |>
  dplyr::mutate(lr = sub("CD96_", "CD96–", interaction),
                target = factor(target, levels = rev(rord)))
ord <- aml |> dplyr::filter(target == "AML (leukemic)") |>
  dplyr::group_by(sample) |> dplyr::summarise(p = sum(prob), .groups = "drop") |>
  dplyr::arrange(p) |> dplyr::pull(sample)
aml$sample <- factor(aml$sample, levels = ord)
aml$sig <- ifelse(!is.na(aml$pval) & aml$pval <= 0.05, "*", "")
pC <- ggplot(aml, aes(target, sample, fill = prob)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = sig), size = 3.2, vjust = 0.78) +
  facet_wrap(~ lr) +
  scale_fill_gradientn(colours = SEQ, name = "comm\nprob", trans = "sqrt") +
  labs(x = "receiver", y = "patient (CD96-targetable)") +
  theme_pub(8) + ggtitle("Per-Patient CD96-on-AML Signalling  (* p ≤ 0.05)") +
  theme(axis.text.x = element_text(angle = 40, hjust = 1),
        axis.text.y = element_text(size = 6),
        strip.text = element_text(face = "bold"), panel.grid = element_blank())

## ---- D: autocrine signalling vs blast CD96 positivity ----
si <- sd("cd96_sample_info.csv")
auto <- ps |> dplyr::filter(source == "AML (leukemic)", target == "AML (leukemic)") |>
  dplyr::group_by(sample) |>
  dplyr::summarise(prob = sum(prob), sig = any(pval <= 0.05, na.rm = TRUE), .groups = "drop")
dd <- si |> dplyr::left_join(auto, by = "sample") |>
  dplyr::mutate(prob = ifelse(is.na(prob), 0, prob), sig = ifelse(is.na(sig), FALSE, sig),
    Significance = factor(ifelse(sig, "p ≤ 0.05", "n.s."), levels = c("p ≤ 0.05", "n.s.")),
    grp = ifelse(grepl("Healthy", sample), "Healthy BM",
          ifelse(targetable, "CD96-targetable AML", "CD96-cold AML")),
    grp = factor(grp, levels = c("CD96-targetable AML", "CD96-cold AML", "Healthy BM")))
pD <- ggplot(dd, aes(pct_cd96_blasts, prob + 1e-7)) +
  geom_vline(xintercept = 10, linetype = "dashed", colour = "grey50") +
  annotate("text", x = 10.5, y = 2e-3, label = "targetable ≥10%", hjust = 0,
           size = 2.6, colour = "grey40") +
  geom_point(aes(fill = grp, colour = Significance), shape = 21, size = 3.2, stroke = 0.8) +
  scale_y_log10() +
  scale_fill_manual(values = c("CD96-targetable AML" = PAL$CD96,
      "CD96-cold AML" = "grey70", "Healthy BM" = "#2471a3"), name = NULL) +
  scale_colour_manual(values = c("p ≤ 0.05" = "black", "n.s." = "grey80"),
      name = "CellChat permutation") +
  guides(fill = guide_legend(order = 1, override.aes = list(colour = "grey60")),
         colour = guide_legend(order = 2, override.aes = list(fill = "grey85"))) +
  labs(x = "% CD96+ leukemic blasts (per sample)",
       y = "AML→AML signalling strength\n(CellChat comm. prob., +1e-7, log)") +
  theme_pub(9) + ggtitle("CD96 Autocrine Signalling Tracks Blast CD96 Positivity") +
  theme(legend.position = "top", legend.box = "vertical")

## ---- E/F: CD96–nectin checkpoint axis (both directions; CD96 is aberrant on AML AND on T/NK) ----
## Data: build_cd96_axis.py -> cd96_axis_expr.csv (compartment x gene: %pos, mean) +
## cd96_axis_connectivity.csv (AML<->T/NK co-expression connectivity for the CD96-nectin pairs).
e  <- sd("cd96_axis_expr.csv")
cc <- sd("cd96_axis_connectivity.csv")
.COMP <- c("AML (leukemic)", "Naïve CD4 T", "Memory CD4 T", "Treg", "Naïve CD8 T", "Memory CD8 T",
           "GZMK CD8 T", "GZMB CD8 T", "GZMB DNT", "MAIT", "GZMK NK", "GZMB NK", "Proliferating T")
e$compartment <- factor(e$compartment, levels = rev(.COMP))
e$gene <- factor(e$gene, levels = c("CD96", "TIGIT", "CD226", "PVR", "NECTIN1", "NECTIN2"))
e$role <- factor(ifelse(e$role == "receptor", "receptor (on immune / AML)", "nectin ligand"),
                 levels = c("receptor (on immune / AML)", "nectin ligand"))
pE <- ggplot(e, aes(gene, compartment)) +
  geom_point(aes(size = pct, fill = mean), shape = 21, colour = "#3b5870", stroke = 0.3) +
  scale_fill_gradientn(colours = SEQ, name = "mean expr\n(log1p)") +
  scale_size_area(max_size = 6, limits = c(0, 55), breaks = c(5, 20, 40), name = "% positive") +
  facet_grid(~ role, scales = "free_x", space = "free") +
  labs(x = NULL, y = NULL) + ggtitle("CD96–nectin axis: receptors and ligands") +
  theme_pub(9) + theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 8),
    strip.text = element_text(face = "bold", size = 8),
    strip.background = element_rect(fill = "#eef2f6", colour = NA),
    panel.spacing = grid::unit(4, "pt"))
cc$subset <- factor(cc$subset, levels = rev(.COMP[-1]))
pF <- ggplot(cc, aes(score, subset, fill = direction)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.66) +
  scale_fill_manual(values = c("AML → T/NK (nectin → CD96)" = unname(PAL$CD96),
    "T/NK → AML (nectin → CD96)" = "#95a5a6"), name = NULL) +
  labs(x = "co-expression connectivity (% ligand+ × % CD96+ ÷ 100)", y = NULL) +
  ggtitle("Directional CD96-axis signalling") +
  theme_pub(9) + theme(legend.position = c(0.98, 0.5), legend.justification = c(1, 0.5),
    legend.background = element_rect(fill = scales::alpha("white", 0.75), colour = NA))

figS8 <- (pA | pB) / (pC | pD) / (pE | pF) +
  plot_layout(heights = c(1, 1.25, 1.05)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 16))
savefig(figS8, "supplementary/FigureS22_cd96_cellchat", 14, 18)
