# ====================================================================
# fig5.R  |  Figure 5 (CD96 target discovery, panels A-G)
# reproducibility header (auto-generated; edit freely below this block)
# Inputs : fig5A_funnel.csv, fig5A_top15.csv, fig5C_headtohead.csv, fig5D_cyto.csv, fig5D_ls.csv, fig5F_adult_beataml_meta.csv ...
# Outputs: (writes to ../figures or ../source_data — see body)
# Run    : Rscript fig5.R   (or: bash run_all.sh)  |  Deps: see REPRODUCIBILITY_CHECKLIST.md
# Source data for this figure is in __SUBMISSION_PACKAGE/08_Source_Data/.
# ====================================================================

# Figure 5 (R/ggplot2) -- CD96 is the lead surface immunotherapy target.
#   A funnel | B CD96 UMAP | C head-to-head | D cyto+LS dotplot | E coverage summary | F/G bulk validation
source("_style.R")
suppressPackageStartupMessages(library(png))

## ---- A: discovery funnel + top-15 composite ----
fn <- sd("fig5A_funnel.csv")
fn$stage <- factor(fn$stage, levels = rev(fn$stage))   # file order; broadest on top
fn$col <- colorRampPalette(c("#94D2BD", "#0A9396", "#005F73", "#001219"))(nrow(fn))  # heatmap0 teal ramp (Fig 3 scheme)
pA1 <- ggplot(fn, aes(stage, n)) +
  geom_col(aes(fill = col), width = 0.74) +
  scale_fill_identity() +
  geom_text(aes(label = scales::comma(n)), hjust = -0.12, size = 2.4) +
  coord_flip() + scale_y_continuous(expand = expansion(c(0, 0.32))) +
  labs(x = NULL, y = "surface genes") + theme_pub(9)   # commentary -> caption
t15 <- sd("fig5A_top15.csv") |> dplyr::slice(1:15) |>
  dplyr::mutate(lead = gene %in% LEAD_TARGETS,
                col = dplyr::if_else(lead, as.character(LEAD_COLORS[gene]), PAL$grey),  # leads pop; others neutral grey
                gene = forcats::fct_reorder(gene, composite))
pA2 <- ggplot(t15, aes(composite, gene)) +
  geom_segment(aes(x = 0, xend = composite, yend = gene, colour = col), linewidth = 0.7) +
  geom_point(aes(colour = col, size = lead)) +
  scale_colour_identity(guide = "none") +
  scale_size_manual(values = c(`TRUE` = 3.1, `FALSE` = 2.1), guide = "none") +
  scale_y_discrete(labels = function(v) ifelse(v %in% LEAD_TARGETS,
                   sprintf("<b style='color:%s'>%s</b>", LEAD_COLORS[v], v), v)) +
  labs(x = "composite score", y = NULL) + theme_pub(9) +
  theme(axis.text.y = ggtext::element_markdown())      # commentary -> caption
pA <- wrap_elements(pA1 + pA2 + plot_layout(widths = c(0.8, 1)))

## ---- B: CD96 UMAP -- rendered with scanpy (Python), embedded as an image ----
## (UMAP embeddings come from scanpy; per the agreed scanpy exception this panel
## is the scanpy/matplotlib render rather than a ggplot scatter.)
pB <- patchwork::wrap_elements(full = grid::rasterGrob(
  png::readPNG(file.path("..", "assets", "fig6B_leads_umap.png")),
  interpolate = TRUE))

## ---- C: head-to-head — composite candidates vs clinically-pursued AML antigens + FLT3 (+ctrl) ----
# x = normal HSPC/Myeloid leak (toxicity), y = % leukemic positive (efficacy), size = % patients >=20%.
# klass: lead (5 leads, coloured circles) | candidate (gated, grey circles) |
#        clinical (top-6 AML-trial antigens, slate diamonds) | control (FLT3, gold diamond).
hh <- sd("fig5C_headtohead.csv") |>
  # corrected scheme: the retired PAC/PPAC "subset-specific" label is removed. CD7 (a bona fide clinical
  # CAR-T antigen) folds into the clinical comparators; the other former-PPAC markers become additional
  # candidate antigens (grey), matching panel D and the 12-antigen comparison used in Figures 6-7.
  dplyr::mutate(klass = dplyr::case_when(label == "CD7" ~ "clinical",
                                         klass == "ppac" ~ "candidate", TRUE ~ klass)) |>
  dplyr::mutate(
    fill = dplyr::case_when(klass == "lead" ~ as.character(LEAD_COLORS[gene]),
                            klass == "control" ~ "#f1c40f",
                            klass == "clinical" ~ "#5d6d7e",
                            TRUE ~ "#c7ced6"),
    labcol = dplyr::case_when(klass == "lead" ~ as.character(LEAD_COLORS[gene]),
                              klass == "control" ~ "#9a7d0a",
                              klass == "clinical" ~ "#2c3e50",
                              TRUE ~ "#5d6d7e"),
    shp = dplyr::case_when(klass %in% c("clinical", "control") ~ 23, TRUE ~ 21),
    disp = ifelse(do_label, label, ""),
    ptalpha = ifelse(klass == "candidate" & !do_label, 0.5, 1))
pC <- ggplot(hh, aes(leak, scRNA_AML_pct)) +
  annotate("rect", xmin = -2, xmax = 10, ymin = -Inf, ymax = Inf, fill = "#eafaf1", alpha = 0.55) +
  geom_point(aes(size = n_pts_pct, fill = fill, shape = shp, alpha = ptalpha),
             colour = "grey25", stroke = 0.35) +
  scale_fill_identity() + scale_shape_identity() + scale_alpha_identity() +
  scale_size_area(max_size = 8, name = "% patients >=20%+", breaks = c(10, 25, 50, 75)) +
  ggrepel::geom_text_repel(aes(label = disp, colour = labcol), size = 2.7, fontface = "bold",
                           max.overlaps = 40, seed = 1, box.padding = 0.5, point.padding = 0.3,
                           min.segment.length = 0, segment.size = 0.2, segment.colour = "grey70") +
  scale_colour_identity() +
  scale_x_continuous(expand = expansion(mult = c(0.04, 0.16))) +
  labs(x = "% normal HSPC / Myeloid_Pro positive (max)  (toxicity)",
       y = "% leukemic cells positive  (efficacy)") +
  theme_pub(9)                                         # commentary -> caption

## ---- D: coverage dotplot, single-cell THROUGHOUT (one shared colour + size scale). LEFT = cytogenetic
##        subtypes (leukemic cells; MPAL + single-patient subtypes flagged); RIGHT = all 49 leukemic states
##        (favorable -> poor by Cox HR) + normal HSPC/Myeloid baseline. TARGET-bulk validation is shown on
##        the selected targets in panels E/F, so panel D stays one consistent modality. CD96 = crimson. ----
suppressPackageStartupMessages(library(ggdendro))
cyto <- sd("fig5D_cyto.csv"); ls5 <- sd("fig5D_ls.csv")
gset <- c("CD96", "CD9", "SUCNR1", "AMN", "IL1RAP", "LMBR1L", "NRXN2", "FURIN", "UMODL1", "CD7",
          "CD33", "CD123", "CLL-1", "CD38", "CD70", "FLT3", "TNFRSF4", "ABCA7", "ITGAX")
CLIN_TRIAL <- c("CD33", "CD123", "CLL-1", "CD38", "CD70", "CD7", "IL1RAP", "FLT3")   # in AML trials (+FLT3 ctrl)
labcolD <- function(v) ifelse(v == "CD96", "#c0392b",                     # CD96 = crimson hero
                       ifelse(v %in% CLIN_TRIAL, "#e8590c", "#2c3e50"))    # clinical = bright orange; else slate
dispD <- function(v) sprintf("<b style='color:%s'>%s</b>", labcolD(v), v)
# cluster genes by their scRNA prognostic-LS coverage profile (the efficacy readout)
matof <- function(df) { w <- df |> dplyr::select(label, group, pct) |>
    tidyr::pivot_wider(names_from = group, values_from = pct)
  m <- as.matrix(w[, -1]); rownames(m) <- w$label; m[gset[gset %in% rownames(m)], , drop = FALSE] }
lsm <- matof(ls5); go <- rownames(lsm)
hg <- hclust(dist(lsm), "ward.D2"); go <- rownames(lsm)[hg$order]         # cluster GENES
# single-cell cytogenetic subtypes ordered by biology; MPAL (mixed-phenotype) = teal, single-patient = grey
cytoord <- c("RUNX1/RUNX1T1", "CBFB/MYH11", "PML/RARA", "MLLr", "Tri(8)/MLLr", "NUP98/NSD1", "del7q",
             "Tri(8)", "Tri(15)", "CN", "MYB/GATA1", "BCR/ABL", "t(2;3)(p15;q26.2)", "t(7;14)(q21;q32)")
cytoord <- cytoord[cytoord %in% unique(cyto$group)]
cyto_np   <- cyto |> dplyr::distinct(group, n_pt, mpal)
cyto_npt  <- setNames(cyto_np$n_pt, cyto_np$group); cyto_mpal <- setNames(cyto_np$mpal, cyto_np$group)
cyto_xcol <- function(v) sprintf("<b style='color:%s'>%s</b>",
  ifelse(cyto_mpal[v] == 1, "#0E7C7B", ifelse(cyto_npt[v] <= 1, "#9aa0a6", "#2c3e50")), v)  # MPAL teal / 1-pt grey
# ALL 49 leukemic states, ordered favorable -> poor by whole-cohort Cox HR (the builder's file order),
# bracketed by the normal HSPC/Myeloid baseline. x-label colour = significant favorable/poor (FDR<0.05) or n.s.
lspm   <- ls5 |> dplyr::distinct(group, prognosis)
lsprog <- setNames(lspm$prognosis, lspm$group)                           # group -> prognosis (favorable/poor/n.s./normal)
lsord  <- unique(ls5$group)                                              # HSPC, Myeloid_Pro, then 49 LS by HR
n_hr_lo <- ls5 |> dplyr::filter(grepl("^LS_", group)) |> dplyr::distinct(group, hr) |>
  dplyr::summarise(n = sum(hr < 1)) |> dplyr::pull(n)                    # LS on the favorable (HR<1) side
PROGCOL <- c(normal = "#7b3294", favorable = "#0A9396", poor = "#AE2012", `n.s.` = "#9aa0a6")
ls_xcol <- function(v) sprintf("<b style='color:%s'>%s</b>", PROGCOL[lsprog[v]], v)
dseg <- ggdendro::dendro_data(as.dendrogram(hg))$segments
gdend <- ggplot(dseg) +
  geom_segment(aes(x = -y, y = x, xend = -yend, yend = xend), linewidth = 0.3, colour = "#95a5a6") +
  scale_y_continuous(limits = c(0.4, length(go) + 0.6), expand = c(0, 0)) +
  scale_x_continuous(expand = expansion(mult = c(0.06, 0))) + theme_void()
# Both sub-panels are single-cell (leukemic cells), so they SHARE one colour scale (mean log-norm) and
# one size scale (% cells positive) — collected into a single legend. Dot SIZE = % cells positive,
# dot COLOUR = mean log-norm expression.
DMAX <- ceiling(max(c(cyto$mean_expr, ls5$mean_expr)) * 10) / 10
SMAX <- ceiling(max(c(cyto$pct, ls5$pct)) / 10) * 10
dotD <- function(df, xlevels, ttl, ylabs = TRUE, xcol = NULL) {
  df$gy <- factor(df$label, levels = go); df$group <- factor(df$group, levels = xlevels)
  df <- df[!is.na(df$gy) & !is.na(df$group), ]
  xlab <- if (is.null(xcol)) ggplot2::waiver() else function(v) xcol(v)
  xtheme <- if (is.null(xcol)) element_text(angle = 45, hjust = 1, size = 6)
            else ggtext::element_markdown(angle = 45, hjust = 1, size = 6)
  ggplot(df, aes(group, gy)) +
    geom_point(aes(size = pct, fill = mean_expr), shape = 21, colour = "#3b5870", stroke = 0.2) +
    scale_fill_gradientn(colours = SEQ, oob = scales::squish, limits = c(0, DMAX), name = "mean\nlog-norm") +
    scale_size_area(max_size = 4, limits = c(0, SMAX), name = "% cells\npositive", breaks = c(10, 25, 50, 75)) +
    scale_y_discrete(limits = go, labels = function(v) dispD(v), expand = expansion(add = 0.6)) +
    scale_x_discrete(expand = expansion(add = 0.6), labels = xlab) +
    labs(x = NULL, y = NULL, title = ttl) + theme_pub(8.5) +
    theme(axis.text.x = xtheme, plot.title = element_text(size = 8.5, face = "bold"),
          axis.text.y = if (ylabs) ggtext::element_markdown(size = 7.3) else element_blank())
}
pD1 <- dotD(dplyr::filter(cyto, label %in% go), cytoord, "Cytogenetic subtype (scRNA)",
            ylabs = TRUE, xcol = cyto_xcol)
pD2 <- dotD(dplyr::filter(ls5, label %in% go), lsord, "Leukemic states — favorable → poor by Cox HR (scRNA)",
            ylabs = FALSE, xcol = ls_xcol) +
  geom_vline(xintercept = 2.5, colour = "grey55", linewidth = 0.5, linetype = "dashed") +           # normal | leukemic
  geom_vline(xintercept = 2.5 + n_hr_lo, colour = "grey70", linewidth = 0.5)                        # HR<1 | HR>1 (favorable | poor)
pD <- wrap_elements(gdend + pD1 + pD2 + plot_layout(widths = c(0.04, 0.42, 1.55), guides = "collect"))

## ---- E: summary of panel D — % cells positive per antigen across ALL leukemic states (total, grey),
##        and the favorable (teal) / poor (red) prognostic subsets, antigens ordered by MEDIAN poor-
##        prognosis coverage. Dashed line = 25% targetable threshold. Single-cell; bulk in F/G. ----
lsb <- ls5 |> dplyr::filter(grepl("^LS_", group))
ord_box <- lsb |> dplyr::filter(prognosis == "poor") |> dplyr::group_by(label) |>
  dplyr::summarise(m = median(pct), .groups = "drop") |> dplyr::arrange(dplyr::desc(m)) |> dplyr::pull(label)
bxd <- dplyr::bind_rows(lsb |> dplyr::mutate(cls = "total"),
                        lsb |> dplyr::filter(prognosis %in% c("favorable", "poor")) |>
                          dplyr::mutate(cls = as.character(prognosis)))
bxd$cls <- factor(bxd$cls, levels = c("total", "favorable", "poor"))
bxd$label <- factor(bxd$label, levels = ord_box)
hscd <- ls5 |> dplyr::filter(group == "HSPC"); hscd$label <- factor(hscd$label, levels = ord_box)   # normal HSC reference
BOXCOL <- c(total = "#808080", favorable = unname(PROG["favorable"]), poor = unname(PROG["poor"]))
pBox <- ggplot(bxd, aes(label, pct, fill = cls)) +
  geom_hline(yintercept = 25, linetype = "dashed", colour = "grey55", linewidth = 0.4) +
  geom_boxplot(outlier.shape = NA, width = 0.72, alpha = 0.55, linewidth = 0.3,
               position = position_dodge(0.8), colour = "grey35") +
  geom_point(aes(colour = cls), position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8),
             size = 0.7, alpha = 0.6) +
  geom_point(data = hscd, aes(label, pct), shape = 23, fill = "#7b3294", colour = "grey20",
             size = 2, stroke = 0.35, inherit.aes = FALSE,
             position = position_nudge(x = -0.8 / 3)) +                            # normal HSC diamond, aligned to the "total" box
  scale_fill_manual(values = BOXCOL, name = "leukemic states") + scale_colour_manual(values = BOXCOL, guide = "none") +
  scale_x_discrete(labels = function(v) dispD(v)) +
  labs(x = NULL, y = "% cells positive per leukemic state",
       title = "Coverage across leukemic states (single-cell), antigens ordered by median poor-prognosis coverage",
       caption = "purple diamond = normal HSC (toxicity reference)") +
  theme_pub(9) + theme(axis.text.x = ggtext::element_markdown(angle = 40, hjust = 1),
                       legend.position = c(0.995, 0.98), legend.justification = c(1, 1),
                       legend.background = element_rect(fill = "white", colour = NA))

# violin_panel() (shared bulk-expression violin, coloured by group) now lives in _style.R
# so the non-MLL supplement (figS.R) can reuse it.

## ---- F (left): pediatric AML — TARGET bulk violins (4 leads + FLT3 ctrl + 3 clinical) ----
# Pediatric (discovery) cohort sets the shared x-order, by descending median; the adult
# panel reuses it so each gene sits at the same x in both → direct cross-cohort read.
# corrected scheme: drop the retired PPAC markers from the bulk validation; CD7 -> clinical comparator.
# `group` is derived from gene identity when a source CSV lacks it (the two source_data mirrors differ).
LEADS_ALL <- c("CD96", "CD9", "SUCNR1", "IL1RAP", "AMN")
drop_ppac <- function(df) {
  if (!"group" %in% names(df)) df$group <- dplyr::case_when(
    df$label %in% LEADS_ALL ~ "lead", df$label == "FLT3" ~ "control",
    df$label %in% PPAC_SPECIFIC ~ "ppac", TRUE ~ "clinical")
  df |> dplyr::filter(!(group == "ppac" & label != "CD7")) |>
    dplyr::mutate(group = dplyr::if_else(label == "CD7", "clinical", group))
}
# the pediatric (TARGET) and adult (BeatAML) source tables were built on different antigen panels;
# restrict BOTH to the targets present in each so E and F are true cross-cohort twins (same x-order,
# same genes), the direct comparison the panel intends. Genes unique to one cohort are dropped here.
ped_meta   <- drop_ppac(sd("fig5F_meta_all.csv"))
adult_meta <- drop_ppac(sd("fig5F_adult_beataml_meta.csv"))
shared     <- intersect(ped_meta$label, adult_meta$label)
ped_meta   <- ped_meta   |> dplyr::filter(label %in% shared)
adult_meta <- adult_meta |> dplyr::filter(label %in% shared)
gene_ord <- ped_meta |> dplyr::arrange(dplyr::desc(median)) |> dplyr::pull(label)
pF <- violin_panel(drop_ppac(sd("fig5F_violin.csv")) |> dplyr::filter(label %in% shared),
                   ped_meta, gene_ord, "Pediatric AML (TARGET)")

## ---- G (right): adult AML — BeatAML bulk violins (cross-cohort twin of E) ----
# Adult validation cohort: BeatAML (OHSU, Nature 2018), CPM profile, BM blasts ≥80% (n=104).
# Same 8 targets, same log2(CPM+1) scale, same x-order → directly comparable to E.
adult <- drop_ppac(sd("fig5F_adult_beataml_violin.csv")) |> dplyr::filter(label %in% shared) |>
  dplyr::rename(log2cpm = log2)
pE <- violin_panel(adult, adult_meta, gene_ord, "Adult AML (BeatAML)")

# bottom row = cross-cohort bulk twins: pediatric TARGET violin (pF) = panel F (left);
# adult BeatAML violin (pE) = panel G (right). Same 8 targets, same log2-CPM scale.
# panel D (all 49 states) + panel E (coverage summary) each get a full-width band; C sits above them,
# the cross-cohort bulk violins (F pediatric, G adult) below.
fig6 <- (pA | pB) / pC / pD / pBox / (pF | pE) +
  plot_layout(heights = c(1, 0.9, 1.35, 0.85, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 16))
savefig(fig6, "Figure5", 16, 21)
