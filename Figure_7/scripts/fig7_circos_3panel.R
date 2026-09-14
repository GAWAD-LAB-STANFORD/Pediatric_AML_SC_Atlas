# ====================================================================
# fig7_circos_3panel.R  |  Figure 7A (target<->cell-type interaction circos)
# reproducibility header (auto-generated; edit freely below this block)
# Inputs : cellchat_3path_input.csv
# Outputs: (writes to ../figures or ../source_data — see body)
# Run    : Rscript fig7_circos_3panel.R   (or: bash run_all.sh)  |  Deps: see REPRODUCIBILITY_CHECKLIST.md
# Source data for this figure is in __SUBMISSION_PACKAGE/08_Source_Data/.
# ====================================================================

# Figure 7 — three interaction circos (CD96, ITGAX, OX40), one per target, AML-centric.
# Edge weight = co-expression interaction score = (% sender expressing ligand) x
# (% receiver expressing receptor) / 100, where the receptor side takes the MAX over the
# alternative ligands of that axis. CD96 axis = the nectin checkpoint (PVR / NECTIN1 /
# NECTIN2-CD112); ITGAX axis = ICAM1 / FCER2A / C3 / THY1; OX40 = TNFSF4. Autocrine =
# AML->AML self-loop. All three panels share ONE linewidth + loop scale so magnitudes are
# directly comparable across targets (CD96 dominates honestly).
suppressPackageStartupMessages({ library(ggplot2); library(dplyr) })
OUT <- normalizePath(file.path(dirname(dirname(getwd())), "figures", "cellchat"), mustWork = FALSE)
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
d <- read.csv("cellchat_3path_input.csv", stringsAsFactors = FALSE)
names(d)[names(d) == "FCER2"] <- "FCER2A"
GEN <- c("CD96","PVR","NECTIN1","NECTIN2","TNFSF4","TNFRSF4","C3","ICAM1","THY1","FCER2A","ITGAX","ITGB2")
COMPS <- c("AML (leukemic)","HSPC","Myeloid_Pro","Mature myeloid","B-lineage","T cells","NK","Erythroid")
AMLc <- "AML (leukemic)"
pct <- sapply(COMPS, function(cc) sapply(GEN, function(g) 100*mean(d[d$compartment==cc, g] > 0)))
rownames(pct) <- GEN; ncell <- as.numeric(table(factor(d$compartment, levels = COMPS)))

n <- length(COMPS); ang <- pi/2 - 2*pi*(0:(n-1))/n
nd <- data.frame(comp = COMPS, x = cos(ang), y = sin(ang), a = ang, sz = ncell)
nx <- setNames(nd$x, nd$comp); ny <- setNames(nd$y, nd$comp)
nd$lhj <- ifelse(cos(nd$a) > 0.15, 0, ifelse(cos(nd$a) < -0.15, 1, 0.5))
nd$lx <- nd$x * 1.18; nd$ly <- nd$y * 1.18
aml <- nd$comp == AMLc; nd$lx[aml] <- 0; nd$ly[aml] <- nd$y[aml] + 0.40; nd$lhj[aml] <- 0.5

# per-target directed edges (incl AML autocrine); CD96 axis now spans the full nectin family
edges_for <- function(tgt) {
  e <- if (tgt == "CD96")
    data.frame(from = AMLc, to = COMPS, w = pct["CD96", AMLc] *
                 pmax(pct["PVR", ], pct["NECTIN1", ], pct["NECTIN2", ]) / 100)
  else if (tgt == "ITGAX")
    data.frame(from = COMPS, to = AMLc, w = pmax(pct["ICAM1", ], pct["FCER2A", ],
                 pct["C3", ], pct["THY1", ]) * pct["ITGAX", AMLc] / 100)
  else
    data.frame(from = COMPS, to = AMLc, w = pct["TNFSF4", ] * pct["TNFRSF4", AMLc] / 100)
  e[e$w >= 0.2, ]
}
TINFO <- list(
  CD96  = list(col = "#c0392b", sub = "CD96 (AML) → PVR / NECTIN1 / NECTIN2 (CD112)   ·   nectin checkpoint   ·   AML = sender"),
  ITGAX = list(col = "#7f8c8d", sub = "ICAM1 / FCER2A / C3 → ITGAX+ITGB2 (AML)   ·   AML = receiver"),
  TNFRSF4 = list(col = "#16a085", sub = "OX40L (TNFSF4) → TNFRSF4 (AML)   ·   AML = receiver"))

WMAX <- max(sapply(names(TINFO), function(t) max(edges_for(t)$w)))   # shared scale across panels
RLOOP <- 0.16                                            # fixed autocrine-loop radius (same in every panel)

loop_pts <- function(comp, r) {                          # self-loop circle just outside a node
  cx <- nx[comp]*1.15; cy <- ny[comp]*1.15; t <- seq(0, 2*pi, length.out = 80)
  data.frame(x = cx + r*cos(t), y = cy + r*sin(t))
}

mkfig <- function(tgt) {
  E <- edges_for(tgt); col <- TINFO[[tgt]]$col
  auto <- E[E$from == E$to, ]; cross <- E[E$from != E$to, ]
  cross <- cross |> mutate(x0 = nx[from], y0 = ny[from], x1 = nx[to], y1 = ny[to],
                           mx = (x0+x1)/2*1.04, my = (y0+y1)/2*1.04)
  p <- ggplot() +
    geom_point(data = nd, aes(x, y, size = sz), shape = 21, fill = "#34495e", colour = "white", stroke = 0.4) +
    geom_text(data = nd, aes(lx, ly, label = comp, hjust = lhj), size = 2.9)
  if (nrow(cross)) p <- p +
    geom_curve(data = cross, aes(x0, y0, xend = x1, yend = y1, linewidth = w), colour = col,
               curvature = 0.2, alpha = 0.85, lineend = "round",
               arrow = grid::arrow(length = unit(0.16, "cm"), type = "closed")) +
    geom_label(data = cross, aes(mx, my, label = sprintf("%.1f", w)), size = 2.5, colour = col,
               label.size = 0, fill = "white", label.padding = unit(0.5, "pt"))
  if (nrow(auto)) {
    w <- auto$w[1]; lp <- loop_pts(AMLc, RLOOP)          # fixed size; THICKNESS alone carries the weight
    p <- p + geom_path(data = lp, aes(x, y, linewidth = w), colour = col, lineend = "round") +
      geom_label(aes(nx[AMLc] + 0.28, ny[AMLc] + 0.15 + RLOOP, label = sprintf("autocrine %.1f", w)),
                 hjust = 0, size = 2.5, colour = col, label.size = 0, fill = "white",
                 label.padding = unit(0.5,"pt"))
  }
  p + scale_linewidth(limits = c(0, WMAX), range = c(0.4, 3.6), guide = "none") +
    scale_size(range = c(2, 9), guide = "none") +
    coord_equal(clip = "off") + xlim(-1.65, 1.65) + ylim(-1.6, 1.95) +
    labs(title = tgt) +                                  # target name only; commentary -> caption
    theme_void() + theme(plot.title = element_text(face = "bold", size = 13, hjust = 0.5, colour = col))
}

suppressPackageStartupMessages(library(patchwork))
plots <- lapply(names(TINFO), mkfig)
for (i in seq_along(plots))
  ggsave(file.path(OUT, paste0("fig7_circos_", names(TINFO)[i], ".png")), plots[[i]],
         width = 5.6, height = 5.6, dpi = 220, bg = "white", device = ragg::agg_png)
figA <- wrap_plots(plots, nrow = 1)                      # no title/subtitle; commentary -> caption
ggsave(file.path(OUT, "fig7A_circos.png"), figA, width = 16.5, height = 6.2, dpi = 200,
       bg = "white", device = ragg::agg_png)
tryCatch(ggsave(file.path(OUT, "fig7A_circos.pdf"), figA, width = 16.5, height = 6.2,
                bg = "white", device = grDevices::cairo_pdf),
         error = function(e) ggsave(file.path(OUT, "fig7A_circos.pdf"), figA,
                width = 16.5, height = 6.2, bg = "white"))   # vector panel for the fitz assembler
cat(sprintf("wrote fig7A_circos.png (combined) · shared WMAX = %.2f\n", WMAX))
for (tgt in names(TINFO)) { cat("\n== ", tgt, " ==\n", sep=""); print(edges_for(tgt), row.names = FALSE) }
