# ====================================================================
# figS_venn.R  |  Discovered (Fig-5A composite top-15) vs in-trial AML antigens
# DATA-DRIVEN: the composite top-15 is read from fig5A_top15.csv and the trial
# top-10 from fig_aml_trials.csv, so the Venn cannot go stale (it previously hard-
# coded an outdated top-15). Overlap computed by gene symbol. No figure subtitle
# (description lives in the caption). CD96 crimson; IL1RAP sky-blue.
# ====================================================================
source("_style.R")
suppressPackageStartupMessages(library(dplyr))

# --- current strict Fig-5A composite top-15 (gene symbols) ---
top15 <- sd("fig5A_top15.csv") |> arrange(rank) |> head(15) |> pull(gene)
# --- top-10 clinically-pursued AML antigens (ClinicalTrials.gov) + gene symbols for overlap ---
tgene <- c(CD33="CD33", CD123="IL3RA", "CLL-1"="CLEC12A", CD38="CD38", CD7="CD7", CD70="CD70",
           "NKG2D-L"="KLRK1L", "TIM-3"="HAVCR2", "B7-H3"="CD276", IL1RAP="IL1RAP")
trials <- sd("fig_aml_trials.csv") |> arrange(desc(n_trials)) |> head(10) |>
  mutate(gene = tgene[antigen])
ov_genes  <- intersect(top15, trials$gene)                      # overlap by gene
disp_gene <- function(g) { i <- match(g, trials$gene); ifelse(is.na(i), g, trials$label[i]) }  # clinical label if in trials

left_only  <- trials$label[!trials$gene %in% ov_genes]          # trials only
overlap    <- vapply(ov_genes, disp_gene, "")                   # both (clinical label)
right_only <- setdiff(top15, ov_genes)                          # composite only (gene symbols)
nL <- length(left_only); nO <- length(overlap); nR <- length(right_only)

mk <- function(genes, x, ytop, ybot) {
  if (!length(genes)) return(NULL)
  ys <- if (length(genes) == 1) 0 else seq(ytop, ybot, length.out = length(genes))
  data.frame(gene = genes, x = x, y = ys, stringsAsFactors = FALSE)
}
# Labels sit in each circle's wide middle band (compressed vertical range) so the top/
# bottom entries never reach the poles, where the circle narrows and text would spill
# outside the outline. Overlap is spread vertically (robust to 1 or 2 shared genes).
g <- rbind(mk(left_only,  -1.85, 1.25, -1.25),
           mk(overlap,     0.00, 0.42, -0.42),
           mk(right_only,  1.80,  1.35, -1.35))
# colour: CD96 crimson, IL1RAP sky-blue (leads that appear); everything else slate
g$col <- dplyr::case_when(g$gene == "CD96" ~ unname(LEAD_COLORS["CD96"]),
                          grepl("IL1RAP", g$gene) ~ unname(LEAD_COLORS["IL1RAP"]),
                          TRUE ~ "#2c3e50")
g$face <- ifelse(g$gene %in% c("CD96") | grepl("IL1RAP", g$gene), "bold", "plain")

circ <- function(cx, cy, r, n = 200) { t <- seq(0, 2*pi, length.out = n); data.frame(x = cx + r*cos(t), y = cy + r*sin(t)) }
cL <- circ(-0.95, 0, 1.95); cR <- circ(0.95, 0, 1.95)

p <- ggplot() +
  geom_polygon(data = cL, aes(x, y), fill = "#2471a3", alpha = 0.10, colour = "#2471a3", linewidth = 0.7) +
  geom_polygon(data = cR, aes(x, y), fill = "#c0392b", alpha = 0.07, colour = "#c0392b", linewidth = 0.7) +
  geom_text(data = g, aes(x, y, label = gene, colour = col, fontface = face), size = 2.5) +
  scale_colour_identity() +
  annotate("text", x = -1.95, y = 2.28, label = "Clinically pursued in AML\n(top-10, ClinicalTrials.gov)",
           fontface = "bold", size = 3.1, colour = "#1a5276", lineheight = 0.9) +
  annotate("text", x = 1.95, y = 2.28, label = "Figure-5A composite\ntop-15 candidates",
           fontface = "bold", size = 3.1, colour = "#7b241c", lineheight = 0.9) +
  annotate("text", x = -1.85, y = 1.66, label = nL, fontface = "bold", size = 3.4, colour = "#1a5276") +
  annotate("text", x = 0.0,  y = 0.92, label = nO, fontface = "bold", size = 3.4, colour = "#6c3483") +
  annotate("text", x = 1.80, y = 1.66, label = nR, fontface = "bold", size = 3.4, colour = "#7b241c") +
  coord_equal(clip = "off") +
  labs(title = "Composite candidates vs the clinically-pursued AML surface antigens") +
  theme_void(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
        plot.margin = margin(10, 16, 10, 16))

savefig(p, "supplementary/FigureS13_venn_composite_vs_trials", 9, 7.5)
cat(sprintf("Done -> trials-only %d, overlap %d (%s), composite-only %d\n",
            nL, nO, paste(overlap, collapse=","), nR))
