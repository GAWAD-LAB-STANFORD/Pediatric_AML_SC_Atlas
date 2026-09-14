# ====================================================================
# Figure2_clusterprofiler.R  |  CD96 figure pipeline component
# reproducibility header (auto-generated; edit freely below this block)
# Inputs : (reads upstream object / raw matrix — see body)
# Outputs: (writes to ../figures or ../source_data — see body)
# Run    : Rscript Figure2_clusterprofiler.R   (or: bash run_all.sh)  |  Deps: see REPRODUCIBILITY_CHECKLIST.md
# Source data for this figure is in __SUBMISSION_PACKAGE/08_Source_Data/.
# ====================================================================

# Figure 2F — GO:BP over-representation per PPAC (clusterProfiler).
# RE-PATHED for the unified package: reads bundled source_data/legacy/, writes to
# figures/legacy/Figure2/, runs from anywhere (Rscript), and is msigdbr-version robust.
# Original (Sherlock paths) preserved in the working tree.
suppressPackageStartupMessages({
  library(clusterProfiler); library(enrichplot); library(msigdbr)
  library(ggplot2); library(dplyr); library(tidyr); library(org.Hs.eg.db)
})

## ---- portable paths (relative to this script) ----
.a <- commandArgs(FALSE); HERE <- dirname(sub("^--file=", "", grep("^--file=", .a, value = TRUE)))
if (length(HERE) == 0 || HERE == "") HERE <- "."
SD  <- normalizePath(file.path(HERE, "..", "..", "..", "source_data", "legacy"))
FIG <- file.path(HERE, "..", "..", "..", "figures", "legacy", "Figure2")
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)

## ---- MSigDB C5 GO:BP gene sets (handle old + new msigdbr APIs) ----
get_bp <- function() {
  for (call in list(
      function() msigdbr(species = "Homo sapiens", collection = "C5", subcollection = "GO:BP"),
      function() msigdbr(species = "Homo sapiens", category   = "C5", subcategory   = "GO:BP"),
      function() msigdbr(species = "Homo sapiens", category   = "C5", subcategory   = "BP")))
    { r <- tryCatch(call(), error = function(e) NULL); if (!is.null(r) && nrow(r)) return(r) }
  stop("could not load MSigDB C5:BP via msigdbr")
}
BP_t2g <- get_bp() %>% dplyr::select(gs_name, gene_symbol)

## ---- per-PPAC DE tables (Wilcoxon; padj<0.001 & |log2FC|>1) ----
rd <- function(f) read.delim2(file.path(SD, "Figure3", f), sep = ",", header = TRUE)
df_name <- rd("AML_PAC_Rank_gene_wilcoxon.csv")
df_fc   <- rd("AML_PAC_Rank_gene_wilcoxon_logFC.csv")
df_padj <- rd("AML_PAC_Rank_gene_wilcoxon_padj.csv")
AML_list <- c("X0_HSPC","FPAC_1","FPAC_2","Myeloid_Pro","PPAC_1","PPAC_2","PPAC_3","PPAC_4","PPAC_5")

go_dot <- function(genes, title) {
  em <- enricher(genes, TERM2GENE = BP_t2g, pvalueCutoff = 0.5, minGSSize = 1, qvalueCutoff = 0.5)
  if (is.null(em) || nrow(as.data.frame(em)) == 0) return(NULL)
  em@result$Description <- gsub("GOBP_", "", em@result$Description)
  dotplot(em, label_format = 26, title = title) +
    theme(axis.text = element_text(face = "bold", size = 10),
          axis.text.x = element_text(angle = 90)) +
    scale_color_gradient(low = "#0086a8", high = "#d04e00")
}

for (suf in AML_list) {
  df <- cbind(df_name[, grepl(suf, names(df_name))],
              df_fc[,   grepl(suf, names(df_fc))],
              df_padj[, grepl(suf, names(df_padj))])
  colnames(df) <- c("gene_symbol", "log2foldchange", "padj"); df <- data.frame(df)
  df$log2foldchange <- as.numeric(as.character(df$log2foldchange))
  df$padj <- as.numeric(as.character(df$padj))
  sig <- subset(df, padj < 0.001)
  up   <- sig$gene_symbol[sig$log2foldchange >  1]
  down <- sig$gene_symbol[sig$log2foldchange < -1]
  p1 <- go_dot(up,   paste0(suf, " GOBP upregulated"))
  if (!is.null(p1)) { pdf(file.path(FIG, paste0("GOBP_upregulated ", suf, ".pdf")), 5, 6); print(p1); dev.off() }
  p2 <- go_dot(down, paste0(suf, " GOBP downregulated"))
  if (!is.null(p2)) { pdf(file.path(FIG, paste0("GOBP_downregulated ", suf, ".pdf")), 5, 6); print(p2); dev.off() }
}

## ---- Module-15 hub-gene GO panel ----
BP2 <- tidyr::separate(BP_t2g, col = "gs_name", sep = "_", into = c("Term", "ID"), extra = "merge") %>%
  dplyr::select(ID, gene_symbol)
mod <- read.delim2(file.path(SD, "Figure2", "AML_Module_hub_gene_list.csv"), sep = ";", header = TRUE)
g15 <- dplyr::filter(mod, module %in% c("Module15"))$gene_name
em <- enricher(g15, TERM2GENE = BP2, pvalueCutoff = 0.5, minGSSize = 1, maxGSSize = 1000, qvalueCutoff = 0.5)
if (!is.null(em) && nrow(as.data.frame(em))) {
  pdf(file.path(FIG, "Module15_GOBP.pdf"), 5, 6)
  print(dotplot(em, font.size = 8, title = "Module 15", label_format = 30) +
        scale_color_gradient(low = "#0086a8", high = "#d04e00"))
  dev.off()
}
cat("Figure 2 GO dotplots ->", normalizePath(FIG), "\n")
