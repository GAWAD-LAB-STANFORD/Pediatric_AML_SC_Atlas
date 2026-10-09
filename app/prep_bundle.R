# Build the app's data bundle.
# ---------------------------------------------------------------------------
# The app used to read 24 MB of CSV at startup, 23 MB of which was two files:
# the per-gene/per-bin expression table (936,506 rows) and the 96,627-cell
# embedding. Parsing those with read.csv inside WebAssembly is what made the
# browser sit on a spinner, and scanning the 936,506-row table on every plot
# made each gene change slow as well.
#
# This turns both into compact binary: expression becomes a dense uint8 matrix
# indexed [bin, gene] - O(1) lookup instead of a linear scan - and the
# embedding keeps factors instead of repeated strings. Everything is written to
# one .rds, which R deserialises rather than parses.
#
#   Rscript app/prep_bundle.R [source_data] [out.rds]
# ---------------------------------------------------------------------------

build_bundle <- function(DATA) {
  rd <- function(...) read.csv(file.path(DATA, ...), stringsAsFactors = FALSE, check.names = FALSE)

  B <- list(
    screen   = rd("Figure_5", "fig5A_top15.csv"),
    funnel   = rd("Figure_5", "fig5A_funnel.csv"),
    h2h      = rd("Figure_5", "fig5C_headtohead.csv"),
    ls_expr  = rd("Figure_5", "fig5D_ls.csv"),
    cyto     = rd("Figure_5", "fig5D_cyto.csv"),
    meta_ped = rd("Figure_5", "fig5F_meta_all.csv"),
    meta_ad  = rd("Figure_5", "fig5F_adult_beataml_meta.csv"),
    sc_ad    = rd("Figure_S23", "adult_scrna_overall.csv"),
    radar    = rd("Figure_6", "fig6D_radar.csv"),
    tox_hema = rd("Figure_6", "fig6tox_hema.csv"),
    tox_org  = rd("Figure_6", "fig6tox_organ.csv"),
    combos   = rd("Figure_7", "fig7_allpairs.csv"),
    lsspec   = rd("Figure_S16", "funnel_ls_specific.csv"),
    lscomp   = rd("Figure_S16", "funnel_ls_specific_composition.csv"),
    ebins    = rd("Figure_2", "fig2A_bins.csv"))

  # embedding: drop the barcode, factor the labels, round the coordinates
  emb <- rd("Figure_2", "fig2A_umap_cells.csv")
  LS_LEVELS <- paste0("LS_", 1:49)
  B$emb <- data.frame(
    UMAP1       = round(emb$UMAP1, 3),
    UMAP2       = round(emb$UMAP2, 3),
    compartment = factor(emb$compartment),
    cell_type   = factor(emb$cell_type),
    LS          = factor(emb$LS, levels = c("", LS_LEVELS)),
    # which grid square each cell falls in, so expression can be painted onto the
    # cells themselves rather than quantised into squares
    bin         = as.integer(emb$bin))

  # expression: dense uint8 matrix, bins x genes
  gx   <- rd("Figure_2", "fig2A_gene_bin_expr.csv")
  gen  <- sort(unique(gx$gene))
  gmax <- max(gx$mean)
  M    <- matrix(0L, nrow = nrow(B$ebins), ncol = length(gen))
  M[cbind(match(gx$bin, B$ebins$bin), match(gx$gene, gen))] <-
    as.integer(round(gx$mean / gmax * 255))
  stopifnot(!anyNA(M), max(M) <= 255)
  B$gmat  <- matrix(as.raw(M), nrow = nrow(B$ebins), dimnames = list(NULL, gen))
  B$gmax  <- gmax
  B$genes_mapped <- gen
  B
}

if (!interactive() && sys.nframe() == 0L) {
  a   <- commandArgs(TRUE)
  src <- if (length(a) >= 1) a[1] else "source_data"
  out <- if (length(a) >= 2) a[2] else "app/bundle.rds"
  B   <- build_bundle(src)
  # gzip, not xz: xz saves 0.5 MB but decompresses ~5x slower, and startup time
  # is what the browser feels.
  saveRDS(B, out, compress = "gzip")
  cat(sprintf("bundle -> %s (%.1f MB) | %d mapped genes x %d bins | %s cells\n",
              out, file.size(out) / 1e6, length(B$genes_mapped), nrow(B$ebins),
              format(nrow(B$emb), big.mark = ",")))
}
