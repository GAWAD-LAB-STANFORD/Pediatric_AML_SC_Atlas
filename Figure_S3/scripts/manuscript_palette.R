# Shared manuscript palette (ltc, hardcoded hex for portability).
# Sequential = heatmap0; Diverging = heatmap3 (colorblind-safe); prognosis = minou; patients = trio3.
suppressMessages({library(ggplot2); library(scales)})

SEQ  <- c("#001219","#005F73","#0A9396","#94D2BD","#E9D8A6","#EE9B00","#CA6702","#AE2012","#9B2226")  # heatmap0
DIV  <- c("#2c7bb6","#abd9e9","#ffffbf","#fdae61","#d7191c")   # heatmap3, low(blue)->high(red)
PROG <- c(favorable="#0A9396", poor="#AE2012")                  # SEQ endpoints (teal / red) - matches the dot-plot gradient
PAT  <- c(`multi-pt`="#009E73", `few-pt`="#56B4E9", `single-pt`="#E69F00")  # trio3 (Okabe-Ito)
GREY <- "#808080"

scale_fill_seq    <- function(...) scale_fill_gradientn(colours = SEQ, ...)
scale_colour_seq  <- function(...) scale_colour_gradientn(colours = SEQ, ...)
scale_fill_div    <- function(midpoint = 0, ...) scale_fill_gradientn(colours = DIV,
                       rescaler = function(x, to = c(0,1), from = range(x, na.rm = TRUE))
                         scales::rescale_mid(x, to, from, mid = midpoint), ...)
scale_fill_prog   <- function(...) scale_fill_manual(values = PROG, ...)
scale_fill_pat    <- function(...) scale_fill_manual(values = PAT, ...)
