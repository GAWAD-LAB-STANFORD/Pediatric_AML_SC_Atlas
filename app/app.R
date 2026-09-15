# Pediatric AML Atlas Browser
# ---------------------------------------------------------------------------
# Interactive browser over the manuscript Source Data: surface-target ranking,
# per-leukemic-state expression, normal-tissue safety, and combination coverage.
# Everything is driven by the CSVs in ../source_data — no atlas .h5ad required.
#
#   shiny::runApp("app")            # from the repository root
# ---------------------------------------------------------------------------

library(shiny); library(bslib); library(ggplot2); library(DT); library(dplyr)

# --- locate source_data whether run from the repo root or from app/ ---------
DATA <- local({
  for (p in c("source_data", file.path("..", "source_data"),
              file.path(dirname(getwd()), "source_data"))) if (dir.exists(p)) return(normalizePath(p))
  stop("source_data/ not found — run from the repository root or from app/")
})
rd <- function(...) read.csv(file.path(DATA, ...), stringsAsFactors = FALSE, check.names = FALSE)

# --- data ------------------------------------------------------------------
rank_tbl <- rd("Figure_5", "fig5A_top15.csv")
ls_expr  <- rd("Figure_5", "fig5D_ls.csv")
radar    <- rd("Figure_6", "fig6D_radar.csv")
tox_hema <- rd("Figure_6", "fig6tox_hema.csv")
tox_org  <- rd("Figure_6", "fig6tox_organ.csv")
combos   <- rd("Figure_7", "fig7_allpairs.csv")
umap     <- rd("Figure_3", "v3", "monocle_umap_coords.csv")

DEEP <- sort(unique(radar$gene))          # the 12 antigens carried through Fig 6-7
PROG <- c(favorable = "#2F5D70", poor = "#8C1515", n.s. = "#B9C0C7", normal = "#7E57C2")

theme_lab <- function() {
  theme_minimal(base_size = 13) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          plot.title = element_text(face = "bold", size = 13),
          axis.text.x = element_text(angle = 45, hjust = 1))
}
note <- function(txt) div(class = "text-muted", style = "font-size:.85rem;margin-top:.4rem", txt)

# --- ui --------------------------------------------------------------------
ui <- page_navbar(
  title = "Pediatric AML Atlas",
  theme = bs_theme(version = 5, primary = "#2F5D70", base_font = font_google("Inter")),
  sidebar = sidebar(
    width = 300,
    selectizeInput("gene", "Surface target", choices = NULL,
                   options = list(placeholder = "type a gene, e.g. CD96")),
    uiOutput("gene_card"),
    hr(),
    note(paste0(format(nrow(rank_tbl), big.mark = ","),
                " surfaceome genes ranked. ", length(DEEP),
                " carry the full safety and combination panels."))
  ),
  nav_panel("Targets", card(card_header("Composite ranking"), DTOutput("tbl"),
    note("Composite = mean of (% leukemic cells positive) and (% patients with >20% positive). Click a row to select a gene."))),
  nav_panel("Leukemic states", card(card_header(textOutput("ls_title")), plotOutput("ls_plot", height = 420),
    note("Percent of cells positive per leukemic state, coloured by the state's TARGET outcome association. HSPC is the normal-marrow baseline."))),
  nav_panel("Safety", layout_columns(
    card(card_header("Normal marrow"), plotOutput("hema_plot", height = 380)),
    card(card_header("Vital organs (top 15)"), plotOutput("org_plot", height = 380)))),
  nav_panel("Combinations", card(card_header("Coverage vs marrow toxicity"), plotOutput("combo_plot", height = 460),
    note("Each point is a single antigen or an 'X OR Y' pair. Up and left is better: broader patient coverage, lower normal-marrow toxicity."))),
  nav_panel("Atlas", card(card_header("Single-cell atlas"),
    layout_columns(col_widths = c(4, 8),
      radioButtons("umap_fill", "Colour by", c("Cell density" = "n", "Pseudotime" = "pt")),
      plotOutput("umap_plot", height = 460)),
    note("Cells are binned into hexagons so rendering cost scales with the grid, not the cell count.")))
)

# --- server ----------------------------------------------------------------
server <- function(input, output, session) {
  updateSelectizeInput(session, "gene", choices = rank_tbl$gene, selected = "CD96", server = TRUE)
  g <- reactive(if (is.null(input$gene) || !nzchar(input$gene)) "CD96" else input$gene)

  output$gene_card <- renderUI({
    r <- rank_tbl[rank_tbl$gene == g(), ]
    if (!nrow(r)) return(NULL)
    tagList(
      h4(g(), style = "margin-bottom:.2rem"),
      div(class = "text-muted", sprintf("rank %d of %s", r$rank[1], format(nrow(rank_tbl), big.mark = ","))),
      hr(),
      tags$table(class = "table table-sm",
        tags$tbody(
          tags$tr(tags$td("Composite"),      tags$td(sprintf("%.1f", r$composite[1]))),
          tags$tr(tags$td("% leukemic cells"), tags$td(sprintf("%.1f%%", r$scRNA_AML_pct[1]))),
          tags$tr(tags$td("Patients >20%"),  tags$td(sprintf("%d", r$scRNA_n_pts_above20pct[1]))),
          tags$tr(tags$td("Normal HSPC"),    tags$td(sprintf("%.1f%%", r$scRNA_HSPC_pct[1]))),
          tags$tr(tags$td("Myeloid prog."),  tags$td(sprintf("%.1f%%", r$scRNA_Myeloid_pct[1]))))),
      if (!(g() %in% DEEP)) note("No safety or combination data for this gene — it was not carried past the Figure 5 screen."))
  })

  output$tbl <- renderDT({
    d <- rank_tbl[, c("rank", "gene", "composite", "scRNA_AML_pct", "scRNA_n_pts_above20pct",
                      "scRNA_HSPC_pct", "scRNA_Myeloid_pct")]
    names(d) <- c("Rank", "Gene", "Composite", "% AML cells", "Patients >20%", "% HSPC", "% Myeloid")
    datatable(d, selection = "single", rownames = FALSE,
              options = list(pageLength = 15, order = list(list(0, "asc")))) |>
      formatRound(c("Composite", "% AML cells", "% HSPC", "% Myeloid"), 1)
  })
  observeEvent(input$tbl_rows_selected, {
    updateSelectizeInput(session, "gene", selected = rank_tbl$gene[input$tbl_rows_selected])
  })

  output$ls_title <- renderText(sprintf("%s across leukemic states", g()))
  output$ls_plot <- renderPlot({
    d <- ls_expr[ls_expr$gene == g(), ]
    validate(need(nrow(d) > 0, sprintf("%s was not profiled per leukemic state.", g())))
    d$group <- factor(d$group, levels = d$group[order(-d$pct)])
    ggplot(d, aes(group, pct, fill = prognosis)) +
      geom_col() +
      scale_fill_manual(values = PROG, name = "State outcome") +
      labs(x = NULL, y = "% cells positive") + theme_lab()
  })

  output$hema_plot <- renderPlot({
    d <- tox_hema[tox_hema$gene == g(), ]
    validate(need(nrow(d) > 0, sprintf("No normal-marrow profile for %s.", g())))
    d <- d[order(-d$pct), ]; d$cell <- factor(d$cell, levels = rev(d$cell))
    ggplot(d, aes(pct, cell)) +
      geom_col(fill = "#2F5D70") +
      labs(x = "% cells positive", y = NULL) + theme_lab() +
      theme(axis.text.x = element_text(angle = 0))
  })

  output$org_plot <- renderPlot({
    d <- tox_org[tox_org$gene == g(), ]
    validate(need(nrow(d) > 0, sprintf("No vital-organ profile for %s.", g())))
    d <- head(d[order(-d$pct), ], 15); d$organ <- factor(d$organ, levels = rev(d$organ))
    ggplot(d, aes(pct, organ)) +
      geom_col(fill = "#8C1515") +
      labs(x = "% cells positive", y = NULL) + theme_lab() +
      theme(axis.text.x = element_text(angle = 0))
  })

  output$combo_plot <- renderPlot({
    d <- combos
    d$sel <- d$g1 == g() | d$g2 == g()
    ggplot(d, aes(toxicity, across)) +
      geom_point(aes(colour = sel, size = sel), alpha = .85) +
      ggrepel_or_text(d) +
      scale_colour_manual(values = c(`FALSE` = "#B9C0C7", `TRUE` = "#8C1515"), guide = "none") +
      scale_size_manual(values = c(`FALSE` = 2, `TRUE` = 4), guide = "none") +
      labs(x = "Normal-marrow toxicity (%)", y = "% patients with >20% blasts covered") +
      theme_lab() + theme(axis.text.x = element_text(angle = 0))
  })

  output$umap_plot <- renderPlot({
    base <- if (input$umap_fill == "pt") {
      ggplot(umap, aes(UMAP1, UMAP2, z = monocle_pt)) +
        stat_summary_hex(bins = 70, fun = median) +
        scale_fill_viridis_c(name = "Pseudotime", option = "mako")
    } else {
      ggplot(umap, aes(UMAP1, UMAP2)) +
        geom_hex(bins = 70) +
        scale_fill_viridis_c(name = "Cells", option = "mako", trans = "log10")
    }
    base + coord_equal() + theme_lab() + theme(axis.text.x = element_text(angle = 0))
  })
}

# label only the selected/near-frontier combos, without requiring ggrepel
ggrepel_or_text <- function(d) {
  lab <- d[d$sel | d$across >= quantile(d$across, .9), ]
  geom_text(data = lab, aes(label = combo), size = 3, vjust = -1, colour = "#17212B")
}

shinyApp(ui, server)
