# Pediatric AML Target Discovery
# ---------------------------------------------------------------------------
# Screen the surfaceome for marrow-sparing immunotherapy targets, then profile
# any candidate for efficacy, safety, cross-cohort validation and combinations.
# Driven entirely by the manuscript Source Data in ../source_data — no .h5ad.
#
#   shiny::runApp("app")          # from the repository root
# ---------------------------------------------------------------------------

library(shiny); library(bslib); library(ggplot2); library(DT); library(dplyr)

DATA <- local({
  for (p in c("source_data", file.path("..", "source_data"),
              file.path(dirname(getwd()), "source_data"))) if (dir.exists(p)) return(normalizePath(p))
  stop("source_data/ not found — run from the repository root or from app/")
})
rd <- function(...) read.csv(file.path(DATA, ...), stringsAsFactors = FALSE, check.names = FALSE)

# --- data ------------------------------------------------------------------
screen   <- rd("Figure_5", "fig5A_top15.csv")       # 1,704 genes past the marrow gate
funnel   <- rd("Figure_5", "fig5A_funnel.csv")
h2h      <- rd("Figure_5", "fig5C_headtohead.csv")  # 48 genes, with class labels
ls_expr  <- rd("Figure_5", "fig5D_ls.csv")          # 19 genes x leukemic states
cyto     <- rd("Figure_5", "fig5D_cyto.csv")        # 19 genes x cytogenetic subtype
meta_ped <- rd("Figure_5", "fig5F_meta_all.csv")            # TARGET bulk (paediatric)
meta_ad  <- rd("Figure_5", "fig5F_adult_beataml_meta.csv")  # Beat AML (adult)
sc_ad    <- rd("Figure_S23", "adult_scrna_overall.csv")     # adult scRNA
radar    <- rd("Figure_6", "fig6D_radar.csv")
tox_hema <- rd("Figure_6", "fig6tox_hema.csv")
tox_org  <- rd("Figure_6", "fig6tox_organ.csv")
combos   <- rd("Figure_7", "fig7_allpairs.csv")
emb      <- rd("Figure_3", "v3", "monocle_umap_coords.csv")
emb$compartment <- ifelse(grepl("^NORM", emb$group), "Normal", "Leukemic")

PROG  <- c(favorable = "#2F5D70", poor = "#8C1515", n.s. = "#B9C0C7", normal = "#7E57C2")
KLASS <- c(lead = "#8C1515", candidate = "#2F5D70", clinical = "#F6A30C",
           control = "#7E57C2", ppac = "#5B8C5A")
CARD  <- "#8C1515"; TEAL <- "#2F5D70"

theme_lab <- function(rot = 45) {
  theme_minimal(base_size = 13) +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold", size = 13),
          axis.text.x = element_text(angle = rot, hjust = if (rot > 0) 1 else .5))
}
note <- function(...) div(class = "text-muted", style = "font-size:.85rem;margin-top:.5rem", ...)
has  <- function(df, g) nrow(df[df$gene == g, ]) > 0

# --- ui --------------------------------------------------------------------
ui <- page_navbar(
  title = "Pediatric AML Target Discovery",
  theme = bs_theme(version = 5, primary = "#2F5D70", base_font = font_google("Inter")),
  sidebar = sidebar(
    width = 320,
    selectizeInput("gene", "Selected target", choices = NULL,
                   options = list(placeholder = "type a gene, e.g. CD96")),
    uiOutput("gene_card"),
    hr(),
    uiOutput("coverage_note")
  ),

  nav_panel(
    "Discover",
    layout_columns(
      col_widths = c(4, 8),
      card(
        card_header("Screening criteria"),
        sliderInput("f_aml",  "Min % leukemic cells positive", 0, 30, 5, step = 1),
        sliderInput("f_pts",  "Min patients with >20% positive", 0, 17, 2, step = 1),
        sliderInput("f_hspc", "Max % normal HSPC", 0, 10, 10, step = 0.5),
        sliderInput("f_mye",  "Max % myeloid progenitors", 0, 10, 10, step = 0.5),
        hr(),
        sliderInput("w", "Score weighting — breadth ↔ efficacy", 0, 1, 0.5, step = 0.05),
        note("Score = w x (% leukemic cells) + (1-w) x (% patients >20% positive). ",
             "At w = 0.5 this is the manuscript composite."),
        downloadButton("dl", "Download passing genes (CSV)", class = "btn-sm btn-primary mt-2")
      ),
      card(card_header(textOutput("funnel_title")), plotOutput("funnel_plot", height = 260),
           note("The 1,704 screened genes have already passed the ≤10% marrow-sparing gate; ",
                "your criteria narrow them further."))
    ),
    layout_columns(
      col_widths = c(6, 6),
      card(card_header("Efficacy vs normal-marrow expression"), plotOutput("gate_plot", height = 420),
           note("Dashed lines are your gates. Genes passing every criterion are highlighted.")),
      card(card_header("Passing targets"), DTOutput("hits"))
    )
  ),

  nav_panel(
    "Target profile",
    layout_columns(
      col_widths = c(7, 5),
      card(card_header(textOutput("ls_title")), plotOutput("ls_plot", height = 380)),
      card(card_header("Head-to-head with established antigens"), plotOutput("h2h_plot", height = 380),
           note("The 48 antigens profiled side by side. Leads and clinical-stage comparators labelled."))
    ),
    card(card_header("Coverage by cytogenetic subtype"), plotOutput("cyto_plot", height = 330),
         note("A target that works across subtypes is more broadly deployable than one confined to a single lesion."))
  ),

  nav_panel(
    "Safety",
    layout_columns(
      col_widths = c(6, 6),
      card(card_header("Normal marrow"), plotOutput("hema_plot", height = 400)),
      card(card_header("Vital organs (top 15)"), plotOutput("org_plot", height = 400))),
    card(card_header("Five-axis profile"), plotOutput("radar_plot", height = 300),
         note("Each axis is min-max scaled across the twelve antigens carried past the screen. Higher is better."))
  ),

  nav_panel(
    "Validation",
    card(card_header("Does the target hold beyond the discovery cohort?"),
         plotOutput("valid_plot", height = 420),
         note("Percent of patients targetable in paediatric TARGET bulk RNA-seq and adult Beat AML, ",
              "alongside percent of cells positive in an independent adult single-cell cohort."))
  ),

  nav_panel(
    "Combinations",
    card(card_header("Coverage vs marrow toxicity"), plotOutput("combo_plot", height = 440),
         note("Each point is a single antigen or an 'X OR Y' pair. Up and left is better.")),
    card(card_header("Pairs involving the selected target"), DTOutput("combo_tbl"))
  ),

  nav_panel(
    "Atlas",
    card(card_header("Single-cell embedding"),
      layout_columns(col_widths = c(3, 9),
        radioButtons("emb_fill", "Colour by",
                     c("Compartment" = "comp", "Pseudotime" = "pt", "Leukemic fraction" = "frac")),
        plotOutput("emb_plot", height = 460)),
      note("Monocle trajectory embedding. Each of the 57 groups is subsampled to 1,000 cells, ",
           "so apparent cell density reflects the subsampling, not true abundance — ",
           "compartment and pseudotime are the meaningful readouts here."))
  )
)

# --- server ----------------------------------------------------------------
server <- function(input, output, session) {
  updateSelectizeInput(session, "gene", choices = screen$gene, selected = "CD96", server = TRUE)
  g <- reactive(if (is.null(input$gene) || !nzchar(input$gene)) "CD96" else input$gene)

  scored <- reactive({
    d <- screen
    d$score <- input$w * d$scRNA_AML_pct +
      (1 - input$w) * (d$scRNA_n_pts_above20pct / max(d$scRNA_n_pts_above20pct, na.rm = TRUE) * 100)
    d$pass <- d$scRNA_AML_pct >= input$f_aml &
      d$scRNA_n_pts_above20pct >= input$f_pts &
      d$scRNA_HSPC_pct <= input$f_hspc &
      d$scRNA_Myeloid_pct <= input$f_mye
    d[order(-d$score), ]
  })
  hits <- reactive(scored()[scored()$pass, ])

  output$gene_card <- renderUI({
    r <- screen[screen$gene == g(), ]; if (!nrow(r)) return(NULL)
    tagList(
      h4(g(), style = "margin-bottom:.2rem"),
      div(class = "text-muted", sprintf("composite rank %d of %s", r$rank[1],
                                        format(nrow(screen), big.mark = ","))),
      hr(),
      tags$table(class = "table table-sm", tags$tbody(
        tags$tr(tags$td("Composite"),        tags$td(sprintf("%.1f", r$composite[1]))),
        tags$tr(tags$td("% leukemic cells"), tags$td(sprintf("%.1f%%", r$scRNA_AML_pct[1]))),
        tags$tr(tags$td("Patients >20%"),    tags$td(sprintf("%d of 28", r$scRNA_n_pts_above20pct[1]))),
        tags$tr(tags$td("Normal HSPC"),      tags$td(sprintf("%.1f%%", r$scRNA_HSPC_pct[1]))),
        tags$tr(tags$td("Myeloid prog."),    tags$td(sprintf("%.1f%%", r$scRNA_Myeloid_pct[1]))))))
  })

  output$coverage_note <- renderUI({
    avail <- c("per-state" = has(ls_expr, g()), "subtype" = has(cyto, g()),
               "safety" = has(tox_hema, g()), "validation" = has(meta_ped, g()))
    if (all(avail)) note("Full profile available for this target.")
    else note(sprintf("This gene has: %s. Panels without data will say so — only the %d antigens carried past the screen have safety and validation data.",
                      if (any(avail)) paste(names(avail)[avail], collapse = ", ") else "ranking only",
                      length(unique(radar$gene))))
  })

  # ---- Discover ----
  output$funnel_title <- renderText(sprintf("Screening funnel — %s genes pass", format(nrow(hits()), big.mark = ",")))
  output$funnel_plot <- renderPlot({
    d <- rbind(funnel, data.frame(stage = "your criteria", n = nrow(hits())))
    d$stage <- factor(d$stage, levels = rev(d$stage))
    d$mine <- d$stage == "your criteria"
    ggplot(d, aes(n, stage, fill = mine)) +
      geom_col() +
      geom_text(aes(label = format(n, big.mark = ",")), hjust = -0.15, size = 4) +
      scale_fill_manual(values = c(`FALSE` = "#B9C0C7", `TRUE` = CARD), guide = "none") +
      scale_x_continuous(expand = expansion(mult = c(0, .18))) +
      labs(x = NULL, y = NULL) + theme_lab(0)
  })

  output$gate_plot <- renderPlot({
    d <- scored()
    sel <- d[d$gene == g(), ]
    ggplot(d, aes(scRNA_AML_pct, scRNA_HSPC_pct)) +
      geom_point(aes(colour = pass), alpha = .6, size = 1.8) +
      geom_vline(xintercept = input$f_aml, linetype = "dashed", colour = CARD) +
      geom_hline(yintercept = input$f_hspc, linetype = "dashed", colour = CARD) +
      { if (nrow(sel)) geom_point(data = sel, size = 4, shape = 21, fill = "#F6A30C", colour = "black") } +
      { if (nrow(sel)) geom_text(data = sel, aes(label = gene), vjust = -1.2, fontface = "bold") } +
      scale_colour_manual(values = c(`FALSE` = "#D5D9DD", `TRUE` = TEAL), name = "Passes gates") +
      labs(x = "% leukemic cells positive", y = "% normal HSPC positive") + theme_lab(0)
  })

  output$hits <- renderDT({
    d <- hits()[, c("gene", "score", "scRNA_AML_pct", "scRNA_n_pts_above20pct",
                    "scRNA_HSPC_pct", "scRNA_Myeloid_pct", "rank")]
    names(d) <- c("Gene", "Score", "% AML cells", "Patients >20%", "% HSPC", "% Myeloid", "Paper rank")
    datatable(d, selection = "single", rownames = FALSE,
              options = list(pageLength = 12, order = list(list(1, "desc")))) |>
      formatRound(c("Score", "% AML cells", "% HSPC", "% Myeloid"), 1)
  })
  observeEvent(input$hits_rows_selected, {
    updateSelectizeInput(session, "gene", selected = hits()$gene[input$hits_rows_selected])
  })
  output$dl <- downloadHandler(
    filename = function() sprintf("aml_target_screen_%s.csv", Sys.Date()),
    content = function(f) write.csv(hits(), f, row.names = FALSE))

  # ---- Target profile ----
  output$ls_title <- renderText(sprintf("%s across leukemic states", g()))
  output$ls_plot <- renderPlot({
    d <- ls_expr[ls_expr$gene == g(), ]
    validate(need(nrow(d) > 0, sprintf("%s was not profiled per leukemic state (only the 19 antigens taken forward were).", g())))
    d$group <- factor(d$group, levels = d$group[order(-d$pct)])
    ggplot(d, aes(group, pct, fill = prognosis)) + geom_col() +
      scale_fill_manual(values = PROG, name = "State outcome") +
      labs(x = NULL, y = "% cells positive") + theme_lab()
  })

  output$h2h_plot <- renderPlot({
    d <- h2h
    d$is_sel <- d$gene == g()
    ggplot(d, aes(scRNA_AML_pct, scRNA_HSPC_pct, colour = klass)) +
      geom_point(size = 2.6, alpha = .85) +
      geom_text(data = d[d$do_label == "True" | d$is_sel, ], aes(label = label),
                vjust = -1, size = 3, show.legend = FALSE) +
      { if (any(d$is_sel)) geom_point(data = d[d$is_sel, ], size = 5, shape = 21,
                                      fill = NA, colour = "black", stroke = 1.1) } +
      scale_colour_manual(values = KLASS, name = NULL) +
      labs(x = "% leukemic cells positive", y = "% normal HSPC positive") + theme_lab(0)
  })

  output$cyto_plot <- renderPlot({
    d <- cyto[cyto$gene == g(), ]
    validate(need(nrow(d) > 0, sprintf("No cytogenetic breakdown for %s.", g())))
    d <- d[order(-d$pct), ]; d$group <- factor(d$group, levels = d$group)
    ggplot(d, aes(group, pct)) +
      geom_col(fill = TEAL) +
      geom_text(aes(label = sprintf("n=%d", n_pt)), vjust = -0.4, size = 3, colour = "#64707C") +
      labs(x = NULL, y = "% cells positive") +
      scale_y_continuous(expand = expansion(mult = c(0, .12))) + theme_lab()
  })

  # ---- Safety ----
  output$hema_plot <- renderPlot({
    d <- tox_hema[tox_hema$gene == g(), ]
    validate(need(nrow(d) > 0, sprintf("No normal-marrow profile for %s.", g())))
    d <- d[order(-d$pct), ]; d$cell <- factor(d$cell, levels = rev(d$cell))
    ggplot(d, aes(pct, cell)) + geom_col(fill = TEAL) +
      labs(x = "% cells positive", y = NULL) + theme_lab(0)
  })
  output$org_plot <- renderPlot({
    d <- tox_org[tox_org$gene == g(), ]
    validate(need(nrow(d) > 0, sprintf("No vital-organ profile for %s.", g())))
    d <- head(d[order(-d$pct), ], 15); d$organ <- factor(d$organ, levels = rev(d$organ))
    ggplot(d, aes(pct, organ)) + geom_col(fill = CARD) +
      labs(x = "% cells positive", y = NULL) + theme_lab(0)
  })
  output$radar_plot <- renderPlot({
    d <- radar[radar$gene == g(), ]
    validate(need(nrow(d) > 0, sprintf("No five-axis profile for %s.", g())))
    ax <- c("AML_coverage", "Patient_breadth", "Marrow_sparing", "Immune_sparing", "Organ_sparing")
    p <- data.frame(axis = factor(ax, levels = ax), value = as.numeric(d[1, ax]))
    ggplot(p, aes(axis, value)) +
      geom_segment(aes(xend = axis, y = 0, yend = value), colour = "#B9C0C7", linewidth = 1) +
      geom_point(size = 6, colour = CARD) +
      geom_text(aes(label = sprintf("%.0f", value)), colour = "white", size = 3) +
      ylim(0, 100) + labs(x = NULL, y = "score") + theme_lab(20)
  })

  # ---- Validation ----
  output$valid_plot <- renderPlot({
    gg <- g()
    rows <- list()
    a <- meta_ped[meta_ped$gene == gg, ]; b <- meta_ad[meta_ad$gene == gg, ]; c3 <- sc_ad[sc_ad$gene == gg, ]
    if (nrow(a)) rows[[length(rows)+1]] <- data.frame(cohort = "Paediatric TARGET\n(bulk, n=475)", metric = "% patients targetable", value = a$pct_targetable[1])
    if (nrow(b)) rows[[length(rows)+1]] <- data.frame(cohort = sprintf("Adult Beat AML\n(bulk, n=%d)", b$n[1]), metric = "% patients targetable", value = b$pct_targetable[1])
    if (nrow(c3)) rows[[length(rows)+1]] <- data.frame(cohort = sprintf("Adult scRNA\n(%s cells)", format(c3$n_cells[1], big.mark=",")), metric = "% cells positive", value = c3$pct[1])
    validate(need(length(rows) > 0, sprintf("No cross-cohort validation data for %s — only the twelve antigens carried past the screen were validated.", gg)))
    d <- do.call(rbind, rows)
    ggplot(d, aes(cohort, value, fill = metric)) +
      geom_col(width = .6) +
      geom_text(aes(label = sprintf("%.1f%%", value)), vjust = -0.4, size = 4) +
      scale_fill_manual(values = c("% patients targetable" = TEAL, "% cells positive" = CARD), name = NULL) +
      scale_y_continuous(expand = expansion(mult = c(0, .15))) +
      labs(x = NULL, y = NULL) + theme_lab(0)
  })

  # ---- Combinations ----
  output$combo_plot <- renderPlot({
    d <- combos; d$sel <- d$g1 == g() | d$g2 == g()
    lab <- d[d$sel | d$across >= quantile(d$across, .92), ]
    ggplot(d, aes(toxicity, across)) +
      geom_point(aes(colour = sel, size = sel), alpha = .85) +
      geom_text(data = lab, aes(label = combo), size = 3, vjust = -1, colour = "#17212B") +
      scale_colour_manual(values = c(`FALSE` = "#B9C0C7", `TRUE` = CARD), guide = "none") +
      scale_size_manual(values = c(`FALSE` = 2, `TRUE` = 4), guide = "none") +
      labs(x = "Normal-marrow toxicity (%)", y = "% patients with >20% of blasts covered") +
      theme_lab(0)
  })
  output$combo_tbl <- renderDT({
    d <- combos[combos$g1 == g() | combos$g2 == g(), c("combo", "kind", "within", "across", "toxicity")]
    validate(need(nrow(d) > 0, sprintf("%s is not part of the combination analysis.", g())))
    names(d) <- c("Combination", "Kind", "Within-patient depth (%)", "Patient breadth (%)", "Marrow toxicity (%)")
    datatable(d, rownames = FALSE, options = list(pageLength = 8, order = list(list(3, "desc")))) |>
      formatRound(3:5, 1)
  })

  # ---- Atlas ----
  output$emb_plot <- renderPlot({
    p <- switch(input$emb_fill,
      comp = ggplot(emb, aes(UMAP1, UMAP2, colour = compartment)) +
        geom_point(size = .35, alpha = .45) +
        scale_colour_manual(values = c(Leukemic = CARD, Normal = TEAL), name = NULL) +
        guides(colour = guide_legend(override.aes = list(size = 4, alpha = 1))),
      pt = ggplot(emb, aes(UMAP1, UMAP2, z = monocle_pt)) +
        stat_summary_hex(bins = 70, fun = median) +
        scale_fill_viridis_c(name = "Pseudotime", option = "mako"),
      frac = ggplot(emb, aes(UMAP1, UMAP2, z = as.numeric(compartment == "Leukemic"))) +
        stat_summary_hex(bins = 70, fun = mean) +
        scale_fill_gradient(low = TEAL, high = CARD, name = "Leukemic\nfraction", labels = scales::percent))
    p + coord_equal() + theme_lab(0)
  })
}

shinyApp(ui, server)
