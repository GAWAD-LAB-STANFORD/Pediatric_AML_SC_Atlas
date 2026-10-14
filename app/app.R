# Pediatric AML Target Discovery
# ---------------------------------------------------------------------------
# Screen the surfaceome for marrow-sparing immunotherapy targets, then profile
# any candidate for efficacy, safety, cross-cohort validation and combinations.
# Driven entirely by the manuscript Source Data in ../source_data \u2014 no .h5ad.
#
#   shiny::runApp("app")          # from the repository root
# ---------------------------------------------------------------------------

library(shiny); library(bslib); library(ggplot2); library(DT)

# Data arrives as one binary bundle (app/prep_bundle.R). Reading 24 MB of CSV at
# startup was what left the browser on a spinner; the bundle is 2 MB and R
# deserialises it instead of parsing it. Falling back to the raw CSVs keeps the
# app runnable straight from the repository during development.
B <- local({
  for (p in c("bundle.rds", file.path("app", "bundle.rds"), file.path("..", "bundle.rds")))
    if (file.exists(p)) return(readRDS(p))
  prep <- Filter(file.exists, c("prep_bundle.R", file.path("app", "prep_bundle.R")))
  data <- Filter(dir.exists, c("source_data", file.path("..", "source_data"),
                               file.path(dirname(getwd()), "source_data")))
  if (length(prep) && length(data)) { source(prep[1], local = TRUE); return(build_bundle(data[1])) }
  stop("no bundle.rds and no source_data/ - run Rscript app/prep_bundle.R first")
})
screen <- B$screen; funnel <- B$funnel; h2h <- B$h2h; ls_expr <- B$ls_expr
cyto <- B$cyto; meta_ped <- B$meta_ped; meta_ad <- B$meta_ad; sc_ad <- B$sc_ad
radar <- B$radar; tox_hema <- B$tox_hema; tox_org <- B$tox_org; combos <- B$combos
lsspec <- B$lsspec; lscomp <- B$lscomp; emb <- B$emb; ebins <- B$ebins
N_PT <- 28                                        # patients in the discovery cohort
LS_LEVELS <- paste0("LS_", 1:49)

# expression: dense uint8 matrix [bin, gene], so a gene lookup is a column read
GMAT <- B$gmat; GSCALE <- B$gmax / 255
EXPR_GENES <- B$genes_mapped
gene_bins <- function(g) { i <- match(g, EXPR_GENES)
                           if (is.na(i)) NULL else as.integer(GMAT[, i]) * GSCALE }

# Each panel is measured on a different subset of genes. Rather than let someone
# pick a gene and hit an empty panel, the sidebar can narrow the picker to the
# genes a given panel actually covers.
SCOPES <- list(
  "Every screened gene"       = sort(unique(screen$gene)),
  "Mapped on the atlas"       = sort(unique(EXPR_GENES)),
  "Head-to-head panel"        = sort(unique(h2h$gene)),
  "Profiled per leukemic state" = sort(unique(ls_expr$gene)),
  "Full safety and validation"  = sort(unique(radar$gene)))
SCOPE_LABELS <- sprintf("%s (%s)", names(SCOPES),
                        format(lengths(SCOPES), big.mark = ",", trim = TRUE))

PROG  <- c(favorable = "#2F5D70", poor = "#8C1515", n.s. = "#B9C0C7", normal = "#7E57C2")
KLASS <- c(lead = "#8C1515", candidate = "#2F5D70", clinical = "#F6A30C",
           control = "#7E57C2", state_specific = "#5B8C5A")
CARD  <- "#8C1515"; TEAL <- "#2F5D70"

theme_lab <- function(rot = 45) {
  theme_minimal(base_size = 13) +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold", size = 13),
          axis.text.x = element_text(angle = rot, hjust = if (rot > 0) 1 else .5))
}
note <- function(...) div(class = "text-muted", style = "font-size:.85rem;margin-top:.5rem", ...)
ext  <- function(label, url) a(label, href = url, target = "_blank", rel = "noopener",
                               style = "font-weight:600")
gene_links <- function(g) tagList(
  ext("CZ CELLxGENE", paste0("https://cellxgene.cziscience.com/gene-expression?genes=", g)), " \u00b7 ",
  ext("Human Protein Atlas", paste0("https://www.proteinatlas.org/search/", g)), " \u00b7 ",
  ext("GTEx", paste0("https://gtexportal.org/home/gene/", g)))
# Shown in place of a panel this study cannot fill for the selected gene, so the
# answer is a pointer to where the measurement does exist rather than a blank.
elsewhere <- function(g, what, why = NULL)
  div(class = "border rounded p-3", style = "background:#F7F8F9;font-size:.9rem",
      div(strong(g), sprintf(" has no %s in this study.", what)),
      if (!is.null(why)) div(class = "text-muted", style = "margin-top:.25rem", why),
      div(style = "margin-top:.6rem", "Look it up in ", gene_links(g), "."))
has  <- function(df, g) nrow(df[df$gene == g, ]) > 0

# --- ui --------------------------------------------------------------------
ui <- page_navbar(
  id = "nav",
  title = "Pediatric AML Target Discovery",
  header = tags$style(HTML(
    ".navbar .navbar-brand{font-size:1.45rem;font-weight:700}",
    ".navbar .nav-link{font-size:1.18rem;font-weight:600;padding:.5rem 1rem}",
    ".navbar .nav-link.active{font-weight:700}")),
  # No font_google(): webR has no curl, so a downloaded font breaks the WebAssembly build.
  theme = bs_theme(version = 5, primary = "#2F5D70",
                   base_font = font_face(family = "Inter", src = "local('Inter')")),
  sidebar = sidebar(
    width = 300,
    selectizeInput("gene", "Selected target", choices = SCOPES[[1]], selected = "CD96",
                   options = list(placeholder = "type a gene, e.g. CD96", maxOptions = 200)),
    radioButtons("glist", "Limit the list to genes with",
                 choiceNames = SCOPE_LABELS, choiceValues = names(SCOPES),
                 selected = names(SCOPES)[1]),
    uiOutput("gene_card"),
    conditionalPanel(
      "input.nav == 'Discover'",
      hr(),
      selectInput("ctx", "Disease context",
                  c("Whole cohort", sort(unique(lsspec$group)))),
      uiOutput("ctx_note"),
      strong("Screening criteria"),
      sliderInput("f_aml",  "Min % leukemic cells positive", 0, 70, 5, step = 1),
      sliderInput("f_pts",  "Min % of patients positive", 0, 100, 10, step = 5),
      sliderInput("f_hspc", "Max % normal HSPC", 0, 10, 10, step = 0.5),
      sliderInput("f_mye",  "Max % myeloid progenitors", 0, 10, 10, step = 0.5),
      sliderInput("w", "Score weighting \u2014 breadth \u2194 efficacy", 0, 1, 0.5, step = 0.05),
      downloadButton("dl", "Download hits (CSV)", class = "btn-sm btn-primary")
    ),
    hr(),
    uiOutput("coverage_note")
  ),

  nav_panel(
    "Atlas",
    card(
      card_header("Paediatric AML single-cell atlas"),
      layout_columns(
        col_widths = c(3, 9),
        radioButtons("emb_fill", "Colour by",
                     c("Gene expression" = "expr", "Compartment" = "comp", "Cell type" = "type",
                       "Leukemic state" = "state", "Cell density" = "dens")),
        div(uiOutput("emb_miss"), plotOutput("emb_plot", height = 520))),
      note(strong("96,627 cells"), " \u2014 70,108 leukemic and 26,519 normal marrow, across 29 annotated cell types. ",
           "Of the leukemic cells, 57,166 fall in one of the 49 reproducible states (LS_1\u2013LS_49); the other 12,942 sit in ",
           "clusters that failed the Jaccard \u2265 0.35 stability cut and are deliberately left unassigned. ",
           "Expression is the mean log1p(CPM) per grid bin, available for ", textOutput("n_expr_genes", inline = TRUE),
           " of the screened targets.")
    )
  ),

  nav_panel(
    "Discover",
    card(card_header(textOutput("funnel_title")), plotOutput("funnel_plot", height = 300),
         note("Set your criteria in the sidebar. The 1,704 screened genes have already passed the ",
              "\u226410% marrow-sparing gate; your thresholds narrow them further. ",
              "Score = w \u00d7 (% leukemic cells) + (1-w) \u00d7 (% patients >20% positive).")),
    layout_columns(
      col_widths = c(5, 7),
      card(card_header("Efficacy vs normal-marrow expression"), plotOutput("gate_plot", height = 430),
           note("Dashed lines are your gates; the selected target is ringed in amber.")),
      card(card_header("Passing targets"), DTOutput("hits"),
           note("Click a row to profile that gene in the other tabs."))
    )
  ),

  nav_panel(
    "Target profile",
    layout_columns(
      col_widths = c(7, 5),
      card(card_header(textOutput("ls_title")), uiOutput("ls_slot")),
      card(card_header("Head-to-head with established antigens"), plotOutput("h2h_plot", height = 380),
           note("The 48 antigens profiled side by side. Leads and clinical-stage comparators labelled."))
    ),
    card(card_header("Coverage by cytogenetic subtype"), uiOutput("cyto_slot"),
         note("A target that works across subtypes is more broadly deployable than one confined to a single lesion."))
  ),

  nav_panel(
    "Safety",
    card(
      card_header("Browse the top 500 targets by marrow safety"),
      layout_columns(
        col_widths = c(5, 7),
        plotOutput("safety_scatter", height = 400),
        DTOutput("safety_tbl")),
      note("Marrow safety is available for all 1,704 screened genes: the axes are the percentage of ",
           "normal HSPC and myeloid progenitors expressing the target, the two measures the ",
           "\u226410% gate was built on. Lower-left is safer. Click any row to profile that gene. ",
           strong("Deep Census profiling \u2014 the vital-organ and five-axis panels below \u2014 exists only for the "),
           strong("twelve antigens carried past the screen"), ", so it is flagged in the table rather than implied for the rest.")
    ),
    layout_columns(
      col_widths = c(6, 6),
      card(card_header("Normal marrow cell types"), uiOutput("hema_slot")),
      card(card_header("Vital organs (top 15)"), uiOutput("org_slot"))),
    card(card_header("Five-axis profile"), uiOutput("radar_slot"),
         note("Each axis is min-max scaled across the twelve antigens carried past the screen. Higher is better."))
  ),

  nav_panel(
    "Organ check",
    card(
      card_header("Marrow sparing is not safety"),
      layout_columns(
        col_widths = c(6, 6),
        plotOutput("tradeoff_plot", height = 370),
        div(
          p(strong("The two filters disagree."),
            " All 1,704 screened genes carry marrow-sparing data, because that came from this atlas. ",
            "Body-wide organ expression was profiled for only the ", strong("twelve"),
            " antigens carried past the screen, against the CZ CELLxGENE Census."),
          p("Among those twelve the axes come apart. ", strong("CD9"), " scores 92.4 on marrow sparing and ",
            strong("5.1"), " on organ sparing; ", strong("IL1RAP"), " scores 90.3 and 24.7. Both sit in the top ",
            "four of the marrow ranking, and both were disqualified on body-wide expression. A screen that ",
            "stopped at marrow sparing would have nominated them."),
          p(strong("So a marrow rank is not a nomination."),
            " Seven of the top twelve marrow-ranked genes have never been organ-checked here at all: ",
            "AMN, LMBR1L, NRXN2, FURIN, UMODL1, GPR174 and SUSD1. Check one before you trust it."),
          note("Axes are the five-axis profile, min-max scaled across the twelve profiled antigens. ",
               "Higher is safer. The selected gene is ringed in amber when it is one of the twelve.")
        )
      )
    ),
    layout_columns(
      col_widths = c(5, 7),
      card(card_header(textOutput("organ_status_title")), uiOutput("organ_status")),
      card(card_header("Check this gene body-wide, yourself"), uiOutput("census_panel"))
    ),
    card(
      card_header(textOutput("organ_full_title")),
      uiOutput("organ_full_note"),
      layout_columns(col_widths = c(6, 6),
                     DTOutput("organ_full_tbl"),
                     plotOutput("organ_full_plot", height = 760)),
      note("Every normal tissue the Census profile covers, not a top-15 cut. ",
           strong("% positive"), " is the share of normal cells in that tissue with a detected ",
           "transcript; ", strong("n"), " is how many cells that share is computed from, so a high ",
           "percentage on a few thousand cells is weaker evidence than a low one on a million. ",
           "Tissues flagged ", strong("vital"), " are the ones an antigen-directed therapy cannot spare.")
    ),
    card(
      card_header("All 1,704 screened genes \u2014 marrow data here, organ data one click away"),
      DTOutput("organ_tbl"),
      note("Sort on % HSPC or % myeloid to rank by marrow sparing, then open the Census link to check the ",
           "organs yourself. ", strong("Organ screened = yes"), " means this study profiled the gene across 54 ",
           "normal organ cell types; ", strong("no"), " means the marrow numbers are real but body-wide ",
           "expression is unknown here.")
    )
  ),
  nav_panel(
    "Validation",
    card(card_header("Does the target hold beyond the paediatric discovery cohort?"),
         uiOutput("valid_slot"),
         note("Percent of patients targetable in paediatric TARGET bulk RNA-seq and adult Beat AML, ",
              "alongside percent of cells positive in an independent adult single-cell cohort."))
  ),

  nav_panel(
    "Combinations",
    card(card_header("Coverage vs marrow toxicity"), plotOutput("combo_plot", height = 440),
         note("Each point is a single antigen or an 'X OR Y' pair. Up and left is better.")),
    card(card_header("Pairs involving the selected target"), uiOutput("combo_miss"), DTOutput("combo_tbl"))
  )
)

# --- server ----------------------------------------------------------------
server <- function(input, output, session) {
  g <- reactive(if (is.null(input$gene) || !nzchar(input$gene)) "CD96" else input$gene)

  # narrowing the list keeps the current gene when that gene is in the new list
  observeEvent(input$glist, {
    ch  <- SCOPES[[input$glist]]; cur <- isolate(input$gene)
    updateSelectizeInput(session, "gene", choices = ch,
                         selected = if (!is.null(cur) && cur %in% ch) cur else ch[1])
  }, ignoreInit = TRUE)

  # what this study actually measured for the selected gene
  has_map  <- reactive(g() %in% EXPR_GENES)
  has_ls   <- reactive(has(ls_expr,  g()))
  has_cyto <- reactive(has(cyto,     g()))
  has_h2h  <- reactive(has(h2h,      g()))
  has_tox  <- reactive(has(tox_hema, g()))
  has_org  <- reactive(has(tox_org,  g()))
  has_rad  <- reactive(has(radar,    g()))
  has_val  <- reactive(has(meta_ped, g()) || has(meta_ad, g()) || has(sc_ad, g()))
  has_cmb  <- reactive(any(combos$g1 == g() | combos$g2 == g()))

  ONLY19 <- "Only the 19 antigens taken forward from the screen were profiled this way."
  ONLY12 <- "Deep Census profiling was run on the twelve antigens carried past the screen."

  # A gene with no data for a panel gets the note alone - not a note above an
  # empty 400px canvas, which reads as a failure rather than as "not measured".
  slot <- function(ok, miss, id, h) if (ok) plotOutput(id, height = h) else miss
  output$ls_slot    <- renderUI(slot(has_ls(),   uiOutput("ls_miss"),    "ls_plot",    380))
  output$cyto_slot  <- renderUI(slot(has_cyto(), uiOutput("cyto_miss"),  "cyto_plot",  330))
  output$hema_slot  <- renderUI(slot(has_tox(),  uiOutput("hema_miss"),  "hema_plot",  400))
  output$org_slot   <- renderUI(slot(has_org(),  uiOutput("org_miss"),   "org_plot",   400))
  output$radar_slot <- renderUI(slot(has_rad(),  uiOutput("radar_miss"), "radar_plot", 300))
  output$valid_slot <- renderUI(slot(has_val(),  uiOutput("valid_miss"), "valid_plot", 420))

  output$ls_miss    <- renderUI(if (!has_ls())   elsewhere(g(), "per-state expression", ONLY19))
  output$cyto_miss  <- renderUI(if (!has_cyto()) elsewhere(g(), "cytogenetic breakdown", ONLY19))
  output$hema_miss  <- renderUI(if (!has_tox())  elsewhere(g(), "normal-marrow cell-type profile", ONLY12))
  output$org_miss   <- renderUI(if (!has_org())  elsewhere(g(), "vital-organ profile", ONLY12))
  output$radar_miss <- renderUI(if (!has_rad())  elsewhere(g(), "five-axis profile", ONLY12))
  output$valid_miss <- renderUI(if (!has_val())  elsewhere(g(), "cross-cohort validation", ONLY12))
  output$combo_miss <- renderUI(if (!has_cmb())  elsewhere(g(), "combination analysis",
                                 "Pairs were enumerated among the antigens carried past the screen."))
  output$emb_miss   <- renderUI(
    if (identical(input$emb_fill, "expr") && !has_map())
      elsewhere(g(), "expression on the atlas",
                "Genes detected in fewer than 25 cells are left off the map rather than drawn from noise."))

  # unified screening table: whole cohort, or one poor-prognosis context
  base_tbl <- reactive({
    if (identical(input$ctx, "Whole cohort")) {
      data.frame(gene = screen$gene, aml = screen$scRNA_AML_pct,
                 breadth = screen$scRNA_n_pts_above20pct / N_PT * 100,
                 hspc = screen$scRNA_HSPC_pct, mye = screen$scRNA_Myeloid_pct,
                 npt = screen$scRNA_n_pts_above20pct, paper_rank = screen$rank,
                 stringsAsFactors = FALSE)
    } else {
      d <- lsspec[lsspec$group == input$ctx, ]
      data.frame(gene = d$gene, aml = d$cov_pct, breadth = d$rec_pct,
                 hspc = d$HSPC, mye = d$Myeloid, npt = d$n_contrib_patients,
                 paper_rank = match(d$gene, screen$gene), stringsAsFactors = FALSE)
    }
  })
  scored <- reactive({
    d <- base_tbl()
    d$score <- input$w * d$aml + (1 - input$w) * d$breadth
    d$pass <- d$aml >= input$f_aml & d$breadth >= input$f_pts &
      d$hspc <= input$f_hspc & d$mye <= input$f_mye
    d[order(-d$score), ]
  })

  output$ctx_note <- renderUI({
    if (identical(input$ctx, "Whole cohort"))
      return(note("Screening every leukemic cell in the cohort."))
    r <- lscomp[lscomp$group == input$ctx, ]
    if (!nrow(r)) return(NULL)
    warn <- identical(as.character(r$single_patient_flag[1]), "True")
    note(sprintf("%s: %s cells from %d contributing patients; largest single patient contributes %.1f%%.",
                 input$ctx, format(r$n_cells[1], big.mark = ","), r$n_patients_contrib[1],
                 r$dominant_patient_share_pct[1]),
         if (warn) span(style = "color:#8C1515;font-weight:600", " Dominated by one patient \u2014 treat as an artifact risk."))
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
    gg <- g()
    items <- c("Atlas expression map" = has_map(), "Per leukemic state" = has_ls(),
               "Cytogenetic subtype" = has_cyto(), "Head-to-head panel" = has_h2h(),
               "Normal-marrow safety" = has_tox(), "Vital-organ safety" = has_org(),
               "Cross-cohort validation" = has_val(), "Combination analysis" = has_cmb())
    tagList(
      strong("Data for this gene", style = "font-size:.9rem"),
      tags$ul(class = "list-unstyled", style = "font-size:.82rem;margin:.35rem 0 0",
        lapply(names(items), function(n) tags$li(
          span(style = sprintf("color:%s;font-weight:700", if (items[[n]]) "#2F5D70" else "#C2C8CE"),
               if (items[[n]]) "\u2713 " else "\u2014 "), n))),
      div(style = "font-size:.8rem;margin-top:.6rem", "Look up ", strong(gg), " elsewhere: ",
          gene_links(gg)))
  })

  # ---- Discover ----
  output$funnel_title <- renderText(sprintf("%s \u2014 %s of %s genes pass",
      input$ctx, format(nrow(hits()), big.mark = ","), format(nrow(base_tbl()), big.mark = ",")))
  output$funnel_plot <- renderPlot({
    base <- if (identical(input$ctx, "Whole cohort")) funnel else
      data.frame(stage = paste("screened in", input$ctx), n = nrow(base_tbl()))
    d <- rbind(base, data.frame(stage = "your criteria", n = nrow(hits())))
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
    # scored() carries the context-neutral names (aml/hspc), not screen's raw
    # scRNA_* columns - those only exist on the whole-cohort table.
    ggplot(d, aes(aml, hspc)) +
      geom_point(aes(colour = pass), alpha = if (nrow(sel)) .3 else .6, size = 1.8) +
      geom_vline(xintercept = input$f_aml, linetype = "dashed", colour = CARD) +
      geom_hline(yintercept = input$f_hspc, linetype = "dashed", colour = CARD) +
      { if (nrow(sel)) geom_point(data = sel, size = 10, shape = 21, fill = NA,
                                  colour = "#F6A30C", stroke = 1.5) } +
      { if (nrow(sel)) geom_point(data = sel, size = 4.5, shape = 21, fill = "#F6A30C", colour = "black") } +
      { if (nrow(sel)) geom_text(data = sel, aes(label = gene), vjust = -1.2, fontface = "bold") } +
      scale_colour_manual(values = c(`FALSE` = "#D5D9DD", `TRUE` = TEAL), name = "Passes gates") +
      # headroom on the right so the label on the best target is not clipped
      scale_x_continuous(expand = expansion(mult = c(.02, .1))) +
      labs(x = "% leukemic cells positive", y = "% normal HSPC positive") + theme_lab(0)
  })

  output$hits <- renderDT({
    d <- hits()[, c("gene", "score", "aml", "breadth", "npt", "hspc", "mye", "paper_rank")]
    names(d) <- c("Gene", "Score", "% leukemic cells", "% patients", "Patients", "% HSPC", "% Myeloid", "Paper rank")
    datatable(d, selection = "single", rownames = FALSE,
              options = list(pageLength = 12, order = list(list(1, "desc")))) |>
      formatRound(c("Score", "% leukemic cells", "% patients", "% HSPC", "% Myeloid"), 1)
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
    req(has_ls())
    d <- ls_expr[ls_expr$gene == g(), ]
    d$group <- factor(d$group, levels = d$group[order(-d$pct)])
    ggplot(d, aes(group, pct, fill = prognosis)) + geom_col() +
      scale_fill_manual(values = PROG, name = "State outcome") +
      labs(x = NULL, y = "% cells positive") + theme_lab()
  })

  output$h2h_plot <- renderPlot({
    d <- h2h
    d$is_sel <- d$gene == g()
    ggplot(d, aes(scRNA_AML_pct, scRNA_HSPC_pct, colour = klass)) +
      geom_point(size = 2.6, alpha = if (any(d$is_sel)) .4 else .85) +
      geom_text(data = d[d$do_label == "True" | d$is_sel, ], aes(label = label),
                vjust = -1, size = 3, show.legend = FALSE) +
      { if (any(d$is_sel)) geom_point(data = d[d$is_sel, ], size = 10, shape = 21,
                                      fill = NA, colour = "#F6A30C", stroke = 1.5) } +
      { if (any(d$is_sel)) geom_point(data = d[d$is_sel, ], size = 4.5, shape = 21,
                                      fill = "#F6A30C", colour = "black") } +
      scale_colour_manual(values = KLASS, name = NULL) +
      labs(x = "% leukemic cells positive", y = "% normal HSPC positive") + theme_lab(0)
  })

  output$cyto_plot <- renderPlot({
    req(has_cyto())
    d <- cyto[cyto$gene == g(), ]
    d <- d[order(-d$pct), ]; d$group <- factor(d$group, levels = d$group)
    ggplot(d, aes(group, pct)) +
      geom_col(fill = TEAL) +
      geom_text(aes(label = sprintf("n=%d", n_pt)), vjust = -0.4, size = 3, colour = "#64707C") +
      labs(x = NULL, y = "% cells positive") +
      scale_y_continuous(expand = expansion(mult = c(0, .12))) + theme_lab()
  })

  # ---- Safety ----
  top500 <- local({
    d <- screen[order(-screen$composite), ][1:min(500, nrow(screen)), ]
    d$deep <- ifelse(d$gene %in% radar$gene, "Full Census profile", "Marrow only")
    d
  })
  output$safety_scatter <- renderPlot({
    d <- top500; sel <- d[d$gene == g(), ]
    ggplot(d, aes(scRNA_HSPC_pct, scRNA_Myeloid_pct)) +
      geom_point(aes(colour = composite, shape = deep), size = 2.4,
                 alpha = if (nrow(sel)) .4 else .85) +
      { if (nrow(sel)) geom_point(data = sel, size = 11, shape = 21, fill = NA,
                                  colour = "#F6A30C", stroke = 1.5) } +
      { if (nrow(sel)) geom_point(data = sel, size = 5, shape = 21, fill = "#F6A30C", colour = "black") } +
      { if (nrow(sel)) geom_text(data = sel, aes(label = gene), vjust = -1.2, fontface = "bold") } +
      scale_colour_viridis_c(option = "mako", direction = -1, name = "Composite") +
      scale_shape_manual(values = c(`Full Census profile` = 17, `Marrow only` = 16), name = NULL) +
      labs(x = "% normal HSPC positive", y = "% myeloid progenitors positive") +
      theme_lab(0)
  })
  output$safety_tbl <- renderDT({
    d <- top500[, c("rank", "gene", "composite", "scRNA_AML_pct", "scRNA_HSPC_pct", "scRNA_Myeloid_pct", "deep")]
    names(d) <- c("Rank", "Gene", "Composite", "% AML cells", "% HSPC", "% Myeloid", "Safety data")
    datatable(d, selection = "single", rownames = FALSE,
              options = list(pageLength = 10, order = list(list(0, "asc")))) |>
      formatRound(c("Composite", "% AML cells", "% HSPC", "% Myeloid"), 1)
  })
  observeEvent(input$safety_tbl_rows_selected, {
    updateSelectizeInput(session, "gene", selected = top500$gene[input$safety_tbl_rows_selected])
  })
  output$hema_plot <- renderPlot({
    req(has_tox())
    d <- tox_hema[tox_hema$gene == g(), ]
    d <- d[order(-d$pct), ]; d$cell <- factor(d$cell, levels = rev(d$cell))
    ggplot(d, aes(pct, cell)) + geom_col(fill = TEAL) +
      labs(x = "% cells positive", y = NULL) + theme_lab(0)
  })
  output$org_plot <- renderPlot({
    req(has_org())
    d <- tox_org[tox_org$gene == g(), ]
    d <- head(d[order(-d$pct), ], 15); d$organ <- factor(d$organ, levels = rev(d$organ))
    ggplot(d, aes(pct, organ)) + geom_col(fill = CARD) +
      labs(x = "% cells positive", y = NULL) + theme_lab(0)
  })
  output$radar_plot <- renderPlot({
    req(has_rad())
    d <- radar[radar$gene == g(), ]
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
    req(length(rows) > 0)
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
      labs(x = "Normal-marrow toxicity (%)", y = "% patients with >20% of blasts covered",
           subtitle = if (!any(d$sel))
             sprintf("%s is in none of these pairs, so nothing is highlighted. Pairs were enumerated among the antigens carried past the screen.", g())) +
      theme_lab(0)
  })
  output$combo_tbl <- renderDT({
    req(has_cmb())
    d <- combos[combos$g1 == g() | combos$g2 == g(), c("combo", "kind", "within", "across", "toxicity")]
    names(d) <- c("Combination", "Kind", "Within-patient depth (%)", "Patient breadth (%)", "Marrow toxicity (%)")
    datatable(d, rownames = FALSE, options = list(pageLength = 8, order = list(list(3, "desc")))) |>
      formatRound(3:5, 1)
  })

  # ---- Atlas ----
  output$n_expr_genes <- renderText(format(length(EXPR_GENES), big.mark = ","))
  output$emb_plot <- renderPlot({
    p <- switch(input$emb_fill,
      expr = {
        req(has_map())
        # Colour the cells themselves, not the grid squares. The estimate is still
        # the per-bin mean - a bin holds ~16 cells, which is what makes it stable -
        # but painting it onto the real coordinates keeps the embedding at full
        # resolution instead of quantising it into squares, at no cost in data.
        # Cells where the gene was never detected stay grey, so "not detected"
        # cannot be read as "low", and the ramp is square-root so a marker confined
        # to one population does not vanish into a flat wash.
        m <- emb; m$val <- gene_bins(g())[match(m$bin, ebins$bin)]
        m$val[is.na(m$val)] <- 0
        m <- m[order(m$val), ]                     # expressing cells drawn on top
        z <- m$val <= 0
        ggplot() +
          geom_point(data = m[z, ],  aes(UMAP1, UMAP2), colour = "#E3E7EA", size = .28, alpha = .55) +
          geom_point(data = m[!z, ], aes(UMAP1, UMAP2, colour = val), size = .28, alpha = .85) +
          scale_colour_viridis_c(option = "rocket", direction = -1, trans = "sqrt",
                                 name = sprintf("%s\nmean log(CPM+1)", g())) +
          labs(x = "UMAP1", y = "UMAP2")
      },
      comp = ggplot(emb, aes(UMAP1, UMAP2, colour = compartment)) +
        geom_point(size = .25, alpha = .35) +
        scale_colour_manual(values = c(Leukemic = CARD, `Non-leukemic` = TEAL), name = NULL) +
        guides(colour = guide_legend(override.aes = list(size = 4, alpha = 1))),
      type = ggplot(emb, aes(UMAP1, UMAP2, colour = cell_type)) +
        geom_point(size = .25, alpha = .4) +
        scale_colour_manual(values = grDevices::hcl.colors(length(unique(emb$cell_type)), "Spectral"), name = NULL) +
        guides(colour = guide_legend(override.aes = list(size = 3.5, alpha = 1), ncol = 2)),
      state = {
        d <- emb[emb$LS != "", ]; d$LS <- factor(d$LS, levels = LS_LEVELS)
        lab <- aggregate(cbind(UMAP1, UMAP2) ~ LS, d, median)
        ggplot(d, aes(UMAP1, UMAP2, colour = LS)) +
          geom_point(size = .25, alpha = .45) +
          geom_text(data = lab, aes(label = sub("LS_", "", LS)), colour = "black", size = 3, fontface = "bold") +
          scale_colour_manual(values = grDevices::hcl.colors(49, "Spectral"), guide = "none")
      },
      dens = ggplot(emb, aes(UMAP1, UMAP2)) +
        geom_hex(bins = 90) +
        scale_fill_viridis_c(name = "Cells", option = "mako", trans = "log10"))
    p + coord_equal() + theme_lab(0) + theme(panel.grid = element_blank())
  })

  # --- Organ check -----------------------------------------------------------
  # Marrow sparing is screened for all 1,704 genes; body-wide organ expression
  # only for the twelve carried past the screen. The two disagree sharply, so
  # this tab shows the dissociation and hands the reader the Census lookup for
  # every gene the study could not profile itself.
  CENSUS <- function(g) paste0("https://cellxgene.cziscience.com/gene-expression?genes=", g)

  output$tradeoff_plot <- renderPlot({
    d <- radar
    d$sel <- d$gene == g()
    ggplot(d, aes(Marrow_sparing, Organ_sparing)) +
      geom_hline(yintercept = 50, linetype = 2, colour = "#B9C0C7") +
      { if (any(d$sel)) geom_point(data = d[d$sel, ], size = 10, shape = 21, fill = NA,
                                   colour = "#F6A30C", stroke = 1.6) } +
      geom_point(aes(colour = group), size = 3.4) +
      geom_text(aes(label = label), vjust = -1.15, size = 3.2) +
      scale_colour_manual(values = KLASS, name = NULL) +
      labs(x = "Marrow sparing (score)", y = "Organ sparing (score)") +
      ylim(0, 112) + xlim(40, 104) + theme_lab(0)
  })

  output$organ_status_title <- renderText(sprintf("%s \u2014 organ status in this study", g()))

  output$organ_status <- renderUI({
    gg <- g()
    s  <- screen[screen$gene == gg, ]
    if (has_org()) {
      d <- tox_org[tox_org$gene == gg, ]
      d <- head(d[order(-d$pct), ], 5)
      r <- radar[radar$gene == gg, ]
      tagList(
        div(class = "alert alert-success py-2", style = "font-size:.9rem",
            strong("Profiled body-wide."), sprintf(" Organ-sparing score %.0f of 100.", r$Organ_sparing[1])),
        strong("Highest-expressing normal organ cell types"),
        tags$ul(lapply(seq_len(nrow(d)), function(i)
          tags$li(sprintf("%s \u2014 %.1f%% of cells positive", d$organ[i], d$pct[i])))),
        note("Every tissue is listed below; the five-axis panel is on the Safety tab.")
      )
    } else {
      tagList(
        div(class = "alert alert-warning py-2", style = "font-size:.9rem",
            strong("Not profiled body-wide here."),
            " Marrow sparing below is measured; organ expression is unknown in this study."),
        if (nrow(s)) tags$ul(
          tags$li(sprintf("Marrow rank %s of 1,704", s$rank[1])),
          tags$li(sprintf("%.1f%% of normal HSPC positive", s$scRNA_HSPC_pct[1])),
          tags$li(sprintf("%.1f%% of myeloid progenitors positive", s$scRNA_Myeloid_pct[1])),
          tags$li(sprintf("%.1f%% of leukemic cells positive", s$scRNA_AML_pct[1]))),
        note("Use the Census lookup on the right before treating this as a candidate.")
      )
    }
  })

  output$census_panel <- renderUI({
    gg <- g()
    tagList(
      div(style = "font-size:1.08rem;line-height:1.5;margin-bottom:.85rem",
          "The CZ CELLxGENE Census shows ", strong(gg), " across roughly 65 normal human ",
          "tissues, broken down by cell type, from millions of cells. Use it to see whether ",
          strong(gg), " is expressed anywhere an antigen-directed therapy would have to spare."),
      a(sprintf("Open %s in the CZ CELLxGENE Census \u2192", gg), href = CENSUS(gg),
        target = "_blank", rel = "noopener", class = "btn btn-primary btn-lg",
        style = "font-size:1.3rem;font-weight:700;padding:.75rem 1.4rem;white-space:normal"),
      div(style = "font-size:1.05rem;margin-top:.9rem",
          "Opens in a new tab. Also worth checking: ", gene_links(gg), "."),
      hr(),
      div(style = "font-size:1.05rem;font-weight:700", "What disqualifies a target"),
      tags$ul(style = "font-size:.98rem;margin-bottom:.4rem;line-height:1.5",
        tags$li("Appreciable expression on heart, brain, lung, liver, kidney or pancreas parenchyma \u2014 ",
                "the organs an antigen-directed therapy cannot spare."),
        tags$li("Expression on vascular endothelium, which carries the agent everywhere."),
        tags$li("Expression on normal HSPC or myeloid progenitors, which is the marrow axis the screen ",
                "already applies at a \u226410% gate."),
        tags$li("In this study CD9, IL1RAP and ABCA7 were disqualified on body-wide expression despite ",
                "passing the marrow gate.")),
      note("The Census is a cross-platform relative reference across roughly 65 tissues. It was not ",
           "statistically compared with the AML cohort in this study, and it is not a toxicology ",
           "prediction \u2014 it tells you where a transcript is detected in normal human tissue.")
    )
  })

  # The bundle carries all 54 tissues per profiled antigen; the Safety tab plots
  # only the top 15, which hid 39 of them. This shows every one, sortable, with
  # the cell count behind each percentage so a reader can weigh the evidence.
  VITAL <- c("heart", "brain", "central nervous system", "cortex", "spinal cord",
             "lung", "liver", "kidney", "pancreas", "vasculature", "eye",
             "respiratory system", "sensory system")

  output$organ_full_title <- renderText(
    if (has_org()) sprintf("%s across all 54 normal tissues", g())
    else sprintf("%s \u2014 no body-wide profile in this study", g()))

  output$organ_full_note <- renderUI(
    if (!has_org()) elsewhere(g(), "body-wide tissue profile", ONLY12))

  # The table and plot are declared statically in the UI: a DT nested inside a
  # renderUI on a tab that starts hidden initialises before Shiny binds its data
  # and renders "No matching records found".
  organ_full <- reactive({
    if (!has_org()) return(NULL)
    d <- tox_org[tox_org$gene == g(), c("organ", "pct", "mean", "n")]
    d[order(-d$pct), ]
  })

  output$organ_full_tbl <- renderDT({
    d <- organ_full()
    req(!is.null(d))
    d$vital <- ifelse(d$organ %in% VITAL, "vital", "")
    names(d) <- c("Tissue", "% positive", "Mean expr", "n cells", "Flag")
    datatable(d, rownames = FALSE, selection = "none",
              options = list(pageLength = 15, order = list(list(1, "desc")))) |>
      formatRound("% positive", 1) |> formatRound("Mean expr", 3) |>
      formatRound("n cells", 0)
  })

  output$organ_full_plot <- renderPlot({
    d <- organ_full()
    req(!is.null(d))
    d$organ <- factor(d$organ, levels = rev(d$organ))
    d$vital <- d$organ %in% VITAL
    ggplot(d, aes(pct, organ, fill = vital)) +
      geom_col() +
      scale_fill_manual(values = c(`TRUE` = CARD, `FALSE` = "#B9C0C7"),
                        labels = c(`TRUE` = "vital organ", `FALSE` = "other tissue"),
                        name = NULL) +
      labs(x = "% of normal cells positive", y = NULL) +
      theme_lab(0) + theme(legend.position = "top")
  })

  output$organ_tbl <- renderDT({
    d <- screen[order(screen$rank),
                c("rank", "gene", "scRNA_AML_pct", "scRNA_HSPC_pct", "scRNA_Myeloid_pct")]
    d$organ <- ifelse(d$gene %in% radar$gene, "yes", "no")
    d$census <- sprintf('<a href="%s" target="_blank" rel="noopener">Census &rarr;</a>', CENSUS(d$gene))
    names(d) <- c("Marrow rank", "Gene", "% AML cells", "% HSPC", "% Myeloid",
                  "Organ screened", "Look up")
    datatable(d, selection = "single", rownames = FALSE, escape = FALSE,
              options = list(pageLength = 12, order = list(list(0, "asc")))) |>
      formatRound(c("% AML cells", "% HSPC", "% Myeloid"), 1)
  })
  observeEvent(input$organ_tbl_rows_selected, {
    d <- screen[order(screen$rank), ]
    updateSelectizeInput(session, "gene", selected = d$gene[input$organ_tbl_rows_selected])
  })
  outputOptions(output, "organ_full_tbl", suspendWhenHidden = FALSE)
}

shinyApp(ui, server)
