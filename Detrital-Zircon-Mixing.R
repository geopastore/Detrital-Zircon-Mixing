# ============================================================
#  DETRITAL ZIRCON MIXING — SHINY APP  v3.2
#
#  Folder structure:
#    detrital_mixing_app.R      <- this file
#    precompute_proportions.R   <- run once to build props/
#    props/                     <- pre-computed .rds tables
#
#  Dependencies:
#    install.packages(c("shiny","IsoplotR","DT","plotly","bslib","shinyjs"))
# ============================================================

library(shiny)
library(IsoplotR)
library(DT)
library(plotly)
library(bslib)
library(shinyjs)

# ============================================================
#  VALID CONFIGURATIONS
# ============================================================
VALID_CONFIGS <- list(
  list(n = 2, step = 0.1), list(n = 2, step = 1),
  list(n = 3, step = 1),
  list(n = 4, step = 5),   list(n = 4, step = 10),
  list(n = 5, step = 5),   list(n = 5, step = 10),
  list(n = 6, step = 5),   list(n = 6, step = 10)
)
is_valid_config <- function(n, step) {
  any(vapply(VALID_CONFIGS,
             function(cfg) isTRUE(cfg$n == n && cfg$step == step), logical(1)))
}

# ============================================================
#  CORE FUNCTIONS
# ============================================================

prop_filename <- function(n_sources, step) {
  step <- suppressWarnings(as.numeric(step))
  if (is.null(step) || length(step) == 0 || is.na(step)) return(NULL)
  tag <- gsub("\\.", "", sprintf("%.1f", step))
  file.path("props", sprintf("prop_n%d_s%s.rds", n_sources, tag))
}

load_proportions <- function(n_sources, step, progress_cb = NULL) {
  if (!is_valid_config(n_sources, step))
    stop(sprintf("Invalid configuration: n=%d, step=%.1f%%", n_sources, step))
  fname <- prop_filename(n_sources, step)
  if (!is.null(fname) && file.exists(fname)) {
    if (!is.null(progress_cb)) progress_cb("Loading pre-computed proportions ...", 0.05)
    return(readRDS(fname))
  }
  if (!is.null(progress_cb)) progress_cb("Computing proportions (may take a while) ...", 0.02)
  if (!requireNamespace("combinat", quietly = TRUE))
    stop("Package 'combinat' is required for on-the-fly generation.")
  library(combinat)
  total <- round(100 / step)
  vals  <- seq(0L, as.integer(total))
  combs <- expand.grid(rep(list(vals), n_sources))
  valid <- combs[rowSums(combs) == total, , drop = FALSE]
  get_unique_perms <- function(row)
    unique(do.call(rbind, combinat::permn(as.integer(row))))
  mat <- unique(do.call(rbind, lapply(seq_len(nrow(valid)),
                                      function(i) get_unique_perms(valid[i, ]))))
  storage.mode(mat) <- "integer"
  mat
}

# ---- run_mixing -------------------------------------------------------------
# prop columns correspond positionally to sources list order.
# tgt_col label "Final" and mix label "mix" are safe internal names —
# IsoplotR accesses them via dz$Final / dz$mix using exact name matching.
run_mixing <- function(target, sources, prop, step, n_q = 101, progress_cb = NULL) {
  probs        <- seq(0, 1, length.out = n_q)
  src_q        <- lapply(sources, function(s) quantile(s, probs, na.rm = TRUE))
  tgt_col      <- c("Final", as.character(target))
  n_perm       <- nrow(prop)
  KS           <- numeric(n_perm)
  WS           <- numeric(n_perm)
  update_every <- max(1L, as.integer(n_perm / 200))

  for (k in seq_len(n_perm)) {
    if (!is.null(progress_cb) && (k == 1L || k %% update_every == 0L))
      progress_cb(
        sprintf("Testing mix %d / %d ...", k, n_perm),
        0.10 + 0.88 * (k / n_perm)
      )
    p       <- as.integer(prop[k, ])
    mixed   <- unlist(mapply(rep, src_q, p, SIMPLIFY = FALSE), use.names = FALSE)
    # Guard: if ALL sources have 0 reps (shouldn't happen but be safe)
    if (length(mixed) == 0) { KS[k] <- 1; WS[k] <- 1; next }
    mix_col <- c("mix", as.character(mixed))
    max_len <- max(length(tgt_col), length(mix_col))
    mat <- cbind(
      c(tgt_col, rep("NA", max_len - length(tgt_col))),
      c(mix_col, rep("NA", max_len - length(mix_col)))
    )
    dz    <- IsoplotR::read.data(mat, method = "detritals")
    KS[k] <- IsoplotR::diss(dz[["Final"]], dz[["mix"]], method = "KS")
    WS[k] <- IsoplotR::diss(dz[["Final"]], dz[["mix"]], method = "W2")
  }

  iks <- which.min(KS);  iws <- which.min(WS)
  ks1 <- KS <= quantile(KS, 0.01)
  ws1 <- WS <= quantile(WS, 0.01)
  to_pct <- function(row) as.numeric(row) * step

  # Re-quantile best mix to exactly 150 points for compact export
  make_mix_vec <- function(idx) {
    p   <- as.integer(prop[idx, ])
    raw <- unlist(mapply(rep, src_q, p, SIMPLIFY = FALSE), use.names = FALSE)
    as.numeric(quantile(raw, seq(0, 1, length.out = 150), na.rm = TRUE))
  }

  list(
    KS = KS, WS = WS, n_perm = n_perm, step = step,
    best_KS_pct = to_pct(prop[iks, , drop = TRUE]),
    best_WS_pct = to_pct(prop[iws, , drop = TRUE]),
    avg_KS_pct  = colMeans(prop[ks1, , drop = FALSE]) * step,
    avg_WS_pct  = colMeans(prop[ws1, , drop = FALSE]) * step,
    min_KS      = KS[iks],  min_WS = WS[iws],
    mix_KS      = make_mix_vec(iks),
    mix_WS      = make_mix_vec(iws),
    prop        = prop
  )
}

# ---- Results table ----------------------------------------------------------
make_result_table <- function(res, src_names, sink_name) {
  show_avg <- length(src_names) > 3
  dec      <- if (res$step < 1) 2 else 1
  fmt_row  <- function(label, dist, pct_vec) {
    cbind(
      data.frame(Solution = label,
                 Sink     = sink_name,
                 Distance = if (is.na(dist)) NA_real_ else round(dist, 5),
                 stringsAsFactors = FALSE),
      setNames(as.data.frame(t(round(pct_vec, dec))), paste0(src_names, " (%)"))
    )
  }
  rows <- list(
    fmt_row("Best KS",       res$min_KS, res$best_KS_pct),
    fmt_row("Best WS (W2)",  res$min_WS, res$best_WS_pct)
  )
  if (show_avg) {
    rows[[3]] <- fmt_row("Avg top-1% KS", NA, res$avg_KS_pct)
    rows[[4]] <- fmt_row("Avg top-1% WS", NA, res$avg_WS_pct)
  }
  do.call(rbind, rows)
}

# ---- Colour palette ---------------------------------------------------------
COLOUR_POOL <- c(
  "#1f78b4","#33a02c","#ff7f00","#6a3d9a","#b15928",
  "#a6cee3","#b2df8a","#fdbf6f","#cab2d6","#8dd3c7",
  "#a6761d","#666666","#1b9e77","#d95f02","#7570b3",
  "#e7298a","#66a61e","#e6ab02","#377eb8","#4daf4a"
)

# Named vector: sample_label -> colour
# Sink=black, Mix-KS=red solid, Mix-WS=darkred dashed, sources=distinct colours
make_colour_map <- function(sink_name, src_names, seed = 42) {
  set.seed(seed)
  src_cols <- sample(COLOUR_POOL, min(length(src_names), length(COLOUR_POOL)))
  c(setNames("black",   sink_name),
    setNames("red",     "Mix-KS"),
    setNames("darkred", "Mix-WS"),
    setNames(src_cols,  src_names))
}

# ---- CAD plot ---------------------------------------------------------------
# ---- Build IsoplotR detrital matrix -----------------------------------------
build_detrital_mat <- function(target, mix_KS, mix_WS, sources, sink_name) {
  all_cols <- c(
    list(c(sink_name, as.character(target))),
    list(c("Mix-KS",  as.character(mix_KS))),
    if (!is.null(mix_WS)) list(c("Mix-WS", as.character(mix_WS))) else NULL,
    lapply(names(sources), function(nm) c(nm, as.character(sources[[nm]])))
  )
  max_len <- max(sapply(all_cols, length))
  do.call(cbind, lapply(all_cols, function(col)
    c(col, rep("NA", max_len - length(col)))))
}

# ---- Match colours to IsoplotR's internal name order -----------------------
colours_for_dz <- function(dz, colour_map) {
  nm  <- names(dz)
  col <- colour_map[nm]
  col[is.na(col)] <- "#888888"
  unname(col)
}
draw_cad <- function(target, best_mix, sources, score, method_lbl,
                     sink_name, colour_map) {
  mix_label <- paste0("Mix-", method_lbl)
  all_cols  <- c(
    list(c(sink_name, as.character(target))),
    list(c(mix_label, as.character(best_mix))),
    lapply(names(sources), function(nm) c(nm, as.character(sources[[nm]])))
  )
  max_len <- max(sapply(all_cols, length))
  mat     <- do.call(cbind, lapply(all_cols, function(col)
    c(col, rep("NA", max_len - length(col)))))
  dz      <- IsoplotR::read.data(mat, method = "detritals")

  # Match colours to the actual order IsoplotR stored names in
  nm      <- names(dz)
  col_vec <- unname(ifelse(is.na(colour_map[nm]), "#888888", colour_map[nm]))
  lty_vec <- ifelse(nm == mix_label, 2L, 1L)
  lwd_vec <- ifelse(nm %in% c(sink_name, mix_label), 2, 1.2)

  IsoplotR::cad(dz, col = col_vec, lwd = lwd_vec, lty = lty_vec)
  title(main = sprintf("Best %s = %.5f  |  Sink: %s", method_lbl, score, sink_name),
        cex.main = 0.9)
}

# ---- KDE plot (manual, bypasses IsoplotR colour limitation) ----------------
draw_kde <- function(target, mix_KS, mix_WS, sources, sink_name, colour_map,
                     bw = "SJ") {

  # Build named list of all age vectors to plot
  show_ws  <- !isTRUE(all.equal(mix_KS, mix_WS, tolerance = 1e-8))
  all_samp <- c(
    list(target, mix_KS),
    if (show_ws) list(mix_WS) else NULL,
    sources
  )
  all_names <- c(sink_name, "Mix-KS",
                 if (show_ws) "Mix-WS" else NULL,
                 names(sources))

  col_vec <- unname(ifelse(is.na(colour_map[all_names]), "#888888",
                           colour_map[all_names]))
  lty_vec <- ifelse(all_names == "Mix-WS", 2L, 1L)
  lwd_vec <- ifelse(all_names %in% c(sink_name, "Mix-KS", "Mix-WS"), 2, 1.2)

  # Compute densities
  dens <- lapply(all_samp, function(x) {
    x <- as.numeric(x); x <- x[is.finite(x)]
    density(x, bw = bw, from = max(0, min(x) - 50), to = max(x) + 50, n = 512)
  })

  xlim <- range(sapply(dens, function(d) range(d$x)))
  ylim <- c(0, max(sapply(dens, function(d) max(d$y))) * 1.1)

  plot(NA, xlim = xlim, ylim = ylim,
       xlab = "Age (Ma)", ylab = "Density",
       main = "Kernel Density Estimates", las = 1)

  for (i in seq_along(dens))
    lines(dens[[i]], col = col_vec[i], lty = lty_vec[i], lwd = lwd_vec[i])

  # Rug for each sample
  for (i in seq_along(all_samp)) {
    x <- as.numeric(all_samp[[i]]); x <- x[is.finite(x)]
    rug(x, col = col_vec[i], quiet = TRUE)
  }

  legend("topright", legend = all_names, col = col_vec,
         lty = lty_vec, lwd = lwd_vec, cex = 0.75, bty = "n")
}

# ---- MDS plot ---------------------------------------------------------------
draw_mds <- function(target, mix_KS, mix_WS, sources, sink_name, colour_map) {
  show_ws <- !isTRUE(all.equal(mix_KS, mix_WS, tolerance = 1e-8))

  all_cols <- c(
    list(c(sink_name, as.character(target))),
    list(c("Mix-KS",  as.character(mix_KS))),
    if (show_ws) list(c("Mix-WS", as.character(mix_WS))) else NULL,
    lapply(names(sources), function(nm) c(nm, as.character(sources[[nm]])))
  )
  max_len <- max(sapply(all_cols, length))
  mat     <- do.call(cbind, lapply(all_cols, function(col)
    c(col, rep("NA", max_len - length(col)))))
  dz      <- IsoplotR::read.data(mat, method = "detritals")

  nm      <- names(dz)
  col_vec <- unname(ifelse(is.na(colour_map[nm]), "#888888", colour_map[nm]))

  tryCatch(
    IsoplotR::mds(dz, col = col_vec, nnlines = TRUE),
    error = function(e) {
      plot.new()
      text(0.5, 0.5, paste("MDS failed:", conditionMessage(e)),
           col = "red", cex = 0.9)
    }
  )
}

# ---- 3-D plotly (viridis) ---------------------------------------------------
make_3d_plot <- function(res, src_names, i, j, metric) {
  scores <- if (metric == "KS") res$KS else res$WS
  xi     <- res$prop[, i] * res$step
  xj     <- res$prop[, j] * res$step
  plot_ly(x = xi, y = xj, z = scores,
          type = "scatter3d", mode = "markers",
          marker = list(size = 2, color = scores, colorscale = "Viridis",
                        colorbar = list(title = metric, len = 0.55),
                        showscale = TRUE, opacity = 0.75),
          hovertemplate = paste0(src_names[i], ": %{x:.1f}%<br>",
                                 src_names[j], ": %{y:.1f}%<br>",
                                 metric, ": %{z:.5f}<extra></extra>")) |>
    layout(scene = list(
             xaxis  = list(title = paste0(src_names[i], " (%)")),
             yaxis  = list(title = paste0(src_names[j], " (%)")),
             zaxis  = list(title = metric),
             camera = list(eye = list(x = 1.6, y = 1.6, z = 0.8))),
           margin = list(t = 40, b = 0, l = 0, r = 0),
           title  = list(text = sprintf("%s vs %s", src_names[i], src_names[j]),
                         font = list(size = 12), x = 0.5))
}

# ============================================================
#  UI
# ============================================================

ui <- page_navbar(
  title  = tags$span(
    tags$img(src = "logo.png", height = "42px",
             style = "margin-right:2px; vertical-align:middle;"),
    "Detrital Zircon Mixing"),
  theme  = bs_theme(bootswatch = "flatly",
                    base_font  = font_google("Inter"),
                    "navbar-bg" = "#2c3e50"),
  # useShinyjs() must go in header= not as a direct nav child
  header = useShinyjs(),

  # ==============================================================
  # TAB 1
  # ==============================================================
  nav_panel(
    title = tagList(icon("flask"), " Mixing Analysis"),
    layout_sidebar(
      sidebar = sidebar(
        width = 315, open = TRUE,

        tags$h6("1. Load Data", class = "text-primary fw-bold mt-1 mb-1"),
        fileInput("csv_file", NULL, accept = ".csv",
                  placeholder = "Upload CSV file", width = "100%"),
        helpText("Each column = one sample. Row 1 = column headers (sample names)."),

        hr(class = "my-2"),
        tags$h6("2. Sink (target sample)", class = "text-primary fw-bold mb-1"),
        selectInput("sink_col", NULL, choices = NULL, width = "100%"),

        hr(class = "my-2"),
        tags$h6("3. Sources  (2-6 columns)", class = "text-primary fw-bold mb-1"),
        checkboxGroupInput("src_cols", NULL, choices = NULL),

        hr(class = "my-2"),
        tags$h6("4. Grid Precision", class = "text-primary fw-bold mb-1"),
        uiOutput("precision_ui"),

        hr(class = "my-2"),
        fluidRow(
          column(8, actionButton("run_btn",    "Run Mixing",
                                  class = "btn-primary w-100 fw-bold")),
          column(4, actionButton("cancel_btn", "Cancel",
                                  class = "btn-danger w-100 fw-bold"))
        ),
        br(),
        uiOutput("perm_info_ui"),
        uiOutput("status_badge_ui")
      ),

      mainPanel(
        width = 9,

        # Placeholder
        conditionalPanel(
          "!output.results_ready",
          div(
            class = "d-flex flex-column align-items-center px-5 pt-4",
            style = "min-height:75vh;",
            tags$img(src = "logo.png", height = "150px",
                     style = "margin-right:2px; vertical-align:middle;"),
            #tags$h4("Detrital Zircon Mixing",
            #        class = "mt-2 mb-3 fw-bold text-secondary"),
            div(
              style = "max-width:640px; text-align:left;",
              tags$p(tags$b("What this app does:"),
                " Estimates the proportional contribution of up to 6 sediment source
                regions to a detrital sample (the 'sink'), based on U-Pb zircon
                age distribution similarity."),
              tags$p(tags$b("Method:"),
                " All possible mixtures of source distributions are tested at a chosen
                percentage grid (0.1-10%). Each synthetic mixture is compared to the
                sink using Kolmogorov-Smirnov (KS) and Wasserstein (WS / W2) distances.
                The smallest-distance combination is the best-fit provenance model.
                For >3 sources, the average of the top 1% solutions is also reported."),
              tags$p(tags$b("Outputs:"),
                " Proportions table (with sink label), CAD plots (KS & WS), downloadable
                mix vectors (150-point quantile), KDE, MDS, and interactive 3-D
                score surfaces."),
              tags$p(tags$b("How to start:"),
                " Upload a CSV (one sample per column, first row = sample name).
                Select the sink and 2-6 sources, choose a precision, click Run."),
              hr(),
              tags$p(
                style = "font-size:0.82em; color:#adb5bd; text-align:center; margin-top:8px;",
                "Developed by Guido Pastore",
                tags$br(),
                tags$a("guidopastore93@gmail.com",
                       href = "mailto:guidopastore93@gmail.com",
                       style = "color:#adb5bd;"),
                tags$br(),
                tags$a(icon("github"), " GitHub",
                       href = "https://github.com/guidopastore/detrital-mixing",
                       target = "_blank",
                       style = "color:#adb5bd;")
              )
            )
          )
        ),

        # Results
        conditionalPanel(
          "output.results_ready",

          card(
            card_header(class = "bg-light fw-bold",
                        "Mixing Results - proportions in %"),
            DTOutput("result_table"),
            card_footer(
              class = "bg-white pt-2",
              fluidRow(
                column(6, downloadButton("dl_table",
                                         "Download Results Table (CSV)",
                                         class = "btn-sm btn-outline-secondary w-100")),
                column(6, downloadButton("dl_mix_csv",
                                         "Download KS & WS Mix Vectors (CSV)",
                                         class = "btn-sm btn-outline-secondary w-100"))
              )
            )
          ),

          br(),

          card(
            card_header(
              class = "bg-light",
              fluidRow(
                column(8, tags$b("Cumulative Age Distributions (CAD)")),
                column(2, downloadButton("dl_cad_ks", "KS PDF",
                                         class = "btn-sm btn-outline-primary")),
                column(2, downloadButton("dl_cad_ws", "WS PDF",
                                         class = "btn-sm btn-outline-danger"))
              )
            ),
            fluidRow(
              column(6,
                tags$p("Best KS solution",
                       class = "text-center text-primary fw-semibold mt-2 mb-0"),
                plotOutput("cad_ks", height = "420px")),
              column(6,
                tags$p("Best WS (W2) solution",
                       class = "text-center text-danger fw-semibold mt-2 mb-0"),
                plotOutput("cad_ws", height = "420px"))
            )
          )
        )
      )
    )
  ),

  # ==============================================================
  # TAB 2 — KDE & MDS
  # ==============================================================
  nav_panel(
    title = tagList(icon("chart-area"), " KDE & MDS"),
    conditionalPanel(
      "!output.results_ready",
      div(class = "d-flex flex-column align-items-center justify-content-center",
          style = "height:60vh; color:#adb5bd;",
          tags$span(style = "font-size:72px;", "\U0001f4ca"), br(),
          h5("Run the Mixing Analysis first (Tab 1).", class = "mt-2 text-muted"))
    ),
    conditionalPanel(
      "output.results_ready",
      div(
        class = "px-3 pt-3",
        tags$p("Sink = black | Mix-KS = red | Mix-WS = darkred dashed | Sources = colours",
               class = "text-muted small mb-3"),
        fluidRow(
          column(6, card(
            card_header(class = "bg-light fw-bold", "Kernel Density Estimates"),
            plotOutput("kde_plot", height = "500px"))),
            column(6,
                 card(
                   card_header(class = "bg-light fw-bold", "MDS — KS distance"),
                   plotOutput("mds_ks_plot", height = "400px")),
                 br(),
                 card(
                   card_header(class = "bg-light fw-bold", "MDS — WS (W2) distance"),
                   plotOutput("mds_ws_plot", height = "400px")))
        )
      )
    )
  ),

  # ==============================================================
  # TAB 3 — 3-D Score Surfaces
  # ==============================================================
  nav_panel(
    title = tagList(icon("cube"), " 3-D Score Surfaces"),
    conditionalPanel(
      "!output.results_ready",
      div(class = "d-flex flex-column align-items-center justify-content-center",
          style = "height:60vh; color:#adb5bd;",
          tags$span(style = "font-size:72px;", "\U0001f4e6"), br(),
          h5("Run the Mixing Analysis first (Tab 1).", class = "mt-2 text-muted"))
    ),
    conditionalPanel(
      "output.results_ready",
      div(
        class = "px-3 pt-3",
        tags$p("X/Y = source proportions (%). Z & colour (Viridis) = distance score.
                Drag to rotate, scroll to zoom.",
               class = "text-muted small mb-2"),
        fluidRow(
          column(6, tags$h5("KS distance",    class = "text-center text-primary")),
          column(6, tags$h5("WS (W2) distance", class = "text-center text-danger"))
        ),
        uiOutput("plots_3d_ui")
      )
    )
  )
)

# ============================================================
#  SERVER
# ============================================================

server <- function(input, output, session) {

  rv <- reactiveValues(
    data       = NULL, result    = NULL,
    sources    = NULL, target    = NULL,
    sink_name  = NULL, src_names = NULL,
    colour_map = NULL, running   = FALSE,
    cancelled  = FALSE
  )

  # ---- Load CSV ---------------------------------------------------------------
  observeEvent(input$csv_file, {
    req(input$csv_file)
    tryCatch({
      df <- read.csv(input$csv_file$datapath, header = TRUE,
                     stringsAsFactors = FALSE, check.names = FALSE)
      rv$data      <- df;  rv$result    <- NULL
      rv$sources   <- NULL; rv$target   <- NULL
      rv$sink_name <- NULL; rv$src_names <- NULL
      rv$colour_map <- NULL; rv$running  <- FALSE; rv$cancelled <- FALSE
      cols <- colnames(df)
      updateSelectInput(session, "sink_col", choices = cols, selected = cols[1])
      updateCheckboxGroupInput(session, "src_cols", choices = cols,
                               selected = if (length(cols) > 1) cols[-1] else character(0))
    }, error = function(e)
      showNotification(paste("CSV error:", e$message), type = "error", duration = 8))
  })

  # ---- Cancel -----------------------------------------------------------------
  observeEvent(input$cancel_btn, {
    if (isTRUE(rv$running)) {
      rv$cancelled <- TRUE
      showNotification("Cancelling after current iteration...", type = "warning", duration = 4)
    }
  })

  # ---- Precision UI -----------------------------------------------------------
  # Default for 2 sources is NOW 1% (not 0.1%)
  output$precision_ui <- renderUI({
    n <- length(input$src_cols)
    if (n < 2) return(helpText("Select at least 2 sources first."))
    if (n == 2) {
      choices  <- c("1%  (~101 rows)" = "1", "0.1%  (~1001 rows, slow)" = "0.1")
      selected <- "1"            # <-- default is 1%, 0.1% available but not default
      note     <- "1% recommended; 0.1% for higher precision."
    } else if (n == 3) {
      choices  <- c("1%  (~5151 rows)" = "1")
      selected <- "1";  note <- "1% grid for 3 sources."
    } else {
      choices  <- c("5%  (recommended)" = "5", "10%  (fast)" = "10")
      selected <- "5"
      note     <- sprintf("For %d sources at 5%% - top-1%% average also shown.", n)
    }
    # Mark missing .rds files (only modify label names, never values)
    for (v in names(choices)) {
      fname <- prop_filename(n, as.numeric(v))
      if (is.null(fname) || !file.exists(fname))
        names(choices)[choices == v] <- paste0(names(choices)[choices == v], " [will compute]")
    }
    tagList(radioButtons("step", NULL, choices = choices, selected = selected),
            helpText(note))
  })

  # ---- Permutation count info -------------------------------------------------
  output$perm_info_ui <- renderUI({
    n <- length(input$src_cols)
    if (n < 2 || is.null(input$step) || nchar(trimws(input$step)) == 0) return(NULL)
    step <- suppressWarnings(as.numeric(trimws(input$step)))
    if (is.na(step) || step <= 0 || !is_valid_config(n, step)) return(NULL)
    fname <- prop_filename(n, step)
    if (!is.null(fname) && file.exists(fname)) {
      n_rows <- nrow(readRDS(fname))
      badge  <- tags$span(class = "badge bg-success ms-1", "cached")
    } else {
      total  <- round(100 / step)
      n_rows <- choose(total + n - 1, n - 1) * factorial(n)
      badge  <- tags$span(class = "badge bg-warning text-dark ms-1", "will compute")
    }
    div(class = "text-muted small mt-1 mb-2",
        sprintf("~%s permutations  ", formatC(n_rows, format = "fg", big.mark = ",")), badge)
  })

  # ---- Status badge -----------------------------------------------------------
  output$status_badge_ui <- renderUI({
    if (isTRUE(rv$running))
      div(class = "alert alert-info p-2 small mt-2 mb-0", "Running...")
    else if (isTRUE(rv$cancelled))
      div(class = "alert alert-warning p-2 small mt-2 mb-0", "Cancelled.")
    else if (!is.null(rv$result))
      div(class = "alert alert-success p-2 small mt-2 mb-0",
          sprintf("Done - %s permutations tested.",
                  formatC(rv$result$n_perm, format = "fg", big.mark = ",")))
  })

  # ---- Run mixing -------------------------------------------------------------
  observeEvent(input$run_btn, {
    req(rv$data, input$sink_col, input$src_cols)

    sink <- input$sink_col
    # Force source order to match CSV column order (not checkbox tick order)
    all_cols_csv <- colnames(rv$data)
    srcs <- all_cols_csv[all_cols_csv %in% input$src_cols]

    if (is.null(input$step) || nchar(trimws(input$step)) == 0) {
      showNotification("Please select a grid precision.", type = "warning"); return()
    }
    step <- suppressWarnings(as.numeric(trimws(input$step)))

    if (length(srcs) < 2)  { showNotification("Select >= 2 sources.", type = "warning"); return() }
    if (length(srcs) > 6)  { showNotification("Maximum 6 sources.",   type = "warning"); return() }
    if (sink %in% srcs)    { showNotification("Sink must differ from sources.", type = "warning"); return() }
    if (is.na(step) || step <= 0 || !is_valid_config(length(srcs), step)) {
      showNotification("Invalid or unsupported (n, step) combination.", type = "error"); return()
    }

    target  <- { x <- as.numeric(rv$data[[sink]]); x[!is.na(x)] }
    sources <- lapply(srcs, function(nm) { x <- as.numeric(rv$data[[nm]]); x[!is.na(x)] })
    names(sources) <- srcs

    if (length(target) == 0 || any(vapply(sources, length, integer(1)) == 0)) {
      showNotification("One or more columns are entirely NA.", type = "error"); return()
    }

    rv$running    <- TRUE
    rv$cancelled  <- FALSE
    rv$result     <- NULL
    rv$colour_map <- make_colour_map(sink, srcs)

    withProgress(message = "Starting...", value = 0, {
      prop <- tryCatch(
        load_proportions(length(srcs), step,
                         progress_cb = function(msg, frac) incProgress(frac, message = msg)),
        error = function(e) {
          showNotification(paste("Proportion error:", e$message), type = "error", duration = 10)
          NULL
        })
      if (is.null(prop)) { rv$running <- FALSE; return() }

      incProgress(0.10, message = sprintf("Loaded %s permutations. Running mixing...",
                                          formatC(nrow(prop), format = "fg", big.mark = ",")))
      result <- tryCatch(
        run_mixing(target, sources, prop, step,
                   progress_cb = function(msg, frac) incProgress(frac - 0.10, message = msg)),
        error = function(e) {
          showNotification(paste("Mixing error:", e$message), type = "error", duration = 10)
          NULL
        })
    })

    rv$running <- FALSE
    if (is.null(result) || isTRUE(rv$cancelled)) return()

    rv$result    <- result
    rv$sources   <- sources
    rv$target    <- target
    rv$sink_name <- sink
    rv$src_names <- srcs
  })

  # ---- results_ready ----------------------------------------------------------
  output$results_ready <- reactive({ !is.null(rv$result) })
  outputOptions(output, "results_ready", suspendWhenHidden = FALSE)

  # ---- Results table ----------------------------------------------------------
  output$result_table <- renderDT({
    req(rv$result)
    tbl      <- make_result_table(rv$result, rv$src_names, rv$sink_name)
    pct_cols <- paste0(rv$src_names, " (%)")
    datatable(
      tbl, rownames = FALSE, selection = "none",
      options = list(dom = "t", pageLength = 4, ordering = FALSE,
                     columnDefs = list(list(className = "dt-center",
                                           targets = seq_len(ncol(tbl)) - 1L))),
      class = "compact stripe hover cell-border"
    ) |>
      formatStyle("Solution", fontWeight = "bold", textAlign = "left") |>
      formatStyle("Sink",     fontWeight = "bold", color = "black") |>
      formatStyle("Distance",
                  background = styleColorBar(c(0, 0.5), "#cce5ff"),
                  backgroundSize = "95% 65%", backgroundRepeat = "no-repeat",
                  backgroundPosition = "center") |>
      formatStyle(pct_cols,
                  background = styleColorBar(c(0, 100), "#d4edda"),
                  backgroundSize = "95% 65%", backgroundRepeat = "no-repeat",
                  backgroundPosition = "center")
  })

  # ---- CAD plots --------------------------------------------------------------
  output$cad_ks <- renderPlot({
    req(rv$result, rv$colour_map)
    draw_cad(rv$target, rv$result$mix_KS, rv$sources,
             rv$result$min_KS, "KS", rv$sink_name, rv$colour_map)
  }, bg = "white")

  output$cad_ws <- renderPlot({
    req(rv$result, rv$colour_map)
    draw_cad(rv$target, rv$result$mix_WS, rv$sources,
             rv$result$min_WS, "WS", rv$sink_name, rv$colour_map)
  }, bg = "white")

  # ---- CAD PDF downloads ------------------------------------------------------
  output$dl_cad_ks <- downloadHandler(
    filename = function() sprintf("CAD_KS_%s_%s.pdf", rv$sink_name, Sys.Date()),
    content  = function(file) {
      pdf(file, width = 8, height = 6)
      draw_cad(rv$target, rv$result$mix_KS, rv$sources,
               rv$result$min_KS, "KS", rv$sink_name, rv$colour_map)
      dev.off()
    })

  output$dl_cad_ws <- downloadHandler(
    filename = function() sprintf("CAD_WS_%s_%s.pdf", rv$sink_name, Sys.Date()),
    content  = function(file) {
      pdf(file, width = 8, height = 6)
      draw_cad(rv$target, rv$result$mix_WS, rv$sources,
               rv$result$min_WS, "WS", rv$sink_name, rv$colour_map)
      dev.off()
    })

  # ---- KDE plot ---------------------------------------------------------------
  output$kde_plot <- renderPlot({
    req(rv$result, rv$colour_map)
    #ks_eq_ws <- isTRUE(all.equal(rv$result$mix_KS, rv$result$mix_WS, tolerance = 1e-8))
    mat <- build_detrital_mat(rv$target, rv$result$mix_KS,
                              #if (ks_eq_ws) NULL else 
                              rv$result$mix_WS,
                              rv$sources, rv$sink_name)
    dz        <- IsoplotR::read.data(mat, method = "detritals")
    col_vec   <- colours_for_dz(dz, rv$colour_map)
    named_col <- setNames(col_vec, names(dz))
    IsoplotR::kde(dz, col = named_col, rug = TRUE, show.hist = FALSE, normalise = FALSE)
  }, bg = "white")

  # ---- MDS plots (KS distance | WS distance) ----------------------------------
  make_mds_plot <- function(mix_vec, mix_label, dist_method) {
    all_cols <- c(
      list(c(rv$sink_name, as.character(rv$target))),
      list(c(mix_label,    as.character(mix_vec))),
      lapply(names(rv$sources), function(nm) c(nm, as.character(rv$sources[[nm]])))
    )
    max_len <- max(sapply(all_cols, length))
    mat     <- do.call(cbind, lapply(all_cols, function(col)
      c(col, rep("NA", max_len - length(col)))))
    dz      <- IsoplotR::read.data(mat, method = "detritals")
    col_vec <- colours_for_dz(dz, rv$colour_map)
    tryCatch(
      IsoplotR::mds(dz, col = col_vec, nnlines = TRUE, method = dist_method),
      error = function(e) {
        plot.new()
        text(0.5, 0.5, paste("MDS failed:", conditionMessage(e)), col = "red", cex = 0.9)
      }
    )
    title(main = paste("MDS —", dist_method, "distance"), cex.main = 0.95)
  }
  
  output$mds_ks_plot <- renderPlot({
    req(rv$result, rv$colour_map)
    make_mds_plot(rv$result$mix_KS, "Mix-KS", "KS")
  }, bg = "white")
  
  output$mds_ws_plot <- renderPlot({
    req(rv$result, rv$colour_map)
    make_mds_plot(rv$result$mix_WS, "Mix-WS", "W2")
  }, bg = "white")

  # ---- 3-D plots --------------------------------------------------------------
  output$plots_3d_ui <- renderUI({
    req(rv$result, rv$src_names)
    pairs <- combn(seq_len(length(rv$src_names)), 2, simplify = FALSE)
    rows  <- lapply(seq_along(pairs), function(k)
      fluidRow(class = "mb-2",
               column(6, plotlyOutput(paste0("p3d_KS_", k), height = "420px")),
               column(6, plotlyOutput(paste0("p3d_WS_", k), height = "420px"))))
    do.call(tagList, rows)
  })

  observe({
    req(rv$result, rv$src_names)
    if (ncol(rv$result$prop) != length(rv$src_names)) return()
    pairs <- combn(seq_len(length(rv$src_names)), 2, simplify = FALSE)
    for (k in seq_along(pairs)) {
      local({
        kk <- k; ii <- pairs[[kk]][1]; jj <- pairs[[kk]][2]
        output[[paste0("p3d_KS_", kk)]] <- renderPlotly(
          make_3d_plot(rv$result, rv$src_names, ii, jj, "KS"))
        output[[paste0("p3d_WS_", kk)]] <- renderPlotly(
          make_3d_plot(rv$result, rv$src_names, ii, jj, "WS"))
      })
    }
  })

  # ---- Downloads --------------------------------------------------------------
  output$dl_table <- downloadHandler(
    filename = function() sprintf("mixing_results_%s.csv", Sys.Date()),
    content  = function(f)
      write.csv(make_result_table(rv$result, rv$src_names, rv$sink_name), f, row.names = FALSE))

  output$dl_mix_csv <- downloadHandler(
    filename = function() sprintf("mix_vectors_%s_%s.csv", rv$sink_name, Sys.Date()),
    content  = function(f)
      write.csv(data.frame(KS_mix = rv$result$mix_KS,
                           WS_mix = rv$result$mix_WS), f, row.names = FALSE))
}

# ============================================================
shinyApp(ui, server)
