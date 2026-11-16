## ─────────────────────────────────────────────────────────────────────────────
## server/smlm/21_density_normalisation.R  -  labelling-density profile and thinning
## ─────────────────────────────────────────────────────────────────────────────
## Part of smlm_server(). Not a standalone script: its statements are evaluated in
## the server function's environment, in file-name order (server/parts.R).
##
## Measuring each file's density, choosing the common target, the
## keep-probabilities, and thinning one file's units. The thinning changes the
## data of every batch panel.
##
## INTERFACE (generated: Rscript tools/server_map.R --headers; do not edit by hand)
##   DEFINES for other parts:
##     apply_density_thinning, build_density_map_fig, build_density_plot_fig,
##     build_density_sweep_plot_fig, density_p_map
##   READS from other parts:
##     01_view_upload: smlm_file_store
##     20_dataset_groups: apply_dataset_groups, file_group_map,
##                        series_palette
##     27_fiducial_qc: qc_exclude_from
## END INTERFACE
## ─────────────────────────────────────────────────────────────────────────────

## ══════════════════════════════════════════════════════════════════════════
## LABELLING-DENSITY NORMALISATION  (opt-in, default OFF)
## ══════════════════════════════════════════════════════════════════════════
## GOAL: make clustering read-outs comparable between acquisitions whose
## LABELLING DENSITY differs for reasons that are not biology - transient
## transfection giving each cell an undefined plasmid load, different expression
## level, run-to-run labelling efficiency.
##
## Two independent opt-ins, both off by default because most experiments do not
## need them (a stable line, or one acquisition, has no such confound):
##   1. `density_profile_on`  - just MEASURE and report per-file density.
##   2. `density_thin_on`     - additionally THIN every file to one common
##                              density before any analysis sees it.
##
## The maths, the choice of thinning target and the trap in using NN distances
## as that target are all documented at the "8b. LABELLING-DENSITY
## NORMALISATION" banner in R/smlm_support.R. Read that before changing anything
## here. In one line: g(r) and K(r) are exactly invariant under independent
## thinning, so thinning removes a density artefact without removing real
## clustering - but the TARGET must be the global density n/support, which
## clustering cannot inflate, not the NN-based local intensity, which it can
## (measured 4.25x on a clustered pattern).
density_profile_rv <- reactiveVal(NULL)   # one row per file (see compute_unit_density)
density_sweep_rv   <- reactiveVal(NULL)   # support vs closing radius, per file
density_local_rv   <- reactiveVal(NULL)   # long df(dataset, local) - per-UNIT density

# Which unit the density is measured and thinned on. De-blinked molecules
# whenever de-blinking is on: a trace ID is ONE BLINK, so trace-ID density
# mixes expression level with photophysics and laser power and is not
# comparable between runs.
density_unit <- reactive({
  u <- input$dens_norm_unit %||% SMLM_DEFAULTS$dens_norm_unit
  if (u != "auto") return(u)
  if (isTRUE(input$enable_clustering)) "cluster" else "tid"
})

# Master gates. Thinning implies profiling (it needs the per-file densities).
density_profile_on <- reactive(isTRUE(input$dens_norm_enable))
density_thin_on    <- reactive(isTRUE(input$dens_norm_enable) && isTRUE(input$dens_thin_enable))

## ── Measure: build the per-file density profile ───────────────────────────
observeEvent(input$run_density, {
  run_set <- smlm_settings_from_input(input)   # one settings snapshot for this run
  store <- smlm_file_store()
  if (!length(store)) { showNotification("No files loaded.", type = "warning"); return() }

  unit    <- density_unit()
  use_3d  <- isTRUE(input$dens_norm_3d)
  smode   <- input$dens_support_mode %||% "grid"
  h_nm    <- max(10,  as.numeric(input$dens_grid_h     %||% 100))
  close_n <- max(0,   as.numeric(input$dens_close_nm   %||% 200))
  kk      <- max(1L,  as.integer(input$dens_nn_k       %||% 3L))
  kde_bw  <- max(10,  as.numeric(input$dens_kde_bw     %||% 150))
  kde_fr  <- min(1, max(0.5, as.numeric(input$dens_kde_frac %||% 0.99)))
  needs_cluster <- unit == "cluster"

  auto_close <- isTRUE(input$dens_close_auto) && smode == "grid"
  want_sweep <- auto_close || (isTRUE(input$dens_show_sweep) && smode == "grid")

  prepped <- list()
  out <- withProgress(message = "Measuring labelling density...", value = 0, {
    nf <- length(store)
    # Filter (and de-blink) each file. NOT thinned: this is the pre-thinning
    # reference the whole panel is built on.
    for (i in seq_along(store)) {
      fname <- names(store)[[i]]
      incProgress(0.6 / nf, detail = fname)
      df <- qc_exclude_from(store[[fname]], fname)
      if (!is.data.frame(df) || !nrow(df)) next
      df <- tryCatch(apply_meta_filters_batch(df, run_set), error = function(e) NULL)
      if (is.null(df) || !nrow(df)) next
      if (needs_cluster) {
        df <- tryCatch(apply_clustering_pipeline(df, run_set), error = function(e) df)
      }
      prepped[[fname]] <- df
    }
    density_profile_files(prepped, unit = unit, use_3d = use_3d, support_mode = smode,
                          h_nm = h_nm, close_nm = close_n, auto_close = auto_close,
                          want_sweep = want_sweep, k = kk,
                          support_from = input$dens_support_from %||% "units",
                          bw_nm = kde_bw, frac = kde_fr,
                          progress = function(frac, detail) incProgress(frac, detail = detail))
  })
  res <- out$profile
  if (is.null(res) || !nrow(res)) {
    showNotification("Could not measure density on any file (need >=3 units per file).",
                     type = "warning")
    density_profile_rv(NULL); return()
  }
  density_profile_rv(res)
  density_sweep_rv(out$sweeps)
  density_local_rv(out$locals)
  rr <- range(res$density, na.rm = TRUE)
  showNotification(sprintf(
    "Density measured on %d file(s): %.1f - %.1f units/%s (%.1fx spread).",
    nrow(res), rr[1], rr[2], if (use_3d) "um3" else "um2",
    if (rr[1] > 0) rr[2] / rr[1] else NA_real_),
    type = "message", duration = 5)
})

## ── Target: the one density every file is thinned to ──────────────────────
# Returned in the same units as density_profile_rv()$density. The target must be
# REACHABLE by thinning alone: nothing can be up-sampled, so a file already below
# the target is left untouched and flagged rather than silently misrepresented.
density_target <- reactive({
  prof <- density_profile_rv()
  if (is.null(prof) || !nrow(prof)) return(NULL)
  mode <- input$dens_target_mode %||% "min"
  if (mode == "manual") {
    v <- suppressWarnings(as.numeric(input$dens_target_manual %||% NA))
    return(if (is.finite(v) && v > 0) v else NULL)
  }
  if (mode == "control" && (input$analysis_by %||% "file") == "group") {
    gm <- file_group_map()
    ctl <- input$dens_target_group %||% ""
    if (nzchar(ctl) && length(gm)) {
      fs <- names(gm)[gm == ctl]
      v  <- stats::median(prof$density[prof$dataset %in% fs], na.rm = TRUE)
      if (is.finite(v) && v > 0) return(v)
    }
    return(NULL)
  }
  # "min" (default): the lowest density in the comparison set. Requires no
  # up-sampling of any file, so it is the only target that is always reachable.
  scope <- prof
  if ((input$analysis_by %||% "file") == "group") {
    gm <- file_group_map()
    if (length(gm)) scope <- prof[prof$dataset %in% names(gm), , drop = FALSE]
  }
  if (!nrow(scope)) return(NULL)
  v <- min(scope$density, na.rm = TRUE)
  if (is.finite(v) && v > 0) v else NULL
})

# Per-file keep-probability p = target / observed, capped at 1.
density_p_map <- reactive({
  prof <- density_profile_rv(); tgt <- density_target()
  if (is.null(prof) || is.null(tgt)) return(NULL)
  p <- pmin(1, tgt / prof$density)
  stats::setNames(p, prof$dataset)
})

## ── The gate: thin one file's localizations ───────────────────────────────
# Mirrors qc_exclude_from's contract - call it on a df together with the file it
# came from, and it returns the df unchanged whenever the feature is off. Must be
# applied AFTER the metadata filters and BEFORE clustering, because making the
# fixed-(eps, minPts) DBSCAN comparable is the entire point.
#
# The seed is per-file and fixed, so a given file always thins to the same subset
# (results are reproducible across re-runs), and thin_units_to_p saves/restores
# the global RNG stream so no other randomised analysis is perturbed.
# `settings`: the settings object the caller's chain runs with (prepared_df
# passes its own), so a molecule-level thinning de-blinks with the same values.
apply_density_thinning <- function(df, fname, settings = input) {
  if (!density_thin_on() || is.null(df) || !nrow(df)) return(df)
  pm <- density_p_map()
  if (is.null(pm) || !fname %in% names(pm)) return(df)
  p <- pm[[fname]]
  if (!is.finite(p) || p >= 1) return(df)
  unit <- density_unit()
  # Molecule-level thinning needs the de-blinked IDs to exist. De-blinking is a
  # deterministic DBSCAN, so deriving it here and again downstream on the
  # surviving points gives the same molecule assignment.
  if (unit == "cluster" && !"cluster" %in% names(df)) {
    df <- tryCatch(apply_clustering_pipeline(df, settings), error = function(e) df)
  }
  seed <- as.integer(sum(utf8ToInt(fname)) %% 100000L) + 1L
  thin_units_to_p(df, unit = unit, p = p, seed = seed)
}

## ── Readouts ──────────────────────────────────────────────────────────────
output$dens_norm_active <- reactive(density_thin_on())
outputOptions(output, "dens_norm_active", suspendWhenHidden = FALSE)

output$density_status <- renderUI({
  if (!density_profile_on())
    return(helpText("Off - analyses run on the data as acquired."))
  prof <- density_profile_rv()
  if (is.null(prof) || !nrow(prof))
    return(helpText("Press 'Measure density' to build the per-file profile."))
  un  <- if (isTRUE(input$dens_norm_3d)) "units/µm³" else "units/µm²"
  ul  <- if (prof$unit[1] == "cluster") "de-blinked molecules" else "trace IDs"
  rr  <- range(prof$density, na.rm = TRUE)
  tgt <- density_target()
  spread <- if (rr[1] > 0) rr[2] / rr[1] else NA_real_

  # Headline: the spread IS the confound. Everything else is a flag, shown only
  # when it fires, so the panel stays quiet when the measurement is clean.
  lines <- list(tags$div(style = "font-weight:600;", sprintf(
    "%d file(s) · %.1f–%.1f %s · %.1f× spread · unit = %s",
    nrow(prof), rr[1], rr[2], un, spread, ul)))
  # From the PROFILE, so the description matches the measurement.
  smode_now <- (if ("run_mode" %in% names(prof)) prof$run_mode[1] else NULL) %||%
               (input$dens_support_mode %||% "kde")
  if (smode_now == "grid")
    lines <- c(lines, list(tags$div(style = "color:#777;", sprintf(
      "Support: %g nm grid, gap fill %g nm (%s) - dark regions wider than ~%g nm stay excluded.",
      prof$h_nm[1], prof$close_nm[1],
      if (isTRUE(prof$close_auto[1])) "auto" else "manual", 2 * prof$close_nm[1]))))
  if (smode_now == "kde")
    lines <- c(lines, list(tags$div(style = "color:#777;", sprintf(
      "Support: smoothed density, region enclosing %.0f%% of units, %g nm bandwidth on a %g nm grid - dark regions wider than ~%g nm fall below the level.",
      100 * ((if ("run_frac" %in% names(prof)) prof$run_frac[1] else NULL) %||% 0.99),
      (if ("run_bw" %in% names(prof)) prof$run_bw[1] else NULL) %||% 150,
      prof$h_nm[1],
      6 * ((if ("run_bw" %in% names(prof)) prof$run_bw[1] else NULL) %||% 150)))))

  warn <- function(txt) tags$div(style = "color:#b26a00;", txt)
  if (is.finite(spread) && spread < 1.5)
    lines <- c(lines, list(tags$div(style = "color:#1a8a1a;",
      "Already comparable (< 1.5×) - thinning would only cost points.")))

  # THE failure mode of the support estimate: with too small a gap-fill radius the
  # occupancy grid outlines the STRUCTURES instead of the cell, so the density
  # denominator absorbs the clustering it is meant to be independent of - and by a
  # different amount per file. Measured: 200 nm gave a 3.13 ratio where truth was
  # 4.0; the plateau gave 3.82.
  if (smode_now %in% c("grid", "kde")) {
    fb <- prof$frac_of_bbox
    if (any(is.finite(fb) & fb < 0.5))
      lines <- c(lines, list(warn(sprintf(
        "Support is only %.0f–%.0f%% of the frame. If your sample fills more than that, the smoothing is too tight and the support is outlining your structures, not the cell - check the support map.",
        100 * min(fb, na.rm = TRUE), 100 * max(fb, na.rm = TRUE)))))
    if (any(is.finite(fb) & fb > 0.97))
      lines <- c(lines, list(warn(
        "Support has reached the bounding box - dark regions are being counted as sample. Tighten the smoothing (lower bandwidth, or a lower enclosed fraction).")))
    if (sum(is.finite(fb)) >= 2 && diff(range(fb, na.rm = TRUE)) > 0.25)
      lines <- c(lines, list(warn(
        "Support fraction differs by > 25 points between files, so the denominators are not the same kind of region.")))
  }
  if (density_thin_on()) {
    pm <- density_p_map()
    if (is.null(pm) || is.null(tgt)) {
      lines <- c(lines, list(warn("Thinning on but no target resolved - check the target setting.")))
    } else {
      worst <- if (length(pm)) min(pm, na.rm = TRUE) else NA_real_
      lines <- c(lines, list(tags$div(style = "font-weight:600; color:#c03030;", sprintf(
        "Thinning ON → %.1f %s · %d/%d thinned · densest keeps %.0f%%",
        tgt, un, sum(pm < 1), length(pm), 100 * worst))))
      if (is.finite(worst) && worst < 0.2)
        lines <- c(lines, list(warn(sprintf(
          "Only %.0f%% of the densest file kept. g(r) survives this, but count-based read-outs lose precision and small clusters can fall below minPts.",
          100 * worst))))
    }
  }
  do.call(tags$div, c(list(style = "font-size:0.85em; line-height:1.5;"), lines))
})

output$density_table <- DT::renderDataTable({
  prof <- density_profile_rv()
  shiny::validate(shiny::need(!is.null(prof) && nrow(prof) > 0,
                              "Press 'Measure density'."))
  pm  <- density_p_map()
  d3  <- isTRUE(input$dens_norm_3d)
  out <- data.frame(
    File           = prof$dataset,
    Units          = prof$n_units,
    `Support (um2/3)` = signif(prof$support_nm / 1000^prof$dims, 4),
    `Support/bbox` = round(prof$frac_of_bbox, 3),
    Density        = round(prof$density, 2),
    `Density (NN)` = round(prof$density_nn, 2),
    `Median NN (nm)` = round(prof$nn_median, 1),
    `Keep p`       = if (!is.null(pm)) round(unname(pm[prof$dataset]), 3) else NA_real_,
    check.names    = FALSE, stringsAsFactors = FALSE)
  DT::datatable(out, rownames = FALSE, options = list(dom = "t", paging = FALSE,
                scrollX = TRUE), class = "compact stripe") %>%
    DT::formatStyle("Keep p", color = DT::styleInterval(c(0.2), c("#b26a00", "#333")))
})

## Per-file density DISTRIBUTION: one horizontal violin per dataset over the
## per-unit local density, with that file's global n/support marked on it.
##
## Why a distribution and not one dot per file: the global number says how much
## label there is, the shape says how it is arranged. A cell with real aggregates
## is heterogeneous (sparse background mode + high-density tail); a uniformly
## labelled cell at the same mean is narrow and unimodal. So the width here is
## itself clustering information - and it is also the visual reason the mean of
## this distribution must not be the thinning target (see nn_local_intensity).
##
## Series are identified by the y-axis tick labels, so no legend is needed - which
## is what keeps the panel readable with several files.
##
## AXIS: the x range is set EXPLICITLY from a robust quantile span of the data
## (plus the global-density markers), never left to plotly's autorange - a violin's
## KDE tails run well past the data and, on a log axis, autorange then pads out to
## whole decades, which is what left the panel long and mostly empty. Log is
## offered because local density is a ratio-scale quantity that can span decades
## inside one cell, and density comparisons are multiplicative (a 4x difference
## should look the same anywhere on the axis) - but it is a toggle, since a sample
## with a narrow spread reads better linear.
build_density_plot_fig <- function() {
  prof <- density_profile_rv()
  shiny::validate(shiny::need(!is.null(prof) && nrow(prof) > 0, "Press 'Measure density'."))
  loc  <- density_local_rv()
  grouped <- (input$analysis_by %||% "file") == "group"
  # dims from the profile: the 2D/3D switch is a compute parameter.
  un   <- if (isTRUE((prof$dims[1] %||% 2L) == 3L)) "units/µm³" else "units/µm²"

  pr <- apply_dataset_groups(prof)
  shiny::validate(shiny::need(nrow(pr) > 0, "No files assigned to groups (group mode)."))
  lc <- if (!is.null(loc)) apply_dataset_groups(loc) else NULL

  # Order rows by the group/file median so the panel reads bottom-to-top.
  ord <- stats::aggregate(density ~ dataset, data = pr, FUN = stats::median)
  ord <- ord$dataset[order(ord$density)]
  pal <- series_palette(ord)

  use_log <- isTRUE(input$dens_plot_log)
  ax      <- density_axis_range(
    if (!is.null(lc)) lc$local else numeric(0), pr$density, use_log = use_log)
  xrange  <- ax$range
  clip    <- ax$clip

  fig <- plotly::plot_ly()
  for (g in ord) {
    col <- unname(pal[g])
    if (!is.null(lc) && any(lc$dataset == g)) {
      fig <- fig %>% plotly::add_trace(
        x = clip(lc$local[lc$dataset == g]), y = rep(g, sum(lc$dataset == g)),
        type = "violin", orientation = "h", name = g, legendgroup = g,
        showlegend = FALSE, side = "positive", points = FALSE,
        line = list(color = col, width = 1),
        fillcolor = sprintf("rgba(%s,0.35)", paste(grDevices::col2rgb(col), collapse = ",")),
        hoverinfo = "x", meanline = list(visible = FALSE), box = list(visible = FALSE))
    }
    sub <- pr[pr$dataset == g, , drop = FALSE]
    fig <- fig %>% plotly::add_trace(
      x = sub$density, y = rep(g, nrow(sub)), type = "scatter", mode = "markers",
      name = g, legendgroup = g, showlegend = FALSE,
      marker = list(color = col, size = 11, symbol = "diamond",
                    line = list(color = "#222", width = 1)),
      text = sprintf("%s<br><b>%.1f %s</b> (global n/support)<br>%d units · support/bbox %.2f",
                     sub$dataset, sub$density, un, sub$n_units, sub$frac_of_bbox),
      hoverinfo = "text")
  }

  tgt <- density_target()
  shp <- list(); ann <- list()
  if (!is.null(tgt) && density_thin_on()) {
    shp <- list(list(type = "line", x0 = tgt, x1 = tgt, y0 = 0, y1 = 1, yref = "paper",
                     line = list(color = "#c03030", dash = "dash", width = 2)))
    # An annotation's x is in AXIS coordinates, which on a log axis means log10 of
    # the value - a shape's x0/x1 are not. Getting this wrong parks the label off
    # frame.
    ann <- list(list(x = if (use_log) log10(max(tgt, 1e-9)) else tgt,
                     y = 1, yref = "paper",
                     text = " target", showarrow = FALSE, xanchor = "left",
                     font = list(color = "#c03030", size = 11)))
  }
  fig %>% plotly::layout(
    shapes = shp, annotations = ann, showlegend = FALSE,
    xaxis = list(title = sprintf("Density (%s%s)", un, if (use_log) ", log" else ""),
                 type = if (use_log) "log" else "linear",
                 range = xrange, zeroline = FALSE),
    yaxis = list(title = "", categoryorder = "array", categoryarray = ord,
                 automargin = TRUE),
    margin = list(l = 10, r = 24, t = 28, b = 44),
    title = list(text = "◆ global density   ▏shaded = spread of per-unit local density",
                 font = list(size = 11, color = "#666"), x = 0),
    paper_bgcolor = "white", plot_bgcolor = "white")
}
output$density_plot <- renderPlotly(build_density_plot_fig())

## Support vs gap-fill radius - the diagnostic that makes the closing radius a
## choice rather than a guess. See support_closing_sweep for how to read it.
build_density_sweep_plot_fig <- function() {
  sw <- density_sweep_rv()
  shiny::validate(shiny::need(!is.null(sw) && nrow(sw) > 0,
    "Tick 'Show gap-fill sweep' and press 'Measure density'."))
  pal  <- series_palette(unique(sw$dataset))
  # Mark the radius ACTUALLY used. With auto on that is the resolved plateau, not
  # the (hidden) manual input, which would otherwise draw a misleading marker.
  prof <- density_profile_rv()
  cur  <- if (!is.null(prof) && nrow(prof) && is.finite(prof$close_nm[1])) prof$close_nm[1]
          else max(0, as.numeric(input$dens_close_nm %||% 200))
  fig <- plotly::plot_ly()
  for (nm in unique(sw$dataset)) {
    s <- sw[sw$dataset == nm, ]; s <- s[order(s$close_nm), ]
    fig <- fig %>% plotly::add_trace(
      x = s$close_nm, y = s$frac_of_bbox, type = "scatter", mode = "lines+markers",
      name = nm, line = list(color = unname(pal[nm]), width = 2),
      marker = list(color = unname(pal[nm]), size = 6),
      text = sprintf("%s<br>fill %.0f nm<br>support %.3g nm^%d (%.0f%% of bbox)<br>step %+.1f%%",
                     nm, s$close_nm, s$support, if (isTRUE(input$dens_norm_3d)) 3 else 2,
                     100 * s$frac_of_bbox, 100 * s$rel_step),
      hoverinfo = "text")
  }
  fig %>% plotly::layout(
    shapes = list(list(type = "line", x0 = cur, x1 = cur, y0 = 0, y1 = 1, yref = "paper",
                       line = list(color = "#c03030", dash = "dot", width = 2))),
    annotations = list(list(x = cur, y = 1, yref = "paper",
                            text = sprintf(" in use: %g nm", cur),
                            showarrow = FALSE, xanchor = "left",
                            font = list(color = "#c03030", size = 10))),
    xaxis = list(title = "Gap fill radius (nm)"),
    yaxis = list(title = "Support / frame", range = c(0, 1.04)),
    # The curve rises left-to-right, so the top-left of the panel is always empty -
    # putting the legend there instead of below the axis buys back the vertical
    # space a horizontal legend was taking from the plot itself.
    margin = list(l = 58, r = 24, t = 26, b = 46),
    title = list(text = "Pick a radius where the curve has flattened",
                 font = list(size = 11, color = "#666"), x = 0),
    legend = list(x = 0.02, y = 0.98, xanchor = "left", yanchor = "top",
                  font = list(size = 10), bgcolor = "rgba(255,255,255,0.75)",
                  bordercolor = "#ddd", borderwidth = 1),
    paper_bgcolor = "white", plot_bgcolor = "white")
}
output$density_sweep_plot <- renderPlotly(build_density_sweep_plot_fig())

## ── Support map: what the density denominator actually counted ─────────────
## GOAL: draw one dataset in 2D with the estimated support outlined in GREEN and
## everything inside the imaged frame that was DISREGARDED shaded in RED.
##
## This is the panel that makes the support estimate falsifiable. The numbers say
## "support is 30% of the frame"; only the picture says whether that 30% is the
## cell or is the outline of your own structures - which is the difference between
## a density and a disguised clustering measure.
##
## Recomputed on demand for the SELECTED file only (masks for every file would be
## needlessly heavy), using the same filter chain and the same resolved gap-fill
## radius the profile used, so what is drawn is what was measured.
observe({
  fs <- names(smlm_file_store())
  sel <- isolate(input$dens_map_file)
  updateSelectInput(session, "dens_map_file", choices = fs,
                    selected = if (!is.null(sel) && sel %in% fs) sel else
                               (isolate(input$active_file) %||% fs[1]))
})

build_density_map_fig <- function() {
  prof <- density_profile_rv()
  shiny::validate(shiny::need(!is.null(prof) && nrow(prof) > 0, "Press 'Measure density'."))
  smode <- input$dens_support_mode %||% "grid"
  store <- smlm_file_store()
  fname <- input$dens_map_file %||% (input$active_file %||% names(store)[1])
  shiny::validate(shiny::need(!is.null(fname) && fname %in% names(store), "Pick a dataset."))

  use_3d <- isTRUE(input$dens_norm_3d)
  unit   <- density_unit()
  h_nm   <- max(10, as.numeric(input$dens_grid_h %||% 100))
  # The radius the PROFILE resolved (auto plateau or manual), not the raw input -
  # otherwise the picture could disagree with the measured number.
  close_n <- if ("close_nm" %in% names(prof) && is.finite(prof$close_nm[1]))
    prof$close_nm[1] else max(0, as.numeric(input$dens_close_nm %||% 200))

  df <- qc_exclude_from(store[[fname]], fname)
  df <- tryCatch(apply_meta_filters_batch(df, as_smlm_settings(input)), error = function(e) NULL)
  shiny::validate(shiny::need(!is.null(df) && nrow(df) > 0, "No localizations after filtering."))
  if (unit == "cluster")
    df <- tryCatch(apply_clustering_pipeline(df, as_smlm_settings(input)), error = function(e) df)

  keep <- if (use_3d) c("x", "y", "z") else c("x", "y")
  shiny::validate(shiny::need(all(keep %in% names(df)), "Missing coordinates."))
  sp <- df[, keep, drop = FALSE]
  if (!use_3d) sp$z <- 0
  sp <- sp[stats::complete.cases(sp), , drop = FALSE]
  shiny::validate(shiny::need(nrow(sp) >= 3, "Too few localizations."))

  sup <- estimate_support(sp, use_3d = use_3d, mode = smode, h_nm = h_nm,
                          close_nm = close_n, return_mask = TRUE,
                          bw_nm = max(10, as.numeric(input$dens_kde_bw %||% 150)),
                          frac = min(1, max(0.5, as.numeric(input$dens_kde_frac %||% 0.99))))
  shiny::validate(shiny::need(!is.null(sup), "Could not estimate the support."))

  win  <- if (!is.null(sup$win)) sup$win else apply(as.matrix(sp[, c("x","y")]), 2, range)
  npts <- nrow(sp)
  sp_d <- if (npts > 30000L) sp[seq(1, npts, length.out = 30000L), , drop = FALSE] else sp

  # Red = disregarded, but only INSIDE the data bounding box. Two reasons: the grid's
  # padding ring exists solely to give the dilation room, and a localization table
  # carries no record of the real scan field - so the widest defensible "frame" is
  # the extent of the data itself. Consequence, stated in the title: red marks
  # interior dark regions and concave margins, never area outside the data.
  fig <- plotly::plot_ly()
  shp <- list()

  # Grid and KDE both produce a cell mask, so they share one renderer; only hull
  # and bbox need the polygon path. `sup$h` is the cell size the estimator
  # actually used, which is NOT necessarily the h_nm input.
  if (smode %in% c("grid", "kde")) {
    shiny::validate(shiny::need(!is.null(sup$mask), "Could not build the support mask."))
    m2 <- mask_project_2d(sup$mask)
    d  <- dim(m2)
    cell <- sup$h %||% h_nm
    x0 <- sup$origin[1]; y0 <- sup$origin[2]
    # Cell CENTRES for the heatmap (plotly places a heatmap cell on its centre);
    # the outline below uses cell EDGES. Mixing the two shifts the overlay by h/2.
    xc <- x0 + (seq_len(d[1]) - 0.5) * cell
    yc <- y0 + (seq_len(d[2]) - 0.5) * cell
    in_frame <- outer(xc >= win[1, 1] & xc <= win[2, 1],
                      yc >= win[1, 2] & yc <= win[2, 2], "&")
    zred <- ifelse(!m2 & in_frame, 1, NA_real_)
    fig <- fig %>% plotly::add_heatmap(
      x = xc, y = yc, z = t(zred), showscale = FALSE, hoverinfo = "skip",
      colorscale = list(list(0, "rgba(200,48,48,0.22)"), list(1, "rgba(200,48,48,0.22)")))
    outl <- mask_outline_segments(m2, c(x0, y0), cell)
  } else {
    # Hull and bbox have no interior structure, so they are drawn exactly as
    # POLYGONS rather than discretised: a red rectangle over the whole data box,
    # then the support polygon filled opaque white on top of it. Points and the
    # green outline are traces, which plotly draws above all shapes.
    rect_path <- sprintf("M %f,%f L %f,%f L %f,%f L %f,%f Z",
                         win[1,1], win[1,2], win[2,1], win[1,2],
                         win[2,1], win[2,2], win[1,1], win[2,2])
    if (smode == "hull") {
      hidx <- tryCatch(grDevices::chull(sp$x, sp$y), error = function(e) integer(0))
      shiny::validate(shiny::need(length(hidx) >= 3, "Could not build a convex hull."))
      hx <- sp$x[hidx]; hy <- sp$y[hidx]
      hull_path <- paste0("M ", paste(sprintf("%f,%f", hx, hy), collapse = " L "), " Z")
      shp <- list(
        list(type = "path", path = rect_path, fillcolor = "rgba(200,48,48,0.22)",
             line = list(width = 0), layer = "below"),
        list(type = "path", path = hull_path, fillcolor = "#ffffff",
             line = list(width = 0), layer = "below"))
      outl <- list(x = c(hx, hx[1]), y = c(hy, hy[1]))
    } else {
      # bbox: the support IS the frame, so nothing is disregarded.
      shp  <- list()
      outl <- list(x = c(win[1,1], win[2,1], win[2,1], win[1,1], win[1,1]),
                   y = c(win[1,2], win[1,2], win[2,2], win[2,2], win[1,2]))
    }
  }

  fig <- fig %>% plotly::add_trace(
    x = sp_d$x, y = sp_d$y, type = "scattergl", mode = "markers",
    marker = list(size = 2, color = "rgba(40,40,40,0.55)"),
    name = "localizations", hoverinfo = "skip", showlegend = FALSE)
  if (length(outl$x))
    fig <- fig %>% plotly::add_trace(
      x = outl$x, y = outl$y, type = "scatter", mode = "lines",
      line = list(color = "#1a8a1a", width = 1.6), connectgaps = FALSE,
      name = "counted as sample", hoverinfo = "skip", showlegend = FALSE)
  # Frame outline, so "how much of the field was disregarded" is readable.
  fig <- fig %>% plotly::add_trace(
    x = c(win[1,1], win[2,1], win[2,1], win[1,1], win[1,1]),
    y = c(win[1,2], win[1,2], win[2,2], win[2,2], win[1,2]),
    type = "scatter", mode = "lines",
    line = list(color = "#999", width = 1, dash = "dot"),
    name = "frame", hoverinfo = "skip", showlegend = FALSE)

  mode_lbl <- switch(smode,
    kde  = sprintf("smoothed density, %.0f%% of units enclosed, %g nm bandwidth",
                   100 * (sup$frac %||% 0.99), sup$bw %||% NA_real_),
    grid = sprintf("%g nm grid, %g nm gap fill (dark regions wider than ~%g nm stay excluded)",
                   h_nm, close_n, 2 * close_n),
    hull = "convex hull - bridges every concavity and dark region by construction",
    "bounding box - nothing is excluded")
  ttl <- sprintf(
    "%s - <span style='color:#1a8a1a'>green = counted as sample (%.0f%% of the data's bounding box)</span>, <span style='color:#c03030'>red = disregarded</span> · %s%s",
    fname, 100 * sup$frac_of_bbox, mode_lbl,
    if (use_3d && smode == "grid") " · plan view of a 3-D support (a column counts if any z slice does)" else "")
  fig %>% plotly::layout(
    shapes = shp,
    title = list(text = ttl, font = list(size = 11, color = "#555"), x = 0),
    xaxis = list(title = "x (nm)", zeroline = FALSE, constrain = "domain"),
    # Equal aspect: 1 nm in x is 1 nm in y, or the support's shape is a lie.
    yaxis = list(title = "y (nm)", zeroline = FALSE,
                 scaleanchor = "x", scaleratio = 1),
    margin = list(l = 60, r = 20, t = 34, b = 44),
    paper_bgcolor = "white", plot_bgcolor = "white") %>%
    plotly::config(doubleClick = "reset")
}
output$density_map <- renderPlotly(build_density_map_fig())

output$download_density <- downloadHandler(
  filename = function() paste0("labelling_density_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
  content = function(file) {
    prof <- density_profile_rv()
    shiny::validate(shiny::need(!is.null(prof) && nrow(prof) > 0, "No results."))
    pm <- density_p_map()
    if (!is.null(pm)) prof$keep_p <- unname(pm[prof$dataset])
    prof$target_density <- density_target() %||% NA_real_
    write.csv(prof, file, row.names = FALSE)
  }
)

# Keep the control-group selector populated from the group names.
observe({
  gm <- file_group_map()
  gs <- unique(unname(gm))
  updateSelectInput(session, "dens_target_group", choices = gs,
                    selected = isolate(input$dens_target_group) %||% gs[1])
})
