## ─────────────────────────────────────────────────────────────────────────────
## server/tracking/16_tsa_cycles.R  -  panel: time-signature analysis (cycle classes)
## ─────────────────────────────────────────────────────────────────────────────
## Part of tracking_server(). Not a standalone script: its statements are evaluated in
## the server function's environment, in file-name order (server/parts.R).
##
## Integrated versus dark steps, before/after cleaning, cycle progression, and
## the MSD downloads.
##
## INTERFACE (generated: Rscript tools/server_map.R --headers; do not edit by hand)
##   DEFINES for other parts:
##     cycle_class_bundle
##   READS from other parts:
##     01_state_sigma: tsa_ran
##     04_datasets_upload: trk_raw_original
##     06_filter_pipeline: trk_spatial_filtered
##     07_cycle_time_sequence: trk_cluster_cycle, trk_cycle_time,
##                             trk_hw_cycle_time
##     11_data_frames: trk
##     13_dt_msd_clusters: dt_vectors, msd_export
## END INTERFACE
## ─────────────────────────────────────────────────────────────────────────────

# ==========================================================================
# cycle_class_bundle() - INTEGRATED vs DARK CYCLE CLASSIFICATION (TSA core)
# ==========================================================================
#
# GOAL
#   Reconstruct WHAT THE MICROSCOPE WAS DOING between consecutive
#   localizations. MINFLUX runs on a fixed acquisition cycle, but it does not
#   report a position every cycle - if the molecule was dark (blinking) or too
#   few photons arrived, that cycle yields nothing and the next localization
#   appears one or more cycles later. This turns the raw timestamps into that
#   cycle-level story, which is the foundation every Δt-aware estimator uses.
#
#   Adapted from Bela Vogler's MINFLUX Time-Signature Analysis (MFX-TSA),
#   Vogler et al. 2025 (Commun. Biol.); reference implementation
#   Toolbox_for_MINFLUX_SPT (cf. its `print_cycle_stats`).
#
# THE MATH (deliberately simple, and that is the point)
#   For each consecutive interval within a track:
#
#        n = round(Δt / t_cycle)          cycles spanned by the step
#
#        n = 1  -> INTEGRATED: localized on the very next cycle, no gap.
#        n >= 2 -> (n − 1) DARK cycles passed before the molecule reappeared.
#
#   Everything downstream is built on this one integer:
#     - integrated fraction = the headline data-quality score
#     - the gap-split threshold is expressed in dark cycles
#     - the MSD's cumulative cycle index is cumsum(n)
#     - CVE, JDD and directionality all filter on n == 1
#
# DIFFERENCE FROM VOGLER'S IMPLEMENTATION (stated honestly)
#   His TSA SIMULATES the hardware timing (dwell time, pattern vertices, EOD
#   and software overheads) and fits those parameters to the data. We infer
#   t_cycle from the Δt distribution instead - no hardware simulation. That is
#   why the app's cycle time is a data-driven estimate (cluster cycle, or an
#   optional regression) rather than a fitted overhead model, and why the
#   sequence-JSON import exists as a cross-check.
#
#   IMPORTANT: the MODAL Δt is a poor estimate of t_cycle when dark cycles are
#   common, because the mode can land on a 2-cycle gap. The cluster-cycle
#   estimator (smallest robust Δt cluster) is used instead.
#
# OBJECTS IT CALLS (inputs)
#   dt_vectors() (per-track Δt), trk_cycle_time(), trk_cluster_cycle()
#   inputs: trk_tsa_run (Run trigger), trk_cycle_dark_max
#
# RETURNS the per-interval cycle counts plus the summary statistics the TSA
# panel displays: integrated fraction, dark-cycle distribution, per-track
# cycle progression, and the before/after comparison for gap splitting.
#
# CAVEAT: descriptive, not inferential. If t_cycle is wrong, EVERY cycle count
# is wrong - check the data-quality verdict and the sequence JSON before
# trusting anything downstream.
# ==========================================================================
cycle_class_bundle <- eventReactive(input$trk_tsa_run, {
  req(isTRUE(input$trk_show_tsa))
  withProgress(message = "Time Signature Analysis…", value = 0.5, {
  t_c <- trk_cycle_time()   # single source of truth
  dts <- dt_vectors()
  vec <- unlist(dts, use.names = FALSE)
  vec <- vec[is.finite(vec) & vec > 0]
  shiny::validate(shiny::need(length(vec) > 0, "No Δt data for cycle classification."))
  shiny::validate(shiny::need(is.finite(t_c) && t_c > 0, "Could not determine cycle time."))
  hw_used <- { hw <- trk_hw_cycle_time(); is.finite(hw) && hw > 0 }

  n_cyc      <- pmax(1L, as.integer(round(vec / t_c)))
  n_total    <- length(n_cyc)
  n_integ    <- sum(n_cyc == 1L)
  n_dark_int <- sum(n_cyc >= 2L)           # intervals containing ≥1 dark cycle
  dark_cycles_total <- sum(n_cyc - 1L)     # total dark cycles across all intervals
  frac_integ <- n_integ / n_total

  max_disp <- max(5L, as.integer(input$trk_cycle_dark_max %||% 20L))
  tab      <- tabulate(pmin(n_cyc, max_disp), nbins = max_disp)
  list(t_cycle = t_c, n_cyc = n_cyc, counts = tab, max_disp = max_disp,
       n_total = n_total, n_integ = n_integ, n_dark_int = n_dark_int,
       dark_cycles_total = dark_cycles_total, frac_integ = frac_integ,
       hw_used = hw_used)
  })  # withProgress
})

output$trk_cycle_class_plot <- renderPlotly({
  cc <- cycle_class_bundle()
  dark <- isTRUE(isolate(input$trk_dark_mode))
  bg_col <- if (dark) "#000000" else "#FFFFFF"
  fg_col <- if (dark) "#FFFFFF" else "#333333"

  x_lab  <- seq_len(cc$max_disp)
  frac   <- cc$counts / cc$n_total
  # Integrated (n=1) green, Dark (n>=2) amber→red gradient
  cols <- ifelse(x_lab == 1L, "#2CA02C",
            grDevices::colorRampPalette(c("#E0A030", "#C03030"))(cc$max_disp)[x_lab])
  x_txt <- as.character(x_lab); x_txt[cc$max_disp] <- paste0("≥", cc$max_disp)

  plot_ly(x = x_lab, y = frac, type = "bar",
          marker = list(color = cols),
          text = sprintf("%d cycle%s: %.1f%%", x_lab, ifelse(x_lab==1,"","s"), frac*100),
          hovertemplate = "%{text}<extra></extra>") %>%
    layout(
      xaxis = list(title = "Cycles per interval (1 = Integrated, ≥2 = Dark)",
                   tickmode = "array", tickvals = x_lab, ticktext = x_txt,
                   color = fg_col),
      yaxis = list(title = "Fraction of intervals", color = fg_col),
      paper_bgcolor = bg_col, plot_bgcolor = bg_col,
      font = list(color = fg_col), showlegend = FALSE,
      bargap = 0.1, uirevision = "cycle_class"
    )
})

output$trk_cycle_class_summary <- renderUI({
  cc <- cycle_class_bundle()
  htmltools::div(
    style = "font-size:0.87em; color:#444; margin:6px 0;",
    htmltools::tags$b("Cycle classification "),
    sprintf("(cycle time = %.4g ms, %s):", cc$t_cycle,
            if (isTRUE(cc$hw_used)) "snapped to hardware model" else "from modal Δt"),
    htmltools::br(),
    sprintf("Integrated (1 cycle): %.1f%% of %s intervals  |  with Dark cycles (≥2): %.1f%%",
            cc$frac_integ * 100, format(cc$n_total, big.mark = ","),
            (1 - cc$frac_integ) * 100),
    htmltools::br(),
    sprintf("Total Dark cycles: %s  |  mean cycles per interval: %.2f",
            format(cc$dark_cycles_total, big.mark = ","),
            mean(cc$n_cyc)),
    htmltools::br(),
    htmltools::tags$span(
      style = "color:#888; font-size:0.9em;",
      if (cc$frac_integ > 0.9)
        "Excellent tracking quality - almost every cycle yielded a localization."
      else if (cc$frac_integ > 0.7)
        "Good tracking - most cycles localized, modest dark fraction."
      else
        "Frequent dropouts. Splitting tracks at gaps (Display tab) keeps long dark gaps out of the analyses."
    )
  )
})


# Raw vs analyzed cycle-stats: the NET effect of all track cleaning (gap-split +
# end-trim + quality gate + metadata + length/duration) on the Integrated/Dark
# balance. "Raw" = the uploaded data; "Analyzed" = what the analyses actually use.
# They differ whenever any cleaning is active; identical means nothing is cleaning.
output$trk_cycle_class_beforeafter <- renderUI({
  req(isTRUE(input$trk_show_tsa), tsa_ran())
  t_c    <- trk_cycle_time()
  before <- compute_cycle_stats(tryCatch(trk_raw_original(),     error = function(e) NULL), t_c)
  after  <- compute_cycle_stats(tryCatch(trk_spatial_filtered(), error = function(e) NULL), t_c)
  if (is.null(before)) return(NULL)
  any_clean <- isTRUE(input$trk_split_gaps_enable) || isTRUE(input$trk_endtrim_enable) ||
               isTRUE(input$trk_quality_enable)

  cell <- function(x) htmltools::tags$td(x, style = "padding:2px 10px; text-align:right;")
  row  <- function(lbl, b, a) htmltools::tags$tr(
    htmltools::tags$td(lbl, style = "padding:2px 10px;"), cell(b), cell(a))
  fmt  <- function(v, pct = FALSE) if (is.null(v)) "-" else
    if (pct) sprintf("%.1f%%", v * 100) else format(round(v, 2), big.mark = ",")

  htmltools::div(style = "font-size:0.86em; color:#444; margin:8px 0;",
    htmltools::tags$b("Cycle stats - raw vs analyzed (after cleaning)"),
    if (!any_clean) htmltools::tags$span(style = "color:#888;",
      " (no cleaning active - enable Split / Trim / Quality gate in the Display tab to see the effect)"),
    htmltools::tags$table(style = "border-collapse:collapse; margin-top:4px;",
      htmltools::tags$thead(htmltools::tags$tr(
        htmltools::tags$th("", style = "padding:2px 10px; border-bottom:1px solid #ccc;"),
        htmltools::tags$th("Raw", style = "padding:2px 10px; border-bottom:1px solid #ccc;"),
        htmltools::tags$th("Analyzed", style = "padding:2px 10px; border-bottom:1px solid #ccc;"))),
      htmltools::tags$tbody(
        row("Tracks / segments", fmt(before$n_tracks), fmt(after$n_tracks)),
        row("Integrated (1 cycle)", fmt(before$frac_integ, TRUE), fmt(after$frac_integ, TRUE)),
        row("With Dark cycles (≥2)", fmt(1 - before$frac_integ, TRUE), fmt(1 - after$frac_integ, TRUE)),
        row("Total Dark cycles", fmt(before$dark_total), fmt(after$dark_total)),
        row("Mean cycles / interval", fmt(before$mean_cyc), fmt(after$mean_cyc)),
        row("Max cycles in one gap", fmt(before$max_cyc), fmt(after$max_cyc)))),
    htmltools::tags$p(style = "color:#888; font-size:0.95em; margin-top:4px;",
      "Cleaning (split at dark gaps, trim degraded ends, drop dark-heavy tracks) should raise ",
      "the Integrated fraction, cut total Dark cycles, and shrink the max gap - the analyzed ",
      "intervals are clean, near-consecutive localizations."))
})

# Per-track cycle progression (plot_cycle_progression_for_ID rewrite to R)
# Track IDs are served on demand (server-side selectize): a dataset can hold
# thousands of tracks, and sending them all to the browser makes the selector slow.
observe({
  req(isTRUE(input$trk_show_tsa))
  df <- tryCatch(trk(), error = function(e) NULL)
  req(!is.null(df), nrow(df) > 0)
  counts <- table(df$tid)
  tids   <- names(counts[counts >= 5L])
  req(length(tids) > 0)
  tids <- tids[order(as.numeric(counts[tids]), decreasing = TRUE)]   # longest tracks first
  cur  <- isolate(input$trk_cycle_prog_tid)
  updateSelectizeInput(session, "trk_cycle_prog_tid", choices = tids,
                       selected = if (isTRUE(cur %in% tids)) cur else tids[1], server = TRUE)
})

output$trk_cycle_prog_plot <- renderPlotly({
  req(isTRUE(input$trk_show_tsa), tsa_ran())
  req(nzchar(input$trk_cycle_prog_tid))
  df  <- trk(); req(nrow(df) > 0)
  t_c <- trk_cycle_time()
  if (!is.finite(t_c) || t_c <= 0) t_c <- 1
  tb <- df[df$tid == input$trk_cycle_prog_tid, , drop = FALSE]
  tb <- tb[order(tb$tim), ]
  shiny::validate(shiny::need(nrow(tb) >= 2, "Track too short."))

  # Cluster (atomic) cycle + its SW offset for the cycle-count mapping
  cc <- tryCatch(trk_cluster_cycle(), error = function(e) NULL)
  sw <- if (!is.null(cc) && is.finite(cc$sw)) cc$sw else 0

  dt_v  <- diff(tb$tim)
  n_cyc <- pmax(1L, as.integer(round((dt_v - sw) / t_c)))  # total cycles to next loc
  dark  <- n_cyc - 1L                                       # dark cycles (0 = Integrated)
  step  <- seq_along(n_cyc)
  is_int <- dark == 0L

  dm <- isTRUE(isolate(input$trk_dark_mode))
  bg_col <- if (dm) "#000000" else "#FFFFFF"
  fg_col <- if (dm) "#FFFFFF" else "#333333"
  hov <- function(i) sprintf("step %d: %d cycle%s = %s (%.3g ms)",
            step[i], n_cyc[i], ifelse(n_cyc[i]==1,"","s"),
            ifelse(is_int[i], "Integrated", sprintf("%d Dark", dark[i])), dt_v[i])

  plot_ly() %>%
    add_trace(x = step, y = dark, type = "scatter", mode = "lines",
              line = list(color = "rgba(120,120,120,0.4)", width = 1),
              showlegend = FALSE, hoverinfo = "skip") %>%
    add_trace(x = step[is_int], y = dark[is_int], type = "scatter", mode = "markers",
              marker = list(color = "#2070D0", size = 6),
              name = "Integrated (0 dark)",
              text = vapply(which(is_int), hov, character(1)),
              hovertemplate = "%{text}<extra></extra>") %>%
    add_trace(x = step[!is_int], y = dark[!is_int], type = "scatter", mode = "markers",
              marker = list(color = "#C03030", size = 6),
              name = "Dark cycles present",
              text = vapply(which(!is_int), hov, character(1)),
              hovertemplate = "%{text}<extra></extra>") %>%
    layout(
      xaxis = list(title = "Localization index along track", color = fg_col),
      yaxis = list(title = "Dark cycles before this localization", color = fg_col,
                   rangemode = "tozero"),
      paper_bgcolor = bg_col, plot_bgcolor = bg_col,
      font = list(color = fg_col),
      legend = list(orientation = "h", x = 0, y = 1.12, font = list(color = fg_col)),
      uirevision = "cycle_prog"
    )
})

# Plain CSV: tidy columns for Excel / Prism
output$trk_dl_msd_csv <- downloadHandler(
  filename = function() {
    ex <- msd_export(); mode <- ex$params$mode
    paste0("msd_", mode %||% "time", ".csv")
  },
  content = function(file) {
    ex <- msd_export()
    utils::write.csv(ex$data, file, row.names = FALSE)
  }
)

# ZIP: includes msd_data.csv + params.txt (easy to read in Excel/Notepad/Prism)
output$trk_dl_msd_zip <- downloadHandler(
  filename = function() {
    ex <- msd_export(); mode <- ex$params$mode
    paste0("msd_", mode %||% "time", "_export.zip")
  },
  content = function(file) {
    ex <- msd_export()
    td <- tempdir()
    f_csv <- file.path(td, "msd_data.csv")
    f_txt <- file.path(td, "params.txt")
    
    # data
    utils::write.csv(ex$data, f_csv, row.names = FALSE)
    
    # human-readable params
    lines <- c(
      "MSD export parameters:",
      sprintf("mode: %s", ex$params$mode),
      sprintf("time_unit: %s", ex$params$time_unit),
      sprintf("space_unit: %s", ex$params$space_unit),
      sprintf("tau_value: %s", ex$params$tau_value),
      sprintf("msd_min_pts: %s", ex$params$msd_min_pts),
      sprintf("msd_weighting: %s", ex$params$msd_weighting),
      sprintf("msd_fit_method: %s", ex$params$msd_fit_method),
      sprintf("msd_tau_max: %s", ex$params$msd_tau_max),
      sprintf("msd_tau_bin_used_ms: %s", ex$params$msd_tau_bin_used),
      sprintf("msd_tau_bin_auto: %s", ex$params$msd_tau_bin_auto),
      sprintf("msd_fit_min_lag_cycles: %s", ex$params$msd_fit_min_lag_cycles),
      sprintf("msd_fit_min_lag_ms: %s", ex$params$msd_fit_min_lag_ms),
      sprintf("msd_ndim: %s", ex$params$msd_ndim),
      sprintf("msd_formula: %s", ex$params$msd_formula),
      sprintf("dt0_ms (median dt): %s", ex$params$dt0_ms)
    )
    writeLines(lines, f_txt)
    
    # zip it
    zip::zipr(zipfile = file, files = c(f_csv, f_txt), root = td)
  }
)
