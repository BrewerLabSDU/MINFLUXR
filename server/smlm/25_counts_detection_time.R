## ─────────────────────────────────────────────────────────────────────────────
## server/smlm/25_counts_detection_time.R  -  panels: traces per cluster, detection over time
## ─────────────────────────────────────────────────────────────────────────────
## Part of smlm_server(). Not a standalone script: its statements are evaluated in
## the server function's environment, in file-name order (server/parts.R).
##
## Traces per cluster, Median Cluster detection over time and trace detection
## over time.
##
## INTERFACE (generated: Rscript tools/server_map.R --headers; do not edit by hand)
##   DEFINES for other parts:
##     build_det_time_plot_fig, build_tid_count_hist_fig,
##     build_tid_time_plot_fig, tid_count_results
##   READS from other parts:
##     01_view_upload: smlm_file_store
##     20_dataset_groups: apply_dataset_groups, series_palette
##     21_density_normalisation: apply_density_thinning
##     23_prepared_cache: batch_filters, cluster2_label, cluster2_on,
##                        deblink_label, prepared_df
##     27_fiducial_qc: qc_exclude_from
## END INTERFACE
## ─────────────────────────────────────────────────────────────────────────────

# ------ Traces per Cluster -------------------------------------------------
tid_count_results <- reactiveVal(NULL)
tid_count_params  <- reactiveVal(NULL)   # compute params, frozen at Run

observeEvent(input$run_tid_count, {
  store <- smlm_file_store()
  if (!length(store)) {
    showNotification("No files loaded.", type = "warning"); return()
  }

  # The GROUPING is ALWAYS the median cluster (cluster2). The UNIT counted
  # inside each median cluster is either raw trace IDs ("tid") or primary
  # de-blinked molecules ("pcid" - needs de-blinking enabled).
  unit      <- input$tid_count_unit %||% "tid"
  count_col <- if (unit == "pcid") "cluster" else "tid"

  all_counts <- list()
  withProgress(message = "Counting units per median cluster...", value = 0, {
    n_files <- length(store)
    for (i in seq_along(store)) {
      fname <- names(store)[[i]]
      incProgress(1 / n_files, detail = fname)
      df <- prepared_df(fname)   # cached filter -> thin -> cluster chain
      if (is.null(df) || !nrow(df)) next
      # Need median clusters (cluster2) to group by.
      if (!"cluster2" %in% names(df) || all(is.na(df$cluster2) | df$cluster2 == 0L)) next
      # Counting de-blinked molecules needs the de-blink (primary `cluster`) column present.
      if (count_col == "cluster" &&
          (!"cluster" %in% names(df) || all(is.na(df$cluster) | df$cluster == 0L))) next

      # Split after preparation, so a 2-colour split counts as channels too.
      # Every Median Cluster already belongs to exactly one channel (the
      # clustering ran inside channels), so this partitions the clusters
      # without changing any count.
      parts <- if (isTRUE(input$tid_count_by_channel)) channel_parts(df, fname)
               else list(list(df = df, tag = fname))

      for (pt in parts) {
        sub <- pt$df %>% dplyr::filter(!is.na(cluster2) & cluster2 > 0L)
        # Legacy guard: de-blinking assigns every trace, so cluster 0 should not occur.
        if (count_col == "cluster") sub <- sub %>% dplyr::filter(!is.na(cluster) & cluster > 0L)
        if (!nrow(sub)) next
        cnt <- sub %>%
          dplyr::group_by(cluster2) %>%
          dplyr::summarise(n_tids = dplyr::n_distinct(.data[[count_col]]), .groups = "drop") %>%
          dplyr::mutate(dataset = pt$tag)
        all_counts[[pt$tag]] <- cnt
      }
    }
  })
  counts <- do.call(rbind, Filter(Negate(is.null), all_counts))

  if (is.null(counts) || !nrow(counts)) {
    hint <- if (count_col == "cluster")
      sprintf("Enable BOTH %s and %s in the Display tab, then press Plot for each file.",
              deblink_label(), cluster2_label())
    else
      sprintf("Enable %s in the Display tab, then press Plot for each file.", cluster2_label())
    showNotification(paste0("No median clusters found. ", hint), type = "warning")
    tid_count_results(NULL)
    return()
  }
  tid_count_results(counts)
  # Compute parameters frozen at Run: reading them live relabelled an existing
  # plot without recomputing it, so the axis could describe a different unit
  # than the data on screen.
  tid_count_params(list(unit = unit))
  unit_lbl <- if (count_col == "cluster") "de-blinked molecules" else "trace IDs"
  showNotification(
    sprintf("Counted %s in %d median clusters across %d file(s) (median = %g).",
            unit_lbl, nrow(counts), length(unique(counts$dataset)),
            stats::median(counts$n_tids)),
    type = "message", duration = 4)
})

build_tid_count_hist_fig <- function() {
  counts <- tid_count_results()
  shiny::validate(shiny::need(!is.null(counts) && nrow(counts) > 0,
                              "Press 'Run' to count traces per cluster."))
  counts <- apply_dataset_groups(counts)
  shiny::validate(shiny::need(nrow(counts) > 0, "No files assigned to groups (group mode)."))

  bins      <- max(5L, as.integer(input$tid_count_bins %||% 30L))
  density   <- identical(input$tid_count_yaxis %||% "count", "density")
  unit      <- (tid_count_params()$unit) %||% (input$tid_count_unit %||% "tid")
  # Grouping is always the median cluster; the unit is what we counted inside it.
  unit_lbl  <- if (unit == "pcid") "De-blinked molecules" else "Trace IDs"
  aggregate <- isTRUE(input$tid_count_aggregate)

  med_n   <- stats::median(counts$n_tids)
  grouped <- (input$analysis_by %||% "file") == "group" && !aggregate
  # y-extent for per-group median line-traces (approx histogram top)
  rng_c <- range(counts$n_tids, na.rm = TRUE)
  if (!all(is.finite(rng_c)) || diff(rng_c) <= 0) rng_c <- c(rng_c[1], rng_c[1] + 1)
  brk_c <- seq(rng_c[1], rng_c[2], length.out = bins + 1L)
  ytop_c <- suppressWarnings(max(vapply(
    if (grouped) split(counts$n_tids, counts$dataset) else list(counts$n_tids),
    function(x) { h <- graphics::hist(x, breaks = brk_c, plot = FALSE, include.lowest = TRUE)
                  max(if (density) h$density else h$counts) }, numeric(1)), na.rm = TRUE))
  if (!is.finite(ytop_c) || ytop_c <= 0) ytop_c <- 1
  fig   <- plotly::plot_ly()

  if (aggregate) {
    fig <- fig %>% plotly::add_histogram(
      x        = counts$n_tids,
      name     = sprintf("All datasets (N clusters = %d)", nrow(counts)),
      nbinsx   = bins,
      histnorm = if (density) "probability density" else "",
      marker   = list(color = "#4C8FBF99", line = list(color = "white", width = 0.5)),
      opacity  = 0.75
    )
  } else {
    datasets <- if ("dataset" %in% names(counts)) unique(counts$dataset) else NULL
    n_ds     <- length(datasets)
    pal <- if (n_ds >= 1) series_palette(datasets) else NULL

    if (!is.null(pal)) {
      for (ds in datasets) {
        sub <- counts[counts$dataset == ds, , drop = FALSE]
        fig <- fig %>% plotly::add_histogram(
          x        = sub$n_tids,
          name     = ds,
          legendgroup = ds,
          nbinsx   = bins,
          histnorm = if (density) "probability density" else "",
          marker   = list(color = paste0(pal[ds], "99"),
                          line  = list(color = pal[ds], width = 0.5)),
          opacity  = 0.75
        )
      }
      # per-group median line-traces (toggle away with their group)
      if (grouped) {
        for (ds in datasets) {
          m <- stats::median(counts$n_tids[counts$dataset == ds], na.rm = TRUE)
          fig <- fig %>% plotly::add_trace(x = c(m, m), y = c(0, ytop_c),
            type = "scatter", mode = "lines", legendgroup = ds, showlegend = FALSE,
            hoverinfo = "skip", line = list(color = pal[ds], width = 2, dash = "dash"))
        }
      }
    } else {
      fig <- fig %>% plotly::add_histogram(
        x        = counts$n_tids,
        nbinsx   = bins,
        histnorm = if (density) "probability density" else "",
        marker   = list(color = "#4C8FBF", line = list(color = "white", width = 0.5)),
        name     = "Clusters"
      )
    }
  }

  # overall median line/label only when not grouped (grouped shows per-group lines)
  ov_shapes <- if (!grouped) list(list(
      type = "line", x0 = med_n, x1 = med_n, y0 = 0, y1 = 1, yref = "paper",
      line = list(color = "#E53935", dash = "dash", width = 1.5))) else list()
  ov_annot  <- if (!grouped) list(list(
      x = med_n, y = 1, yref = "paper", xanchor = "left",
      text = paste0(" median = ", round(med_n, 1)), showarrow = FALSE,
      font = list(color = "#E53935", size = 12))) else list()
  fig %>%
    plotly::layout(
      barmode = "overlay",
      shapes = ov_shapes,
      annotations = ov_annot,
      xaxis = list(title = paste0(unit_lbl, " per Median Cluster")),
      yaxis = list(title = if (density) "Density" else "Number of median clusters"),
      legend = list(orientation = "h", y = -0.18, groupclick = "togglegroup"),
      paper_bgcolor = "white",
      plot_bgcolor  = "white"
    )
}
output$tid_count_hist <- renderPlotly(build_tid_count_hist_fig())

output$download_tid_count <- downloadHandler(
  filename = function()
    paste0("tid_per_cluster_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
  content = function(file) {
    counts <- tid_count_results()
    shiny::validate(shiny::need(!is.null(counts) && nrow(counts) > 0, "No results."))
    write.csv(counts, file, row.names = FALSE)
  }
)

## ── Cluster2 Detection over Time ──────────────────────────────────────────
# Each TID is one detection event placed at its median acquisition time
# (tim_med), tagged in_c2 = whether it belongs to a Cluster2 structure. The
# plot contrasts structure vs background detection over the acquisition.
det_time_results <- reactiveVal(NULL)

observeEvent(input$run_det_time, {
  store <- smlm_file_store()
  if (!length(store)) { showNotification("No files loaded.", type = "warning"); return() }
  if (!cluster2_on()) {
    showNotification(sprintf("Enable %s in the Display tab first.", cluster2_label()),
                     type = "warning"); return()
  }
  if (!isTRUE(input$det_time_all_files)) {
    af <- input$active_file
    if (!is.null(af) && af %in% names(store)) store <- store[af]
  }

  all_tid <- list()
  withProgress(message = "Computing detection over time...", value = 0, {
    n_files <- length(store)
    for (i in seq_along(store)) {
      fname <- names(store)[[i]]
      incProgress(1 / n_files, detail = fname)
      df <- prepared_df(fname)   # cached filter -> thin -> cluster chain
      if (is.null(df) || !nrow(df) || !"cluster2" %in% names(df)) next

      # Per channel when asked. This is the panel the multiplex boundary
      # question is read from: with a TIM-split build, the detection rate per
      # channel is what shows whether a channel is losing traces over its own
      # window rather than over the acquisition as a whole.
      parts <- if (isTRUE(input$det_time_by_channel)) channel_parts(df, fname)
               else list(list(df = df, tag = fname))

      for (pt in parts) {
        per_tid <- pt$df %>%
          dplyr::group_by(tid) %>%
          dplyr::summarise(
            tim_med = stats::median(tim, na.rm = TRUE),
            in_c2   = any(!is.na(cluster2) & cluster2 > 0L),
            .groups = "drop"
          ) %>%
          dplyr::filter(is.finite(tim_med)) %>%
          dplyr::mutate(dataset = pt$tag)
        if (nrow(per_tid)) all_tid[[pt$tag]] <- per_tid
      }
    }
  })
  res <- do.call(rbind, Filter(Negate(is.null), all_tid))
  if (is.null(res) || !nrow(res)) {
    showNotification(sprintf("No traces with a tim column and Median cluster found. Enable %s and press Plot for each file.", cluster2_label()),
                     type = "warning")
    det_time_results(NULL); return()
  }
  det_time_results(res)
  showNotification(
    sprintf("Detection over time: %d traces across %d file(s); %.0f%% in Median Cluster",
            nrow(res), length(unique(res$dataset)), 100 * mean(res$in_c2)),
    type = "message", duration = 4)
})

build_det_time_plot_fig <- function() {
  res <- det_time_results()
  shiny::validate(shiny::need(!is.null(res) && nrow(res) > 0,
                              "Press 'Run' to compute detection over time."))
  res <- apply_dataset_groups(res)
  shiny::validate(shiny::need(nrow(res) > 0, "No files assigned to groups (group mode)."))
  grouped <- (input$analysis_by %||% "file") == "group"

  bins <- max(5L, as.integer(input$det_time_bins %||% 40L))
  mode <- input$det_time_yaxis %||% "count"
  normalized <- mode %in% c("norm", "norm_cumulative")

  rng <- range(res$tim_med, na.rm = TRUE)
  if (!all(is.finite(rng)) || diff(rng) <= 0) rng <- c(rng[1], rng[1] + 1)
  breaks  <- seq(rng[1], rng[2], length.out = bins + 1L)
  ctr     <- (breaks[-1] + breaks[-(bins + 1L)]) / 2
  res$bin <- pmin(pmax(findInterval(res$tim_med, breaks, rightmost.closed = TRUE), 1L), bins)

  # per-subset y-values (respecting the selected mode) + KS gap
  yvals <- function(cnt_c2, cnt_bg) {
    cum_c2 <- cumsum(cnt_c2); cum_bg <- cumsum(cnt_bg)
    tot_c2 <- sum(cnt_c2);    tot_bg <- sum(cnt_bg)
    sc2 <- if (tot_c2 > 0) tot_c2 else 1; sbg <- if (tot_bg > 0) tot_bg else 1
    ks  <- if (tot_c2 > 0 && tot_bg > 0) max(abs(cum_c2 / tot_c2 - cum_bg / tot_bg)) else NA_real_
    v <- switch(mode,
      cumulative      = list(c2 = cum_c2,      bg = cum_bg,      yt = "Cumulative traces detected"),
      norm            = list(c2 = cnt_c2 / sc2, bg = cnt_bg / sbg, yt = "Detection fraction per bin"),
      norm_cumulative = list(c2 = cum_c2 / sc2, bg = cum_bg / sbg, yt = "Cumulative detection fraction"),
                        list(c2 = cnt_c2,      bg = cnt_bg,      yt = "Traces detected per bin"))
    c(v, list(ks = ks))
  }

  fig    <- plotly::plot_ly()
  ytitle <- "Traces detected per bin"
  ttl    <- NULL

  if (!grouped) {
    cnt_c2 <- tabulate(res$bin[res$in_c2],  nbins = bins)
    cnt_bg <- tabulate(res$bin[!res$in_c2], nbins = bins)
    yy <- yvals(cnt_c2, cnt_bg); ytitle <- yy$yt
    fig <- fig %>%
      plotly::add_trace(x = ctr, y = yy$c2, type = "scatter", mode = "lines+markers",
                        line = list(color = "#2C7FB8"), marker = list(color = "#2C7FB8"), name = "In Median Cluster") %>%
      plotly::add_trace(x = ctr, y = yy$bg, type = "scatter", mode = "lines+markers",
                        line = list(color = "#9E9E9E"), marker = list(color = "#9E9E9E"), name = "Unclustered")
    if (normalized && is.finite(yy$ks))
      ttl <- list(text = sprintf("Detection-timing discrepancy (KS gap) = %.3f", yy$ks), font = list(size = 13))
  } else {
    gs <- unique(res$dataset); ng <- length(gs)
    pal <- series_palette(gs)
    for (g in gs) {
      sub    <- res[res$dataset == g, , drop = FALSE]
      cnt_c2 <- tabulate(sub$bin[sub$in_c2],  nbins = bins)
      cnt_bg <- tabulate(sub$bin[!sub$in_c2], nbins = bins)
      yy <- yvals(cnt_c2, cnt_bg); ytitle <- yy$yt
      kslab <- if (normalized && is.finite(yy$ks)) sprintf(" (KS %.2f)", yy$ks) else ""
      fig <- fig %>%
        plotly::add_trace(x = ctr, y = yy$c2, type = "scatter", mode = "lines",
                          line = list(color = pal[g], width = 2),
                          name = paste0(g, " - In Median Cluster", kslab)) %>%
        plotly::add_trace(x = ctr, y = yy$bg, type = "scatter", mode = "lines",
                          line = list(color = pal[g], width = 2, dash = "dash"),
                          name = paste0(g, " - Unclustered"))
    }
  }
  fig %>% plotly::layout(
    title  = ttl,
    xaxis  = list(title = "Acquisition time (tim)"),
    yaxis  = list(title = ytitle),
    legend = list(orientation = "h")
  )
}
output$det_time_plot <- renderPlotly(build_det_time_plot_fig())

output$download_det_time <- downloadHandler(
  filename = function()
    paste0("cluster2_detection_over_time_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
  content = function(file) {
    res <- det_time_results()
    shiny::validate(shiny::need(!is.null(res) && nrow(res) > 0, "No results."))
    write.csv(res, file, row.names = FALSE)
  }
)

# ==========================================================================
# Trace-ID detection over time - unique traces vs acquisition time
# ==========================================================================
# GOAL: show how many UNIQUE TRACE IDs are detected as a function of
# acquisition time. A trace ID is one continuous detection of a blinking
# fluorophore, so this is the detection-event RATE over the run - a direct QC
# of whether imaging stayed stable (a declining rate flags photobleaching,
# exhausted blinking, drift, or label/sample depletion).
#
# DIFFERS from the Median-Cluster detection plot above: no clustering is
# required (works on raw traces), and it draws ONE SERIES PER FILE - or per
# group in group mode - rather than a clustered/unclustered split.
#
# Each trace is counted ONCE, at its median tim, so a long trace with many
# localizations does not dominate the early bins. Per-localization metadata
# filters apply (via apply_meta_filters_batch); the trace-length filter is
# deliberately NOT applied, because short blinking events are real detections.
#
# Each dataset has its OWN tim timeline, so each series is drawn only across
# its own acquisition window (see the plot): the line STOPS where that
# dataset's tim ends, rather than dropping to a misleading zero.
#
# RESULT stored: data.frame(tid, tim_med, dataset) - one row per trace.
tid_time_results <- reactiveVal(NULL)

observeEvent(input$run_tid_time, {
  run_set <- smlm_settings_from_input(input)   # one settings snapshot for this run
  store <- smlm_file_store()
  if (!length(store)) { showNotification("No files loaded.", type = "warning"); return() }
  # Active-file-only unless "all files" is ticked.
  if (!isTRUE(input$tid_time_all_files)) {
    af <- input$active_file
    if (!is.null(af) && af %in% names(store)) store <- store[af]
  }

  # Trace-ID length prefilter (localizations per trace), the SAME range applied
  # to every dataset. len_max = 0 means "no upper limit". Read once here.
  len_min <- max(1L, suppressWarnings(as.integer(input$tid_time_len_min %||% 1L)))
  len_max <- suppressWarnings(as.integer(input$tid_time_len_max %||% 0L))
  if (is.na(len_min)) len_min <- 1L
  if (is.na(len_max)) len_max <- 0L

  all_tid <- list()
  withProgress(message = "Computing trace detection over time...", value = 0, {
    n_files <- length(store)
    for (i in seq_along(store)) {
      fname <- names(store)[[i]]
      incProgress(1 / n_files, detail = fname)
      df <- qc_exclude_from(store[[fname]], fname)
      if (!is.data.frame(df) || !nrow(df)) next
      # Row-level metadata filters only (efo/cfr/dcr/...); NOT spread/length.
      df <- tryCatch(batch_filters(df, run_set), error = function(e) NULL)
      if (is.null(df) || !nrow(df) || !all(c("tid", "tim") %in% names(df))) next
      df <- apply_density_thinning(df, fname, run_set)   # opt-in labelling-density match
      if (is.null(df) || !nrow(df)) next
      # Trace-length prefilter: keep only trace IDs whose localization count is
      # in [len_min, len_max]. Counts are per-file (a tid is within one file).
      tl <- table(df$tid)
      keep_tids <- names(tl)[tl >= len_min & (len_max <= 0L | tl <= len_max)]
      df <- df[df$tid %in% keep_tids, , drop = FALSE]
      if (!nrow(df)) next
      # One row per trace, placed at its median acquisition time.
      per_tid <- df %>%
        dplyr::group_by(tid) %>%
        dplyr::summarise(tim_med = stats::median(tim, na.rm = TRUE), .groups = "drop") %>%
        dplyr::filter(is.finite(tim_med)) %>%
        dplyr::mutate(dataset = fname)
      if (nrow(per_tid)) all_tid[[fname]] <- per_tid
    }
  })
  res <- do.call(rbind, Filter(Negate(is.null), all_tid))
  if (is.null(res) || !nrow(res)) {
    showNotification("No traces with a tim column found. Check the metadata filters.",
                     type = "warning")
    tid_time_results(NULL); return()
  }
  tid_time_results(res)
  showNotification(
    sprintf("Trace detection over time: %d traces across %d file(s).",
            nrow(res), length(unique(res$dataset))),
    type = "message", duration = 4)
})

build_tid_time_plot_fig <- function() {
  res <- tid_time_results()
  shiny::validate(shiny::need(!is.null(res) && nrow(res) > 0,
                              "Press 'Run' to compute trace detection over time."))
  # Relabel file -> group when in group mode (drops unassigned files).
  res <- apply_dataset_groups(res)
  shiny::validate(shiny::need(nrow(res) > 0, "No files assigned to groups (group mode)."))

  bins <- max(5L, as.integer(input$tid_time_bins %||% 40L))
  mode <- input$tid_time_yaxis %||% "count"

  # Shared time axis (union of all datasets) so files/groups line up on x.
  rng <- range(res$tim_med, na.rm = TRUE)
  if (!all(is.finite(rng)) || diff(rng) <= 0) rng <- c(rng[1], rng[1] + 1)
  breaks  <- seq(rng[1], rng[2], length.out = bins + 1L)
  ctr     <- (breaks[-1] + breaks[-(bins + 1L)]) / 2
  res$bin <- pmin(pmax(findInterval(res$tim_med, breaks, rightmost.closed = TRUE), 1L), bins)

  # Per-series y-values honouring the selected mode. `active` marks the bins
  # that fall inside THIS dataset's own tim window; bins outside it are set to
  # NA so the line STOPS at the dataset's end instead of dropping to zero
  # (zero would wrongly read as "detected nothing" rather than "acquisition
  # over"). A genuine empty bin WITHIN the window stays 0 - a real dip.
  yvals <- function(cnt, active) {
    tot <- sum(cnt); sc <- if (tot > 0) tot else 1
    y <- switch(mode,
      cumulative      = cumsum(cnt),
      norm            = cnt / sc,
      norm_cumulative = cumsum(cnt) / sc,
                        cnt)
    y[!active] <- NA_real_
    yt <- switch(mode,
      cumulative      = "Cumulative unique traces",
      norm            = "Detection fraction per bin",
      norm_cumulative = "Cumulative detection fraction",
                        "Unique traces per bin")
    list(y = y, yt = yt)
  }

  series <- unique(res$dataset)
  pal    <- series_palette(series)
  fig    <- plotly::plot_ly()
  ytitle <- "Unique traces per bin"
  for (s in series) {
    tm  <- res$tim_med[res$dataset == s]
    cnt <- tabulate(res$bin[res$dataset == s], nbins = bins)
    # Bins spanned by this dataset's own [min, max] tim.
    b_lo   <- pmin(pmax(findInterval(min(tm), breaks, rightmost.closed = TRUE), 1L), bins)
    b_hi   <- pmin(pmax(findInterval(max(tm), breaks, rightmost.closed = TRUE), 1L), bins)
    active <- seq_len(bins) >= b_lo & seq_len(bins) <= b_hi
    yy  <- yvals(cnt, active); ytitle <- yy$yt
    fig <- fig %>% plotly::add_trace(
      x = ctr, y = yy$y, type = "scatter", mode = "lines+markers",
      line = list(color = pal[s], width = 2), marker = list(color = pal[s], size = 5),
      connectgaps = FALSE, name = s)
  }
  fig %>% plotly::layout(
    xaxis  = list(title = "Acquisition time (tim)"),
    yaxis  = list(title = ytitle),
    legend = list(orientation = "h")
  )
}
output$tid_time_plot <- renderPlotly(build_tid_time_plot_fig())

output$download_tid_time <- downloadHandler(
  filename = function()
    paste0("trace_detection_over_time_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
  content = function(file) {
    res <- tid_time_results()
    shiny::validate(shiny::need(!is.null(res) && nrow(res) > 0, "No results."))
    write.csv(res, file, row.names = FALSE)
  }
)
