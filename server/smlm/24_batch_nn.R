## ─────────────────────────────────────────────────────────────────────────────
## server/smlm/24_batch_nn.R  -  panel: batch nearest-neighbour histogram
## ─────────────────────────────────────────────────────────────────────────────
## Part of smlm_server(). Not a standalone script: its statements are evaluated in
## the server function's environment, in file-name order (server/parts.R).
##
## Run handler, filtered view, figure, statistics and download.
##
## INTERFACE (generated: Rscript tools/server_map.R --headers; do not edit by hand)
##   DEFINES for other parts:
##     build_nn_batch_hist_fig, nn_batch_params_live, nn_batch_results
##   READS from other parts:
##     01_view_upload: smlm_file_store, uploaded_df_rv
##     02_multiplex: mplex_split_active
##     20_dataset_groups: apply_dataset_groups, file_group_map, file_palette
##     21_density_normalisation: apply_density_thinning
##     23_prepared_cache: resolve_unit
##     27_fiducial_qc: qc_exclude_from
## END INTERFACE
## ─────────────────────────────────────────────────────────────────────────────

# ------ Batch NN Distance Histogram ----------------------------------------
nn_batch_results <- reactiveVal(NULL)
nn_batch_params  <- reactiveVal(NULL)   # RESOLVED compute params, frozen at Run
nn_batch_params_live <- reactiveVal(NULL)  # raw selector values, for the stale badge

# Two-color group filter - only visible when two-color is enabled
output$nn_batch_group_ui <- renderUI({
  if (!isTRUE(input$enable_two_color)) return(NULL)
  selectInput("nn_batch_group", "Analysis group",
              choices = c("All"               = "all",
                          "Within Group 1"    = "within_1",
                          "Within Group 2"    = "within_2",
                          "Between Ch1 → Ch2" = "between_12",
                          "Between Ch2 → Ch1" = "between_21"),
              selected = "all")
})

observeEvent(input$run_batch_nn, {
  run_set <- smlm_settings_from_input(input)   # one settings snapshot for this run
  unit  <- resolve_unit(input$nn_batch_unit, "pcid", "tid", "Batch NN")
  k     <- max(1L, as.integer(input$nn_batch_k %||% 1L))
  gf    <- input$nn_batch_group %||% "all"
  group_filter <- switch(gf, "within_1" = 1L, "within_2" = 2L, NA_integer_)
  cross_groups <- switch(gf,
    "between_12" = c(1L, 2L), "between_21" = c(2L, 1L), NULL)

  store <- smlm_file_store()
  if (!length(store)) {
    showNotification("No files loaded.", type = "warning"); return()
  }

  # ── Multiplex NN modes ──────────────────────────────────────────────────
  adf      <- uploaded_df_rv()
  mplex_on <- !is.null(adf) && "mplex_label" %in% names(adf)
  mmode    <- if (mplex_on) (input$nn_batch_mplex_mode %||% "split") else "off"
  chx      <- input$nn_batch_chx %||% ""
  chy      <- input$nn_batch_chy %||% ""
  mplex_groups <- FALSE   # analysis_group derived from channel membership (between)

  if (mmode == "split") {
    cs <- mplex_split_active()
    if (length(cs)) store <- cs
  } else if (mmode == "within") {
    if (!nzchar(chx)) { showNotification("Pick Channel X.", type = "warning"); return() }
    sub <- adf[!is.na(adf$mplex_label) & adf$mplex_label == chx, , drop = FALSE]
    if (!nrow(sub)) { showNotification("Channel X has no data.", type = "warning"); return() }
    store <- stats::setNames(list(sub), chx)
    group_filter <- NA_integer_; cross_groups <- NULL
  } else if (mmode == "between") {
    if (!nzchar(chx) || !nzchar(chy) || identical(chx, chy)) {
      showNotification("Pick two distinct channels (X and Y).", type = "warning"); return()
    }
    store <- stats::setNames(list(adf), paste0(chx, " → ", chy))
    group_filter <- NA_integer_; cross_groups <- c(1L, 2L)
    mplex_groups <- TRUE
  }

  results   <- list()
  n_skipped <- 0L
  withProgress(message = "Running batch NN...", value = 0, {
    n_files <- length(store)
    for (i in seq_along(store)) {
      fname <- names(store)[[i]]
      incProgress(1 / n_files, detail = fname)
      df <- qc_exclude_from(store[[fname]], fname)
      if (!is.data.frame(df) || !nrow(df)) { n_skipped <- n_skipped + 1L; next }

      df <- tryCatch(apply_meta_filters_batch(df, run_set), error = function(e) NULL)
      if (is.null(df) || !nrow(df)) { n_skipped <- n_skipped + 1L; next }
      df <- apply_density_thinning(df, fname, run_set)   # opt-in labelling-density match
      if (is.null(df) || !nrow(df)) { n_skipped <- n_skipped + 1L; next }

      got <- batch_nn_for_file(df, fname, run_set, unit = unit, k = k,
                               group_filter = group_filter, cross_groups = cross_groups,
                               mplex_pair = if (mplex_groups) c(chx, chy) else NULL,
                               by_channel = isTRUE(input$analysis_by_channel %||% FALSE))
      if (!length(got)) { n_skipped <- n_skipped + 1L; next }
      results[names(got)] <- got
    }
  })

  if (!length(results)) {
    hint <- switch(unit,
      pcid          = "Check that Primary point DBSCAN Clustering is enabled and eps/minPts produce clusters.",
      cluster2_tid  = "Check that Median Clustering is enabled.",
      cluster2_pcid = "Check that both clustering levels are enabled.",
      "Check filters - all rows may have been filtered out."
    )
    showNotification(paste0("No distances computed. ", hint), type = "warning")
    nn_batch_results(NULL)
    return()
  }
  nn_batch_results(do.call(rbind, results))
  # Snapshot the compute parameters so the plot labels the data it is actually
  # showing. Reading them live meant changing the Unit selector relabelled an
  # old plot without recomputing it - and with unit = "auto" the axis read
  # "[auto]" because the label switch has no such case.
  nn_batch_params(list(unit = unit, k = k))
  nn_batch_params_live(list(unit = input$nn_batch_unit, k = input$nn_batch_k))
  msg <- sprintf("NN computed across %d of %d dataset(s).", length(results), length(store))
  if (n_skipped > 0)
    msg <- paste0(msg, sprintf(" %d file(s) yielded no distances (filtered out or no clusters).", n_skipped))
  showNotification(msg, type = "message", duration = 5)
})

# Shared reactive: stored results filtered by the distance range cutoff.
# Both the histogram and the stats table read from this so they stay in sync.
nn_batch_filtered <- reactive({
  res <- nn_batch_results()
  if (is.null(res) || !nrow(res)) return(NULL)
  dr <- input$nn_batch_dist_range %||% c(0, Inf)
  lo <- if (is.finite(dr[1])) dr[1] else 0
  hi <- if (is.finite(dr[2])) dr[2] else Inf
  res <- res[res$distance >= lo & res$distance <= hi, , drop = FALSE]
  # Keep the source FILE name alongside the (possibly group-relabelled) dataset,
  # so the plot can draw one series per file and colour it by its group.
  if (nrow(res)) res$file <- res$dataset
  apply_dataset_groups(res)   # relabel dataset -> group in group mode
})

build_nn_batch_hist_fig <- function() {
  res <- nn_batch_filtered()
  shiny::validate(shiny::need(!is.null(res) && nrow(res) > 0,
                              "Press 'Run Batch NN' to generate histogram (or widen the distance filter)."))

  bins       <- max(5L, as.integer(input$nn_batch_bins %||% 60L))
  density    <- identical(input$nn_batch_yaxis %||% "count", "density")
  error_type <- input$nn_batch_error %||% "iqr"
  aggregate  <- isTRUE(input$nn_batch_aggregate)
  col_g1     <- input$color_group1 %||% "#1f77b4"
  col_g2     <- input$color_group2 %||% "#d62728"
  group_col  <- switch(input$nn_batch_group %||% "all",
    within_1   = col_g1,
    within_2   = col_g2,
    between_12 = col_g1,
    between_21 = col_g2,
    "#2166ac"
  )

  # shared bins across all series so overlaid histograms + indicator lines align
  rng_all <- range(res$distance, na.rm = TRUE)
  if (!all(is.finite(rng_all)) || diff(rng_all) <= 0) rng_all <- c(rng_all[1], rng_all[1] + 1)
  bsize   <- diff(rng_all) / bins
  brk     <- seq(rng_all[1], rng_all[2], length.out = bins + 1L)
  xbins   <- list(start = rng_all[1], end = rng_all[2], size = bsize)
  bin_hts <- function(x) {
    h <- graphics::hist(x, breaks = brk, plot = FALSE, include.lowest = TRUE)
    if (density) h$density else h$counts
  }
  subsets <- if (aggregate) list(res$distance) else split(res$distance, res$dataset)
  ytop <- suppressWarnings(max(vapply(subsets,
            function(x) if (length(x)) max(bin_hts(x)) else 0, numeric(1)), na.rm = TRUE))
  if (!is.finite(ytop) || ytop <= 0) ytop <- 1

  # median + IQR/MAD as vertical LINE TRACES sharing the histogram's legendgroup
  # (so they toggle away with the group when its legend entry is clicked)
  add_spread <- function(fig, d, col, lg) {
    med <- stats::median(d)
    sv  <- switch(error_type,
      iqr = c(stats::quantile(d, 0.25, names = FALSE), stats::quantile(d, 0.75, names = FALSE)),
      mad = { mv <- stats::mad(d, constant = 1); c(max(0, med - mv), med + mv) },
      NULL)
    fig <- fig %>% plotly::add_trace(x = c(med, med), y = c(0, ytop),
      type = "scatter", mode = "lines", legendgroup = lg, showlegend = FALSE,
      hoverinfo = "skip", opacity = 0.9, line = list(color = col, width = 2))
    for (sv_i in sv) fig <- fig %>% plotly::add_trace(x = c(sv_i, sv_i), y = c(0, ytop),
      type = "scatter", mode = "lines", legendgroup = lg, showlegend = FALSE,
      hoverinfo = "skip", opacity = 0.65, line = list(color = col, width = 1.2, dash = "dot"))
    fig
  }

  fig <- plotly::plot_ly()

  if (aggregate) {
    fig <- fig %>% plotly::add_histogram(
      x        = res$distance,
      name     = sprintf("All datasets (N = %d)", nrow(res)),
      legendgroup = "all",
      marker   = list(color = paste0(group_col, "99")),
      xbins    = xbins,
      histnorm = if (density) "probability density" else "",
      opacity  = 0.75
    )
    fig <- add_spread(fig, res$distance, group_col, "all")
  } else {
    grouped_mode <- (input$analysis_by %||% "file") == "group"
    # ONE SERIES PER FILE in both modes. Colour encodes what distinguishes them:
    #   file mode  - a categorical colour per file (the old code shaded ONE blue
    #                by alpha, which made every file a shade of blue);
    #   group mode - the hue is the group and the shade is the file within it, so
    #                the conditions separate at a glance while replicates stay
    #                individually visible.
    files <- unique(res$file %||% res$dataset)
    fpal  <- file_palette(files)
    # Draw grouped by condition so the legend lists a group's files together.
    if (grouped_mode) files <- files[order(match(files, names(file_group_map())))]

    for (fl in files) {
      sub <- res[(res$file %||% res$dataset) == fl, , drop = FALSE]
      if (!nrow(sub)) next
      col_base <- unname(fpal[fl])
      if (is.na(col_base)) col_base <- group_col
      lgroup   <- if (grouped_mode) sub$dataset[[1L]] else fl
      lbl      <- if (grouped_mode) sprintf("%s - %s", sub$dataset[[1L]], fl)
                  else if ("tag" %in% names(sub)) sub$tag[[1L]] else fl
      fig <- fig %>% plotly::add_histogram(
        x        = sub$distance,
        name     = lbl,
        legendgroup = lgroup,
        marker   = list(color = paste0(col_base, "99")),
        xbins    = xbins,
        histnorm = if (density) "probability density" else "",
        opacity  = 0.75
      )
      fig <- add_spread(fig, sub$distance, col_base, lgroup)
    }
  }

  indicator_note <- switch(error_type,
    iqr  = "- median  ··· Q25 / Q75",
    mad  = "- median  ··· ±MAD",
    "- median"
  )
  # From the frozen snapshot (see nn_batch_params), never the live selector.
  npar <- nn_batch_params()
  u_used <- npar$unit %||% "tid"
  k_used <- npar$k    %||% 1L
  unit_lbl <- switch(u_used,
    tid           = "TID centroid",
    pcid          = "De-blinked molecule centroid",
    cluster2_tid  = "TID-within-Median-Cluster",
    cluster2_pcid = "De-blinked molecule centroid within Median Cluster",
    u_used
  )
  fig %>% plotly::layout(
    barmode = "overlay",
    xaxis   = list(title = sprintf("k=%d NN distance (nm)  [%s]",
                                   as.integer(k_used), unit_lbl)),
    yaxis   = list(title = if (density) "Density" else "Count"),
    legend  = list(orientation = "h", y = -0.18, groupclick = "togglegroup"),
    annotations = list(list(
      xref = "paper", yref = "paper", x = 0.01, y = 0.99,
      xanchor = "left", yanchor = "top", showarrow = FALSE,
      text = indicator_note, font = list(size = 11, color = "grey40")
    )),
    paper_bgcolor = "white",
    plot_bgcolor  = "white"
  )
}
output$nn_batch_hist <- renderPlotly(build_nn_batch_hist_fig())

output$nn_batch_stats <- renderTable({
  res        <- nn_batch_filtered()
  error_type <- input$nn_batch_error %||% "iqr"
  aggregate  <- isTRUE(input$nn_batch_aggregate)
  if (is.null(res) || !nrow(res)) return(NULL)

  make_stats_row <- function(d, file_lbl, ds_lbl) {
    med <- stats::median(d)
    spread <- switch(error_type,
      iqr = {
        q25 <- stats::quantile(d, 0.25, names = FALSE)
        q75 <- stats::quantile(d, 0.75, names = FALSE)
        data.frame(Q25 = round(q25, 1), Q75 = round(q75, 1),
                   IQR = round(q75 - q25, 1), check.names = FALSE)
      },
      mad = {
        mv <- stats::mad(d, constant = 1)
        data.frame(`Median - MAD` = round(med - mv, 1),
                   `Median + MAD` = round(med + mv, 1),
                   MAD = round(mv, 1), check.names = FALSE)
      },
      data.frame()
    )
    cbind(data.frame(File = file_lbl, Dataset = ds_lbl, N = length(d),
                     `Median (nm)` = round(med, 1), check.names = FALSE), spread)
  }

  if (aggregate) {
    return(make_stats_row(res$distance, "All", "Aggregated"))
  }

  datasets <- unique(res$dataset)
  rows <- lapply(datasets, function(ds) {
    sub <- res[res$dataset == ds, , drop = FALSE]
    lbl <- if ("tag" %in% names(sub)) sub$tag[[1L]] else ds
    make_stats_row(sub$distance, ds, lbl)
  })
  do.call(rbind, rows)
}, striped = TRUE, hover = TRUE, spacing = "s", digits = 1)

output$download_batch_nn <- downloadHandler(
  filename = function()
    paste0("batch_nn_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
  content = function(file) {
    res <- isolate(nn_batch_filtered())
    shiny::validate(shiny::need(!is.null(res) && nrow(res) > 0, "No NN results."))
    keep <- intersect(c("dataset", "tag", "distance"), names(res))
    write.csv(res[, keep, drop = FALSE], file, row.names = FALSE)
  }
)
