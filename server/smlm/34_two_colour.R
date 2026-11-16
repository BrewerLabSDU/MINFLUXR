# ── Two-Color: Cross-Channel NN ──────────────────────────────────────────
cross_nn_results <- reactiveVal(NULL)

observeEvent(input$run_cross_nn, {
  run_set <- smlm_settings_from_input(input)   # one settings snapshot for this run
  store <- smlm_file_store()
  if (!length(store)) { showNotification("No datasets loaded.", type = "warning"); return() }

  unit <- input$cross_nn_unit %||% "tid"
  k    <- max(1L, as.integer(input$cross_nn_k %||% 1L))
  needs_cluster  <- unit %in% c("cluster", "cluster2_cluster")
  needs_cluster2 <- unit %in% c("cluster2_tid", "cluster2_cluster", "cluster2_centroid")
  base_unit <- if (unit == "cluster2_centroid") "cluster2"
               else if (unit %in% c("tid", "cluster2_tid")) "tid" else "cluster"

  two_color_mode <- isTRUE(input$enable_two_color)
  mplex_pair     <- mplex_cross_pair()

  results <- list()
  withProgress(message = "Cross-channel NN...", value = 0, {
    for (i in seq_along(store)) {
      incProgress(1 / length(store))
      nm <- names(store)[i]
      df <- apply_meta_filters_batch(store[[nm]], run_set)
      if (!nrow(df)) next

      # A chosen multiplex pair is collapsed to analysis_group 1/2 FIRST;
      # resolve_channels() then preserves that choice rather than re-deriving
      # it from channel order. cluster2 units need the chain run even when the
      # counted unit is the raw trace, so it is not conditional here.
      if (!is.null(mplex_pair))
        df <- assign_analysis_group_from_mplex(df, mplex_pair[1], mplex_pair[2])
      if (nrow(df)) df <- apply_clustering_pipeline(df, run_set)
      if (!nrow(df)) next
      if (!all(c(1L, 2L) %in% df$analysis_group)) next

      # ── MEDIAN-CLUSTER MODES RESTRICT, THEY DO NOT MIX ──────────────────
      # A Median Cluster belongs to ONE channel by design - the 2-colour
      # pipeline clusters each group separately, and that separation is the
      # point of the 2-colour split. So these modes do NOT pair molecules
      # inside a shared structure (no such structure exists, and building one
      # would merge the two fluorophore populations). They keep every distance
      # strictly between channels and drop the unclustered background from
      # BOTH sides, so the question asked is "how far is a group-1 STRUCTURE
      # from group-2 material" rather than "...from any stray group-2 point".
      if (needs_cluster2) {
        if (!"cluster2" %in% names(df) || !any(!is.na(df$cluster2) & df$cluster2 > 0L)) next
        df <- dplyr::filter(df, !is.na(cluster2) & cluster2 > 0L)
        if (!nrow(df) || !all(c(1L, 2L) %in% df$analysis_group)) next
      }
      d <- compute_cross_nn(df, unit = base_unit, k = k)
      if (!is.null(d) && nrow(d)) results[[nm]] <- cbind(dataset = nm, d)
    }
  })

  if (!length(results)) {
    hint <- if (needs_cluster2)
      sprintf(" One or both channels have no Median Cluster structures - check that %s produces clusters (eps/MinPts).", cluster2_label())
    else if (needs_cluster)
      sprintf(" The molecule units need %s to be enabled.", deblink_label())
    else ""
    showNotification(paste0(
      "No cross-channel data found. Check that both channels are present and analysis_group is set correctly.",
      hint), type = "warning")
    cross_nn_results(NULL)
    return()
  }
  cross_nn_results(do.call(rbind, results))
})

build_cross_nn_hist_fig <- function() {
  res <- cross_nn_results()
  shiny::validate(shiny::need(!is.null(res) && nrow(res) > 0, "Run cross-channel NN first."))

  aggregate   <- isTRUE(input$cross_nn_aggregate)
  hist_type   <- input$cross_nn_hist_type %||% "density"
  nbins       <- as.integer(input$cross_nn_bins %||% 40L)
  histnorm    <- if (hist_type == "density") "probability density" else ""
  col_g1      <- input$color_group1 %||% "#1f77b4"
  col_g2      <- input$color_group2 %||% "#d62728"
  grp_colors  <- c("Group 1 → Group 2" = col_g1, "Group 2 → Group 1" = col_g2)

  fig <- plotly::plot_ly()

  if (aggregate) {
    for (dir in c("Group 1 → Group 2", "Group 2 → Group 1")) {
      d <- res$distance[res$direction == dir]
      if (!length(d)) next
      fig <- fig %>% plotly::add_histogram(
        x = d, name = dir, histnorm = histnorm, nbinsx = nbins,
        marker = list(color = grp_colors[dir], opacity = 0.7),
        autobinx = FALSE
      )
    }
  } else {
    datasets <- unique(res$dataset)
    n_ds     <- length(datasets)
    alphas   <- if (n_ds == 1L) "cc"
                else sprintf("%02x", as.integer(seq(204, 80, length.out = n_ds)))
    names(alphas) <- datasets
    for (nm in datasets) {
      for (dir in c("Group 1 → Group 2", "Group 2 → Group 1")) {
        d <- res$distance[res$dataset == nm & res$direction == dir]
        if (!length(d)) next
        base_col <- grp_colors[dir]
        fig <- fig %>% plotly::add_histogram(
          x = d, name = sprintf("%s | %s", nm, dir),
          histnorm = histnorm, nbinsx = nbins,
          marker = list(color = paste0(base_col, alphas[nm]), opacity = 0.75),
          autobinx = FALSE
        )
      }
    }
  }

  # Median marker per drawn series. The cross-NN distribution is strongly
  # right-skewed, so the median (not the mode of the tallest bin) is the number
  # that summarises it.
  # NOTE when not aggregating: one line per (dataset x direction) would be
  # unreadable on many files, so the line marks the POOLED median of each
  # direction over the files shown - which is not the median of the per-file
  # medians, and is weighted by how many points each file contributed.
  if (isTRUE(input$cross_nn_median_lines %||% TRUE)) {
    meds <- vapply(c("Group 1 → Group 2", "Group 2 → Group 1"), function(dir)
      stats::median(res$distance[res$direction == dir], na.rm = TRUE), numeric(1))
    fig <- add_median_vlines(fig, meds, grp_colors, unit = "nm")
  }

  fig %>% plotly::layout(
    barmode = "overlay",
    xaxis = list(title = "Distance (nm)"),
    yaxis = list(title = if (hist_type == "density") "Density" else "Count"),
    title = "Cross-Channel Nearest-Neighbour Distances",
    paper_bgcolor = "white", plot_bgcolor = "white"
  )
}
output$cross_nn_hist <- renderPlotly(build_cross_nn_hist_fig())

build_coloc_curve_fig <- function() {
  res <- cross_nn_results()
  shiny::validate(shiny::need(!is.null(res) && nrow(res) > 0, "Run cross-channel NN first."))

  aggregate   <- isTRUE(input$cross_nn_aggregate)
  col_ch1     <- input$color_group1 %||% "#1f77b4"
  col_ch2     <- input$color_group2 %||% "#d62728"
  ch_colors   <- c("Group 1 → Group 2" = col_ch1, "Group 2 → Group 1" = col_ch2)

  fig <- plotly::plot_ly()

  add_ecdf_trace <- function(d, nm, col, dash = "solid", show_leg = TRUE) {
    if (!length(d)) return(invisible(NULL))
    ord  <- sort(d)
    frac <- seq_along(ord) / length(ord)
    fig <<- fig %>% plotly::add_trace(
      x = ord, y = frac * 100,
      type = "scatter", mode = "lines", name = nm,
      showlegend = show_leg,
      line = list(color = col, width = 2, dash = dash))
  }

  add_null_trace <- function(sub, col) {
    if (!all(c("n_ref", "area") %in% names(sub))) return(invisible(NULL))
    n_ref <- sub$n_ref[1L]; area <- sub$area[1L]
    if (!is.finite(area) || area <= 0 || n_ref < 1L) return(invisible(NULL))
    r_seq <- seq(0, max(sub$distance, na.rm = TRUE), length.out = 200L)
    null_frac <- (1 - pmax(0, 1 - pi * r_seq^2 / area)^n_ref) * 100
    fig <<- fig %>% plotly::add_trace(
      x = r_seq, y = null_frac,
      type = "scatter", mode = "lines", name = "CSR null",
      showlegend = FALSE,
      line = list(color = col, width = 1, dash = "dash"),
      opacity = 0.45)
  }

  if (aggregate) {
    for (dir in c("Group 1 → Group 2", "Group 2 → Group 1")) {
      sub <- res[res$direction == dir, ]
      add_ecdf_trace(sub$distance, dir, ch_colors[dir])
      add_null_trace(sub, ch_colors[dir])
    }
  } else {
    datasets <- unique(res$dataset)
    pal <- if (length(datasets) <= 8) RColorBrewer::brewer.pal(max(3, length(datasets)), "Set2")[seq_along(datasets)]
           else scales::hue_pal()(length(datasets))
    names(pal) <- datasets
    for (nm in datasets) {
      for (dir in c("Group 1 → Group 2", "Group 2 → Group 1")) {
        sub  <- res[res$dataset == nm & res$direction == dir, ]
        dash <- if (dir == "Group 1 → Group 2") "solid" else "dot"
        add_ecdf_trace(sub$distance, sprintf("%s | %s", nm, dir), pal[nm], dash)
        add_null_trace(sub, pal[nm])
      }
    }
  }

  fig %>% plotly::layout(
    xaxis = list(title = "Distance threshold (nm)"),
    yaxis = list(title = "% molecules with neighbour within d", range = c(0, 100)),
    title = "Co-localization Fraction  (solid = observed, dashed = CSR null)",
    hovermode = "x unified",
    paper_bgcolor = "white", plot_bgcolor = "white"
  )
}
output$coloc_curve <- renderPlotly(build_coloc_curve_fig())

output$download_cross_nn <- downloadHandler(
  filename = function() paste0("cross_nn_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
  content = function(file) {
    res <- cross_nn_results()
    shiny::validate(shiny::need(!is.null(res) && nrow(res) > 0, "No results."))
    write.csv(res, file, row.names = FALSE)
  }
)

# ── Two-Color: Bivariate Ripley's L12 ────────────────────────────────────
ripley12_results     <- reactiveVal(NULL)
ripley12_env_results <- reactiveVal(NULL)  # Monte-Carlo CSR envelope (per dataset)

observeEvent(input$run_ripley12, {
  run_set <- smlm_settings_from_input(input)   # one settings snapshot for this run
  store <- smlm_file_store()
  if (!length(store)) { showNotification("No datasets loaded.", type = "warning"); return() }

  unit      <- input$ripley12_unit %||% "tid"
  use_3d    <- isTRUE(input$ripley12_use_3d)
  r_max_raw <- as.numeric(input$ripley12_r_max %||% 0)
  r_steps   <- as.integer(input$ripley12_r_steps %||% 60L)
  want_env  <- isTRUE(input$ripley12_envelope)
  nsim      <- max(19L, min(199L, as.integer(input$ripley12_nsim %||% 99L)))
  env_file_cap <- 6L
  needs_cluster <- unit == "cluster"

  # Prep every file FIRST, so one gap-fill radius can be resolved for the whole
  # set before any curve is computed (see spatial_resolve_close()).
  prepped <- Filter(Negate(is.null), lapply(seq_along(store), function(i) {
    nm <- names(store)[i]
    df <- prep_tc_df(store[[nm]], run_set, needs_cluster)
    if (is.null(df)) return(NULL)
    list(df = df, tag = nm)
  }))
  if (!length(prepped)) {
    showNotification("No file has both channels after filtering.", type = "warning"); return()
  }
  cc <- spatial_gate_close(prepped, use_3d)
  common_close <- cc$close_nm

  results  <- list()
  env_list <- list()
  withProgress(message = "Bivariate Ripley L12...", value = 0, {
    for (i in seq_along(prepped)) {
      incProgress(1 / length(prepped))
      nm <- prepped[[i]]$tag
      df <- prepped[[i]]$df

      r_seq <- if (r_max_raw > 0) seq(0, r_max_raw, length.out = r_steps) else NULL
      pts1 <- get_tc_pts(df, 1L, unit); pts2 <- get_tc_pts(df, 2L, unit)
      if (!nrow(pts1) || !nrow(pts2)) next
      sup <- spatial_support_of(df, use_3d, tag = nm, close_nm = common_close)
      out <- cross_curve_for_points(pts1, pts2,
                                    estimator = "L12", r_seq = r_seq, use_3d = use_3d,
                                    support = sup,
                                    edge_correct = isTRUE(input$spatial_edge_correct),
                                    envelope = want_env && length(env_list) < env_file_cap,
                                    nsim = nsim, null_mode = input$coloc_null_mode %||% "shift",
                                    tag = nm)
      if (is.null(out)) next
      results[[nm]] <- out$curve
      if (!is.null(out$env)) env_list[[nm]] <- out$env
    }
  })

  if (!length(results)) {
    showNotification("No results. Ensure both channels are present.", type = "warning")
    ripley12_results(NULL)
    ripley12_env_results(NULL)
    return()
  }
  ripley12_results(do.call(rbind, results))
  ripley12_env_results(if (length(env_list)) do.call(rbind, env_list) else NULL)
  spatial_report_windows(vapply(prepped, function(it) as.character(it$tag), character(1)))
})

build_ripley12_plot_fig <- function() {
  res <- ripley12_results()
  shiny::validate(shiny::need(!is.null(res) && nrow(res) > 0, "Run Bivariate Ripley L12 first."))

  aggregate <- isTRUE(input$ripley12_aggregate)
  do_smooth <- isTRUE(input$ripley12_smooth)
  span_val  <- as.numeric(input$ripley12_smooth_span %||% 0.2)
  env_res   <- ripley12_env_results()
  smooth_y  <- function(sub) {
    if (!do_smooth || nrow(sub) < 5L) return(sub$Lr_minus_r)
    tryCatch(predict(loess(Lr_minus_r ~ r, data = sub, span = span_val)),
             error = function(e) sub$Lr_minus_r)
  }
  hex2rgba <- function(hex, a) {
    cc <- grDevices::col2rgb(hex); sprintf("rgba(%d,%d,%d,%.2f)", cc[1], cc[2], cc[3], a)
  }
  fig <- plotly::plot_ly()

  if (aggregate) {
    agg <- res %>%
      dplyr::group_by(r) %>%
      dplyr::summarise(Lr_minus_r = mean(Lr_minus_r, na.rm = TRUE), .groups = "drop")
    fig <- fig %>% plotly::add_trace(
      x = agg$r, y = smooth_y(agg),
      type = "scatter", mode = "lines",
      name = "All datasets (mean)",
      line = list(color = "#7B2D8B", width = 2)
    )
  } else {
    datasets <- unique(res$dataset)
    pal <- if (length(datasets) <= 8) RColorBrewer::brewer.pal(max(3, length(datasets)), "Set2")[seq_along(datasets)]
           else scales::hue_pal()(length(datasets))
    names(pal) <- datasets
    for (nm in datasets) {
      sub <- res[res$dataset == nm, ]
      if (!is.null(env_res) && nm %in% env_res$dataset) {
        e <- env_res[env_res$dataset == nm, ]; e <- e[order(e$r), ]
        fig <- fig %>%
          plotly::add_trace(x = e$r, y = e$hi, type = "scatter", mode = "lines",
            line = list(width = 0), legendgroup = nm, showlegend = FALSE, hoverinfo = "skip") %>%
          plotly::add_trace(x = e$r, y = e$lo, type = "scatter", mode = "lines", fill = "tonexty",
            fillcolor = hex2rgba(unname(pal[nm]), 0.15), line = list(width = 0),
            legendgroup = nm, name = paste0(nm, " CSR 95% envelope"),
            showlegend = TRUE, hoverinfo = "skip") %>%
          plotly::add_trace(x = e$r, y = e$null_mean, type = "scatter", mode = "lines",
            line = list(color = unname(pal[nm]), width = 1, dash = "dot"),
            legendgroup = nm, showlegend = FALSE, hoverinfo = "skip")
        fig <- .add_global_band(fig, e, unname(pal[nm]), nm)
      }
      fig <- fig %>% plotly::add_trace(
        x = sub$r, y = smooth_y(sub),
        type = "scatter", mode = "lines", name = nm, legendgroup = nm,
        line = list(color = pal[nm], width = 2)
      )
    }
  }

  csr_label <- if (!is.null(env_res) && !aggregate)
    "theoretical CSR (uncorrected - see envelope)" else "Random (L−r = 0)"
  fig %>%
    plotly::add_trace(
      x = range(res$r), y = c(0, 0),
      type = "scatter", mode = "lines",
      name = csr_label,
      line = list(color = "black", width = 1, dash = "dash"),
      showlegend = TRUE
    ) %>%
    plotly::layout(
      xaxis = list(title = "r (nm)"),
      yaxis = list(title = "L₁₂(r) − r (nm)"),
      annotations = .global_p_annot(env_res),
      title = "Bivariate Ripley’s L₁₂(r)",
      hovermode = "x unified",
      paper_bgcolor = "white", plot_bgcolor = "white"
    )
}
output$ripley12_plot <- renderPlotly(build_ripley12_plot_fig())

output$download_ripley12 <- downloadHandler(
  filename = function() paste0("ripley_L12_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
  content = function(file) {
    res <- ripley12_results()
    shiny::validate(shiny::need(!is.null(res) && nrow(res) > 0, "No results."))
    write.csv(res, file, row.names = FALSE)
  }
)

# ── shared helper: build per-group point sets for two-color analyses ─────

# ── shared helper: prepare data for one file in two-color batch ───────────
#   settings : the run's settings snapshot (smlm_settings_from_input)
prep_tc_df <- function(raw, settings, needs_cluster) {
  df <- apply_meta_filters_batch(raw, settings)
  if (!nrow(df)) return(NULL)
  mplex_pair <- mplex_cross_pair()
  if (!is.null(mplex_pair))
    df <- assign_analysis_group_from_mplex(df, mplex_pair[1], mplex_pair[2])
  if (nrow(df)) df <- apply_clustering_pipeline(df, settings)
  if (!nrow(df) || !all(c(1L, 2L) %in% df$analysis_group)) return(NULL)
  df
}

# ── Pair Cross-Correlation g₁₂(r) ────────────────────────────────────────
pcf12_results     <- reactiveVal(NULL)
pcf12_env_results <- reactiveVal(NULL)  # Monte-Carlo CSR envelope (per dataset)

observeEvent(input$run_pcf12, {
  run_set <- smlm_settings_from_input(input)   # one settings snapshot for this run
  store <- smlm_file_store()
  if (!length(store)) { showNotification("No datasets loaded.", type = "warning"); return() }

  unit          <- input$pcf12_unit %||% "tid"
  use_3d        <- isTRUE(input$pcf12_use_3d)
  r_max_raw     <- as.numeric(input$pcf12_r_max %||% 0)
  r_steps       <- as.integer(input$pcf12_r_steps %||% 60L)
  want_env      <- isTRUE(input$pcf12_envelope)
  nsim          <- max(19L, min(199L, as.integer(input$pcf12_nsim %||% 99L)))
  env_file_cap  <- 6L
  needs_cluster <- unit == "cluster"

  # Prep first so one gap-fill radius covers the whole set (spatial_resolve_close()).
  prepped <- Filter(Negate(is.null), lapply(seq_along(store), function(i) {
    nm <- names(store)[i]
    df <- prep_tc_df(store[[nm]], run_set, needs_cluster)
    if (is.null(df)) return(NULL)
    list(df = df, tag = nm)
  }))
  if (!length(prepped)) {
    showNotification("No file usable after filtering.", type = "warning"); return()
  }
  cc <- spatial_gate_close(prepped, use_3d)
  common_close <- cc$close_nm

  results  <- list()
  env_list <- list()
  withProgress(message = "Pair cross-correlation g₁₂...", value = 0, {
    for (i in seq_along(prepped)) {
      incProgress(1 / length(prepped))
      nm <- prepped[[i]]$tag
      df <- prepped[[i]]$df
      r_seq <- if (r_max_raw > 0) seq(0, r_max_raw, length.out = r_steps) else NULL
      pts1 <- get_tc_pts(df, 1L, unit); pts2 <- get_tc_pts(df, 2L, unit)
      if (!nrow(pts1) || !nrow(pts2)) next
      sup   <- spatial_support_of(df, use_3d, tag = nm, close_nm = common_close)
      out <- cross_curve_for_points(pts1, pts2,
                                    estimator = "g12", r_seq = r_seq, use_3d = use_3d,
                                    support = sup,
                                    edge_correct = isTRUE(input$spatial_edge_correct),
                                    envelope = want_env && length(env_list) < env_file_cap,
                                    nsim = nsim, null_mode = input$coloc_null_mode %||% "shift",
                                    tag = nm)
      if (is.null(out)) next
      results[[nm]] <- out$curve
      if (!is.null(out$env)) env_list[[nm]] <- out$env
    }
  })

  if (!length(results)) {
    showNotification("No results. Check that both channels are present.", type = "warning")
    pcf12_results(NULL); pcf12_env_results(NULL); return()
  }
  pcf12_results(do.call(rbind, results))
  pcf12_env_results(if (length(env_list)) do.call(rbind, env_list) else NULL)
  spatial_report_windows(vapply(prepped, function(it) as.character(it$tag), character(1)))
})

build_pcf12_plot_fig <- function() {
  res <- pcf12_results()
  shiny::validate(shiny::need(!is.null(res) && nrow(res) > 0, "Run g₁₂ first."))

  aggregate <- isTRUE(input$pcf12_aggregate)
  env_res   <- pcf12_env_results()
  hex2rgba  <- function(hex, a) {
    cc <- grDevices::col2rgb(hex); sprintf("rgba(%d,%d,%d,%.2f)", cc[1], cc[2], cc[3], a)
  }
  fig <- plotly::plot_ly()

  if (aggregate) {
    agg <- res %>%
      dplyr::group_by(r) %>%
      dplyr::summarise(g12 = mean(g12, na.rm = TRUE), .groups = "drop")
    fig <- fig %>% plotly::add_trace(x = agg$r, y = agg$g12,
      type = "scatter", mode = "lines", name = "All datasets (mean)",
      line = list(color = "#7B2D8B", width = 2))
  } else {
    datasets <- unique(res$dataset)
    pal <- if (length(datasets) <= 8) RColorBrewer::brewer.pal(max(3, length(datasets)), "Set2")[seq_along(datasets)]
           else scales::hue_pal()(length(datasets))
    names(pal) <- datasets
    for (nm in datasets) {
      sub <- res[res$dataset == nm, ]
      if (!is.null(env_res) && nm %in% env_res$dataset) {
        e <- env_res[env_res$dataset == nm, ]; e <- e[order(e$r), ]
        fig <- fig %>%
          plotly::add_trace(x = e$r, y = e$hi, type = "scatter", mode = "lines",
            line = list(width = 0), legendgroup = nm, showlegend = FALSE, hoverinfo = "skip") %>%
          plotly::add_trace(x = e$r, y = e$lo, type = "scatter", mode = "lines", fill = "tonexty",
            fillcolor = hex2rgba(unname(pal[nm]), 0.15), line = list(width = 0),
            legendgroup = nm, name = paste0(nm, " CSR 95% envelope"),
            showlegend = TRUE, hoverinfo = "skip") %>%
          plotly::add_trace(x = e$r, y = e$null_mean, type = "scatter", mode = "lines",
            line = list(color = unname(pal[nm]), width = 1, dash = "dot"),
            legendgroup = nm, showlegend = FALSE, hoverinfo = "skip")
        fig <- .add_global_band(fig, e, unname(pal[nm]), nm)
      }
      fig <- fig %>% plotly::add_trace(x = sub$r, y = sub$g12,
        type = "scatter", mode = "lines", name = nm, legendgroup = nm,
        line = list(color = pal[nm], width = 2))
    }
  }

  fig %>%
    plotly::add_trace(x = range(res$r, na.rm = TRUE), y = c(1, 1),
      type = "scatter", mode = "lines", name = "CSR (g = 1)",
      line = list(color = "black", width = 1, dash = "dash"), showlegend = TRUE) %>%
    plotly::layout(
      xaxis = list(title = "r (nm)"),
      yaxis = list(title = "g₁₂(r)"),
      annotations = .global_p_annot(env_res),
      title = "Pair Cross-Correlation Function g₁₂(r)",
      hovermode = "x unified",
      paper_bgcolor = "white", plot_bgcolor = "white"
    )
}
output$pcf12_plot <- renderPlotly(build_pcf12_plot_fig())

output$download_pcf12 <- downloadHandler(
  filename = function() paste0("pcf12_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
  content  = function(file) {
    res <- pcf12_results()
    shiny::validate(shiny::need(!is.null(res) && nrow(res) > 0, "No results."))
    write.csv(res, file, row.names = FALSE)
  }
)

# ── 2D KDE Density Co-localization Map (active file only) ───────────────
kde_coloc_results <- reactiveVal(NULL)

output$kde_plane_ui <- renderUI({
  nm  <- input$active_file %||% ""
  raw <- if (nzchar(nm)) smlm_file_store()[[nm]] else NULL
  has_z <- !is.null(raw) && "z" %in% names(raw) && any(is.finite(raw$z))
  choices <- c("XY plane" = "xy")
  if (has_z) choices <- c(choices, "XZ plane" = "xz", "YZ plane" = "yz")
  selectInput("kde_plane", "Projection plane", choices = choices, selected = "xy")
})

observeEvent(input$run_kde_coloc, {
  nm <- input$active_file %||% ""
  if (!nzchar(nm)) { showNotification("No active file selected.", type = "warning"); return() }
  store <- smlm_file_store()
  raw   <- store[[nm]]
  if (is.null(raw) || !nrow(raw)) { showNotification("Active file has no data.", type = "warning"); return() }

  unit          <- input$kde_unit %||% "tid"
  bw_raw        <- as.numeric(input$kde_bandwidth %||% 0)
  n_grid        <- as.integer(input$kde_n_grid %||% 100L)
  do_sweep      <- isTRUE(input$kde_bw_sweep)
  costes_n      <- max(0L, as.integer(input$kde_costes_nsim %||% 0L))
  needs_cluster <- unit == "cluster"

  df <- prep_tc_df(raw, smlm_settings_from_input(input), needs_cluster)
  if (is.null(df)) {
    showNotification("Could not assign two-color groups. Check 2-color settings.", type = "warning")
    kde_coloc_results(NULL); return()
  }

  pts1 <- get_tc_pts(df, 1L, unit)
  pts2 <- get_tc_pts(df, 2L, unit)
  if (!nrow(pts1) || !nrow(pts2)) {
    showNotification("One or both channels empty after filtering.", type = "warning")
    kde_coloc_results(NULL); return()
  }

  plane <- input$kde_plane %||% "xy"
  ax <- switch(plane,
    xz = c("x", "z"), yz = c("y", "z"), c("x", "y"))
  has_ax <- function(pts) all(ax %in% names(pts)) &&
    any(is.finite(pts[[ax[1]]])) && any(is.finite(pts[[ax[2]]]))
  if (!has_ax(pts1) || !has_ax(pts2)) {
    showNotification("Selected projection plane requires z coordinates which are absent.", type = "warning")
    kde_coloc_results(NULL); return()
  }

  bw  <- if (bw_raw > 0) bw_raw else NULL

  # Bandwidth sweep grid: log-spaced fractions of the field span (0.4%–20%),
  # so the bandwidth-dependence of Pearson r is visible.
  bw_sweep <- NULL
  if (do_sweep) {
    a1 <- c(pts1[[ax[1]]], pts2[[ax[1]]]); a2 <- c(pts1[[ax[2]]], pts2[[ax[2]]])
    span <- min(diff(range(a1, na.rm = TRUE)), diff(range(a2, na.rm = TRUE)))
    if (is.finite(span) && span > 0)
      bw_sweep <- span * exp(seq(log(0.004), log(0.2), length.out = 12))
  }

  # Single-file panel (active file only), so there is no set to share a gap-fill
  # radius with and the per-file plateau is the right choice here.
  kde_sup <- spatial_support_of(df, FALSE)
  res <- withProgress(message = "KDE co-localization...", value = 0.3,
    compute_kde_coloc(pts1[[ax[1]]], pts1[[ax[2]]],
                      pts2[[ax[1]]], pts2[[ax[2]]], bw, n_grid,
                      bw_sweep = bw_sweep, costes_nsim = costes_n,
                      null_mode = input$kde_null_mode %||% "shift",
                      support = kde_sup))
  if (is.null(res)) {
    showNotification("KDE computation failed.", type = "warning")
    kde_coloc_results(NULL); return()
  }
  res$dataset <- nm
  res$ax      <- ax
  # Whether the randomisation ran inside the sample footprint. With edge correction
  # on and the footprint switch off, the window passed in is still the bounding
  # box, so the switch (and an accepted footprint) is what decides the reading.
  res$footprint_null <- isTRUE(input$spatial_use_support) && !is.null(kde_sup) &&
                        !isTRUE(kde_sup$rejected_footprint)
  kde_coloc_results(res)
})

build_kde_coloc_plot_fig <- function() {
  r <- kde_coloc_results()
  shiny::validate(shiny::need(!is.null(r), "Run KDE first."))

  ax     <- r$ax %||% c("x", "y")
  xlabel <- paste0(ax[1], " (nm)")
  ylabel <- paste0(ax[2], " (nm)")

  style   <- input$kde_style    %||% "lines"
  cpair   <- input$kde_colorpair %||% "redblue"
  opacity <- as.numeric(input$kde_opacity %||% 0.75)

  color_pairs <- list(
    redblue      = list(ch1 = list(c(0,"white"), c(1,"#e41a1c")),
                        ch2 = list(c(0,"white"), c(1,"#377eb8"))),
    orangeteal   = list(ch1 = list(c(0,"white"), c(1,"#f46d43")),
                        ch2 = list(c(0,"white"), c(1,"#1a9988"))),
    magentagreen = list(ch1 = list(c(0,"white"), c(1,"#d01c8b")),
                        ch2 = list(c(0,"white"), c(1,"#4dac26"))),
    yellowpurple = list(ch1 = list(c(0,"white"), c(1,"#fc8d59")),
                        ch2 = list(c(0,"white"), c(1,"#542788")))
  )
  pal <- color_pairs[[cpair]] %||% color_pairs$redblue

  op <- if (style == "lines") 1.0 else opacity

  plotly::plot_ly() %>%
    plotly::add_contour(x = r$x, y = r$y, z = t(r$z1),
      name = "Ch1", opacity = op,
      colorscale = pal$ch1,
      showscale = FALSE, contours = list(coloring = style),
      line = list(width = if (style == "lines") 1.5 else 0.5)) %>%
    plotly::add_contour(x = r$x, y = r$y, z = t(r$z2),
      name = "Ch2", opacity = op,
      colorscale = pal$ch2,
      showscale = FALSE, contours = list(coloring = style),
      line = list(width = if (style == "lines") 1.5 else 0.5)) %>%
    plotly::layout(
      xaxis = list(title = xlabel),
      yaxis = list(title = ylabel, scaleanchor = "x", scaleratio = 1),
      title = sprintf("KDE - %s  [%s]  (Ch1 / Ch2)",
                      r$dataset, paste(ax, collapse = "/")),
      paper_bgcolor = "white", plot_bgcolor = "white"
  )
}
output$kde_coloc_plot <- renderPlotly(build_kde_coloc_plot_fig())

output$kde_pearson_display <- renderUI({
  r <- kde_coloc_results()
  if (is.null(r)) return(NULL)
  costes_line <- if (!is.null(r$costes_p)) {
    sig <- if (r$costes_p < 0.05) "significant co-localization" else "not significant"
    tagList(
      tags$p(sprintf("Spatial-randomization p = %.3f  (%s; %d randomizations, null r = %.3f ± %.3f)",
                     r$costes_p, sig, r$costes_nsim, r$costes_null_mean, r$costes_null_sd),
             style = "font-size:13px; margin:2px 0; color:#555;"),
      tags$p(sprintf("Null: %s", r$null_label %||% "-"),
             style = "font-size:12px; margin:0 0 2px; color:#888;"),
      # Measured false-positive rates (simulation, see compute_kde_coloc header):
      # signal spread over the whole window -> shift 3% / csr 3%; signal confined to
      # one cell inside a larger window -> shift 53% / csr 93%, vs the 5% a calibrated
      # test should give over the bounding box. Drawn inside the footprint instead,
      # the test is conservative, so the caveat depends on which window was used.
      if (identical(r$null_mode, "csr")) tags$p(
        HTML("&#9888; The uniform-CSR null tests against complete spatial randomness, so it largely reports &ldquo;both channels are non-uniform&rdquo; and rejects far too readily. Prefer the <b>random-shift</b> null."),
        style = "font-size:12px; margin:0; color:#b26a00;")
      else tags$p(
        HTML(if (isTRUE(r$footprint_null))
               "Null drawn inside the sample footprint, so the test is conservative: p &gt; 0.05 can be quoted as not significant. Take a positive claim from g&#8321;&#8322;(r)."
             else
               "&#9888; Null drawn over the bounding box, so for a structure filling part of the view this rejects too readily. Switch the footprint window on, or confirm with g&#8321;&#8322;(r)."),
        style = paste0("font-size:12px; margin:0; color:",
                       if (isTRUE(r$footprint_null)) "#555;" else "#b26a00;")))
  } else NULL
  tags$div(
    style = "margin-top:6px;",
    tags$b("Density-density Pearson r:"),
    tags$p(sprintf("%.3f  (bandwidth = %.1f nm)", r$pearson_r, r$bandwidth),
           style = "font-size:13px; margin:2px 0;"),
    costes_line,
    tags$p(HTML("Pearson r is bandwidth-dependent. For a bandwidth-free, scale-resolved colocalization measure use Pair Cross-Correlation g₁₂(r)."),
           style = "font-size:12px; margin:4px 0 0; color:#888; font-style:italic;")
  )
})

output$kde_bw_sweep_plot <- renderPlotly({
  r <- kde_coloc_results()
  shiny::validate(shiny::need(!is.null(r) && !is.null(r$bw_sweep),
                              "Enable 'Bandwidth sweep' and run KDE to see r(bandwidth)."))
  sw <- r$bw_sweep
  plotly::plot_ly() %>%
    plotly::add_trace(x = sw$bw, y = sw$pearson_r, type = "scatter", mode = "lines+markers",
      line = list(color = "#7B2D8B", width = 2), marker = list(size = 5),
      name = "Pearson r") %>%
    plotly::add_trace(x = r$bandwidth, y = r$pearson_r, type = "scatter", mode = "markers",
      marker = list(color = "#e41a1c", size = 10, symbol = "diamond"),
      name = "chosen bandwidth") %>%
    plotly::layout(
      xaxis = list(title = "KDE bandwidth (nm)"),
      yaxis = list(title = "Density-density Pearson r"),
      title = "Bandwidth dependence of Pearson r",
      hovermode = "x unified",
      paper_bgcolor = "white", plot_bgcolor = "white")
})

output$download_kde_coloc <- downloadHandler(
  filename = function() paste0("kde_coloc_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
  content  = function(file) {
    r <- kde_coloc_results()
    shiny::validate(shiny::need(!is.null(r), "No results."))
    tbl <- expand.grid(x = r$x, y = r$y) %>%
      dplyr::mutate(dataset     = r$dataset,
                    density_ch1 = as.vector(t(r$z1)),
                    density_ch2 = as.vector(t(r$z2)),
                    pearson_r   = r$pearson_r)
    write.csv(tbl, file, row.names = FALSE)
  }
)
