## ─────────────────────────────────────────────────────────────────────────────
## server/smlm/31_ripley_pcf.R  -  panels: Ripley L and pair correlation g(r)
## ─────────────────────────────────────────────────────────────────────────────
## Part of smlm_server(). Not a standalone script: its statements are evaluated in
## the server function's environment, in file-name order (server/parts.R).
##
## Run handlers with the CSR envelope, the correlation length, figures,
## summary and downloads.
##
## INTERFACE (generated: Rscript tools/server_map.R --headers; do not edit by hand)
##   DEFINES for other parts:
##     build_pcf_plot_fig, build_ripley_plot_fig, pcf_params, ripley_params
##   READS from other parts:
##     01_view_upload: smlm_file_store
##     17_plottable_data: plottable_data
##     20_dataset_groups: file_group_map, file_palette, series_palette
##     22_spatial_window: spatial_gate_close, spatial_report_windows,
##                        spatial_support_of
##     23_prepared_cache: prepared_df, resolve_unit
##     32_hotspots: spatial_outline_pts
## END INTERFACE
## ─────────────────────────────────────────────────────────────────────────────

## ── Ripley's L(r) ─────────────────────────────────────────────────────────
## ── Global (whole-curve) CSR test: band + caption ──────────────────────────


ripley_results     <- reactiveVal(NULL)
ripley_params      <- reactiveVal(NULL)   # compute settings, frozen at Run
ripley_env_results <- reactiveVal(NULL)  # Monte-Carlo CSR envelope (per dataset)

observeEvent(input$run_ripley, {
  unit      <- resolve_unit(input$ripley_unit, "cluster", "tid", "Ripley's L")
  use_3d    <- isTRUE(input$ripley_use_3d)
  r_max_in  <- as.numeric(input$ripley_r_max  %||% 0)
  r_steps   <- max(10L, as.integer(input$ripley_r_steps %||% 60L))
  all_files <- isTRUE(input$ripley_all_files)
  want_env  <- isTRUE(input$ripley_envelope)
  nsim      <- max(19L, min(199L, as.integer(input$ripley_nsim %||% 99L)))
  env_file_cap <- 6L  # cap envelope work in all-files mode (nsim × files is O(n²)-heavy)

  needs_cluster <- unit == "cluster"

  compute_ripley_for_df <- function(df, tag, r_seq_use, do_env, close_nm = NULL,
                                    win_df = NULL, win_tag = NULL) {
    pts <- unit_centroids(df, unit = unit, use_3d = use_3d)
    if (is.null(pts) || nrow(pts) < 2L) return(NULL)
    # THE WINDOW COMES FROM win_df, NOT FROM pts. Normally they are the same
    # frame. They differ when this file has been split into channels: there
    # win_df is the pooled file, so both channels are measured in one window
    # (see split_df_by_channel). Deriving it from the channel's own points
    # would shrink the window around the sparser channel and deflate its curve.
    wdf  <- if (is.data.frame(win_df) && nrow(win_df)) win_df else df
    wpts <- if (identical(wdf, df)) pts else unit_centroids(wdf, unit = unit, use_3d = use_3d)
    if (is.null(wpts) || !nrow(wpts)) wpts <- pts
    axes <- if (use_3d) c("x", "y", "z") else c("x", "y")
    win  <- apply(as.matrix(wpts[, axes, drop = FALSE]), 2, range)
    # Support window (opt-in), keyed on win_tag so channels of one file share a
    # single cached support. Same shared choice of outline points as the
    # hotspot panel, so one cell gets one area whichever panel measured it.
    sup  <- spatial_support_of(spatial_outline_pts(wdf, wpts), use_3d,
                               tag = win_tag %||% tag, close_nm = close_nm)
    spatial_curve_for_points(pts, win, estimator = "L", r_seq = r_seq_use,
                             use_3d = use_3d, support = sup,
                             edge_correct = isTRUE(input$spatial_edge_correct),
                             envelope = do_env, nsim = nsim, progress = env_prog, tag = tag)
  }

  get_df_list <- function() {
    if (all_files) {
      store <- smlm_file_store()
      shiny::validate(shiny::need(length(store) > 0, "No files loaded."))
      # Prepare each FILE once, then split it into channels if asked. Splitting
      # after preparation (rather than handing the pipeline pre-split raw data,
      # as this used to) is what lets both channels keep the pooled file as
      # their shared observation window - and it means the split now covers a
      # 2-colour DCR split as well as a multiplex build.
      items <- lapply(names(store), function(nm) {
        # Cached chain. NOTE this also brings run_ripley in line with every other
        # panel by applying qc_exclude_from, which it previously skipped (a no-op
        # unless 'Exclude flagged traces' is on).
        df <- prepared_df(nm, needs_cluster = needs_cluster)
        if (is.null(df)) df <- store[[nm]][0, , drop = FALSE]
        list(df = df, tag = nm)
      })
      if (isTRUE(input$ripley_by_channel))
        items <- unlist(lapply(items, function(it) split_df_by_channel(it$df, it$tag)),
                        recursive = FALSE)
      else
        items <- lapply(items, function(it) c(it, list(win_df = it$df, win_tag = it$tag)))
      items
    } else {
      df <- plottable_data()
      shiny::validate(shiny::need(!is.null(df) && nrow(df) > 0, "Press Plot first."))
      nm <- input$active_file %||% "Active file"
      if (isTRUE(input$ripley_by_channel)) split_df_by_channel(df, nm)
      else list(list(df = df, tag = nm, win_df = df, win_tag = nm))
    }
  }

  df_list <- get_df_list()

  # One gap-fill radius for the whole set (see spatial_resolve_close()): a
  # per-file radius makes the windows - and therefore K, which scales with the
  # window area - incomparable between files. Never blocks.
  # The sweep runs on the same points the window will be drawn from.
  cc <- spatial_gate_close(lapply(df_list, function(it) {
    u <- tryCatch(unit_centroids(it$df, unit = unit, use_3d = use_3d), error = function(e) NULL)
    list(df = spatial_outline_pts(it$df, u), tag = it$tag)
  }), use_3d)
  common_close <- cc$close_nm

  if (want_env && length(df_list) > env_file_cap) {
    showNotification(sprintf(
      "CSR envelope computed for the first %d of %d files only (cost is nsim × files).",
      env_file_cap, length(df_list)), type = "message")
  }

  # The envelope dominates the runtime, so the bar advances per SIMULATION.
  n_env_files <- if (want_env) min(length(df_list), env_file_cap) else 0L
  env_prog <- if (want_env)
    function(frac) incProgress(frac / max(1L, n_env_files)) else NULL

  run_all <- function() {
    lapply(seq_along(df_list), function(i) {
      item      <- df_list[[i]]
      r_seq_use <- if (r_max_in > 0) seq(0, r_max_in, length.out = r_steps + 1L)[-1] else NULL
      do_env    <- want_env && i <= env_file_cap
      compute_ripley_for_df(item$df, item$tag, r_seq_use, do_env, common_close,
                            win_df = item$win_df, win_tag = item$win_tag)
    })
  }
  all_out <- if (want_env) {
    withProgress(message = sprintf("Ripley L: %d CSR simulations x %d file(s)...",
                                   nsim, n_env_files), value = 0, run_all())
  } else {
    # Without an envelope there is no per-simulation tick, but building the
    # windows still takes seconds - show a bar rather than appear frozen.
    withProgress(message = sprintf("Ripley L: %d file(s)...", length(df_list)),
                 value = 0.5, run_all())
  }
  spatial_report_windows(vapply(df_list, function(it) as.character(it$tag), character(1)))

  all_out    <- Filter(Negate(is.null), all_out)
  all_curves <- do.call(rbind, lapply(all_out, `[[`, "curve"))
  shiny::validate(shiny::need(!is.null(all_curves) && nrow(all_curves) > 0, "Ripley computation returned no data."))
  ripley_results(all_curves)
  ripley_params(list(unit = input$ripley_unit, `r max` = input$ripley_r_max,
                     `r steps` = input$ripley_r_steps, `3-D` = use_3d,
                     `all files` = all_files, envelope = want_env,
                     simulations = input$ripley_nsim,
                     `edge correction` = isTRUE(input$spatial_edge_correct),
                     `support window` = isTRUE(input$spatial_use_support)))

  envs <- Filter(Negate(is.null), lapply(all_out, `[[`, "env"))
  ripley_env_results(if (length(envs)) do.call(rbind, envs) else NULL)

})

build_ripley_plot_fig <- function() {
  res <- ripley_results()
  shiny::validate(shiny::need(!is.null(res) && nrow(res) > 0, "Run Ripley's L first."))

  do_smooth <- isTRUE(input$ripley_smooth)
  span_val  <- as.numeric(input$ripley_smooth_span %||% 0.2)
  smooth_vec <- function(r, y) {
    if (!do_smooth || length(y) < 5L) return(y)
    tryCatch(stats::predict(stats::loess(y ~ r, data = data.frame(r = r, y = y), span = span_val)),
             error = function(e) y)
  }
  hex2rgba <- function(hex, a) {
    cc <- grDevices::col2rgb(hex); sprintf("rgba(%d,%d,%d,%.2f)", cc[1], cc[2], cc[3], a)
  }
  grouped <- (input$analysis_by %||% "file") == "group"
  env_res <- ripley_env_results()
  # Declared up front because the axis code below reads them in both branches.
  rmin_c <- NA_real_; rmax_c <- NA_real_
  show_files <- isTRUE(input$ripley_show_files)
  fig <- plotly::plot_ly()

  if (!grouped) {
    datasets <- unique(res$dataset)
    pal <- series_palette(datasets)
    for (nm in datasets) {
      sub <- res[res$dataset == nm, ]
      # CSR envelope band (drawn first so the observed curve sits on top)
      if (!is.null(env_res) && nm %in% env_res$dataset) {
        # The BAND MUST GET THE SAME SMOOTHING AS THE CURVE. Comparing a smoothed
        # observed curve against a raw envelope compares two different statistics:
        # smoothing damps the curve's excursions but left the band at its full
        # unsmoothed width, so the comparison silently became conservative (and the
        # band looked jagged next to a smooth curve). Smoothing the quantiles is a
        # close approximation to taking quantiles of smoothed simulations - the two
        # differ only in how a monotone smoother interacts with the quantile, which
        # is negligible here and keeps the band consistent with whatever smoothing
        # is displayed (smoothing is a live display control, the envelope is not).
        e <- env_res[env_res$dataset == nm, ]
        e <- e[order(e$r), ]
        e_hi <- smooth_vec(e$r, e$hi); e_lo <- smooth_vec(e$r, e$lo)
        e_nm <- smooth_vec(e$r, e$null_mean)
        fig <- fig %>%
          plotly::add_trace(x = e$r, y = e_hi, type = "scatter", mode = "lines",
            line = list(width = 0), legendgroup = nm, showlegend = FALSE, hoverinfo = "skip") %>%
          plotly::add_trace(x = e$r, y = e_lo, type = "scatter", mode = "lines", fill = "tonexty",
            fillcolor = hex2rgba(unname(pal[nm]), 0.15), line = list(width = 0),
            legendgroup = nm, name = paste0(nm, " CSR 95% envelope"),
            showlegend = TRUE, hoverinfo = "skip") %>%
          plotly::add_trace(x = e$r, y = e_nm, type = "scatter", mode = "lines",
            line = list(color = unname(pal[nm]), width = 1, dash = "dot"),
            legendgroup = nm, showlegend = FALSE, hoverinfo = "skip")
        fig <- .add_global_band(fig, e, unname(pal[nm]), nm, smooth_vec)
      }
      fig <- fig %>% plotly::add_trace(
        x = sub$r, y = smooth_vec(sub$r, sub$Lr_minus_r),
        type = "scatter", mode = "lines", name = nm, legendgroup = nm,
        line = list(color = unname(pal[nm]), width = 2), hoverinfo = "skip")
    }
  } else {
    gm <- file_group_map()
    shiny::validate(shiny::need(length(gm) > 0, "No files assigned to groups (group mode)."))
    res <- res[res$dataset %in% names(gm), , drop = FALSE]
    shiny::validate(shiny::need(nrow(res) > 0, "No group files in the Ripley results."))
    res$group <- unname(gm[res$dataset])
    pal  <- series_palette(unique(res$group))
    fpal <- file_palette(unique(res$dataset))
    band <- input$ripley_group_band %||% "iqr"

    # common r-grid up to the smallest per-file r max (avoid extrapolation),
    # then average each group's per-file curves + between-file spread band.
    file_maxr <- tapply(res$r, res$dataset, function(v) max(v, na.rm = TRUE))
    rmax_c <- suppressWarnings(min(file_maxr, na.rm = TRUE))
    # Common grid spans only the range EVERY file actually measured. It used to
    # start at r = 0 while the stored curves start at r_max/steps, so approx(rule=2)
    # flat-extrapolated the first real value back to r = 0 - a false non-zero point
    # at a radius where L(r)-r is 0 by definition, injected into every group curve
    # and its band.
    file_minr <- tapply(res$r, res$dataset, function(v) min(v, na.rm = TRUE))
    rmin_c <- suppressWarnings(max(file_minr, na.rm = TRUE))
    if (!is.finite(rmin_c) || rmin_c < 0) rmin_c <- 0
    # Grid density from the STORED curve, not the live control: r steps is a
    # compute parameter and only takes effect on Run.
    steps  <- max(10L, length(unique(res$r)))
    shiny::validate(shiny::need(is.finite(rmax_c) && rmax_c > 0, "Could not form a common r-grid."))
    rr <- seq(rmin_c, rmax_c, length.out = steps)

    for (g in unique(res$group)) {
      files <- unique(res$dataset[res$group == g])
      mat <- vapply(files, function(f) {
        d <- res[res$dataset == f, ]
        stats::approx(d$r, d$Lr_minus_r, xout = rr, rule = 2)$y
      }, numeric(length(rr)))
      if (is.null(dim(mat))) mat <- matrix(mat, ncol = 1)
      meanc <- smooth_vec(rr, rowMeans(mat, na.rm = TRUE))
      col   <- unname(pal[g]); nfl <- length(files)

      # ── BAND MUST MATCH THE DRAWN MEAN ──────────────────────────────────
      # Three faults were fixed here together:
      #  1. the mean was smoothed and the band was not, so a smooth line sat in a
      #     jagged band (measured: band edge 16-18x the roughness of the mean).
      #     Same argument as the CSR envelope above - a band and a curve that got
      #     different treatment are not comparable.
      #  2. the SD band was built around the RAW rowMeans while the SMOOTHED mean
      #     was drawn, so the line was not centred in its own band (up to 22 nm
      #     off). It is now centred on the curve actually shown.
      #  3. quartiles need >= 4 values. On 2 files the IQR band is the middle half
      #     of a 2-point range and contains NEITHER file (measured 0% coverage),
      #     which draws a reassuringly tight band over maximal disagreement. Under
      #     4 files the full min-max range is drawn instead, and the legend says so.
      # ── CSR ENVELOPE IN GROUP MODE ──────────────────────────────────────
      # This used to be omitted entirely: envelopes were computed per file but
      # the grouped branch never drew them, so ticking "Show CSR envelope" did
      # nothing as soon as you switched to groups - and the between-file spread
      # band was easy to mistake for it.
      #
      # The group curve is the mean of its files' curves, so its CSR reference
      # is the mean of those files' envelopes, interpolated onto the same grid
      # and smoothed the same way. Drawn as a dashed outline rather than a fill,
      # so it cannot be confused with the solid spread band: they answer
      # different questions - the spread band is "how much do my cells differ",
      # the envelope is "what would randomness look like".
      env_n <- 0L
      if (!is.null(env_res)) {
        ef <- intersect(files, unique(env_res$dataset))
        if (length(ef)) {
          elo <- vapply(ef, function(f) {
            e <- env_res[env_res$dataset == f, ]; e <- e[order(e$r), ]
            stats::approx(e$r, e$lo, xout = rr, rule = 2)$y }, numeric(length(rr)))
          ehi <- vapply(ef, function(f) {
            e <- env_res[env_res$dataset == f, ]; e <- e[order(e$r), ]
            stats::approx(e$r, e$hi, xout = rr, rule = 2)$y }, numeric(length(rr)))
          if (is.null(dim(elo))) { elo <- matrix(elo, ncol = 1); ehi <- matrix(ehi, ncol = 1) }
          env_n <- length(ef)
          e_lo <- smooth_vec(rr, rowMeans(elo, na.rm = TRUE))
          e_hi <- smooth_vec(rr, rowMeans(ehi, na.rm = TRUE))
          # TWO lines per group - the upper and lower bound of the same band -
          # so N groups draw 2N dotted lines. They are given their own legend
          # entry and their own legendgroup, so it is obvious what they are and
          # they can be switched off independently of the group curve.
          fig <- fig %>% plotly::add_trace(
            x = rr, y = e_hi, type = "scatter", mode = "lines",
            line = list(color = col, width = 1, dash = "dot"),
            name = sprintf("%s – CSR 95%% band", g),
            legendgroup = paste0(g, "__csr"), showlegend = TRUE, hoverinfo = "skip")
          fig <- fig %>% plotly::add_trace(
            x = rr, y = e_lo, type = "scatter", mode = "lines",
            line = list(color = col, width = 1, dash = "dot"),
            legendgroup = paste0(g, "__csr"), showlegend = FALSE, hoverinfo = "skip")
        }
      }

      band_lbl <- ""
      if (band != "none" && nfl >= 2) {
        if (band == "iqr" && nfl >= 4L) {
          lo <- smooth_vec(rr, apply(mat, 1, stats::quantile, 0.25, na.rm = TRUE))
          hi <- smooth_vec(rr, apply(mat, 1, stats::quantile, 0.75, na.rm = TRUE))
          band_lbl <- ", IQR"
        } else if (band == "iqr") {
          lo <- smooth_vec(rr, apply(mat, 1, min, na.rm = TRUE))
          hi <- smooth_vec(rr, apply(mat, 1, max, na.rm = TRUE))
          band_lbl <- ", range"
        } else {
          sdv <- smooth_vec(rr, apply(mat, 1, stats::sd, na.rm = TRUE))
          lo <- meanc - sdv; hi <- meanc + sdv
          band_lbl <- ", ±1 SD"
        }
        fig <- fig %>%
          plotly::add_trace(x = rr, y = hi, type = "scatter", mode = "lines",
            line = list(width = 0), legendgroup = g, showlegend = FALSE, hoverinfo = "skip") %>%
          plotly::add_trace(x = rr, y = lo, type = "scatter", mode = "lines", fill = "tonexty",
            fillcolor = hex2rgba(col, 0.18), line = list(width = 0),
            legendgroup = g, showlegend = FALSE, hoverinfo = "skip")
      }
      # Individual files as shades of the group colour, drawn thin and behind the
      # group mean, so a group difference can be judged against its replicates.
      if (show_files && nfl >= 2) {
        for (f in files) {
          d <- res[res$dataset == f, ]; d <- d[order(d$r), ]
          fc <- unname(fpal[f]); if (is.na(fc)) fc <- col
          fig <- fig %>% plotly::add_trace(
            x = d$r, y = d$Lr_minus_r, type = "scatter", mode = "lines",
            name = f, legendgroup = g, showlegend = FALSE,
            line = list(color = fc, width = 1), opacity = 0.75,
            hoverinfo = "skip")
        }
      }
      # Say how many files the envelope averages. With more files than the
      # envelope cap it covers only some of them, and a reader has to know that
      # before reading the dotted lines as the group's null.
      env_lbl <- if (env_n > 0L && env_n < nfl) sprintf(", CSR n=%d", env_n)
                 else if (env_n > 0L) ", CSR" else ""
      fig <- fig %>% plotly::add_trace(x = rr, y = meanc, type = "scatter", mode = "lines",
        name = sprintf("%s (n=%d%s%s)", g, nfl, band_lbl, env_lbl), legendgroup = g,
        line = list(color = col, width = 2.5), hoverinfo = "skip")
    }
  }

  # ── The x extent that is actually DRAWN ───────────────────────────────
  # In group mode the mean and its band live on a common grid capped at the
  # SMALLEST per-file r max, so that no file is extrapolated. The reference line
  # used to be drawn over range(res$r) instead - the range across ALL files - so
  # one file with a larger r max stretched it far past every curve, plotly
  # autoranged to the line, and the plot ended in a wide empty strip. The
  # reference line now spans exactly what is plotted, and nothing else reaches
  # beyond it unless the individual per-file curves are switched on, which draw
  # at their own full range and legitimately extend the axis.
  x_drawn <- if (grouped) {
    c(rmin_c, if (isTRUE(show_files)) max(res$r, na.rm = TRUE) else rmax_c)
  } else range(res$r, na.rm = TRUE)

  csr_label <- if (!is.null(env_res) && !grouped)
    "theoretical CSR (uncorrected - see envelope)" else "CSR (L−r = 0)"
  fig %>%
    plotly::add_trace(
      x = x_drawn, y = c(0, 0), type = "scatter", mode = "lines",
      name = csr_label,
      line = list(color = "black", width = 1, dash = "dash"), showlegend = TRUE,
      hoverinfo = "skip") %>%
    plotly::layout(
      xaxis  = list(title = .rmax_axis_title(grouped, res, rmax_c, show_files)),
      yaxis  = list(title = "L(r) − r (nm)"),
      annotations = .global_p_annot(env_res),
      # Hover deliberately OFF on this plot: the curve is read as a shape against
      # the envelope, and a unified tooltip listing every series' L−r at each r
      # was noise that also covered the plot.
      hovermode = FALSE,
      legend = list(groupclick = "togglegroup"),
      paper_bgcolor = "white", plot_bgcolor  = "white")
}
output$ripley_plot <- renderPlotly(build_ripley_plot_fig())

output$download_ripley <- downloadHandler(
  filename = function() paste0("ripley_L_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
  content = function(file) {
    res <- ripley_results()
    shiny::validate(shiny::need(!is.null(res) && nrow(res) > 0, "No results."))
    write.csv(res, file, row.names = FALSE)
  }
)

## ── Univariate pair-correlation g(r) ──────────────────────────────────────
## GOAL: the scale-LOCAL clustering statistic, and the recommended endpoint for
## comparing conditions whose labelling density differs.
##
## WHY IT IS THE RIGHT ENDPOINT HERE (see also the LABELLING-DENSITY block):
##   * g is normalised by the intensity SQUARED, so it is invariant to how
##     densely the sample happens to be labelled - verified invariant to within
##     1-2% under thinning down to 10% of the units.
##   * being a derivative rather than a cumulative statistic, an aggregate of
##     radius R shows up as a PEAK near R instead of the broad rising ramp L(r)
##     gives, and a short-r blinking/label artefact stays in the first few bins
##     instead of contaminating every larger r.
## Runs on de-blinked units, same as Ripley's L.
pcf_results     <- reactiveVal(NULL)
pcf_params      <- reactiveVal(NULL)   # compute settings, frozen at Run
pcf_env_results <- reactiveVal(NULL)

observeEvent(input$run_pcf, {
  unit      <- resolve_unit(input$pcf_unit, "cluster", "tid", "g(r)")
  use_3d    <- isTRUE(input$pcf_use_3d)
  r_max_in  <- as.numeric(input$pcf_r_max %||% 0)
  r_steps   <- max(10L, as.integer(input$pcf_r_steps %||% 50L))
  all_files <- isTRUE(input$pcf_all_files)
  want_env  <- isTRUE(input$pcf_envelope)
  nsim      <- max(19L, min(199L, as.integer(input$pcf_nsim %||% 99L)))
  env_file_cap <- 6L

  needs_cluster <- unit == "cluster"

  compute_for_df <- function(df, tag, r_seq_use, do_env, close_nm = NULL,
                             win_df = NULL, win_tag = NULL) {
    pts <- unit_centroids(df, unit = unit, use_3d = use_3d)
    if (is.null(pts) || nrow(pts) < 3L) return(NULL)
    # Window from win_df, not from pts - see compute_ripley_for_df(). When this
    # file is split into channels, win_df is the pooled file so both channels
    # are measured in one window; g(r) divides by density, so a per-channel
    # window would invent a difference out of the count imbalance alone.
    wdf  <- if (is.data.frame(win_df) && nrow(win_df)) win_df else df
    wpts <- if (identical(wdf, df)) pts else unit_centroids(wdf, unit = unit, use_3d = use_3d)
    if (is.null(wpts) || !nrow(wpts)) wpts <- pts
    # Freeze the window so the curve and its envelope share one window and the
    # (uncorrected) edge bias cancels between them.
    axes <- if (use_3d) c("x", "y", "z") else c("x", "y")
    win  <- apply(as.matrix(wpts[, axes, drop = FALSE]), 2, range)
    # Same shared choice of outline points as the hotspot panel, so one cell
    # gets one area whichever panel measured it.
    sup  <- spatial_support_of(spatial_outline_pts(wdf, wpts), use_3d,
                               tag = win_tag %||% tag, close_nm = close_nm)
    spatial_curve_for_points(pts, win, estimator = "g", r_seq = r_seq_use,
                             use_3d = use_3d, support = sup,
                             edge_correct = isTRUE(input$spatial_edge_correct),
                             envelope = do_env, nsim = nsim, progress = env_prog, tag = tag)
  }

  get_df_list <- function() {
    if (all_files) {
      store <- smlm_file_store()
      shiny::validate(shiny::need(length(store) > 0, "No files loaded."))
      # Prepare each FILE once, then split into channels - so both channels
      # share the pooled file as their observation window.
      items <- lapply(names(store), function(nm) {
        df <- prepared_df(nm, needs_cluster = needs_cluster)
        if (is.null(df)) df <- store[[nm]][0, , drop = FALSE]
        list(df = df, tag = nm)
      })
      if (isTRUE(input$pcf_by_channel))
        unlist(lapply(items, function(it) split_df_by_channel(it$df, it$tag)),
               recursive = FALSE)
      else
        lapply(items, function(it) c(it, list(win_df = it$df, win_tag = it$tag)))
    } else {
      df <- plottable_data()
      shiny::validate(shiny::need(!is.null(df) && nrow(df) > 0, "Press Plot first."))
      nm <- input$active_file %||% "Active file"
      if (isTRUE(input$pcf_by_channel)) split_df_by_channel(df, nm)
      else list(list(df = df, tag = nm, win_df = df, win_tag = nm))
    }
  }

  df_list <- get_df_list()

  # One gap-fill radius for the whole set (see spatial_resolve_close()), swept on
  # the same points the window will be drawn from.
  cc <- spatial_gate_close(lapply(df_list, function(it) {
    u <- tryCatch(unit_centroids(it$df, unit = unit, use_3d = use_3d),
                  error = function(e) NULL)
    list(df = spatial_outline_pts(it$df, u), tag = it$tag)
  }), use_3d)
  common_close <- cc$close_nm

  if (want_env && length(df_list) > env_file_cap)
    showNotification(sprintf(
      "CSR envelope computed for the first %d of %d files only (cost is simulations x files).",
      env_file_cap, length(df_list)), type = "message")

  n_env_files <- if (want_env) min(length(df_list), env_file_cap) else 0L
  env_prog <- if (want_env)
    function(frac) incProgress(frac / max(1L, n_env_files)) else NULL

  run_all <- function() lapply(seq_along(df_list), function(i) {
    r_seq_use <- if (r_max_in > 0) seq(0, r_max_in, length.out = r_steps) else NULL
    compute_for_df(df_list[[i]]$df, df_list[[i]]$tag, r_seq_use,
                   want_env && i <= env_file_cap, common_close,
                   win_df = df_list[[i]]$win_df, win_tag = df_list[[i]]$win_tag)
  })
  all_out <- if (want_env) withProgress(
               message = sprintf("g(r): %d CSR simulations x %d file(s)...", nsim, n_env_files),
               value = 0, run_all())
             else withProgress(message = sprintf("g(r): %d file(s)...", length(df_list)),
                               value = 0.5, run_all())
  spatial_report_windows(vapply(df_list, function(it) as.character(it$tag), character(1)))

  all_out <- Filter(Negate(is.null), all_out)
  curves  <- do.call(rbind, lapply(all_out, `[[`, "curve"))
  shiny::validate(shiny::need(!is.null(curves) && nrow(curves) > 0,
                             "g(r) returned no data (need >=3 units per file)."))
  pcf_results(curves)
  pcf_params(list(unit = input$pcf_unit, `r max` = input$pcf_r_max,
                  `r steps` = input$pcf_r_steps, `3-D` = use_3d,
                  `all files` = all_files, envelope = want_env,
                  simulations = input$pcf_nsim,
                  `edge correction` = isTRUE(input$spatial_edge_correct),
                  `support window` = isTRUE(input$spatial_use_support)))
  envs <- Filter(Negate(is.null), lapply(all_out, `[[`, "env"))
  pcf_env_results(if (length(envs)) do.call(rbind, envs) else NULL)
})

# ── Correlation length xi, one row per curve ──────────────────────────────
# Fits each curve the plot actually draws: per FILE in file mode, and per GROUP
# in group mode - where the group curve is the mean of its files' curves on a
# common r-grid, exactly as plotted, so the number matches the line on screen.
#
# Per-file rows are kept in group mode too. A single xi per group hides whether
# the group is tight or whether one cell is carrying it, and that spread is the
# thing a test between conditions has to work with.
pcf_xi_results <- reactive({
  req(isTRUE(input$pcf_fit_xi))
  res <- pcf_results()
  if (is.null(res) || !nrow(res)) return(NULL)
  r_min <- suppressWarnings(as.numeric(input$pcf_fit_rmin))
  if (!length(r_min) || !is.finite(r_min) || r_min < 0) r_min <- 0

  per_file <- do.call(rbind, lapply(split(res, res$dataset), function(d) {
    d <- d[order(d$r), ]
    f <- fit_correlation_length(d$r, d$g, r_min = r_min)
    cbind(data.frame(curve = d$dataset[1], level = "file",
                     stringsAsFactors = FALSE), f)
  }))

  out <- per_file
  if ((input$analysis_by %||% "file") == "group") {
    gm <- file_group_map()
    if (length(gm)) {
      res$group <- unname(gm[res$dataset])
      res <- res[!is.na(res$group), , drop = FALSE]
      if (nrow(res)) {
        # Same common-grid averaging the plot uses - never pool coordinates.
        per_grp <- do.call(rbind, lapply(split(res, res$group), function(d) {
          files <- unique(d$dataset)
          rmin_c <- max(tapply(d$r, d$dataset, min), na.rm = TRUE)
          rmax_c <- min(tapply(d$r, d$dataset, max), na.rm = TRUE)
          if (!is.finite(rmin_c) || !is.finite(rmax_c) || rmax_c <= rmin_c) return(NULL)
          rr <- seq(rmin_c, rmax_c, length.out = max(10L, length(unique(d$r))))
          mat <- vapply(files, function(f) {
            dd <- d[d$dataset == f, ]; dd <- dd[order(dd$r), ]
            stats::approx(dd$r, dd$g, xout = rr, rule = 2)$y
          }, numeric(length(rr)))
          if (is.null(dim(mat))) mat <- matrix(mat, ncol = 1)
          f <- fit_correlation_length(rr, rowMeans(mat, na.rm = TRUE), r_min = r_min)
          cbind(data.frame(curve = d$group[1], level = "group",
                           stringsAsFactors = FALSE), f)
        }))
        if (!is.null(per_grp) && nrow(per_grp))
          out <- rbind(per_grp, per_file[per_file$curve %in% names(gm), , drop = FALSE])
      }
    }
  }
  out
})

output$pcf_xi_table <- renderTable({
  x <- pcf_xi_results()
  shiny::validate(shiny::need(!is.null(x) && nrow(x), "Run g(r) first."))
  data.frame(
    Curve   = x$curve,
    Level   = x$level,
    `xi (nm)` = ifelse(is.na(x$xi), "-", sprintf("%.0f", x$xi)),
    `± SE`  = ifelse(is.na(x$xi_se), "-", sprintf("%.0f", x$xi_se)),
    A       = ifelse(is.na(x$amplitude), "-", sprintf("%.2f", x$amplitude)),
    `R²`    = ifelse(is.na(x$r2), "-", sprintf("%.3f", x$r2)),
    `fit range (nm)` = sprintf("%.0f-%.0f", x$r_min, x$r_max),
    n       = x$n_used,
    check.names = FALSE, stringsAsFactors = FALSE)
}, striped = TRUE, spacing = "xs", width = "100%")

output$pcf_xi_note <- renderUI({
  x <- pcf_xi_results(); if (is.null(x) || !nrow(x)) return(NULL)
  msgs <- list()
  if (any(is.na(x$xi)))
    msgs <- c(msgs, sprintf("%d curve(s) show no decaying excess over the fit range - no correlation length is reported for them, rather than a number fitted to noise.",
                            sum(is.na(x$xi))))
  if (any(isTRUE(x$below_floor)) || any(x$below_floor %in% TRUE))
    msgs <- c(msgs, "At least one ξ is SMALLER than the fit's starting radius, so the fit began past the decay it was measuring - lower 'Fit from r' for those curves before quoting the value.")
  if (any(grepl("log-linear", x$method %||% "")))
    msgs <- c(msgs, "One or more fits fell back to the log-linear estimator (nls did not converge); that estimator over-estimates ξ on noisy tails.")
  if (!length(msgs)) return(NULL)
  tags$div(style = "font-size:0.82em; color:#7a5c00; margin-top:4px;",
           lapply(msgs, function(m) tags$div(m)))
})

build_pcf_plot_fig <- function() {
  res <- pcf_results()
  shiny::validate(shiny::need(!is.null(res) && nrow(res) > 0, "Run g(r) first."))

  do_smooth <- isTRUE(input$pcf_smooth)
  span_val  <- as.numeric(input$pcf_smooth_span %||% 0.25)
  smooth_vec <- function(r, y) {
    ok <- is.finite(y)
    if (!do_smooth || sum(ok) < 5L) return(y)
    tryCatch(stats::predict(stats::loess(y ~ r, data = data.frame(r = r, y = y),
                                         span = span_val, na.action = stats::na.exclude),
                            newdata = data.frame(r = r)),
             error = function(e) y)
  }
  hex2rgba <- function(hex, a) {
    cc <- grDevices::col2rgb(hex); sprintf("rgba(%d,%d,%d,%.2f)", cc[1], cc[2], cc[3], a)
  }
  grouped <- (input$analysis_by %||% "file") == "group"
  env_res <- pcf_env_results()
  # Declared up front because the axis code below reads them in both branches.
  rmin_c <- NA_real_; rmax_c <- NA_real_
  show_files <- isTRUE(input$pcf_show_files)
  fig <- plotly::plot_ly()

  if (!grouped) {
    datasets <- unique(res$dataset)
    pal <- series_palette(datasets)
    for (nm in datasets) {
      sub <- res[res$dataset == nm, ]; sub <- sub[order(sub$r), ]
      if (!is.null(env_res) && nm %in% env_res$dataset) {
        # The BAND MUST GET THE SAME SMOOTHING AS THE CURVE. Comparing a smoothed
        # observed curve against a raw envelope compares two different statistics:
        # smoothing damps the curve's excursions but left the band at its full
        # unsmoothed width, so the comparison silently became conservative (and the
        # band looked jagged next to a smooth curve). Smoothing the quantiles is a
        # close approximation to taking quantiles of smoothed simulations - the two
        # differ only in how a monotone smoother interacts with the quantile, which
        # is negligible here and keeps the band consistent with whatever smoothing
        # is displayed (smoothing is a live display control, the envelope is not).
        e <- env_res[env_res$dataset == nm, ]; e <- e[order(e$r), ]
        e_hi <- smooth_vec(e$r, e$hi); e_lo <- smooth_vec(e$r, e$lo)
        fig <- fig %>%
          plotly::add_trace(x = e$r, y = e_hi, type = "scatter", mode = "lines",
            line = list(width = 0), legendgroup = nm, showlegend = FALSE, hoverinfo = "skip") %>%
          plotly::add_trace(x = e$r, y = e_lo, type = "scatter", mode = "lines", fill = "tonexty",
            fillcolor = hex2rgba(unname(pal[nm]), 0.15), line = list(width = 0),
            legendgroup = nm, name = paste0(nm, " CSR 95% envelope"), hoverinfo = "skip")
        fig <- .add_global_band(fig, e, unname(pal[nm]), nm, smooth_vec)
      }
      fig <- fig %>% plotly::add_trace(
        x = sub$r, y = smooth_vec(sub$r, sub$g), type = "scatter", mode = "lines",
        name = sprintf("%s (n=%d)", nm, sub$n_units[1]), legendgroup = nm,
        line = list(color = unname(pal[nm]), width = 2), hoverinfo = "skip")
    }
  } else {
    gm <- file_group_map()
    shiny::validate(shiny::need(length(gm) > 0, "No files assigned to groups (group mode)."))
    res <- res[res$dataset %in% names(gm), , drop = FALSE]
    shiny::validate(shiny::need(nrow(res) > 0, "No group files in the g(r) results."))
    res$group <- unname(gm[res$dataset])
    pal  <- series_palette(unique(res$group))
    fpal <- file_palette(unique(res$dataset))
    band <- input$pcf_group_band %||% "iqr"

    # Same contract as the Ripley group mode: interpolate each FILE's curve onto a
    # common r-grid and average within the group. Coordinates are never pooled -
    # pooling would mix windows and densities and invent structure.
    file_maxr <- tapply(res$r, res$dataset, function(v) max(v, na.rm = TRUE))
    rmax_c <- suppressWarnings(min(file_maxr, na.rm = TRUE))
    # Common grid spans only the range EVERY file actually measured. It used to
    # start at r = 0 while the stored curves start at r_max/steps, so approx(rule=2)
    # flat-extrapolated the first real value back to r = 0 - a false non-zero point
    # at a radius where L(r)-r is 0 by definition, injected into every group curve
    # and its band.
    file_minr <- tapply(res$r, res$dataset, function(v) min(v, na.rm = TRUE))
    rmin_c <- suppressWarnings(max(file_minr, na.rm = TRUE))
    if (!is.finite(rmin_c) || rmin_c < 0) rmin_c <- 0
    shiny::validate(shiny::need(is.finite(rmax_c) && rmax_c > 0, "Could not form a common r-grid."))
    # See the Ripley note: grid density comes from the stored curve.
    rr <- seq(rmin_c, rmax_c, length.out = max(10L, length(unique(res$r))))

    for (g in unique(res$group)) {
      files <- unique(res$dataset[res$group == g])
      mat <- vapply(files, function(f) {
        d <- res[res$dataset == f, ]; d <- d[order(d$r), ]
        stats::approx(d$r, d$g, xout = rr, rule = 2)$y
      }, numeric(length(rr)))
      if (is.null(dim(mat))) mat <- matrix(mat, ncol = 1)
      meanc <- smooth_vec(rr, rowMeans(mat, na.rm = TRUE))
      col <- unname(pal[g]); nfl <- length(files)
      # ── BAND MUST MATCH THE DRAWN MEAN ──────────────────────────────────
      # Three faults were fixed here together:
      #  1. the mean was smoothed and the band was not, so a smooth line sat in a
      #     jagged band (measured: band edge 16-18x the roughness of the mean).
      #     Same argument as the CSR envelope above - a band and a curve that got
      #     different treatment are not comparable.
      #  2. the SD band was built around the RAW rowMeans while the SMOOTHED mean
      #     was drawn, so the line was not centred in its own band (up to 22 nm
      #     off). It is now centred on the curve actually shown.
      #  3. quartiles need >= 4 values. On 2 files the IQR band is the middle half
      #     of a 2-point range and contains NEITHER file (measured 0% coverage),
      #     which draws a reassuringly tight band over maximal disagreement. Under
      #     4 files the full min-max range is drawn instead, and the legend says so.
      band_lbl <- ""
      if (band != "none" && nfl >= 2) {
        if (band == "iqr" && nfl >= 4L) {
          lo <- smooth_vec(rr, apply(mat, 1, stats::quantile, 0.25, na.rm = TRUE))
          hi <- smooth_vec(rr, apply(mat, 1, stats::quantile, 0.75, na.rm = TRUE))
          band_lbl <- ", IQR"
        } else if (band == "iqr") {
          lo <- smooth_vec(rr, apply(mat, 1, min, na.rm = TRUE))
          hi <- smooth_vec(rr, apply(mat, 1, max, na.rm = TRUE))
          band_lbl <- ", range"
        } else {
          sdv <- smooth_vec(rr, apply(mat, 1, stats::sd, na.rm = TRUE))
          lo <- meanc - sdv; hi <- meanc + sdv
          band_lbl <- ", ±1 SD"
        }
        fig <- fig %>%
          plotly::add_trace(x = rr, y = hi, type = "scatter", mode = "lines",
            line = list(width = 0), legendgroup = g, showlegend = FALSE, hoverinfo = "skip") %>%
          plotly::add_trace(x = rr, y = lo, type = "scatter", mode = "lines", fill = "tonexty",
            fillcolor = hex2rgba(col, 0.18), line = list(width = 0),
            legendgroup = g, showlegend = FALSE, hoverinfo = "skip")
      }
      # Individual files as shades of the group colour, behind the group mean.
      if (show_files && nfl >= 2) {
        for (f in files) {
          d <- res[res$dataset == f, ]; d <- d[order(d$r), ]
          fc <- unname(fpal[f]); if (is.na(fc)) fc <- col
          fig <- fig %>% plotly::add_trace(
            x = d$r, y = d$g, type = "scatter", mode = "lines",
            name = f, legendgroup = g, showlegend = FALSE,
            line = list(color = fc, width = 1), opacity = 0.75, hoverinfo = "skip")
        }
      }
      # CSR envelope for the group - same reasoning as the Ripley panel: the
      # group curve is the mean of its files' curves, so its null is the mean of
      # those files' envelopes. Dotted, so it reads as a reference rather than
      # as the solid between-file spread band.
      env_n <- 0L
      if (!is.null(env_res)) {
        ef <- intersect(files, unique(env_res$dataset))
        if (length(ef)) {
          elo <- vapply(ef, function(f) {
            e <- env_res[env_res$dataset == f, ]; e <- e[order(e$r), ]
            stats::approx(e$r, e$lo, xout = rr, rule = 2)$y }, numeric(length(rr)))
          ehi <- vapply(ef, function(f) {
            e <- env_res[env_res$dataset == f, ]; e <- e[order(e$r), ]
            stats::approx(e$r, e$hi, xout = rr, rule = 2)$y }, numeric(length(rr)))
          if (is.null(dim(elo))) { elo <- matrix(elo, ncol = 1); ehi <- matrix(ehi, ncol = 1) }
          env_n <- length(ef)
          # Two lines per group: the upper and lower bound of one band. Own
          # legend entry and legendgroup so they are identifiable and can be
          # toggled off without hiding the group curve.
          fig <- fig %>% plotly::add_trace(
            x = rr, y = smooth_vec(rr, rowMeans(ehi, na.rm = TRUE)),
            type = "scatter", mode = "lines",
            line = list(color = col, width = 1, dash = "dot"),
            name = sprintf("%s – CSR 95%% band", g),
            legendgroup = paste0(g, "__csr"), showlegend = TRUE, hoverinfo = "skip")
          fig <- fig %>% plotly::add_trace(
            x = rr, y = smooth_vec(rr, rowMeans(elo, na.rm = TRUE)),
            type = "scatter", mode = "lines",
            line = list(color = col, width = 1, dash = "dot"),
            legendgroup = paste0(g, "__csr"), showlegend = FALSE, hoverinfo = "skip")
        }
      }
      env_lbl <- if (env_n > 0L && env_n < nfl) sprintf(", CSR n=%d", env_n)
                 else if (env_n > 0L) ", CSR" else ""
      fig <- fig %>% plotly::add_trace(x = rr, y = meanc, type = "scatter", mode = "lines",
        name = sprintf("%s (n=%d%s%s)", g, nfl, band_lbl, env_lbl), legendgroup = g,
        line = list(color = col, width = 2.5), hoverinfo = "skip")
    }
  }

  # Same as the Ripley plot: the reference line spans what is actually drawn, so
  # a file with a larger r max cannot stretch the axis past every curve.
  x_drawn <- if (grouped) {
    c(rmin_c, if (isTRUE(show_files)) max(res$r, na.rm = TRUE) else rmax_c)
  } else range(res$r, na.rm = TRUE)

  fig %>%
    plotly::add_trace(x = x_drawn, y = c(1, 1), type = "scatter",
      mode = "lines", name = "CSR (g = 1)",
      line = list(color = "black", width = 1, dash = "dash"), hoverinfo = "skip") %>%
    plotly::layout(
      xaxis = list(title = .rmax_axis_title(grouped, res, rmax_c, show_files)),
      yaxis = list(title = "g(r)"),
      annotations = .global_p_annot(env_res),
      # Hover off, as on the Ripley plot: the curve is read as a shape against the
      # envelope, and a tooltip per series at every r covered the plot.
      hovermode = FALSE,
      legend = list(groupclick = "togglegroup"),
      paper_bgcolor = "white", plot_bgcolor = "white")
}
output$pcf_plot <- renderPlotly(build_pcf_plot_fig())

## Summary numbers - the actual endpoint you would compare between conditions.
## Reports the MEAN g over a user-set window rather than only the peak: a max
## over ~50 bins is upward-biased at low unit counts, so peak g is not safely
## comparable between files of different density even though g itself is.
## Per-FILE scalars first, then the group summary - never scalars from pooled rows.
##
## WHY THIS ORDER MATTERS. A group's rows are several files' curves concatenated, so
## reducing them together makes the CELL stop being the unit of replication: the peak
## becomes a max over files, and a denser file contributes more rows and so more
## weight. Reducing each file first and then averaging the scalars keeps one vote per
## cell, which is the standard unit of replication in SMLM, and it yields a spread
## across cells - the only thing that says whether a between-condition difference is
## bigger than the between-cell variability.
##
## Measured on the four real a real imaging file acquisitions (mean g over 30-300 nm, common r grid):
## cond01 1.24 / 1.50, cond02 2.31 / 4.38. Between-condition 2.44x against a
## within-condition spread of 1.21x and 1.89x, and no overlap - so averaging is
## justified there. But at r = 20 nm the two cond02 files read 5.03 and 14.32, so the
## MEAN CURVE (9.68) describes neither: read strengths and sizes off the per-file
## rows, and use the group row only for the comparison.
output$pcf_summary <- renderUI({
  res <- pcf_results()
  if (is.null(res) || !nrow(res)) return(NULL)
  r1 <- as.numeric(input$pcf_win_lo %||% 0)
  r2 <- as.numeric(input$pcf_win_hi %||% 0)
  if (!is.finite(r1) || !is.finite(r2) || r2 <= r1) {
    rmx <- max(res$r, na.rm = TRUE); r1 <- 0.05 * rmx; r2 <- 0.4 * rmx
  }
  grouped <- (input$analysis_by %||% "file") == "group"

  # One row per FILE, always computed on the file's own curve.
  per_file <- do.call(rbind, lapply(split(res, res$dataset), function(d) {
    d  <- d[order(d$r), ]
    ok <- is.finite(d$g) & d$r >= r1 & d$r <= r2
    pk <- if (any(is.finite(d$g))) which.max(replace(d$g, !is.finite(d$g), -Inf)) else NA_integer_
    # Size readout: first radius AFTER the peak where g returns to ~1. NA when the
    # curve has not settled inside the r range, which is itself informative -
    # the structure is larger than the window can resolve.
    rz <- NA_real_
    if (!is.na(pk)) {
      ti <- which(seq_along(d$r) > pk & is.finite(d$g) & d$g <= 1.05)
      if (length(ti)) rz <- d$r[ti[1]]
    }
    data.frame(file = d$dataset[1],
               g_win = if (any(ok)) mean(d$g[ok]) else NA_real_,
               excess = if (any(ok)) sum((d$g[ok] - 1) * c(diff(d$r[ok]), 0)) else NA_real_,
               r_settle = rz, rmax = max(d$r, na.rm = TRUE),
               stringsAsFactors = FALSE)
  }))
  gm <- if (grouped) file_group_map() else character(0)
  per_file$group <- if (grouped) unname(gm[per_file$file]) else per_file$file
  per_file <- per_file[!is.na(per_file$group), , drop = FALSE]
  if (!nrow(per_file)) return(NULL)

  num <- function(v, dg = 2) if (!is.finite(v)) "-" else formatC(v, format = "f", digits = dg)
  hdr <- function(cols) tags$tr(lapply(cols, function(h)
    tags$th(style = "text-align:left; padding:2px 12px 2px 0; border-bottom:1px solid #ddd;", h)))
  row_td <- function(...) tags$tr(lapply(list(...), function(v)
    tags$td(style = "padding:1px 12px 1px 0;", v)))

  file_tbl <- tags$table(style = "border-collapse:collapse;",
    hdr(c("File", if (grouped) "Group" else NULL,
          sprintf("mean g (%.0f-%.0f nm)", r1, r2), "∫(g-1)dr (nm)",
          "r where g→1 (nm)", "r max")),
    lapply(seq_len(nrow(per_file)), function(i) {
      args <- list(per_file$file[i])
      if (grouped) args <- c(args, list(per_file$group[i]))
      args <- c(args, list(num(per_file$g_win[i]),
                           num(per_file$excess[i], 0),
                           if (is.finite(per_file$r_settle[i])) num(per_file$r_settle[i], 0)
                           else sprintf("> %.0f", per_file$rmax[i]),
                           num(per_file$rmax[i], 0)))
      do.call(row_td, args)
    }))

  out <- list(tags$div(style = "font-weight:600; margin:2px 0 3px;",
                       "Per file (one row per cell - this is the unit of replication)"),
              file_tbl)

  if (grouped) {
    ag <- do.call(rbind, lapply(split(per_file, per_file$group), function(d) data.frame(
      group = d$group[1], n = nrow(d),
      mean_g = mean(d$g_win, na.rm = TRUE),
      sd_g   = if (nrow(d) > 1) stats::sd(d$g_win, na.rm = TRUE) else NA_real_,
      lo = min(d$g_win, na.rm = TRUE), hi = max(d$g_win, na.rm = TRUE),
      stringsAsFactors = FALSE)))
    grp_tbl <- tags$table(style = "border-collapse:collapse;",
      hdr(c("Group", "n cells", sprintf("mean g (%.0f-%.0f nm)", r1, r2),
            "SD across cells", "range")),
      lapply(seq_len(nrow(ag)), function(i) row_td(
        ag$group[i], ag$n[i], num(ag$mean_g[i]),
        if (is.finite(ag$sd_g[i])) num(ag$sd_g[i]) else "-",
        sprintf("%s - %s", num(ag$lo[i]), num(ag$hi[i])))))
    # The decisive comparison: is the between-group gap bigger than the within-group
    # spread? With no overlap it is, and that does not need a test to state.
    verdict <- NULL
    if (nrow(ag) == 2 && all(is.finite(c(ag$lo, ag$hi)))) {
      sep <- ag$lo[which.max(ag$mean_g)] > ag$hi[which.min(ag$mean_g)]
      ratio <- max(ag$mean_g) / min(ag$mean_g)
      verdict <- tags$p(style = paste0("margin:5px 0 0;",
                          if (sep) "color:#1a8a1a;" else "color:#b26a00;"),
        sprintf("Group means differ %.2f×; per-cell ranges %s.", ratio,
                if (sep) "do NOT overlap, so the difference exceeds between-cell variability"
                else "OVERLAP, so this difference is within between-cell variability"),
        if (min(ag$n) < 3) " With fewer than 3 cells per group this is descriptive only - no test is meaningful." else "")
    }
    out <- c(out, list(tags$div(style = "font-weight:600; margin:10px 0 3px;",
                                "Per group (mean of the per-cell values)"),
                      grp_tbl, verdict))
  }

  out <- c(out, list(tags$p(style = "margin:8px 0 0; color:#7a5c00;",
    "g = 1 is complete randomness. Read STRENGTH from mean g and SIZE from the radius where g returns to 1. ",
    "A value shown as '> r max' means the structure is larger than this window can resolve - set a fixed r max so files are comparable. ",
    "Do not read a size off a group mean curve: if the cells differ, the average describes none of them.")))
  do.call(tags$div, c(list(style = "font-size:0.85em;"), out))
})

output$download_pcf <- downloadHandler(
  filename = function() paste0("pcf_g_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
  content = function(file) {
    res <- pcf_results()
    shiny::validate(shiny::need(!is.null(res) && nrow(res) > 0, "No results."))
    write.csv(res, file, row.names = FALSE)
  }
)
