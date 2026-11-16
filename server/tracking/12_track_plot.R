## ─────────────────────────────────────────────────────────────────────────────
## server/tracking/12_track_plot.R  -  track plot
## ─────────────────────────────────────────────────────────────────────────────
## Part of tracking_server(). Not a standalone script: its statements are evaluated in
## the server function's environment, in file-name order (server/parts.R).
##
## output$trk_track2D, the proxy draw functions of the play loop, and the
## dwell heatmap store.
##
## INTERFACE (generated: Rscript tools/server_map.R --headers; do not edit by hand)
##   DEFINES for other parts:
##     draw_frame, draw_tid_point
##   READS from other parts:
##     01_state_sigma: rv, slider_guard, trace_map, trk_rechunk_bbox
##     03_playback_state_zcal: is_seq, loop, trk_z_applied
##     09_overlay_settings_io: trk_image_overlay
##     10_triggers: layout_trigger, plot_trigger, trk_axis_x_span,
##                  trk_data_key, trk_locked_ranges, trk_plot_ready, use_3d
##     11_data_frames: frames_store, seq_data, trk
## END INTERFACE
## ─────────────────────────────────────────────────────────────────────────────

# ---- Initial render ----
output$trk_track2D <- renderPlotly({
  req(!is.na(layout_trigger()))
  # Record the EFFECTIVE view-mode / z-stretch / chunk this render used (with the
  # same defaults the render applies). The view_mode / z_scale / chunk_sel
  # observers compare against these and skip firing when unchanged - so their
  # INPUT CREATION after an upload (value == the just-rendered default) never
  # triggers a redundant rebuild, without any fragile time window.
  loop$rendered_view  <- isolate(input$trk_view_mode  %||% "2d")
  loop$rendered_zs    <- isolate(trk_z_applied())
  loop$rendered_chunk <- isolate(input$trk_chunk_sel  %||% "__all__")
  pre  <- isolate(frames_store()); req(is.list(pre), length(pre$frames) > 0)
  f1    <- pre$frames[[1]]
  tids  <- pre$tids
  cmap  <- pre$col_map
  sp    <- pre$sp
  use3d <- isTRUE(pre$has_3d)

  # 2D viewing plane the frames were built in (XY / XZ / YZ). The coordinates are
  # already remapped in frames_store, so here we only derive axis labels and decide
  # whether the XY-physical overlays (dwell heatmap, background image) still apply.
  plane     <- pre$plane %||% "XY"
  loop$rendered_plane <- plane
  plane_lab <- switch(plane,
    XZ = list(x = "X (nm)", y = "Z (nm)"),
    YZ = list(x = "Y (nm)", y = "Z (nm)"),
    list(x = "X (nm)", y = "Y (nm)"))
  plane_xy  <- identical(plane, "XY")

  fig <- plot_ly(source = "track")

  # Faint background of ALL raw tracks (pre$sp_raw, when the Display toggle is on).
  # Drawn first so it sits behind everything; consolidated into one trace.
  draw_raw_faint <- function(fig) {
    spr <- pre$sp_raw
    if (is.null(spr) || !length(spr)) return(fig)
    dk  <- isTRUE(isolate(input$trk_dark_mode))
    col <- if (dk) "rgba(255,255,255,0.10)" else "rgba(120,120,120,0.12)"
    nms <- names(spr)
    rx  <- unlist(lapply(nms, function(nm) c(spr[[nm]]$x, NA_real_)), use.names = FALSE)
    ry  <- unlist(lapply(nms, function(nm) c(spr[[nm]]$y, NA_real_)), use.names = FALSE)
    tr  <- list(x = rx, y = ry, type = if (use3d) "scatter3d" else "scattergl", mode = "lines",
                line = list(width = 1, color = col), hoverinfo = "skip",
                showlegend = TRUE, name = "All raw tracks")
    if (use3d) tr$z <- unlist(lapply(nms, function(nm) c(spr[[nm]]$z, NA_real_)), use.names = FALSE)
    do.call(add_trace, c(list(fig), tr))
  }

  # Dwell-time heatmap (2D only) - shared by FULL-tracks and PLAYBACK modes so the
  # toggle works in both. heatmap_store() is computed from trk() on plot_trigger
  # (not playback-specific), so it has data in either mode. Added as a background
  # layer (before the tracks) so paths draw on top.
  # The dwell-time heatmap is computed in real XY space, so it is only meaningful in
  # the XY plane - suppress it (like the background image) in the XZ/YZ side views.
  hm_data <- if (!use3d && plane_xy && isTRUE(isolate(input$trk_show_heatmap)))
               tryCatch(isolate(heatmap_store()), error = function(e) NULL) else NULL
  draw_heatmap <- function(fig) {
    if (is.null(hm_data)) return(fig)
    z_proc <- apply_heatmap_contrast(hm_data$z,
                clip = isolate(input$trk_hm_clip_pct), gamma = isolate(input$trk_hm_gamma))
    fig %>% add_heatmap(x = hm_data$x, y = hm_data$y, z = z_proc, zmin = 0, zmax = 1,
                        showscale = TRUE, opacity = isolate(input$trk_hm_opacity),
                        colorscale = hm_colorscale(isolate(input$trk_hm_colors)),
                        hoverinfo = "skip")
  }

  # ── Full-tracks static mode ───────────────────────────────────────────────
  if (identical(isolate(input$trk_view_mode_main), "full")) {
    connect  <- isTRUE(isolate(input$trk_full_connect))
    sp_all   <- if (!is.null(pre$sp_static)) pre$sp_static else sp
    all_tids <- names(sp_all)
    dark     <- isTRUE(isolate(input$trk_dark_mode))
    bg_col   <- if (dark) "#000000" else "#FFFFFF"
    fg_col   <- if (dark) "#FFFFFF" else "#63666A"
    cmode_st <- isolate(input$trk_color_mode %||% "tid")
    pcol_st  <- switch(cmode_st, dark = "col_dark", speed = "col_speed", NULL)
    color_by_dark <- !is.null(pcol_st) &&
                     any(vapply(sp_all, function(p) pcol_st %in% names(p), logical(1)))
    dark_line_col <- if (dark) "rgba(200,200,200,0.30)" else "rgba(120,120,120,0.30)"
    # Integration-status display filter: keep all / integrated (n_dark==0) / dark (>=1).
    dark_show <- isolate(input$trk_dark_show %||% "all")
    dark_keep <- function(pts) {
      if (dark_show == "all" || !("n_dark" %in% names(pts))) return(pts)
      m <- if (dark_show == "integrated") pts$n_dark == 0L else pts$n_dark >= 1L
      pts[m, , drop = FALSE]
    }

    fig <- draw_heatmap(fig)     # dwell-time heatmap background (2D only; no-op otherwise)
    fig <- draw_raw_faint(fig)   # all-raw-tracks background (behind the filtered tracks)

    if (connect) {
      # One line+marker trace per TID (allows per-TID color on lines + hover)
      tr_type <- if (use3d) "scatter3d" else "scattergl"
      for (i in seq_along(all_tids)) {
        nm  <- all_tids[i]
        pts <- sp_all[[nm]]
        if (!nrow(pts)) next
        pts <- pts[order(pts$t), ]
        pts <- dark_keep(pts)
        if (!nrow(pts)) next
        col_i    <- as.character(cmap[[nm]])
        tim_col  <- if ("tim" %in% names(pts)) pts$tim else pts$t
        hover_z  <- if (use3d && "z" %in% names(pts))
                      sprintf("<br>z: %.1f nm", pts$z) else
                      rep("", nrow(pts))
        hover_d  <- if ("n_dark" %in% names(pts))
                      sprintf("<br>dark cycles: %d", pts$n_dark) else rep("", nrow(pts))
        hover_txt <- sprintf("TID: %s<br>t: %.1f ms<br>x: %.1f nm<br>y: %.1f nm%s%s",
                             nm, tim_col, pts$x, pts$y, hover_z, hover_d)
        mk_col <- if (color_by_dark && pcol_st %in% names(pts)) as.character(pts[[pcol_st]]) else col_i
        ln_col <- if (color_by_dark) dark_line_col else col_i
        tr_args <- list(
          x = pts$x, y = pts$y,
          type = tr_type, mode = "lines+markers",
          line   = list(width = 1, color = ln_col),
          marker = list(size = isolate(input$trk_pt_size), color = mk_col),
          text = hover_txt, hovertemplate = "%{text}<extra></extra>",
          showlegend = (i == 1L && !color_by_dark), legendgroup = "Localizations",
          name = if (i == 1L) "Localizations" else paste0("tid_", nm)
        )
        if (use3d) tr_args$z <- pts$z
        fig <- do.call(add_trace, c(list(fig), tr_args))
      }
      # Connect mode draws one trace per TID with hex colours, so it gets no automatic
      # colourbar - append the scale explicitly.
      if (color_by_dark) {
        a1 <- sp_all[[all_tids[1]]]
        fig <- add_colorbar_trace(fig, cmode_st, fg_col,
                                  vmax = suppressWarnings(as.numeric(pre$speed_max)),
                                  use3d = use3d,
                                  ax = a1$x[1], ay = a1$y[1],
                                  az = if (use3d && "z" %in% names(a1)) a1$z[1] else 0)
      }
    } else {
      # Single fast scattergl trace with per-point colors + hover text
      tr_type <- if (use3d) "scatter3d" else "scattergl"
      x_v <- c(); y_v <- c(); z_v <- c(); col_v <- c(); nd_v <- c(); vd_v <- c(); txt_v <- c()
      for (nm in all_tids) {
        pts   <- sp_all[[nm]]
        if (!nrow(pts)) next
        pts <- dark_keep(pts)
        if (!nrow(pts)) next
        col_i   <- as.character(cmap[[nm]])
        tim_col <- if ("tim" %in% names(pts)) pts$tim else pts$t
        hover_z <- if (use3d && "z" %in% names(pts))
                     sprintf("<br>z: %.1f nm", pts$z) else
                     rep("", nrow(pts))
        hover_d <- if ("n_dark" %in% names(pts))
                     sprintf("<br>dark cycles: %d", pts$n_dark) else rep("", nrow(pts))
        hover_v <- if ("v_drift" %in% names(pts))
                     ifelse(is.finite(pts$v_drift),
                            sprintf("<br>drift speed: %.2f µm/s", pts$v_drift),
                            "<br>drift speed: n/a") else rep("", nrow(pts))
        txt_v <- c(txt_v, sprintf("TID: %s<br>t: %.1f ms<br>x: %.1f nm<br>y: %.1f nm%s%s%s",
                                  nm, tim_col, pts$x, pts$y, hover_z, hover_d, hover_v))
        x_v   <- c(x_v,   pts$x)
        y_v   <- c(y_v,   pts$y)
        if (use3d) z_v <- c(z_v, pts$z)
        col_v <- c(col_v, rep(col_i, nrow(pts)))
        if ("n_dark" %in% names(pts))  nd_v <- c(nd_v, pts$n_dark)
        if ("v_drift" %in% names(pts)) vd_v <- c(vd_v, pts$v_drift)
      }
      # Dark-cycle colouring uses a numeric colour array + colorscale so plotly draws
      # a self-documenting colorbar (0 = integrated/blue → deeper red = more dark cycles).
      if (identical(cmode_st, "dark") && color_by_dark && length(nd_v) == length(x_v)) {
        mk <- list(size = isolate(input$trk_pt_size), color = nd_v,
                   colorscale = dark_cycle_colorscale(), cmin = 0, cmax = DARK_NMAX,
                   showscale = TRUE,
                   colorbar = list(title = list(text = "dark<br>cycles", font = list(color = fg_col)),
                                   tickmode = "array", tickvals = 0:DARK_NMAX,
                                   ticktext = c(as.character(0:(DARK_NMAX - 1L)), paste0(DARK_NMAX, "+")),
                                   tickfont = list(color = fg_col), len = 0.5, thickness = 14))
        show_lg <- FALSE
      } else if (identical(cmode_st, "speed") && color_by_dark && length(vd_v) == length(x_v)) {
        vmax_st <- suppressWarnings(as.numeric(pre$speed_max))
        if (!is.finite(vmax_st) || vmax_st <= 0)
          vmax_st <- suppressWarnings(stats::quantile(vd_v, 0.95, na.rm = TRUE))
        if (!is.finite(vmax_st) || vmax_st <= 0) vmax_st <- 1
        mk <- list(size = isolate(input$trk_pt_size), color = vd_v,
                   colorscale = speed_colorscale(), cmin = 0, cmax = vmax_st,
                   showscale = TRUE,
                   colorbar = list(title = list(text = "drift speed<br>(µm/s)",
                                               font = list(color = fg_col)),
                                   tickfont = list(color = fg_col), len = 0.5, thickness = 14))
        show_lg <- FALSE
      } else {
        mk <- list(size = isolate(input$trk_pt_size), color = col_v)
        show_lg <- TRUE
      }
      tr_args <- list(
        x = x_v, y = y_v,
        type = tr_type, mode = "markers",
        marker = mk,
        text = txt_v, hovertemplate = "%{text}<extra></extra>",
        showlegend = show_lg, name = "Localizations", legendgroup = "Localizations"
      )
      if (use3d) tr_args$z <- z_v
      fig <- do.call(add_trace, c(list(fig), tr_args))
    }

    # Axis layout (reuse shared logic below)
    all_x <- unlist(lapply(sp_all, `[[`, "x"), use.names = FALSE)
    all_y <- unlist(lapply(sp_all, `[[`, "y"), use.names = FALSE)
    all_x <- all_x[is.finite(all_x)]; all_y <- all_y[is.finite(all_y)]
    hide_grid <- isTRUE(isolate(input$trk_hide_grid))
    grid_col  <- if (dark) "rgba(255,255,255,0.18)" else "rgba(0,0,0,0.18)"
    axis_base <- list(showgrid = !hide_grid, zeroline = !hide_grid,
                      showticklabels = !hide_grid,
                      gridcolor = grid_col, color = fg_col)
    ovl <- if (plane_xy) tryCatch(isolate(trk_image_overlay()), error = function(e) NULL) else NULL

    if (use3d) {
      all_z <- unlist(lapply(sp_all, `[[`, "z"), use.names = FALSE)
      lim3  <- get_equal_limits(all_x, all_y, all_z)
      ax3   <- c(list(backgroundcolor = bg_col, autorange = FALSE), axis_base)
      fig   <- fig %>% layout(
        paper_bgcolor = bg_col,
        showlegend = TRUE,                      # always show the "Localizations" toggle
        legend = list(font = list(color = fg_col), groupclick = "togglegroup"),
        scene = list(
          bgcolor = bg_col,
          xaxis = c(list(title = if (hide_grid) "" else "X (nm)", range = lim3$x), ax3),
          yaxis = c(list(title = if (hide_grid) "" else "Y (nm)", range = lim3$y), ax3),
          zaxis = c(list(title = if (hide_grid) "" else "Z (nm)", range = lim3$z), ax3)
        )
      )
    } else {
      xspan <- if (length(all_x) > 1) diff(range(all_x)) else 1000
      yspan <- if (length(all_y) > 1) diff(range(all_y)) else 1000
      half  <- max(xspan, yspan) / 2 * 1.04
      xctr  <- if (length(all_x) > 1) mean(range(all_x)) else 0
      yctr  <- if (length(all_y) > 1) mean(range(all_y)) else 0
      xr    <- c(xctr - half, xctr + half)
      yr    <- c(yctr - half, yctr + half)
      trk_locked_ranges(list(x = xr, y = yr))
      fig <- fig %>% layout(
        paper_bgcolor = bg_col, plot_bgcolor = bg_col,
        xaxis = c(list(title = if (hide_grid) "" else plane_lab$x, range = xr), axis_base),
        yaxis = c(list(title = if (hide_grid) "" else plane_lab$y, range = yr,
                       scaleanchor = "x", scaleratio = 1), axis_base),
        images = if (!is.null(ovl)) list(ovl) else list(),
        # showlegend=TRUE forces the legend to appear even with a single track trace
        # (plotly hides it otherwise, and the heatmap adds only a colorbar) so the
        # "Localizations" entry is always clickable to reveal the dwell-time heatmap.
        # groupclick=togglegroup so one click hides ALL track traces (connect mode
        # draws one trace per TID).
        showlegend = TRUE,
        legend = list(font = list(color = fg_col), groupclick = "togglegroup"),
        uirevision = isolate(trk_data_key())
      )
      # Diagnostic markers (same inline logic as the main render path)
      diag <- isolate(rv$trk_diag)
      if (!is.null(diag) && isTRUE(isolate(input$trk_show_diag))) {
        if (!is.null(diag$data_clusters) && nrow(diag$data_clusters) > 0)
          fig <- fig %>% plotly::add_trace(
            x = diag$data_clusters$x, y = diag$data_clusters$y,
            type = "scatter", mode = "markers",
            marker = list(symbol = "cross", size = 16, color = "red",
                          line = list(color = "red", width = 3)),
            name = "MINFLUX clusters", showlegend = TRUE, hoverinfo = "skip")
        if (!is.null(diag$img_blobs) && nrow(diag$img_blobs) > 0)
          fig <- fig %>% plotly::add_trace(
            x = diag$img_blobs$x, y = diag$img_blobs$y,
            type = "scatter", mode = "markers",
            marker = list(symbol = "circle-open", size = 20, color = "#00FF80",
                          line = list(color = "#00FF80", width = 3)),
            name = "Confocal blobs", showlegend = TRUE, hoverinfo = "skip")
      }
      # BAKE the scale bar into the figure (paper-coord shapes), exactly like the
      # SMLM tab - reading inputs ISOLATED so toggling the checkbox does NOT invalidate
      # this render (previously trk_push_scale_bar()'s non-isolated input reads made the
      # plot depend on the checkbox → a rebuild that wiped the proxy-added bar). The
      # checkbox observer still proxy-pushes it live; a genuine rebuild re-bakes it here
      # so it never vanishes. Width uses the rendered x-span; zoom updates via the
      # plotly_relayout observer. trk_axis_x_span is seeded so the live proxy matches.
      trk_axis_x_span(diff(xr))
      if (isTRUE(isolate(input$trk_show_scale_bar))) {
        sb_len   <- suppressWarnings(as.numeric(isolate(input$trk_scale_bar_length)))
        sb_thick <- suppressWarnings(as.numeric(isolate(input$trk_scale_bar_thickness)))
        sb_font  <- suppressWarnings(as.numeric(isolate(input$trk_scale_font_size)))
        sb <- make_scale_bar(sb_len, diff(xr), sb_thick, dark = dark,
                             scale_font_size = if (is.finite(sb_font) && sb_font > 0) sb_font else 18)
        if (!is.null(sb))
          fig <- fig %>% plotly::layout(shapes = sb$shapes, annotations = sb$annotations)
      }
    }
    trk_plot_ready(TRUE)
    return(fig %>% plotly::config(scrollZoom = TRUE, displayModeBar = TRUE) %>%
             plotly::event_register('plotly_relayout'))
  }
  # ── End full-tracks static mode ───────────────────────────────────────────

  # ---- Heatmap (optional, 2D only) - shared helper (see draw_heatmap above) ----
  fig    <- draw_heatmap(fig)
  n_heat <- if (!is.null(hm_data)) 1L else 0L

  # ---- Faint static traces (optional) ----
  # Consolidated into ONE trace with NA separators instead of one per TID.
  # This drops per-TID colouring (opaciOkaty 0.15 grey is visually equivalent),
  # but reduces WebGL trace count from O(N_TIDs) to O(1), which dramatically
  # shortens the initial render time and eliminates the proxy-during-init blank.
  show_faint <- isTRUE(isolate(input$trk_show_full_tracks_faint))
  if (show_faint) {
    # Faint background = the full (post-filtered, post-split) paths of exactly the
    # TIDs being ANIMATED (chunk-scoped), so it never shows tracks that aren't
    # playing. With "All TIDs" selected this is every filtered track.
    #
    # pre$tids is windowed to the playback window [tmin, tmax], but sequential
    # mode plays EVERY chunk TID (seq_data$tids), not just the windowed subset.
    # Use the full played set in sequential mode so every coloured trail that
    # runs has a faint full-track underneath it.
    faint_tids <- pre$tids
    if (isolate(is_seq())) {
      sd_faint <- isolate(seq_data())
      if (is.list(sd_faint) && length(sd_faint$tids)) faint_tids <- sd_faint$tids
    }
    sp_st    <- if (!is.null(pre$sp_static))
                  pre$sp_static[intersect(faint_tids, names(pre$sp_static))] else sp
    all_tids <- names(sp_st)
    static_x <- unlist(lapply(all_tids, function(nm) c(sp_st[[nm]]$x, NA_real_)), use.names = FALSE)
    static_y <- unlist(lapply(all_tids, function(nm) c(sp_st[[nm]]$y, NA_real_)), use.names = FALSE)
    static_z <- if (use3d) unlist(lapply(all_tids, function(nm) c(sp_st[[nm]]$z, NA_real_)), use.names = FALSE) else NULL
    tr_args <- list(
      x = static_x, y = static_y,
      type = if (use3d) "scatter3d" else "scattergl", mode = "lines",
      line = list(width = 1, color = if (isTRUE(isolate(input$trk_dark_mode))) {
        if (use3d) "rgba(255,255,255,0.55)" else "rgba(255,255,255,0.25)"
      } else {
        if (use3d) "rgba(80,80,80,0.55)" else "rgba(80,80,80,0.15)"
      }),
      hoverinfo = "skip", showlegend = TRUE, name = "Full tracks"
    )
    if (use3d) tr_args$z <- static_z
    fig <- do.call(add_trace, c(list(fig), tr_args))
  }
  n_static <- if (show_faint) 1L else 0L

  # ---- Faint ALL-RAW-tracks background (optional; one consolidated trace) ----
  has_raw_bg <- !is.null(pre$sp_raw) && length(pre$sp_raw) > 0
  if (has_raw_bg) fig <- draw_raw_faint(fig)
  n_raw <- if (has_raw_bg) 1L else 0L

  # IMPORTANT: total background traces BEFORE dynamic ones
  n_bg <- n_heat + n_static + n_raw

  # ---- Dynamic setup ----
  show_trails <- isTRUE(isolate(input$trk_show_moving_trails))

  if (isolate(is_seq())) {
    # Sequential playback = the moving trail + a leading marker sweeping through
    # each TID in chronological order. The trail IS the playback here (it's how
    # you watch a TID light up), so it is always present regardless of the
    # simultaneous-mode "Animate moving trails" toggle. The full-track context
    # is drawn faintly underneath (see the faint-tracks block above).
    first_tid <- tids[1]; col1 <- as.character(cmap[[first_tid]])
    tr_type   <- if (use3d) "scatter3d" else "scattergl"
    trail_args  <- list(x = numeric(0), y = numeric(0), type = tr_type, mode = "lines",
                        line = list(width = 3, color = col1), hoverinfo = "skip",
                        showlegend = TRUE, legendgroup = "Localizations", name = "Localizations")
    marker_args <- list(x = numeric(0), y = numeric(0), type = tr_type, mode = "markers",
                        marker = list(size = isolate(input$trk_pt_size), color = col1),
                        hoverinfo = "skip", showlegend = FALSE,
                        legendgroup = "Localizations", name = "marker_seq")
    if (use3d) { trail_args$z <- numeric(0); marker_args$z <- numeric(0) }
    fig <- do.call(add_trace, c(list(fig), trail_args))
    fig <- do.call(add_trace, c(list(fig), marker_args))

    trace_map(list(
      trails_mode = "sequential",
      trail_start = n_bg,
      n_tid       = 0L,
      markers_idx = n_bg + 1L,
      static_idx  = if (show_faint) n_heat else NA_integer_
    ))

  } else {
    if (show_trails) {
      # One line trace per TID, allows per-TID color without per-point tricks.
      # Markers are hidden while trails are active.
      tr_type <- if (use3d) "scatter3d" else "scattergl"
      for (i in seq_along(tids)) {
        col_i   <- as.character(cmap[[tids[i]]])
        tr_args <- list(x = numeric(0), y = numeric(0), type = tr_type, mode = "lines",
                        line = list(width = 3, color = col_i), hoverinfo = "skip",
                        showlegend = (i == 1L), legendgroup = "Localizations",
                        name = if (i == 1L) "Localizations" else paste0("trail_", i))
        if (use3d) tr_args$z <- numeric(0)
        fig <- do.call(add_trace, c(list(fig), tr_args))
      }
      trail_mode  <- "per_tid"
      trail_start <- n_bg
      markers_idx <- n_bg + length(tids)
    } else {
      trail_mode  <- "none"
      trail_start <- NA_integer_
      markers_idx <- n_bg
    }

    # markers - empty when trails are shown (lines replace dots)
    mk_type <- if (use3d) "scatter3d" else "scattergl"
    mk_args <- list(
      x = if (show_trails) numeric(0) else f1$markers_x,
      y = if (show_trails) numeric(0) else f1$markers_y,
      type = mk_type, mode = "markers",
      marker = list(size = isolate(input$trk_pt_size),
                    color = if (show_trails) "rgba(0,0,0,0)" else as.character(f1$marker_cols)),
      hoverinfo = "skip",
      showlegend = !show_trails,
      legendgroup = "Localizations",
      name = "Localizations"
    )
    if (use3d) mk_args$z <- if (show_trails) numeric(0) else (f1$markers_z %||% numeric(0))
    fig <- do.call(add_trace, c(list(fig), mk_args))

    trace_map(list(
      trails_mode = trail_mode,
      trail_start = trail_start,
      n_tid       = if (identical(trail_mode, "per_tid")) length(tids) else 0L,
      markers_idx = markers_idx,
      static_idx  = if (show_faint) n_heat else NA_integer_
    ))

    # Colour scale for playback. Added AFTER trace_map is set so the extra trace can
    # never disturb trail_start / markers_idx (draw_frame restyles by index).
    # fg_col is local to the static branch (which returns above), so derive it here.
    cmode_pb <- isolate(input$trk_color_mode %||% "tid")
    if (cmode_pb %in% c("dark", "speed")) {
      fg_pb <- if (isTRUE(isolate(input$trk_dark_mode))) "#FFFFFF" else "#63666A"
      # Anchor on a real localization (frame 1, else the first track) - an all-NA
      # trace is dropped by plotly and would draw no colourbar.
      a1x <- f1$markers_x[1]; a1y <- f1$markers_y[1]
      a1z <- if (use3d) (f1$markers_z %||% numeric(0))[1] else 0
      if (!length(a1x) || !is.finite(a1x %||% NA_real_)) {
        fb  <- if (length(tids)) sp[[tids[1]]] else NULL
        a1x <- if (!is.null(fb)) fb$x[1] else 0
        a1y <- if (!is.null(fb)) fb$y[1] else 0
        a1z <- if (!is.null(fb) && use3d && "z" %in% names(fb)) fb$z[1] else 0
      }
      fig <- add_colorbar_trace(fig, cmode_pb, fg_pb,
                                vmax = suppressWarnings(as.numeric(pre$speed_max)),
                                use3d = use3d, ax = a1x, ay = a1y, az = a1z)
    }
  }
  
  # Bounding box for axis locking.  After Filter & Rechunk, use the stored
  # full spatial extent of all remaining TIDs so the view stays stable while
  # playing through individual TIDs across chunks.  Without a rechunk bbox,
  # fall back to frame 1's markers (or the whole current chunk if frame 1 is empty).
  rechunk_bb <- isolate(trk_rechunk_bbox())
  if (!is.null(rechunk_bb)) {
    bbox_x <- rechunk_bb$x
    bbox_y <- rechunk_bb$y
  } else if (isolate(is_seq())) {
    # Sequential plays one TID at a time, so frame 1 is just the earliest TID's
    # first bin - using it would zoom the view into a tiny region and push the
    # faint full-track context off-screen. Use the full spatial extent of every
    # PLAYED track (seq_data$tids, matching the faint overlay) so every track
    # stays visible while stepping through TIDs.
    seq_tids <- pre$tids
    sd_bb <- isolate(seq_data())
    if (is.list(sd_bb) && length(sd_bb$tids)) seq_tids <- sd_bb$tids
    src <- if (!is.null(pre$sp_static))
             pre$sp_static[intersect(seq_tids, names(pre$sp_static))] else sp
    bbox_x <- unlist(lapply(src, `[[`, "x"), use.names = FALSE)
    bbox_y <- unlist(lapply(src, `[[`, "y"), use.names = FALSE)
    bbox_x <- bbox_x[is.finite(bbox_x)]; bbox_y <- bbox_y[is.finite(bbox_y)]
    if (!length(bbox_x) || !length(bbox_y)) {
      bbox_x <- unlist(lapply(sp, `[[`, "x"), use.names = FALSE)
      bbox_y <- unlist(lapply(sp, `[[`, "y"), use.names = FALSE)
      bbox_x <- bbox_x[is.finite(bbox_x)]; bbox_y <- bbox_y[is.finite(bbox_y)]
    }
  } else {
    bbox_x <- f1$markers_x[is.finite(f1$markers_x)]
    bbox_y <- f1$markers_y[is.finite(f1$markers_y)]
    if (!length(bbox_x) || !length(bbox_y)) {
      bbox_x <- unlist(lapply(sp, `[[`, "x"), use.names = FALSE)
      bbox_y <- unlist(lapply(sp, `[[`, "y"), use.names = FALSE)
      bbox_x <- bbox_x[is.finite(bbox_x)]; bbox_y <- bbox_y[is.finite(bbox_y)]
    }
  }
  have_bbox <- length(bbox_x) > 0 && length(bbox_y) > 0
  xr <- yr <- NULL
  if (have_bbox) {
    xr   <- range(bbox_x); yr <- range(bbox_y)
    xpad <- max(1, diff(xr) * 0.04); ypad <- max(1, diff(yr) * 0.04)
    if (!use3d) {
      # 2D only - locked ranges and scale bar depend on these
      trk_axis_x_span(diff(xr) + 2 * xpad)
      trk_locked_ranges(list(x = c(xr[1]-xpad, xr[2]+xpad),
                             y = c(yr[1]-ypad, yr[2]+ypad)))
    }
  }

  ovl       <- if (plane_xy) tryCatch(isolate(trk_image_overlay()), error = function(e) NULL) else NULL
  hide_grid <- isTRUE(isolate(input$trk_hide_grid))
  dark      <- isTRUE(isolate(input$trk_dark_mode))
  bg_col    <- if (dark) "#000000" else "#FFFFFF"
  fg_col    <- if (dark) "#FFFFFF" else "#63666A"
  grid_col  <- if (dark) "rgba(255,255,255,0.18)" else "rgba(0,0,0,0.18)"
  axis_grid <- list(
    showgrid       = !hide_grid, zeroline   = !hide_grid,
    showticklabels = !hide_grid,
    gridcolor  = grid_col, color = fg_col
  )
  if (use3d) {
    # 3D scene layout - uses scene instead of xaxis/yaxis
    ax3 <- list(color = fg_col, gridcolor = grid_col,
                showgrid = !hide_grid, zeroline = !hide_grid,
                showticklabels = !hide_grid, backgroundcolor = bg_col)

    # Equal axis ranges: use the same helper as the SMLM 3D plot.
    # get_equal_limits() finds the largest span across X, Y, Z and applies it
    # to all three axes (centred on each axis's own midpoint) - identical nm/tick.
    all_x <- unlist(lapply(sp, `[[`, "x"), use.names = FALSE)
    all_y <- unlist(lapply(sp, `[[`, "y"), use.names = FALSE)
    all_z <- unlist(lapply(sp, `[[`, "z"), use.names = FALSE)
    lim3  <- get_equal_limits(all_x, all_y, all_z)

    # autorange = FALSE is required - Plotly 3D ignores range unless autorange is
    # explicitly disabled.
    ax3_x <- c(list(title = if (hide_grid) "" else "X (nm)", range = lim3$x, autorange = FALSE), ax3)
    ax3_y <- c(list(title = if (hide_grid) "" else "Y (nm)", range = lim3$y, autorange = FALSE), ax3)
    ax3_z <- c(list(title = if (hide_grid) "" else "Z (nm)", range = lim3$z, autorange = FALSE), ax3)

    lay <- list(
      title        = list(text = "t = frame 1", font = list(color = fg_col)),
      paper_bgcolor = bg_col,
      font          = list(color = fg_col),
      scene = list(
        bgcolor     = bg_col,
        xaxis       = ax3_x,
        yaxis       = ax3_y,
        zaxis       = ax3_z,
        aspectmode  = "manual",
        aspectratio = list(x = 1, y = 1, z = 1)
      ),
      uirevision  = isolate(trk_data_key()),
      transition  = list(duration = 0)
    )
  } else {
    lay <- list(
      # No frame text in the title: updating it would need a relayout, which re-solves
      # scaleanchor and jitters the FOV. The live frame readout lives in trk_frame_meta.
      title = list(text = "", font = list(color = fg_col)),
      paper_bgcolor = bg_col, plot_bgcolor = bg_col,
      font = list(color = fg_col),
      xaxis = c(
        if (have_bbox) list(title = if (hide_grid) "" else plane_lab$x,
                            range = c(xr[1]-xpad, xr[2]+xpad), autorange = FALSE)
        else           list(title = if (hide_grid) "" else plane_lab$x),
        axis_grid
      ),
      yaxis = c(
        if (have_bbox) list(title = if (hide_grid) "" else plane_lab$y,
                            range = c(yr[1]-ypad, yr[2]+ypad), autorange = FALSE,
                            scaleanchor = "x", scaleratio = 1)
        else           list(title = if (hide_grid) "" else plane_lab$y,
                            scaleanchor = "x", scaleratio = 1),
        axis_grid
      ),
      uirevision = isolate(trk_data_key()),
      transition = list(duration = 0)
    )
    if (!is.null(ovl)) lay$images <- list(ovl)
    # ---- Scale bar (2D only) ----
    if (have_bbox && isTRUE(isolate(input$trk_show_scale_bar))) {
      sb_len   <- suppressWarnings(as.numeric(isolate(input$trk_scale_bar_length)))
      sb_thick <- suppressWarnings(as.numeric(isolate(input$trk_scale_bar_thickness)))
      sb_font  <- suppressWarnings(as.numeric(isolate(input$trk_scale_font_size)))
      sb <- make_scale_bar(sb_len, isolate(trk_axis_x_span()), sb_thick, dark = dark,
                           scale_font_size = if (is.finite(sb_font) && sb_font > 0) sb_font else 18)
      if (!is.null(sb)) {
        lay$shapes      <- sb$shapes
        lay$annotations <- sb$annotations
      }
    }
  }
  # ---- Diagnostic: detected cluster centres (red +) and confocal blobs (green o) ----
  # Use type="scatter" (SVG), not scattergl - scattergl lacks circle-open / cross symbols.
  diag <- isolate(rv$trk_diag)
  if (!use3d && isTRUE(isolate(input$trk_show_diag)) && !is.null(diag)) {
    if (!is.null(diag$data_clusters) && nrow(diag$data_clusters) > 0) {
      fig <- fig %>% plotly::add_trace(
        x = diag$data_clusters$x, y = diag$data_clusters$y,
        type = "scatter", mode = "markers",
        marker = list(symbol = "cross", size = 16, color = "red",
                      line = list(color = "red", width = 3)),
        name = "MINFLUX clusters", showlegend = TRUE, hoverinfo = "skip"
      )
    }
    if (!is.null(diag$img_blobs) && nrow(diag$img_blobs) > 0) {
      fig <- fig %>% plotly::add_trace(
        x = diag$img_blobs$x, y = diag$img_blobs$y,
        type = "scatter", mode = "markers",
        marker = list(symbol = "circle-open", size = 20, color = "#00FF80",
                      line = list(color = "#00FF80", width = 3)),
        name = "Confocal blobs", showlegend = TRUE, hoverinfo = "skip"
      )
    }
  }

  trk_plot_ready(TRUE)
  do.call(plotly::layout, c(list(fig), lay)) %>% plotly::event_register('plotly_relayout')
})

# After every layout rebuild, restore the locked axis ranges via proxy (2D only).
# lr is read INSIDE the delay so it captures the value written by renderPlotly
# (which runs in the same flush cycle as layout_trigger but finishes before the
# 120ms callback fires).  Reading it outside would capture the stale pre-render value.
observeEvent(layout_trigger(), {
  if (isTRUE(isolate(use_3d()))) return()   # 3D camera manages its own view
  shinyjs::delay(120, {
    lr <- isolate(trk_locked_ranges())
    if (!is.null(lr)) {
      plotlyProxyInvoke(
        plotlyProxy("trk_track2D", session), "relayout",
        list(`xaxis.range[0]` = lr$x[1], `xaxis.range[1]` = lr$x[2],
             `yaxis.range[0]` = lr$y[1], `yaxis.range[1]` = lr$y[2],
             `xaxis.autorange` = FALSE,   `yaxis.autorange` = FALSE)
      )
    }
  })
}, ignoreInit = TRUE)


# ---- Draw helpers ----
draw_frame <- function(k, force = FALSE) {
  if (identical(isolate(input$trk_view_mode_main), "full")) {
    return(invisible(NULL))
  }
  # Pull current prebuilt frames and mapping safely
  # Use if/return guards (not req())  -  req() throws inside later() callbacks
  pre <- isolate(frames_store())
  if (!is.list(pre) || !length(pre$frames)) { return(invisible(NULL)) }
  tm  <- isolate(trace_map())
  if (is.null(tm)) { return(invisible(NULL)) }
  pt_size <- suppressWarnings(isolate(input$trk_pt_size))

  # Respect play/pause unless explicitly forced
  if (!force && !(isTRUE(loop$gplaying) || isTRUE(loop$playing))) { return(invisible(NULL)) }
  
  # Clamp k and extract frame
  k  <- max(1L, min(as.integer(k), length(pre$frames)))
  fr <- pre$frames[[k]]
  
  proxy <- plotlyProxy("trk_track2D", session)

  is3d <- isTRUE(pre$has_3d)

  # ----- Trails (per-TID line traces, batched into ONE proxy call) -----
  if (identical(tm$trails_mode, "per_tid")) {
    tpt <- if (!is.null(fr$trails_per_tid)) fr$trails_per_tid else list()
    n   <- tm$n_tid
    xs  <- vector("list", n); ys <- vector("list", n); zs <- vector("list", n)
    # Per-point colour mode: recolour each trail LINE from its leading point, so the
    # trails carry the same colour meaning the markers would (blue→red speed, or the
    # dark-cycle gradient) and no dots are needed on top.
    per_pt <- isTRUE(pre$color_by_dark)
    lc <- if (per_pt) vector("list", n) else NULL
    for (i in seq_len(n)) {
      td      <- if (i <= length(tpt)) tpt[[i]] else NULL
      xs[[i]] <- if (!is.null(td) && length(td$x)) td$x else NA_real_
      ys[[i]] <- if (!is.null(td) && length(td$y)) td$y else NA_real_
      if (is3d) zs[[i]] <- if (!is.null(td) && length(td$z)) td$z else NA_real_
      if (per_pt) {
        # Empty trail -> transparent (its x/y are NA so nothing is drawn anyway).
        lc[[i]] <- if (!is.null(td) && !is.null(td$col) && nzchar(td$col)) td$col else "rgba(0,0,0,0)"
      }
    }
    restyle_data <- list(x = xs, y = ys)
    if (is3d) restyle_data$z <- zs
    if (per_pt) restyle_data$`line.color` <- lc
    plotlyProxyInvoke(proxy, "restyle", restyle_data,
      seq(tm$trail_start, tm$trail_start + n - 1L))
    # Markers stay empty while trails are shown - the trails themselves now carry the
    # per-point colour, so overlaying dots would just double-draw the same information.
    empty_mk <- list(x = list(NA_real_), y = list(NA_real_))
    if (is3d) empty_mk$z <- list(NA_real_)
    plotlyProxyInvoke(proxy, "restyle", empty_mk, tm$markers_idx)
  } else {
    # ----- Markers (no trails active) -----
    if (length(fr$markers_x)) {
      mk_data <- list(
        x = list(fr$markers_x),
        y = list(fr$markers_y),
        `marker.size`  = if (is.finite(pt_size)) pt_size else 6,
        `marker.color` = list(as.character(fr$marker_cols))
      )
      if (is3d) mk_data$z <- list(fr$markers_z %||% numeric(0))
      plotlyProxyInvoke(proxy, "restyle", mk_data, tm$markers_idx)
    } else {
      empty_data <- list(x = list(NA_real_), y = list(NA_real_),
                         `marker.size` = if (is.finite(pt_size)) pt_size else 6)
      if (is3d) empty_data$z <- list(NA_real_)
      plotlyProxyInvoke(proxy, "restyle", empty_data, tm$markers_idx)
    }
  }

  # ----- NO RELAYOUT HERE, EVER -----
  # draw_frame is RESTYLE-ONLY in every path. Any relayout - even a title-only one -
  # makes plotly re-solve the scaleanchor aspect constraint and visibly rescales the
  # FOV (the sequential loop skips relayout for the same reason; see
  # run_sequential_loop_plain). Gating it on "not playing" was not enough: draw_frame
  # is also driven by the slider observer and by stale slider echoes that arrive after
  # playback ends (slider_guard is reset to 0 there), so the flags can read "paused"
  # while frames are still advancing - which produced a per-frame FOV jitter after
  # stop-then-replay. The frame/marker readout now lives in trk_frame_meta (driven by
  # tick()), which costs the plot nothing. uirevision keeps the view - and any user
  # zoom - stable across the restyles.

  # ----- Slider sync (throttled to every 5 frames during playback) -----
  if (!is.null(isolate(input$trk_frame_idx))) {
    if (!isTRUE(loop$gplaying) && !isTRUE(loop$playing) || (k %% 5L) == 1L) {
      slider_guard(isolate(slider_guard()) + 1L)
      updateSliderInput(session, "trk_frame_idx", value = k, max = length(pre$frames))
    }
  }
  
  invisible(k)
}


# Draw one step of the CURRENT TID by point index (1..nrow(tb))
draw_tid_point <- function(tid, i) {
  # Read reactives safely (no dependencies)
  pre <- isolate(frames_store()); if (is.null(pre)) return(invisible(NULL))
  tm  <- isolate(trace_map());    if (is.null(tm))  return(invisible(NULL))
  sd  <- isolate(seq_data());     if (is.null(sd) || !(tid %in% sd$tids)) return(invisible(NULL))
  
  tb <- sd$sp[[tid]]; if (is.null(tb) || !nrow(tb)) return(invisible(NULL))
  n  <- nrow(tb)
  i  <- max(1L, min(as.integer(i), n))
  
  tnow  <- tb$t[i]
  is3d  <- isTRUE(pre$has_3d)
  has_z_col <- is3d && "z" %in% names(tb)
  w <- which(tb$t >= (tnow - pre$trail_len) & tb$t <= tnow)
  if (length(w) > pre$trail_nmax) w <- tail(w, pre$trail_nmax)
  xx <- if (length(w)) c(tb$x[w], NA_real_) else numeric(0)
  yy <- if (length(w)) c(tb$y[w], NA_real_) else numeric(0)
  zz <- if (has_z_col && length(w)) c(tb$z[w], NA_real_) else numeric(0)

  col <- as.character(pre$col_map[[tid]])
  proxy <- plotlyProxy("trk_track2D", session)

  # Trail (single trace in sequential mode) - only when the trail trace exists
  # (i.e. "Animate moving trails" is on). Otherwise it's a plain moving point.
  if (!is.na(tm$trail_start)) {
    trail_data <- list(x = list(xx), y = list(yy), `line.color` = list(col))
    if (is3d) trail_data$z <- list(zz)
    plotlyProxyInvoke(proxy, "restyle", trail_data, tm$trail_start)
  }

  # Marker
  mk_data <- list(x = list(tb$x[i]), y = list(tb$y[i]),
                  `marker.size`  = isolate(input$trk_pt_size),
                  `marker.color` = list(col))
  if (has_z_col) mk_data$z <- list(tb$z[i])
  plotlyProxyInvoke(proxy, "restyle", mk_data, tm$markers_idx)
  
  plotlyProxyInvoke(proxy, "relayout",
                    list(title = paste0("tid = ", tid,
                                        " (", which(pre$tids == tid), "/", length(pre$tids), ")",
                                        ", idx = ", i,
                                        ", t = ", signif(tnow, 6))))
  
  # Keep the slider in sync (ok to use isolate here)
  if (!is.null(isolate(input$trk_frame_idx))) {
    slider_guard(isolate(slider_guard()) + 1L)
    updateSliderInput(session, "trk_frame_idx", value = i, max = n)
  }
}


# ── heatmap_store() - dwell-time heatmap ────────────────────────────────────
# GOAL: a 2D map of how much TIME molecules spent in each region of the field,
# rather than how many localizations landed there. Time-weighting matters
# because a slow/trapped molecule produces many localizations in one spot; a
# raw count map would exaggerate it purely as a sampling artefact.
# Bin count is user-set, or auto-scaled from the dataset size.
heatmap_store <- eventReactive(plot_trigger(), {
  df <- isolate(trk()); req(nrow(df) > 0)
  # Auto bin count keeps the map readable across very different dataset sizes.
  bins <- if (isTruthy(isolate(input$trk_hm_bins))) isolate(input$trk_hm_bins) else auto_hm_bins(nrow(df))
  build_dwell_heatmap(df[, c("x","y","tid","tim")], nbins = bins)   # needs tim for dwell weighting
}, ignoreInit = TRUE)
