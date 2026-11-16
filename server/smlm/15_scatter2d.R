## ─────────────────────────────────────────────────────────────────────────────
## server/smlm/15_scatter2d.R  -  main 2D plot
## ─────────────────────────────────────────────────────────────────────────────
## Part of smlm_server(). Not a standalone script: its statements are evaluated in
## the server function's environment, in file-name order (server/parts.R).
##
## output$scatter2D, the live scale bar and the 2D crop-to-view.
##
## INTERFACE (generated: Rscript tools/server_map.R --headers; do not edit by hand)
##   DEFINES for other parts:
##     (nothing)
##   READS from other parts:
##     01_view_upload: smlm_2d_crop_rv, smlm_axis_x_span,
##                     smlm_default_x_span, smlm_push_scale_bar,
##                     smlm_zoom_ranges
##     02_multiplex: .mplex_clu_memo
##     03_dataset_zcal: uploaded_df
##     10_triggers_stores: plot_trigger, smlm_diag
##     13_image_overlay: image_overlay
##     14_gt_estimator: add_gt_overlay_2d
##     17_plottable_data: plottable_data
##     18_nn_tables: secondary_nn_lines
## END INTERFACE
## ─────────────────────────────────────────────────────────────────────────────

# What the 2D view last showed: the multiplex point-set key and view plane of
# the last render, and the view plane a zoom was made in. Used to keep the
# zoom across a re-render of the same data (see the layout step below).
.smlm2d_view <- new.env()
observeEvent(smlm_zoom_ranges(), {
  .smlm2d_view$zoom_vp <- .smlm2d_view$vp %||% (input$view_plane_2d %||% "XY")
}, ignoreNULL = FALSE)

output$scatter2D <- renderPlotly({
  req(input$mode3d == "2D")
  if (is.null(input$file) || is.null(plot_trigger())) {
    return(plotly::plot_ly(source = "smlm2d", type = "scatter", mode = "markers",
                           x = numeric(0), y = numeric(0)) %>%
             plotly::event_register("plotly_relayout"))
  }

  df <- plottable_data()
  shiny::validate(shiny::need(nrow(df) > 0, "No points to display - please check filters or upload a valid file."))

  # Remap x/y columns for XZ or YZ viewing planes
  vp <- input$view_plane_2d %||% "XY"
  axis_labels <- list(x = "X (nm)", y = "Y (nm)")
  if (vp == "XZ" && "z" %in% names(df)) {
    df <- df %>% dplyr::mutate(y = z)
    axis_labels <- list(x = "X (nm)", y = "Z (nm)")
  } else if (vp == "YZ" && "z" %in% names(df)) {
    df <- df %>% dplyr::mutate(x = y, y = z)
    axis_labels <- list(x = "Y (nm)", y = "Z (nm)")
  }

  df_raw <- uploaded_df()
  # Ensure x/y exist on raw
  df_raw <- df_raw %>% dplyr::mutate(x = .data[["x"]], y = .data[["y"]])

  # Base limits from the data
  lim <- get_equal_limits_2d(df$x, df$y)
  
  # Try to build an overlay (may be NULL)
  ovl <- image_overlay()
  
  # Expand limits to include image bbox
  lim <- expand_limits_with_image(lim, ovl, invert_y = isTRUE(input$img_invert_y))
  
  # -------------------------
  # 0) Start empty figure
  # -------------------------
  fig    <- plotly::plot_ly(source = "smlm2d")
  fg_col <- if (isTRUE(input$dark_mode)) "#FFFFFF" else "#63666A"

  # -------------------------
  # 1) HEATMAP (below points)
  # -------------------------
  if (isTRUE(input$show_heatmap)) {
    df_used <- heatmap_input_df(
      df,
      use_raw    = isTRUE(input$heatmap_use_raw),
      clustering = isTRUE(input$enable_clustering)
    )
    lim_hm <- get_equal_limits_2d(df_used$x, df_used$y)
    ovl_for_hm <- image_overlay()
    lim_hm <- expand_limits_with_image(lim_hm, ovl_for_hm, invert_y = isTRUE(input$img_invert_y))
    
    est <- compute_kde(
      df_used,
      lim = lim_hm,
      n   = if (is.null(input$heatmap_n)) 200 else input$heatmap_n,
      bw  = if (is.null(input$heatmap_bw_abs)) NA_real_ else input$heatmap_bw_abs,
      bw_factor = if (is.null(input$heatmap_bw_factor)) 1 else input$heatmap_bw_factor
    )
    
    if (!is.null(est)) {
      dens <- est$z
      floor_p <- if (is.null(input$heatmap_floor)) 0 else input$heatmap_floor
      if (is.finite(floor_p) && floor_p > 0) {
        thr <- stats::quantile(dens, probs = floor_p, na.rm = TRUE)
        dens[dens < thr] <- 0
        if (max(dens, na.rm = TRUE) > 0) dens <- dens / max(dens, na.rm = TRUE)
      }
      dens[!is.finite(dens)] <- 0
      
      clip_p <- if (is.null(input$heatmap_clip)) 0.99 else input$heatmap_clip
      cap    <- stats::quantile(dens, probs = clip_p, na.rm = TRUE)
      dens   <- pmin(dens, cap)
      
      if (max(dens, na.rm = TRUE) > 0) dens <- dens / max(dens, na.rm = TRUE)
      
      if (isTRUE(input$heatmap_log)) {
        dens <- log1p(dens)
        dens <- dens / max(dens, na.rm = TRUE)
      }
      
      gm <- if (is.null(input$heatmap_gamma)) 1 else input$heatmap_gamma
      if (is.finite(gm) && gm > 0) {
        dens <- dens ^ gm
        if (max(dens, na.rm = TRUE) > 0) dens <- dens / max(dens, na.rm = TRUE)
      }
      
      fig <- fig %>% plotly::add_heatmap(
        x = est$x,
        y = est$y,
        z = t(dens),
        zmin = 0, zmax = 1,
        colors = palette_colors(input$heatmap_palette),
        showscale = isTRUE(input$heatmap_showlegend),
        colorbar = list(
          titlefont = list(color = fg_col),
          tickfont  = list(color = fg_col)
        ),
        opacity = if (is.null(input$heatmap_opacity)) 0.35 else input$heatmap_opacity,
        hoverinfo = "skip",
        inherit = FALSE
      )
      
      # lock ranges to heatmap grid
      lim <- lim_hm
    }
  }
  
  # -------------------------
  # 2) BACKGROUND IMAGE (below points)
  # -------------------------
  if (!is.null(ovl)) {
    fig <- fig %>% layout(images = list(ovl))
  }
  
  # -----------------------------------
  # 3) POINTS (build pts + marker once)
  # -----------------------------------
  # --- build the points once (same idea as 3D) ---
  disp <- smlm_display_from_input(input)
  pts <- build_pts_2d(df, disp, deblinked = isTRUE(input$enable_clustering))
  pts$color <- build_colors_for_pts(pts, disp, default_palette_col = "color")
  cmap <- color_mapping_for_pts(pts, disp)
  
  # bail out early if empty; still keep overlay + axes
  if (!nrow(pts)) {
    if (!is.null(ovl)) fig <- fig %>% layout(images = list(ovl))
    return(
      fig %>% layout(
        title = "2D Scatter Plot",
        xaxis = list(title = axis_labels$x, range = lim$x, scaleanchor = "y", scaleratio = 1),
        yaxis = list(title = axis_labels$y, range = lim$y),
        paper_bgcolor = if (isTRUE(input$dark_mode)) "#1e1e1e" else "white",
        plot_bgcolor  = if (isTRUE(input$dark_mode)) "#1e1e1e" else "white"
      )
    )
  }
  
  # --- draw the points once ---
  mplex_draw <- identical(input$color_mode %||% "default", "multiplex") &&
                "mplex_label" %in% names(pts) && any(!is.na(pts$mplex_label))
  if (mplex_draw) {
    # one trace per channel so each gets a named, colored legend entry
    ptso  <- pts[order(pts$mplex_channel), , drop = FALSE]
    chans <- unique(ptso$mplex_channel[!is.na(ptso$mplex_channel)])
    for (chv in chans) {
      sub <- ptso[!is.na(ptso$mplex_channel) & ptso$mplex_channel == chv, , drop = FALSE]
      if (!nrow(sub)) next
      # WebGL: a multiplex is usually many channels of full acquisitions, and
      # SVG redraws of tens of thousands of markers are what makes each
      # re-render (e.g. an alignment nudge) feel slow. Plain circle markers
      # only, which scattergl supports.
      fig <- fig %>% plotly::add_trace(
        data = sub, x = ~x, y = ~y,
        type = "scattergl", mode = "markers",
        name = sub$mplex_label[1],
        marker = list(size = input$median_point_size %||% 4, color = sub$mplex_color[1]),
        text = ~hover_text,
        hovertemplate = "%{text}<extra></extra>",
        hoverinfo = "text",
        showlegend = TRUE, inherit = FALSE
      )
    }
  } else if (isTRUE(cmap$showscale) && identical(cmap$kind, "scale")) {
    fig <- fig %>% plotly::add_trace(
      data = pts,
      x = ~x, y = ~y,
      type = "scatter", mode = "markers",
      name = if (isTRUE(input$show_medians)) "Localization Medians" else "Localizations",
      marker = list(
        size = input$median_point_size %||% 4,
        color = cmap$colorvar,
        colorscale = cmap$colorscale,
        cmin = cmap$cmin, cmax = cmap$cmax,
        showscale = TRUE,
        colorbar = list(
          title     = cmap$title,
          titlefont = list(color = fg_col),
          tickfont  = list(color = fg_col)
        )
      ),
      text = ~hover_text,
      hovertemplate = "%{text}<extra></extra>",
      hoverinfo = "text",
      showlegend = TRUE, inherit = FALSE
    )
  } else {
    hex <- if (!is.null(cmap$hex)) cmap$hex else rep("#bdbdbd", nrow(pts))
    fig <- fig %>% plotly::add_trace(
      data = pts,
      x = ~x, y = ~y,
      type = "scatter", mode = "markers",
      name = if (isTRUE(input$show_medians)) "Localization Medians" else "Localizations",
      marker = list(size = input$median_point_size %||% 4, color = hex),
      text = ~hover_text,
      hovertemplate = "%{text}<extra></extra>",
      hoverinfo = "text",
      showlegend = TRUE, inherit = FALSE
    )
  }
  
  # ------------------------------------
  # 5) CLUSTER2 NN / FITS (unchanged)
  # ------------------------------------
  if (isTRUE(input$show_cluster2_nn_lines) && isTRUE(input$enable_cluster2) && !is.null(df$cluster2)) {
    lines <- secondary_nn_lines()
    if (length(lines)) {
      seg_x <- c(); seg_y <- c(); lab_x <- c(); lab_y <- c(); labs <- c()
      for (ln in lines) {
        seg_x <- c(seg_x, ln$x, NA_real_)
        seg_y <- c(seg_y, ln$y, NA_real_)
        lab_x <- c(lab_x, ln$midpoint[1]); lab_y <- c(lab_y, ln$midpoint[2]); labs <- c(labs, ln$text)
      }
      fig <- fig %>% plotly::add_trace(
        x = seg_x, y = seg_y, type = "scatter", mode = "lines",
        line = list(color = "green", width = 2),
        hoverinfo = "none", showlegend = FALSE, inherit = FALSE
      )
      if (length(labs) && isTRUE(input$show_c2_nn_labels)) {
        nn_lfs <- suppressWarnings(as.numeric(input$c2_label_font_size %||% 12))
        if (!is.finite(nn_lfs) || nn_lfs < 1) nn_lfs <- 12
        fig <- fig %>% plotly::add_trace(
          x = lab_x, y = lab_y, type = "scatter", mode = "text",
          text = labs, textfont = list(color = "green", size = nn_lfs),
          hoverinfo = "none", showlegend = FALSE, inherit = FALSE
        )
      }
    }
  }


  # ------------------------------------
  # 6) GENERAL NN + 2-color NN (your blocks)
  # ------------------------------------
  if (isTRUE(input$show_nn_lines) && nrow(df) > 1) {
    nn <- build_nn_segments(df, smlm_nn_options(input), is3d = FALSE)
    if (!is.null(nn)) {
      ln_col <- if (isTRUE(input$dark_mode)) "white" else "black"
      # With channels, colour each line by the channels it joins: a
      # within-channel line takes that channel's colour, a cross-channel line
      # is magenta. Without channels it is one plain series, as before.
      ed <- nn$edges
      if (!is.null(ed) && nrow(ed) && length(unique(stats::na.omit(ed$ch1))) +
          length(unique(stats::na.omit(ed$ch2))) > 2L) {
        for (cc in unique(ed$color)) {
          sel <- which(ed$color == cc)
          fig <- fig %>% plotly::add_trace(
            x = as.numeric(t(cbind(ed$x1[sel], ed$x2[sel], NA))),
            y = as.numeric(t(cbind(ed$y1[sel], ed$y2[sel], NA))),
            type = "scatter", mode = "lines",
            line = list(color = cc, width = 2),
            hoverinfo = "none", showlegend = FALSE, inherit = FALSE
          )
        }
      } else {
      fig <- fig %>% plotly::add_trace(
        x = nn$seg_x, y = nn$seg_y,
        type = "scatter", mode = "lines",
        line = list(color = ln_col, width = 1),
        hoverinfo = "none", showlegend = FALSE, inherit = FALSE
      )
      }
      if (isTRUE(input$show_nn_labels) && length(nn$labs)) {
        lab_col <- if (isTRUE(input$dark_mode)) "orange" else "red"
        nn_lfs <- suppressWarnings(as.numeric(input$nn_label_font_size %||% 12))
        if (!is.finite(nn_lfs) || nn_lfs < 1) nn_lfs <- 12
        fig <- fig %>% plotly::add_trace(
          x = nn$lab_x, y = nn$lab_y,
          type = "scatter", mode = "text",
          text = nn$labs, textfont = list(color = lab_col, size = nn_lfs),
          hoverinfo = "none", showlegend = FALSE, inherit = FALSE
        )
      }
    }
  }
  
  # The separate 2-colour NN block that used to live here is gone. NN lines are
  # drawn once, by the shared panel above, which now colours each line by the
  # channels it joins - so enabling channels ADAPTS the Nearest Neighbour
  # controls rather than adding a parallel set of them.

  # ------------------------------------
  # 7) AXES + THEMING (equal aspect)
  # ------------------------------------
  bg_col   <- if (isTRUE(input$dark_mode)) "#000000" else "#FFFFFF"
  fg_col   <- if (isTRUE(input$dark_mode)) "#FFFFFF" else "#63666A"
  show_g   <- !isTRUE(input$hide_grid)
  grid_col <- if (isTRUE(input$dark_mode)) "rgba(255,255,255,0.18)" else "rgba(0,0,0,0.18)"
  zero_col <- if (isTRUE(input$dark_mode)) "rgba(255,255,255,0.25)" else "rgba(0,0,0,0.25)"
  
  # Multiplex alignment: a re-render of the SAME point set (same rows, same
  # settings; only a channel shift or a display option changed) keeps the
  # user's zoom, so one can zoom onto a structure and nudge. Any other change
  # of the point set, or a different view plane, falls back to the full extent.
  rows_key <- .mplex_clu_memo$key
  zr <- isolate(smlm_zoom_ranges())
  if (!is.null(rows_key) && identical(rows_key, .smlm2d_view$rows_key) &&
      identical(vp, .smlm2d_view$zoom_vp) && !is.null(zr))
    lim <- list(x = zr$x, y = zr$y)
  .smlm2d_view$rows_key <- rows_key
  .smlm2d_view$vp       <- vp

  fig <- fig %>%
    plotly::layout(
      title = "2D Scatter Plot",
      paper_bgcolor = bg_col, plot_bgcolor = bg_col,
      xaxis = list(
        title = axis_labels$x,
        range = lim$x,
        showgrid = show_g, gridcolor = grid_col, gridwidth = 0.4,
        zeroline = show_g, zerolinecolor = zero_col, zerolinewidth = 0.6,
        showticklabels = show_g, tickfont = list(color = fg_col),
        color = fg_col
      ),
      yaxis = list(
        title = axis_labels$y,
        range = lim$y,
        showgrid = show_g, gridcolor = grid_col, gridwidth = 0.4,
        zeroline = show_g, zerolinecolor = zero_col, zerolinewidth = 0.6,
        showticklabels = show_g, tickfont = list(color = fg_col),
        color = fg_col,
        scaleanchor = "x", scaleratio = 1
      ),
      showlegend = TRUE,
      font = list(color = fg_col),
      legend = list(
        font = list(color = fg_col),
        orientation = "h",
        x = 0, y = 0.97,
        xanchor = "left", yanchor = "bottom"
      )
    )

  smlm_axis_x_span(diff(lim$x))

  # --- 2D Cluster2 convex hulls (optional) ---
  if (!is.null(input$show_cluster2_shape_2D) &&
      isTRUE(input$show_cluster2_shape_2D) &&
      isTRUE(input$enable_cluster2) &&
      !is.null(df$cluster2)) {
    
    hull_opacity <- if (!is.null(input$hull_opacity)) input$hull_opacity else 0.2
    
    for (cl2 in unique(df$cluster2[df$cluster2 != 0])) {
      group_df <- df %>% dplyr::filter(cluster2 == cl2)
      if (nrow(group_df) >= 3) {
        hidx <- tryCatch(grDevices::chull(group_df$x, group_df$y), error = function(e) integer(0))
        if (length(hidx) >= 3) {
          hull_xy <- group_df[c(hidx, hidx[1]), c("x","y")]
          base_col <- if ("cluster2_color" %in% names(group_df) && !is.na(group_df$cluster2_color[1]))
            group_df$cluster2_color[1] else "#00BFFF"
          fill_col <- grDevices::adjustcolor(base_col, alpha.f = hull_opacity)
          
          fig <- fig %>%
            plotly::add_trace(
              x = hull_xy$x, y = hull_xy$y,
              type = "scatter", mode = "lines",
              fill = "toself",
              fillcolor = fill_col,
              line = list(color = base_col, width = 2),
              name = paste("Median Cluster", cl2),
              showlegend = TRUE,
              inherit = FALSE
            )
        }
      }
    }
  }

  # ---- Scale bar ----
  if (isTRUE(input$show_scale_bar) && input$mode3d == "2D") {
    sb_len   <- suppressWarnings(as.numeric(input$scale_bar_length))
    sb_thick <- suppressWarnings(as.numeric(input$scale_bar_thickness))
    sb_font  <- suppressWarnings(as.numeric(input$scale_font_size))
    sb <- make_scale_bar(sb_len, diff(lim$x), sb_thick, dark = isTRUE(input$dark_mode),
                         scale_font_size = if (is.finite(sb_font) && sb_font > 0) sb_font else 18)
    if (!is.null(sb))
      fig <- fig %>% plotly::layout(shapes = sb$shapes,
                                    annotations = sb$annotations)
  }

  # ---- Diagnostic: detected cluster centres (red +) and confocal blobs (green ○) ----
  # type="scatter" (SVG) - scattergl lacks circle-open / cross symbol support.
  diag <- smlm_diag()   # non-isolated: toggling show_diag or new alignment re-renders
  if (isTRUE(input$show_diag) && !is.null(diag)) {
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

  # ---- Estimated ground-truth candidates (hollow circles, green->red) ----
  fig <- add_gt_overlay_2d(fig, vp)

  fig %>% plotly::event_register('plotly_relayout')
})

observe({
  # Only read the relayout event while the SMLM tab is on screen, so plotly's
  # deferred "source not registered" warning is never scheduled.
  if (!identical(input$main_navbar, "tab_smlm")) return()
  req(input$file, plot_trigger(), identical(input$mode3d, "2D"))
  ed <- suppressWarnings(plotly::event_data("plotly_relayout", source = "smlm2d"))
  if (!is.list(ed)) return()

  x0 <- suppressWarnings(as.numeric(ed[["xaxis.range[0]"]]))
  x1 <- suppressWarnings(as.numeric(ed[["xaxis.range[1]"]]))
  if (length(x0) == 1 && length(x1) == 1 && is.finite(x0) && is.finite(x1) && x1 > x0) {
    smlm_axis_x_span(x1 - x0)
    y0 <- suppressWarnings(as.numeric(ed[["yaxis.range[0]"]]))
    y1 <- suppressWarnings(as.numeric(ed[["yaxis.range[1]"]]))
    if (length(y0) == 1 && length(y1) == 1 && is.finite(y0) && is.finite(y1) && y1 > y0)
      smlm_zoom_ranges(list(x = c(x0, x1), y = c(y0, y1)))
  } else if (isTRUE(ed[["xaxis.autorange"]])) {
    span <- smlm_default_x_span()
    if (length(span) == 1 && is.finite(span) && span > 0) smlm_axis_x_span(span)
    smlm_zoom_ranges(NULL)
  } else {
    return()
  }

  if (isTRUE(input$show_scale_bar)) smlm_push_scale_bar()
})

observeEvent(list(input$show_scale_bar, input$scale_bar_length,
                  input$scale_bar_thickness, input$scale_font_size,
                  input$dark_mode, input$mode3d), {
  smlm_push_scale_bar()
}, ignoreInit = TRUE)

# ---- 2D crop-to-view → feeds the existing box trimmer -----------------------
# When the user zooms in 2D and presses "Crop to view", we enable the box
# trimmer and push the current axis ranges into trim_x_range / trim_y_range.
# Because those inputs are created by renderUI (trim_box_ui) which only renders
# after enable_trimmer=TRUE + trim_mode="box", we use a short later() delay to
# let the DOM catch up before sending the update.  The crop then persists
# automatically when the user switches to 3D.
observeEvent(input$smlm_crop_to_view, {
  zr <- smlm_zoom_ranges()
  if (is.null(zr)) {
    showNotification("Zoom into a region first, then press Crop to view.",
                     type = "warning", duration = 5)
    return()
  }
  smlm_2d_crop_rv(zr)
  plot_trigger(Sys.time())
  showNotification(
    sprintf("Crop active: x [%.0f, %.0f] nm  |  y [%.0f, %.0f] nm (retained in 3D).",
            zr$x[1], zr$x[2], zr$y[1], zr$y[2]),
    type = "message", duration = 6
  )
}, ignoreInit = TRUE)

observeEvent(input$smlm_clear_crop, {
  smlm_2d_crop_rv(NULL)
  plot_trigger(Sys.time())
  showNotification("Spatial crop cleared.", type = "message", duration = 3)
}, ignoreInit = TRUE)

spatial_info_ui <- reactive({
  crop <- smlm_2d_crop_rv()
  if (is.null(crop))
    tags$p(style = "font-size:0.82em; color:#888; margin:2px 0;",
           "No 2D crop active.")
  else
    tags$p(style = "font-size:0.82em; color:#2070a0; margin:2px 0;",
           sprintf("Active: x [%.0f, %.0f]  |  y [%.0f, %.0f] nm",
                   crop$x[1], crop$x[2], crop$y[1], crop$y[2]))
})
output$smlm_spatial_info    <- renderUI({ spatial_info_ui() })
output$smlm_spatial_info_3d <- renderUI({ spatial_info_ui() })

observeEvent(input$smlm_clear_crop_3d, {
  smlm_2d_crop_rv(NULL)
  plot_trigger(Sys.time())
  showNotification("Spatial crop cleared.", type = "message", duration = 3)
}, ignoreInit = TRUE)
