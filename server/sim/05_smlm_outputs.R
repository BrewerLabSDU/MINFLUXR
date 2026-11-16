## ─────────────────────────────────────────────────────────────────────────────
## server/sim/05_smlm_outputs.R  -  SMLM simulator outputs
## ─────────────────────────────────────────────────────────────────────────────
## Part of sim_server(). Not a standalone script: its statements are evaluated in
## the server function's environment, in file-name order (server/parts.R).
##
## The plot, summaries, precision budget, comparison views, confocal view,
## acquisition plots and CSV exports.
##
## INTERFACE (generated: Rscript tools/server_map.R --headers; do not edit by hand)
##   DEFINES for other parts:
##     crlb_inputs
##   READS from other parts:
##     01_view_scalebars: sim_axis_x_span, sim_diff_full, sim_diff_ready,
##                        sim_diff_x_span, sim_full_x_span, sim_plot_ready,
##                        sim_push_diff_scale_bar, sim_push_scale_bar,
##                        sim_use_fov
##     02_smlm_simulation: sim_acq_live, sim_result
##     03_real_run_calibration: sim_real_rv
##     04_compare_strategies: compare_result
## END INTERFACE
## ─────────────────────────────────────────────────────────────────────────────

# ── Phase 1 plot ─────────────────────────────────────────────────────────────
output$sim_plot <- renderPlotly({
  res <- sim_result()
  shiny::validate(
    shiny::need(!is.null(res), "Click 'Generate' to run the simulation."))

  dark    <- isTRUE(input$sim_dark_mode)
  bg_col  <- if (dark) "#111111" else "#FFFFFF"
  fg_col  <- if (dark) "#DDDDDD" else "#444444"
  use_3d  <- isTRUE(res$use_3d)
  tr_type <- if (use_3d) "scatter3d" else "scattergl"
  pt_size <- max(1L, as.integer(input$sim_point_size))
  gc      <- if (dark) "rgba(255,255,255,0.1)" else "rgba(0,0,0,0.1)"
  show_grid <- !isTRUE(input$sim_hide_grid)

  # Field of view - display window only, nothing is cropped. Axes are pinned to
  # a real-camera-style FOV; the 3D aspect ratio reflects the true X:Y:Z
  # proportions so Z appears compressed (4000×4000×1000). Always on in
  # Acquisition mode; in Normal mode it's opt-in via "Pin axes to a FOV (ROI)".
  use_fov <- sim_use_fov()
  fov_x   <- max(1, as.numeric(input$sim_fov_x %||% 4000))
  fov_y   <- max(1, as.numeric(input$sim_fov_y %||% 4000))
  fov_z   <- max(1, as.numeric(input$sim_fov_z %||% 1000))

  plt <- plot_ly(source = "sim2d")

  # Helper: add a plain marker layer
  add_layer <- function(p, df, symbol, size, color, name, htxt = NULL) {
    if (nrow(df) == 0) return(p)
    z_val <- if (use_3d && "z" %in% names(df)) df$z else NULL
    args  <- list(
      x = df$x, y = df$y,
      type = tr_type, mode = "markers",
      marker = list(symbol = symbol, size = size, color = color,
                    line = list(color = color, width = 1)),
      name = name,
      hoverinfo = if (is.null(htxt)) "skip" else "text",
      text = htxt
    )
    if (!is.null(z_val)) args$z <- z_val
    do.call(add_trace, c(list(p), args))
  }

  # Build hover text helpers
  hover_xy <- function(df, prefix = "") {
    z_part <- if (use_3d && "z" %in% names(df))
                paste0("<br>z: ", round(df$z, 1), " nm")
              else rep("", nrow(df))
    paste0(prefix,
           "x: ", round(df$x, 1), " nm",
           "<br>y: ", round(df$y, 1), " nm",
           z_part)
  }

  # GT layer - bright distinct markers so they stand out in both 2D and 3D.
  # In image mode these are the intensity-weighted sampled molecule positions.
  if (isTRUE(input$sim_show_gt) && nrow(res$gt) > 0) {
    gt_clean <- res$gt[, intersect(c("x","y","z"), names(res$gt)), drop = FALSE]
    gt_label <- if (identical(input$sim_shape, "image"))
                  "GT (sampled positions)" else "Ground truth"
    ht <- hover_xy(gt_clean, "GT<br>")
    plt <- add_layer(plt, gt_clean,
                     if (use_3d) "circle-open" else "cross",
                     round(pt_size * 2),
                     if (dark) "#FFEE44" else "#CC2200",
                     gt_label, ht)
  }


  # Localization layer - coloured by TID, Z, time, or local density.
  if (isTRUE(input$sim_show_locs) && nrow(res$locs) > 0) {
    cmode <- input$sim_color_mode %||% "tid"
    if (identical(cmode, "z") && !(use_3d && "z" %in% names(res$locs)))
      cmode <- "tid"                                  # no depth in 2D \u2192 TID
    z_val <- if (use_3d && "z" %in% names(res$locs)) res$locs$z else NULL
    z_ht  <- if (!is.null(z_val))
                paste0("<br>z: ", round(z_val, 1), " nm")
             else rep("", nrow(res$locs))
    ht <- paste0("TID: ", res$locs$tid,
                 "<br>x: ", round(res$locs$x, 1), " nm",
                 "<br>y: ", round(res$locs$y, 1), " nm",
                 z_ht)

    if (identical(cmode, "tid")) {
      pal    <- tid_color_map(unique(res$locs$tid))
      marker <- list(size = pt_size,
                     color = unname(pal[res$locs$tid]), opacity = 0.85)
    } else {
      cv <- switch(cmode,
        z       = list(v = res$locs$z,   title = "Depth z (nm)"),
        tim     = list(v = res$locs$tim, title = "Time (s)"),
        locdens = list(v = sim_local_density(res$locs,
                            max(1, as.numeric(input$sim_dens_radius_nm %||% 500))),
                       title = "Local density"))
      v   <- as.numeric(cv$v)
      rng <- range(v[is.finite(v)])
      if (!all(is.finite(rng))) rng <- c(0, 1)
      marker <- list(
        size = pt_size, color = v, opacity = 0.85,
        colorscale = plotly_colorscale(
          input$sim_numeric_palette %||% "Spectral",
          reverse = isTRUE(input$sim_palette_reverse)),
        cmin = rng[1], cmax = rng[2], showscale = TRUE,
        colorbar = list(title = cv$title,
                        titlefont = list(color = fg_col),
                        tickfont  = list(color = fg_col)))
    }
    args <- list(
      x = res$locs$x, y = res$locs$y,
      type = tr_type, mode = "markers",
      marker = marker,
      name = sprintf("Localizations (n\u202F=\u202F%d)", nrow(res$locs)),
      hoverinfo = "text", text = ht
    )
    if (!is.null(z_val)) args$z <- z_val
    plt <- do.call(add_trace, c(list(plt), args))
  }

  # Lines: GT molecule position → median of each TID.
  #
  # One line per TID, anchored at the true GT molecule the TID came from.
  # For single-dye labels one short line per molecule shows the small
  # link offset; for polyclonal + secondary many longer lines fan out
  # from the same GT point, reflecting the greater linker distance and
  # the number of fluorophores.  Works in all shape modes including image.
  if (isTRUE(input$sim_show_gt_lines) &&
      nrow(res$gt) > 0 && nrow(res$fluoro) > 0 && nrow(res$locs) > 0 &&
      "mol_id" %in% names(res$fluoro)) {
    mol_ids  <- res$fluoro$mol_id
    has_gz   <- use_3d && "z" %in% names(res$gt)
    has_lz   <- use_3d && "z" %in% names(res$locs)
    seg_col  <- if (dark) "rgba(255,120,80,0.75)" else "rgba(190,50,30,0.75)"
    n_fluoro <- nrow(res$fluoro)

    # One segment per TID: GT[mol_id] → median(locs for that TID)
    segs <- do.call(rbind, lapply(seq_len(n_fluoro), function(k) {
      lk <- res$locs[res$locs$tid == as.character(k), , drop = FALSE]
      if (nrow(lk) == 0) return(NULL)
      gi <- mol_ids[k]
      if (is.na(gi) || gi < 1L || gi > nrow(res$gt)) return(NULL)
      data.frame(
        x0 = res$gt$x[gi],  y0 = res$gt$y[gi],
        x1 = median(lk$x),  y1 = median(lk$y),
        z0 = if (has_gz) res$gt$z[gi]  else NA_real_,
        z1 = if (has_lz) median(lk$z)  else NA_real_
      )
    }))

    if (!is.null(segs) && nrow(segs) > 0) {
      if (use_3d) {
        xs <- c(rbind(segs$x0, segs$x1, NA))
        ys <- c(rbind(segs$y0, segs$y1, NA))
        zs <- c(rbind(segs$z0, segs$z1, NA))
        plt <- do.call(add_trace,
                       list(plt, x = xs, y = ys, z = zs,
                            type = "scatter3d", mode = "lines",
                            line = list(color = seg_col, width = 2),
                            name = "GT \u2192 TID median",
                            hoverinfo = "skip", showlegend = TRUE))
      } else {
        plt <- plt %>%
          add_segments(x = ~x0, xend = ~x1, y = ~y0, yend = ~y1,
                       data = segs, inherit = FALSE,
                       line = list(color = seg_col, width = 1.5),
                       name = "GT \u2192 TID median",
                       hoverinfo = "skip", showlegend = TRUE)
      }
    }
  }

  title_txt <- sprintf("Simulated localizations%s",
                       if (use_3d) " (3D)" else "")

  if (use_3d) {
    ref_pts <- if (nrow(res$locs) > 0) res$locs else res$gt
    if (use_fov) {
      # FOV pins each axis; aspect ratio = true X:Y:Z so Z is compressed.
      xr <- c(-fov_x / 2, fov_x / 2)
      yr <- c(-fov_y / 2, fov_y / 2)
      zr <- c(-fov_z / 2, fov_z / 2)
      aspect <- list(x = 1, y = fov_y / fov_x, z = fov_z / fov_x)
    } else {
      # Equal-scale axes: Z gets the same span as XY so the plot looks cubic.
      # Localizations remain at their actual Z positions near the centre.
      x_rng <- range(ref_pts$x, na.rm = TRUE)
      y_rng <- range(ref_pts$y, na.rm = TRUE)
      half  <- max(diff(x_rng), diff(y_rng)) / 2 * 1.05   # 5% padding
      if (half < 1) half <- 1
      x_ctr <- mean(x_rng); y_ctr <- mean(y_rng)
      z_ctr <- if ("z" %in% names(ref_pts))
                  median(ref_pts$z, na.rm = TRUE) else 0
      xr <- c(x_ctr - half, x_ctr + half)
      yr <- c(y_ctr - half, y_ctr + half)
      zr <- c(z_ctr - half, z_ctr + half)
      aspect <- list(x = 1, y = 1, z = 1)
    }

    plt <- plt %>% layout(
      paper_bgcolor = bg_col,
      scene = list(
        bgcolor = bg_col,
        xaxis = list(title = "X (nm)", color = fg_col, range = xr,
                     showgrid = show_grid, zeroline = show_grid),
        yaxis = list(title = "Y (nm)", color = fg_col, range = yr,
                     showgrid = show_grid, zeroline = show_grid),
        zaxis = list(title = "Z (nm)", color = fg_col, range = zr,
                     showgrid = show_grid, zeroline = show_grid),
        aspectmode  = "manual",
        aspectratio = aspect
      ),
      font = list(color = fg_col), showlegend = TRUE,
      legend = list(x = 0.01, y = 0.99, xanchor = "left", yanchor = "top",
                    bgcolor = if (dark) "rgba(17,17,17,0.6)"
                              else "rgba(255,255,255,0.6)",
                    font = list(size = 11)),
      title = list(text = title_txt, font = list(color = fg_col))
    )
  } else {
    img_layout <- if (!is.null(res$img_b64)) {
      list(list(source  = res$img_b64,
                xref    = "x", yref = "y",
                x       = -res$img_w_nm / 2,
                y       = res$img_h_nm / 2,
                sizex   = res$img_w_nm, sizey = res$img_h_nm,
                sizing  = "stretch", opacity = 0.25, layer = "below"))
    } else NULL

    xaxis <- list(title = "X (nm)", color = fg_col,
                  scaleanchor = "y", scaleratio = 1, gridcolor = gc,
                  showgrid = show_grid, zeroline = show_grid)
    yaxis <- list(title = "Y (nm)", color = fg_col, gridcolor = gc,
                  showgrid = show_grid, zeroline = show_grid)
    if (use_fov) {
      xaxis$range <- c(-fov_x / 2, fov_x / 2)
      yaxis$range <- c(-fov_y / 2, fov_y / 2)
    }

    # Scale bar (2D only) - reuses the shared SMLM make_scale_bar helper.
    sb <- NULL
    if (isTRUE(input$sim_show_scale_bar)) {
      x_span <- if (use_fov) fov_x else {
        ref <- if (nrow(res$locs) > 0) res$locs else res$gt
        diff(range(ref$x, na.rm = TRUE))
      }
      sb <- make_scale_bar(
        suppressWarnings(as.numeric(input$sim_scale_bar_length)),
        x_span,
        suppressWarnings(as.numeric(input$sim_scale_bar_thickness)),
        dark = dark,
        scale_font_size = suppressWarnings(
          as.numeric(input$sim_scale_font_size %||% 18)))
    }

    lay_args <- list(plt,
      paper_bgcolor = bg_col, plot_bgcolor = bg_col,
      font   = list(color = fg_col),
      xaxis  = xaxis, yaxis = yaxis,
      showlegend = TRUE,
      legend = list(x = 0.01, y = 0.99, xanchor = "left", yanchor = "top",
                    bgcolor = if (dark) "rgba(17,17,17,0.6)"
                              else "rgba(255,255,255,0.6)",
                    font = list(size = 11)),
      title  = list(text = title_txt, font = list(color = fg_col)),
      images = img_layout,
      uirevision = "sim"
    )
    if (!is.null(sb)) {
      lay_args$shapes      <- sb$shapes
      lay_args$annotations <- sb$annotations
    }
    plt <- do.call(layout, lay_args)
  }
  plt <- plt %>% config(displaylogo = FALSE,
                 toImageButtonOptions = list(format = "svg",
                                             filename = "sim_localisations"))
  # Register the relayout event so the "sim2d" source is recognised in both
  # 2D and 3D (avoids the "source not registered" console warning). Flag it as
  # rendered so the relayout observer only reads the event after registration.
  sim_plot_ready(TRUE)
  plt %>% plotly::event_register("plotly_relayout")
})

# Capture the visible x-span on zoom/pan and re-push the scale bar so it
# tracks the current view (matches the SMLM tab's behaviour).
observe({
  # Only read the relayout event when the plot has rendered + registered AND is
  # on screen (sim tab + MINFLUX subtab + 2D); otherwise plotly warns "source
  # not registered". sim_plot_ready stays FALSE until the first Generate render.
  if (!isTRUE(sim_plot_ready()) ||
      !identical(input$main_navbar, "tab_sim") ||
      !identical(input$sim_main_tabs, "MINFLUX") ||
      !identical(input$sim_view_mode, "2d")) return()
  ed <- suppressWarnings(plotly::event_data("plotly_relayout", source = "sim2d"))
  if (!is.list(ed)) return()
  x0 <- suppressWarnings(as.numeric(ed[["xaxis.range[0]"]]))
  x1 <- suppressWarnings(as.numeric(ed[["xaxis.range[1]"]]))
  if (length(x0) == 1 && length(x1) == 1 &&
      is.finite(x0) && is.finite(x1) && x1 > x0) {
    sim_axis_x_span(x1 - x0)
  } else if (isTRUE(ed[["xaxis.autorange"]])) {
    span <- sim_full_x_span()
    if (length(span) == 1 && is.finite(span) && span > 0) sim_axis_x_span(span)
  } else return()
  if (isTRUE(input$sim_show_scale_bar)) sim_push_scale_bar()
})

# Confocal view: capture its own zoom and re-push its (white) scale bar.
observe({
  if (!isTRUE(sim_diff_ready()) ||
      !identical(input$main_navbar, "tab_sim") ||
      !identical(input$sim_main_tabs, "Confocal view")) return()
  ed <- suppressWarnings(plotly::event_data("plotly_relayout", source = "sim_diff"))
  if (!is.list(ed)) return()
  x0 <- suppressWarnings(as.numeric(ed[["xaxis.range[0]"]]))
  x1 <- suppressWarnings(as.numeric(ed[["xaxis.range[1]"]]))
  if (length(x0) == 1 && length(x1) == 1 &&
      is.finite(x0) && is.finite(x1) && x1 > x0) {
    sim_diff_x_span(x1 - x0)
  } else if (isTRUE(ed[["xaxis.autorange"]])) {
    sim_diff_x_span(NULL)
  } else return()
  if (isTRUE(input$sim_show_scale_bar)) sim_push_diff_scale_bar()
})

observeEvent(list(input$sim_show_scale_bar, input$sim_scale_bar_length,
                  input$sim_scale_bar_thickness, input$sim_scale_font_size,
                  input$sim_dark_mode, input$sim_view_mode,
                  input$sim_mode, input$sim_view_fov, input$sim_fov_x), {
  sim_push_scale_bar()
  sim_push_diff_scale_bar()
}, ignoreInit = TRUE)

# ── Phase 1 summary ──────────────────────────────────────────────────────────
output$sim_summary <- renderUI({
  res <- sim_result()
  if (is.null(res)) return(NULL)
  nn_dist <- if (nrow(res$locs) >= 2) {
    idx  <- if (nrow(res$locs) > 500) sample(nrow(res$locs), 500)
            else seq_len(nrow(res$locs))
    cols <- intersect(c("x","y","z"), names(res$locs))
    pts  <- as.matrix(res$locs[idx, cols])
    d    <- as.matrix(dist(pts)); diag(d) <- Inf
    mean(apply(d, 1, min), na.rm = TRUE)
  } else NA_real_

  row <- function(lbl, val)
    tags$tr(tags$td(style = "padding:2px 8px;", lbl),
            tags$td(style = "padding:2px 8px;", tags$b(val)))

  # Steric warning block
  steric_ui <- if (!is.null(res$steric_warn)) {
    alert_cls <- if (identical(res$steric_warn$level, "danger"))
                   "alert alert-danger" else "alert alert-warning"
    tags$div(class = alert_cls,
             style = "padding:6px 10px; font-size:0.84em; margin-top:6px;",
             tags$b(if (identical(res$steric_warn$level, "danger"))
                      "\u26A0 Severe steric exclusion" else "\u26A0 Steric exclusion warning"),
             tags$br(), res$steric_warn$msg)
  } else NULL

  # Blinking TID inflation row (shown only when blinking is active)
  blink_row <- if (isTRUE(res$blink_active) && !is.null(res$n_fluoro_base)) {
    row("Blinking TID inflation",
        sprintf("%d base \u2192 %d apparent  (\u00D7%.1f, \u03BB\u202F=\u202F%.2f)",
                res$n_fluoro_base, res$n_fluoro,
                res$n_fluoro / max(1, res$n_fluoro_base),
                res$blink_lambda_used))
  } else NULL

  tagList(
    steric_ui,
    tags$div(
      class = "well well-sm", style = "font-size:0.88em; margin-top:10px;",
      tags$b("Simulation summary"), tags$br(),
      tags$table(style = "width:100%; border-collapse:collapse;",
        row("Molecules placed",  res$n_prot),
        row("Labeling strategy", res$lbl_note),
        row("Fluorophores (trace IDs)",
            sprintf("%d  (~%.1f/molecule)", res$n_fluoro,
                    res$n_fluoro / max(1, res$n_prot))),
        blink_row,
        row("Linker distance",   sprintf("~%.0f nm", res$link_nm)),
        if (!is.na(res$footprint_nm))
          row("Label footprint (steric check)", sprintf("~%.0f nm", res$footprint_nm)),
        row("Total localizations", res$n_locs),
        row("Avg locs / trace ID",
            sprintf("%.1f", res$n_locs / max(1, res$n_fluoro))),
        row("Precision XY \u03C3",    sprintf("%.1f nm", res$sigma_loc_xy)),
        if (isTRUE(res$use_3d))
          row("Precision Z \u03C3",   sprintf("%.1f nm", res$sigma_loc_z)),
        if (is.finite(nn_dist))
          row("Mean nearest-neighbour", sprintf("%.1f nm", nn_dist)),
        if (!is.null(res$rmse_vals) && length(res$rmse_vals) >= 1)
          row("Mean TID displacement error",
              sprintf("%.1f nm (SD %.1f, n\u202F=\u202F%d)",
                      mean(res$rmse_vals), sd(res$rmse_vals),
                      length(res$rmse_vals))),
        if (!is.null(res$rmse_vals) && length(res$rmse_vals) >= 1)
          row("Theoretical combined \u03C3  (\u221A(link\u00B2 + \u03C3\u00B2_xy))",
              sprintf("%.1f nm", sqrt(res$link_nm^2 + res$sigma_loc_xy^2)))
      )
    )
  )
})

# ── Labeling strategy info card (live, no Generate needed) ─────────────────
output$sim_label_info <- renderUI({
  key <- input$sim_label_type %||% "Custom"
  if (identical(key, "Custom")) return(NULL)
  p <- LABEL_PRESETS[[key]]
  if (is.null(p) || is.na(p$n_fluor)) return(NULL)
  tags$div(
    class = "alert alert-info",
    style = "padding:6px 10px; font-size:0.84em; margin-top:4px;",
    tags$b(sprintf("%d fluorophore(s) per molecule", p$n_fluor)),
    sprintf(" - linker ~%d nm", round(p$link_nm)),
    tags$br(), p$note
  )
})

# ── CRLB precision inputs ────────────────────────────────────────────────────
# GOAL: one place that reads the CRLB-mode inputs, so Generate and the live
# preview can never disagree. Fallbacks match the sim_ui.R defaults.
crlb_inputs <- reactive({
  # Cleared numericInputs arrive as NA; fall back to the default instead.
  num <- function(x, d) {
    x <- suppressWarnings(as.numeric(x))
    if (length(x) != 1 || !is.finite(x)) d else x
  }
  L_nm   <- max(1, num(input$sim_crlb_L,       40))    # physical L
  n_ph   <- max(1, num(input$sim_crlb_photons, 80))    # total photons
  bg     <- max(0, min(0.99, num(input$sim_crlb_bg, 0.05)))
  k      <- max(0, num(input$sim_crlb_k,       1))     # × bound
  floor0 <- max(0, num(input$sim_crlb_floor,   6.5))   # nm, quadrature
  s <- smlm_sigma_loc(L_nm, n_ph, bg, k, floor0)
  c(s, list(floor = floor0))
})

# ── CRLB precision preview (live - no Generate needed) ──────────────────────
output$sim_crlb_preview <- renderUI({
  req(identical(input$sim_prec_mode %||% "manual", "crlb"))
  s       <- crlb_inputs()
  sig_xy  <- max(0.1, s$total)
  use_3d  <- identical(input$sim_view_mode, "3d")
  z_fac   <- max(0.1, as.numeric(input$sim_crlb_z_factor %||% 0.7))
  sig_z   <- sig_xy * z_fac
  # Share of the variance the floor owns: tells the user whether more photons help.
  floor_share <- if (s$total > 0) s$floor^2 / s$total^2 else 0
  tags$div(
    class = "alert alert-success",
    style = "padding:6px 10px; font-size:0.84em; margin-top:4px;",
    tags$b(sprintf("\u03C3_xy \u2248 %.2f nm", sig_xy)),
    if (use_3d) sprintf("   |   \u03C3_z \u2248 %.2f nm", sig_z),
    tags$br(),
    tags$small(sprintf(
      "Ideal bound %.2f nm, \u00D7k = %.2f nm, floor %.2f nm (%.0f%% of variance)",
      s$crlb, s$photon, s$floor, 100 * floor_share))
  )
})

# ── Phase 1 per-molecule TID count summary ───────────────────────────────────
output$sim_tid_summary <- renderUI({
  res <- sim_result()
  if (is.null(res) || !("mol_id" %in% names(res$fluoro)) ||
      nrow(res$fluoro) == 0) return(NULL)

  # Count detected TIDs per molecule (background TIDs have mol_id = NA)
  real_tids  <- res$fluoro[!is.na(res$fluoro$mol_id), , drop = FALSE]
  if (nrow(real_tids) == 0) return(NULL)
  tids_per_mol <- as.integer(table(real_tids$mol_id))
  # Molecules that lost all fluorophores to detection efficiency = 0 TIDs
  n_zero  <- res$n_prot - length(tids_per_mol)
  counts  <- c(tids_per_mol, rep(0L, max(0L, n_zero)))

  row <- function(lbl, val)
    tags$tr(tags$td(style = "padding:2px 8px;", lbl),
            tags$td(style = "padding:2px 8px;", tags$b(val)))

  tagList(tags$div(
    class = "well well-sm",
    style = "font-size:0.88em; margin-top:6px;",
    tags$b("TIDs per molecule"),
    tags$table(style = "width:100%; border-collapse:collapse;",
      row("Mean",    sprintf("%.2f", mean(counts))),
      row("SD",      sprintf("%.2f", sd(counts))),
      row("Median",  sprintf("%d",   as.integer(median(counts)))),
      row("Range",   sprintf("%d \u2013 %d", min(counts), max(counts))),
      row("Molecules with 0 TIDs",
          sprintf("%d (%.0f%%)", n_zero,
                  100 * n_zero / max(1, res$n_prot))),
      if (res$n_bg > 0)
        row("Background TIDs added", res$n_bg)
    )
  ))
})

# ── Phase 1 precision budget ─────────────────────────────────────────────────
# Bar chart breaking down: linker error | localization precision | combined
# (+ observed RMSE when GT data are available).
output$sim_precision_budget <- renderUI({
  res <- sim_result()
  if (is.null(res)) return(NULL)

  # ── Error budget, all terms as RADIAL RMS displacement (nm) ───────────────
  # Every component must be expressed in the SAME quantity before it can go in
  # quadrature.  The observed value this is compared against is a radial
  # distance |trace median - GT|, so radial RMS is the common currency.
  #
  #   protein   epitope offset: a FIXED magnitude r        -> r^2
  #   linker    |N(0, s)| with s = link*sqrt(pi/2), so the MEAN magnitude is
  #             link_nm but E[|d|^2] = s^2 = (pi/2)*link^2.  Using link^2 here
  #             (as this panel used to) understates the term by 20 %.
  #   dye       per-AXIS sigma -> radial variance is D * dye_spread^2
  #   loc       per-AXIS sigma, and the TID position is the MEDIAN of n_loc
  #             localizations, so the localization term is DIVIDED DOWN by the
  #             trace length.  Omitting that made "observed vs theoretical"
  #             compare two different quantities.
  #             The textbook Var(median) = (pi/2)*sigma^2/n is ASYMPTOTIC and
  #             runs +17 % at n=3, +9 % at n=5 - exactly the short traces
  #             MINFLUX produces.  (pi/2)/(n + pi/2 - 1) holds to ~2 % from
  #             n=3 up.  It smooths over the odd/even-n zigzag in the exact
  #             value, which is right here because Poisson trace lengths mix
  #             both parities anyway.
  D       <- if (isTRUE(res$use_3d)) 3L else 2L
  link    <- res$link_nm
  dspread <- res$dye_spread_nm %||% 0
  prot_r  <- res$protein_r_nm  %||% 0
  sig_xy  <- res$sigma_loc_xy
  sig_z   <- res$sigma_loc_z %||% sig_xy

  # Localizations per REAL trace (background TIDs excluded).  The WHOLE vector
  # is kept, not just its median: trace lengths are Poisson/geometric and the
  # median-variance factor is convex, so evaluating it at the median length
  # under-states the ensemble mean by 5-7 % (Jensen).  Averaging the factor
  # over the realized lengths removes that exactly.
  n_vec     <- numeric(0)
  n_loc_med <- 1
  if (!is.null(res$locs) && nrow(res$locs) > 0 && !is.null(res$fluoro)) {
    tid_i <- suppressWarnings(as.integer(res$locs$tid))
    real  <- !is.na(tid_i) & tid_i >= 1L & tid_i <= nrow(res$fluoro)
    if (any(real)) {
      n_vec     <- as.numeric(table(tid_i[real]))
      n_loc_med <- max(1, stats::median(n_vec))
    }
  }

  # Variance of the median of n Gaussian draws, per unit sigma^2.
  med_var_f <- function(n) (pi / 2) / (pmax(1, n) + pi / 2 - 1)
  var_loc1  <- if (D == 3L) 2 * sig_xy^2 + sig_z^2 else 2 * sig_xy^2
  v_prot    <- prot_r^2
  v_link    <- (pi / 2) * link^2
  v_dye     <- D * dspread^2
  v_loc     <- var_loc1 * if (length(n_vec)) mean(med_var_f(n_vec))
                          else med_var_f(n_loc_med)

  theory <- sqrt(v_prot + v_link + v_dye + v_loc)
  # What NO amount of acquisition time removes: the per-molecule label offset.
  # dye_spread is excluded because it averages once a molecule's TIDs are
  # merged; protein + linker are drawn once per molecule and never do.
  lab_floor <- sqrt(v_prot + v_link)

  # Observed as an RMS too (it was a mean of magnitudes, which is a different
  # statistic and always smaller), so the two are directly comparable.
  obs <- if (!is.null(res$rmse_vals) && length(res$rmse_vals) >= 1)
           sqrt(mean(res$rmse_vals^2)) else NA_real_

  # Per-MOLECULE error after merging that molecule's TIDs — the quantity that
  # decides whether two molecules can be told apart.  dye_spread averages over
  # the dyes detected, sigma_loc over all their localizations.
  n_dye_det <- max(1, res$n_fluor * (res$eff_pct %||% 100) / 100)
  n_tid_mol <- n_dye_det * (1 + (res$blink_lambda %||% 0))
  n_loc_mean <- if (length(n_vec)) mean(n_vec) else n_loc_med
  eps_mol   <- sqrt(v_prot + v_link + D * dspread^2 / n_dye_det +
                    var_loc1 * med_var_f(n_tid_mol * n_loc_mean))

  # Two molecules are separable at ~k x their position error. k = 2 is kept
  # for continuity with previous numbers; strictly, two INDEPENDENT estimates
  # each of error eps separate at 2*sqrt(2)*eps.
  resolve_limit <- 2 * eps_mol
  max_val       <- max(theory, obs, 1, na.rm = TRUE) * 1.05

  make_bar <- function(label, value, color) {
    if (is.na(value) || value < 0) return(NULL)
    w_pct <- sprintf("%.1f%%", 100 * min(value, max_val) / max_val)
    tags$div(style = "margin-bottom:7px;",
      tags$div(
        style = paste0("font-size:0.80em; margin-bottom:2px;",
                       " display:flex; justify-content:space-between;"),
        tags$span(label),
        tags$b(sprintf("%.1f nm", value))
      ),
      tags$div(
        style = "background:#e8e8e8; border-radius:3px; height:11px;",
        tags$div(style = sprintf(
          "background:%s; width:%s; height:11px; border-radius:3px;",
          color, w_pct))
      )
    )
  }

  tags$div(
    class = "well well-sm",
    style = "font-size:0.88em; margin-top:8px;",
    tags$b("Precision Budget"),
    tags$br(),
    tags$small(sprintf(
      "Radial RMS displacement of one TID from its molecule (%dD, %.0f locs/trace):",
      D, n_loc_med)),
    tags$br(), tags$br(),
    if (prot_r > 0) make_bar("Protein size (epitope offset)", prot_r, "#be185d"),
    make_bar("Linker / label error  (\u221A(\u03C0/2)\u00B7link)",
             sqrt(v_link), "#e07b39"),
    if (dspread > 0)
      make_bar("Dye spread in complex", sqrt(v_dye), "#0891b2"),
    make_bar("Localization  (\u03c3 averaged over the trace)",
             sqrt(v_loc), "#3b82c4"),
    make_bar("Combined - theoretical", theory, "#7c3aed"),
    if (!is.na(obs))
      make_bar(sprintf("Observed TID displacement (RMS, n\u202F=\u202F%d)",
                       length(res$rmse_vals)),
               obs, "#16a34a"),
    tags$small(style = "color:#888;",
      "Observed and theoretical should agree; both are radial RMS."),
    # The label floor is the headline: it is what the acquisition CANNOT fix.
    tags$div(
      style = paste0("margin-top:8px; padding:5px 8px; border-radius:4px;",
                     " font-size:0.82em; background:#f1f5f9; color:#334155;"),
      tags$b(sprintf("Label floor: %.1f nm", lab_floor)),
      tags$br(),
      "Protein offset + linker are drawn once per molecule, so no extra ",
      "localizations or acquisition time reduce them."),
    if (!is.na(res$gt_nnd)) {
      resolvable <- res$gt_nnd >= resolve_limit
      never      <- res$gt_nnd < 2 * lab_floor
      tags$div(
        style = sprintf(
          paste0("margin-top:8px; padding:5px 8px; border-radius:4px;",
                 " font-size:0.82em; background:%s; color:%s;"),
          if (resolvable) "#dcfce7" else "#fef9c3",
          if (resolvable) "#166534" else "#854d0e"),
        tags$b(if (resolvable) "\u2714 Structure resolvable"
               else if (never) "\u2716 Not resolvable with this label"
               else "\u26A0 Below resolution limit"),
        tags$br(),
        sprintf(paste0("GT spacing \u2248 %.0f nm vs 2\u00D7 per-molecule error",
                       " (%.1f nm) = %.0f nm.  %s"),
          res$gt_nnd, eps_mol, resolve_limit,
          if (resolvable)
            "Individual molecules can be resolved."
          else if (never)
            paste0("The spacing is below 2\u00D7 the label floor, so this ",
                   "cannot be fixed by imaging longer - a shorter linker ",
                   "is required.")
          else
            paste0("Currently unresolved, but the label floor still allows ",
                   "it - more localizations per molecule would close the gap."))
      )
    }
  )
})

# ── Compare strategies: helper to build a stripped-down scatter plot ─────────
make_cmp_plot <- function(gt, locs, use_3d, title,
                           loc_color = "rgba(220,38,38,0.65)",
                           dark = FALSE) {
  bg_col <- if (dark) "#111111" else "#FFFFFF"
  fg_col <- if (dark) "#DDDDDD" else "#444444"
  tr_type <- if (use_3d) "scatter3d" else "scattergl"
  show_grid <- !isTRUE(input$sim_hide_grid)

  # Compute isotropic axis ranges from GT so both panels share identical scales
  # and 1 nm maps to the same pixel distance on every axis.
  has_z  <- use_3d && "z" %in% names(gt)
  cx     <- mean(range(gt$x, na.rm = TRUE))
  cy     <- mean(range(gt$y, na.rm = TRUE))
  cz     <- if (has_z) mean(range(gt$z, na.rm = TRUE)) else 0
  spans  <- c(diff(range(gt$x, na.rm = TRUE)),
              diff(range(gt$y, na.rm = TRUE)),
              if (has_z) diff(range(gt$z, na.rm = TRUE)) else 0)
  half   <- max(spans, 1) / 2 * 1.15   # 15 % padding, same on all axes
  xrng   <- cx + c(-half, half)
  yrng   <- cy + c(-half, half)
  zrng   <- cz + c(-half, half)

  gt_col <- if (dark) "#FFFFFF" else "#000000"

  # Draw localizations first so GT markers render on top of them.
  plt <- plot_ly() %>%
    add_trace(type = tr_type,
              x = locs$x, y = locs$y,
              z = if (use_3d && "z" %in% names(locs)) locs$z else NULL,
              mode = "markers",
              marker = list(size = if (use_3d) 4 else 7,
                            color = loc_color,
                            opacity = 0.75),
              name = "Localizations",
              hoverinfo = "skip")

  # GT layer drawn on top so it is never hidden by localizations.
  # All shapes use open circle markers (no connecting lines).
  plt <- add_trace(plt,
            type = tr_type,
            x = gt$x, y = gt$y,
            z = if (has_z) gt$z else NULL,
            mode = "markers",
            marker = list(size   = if (use_3d) 10 else 14,
                          color  = gt_col,
                          symbol = "circle-open"),
            name = "GT molecules",
            hoverinfo = "skip")

  plt <- plt %>%
    layout(
      title = list(text = title, font = list(size = 12, color = fg_col)),
      paper_bgcolor = bg_col, plot_bgcolor = bg_col,
      xaxis = list(color = fg_col, title = "x (nm)",
                   scaleanchor = "y", scaleratio = 1,
                   range = xrng, constrain = "domain", showgrid = show_grid, zeroline = show_grid),
      yaxis = list(color = fg_col, title = "y (nm)",
                   range = yrng, constrain = "domain", showgrid = show_grid, zeroline = show_grid),
      showlegend = TRUE,
      legend = list(font = list(size = 10, color = fg_col))
    )
  if (use_3d) plt <- plt %>%
    layout(scene = list(
      aspectmode  = "manual",
      aspectratio = list(x = 1, y = 1, z = 1),
      xaxis = list(title = "x (nm)", color = fg_col, range = xrng,
                   showgrid = show_grid, zeroline = show_grid),
      yaxis = list(title = "y (nm)", color = fg_col, range = yrng,
                   showgrid = show_grid, zeroline = show_grid),
      zaxis = list(title = "z (nm)", color = fg_col, range = zrng,
                   showgrid = show_grid, zeroline = show_grid),
      bgcolor = bg_col))
  plt
}

output$sim_cmp_plot_a <- renderPlotly({
  cmp <- compare_result()
  shiny::validate(shiny::need(!is.null(cmp), "Click Compare."))
  make_cmp_plot(cmp$gt, cmp$a$locs, cmp$use_3d,
                "Strategy A", "rgba(220,38,38,0.75)",
                dark = isTRUE(input$sim_dark_mode))
})

output$sim_cmp_plot_b <- renderPlotly({
  cmp <- compare_result()
  shiny::validate(shiny::need(!is.null(cmp), "Click Compare."))
  make_cmp_plot(cmp$gt, cmp$b$locs, cmp$use_3d,
                "Strategy B", "rgba(37,99,235,0.75)",
                dark = isTRUE(input$sim_dark_mode))
})

output$sim_cmp_table <- renderUI({
  cmp <- compare_result()
  if (is.null(cmp)) return(NULL)

  sig <- cmp$sigma_loc_xy
  a   <- cmp$a;  b <- cmp$b

  pr_a     <- a$protein_r_nm %||% 0
  pr_b     <- b$protein_r_nm %||% 0
  theory_a <- sqrt(pr_a^2 + a$link_nm^2 + sig^2)
  theory_b <- sqrt(pr_b^2 + b$link_nm^2 + sig^2)
  obs_a    <- if (!is.null(a$rmse_vals) && length(a$rmse_vals) >= 1)
                sprintf("%.1f nm (n=%d)", mean(a$rmse_vals), length(a$rmse_vals))
              else "-"
  obs_b    <- if (!is.null(b$rmse_vals) && length(b$rmse_vals) >= 1)
                sprintf("%.1f nm (n=%d)", mean(b$rmse_vals), length(b$rmse_vals))
              else "-"

  cell_style <- "padding:4px 8px; white-space:normal; word-break:break-word;"
  hdr  <- function(...) tags$th(style = paste0(cell_style, " text-align:left;"), ...)
  cel  <- function(...) tags$td(style = cell_style, ...)
  bcel <- function(...) tags$td(style = cell_style, tags$b(...))
  row2 <- function(lbl, va, vb)
    tags$tr(cel(lbl), bcel(va), bcel(vb))

  tags$div(
    class = "well well-sm",
    style = "font-size:0.86em; margin-top:10px; overflow-x:auto;",
    tags$b("Strategy comparison"),
    tags$table(
      style = "width:100%; border-collapse:collapse; margin-top:6px; table-layout:fixed;",
      tags$colgroup(
        tags$col(style = "width:36%;"),
        tags$col(style = "width:32%;"),
        tags$col(style = "width:32%;")
      ),
      tags$thead(tags$tr(
        hdr("Metric"),
        hdr(tags$span(style = "color:#dc2626;", "Strategy A")),
        hdr(tags$span(style = "color:#2563eb;", "Strategy B"))
      )),
      tags$tbody(
        row2("Label",
             a$lbl_note,
             b$lbl_note),
        row2("Fluorophores / molecule (mean)",
             sprintf("%.1f", a$n_fluoro / max(1, nrow(cmp$gt))),
             sprintf("%.1f", b$n_fluoro / max(1, nrow(cmp$gt)))),
        row2("Total TIDs",
             a$n_fluoro, b$n_fluoro),
        row2("Linker distance",
             sprintf("%.0f nm", a$link_nm),
             sprintf("%.0f nm", b$link_nm)),
        row2("Label footprint",
             if (!is.na(a$footprint_nm)) sprintf("~%.0f nm", a$footprint_nm) else "-",
             if (!is.na(b$footprint_nm)) sprintf("~%.0f nm", b$footprint_nm) else "-"),
        if (pr_a > 0 || pr_b > 0)
          row2("Protein radius",
               sprintf("%.1f nm", pr_a),
               sprintf("%.1f nm", pr_b)),
        row2("Combined σ",
             sprintf("%.1f nm", theory_a),
             sprintf("%.1f nm", theory_b)),
        row2("Observed TID displacement",
             obs_a, obs_b)
      )
    )
  )
})

# ── Confocal view: live PSF preview  ────────────────────
output$sim_diff_psf_preview <- renderUI({
  na      <- max(0.5, min(1.5, as.numeric(input$sim_diff_NA %||% 1.4)))
  au      <- max(0.1, as.numeric(input$sim_diff_AU %||% 0.83))
  sig     <- compute_confocal_sigma(640, 670, na, au)
  fwhm    <- round(sig * 2.355)
  au_nm   <- round(1.22 * 670 / na)
  ee      <- confocal_collection_efficiency(au)
  tags$div(
    class = "alert alert-info",
    style = "padding:6px 10px; font-size:0.84em; margin-top:4px;",
    tags$b(sprintf("\u03C3 \u2248 %.1f nm  (FWHM \u2248 %d nm)  \u00B7  signal \u2248 %.0f%%",
                   sig, fwhm, 100 * ee)),
    tags$br(),
    tags$small(sprintf(
      "1 AU = %d nm  \u00B7  pinhole radius = %.0f nm  \u00B7  closing sharpens + sections but dims",
      au_nm, au * au_nm / 2))
  )
})

# ── Confocal view plot ───────────────────────────────────────────────────────
output$sim_diff_plot <- renderPlotly({
  res <- sim_result()
  shiny::validate(
    shiny::need(!is.null(res) && nrow(res$gt) > 0,
                "Click 'Generate' first to build the GT positions."))

  dark   <- isTRUE(input$sim_dark_mode)
  bg_col <- if (dark) "#111111" else "#FFFFFF"
  fg_col <- if (dark) "#DDDDDD" else "#444444"
  gc     <- if (dark) "rgba(255,255,255,0.1)" else "rgba(0,0,0,0.1)"

  na      <- max(0.5, min(1.5, as.numeric(input$sim_diff_NA %||% 1.4)))
  au      <- max(0.1, as.numeric(input$sim_diff_AU %||% 0.83))
  psf_sig <- compute_confocal_sigma(640, 670, na, au)
  n_ph    <- max(1, as.numeric(input$sim_diff_photons %||% 200))
  bg_ph   <- max(0, as.numeric(input$sim_diff_bg      %||% 5))

  # Pinhole physics: the collected signal is scaled by the encircled energy
  # (closing the pinhole → dimmer, noisier via Poisson), while the out-of-focus
  # haze grows with the pinhole AREA (~AU², opening → washed-out background).
  ee       <- confocal_collection_efficiency(au)
  n_ph_eff <- n_ph * ee                 # signal photons actually collected
  bg_eff   <- bg_ph * au^2              # out-of-focus background admitted

  # Render from the UNIQUE fluorophore pool (one physical dye = one PSF), not
  # the over-counted MINFLUX detection list - re-detecting a molecule over the
  # acquisition does not make its confocal spot brighter.
  #   • Each fluorophore = one PSF at its displaced position, n_ph photons.
  #     Multi-dye labels (e.g. polyclonal + secondary) still stack into a
  #     brighter spot because they are distinct dyes.
  #   • Background TIDs (non-specific binding) render as single n_ph spots.
  fp     <- if (!is.null(res$fluoro_pool)) res$fluoro_pool else res$fluoro
  src_x  <- if (nrow(fp) > 0) fp$x else numeric(0)
  src_y  <- if (nrow(fp) > 0) fp$y else numeric(0)
  src_ph <- rep(n_ph_eff, length(src_x))

  if (!is.null(res$bg_pos) && nrow(res$bg_pos) > 0) {
    src_x  <- c(src_x,  res$bg_pos$x)
    src_y  <- c(src_y,  res$bg_pos$y)
    src_ph <- c(src_ph, rep(n_ph_eff, nrow(res$bg_pos)))
  }

  # Cap at 5 000 sources total to keep rendering fast (sample uniformly)
  if (length(src_x) > 5000L) {
    idx    <- sample(length(src_x), 5000L)
    src_x  <- src_x[idx]; src_y <- src_y[idx]; src_ph <- src_ph[idx]
  }

  shiny::validate(shiny::need(length(src_x) > 0,
    "No fluorophores to render. Adjust labeling efficiency or generate first."))

  px_user  <- max(1, as.numeric(input$sim_diff_px %||% 50))
  # Pin the confocal window to the FOV when it's active so the confocal spans
  # the same ROI as the MINFLUX view rather than cropping to the sources.
  use_fov  <- sim_use_fov()
  fov_ext  <- if (use_fov) {
                fx <- max(1, as.numeric(input$sim_fov_x %||% 4000))
                fy <- max(1, as.numeric(input$sim_fov_y %||% 4000))
                c(-fx / 2, fx / 2, -fy / 2, fy / 2)
              } else NULL
  img_data <- render_diffraction_image(src_x, src_y, psf_sig, n_ph_eff, bg_eff,
                                        px_nm_user = px_user,
                                        photons_vec = src_ph,
                                        extent = fov_ext)
  shiny::validate(shiny::need(!is.null(img_data), "Rendering failed."))

  # FWHM = 2.355 × σ: optical resolution - apparent width of a point source
  # and the minimum resolvable separation between two identical emitters.
  fwhm_nm   <- round(psf_sig * 2.355)
  title_txt <- sprintf(
    "Confocal view  (\u03C3 \u2248 %.0f nm, FWHM \u2248 %d nm, signal %.0f%% @ %.2f AU)",
    psf_sig, fwhm_nm, 100 * ee, au)

  # Hot Red LUT (ImageJ "Red Hot") - fixed regardless of dark mode; the
  # confocal image is always rendered on black. Low → high:
  # black → red → orange → yellow → white.
  cscale <- list(
    list(0.00, "rgb(0,0,0)"),
    list(0.25, "rgb(128,0,0)"),
    list(0.45, "rgb(255,0,0)"),
    list(0.65, "rgb(255,140,0)"),
    list(0.85, "rgb(255,255,0)"),
    list(1.00, "rgb(255,255,255)"))

  # The confocal image is inherently a black-background picture, so render it
  # on black with no gridlines and white axis text regardless of dark mode.
  cf_bg <- "#000000"; cf_fg <- "#DDDDDD"

  # Scale bar - same controls as the MINFLUX view, forced white so it stands
  # out on the black image, measured against the rendered image width.  The
  # full image span is stored so the zoom observer can re-scale the bar.
  x_span <- diff(range(img_data$xs, na.rm = TRUE))
  sim_diff_full(x_span)
  sb <- NULL
  if (isTRUE(input$sim_show_scale_bar)) {
    sb <- make_scale_bar(
      suppressWarnings(as.numeric(input$sim_scale_bar_length)),
      x_span,
      suppressWarnings(as.numeric(input$sim_scale_bar_thickness)),
      dark = TRUE,
      scale_font_size = suppressWarnings(
        as.numeric(input$sim_scale_font_size %||% 18)))
  }

  lay_args <- list(
    paper_bgcolor = cf_bg, plot_bgcolor = cf_bg,
    font  = list(color = cf_fg),
    xaxis = list(title = "X (nm)", color = cf_fg,
                 scaleanchor = "y", scaleratio = 1,
                 showgrid = FALSE, zeroline = FALSE),
    yaxis = list(title = "Y (nm)", color = cf_fg,
                 showgrid = FALSE, zeroline = FALSE),
    title = list(text = title_txt, font = list(color = cf_fg)))
  if (!is.null(sb)) {
    lay_args$shapes      <- sb$shapes
    lay_args$annotations <- sb$annotations
  }

  fig <- plot_ly(source = "sim_diff") %>%
    add_trace(
      type       = "heatmap",
      x          = img_data$xs, y = img_data$ys, z = img_data$img,
      colorscale = cscale,
      # Auto-contrast the display (like ImageJ auto B&C) so it looks like a real
      # confocal image. The pinhole still shows through the DATA - closing it
      # collects fewer photons (grainier via Poisson) and tightens the PSF;
      # opening it adds out-of-focus background - not through raw brightness.
      showscale  = FALSE, hoverinfo = "skip",
      name       = "Diffraction-limited"
    )
  sim_diff_ready(TRUE)
  do.call(layout, c(list(fig), lay_args)) %>%
    plotly::event_register("plotly_relayout") %>%
    config(displaylogo = FALSE,
           toImageButtonOptions = list(format = "svg",
                                       filename = "diffraction_limited"))
})

# ── Phase 1 CSV export ───────────────────────────────────────────────────────
# Each row is one localisation.  When GT data exist, mol_id and the GT
# coordinates (gt_x, gt_y, gt_z) are appended so every loc can be traced
# back to its source molecule.  Background TIDs get NA in these columns.
output$sim_dl_csv <- downloadHandler(
  filename = function()
    paste0("sim_localisations_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
  content = function(file) {
    res <- sim_result()
    shiny::validate(shiny::need(!is.null(res) && nrow(res$locs) > 0,
                                "No localizations to export."))
    out <- res$locs[, intersect(c("x", "y", "z", "tid", "tim"),
                                names(res$locs)), drop = FALSE]

    if ("mol_id" %in% names(res$fluoro) && nrow(res$gt) > 0) {
      n_real  <- nrow(res$fluoro)
      tid_idx <- suppressWarnings(as.integer(out$tid))
      valid   <- !is.na(tid_idx) & tid_idx >= 1L & tid_idx <= n_real
      mol_ids <- rep(NA_integer_, nrow(out))
      mol_ids[valid] <- res$fluoro$mol_id[tid_idx[valid]]
      out$mol_id <- mol_ids
      out$gt_x   <- ifelse(is.na(mol_ids), NA_real_, res$gt$x[mol_ids])
      out$gt_y   <- ifelse(is.na(mol_ids), NA_real_, res$gt$y[mol_ids])
      if ("z" %in% names(res$gt))
        out$gt_z <- ifelse(is.na(mol_ids), NA_real_, res$gt$z[mol_ids])
    }

    write.csv(out, file, row.names = FALSE)
  }
)

# ── Acquisition timeline: resolution-over-time plot ──────────────────────────
output$sim_acq_plot <- renderPlotly({
  shiny::validate(shiny::need(
    identical(input$sim_mode %||% "normal", "acquisition"),
    "Switch to Acquisition mode (top of the sidebar) to model an acquisition timeline."))
  a <- sim_acq_live()
  shiny::validate(shiny::need(
    !is.null(a),
    "Click 'Generate' to build the acquisition timeline."))
  dark   <- isTRUE(input$sim_dark_mode)
  bg_col <- if (dark) "#111111" else "#FFFFFF"
  fg_col <- if (dark) "#DDDDDD" else "#444444"
  gc     <- if (dark) "rgba(255,255,255,0.1)" else "rgba(0,0,0,0.1)"

  # Work in minutes; lead with the detection rate (traces/hour) since that is
  # the curve directly comparable to a real run (early peak → constant floor).
  rt_min  <- a$rate_t / 60
  t_min   <- a$t / 60
  # % resolved is relative to the LABELED molecules, so the curve plateaus at
  # the fraction actually reachable in this run (reaches 100 % when all
  # labeled molecules are found).
  cum_pct <- 100 * a$cum_unique / max(1, a$n_mol_resolvable)
  rmax    <- max(a$rate_v, 1, na.rm = TRUE)

  p <- plot_ly()
  if (!is.null(a$rate_t))
    p <- p %>% add_trace(x = rt_min, y = a$rate_v, type = "scatter",
                         mode = "lines", fill = "tozeroy",
                         fillcolor = "rgba(224,123,57,0.22)",
                         line = list(color = "#e07b39", width = 1.5),
                         name = "Detection rate", hoverinfo = "x+y")
  # Optional decomposition of the PREDICTED rate into its structure (decaying +
  # bleaching) and constant background components - makes explicit that the
  # orange total = structure + background, and the floor is background.
  # Convention: solid = model, dotted = measured (real run, below).
  if (isTRUE(input$sim_real_show_split) && !is.null(a$rate_struct)) {
    p <- p %>%
      add_trace(x = rt_min, y = a$rate_struct, type = "scatter",
                mode = "lines",
                line = list(color = "#2a9d8f", width = 1.5),
                name = "Structure (model)", hoverinfo = "x+y") %>%
      add_segments(x = 0, xend = max(rt_min),
                   y = a$bg_rate, yend = a$bg_rate,
                   line = list(color = "#b07aa1", width = 1.5),
                   name = "Background (model)", inherit = FALSE)
  }
  # Late-phase floor line (median of the last quarter of the run).
  if (!is.null(a$rate_floor) && is.finite(a$rate_floor))
    p <- p %>% add_segments(x = 0, xend = max(rt_min),
                            y = a$rate_floor, yend = a$rate_floor,
                            line = list(color = "#e07b39", width = 1,
                                        dash = "dot"),
                            name = sprintf("floor ~%.0f/hr", a$rate_floor),
                            inherit = FALSE)
  # Cumulative molecules resolved (%) on the secondary axis.
  p <- p %>% add_trace(x = t_min, y = cum_pct, type = "scatter",
                       mode = "lines", yaxis = "y2",
                       line = list(color = "#3b82c4", width = 2),
                       name = "Molecules resolved", hoverinfo = "x+y")
  # Measured detection rate from a loaded real run (calibration overlay).
  # Optionally decomposed into the spatially-split structure / background.
  real <- sim_real_rv()
  if (!is.null(real)) {
    if (isTRUE(input$sim_real_show_split) && !is.null(real$rate_v_struct)) {
      p <- p %>%
        add_trace(x = real$rate_t / 60, y = real$rate_v_struct,
                  type = "scatter", mode = "lines",
                  line = list(color = "#2a9d8f", width = 1.5, dash = "dot"),
                  name = "Measured structure", hoverinfo = "x+y") %>%
        add_trace(x = real$rate_t / 60, y = real$rate_v_bg,
                  type = "scatter", mode = "lines",
                  line = list(color = "#b07aa1", width = 1.5, dash = "dot"),
                  name = "Measured background", hoverinfo = "x+y")
    } else {
      p <- p %>% add_trace(x = real$rate_t / 60, y = real$rate_v,
                           type = "scatter", mode = "lines",
                           line = list(color = "#777777", width = 1.5,
                                       dash = "dot"),
                           name = "Measured (real run)", hoverinfo = "x+y")
    }
  }
  # Resolution marker only when the threshold is actually reached.
  if (!is.null(a$res_time) && is.finite(a$res_time))
    p <- p %>% add_segments(x = a$res_time / 60, xend = a$res_time / 60,
                            y = 0, yend = rmax,
                            line = list(color = "#16a34a", width = 1.5,
                                        dash = "dash"),
                            name = sprintf("%.0f%% resolved", 100 * a$resolve_q),
                            inherit = FALSE)
  p %>%
    layout(
      paper_bgcolor = bg_col, plot_bgcolor = bg_col,
      font   = list(color = fg_col),
      margin = list(l = 60, r = 95, t = 50, b = 110),
      xaxis  = list(title = "Acquisition time (min)", color = fg_col,
                    gridcolor = gc),
      yaxis  = list(title = "Detection rate (traces/hour)", color = "#e07b39",
                    gridcolor = gc, rangemode = "tozero", automargin = TRUE),
      yaxis2 = list(title = list(text = "Molecules resolved (%)",
                                 standoff = 22),
                    color = "#3b82c4", overlaying = "y", side = "right",
                    range = c(0, 105), showgrid = FALSE, automargin = TRUE),
      title  = list(text = "Detection rate & resolution over time",
                    font = list(color = fg_col), x = 0.02, xanchor = "left",
                    y = 0.97, yanchor = "top"),
      showlegend = TRUE,
      legend = list(orientation = "h", x = 1, xanchor = "right",
                    y = -0.28, yanchor = "top", font = list(size = 10))
    ) %>%
    config(displaylogo = FALSE)
})

# ── Acquisition timeline summary ─────────────────────────────────────────────
output$sim_acq_summary <- renderUI({
  a <- sim_acq_live()
  if (is.null(a)) return(NULL)

  fmt_t <- function(s) if (s >= 60) sprintf("%.1f min", s / 60)
                       else          sprintf("%.1f s", s)
  reached <- !is.null(a$res_time) && is.finite(a$res_time)
  capped  <- !is.null(a$struct_cap) && a$struct_cap > 0 &&
             a$n_real_tid >= a$struct_cap

  row <- function(lbl, val)
    tags$tr(tags$td(style = "padding:2px 8px;", lbl),
            tags$td(style = "padding:2px 8px;", tags$b(val)))

  tags$div(
    class = "well well-sm", style = "font-size:0.88em; margin-top:10px;",
    tags$b("Acquisition timeline"), tags$br(),
    tags$table(style = "width:100%; border-collapse:collapse;",
      row("Total acquisition time", fmt_t(a$T_sec)),
      row("Structure detection rate",
          sprintf("%.0f traces/hour  (mean %.0f)", a$struct_rate, a$rate_mean)),
      row("Structure trace IDs",
          if (capped) sprintf("%d  (photon budget reached)", a$n_real_tid)
          else        sprintf("%d", a$n_real_tid)),
      row("Molecules resolved",
          sprintf("%d / %d labeled  (%.0f%%; %d of %d total)",
                  a$n_mol_det, a$n_mol_resolvable, 100 * a$resolved_frac,
                  a$n_mol_det, a$n_mol_total)),
      row(sprintf("Time to %.0f%% resolved", 100 * a$resolve_q),
          if (reached) fmt_t(a$res_time)
          else sprintf("not reached (%.0f%% max)", 100 * a$resolved_frac)),
      if (reached)
        row("Time after resolved",
            sprintf("%s  (%.0f%% of run, mostly re-detecting)",
                    fmt_t(max(0, a$T_sec - a$res_time)),
                    100 * max(0, a$T_sec - a$res_time) / max(1e-9, a$T_sec))),
      row("Detections per molecule",
          sprintf("%.1f× on average  (%.0f%% of traces are repeats)",
                  a$det_per_mol, 100 * a$overcount_frac)),
      if (!is.null(a$n_bg) && a$n_bg > 0)
        row("Background TIDs",
            sprintf("%d  (@ %.0f/hr)", a$n_bg, a$bg_rate))
    ),
    tags$small(style = "color:#888;",
      if (reached)
        paste("Resolved early, then later traces re-detect the same",
              "molecules - the over-counting you group as repeat blinks.",
              "Increase the rate or duration to resolve more; decrease them",
              "(or shorten the run) to see partial resolution.")
      else
        paste(sprintf("Only %.0f%% of the structure resolved.",
                      100 * a$resolved_frac),
              "Increase the detection rate or duration (or improve labeling)",
              "to reach the resolve threshold."))
  )
})

# ── Phase 1 GT positions export ──────────────────────────────────────────────
output$sim_dl_gt_csv <- downloadHandler(
  filename = function()
    paste0("sim_gt_positions_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
  content = function(file) {
    res <- sim_result()
    shiny::validate(shiny::need(!is.null(res) && nrow(res$gt) > 0,
                                "No ground-truth positions to export."))
    cols <- intersect(c("x","y","z"), names(res$gt))
    gt_out <- res$gt[, cols, drop = FALSE]
    gt_out$mol_id <- seq_len(nrow(gt_out))
    write.csv(gt_out, file, row.names = FALSE)
  }
)
