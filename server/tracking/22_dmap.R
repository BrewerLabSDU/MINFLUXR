# ==========================================================================
# dmap_bundle() - SPATIALLY RESOLVED DIFFUSION MAP
# ==========================================================================
#
# GOAL
#   Every other estimator here collapses the data to ONE number (or a few
#   states) for the whole field of view. This asks instead: does D vary with
#   LOCATION? It bins the field into cells and estimates a local D in each,
#   revealing slow zones (membranes, walls, condensates) vs fast interior.
#
# THE MATH - a SINGLE-LAG estimator, not a fitted MSD curve
#   Rearranging <r²> = 2dDΔt for one step gives the local diffusivity directly:
#
#        D_local = r² / (2·d·Δt)        (2d = 4 in 2D, 6 in 3D)
#
#   Each step is assigned to the cell containing its MIDPOINT, and the cell's D
#   is the MEAN of its steps' D, with the noise floor subtracted once per cell:
#
#        D_cell = max(0, mean(r²/(4Δt)) − σ²·mean(1/Δt)) / (1 − 2R·t_c·mean(1/Δt))
#
#   (the denominator is the motion-blur correction added 2026-09-17; without it
#   the map read −30% on blurred simulated tracks. Estimator: est_dmap.R)
#
#   WHY THE MEAN AND NOT THE MEDIAN. The relation <r²> = 4DΔt + 4σ² holds for
#   the MEAN only. A single 2D step's r² is exponentially distributed, so its
#   median sits at ln2 ≈ 0.69 of the mean. The old median form, with the noise
#   floor then subtracted per step and clamped at 0, measured (σ 5 nm, Δt 0.2 ms)
#   −50% at D = 0.2 µm²/s and read exactly 0 at D ≤ 0.05. The mean form is
#   unbiased at every D, whatever the step distribution. The price is
#   robustness: one mis-linked jump inflates its cell, so the per-cell standard
#   error is returned alongside D and upstream track cleaning carries the load.
#   Subtracting per cell (not per step) avoids clamping each step at 0, which
#   would rectify the noise and bias slow cells upward.
#
#   WHY SINGLE-LAG RATHER THAN A PER-CELL MSD FIT: a cell contains few steps.
#   Fitting a multi-parameter MSD curve per cell would be badly under-
#   constrained, whereas the single-lag form needs only first moments. It is
#   the short-lag limit of the covariance-based estimator (Vestergaard,
#   Blainey & Flyvbjerg 2014, Phys. Rev. E 89:022726); MINFLUX-specific
#   validity in Vogler et al. 2025 (Commun. Biol. 8:1573).
#   [Previously mis-cited to Ostersehlt 2022 - a DNA-PAINT imaging paper with
#   no diffusion mapping. Corrected.]
#
#   NOISE FLOOR. An immobile molecule still yields apparent D = σ²/Δt (header
#   §2b eq. 3), so without correction every cell reads as diffusing. Subtracting
#   that offset is what lets genuinely immobile regions show as ~0.
#
# THE MINFLUX ADAPTATION
#   "Integrated steps only" drops steps spanning dark cycles. Two reasons: a
#   long-gap step has a large Δt whose midpoint is not a meaningful location
#   (the molecule wandered between the two ends), and it would smear the map.
#
# OBJECTS IT CALLS (inputs)
#   trk(), trk_cycle_time(), trk_sigma_loc_eff()
#   inputs: trk_dmap_run (Run trigger), trk_dmap_bins, trk_dmap_min_pts,
#           trk_dmap_integrated_only, trk_dmap_sub_noise
#
# STEP BY STEP
#   1. Build per-step r², Δt and midpoint (x,y) within each track.
#   2. Optionally keep only single-cycle (integrated) steps.
#   3. Bin midpoints onto a regular grid.
#   4. Per cell: D = [mean(r²/(4Δt)) − σ²·mean(1/Δt)]₊ / (1 − 2R·t_c·mean(1/Δt)), plus its SE.
#   5. Blank cells with fewer than `min_pts` steps - too few to trust.
#
# CAVEAT: single-lag D is an INSTANTANEOUS estimate; it says nothing about
# confinement or anomaly. And where sampling is spatially uneven, the map's
# blank/noisy regions reflect coverage, not biology.
# ==========================================================================
dmap_bundle <- eventReactive(input$trk_dmap_run, {
  req(isTRUE(input$trk_show_dmap))
  df <- trk(); req(nrow(df) > 0)
  t_c_dmap <- tryCatch(trk_cycle_time(), error = function(e) NA_real_)
  # Motion blur needs the cycle time. Without it the fit would silently run with no
  # blur correction, which reads D about a third low at MINFLUX cycle times, so the
  # panel asks for the calibration instead of guessing.
  shiny::validate(shiny::need(is.finite(t_c_dmap) && t_c_dmap > 0,
    "Cycle time unknown - run Data Quality & Calibration first (the blur correction needs it)."))
  # σ_loc floor is optional (toggle); the app-wide value is used when set.
  sig_loc  <- tryCatch(trk_sigma_loc_eff(), error = function(e) NA_real_)
  sig_use  <- if (isTRUE(input$trk_dmap_sub_noise) && is.finite(sig_loc) && sig_loc > 0) sig_loc else NA_real_
  # Estimator lives in est_dmap.R (pure, tested); blur offset as in the MSD fits.
  withProgress(message = "Computing D-map…", value = 0.1, {
    res <- dmap_compute(df, n_bins = input$trk_dmap_bins, min_pts = input$trk_dmap_min_pts,
                        t_c = t_c_dmap, integrated_only = isTRUE(input$trk_dmap_integrated_only),
                        sigma_loc = sig_use,
                        R_abs = trk_blur()$R_abs)
  })
  shiny::validate(shiny::need(is.null(res$error), res$error %||% ""))
  res
})

output$trk_dmap_plot <- renderPlotly({
  bd     <- dmap_bundle()
  dark   <- isTRUE(isolate(input$trk_dark_mode))
  bg_col <- if (dark) "#000000" else "#FFFFFF"
  fg_col <- if (dark) "#FFFFFF" else "#333333"

  # D_mat[ix, iy]: ix = x-column, iy = y-row
  # Plotly heatmap z[row, col] => transpose so y is rows, x is cols
  z     <- t(bd$D_mat)
  z_log <- log10(z); z_log[!is.finite(z_log)] <- NA_real_

  finite_vals <- z_log[is.finite(z_log)]
  shiny::validate(shiny::need(length(finite_vals) > 0,
    "No cells met the minimum step threshold - lower 'Min steps per cell' or reduce grid resolution."))

  # Clip to 2nd–98th percentile so one extreme cell doesn't dominate the scale
  zmin <- floor(quantile(finite_vals, 0.02, na.rm = TRUE))
  zmax <- ceiling(quantile(finite_vals, 0.98, na.rm = TRUE))
  if (zmin >= zmax) { zmin <- zmin - 1; zmax <- zmax + 1 }

  cs <- list(
    list(0,    "#08306b"),
    list(0.25, "#2171b5"),
    list(0.5,  "#f7f7f7"),
    list(0.75, "#f4a582"),
    list(1,    "#b2182b")
  )

  # Build hover text matrix (D in nm²/ms + step count)
  n_row <- nrow(z_log); n_col <- ncol(z_log)
  N_t   <- t(bd$N_mat)
  hover <- matrix("", n_row, n_col)
  for (r in seq_len(n_row))
    for (c in seq_len(n_col)) {
      d  <- bd$D_mat[c, r]
      se <- bd$SE_mat[c, r]
      n  <- bd$N_mat[c, r]
      hover[r, c] <- if (is.finite(d))
        sprintf("x: %.0f nm<br>y: %.0f nm<br>D: %.3g \u00B1 %.2g nm\u00B2/ms (SE)<br>n steps: %d",
                bd$x_centers[c], bd$y_centers[r], d, se, n)
      else sprintf("x: %.0f nm<br>y: %.0f nm<br>(no data)", bd$x_centers[c], bd$y_centers[r])
    }

  fig <- plot_ly(
    x = bd$x_centers, y = bd$y_centers, z = z_log,
    type = "heatmap",
    colorscale = cs, zmin = zmin, zmax = zmax,
    text = hover, hovertemplate = "%{text}<extra></extra>",
    colorbar = list(
      title = list(text = "log\u2081\u2080 D<br>(nm\u00B2/ms)"),
      tickvals = seq(zmin, zmax, by = 1),
      ticktext = paste0("10^", seq(zmin, zmax, by = 1)),
      len = 0.6
    )
  )

  # Optional: overlay step count as text annotations
  if (isTRUE(isolate(input$trk_dmap_show_n))) {
    ann <- list()
    for (r in seq_len(n_row))
      for (c in seq_len(n_col))
        if (bd$N_mat[c, r] > 0)
          ann <- c(ann, list(list(
            x = bd$x_centers[c], y = bd$y_centers[r],
            text = as.character(bd$N_mat[c, r]),
            showarrow = FALSE, xref = "x", yref = "y",
            font = list(size = 8, color = "rgba(0,0,0,0.6)")
          )))
  } else { ann <- list() }

  fig %>% layout(
    xaxis = list(title = "x (nm)", color = fg_col,
                 scaleanchor = "y", scaleratio = 1),
    yaxis = list(title = "y (nm)", color = fg_col),
    paper_bgcolor = bg_col, plot_bgcolor = bg_col,
    font = list(color = fg_col),
    annotations = ann,
    uirevision = "dmap"
  )
})

output$trk_dmap_summary <- renderUI({
  bd <- dmap_bundle()
  D_vals <- bd$D_mat[is.finite(bd$D_mat) & bd$D_mat > 0]
  if (!length(D_vals)) return(NULL)
  corr_lbl <- paste0(
    if (isTRUE(bd$integrated_only)) "integrated steps only" else "all steps",
    "  |  ",
    if (isTRUE(bd$noise_subtracted))
      sprintf("noise floor subtracted (\u03C3_loc %.2f nm)", bd$sigma_loc_nm)
    else "no noise-floor subtraction")
  tags$div(style = "font-size:0.88em; margin:4px 0;",
    sprintf(
      "Filled cells: %d / %d  |  D: %.3g \u2013 %.3g nm\u00B2/ms  |  median D: %.3g nm\u00B2/ms  |  cell size: %.0f \u00D7 %.0f nm  |  median \u0394t: %.2g ms  |  total steps: %d",
      length(D_vals), bd$n_bins^2,
      min(D_vals), max(D_vals), median(D_vals),
      bd$cell_w_nm, bd$cell_h_nm,
      bd$dt_med_ms, bd$n_steps
    ),
    tags$br(),
    tags$span(style = "color:#888;", corr_lbl)
  )
})

output$trk_dl_dmap_csv <- downloadHandler(
  filename = function() paste0("dmap_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
  content  = function(file) {
    bd  <- dmap_bundle()
    rows <- lapply(seq_len(bd$n_bins), function(ix)
      lapply(seq_len(bd$n_bins), function(iy)
        data.frame(x_nm = bd$x_centers[ix], y_nm = bd$y_centers[iy],
                   D_nm2_per_ms = bd$D_mat[ix, iy],
                   D_se_nm2_per_ms = bd$SE_mat[ix, iy],
                   n_steps = bd$N_mat[ix, iy],
                   stringsAsFactors = FALSE)
      )
    )
    out <- do.call(rbind, do.call(c, rows))
    out <- out[is.finite(out$D_nm2_per_ms), ]
    write.csv(out, file, row.names = FALSE)
  }
)
