## ─────────────────────────────────────────────────────────────────────────────
## server/tracking/19_directionality.R  -  panel: directionality
## ─────────────────────────────────────────────────────────────────────────────
## Part of tracking_server(). Not a standalone script: its statements are evaluated in
## the server function's environment, in file-name order (server/parts.R).
##
## Turning angles and the velocity autocorrelation.
##
## INTERFACE (generated: Rscript tools/server_map.R --headers; do not edit by hand)
##   DEFINES for other parts:
##     (nothing)
##   READS from other parts:
##     01_state_sigma: trk_sigma_loc_eff
##     07_cycle_time_sequence: trk_cycle_time
##     11_data_frames: trk
## END INTERFACE
## ─────────────────────────────────────────────────────────────────────────────

# ============================================================================
# dir_bundle() - DIRECTIONALITY (is the motion persistent, random, or reversing?)
# ============================================================================
#
# GOAL
#   Free diffusion has no memory of direction: each step is independent, so the
#   angle between consecutive steps is completely random. Any DEVIATION from
#   that is informative - persistence points to directed transport, systematic
#   reversal points to confinement or tethering.
#
# THE MATH
#   (a) TURN ANGLE between consecutive displacement vectors:
#
#          cos θ = (Δr_i · Δr_{i+1}) / (|Δr_i| · |Δr_{i+1}|),   θ = acos(...)
#
#       Reading the distribution of θ:
#          flat/uniform  -> free Brownian (no preferred direction)
#          peak near 0°  -> persistent / directed motion
#          peak near 180 -> anti-persistent: confinement, tethering, OR noise
#
#       Crucially this null is Δt-INDEPENDENT: the angle between two
#       independent isotropic vectors is uniform no matter how long each step
#       lasted. So the turn-angle test needs no Δt normalisation at all.
#
#   (b) VELOCITY AUTOCORRELATION over unit step vectors v̂ = Δr/|Δr|:
#
#          C(τ) = < v̂(t) · v̂(t + τ) >
#
#       C > 0 and slowly decaying = persistent motion; C ~ 0 = free diffusion.
#
#   (c) THE NOISE FLOOR (why a negative first lag is EXPECTED, not a finding).
#       A pure-noise step has RMS length sqrt(2d)·σ ( = 2σ in 2D), because each
#       of the two localizations contributes σ per axis. More importantly the
#       SHARED localization error at the midpoint enters the two neighbouring
#       steps with OPPOSITE sign - the very same effect the CVE exploits. That
#       makes lag-1 anti-correlated even for a perfectly free molecule:
#
#          C(first lag) ≈ − σ² / MSD(one cycle)
#
#       So a small negative first lag is the signature of measurement noise.
#       Only a value substantially more negative than the following lags is
#       evidence of real anti-persistence. The fraction of steps below the
#       noise floor is reported so you can judge how much is noise.
#
# THE MINFLUX ADAPTATION
#   The autocorrelation lag must be ELAPSED TIME, not step index - on
#   irregularly sampled data a step-index lag would draw unevenly spaced points
#   as if evenly spaced and smear the decay. So pairs are binned by the
#   difference of their CUMULATIVE CYCLE INDEX, and the axis is cycles × t_c.
#   Steps spanning more than `trk_dir_max_dark` dark cycles are dropped
#   entirely: their direction is unreliable and a long gap may even bridge two
#   different molecules.
#
# OBJECTS IT CALLS (inputs)
#   trk(), trk_cycle_time(), trk_sigma_loc_eff()
#   inputs: trk_dir_run (Run trigger), trk_dir_min_pts, trk_dir_max_dark
#
# STEP BY STEP
#   1. Split into tracks; per track build steps, their Δt and cycle counts.
#   2. Mark "usable" steps: non-zero length AND <= max_dark dark cycles.
#   3. Turn angles: only between two CONSECUTIVE usable steps (they share a
#      localization, so the angle is meaningful).
#   4. Autocorrelation: for each pair of usable steps in a track, bin the dot
#      product of their unit vectors by the elapsed cycle lag between them.
#   5. Count steps below the noise floor to quantify how much of the signal is
#      localization error rather than motion.
#
# RETURNS a list:
#   angles ............. all turn angles, degrees (pooled across tracks)
#   ac[lag], lags ...... autocorrelation per integer cycle lag
#   lag_time, have_cyc . the lag axis in ms (falls back to step index if the
#                        cycle time is unknown - and says so on the plot)
#   t_cycle, max_dark, r_noise, sigma_loc
#   frac_below_noise ... fraction of steps shorter than the noise floor
#   frac_dropped_dark .. fraction excluded for spanning too many dark cycles
#
# CAVEAT: single short tracks carry little power - these are pooled statistics.
# ============================================================================
dir_bundle <- eventReactive(input$trk_dir_run, {
  req(isTRUE(input$trk_show_dir))
  df <- trk(); req(nrow(df) > 0)
  min_pts  <- max(4L, as.integer(input$trk_dir_min_pts))
  max_dark <- max(0L, as.integer(input$trk_dir_max_dark %||% 2L))
  # Cycle time for the lag axis and the app-wide sigma_loc for the noise floor;
  # the computation is dir_compute() in R/est_directionality.R.
  t_c     <- tryCatch(trk_cycle_time(), error = function(e) NA_real_)
  sig_dir <- tryCatch(trk_sigma_loc_eff(), error = function(e) NA_real_)
  withProgress(message = "Computing directionality…", value = 0, {
    res <- dir_compute(df, t_c, sig_dir, min_pts, max_dark, max_lag = 20L, progress = incProgress)
    shiny::validate(shiny::need(!is.null(res),
      "Too few usable steps - lower 'Min points per TID', raise 'Max dark cycles', or relax filters."))
    res
  })  # withProgress
})

output$trk_dir_angle_plot <- renderPlotly({
  bd   <- dir_bundle()
  dark <- isTRUE(isolate(input$trk_dark_mode))
  bg_col <- if (dark) "#000000" else "#FFFFFF"
  fg_col <- if (dark) "#FFFFFF" else "#333333"

  h <- hist(bd$angles, breaks = seq(0, 180, by = 10), plot = FALSE)
  plot_ly() %>%
    add_bars(
      x = h$mids, y = h$density,
      marker = list(color = "rgba(80,140,210,0.65)",
                    line  = list(color = "rgba(50,100,170,0.8)", width = 0.5))
    ) %>%
    layout(
      xaxis = list(title = "Turn angle (\u00B0)", range = c(0, 180),
                   tickvals = seq(0, 180, 30), color = fg_col),
      yaxis = list(title = "Density", color = fg_col),
      paper_bgcolor = bg_col, plot_bgcolor = bg_col,
      font = list(color = fg_col), showlegend = FALSE,
      uirevision = "dir_angle"
    )
})

output$trk_dir_autocorr_plot <- renderPlotly({
  bd   <- dir_bundle()
  dark <- isTRUE(isolate(input$trk_dark_mode))
  bg_col <- if (dark) "#000000" else "#FFFFFF"
  fg_col <- if (dark) "#FFFFFF" else "#333333"

  # x-axis is real elapsed time (cycle lag × cycle time) when the cycle time is
  # known; step-index only as a labelled fallback.
  xv    <- if (isTRUE(bd$have_cyc)) bd$lag_time else bd$lags
  xlab  <- if (isTRUE(bd$have_cyc)) "Time lag (ms)" else "Lag (steps - cycle time unknown)"
  plot_ly() %>%
    add_lines(x = xv, y = bd$ac,
              line = list(color = "#D06030", width = 2)) %>%
    add_markers(x = xv, y = bd$ac,
                marker = list(color = "#D06030", size = 5)) %>%
    add_lines(x = range(xv, na.rm = TRUE), y = c(0, 0),
              line = list(color = fg_col, width = 1, dash = "dot"),
              showlegend = FALSE) %>%
    layout(
      xaxis  = list(title = xlab, color = fg_col),
      yaxis  = list(title = "Velocity autocorrelation", range = c(-1.05, 1.05),
                    color = fg_col),
      paper_bgcolor = bg_col, plot_bgcolor = bg_col,
      font = list(color = fg_col), showlegend = FALSE,
      uirevision = "dir_ac"
    )
})

output$trk_dir_summary <- renderUI({
  bd <- dir_bundle()
  med_ang <- median(bd$angles, na.rm = TRUE)
  pct_fwd  <- mean(bd$angles <  60, na.rm = TRUE) * 100
  pct_back <- mean(bd$angles > 120, na.rm = TRUE) * 100
  ac1 <- if (length(bd$ac) >= 1 && is.finite(bd$ac[1])) round(bd$ac[1], 3) else NA
  ac1_lbl <- if (isTRUE(bd$have_cyc)) sprintf("AC(%.0f ms)", bd$t_cycle) else "AC(lag 1 step)"
  dark_note <- if (is.finite(bd$frac_dropped_dark %||% NA_real_) && (bd$frac_dropped_dark %||% 0) > 0)
    tags$span(style = "color:#888;", tags$br(),
      sprintf("%.0f%% of steps span more than %d dark cycle(s) and were left out.",
              100 * bd$frac_dropped_dark, bd$max_dark))
  else NULL
  noise_note <- if (is.finite(bd$frac_below_noise %||% NA_real_))
    tags$span(style = "color:#b26a00;", tags$br(),
      sprintf("⚠ %.0f%% of steps are below the localization-noise floor (2·", 100 * bd$frac_below_noise),
      sigloc(),
      sprintf(" = %.1f nm) - their directions are noise-dominated and inflate the backward (180°) peak / lag-1 anti-correlation.",
              bd$r_noise))
  else
    tags$span(style = "color:#888;", tags$br(),
      "Set ", sigloc(), " (run the MSD, or enter it in Acquisition Parameters) to see how much of the turn-angle signal is localization noise.")
  tags$div(style = "font-size:0.88em; margin:4px 0;",
    tags$b(sprintf("n = %d turn angles", length(bd$angles))),
    tags$span(sprintf(
      "  |  median angle = %.1f\u00B0  |  forward (<60\u00B0): %.0f%%  |  backward (>120\u00B0): %.0f%%  |  %s = %s",
      med_ang, pct_fwd, pct_back, ac1_lbl, if (is.na(ac1)) "-" else ac1
    )),
    dark_note,
    noise_note
  )
})

output$trk_dl_dir_csv <- downloadHandler(
  filename = function() paste0("directionality_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
  content  = function(file) {
    bd <- dir_bundle()
    # Two sheets packed into one CSV with a blank-row separator
    lag_ms <- if (isTRUE(bd$have_cyc)) bd$lag_time else rep(NA_real_, length(bd$lags))
    ac_df  <- data.frame(lag_cycles = bd$lags, lag_ms = lag_ms, velocity_autocorr = bd$ac)
    ang_df <- data.frame(turn_angle_deg = bd$angles)
    out <- rbind(
      cbind(sheet = "autocorrelation", ac_df,
            data.frame(turn_angle_deg = rep(NA, nrow(ac_df)))),
      data.frame(sheet = NA, lag_cycles = NA, lag_ms = NA, velocity_autocorr = NA,
                 turn_angle_deg = NA),
      cbind(sheet = "angles", lag_cycles = NA, lag_ms = NA, velocity_autocorr = NA,
            turn_angle_deg = ang_df$turn_angle_deg)
    )
    write.csv(out, file, row.names = FALSE)
  }
)
