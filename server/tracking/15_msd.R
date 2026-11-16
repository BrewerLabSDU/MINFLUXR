# --- Diffusion limit (D_max) panel -----------------------------------------
output$trk_dmax_panel <- renderUI({
  req(isTRUE(input$trk_show_dmax))
  df  <- trk(); req(nrow(df) > 0)
  ub  <- tryCatch(msd_unified_bundle(), error = function(e) NULL)

  L    <- suppressWarnings(as.numeric(input$trk_L_tcp))
  t_hw <- suppressWarnings(as.numeric(input$trk_dmax_t_hw))
  if (!is.finite(L)    || L    <= 0) L    <- DMAX_L_TCP_DEFAULT
  if (!is.finite(t_hw) || t_hw <  0) t_hw <- 0.05

  de <- dmax_estimate(df, L, t_hw)            # shared helper (single source of truth)
  efo_med      <- de$efo_med;      n_pl         <- de$n_pl
  t_loc_eco    <- de$t_loc_eco;    t_loc_dt     <- de$t_loc_dt
  t_loc        <- de$t_loc
  D_max_simple <- de$D_max_simple; D_max_vogler <- de$D_max_vogler

  # --- Measured D from the Unified MSD (eDIFF, ensemble fit) -----------------
  D_meas <- if (!is.null(ub)) {
    fr <- tryCatch(msd_unified_fit(ub), error = function(e) NULL)
    if (!is.null(fr) && is.finite(fr$D)) fr$D / 1e3 else NA_real_  # nm2/ms -> um2/s
  } else NA_real_

  # --- Build output ---------------------------------------------------------
  fmt <- function(x, digits = 3)
    if (is.finite(x)) sprintf("%.*g", digits, x) else "-"

  ratio_simple <- if (is.finite(D_max_simple) && D_max_simple > 0 && is.finite(D_meas))
    D_meas / D_max_simple else NA_real_
  ratio_vogler <- if (is.finite(D_max_vogler) && D_max_vogler > 0 && is.finite(D_meas))
    D_meas / D_max_vogler else NA_real_

  flag_col <- function(ratio) {
    if (!is.finite(ratio)) return("#444")
    if (ratio > 0.8)  return("#CC2200")   # red: likely affected
    if (ratio > 0.5)  return("#BB7700")   # amber: caution
    return("#228833")                      # green: safe
  }
  # Short verdicts - the mechanism (the loop loses lock above D_max, so D reads low)
  # is in README "Fast-Diffusion Ceiling" and the panel footnote below.
  flag_msg <- function(ratio) {
    if (!is.finite(ratio)) return("")
    if (ratio > 0.8)  return(" ⚠ at the ceiling - D is likely under-estimated")
    if (ratio > 0.5)  return(" ⚠ approaching the ceiling")
    return(" ✓ safe")
  }

  use_ratio <- if (is.finite(ratio_vogler)) ratio_vogler else ratio_simple
  main_dmax <- if (is.finite(D_max_vogler)) D_max_vogler else D_max_simple

  htmltools::div(
    style = "font-size:0.88em; color:#333; padding:6px 0;",

    htmltools::tags$table(
      style = "width:100%; border-collapse:collapse;",
      htmltools::tags$tr(
        htmltools::tags$td(style="padding:3px 6px; color:#666;", "L_TCP (final iteration)"),
        htmltools::tags$td(style="padding:3px 6px; font-weight:bold;", sprintf("%g nm", L),
          # D_max goes as L^2, so the wrong iteration moves the ceiling by an order
          # of magnitude silently. Only flag clearly imaging-range diameters: fine
          # patterns are legitimate for slow targets, so the threshold sits well
          # below the 100-300 nm a tracking iteration normally uses.
          if (is.finite(L) && L < 60) htmltools::span(style = "font-weight:normal; color:#B8860B;",
            htmltools::HTML(" &#9888; tight for tracking. Confirm this is the reiterated iteration."))
          else NULL)
      ),
      htmltools::tags$tr(
        htmltools::tags$td(style="padding:3px 6px; color:#666;", "t_loc from ECO/EFO"),
        htmltools::tags$td(style="padding:3px 6px;",
          if(is.finite(t_loc_eco)) sprintf("%.3f ms", t_loc_eco) else "- (needs ECO & EFO)")
      ),
      htmltools::tags$tr(
        htmltools::tags$td(style="padding:3px 6px; color:#666;", "Median Δt within tracks"),
        htmltools::tags$td(style="padding:3px 6px;",
          if(is.finite(t_loc_dt)) sprintf("%.3f ms", t_loc_dt) else "-")
      ),
      htmltools::tags$tr(
        htmltools::tags$td(style="padding:3px 6px; color:#666;", "Median EFO (ψ_MFX)"),
        htmltools::tags$td(style="padding:3px 6px;",
          if(is.finite(efo_med)) sprintf("%.0f Hz", efo_med) else "- (no EFO column)")
      ),
      htmltools::tags$tr(
        htmltools::tags$td(style="padding:3px 6px; color:#666;", "Median ECO (N_PL, photons/loc)"),
        htmltools::tags$td(style="padding:3px 6px;",
          if(is.finite(n_pl)) sprintf("%.1f", n_pl) else "- (no ECO column)")
      ),
      htmltools::tags$tr(style="border-top:1px solid #ddd;",
        htmltools::tags$td(style="padding:4px 6px; color:#666;", htmltools::HTML("D<sub>max</sub> simple (L_TCP² / 16·t_loc)")),
        htmltools::tags$td(style="padding:4px 6px;", if(is.finite(D_max_simple)) sprintf("%.3g µm²/s", D_max_simple) else "-")
      ),
      htmltools::tags$tr(
        htmltools::tags$td(style="padding:4px 6px; color:#666;", htmltools::HTML("D<sub>max</sub> Vogler (EFO/N_PL)")),
        htmltools::tags$td(style="padding:4px 6px;", if(is.finite(D_max_vogler)) sprintf("%.3g µm²/s", D_max_vogler) else "- (needs ECO & EFO)")
      ),
      htmltools::tags$tr(style="border-top:1px solid #ddd;",
        htmltools::tags$td(style="padding:4px 6px; color:#666;", "Measured D (MSD fit)"),
        htmltools::tags$td(style="padding:4px 6px;", if(is.finite(D_meas)) sprintf("%.3g µm²/s", D_meas) else "- (run MSD first)")
      ),
      if (is.finite(use_ratio)) htmltools::tags$tr(
        htmltools::tags$td(style="padding:4px 6px; color:#666;", htmltools::HTML("D / D<sub>max</sub>")),
        htmltools::tags$td(
          style = sprintf("padding:4px 6px; font-weight:bold; color:%s;", flag_col(use_ratio)),
          htmltools::HTML(sprintf("%.4f%s", use_ratio, flag_msg(use_ratio)))
        )
      )
    ),

    htmltools::tags$div(
      style = "margin-top:8px; font-size:0.82em; color:#666;",
      htmltools::HTML("L<sub>TCP</sub>²/(16·t_loc) &nbsp;·&nbsp; Vogler: (L<sub>TCP</sub>²/16)·ψ/(N_PL + ψ·t_hw), "),
      "with L the final iteration reiterated during tracking. The ceiling scales as L². ",
      htmltools::tags$b("Theoretical ceiling"), " (Vogler 2025 report ≈10–11 µm²/s); tracking is only ",
      "demonstrated to ≈2.5 µm²/s, so treat D above ~2–3 µm²/s as suspect even on a green flag."
    )
  )
})

# Per-track D histogram (log10 or linear scale)
output$trk_msd_hist_plot <- renderPlotly({
  req(isTRUE(input$trk_msd_show_hist))
  fits <- msd_pertk_fits()
  shiny::validate(shiny::need(!is.null(fits) && nrow(fits) >= 3L,
    "Need at least 3 successfully fitted tracks. Try lowering Min points per TID."))

  dark    <- isTRUE(isolate(input$trk_dark_mode))
  bg_col  <- if (dark) "#000000" else "#FFFFFF"
  txt_col <- if (dark) "#CCCCCC" else "#333333"
  use_log <- isTRUE(input$trk_msd_hist_log %||% TRUE)

  # Exclude D ≤ 0 (noise-clamped, unphysical); keep only positive fits
  D_vals <- fits$D_um2_s[is.finite(fits$D_um2_s) & fits$D_um2_s > 0]
  shiny::validate(shiny::need(length(D_vals) >= 3L, "No tracks with positive D after fitting."))

  if (use_log) {
    x_data  <- log10(D_vals)
    x_title <- "log10(D) [µm²/s]"
    ann_txt <- sprintf("n = %d tracks  |  median log10(D) = %.2f  |  median D = %.3g µm²/s",
                       length(D_vals), median(x_data), median(D_vals))
  } else {
    x_data  <- D_vals
    x_title <- "D [µm²/s]"
    ann_txt <- sprintf("n = %d tracks  |  median D = %.3g µm²/s",
                       length(D_vals), median(D_vals))
  }

  plot_ly() %>%
    add_trace(x = x_data, type = "histogram", nbinsx = 35,
              marker = list(color = "rgba(48,96,192,0.75)",
                            line = list(color = "rgba(30,60,140,0.9)", width = 0.5))) %>%
    layout(
      showlegend = FALSE,
      plot_bgcolor  = bg_col, paper_bgcolor = bg_col,
      font = list(color = txt_col),
      xaxis = list(title = x_title, color = txt_col,
                   gridcolor = if (dark) "rgba(255,255,255,0.15)" else "rgba(0,0,0,0.12)"),
      yaxis = list(title = "Count", color = txt_col,
                   gridcolor = if (dark) "rgba(255,255,255,0.15)" else "rgba(0,0,0,0.12)"),
      bargap = 0.05,
      annotations = list(list(
        text = ann_txt,
        xref = "paper", yref = "paper", x = 0.5, y = 1.06,
        xanchor = "center", showarrow = FALSE,
        font = list(size = 11, color = txt_col)
      ))
    ) %>%
    config(displaylogo = FALSE,
           toImageButtonOptions = list(format = "svg", filename = "msd_pertk_hist"))
})
output$trk_dl_msd_pertk_csv <- downloadHandler(
  filename = function() paste0("msd_pertk_fits_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
  content  = function(file) write.csv(msd_pertk_fits(), file, row.names = FALSE)
)

# MSD cluster scatter: D (log10 or linear) vs CI, colored by cluster
output$trk_msd_cluster_plot <- renderPlotly({
  req(isTRUE(input$trk_show_msd), isTRUE(input$trk_msd_cluster))
  cl <- tryCatch(msd_cluster_rv(), error = function(e) NULL)
  shiny::validate(shiny::need(!is.null(cl),
    "Clustering not ready. Need at least 3 fitted tracks and k \u2264 n tracks."))
  fits    <- cl$fits
  pal     <- cl$palette
  use_log <- isTRUE(input$trk_msd_hist_log %||% TRUE)

  # Exclude D \u2264 0 as in the histogram
  ok_d <- is.finite(fits$D_um2_s) & fits$D_um2_s > 0
  fits <- fits[ok_d, , drop = FALSE]
  shiny::validate(shiny::need(nrow(fits) >= 3L, "No tracks with positive D."))

  x_data  <- if (use_log) log10(fits$D_um2_s) else fits$D_um2_s
  x_title <- if (use_log) "log10(D) [\u00B5m\u00B2/s]" else "D [\u00B5m\u00B2/s]"
  x_min   <- min(x_data, na.rm = TRUE); x_max <- max(x_data, na.rm = TRUE)

  plt <- plot_ly()
  for (c_lbl in seq_len(cl$n_clusters)) {
    idx <- fits$cluster == c_lbl
    if (!any(idx)) next
    hover_x <- if (use_log)
      paste0("log10(D): ", round(x_data[idx], 3))
    else
      paste0("D: ", round(fits$D_um2_s[idx], 6), " \u00B5m\u00B2/s")
    plt <- plt %>% add_trace(
      x    = x_data[idx],
      y    = fits$conf_index[idx],
      type = "scattergl", mode = "markers",
      marker = list(color = pal[c_lbl], size = 6, opacity = 0.75,
                    line = list(color = "#ffffff", width = 0.4)),
      name = sprintf("Cluster %d (n=%d)", c_lbl, sum(idx)),
      text = paste0("TID: ", fits$tid[idx], "<br>", hover_x,
                    "<br>D: ", round(fits$D_um2_s[idx], 6), " \u00B5m\u00B2/s",
                    "<br>CI: ", round(fits$conf_index[idx], 3)),
      hoverinfo = "text"
    )
  }
  plt %>%
    layout(
      xaxis = list(title = x_title),
      yaxis = list(title = "Confinement index (CI)", type = "log"),
      showlegend = TRUE,
      shapes = list(
        list(type = "line", x0 = x_min, x1 = x_max,
             y0 = 1, y1 = 1, xref = "x", yref = "y",
             line = list(color = "rgba(80,80,80,0.35)", dash = "dash", width = 1))
      ),
      annotations = list(
        list(x = 0.5, y = 1.02, xref = "paper", yref = "paper", showarrow = FALSE,
             text = "CI = 1: free diffusion boundary",
             font = list(size = 10, color = "rgba(100,100,100,0.6)"),
             xanchor = "center", yanchor = "bottom")
      )
    ) %>%
    config(displaylogo = FALSE,
           toImageButtonOptions = list(format = "svg", filename = "msd_cluster_scatter"))
})

# ---- ms hint for the standard MSD cycle-based min-lag input -----------------
# Input is in cycles; show the ms equivalent at the current cycle time.
output$trk_msd_min_lag_cycle_hint <- renderUI({
  t_c <- trk_cycle_time()
  v   <- suppressWarnings(as.numeric(input$trk_msd_min_lag_cyc))
  if (!is.finite(t_c) || t_c <= 0 || !is.finite(v) || v <= 0) return(NULL)
  tags$p(style = "font-size:0.81em; color:#2070a0; margin:2px 0;",
         sprintf("\u2248 %.4g ms  (at %.4g ms/cycle)", v * t_c, t_c))
})

# Sync button: make the standard MSD reproduce the Advanced (cycle) MSD as
# closely as possible.  Auto-detects the cycle time from the modal \u0394t (same
# logic as the Advanced "Auto" button), then aligns every parameter that
# affects the result: \u03c4-bin, fit-from-lag, \u03c4-max, weighting mode, and OLSF.
observeEvent(input$trk_msd_sync_cycle, {
  # 1. Cycle time from the single source of truth (hardware \u2192 manual \u2192 modal \u0394t)
  t_c <- trk_cycle_time()
  if (!is.finite(t_c) || t_c <= 0) {
    showNotification("Could not determine cycle time (no \u0394t data). Set it manually or via the hardware model.",
                     type = "warning", duration = 6)
    return()
  }

  max_lag  <- suppressWarnings(as.integer(input$trk_msd_cycle_max_lag))
  if (!is.finite(max_lag) || max_lag < 5L) max_lag <- 50L
  min_lag  <- suppressWarnings(as.integer(input$trk_msd_cycle_min_lag))
  if (!is.finite(min_lag) || min_lag < 1L) min_lag <- 4L

  # 2. Push the same cycle time to the Advanced panel so both use one value.
  updateNumericInput(session, "trk_msd_cycle_ms",  value = round(t_c, 4))
  # 3. Align standard-MSD binning & fit range (fit window in cycles).
  updateNumericInput(session, "trk_msd_tau_bin",     value = round(t_c, 4))
  updateNumericInput(session, "trk_msd_tau_max",     value = round(max_lag * t_c, 4))
  updateNumericInput(session, "trk_msd_min_lag_cyc", value = min_lag)

  showNotification(
    sprintf(paste0("Synced to cycle = %.4g ms:  \u03c4-bin = 1 cycle  |  \u03c4-max = %d cycles  |  ",
                   "fit from lag %d (%.4g ms).  Both use OLSF + pair-weighting. ",
                   "Click 'Plot MSD' to recompute."),
            t_c, max_lag, min_lag, min_lag * t_c),
    type = "message", duration = 9
  )
}, ignoreInit = TRUE)

# ---- Advanced MINFLUX MSD: cycle-based (Vogler et al. 2025) -----------------
# Resolved cycle-time status shown inside the Advanced MSD panel
output$trk_msd_cycle_status <- renderUI({
  req(isTRUE(input$trk_msd_advanced))
  tc  <- trk_cycle_time()
  if (!is.finite(tc) || tc <= 0)
    return(tags$p(style = "font-size:0.82em; color:#888; margin:2px 0;",
                  "Cycle time: not yet determined."))
  hw  <- trk_hw_cycle_time()
  man <- suppressWarnings(as.numeric(input$trk_msd_cycle_ms))
  src <- if (is.finite(hw) && hw > 0) "hardware model"
         else if (is.finite(man) && man > 0) "manual override"
         else "modal \u0394t (data)"
  tags$p(style = "font-size:0.83em; color:#2070a0; margin:2px 0;",
         sprintf("Cycle time in use: %.4g ms (source: %s)", tc, src))
})


# ============================================================================
# msd_unified_bundle() - THE MSD ENGINE
# ============================================================================
#
# GOAL
#   Build the mean-squared-displacement curve, MSD(τ), from which the
#   diffusion coefficient D is fitted. This is the single MSD estimator in the
#   app; every per-track D, the D histogram, the clustering and the CSV export
#   all derive from this one bundle.
#
#   MSD(τ) = < |r(t+τ) − r(t)|² >  averaged over all pairs separated by τ.
#   For free Brownian motion MSD = 2·d·D·τ + 2·d·σ_loc², so the SLOPE gives D
#   and the INTERCEPT gives the localization error.
#
# THE MINFLUX PROBLEM THIS SOLVES
#   MSD needs pairs grouped by a common time lag τ. MINFLUX has no fixed clock:
#   consecutive localizations are separated by integer multiples of the base
#   cycle, with "dark" cycles skipped at random. Binning by STEP INDEX would
#   put physically different Δt into the same bin and corrupt the curve.
#   Fix: bin by CUMULATIVE CYCLE INDEX (Vogler et al. 2025) - see step 3.
#
# OBJECTS IT CALLS (inputs)
#   trk() ................. the analysis-ready dataset (see its header)
#   trk_cycle_time() ...... base cycle time t_c, ms
#   trk_cluster_cycle()$sw  software/settling overhead subtracted before
#                           converting an interval into a cycle count
#   .track_positive_D() ... optional per-track D>0 quality filter
#   inputs: trk_msd_compute (the Run trigger), trk_uni_max_lag,
#           trk_msd_min_pts, trk_msd_ndim, trk_uni_drop_negd,
#           trk_uni_track_weight
#
# STEP BY STEP
#   1. Resolve the cycle time t_c and the max cycle-lag to compute.
#   2. Take trk(), pick 2D/3D, drop non-finite rows, split into one
#      data.frame per track (`sp`).
#   3. For each track, convert its irregular timestamps into a CUMULATIVE
#      CYCLE INDEX:
#         dcyc   = round((Δt − sw) / t_c)   cycles spanned by each step (>= 1)
#         cycidx = cumsum(dcyc)             position on a linearized timeline
#      The lag between two localizations is then the INTEGER difference of
#      their cycidx - exact regardless of how many cycles were skipped.
#      (Per-interval rounding then cumsum is deliberate: it is robust to a
#      small error in t_c, whereas round(total_Δt / t_c) accumulates that
#      error and produces a sawtooth in the curve.)
#   4. For every pair (i, j) within a track, accumulate into three
#      lag × track matrices:
#         tk_dr2  Σ squared displacement    tk_np  pair count
#         tk_dt   Σ real elapsed Δt  (kept so the x-axis can be TRUE mean Δt,
#                                     not an assumed multiple of t_c)
#   5. Per-track TAMSD = tk_dr2 / tk_np. Optionally drop tracks whose TAMSD
#      slope is <= 0 (non-diffusing / artefact).
#   6. Collapse to the ensemble curve, two weightings available:
#         pair-weighted  Σdr² / Σn   (default; smoother, long tracks dominate)
#         track-weighted mean of per-track TAMSDs (each track counts once -
#                        this is Vogler's eMSD convention)
#   7. Compute the IQR band across tracks and drop lags supported by fewer
#      than 10% of the tracks (`keep_lag`) - the long-lag tail is built from
#      very few, highly correlated pairs and is not trustworthy.
#
# RETURNS a list (this is the contract every consumer relies on):
#   lag ........... integer cycle lags 1..max_lag
#   tau ........... TRUE mean elapsed Δt per lag bin (ms)  <- honest x-axis
#   tau_lin ....... lag × t_c + SW (linearized x-axis = Vogler's exact T_LIN)
#   msd ........... ensemble MSD per lag (nm²)
#   npairs ........ pairs contributing to each lag
#   n_trk_lag ..... tracks contributing to each lag
#   keep_lag ...... logical, TRUE where enough tracks support the lag
#   track_msd ..... lag × track matrix of per-track TAMSDs (nm²)
#   tk_names ...... TIDs, in the column order of track_msd
#   msd_q25/q75 ... inter-quartile band across tracks
#   t_cycle, has_z, ndim, track_weighted, n_dropped_negD, n_tracks
#
# NOTE: units are nm² and ms throughout. Consumers convert to µm²/s at the
# point of display (divide nm²/ms by 1e3).
# ============================================================================
# ── Effective MSD divisor ───────────────────────────────────────────────────
# MSD = 2 · d_eff · D · τ. d_eff is NOT simply the number of recorded coordinates:
# it is set by how many directions the molecule can actually move in, and by how
# that motion is foreshortened when projected onto the image plane.
#
#   motion   tracked        d_eff   why
#   volume   3D              3      3 DOF on 3 axes - the plain case
#   volume   2D              2      dropping z leaves 2 independent axes
#   surface  3D              2      2 DOF however many coordinates are recorded
#   surface  2D              1 + ⟨n_z²⟩, geometry-dependent (below)
#
# The last row is the subtle one: a molecule on a curved surface steps in its local
# tangent plane, and the fraction of that step that survives projection depends on
# how the surface is tilted. Averaging over the shape gives ⟨n_z²⟩ = 1/3 for a
# sphere, 1/2 for a cylinder whose axis lies in the image plane, and 1 for a flat
# membrane parallel to it - so d_eff = 4/3, 3/2 and 2. Verified by simulation to
# three decimals. A membrane tilted by θ gives 1 + cos²θ, between the rod and flat
# values; it is not offered separately.
#
# Getting this wrong is a first-order error on D, not a refinement: using the
# recorded dimensionality for a surface-bound molecule under-reports D by 1/3 when
# tracked in 3D, and by 1/3 (sphere) or 1/4 (rod) when tracked in 2D.
# Extended for 1-DOF motion (a motor on a filament). Recorded in 3D the answer is
# just the motion's own dimensionality; recorded in 2D it is reduced by whatever
# fraction of a step points out of the image plane, which is ⟨n_z²⟩ (surface) or
# ⟨u_z²⟩ (filament) for the shape in use:
#   surface  on a sphere    ⟨n_z²⟩ = 1/3 → d_eff = 4/3
#            on a rod       ⟨n_z²⟩ = 1/2 → d_eff = 3/2   (axis in the image plane)
#            on a membrane  ⟨n_z²⟩ = 1   → d_eff = 2     (tilt θ gives 1 + cos²θ)
#   filament in-plane                    → d_eff = 1
# An unspecified shape ("none") is taken to lie in the image plane, the common
# flat-mount case and the least surprising default - it applies no reduction.
# The simulated free-motion nulls understand closed/flat surfaces only. "None" is
# genuinely unbounded, and a filament has no surface to diffuse on, so both fall back
# to the unbounded case - the null then reduces to the infinite-plane one, which is
# the right answer rather than a fudge.
.geom_for_null <- function(g) if (g %in% c("sphere", "rod", "ellipsoid", "membrane")) g else "free"


.msd_d_eff <- function(d_obs) {
  motion <- input$trk_motion %||% "volume"
  if (identical(motion, "volume")) return(as.numeric(d_obs))
  d_mot <- if (identical(motion, "filament")) 1 else 2
  if (d_obs >= 3L) return(d_mot)               # all coordinates recorded: no projection loss
  if (identical(motion, "filament")) return(1) # a track lying in the image plane
  switch(input$trk_geom_type %||% "none",
         sphere = 4 / 3, ellipsoid = 4 / 3, rod = 3 / 2, membrane = 2, 2)
}

msd_unified_bundle <- eventReactive(input$trk_msd_compute, {
  req(isTRUE(input$trk_show_msd), isolate(isTRUE(input$trk_msd_unified)))
  t_c <- isolate(trk_cycle_time())
  cc_sw <- isolate({ cc <- tryCatch(trk_cluster_cycle(), error = function(e) NULL)
                     if (!is.null(cc) && is.finite(cc$sw)) cc$sw else 0 })
  shiny::validate(shiny::need(is.finite(t_c) && t_c > 0, "Could not determine cycle time."))
  # Number of cycle-lags to compute = the requested upper fit limit (so the
  # confined-model plateau is reachable), floored at 60 and capped at 5000 to
  # avoid runaway. Re-click Plot MSD after changing the upper limit to recompute.
  max_lag <- max(60L, min(5000L, as.integer(isolate(input$trk_uni_max_lag) %||% 50L)))
  min_pts <- isolate(max(5L, as.integer(input$trk_msd_min_pts)))

  df <- isolate(trk()); req(nrow(df) > 0)
  # Dimensionality: honour the user's 2D/3D selector (auto = 3D only if z varies).
  ndim_sel   <- isolate(input$trk_msd_ndim %||% "auto")
  has_z_auto <- .has_real_z(df)
  has_z <- if (identical(ndim_sel, "2")) FALSE
           else if (identical(ndim_sel, "3")) ("z" %in% names(df))
           else has_z_auto
  df <- df[is.finite(df$x) & is.finite(df$y) & is.finite(df$tim) &
             if (has_z) is.finite(df$z) else TRUE, ]
  # Per-track sums at every cycle lag (msd_track_sums, est_msd.R): the lag is the
  # difference in cumulative cycle index (Bela's TOT), robust to cycle-time error,
  # and the real dt is kept for the mean-dt-per-bin axis. Vectorised over index
  # offsets; it replaced a loop over every localization that took about a minute
  # on a single cell.
  sums <- withProgress(message = "Computing MSD…", value = 0.5,
                       msd_track_sums(df, t_c, cc_sw, max_lag, min_pts, has_z))
  tk_names <- sums$tk_names
  tk_dr2 <- sums$dr2; tk_np <- sums$np; tk_dt <- sums$dt

  track_msd <- ifelse(tk_np > 0L, tk_dr2 / tk_np, NA_real_)   # per-track TAMSD (nm²)

  # Drop non-positive-D tracks before the ensemble. OFF by default to match
  # Bela's eDIFF, which keeps EVERY track (his D≤0 filter applies only to the
  # per-track tDIFF, which we still enforce there). Toggle on for the previous
  # QC behaviour.
  n_dropped_negD <- 0L
  if (isTRUE(input$trk_uni_drop_negd)) {
    posD <- .track_positive_D(track_msd, seq_len(max_lag) * t_c)
    n_dropped_negD <- if (any(posD) && !all(posD)) sum(!posD) else 0L
    if (any(posD) && !all(posD)) {
      track_msd <- track_msd[, posD, drop = FALSE]
      tk_dr2    <- tk_dr2[, posD, drop = FALSE]
      tk_np     <- tk_np[, posD, drop = FALSE]
      tk_dt     <- tk_dt[, posD, drop = FALSE]
      tk_names  <- tk_names[posD]
    }
  }
  colnames(track_msd) <- tk_names

  np_lag    <- rowSums(tk_np)
  dr2_lag   <- rowSums(tk_dr2)
  dt_lag    <- rowSums(tk_dt)
  # Ensemble weighting: pair-weighted (default, smoother) or track-weighted
  # (Vogler/Bela eMSD = mean of per-track TAMSDs).
  track_weight <- isolate(isTRUE(input$trk_uni_track_weight))
  msd_pair  <- if (track_weight) rowMeans(track_msd, na.rm = TRUE)
               else ifelse(np_lag > 0, dr2_lag / np_lag, NA_real_)
  tau_real  <- ifelse(np_lag > 0, dt_lag  / np_lag, seq_len(max_lag) * t_c)  # true mean Δt
  n_trk_lag <- rowSums(!is.na(track_msd))
  msd_q25   <- apply(track_msd, 1L, quantile, probs = 0.25, na.rm = TRUE)
  msd_q75   <- apply(track_msd, 1L, quantile, probs = 0.75, na.rm = TRUE)
  # Median of the same per-track TAMSDs, so the plotted band has a central line
  # that is commensurable with it (the eMSD is a mean of a skewed distribution).
  msd_med   <- apply(track_msd, 1L, stats::median, na.rm = TRUE)

  min_trk <- max(2L, floor(0.10 * max(n_trk_lag)))
  keep_lag <- n_trk_lag >= min_trk
  # Linearized (Vogler/MinSpt) time axis = k·t_c + SW, matching his T_LIN exactly
  # (simulate_cycles(k) = k·t_c + SW_OVERHEAD; get_cycle_time subtracts the SW).
  # The +SW is a CONSTANT offset, so it does NOT change the fitted slope / D - it
  # only shifts the curve right by SW and aligns the OLSF intercept (σ_loc) with his.
  list(lag = seq_len(max_lag), tau = tau_real, tau_lin = seq_len(max_lag) * t_c + cc_sw,
       msd = msd_pair,
       npairs = np_lag, n_trk_lag = n_trk_lag, keep_lag = keep_lag,
       track_msd = track_msd, tk_names = tk_names,
       msd_q25 = msd_q25, msd_q75 = msd_q75, msd_med = msd_med,
       t_cycle = t_c, has_z = has_z, ndim = if (has_z) 3L else 2L,
       # MSD = 2*d_eff*D*tau. d_eff follows the MOTION and its projection, not the
       # recorded coordinate count -- see .msd_d_eff.
       d_eff = .msd_d_eff(if (has_z) 3L else 2L),
       coef_eff = 2 * .msd_d_eff(if (has_z) 3L else 2L),
       track_weighted = track_weight, n_dropped_negD = n_dropped_negD,
       n_tracks = sum(colSums(tk_np) > 0))
})

# ── Free-diffusion-on-a-geometry MSD reference ──────────────────────────────
# On a closed cell the MSD of a FREE molecule bends and eventually saturates, so a
# bend in the measured curve is not by itself evidence of trapping. This simulates
# free diffusion on the chosen shape at the fitted D and returns the curve it would
# produce, on the same lag grid and the same observed coordinates as the data.
#
# The geometry is taken from the Displacement vs Diffusion panel rather than
# duplicated here - one cell definition, used by both tests.
#
# D fed in is the FITTED D from this panel, which is already corrected for motion
# degrees of freedom and projection (see .msd_d_eff), i.e. a true surface D. That is
# why no projection re-calibration is applied, unlike the net-displacement null which
# receives a raw per-axis D.
msd_geom_null <- reactive({
  req(isTRUE(input$trk_msd_geom))
  ub <- tryCatch(msd_unified_bundle(), error = function(e) NULL); req(!is.null(ub))
  fr <- tryCatch(msd_unified_fit(ub), error = function(e) NULL)
  req(!is.null(fr), is.finite(fr$D), fr$D > 0)
  geom <- .geom_for_null(input$trk_geom_type %||% "none")
  params <- list(
    R    = max(1, as.numeric(input$trk_disp_R    %||% 500)),
    r    = max(1, as.numeric(input$trk_disp_r    %||% 400)),
    Lcyl = max(1, as.numeric(input$trk_disp_L    %||% 2000)),
    a    = max(1, as.numeric(input$trk_disp_a    %||% 500)),
    b    = max(1, as.numeric(input$trk_disp_b    %||% 500)),
    c    = max(1, as.numeric(input$trk_disp_c    %||% 900)),
    mtilt = max(0, as.numeric(input$trk_disp_mtilt %||% 0)))
  # Anything that is not volume-filling sits on the boundary (same mapping the
  # displacement null uses), so filament and surface both simulate on the surface.
  mode <- if (identical(input$trk_motion %||% "volume", "volume")) "volume" else "surface"
  keep <- ub$keep_lag & is.finite(ub$msd)
  lags <- ub$lag[keep]
  if (!length(lags)) return(NULL)
  # Lag grid: log-spaced across the whole range so the overlay reads well on any
  # scale, PLUS an even spread inside the fit window. The fit-window points are what
  # the D-shortfall estimate below is regressed on, and a purely log grid leaves that
  # regression resting on a handful of long lags, which made it noisy.
  lo_w <- max(1L, as.integer(input$trk_uni_min_lag %||% 4L))
  hi_w <- max(lo_w + 2L, as.integer(input$trk_uni_max_lag %||% 50L))
  hi_w <- min(hi_w, max(lags))
  grid <- unique(c(
    round(exp(seq(log(1), log(max(lags)), length.out = 40))),
    round(seq(lo_w, hi_w, length.out = min(30L, max(3L, hi_w - lo_w + 1L))))))
  grid <- sort(grid[grid >= 1L])
  # Cost scales with the LARGEST lag, not the number of points, so keep the ensemble
  # large - the long-lag points are the noisy ones and they carry the signal.
  M <- 4000L
  res <- withProgress(message = "Simulating free diffusion on the geometry…", value = 0.4,
    sim_free_geom_msd(grid, D_surf = fr$D, geom = geom, params = params,
                      cyc_ms = ub$t_cycle, proj2d = (ub$ndim == 2L),
                      mode = mode, M = M))
  if (is.null(res) || !nrow(res)) return(NULL)
  res$tau_lin <- res$lag * ub$t_cycle
  attr(res, "geom") <- geom; attr(res, "mode") <- mode
  attr(res, "dims") <- .disp_geom_dims(geom, params)
  attr(res, "D") <- fr$D
  attr(res, "coef") <- ub$coef_eff %||% (2 * ub$ndim)
  res
})

# Fit helper shared by the unified plot + estimates (OLSF + Michalet p_min)
# The MSD panel's fit settings, for computing the same upper-regime D anywhere
# (msd_headline, est_msd.R): Data Quality's headline D and Compare Cells.
.msd_settings <- function(df) {
  ndim_sel <- input$trk_msd_ndim %||% "auto"
  has_z <- if (identical(ndim_sel, "2")) FALSE
           else if (identical(ndim_sel, "3")) ("z" %in% names(df))
           else .has_real_z(df)
  lo <- max(1L, as.integer(input$trk_uni_min_lag %||% 5L))
  hi <- max(lo + 2L, as.integer(input$trk_uni_max_lag %||% 50L))
  list(window = c(lo, hi), track_weight = isTRUE(input$trk_uni_track_weight %||% TRUE),
       iterative = isTRUE(input$trk_uni_iterative %||% TRUE),
       linear_tau = isTRUE(input$trk_uni_linear_tau %||% TRUE),
       min_pts = max(5L, as.integer(input$trk_msd_min_pts %||% 50L)),
       R = tryCatch(trk_blur()$R_msd, error = function(e) 0),
       coef = c(d2 = 2 * .msd_d_eff(2L), d3 = 2 * .msd_d_eff(3L)), has_z = has_z)
}

msd_unified_fit <- function(ub) {
  min_lag <- max(1L, as.integer(isolate(input$trk_uni_min_lag) %||% 4L))
  # Upper fit limit (Vogler fit_limits, default 50). Read WITHOUT isolate so the
  # plot/estimates that call this re-render when the user changes it.
  max_lag <- max(min_lag + 2L, as.integer(input$trk_uni_max_lag %||% 50L))
  tau_use <- if (isTRUE(isolate(input$trk_uni_linear_tau))) ub$tau_lin else ub$tau
  ok  <- is.finite(tau_use) & is.finite(ub$msd) & ub$keep_lag &
         ub$lag >= min_lag & ub$lag <= max_lag
  if (sum(ok) < 3L) return(NULL)
  if (isTRUE(input$trk_uni_iterative)) {
    # Bela's exact iterative loop - converges on the lowest, least-noisy lags.
    msd_fit_olsf_iterative(tau_use[ok], ub$msd[ok], ndim = ub$ndim, R = trk_blur()$R_msd,
                           dt_cycle = ub$t_cycle, N = sum(ok),
                           coef_override = ub$coef_eff)
  } else {
    # Robust full-window-anchored Michalet (default).
    msd_fit_olsf(tau_use[ok], ub$msd[ok], ub$npairs[ok],
                 k = sum(ok), ndim = ub$ndim,
                 auto_pmin = TRUE, R = trk_blur()$R_msd, Npts = 1e4, dt_cycle = ub$t_cycle,
                 coef_override = ub$coef_eff)
  }
}

# .msd_model_fit() (linear / confined / directed / anomalous) lives in est_msd.R.

# Ensemble-MSD model fit for the unified plot (selectable model). Non-linear
# models fit the whole [min_lag, max_lag] window (they need curvature/plateau).
msd_unified_modelfit <- function(ub) {
  model <- input$trk_uni_model %||% "linear"
  if (identical(model, "linear")) return(msd_unified_fit(ub))
  min_lag <- max(1L, as.integer(input$trk_uni_min_lag %||% 4L))
  max_lag <- max(min_lag + 2L, as.integer(input$trk_uni_max_lag %||% 50L))
  tau_use <- if (isTRUE(input$trk_uni_linear_tau)) ub$tau_lin else ub$tau
  ok <- is.finite(tau_use) & is.finite(ub$msd) & ub$keep_lag &
        ub$lag >= min_lag & ub$lag <= max_lag
  if (sum(ok) < 4L) return(NULL)
  .msd_model_fit(tau_use[ok], ub$msd[ok], ub$npairs[ok], ub$ndim, model,
                 isTRUE(input$trk_uni_iterative), ub$t_cycle,
                 coef_override = ub$coef_eff, R = trk_blur()$R_msd)
}

# Cached per-track analysis from the UNIFIED MSD - the single source for the
# time-averaged D (tDIFF), the per-track D histogram, the D/CI clustering and the
# per-track CSV export. Runs the same per-track OLSF as the ensemble (iterative
# when toggled, else anchored) on every track's TAMSD over the fit window, and
# computes a confinement index (CI) from the same windowed curve. Cached so the
# heavy loop only re-runs when the bundle or the window/method changes.
# Returns list(tDIFF µm²/s, n_tD, per_D, pertrack = data.frame).
#
# WHY TWO DIFFUSION NUMBERS - eDIFF vs tDIFF (they answer different questions)
#   eDIFF ("ensemble")   fit ONE line to the pooled ensemble MSD curve.
#                        -> the population's typical D. Computed by
#                           msd_unified_fit() on msd_unified_bundle()$msd.
#   tDIFF ("time-average") fit EACH track separately, then average the D values.
#                        -> the average individual molecule. Computed here.
#
#   For a homogeneous, ergodic sample the two agree. When they DISAGREE that is
#   itself the result - it means the population is heterogeneous or
#   non-ergodic (a few fast molecules pull the ensemble fit, while the typical
#   molecule is slower, or vice versa). The ratio is reported as
#
#        ergodicity = <D>_time-average / D_ensemble      (1.0 = ergodic)
#
#   Only tracks with D > 0 enter tDIFF (a fitted negative slope is noise, not a
#   molecule diffusing backwards).
#
#   CONFINEMENT INDEX (CI): the ratio of the MSD curve's late-window rise to
#   its early-window rise. Free diffusion rises linearly so CI ~ 1; a confined
#   molecule's MSD flattens towards a plateau, so CI < 1.
trk_unified_tdiff <- reactive({
  ub <- tryCatch(msd_unified_bundle(), error = function(e) NULL)
  if (is.null(ub) || is.null(ub$track_msd) || ncol(ub$track_msd) == 0L) return(NULL)
  min_lag   <- max(1L, as.integer(input$trk_uni_min_lag %||% 4L))
  max_lag_e <- max(min_lag + 2L, as.integer(input$trk_uni_max_lag %||% 50L))
  use_iter  <- isTRUE(input$trk_uni_iterative)
  tau_sel   <- if (isTRUE(input$trk_uni_linear_tau)) ub$tau_lin else ub$tau
  win   <- ub$lag >= min_lag & ub$lag <= max_lag_e & is.finite(tau_sel)
  tau_w <- tau_sel[win]; t_c_e <- ub$t_cycle
  tids  <- if (!is.null(ub$tk_names)) ub$tk_names else as.character(seq_len(ncol(ub$track_msd)))
  rows <- lapply(seq_len(ncol(ub$track_msd)), function(j) {
    col <- ub$track_msd[, j]
    mv  <- col[win]
    okc <- is.finite(mv); if (sum(okc) < 3L) return(NULL)
    tt <- tau_w[okc]; yy <- mv[okc]; N_tr <- sum(is.finite(col))
    fr_tr <- if (use_iter)
      msd_fit_olsf_iterative(tt, yy, ndim = ub$ndim, R = trk_blur()$R_msd, dt_cycle = t_c_e, N = N_tr,
                             coef_override = ub$coef_eff)
    else
      msd_fit_olsf(tt, yy, rep(1L, length(tt)), k = length(tt), ndim = ub$ndim,
                   auto_pmin = TRUE, R = trk_blur()$R_msd, Npts = N_tr, dt_cycle = t_c_e,
                   coef_override = ub$coef_eff)
    if (is.null(fr_tr) || !is.finite(fr_tr$D) || fr_tr$D <= 0) return(NULL)
    # Confinement index: actual MSD growth vs the free-diffusion expectation,
    # using the median of the first/last 15% of the windowed bins (robust to a
    # noisy endpoint). CI < 1 confined, ~1 Brownian, > 1 directed.
    n_ci <- length(tt); n_edge <- max(2L, floor(0.15 * n_ci))
    tf <- median(tt[seq_len(n_edge)]); tl <- median(tt[seq(n_ci - n_edge + 1L, n_ci)])
    mf <- median(yy[seq_len(n_edge)]); ml <- median(yy[seq(n_ci - n_edge + 1L, n_ci)])
    ci <- if (is.finite(tf) && tf > 0 && is.finite(mf) && mf > 0 &&
              is.finite(tl) && tl > tf && is.finite(ml)) (ml / mf) / (tl / tf) else NA_real_
    data.frame(tid = tids[j], D_um2_s = fr_tr$D / 1e3, alpha = 1,
               sigma_loc_nm = fr_tr$sigma_loc, conf_index = ci,
               n_bins = sum(okc), stringsAsFactors = FALSE)
  })
  pertrack <- do.call(rbind, Filter(Negate(is.null), rows))
  per_D <- if (!is.null(pertrack)) pertrack$D_um2_s else numeric(0)
  ci <- .mean_ci(per_D)   # track-level SD + 95% CI on the per-track-mean D
  list(tDIFF = ci$mean, n_tD = ci$n, sd_tD = ci$sd, ci_lo = ci$lo, ci_hi = ci$hi,
       med_tD = ci$median, per_D = per_D, pertrack = pertrack)
})

# Per-cluster diffusion summary (cached): eDIFF from the cluster eMSD fit (same
# window + method as the overall ensemble) and tDIFF from the cluster's per-track
# D's. Used by BOTH the unified plot legend and the estimates breakdown so the
# cluster eMSD is fit only once. Returns a data.frame(cluster, eDIFF, tDIFF, n,
# color) or NULL.
trk_unified_cluster_diff <- reactive({
  if (!isTRUE(input$trk_msd_cluster)) return(NULL)
  ub <- tryCatch(msd_unified_bundle(), error = function(e) NULL)
  cl <- tryCatch(msd_cluster_rv(),    error = function(e) NULL)
  td <- tryCatch(trk_unified_tdiff(), error = function(e) NULL)
  if (is.null(ub) || is.null(cl) || is.null(ub$track_msd)) return(NULL)
  tids_u    <- if (!is.null(ub$tk_names)) ub$tk_names else as.character(seq_len(ncol(ub$track_msd)))
  cids      <- cl$cluster_ids[tids_u]
  min_lag   <- max(1L, as.integer(input$trk_uni_min_lag %||% 4L))
  max_lag_c <- max(min_lag + 2L, as.integer(input$trk_uni_max_lag %||% 50L))
  tau_sel   <- if (isTRUE(input$trk_uni_linear_tau)) ub$tau_lin else ub$tau
  use_iter  <- isTRUE(input$trk_uni_iterative)
  pt <- if (!is.null(td)) td$pertrack else NULL
  do.call(rbind, lapply(seq_len(cl$n_clusters), function(clk) {
    cols_clk <- which(!is.na(cids) & cids == clk)
    if (!length(cols_clk)) return(NULL)
    emsd <- rowMeans(ub$track_msd[, cols_clk, drop = FALSE], na.rm = TRUE)  # nm²
    winc <- ub$lag >= min_lag & ub$lag <= max_lag_c & is.finite(tau_sel) & is.finite(emsd)
    eD <- if (sum(winc) >= 3L) {
      frc <- tryCatch(
        if (use_iter) msd_fit_olsf_iterative(tau_sel[winc], emsd[winc], ndim = ub$ndim,
                                             R = trk_blur()$R_msd, dt_cycle = ub$t_cycle, N = sum(winc),
                                             coef_override = ub$coef_eff)
        else msd_fit_olsf(tau_sel[winc], emsd[winc], rep(1L, sum(winc)), k = sum(winc),
                          ndim = ub$ndim, auto_pmin = TRUE, R = trk_blur()$R_msd, Npts = 1e4,
                          dt_cycle = ub$t_cycle, coef_override = ub$coef_eff),
        error = function(e) NULL)
      if (!is.null(frc) && is.finite(frc$D)) frc$D / 1e3 else NA_real_
    } else NA_real_
    tD <- if (!is.null(pt)) {
      dv <- pt$D_um2_s[!is.na(cl$cluster_ids[pt$tid]) & cl$cluster_ids[pt$tid] == clk]
      if (length(dv)) mean(dv, na.rm = TRUE) else NA_real_
    } else NA_real_
    data.frame(cluster = clk, eDIFF = eD, tDIFF = tD, n = length(cols_clk),
               color = cl$palette[clk], stringsAsFactors = FALSE)
  }))
})

output$trk_msd_unified_plot <- renderPlotly({
  req(isTRUE(input$trk_msd_unified))
  ub <- msd_unified_bundle()
  shiny::validate(shiny::need(!is.null(ub), "Click 'Plot MSD' to compute."))
  dark <- isTRUE(isolate(input$trk_dark_mode))
  bg <- if (dark) "#000000" else "#FFFFFF"; fg <- if (dark) "#FFFFFF" else "#333333"
  show_band <- isTRUE(input$trk_uni_show_band); show_ntau <- isTRUE(input$trk_uni_show_ntau)
  min_lag <- max(1L, as.integer(input$trk_uni_min_lag %||% 4L))
  linear_tau <- isTRUE(input$trk_uni_linear_tau)
  tau_full <- if (linear_tau) ub$tau_lin else ub$tau

  ok <- is.finite(tau_full) & is.finite(ub$msd) & ub$keep_lag
  tau <- tau_full[ok]; msd <- ub$msd[ok] / 1e6   # µm²
  fr  <- msd_unified_modelfit(ub)
  # Regime-cutoff x-position: exact on the linearized axis; the real τ at the
  # cutoff lag on the mean-Δt axis.
  regime_x <- if (linear_tau) min_lag * ub$t_cycle else {
    idx <- which(ub$lag == min_lag)
    if (length(idx) && is.finite(ub$tau[idx])) ub$tau[idx] else min_lag * ub$t_cycle
  }

  fig <- plot_ly()

  # Cluster map (tid -> cluster / colour) when clustering is enabled.
  cl <- if (isTRUE(input$trk_msd_cluster)) tryCatch(msd_cluster_rv(), error = function(e) NULL) else NULL
  tids_u <- if (!is.null(ub$tk_names)) ub$tk_names else as.character(seq_len(ncol(ub$track_msd)))

  # Per-track TAMSD spaghetti (optional, drawn first = behind). Coloured by
  # cluster when clustering is on, else faint grey. Sampled for performance.
  if (isTRUE(input$trk_uni_show_spag) && !is.null(ub$track_msd) && ncol(ub$track_msd) > 0) {
    tm <- ub$track_msd
    ncols <- ncol(tm); max_lines <- 200L
    set.seed(42L)
    idx_s <- if (ncols > max_lines) sort(sample(ncols, max_lines)) else seq_len(ncols)
    if (!is.null(cl)) {
      cids <- cl$cluster_ids[tids_u]      # cluster per column (NA if not fit)
      for (clk in seq_len(cl$n_clusters)) {
        sx <- numeric(0); sy <- numeric(0); ncl <- 0L
        for (j in idx_s) {
          if (is.na(cids[j]) || cids[j] != clk) next
          yv <- tm[ok, j] / 1e6; okj <- is.finite(yv)
          if (sum(okj) < 2L) next
          sx <- c(sx, tau[okj], NA_real_); sy <- c(sy, yv[okj], NA_real_); ncl <- ncl + 1L
        }
        if (length(sx))
          fig <- fig %>% add_trace(x = sx, y = sy, type = "scattergl", mode = "lines",
            line = list(color = hex_to_rgba(cl$palette[clk], 0.22), width = 1),
            name = sprintf("Cluster %d tracks (n=%d)", clk, ncl),
            hoverinfo = "skip", showlegend = TRUE)
      }
    } else {
      sx <- numeric(0); sy <- numeric(0)
      for (j in idx_s) {
        yv <- tm[ok, j] / 1e6; okj <- is.finite(yv)
        if (sum(okj) < 2L) next
        sx <- c(sx, tau[okj], NA_real_); sy <- c(sy, yv[okj], NA_real_)
      }
      if (length(sx))
        fig <- fig %>% add_trace(x = sx, y = sy, type = "scattergl", mode = "lines",
          line = list(color = "rgba(140,140,140,0.22)", width = 1),
          name = "Per-track TAMSD", hoverinfo = "skip", showlegend = TRUE)
    }
  }

  if (show_band) {
    q25 <- ub$msd_q25[ok] / 1e6; q75 <- ub$msd_q75[ok] / 1e6
    okq <- is.finite(q25) & is.finite(q75)
    if (any(okq)) fig <- fig %>% add_trace(
      x = c(tau[okq], rev(tau[okq])), y = c(q75[okq], rev(q25[okq])),
      # mode = "lines" with a transparent line, NOT mode = "none": passing a line
      # object with mode="none" makes plotly print "A line object has been specified,
      # but lines is not in the mode / Adding lines to the mode..." to the console.
      type = "scatter", mode = "lines", fill = "toself",
      fillcolor = "rgba(70,130,200,0.18)",
      line = list(color = "rgba(0,0,0,0)", width = 0),
      # Named for what it is: the spread of per-track TAMSDs, NOT a confidence
      # band on the plotted ensemble value. The line is a MEAN of that same
      # distribution, and per-track TAMSD is strongly right-skewed, so the mean
      # sits high in the band - near Q3 - with Q1 far below it.
      hoverinfo = "skip", name = "IQR across tracks", showlegend = TRUE)
  }
  fig <- fig %>% add_trace(x = tau, y = msd, type = "scatter", mode = "lines+markers",
      marker = list(size = 5, color = "steelblue"), line = list(color = "steelblue", width = 2),
      # Name the curve for what it actually is: eMSD = mean of the per-track TAMSDs
      # (track-weighted), otherwise the pair-weighted ensemble average.
      name = if (isTRUE(ub$track_weighted)) "eMSD (track-weighted)"
             else "Ensemble MSD (pair-weighted)")
  # Median of the same per-track TAMSDs. The eMSD above is a MEAN of a strongly
  # right-skewed distribution, so it sits high inside the IQR band - near Q3, with
  # Q1 far below. Without the median next to it the line looks like it is sitting
  # outside its own band. Legend-toggleable, like every other trace.
  if (show_band && !is.null(ub$msd_med)) {
    mm <- ub$msd_med[ok] / 1e6
    okm <- is.finite(mm)
    if (any(okm)) fig <- fig %>% add_trace(
      x = tau[okm], y = mm[okm], type = "scatter", mode = "lines",
      line = list(color = "rgba(70,130,200,0.95)", width = 1.6, dash = "dot"),
      name = "median track")
  }
  # Free-diffusion-on-the-cell reference. Drawn in the same green as the geometry
  # band in Displacement vs Diffusion, since it is the same idea applied to the MSD.
  gm <- tryCatch(msd_geom_null(), error = function(e) NULL)
  if (!is.null(gm) && nrow(gm)) {
    gx <- if (linear_tau) gm$tau_lin else gm$tau
    fig <- fig %>% add_trace(
      x = gx, y = gm$msd / 1e6, type = "scatter", mode = "lines",
      line = list(color = "#1B9E5A", width = 2, dash = "dash"),
      name = sprintf("free on %s (%s)", attr(gm, "geom"), attr(gm, "mode")),
      hovertemplate = "free-on-geometry<br>τ %{x:.4g} ms<br>MSD %{y:.4g} µm²<extra></extra>")
  }
  if (!is.null(fr)) {
    # Model fit curve (from fr$tau_fit / fr$msd_fit; same for every model).
    fx <- if (!is.null(fr$tau_fit)) fr$tau_fit else seq(regime_x, max(tau), length.out = 100)
    fy <- if (!is.null(fr$msd_fit)) pmax(0, fr$msd_fit) / 1e6 else {
      coef_v <- 2L * ub$ndim
      pmax(0, coef_v * fr$D * fx + coef_v * fr$sigma_loc^2) / 1e6
    }
    fit_name <- switch(fr$model %||% "linear",
      confined  = if (isTRUE(fr$degenerate))
                    sprintf("Confined: no plateau in window (≈linear, D=%.3g)", fr$D / 1e3)
                  else sprintf("Confined: D=%.3g µm²/s, L_c=%.0f nm", fr$D / 1e3, fr$L_conf),
      directed  = sprintf("Directed: D=%.3g µm²/s, v=%.3g µm/s", fr$D / 1e3, fr$v_um_s),
      anomalous = sprintf("Anomalous: D=%.3g, α=%.2f", fr$D / 1e3, fr$alpha),
      sprintf("OLSF D=%.3g µm²/s (%d pts)", fr$D / 1e3, fr$n_points %||% NA))
    # Draw the fit SOLID only across the lags it was fitted to, and faint beyond.
    # OLSF picks an optimal (usually small) number of low lags, but the fitted
    # curve was being drawn across the whole axis - on a bending MSD that puts the
    # line far above the data at long lags and reads as a failed fit rather than an
    # extrapolation. Splitting the trace makes the fitted range explicit.
    # OLSF fits the first n_points of the WINDOWED curve, i.e. lags starting at the
    # lower cutoff (regime_x in tau), so the fitted range runs from there. The
    # non-linear models use the whole window, so they fall back to the upper limit.
    fit_hi <- {
      v_lag  <- suppressWarnings(as.numeric(input$trk_uni_max_lag))
      tau_hi <- if (is.finite(v_lag) && v_lag > 0) v_lag * ub$t_cycle else max(fx, na.rm = TRUE)
      n_used <- suppressWarnings(as.integer(fr$n_points %||% NA))
      # The fitted range ends at the n_used-th point OF THE WINDOWED CURVE, so read
      # it off the tau grid itself. Adding n_used * t_cycle to the cutoff is only
      # right on the linearized axis; on the mean-Delta-t axis the bins are spaced
      # wider than one cycle (they carry dark cycles), so that arithmetic lands
      # short and the solid segment renders far shorter than the range it stands for.
      tw <- sort(tau[is.finite(tau)])
      if (is.finite(n_used) && n_used > 0 && length(tw))
        min(tau_hi, tw[min(n_used, length(tw))])
      else tau_hi
    }
    in_fit <- is.finite(fx) & fx <= fit_hi
    # One curve, two segments: solid-dashed over the lags actually fitted, lighter
    # dotted beyond. Extrapolating an OLSF fit past its optimal point count is what
    # the method warns against, so the split stays - but both halves share a legend
    # group, so they read (and toggle) as a single trace.
    if (any(in_fit)) fig <- fig %>% add_trace(
      x = fx[in_fit], y = fy[in_fit], type = "scatter", mode = "lines",
      legendgroup = "msdfit",
      line = list(color = "#E03020", width = 2.5, dash = "dash"), name = fit_name)
    if (any(!in_fit)) fig <- fig %>% add_trace(
      x = fx[!in_fit], y = fy[!in_fit], type = "scatter", mode = "lines",
      legendgroup = "msdfit", showlegend = FALSE,
      line = list(color = "rgba(224,48,32,0.55)", width = 2, dash = "dot"),
      name = "fit, extrapolated",
      hovertemplate = "extrapolated beyond fitted lags<br>tau %{x:.4g} ms<br>MSD %{y:.4g} um2<extra></extra>")
  }

  # Per-cluster eMSD: track-weighted mean TAMSD of each cluster's tracks, labelled
  # with the cluster eDIFF + tDIFF (from the cached trk_unified_cluster_diff).
  if (!is.null(cl) && !is.null(ub$track_msd)) {
    cdiff <- tryCatch(trk_unified_cluster_diff(), error = function(e) NULL)
    cids  <- cl$cluster_ids[tids_u]
    for (clk in seq_len(cl$n_clusters)) {
      cols_clk <- which(!is.na(cids) & cids == clk)
      if (length(cols_clk) < 1L) next
      emsd_clk <- rowMeans(ub$track_msd[, cols_clk, drop = FALSE], na.rm = TRUE)[ok] / 1e6
      okc <- is.finite(emsd_clk); if (sum(okc) < 2L) next
      eD <- tD <- NA_real_
      if (!is.null(cdiff)) { r <- cdiff[cdiff$cluster == clk, ]; if (nrow(r)) { eD <- r$eDIFF[1]; tD <- r$tDIFF[1] } }
      nm <- if (is.finite(eD) || is.finite(tD))
              sprintf("Cluster %d eMSD (eDIFF %s, tDIFF %s)", clk,
                      if (is.finite(eD)) sprintf("%.3g", eD) else "-",
                      if (is.finite(tD)) sprintf("%.3g", tD) else "-")
            else sprintf("Cluster %d eMSD", clk)
      fig <- fig %>% add_trace(x = tau[okc], y = emsd_clk[okc], type = "scatter", mode = "lines",
        line = list(color = cl$palette[clk], width = 2.8), name = nm)
    }
  }
  regime <- regime_x
  # N(τ) on the right axis - same warm style as the standard plot.
  show_n2 <- show_ntau && diff(range(ub$n_trk_lag[ok])) > 0
  if (show_n2)
    fig <- fig %>% add_trace(x = tau, y = ub$n_trk_lag[ok], type = "scatter", mode = "lines",
      line = list(color = "rgba(180,100,40,0.70)", width = 1.5, dash = "dot"),
      yaxis = "y2", name = "N(τ) tracks")
  # Lower-regime shading (standard-plot style): shade τ < cutoff as the
  # noise / linker-dominated zone excluded from the fit.
  regime_shapes <- if (is.finite(regime) && regime > 0 && regime < max(tau))
    list(list(type = "rect", xref = "x", yref = "paper",
              x0 = 0, x1 = regime, y0 = 0, y1 = 1,
              fillcolor = "rgba(220,160,80,0.10)",
              line = list(color = "rgba(200,140,60,0.40)", width = 1, dash = "dot")))
    else list()
  lay <- list(xaxis = list(title = if (linear_tau)
                             "τ (ms, linearized = lag × cycle)" else
                             "τ (ms, true mean Δt per bin)", color = fg),
              yaxis = list(title = "MSD (µm²)", color = fg, rangemode = "tozero"),
              paper_bgcolor = bg, plot_bgcolor = bg, font = list(color = fg),
              # Legend BELOW the plot: the right-hand side is occupied by the N(τ)
              # axis, and plotly's default right-hand legend sits on top of it.
              # y is in paper coordinates (a fraction of the plotting area), so it
              # is derived from the element height rather than hard-coded: the gap
              # below the axes then stays put when TRK_MSD_PLOT_PX changes.
              legend = list(font = list(color = fg, size = 10),
                            orientation = "h", x = 0, xanchor = "left",
                            y = -TRK_MSD_LEGEND_GAP /
                                 max(60, TRK_MSD_PLOT_PX - TRK_MSD_MARGIN_T - TRK_MSD_MARGIN_B),
                            yanchor = "top"),
              margin = list(l = 55, r = 60, t = TRK_MSD_MARGIN_T, b = TRK_MSD_MARGIN_B),
              shapes = regime_shapes, uirevision = "msd_unified")
  if (is.finite(regime) && regime > 0)
    lay$annotations <- list(list(
      x = regime * 0.5, y = 0.97, xref = "x", yref = "paper",
      text = "lower\nregime", showarrow = FALSE,
      font = list(size = 9, color = "rgba(180,120,40,0.70)"),
      xanchor = "center", yanchor = "top"))
  if (show_n2)
    lay$yaxis2 <- list(overlaying = "y", side = "right", title = "N(τ) tracks",
                       showgrid = FALSE, rangemode = "tozero",
                       tickfont = list(color = "rgba(180,100,40,0.85)"),
                       titlefont = list(color = "rgba(180,100,40,0.85)"))
  do.call(plotly::layout, c(list(fig), lay))
})

# ── Cycle-binning diagnostic ────────────────────────────────────────────────
# Re-inspects whether displacements are being binned correctly by cycle. Reads
# the existing bundle (no recompute). Two signals per cycle-lag:
#
#   • BLUE, pairs per lag (left, log) - how much data backs each lag. Jagged at low
#     lag because dark gaps make some cycle-lags reachable more often than others;
#     then falls away, because a track of N localizations offers only ~(N−k) pairs at
#     lag k and short tracks stop contributing at all. Where the count collapses the
#     MSD goes noisy/saw-toothed - hence the "counts drop past N×" marker.
#
#   • RED, real elapsed Δt / linearized axis (right) - a rise to a PLATEAU is the
#     EXPECTED, benign shape, not a fault. Per step, real elapsed = dcyc·t_c + SW,
#     so a pair spanning k cycles in m steps takes k·t_c + m·SW, while the linearized
#     axis (Vogler's T_LIN) charges SW just once: k·t_c + SW. Hence
#
#         ratio(k) = (k·t_c + m·SW) / (k·t_c + SW)
#
#     which is exactly 1 at lag 1 (m = 1) and climbs to a ceiling
#
#         plateau = 1 + (m/k)·SW/t_c ,   m/k → fraction of cycles that start a step
#
#     i.e. ≤ 1 + SW/t_c, hit only if every step is integrated. Verified in silico
#     (t_c = 0.5, SW = 0.08, 35% dark: predicted 1.1041, observed 1.1033). So the
#     PLATEAU HEIGHT reads out the overhead, and its excess over 1 is the absolute-D
#     difference between the two axis conventions (it cancels in A-vs-B comparisons).
#     What would actually be wrong: a ratio that keeps climbing without flattening,
#     or one that saw-tooths - that is a mis-estimated t_c, not overhead.
trk_msd_bindiag <- reactive({
  ub <- msd_unified_bundle(); req(!is.null(ub))
  lag   <- ub$lag; tc <- ub$t_cycle
  npair <- ub$npairs; ntrk <- ub$n_trk_lag
  keep  <- ub$keep_lag & is.finite(npair) & npair > 0
  # Compare real elapsed against the ACTUAL linearized axis (ub$tau_lin = k·t_c + SW,
  # Vogler's T_LIN), not the bare k·t_c - that's the axis the D is fit on. The two
  # coincide at lag 1 (both carry one SW) and diverge as real time accumulates an
  # extra SW per step; that divergence is the real-vs-linearized D difference.
  ratio <- ifelse(npair > 0 & is.finite(ub$tau) & is.finite(ub$tau_lin) & ub$tau_lin > 0,
                  ub$tau / ub$tau_lin, NA_real_)
  # Trust boundary: last kept lag whose pair count is still ≥ 20% of the peak.
  pk <- if (any(keep)) max(npair[keep], na.rm = TRUE) else NA_real_
  trust_lag <- if (is.finite(pk)) {
    good <- which(keep & npair >= 0.20 * pk)
    if (length(good)) max(lag[good]) else NA_integer_
  } else NA_integer_
  # Observed ceiling of the red curve: median ratio over the top decile of kept lags.
  # Its excess over 1 is (steps per cycle)·SW/t_c - the overhead read-out.
  kl  <- lag[keep & is.finite(ratio)]
  plateau <- if (length(kl) >= 10L)
    stats::median(ratio[keep & is.finite(ratio) &
                        lag >= stats::quantile(kl, 0.90, names = FALSE)], na.rm = TRUE)
  else NA_real_
 
  # UI says so. Use the TSA / cycle-time panel to actually determine t_c.
  sw_eff  <- if (length(ub$tau_lin) && is.finite(ub$tau_lin[1]) && is.finite(tc))
               ub$tau_lin[1] - tc else NA_real_
  ceiling_ratio <- if (is.finite(sw_eff) && is.finite(tc) && tc > 0)
                     1 + sw_eff / tc else NA_real_
  tc_flag <- if (!is.finite(plateau) || !is.finite(ceiling_ratio)) NA_character_
             else if (plateau < 0.995) "too_large"
             else if (plateau > ceiling_ratio * 1.005) "too_small"
             else "ok"
  list(lag = lag, npair = npair, ntrk = ntrk, ratio = ratio, keep = keep,
       tc = tc, trust_lag = trust_lag, plateau = plateau,
       sw_eff = sw_eff, ceiling_ratio = ceiling_ratio, tc_flag = tc_flag,
       min_lag = max(1L, as.integer(input$trk_uni_min_lag %||% 4L)),
       max_lag = max(3L, as.integer(input$trk_uni_max_lag %||% 50L)))
})

output$trk_msd_bindiag_plot <- renderPlotly({
  req(isTRUE(input$trk_msd_unified), isTRUE(input$trk_msd_show_bindiag))
  d <- tryCatch(trk_msd_bindiag(), error = function(e) NULL)
  shiny::validate(shiny::need(!is.null(d), "Click 'Plot MSD' first."))
  dark <- isTRUE(isolate(input$trk_dark_mode))
  bg <- if (dark) "#000000" else "#FFFFFF"; fg <- if (dark) "#FFFFFF" else "#333333"
  k  <- d$keep
  fig <- plot_ly() %>%
    add_trace(x = d$lag[k], y = d$npair[k], type = "scatter", mode = "lines+markers",
              line = list(color = "rgba(70,110,180,0.9)", width = 2),
              marker = list(size = 3, color = "rgba(70,110,180,0.9)"),
              name = "pairs / lag", yaxis = "y") %>%
    add_trace(x = d$lag[k], y = d$ratio[k], type = "scatter", mode = "lines",
              line = list(color = "#C0392B", width = 2),
              name = "real Δt / linearized axis", yaxis = "y2")
  # Reference: real Δt == the linearized axis (ratio = 1) on the right axis.
  shp <- list(list(type = "line", x0 = min(d$lag[k]), x1 = max(d$lag[k]),
                   y0 = 1, y1 = 1, yref = "y2",
                   line = list(color = if (dark) "rgba(200,200,200,0.4)" else "rgba(90,90,90,0.4)",
                               width = 1, dash = "dot")))
  # Shade the fit window, and mark the trust boundary.
  shp[[length(shp)+1L]] <- list(type = "rect", x0 = d$min_lag, x1 = d$max_lag,
    y0 = 0, y1 = 1, yref = "paper", fillcolor = "rgba(120,160,120,0.10)",
    line = list(width = 0), layer = "below")
  ann <- list()
  if (is.finite(d$trust_lag)) {
    shp[[length(shp)+1L]] <- list(type = "line", x0 = d$trust_lag, x1 = d$trust_lag,
      y0 = 0, y1 = 1, yref = "paper",
      line = list(color = "#E67E22", width = 1.5, dash = "dash"))
    ann[[1]] <- list(x = d$trust_lag, y = 1, yref = "paper", yanchor = "bottom",
      text = sprintf("counts drop past %d×", d$trust_lag), showarrow = FALSE,
      font = list(size = 9, color = "#E67E22"))
  }
  fig %>% layout(
    xaxis = list(title = "cycle lag", color = fg),
    yaxis = list(title = "pairs / lag", type = "log", color = fg,
                 gridcolor = if (dark) "rgba(255,255,255,0.12)" else "rgba(0,0,0,0.10)"),
    yaxis2 = list(title = "real Δt / linearized axis", overlaying = "y", side = "right",
                  color = "#C0392B", zeroline = FALSE),
    paper_bgcolor = bg, plot_bgcolor = bg, font = list(color = fg),
    legend = list(orientation = "v", x = 1.08, xanchor = "left", y = 1,
                  font = list(color = fg)),
    margin = list(r = 10),
    shapes = shp, annotations = ann,
    uirevision = "msd_bindiag"
  )
})

output$trk_msd_bindiag_note <- renderUI({
  req(isTRUE(input$trk_msd_unified), isTRUE(input$trk_msd_show_bindiag))
  d <- tryCatch(trk_msd_bindiag(), error = function(e) NULL); req(!is.null(d))
  k <- d$keep & d$lag >= d$min_lag & d$lag <= d$max_lag
  mean_ratio <- if (any(k)) mean(d$ratio[k], na.rm = TRUE) else NA_real_
  bias_pct   <- if (is.finite(mean_ratio)) 100 * (mean_ratio - 1) else NA_real_
  linear_tau <- isTRUE(input$trk_uni_linear_tau)
  tags$div(style = "font-size:0.82em; color:#555; margin:4px 0;",
    if (is.finite(mean_ratio))
      tags$div(
        tags$b("Linearized vs real-Δt axis: "),
        sprintf("real Δt averages %+.1f%% vs the linearized axis (lag × cycle + SW) over the fit window. ", bias_pct),
        if (is.finite(bias_pct) && abs(bias_pct) > 1)
          tags$span(
            "This shifts absolute D by about that much and cancels in comparisons between datasets. ",
            if (linear_tau)
              "Linearized axis in use, the published convention (Vogler et al. 2025)."
            else
              "Real-\u0394t axis in use, as a cross-check; the linearized axis is the published convention.")
        else tags$span(style = "color:#3a7d3a;", "The two axes agree here.")),
    if (is.finite(d$trust_lag))
      tags$div(style = "margin-top:2px;",
        tags$b("Pair-count trust boundary: "),
        sprintf("counts stay \u2265 20%% of peak up to lag %d (\u2248 %.3g ms); keep the upper fit limit at or below it.",
                d$trust_lag, d$trust_lag * d$tc)),
    tags$div(style = "color:#888; margin-top:2px;",
      tags$b("Red curve: "),
      "real time over lag \u00d7 cycle. It should rise from 1 to a plateau between 1 and 1 + SW/t_c. ",
      if (is.finite(d$plateau))
        tags$span("Plateau here: ", tags$b(sprintf("%.3f", d$plateau)),
                  sprintf(" (+%.1f%%)", 100 * (d$plateau - 1)),
                  if (is.finite(d$ceiling_ratio))
                    sprintf(", band 1.000\u2013%.3f (SW \u2248 %.3g ms). ",
                            d$ceiling_ratio, d$sw_eff) else ". ") else NULL,
      tags$b("Cycle-time check: "),
      if (identical(d$tc_flag, "too_large"))
        tags$span(style = "color:#c03030;",
          sprintf("\u26a0 %.3f is below 1, so the cycle time is too large. Re-check it in TSA.",
                  d$plateau))
      else if (identical(d$tc_flag, "too_small"))
        tags$span(style = "color:#c03030;",
          sprintf("\u26a0 %.3f is above the %.3f ceiling, so the cycle time is too small. Re-check it in TSA.",
                  d$plateau, d$ceiling_ratio))
      else if (identical(d$tc_flag, "ok"))
        tags$span(style = "color:#3a7d3a;", "Within the band, consistent with the cycle time in use.")
      else NULL,
      " A mildly over-estimated cycle time can still land inside the band.",
      tags$br(),
      tags$b("Blue curve: "),
      "pairs per lag (log axis). Keep the upper fit limit where it is well populated.")
  )
})

# Cache the unified MSD's σ_loc (the LINEAR OLSF intercept - never a non-linear
# model's offset, which is not separable) as the FALLBACK noise floor, used only when
# neither a manual value nor the CVE estimate is available.
#
# This intercept is an extrapolation of the fit window back to τ = 0, so it absorbs
# anything that is not linear-in-τ at short lags, not just localization noise. On
# confined or otherwise bending data it therefore reads HIGH - simulated confined
# tracks with a true σ_loc of 2.4 nm returned 3.7 nm from the linear intercept. The
# CVE estimate is model-free (it reads the anti-correlation of consecutive steps,
# which is a pure noise signature) and is preferred wherever available; that ordering
# is enforced in trk_sigma_loc_eff.
#
# Require strictly > 0: a fit whose intercept collapses to exactly zero has not
# measured anything, and caching 0 would silently disable every noise correction.
observeEvent(msd_unified_bundle(), {
  ub <- tryCatch(msd_unified_bundle(), error = function(e) NULL); if (is.null(ub)) return()
  fr <- tryCatch(msd_unified_fit(ub), error = function(e) NULL)
  if (!is.null(fr) && is.finite(fr$sigma_loc) && fr$sigma_loc > 0)
    trk_sigma_loc_rv(fr$sigma_loc)
}, ignoreInit = TRUE)

output$trk_sigma_loc_status <- renderUI({
  eff <- trk_sigma_loc_eff()
  txt <- if (is.finite(eff))
           htmltools::tagList("Effective ", sigloc(),
             sprintf(" = %.2f nm (%s). Used as the noise floor everywhere.",
                     eff, trk_sigma_loc_src()))
         else htmltools::tagList(sigloc(),
                     " not set - estimate it in Analysis → Data Quality & Calibration ",
                     "(no MSD needed), or enter it manually. Noise-floor corrections stay ",
                     "inactive until then.")
  htmltools::div(style = "font-size:0.8em; color:#666; margin-top:-4px;", txt)
})

# eDIFF 95% CI by CLUSTER-BOOTSTRAP over tracks: resample track TAMSD columns with
# replacement, rebuild the (track-weighted) ensemble MSD, refit - the honest CI
# given within-track correlation. Gated behind a checkbox (B refits are slower);
# recomputes only when Plot MSD is pressed with the box ticked.
trk_unified_boot <- eventReactive(input$trk_msd_compute, {
  req(isTRUE(input$trk_uni_bootstrap))
  ub <- tryCatch(msd_unified_bundle(), error = function(e) NULL); req(!is.null(ub))
  tm <- ub$track_msd; req(!is.null(tm), ncol(tm) >= 5L)
  B  <- 300L; nc <- ncol(tm)
  boot_D <- withProgress(message = "Bootstrapping eDIFF 95% CI…", value = 0, {
    vapply(seq_len(B), function(i) {
      if (i %% 30L == 0L) incProgress(30 / B)
      ub2 <- ub
      ub2$msd <- rowMeans(tm[, sample.int(nc, nc, replace = TRUE), drop = FALSE], na.rm = TRUE)
      fr <- tryCatch(msd_unified_fit(ub2), error = function(e) NULL)
      if (is.null(fr) || !is.finite(fr$D)) NA_real_ else fr$D / 1e3   # µm²/s
    }, numeric(1))
  })
  q <- stats::quantile(boot_D, c(0.025, 0.975), names = FALSE, na.rm = TRUE)
  list(lo = q[1], hi = q[2], B = sum(is.finite(boot_D)))
})

output$trk_msd_unified_estimates <- renderUI({
  req(isTRUE(input$trk_msd_unified))
  ub <- tryCatch(msd_unified_bundle(), error = function(e) NULL); if (is.null(ub)) return(NULL)
  fr <- msd_unified_modelfit(ub)
  if (is.null(fr)) return(htmltools::div(style="color:#888;font-size:0.86em;", "Fit not available."))
  # eDIFF = fr$D is the ensemble-fit D. tDIFF = the time-averaged D (mean of
  # per-track OLSF D, D>0) - computed in the cached trk_unified_tdiff reactive so
  # the heavy per-track fit doesn't rerun on every estimates refresh.
  model <- fr$model %||% "linear"
  linear_tau <- isTRUE(input$trk_uni_linear_tau)
  td    <- trk_unified_tdiff()
  tDIFF <- if (!is.null(td)) td$tDIFF else NA_real_
  n_tD  <- if (!is.null(td)) td$n_tD  else 0L
  tD_ci <- if (!is.null(td) && is.finite(td$ci_lo))
             sprintf(" [95%% CI %.3g–%.3g, SD %.3g]", td$ci_lo, td$ci_hi, td$sd_tD) else ""
  eb    <- tryCatch(trk_unified_boot(), error = function(e) NULL)
  eD_ci <- if (!is.null(eb) && is.finite(eb$lo)) sprintf(" (95%% CI %.3g–%.3g)", eb$lo, eb$hi) else ""
  sig_eff <- tryCatch(trk_sigma_loc_eff(), error = function(e) NA_real_)   # noise floor for sanity flags
  ergo  <- if (is.finite(tDIFF) && is.finite(fr$D) && fr$D > 0) tDIFF / (fr$D / 1e3) else NA_real_
  # Say which MSD-to-D relation is in force whenever it is not the plain one, since
  # it scales D directly.
  dof_note <- {
    de <- ub$d_eff %||% ub$ndim
    if (abs(de - ub$ndim) > 1e-9)
      htmltools::div(style = "color:#555; font-size:0.86em; margin-top:2px;",
        htmltools::HTML(sprintf(
          "Restricted motion: MSD = %.4g·D·τ (d<sub>eff</sub> = %.4g) on %dD coordinates%s. Reading it as %dD motion would report D a factor %.3g too low.",
          ub$coef_eff, de, ub$ndim,
          if (ub$ndim == 2L) sprintf(", projected from a %s", .geom_label(input$trk_geom_type)) else "",
          ub$ndim, de / ub$ndim)))
    else NULL
  }
  extra <- switch(model,
    confined  = if (isTRUE(fr$degenerate))
                  "  |  \u26a0 no plateau in the window, so L_c is unconstrained. Raise the upper limit, or the population is not confined on this timescale."
                else sprintf("  |  confinement size L_c = %.0f nm%s  ·  plateau = %.3g µm²  ·  τ_c = %.3g ms",
                             fr$L_conf,
                             if (is.finite(sig_eff) && is.finite(fr$L_conf) && fr$L_conf < 2 * sig_eff)
                               # extra[] flows through sprintf into an escaped string, so no
                               # HTML here - use plain 2σ notation instead of an underscore.
                               sprintf(" ⚠ below precision (2σ = %.0f nm - unresolved)", 2 * sig_eff) else "",
                             fr$plateau_um2, fr$tau_c_ms),
    directed  = {
      vnoise <- if (is.finite(sig_eff) && is.finite(ub$t_cycle) && ub$t_cycle > 0)
                  sqrt(2 * ub$ndim) * sig_eff / ub$t_cycle else NA_real_   # nm/ms = µm/s
      sprintf("  |  velocity v = %.3g µm/s%s%s", fr$v_um_s,
              if (is.finite(vnoise) && is.finite(fr$v_um_s) && fr$v_um_s < vnoise)
                sprintf(" ⚠ at/below noise floor (%.3g µm/s - consistent with noise)", vnoise) else "",
              if (isTRUE(fr$degenerate)) " ⚠ fit did not converge - parameters not identifiable over this window" else "")
    },
    anomalous = sprintf("  |  α = %.2f (%s)%s", fr$alpha,
                        if (fr$alpha < 0.9) "sub-diffusive / confined"
                        else if (fr$alpha > 1.1) "super-diffusive / directed" else "≈ Brownian",
                        if (isTRUE(fr$degenerate)) " ⚠ fit did not converge - parameters not identifiable over this window" else ""),
    "")
  method_lbl <- switch(model,
    linear = if (identical(fr$method, "OLSF-iterative")) "iterative OLSF" else "anchored Michalet",
    confined = "confined model", directed = "directed model", anomalous = "anomalous model", "model")
  D_lbl <- if (identical(model, "linear")) "eDIFF (ensemble)" else "D (local, ensemble)"
  # How much is the CELL bending the curve inside the fit window? Compares the
  # simulated free-on-geometry MSD to a straight line at the same D over the same
  # lags: a large shortfall means the fitted D is being pulled down by the geometry,
  # not by biology, and the window should be shortened (or the bend interpreted as
  # the cell).
  geom_note <- {
    gm <- tryCatch(msd_geom_null(), error = function(e) NULL)
    if (is.null(gm) || !nrow(gm)) NULL else {
      lo <- max(1L, as.integer(input$trk_uni_min_lag %||% 4L))
      hi <- max(lo + 2L, as.integer(input$trk_uni_max_lag %||% 50L))
      w  <- gm$lag >= lo & gm$lag <= hi
      # Fit the SELECTED model to the simulated free-on-geometry curve, over the same
      # window, and report what it returns. Nothing is confined, drifting or anomalous
      # in that simulation, so whatever comes back is the shape alone - the baseline
      # the measured parameter has to beat. This generalises the D-shortfall check to
      # every model: on a closed shape the confined model reports a spurious corral
      # and the anomalous model a sub-diffusive exponent, purely from the ceiling.
      fg <- if (sum(w) < 4L) NULL else
        tryCatch(.msd_model_fit(gm$tau[w], gm$msd[w], rep(1000L, sum(w)), ub$ndim,
                                model, FALSE, ub$t_cycle, coef_override = attr(gm, "coef"), R = trk_blur()$R_msd),
                 error = function(e) NULL)
      if (is.null(fg)) NULL else {
        Dtrue <- attr(gm, "D")
        rat   <- if (is.finite(fg$D) && Dtrue > 0) fg$D / Dtrue else NA_real_
        # Model-specific line: what does pure geometry hand back?
        say <- switch(model,
          confined = if (is.finite(fg$L_conf))
              sprintf("would itself fit as a %s of L_c &asymp; %.0f nm%s",
                      if (isTRUE(fg$degenerate)) "non-converged corral" else "corral",
                      fg$L_conf,
                      if (is.finite(fr$L_conf))
                        sprintf(", against your fitted %.0f nm", fr$L_conf) else "")
            else NULL,
          directed = if (is.finite(fg$v_um_s))
              sprintf("would itself fit as v &asymp; %.3g &micro;m/s%s", fg$v_um_s,
                      if (is.finite(fr$v_um_s)) sprintf(", against your fitted %.3g", fr$v_um_s) else "")
            else NULL,
          anomalous = if (is.finite(fg$alpha))
              sprintf("would itself fit as &alpha; &asymp; %.2f%s", fg$alpha,
                      if (is.finite(fr$alpha)) sprintf(", against your fitted %.2f", fr$alpha) else "")
            else NULL,
          if (is.finite(rat)) sprintf("fitted the same way returns about %.0f%% of its true D",
                                      5 * round(100 * rat / 5)) else NULL)
        if (is.null(say)) NULL else {
          # Warn when the shape alone already explains most of what was measured.
          bad <- switch(model,
            confined  = is.finite(fg$L_conf) && is.finite(fr$L_conf) && fg$L_conf < 2 * fr$L_conf,
            anomalous = is.finite(fg$alpha) && fg$alpha < 0.95,
            directed  = FALSE,
            is.finite(rat) && rat < 0.95)
          col <- if (isTRUE(bad)) "#b26a00" else "#3a7d3a"
          htmltools::div(style = sprintf("color:%s; font-size:0.86em; margin-top:2px;", col),
            htmltools::HTML(sprintf(
              "Simulated free diffusion on this %s (%s, %s) %s, %s",
              attr(gm, "geom"), attr(gm, "dims"), attr(gm, "mode"), say,
              if (isTRUE(bad))
                "geometry, not biology: shorten the upper fit limit."
              else "so the shape does not explain your value over this window.")))
        }
      }
    }
  }
  htmltools::div(style = "font-size:0.87em; color:#444; margin:6px 0;",
    dof_note, geom_note,
    htmltools::tags$b("MSD: "),
    sprintf("%s = %.4g%s  |  tDIFF (per-track mean, n=%d) = %.4g µm²/s%s  |  ",
            D_lbl, fr$D / 1e3, eD_ci, n_tD, tDIFF, tD_ci),
    sigloc(),
    # The offset is only separable in the LINEAR model. In the non-linear models it
    # competes with the plateau / exponent for the same part of the curve and is
    # bounded below at zero, so it routinely collapses to exactly 0 - which is a
    # failure to measure it, not a measurement of zero. Say so rather than print 0.00.
    if (identical(model, "linear") || (is.finite(fr$sigma_loc) && fr$sigma_loc > 0.01))
      sprintf(" = %.2f nm", fr$sigma_loc)
    else " not separable in this model",
    sprintf("  |  %s  |  ergodicity = %.2f%s", method_lbl, ergo, extra),
    htmltools::br(),
    htmltools::tags$span(style = "color:#888; font-size:0.9em;",
      sprintf("n = %d tracks%s  |  %s  |  cycle = %.4g ms  |  τ axis = %s",
              ub$n_tracks,
              if (!is.null(ub$n_dropped_negD) && ub$n_dropped_negD > 0)
                sprintf(" (−%d D≤0)", ub$n_dropped_negD) else "",
              if (isTRUE(ub$track_weighted)) "track-weighted (eMSD)" else "pair-weighted",
              ub$t_cycle,
              if (linear_tau) "linearized lag×cycle" else "real mean Δt (cycle-time-independent)")),
    # Confinement honesty. Bela Vogler's MINFLUX MSD workflow (and Vogler 2025) is built for
    # FREE Brownian diffusion; anomalous/confined MSD is explicitly "under development / implement
    # on your own" in his toolbox. So when fitting Linear, remind the user that in a confined
    # space (e.g. E. coli) the reported D is only an apparent/local D, and point to the CI +
    # Confined/Anomalous models - flagged as beyond-standard-for-MINFLUX.
    # One line, not a lecture: the reasoning is in the panel's Method box.
    if (identical(model, "linear"))
      htmltools::div(style = "color:#8a6d00; font-size:0.82em; margin-top:4px;",
        htmltools::HTML("&#9432; "),
        htmltools::tags$b("Free-Brownian fit"),
        " - under confinement this D is apparent/local. Cross-check the confinement index."),
    # Per-cluster eDIFF / tDIFF breakdown (when clustering is enabled).
    local({
      cdiff <- tryCatch(trk_unified_cluster_diff(), error = function(e) NULL)
      if (is.null(cdiff) || !nrow(cdiff)) return(NULL)
      rows <- lapply(seq_len(nrow(cdiff)), function(i) {
        r <- cdiff[i, ]
        htmltools::tags$div(style = sprintf("color:%s; font-size:0.85em; margin-left:8px;", r$color),
          sprintf("● Cluster %d: eDIFF = %s, tDIFF = %s µm²/s  (n = %d)", r$cluster,
                  if (is.finite(r$eDIFF)) sprintf("%.3g", r$eDIFF) else "-",
                  if (is.finite(r$tDIFF)) sprintf("%.3g", r$tDIFF) else "-", r$n))
      })
      htmltools::tagList(htmltools::br(),
        htmltools::tags$b(style = "font-size:0.85em;", "Per-cluster diffusion:"), rows)
    }))
})
