## ─────────────────────────────────────────────────────────────────────────────
## server/tracking/20_speed_transport.R  -  panel: speed and transport
## ─────────────────────────────────────────────────────────────────────────────
## Part of tracking_server(). Not a standalone script: its statements are evaluated in
## the server function's environment, in file-name order (server/parts.R).
##
## Per-track drift velocity and local alpha (R/est_transport.R), the
## population export summary, plots and downloads.
##
## INTERFACE (generated: Rscript tools/server_map.R --headers; do not edit by hand)
##   DEFINES for other parts:
##     (nothing)
##   READS from other parts:
##     01_state_sigma: speed_ran, trk_sigma_loc_eff
##     11_data_frames: trk
## END INTERFACE
## ─────────────────────────────────────────────────────────────────────────────

# Effective σ_loc for the speed/transport noise floor: the user's field if set
# (> 0), else the app-wide value (unified MSD / manual), else 0. Keeps the drift
# noise floor √2σ/W and the α de-bias consistent with the rest of the app.
trk_speed_sigma_eff <- reactive({
  v <- suppressWarnings(as.numeric(input$trk_speed_sigma_nm %||% 0))
  if (is.finite(v) && v > 0) return(v)
  s <- tryCatch(trk_sigma_loc_eff(), error = function(e) NA_real_)
  if (is.finite(s) && s > 0) s else 0
})

# ==========================================================================
# trk_export_summary() - DIRECTED TRANSPORT ("export") DETECTION, per track
# ==========================================================================
#
# GOAL
#   Detect when a molecule stops diffusing and starts being CARRIED - a
#   sustained, directional excursion (e.g. motor-driven transport, or a
#   molecule leaving a compartment). Produces one summary row per track: did
#   it export, how fast, and WHERE along the track it happened.
#
# THE CENTRAL DIFFICULTY: DIFFUSION LOOKS FAST TOO
#   You cannot detect transport from speed. From <r²> = 2dDΔt, the per-step
#   speed of a purely diffusing molecule is
#
#        v_step = r/Δt ~ sqrt(2dD/Δt)    ->  DIVERGES as Δt -> 0
#
#   so a freely diffusing molecule shows an arbitrarily high instantaneous
#   "speed" if you sample it finely enough (header §2b: Brownian motion has no
#   instantaneous velocity). Any threshold on path speed therefore measures
#   your sampling rate, not transport.
#
#   THE FIX - NET DISPLACEMENT OVER A WINDOW. Diffusion is unbiased: over a
#   window W its net displacement grows only as sqrt(W), because the steps
#   cancel. Directed motion is biased: its net displacement grows LINEARLY
#   with W. So the drift velocity
#
#        v_drift(t) = |r(t+W) − r(t)| / W
#
#   stays small for diffusion but stays high for transport, and the two
#   separate more strongly the longer the window. That is what is thresholded.
#
# WHAT COUNTS AS AN EXPORT (three guards against false positives)
#   1. THRESHOLD = max(user or auto threshold, localization-noise floor).
#      Auto = median + 3·MAD (robust, so a few fast tracks don't set the bar).
#      Noise floor = sqrt(2)·σ_loc / W - the apparent drift velocity produced
#      by localization error alone over that window.
#   2. SUSTAINED RUN of >= min_run consecutive windows above threshold. A
#      single spike is noise; transport persists.
#   3. Tracks too short for a full window are skipped rather than guessed at.
#
#   onset_frac = where the export starts along the track (0 = start, 1 = end).
#   This is the diagnostic that reveals the "exports cluster at track ends"
#   pattern - which is usually an ARTEFACT (photobleaching / loss of lock at
#   the end of a track), not biology. Check this histogram before believing
#   an export result.
#
# OBJECTS IT CALLS (inputs)
#   trk(), trk_speed_sigma_eff() (σ_loc for the noise floor),
#   trk_drift_velocity_one() (the per-track drift-velocity curve)
#   inputs: trk_speed_drift_win (W), trk_speed_drift_minpts,
#           trk_speed_min_pts, threshold + min-run settings
#
# RETURNS a data.frame, one row per qualifying track:
#   tid, has_export, v_max / v_export (µm/s), export_pos_frac (onset),
#   run length, and the threshold actually used.
#
# CAVEAT: this reports SUSTAINED DIRECTIONAL BIAS. Confirm against the MSD
# directed model (which fits a velocity) and the HMM/per-state MSD before
# calling it active transport - confinement release and flow can look similar.
# ==========================================================================
trk_export_summary <- eventReactive(input$trk_speed_run, {
  req(isTRUE(input$trk_show_speed))
  df <- trk(); req(nrow(df) > 0)
  min_pts   <- max(5L, as.integer(input$trk_speed_min_pts %||% 20L))
  smooth    <- as.integer(input$trk_speed_smooth %||% 5L)
  drift_win <- as.numeric(input$trk_speed_drift_win %||% 200)
  drift_mp  <- as.integer(input$trk_speed_drift_minpts %||% 5L)
  sigma_nm  <- trk_speed_sigma_eff()
  onset_thr <- as.numeric(input$trk_speed_onset_thr %||% 0)
  fold      <- as.numeric(input$trk_speed_fold %||% 3)
  min_run   <- max(1L, as.integer(input$trk_speed_min_run %||% 3L))
  guard_win <- max(0, as.numeric(input$trk_speed_edge_guard %||% 1))
  a_win     <- as.integer(input$trk_speed_alpha_win %||% 21L)
  a_maxlag  <- as.integer(input$trk_speed_alpha_maxlag %||% 6L)
  vfloor    <- if (is.finite(sigma_nm) && sigma_nm > 0 && drift_win > 0) sqrt(2) * sigma_nm / drift_win else 0

  counts <- table(df$tid); tids <- names(counts[counts >= min_pts])
  shiny::validate(shiny::need(length(tids) > 0, "No tracks with enough points - lower 'Min points per TID'."))

  rows <- withProgress(message = "Speed & transport…", value = 0,
    lapply(seq_along(tids), function(ti_sp) {
    if (ti_sp %% 25L == 0L) incProgress(25 / length(tids), detail = "per-track drift/α")
    t_i <- tids[ti_sp]
    d  <- df[df$tid == t_i, , drop = FALSE]
    dv <- trk_drift_velocity_one(d, drift_win, drift_mp, smooth)
    if (is.null(dv) || nrow(dv) < max(2L, min_run)) return(NULL)
    al <- trk_local_alpha_one(d, a_win, a_maxlag, sigma_nm)
    v  <- dv$v_drift_med
    thr_eff <- trk_export_thr(v, onset_thr, fold, vfloor)
    # Edge guard: the first/last window of a track is unreliable (initial
    # localization jitter at track start, collapsing span at the end), which
    # otherwise piles false "export at frac 0/1" events at the edges. Restrict
    # detection to interior windows >= guard_win·W from each end.
    tmin <- min(dv$t_ms); tmax <- max(dv$t_ms); t_span <- tmax - tmin
    guard    <- guard_win * drift_win
    interior <- dv$t_ms >= (tmin + guard) & dv$t_ms <= (tmax - guard)
    n_int    <- sum(interior)
    above <- (v >= thr_eff) & interior
    r <- rle(above); ends <- cumsum(r$lengths); starts <- ends - r$lengths + 1L
    seg_ok <- r$values & r$lengths >= min_run
    has_exp <- isTRUE(any(seg_ok))          # NA-proof
    n_events <- sum(seg_ok, na.rm = TRUE)   # number of directed episodes (hops/exports)
    # "Sits still" fraction: interior windows whose net motion is at/below the
    # localization-noise floor - i.e. confined/immobile at this timescale. High
    # immobile_frac + several n_events = the hop-diffusion (dwell → hop → dwell) look.
    immobile_frac <- if (n_int > 0) sum(interior & v <= max(vfloor, 0), na.rm = TRUE) / n_int else NA_real_
    onset_ms <- NA_real_; onset_frac <- NA_real_; peak_frac <- NA_real_
    export_v <- NA_real_; export_dur <- NA_real_
    # Per-event corroboration (artifact vs signal). A GENUINE directed excursion is
    # made of INTEGRATED localizations with normal photon counts - the molecule is
    # present and tracked throughout. A tracker SEARCH after losing the molecule has
    # the same gross signature (fast, outward, track ends) but its localizations are
    # dark-cycle-dominated and photon-starved. These fields come from dark cycles and
    # photons, NOT from drift velocity, so they do not restate what defined the event.
    ev_integ <- NA_real_; ev_eco_ratio <- NA_real_; ev_cfr <- NA_real_
    ev_nloc  <- NA_integer_; ev_at_end <- NA; ev_support <- NA_character_
    t0_abs <- suppressWarnings(min(d$tim, na.rm = TRUE))
    if (has_exp) {
      first    <- which(seg_ok)[1]
      onset_ms <- dv$t_ms[starts[first]]
      onset_frac <- if (is.finite(t_span) && t_span > 0) (onset_ms - tmin) / t_span else NA_real_
      keep <- unlist(mapply(function(s, e) s:e, starts[seg_ok], ends[seg_ok], SIMPLIFY = FALSE))
      # Export POSITION = the strongest directed window (peak drift velocity),
      # a far less start-biased locator than the first threshold crossing.
      peak_idx  <- keep[which.max(v[keep])]
      peak_frac <- if (is.finite(t_span) && t_span > 0) (dv$t_ms[peak_idx] - tmin) / t_span else NA_real_
      export_v   <- stats::median(v[keep])
      export_dur <- if (n_int > 0) sum(r$lengths[seg_ok]) / n_int * t_span else NA_real_

      # Corroborate the STRONGEST event (the one export_v / export_pos_frac report):
      # find the run containing the peak window, widen by W/2 to cover the
      # localizations the windows were actually built from, then score that span.
      wrun <- which(seg_ok)
      prun <- wrun[which(peak_idx >= starts[wrun] & peak_idx <= ends[wrun])[1]]
      if (!is.na(prun)) {
        ev_lo <- dv$t_ms[starts[prun]] - drift_win / 2
        ev_hi <- dv$t_ms[ends[prun]]   + drift_win / 2
        rel   <- d$tim - t0_abs
        insid <- is.finite(rel) & rel >= ev_lo & rel <= ev_hi
        ev_nloc <- sum(insid)
        if (ev_nloc > 0) {
          if ("n_dark" %in% names(d))
            ev_integ <- mean(d$n_dark[insid] == 0L, na.rm = TRUE)
          if ("eco" %in% names(d)) {
            e_ev <- stats::median(d$eco[insid], na.rm = TRUE)
            e_tr <- stats::median(d$eco,        na.rm = TRUE)   # track's own baseline
            if (is.finite(e_ev) && is.finite(e_tr) && e_tr > 0) ev_eco_ratio <- e_ev / e_tr
          }
          if ("cfr" %in% names(d)) ev_cfr <- stats::median(d$cfr[insid], na.rm = TRUE)
        }
        ev_at_end <- isTRUE((tmax - dv$t_ms[ends[prun]]) <= drift_win)
        ev_support <-
          if (is.finite(ev_integ) && ev_integ < 0.5) "dark-dominated"
          else if (is.finite(ev_eco_ratio) && ev_eco_ratio < 0.6) "low-photon"
          else if (is.finite(ev_integ) && ev_integ >= 0.5) "supported"
          else "unscored"
      }
    }
    # Track-level TAIL quality: integrated fraction over the final W of the track.
    # Directly answers "is the outward tail the molecule, or the tracker hunting?"
    tail_integ <- NA_real_
    if ("n_dark" %in% names(d) && is.finite(t0_abs)) {
      rel <- d$tim - t0_abs
      rmx <- suppressWarnings(max(rel, na.rm = TRUE))
      tl  <- is.finite(rel) & rel >= (rmx - drift_win)
      if (any(tl)) tail_integ <- mean(d$n_dark[tl] == 0L, na.rm = TRUE)
    }
    data.frame(
      tid = t_i, n_loc = nrow(d), track_dur_ms = t_span,
      v_drift_peak = max(v), v_drift_median = stats::median(v),
      directed_frac = if (n_int > 0) sum(above, na.rm = TRUE) / n_int else 0,
      immobile_frac = immobile_frac,
      alpha_max = if (!is.null(al)) max(al$alpha) else NA_real_,
      has_export = has_exp, n_events = n_events,
      export_onset_ms = onset_ms, export_onset_frac = onset_frac,
      export_pos_frac = peak_frac,
      export_v = export_v, export_dur_ms = export_dur,
      # --- corroboration for the strongest event (independent of drift velocity) ---
      event_support     = ev_support,      # supported | dark-dominated | low-photon
      event_integ_frac  = ev_integ,        # 1 = every loc in the event was integrated
      event_eco_ratio   = ev_eco_ratio,    # event photons / this track's own median
      event_cfr_med     = ev_cfr,
      event_n_loc       = ev_nloc,
      event_at_track_end = ev_at_end,
      tail_integ_frac   = tail_integ,      # integrated fraction in the track's final W
      stringsAsFactors = FALSE
    )
  }))
  out <- do.call(rbind, rows)
  shiny::validate(shiny::need(!is.null(out) && nrow(out) > 0, "No tracks produced a drift-velocity curve."))
  out
})

output$trk_speed_tid_ui <- renderUI({
  req(isTRUE(input$trk_show_speed), speed_ran())
  s <- tryCatch(trk_export_summary(), error = function(e) NULL)
  req(!is.null(s), nrow(s) > 0)
  # Rank by segment strength so the picker doubles as a follow-up list:
  # flagged tracks first, then by peak drift velocity (descending).
  s <- s[order(-as.integer(s$has_export %in% TRUE), -s$v_drift_peak), , drop = FALSE]
  labs <- sprintf("%s%s  (peak %.2g µm/s)",
                  ifelse(s$has_export %in% TRUE, "★ ", "  "),
                  s$tid, s$v_drift_peak)
  choices <- stats::setNames(s$tid, labs)
  sel <- if (!is.null(input$trk_speed_tid) && input$trk_speed_tid %in% s$tid)
           input$trk_speed_tid else s$tid[1]
  selectInput("trk_speed_tid", "Track (TID) - ★ = candidate directed segment, ranked by peak drift",
              choices = choices, selected = sel)
})

output$trk_speed_plot <- renderPlotly({
  req(isTRUE(input$trk_show_speed), speed_ran())
  df <- trk(); req(nrow(df) > 0)
  dark <- isTRUE(isolate(input$trk_dark_mode))
  bg <- if (dark) "#000000" else "#FFFFFF"; fg <- if (dark) "#FFFFFF" else "#333333"
  min_pts <- max(5L, as.integer(input$trk_speed_min_pts %||% 20L))
  smooth  <- as.integer(input$trk_speed_smooth %||% 5L)
  counts  <- table(df$tid); tids <- names(counts[counts >= min_pts])
  shiny::validate(shiny::need(length(tids) > 0, "No tracks with enough points - lower 'Min points per TID'."))
  tid_sel <- input$trk_speed_tid
  if (is.null(tid_sel) || !(tid_sel %in% tids)) tid_sel <- tids[1]

  # New detectors (default on): rolling drift velocity + local exponent α(t).
  show_path  <- isTRUE(input$trk_speed_show_path  %||% TRUE)
  show_drift <- isTRUE(input$trk_speed_show_drift %||% TRUE)
  show_alpha <- isTRUE(input$trk_speed_show_alpha %||% TRUE)
  drift_win  <- as.numeric(input$trk_speed_drift_win %||% 200)
  drift_mp   <- as.integer(input$trk_speed_drift_minpts %||% 5L)
  sigma_nm   <- trk_speed_sigma_eff()
  a_win      <- as.integer(input$trk_speed_alpha_win %||% 21L)
  a_maxlag   <- as.integer(input$trk_speed_alpha_maxlag %||% 6L)
  onset_thr  <- as.numeric(input$trk_speed_onset_thr %||% 0)
  fold       <- as.numeric(input$trk_speed_fold %||% 3)

  fig <- plot_ly()

  dsel <- df[df$tid == tid_sel, , drop = FALSE]
  spd  <- trk_speed_one(dsel, smooth)
  dv   <- trk_drift_velocity_one(dsel, drift_win, drift_mp, smooth)
  al   <- trk_local_alpha_one(dsel, a_win, a_maxlag, sigma_nm)
  shiny::validate(shiny::need(!is.null(spd) && nrow(spd) >= 2, "Selected track has too few valid steps."))

  # Path speed (per-step) - kept for reference; it is dominated by diffusion+noise
  if (show_path) {
    fig <- fig %>%
      add_trace(x = spd$t_ms, y = spd$speed_um_s, type = "scattergl", mode = "markers",
                marker = list(size = 4, color = "rgba(140,140,140,0.30)"),
                name = "path speed (raw)", hoverinfo = "x+y") %>%
      add_trace(x = spd$t_ms, y = spd$speed_med, type = "scattergl", mode = "lines",
                line = list(color = "rgba(120,120,120,0.85)", width = 1.5, dash = "dot"),
                name = "path speed (median)")
  }

  # Rolling drift (transport) velocity - the segment-detecting curve
  if (show_drift && !is.null(dv) && nrow(dv) >= 2) {
    fig <- fig %>%
      add_trace(x = dv$t_ms, y = dv$v_drift_med, type = "scattergl", mode = "lines",
                line = list(color = "#1B9E5A", width = 3), name = "drift velocity")
    vfloor <- if (is.finite(sigma_nm) && sigma_nm > 0 && is.finite(drift_win) && drift_win > 0)
                sqrt(2) * sigma_nm / drift_win else 0
    # Localization-noise floor √2·σ/W (µm/s): drift bumps below this are noise.
    if (vfloor > 0) {
      xr <- range(dv$t_ms)
      fig <- fig %>% add_trace(x = xr, y = c(vfloor, vfloor), type = "scatter", mode = "lines",
        line = list(color = "rgba(200,120,0,0.9)", width = 1.2, dash = "dash"),
        name = "noise floor √2σ/W", hoverinfo = "name")
    }
    # Onset overlay: same threshold as the population classifier (fold × baseline,
    # or absolute if set), highlight windows above it.
    thr <- trk_export_thr(dv$v_drift_med, onset_thr, fold, vfloor)
    if (is.finite(thr) && thr > 0) {
      xr <- range(dv$t_ms)
      above <- dv$v_drift_med >= thr
      fig <- fig %>% add_trace(x = xr, y = c(thr, thr), type = "scatter", mode = "lines",
        line = list(color = "rgba(224,48,32,0.85)", width = 1.2, dash = "dash"),
        name = "onset threshold", hoverinfo = "name")
      if (any(above))
        fig <- fig %>% add_trace(x = dv$t_ms[above], y = dv$v_drift_med[above],
          type = "scattergl", mode = "markers",
          marker = list(size = 6, color = "#E03020"),
          name = "above threshold", hoverinfo = "x+y")
    }
  }

  # Local exponent α(t) on a secondary axis (0..2): 1 = Brownian, 2 = directed
  if (show_alpha && !is.null(al) && nrow(al) >= 2) {
    fig <- fig %>% add_trace(x = al$t_ms, y = al$alpha, type = "scattergl", mode = "lines",
      line = list(color = "#7B4FBE", width = 2), name = "local α (right axis)",
      yaxis = "y2", hoverinfo = "x+y")
  }

  # α reference lines (1 = Brownian, 2 = ballistic) on the secondary axis
  a_shapes <- list()
  if (show_alpha && !is.null(al) && nrow(al) >= 2) {
    for (av in c(1, 2))
      a_shapes[[length(a_shapes) + 1L]] <- list(type = "line", xref = "paper",
        x0 = 0, x1 = 1, yref = "y2", y0 = av, y1 = av,
        line = list(color = "rgba(123,79,190,0.35)", width = 1, dash = "dot"))
  }

  fig %>% layout(
    xaxis = list(title = "time within track (ms)", color = fg),
    yaxis = list(title = "velocity (µm/s)", color = fg, rangemode = "tozero"),
    yaxis2 = list(title = "local α", overlaying = "y", side = "right",
                  range = c(0, 2.1), color = "#7B4FBE", showgrid = FALSE),
    shapes = a_shapes,
    paper_bgcolor = bg, plot_bgcolor = bg, font = list(color = fg),
    legend = list(font = list(color = fg, size = 10), orientation = "h", y = -0.2),
    title = list(text = sprintf("TID %s", tid_sel), font = list(color = fg, size = 12)),
    uirevision = "speed")
})

# ---- Population view: quantify export across all tracks -----------------------
output$trk_export_summary_txt <- renderUI({
  req(isTRUE(input$trk_show_speed))
  s <- trk_export_summary()
  n  <- nrow(s); ne <- sum(s$has_export, na.rm = TRUE)
  med_peak <- stats::median(s$v_drift_peak, na.rm = TRUE)
  exp_rows <- s[s$has_export %in% TRUE, , drop = FALSE]
  med_onset <- if (nrow(exp_rows)) stats::median(exp_rows$export_pos_frac, na.rm = TRUE) else NA_real_
  med_expv  <- if (nrow(exp_rows)) stats::median(exp_rows$export_v, na.rm = TRUE) else NA_real_
  n_multi   <- if (nrow(exp_rows)) sum(exp_rows$n_events >= 2, na.rm = TRUE) else 0
  # Corroboration breakdown: how many candidate events are actually supported by
  # evidence independent of the drift threshold that found them.
  sup  <- if (nrow(exp_rows)) sum(exp_rows$event_support %in% "supported")      else 0L
  dkd  <- if (nrow(exp_rows)) sum(exp_rows$event_support %in% "dark-dominated") else 0L
  lowp <- if (nrow(exp_rows)) sum(exp_rows$event_support %in% "low-photon")     else 0L
  med_ti <- stats::median(s$tail_integ_frac, na.rm = TRUE)
  htmltools::div(style = "font-size:0.9em; line-height:1.5;",
    htmltools::HTML(sprintf(
      "<b>%d</b> tracks analysed &middot; <b>%d (%.0f%%)</b> with a candidate directed segment%s<br>%s%smedian peak drift = <b>%.3g µm/s</b>",
      n, ne, if (n > 0) 100 * ne / n else 0,
      if (n_multi > 0) sprintf(" &middot; <b>%d</b> with multiple episodes (hop/multi-step)", n_multi) else "",
      if (is.finite(med_expv))  sprintf("median segment velocity = <b>%.3g µm/s</b> &middot; ", med_expv) else "",
      if (is.finite(med_onset)) sprintf("median segment at <b>%.0f%%</b> of track length &middot; ", 100 * med_onset) else "",
      med_peak)),
    # Corroboration = evidence from dark cycles and photons, i.e. NOT from the drift
    # threshold that found the segment. The circularity argument itself is in the
    # panel's Method box, not reprinted on every run.
    if (ne > 0) htmltools::HTML(sprintf(
      "<div style='margin-top:6px;'>Corroboration (dark cycles &amp; photons, independent of drift): <b style='color:#1B9E5A;'>%d supported</b> &middot; <b style='color:#C0392B;'>%d dark-dominated</b> &middot; <b style='color:#B8860B;'>%d low-photon</b></div>",
      sup, dkd, lowp)),
    if (is.finite(med_ti)) htmltools::HTML(sprintf(
      "<div style='margin-top:4px; color:#555;'>Median tail integrated fraction = <b>%.2f</b>%s</div>",
      med_ti, if (med_ti < 0.5) " - track ends are search-dominated" else "")),
    htmltools::HTML(
      "<div style='margin-top:6px; color:#B8860B; font-size:0.95em;'>&#9888; Segment velocity is <b>not independent evidence</b>: the flag that found it was a speed. Judge on corroboration and a matched negative control.</div>")
  )
})

# Distribution of a per-track metric (the "histogram of speeds" across tracks)
output$trk_export_hist <- renderPlotly({
  req(isTRUE(input$trk_show_speed))
  s <- trk_export_summary()
  dark <- isTRUE(isolate(input$trk_dark_mode))
  bg <- if (dark) "#000000" else "#FFFFFF"; fg <- if (dark) "#FFFFFF" else "#333333"
  metric <- input$trk_export_hist_metric %||% "v_drift_peak"
  lab <- switch(metric,
    v_drift_peak     = "peak drift velocity (µm/s)",
    export_v         = "segment velocity (µm/s)",
    directed_frac    = "fraction of track time directed",
    immobile_frac    = "fraction of track time immobile/confined",
    n_events         = "number of transport episodes per track",
    alpha_max        = "max local α",
    event_integ_frac = "integrated fraction within the event (1 = all integrated)",
    event_eco_ratio  = "event photons / track median photons",
    tail_integ_frac  = "integrated fraction in the track's final window",
    "value")
  # Event-scoped metrics exist only for flagged tracks; the rest span all tracks.
  event_only <- c("export_v", "event_integ_frac", "event_eco_ratio")
  vals_all <- s[[metric]]
  exp_mask <- s$has_export %in% TRUE
  fig <- plot_ly()
  if (!(metric %in% event_only)) {
    # unflagged vs flagged stacked so the split is visible
    ve <- vals_all[exp_mask];  vn <- vals_all[!exp_mask]
    ve <- ve[is.finite(ve)];   vn <- vn[is.finite(vn)]
    if (length(vn)) fig <- fig %>% add_histogram(x = vn, name = "unflagged",
      marker = list(color = "rgba(150,150,150,0.55)"), nbinsx = 30)
    if (length(ve)) fig <- fig %>% add_histogram(x = ve, name = "flagged",
      marker = list(color = "rgba(27,158,90,0.75)"), nbinsx = 30)
    fig <- fig %>% layout(barmode = "stack")
  } else {
    v <- vals_all[exp_mask]; v <- v[is.finite(v)]
    shiny::validate(shiny::need(length(v) > 0, "No flagged tracks yet - loosen the threshold / min-run."))
    fig <- fig %>% add_histogram(x = v, name = "flagged",
      marker = list(color = "rgba(27,158,90,0.75)"), nbinsx = 30)
  }
  fig %>% layout(
    xaxis = list(title = lab, color = fg),
    yaxis = list(title = "tracks", color = fg),
    paper_bgcolor = bg, plot_bgcolor = bg, font = list(color = fg),
    legend = list(font = list(color = fg, size = 10), orientation = "h", y = -0.25),
    title = list(text = "Per-track distribution", font = list(color = fg, size = 12)),
    uirevision = "exporthist")
})

# Where along the track is transport strongest? (0 = start, 1 = end)
# Uses the peak-drift-velocity position (not first crossing, which is start-biased).
# ── Is a late-track excursion real, or the tracker searching? ────────────────
# WHERE an event sits along the track cannot answer that on its own: a molecule
# leaving and a tracker losing its molecule both produce fast outward motion at the
# end. What separates them is whether the localizations during the event were still
# being found on schedule. So plot position (x) against the event's integrated
# fraction (y): a molecule that is genuinely present keeps being localized and sits
# high; a tracker searching accumulates dark cycles and sinks, typically at the far
# right where tracks end. y is independent of the drift threshold that found the
# event, which is what makes it evidence rather than restatement.
output$trk_export_corrob <- renderPlotly({
  req(isTRUE(input$trk_show_speed))
  s <- trk_export_summary()
  dark <- isTRUE(isolate(input$trk_dark_mode))
  bg <- if (dark) "#000000" else "#FFFFFF"; fg <- if (dark) "#FFFFFF" else "#333333"
  e <- s[s$has_export %in% TRUE, , drop = FALSE]
  e <- e[is.finite(e$export_onset_frac) & is.finite(e$event_integ_frac), , drop = FALSE]
  shiny::validate(shiny::need(nrow(e) > 0,
    "No candidate segments yet - press 'Run speed & transport', or relax the threshold."))
  pal <- c(supported = "#1B9E5A", `dark-dominated` = "#C0392B", `low-photon` = "#B8860B")
  fig <- plot_ly()
  for (lv in names(pal)) {
    i <- which(e$event_support %in% lv)
    if (!length(i)) next
    fig <- fig %>% add_trace(
      x = e$export_onset_frac[i], y = e$event_integ_frac[i],
      type = "scatter", mode = "markers",
      marker = list(size = 8, color = pal[[lv]], opacity = 0.75,
                    line = list(width = 0.4, color = "rgba(0,0,0,0.35)")),
      name = sprintf("%s (%d)", lv, length(i)),
      # hoverinfo="text" rather than a hovertemplate: no template parsing, so the
      # tooltip cannot fall back to printing the placeholder itself.
      text = sprintf("TID %s<br>starts at %.0f%% of track<br>integrated fraction %.2f<br>photons %.2f x track median<br>%s",
                     as.character(e$tid[i]), 100 * e$export_onset_frac[i],
                     e$event_integ_frac[i], e$event_eco_ratio[i],
                     as.character(e$event_support[i])),
      hoverinfo = "text")
  }
  fig %>% layout(
    xaxis = list(title = "where the segment starts along the track (0 = start, 1 = end)",
                 range = c(-0.02, 1.02), color = fg),
    yaxis = list(title = "integrated fraction during the segment",
                 range = c(-0.02, 1.02), color = fg),
    shapes = list(list(type = "line", x0 = -0.02, x1 = 1.02, y0 = 0.5, y1 = 0.5,
                       line = list(color = "rgba(120,120,120,0.45)", dash = "dot", width = 1))),
    annotations = list(list(x = 0.02, y = 0.52, xref = "x", yref = "y",
      text = "above: molecule still being found · below: dark-cycle dominated",
      showarrow = FALSE, xanchor = "left", yanchor = "bottom",
      font = list(size = 9, color = "rgba(120,120,120,0.9)"))),
    paper_bgcolor = bg, plot_bgcolor = bg, font = list(color = fg),
    legend = list(font = list(color = fg, size = 10), orientation = "h",
                  x = 0, xanchor = "left", y = -0.24, yanchor = "top"),
    margin = list(l = 60, r = 20, t = 10, b = 90),
    uirevision = "export_corrob")
})

output$trk_export_onset_hist <- renderPlotly({
  req(isTRUE(input$trk_show_speed))
  s <- trk_export_summary()
  dark <- isTRUE(isolate(input$trk_dark_mode))
  bg <- if (dark) "#000000" else "#FFFFFF"; fg <- if (dark) "#FFFFFF" else "#333333"
  v <- s$export_pos_frac[s$has_export %in% TRUE]; v <- v[is.finite(v)]
  shiny::validate(shiny::need(length(v) > 0, "No flagged tracks yet - loosen the threshold / min-run."))
  plot_ly(x = v, type = "histogram", nbinsx = 20,
          marker = list(color = "rgba(123,79,190,0.75)")) %>%
    layout(
      xaxis = list(title = "segment position (0 = start … 1 = end of track)",
                   range = c(0, 1), color = fg),
      yaxis = list(title = "flagged tracks", color = fg),
      paper_bgcolor = bg, plot_bgcolor = bg, font = list(color = fg),
      # Neutral wording: a late pile-up is equally "molecule left" and "tracker lost it".
      title = list(text = sprintf("Segment location (n=%d) - pile-up near 1 is ambiguous: check corroboration", length(v)),
                   font = list(color = fg, size = 12)),
      uirevision = "onsethist")
})

# Per-track summary table download (the quantification deliverable)
output$trk_dl_export_csv <- downloadHandler(
  filename = function() paste0("export_summary_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
  content  = function(file) write.csv(trk_export_summary(), file, row.names = FALSE)
)

output$trk_dl_speed_csv <- downloadHandler(
  filename = function() paste0("speed_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
  content  = function(file) {
    df <- tryCatch(trk(), error = function(e) NULL); req(!is.null(df))
    min_pts   <- max(5L, as.integer(input$trk_speed_min_pts %||% 20L))
    smooth    <- as.integer(input$trk_speed_smooth %||% 5L)
    drift_win <- as.numeric(input$trk_speed_drift_win %||% 200)
    drift_mp  <- as.integer(input$trk_speed_drift_minpts %||% 5L)
    a_win     <- as.integer(input$trk_speed_alpha_win %||% 21L)
    a_maxlag  <- as.integer(input$trk_speed_alpha_maxlag %||% 6L)
    sigma_nm  <- trk_speed_sigma_eff()
    counts  <- table(df$tid); tids <- names(counts[counts >= min_pts])
    # Long/tidy format: all three metrics stacked (they have different time bases).
    out <- do.call(rbind, lapply(tids, function(t_i) {
      d <- df[df$tid == t_i, , drop = FALSE]
      sp <- trk_speed_one(d, smooth)
      dv <- trk_drift_velocity_one(d, drift_win, drift_mp, smooth)
      al <- trk_local_alpha_one(d, a_win, a_maxlag, sigma_nm)
      rows <- list()
      if (!is.null(sp)) rows[[length(rows)+1L]] <- data.frame(tid = t_i, metric = "path_speed_um_s",
        t_ms = sp$t_ms, value = sp$speed_um_s, value_smoothed = sp$speed_med)
      if (!is.null(dv)) rows[[length(rows)+1L]] <- data.frame(tid = t_i, metric = "drift_velocity_um_s",
        t_ms = dv$t_ms, value = dv$v_drift, value_smoothed = dv$v_drift_med)
      if (!is.null(al)) rows[[length(rows)+1L]] <- data.frame(tid = t_i, metric = "local_alpha",
        t_ms = al$t_ms, value = al$alpha, value_smoothed = NA_real_)
      if (!length(rows)) return(NULL)
      do.call(rbind, rows)
    }))
    write.csv(out, file, row.names = FALSE)
  }
)
