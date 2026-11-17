## ─────────────────────────────────────────────────────────────────────────────
## server/tracking/18_confinement_displacement.R  -  panels: confinement index; displacement versus diffusion
## ─────────────────────────────────────────────────────────────────────────────
## Part of tracking_server(). Not a standalone script: its statements are evaluated in
## the server function's environment, in file-name order (server/parts.R).
##
## The MSD confinement index with its free-diffusion null, and net
## displacement against the free and geometry-aware nulls.
##
## INTERFACE (generated: Rscript tools/server_map.R --headers; do not edit by hand)
##   DEFINES for other parts:
##     (nothing)
##   READS from other parts:
##     01_state_sigma: .sigma_stamp, trk_sigma_loc_eff, trk_sigma_loc_src
##     02_data_quality_cve: trk_cve_bundle
##     07_cycle_time_sequence: trk_blur, trk_cycle_time
##     11_data_frames: trk
##     15_msd: .geom_for_null, msd_unified_bundle, msd_unified_modelfit,
##             trk_unified_tdiff
## END INTERFACE
## ─────────────────────────────────────────────────────────────────────────────

# ============================================================================
# conf_ci_bundle() - PER-TRACK CONFINEMENT FROM THE MSD SHAPE (CI + simulated null)
# ============================================================================
# The most sensitive per-track confinement test in the app. Uses the confinement
# index already computed per track by the MSD:
#
#     CI = (MSD_late / MSD_early) / (tau_late / tau_early)
#
# WHY THIS ONE. It compares a track's own SHORT lags against its own LONG lags, so
# it needs no external D and cannot be biased by D varying between molecules - the
# confound that inflates the Displacement-vs-Diffusion split. It is not circular
# (short and long lags are different information) and it inherits the cycle-index
# binning, so irregular Δt is already handled.
#
# WHY IT NEEDS A SIMULATED NULL. For ideal free motion CI = 1, but the OBSERVED free
# value is not: the 2·d·σ_loc² intercept sits on every MSD point and lifts the EARLY
# window more than the late one, pushing free CI below 1. How far below depends on
# σ_loc, D, cycle time, fit window and track length, so there is no universal
# threshold - the null is simulated at the data's own conditions and matched per
# track on localization count. Simulated free CI medians run ≈0.67–0.77 with the 5th
# percentile moving 0.28 → 0.63 from N=50 to N=800, which is why matching on N matters.
#
# VALIDATION (in silico, MINFLUX-like: t_c 0.191 ms, σ_loc 10 nm, D 0.2 µm²/s, N=200):
#   free → 7.3% confined / 89.2% free-like / 3.5% directed  (target 5/90/5)
#   confined L=100 nm → 89.5% detected;  L=250 nm → 25.0% detected
#   directed v=0.4 µm/s → NOT detected (6.2%) - see the caveat below.
#
# LIMITS. Detection is set by how far a free molecule would travel in the same time
# relative to the confinement size: the resolution floor √(2·d·D·T) is reported so a
# null result can be read as "no confinement below this scale" rather than "free".
# CI is weak for DIRECTED motion (a modest drift barely lifts the late/early ratio) -
# use Directionality or Speed & Transport for transport.
conf_ci_bundle <- eventReactive(input$trk_conf_ci_run, {
  td <- tryCatch(trk_unified_tdiff(), error = function(e) NULL)
  shiny::validate(shiny::need(!is.null(td) && !is.null(td$pertrack) && nrow(td$pertrack) > 0,
    "Run the MSD first (Plot MSD) - the confinement index comes from the per-track MSD."))
  pt <- td$pertrack
  pt <- pt[is.finite(pt$conf_index), , drop = FALSE]
  shiny::validate(shiny::need(nrow(pt) >= 5L, "Too few tracks with a usable confinement index."))

  df  <- trk(); t_c <- tryCatch(trk_cycle_time(), error = function(e) NA_real_)
  shiny::validate(shiny::need(is.finite(t_c) && t_c > 0, "Cycle time unknown."))
  ub  <- tryCatch(msd_unified_bundle(), error = function(e) NULL)
  ndim <- if (!is.null(ub)) ub$ndim else 2L
  sig <- tryCatch(trk_sigma_loc_eff(), error = function(e) NA_real_)
  if (!is.finite(sig) || sig < 0) sig <- 0
  b   <- tryCatch(trk_cve_bundle(), error = function(e) NULL)
  D_use <- if (!is.null(b) && is.finite(b$D_cve) && b$D_cve > 0) b$D_cve else {
    dm <- stats::median(pt$D_um2_s[is.finite(pt$D_um2_s) & pt$D_um2_s > 0], na.rm = TRUE)
    if (is.finite(dm)) dm * 1e3 else NA_real_
  }
  shiny::validate(shiny::need(is.finite(D_use) && D_use > 0,
    "Need a D - run Data Quality & Calibration, or the MSD."))

  # Per-track localization count and the dataset's real dark-gap distribution, so the
  # null tracks are sampled exactly like the data.
  npts <- table(df$tid)
  pt$n_loc <- as.integer(npts[as.character(pt$tid)])
  pt <- pt[is.finite(pt$n_loc) & pt$n_loc >= 5L, , drop = FALSE]
  shiny::validate(shiny::need(nrow(pt) >= 5L, "Too few tracks after matching localization counts."))
  d2 <- df[order(df$tid, df$tim), ]
  same <- c(FALSE, d2$tid[-1] == d2$tid[-nrow(d2)])
  kk <- pmax(1L, as.integer(round(c(NA, diff(d2$tim))[same] / t_c)))
  kk <- kk[is.finite(kk) & kk >= 1L]
  gap_p <- if (length(kk) > 50L) { tb <- tabulate(pmin(kk, 15L), nbins = 15L); tb / sum(tb) }
           else c(0.55, 0.22, 0.09, 0.05, 0.04, 0.05)

  min_lag <- max(1L, as.integer(input$trk_uni_min_lag %||% TRK_DEFAULTS$trk_uni_min_lag))
  max_lag <- max(min_lag + 2L, as.integer(input$trk_uni_max_lag %||% TRK_DEFAULTS$trk_uni_max_lag))
  Ns <- pt$n_loc
  grid <- if (length(unique(Ns)) <= 8L) sort(unique(Ns))
          else sort(unique(pmax(5L, round(stats::quantile(Ns, seq(0, 1, length.out = 8L),
                                                          names = FALSE)))))
  nsim <- max(50L, min(500L, as.integer(input$trk_conf_ci_nsim %||% TRK_DEFAULTS$trk_conf_ci_nsim)))

  # ── D grid ────────────────────────────────────────────────────────────────
  # The free CI depends on x = σ_loc²/(D·t_c), so the null is matched on each
  # track's OWN D as well as its length - otherwise the slow half of a
  # heterogeneous population is measured against a null built for the fast half and
  # comes out "confined". The grid collapses to a single point (the cheap path, and
  # identical to matching on N alone) when that dependence cannot bite: either the
  # noise is negligible at the fastest tracks, or the per-track D barely varies.
  Dtr <- pt$D_um2_s * 1e3                       # per-track D, nm²/ms
  ok_D <- is.finite(Dtr) & Dtr > 0
  Dtr[!ok_D] <- D_use
  x_of <- function(D) if (sig > 0) sig^2 / (D * t_c) else 0
  D_rng <- if (any(ok_D)) stats::quantile(Dtr[ok_D], c(0.05, 0.95), names = FALSE)
           else c(D_use, D_use)
  needs_D_grid <- sig > 0 && is.finite(D_rng[1]) && D_rng[1] > 0 &&
                  x_of(D_rng[2]) > 0.3 &&        # noise matters somewhere
                  (D_rng[2] / D_rng[1]) > 2      # and D actually varies
  D_grid <- if (needs_D_grid)
              exp(seq(log(D_rng[1]), log(D_rng[2]), length.out = 4L))
            else D_use
  ncell <- length(grid) * length(D_grid)
  # Noise for the null: when σ_loc comes from the photon model, resample the
  # per-localization values (the MSD offset is set by mean σ², which a long-tailed
  # photon count pushes above median σ²). A manual or CVE/MSD σ stays static.
  sig_null <- if (identical(trk_sigma_loc_src(), "photon budget") && !is.null(b) &&
                  sum(is.finite(b$sig_phot)) >= 50L) b$sig_phot[is.finite(b$sig_phot)]
              else sig
  # The simulation is ~85% of the wait, so it owns 0 -> 0.85 of the bar and the
  # per-track matching that follows owns the rest. sim_free_ci_null() reports a
  # CUMULATIVE fraction, so it is set, not added: incrementing would run the bar
  # to 100% at the end of the simulation and leave the matching looking hung.
  nul <- withProgress(message = "Confinement: simulating the free null…", value = 0,
    sim_free_ci_null(grid, D_grid = D_grid, n_sim = nsim, t_c = t_c, sigma = sig_null,
                     ndim = ndim, gap_p = gap_p, min_lag = min_lag, max_lag = max_lag,
                     progress = function(f) setProgress(0.85 * f,
                       detail = sprintf("%d of %d (length, D) cells", round(f * ncell), ncell)),
                     R_blur = trk_blur()$R_msd))   # same blur coefficient as the MSD fits
  # nearest cell in (log N, log D)
  idx <- withProgress(message = "Confinement: matching tracks to the null…", value = 0.85, {
    ln <- log(nul$n); lD <- log(nul$D)
    vapply(seq_len(nrow(pt)), function(i) {
      which.min((ln - log(Ns[i]))^2 + (lD - log(Dtr[i]))^2)
    }, integer(1))
  })
  pt$null_D  <- nul$D[idx]
  pt$ci_lo05 <- nul$lo05[idx]; pt$ci_med <- nul$med[idx]; pt$ci_hi95 <- nul$hi95[idx]
  pt$ci_class <- ifelse(pt$conf_index < pt$ci_lo05, "confined",
                 ifelse(pt$conf_index > pt$ci_hi95, "directed", "free-like"))
  # Resolution floor: how far a free molecule travels over the median track duration.
  durs <- tapply(df$tim, df$tid, function(v) diff(range(v)))
  T_med <- stats::median(as.numeric(durs), na.rm = TRUE)
  attr(pt, "resolv") <- sqrt(2 * ndim * D_use * T_med)
  attr(pt, "D_use") <- D_use; attr(pt, "sigma") <- sig
  attr(pt, "T_med") <- T_med; attr(pt, "null") <- nul; attr(pt, "nsim") <- nsim
  attr(pt, "D_matched") <- length(D_grid) > 1L
  attr(pt, "x_noise") <- x_of(D_use)
  attr(pt, "sigma_info") <- .sigma_stamp(sig, note = if (length(sig_null) > 1L)
    "the null resamples the per-localization photon values" else NULL)
  pt
})

output$trk_conf_ci_summary <- renderUI({
  req(isTRUE(input$trk_show_conf))
  pt <- tryCatch(conf_ci_bundle(), error = function(e) NULL)
  if (is.null(pt) || !nrow(pt)) return(NULL)
  f_c <- mean(pt$ci_class == "confined"); f_f <- mean(pt$ci_class == "free-like")
  f_d <- mean(pt$ci_class == "directed")
  rv2 <- attr(pt, "resolv")
  tags$div(style = "font-size:0.9em; margin:6px 0; padding:6px 8px; border-left:3px solid #1a8a1a; background:rgba(26,138,26,0.05);",
    tags$b("MSD confinement index vs a simulated free null"),
    tags$span(style = "color:#666; font-weight:normal;", " - the primary call on confinement"), tags$br(),
    tags$b(style = "color:#1a8a1a; font-size:1.02em;",
      sprintf("Confined %.0f%%  ·  free-like %.0f%%  ·  directed %.0f%%",
              100 * f_c, 100 * f_f, 100 * f_d)),
    tags$div(style = "color:#666; margin-top:2px;",
      sprintf("n = %d tracks · null: σ_loc = %.2f nm, D = %.3g µm²/s, %d free tracks per grid point, matched on localization count%s",
              nrow(pt), attr(pt, "sigma"), attr(pt, "D_use") / 1e3, attr(pt, "nsim"),
              if (isTRUE(attr(pt, "D_matched"))) " and per-track D" else "")),
    tags$div(style = "color:#666; margin-top:2px;",
      htmltools::HTML(sprintf(
        "Resolution floor <b>%.0f nm</b> in the median track (%.1f ms) - confinement larger than this is invisible.",
        rv2, attr(pt, "T_med")))),
    tags$div(style = "color:#888; margin-top:2px;",
      "~5% per tail by construction; weak for directed motion - use Directionality or Speed & Transport."))
})

output$trk_dl_conf_ci_csv <- downloadHandler(
  filename = function() paste0("confinement_index_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
  content  = function(file) write.csv(conf_ci_bundle(), file, row.names = FALSE)
)

# CI against localization count, with the simulated free band behind it. Plotting
# against N makes the N-matching visible: the band narrows as tracks get longer,
# which is exactly why a single fixed threshold would be wrong.
output$trk_conf_ci_plot <- renderPlotly({
  req(isTRUE(input$trk_show_conf))
  pt <- tryCatch(conf_ci_bundle(), error = function(e) NULL)
  shiny::validate(shiny::need(!is.null(pt) && nrow(pt) > 0,
    "Run the MSD, then click 'Run confinement index'."))
  dark <- isTRUE(isolate(input$trk_dark_mode))
  bg <- if (dark) "#000000" else "#FFFFFF"; fg <- if (dark) "#FFFFFF" else "#333333"
  nul <- attr(pt, "null")
  # When the null is matched on D too there is a band per D value; show the one
  # nearest the ensemble D so the plot stays a single readable band.
  if (!is.null(nul) && length(unique(nul$D)) > 1L) {
    Dref <- unique(nul$D)[which.min(abs(log(unique(nul$D)) - log(attr(pt, "D_use"))))]
    nul <- nul[nul$D == Dref, , drop = FALSE]
  }
  cols <- c(confined = "#2c7fb8", `free-like` = "#9e9e9e", directed = "#e31a1c")
  fig <- plot_ly()
  if (!is.null(nul) && nrow(nul) > 1) {
    o <- order(nul$n)
    fig <- fig %>%
      add_trace(x = nul$n[o], y = nul$hi95[o], type = "scatter", mode = "lines",
                line = list(color = "rgba(120,120,120,0)"), showlegend = FALSE,
                hoverinfo = "skip") %>%
      add_trace(x = nul$n[o], y = nul$lo05[o], type = "scatter", mode = "lines",
                fill = "tonexty", fillcolor = "rgba(150,150,150,0.18)",
                line = list(color = "rgba(120,120,120,0)"),
                name = if (isTRUE(attr(pt, "D_matched")))
                         "free 5–95% (at the ensemble D)" else "free 5–95% (simulated)",
                hoverinfo = "skip") %>%
      add_trace(x = nul$n[o], y = nul$med[o], type = "scatter", mode = "lines",
                line = list(color = "#888", dash = "dash", width = 2),
                name = "free median", hoverinfo = "skip")
  }
  for (cl in names(cols)) {
    idx <- pt$ci_class == cl
    if (!any(idx)) next
    fig <- fig %>% add_trace(
      x = pt$n_loc[idx], y = pt$conf_index[idx], type = "scatter", mode = "markers",
      marker = list(size = 6, color = cols[[cl]], opacity = 0.75,
                    line = list(width = 0.3, color = "rgba(0,0,0,0.3)")),
      name = sprintf("%s (%d)", cl, sum(idx)),
      text = sprintf("TID %s<br>%d locs · CI %.3f<br>%s",
                     pt$tid[idx], pt$n_loc[idx], pt$conf_index[idx], cl),
      hovertemplate = "%{text}<extra></extra>")
  }
  fig %>% layout(
    xaxis = list(title = "localizations in track", color = fg, type = "log"),
    yaxis = list(title = "confinement index (CI)", color = fg),
    paper_bgcolor = bg, plot_bgcolor = bg, font = list(color = fg),
    legend = list(font = list(color = fg, size = 10), orientation = "h",
                  x = 0, xanchor = "left", y = -0.34, yanchor = "top"),
    margin = list(l = 55, r = 20, t = 10, b = 110),
    uirevision = "conf_ci")
})


# ============================================================================
# trk_disp_bundle() - NET DISPLACEMENT vs DIFFUSION  (Step 1 of the geometry test)
# ----------------------------------------------------------------------------
# Per track: the NET start→end displacement R_net = |r(T) − r(0)|, the max excursion
# from the start, the (sampling-dependent) path length, and the duration T. Each is
# compared to the FREE-diffusion expectation at the ensemble D (from the CVE panel):
# for 2D Brownian motion ⟨R_net²⟩ = 4·D·T, and R_net itself is Rayleigh-distributed
# (scale √(2DT)), so a single track has a wide spread even when free.
#
# WHY 2D + WHY THIS IS ONLY STEP 1: displacement is measured in the 2D image plane
# (x, y). This is the INFINITE-PLANE null - it does NOT account for a finite/curved
# surface, which itself caps how far a FREE molecule can get. So "sub-diffusive" here
# means "confined OR geometrically constrained"; Step 2 (simulate free diffusion ON a
# chosen geometry at this D) separates the two. Path length is reported but NEVER used
# to classify - like path speed it scales with sampling rate, not the molecule.
#
# NOISE: endpoints carry σ_loc, so ⟨R_net²⟩ gains 4σ² (2 endpoints × 2 axes). That
# (and the extreme-value analogue for Rmax) is subtracted when σ_loc is set.
# ============================================================================
trk_disp_bundle <- eventReactive(input$trk_disp_run, {
  req(isTRUE(input$trk_show_disp))
  df <- trk(); req(nrow(df) > 0)
  min_pts <- max(3L, as.integer(input$trk_disp_min_pts %||% TRK_DEFAULTS$trk_disp_min_pts))

  # Dimensionality: net displacement (and the free null) in 2D or 3D. "Auto" uses z when
  # it carries real signal. 2D data → 2D; 3D data → 3D (Rayleigh vs Maxwell null). The
  # geometry null (Step 2) then projects (2D) or stays 3D to match.
  z_has <- .has_real_z(df)
  dmode <- input$trk_disp_dims %||% TRK_DEFAULTS$trk_disp_dims
  use3d <- (identical(dmode, "auto") && z_has) || identical(dmode, "3d")
  if (use3d && !z_has) use3d <- FALSE                       # asked for 3D but no z signal
  d_obs <- if (use3d) 3L else 2L
  # ── Degrees of freedom of the MOTION vs of the OBSERVATION ─────────────────
  # d_obs is how many coordinates are recorded; d_mot is how many the molecule can
  # actually move in. They differ for a molecule confined to a SURFACE (a membrane,
  # a cell envelope) but tracked in 3D: 2 DOF observed on 3 axes.
  # Two consequences, both handled here:
  #   1. A per-axis estimator (CVE, or MSD with ndim = 3) spreads 2-DOF motion over
  #      3 axes, so it reports D_meas = D_true · d_mot/d_obs - exactly 2/3 for a
  #      surface in 3D (simulated: 0.675). Undo it: D_mot = D_meas · d_obs/d_mot.
  #   2. R_net² follows χ² with d_mot degrees of freedom, not d_obs.
  # Getting (1) right but not (2) leaves the MEAN correct while the BAND is the wrong
  # shape: on simulated free-on-surface data the old d_obs-only form gave
  # 11.4% sub / 81.4% free / 7.2% super instead of 5/90/5; with d_mot it gives
  # 5.2/90.0/4.8. When d_mot == d_obs (2D data, or true volume diffusion) both terms
  # cancel and this reduces exactly to the previous expression.
  # Degrees of freedom of the motion, from the shared Sample geometry setting.
  d_mot <- {
    mo <- input$trk_motion %||% TRK_DEFAULTS$trk_motion
    if (identical(mo, "filament")) 1L
    else if (identical(mo, "surface")) min(2L, d_obs)
    else d_obs
  }

  keep <- is.finite(df$x) & is.finite(df$y) & is.finite(df$tim) &
          (if (use3d) is.finite(df$z) else TRUE)
  df <- df[keep, ]

  # ── Which D sets the free null? ────────────────────────────────────────────
  # The null is √(2·d·D·T), so D directly sets the verdict: a larger D raises the band
  # and flags MORE tracks sub-diffusive. It must be the SHORT-TIME (local) mobility.
  #   CVE      - short-lag, noise-corrected (Vestergaard), model-free. DEFAULT, because
  #              it is measured before any confinement bends the MSD.
  #   MSD fit  - fitted over the lag window. If the population IS confined the MSD has
  #              already bent down inside that window, so this D is DEFLATED and the
  #              test loses power (it under-calls confinement). Offered for comparison.
  #   Manual   - a literature/known value, or for sensitivity testing.
  dsrc  <- input$trk_disp_dsrc %||% TRK_DEFAULTS$trk_disp_dsrc
  b     <- tryCatch(trk_cve_bundle(), error = function(e) NULL)
  D_cve <- if (!is.null(b)) b$D_cve else NA_real_
  D_msd <- {
    ub <- tryCatch(msd_unified_bundle(), error = function(e) NULL)
    fr <- if (!is.null(ub)) tryCatch(msd_unified_modelfit(ub), error = function(e) NULL) else NULL
    if (!is.null(fr) && is.finite(fr$D)) fr$D else NA_real_
  }
  D_man <- suppressWarnings(as.numeric(input$trk_disp_dman %||% TRK_DEFAULTS$trk_disp_dman)) * 1e3  # µm²/s → nm²/ms
  D_use <- switch(dsrc, msd = D_msd, manual = D_man, D_cve)
  d_lab <- switch(dsrc, msd = "MSD fit", manual = "manual", "CVE")
  # In-plane (motion) D. CVE and the MSD fit are per-axis over d_obs recorded axes, so
  # for 2-DOF motion tracked in 3D they read low by d_mot/d_obs and are rescaled here.
  # A MANUAL D is taken to be the true motion D already, so it is passed through.
  D_mot <- if (identical(dsrc, "manual")) D_use else D_use * d_obs / d_mot
  shiny::validate(shiny::need(is.finite(D_use) && D_use > 0, switch(dsrc,
    msd    = "No MSD D yet - press Plot MSD in the MSD panel (or switch the D source to CVE).",
    manual = "Enter a positive D (µm²/s) for the manual D source.",
    "Need an ensemble D first - run Analysis → Data Quality & Calibration (Estimate σloc & quality).")))
  sig    <- tryCatch(trk_sigma_loc_eff(), error = function(e) NA_real_)
  # Endpoint noise: 2 endpoints × d_obs axes × σ².
  noise2 <- if (is.finite(sig) && sig > 0) 2 * d_obs * sig^2 else 0

  sp <- split(df, df$tid)
  rows <- withProgress(message = "Displacement vs diffusion…", value = 0,
    lapply(seq_along(sp), function(i) {
      if (i %% 100L == 0L) incProgress(100 / length(sp))
      tb <- sp[[i]]; if (nrow(tb) < min_pts) return(NULL)
      tb <- tb[order(tb$tim), ]; ne <- nrow(tb)
      T_ms <- tb$tim[ne] - tb$tim[1]
      if (!is.finite(T_ms) || T_ms <= 0) return(NULL)
      dz2 <- if (use3d) (tb$z[ne] - tb$z[1])^2 else 0
      Rnet <- sqrt(max(0, (tb$x[ne] - tb$x[1])^2 + (tb$y[ne] - tb$y[1])^2 + dz2 - noise2))
      ex2  <- (tb$x - tb$x[1])^2 + (tb$y - tb$y[1])^2 + if (use3d) (tb$z - tb$z[1])^2 else 0
      Rmax <- sqrt(max(0, max(ex2) - noise2))
      pl   <- sum(sqrt(diff(tb$x)^2 + diff(tb$y)^2 + if (use3d) diff(tb$z)^2 else 0))
      # ⟨R_net²⟩ = 2·d_mot·D_mot·T. With d_mot = d_obs this is √(4DT) in 2D and
      # √(6DT) in 3D; for a surface tracked in 3D it becomes √(4·D_mot·T).
      free_rms <- sqrt(2 * d_mot * D_mot * T_ms)
      data.frame(tid = names(sp)[i], n_pts = ne, dur_ms = T_ms,
                 Rnet_nm = Rnet, Rmax_nm = Rmax, path_nm = pl,
                 free_rms_nm = free_rms,
                 rho = if (free_rms > 0) Rnet / free_rms else NA_real_,
                 stringsAsFactors = FALSE)
    }))
  rows <- rows[!vapply(rows, is.null, logical(1))]
  shiny::validate(shiny::need(length(rows) > 0, "Too few tracks - lower 'Min points per TID'."))
  res <- do.call(rbind, rows)

  # Per-track free (infinite-plane) 5/95 band and two-tailed classification. R_net² for
  # free diffusion is 2·D_mot·T·χ²(d_mot) - the χ² degrees of freedom follow the MOTION,
  # not the number of recorded coordinates (see the d_mot note above).
  sc2 <- 2 * D_mot * res$dur_ms
  res$free_lo05 <- sqrt(sc2 * stats::qchisq(0.05, d_mot))
  res$free_hi95 <- sqrt(sc2 * stats::qchisq(0.95, d_mot))
  res$class <- ifelse(res$Rnet_nm < res$free_lo05, "sub-diffusive",
               ifelse(res$Rnet_nm > res$free_hi95, "super-diffusive", "free-consistent"))
  # Self-contained columns so the CSV carries everything Step 2 needs.
  res$D_ens_um2_s  <- D_use / 1e3
  res$D_motion_um2_s <- D_mot / 1e3
  res$sigma_loc_nm <- if (is.finite(sig)) sig else NA_real_
  res$dims         <- d_obs
  res$motion_dof   <- d_mot
  res$D_source     <- d_lab
  attr(res, "D_use") <- D_use; attr(res, "sigma") <- sig; attr(res, "d_obs") <- d_obs
  attr(res, "d_mot") <- d_mot; attr(res, "D_mot") <- D_mot
  attr(res, "d_lab") <- d_lab
  attr(res, "D_cve") <- D_cve; attr(res, "D_msd") <- D_msd
  attr(res, "sigma_info") <- .sigma_stamp(sig)
  res
})

# Echo the shared geometry inside the displacement panel, so the Step-2 band states
# which shape it is using without owning a second copy of the control.
output$trk_disp_geom_echo <- renderUI({
  g <- input$trk_geom_type %||% TRK_DEFAULTS$trk_geom_type
  params <- list(R = input$trk_disp_R, r = input$trk_disp_r, Lcyl = input$trk_disp_L,
                 a = input$trk_disp_a, b = input$trk_disp_b, c = input$trk_disp_c,
                 mtilt = input$trk_disp_mtilt)
  dims <- tryCatch(.disp_geom_dims(.geom_for_null(g), params), error = function(e) "")
  htmltools::div(style = "font-size:0.85em; margin:2px 0 4px;",
    htmltools::tags$b("Geometry: "), .geom_label(g),
    if (nzchar(dims)) htmltools::span(style = "color:#666;", sprintf("  (%s)", dims)) else NULL,
    if (identical(.geom_for_null(g), "free"))
      htmltools::div(style = "color:#b26a00;",
        "No bounding shape, so this band matches the unbounded null. Set a shape in Sample geometry to test a closed cell.")
    else NULL)
})

# ── Step 2: GEOMETRY-AWARE free null (free diffusion ON the chosen geometry) ───────
# Replaces the infinite-plane band with a band simulated from free diffusion on the
# selected shape (sphere / rod / ellipsoid / membrane / free - the Diffusion Simulator's
# geometries), so "sub-diffusive" can be split into TRULY confined vs merely
# geometrically constrained. The measured 2D D is a PROJECTED D (< the true surface D on a
# curved shape), so sim_free_geom_rnet empirically calibrates that per-geometry factor and
# drives the null at the matching surface D → both bands agree at short T and the
# closed-surface band diverges only where finite size + curvature bite. Same physics as
# the Diffusion Simulator.
trk_disp_geom_null <- reactive({
  if (!isTRUE(input$trk_disp_geom)) return(NULL)
  res <- trk_disp_bundle(); req(nrow(res) > 0)
  D_use <- attr(res, "D_use")
  geom  <- .geom_for_null(input$trk_geom_type %||% TRK_DEFAULTS$trk_geom_type)
  params <- list(
    R    = max(1, as.numeric(input$trk_disp_R    %||% TRK_DEFAULTS$trk_disp_R)),
    r    = max(1, as.numeric(input$trk_disp_r    %||% TRK_DEFAULTS$trk_disp_r)),
    Lcyl = max(1, as.numeric(input$trk_disp_L    %||% TRK_DEFAULTS$trk_disp_L)),
    a    = max(1, as.numeric(input$trk_disp_a    %||% TRK_DEFAULTS$trk_disp_a)),
    b    = max(1, as.numeric(input$trk_disp_b    %||% TRK_DEFAULTS$trk_disp_b)),
    c    = max(1, as.numeric(input$trk_disp_c    %||% TRK_DEFAULTS$trk_disp_c)),
    mtilt = max(0, as.numeric(input$trk_disp_mtilt %||% TRK_DEFAULTS$trk_disp_mtilt)))
  cyc   <- tryCatch(isolate(trk_cycle_time()), error = function(e) NA_real_)
  if (!is.finite(cyc) || cyc <= 0) cyc <- stats::median(diff(sort(res$dur_ms))) # crude fallback
  if (!is.finite(cyc) || cyc <= 0) return(NULL)
  # Match the real DURATION distribution (R_net depends on T) over its FULL range.
  # This used to be truncated at the 98th percentile to bound the per-cycle loop, which
  # left the longest tracks - the ones carrying the most information about confinement -
  # judged against a flat-extrapolated band. Cost scales with the number of simulation
  # steps (longest duration / step), so instead of dropping those tracks the SIMULATION
  # STEP is coarsened: Brownian statistics are step-size invariant, so a larger step with
  # SD = √(2·D·dt) gives the same R_net distribution provided the step stays small
  # compared with the geometry. Measured on a rod (r = 400 nm, D = 200 nm²/ms): stepping
  # at 1×/5×/10×/25× the cycle gave median R_net 292/305/306/302 nm and p95 753/782/787/786
  # while the run time fell from 2.6 s to 0.09 s.
  # Simulate on a LOG-SPACED GRID spanning the observed durations, with replicates at
  # each point, rather than resampling the empirical durations. Resampling plus
  # equal-count binning put the last bin's MEDIAN far below the longest track, so the
  # band visibly stopped short of the data and the few longest tracks - the most
  # informative ones - were classified against a flat extrapolation. On a grid the
  # band is defined at the maximum duration by construction.
  dv <- res$dur_ms[is.finite(res$dur_ms) & res$dur_ms > 0]
  if (!length(dv)) return(NULL)
  d_rng <- range(dv)
  if (d_rng[1] <= 0 || !is.finite(d_rng[2])) return(NULL)
  n_grid_t <- 24L; n_rep <- 200L
  grid_t <- if (d_rng[2] / d_rng[1] > 1.05)
              exp(seq(log(d_rng[1]), log(d_rng[2]), length.out = n_grid_t))
            else rep(d_rng[2], 1L)
  durs  <- rep(grid_t, each = n_rep)
  T_max <- max(durs, na.rm = TRUE)
  step_budget <- 4000L
  dt_sim <- max(cyc, T_max / step_budget)
  # Keep the step well below the smallest geometry dimension, or the tangent-step-and-
  # reproject operator loses accuracy. Unbounded shapes have no such limit.
  gmin <- switch(geom, sphere = params$R, rod = params$r,
                 ellipsoid = min(params$a, params$b, params$c), Inf)
  if (is.finite(gmin) && gmin > 0 && is.finite(D_use) && D_use > 0) {
    dt_geom <- (0.05 * gmin)^2 / (2 * D_use)     # step SD ≤ 5% of the smallest dimension
    dt_sim  <- min(dt_sim, max(cyc, dt_geom))
  }
  # 2D data → simulate the 3D geometry and PROJECT (the cell is 3D even if tracking is 2D);
  # 3D data → keep all three axes.
  proj2d <- (attr(res, "d_obs") %||% 2L) == 2L
  # Same "Molecule moves" choice that sets the Step-1 degrees of freedom.
  # "surface" and "filament" both mean the molecule is on the boundary; only a
  # volume-filling molecule reflects around inside it.
  gmode  <- if (identical(input$trk_motion %||% TRK_DEFAULTS$trk_motion, "volume")) "volume" else "surface"
  rnet <- withProgress(message = "Simulating free-on-geometry null…", value = 0.4,
                       sim_free_geom_rnet(durs, D_use, geom, params, dt_sim, proj2d = proj2d,
                                          mode = gmode))
  if (!length(rnet)) return(NULL)
  D_surf <- attr(rnet, "D_surf") %||% D_use
  # One 5/50/95 summary per grid duration - the band therefore spans the full range of
  # observed track durations, ending exactly at the longest track.
  rs  <- as.numeric(rnet)
  gid <- rep(seq_along(grid_t), each = n_rep)[seq_along(rs)]
  bt  <- grid_t[sort(unique(gid))]
  blo <- as.numeric(tapply(rs, gid, stats::quantile, 0.05))
  bmd <- as.numeric(tapply(rs, gid, stats::quantile, 0.50))
  bhi <- as.numeric(tapply(rs, gid, stats::quantile, 0.95))
  # Classify each REAL track against the geometry null at its own duration.
  lo_r <- stats::approx(bt, blo, xout = res$dur_ms, rule = 2)$y
  hi_r <- stats::approx(bt, bhi, xout = res$dur_ms, rule = 2)$y
  gclass <- ifelse(res$Rnet_nm < lo_r, "confined",
            ifelse(res$Rnet_nm > hi_r, "directed", "free-on-geometry"))
  list(bt = bt, blo = blo, bmd = bmd, bhi = bhi, gclass = gclass, D_surf = D_surf,
       geom = geom, dims = .disp_geom_dims(geom, params), cyc = cyc, gmode = gmode)
})


output$trk_disp_plot <- renderPlotly({
  req(isTRUE(input$trk_show_disp))
  res <- trk_disp_bundle(); req(nrow(res) > 0)
  dark <- isTRUE(isolate(input$trk_dark_mode))
  bg <- if (dark) "#000000" else "#FFFFFF"; fg <- if (dark) "#FFFFFF" else "#333333"
  D_use <- attr(res, "D_use"); d_obs <- attr(res, "d_obs") %||% 2L
  # Band uses the MOTION dof / motion D (see trk_disp_bundle); identical to the old
  # expression whenever d_mot == d_obs.
  d_mot <- attr(res, "d_mot") %||% d_obs
  D_mot <- attr(res, "D_mot") %||% D_use
  tg  <- seq(max(1e-6, min(res$dur_ms)), max(res$dur_ms), length.out = 120)
  sc2 <- 2 * D_mot * tg   # R_net² = 2·D_mot·T·χ²(d_mot)
  med <- sqrt(sc2 * stats::qchisq(0.5, d_mot))
  lo  <- sqrt(sc2 * stats::qchisq(0.05, d_mot))
  hi  <- sqrt(sc2 * stats::qchisq(0.95, d_mot))
  gn  <- trk_disp_geom_null()   # NULL unless the Step-2 geometry option is on
  pt_class <- if (!is.null(gn)) gn$gclass else res$class
  fig <- plot_ly() %>%
    add_trace(x = tg, y = hi, type = "scatter", mode = "lines",
              line = list(color = "rgba(120,120,120,0)"), showlegend = FALSE, hoverinfo = "skip") %>%
    add_trace(x = tg, y = lo, type = "scatter", mode = "lines", fill = "tonexty",
              fillcolor = "rgba(150,150,150,0.15)", line = list(color = "rgba(120,120,120,0)"),
              name = "free 5–95% (unbounded)", hoverinfo = "skip") %>%
    add_trace(x = tg, y = med, type = "scatter", mode = "lines",
              line = list(color = "#888", dash = "dash", width = 2),
              name = "free median (unbounded)", hoverinfo = "skip")
  if (!is.null(gn)) {
    # Free-on-E.coli band (green): the geometry-aware null.
    fig <- fig %>%
      add_trace(x = gn$bt, y = gn$bhi, type = "scatter", mode = "lines",
                line = list(color = "rgba(27,158,90,0)"), showlegend = FALSE, hoverinfo = "skip") %>%
      add_trace(x = gn$bt, y = gn$blo, type = "scatter", mode = "lines", fill = "tonexty",
                fillcolor = "rgba(27,158,90,0.18)", line = list(color = "rgba(27,158,90,0)"),
                name = "free 5–95% (on geometry)", hoverinfo = "skip") %>%
      add_trace(x = gn$bt, y = gn$bmd, type = "scatter", mode = "lines",
                line = list(color = "#1B9E5A", width = 2),
                name = "free median (on geometry)", hoverinfo = "skip")
  }
  fig <- fig %>%
    add_trace(x = res$dur_ms, y = res$Rnet_nm, type = "scatter", mode = "markers",
              marker = list(size = 6, color = .disp_class_col(pt_class),
                            line = list(width = 0.3, color = "rgba(0,0,0,0.3)")),
              text = sprintf("TID %s<br>dur %.1f ms · %d pts<br>R_net %.0f nm · ρ %.2f<br>%s",
                             res$tid, res$dur_ms, res$n_pts, res$Rnet_nm, res$rho, pt_class),
              hovertemplate = "%{text}<extra></extra>", name = "tracks")
  fig %>%
    layout(
      xaxis = list(title = "track duration T (ms)", color = fg, zeroline = FALSE),
      yaxis = list(title = "net displacement R_net (nm)", color = fg, zeroline = FALSE),
      paper_bgcolor = bg, plot_bgcolor = bg, font = list(color = fg),
      legend = list(font = list(color = fg, size = 10), orientation = "h", y = -0.2),
      title = list(text = "Net displacement vs free-diffusion expectation",
                   font = list(color = fg, size = 12)),
      uirevision = "disp")
})

output$trk_disp_summary <- renderUI({
  req(isTRUE(input$trk_show_disp))
  res <- trk_disp_bundle(); req(nrow(res) > 0)
  D_use <- attr(res, "D_use"); sig <- attr(res, "sigma")
  n <- nrow(res)
  f_sub <- mean(res$class == "sub-diffusive"); f_free <- mean(res$class == "free-consistent")
  f_sup <- mean(res$class == "super-diffusive")
  T_med <- stats::median(res$dur_ms)
  d_obs <- attr(res, "d_obs") %||% 2L
  resolv <- sqrt(2 * d_obs * D_use * T_med)   # RMS free extent; confinement below it is invisible
  htmltools::div(style = "font-size:0.9em; line-height:1.55;",
    htmltools::HTML(sprintf(
      "<b>%d</b> tracks · <b>%dD</b> net displacement · D = <b>%.3g µm²/s</b> (%s)%s · median duration <b>%.1f ms</b>.",
      n, attr(res, "d_obs") %||% 2L, D_use / 1e3, attr(res, "d_lab") %||% "CVE",
      if (is.finite(sig)) sprintf(" · σ<sub>loc</sub> = %.2f nm", sig) else "", T_med)),
    {
      # Only worth stating when the motion dof differs from the recorded coordinates.
      dm <- attr(res, "d_mot") %||% d_obs; Dm <- attr(res, "D_mot") %||% D_use
      if (dm != d_obs)
        htmltools::div(style = "color:#555; margin-top:3px;",
          htmltools::HTML(sprintf(
            "Surface motion: %d DOF on %d axes &rarr; in-plane D = <b>%.3g µm²/s</b>, band uses χ²(%d).",
            dm, d_obs, Dm / 1e3, dm)))
      else NULL
    },
    {
      # Cross-check: if the two estimators disagree, the verdict is D-sensitive.
      Dc <- attr(res, "D_cve"); Dm <- attr(res, "D_msd")
      if (is.finite(Dc) && is.finite(Dm) && Dc > 0 && Dm > 0 && abs(log(Dm / Dc)) > log(1.3))
        htmltools::div(style = "color:#B8860B; margin-top:3px; font-size:0.93em;",
          htmltools::HTML(sprintf(
            "&#9888; CVE D = %.3g vs MSD D = %.3g \u00b5m\u00b2/s (%.1f\u00d7 apart), so the split below depends on which D is used.",
            Dc / 1e3, Dm / 1e3, max(Dm, Dc) / min(Dm, Dc))))
      else NULL
    },
    htmltools::br(),
    htmltools::HTML(sprintf(
      "Infinite-plane null: <b style='color:#2c7fb8;'>%.0f%% sub-diffusive</b> · <b style='color:#9e9e9e;'>%.0f%% free</b> · <b style='color:#e31a1c;'>%.0f%% super-diffusive</b> (~5%% per tail by construction)",
      100 * f_sub, 100 * f_free, 100 * f_sup)),
    {
      # This panel corroborates the confinement index; it does not compete with it.
      # The two ask different questions of the same track - the CI reads the shape of
      # the MSD curve across lags, this reads how far the track actually travelled -
      # so their percentages are not meant to match. What is informative is how far
      # they agree track by track, which is only shown once both have been run.
      pt <- tryCatch(conf_ci_bundle(), error = function(e) NULL)
      if (is.null(pt) || !nrow(pt) || !all(c("tid", "ci_class") %in% names(pt)))
        htmltools::div(style = "color:#555; margin-top:4px;",
          htmltools::HTML("Run the <b>confinement index</b>, the primary call, to see how far the two agree."))
      else {
        m <- merge(res[, c("tid", "class")], pt[, c("tid", "ci_class")], by = "tid")
        if (!nrow(m)) NULL else {
          ci_conf <- m$ci_class == "confined"
          both    <- sum(ci_conf & m$class == "sub-diffusive")
          htmltools::div(style = "color:#555; margin-top:4px;",
            htmltools::HTML(sprintf(
              "Of the <b>%d</b> tracks the confinement index calls confined, <b>%d</b> (%.0f%%) also fall below this band (%d tracks scored by both).",
              sum(ci_conf), both,
              if (sum(ci_conf)) 100 * both / sum(ci_conf) else 0, nrow(m))))
        }
      }
    },
    htmltools::div(style = "color:#555; margin-top:5px;",
      htmltools::HTML(sprintf(
        "Resolution floor: <b>%.0f nm</b> RMS free extent. Confinement on a larger scale is invisible.", resolv))),
    {
      # ONE ensemble D sets the band for every track, so a population with a wide
      # spread of per-track D reads as sub-diffusive at the slow end even when every
      # molecule is free. Quantified on simulated all-free populations: the
      # sub-diffusive call runs ~5.5% at no spread, ~6.9% at sd(ln D)=0.6 and ~12%
      # at sd(ln D)=1.0. Warn using the per-track D spread when the MSD has been run.
      td <- tryCatch(trk_unified_tdiff(), error = function(e) NULL)
      Dv <- if (!is.null(td) && !is.null(td$pertrack)) td$pertrack$D_um2_s else NULL
      Dv <- Dv[is.finite(Dv) & Dv > 0]
      sdl <- if (length(Dv) >= 20L) stats::sd(log10(Dv)) else NA_real_
      if (is.finite(sdl) && sdl > 0.35)
        htmltools::div(style = "color:#B8860B; margin-top:5px; font-size:0.95em;",
          htmltools::HTML(sprintf(
            "&#9888; Wide per-track D (sd log<sub>10</sub>D = %.2f, ~%.0f\u00d7 across the 5\u201395%% range): slow free tracks fall below the band, so the sub-diffusive %% is an upper bound.",
            sdl, 10^(2 * 1.645 * sdl))))
      else NULL
    },
    {
      gn <- trk_disp_geom_null()
      if (is.null(gn))
        htmltools::div(style = "color:#B8860B; margin-top:5px; font-size:0.95em;",
          htmltools::HTML("&#9888; Infinite-plane null: free diffusion on a curved cell already reads mildly sub-diffusive. <b>Step 2</b> separates the two. A low directed % does not rule out transport."))
      else {
        gc <- gn$gclass
        fc <- mean(gc == "confined"); fg2 <- mean(gc == "free-on-geometry"); fd <- mean(gc == "directed")
        geom_lab <- c(sphere = "sphere", rod = "spherocylinder (rod)", ellipsoid = "ellipsoid",
                      membrane = "membrane plane", free = "free (unbounded)")[gn$geom]
        unbounded <- gn$geom %in% c("free", "membrane")
        if (!unbounded)
          geom_lab <- paste0(geom_lab, if (identical(gn$gmode, "volume"))
                             " - inside the volume" else " - on the surface")
        htmltools::tagList(
          htmltools::div(style = "margin-top:7px; padding-top:5px; border-top:1px solid #ddd;",
            htmltools::HTML(sprintf(
              "<b>Step 2, free diffusion on %s</b> (%s; surface D = %.3g \u00b5m\u00b2/s): <b style='color:#2c7fb8;'>%.0f%% confined</b> · <b style='color:#1B9E5A;'>%.0f%% free-on-geometry</b> · <b style='color:#e31a1c;'>%.0f%% directed</b>.",
              geom_lab, gn$dims, gn$D_surf / 1e3, 100 * fc, 100 * fg2, 100 * fd))),
          htmltools::div(style = "color:#555; margin-top:4px;",
            htmltools::HTML(
              if (unbounded)
                sprintf("Unbounded geometry, so this matches the infinite-plane null (%.0f%% confined). Pick a closed surface to test finite size.", 100 * fc)
              else sprintf(
                "Of the %.0f%% called sub-diffusive on the plane, geometry keeps <b>%.0f%% as truly confined</b>%s",
                100 * f_sub, 100 * fc,
                if (fc <= 0.08) ", at the ~5% free baseline: consistent with free surface diffusion."
                else " - the excess over the ~5% free baseline."))),
          htmltools::div(style = "color:#555; margin-top:4px; font-size:0.94em;",
            htmltools::HTML(
              if ((attr(res, "d_obs") %||% 2L) == 2L)
                "2D data: the geometry is simulated in 3D and projected to the image plane."
              else "3D data: the geometry null is evaluated in 3D, no projection.")),
          htmltools::div(style = "color:#B8860B; margin-top:4px; font-size:0.94em;",
            htmltools::HTML("&#9888; The null reflects at the surface but does not model dark-cycle survivorship."))
        )
      }
    }
  )
})

output$trk_dl_disp_csv <- downloadHandler(
  filename = function() paste0("displacement_vs_D_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
  content  = function(file) write.csv(trk_disp_bundle(), file, row.names = FALSE)
)
