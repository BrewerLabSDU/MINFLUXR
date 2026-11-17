## ─────────────────────────────────────────────────────────────────────────────
## server/tracking/21_hmm_state_msd.R  -  panels: diffusive-state HMM; per-state MSD
## ─────────────────────────────────────────────────────────────────────────────
## Part of tracking_server(). Not a standalone script: its statements are evaluated in
## the server function's environment, in file-name order (server/parts.R).
##
## hmm_bundle (R/est_hmm.R) with its readouts, and the MSD of each decoded
## state.
##
## INTERFACE (generated: Rscript tools/server_map.R --headers; do not edit by hand)
##   DEFINES for other parts:
##     (nothing)
##   READS from other parts:
##     01_state_sigma: .sigma_stamp, trk_sigma_loc_eff, trk_sigma_loc_src
##     02_data_quality_cve: trk_cve_bundle
##     06_filter_pipeline: .a_in_use
##     07_cycle_time_sequence: trk_blur, trk_cluster_cycle, trk_cycle_time
##     11_data_frames: trk
##     15_msd: .msd_d_eff
##     17_jdd: jdd_bundle
## END INTERFACE
## ─────────────────────────────────────────────────────────────────────────────

# ============================================================================
# hmm_bundle() - DIFFUSION-STATE SEGMENTATION (continuous-time HMM, est_hmm.R)
# ============================================================================
#
# GOAL
#   Answer "does this population switch between diffusive states, and how fast?"
#   with K states (D_k) and switching RATES estimated by maximum likelihood over
#   all tracks at once. The estimator itself lives in est_hmm.R (pure functions,
#   tested against simulated ground truth); this reactive only gathers inputs,
#   calls it and packages the result for the plots.
#
# THE MODEL (full statement in the est_hmm.R header)
#   Hidden state per step, a continuous-time Markov chain with generator Q, so a
#   step of any Δt has transition matrix expm(Q·Δt): dark cycles are handled
#   exactly. Emission: r² over x,y is Gamma(1, 2·v) with the per-axis variance
#       v = 2·D_s·(Δt − 2·R·t_c) + σ_j² + σ_{j+1}²
#   i.e. diffusion, motion blur (trk_blur()) and the noise of both localizations.
#   D is therefore noise- and blur-corrected by construction.
#
# WHAT CHANGED AND WHY (2026-09-17)
#   The old version fitted a per-step matrix by Baum-Welch and counted switches
#   along the Viterbi path. Path counting absorbs short dwells: on simulated data
#   with dwells of 10 / 20 ms it read 15.6 / 30.2 ms (+50%), while this MLE read
#   9.9 / 18.8 ms. The Viterbi path is still computed, but ONLY to colour tracks
#   and to feed the per-state MSD.
#
# OBJECTS IT CALLS (inputs)
#   trk() ................. analysis-ready data (tid, x, y, tim in ms, eco)
#   trk_cycle_time() ...... t_c, for the blur offset and starting rates
#   trk_sigma_loc_eff() / trk_sigma_loc_src() ... the noise model: per-localization
#                           photon σ when the photon model is active, else one σ
#   trk_cve_bundle() ...... optional, supplies the emitter offset x for the photon σ
#   inputs: trk_hmm_run, trk_hmm_n_states, trk_hmm_min_pts, trk_hmm_n_ex,
#           trk_hmm_boot (bootstrap resamples, 0 = off)
#
# RETURNS a list (fields read by the plots, the summary and seg_msd_rv):
#   K, D_vals[K] (nm²/ms), fracs[K] (stationary occupancy), P_cycle (K×K
#   transition probabilities over one cycle), Qrate (per ms), dwell_ms[K],
#   fit (the full est_hmm result incl. CIs), boot (bootstrap CIs or NULL),
#   seg_list (per-track data.frame(tid, x, y, tim, state)), n_tracks, n_steps,
#   dt_med_ms, t_cycle_ms, sigma_loc_nm, sigma_mode, D_noise, noise_state[K],
#   resolv (per-state flags from the validated range), bic_other (BIC at K±1)
# ============================================================================
hmm_bundle <- eventReactive(input$trk_hmm_run, {
  req(isTRUE(input$trk_show_hmm))
  df <- trk(); req(nrow(df) > 0)
  K       <- max(2L, min(4L, as.integer(input$trk_hmm_n_states)))
  min_pts <- max(4L, as.integer(input$trk_hmm_min_pts))
  B       <- max(0L, min(200L, as.integer(input$trk_hmm_boot %||% TRK_DEFAULTS$trk_hmm_boot)))
  t_c     <- tryCatch(trk_cycle_time(), error = function(e) NA_real_)
  shiny::validate(shiny::need(is.finite(t_c) && t_c > 0,
    "Cycle time unknown - run Data Quality & Calibration first."))
  R_abs   <- trk_blur()$R_abs                         # blur offset in ms, shared (est_blur.R)

  # ── Noise model ──────────────────────────────────────────────────────────
  # The emission needs each step's noise variance. With the photon model active,
  # every localization gets its own σ from its photon count (the same model the
  # rest of the tab uses); otherwise one σ for all. Without any σ the slow state
  # would absorb the noise and read too fast, so that case is refused.
  sig_eff <- tryCatch(trk_sigma_loc_eff(), error = function(e) NA_real_)
  src     <- tryCatch(trk_sigma_loc_src(), error = function(e) "not set")
  # If Data Quality & Calibration has not published a σ yet, derive the photon-budget
  # value here from the same model rather than refusing: every input it needs (ECO,
  # L_TCP, the prefactor) is already on this tab.
  if (!is.finite(sig_eff) || sig_eff <= 0) {
    if ("eco" %in% names(df) && any(is.finite(df$eco) & df$eco > 0)) {
      L0 <- suppressWarnings(as.numeric(input$trk_L_tcp %||% DMAX_L_TCP_DEFAULT))
      if (!is.finite(L0) || L0 <= 0) L0 <- DMAX_L_TCP_DEFAULT
      n0 <- suppressWarnings(as.numeric(input$trk_beam_fwhm %||% TCP_BEAM_FWHM_NM))
      if (!is.finite(n0) || n0 <= 0) n0 <- TCP_BEAM_FWHM_NM
      a0 <- .a_in_use()
      ok0 <- is.finite(df$eco) & df$eco > 0
      sig_eff <- stats::median(.crlb_doughnut(L0, n0, df$eco[ok0], pref = a0,
                                              pattern = input$trk_tcp_pattern %||% TCP_PATTERN_DEFAULT),
                               na.rm = TRUE)
      src <- "photon budget"
    }
  }
  shiny::validate(shiny::need(is.finite(sig_eff) && sig_eff > 0,
    "Localization precision unknown: the file carries no photon counts, so enter σ_loc on this tab."))
  d_in <- df[, intersect(c("tid", "x", "y", "tim", "eco"), names(df)), drop = FALSE]

  # ── Step budget ──────────────────────────────────────────────────────────
  # Fit time scales with the number of steps, and the optimiser evaluates the
  # likelihood over every step many times. The spread of D across molecules, not
  # the step count, is what limits how precisely D is known, so beyond a budget
  # the extra steps buy time rather than certainty. Whole TRACKS are drawn, never
  # parts of them, so no track is cut in half.
  cap <- suppressWarnings(as.integer(input$trk_hmm_max_steps %||% TRK_DEFAULTS$trk_hmm_max_steps))
  n_sub_tracks <- NA_integer_
  if (is.finite(cap) && cap > 0) {
    per_tid <- table(d_in$tid)
    if (sum(per_tid - 1L) > cap) {
      set.seed(1L)                                     # same subset on a re-run
      shuf <- sample(names(per_tid))
      keep <- shuf[cumsum(per_tid[shuf] - 1L) <= cap]
      if (length(keep) >= 2L) {
        d_in <- d_in[d_in$tid %in% keep, , drop = FALSE]
        n_sub_tracks <- length(per_tid)
      }
    }
  }

  sigma_mode <- "one value for all localizations"
  if (identical(src, "photon budget") && "eco" %in% names(d_in)) {
    cve <- tryCatch(trk_cve_bundle(), error = function(e) NULL)
    L_tcp  <- suppressWarnings(as.numeric(input$trk_L_tcp %||% DMAX_L_TCP_DEFAULT))
    if (!is.finite(L_tcp) || L_tcp <= 0) L_tcp <- DMAX_L_TCP_DEFAULT
    nu     <- suppressWarnings(as.numeric(input$trk_beam_fwhm %||% TCP_BEAM_FWHM_NM))
    if (!is.finite(nu) || nu <= 0) nu <- TCP_BEAM_FWHM_NM
    a_pref <- .a_in_use()
    D_off  <- if (!is.null(cve) && is.finite(cve$D_cve) && cve$D_cve > 0) cve$D_cve else 0
    okc    <- is.finite(d_in$eco) & d_in$eco > 0
    d_in$sig <- NA_real_                               # rows without eco fall back to sig_eff
    if (any(okc)) d_in$sig[okc] <- .crlb_doughnut(L_tcp, nu, d_in$eco[okc], D = D_off, t_loc = t_c,
                                                  pref = a_pref,
                                                  pattern = input$trk_tcp_pattern %||% TCP_PATTERN_DEFAULT)
    sigma_mode <- "per localization, from its photon count"
  }

  pr <- hmm_ct_prepare(d_in, ndim = 2L, sigma = sig_eff, min_pts = min_pts)
  shiny::validate(shiny::need(!is.null(pr) && pr$n_steps >= 50L,
    "Too few valid steps - lower 'Min localizations per TID'."))

  fit <- withProgress(message = "Fitting continuous-time HMM…", value = 0.2, {
    f <- hmm_ct_fit(pr, K = K, R_abs = R_abs, t_c = t_c)
    incProgress(0.5)
    f
  })
  shiny::validate(shiny::need(!is.null(fit), "The HMM fit failed - try fewer states."))

  # Neighbouring model sizes, so the summary can show whether the data support K.
  # Two further fits, so it is opt-in: it roughly triples the run time.
  bic_other <- if (!isTRUE(input$trk_hmm_bic)) NULL else
    withProgress(message = "Comparing state counts (BIC)…", value = 0.3, {
    ks <- setdiff(c(K - 1L, K + 1L), 1L); ks <- ks[ks <= 4L]
    stats::setNames(vapply(ks, function(k) {
      f2 <- tryCatch(hmm_ct_fit(pr, K = k, R_abs = R_abs, t_c = t_c, hessian = FALSE,
                                n_starts = 2L), error = function(e) NULL)
      if (is.null(f2)) NA_real_ else f2$BIC
    }, numeric(1)), ks)
  })

  boot <- if (B >= 10L) withProgress(message = sprintf("Bootstrap over tracks (%d refits)…", B), value = 0,
    hmm_ct_bootstrap(pr$df, fit, B = B, ndim = 2L, sigma = sig_eff, t_c = t_c,
                     min_pts = min_pts, progress = function(p) setProgress(p)))
    else NULL

  # Display-only decode under the fitted model.
  step_states <- hmm_ct_viterbi(pr, fit)
  lab <- hmm_ct_states_by_loc(pr, step_states)
  lab <- lab[is.finite(lab$state), , drop = FALSE]
  seg_list <- lapply(split(lab, lab$tid), function(d)
    data.frame(tid = d$tid, x = d$x, y = d$y, tim = d$tim, state = as.integer(d$state),
               stringsAsFactors = FALSE))
  # Longest tracks first, so the example panel shows tracks with room to switch.
  seg_list <- seg_list[order(-vapply(seg_list, nrow, integer(1)))]

  dt_med  <- stats::median(pr$DT, na.rm = TRUE)
  D_noise <- sig_eff^2 / dt_med                        # D equivalent of the noise per step
  # Validated range (tests/validation/hmm_recovery_sweep.R, HMM_RESOLV_* in packages.R).
  exp_sw  <- fit$expected_switches
  resolv <- data.frame(
    state        = seq_len(K),
    dwell_cycles = fit$dwell_ms / t_c,
    too_short    = fit$dwell_ms / t_c < HMM_RESOLV_MIN_DWELL_CYC,
    few_switches = rowSums(exp_sw, na.rm = TRUE) < HMM_RESOLV_MIN_SWITCHES,
    stringsAsFactors = FALSE)
  # Two states with the same D differ only in how long they last: the model is
  # describing a non-exponential dwell distribution of ONE diffusivity, not a
  # separate state. Seen on the tracking test file at K = 3 (D 0.026 / 0.026 / 0.49 µm²/s,
  # dwells 17 / 0.26 / 0.17 ms). Flag adjacent states within 5% of each other.
  same_D <- rep(FALSE, K)
  for (k in seq_len(K - 1L)) if (abs(fit$D[k + 1L] / fit$D[k] - 1) < 0.05) same_D[c(k, k + 1L)] <- TRUE
  resolv$same_D <- same_D

  structure(list(K = K, D_vals = fit$D, fracs = fit$pi,
       P_cycle = .expm_small(fit$Q * t_c), Qrate = fit$Q, dwell_ms = fit$dwell_ms,
       fit = fit, boot = boot, bic_other = bic_other, n_tracks_all = n_sub_tracks,
       seg_list = seg_list, n_tracks = fit$n_tracks, n_steps = fit$n_steps,
       dt_med_ms = dt_med, t_cycle_ms = t_c, R_abs = R_abs,
       sigma_loc_nm = sig_eff, sigma_mode = sigma_mode,
       D_noise = D_noise, noise_state = fit$D <= D_noise, resolv = resolv,
       norm_dt = TRUE),
    # src is the app-wide source, or "photon budget" computed in this panel
    sigma_info = .sigma_stamp(sig_eff, src, note = sigma_mode))
})

output$trk_hmm_params_plot <- renderPlotly({
  bd   <- hmm_bundle()
  dark <- isTRUE(isolate(input$trk_dark_mode))
  bg_col <- if (dark) "#000000" else "#FFFFFF"
  fg_col <- if (dark) "#FFFFFF" else "#333333"
  state_cols <- c("#C03030","#D08020","#2080C0","#7030A0")[seq_len(bd$K)]
  D_um <- bd$D_vals / 1e3                               # µm²/s
  ci   <- if (!is.null(bd$boot)) bd$boot$D else bd$fit$ci$D   # bootstrap CI preferred
  lo   <- if (!is.null(ci)) ci[1, ] / 1e3 else rep(NA_real_, bd$K)
  hi   <- if (!is.null(ci)) ci[2, ] / 1e3 else rep(NA_real_, bd$K)
  is_noise <- bd$noise_state %||% rep(FALSE, bd$K)
  labels <- sprintf("State %d\n%.0f%% of time", seq_len(bd$K), bd$fracs * 100)
  labels <- ifelse(is_noise, paste0(labels, "\n⚠ at noise level"), labels)
  bar_cols <- ifelse(is_noise, "#9e9e9e", state_cols)
  p <- plot_ly() %>%
    add_bars(
      x = D_um, y = labels, orientation = "h",
      marker = list(color = bar_cols),
      error_x = if (all(is.finite(c(lo, hi))))
                  list(type = "data", symmetric = FALSE, array = hi - D_um,
                       arrayminus = D_um - lo, color = fg_col) else NULL,
      text = sprintf("D = %.3g µm²/s  (95%% CI %.3g to %.3g)", D_um, lo, hi),
      hovertemplate = "%{text}<extra></extra>")
  p %>% layout(
    xaxis = list(title = "D (µm²/s)", color = fg_col),
    yaxis = list(title = "", autorange = "reversed", color = fg_col),
    paper_bgcolor = bg_col, plot_bgcolor = bg_col,
    font = list(color = fg_col), showlegend = FALSE, uirevision = "hmm_params")
})

output$trk_hmm_trans_plot <- renderPlotly({
  bd   <- hmm_bundle()
  dark <- isTRUE(isolate(input$trk_dark_mode))
  bg_col <- if (dark) "#000000" else "#FFFFFF"
  fg_col <- if (dark) "#FFFFFF" else "#333333"
  K      <- bd$K
  labels <- paste0("S", seq_len(K))
  P      <- bd$P_cycle                                  # expm(Q·t_c): one cycle
  plot_ly(
    x = labels, y = rev(labels),
    z = P[rev(seq_len(K)), , drop = FALSE],
    type = "heatmap",
    colorscale = list(list(0, "#08306b"), list(0.5, "#4393c3"), list(1, "#ffffff")),
    zmin = 0, zmax = 1,
    text  = matrix(sprintf("%.3f", P[rev(seq_len(K)), ]), K, K),
    texttemplate = "%{text}",
    hovertemplate = "From %{y} → To %{x} within one cycle: %{text}<extra></extra>"
  ) %>% layout(
    xaxis = list(title = "To state",   color = fg_col, side = "top"),
    yaxis = list(title = "From state", color = fg_col),
    paper_bgcolor = bg_col, plot_bgcolor = bg_col,
    font = list(color = fg_col), uirevision = "hmm_trans")
})

output$trk_hmm_example_plot <- renderPlotly({
  bd     <- hmm_bundle()
  n_ex   <- max(1L, as.integer(isolate(input$trk_hmm_n_ex)))
  dark   <- isTRUE(isolate(input$trk_dark_mode))
  bg_col <- if (dark) "#000000" else "#FFFFFF"
  fg_col <- if (dark) "#FFFFFF" else "#333333"
  state_cols <- c("#C03030","#D08020","#2080C0","#7030A0")[seq_len(bd$K)]

  ex <- head(bd$seg_list, n_ex)
  # Lay out tracks in a grid - compute a consistent offset spacing from
  # the 90th percentile of Rmax across example tracks so the grid never
  # overlaps but also isn't wastefully spread out.
  spans <- sapply(ex, function(seg) {
    cx <- mean(seg$x); cy <- mean(seg$y)
    r  <- sqrt((seg$x - cx)^2 + (seg$y - cy)^2)
    max(r, na.rm = TRUE)
  })
  spacing <- max(50, quantile(spans, 0.90, na.rm = TRUE) * 2.5)
  n_cols  <- max(1L, ceiling(sqrt(length(ex))))

  fig <- plot_ly()
  is_first <- TRUE
  for (i in seq_along(ex)) {
    seg   <- ex[[i]]
    col_i <- ((i - 1L) %% n_cols)
    row_i <- ((i - 1L) %/% n_cols)
    cx <- mean(seg$x); cy <- mean(seg$y)
    xc <- seg$x - cx + col_i * spacing
    yc <- seg$y - cy - row_i * spacing   # negative so row 0 is at top
    # Grey track line
    fig <- fig %>% add_lines(
      x = xc, y = yc,
      line = list(color = "rgba(150,150,150,0.3)", width = 1),
      hoverinfo = "skip", showlegend = FALSE
    )
    # Points coloured by state
    for (k in seq_len(bd$K)) {
      idx <- seg$state == k
      if (!any(idx)) next
      fig <- fig %>% add_markers(
        x = xc[idx], y = yc[idx],
        marker = list(size = 4, color = state_cols[k], opacity = 0.8),
        name = paste0("State ", k),
        legendgroup = paste0("state", k),
        showlegend = is_first,
        hovertemplate = paste0("TID: ", seg$tid[1L], "  State ", k, "<extra></extra>")
      )
    }
    is_first <- FALSE
  }
  fig %>% layout(
    xaxis = list(title = "position (nm, centred + offset)", color = fg_col,
                 zeroline = FALSE, showticklabels = FALSE),
    yaxis = list(title = "", color = fg_col, zeroline = FALSE,
                 showticklabels = FALSE, scaleanchor = "x", scaleratio = 1),
    paper_bgcolor = bg_col, plot_bgcolor = bg_col,
    font = list(color = fg_col),
    uirevision = "hmm_ex"
  )
})

output$trk_hmm_summary <- renderUI({
  bd  <- hmm_bundle()
  K   <- bd$K; fit <- bd$fit
  fmt_time <- function(ms) if (!is.finite(ms)) "-"
                           else if (ms >= 1000) sprintf("%.2f s", ms / 1000)
                           else if (ms >= 1) sprintf("%.1f ms", ms)
                           else sprintf("%.2f ms", ms)
  ci_src <- if (!is.null(bd$boot)) bd$boot else fit$ci   # bootstrap preferred
  ci_lbl <- if (!is.null(bd$boot)) sprintf("bootstrap over tracks, %d refits", bd$boot$n_ok)
            else "from the likelihood curvature; too narrow when molecules differ"
  rng <- function(m, k, scale = 1, f = function(v) sprintf("%.3g", v))
    if (is.null(m)) "" else sprintf(" (%s to %s)", f(m[1, k] * scale), f(m[2, k] * scale))
  rs  <- bd$resolv
  rows <- lapply(seq_len(K), function(k) {
    flag <- c(if (isTRUE(bd$noise_state[k])) "D at the noise level",
              if (isTRUE(rs$too_short[k])) "dwell too short to resolve reliably",
              if (isTRUE(rs$few_switches[k])) "too few switches to constrain the rate",
              if (isTRUE(rs$same_D[k])) "same D as a neighbouring state (one diffusivity with non-exponential dwells)")
    tags$tr(
      tags$td(sprintf("State %d", k), style = "padding:2px 8px;"),
      tags$td(sprintf("%.0f%%", bd$fracs[k] * 100), style = "padding:2px 8px;"),
      tags$td(sprintf("%.3g", bd$D_vals[k] / 1e3), rng(ci_src$D, k, 1e-3), style = "padding:2px 8px;"),
      tags$td(fmt_time(bd$dwell_ms[k]), rng(ci_src$dwell_ms, k, 1, fmt_time),
              sprintf(" ≈ %.0f cycles", rs$dwell_cycles[k]), style = "padding:2px 8px;"),
      tags$td(if (length(flag)) tags$span(paste("⚠", paste(flag, collapse = "; ")),
                                          style = "color:#b26a00;") else "",
              style = "padding:2px 8px;"))
  })
  rate_txt <- {
    off <- character(0)
    for (a in seq_len(K)) for (b in seq_len(K)) if (a != b)
      off <- c(off, sprintf("S%d→S%d %.3g s⁻¹%s", a, b, bd$Qrate[a, b] * 1000,
                            if (!is.null(ci_src$rates_lo))
                              sprintf(" (%.3g to %.3g)", ci_src$rates_lo[a, b], ci_src$rates_hi[a, b])
                            else ""))
    paste(off, collapse = "  ·  ")
  }
  bic_txt <- {
    bo <- bd$bic_other
    if (length(bo) && any(is.finite(bo))) {
      all_bic <- c(stats::setNames(fit$BIC, K), bo)
      best <- as.integer(names(all_bic)[which.min(all_bic)])
      sprintf("BIC prefers %d states (%s).", best,
              paste(sprintf("K=%s: %+.0f", names(all_bic), all_bic - min(all_bic, na.rm = TRUE)),
                    collapse = ", "))
    } else ""
  }
  # JDD cross-check: at short lags the JDD component weights estimate the same
  # occupancies (per step rather than per unit time), with no segmentation at all.
  jdd_txt <- {
    jb <- tryCatch(jdd_bundle(), error = function(e) NULL)
    if (!is.null(jb) && isTRUE(jb$n_comp == K) && length(jb$ws) == K) {
      o <- order(jb$D_vals)
      sprintf("JDD (%d components): D %s µm²/s, weights %s.", K,
              paste(sprintf("%.3g", jb$D_vals[o] / 1e3), collapse = " / "),
              paste(sprintf("%.0f%%", 100 * jb$ws[o]), collapse = " / "))
    } else NULL
  }
  tags$div(style = "font-size:0.88em; margin:4px 0;",
    tags$p(style = "color:#888; margin:2px 0;",
      sprintf("n = %d tracks%s, %d steps  |  cycle %.3g ms  |  σ_loc %s (median %.2f nm)  |  log-likelihood %.1f",
              bd$n_tracks,
              if (is.finite(bd$n_tracks_all %||% NA)) sprintf(" of %d, drawn at random to meet the step budget",
                                                              bd$n_tracks_all) else "",
              bd$n_steps, bd$t_cycle_ms, bd$sigma_mode, bd$sigma_loc_nm, fit$loglik)),
    tags$table(style = "border-collapse:collapse;",
      tags$thead(tags$tr(lapply(
        c("", "Occupancy", "D (µm²/s, 95% CI)", "Mean dwell (95% CI)", ""),
        function(h) tags$th(h, style = "padding:2px 8px; border-bottom:1px solid #ccc;")))),
      tags$tbody(rows)),
    tags$p(style = "color:#444; margin:6px 0 2px;", tags$b("Switching rates: "), rate_txt),
    tags$p(style = "color:#888; font-size:0.9em; margin:2px 0;",
      sprintf("CIs %s. %s", ci_lbl, bic_txt)),
    if (!is.null(jdd_txt))
      tags$p(style = "color:#888; font-size:0.9em; margin:2px 0;", tags$b("Cross-check. "), jdd_txt),
    if (!isTRUE(fit$converged))
      tags$p(style = "color:#b26a00; margin:2px 0;", "The optimiser did not report convergence; treat the estimates with caution."))
})

output$trk_dl_hmm_csv <- downloadHandler(
  filename = function() paste0("hmm_states_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
  content  = function(file) {
    bd  <- hmm_bundle()
    out <- do.call(rbind, bd$seg_list)
    write.csv(out, file, row.names = FALSE)
  }
)

# ==========================================================================
# seg_msd_rv() - PER-STATE MSD (one MSD curve per HMM diffusive state)
# ==========================================================================
#
# GOAL
#   Fix the biggest weakness of a whole-track MSD: AVERAGING OVER PHASES.
#   If a molecule spends part of its track confined at a wall and part of it
#   being transported, one MSD curve mixes both and describes neither. The
#   fitted "D" is then an artefact of the mixing ratio, not a property of the
#   molecule.
#
#   This splits every track by the HMM's state labels FIRST, then builds a
#   separate MSD for each state. Now the physics-specific models finally mean
#   something:
#       Confined model on the slow/wall state  -> a real confinement size L
#       Directed model on the fast state       -> a real velocity v
#
# THE MATH
#   Identical to the main MSD engine - same cumulative cycle-index binning,
#   same fit window, same model set - but pairs are pooled PER STATE instead of
#   per dataset. Only pairs where BOTH localizations carry the same state label
#   contribute, so a curve never mixes phases.
#
#   Being an ensemble over many short same-state runs (rather than a few long
#   tracks) is actually an advantage here: it is the vbSPT-style trick of
#   pooling statistical power across many short segments.
#
# OBJECTS IT CALLS (inputs)
#   hmm_bundle()$seg_list ... per-localization state labels - THE dependency;
#                             this analysis cannot run until the HMM has run
#   trk_cycle_time(), trk_cluster_cycle()$sw
#   .msd_model_fit() ........ the shared fitter, so the model/window/linear-τ/
#                             iterative settings are INHERITED from the main
#                             MSD panel rather than duplicated here
#   inputs: trk_seg_msd_run (Run trigger), trk_seg_msd_min_pts,
#           trk_msd_ndim, trk_uni_max_lag (+ the inherited fit settings)
#
# STEP BY STEP
#   1. Take the HMM's Viterbi labels; states are already ordered slow -> fast.
#   2. Within each track, run-length-split into CONTIGUOUS same-state runs and
#      discard runs shorter than `min_pts` (too few pairs to be meaningful).
#   3. For each run, accumulate displacement pairs into that STATE's bins,
#      using the same cumulative cycle index as the main MSD.
#   4. Per state: ensemble MSD = Σdr² / Σn per lag.
#   5. Fit each state's curve with the inherited model and window.
#   6. Cross-check: the fitted D per state should broadly agree with the HMM's
#      own D_k = σ_k²/2. Large disagreement means the state assignment is
#      unreliable - treat both with suspicion.
#
# RETURNS a list with, per state: the MSD curve (lag, tau, msd, npairs), the
# fit result and its model-specific extras (L for confined, v for directed,
# α for anomalous), segment and point counts.
#
# CAVEAT: everything here is downstream of the HMM. A poor state decode makes
# these curves meaningless no matter how good the fit looks. Confined/anomalous
# MSD on MINFLUX is also flagged exploratory (see the MSD model note).
# ==========================================================================
SEG_MSD_PAL <- c("#3060C0","#C03030","#2CA02C","#FF7F0E","#9467BD","#8C564B")

seg_msd_rv <- eventReactive(input$trk_seg_msd_run, {
  req(isTRUE(input$trk_show_seg_msd))
  bd <- tryCatch(hmm_bundle(), error = function(e) NULL)
  shiny::validate(shiny::need(
    !is.null(bd) && !is.null(bd$seg_list) && length(bd$seg_list) > 0,
    "Enable the Diffusion State HMM (above) first - its Viterbi states define the segments."))
  K   <- bd$K
  t_c <- trk_cycle_time(); req(is.finite(t_c) && t_c > 0)
  cc  <- tryCatch(trk_cluster_cycle(), error = function(e) NULL)
  sw  <- if (!is.null(cc) && is.finite(cc$sw)) cc$sw else 0
  min_seg <- max(5L, as.integer(input$trk_seg_msd_min_pts %||% TRK_DEFAULTS$trk_seg_msd_min_pts))
  has_z <- any(vapply(bd$seg_list, function(d) "z" %in% names(d) && any(is.finite(d$z)), logical(1)))
  ndim_sel <- input$trk_msd_ndim %||% TRK_DEFAULTS$trk_msd_ndim
  ndim <- if (identical(ndim_sel, "3")) 3L else if (identical(ndim_sel, "2")) 2L else if (has_z) 3L else 2L
  max_lag <- max(60L, min(5000L, as.integer(input$trk_uni_max_lag %||% TRK_DEFAULTS$trk_uni_max_lag)))

  dr2 <- matrix(0,  max_lag, K); np <- matrix(0L, max_lag, K); dts <- matrix(0, max_lag, K)
  n_seg <- integer(K); n_pts <- integer(K)

  withProgress(message = "Per-segment MSD…", value = 0, {
    for (di in seq_along(bd$seg_list)) {
      incProgress(1 / length(bd$seg_list))
      d <- bd$seg_list[[di]]; d <- d[order(d$tim), ]
      st <- d$state; n <- nrow(d); if (n < 2L) next
      runs <- rle(st); ends <- cumsum(runs$lengths); starts <- ends - runs$lengths + 1L
      for (r in seq_along(runs$lengths)) {
        sstate <- runs$values[r]
        if (is.na(sstate) || sstate < 1L || sstate > K) next
        idx <- starts[r]:ends[r]; if (length(idx) < min_seg) next
        tb <- d[idx, , drop = FALSE]
        n_seg[sstate] <- n_seg[sstate] + 1L; n_pts[sstate] <- n_pts[sstate] + length(idx)
        dcyc <- pmax(1L, as.integer(round((diff(tb$tim) - sw) / t_c)))
        cyc  <- c(0L, cumsum(dcyc)); nn <- nrow(tb)
        for (i in seq_len(nn - 1L)) {
          dt_v <- tb$tim[seq.int(i + 1L, nn)] - tb$tim[i]
          lag  <- cyc[seq.int(i + 1L, nn)] - cyc[i]
          okk  <- lag >= 1L & lag <= max_lag
          if (!any(okk)) next
          js <- seq.int(i + 1L, nn)[okk]; b <- lag[okk]; dtk <- dt_v[okk]
          r2 <- (tb$x[js] - tb$x[i])^2 + (tb$y[js] - tb$y[i])^2
          if (has_z && "z" %in% names(tb)) r2 <- r2 + (tb$z[js] - tb$z[i])^2
          rs <- rowsum(r2, b, reorder = FALSE); ix <- as.integer(rownames(rs))
          dr2[ix, sstate] <- dr2[ix, sstate] + rs[, 1L]
          rt <- rowsum(dtk, b, reorder = FALSE); it <- as.integer(rownames(rt))
          dts[it, sstate] <- dts[it, sstate] + rt[, 1L]
          np[, sstate] <- np[, sstate] + tabulate(b, nbins = max_lag)
        }
      }
    }
  })

  model    <- input$trk_uni_model %||% TRK_DEFAULTS$trk_uni_model
  min_lag  <- max(1L, as.integer(input$trk_uni_min_lag %||% TRK_DEFAULTS$trk_uni_min_lag))
  max_fit  <- max(min_lag + 2L, as.integer(input$trk_uni_max_lag %||% TRK_DEFAULTS$trk_uni_max_lag))
  use_iter <- isTRUE(input$trk_uni_iterative)
  linear_tau <- isTRUE(input$trk_uni_linear_tau)
  lag_v <- seq_len(max_lag)
  states <- Filter(Negate(is.null), lapply(seq_len(K), function(k) {
    npk <- np[, k]; if (sum(npk) == 0L) return(NULL)
    msd <- ifelse(npk > 0L, dr2[, k] / npk, NA_real_)
    tau_use <- if (linear_tau) lag_v * t_c else ifelse(npk > 0L, dts[, k] / npk, lag_v * t_c)
    okf <- is.finite(tau_use) & is.finite(msd) & lag_v >= min_lag & lag_v <= max_fit & npk > 0L
    fr  <- if (sum(okf) >= 3L)
             # Same motion/projection correction as the main MSD panel, so a
             # per-state D is on the same footing as the ensemble D.
             tryCatch(.msd_model_fit(tau_use[okf], msd[okf], npk[okf], ndim, model, use_iter, t_c,
                                     coef_override = 2 * .msd_d_eff(ndim), R = trk_blur()$R_msd),
                      error = function(e) NULL) else NULL
    list(state = k, n_seg = n_seg[k], n_pts = n_pts[k], D_hmm = bd$D_vals[k] / 1e3,
         lag = lag_v, tau = tau_use, msd = msd, npairs = npk, fr = fr)
  }))
  shiny::validate(shiny::need(length(states) > 0, "No segments long enough - lower 'Min points per segment' or the HMM min points."))
  list(states = states, K = K, ndim = ndim, t_cycle = t_c, model = model, linear_tau = linear_tau)
})

output$trk_seg_msd_plot <- renderPlotly({
  req(isTRUE(input$trk_show_seg_msd))
  sm <- seg_msd_rv()
  dark <- isTRUE(isolate(input$trk_dark_mode))
  bg <- if (dark) "#000000" else "#FFFFFF"; fg <- if (dark) "#FFFFFF" else "#333333"
  fig <- plot_ly()
  for (st in sm$states) {
    ok <- is.finite(st$tau) & is.finite(st$msd) & st$npairs > 0L
    if (sum(ok) < 2L) next
    col <- SEG_MSD_PAL[((st$state - 1L) %% length(SEG_MSD_PAL)) + 1L]
    fig <- fig %>% add_trace(x = st$tau[ok], y = st$msd[ok] / 1e6, type = "scatter",
      mode = "lines+markers", marker = list(size = 4, color = col),
      line = list(color = col, width = 2), name = sprintf("State %d eMSD (%d seg)", st$state, st$n_seg))
    fr <- st$fr
    if (!is.null(fr)) {
      lbl <- switch(fr$model %||% "linear",
        confined  = if (isTRUE(fr$degenerate)) sprintf("State %d: ≈linear (no plateau)", st$state)
                    else sprintf("State %d confined: D=%.3g, L_c=%.0f nm", st$state, fr$D / 1e3, fr$L_conf),
        directed  = sprintf("State %d directed: D=%.3g, v=%.3g µm/s", st$state, fr$D / 1e3, fr$v_um_s),
        anomalous = sprintf("State %d: D=%.3g, α=%.2f", st$state, fr$D / 1e3, fr$alpha),
        sprintf("State %d linear: D=%.3g µm²/s", st$state, fr$D / 1e3))
      fig <- fig %>% add_trace(x = fr$tau_fit, y = pmax(0, fr$msd_fit) / 1e6, type = "scatter",
        mode = "lines", line = list(color = col, width = 2.5, dash = "dash"), name = lbl)
    }
  }
  fig %>% layout(
    xaxis = list(title = if (sm$linear_tau) "τ (ms, lag × cycle)" else "τ (ms, mean Δt)", color = fg),
    yaxis = list(title = "MSD (µm²)", color = fg, rangemode = "tozero"),
    paper_bgcolor = bg, plot_bgcolor = bg, font = list(color = fg),
    legend = list(font = list(color = fg, size = 10)),
    margin = list(l = 55, r = 10, t = 10, b = 45), uirevision = "segmsd")
})

output$trk_seg_msd_summary <- renderUI({
  req(isTRUE(input$trk_show_seg_msd))
  sm <- tryCatch(seg_msd_rv(), error = function(e) NULL); if (is.null(sm)) return(NULL)
  rows <- lapply(sm$states, function(st) {
    fr <- st$fr
    param <- if (is.null(fr)) "-" else switch(fr$model %||% "linear",
      confined  = if (isTRUE(fr$degenerate)) "no plateau (≈linear)" else sprintf("L_c = %.0f nm", fr$L_conf),
      directed  = sprintf("v = %.3g µm/s", fr$v_um_s),
      anomalous = sprintf("α = %.2f", fr$alpha), "-")
    D_msd <- if (!is.null(fr) && is.finite(fr$D)) sprintf("%.3g", fr$D / 1e3) else "-"
    col <- SEG_MSD_PAL[((st$state - 1L) %% length(SEG_MSD_PAL)) + 1L]
    htmltools::tags$tr(
      htmltools::tags$td(htmltools::tags$span(style = sprintf("color:%s;font-weight:bold;", col), sprintf("State %d", st$state))),
      htmltools::tags$td(st$n_seg), htmltools::tags$td(format(st$n_pts, big.mark = ",")),
      htmltools::tags$td(D_msd), htmltools::tags$td(param),
      htmltools::tags$td(sprintf("%.3g", st$D_hmm)))
  })
  htmltools::div(style = "font-size:0.85em; margin-top:6px;",
    htmltools::tags$b(sprintf("Per-state MSD (%s model):", sm$model)),
    htmltools::tags$table(class = "table table-condensed", style = "width:auto; margin-top:4px;",
      htmltools::tags$thead(htmltools::tags$tr(lapply(
        c("State", "segments", "points", "D (MSD, µm²/s)", "model param", "D (HMM, µm²/s)"),
        htmltools::tags$th))),
      htmltools::tags$tbody(rows)),
    htmltools::tags$div(style = "color:#888; font-size:0.92em;",
      "States are the HMM diffusion states (sorted slow→fast). Switch the MSD ",
      "'Fit model' to Confined (for a slow/wall state → size L) or Directed (fast/exported → velocity v)."))
})

output$trk_dl_seg_msd_csv <- downloadHandler(
  filename = function() paste0("segment_msd_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
  content  = function(file) {
    sm <- seg_msd_rv()
    out <- do.call(rbind, lapply(sm$states, function(st) {
      ok <- st$npairs > 0L
      data.frame(state = st$state, cycle_lag = st$lag[ok], tau_ms = st$tau[ok],
                 MSD_um2 = st$msd[ok] / 1e6, n_pairs = st$npairs[ok])
    }))
    write.csv(out, file, row.names = FALSE)
  }
)
