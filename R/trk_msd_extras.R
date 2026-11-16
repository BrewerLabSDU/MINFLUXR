## ─────────────────────────────────────────────────────────────────────────────
## trk_msd_extras.R  -  tracking: MSD fit window, nls model fit and D_max
## ─────────────────────────────────────────────────────────────────────────────
##
## PURPOSE
##   The fit-window mask shared by the MSD views, the nls MSD fit used by the
##   MSD cluster models, and the trackable-diffusion ceiling D_max.
##
##   Plain functions and constants with no Shiny state, moved unchanged out of
##   the tracking_server() closure (server_tracking.R). The server gathers the
##   inputs and calls these.
##
## UNITS: nm and ms; D in nm^2/ms unless a name says um2_s.
## ─────────────────────────────────────────────────────────────────────────────

# MSD fit core (blur R from trk_blur(), .olsf_sig2, .michalet_p, msd_fit_olsf,
# msd_fit_olsf_iterative) lives in est_msd.R since 2026-09-17.

msd_fit_result <- function(tau, msd, k, anomalous = TRUE, ndim = 2L) {
  df <- data.frame(tau = tau[seq_len(k)], msd = msd[seq_len(k)])
  df <- df[is.finite(df$tau) & is.finite(df$msd) & df$tau > 0, ]
  if (nrow(df) < 3L) return(NULL)
  coef_val <- 2L * ndim   # 4 for 2D, 6 for 3D

  D0     <- max(1e-6, (max(df$msd) - min(df$msd)) / (coef_val * max(df$tau)))
  sig0   <- max(0, sqrt(abs(df$msd[1]) / coef_val))
  alpha0 <- 1.0
  alpha_se <- NA_real_   # parametric SE on alpha from the ensemble fit (set below)

  fm_anom <- as.formula(sprintf("msd ~ %d*D*tau^alpha + %d*sig2", coef_val, coef_val))
  fm_lin  <- as.formula(sprintf("msd ~ %d*D*tau + %d*sig2",       coef_val, coef_val))

  if (anomalous) {
    fit <- tryCatch(
      nls(fm_anom, data = df,
          start   = list(D = D0, alpha = alpha0, sig2 = sig0^2),
          lower   = list(D = 0,  alpha = 0.05,   sig2 = 0),
          upper   = list(D = Inf, alpha = 2.5,   sig2 = Inf),
          algorithm = "port",
          control = nls.control(maxiter = 200, warnOnly = TRUE)),
      error = function(e) NULL
    )
    if (is.null(fit)) {
      fit <- tryCatch(
        nls(fm_lin, data = df,
            start = list(D = D0, sig2 = sig0^2),
            lower = list(D = 0, sig2 = 0),
            algorithm = "port",
            control = nls.control(maxiter = 200, warnOnly = TRUE)),
        error = function(e) NULL
      )
      if (is.null(fit)) return(NULL)
      cf <- coef(fit)
      D_est <- unname(cf["D"]); alpha_est <- 1.0
      sig_est <- sqrt(max(0, unname(cf["sig2"]))); forced_linear <- TRUE
    } else {
      cf <- coef(fit)
      D_est <- unname(cf["D"]); alpha_est <- unname(cf["alpha"])
      sig_est <- sqrt(max(0, unname(cf["sig2"]))); forced_linear <- FALSE
      # Parametric SE on alpha from the fit covariance. NB this treats the averaged
      # MSD points as independent, so it is a LOWER bound on the true uncertainty
      # (MSD lags are correlated) - surfaced in the UI with that caveat.
      alpha_se <- tryCatch(unname(sqrt(diag(vcov(fit)))["alpha"]), error = function(e) NA_real_)
    }
  } else {
    fit <- tryCatch(
      nls(fm_lin, data = df,
          start = list(D = D0, sig2 = sig0^2),
          lower = list(D = 0, sig2 = 0),
          algorithm = "port",
          control = nls.control(maxiter = 200, warnOnly = TRUE)),
      error = function(e) NULL
    )
    if (is.null(fit)) return(NULL)
    cf <- coef(fit)
    D_est <- unname(cf["D"]); alpha_est <- 1.0
    sig_est <- sqrt(max(0, unname(cf["sig2"]))); forced_linear <- FALSE
  }

  tau_seq <- seq(min(df$tau), max(df$tau), length.out = 200)
  yfit    <- pmax(0, coef_val * D_est * tau_seq^alpha_est + coef_val * sig_est^2)

  list(D = D_est, alpha = alpha_est, alpha_se = alpha_se, sigma_loc = sig_est,
       tau_fit = tau_seq, msd_fit = yfit,
       forced_linear = if (anomalous) forced_linear else FALSE,
       ndim = ndim)
}

# Shared fit-window selector - applied consistently by the ensemble estimates,
# the plotted fit line, and the per-track fits so they always agree.
# Lower bound: min_lag_ms.  Upper bound: fit_max_ms (absolute) when > 0, else
# fit_pct % of the remaining points.  Returns a logical mask over `tau`.
msd_fit_mask <- function(tau, min_lag_ms = 0, fit_max_ms = 0, fit_pct = 30) {
  keep <- is.finite(tau)
  if (min_lag_ms > 0) keep <- keep & tau >= min_lag_ms
  idx  <- which(keep)
  if (!length(idx)) return(keep)
  if (fit_max_ms > 0) {
    keep <- keep & tau <= fit_max_ms
  } else {
    k <- max(3L, floor((fit_pct / 100) * length(idx)))
    cutoff <- sort(tau[idx])[min(k, length(idx))]
    keep <- keep & tau <= cutoff
  }
  keep
}

# Shared trackable-diffusion ceiling (D_max) estimate, backing the D_max panel.
# Above ~D_max the MINFLUX localization loop can no longer keep the doughnut centred,
# so D is systematically UNDER-estimated (Vogler et al. 2025; SimuFLUX, Marin & Ries 2025).
#
# WHICH L THIS IS. L_TCP is the diameter of the FINAL iteration - the one that is
# reiterated continuously while tracking - because that is the pattern the particle
# has to stay inside between consecutive fixes. The coarse search iterations that
# precede lock-on are irrelevant here, and so is the much tighter final iteration an
# IMAGING sequence would use: SPT deliberately stops at a larger L, since precision
# improves only linearly as L shrinks (sigma_CRLB ~ L/(4*sqrt(N))) while the escape
# risk grows quadratically (D_max ~ L^2). 150 nm is the tracking-iteration diameter
# in Vogler et al. 2025 and matches the simulator's L default (sim_server.R).
#
# Pure function of the localization table + L_TCP (nm) & hardware dead-time t_hw (ms).
# Returns D_max in um^2/s (Vogler form preferred, simple form fallback).
dmax_estimate <- function(df, L = DMAX_L_TCP_DEFAULT, t_hw = 0.05) {
  if (is.null(df) || !nrow(df) || !all(c("tid","tim") %in% names(df))) return(NULL)
  if (!is.finite(L)    || L    <= 0) L    <- DMAX_L_TCP_DEFAULT
  if (!is.finite(t_hw) || t_hw <  0) t_hw <- 0.05
  efo_med <- if ("efo" %in% names(df)) { v <- df$efo[is.finite(df$efo) & df$efo > 0]; if (length(v)) median(v) else NA_real_ } else NA_real_
  eco_med <- if ("eco" %in% names(df)) { v <- df$eco[is.finite(df$eco) & df$eco > 0]; if (length(v)) median(v) else NA_real_ } else NA_real_
  n_pl <- eco_med   # ECO is the photon count per localization (= N_PL directly)
  t_loc_eco <- if (is.finite(eco_med) && is.finite(efo_med) && efo_med > 0)
    (eco_med / efo_med) * 1000 else NA_real_   # ms to collect ECO photons at rate EFO
  # Median within-track dt fallback (fast integer-coded, not slow character tapply).
  tf  <- match(df$tid, unique(df$tid))
  ord <- order(tf, df$tim, method = "radix")
  tf2 <- tf[ord]; tim2 <- df$tim[ord]; n <- length(tf2)
  same <- c(FALSE, tf2[-1] == tf2[-n])
  dtd  <- c(NA_real_, diff(tim2))[same]; dtd <- dtd[is.finite(dtd) & dtd > 0]
  t_loc_dt <- if (length(dtd)) median(dtd) else NA_real_
  t_loc <- if (is.finite(t_loc_eco)) t_loc_eco else t_loc_dt
  D_max_simple <- if (is.finite(t_loc) && t_loc > 0) (L^2 / (16 * t_loc)) / 1e3 else NA_real_
  D_max_vogler <- if (is.finite(efo_med) && is.finite(n_pl) && is.finite(t_hw)) {
    psi <- efo_med; t_hw_s <- t_hw / 1000
    ((L^2 / 16) * psi / (n_pl + psi * t_hw_s)) / 1e6
  } else NA_real_
  list(efo_med = efo_med, n_pl = n_pl, t_loc_eco = t_loc_eco, t_loc_dt = t_loc_dt,
       t_loc = t_loc, D_max_simple = D_max_simple, D_max_vogler = D_max_vogler,
       D_max = if (is.finite(D_max_vogler)) D_max_vogler else D_max_simple)
}

