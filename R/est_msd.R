## ─────────────────────────────────────────────────────────────────────────────
## est_msd.R  -  MSD model fitting (OLSF linear fit and the non-linear models)
## ─────────────────────────────────────────────────────────────────────────────
##
## PURPOSE
##   The fitting core behind the Tracking tab's MSD panels, as pure functions so
##   tests/ can check them against simulated ground truth. Moved verbatim out of
##   tracking_server() on 2026-09-17 (publication audit C1); the server still
##   builds the MSD curves (msd_unified_bundle) and calls these to fit them.
##
## WHAT LIVES HERE
##   (blur R: est_blur.R; callers pass the sequence-derived value) ... shared by the MSD fits, the
##                              confinement null and the HMM
##   .olsf_sig2() ............. intercept -> σ² with the blur offset added back
##   .michalet_p() ............ Michalet & Berglund optimal point count
##   msd_fit_olsf() ........... anchored OLSF (app default)
##   msd_fit_olsf_iterative() . MinSpt's exact iterative OLSF
##   .msd_model_fit() ......... linear / confined / directed / anomalous dispatcher
##
## UNITS: τ in ms, MSD in nm², D in nm²/ms (÷1000 → µm²/s), σ in nm.
##
## ATTRIBUTION: the OLSF workflow is an R port of MinSpt (Vogler et al. 2025), which is
## GPL-3.0; that is why MINFLUXR as a whole is GPL-3.0. See THIRD-PARTY.md section 1.
## ─────────────────────────────────────────────────────────────────────────────

# --- Shared MSD fit helpers --------------------------------------------------
# ndim = 2 (2D): MSD = 4*D*tau^alpha + 4*sigma^2
# ndim = 3 (3D): MSD = 6*D*tau^alpha + 6*sigma^2
# Units throughout: tau in ms, MSD in nm², D in nm²/ms.
# Display conversion: D [µm²/s] = D [nm²/ms] / 1e3   (1 nm²/ms = 0.001 µm²/s)
#                     MSD [µm²] = MSD [nm²]   / 1e6

# OLSF - faithful port of Vogler's OptimizedLeastSquaresFitting (Michalet &
# Berglund). UNWEIGHTED linear fit MSD = a + b·τ with SEPARATE optimal point
# counts for the slope (→ D) and the intercept (→ σ²), iterated to a fixed
# point.  Edge cases and formulas match his olsf_base_math / olsf_fit_utilities:
#   x   = σ²/(D·dt) − 2R                 (reduced localization error)
#   a   = 2d·σ² − 4d·R·D·dt   (a<0 ⟺ x<0 → p = 2)
#   b   = 2d·D·dt            (b<0 → p = L)
#   f_b = 2 + 1.35·x^0.6 ,  L_b = 0.8 + 0.564·N            (slope)
#   f_a = 2 + 1.6·x^0.51 ,  L_a = 3 + (4.5·N^0.4 − 8.5)^1.2 (intercept)
#   p   = f·L / (f³ + L³)^(1/3)
.michalet_p <- function(D, sig2, dt, N, R, d, which) {
  if (!is.finite(D) || D <= 0 || !is.finite(dt) || dt <= 0) return(NA_real_)
  a <- 2*d*sig2 - 4*d*R*D*dt
  b <- 2*d*D*dt
  if (which == "slope") {
    if (a < 0) return(2); if (b < 0) return(0.8 + 0.564*N)
    x <- sig2/(D*dt) - 2*R; if (x < 0) x <- 0
    fb <- 2 + 1.35*x^0.6; Lb <- 0.8 + 0.564*N
    fb*Lb/((fb^3 + Lb^3)^(1/3))
  } else {
    if (a < 0) return(2); if (b < 0) return(3 + (4.5*N^0.4 - 8.5)^1.2)
    x <- sig2/(D*dt) - 2*R; if (x < 0) x <- 0
    fa <- 2 + 1.6*x^0.51; La <- 3 + (4.5*N^0.4 - 8.5)^1.2
    fa*La/((fa^3 + La^3)^(1/3))
  }
}

# Motion-blur coefficient R for the MSD fits. A MINFLUX localization averages the
# molecule's path over the cycle, which LOWERS every lag's MSD by the same offset
# (Michalet & Berglund 2012):
#     MSD(τ) = coef·D·τ + coef·σ² − 2·coef·R·D·dT        (dT = cycle time)
# so the intercept is coef·(σ² − 2·R·D·dT), not coef·σ². Reading σ² straight off
# the intercept (R = 0) under-reports σ_loc: at D = 0.215 µm²/s, t_c = 0.192 ms and
# σ = 5 nm it returns 3.3 nm (simulated, full-cycle window). D is unaffected.
# R here is relative to the cycle time dT. It comes from the sequence timing
# (blur_R_step() / blur_offsets() in est_blur.R, R_msd = R_step * dt / t_c); the
# server passes it at every call. MinSpt passes a constant 1/6.2, which on MINFLUX
# tracking timings is ~50% above what the photon-collection timeline gives (see
# est_blur.R); the default below is the fallback for the standard 2D sequence.
MSD_R_FALLBACK <- 0.106

# Convert a linear fit's intercept to σ² by adding back the blur offset.
# GOAL: one place for the intercept → σ² relation so every OLSF variant agrees.
# intercept_over_coef = cf[1]/coef; D = slope/coef (nm²/ms); dT = cycle time (ms).
.olsf_sig2 <- function(intercept_over_coef, D, R, dT) {
  blur <- if (is.finite(R) && R > 0 && is.finite(D) && D > 0 && is.finite(dT) && dT > 0)
            2 * R * D * dT else 0          # blur offset is only defined for a positive slope
  intercept_over_coef + blur
}

# auto_pmin = FALSE: fit the first k points (caller-chosen window).
# auto_pmin = TRUE : two-stage iterative Michalet optimal point selection.
# `dt_cycle` is the cycle time (for x and the blur offset); `Npts` the trajectory
# length for L(N). `R` is the motion-blur coefficient relative to dT (blur_offsets()$R_msd).
# Returns sigma_loc blur-corrected, and `intercept` (nm², the fitted offset) which is
# what the drawn line uses, since the data sit on the blurred line.
# coef_override: the MSD-vs-D coefficient, when it is NOT simply 2*ndim. A molecule
# confined to a SURFACE moves in 2 degrees of freedom however many coordinates are
# recorded, and a 2D projection of a curved surface is foreshortened by a
# geometry-dependent factor, so the divisor is not the recorded dimensionality (see
# .msd_d_eff). `ndim` still carries the GEOMETRIC dimension, because Michalet's
# optimal point count is derived for real spatial dimensions.
msd_fit_olsf <- function(tau, msd, npairs, k, ndim = 2L,
                         auto_pmin = FALSE, R = 0, Npts = NA_real_,
                         dt_cycle = NA_real_, coef_override = NULL) {
  coef_val <- if (is.null(coef_override)) 2L * ndim else coef_override
  ok_all <- is.finite(tau) & is.finite(msd) & tau > 0
  tau <- tau[ok_all]; msd <- msd[ok_all]
  n_avail <- length(tau)
  if (n_avail < 3L) return(NULL)
  # Cycle time for the Michalet x and the blur offset (fallback: the lag spacing).
  dt <- if (is.finite(dt_cycle) && dt_cycle > 0) dt_cycle
        else if (n_avail >= 2) stats::median(diff(tau)) else tau[1]

  # Unweighted linear fit on the first p points → D, σ² (his curve_fit).
  # sig2 is the blur-corrected σ²; icpt keeps the raw fitted offset for drawing.
  fit_p <- function(p) {
    p  <- max(2L, min(as.integer(p), n_avail))
    s <- tau[seq_len(p)]; y <- msd[seq_len(p)]
    ok <- is.finite(s) & is.finite(y)
    if (sum(ok) < 2L) return(NULL)
    cf <- tryCatch(coef(lm(y[ok] ~ s[ok])), error = function(e) NULL)
    if (is.null(cf) || any(!is.finite(cf))) return(NULL)
    D_p <- cf[2] / coef_val
    list(D = D_p, sig2 = .olsf_sig2(cf[1] / coef_val, D_p, R, dt),
         icpt = unname(cf[1]), p = p)
  }

  if (!auto_pmin) {
    r <- fit_p(max(3L, min(as.integer(k), n_avail)))
    if (is.null(r)) return(NULL)
    D_est <- max(0, r$D); sig2_est <- max(0, r$sig2); p_slope <- r$p; p_int <- r$p
    icpt_est <- r$icpt
  } else {
    N  <- if (is.finite(Npts) && Npts > 0) Npts else n_avail
    d  <- ndim
    # Anchor the reduced localization error x on a STABLE fit over the whole
    # [min_lag, fit_max] window.  Michalet's iteration is otherwise init-
    # dependent: a few-point fit gives a too-small σ² (hence too-small x) and
    # collapses p to the starting guess.  The full-window σ² is the stable
    # baseline that makes the optimal point count well-defined and reproducible.
    r_anc   <- fit_p(n_avail)
    if (is.null(r_anc)) return(NULL)
    D_anc   <- max(1e-9, r_anc$D)
    sig2_anc<- max(0, r_anc$sig2)
    # Separate optimal point counts for the slope (→ D) and intercept (→ σ²),
    # per Michalet & Berglund (Vogler's p_b / p_a), with his +1 shift.
    p_b <- .michalet_p(D_anc, sig2_anc, dt, N, R, d, "slope")
    p_a <- .michalet_p(D_anc, sig2_anc, dt, N, R, d, "intercept")
    p_slope <- if (is.finite(p_b)) max(3L, min(n_avail, as.integer(round(p_b)) + 1L)) else n_avail
    p_int   <- if (is.finite(p_a)) max(3L, min(n_avail, as.integer(round(p_a)) + 1L)) else n_avail
    r_s <- fit_p(p_slope); r_i <- fit_p(p_int)
    D_est    <- if (!is.null(r_s)) max(0, r_s$D)    else max(0, D_anc)
    sig2_est <- if (!is.null(r_i)) max(0, r_i$sig2) else sig2_anc
    icpt_est <- if (!is.null(r_i)) r_i$icpt else r_anc$icpt
  }
  if (!is.finite(D_est)) return(NULL)
  if (!is.finite(sig2_est)) sig2_est <- 0
  if (!is.finite(icpt_est)) icpt_est <- coef_val * sig2_est
  sig_est <- sqrt(sig2_est)
  tau_seq <- seq(min(tau), max(tau), length.out = 200)
  # Draw the line the data actually follow: slope from D, offset = fitted intercept
  # (blurred), so the overlay does not float above the points by the blur term.
  yfit    <- pmax(0, coef_val * D_est * tau_seq + icpt_est)

  list(D = D_est, alpha = 1.0, sigma_loc = sig_est, intercept = icpt_est,
       tau_fit = tau_seq, msd_fit = yfit, n_points = p_slope,
       n_points_intercept = p_int, R_blur = R,
       forced_linear = TRUE, ndim = ndim, method = "OLSF")
}

# Bela's EXACT iterative OLSF (optimized_least_squares_fitting.py). Unlike the
# anchored variant above, it does NOT pre-fit the whole window - it starts at a
# few points and iterates the optimal point count until it repeats, separately
# for slope (→ D) and intercept (→ σ²). This deliberately converges on the
# LOWEST lags (least error-prone), reproducing his numbers. `tau` is the lag
# axis (use the linearized lag×cycle axis for a faithful match); D = slope/(2d).
# coef_override as in msd_fit_olsf.
msd_fit_olsf_iterative <- function(tau, msd, ndim = 2L, R = 0, coef_override = NULL,
                                   dt_cycle = NA_real_, N = NA_real_,
                                   max_attempts = 1000L) {
  coef_val <- if (is.null(coef_override)) 2L * ndim else coef_override
  ok <- is.finite(tau) & is.finite(msd) & tau > 0
  tau <- tau[ok]; msd <- msd[ok]
  n_avail <- length(tau)
  if (n_avail < 3L) return(NULL)
  d  <- ndim
  if (!is.finite(N) || N <= 0) N <- n_avail
  if (!is.finite(dt_cycle) || dt_cycle <= 0)
    dt_cycle <- if (n_avail >= 2) stats::median(diff(tau)) else tau[1]

  # sig2 is blur-corrected (his curve_fit fits σ² inside a = 2dσ² − 4dRDdT);
  # icpt is the raw fitted offset, used to draw the line through the data.
  fit_p <- function(p) {
    p <- max(2L, min(as.integer(p), n_avail))
    s <- tau[seq_len(p)]; y <- msd[seq_len(p)]
    cf <- tryCatch(stats::coef(stats::lm(y ~ s)), error = function(e) NULL)
    if (is.null(cf) || any(!is.finite(cf))) return(NULL)
    D_p <- cf[2] / coef_val
    list(D = D_p, sig2 = .olsf_sig2(cf[1] / coef_val, D_p, R, dt_cycle),
         icpt = unname(cf[1]), p = p)
  }

  p0 <- if (n_avail %/% 10L > 3L) n_avail %/% 10L else 3L
  nps_s <- p0; nps_i <- p0
  hist_s <- p0; hist_i <- p0
  s_done <- FALSE; i_done <- FALSE
  D_opt <- NA_real_; sig2_opt <- NA_real_
  p_s_fin <- p0; p_i_fin <- p0

  for (it in seq_len(max_attempts)) {
    rs <- if (!s_done) fit_p(nps_s) else NULL
    ri <- if (!i_done) fit_p(nps_i) else NULL
    if ((!s_done && is.null(rs)) || (!i_done && is.null(ri))) break
    # Cross-coupled (his exact wiring): σ² from the INTERCEPT fit, D from SLOPE.
    D_cur    <- if (!s_done) rs$D    else D_opt
    sig2_cur <- if (!i_done) ri$sig2 else sig2_opt
    if (!s_done) {
      pb <- .michalet_p(max(1e-9, D_cur), max(0, sig2_cur), dt_cycle, N, R, d, "slope")
      new_s <- if (is.finite(pb)) as.integer(round(pb)) + 1L else n_avail
      if (new_s %in% hist_s) { D_opt <- rs$D; p_s_fin <- nps_s; s_done <- TRUE }
      else { hist_s <- c(hist_s, new_s); nps_s <- max(2L, min(n_avail, new_s)) }
    }
    if (!i_done) {
      pa <- .michalet_p(max(1e-9, D_cur), max(0, sig2_cur), dt_cycle, N, R, d, "intercept")
      new_i <- if (is.finite(pa)) as.integer(round(pa)) + 1L else n_avail
      if (new_i %in% hist_i) { sig2_opt <- ri$sig2; p_i_fin <- nps_i; i_done <- TRUE }
      else { hist_i <- c(hist_i, new_i); nps_i <- max(2L, min(n_avail, new_i)) }
    }
    if (s_done && i_done) break
  }
  if (!is.finite(D_opt))    { r <- fit_p(nps_s); if (!is.null(r)) D_opt <- r$D }
  if (!is.finite(sig2_opt)) { r <- fit_p(nps_i); if (!is.null(r)) sig2_opt <- r$sig2 }
  if (!is.finite(D_opt)) return(NULL)
  D_est <- max(0, D_opt); sig2_est <- max(0, if (is.finite(sig2_opt)) sig2_opt else 0)
  sig_est <- sqrt(sig2_est)
  # Intercept of the converged intercept fit (blurred offset) for the drawn line.
  icpt_est <- { r <- fit_p(p_i_fin); if (!is.null(r)) r$icpt else coef_val * sig2_est }
  tau_seq <- seq(min(tau), max(tau), length.out = 200)
  yfit <- pmax(0, coef_val * D_est * tau_seq + icpt_est)
  list(D = D_est, alpha = 1.0, sigma_loc = sig_est, intercept = icpt_est,
       tau_fit = tau_seq, msd_fit = yfit, n_points = p_s_fin,
       n_points_intercept = p_i_fin, R_blur = R,
       forced_linear = TRUE, ndim = ndim, method = "OLSF-iterative")
}

# Ensemble-MSD model fit with a selectable model. "linear" = the OLSF above
# (Brownian, D from the slope). "confined" = a plateauing corral model (→ local
# D + confinement size L). "directed" = Brownian + drift (→ D + velocity v).
# "anomalous" = power law (→ D + exponent α). The non-linear models fit the WHOLE
# [min_lag, max_lag] window (they need the curvature/plateau, not just low lags).
# Returns a list with $model, $D (nm²/ms), $sigma_loc (nm), $tau_fit, $msd_fit
# (nm²) and model-specific extras (L_conf, tau_c_ms, plateau_um2, v_um_s, alpha).
# Pure model-fit core, reused by the unified plot AND the per-segment analysis.
# s = τ (ms, already windowed), y = MSD (nm²). model ∈ linear/confined/directed/
# anomalous. Returns $model, $D (nm²/ms), $sigma_loc (nm), $tau_fit, $msd_fit (nm²)
# + extras (L_conf, tau_c_ms, plateau_um2, degenerate / v_um_s / alpha).
.msd_model_fit <- function(s, y, npairs, ndim, model, use_iter, t_cycle,
                           coef_override = NULL, R = MSD_R_FALLBACK) {
  coef_v <- if (is.null(coef_override)) 2L * ndim else coef_override
  if (identical(model, "linear")) {
    if (length(s) < 3L) return(NULL)
    return(if (isTRUE(use_iter))
             msd_fit_olsf_iterative(s, y, ndim = ndim, R = R, dt_cycle = t_cycle,
                                    N = length(s), coef_override = coef_override)
           else
             msd_fit_olsf(s, y, npairs, k = length(s), ndim = ndim,
                          auto_pmin = TRUE, R = R, Npts = 1e4, dt_cycle = t_cycle,
                          coef_override = coef_override))
  }
  if (length(s) < 4L) return(NULL)
  # `warnOnly = TRUE` makes nls RETURN a fit instead of erroring when it cannot
  # converge - but it still prints "Convergence failure: singular convergence (7)"
  # to the console, once per attempt, which is noise the user cannot act on (each
  # cluster × each of the three non-linear models is an attempt). The failure is not
  # ignored: it is recorded on the fit and folded into the `degenerate` flag below,
  # which is what the UI already reports. Singular convergence means the parameters
  # are not identifiable over this window - usually no plateau has been reached -
  # so the coefficients that come back must not be presented as a real fit.
  nls_fit <- function(formula, start, lower, upper) {
    f <- tryCatch(suppressWarnings(
      stats::nls(formula, data = data.frame(s = s, y = y), start = start,
                 lower = lower, upper = upper, algorithm = "port",
                 control = stats::nls.control(maxiter = 200, warnOnly = TRUE))),
      error = function(e) NULL)
    if (is.null(f)) return(NULL)
    attr(f, "converged") <- isTRUE(tryCatch(f$convInfo$isConv, error = function(e) NA))
    f
  }
  D0   <- max(1e-6, (max(y) - min(y)) / (coef_v * max(s)))
  C0   <- max(0, min(y, na.rm = TRUE))
  tseq <- seq(min(s), max(s), length.out = 200)
  if (identical(model, "confined")) {
    fit <- nls_fit(y ~ A * (1 - exp(-s / tc)) + C,
                   start = list(A = max(y) * 0.9, tc = max(s) / 3, C = C0),
                   lower = list(A = 0, tc = 1e-4, C = 0),
                   upper = list(A = Inf, tc = Inf, C = Inf))
    if (is.null(fit)) return(NULL)
    cf <- coef(fit); A <- unname(cf["A"]); tc <- unname(cf["tc"]); C <- unname(cf["C"])
    # A non-converged fit is degenerate too: singular convergence means the plateau
    # was never resolved, which is exactly what this flag already communicates.
    degenerate <- (tc > 2 * max(s)) || (A > 3 * max(y, na.rm = TRUE)) ||
                  !isTRUE(attr(fit, "converged"))
    # Confinement size L_c from the fitted plateau A. For a particle uniformly
    # exploring a d-dimensional box of side L_c, each axis contributes
    # Var(Δx) = L_c²/6, so the steady-state MSD is A = d·L_c²/6 and
    #     L_c = sqrt(6·A/d)
    # which gives the familiar sqrt(3A) in 2D and sqrt(2A) in 3D. Verified against
    # reflecting-boundary simulations (d=2 and d=3, L_c = 140/200/300 nm, recovered
    # to <1%). Using the 2D form in 3D would inflate L_c by sqrt(3/2) = 1.22×.
    list(model = "confined", D = A / (coef_v * tc), sigma_loc = sqrt(max(0, C / coef_v)),
         L_conf = sqrt(6 * A / ndim), tau_c_ms = tc, plateau_um2 = A / 1e6, degenerate = degenerate,
         tau_fit = tseq, msd_fit = A * (1 - exp(-tseq / tc)) + C, ndim = ndim, method = "confined")
  } else if (identical(model, "directed")) {
    fit <- nls_fit(y ~ coef_v * D * s + (v * s)^2 + C,
                   start = list(D = D0, v = sqrt(D0) / max(s), C = C0),
                   lower = list(D = 0, v = 0, C = 0),
                   upper = list(D = Inf, v = Inf, C = Inf))
    if (is.null(fit)) return(NULL)
    cf <- coef(fit); D <- unname(cf["D"]); v <- unname(cf["v"]); C <- unname(cf["C"])
    list(model = "directed", D = D, sigma_loc = sqrt(max(0, C / coef_v)), v_um_s = v,
         degenerate = !isTRUE(attr(fit, "converged")),
         tau_fit = tseq, msd_fit = coef_v * D * tseq + (v * tseq)^2 + C, ndim = ndim, method = "directed")
  } else if (identical(model, "anomalous")) {
    fit <- nls_fit(y ~ coef_v * D * s^a + C,
                   start = list(D = D0, a = 1, C = C0),
                   lower = list(D = 0, a = 0.05, C = 0),
                   upper = list(D = Inf, a = 2.5, C = Inf))
    if (is.null(fit)) return(NULL)
    cf <- coef(fit); D <- unname(cf["D"]); a <- unname(cf["a"]); C <- unname(cf["C"])
    list(model = "anomalous", D = D, alpha = a, sigma_loc = sqrt(max(0, C / coef_v)),
         degenerate = !isTRUE(attr(fit, "converged")),
         tau_fit = tseq, msd_fit = coef_v * D * tseq^a + C, ndim = ndim, method = "anomalous")
  } else NULL
}

# ── msd_track_sums() ──────────────────────────────────────────────────────────
# GOAL: the per-track sums behind the cycle-binned MSD. The lag of a pair is the
# difference of the two localizations' cumulative cycle index (Vogler's TOT),
# round((dt - SW) / t_c) summed along the track, so dark cycles count as elapsed
# cycles. Vectorised over index offsets: a pair m localizations apart spans at
# least m cycles, so offsets 1..max_lag find every pair with lag <= max_lag.
# df      : localizations (tid, x, y, tim in ms, optional z)
# min_pts : tracks with fewer localizations contribute nothing (zero columns)
# RETURNS list(tk_names, dr2, np, dt) - max_lag x n_tracks matrices of summed
#   squared displacement (nm^2), pair counts and summed real dt (ms) per cycle lag;
#   columns follow the sorted track ids, as split() orders them.
msd_track_sums <- function(df, t_c, sw = 0, max_lag = 50L, min_pts = 5L, has_z = FALSE) {
  df <- df[order(df$tid, df$tim), ]
  tk_names <- levels(factor(df$tid))
  nt  <- length(tk_names)
  dr2 <- matrix(0, max_lag, nt); np <- matrix(0L, max_lag, nt); dts <- matrix(0, max_lag, nt)
  out <- list(tk_names = tk_names, dr2 = dr2, np = np, dt = dts)
  if (!nt || !is.finite(t_c) || t_c <= 0) return(out)
  tf   <- match(df$tid, tk_names)
  keep <- (tabulate(tf, nt) >= min_pts)[tf]               # short tracks contribute nothing
  d <- df[keep, , drop = FALSE]; tfd <- tf[keep]; n <- nrow(d)
  if (n < 2L) return(out)
  newtrk <- c(TRUE, tfd[-1] != tfd[-n])
  dcyc <- pmax(1L, as.integer(round((c(0, diff(d$tim)) - sw) / t_c)))
  dcyc[newtrk] <- 0L
  cyc  <- cumsum(dcyc); cyc <- cyc - cyc[newtrk][cumsum(newtrk)]   # cycle index within track
  nb   <- max_lag * nt
  for (m in seq_len(min(max_lag, n - 1L))) {
    i  <- seq_len(n - m); j <- i + m
    ok <- tfd[i] == tfd[j]
    if (!any(ok)) next
    i <- i[ok]; j <- j[ok]
    lag <- cyc[j] - cyc[i]
    ok  <- lag >= 1L & lag <= max_lag
    if (!any(ok)) next
    i <- i[ok]; j <- j[ok]; lag <- lag[ok]
    v <- (d$x[j] - d$x[i])^2 + (d$y[j] - d$y[i])^2
    if (has_z) v <- v + (d$z[j] - d$z[i])^2
    key <- lag + (tfd[i] - 1L) * max_lag                   # column-major cell of (lag, track)
    s <- rowsum(v, key); ix <- as.integer(rownames(s))
    dr2[ix] <- dr2[ix] + s[, 1L]
    s <- rowsum(d$tim[j] - d$tim[i], key); ix <- as.integer(rownames(s))
    dts[ix] <- dts[ix] + s[, 1L]
    np <- np + tabulate(key, nbins = nb)
  }
  list(tk_names = tk_names, dr2 = dr2, np = np, dt = dts)
}

# ── msd_ensemble() ────────────────────────────────────────────────────────────
# GOAL: the ensemble MSD curve from msd_track_sums(): track-weighted (the mean of
# per-track TAMSDs, Vogler's eMSD) or pair-weighted, with the real mean dt per lag,
# the linearised axis k * t_c + SW, and which lags enough tracks reach.
msd_ensemble <- function(sums, t_c, sw = 0, track_weight = TRUE) {
  L <- nrow(sums$dr2)
  track_msd <- ifelse(sums$np > 0L, sums$dr2 / sums$np, NA_real_)
  np_lag <- rowSums(sums$np)
  msd <- if (track_weight) rowMeans(track_msd, na.rm = TRUE)
         else ifelse(np_lag > 0, rowSums(sums$dr2) / np_lag, NA_real_)
  msd[!is.finite(msd)] <- NA_real_
  n_trk_lag <- rowSums(!is.na(track_msd))
  keep_lag  <- n_trk_lag >= max(2L, floor(0.10 * max(c(0, n_trk_lag))))
  list(lag = seq_len(L), msd = msd, npairs = np_lag, n_trk_lag = n_trk_lag, keep_lag = keep_lag,
       tau = ifelse(np_lag > 0, rowSums(sums$dt) / np_lag, seq_len(L) * t_c),
       tau_lin = seq_len(L) * t_c + sw, track_msd = track_msd)
}

# ── msd_headline() ────────────────────────────────────────────────────────────
# GOAL: D from the upper lag regime of the ensemble MSD, similar to the estimator Vogler et
# al. 2025. A similar method where the MSD panel runs (OLSF over
# the fit window, iterative or anchored). Short lags are left out because they carry
# non-diffusive motion on top of noise (probe wobble, tracker wander).
# The OLSF fits only the first p lags of the window, p from Michalet's optimal
# point count (often 4-6), so the lags actually fitted are returned beside the
# window: on a bending MSD, D depends on which lags those are.
# window : c(lower, upper) fit window in cycles
# R      : motion-blur coefficient per step; coef : 2 x effective dimensionality
# RETURNS list(D (nm^2/ms), sigma (nm), n_points, window, fit_lags = first and
#         last lag fitted, fit_tau = their lag times in ms) or NULL.
msd_headline <- function(df, t_c, sw = 0, window = c(5L, 50L), track_weight = TRUE,
                         iterative = TRUE, linear_tau = TRUE, R = 0, coef = 4,
                         min_pts = 50L, has_z = FALSE) {
  lo <- max(1L, as.integer(window[1])); hi <- max(lo + 2L, as.integer(window[2]))
  df <- df[is.finite(df$x) & is.finite(df$y) & is.finite(df$tim) &
             (if (has_z) is.finite(df$z) else TRUE), , drop = FALSE]
  sums <- msd_track_sums(df, t_c, sw, max(60L, min(5000L, hi)), min_pts, has_z)
  ens  <- msd_ensemble(sums, t_c, sw, track_weight)
  tau  <- if (linear_tau) ens$tau_lin else ens$tau
  ok   <- is.finite(tau) & is.finite(ens$msd) & ens$keep_lag & ens$lag >= lo & ens$lag <= hi
  if (sum(ok) < 3L) return(NULL)
  ndim <- if (has_z) 3L else 2L
  fr <- if (iterative)
    msd_fit_olsf_iterative(tau[ok], ens$msd[ok], ndim = ndim, R = R, dt_cycle = t_c,
                           N = sum(ok), coef_override = coef)
  else
    msd_fit_olsf(tau[ok], ens$msd[ok], ens$npairs[ok], k = sum(ok), ndim = ndim,
                 auto_pmin = TRUE, R = R, Npts = 1e4, dt_cycle = t_c, coef_override = coef)
  if (is.null(fr) || !is.finite(fr$D)) return(NULL)
  # the lags the slope fit used: the first n_points of the usable lags in the window
  np    <- min(sum(ok), as.integer(fr$n_points %||% sum(ok)))
  used  <- which(ok)[seq_len(max(1L, np))]
  list(D = fr$D, sigma = fr$sigma_loc, n_points = fr$n_points %||% NA_integer_, window = c(lo, hi),
       fit_lags = range(ens$lag[used]), fit_tau = range(tau[used]))
}

# ── msd_unified_curves() ──────────────────────────────────────────────────────
# GOAL: everything the MSD panel plots and fits, from the per-track sums of
# msd_track_sums(): per-track TAMSD, the ensemble curve (track- or pair-weighted),
# its quartile band and median, the real and the linearised lag axes, and which
# lags enough tracks reach. Moved unchanged out of msd_unified_bundle()
# (server/tracking/15_msd.R).
# sums         : result of msd_track_sums()
# t_c, cc_sw   : cycle time and per-step overhead (ms)
# max_lag      : number of cycle lags in `sums`
# has_z        : whether z entered the squared displacements
# drop_negd    : drop tracks whose TAMSD has a non-positive slope before averaging
# track_weight : TRUE = mean of per-track TAMSDs (eMSD), FALSE = pair-weighted
# d_eff        : effective dimensionality of the motion as observed (MSD = 2 d_eff D tau)
# RETURNS list(lag, tau, tau_lin, msd, npairs, n_trk_lag, keep_lag, track_msd,
#   tk_names, msd_q25, msd_q75, msd_med, t_cycle, has_z, ndim, d_eff, coef_eff,
#   track_weighted, n_dropped_negD, n_tracks).
msd_unified_curves <- function(sums, t_c, cc_sw = 0, max_lag, has_z = FALSE, drop_negd = FALSE,
                               track_weight = FALSE, d_eff = if (has_z) 3 else 2) {
  tk_names <- sums$tk_names
  tk_dr2 <- sums$dr2; tk_np <- sums$np; tk_dt <- sums$dt

  track_msd <- ifelse(tk_np > 0L, tk_dr2 / tk_np, NA_real_)   # per-track TAMSD (nm²)

  # Drop non-positive-D tracks before the ensemble. OFF by default to match
  # Bela's eDIFF, which keeps EVERY track (his D≤0 filter applies only to the
  # per-track tDIFF, which we still enforce there). Toggle on for the previous
  # QC behaviour.
  n_dropped_negD <- 0L
  if (isTRUE(drop_negd)) {
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
       # recorded coordinate count ; the caller supplies it.
       d_eff = d_eff,
       coef_eff = 2 * d_eff,
       track_weighted = track_weight, n_dropped_negD = n_dropped_negD,
       n_tracks = sum(colSums(tk_np) > 0))
}

# ── msd_pertrack_fits() ───────────────────────────────────────────────────────
# GOAL: fit every track's own TAMSD over the fit window with the same OLSF as the
# ensemble, and summarise the per-track D (the time-averaged D, tDIFF). Also the
# confinement index of each track: late-window rise over early-window rise of its
# MSD, relative to the same ratio of lag times (< 1 confined, ~1 Brownian, > 1
# directed). Only tracks with D > 0 enter. Moved unchanged out of
# trk_unified_tdiff() (server/tracking/15_msd.R).
# ub         : result of msd_unified_curves()
# min_lag, max_lag_e : fit window in cycles
# use_iter   : iterative OLSF (TRUE) or anchored OLSF with Michalet's p_min
# linear_tau : fit against the linearised lag axis (k t_c + SW) or the real mean dt
# R_msd      : motion-blur coefficient per step
# RETURNS list(tDIFF (um2/s), n_tD, sd_tD, ci_lo, ci_hi, med_tD, per_D, pertrack).
msd_pertrack_fits <- function(ub, min_lag = 5L, max_lag_e = 50L, use_iter = TRUE, linear_tau = TRUE,
                              R_msd = 0) {
  tau_sel   <- if (isTRUE(linear_tau)) ub$tau_lin else ub$tau
  win   <- ub$lag >= min_lag & ub$lag <= max_lag_e & is.finite(tau_sel)
  tau_w <- tau_sel[win]; t_c_e <- ub$t_cycle
  tids  <- if (!is.null(ub$tk_names)) ub$tk_names else as.character(seq_len(ncol(ub$track_msd)))
  rows <- lapply(seq_len(ncol(ub$track_msd)), function(j) {
    col <- ub$track_msd[, j]
    mv  <- col[win]
    okc <- is.finite(mv); if (sum(okc) < 3L) return(NULL)
    tt <- tau_w[okc]; yy <- mv[okc]; N_tr <- sum(is.finite(col))
    fr_tr <- if (use_iter)
      msd_fit_olsf_iterative(tt, yy, ndim = ub$ndim, R = R_msd, dt_cycle = t_c_e, N = N_tr,
                             coef_override = ub$coef_eff)
    else
      msd_fit_olsf(tt, yy, rep(1L, length(tt)), k = length(tt), ndim = ub$ndim,
                   auto_pmin = TRUE, R = R_msd, Npts = N_tr, dt_cycle = t_c_e,
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
}
