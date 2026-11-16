## ─────────────────────────────────────────────────────────────────────────────
## sim_acquisition.R  -  simulator: detection-rate model and calibration from a real run
## ─────────────────────────────────────────────────────────────────────────────
##
## PURPOSE
##   The acquisition model of the Simulator tab: the shape of the trace-detection rate over an acquisition, sampling detection times from it, the resolution-over-time timeline, and the analysis and quality gate of a real MINFLUX run used to calibrate those.
##
##   Plain functions with no Shiny state: they take numbers and tables and return
##   numbers and tables, so they can be called and tested without a running app.
##   Moved unchanged out of the sim_server() closure (sim_server.R); the server
##   gathers the inputs and calls these.
##
## RANDOM NUMBERS: the generators draw from R's global RNG in a fixed order. The
##   golden tests (tests/testthat/test-golden-sim.R) pin that order for a fixed seed.
## UNITS: nm, seconds for SMLM acquisition time, ms for tracking.
## ─────────────────────────────────────────────────────────────────────────────

# ── Detection-rate SHAPE models ─────────────────────────────────────────────
# The instantaneous structure detection rate is  struct_rate · g(t),  where
# g(t) ∈ [0, 1] is a normalized SHAPE and struct_rate is the CHARACTERISTIC
# (reference) rate - the peak for a decaying run, the plateau for a rising or
# constant one. Three shapes cover the profiles seen in real MINFLUX runs:
#
#   "decay"    g(t) = level + (1-level)·e^{-t/τ}     starts at 1, → level
#              Blinking dyes: bright early, bleach toward a floor. `level` =
#              late/peak ratio (0 = decays to zero, 1 = no decay).
#   "constant" g(t) = 1                              flat (variation is Poisson)
#   "rise"     g(t) = 1 − (1-level)·e^{-t/τ}         starts at `level`, → 1
#              Activation-limited: almost no detection early, then rises to a
#              cap. `level` = initial/plateau ratio (0 = starts from zero).
#
# `t` and `tau` must share units (both seconds, or both minutes). level clamped
# to [0,1]. The decay case is identical to the old hard-coded formula, so
# existing calibrations are unchanged.
acq_shape <- function(t, model, tau, level) {
  tau <- max(1e-6, tau); lv <- max(0, min(1, level))
  e   <- exp(-t / tau)
  switch(model %||% "decay",
    constant = rep(1, length(t)),
    rise     = 1 - (1 - lv) * e,
    lv + (1 - lv) * e)                 # decay (default)
}

# Closed-form ∫₀ᵀ g(t) dt  (T, tau in the SAME units). Used for expected counts.
acq_shape_integral <- function(T, model, tau, level) {
  tau <- max(1e-6, tau); lv <- max(0, min(1, level))
  ee  <- tau * (1 - exp(-T / tau))     # ∫₀ᵀ e^{-t/τ} dt
  switch(model %||% "decay",
    constant = T,
    rise     = T - (1 - lv) * ee,
    lv * T + (1 - lv) * ee)            # decay (default)
}

# Sample n detection times on [0, T] from the (front- or back-loaded) rate
# shape, by inverse-CDF sampling of the density g(t). model/level select the
# profile (see acq_shape); tau_sec is the rate-change time constant.
sim_sample_detection_times <- function(n, T_sec, tau_sec, floor_frac,
                                       model = "decay") {
  if (n <= 0) return(numeric(0))
  if (!is.finite(T_sec) || T_sec <= 0) return(rep(0, n))
  ng   <- 2000L
  tg   <- seq(0, T_sec, length.out = ng)
  dens <- pmax(1e-9, acq_shape(tg, model, tau_sec, floor_frac))
  cdf  <- cumsum(dens); cdf <- cdf / cdf[ng]
  approx(cdf, tg, runif(n), rule = 2, ties = "ordered")$y
}

# Expected number of structure detections over [0, T] for an ABSOLUTE rate
# r(t) = struct_rate · g(t)  [traces/min]:  N_T = struct_rate · ∫₀ᵀ g dt.
# (T, τ in min.) For decay this is the old floor·T + (1-floor)·τ·(1-e^{-T/τ}).
acq_expected_count <- function(struct_rate, T_sec, decay_min, floor_frac,
                               model = "decay") {
  T_min <- T_sec / 60
  struct_rate * acq_shape_integral(T_min, model, max(1e-6, decay_min), floor_frac)
}

# Acquisition timeline from a fixed fluorophore pool + live acquisition
# parameters.  Operates only on molecule IDs (no positions/localizations), so
# it is cheap enough to recompute live as rate / duration / background change
# without re-running the whole structure.  pool_mol = mol_id of every unique
# fluorophore (length = pool size); n_mol_total = molecules placed (for
# display).  Returns the same shape of list the acquisition plot/summary use.
compute_acq_timeline <- function(pool_mol, n_mol_total, struct_rate, T_sec,
                                 blinks, decay_min, floor_frac, bg_rate,
                                 resolve_q, crowd_m = NULL, model = "decay") {
  n_pool <- length(pool_mol)
  if (n_pool == 0 || !is.finite(T_sec) || T_sec <= 0 || struct_rate <= 0)
    return(NULL)
  n_mol_resolvable <- length(unique(pool_mol))

  # Deterministic EXPECTED timeline (no per-realization jitter), so the curves
  # and the time-to-resolve are stable as the user drags rate / duration.
  # struct_rate / bg_rate are TRACES PER HOUR (easier to match real settings);
  # converted to per-minute internally.  Detections follow an ABSOLUTE rate
  # r(t) = X·(floor + (1-floor)·e^{-t/τ}), so how fast the structure resolves
  # depends on the rate and structure - NOT on the chosen duration.
  T_min     <- T_sec / 60
  tau_min   <- max(1e-6, decay_min)
  rate_min  <- struct_rate / 60                    # per hour → per minute
  bg_min    <- bg_rate / 60

  # Photon budget: each fluorophore yields ~blinks detections before bleaching.
  # This caps the TOTAL TIDs (over-counting) only - it does NOT gate resolution
  # (a molecule is resolved by its FIRST detection; the extra blinks just pile
  # more TIDs on top).
  struct_cap <- n_pool * max(1, blinks)
  f_m        <- as.numeric(tapply(rep(1, n_pool), pool_mol, sum))  # dyes / molecule

  ng     <- 240L
  tg_min <- seq(0, T_min, length.out = ng)
  tg_sec <- tg_min * 60
  # Cumulative fluorophore-detection attempts from the absolute rate:
  # ∫₀ᵗ g(s) ds for the chosen shape (decay / constant / rise).
  cum_R <- rate_min * vapply(tg_min,
             function(tt) acq_shape_integral(tt, model, tau_min, floor_frac),
             numeric(1))
  # Photobleaching consumes the photon budget GRADUALLY (fluorophores bleach at
  # different times), so the total saturates smoothly toward struct_cap and the
  # structure rate tapers off - no hard cliff when the budget is reached.
  bleach <- exp(-cum_R / struct_cap)          # fraction of budget still emitting
  N_t    <- struct_cap * (1 - bleach)         # total recorded TIDs (over-counting)
  n_s    <- N_t[ng]

  # Expected distinct molecules resolved by time t: a molecule is found once
  # ANY of its f_m fluorophores is first detected.  Driven by cum_R (rate ×
  # time × pool) - INDEPENDENT of the per-fluorophore blink budget, so
  # "detections before bleaching" never shifts the road to X% resolved.
  # crowd_m (≤ 1, per molecule, aligned to the sorted unique mol_ids) slows
  # discovery where emitters crowd the donut: a crowded molecule's effective
  # detection rate is scaled down, so dense structures resolve slower.
  cw <- if (is.null(crowd_m) || length(crowd_m) != length(f_m))
          rep(1, length(f_m)) else crowd_m
  cum_unique <- vapply(cum_R,
                       function(c) sum(1 - exp(-f_m * cw * c / n_pool)),
                       numeric(1))
  n_mol_det  <- cum_unique[ng]
  denom      <- max(1, n_mol_resolvable)

  # Instantaneous detection rate (traces/HOUR): structure (the chosen shape g(t)
  # × the photobleaching taper) + constant background floor.
  rate_struct <- struct_rate * acq_shape(tg_min, model, tau_min, floor_frac) * bleach
  rate_v <- rate_struct + bg_rate
  n_bg   <- max(0L, min(200000L, as.integer(round(bg_min * T_min))))

  q        <- max(0.5, min(1, resolve_q))
  res_time <- if (n_mol_det >= q * denom)
                approx(cum_unique, tg_sec, q * denom, ties = "ordered")$y
              else NA_real_

  list(
    T_sec = T_sec, struct_rate = struct_rate, bg_rate = bg_rate,
    n_pool = n_pool, n_mol_total = n_mol_total,
    n_mol_resolvable = n_mol_resolvable, struct_cap = struct_cap, n_bg = n_bg,
    t = tg_sec, cum_unique = cum_unique, rate_t = tg_sec, rate_v = rate_v,
    rate_struct = rate_struct,          # structure-only component of rate_v
    rate_peak  = max(rate_v),
    rate_floor = stats::median(tail(rate_v, max(1L, length(rate_v) %/% 4L))),
    rate_mean  = n_s / max(1e-9, T_min) * 60,      # traces/hour
    n_mol_det = round(n_mol_det), resolved_frac = n_mol_det / denom,
    n_real_tid = round(n_s), det_per_mol = n_s / max(1, n_mol_det),
    resolve_q = q, res_time = res_time,
    overcount_frac = 1 - n_mol_det / max(1, n_s))
}

# Calibrate the acquisition model from a REAL MINFLUX run: measure the
# detection-rate profile (peak / mean / floor / decay τ), the run duration and
# the FOV extent from the data.  Works on the first detection time of each
# unique trace (tid).  df = canonical MINFLUX frame (tid, tim [s], x, y[, z] nm).
analyze_real_run <- function(df) {
  if (!all(c("tid", "tim", "x", "y") %in% names(df))) return(NULL)
  df <- df[is.finite(df$tim) & is.finite(df$x) & is.finite(df$y), , drop = FALSE]
  if (nrow(df) < 10) return(NULL)

  # First detection time + median position of each trace (same tid order).
  t0_tid <- as.numeric(tapply(df$tim, df$tid, min))
  tx     <- as.numeric(tapply(df$x,   df$tid, median))
  ty     <- as.numeric(tapply(df$y,   df$tid, median))
  t0_tid <- t0_tid - min(t0_tid)                    # seconds from start
  n_tr   <- length(t0_tid)
  T_sec  <- max(t0_tid)
  if (!is.finite(T_sec) || T_sec <= 0) return(NULL)

  # Binned detection rate (traces / hour).
  nb     <- max(8L, min(60L, n_tr %/% 10L))
  br     <- seq(0, T_sec, length.out = nb + 1L)
  bin_hr <- max(1e-9, (br[2] - br[1]) / 3600)
  h      <- hist(t0_tid, breaks = br, plot = FALSE)
  rate_t <- h$mids                                  # seconds
  rate_v <- h$counts / bin_hr

  R0     <- mean(head(rate_v, max(1L, nb %/% 10L)))           # early rate
  Rf     <- stats::median(tail(rate_v, max(1L, nb %/% 4L)))   # late rate

  R_mean <- n_tr / (T_sec / 3600)

  # τ from a log-linear fit of |rate − late_rate| toward 0 - works whether the
  # rate falls (dye bleaching) OR rises (e.g. buffer activating over the run).
  tau_min <- 15
  exc     <- abs(rate_v - Rf)
  ok      <- exc > 0.05 * max(exc, 1e-9) & rate_t < T_sec * 0.8
  if (sum(ok) >= 3) {
    fit <- tryCatch(stats::lm(log(exc[ok]) ~ rate_t[ok]), error = function(e) NULL)
    if (!is.null(fit)) {
      slope <- stats::coef(fit)[[2]]
      if (is.finite(slope) && slope < 0) tau_min <- max(0.5, (-1 / slope) / 60)
    }
  }
  # If the rate visibly changed over the run, τ can't exceed it - a noisy
  # (especially rising) profile can fit an absurd τ ≫ run, which would make the
  # simulated change look flat.  Cap at the run length.
  tau_min <- min(tau_min, T_sec / 60)

  # ── Choose the detection-rate MODEL that best fits the profile ─────────────
  # All three are the same family  A + B·e^{-t/τ}:  decay B>0 (high→low),
  # rise B<0 (low→high), constant B≈0. Compare the exponential fit's residual
  # to a flat line; if it barely helps, or the endpoints hardly differ, call it
  # constant. Otherwise the SIGN of (early − late) picks decay vs rise, and the
  # level/reference are set so struct_rate is the peak (decay) or plateau (rise).
  okv       <- is.finite(rate_v)
  rt        <- rate_t[okv]; rv <- rate_v[okv]
  rss_const <- sum((rv - mean(rv))^2)
  pred_exp  <- Rf + (R0 - Rf) * exp(-rt / (tau_min * 60))
  rss_exp   <- sum((rv - pred_exp)^2)
  rel_chg   <- abs(R0 - Rf) / max(1e-9, R_mean)
  # Only fall back to "constant" when the exponential barely helps (< 5% lower
  # residual than a flat line) AND the endpoints differ by < 10% of the mean.
  # (Loosened from 0.9 / 0.15 so a shallow but real decay is not called flat.)
  if (rss_exp >= 0.95 * rss_const || rel_chg < 0.10) {
    det_model <- "constant"; level_frac <- 1;                       ref_rate <- R_mean
  } else if (R0 >= Rf) {
    det_model <- "decay";    level_frac <- max(0, min(1, Rf / max(1e-9, R0))); ref_rate <- R0
  } else {
    det_model <- "rise";     level_frac <- max(0, min(1, R0 / max(1e-9, Rf))); ref_rate <- Rf
  }

  fov_x <- diff(range(df$x, na.rm = TRUE))
  fov_y <- diff(range(df$y, na.rm = TRUE))
  fov_z <- if ("z" %in% names(df)) diff(range(df$z, na.rm = TRUE)) else 0

  # Spatial split: a trace sitting in a locally dense neighbourhood is ON the
  # structure; an isolated trace is background.  Uses the DBSCAN core-point
  # criterion (≥ minPts neighbours within eps) - absolute density, so it works
  # even when the structure itself is densely labeled.  Approximate; the user
  # can override the resulting fraction.  Also bins the two classes separately
  # so the measured structure / background rate curves can be overlaid.
  struct_frac   <- 0.5
  rate_v_struct <- NULL
  rate_v_bg     <- NULL
  if (n_tr >= 20) {
    eps    <- 60          # nm neighbourhood radius
    minPts <- 3L          # core-point threshold
    K      <- min(n_tr, 64L)
    nn     <- RANN::nn2(cbind(tx, ty), cbind(tx, ty), k = K,
                        searchtype = "radius", radius = eps)
    cnt    <- rowSums(nn$nn.idx > 0) - 1L          # neighbours within eps
    is_struct   <- cnt >= minPts
    struct_frac <- max(0, min(1, mean(is_struct)))
    rate_v_struct <- hist(t0_tid[is_struct],  breaks = br, plot = FALSE)$counts / bin_hr
    rate_v_bg     <- hist(t0_tid[!is_struct], breaks = br, plot = FALSE)$counts / bin_hr
  }

  # ── Point-cloud profile ───────────────────────────────────────────────────
  # Localization precision = SD of the repeated final localizations WITHIN a
  # trace (the molecule re-localized while ON), per axis, median over traces
  # with >= 3 localizations.  Photons/loc from the effective count (eco).
  # The CRLB-mode k and σ₀ are fitted from these per-trace variances below.
  #
  # CAVEAT (stated in the UI too): a within-trace SD is not pure localization
  # precision — it also contains any real molecular motion, stage drift and
  # sample wobble during the trace.  On a static target that is small; on
  # anything mobile this is an UPPER BOUND on precision, not an estimate of it.
  nl  <- as.numeric(tapply(df$x, df$tid, length))
  sx  <- as.numeric(tapply(df$x, df$tid, stats::sd))
  sy  <- as.numeric(tapply(df$y, df$tid, stats::sd))
  sz  <- if ("z" %in% names(df))
           as.numeric(tapply(df$z, df$tid, stats::sd))
         else rep(NA_real_, length(nl))
  okp <- nl >= 3 & is.finite(sx) & is.finite(sy)

  # Small-n bias correction.  The sample SD under-estimates sigma at small n,
  # which matters now that the user controls the trace-length gate and can let
  # short traces dominate.  Each trace's s is divided by its OWN factor.
  #
  # NOT c4.  The textbook constant c4(n) = E[s]/sigma un-biases the MEAN of the
  # sampling distribution, but this aggregates with a MEDIAN, and median(s)
  # sits BELOW E[s] at small n.  Dividing by c4 under a median therefore leaves
  # about 2.5 % of the bias in place (measured, short-trace-heavy population).
  # The matching constant for a median aggregator is the median of the same
  # chi distribution:  m(n) = sqrt(qchisq(0.5, n-1) / (n-1)).
  #   n = 3:  c4 = 0.886  but  m = 0.833   (raw median is 7.6 % low)
  # Using m recovers sigma to within ~0.2 % of truth.
  # NOTE: either constant CANCELS exactly in the axis ratios below (both axes
  # share the same n), so the ratios do not depend on this choice.
  sd_median_debias <- function(n) {
    n <- pmax(2, n)
    sqrt(stats::qchisq(0.5, n - 1) / (n - 1))
  }
  cc   <- sd_median_debias(nl)
  sx_c <- sx / cc; sy_c <- sy / cc; sz_c <- sz / cc

  # Median and spread (nMAD) of the per-trace sigma.  The SPREAD is the point:
  # real MINFLUX sigma is heteroscedastic (it tracks photon count), and a
  # single median throws exactly that away.
  med_mad <- function(v) {
    v <- v[is.finite(v)]
    if (!length(v)) return(c(NA_real_, NA_real_))
    c(stats::median(v), 1.4826 * stats::mad(v, constant = 1))
  }
  mx <- med_mad(sx_c[okp]); my <- med_mad(sy_c[okp])
  mz <- med_mad(sz_c[okp & is.finite(sz_c)])
  sigma_x <- mx[1]; sigma_y <- my[1]; sigma_z <- mz[1]

  # Axis ratios, as the MEDIAN OF PER-TRACE RATIOS rather than the ratio of
  # medians: within one trace both axes share the same photon count, so the
  # per-trace ratio already divides the photon dependence out and measures
  # cloud SHAPE alone.  This is what an anisotropic simulated cloud needs.
  okr     <- okp & is.finite(sx) & sx > 0
  ratio_yx <- if (any(okr)) stats::median((sy / sx)[okr]) else NA_real_
  okrz     <- okr & is.finite(sz)
  ratio_zx <- if (any(okrz)) stats::median((sz / sx)[okrz]) else NA_real_

  # Distributions kept whole (not just their medians) — 8.1c resamples these
  # to give every simulated localization its own photon count and every
  # simulated trace its own length.  Capped so a huge run cannot bloat the
  # stored profile; resampling from 200 k is indistinguishable from 5 M.
  CAP  <- 200000L
  subs <- function(v) {
    v <- v[is.finite(v)]
    if (length(v) > CAP) v[sample.int(length(v), CAP)] else v
  }
  eco_vec <- if ("eco" %in% names(df)) subs(as.numeric(df$eco)) else numeric(0)
  len_vec <- as.integer(nl[is.finite(nl)])       # ALL gated traces, not just okp

  photons  <- if ("eco" %in% names(df)) stats::median(df$eco, na.rm = TRUE) else NA_real_

  # ── Precision model fit: σ² = A/N + B ────────────────────────────────────
  # The CRLB-mode model is σ² = (k·σ_CRLB)² + σ₀², and σ_CRLB² ∝ 1/N, so the
  # within-trace variance is LINEAR in 1/N: slope A = photon-limited term,
  # intercept B = σ₀² (drift, vibration, motion). One (σ, N) pair cannot split
  # the two; the spread of photon counts across traces can. Traces are binned
  # by their mean eco and each bin's variance is POOLED (Σ(n−1)s² / Σ(n−1)),
  # which is unbiased at any n, so no median debias is needed here.
  prec_fit <- NULL
  if ("eco" %in% names(df)) {
    # Per-trace mean of 1/N, not 1/mean(N): σ² is linear in 1/N per
    # localization, so a trace's variance averages 1/N (Jensen: 1/mean is low).
    inv_all <- ifelse(is.finite(df$eco) & df$eco > 0, 1 / df$eco, NA_real_)
    inv_tr  <- as.numeric(tapply(inv_all, df$tid, mean, na.rm = TRUE))
    okf <- okp & is.finite(inv_tr) & inv_tr > 0
    if (sum(okf) >= 30) {
      v_tr  <- (sx[okf]^2 + sy[okf]^2) / 2       # per-trace variance, xy-averaged
      dof   <- nl[okf] - 1                       # degrees of freedom per trace
      inv_n <- inv_tr[okf]
      nb    <- max(3L, min(8L, floor(sum(okf) / 20)))   # ≥ ~20 traces per bin
      brk   <- unique(stats::quantile(inv_n, seq(0, 1, length.out = nb + 1)))
      if (length(brk) >= 4) {
        bin  <- cut(inv_n, brk, include.lowest = TRUE)
        sdof <- tapply(dof, bin, sum)
        pv   <- tapply(v_tr * dof, bin, sum) / sdof     # pooled variance per bin
        px   <- tapply(inv_n * dof, bin, sum) / sdof    # matching mean 1/N
        ok_b <- is.finite(pv) & is.finite(px) & sdof > 0
        pv <- pv[ok_b]; px <- px[ok_b]; sdof <- sdof[ok_b]
        # Var(pooled s²) ≈ 2σ⁴/dof, so weight each bin by dof / σ⁴.
        w  <- sdof / pmax(pv, 1e-6)^2
        A <- NA_real_; B <- NA_real_; se_A <- NA_real_; se_B <- NA_real_
        if (length(pv) >= 3) {
          fit <- stats::lm(pv ~ px, weights = w)
          cf  <- summary(fit)$coefficients
          B <- cf[1, 1]; A <- cf[2, 1]; se_B <- cf[1, 2]; se_A <- cf[2, 2]
          # Both terms are variances, so neither may be negative. A negative
          # one means the data cannot support it: drop it and refit the other.
          if (B < 0) {
            B <- 0; se_B <- NA_real_
            A <- sum(w * px * pv) / sum(w * px^2)
          }
          if (A < 0) {
            A <- 0; se_A <- NA_real_
            B <- sum(w * pv) / sum(w)
          }
        }
        # 1/N must span a useful range, or slope and intercept trade off freely.
        span <- if (length(px) >= 2) max(px) / min(px) else 1
        prec_fit <- list(
          A = A, B = B, se_A = se_A, se_B = se_B,
          c_phot   = sqrt(max(0, A)),            # photon term = c_phot / √N, nm
          floor    = sqrt(max(0, B)),            # σ₀, nm
          # SE of σ₀ by the delta method, from SE of B = σ₀².
          se_floor = if (is.finite(se_B) && B > 0) se_B / (2 * sqrt(B)) else NA_real_,
          span = span, n_bins = length(pv), n_traces = sum(okf),
          ok   = is.finite(A) && is.finite(B) && span >= 1.5,
          bin_n = 1 / px, bin_sigma = sqrt(pv))
      }
    }
  }

  list(n_traces = n_tr, duration_hr = T_sec / 3600,
       rate_peak = R0, rate_mean = R_mean, rate_floor = Rf,
       det_model = det_model, level_frac = level_frac, ref_rate = ref_rate,
       tau_min = tau_min, struct_frac = struct_frac,
       fov_x = fov_x, fov_y = fov_y, fov_z = fov_z,
       rate_t = rate_t, rate_v = rate_v,
       rate_v_struct = rate_v_struct, rate_v_bg = rate_v_bg,
       sigma_x = sigma_x, sigma_y = sigma_y, sigma_z = sigma_z,
       sigma_x_mad = mx[2], sigma_y_mad = my[2], sigma_z_mad = mz[2],
       ratio_yx = ratio_yx, ratio_zx = ratio_zx,
       n_sigma_traces = sum(okp),
       debias_med = if (any(okp)) stats::median(cc[okp]) else NA_real_,
       eco_vec = eco_vec, eco_n = if ("eco" %in% names(df)) sum(is.finite(df$eco)) else 0L,
       len_vec = len_vec,
       len_med = if (length(len_vec)) stats::median(len_vec) else NA_real_,
       len_q = if (length(len_vec))
                 stats::quantile(len_vec, c(0.25, 0.75), names = FALSE)
               else c(NA_real_, NA_real_),
       photons = photons, prec_fit = prec_fit)
}

# ── Calibration gate ────────────────────────────────────────────────────────
# GOAL: restrict a real run to the localizations / traces the calibration is
# measured from, and report exactly how much was dropped.
#
# Order is deliberate: EFO and CFR gate individual LOCALIZATIONS first, then
# trace length is recounted on the survivors and the length gate drops whole
# TRACES.  So "trace length" here means USABLE localizations per trace - gating
# hard on EFO and then asking for >= 20 locs is a far stricter request than
# either gate alone, which is why the info box reports both stages.
#
# Non-finite values are KEPT (a gate cannot judge a missing number) and counted
# separately, so leaving a range at its full extent is always a true no-op.
# A NULL range is a no-op too: the column is absent, or the widget is not built.
.rng_ok <- function(r)
  length(r) == 2L && all(is.finite(r)) && r[2] >= r[1]

# A column is usable only if it exists AND carries at least one finite value.
.col_ok <- function(df, nm)
  nm %in% names(df) &&
    any(is.finite(suppressWarnings(as.numeric(df[[nm]]))))

gate_real_run <- function(df, efo_rng = NULL, cfr_rng = NULL, len_rng = NULL) {
  if (is.null(df) || !all(c("tid", "tim", "x", "y") %in% names(df))) return(NULL)
  n_loc_in <- nrow(df)
  n_tr_in  <- length(unique(df$tid))

  # ── Stage 1: localization-level gates (EFO, CFR) ──────────────────────────
  keep    <- rep(TRUE, n_loc_in)
  n_na    <- c(efo = 0L, cfr = 0L)   # values a gate could not judge
  narrowed <- character(0)           # gates the user actually tightened
  for (nm in c("efo", "cfr")) {
    rng <- if (identical(nm, "efo")) efo_rng else cfr_rng
    if (!.rng_ok(rng) || !.col_ok(df, nm)) next
    v    <- suppressWarnings(as.numeric(df[[nm]]))
    fin  <- is.finite(v)
    keep <- keep & (!fin | (v >= rng[1] & v <= rng[2]))
    n_na[nm] <- sum(!fin)
    full <- range(v[fin])
    if (rng[1] > full[1] || rng[2] < full[2]) narrowed <- c(narrowed, toupper(nm))
  }
  g <- df[keep, , drop = FALSE]
  n_loc_mid <- nrow(g)
  n_tr_mid  <- length(unique(g$tid))

  # ── Stage 2: trace-level gate on the SURVIVING localization count ─────────
  if (.rng_ok(len_rng) && n_loc_mid > 0) {
    tl   <- table(g$tid)
    ok   <- names(tl)[tl >= len_rng[1] & tl <= len_rng[2]]
    g    <- g[as.character(g$tid) %in% ok, , drop = FALSE]
    full <- range(as.integer(tl))
    if (len_rng[1] > full[1] || len_rng[2] < full[2])
      narrowed <- c(narrowed, "trace length")
  }

  # Traces with >= 3 localizations are the only ones a within-trace SD can be
  # measured from, so they are the sample size that actually sets sigma.
  n_sig <- if (nrow(g) > 0) sum(table(g$tid) >= 3L) else 0L

  list(df = g, narrowed = narrowed, n_na = n_na,
       n_loc_in = n_loc_in, n_loc_mid = n_loc_mid, n_loc_out = nrow(g),
       n_tr_in  = n_tr_in,  n_tr_mid  = n_tr_mid,
       n_tr_out = length(unique(g$tid)), n_tr_sigma = n_sig)
}

