## ─────────────────────────────────────────────────────────────────────────────
## est_transport.R  -  per-track speed, drift velocity and local anomalous exponent
## ─────────────────────────────────────────────────────────────────────────────
##
## PURPOSE
##   The per-track quantities behind the Speed and Transport panel:
##   instantaneous step speed, rolling drift (net) velocity, the noise-corrected
##   local anomalous exponent alpha(t), and the threshold that flags an
##   export-like excursion.
##
##   Plain functions and constants with no Shiny state, moved unchanged out of
##   the tracking_server() closure (server_tracking.R). The server gathers the
##   inputs and calls these.
##
## UNITS: nm and ms; D in nm^2/ms unless a name says um2_s.
## ─────────────────────────────────────────────────────────────────────────────

# ---- Instantaneous speed (per-step v = |Δr|/Δt) -------------------------
# Per-track speed time-series + rolling median. nm/ms = µm/s directly.
trk_speed_one <- function(tb, smooth = 5L) {
  tb <- tb[order(tb$tim), ]
  t  <- as.numeric(tb$tim)
  has_z <- "z" %in% names(tb) && any(is.finite(tb$z))
  dt <- diff(t)
  dr <- sqrt(diff(tb$x)^2 + diff(tb$y)^2 + (if (has_z) diff(tb$z)^2 else 0))
  ok <- is.finite(dt) & dt > 0 & is.finite(dr)
  if (!any(ok)) return(NULL)
  spd   <- dr[ok] / dt[ok]            # nm/ms = µm/s
  t_rel <- t[-1][ok] - t[1]           # ms relative to track start
  k <- max(1L, as.integer(smooth)); if (k %% 2L == 0L) k <- k + 1L
  med <- if (length(spd) >= k && k > 1L) as.numeric(stats::runmed(spd, k)) else spd
  data.frame(t_ms = t_rel, speed_um_s = spd, speed_med = med)
}

# ---- Rolling drift (transport) velocity --------------------------------------
# v_drift(t) = |net displacement over a time window W| / W.  Unlike per-step
# path speed, this measures COHERENT transport: diffusion's net displacement
# grows only as sqrt(t), so its apparent velocity ~ sqrt(4D/W) shrinks with W,
# while directed transport gives a large, sustained v_drift. Windows are defined
# by TIME (ms), robust to MINFLUX's irregular Δt, with a min-localizations guard.
# nm/ms = µm/s directly. Returns midpoint-time, raw v_drift, and a light median.
# min_span_frac guards the track ends: near the end the forward window can no
# longer extend, so its DURATION collapses toward 0 while displacement noise
# stays - inflating v = dr/span into a spurious spike (which otherwise reads as a
# false "export at the very end"). Requiring the window to span >= frac·W drops
# those unreliable end windows. Cost: the last <frac·W of each track is not
# measured - you cannot estimate a velocity over near-zero time.
trk_drift_velocity_one <- function(tb, win_ms = 200, min_win_pts = 5L, smooth = 5L,
                                   min_span_frac = 0.5) {
  tb <- tb[order(tb$tim), ]
  has_z <- "z" %in% names(tb) && any(is.finite(tb$z))
  # Drop non-finite rows: trk() may carry NA localizations, and the net-displacement
  # endpoints are used directly, so a single NA endpoint would poison v with NA
  # (and crash the downstream if(has_export)). Matches the rest of the tracking code.
  keep <- is.finite(tb$tim) & is.finite(tb$x) & is.finite(tb$y) &
          (if (has_z) is.finite(tb$z) else TRUE)
  tb <- tb[keep, , drop = FALSE]
  t  <- as.numeric(tb$tim)
  n  <- length(t)
  min_win_pts <- max(2L, as.integer(min_win_pts))
  if (n < min_win_pts || !is.finite(win_ms) || win_ms <= 0) return(NULL)
  min_span <- max(0, min_span_frac) * win_ms
  x <- tb$x; y <- tb$y; z <- if (has_z) tb$z else rep(0, n)
  t0 <- t[1]
  out_t <- numeric(n); out_v <- numeric(n); m <- 0L
  j <- 1L
  for (i in seq_len(n)) {
    if (j < i) j <- i
    while (j < n && t[j + 1L] <= t[i] + win_ms) j <- j + 1L   # last pt within [t_i, t_i+W]
    span <- t[j] - t[i]; npts <- j - i + 1L
    if (span >= min_span && span > 0 && npts >= min_win_pts) {
      dr <- sqrt((x[j] - x[i])^2 + (y[j] - y[i])^2 + (z[j] - z[i])^2)
      m <- m + 1L
      out_t[m] <- (t[i] + t[j]) / 2 - t0   # centre of the window
      out_v[m] <- dr / span                # nm/ms = µm/s
    }
  }
  if (m < 1L) return(NULL)
  out_t <- out_t[seq_len(m)]; out_v <- out_v[seq_len(m)]
  k <- max(1L, as.integer(smooth)); if (k %% 2L == 0L) k <- k + 1L
  vmed <- if (length(out_v) >= k && k > 1L) as.numeric(stats::runmed(out_v, k)) else out_v
  data.frame(t_ms = out_t, v_drift = out_v, v_drift_med = vmed)
}

# ---- Local anomalous exponent α(t) -------------------------------------------
# Rolling log-log slope of a short-window MSD: α≈1 Brownian, α→2 ballistic/
# directed. This is the per-track analog of the ensemble TAMSD "shooting upward".
# Windowed by STEP COUNT (needs contiguous lag structure), but the log-log fit uses
# the REAL mean elapsed time per step-lag (τ = mean Δt over that lag), NOT the step
# index - so MINFLUX's irregular Δt / dark cycles don't bias α (a window straddling a
# dark gap would otherwise stretch a "1-step" lag over many cycles and flatten the
# slope). Assigned at the window-centre time. Returns NULL if a track is too short.
# sigma_nm > 0 subtracts the static localization-noise offset (2·d·σ²) from each
# MSD lag before the log-log fit. Without it, noise flattens the short-lag MSD and
# biases α downward (directed motion can read as α≈1 instead of ≈2).
trk_local_alpha_one <- function(tb, win_steps = 21L, max_lag = 6L, sigma_nm = 0) {
  tb <- tb[order(tb$tim), ]
  has_z <- "z" %in% names(tb) && any(is.finite(tb$z))
  keep <- is.finite(tb$tim) & is.finite(tb$x) & is.finite(tb$y) &
          (if (has_z) is.finite(tb$z) else TRUE)
  tb <- tb[keep, , drop = FALSE]
  n  <- nrow(tb)
  win_steps <- max(7L, as.integer(win_steps)); if (win_steps %% 2L == 0L) win_steps <- win_steps + 1L
  if (n < win_steps) return(NULL)
  x <- tb$x; y <- tb$y; z <- if (has_z) tb$z else rep(0, n)
  t <- as.numeric(tb$tim); t0 <- t[1]
  d_dim     <- if (has_z) 3 else 2
  noise_off <- if (is.finite(sigma_nm) && sigma_nm > 0) 2 * d_dim * sigma_nm^2 else 0
  half <- win_steps %/% 2L
  L    <- max(2L, min(as.integer(max_lag), win_steps %/% 3L))
  centres <- (half + 1L):(n - half)
  out_t <- numeric(length(centres)); out_a <- rep(NA_real_, length(centres)); li <- 0L
  lags  <- seq_len(L)
  for (cc in centres) {
    idx <- (cc - half):(cc + half)
    xi <- x[idx]; yi <- y[idx]; zi <- z[idx]; ti <- t[idx]; mlen <- length(idx)
    msd <- rep(NA_real_, L); tau <- rep(NA_real_, L)   # MSD and REAL mean time-lag per step-lag
    for (lg in lags) {
      if (lg >= mlen) next
      dx  <- xi[(lg + 1L):mlen] - xi[1:(mlen - lg)]
      dy  <- yi[(lg + 1L):mlen] - yi[1:(mlen - lg)]
      dz  <- zi[(lg + 1L):mlen] - zi[1:(mlen - lg)]
      dtt <- ti[(lg + 1L):mlen] - ti[1:(mlen - lg)]    # actual elapsed time of each lag-lg pair
      v   <- mean(dx^2 + dy^2 + dz^2) - noise_off
      tt  <- mean(dtt[is.finite(dtt) & dtt > 0])       # mean real time spanned at this step-lag
      if (is.finite(v) && v > 0 && is.finite(tt) && tt > 0) { msd[lg] <- v; tau[lg] <- tt }
    }
    okl <- is.finite(msd) & is.finite(tau)
    li <- li + 1L
    out_t[li] <- t[cc] - t0
    if (sum(okl) >= 2L) {   # slope of log(MSD) vs log(REAL τ) = time-based α
      fit <- suppressWarnings(stats::lm.fit(cbind(1, log(tau[okl])), log(msd[okl])))
      out_a[li] <- unname(fit$coefficients[2])
    }
  }
  ok <- is.finite(out_a)
  if (!any(ok)) return(NULL)
  data.frame(t_ms = out_t[ok], alpha = out_a[ok])
}

# Effective threshold for calling a drift-velocity window "directed":
#   absolute onset_thr if the user set one, else fold × the track's own diffusive
#   baseline (median drift velocity ~ sqrt(4D/W)), always floored at the
#   localization-noise floor √2σ/W. Fold-over-baseline (not median+k·MAD of the
#   signal) is what keeps pure-diffusion tracks from being flagged: validated on
#   synthetic data, fold≈3 gives ~5% false positives with full sensitivity.
trk_export_thr <- function(v, onset_thr, fold, vfloor) {
  base <- stats::median(v, na.rm = TRUE)
  thr  <- if (is.finite(onset_thr) && onset_thr > 0) onset_thr else fold * base
  max(thr, vfloor, na.rm = TRUE)
}

