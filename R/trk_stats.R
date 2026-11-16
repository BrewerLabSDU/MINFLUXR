## ─────────────────────────────────────────────────────────────────────────────
## trk_stats.R  -  tracking: small statistics shared by the panels
## ─────────────────────────────────────────────────────────────────────────────
##
## PURPOSE
##   Track-level mean with confidence interval, bootstrap interval over tracks,
##   the test for real axial data, the positive-slope test used to drop
##   non-diffusing tracks, and per-step cycle statistics.
##
##   Plain functions and constants with no Shiny state, moved unchanged out of
##   the tracking_server() closure (server_tracking.R). The server gathers the
##   inputs and calls these.
##
## UNITS: nm and ms; D in nm^2/ms unless a name says um2_s.
## ─────────────────────────────────────────────────────────────────────────────

# Mean + SD + t-based 95% CI treating each element (a TRACK) as the independent
# unit - the honest CI for track-level statistics, since localizations within a
# track are correlated (the track, not the localization, is the sampling unit).
# Also returns median + quartiles for a robust summary.
.mean_ci <- function(x, conf = 0.95) {
  x <- x[is.finite(x)]; n <- length(x)
  if (n < 1L) return(list(mean = NA_real_, sd = NA_real_, se = NA_real_,
                          lo = NA_real_, hi = NA_real_, n = 0L,
                          median = NA_real_, q1 = NA_real_, q3 = NA_real_))
  m  <- mean(x); md <- stats::median(x)
  qs <- stats::quantile(x, c(0.25, 0.75), names = FALSE, na.rm = TRUE)
  if (n < 2L) return(list(mean = m, sd = NA_real_, se = NA_real_, lo = NA_real_,
                          hi = NA_real_, n = n, median = md, q1 = qs[1], q3 = qs[2]))
  s  <- stats::sd(x); se <- s / sqrt(n)
  tc <- stats::qt(1 - (1 - conf) / 2, df = n - 1)
  list(mean = m, sd = s, se = se, lo = m - tc * se, hi = m + tc * se,
       n = n, median = md, q1 = qs[1], q3 = qs[2])
}

# Percentile bootstrap CI of a statistic over the sampling units (TRACKS): x is a
# per-track vector, resampled with replacement B times. This is the honest way to
# put a CI on an ensemble number, because localizations within a track are
# correlated - the track is the independent unit, not the localization.
.boot_ci <- function(x, statfn = stats::median, B = 2000L, conf = 0.95) {
  x <- x[is.finite(x)]; n <- length(x)
  if (n < 3L) return(c(lo = NA_real_, hi = NA_real_))
  bs <- vapply(seq_len(B), function(i) statfn(x[sample.int(n, n, replace = TRUE)]), numeric(1))
  q  <- stats::quantile(bs, c((1 - conf)/2, 1 - (1 - conf)/2), names = FALSE, na.rm = TRUE)
  c(lo = q[1], hi = q[2])
}

# Does the z column carry REAL axial data, or is it a placeholder?
# 2D MINFLUX exports routinely ship an all-zero (or constant) z column, and 0 is
# finite - so the obvious test `any(is.finite(df$z))` returns TRUE for 2D data and
# every dimension-aware estimator downstream silently switches to 3D. That pools a
# ZERO-VARIANCE axis into the maths: it diluted the pooled σ_loc by √(2/3) (−18%)
# and D_cve by 2/3 (−33%), and it would have the JDD fit a 3D Maxwell law to 2D
# distances. Requiring real spread is the criterion msd_unified_bundle and
# trk_disp_bundle already used; this makes it the single shared definition so the
# dimension decision cannot drift between tools again.
# Tracking's 3D test, kept as a named wrapper because it has many call sites.
# Delegates to the shared predicate so the axial-calibration panel and the
# 2D/3D view toggle can never disagree about whether a dataset has z.
.has_real_z <- function(df, min_sd = 1) {
  has_axial_data(df, min_sd = min_sd, min_n = 11L)
}

# Per-track diffusion-sign filter (Bela's TAMSD D>0 rule). Given a per-track MSD
# matrix (rows = lags/bins, cols = tracks) and the lag-time vector, fit each
# track's MSD ~ a + b·tau by OLS over its finite points and return a logical
# mask of columns whose slope b > 0 (positive diffusion). Tracks with < 2 finite
# points are kept (cannot judge). Used to drop non-diffusing / artefact tracks
# before the ensemble in ALL three MSD pipelines.
.track_positive_D <- function(MSD_mat, tau) {
  if (is.null(dim(MSD_mat))) MSD_mat <- matrix(MSD_mat, ncol = 1L)
  apply(MSD_mat, 2L, function(m) {
    ok <- is.finite(m) & is.finite(tau)
    if (sum(ok) < 2L) return(TRUE)
    x <- tau[ok]; y <- m[ok]
    vx <- sum((x - mean(x))^2)
    if (vx <= 0) return(TRUE)
    b <- sum((x - mean(x)) * (y - mean(y))) / vx
    is.finite(b) && b > 0
  })
}

# Cycle stats on an arbitrary df (Integrated/Dark) at a given cycle time.
compute_cycle_stats <- function(df, t_c) {
  if (is.null(df) || !nrow(df) || !all(c("tid", "tim") %in% names(df)) ||
      !is.finite(t_c) || t_c <= 0) return(NULL)
  tf  <- match(df$tid, unique(df$tid))           # integer codes → radix sort (fast)
  ord <- order(tf, df$tim, method = "radix")
  tf  <- tf[ord]; tim <- df$tim[ord]; n <- length(tf)
  same <- c(FALSE, tf[-1] == tf[-n])             # consecutive pair within same track
  dts  <- c(NA_real_, diff(tim))[same]           # within-track gaps only
  dts <- dts[is.finite(dts) & dts > 0]
  if (!length(dts)) return(NULL)
  n_cyc <- pmax(1L, as.integer(round(dts / t_c)))
  list(n_tracks   = length(unique(tf)),
       n_total    = length(n_cyc),
       frac_integ = mean(n_cyc == 1L),
       dark_total = sum(n_cyc - 1L),
       mean_cyc   = mean(n_cyc),
       max_cyc    = max(n_cyc))
}

