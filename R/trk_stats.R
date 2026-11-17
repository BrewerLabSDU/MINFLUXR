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

# ── cycle_class_counts() ──────────────────────────────────────────────────────
# GOAL: classify every within-track time step by the number of acquisition cycles
# it spans (the time-signature analysis). n_cyc = round(dt / t_c), at least 1;
# 1 = "integrated" (a localization every cycle), k >= 2 = k - 1 dark cycles.
# vec      : positive, finite step durations (ms), all tracks pooled
# t_c      : cycle time (ms)
# max_disp : steps spanning more cycles are counted in the last bin
# RETURNS list(t_cycle, n_cyc, counts (per cycle count, length max_disp),
#   max_disp, n_total, n_integ, n_dark_int, dark_cycles_total, frac_integ).
cycle_class_counts <- function(vec, t_c, max_disp = 20L) {
  n_cyc      <- pmax(1L, as.integer(round(vec / t_c)))
  n_total    <- length(n_cyc)
  n_integ    <- sum(n_cyc == 1L)
  n_dark_int <- sum(n_cyc >= 2L)           # intervals containing ≥1 dark cycle
  dark_cycles_total <- sum(n_cyc - 1L)     # total dark cycles across all intervals
  frac_integ <- n_integ / n_total
  tab      <- tabulate(pmin(n_cyc, max_disp), nbins = max_disp)
  list(t_cycle = t_c, n_cyc = n_cyc, counts = tab, max_disp = max_disp,
       n_total = n_total, n_integ = n_integ, n_dark_int = n_dark_int,
       dark_cycles_total = dark_cycles_total, frac_integ = frac_integ)
}

# ── sigma_used_text() ─────────────────────────────────────────────────────────
# GOAL: one sentence saying which localization precision a panel used and where
# it came from, so a result that depends on sigma (JDD, D-map, change-point, HMM,
# the confinement and displacement nulls, directionality, speed) shows it. The
# sigma is the value the panel read WHEN IT RAN, not the current one.
# info : list(sigma (nm), source, role, used, note, a, a_calibrated, profile) as
#        built by .sigma_stamp() in server/tracking/01_state_sigma.R
#   role         "subtracted" (a noise term is removed) or "noise floor" (used
#                only to flag steps or set a threshold)
#   used         FALSE when the panel's own switch turned the correction off
#   a, a_calibrated, profile : the photon-budget prefactor, whether it was set
#                (typed or from a calibration profile) and the profile's name
# RETURNS a character string, or NULL when info is NULL.
sigma_used_text <- function(info) {
  if (is.null(info)) return(NULL)
  role <- info$role %||% "subtracted"
  if (!isTRUE(info$used))
    return(sprintf("\u03c3_loc not %s: %s.", if (role == "subtracted") "subtracted" else "used",
                   info$note %||% "switched off in this panel"))
  if (!is.finite(info$sigma) || info$sigma <= 0)
    return(sprintf("\u03c3_loc not %s: no value is set. Run Data Quality & Calibration or enter \u03c3_loc.",
                   if (role == "subtracted") "subtracted" else "used"))
  src <- info$source %||% "not set"
  detail <- if (identical(src, "photon budget")) {
    if (isTRUE(info$a_calibrated))
      sprintf("photon budget, prefactor a = %.2f%s", info$a,
              if (!is.null(info$profile) && nzchar(info$profile)) sprintf(" from profile '%s'", info$profile) else "")
    else sprintf("photon budget, prefactor a = %.2f, uncalibrated (ideal estimator; calibrate on an immobile sample)", info$a)
  } else src
  sprintf("\u03c3_loc %s: %.2f nm (%s)%s.", if (role == "subtracted") "subtracted" else "used as the noise floor",
          info$sigma, detail, if (!is.null(info$note) && nzchar(info$note)) paste0("; ", info$note) else "")
}
