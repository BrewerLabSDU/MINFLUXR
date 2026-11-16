## ─────────────────────────────────────────────────────────────────────────────
## est_directionality.R  -  turning angles and velocity autocorrelation
## ─────────────────────────────────────────────────────────────────────────────
##
## PURPOSE
##   The computation behind the Tracking tab's Directionality panel as a plain
##   function. Moved unchanged out of dir_bundle() (server/tracking/
##   19_directionality.R); the server reads the inputs and calls it.
##
## METHOD
##   Turning angle between two consecutive usable steps, cos(a) = (d1 . d2) /
##   (|d1| |d2|). Velocity autocorrelation of unit step directions, binned by
##   ELAPSED TIME in whole acquisition cycles (not by step index, which would
##   treat unevenly spaced MINFLUX localizations as evenly spaced). A step is
##   usable when it has non-zero length and spans at most `max_dark` dark cycles.
##   With localization noise the first lag is expected to be negative
##   (about -sigma^2 / MSD); that is noise, not a finding.
##
## UNITS: nm, ms.
## ─────────────────────────────────────────────────────────────────────────────

# ── dir_compute() ─────────────────────────────────────────────────────────────
# GOAL: turning angles and the direction autocorrelation of a set of tracks.
# df       : data.frame(tid, x, y, tim), tim in ms
# t_c      : cycle time (ms); NA = unknown, lags then count steps instead of cycles
# sig_dir  : sigma_loc (nm) for the noise floor 2*sigma on a 2D step; NA = none
# min_pts  : tracks with fewer localizations are skipped
# max_dark : a usable step spans at most this many dark cycles
# max_lag  : largest autocorrelation lag (cycles, or steps without a cycle time)
# progress : optional function(fraction), called every 50 tracks
# RETURNS list(angles (degrees), ac, lags, lag_time (ms), have_cyc, t_cycle,
#   max_dark, r_noise, sigma_loc, frac_below_noise, frac_dropped_dark), or NULL
#   when fewer than 10 turning angles were found.
dir_compute <- function(df, t_c = NA_real_, sig_dir = NA_real_, min_pts = 4L, max_dark = 2L,
                        max_lag = 20L, progress = NULL) {
  df <- df[is.finite(df$x) & is.finite(df$y) & is.finite(df$tim), ]
  df <- df[order(df$tid, df$tim), ]
  sp <- split(df, df$tid)

  # Cycle time: the autocorrelation lag is measured in ELAPSED TIME, expressed as
  # integer cycle counts (like the MSD), not raw step index. On irregular MINFLUX
  # sampling a "step" is not a fixed Δt, so a step-index lag would misrepresent
  # unevenly-spaced points as evenly spaced. When the cycle time is unknown we fall
  # back to step-index lags and say so.
  have_cyc  <- is.finite(t_c) && t_c > 0

  # Pre-allocate angle storage as a list then unlist once at the end
  angle_list <- vector("list", length(sp))
  ac_sums    <- rep(0,  max_lag)
  ac_n       <- rep(0L, max_lag)
  # Localization-noise floor on a 2D step magnitude: RMS pure-noise step = 2·σ_loc.
  # Steps below it have noise-dominated directions (inflate the 180° turn-angle
  # peak and the lag-1 anti-correlation). Count them to quantify how much of the
  # directionality signal is noise. Uses the app-wide σ_loc.
  r_noise <- if (is.finite(sig_dir) && sig_dir > 0) 2 * sig_dir else NA_real_
  n_steps_tot <- 0L; n_below <- 0L; n_dropped_dark <- 0L

  for (ii in seq_along(sp)) {
    if (is.function(progress) && ii %% 50L == 0L) progress(50 / length(sp))
    tb <- sp[[ii]]
    if (nrow(tb) < min_pts) next
    dx  <- diff(tb$x); dy <- diff(tb$y)
    dtt <- diff(tb$tim)
    if (length(dx) < 2) next

    mag <- sqrt(dx^2 + dy^2)
    # Per-step dark-cycle count. A "usable" step spans ≤ max_dark dark cycles
    # (integrated = 0 dark) and has non-zero length. Steps across long gaps are
    # dropped: the direction over a big gap is unreliable and may bridge a merge
    # with another molecule.
    n_cyc  <- if (have_cyc) pmax(1L, as.integer(round(dtt / t_c))) else rep(1L, length(dtt))
    usable <- mag > 0 & is.finite(n_cyc) & ((n_cyc - 1L) <= max_dark)
    n_steps_tot    <- n_steps_tot + length(mag)
    n_dropped_dark <- n_dropped_dark + sum(mag > 0 & (n_cyc - 1L) > max_dark, na.rm = TRUE)
    if (is.finite(r_noise)) n_below <- n_below + sum(mag < r_noise & usable, na.rm = TRUE)

    # --- Turn angles: only between two CONSECUTIVE usable steps (share a loc) ---
    ok_pair <- usable[-length(usable)] & usable[-1]
    if (any(ok_pair)) {
      dot   <- dx[-length(dx)] * dx[-1] + dy[-length(dy)] * dy[-1]
      denom <- mag[-length(mag)] * mag[-1]
      cos_a <- dot[ok_pair] / denom[ok_pair]
      cos_a <- pmax(-1, pmin(1, cos_a))
      angle_list[[ii]] <- acos(cos_a) * 180 / pi
    }

    # --- Velocity autocorrelation, binned by ELAPSED-TIME (cycle-index) lag ---
    # Each usable step's reference position on the track's cumulative cycle
    # timeline is the cycle index of its start localization. Pairs of usable steps
    # are accumulated into the bin = difference of those cycle indices, so a lag of
    # k means k cycles of real time, however many localizations were skipped.
    cyc_start <- c(0L, cumsum(n_cyc))          # per-loc cumulative cycle index
    cyc_start <- cyc_start[-length(cyc_start)] # start-of-step index (length = #steps)
    ux <- ifelse(usable, dx / mag, NA_real_)
    uy <- ifelse(usable, dy / mag, NA_real_)
    idx_u <- which(usable)
    for (a in idx_u) {
      ca <- cyc_start[a]
      for (b in idx_u[idx_u > a]) {
        lag <- cyc_start[b] - ca
        if (lag < 1L) next
        if (lag > max_lag) break            # cyc_start is increasing → no closer pairs beyond
        ac_sums[lag] <- ac_sums[lag] + ux[a]*ux[b] + uy[a]*uy[b]
        ac_n[lag]    <- ac_n[lag]    + 1L
      }
    }
  }

  all_angles <- unlist(angle_list, use.names = FALSE)

  if (length(all_angles) < 10) return(NULL)        # caller reports "too few usable steps"

  ac_vals  <- ifelse(ac_n > 0, ac_sums / ac_n, NA_real_)
  lag_time <- if (have_cyc) seq_len(max_lag) * t_c else NA_real_   # ms, real elapsed time
  list(angles = all_angles, ac = ac_vals, lags = seq_len(max_lag),
       lag_time = lag_time, have_cyc = have_cyc, t_cycle = t_c, max_dark = max_dark,
       r_noise = r_noise, sigma_loc = if (is.finite(r_noise)) sig_dir else NA_real_,
       frac_below_noise = if (n_steps_tot > 0) n_below / n_steps_tot else NA_real_,
       frac_dropped_dark = if (n_steps_tot > 0) n_dropped_dark / n_steps_tot else NA_real_)
}
