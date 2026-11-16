# Tests for the computations moved out of the tracking server reactives into R/
# (est_directionality.R, est_msd.R additions, est_cve.R additions, trk_stats.R).
# Each is checked against a case with a known answer.
for (f in c("R/est_blur.R", "R/est_hmm.R", "R/est_cve.R", "R/est_msd.R", "R/trk_stats.R",
            "R/est_directionality.R", "R/est_clean.R"))
  source(file.path(.repo_root, f))

# GOAL: tracks with a chosen number of cycles per step, for timing tests.
.timed_tracks <- function(k_per_step, t_c = 0.2, n_tracks = 3L) {
  do.call(rbind, lapply(seq_len(n_tracks), function(i)
    data.frame(tid = paste0("t", i), tim = c(0, cumsum(k_per_step * t_c)),
               x = 0, y = 0, stringsAsFactors = FALSE)))
}

# CLAIM: a step is "integrated" when it spans one cycle; the fraction is counted
# over within-track steps only.
test_that("integrated fraction and continuity follow the cycle counts", {
  k  <- c(1, 1, 2, 1, 4, 1, 1, 3)                       # 5 of 8 steps span one cycle
  df <- .timed_tracks(k)
  expect_equal(cve_integrated_fraction(df, 0.2), 5 / 8)
  expect_equal(cve_continuity_mean(df, 0.2, 0), 5 / 8)  # half-life 0 = integrated fraction
  expect_equal(cve_continuity_mean(df, 0.2, 3), mean(2^(-(k - 1) / 3)))
  # the dataset mean agrees with the per-track scores of the continuity gate
  expect_equal(cve_continuity_mean(df, 0.2, 3), mean(continuity_scores(df, 0.2, 3)))
  expect_true(is.na(cve_integrated_fraction(df[0, ], 0.2)))
  expect_equal(cve_integrated_fraction(df, NA_real_), 1)  # no cycle time: every step counts
})

test_that("cycle classes are counted and capped at the display limit", {
  k <- c(rep(1, 6), 2, 2, 3, 9)
  r <- cycle_class_counts(k * 0.2, 0.2, max_disp = 5L)
  expect_equal(r$counts, c(6L, 2L, 1L, 0L, 1L))         # the 9-cycle step lands in the last bin
  expect_equal(c(r$n_total, r$n_integ, r$n_dark_int), c(10L, 6L, 4L))
  expect_equal(r$dark_cycles_total, sum(k - 1))
  expect_equal(r$frac_integ, 0.6)
})

# CLAIM: straight motion gives turning angles of 0 and a direction autocorrelation
# of 1; an immobile molecule seen through localization noise turns back on itself
# (mean angle above 90 degrees) and has a first-lag autocorrelation near -1/2.
test_that("directionality separates straight motion from pure noise", {
  t_c <- 0.2; n <- 60
  straight <- do.call(rbind, lapply(1:5, function(i)
    data.frame(tid = paste0("s", i), tim = (0:(n - 1)) * t_c, x = (0:(n - 1)) * 10, y = 3 * i)))
  ds <- dir_compute(straight, t_c, NA_real_, min_pts = 4L, max_dark = 2L)
  expect_lt(max(ds$angles), 1e-4)
  expect_equal(ds$ac[1:5], rep(1, 5))
  expect_equal(ds$lag_time[3], 3 * t_c)

  set.seed(4)
  noise <- do.call(rbind, lapply(1:40, function(i)
    data.frame(tid = paste0("n", i), tim = (0:(n - 1)) * t_c, x = rnorm(n, 0, 5), y = rnorm(n, 0, 5))))
  dn <- dir_compute(noise, t_c, 5, min_pts = 4L, max_dark = 2L)
  expect_gt(mean(dn$angles), 100)
  # the correlation of unit directions is weaker than that of the steps (-0.5), but clearly negative
  expect_lt(dn$ac[1], -0.3)
  expect_lt(abs(dn$ac[3]), 0.05)                        # no memory beyond the shared localization
  expect_equal(dn$r_noise, 10)
  expect_null(dir_compute(straight[1:5, ], t_c, NA_real_, min_pts = 4L, max_dark = 2L))
})

# CLAIM: steps across more dark cycles than allowed are left out and counted.
test_that("directionality drops steps across long dark gaps", {
  t_c <- 0.2
  tim <- cumsum(c(0, rep(c(1, 1, 1, 6), 15))) * t_c     # every fourth step spans 6 cycles
  df  <- data.frame(tid = "a", tim = tim, x = seq_along(tim) * 10, y = 0)
  d   <- dir_compute(df, t_c, NA_real_, min_pts = 4L, max_dark = 2L)
  expect_equal(d$frac_dropped_dark, 0.25)
})

# CLAIM: the MSD panel's curve and the headline-D curve are the same ensemble.
# msd_unified_curves() (the panel) and msd_ensemble() (msd_headline) were written
# separately; this keeps them from drifting apart.
test_that("the MSD panel curve equals msd_ensemble()", {
  set.seed(2)
  Q0 <- matrix(c(-1e-9, 1e-9, 1e-9, -1e-9), 2)
  df <- sim_ct_state_tracks(40, 150, D = c(200, 200.0002), Q = Q0, sigma = 4, seed = 2)
  t_c <- 0.192; sw <- 0.011
  sums <- msd_track_sums(df, t_c, sw, 60L, 20L)
  for (tw in c(TRUE, FALSE)) {
    a <- msd_unified_curves(sums, t_c, sw, 60L, track_weight = tw)
    b <- msd_ensemble(sums, t_c, sw, track_weight = tw)
    ok <- is.finite(b$msd)
    expect_equal(a$msd[ok], b$msd[ok])
    expect_false(any(is.finite(a$msd[!ok])))
    expect_equal(a$keep_lag, b$keep_lag)
    expect_equal(a$tau, b$tau); expect_equal(a$tau_lin, b$tau_lin)
    expect_equal(a$npairs, b$npairs)
  }
  expect_equal(msd_unified_curves(sums, t_c, sw, 60L, d_eff = 4 / 3)$coef_eff, 8 / 3)
})

# CLAIM: on one free population the per-track D scatters around the true D and
# the confinement index around 1.
test_that("per-track MSD fits recover D on free diffusion", {
  Q0 <- matrix(c(-1e-9, 1e-9, 1e-9, -1e-9), 2)
  df <- sim_ct_state_tracks(120, 400, D = c(200, 200.0002), Q = Q0, sigma = 3, gap_p = 1, seed = 8)
  t_c <- 0.192; sw <- 0.011
  ub <- msd_unified_curves(msd_track_sums(df, t_c, sw, 60L, 50L), t_c, sw, 60L, track_weight = TRUE)
  pt <- msd_pertrack_fits(ub, 5L, 50L, use_iter = FALSE, linear_tau = TRUE, R_msd = MSD_R_FALLBACK)
  expect_gt(pt$n_tD, 100)
  expect_lt(abs(pt$med_tD / 0.2 - 1), 0.15)             # um2/s
  expect_lt(abs(stats::median(pt$pertrack$conf_index, na.rm = TRUE) - 1), 0.2)
  expect_true(all(pt$per_D > 0))
})
