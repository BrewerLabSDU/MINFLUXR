## ─────────────────────────────────────────────────────────────────────────────
## tests/validation/cpd_dark_cycle_clustering.R
## Does change-point segmentation cut a single-state track where dark cycles cluster?
## ─────────────────────────────────────────────────────────────────────────────
## cpd_segment_tracks() (R/est_cpd.R) runs PELT on log(r^2/dt). With localization
## noise E[r^2/dt] = 2dD + 2d sigma^2/dt, so the noise part depends on dt. This
## script simulates ONE free population and compares the fraction of tracks that
## get cut when (a) there are no dark cycles, (b) dark cycles are spread evenly
## along the track, (c) the same number of dark cycles sits in the last third of
## the track, as when a dye dims before bleaching. Simulated data only.
##
## Result of 2026-11-16 (200 tracks each):
##      D  sigma  noise share   no gaps   spread   late
##     20    3       0.70        0.125    0.090   0.290
##     20    6       0.90        0.105    0.130   0.590
##    215    3       0.18        0.065    0.070   0.065
##    215    6       0.47        0.060    0.045   0.120
##   (D in nm^2/ms; noise share = sigma^2 / (D t_c + sigma^2))
##   At D 20, sigma 6 the cut falls at the step where the dark cycles begin
##   (median 201 of 200) and the two sides read 32 and 0.5 nm^2/ms after the
##   noise and blur correction, for a true D of 20 throughout.
## Reading: evenly spread dark cycles do no harm. Clustered dark cycles produce a
## false "state change" only when localization noise dominates the step, i.e. for
## slow molecules; at D ~ 0.2 um2/s the effect is small.
##
## RUN  Rscript tests/validation/cpd_dark_cycle_clustering.R      (about 1 min)
## ─────────────────────────────────────────────────────────────────────────────
root <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))), "..", ".."))
setwd(root)
suppressMessages(suppressWarnings(source("global.R")))
t_c <- 0.192; sw <- 0.011; R <- 0.106
Q0 <- matrix(c(-1e-9, 1e-9, 1e-9, -1e-9), 2)
# One continuous single-state path per track, with the gap pattern switching at step 200:
# simulate finely (every cycle), then DROP localizations in the last third to create dark cycles.
make <- function(D, sigma, late = TRUE, n_tracks = 200, seed = 1) {
  full <- sim_ct_state_tracks(n_tracks, 430, D = c(D, D * 1.000001), Q = Q0, t_c = t_c, sw = 0,
                              gap_p = 1, sigma = sigma, R_blur = R, seed = seed)
  set.seed(seed + 7)
  do.call(rbind, lapply(split(full, full$tid), function(d) {
    n <- nrow(d); keep <- rep(TRUE, n)
    idx <- if (late) 201:n else 2:n
    p   <- if (late) 0.55 else 0.55 * (n - 200) / (n - 1)      # same number of dropped points overall
    keep[idx] <- stats::runif(length(idx)) > p
    d[keep, ] }))
}
stat <- function(df, sigma) {
  s <- cpd_segment_tracks(df, min_seg = 8L, noise_sig2 = sigma^2, R_abs = R * t_c)
  nst <- tapply(s$n_states, s$tid, max)
  c(cut = mean(nst > 1), mean_segments = mean(nst))
}
out <- NULL
for (D in c(20, 215)) for (sigma in c(3, 6)) {
  nog <- sim_ct_state_tracks(200, 300, D = c(D, D * 1.000001), Q = Q0, t_c = t_c, sw = 0, gap_p = 1, sigma = sigma, R_blur = R, seed = 1)
  out <- rbind(out, data.frame(D = D, sigma = sigma, noise_share = round(sigma^2 / (D * t_c + sigma^2), 2),
     cut_no_gaps = stat(nog, sigma)[1], cut_gaps_spread = stat(make(D, sigma, FALSE), sigma)[1],
     cut_gaps_late = stat(make(D, sigma, TRUE), sigma)[1]))
}
print(out, row.names = FALSE)
# where do the cuts fall, and what D do the two sides read (before noise subtraction)?
df <- make(20, 6, TRUE); s <- cpd_segment_tracks(df, min_seg = 8L, noise_sig2 = 0, R_abs = 0)
two <- s[s$n_states == 2, ]; f <- two[two$seg_idx == 1, ]; g <- two[two$seg_idx == 2, ]
cat("\nD 20, sigma 6, gaps late: tracks with exactly 2 segments:", nrow(f), "of 200\n")
cat("  first segment length, quartiles:", stats::quantile(f$n_steps, c(.25, .5, .75)), " (gaps begin after step 200)\n")
cat("  raw D first / second segment (median):", round(stats::median(f$D_nm2_ms), 1), "/", round(stats::median(g$D_nm2_ms), 1), "  true D = 20 throughout\n")
s2 <- cpd_segment_tracks(df, min_seg = 8L, noise_sig2 = 36, R_abs = R * t_c)
two <- s2[s2$n_states == 2, ]
cat("  after noise + blur correction, first / second:", round(stats::median(two$D_nm2_ms[two$seg_idx == 1]), 1), "/", round(stats::median(two$D_nm2_ms[two$seg_idx == 2]), 1), "\n")
