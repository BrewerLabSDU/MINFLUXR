## ── trk_dmax.R  -  the trackable-diffusion ceiling D_max ──
##
## PURPOSE
##   dmax_estimate(), behind the Tracking tab's D_max panel. A plain function with
##   no Shiny state, moved unchanged out of the tracking_server() closure.
##
## UNITS: nm and ms in; D_max in um^2/s out.

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
