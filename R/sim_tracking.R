## ─────────────────────────────────────────────────────────────────────────────
## sim_tracking.R  -  simulator: tracking-simulation building blocks
## ─────────────────────────────────────────────────────────────────────────────
##
## PURPOSE
##   Pieces of the tracking simulator that hold no session state. For now only the fractional Gaussian noise generator; the simulation itself is still in sim_trk_result (sim_server.R).
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

# ── .sim_fgn() - fractional Gaussian noise (Davies-Harte) ───────────────────
# GOAL: increments for fractional Brownian motion, the standard model of
# ANOMALOUS diffusion. Successive increments are correlated (anti-correlated for
# H < 0.5), which makes MSD grow as t^alpha with alpha = 2H instead of t^1.
# A cell-wall protein in a crowded peptidoglycan mesh is subdiffusive: the the tracking test data
# files measure alpha ~ 0.78.
# n: number of increments. H: Hurst exponent (alpha/2), 0 < H < 1.
# returns: numeric(n), unit variance per increment, correlated.
# Exact via circulant embedding; falls back to white noise at H = 0.5 (Brownian).
.sim_fgn <- function(n, H) {
  if (n < 2L) return(stats::rnorm(n))
  if (abs(H - 0.5) < 1e-6) return(stats::rnorm(n))
  k <- 0:(n - 1)
  g <- 0.5 * (abs(k + 1)^(2*H) - 2 * abs(k)^(2*H) + abs(k - 1)^(2*H))  # fGn autocovariance
  cc  <- c(g, g[n:2])                       # circulant first row, length m = 2n-2
  m   <- length(cc)
  lam <- Re(stats::fft(cc))
  if (any(lam < 0)) lam <- pmax(lam, 0)     # embedding not non-negative: clip (rare)
  half <- m / 2
  w <- complex(length.out = m)
  w[1]        <- sqrt(lam[1] / m) * stats::rnorm(1)
  w[half + 1] <- sqrt(lam[half + 1] / m) * stats::rnorm(1)
  j  <- 2:half
  zr <- stats::rnorm(length(j)); zi <- stats::rnorm(length(j))
  w[j]         <- sqrt(lam[j] / (2 * m)) * complex(real = zr, imaginary = zi)
  w[m - j + 2] <- Conj(w[j])
  Re(stats::fft(w))[1:n]
}

