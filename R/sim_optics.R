## ─────────────────────────────────────────────────────────────────────────────
## sim_optics.R  -  simulator: confocal view and local density
## ─────────────────────────────────────────────────────────────────────────────
##
## PURPOSE
##   The diffraction-limited (confocal) rendering of a simulated sample and the local molecule density used for the crowding model.
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

# ── Confocal PSF computation from optical parameters ────────────────────────
# Gaussian approximation of the lateral confocal PSF sigma.
#
#   σ_ex   = 0.21 × λ_ex / NA               excitation PSF lateral sigma
#   σ_em   = 0.21 × λ_em / NA               emission PSF lateral sigma
#   1 AU   = 1.22 × λ_em / NA  (nm)         one Airy Unit diameter (sample plane)
#   r_ph   = AU × (1 AU_nm) / 2             pinhole radius in nm
#   σ_det  = √(σ_em² + (r_ph / 2.355)²)    detection sigma incl. pinhole size
#   σ_conf = σ_ex × σ_det / √(σ_ex² + σ_det²)   product of two Gaussians
#
# For 640 nm / NA 1.4 / 0.83 AU:  σ ≈ 80 nm, FWHM ≈ 188 nm
# For 560 nm / NA 1.4 / 1.0 AU:   σ ≈ 72 nm, FWHM ≈ 170 nm
# For 488 nm / NA 1.4 / 1.0 AU:   σ ≈ 63 nm, FWHM ≈ 148 nm
compute_confocal_sigma <- function(lambda_ex, lambda_em, na, au) {
  sigma_ex  <- 0.21 * lambda_ex / na
  sigma_em  <- 0.21 * lambda_em / na
  au_nm     <- 1.22 * lambda_em / na        # 1 AU diameter in nm (sample plane)
  ph_radius <- au * au_nm / 2               # pinhole radius in nm
  sigma_det <- sqrt(sigma_em^2 + (ph_radius / 2.355)^2)
  max(1, sigma_ex * sigma_det / sqrt(sigma_ex^2 + sigma_det^2))
}

# Fraction of a point emitter's light that passes the pinhole - the encircled
# energy of the Airy pattern within the pinhole radius. v is the optical radius
# (v = 3.8317 at exactly 1 AU, the first Airy minimum, where ~84 % is enclosed).
# This drives the real confocal trade-off: closing the pinhole (AU < 1) rejects
# light - dimmer, noisier - but tightens the PSF and sections better; opening it
# (AU > 1) barely adds signal yet blurs and lets in out-of-focus haze.
confocal_collection_efficiency <- function(au) {
  v <- 3.8317 * max(0, au)
  max(0.01, min(1, 1 - besselJ(v, 0)^2 - besselJ(v, 1)^2))
}

# ── Diffraction-limited image renderer ──────────────────────────────────────
# Renders each GT molecule as a 2D Gaussian PSF on a pixel grid, sums them,
# adds uniform background, then draws Poisson shot noise.
# Only computes the PSF within a 3-sigma radius per molecule (sparse update)
# so performance stays acceptable for up to ~5 000 molecules.
# Grid is capped to 500 × 500 px.
# px_nm_user: when supplied, that pixel size is used directly (honouring a
# user-chosen camera/scanner pixel size in nm).  If the structure extent
# would push the grid over 500 px at that pixel size, px_nm is silently
# increased just enough to stay within the cap.
render_diffraction_image <- function(x_nm, y_nm, psf_sigma_nm,
                                      photons_per_mol, bg_photons_per_px,
                                      px_nm_user = NULL,
                                      photons_vec = NULL,
                                      extent = NULL) {
  n <- length(x_nm)
  if (n == 0) return(NULL)
  pad  <- 4 * psf_sigma_nm
  if (!is.null(extent) && length(extent) == 4 && all(is.finite(extent))) {
    # Explicit window (pinned to the FOV) so the confocal view spans the same
    # ROI as the MINFLUX view instead of cropping to the data bounding box.
    x0 <- extent[1]; x1 <- extent[2]; y0 <- extent[3]; y1 <- extent[4]
  } else {
    x0 <- min(x_nm) - pad; x1 <- max(x_nm) + pad
    y0 <- min(y_nm) - pad; y1 <- max(y_nm) + pad
  }
  span <- max(x1 - x0, y1 - y0, 1)
  # Use user pixel size if given; auto-size otherwise.
  # Either way, enforce the 500 px cap.
  px_nm <- if (!is.null(px_nm_user) && px_nm_user > 0)
             max(px_nm_user, span / 500)
           else
             max(psf_sigma_nm / 2, span / 500)
  nx    <- max(10L, min(500L, ceiling((x1 - x0) / px_nm)))
  ny    <- max(10L, min(500L, ceiling((y1 - y0) / px_nm)))
  # Keep pixels square at the requested size (don't rescale to fit)
  py_nm <- px_nm
  xs <- x0 + (seq_len(nx) - 0.5) * px_nm
  ys <- y0 + (seq_len(ny) - 0.5) * py_nm
  img  <- matrix(0, nrow = ny, ncol = nx)
  r_px <- ceiling(3 * psf_sigma_nm / min(px_nm, py_nm))
  for (i in seq_len(n)) {
    xi <- floor((x_nm[i] - x0) / px_nm) + 1L
    yi <- floor((y_nm[i] - y0) / py_nm) + 1L
    ix <- max(1L, xi - r_px):min(nx, xi + r_px)
    iy <- max(1L, yi - r_px):min(ny, yi + r_px)
    if (length(ix) == 0 || length(iy) == 0) next
    gx   <- exp(-0.5 * ((xs[ix] - x_nm[i]) / psf_sigma_nm)^2)
    gy   <- exp(-0.5 * ((ys[iy] - y_nm[i]) / psf_sigma_nm)^2)
    ph_i <- if (!is.null(photons_vec)) photons_vec[i] else photons_per_mol
    img[iy, ix] <- img[iy, ix] + ph_i * outer(gy, gx)
  }
  img <- img + bg_photons_per_px
  img <- matrix(rpois(length(img), as.vector(img)), nrow = ny)
  list(img = img, xs = xs, ys = ys, px_nm = px_nm)
}

# Per-localization local density: number of neighbours within radius_nm
# (self excluded). Used by the "Local density" colour mode. k is capped so
# very dense scenes stay fast; counts saturate at k-1 in that regime.
sim_local_density <- function(locs, radius_nm) {
  cols <- intersect(c("x", "y", "z"), names(locs))
  pts  <- as.matrix(locs[, cols, drop = FALSE])
  n    <- nrow(pts)
  if (n < 2) return(rep(0, n))
  k    <- min(n, 256L)
  nn   <- RANN::nn2(pts, pts, k = k)
  rowSums(nn$nn.dists <= radius_nm) - 1L
}

