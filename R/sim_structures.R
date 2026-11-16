## ─────────────────────────────────────────────────────────────────────────────
## sim_structures.R  -  simulator: ground-truth structure generators
## ─────────────────────────────────────────────────────────────────────────────
##
## PURPOSE
##   Molecule positions for each structure the Simulator tab offers (cluster, ring, ellipse, filament, hollow sphere, double ring, arc, stacked oligomer rings) and positions drawn from a user image.
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

# ── Shape generators ────────────────────────────────────────────────────────

gen_cluster <- function(n, radius, distribution = "uniform", use_3d = FALSE) {
  sigma <- radius / 2
  if (use_3d) {
    if (distribution == "gaussian") {
      df <- data.frame(x = rnorm(n, 0, sigma),
                       y = rnorm(n, 0, sigma),
                       z = rnorm(n, 0, sigma))
    } else {
      dirs <- matrix(rnorm(n * 3), ncol = 3)
      lens <- sqrt(rowSums(dirs^2)); lens[lens == 0] <- 1
      dirs <- dirs / lens
      r3   <- radius * runif(n)^(1/3)
      df <- data.frame(x = r3 * dirs[,1], y = r3 * dirs[,2], z = r3 * dirs[,3])
    }
  } else {
    r     <- if (distribution == "gaussian") abs(rnorm(n, 0, sigma))
             else radius * sqrt(runif(n))
    theta <- runif(n, 0, 2 * pi)
    df <- data.frame(x = r * cos(theta), y = r * sin(theta))
  }
  df$ldir_type <- NA_character_
  df
}

gen_ring <- function(n, diameter, ring_width, z_thickness = 0,
                      use_3d = FALSE) {
  radius <- diameter / 2
  theta  <- runif(n, 0, 2 * pi)
  r      <- rnorm(n, mean = radius, sd = max(0.1, ring_width / 4))
  df     <- data.frame(x = r * cos(theta), y = r * sin(theta))
  if (use_3d) df$z <- rnorm(n, 0, max(0.1, z_thickness / 4))
  df$ldir_type <- "radial"
  df$struct_cx <- 0; df$struct_cy <- 0
  if (use_3d) df$struct_cz <- 0
  df
}

gen_ellipse <- function(n, semi_a, semi_b, semi_c = 1,
                         filled = TRUE, use_3d = FALSE) {
  if (use_3d) {
    phi   <- acos(2 * runif(n) - 1)
    theta <- runif(n, 0, 2 * pi)
    r     <- if (filled) runif(n)^(1/3) else rep(1, n)
    df <- data.frame(x = semi_a * r * sin(phi) * cos(theta),
                     y = semi_b * r * sin(phi) * sin(theta),
                     z = semi_c * r * cos(phi))
  } else {
    theta <- runif(n, 0, 2 * pi)
    if (filled) {
      r  <- sqrt(runif(n))
      df <- data.frame(x = semi_a * r * cos(theta), y = semi_b * r * sin(theta))
    } else {
      df <- data.frame(x = semi_a * cos(theta), y = semi_b * sin(theta))
    }
  }
  df$ldir_type <- NA_character_
  df
}

# n molecules spaced at spacing_nm along a straight filament.
# jitter_lat_nm: noise perpendicular to the filament axis (isotropic in perp plane).
# jitter_ax_nm:  noise along the filament axis (arc-length perturbation).
# helix_r / helix_pitch: optional helical arrangement (0 = disabled);
#   default pitch 36 nm mirrors B-DNA major-groove repeat.
gen_filament <- function(n, spacing_nm, angle_deg = 0,
                          jitter_lat_nm = 0, jitter_ax_nm = 0,
                          helix_r = 0, helix_pitch = 36,
                          use_3d = FALSE, z_angle_deg = 0) {
  n      <- max(2L, n)
  s      <- (seq_len(n) - (n + 1) / 2) * spacing_nm   # centred arc-lengths
  if (jitter_ax_nm > 0) s <- s + rnorm(n, 0, jitter_ax_nm)

  xy_rad <- angle_deg  * pi / 180
  z_rad  <- z_angle_deg * pi / 180

  # Filament axis unit vector
  ax <- cos(xy_rad) * cos(z_rad)
  ay <- sin(xy_rad) * cos(z_rad)
  az <- sin(z_rad)

  # Two perpendicular unit vectors (analytical, both already unit length):
  # p1 = lateral perpendicular in XY plane; p2 = cross(axis, p1) → toward Z
  p1x <- -sin(xy_rad); p1y <- cos(xy_rad); p1z <- 0
  p2x <- -sin(z_rad) * cos(xy_rad)
  p2y <- -sin(z_rad) * sin(xy_rad)
  p2z <-  cos(z_rad)

  # Lateral jitter as a random walk along the filament: each molecule's
  # perpendicular offset = previous offset + small increment.  This makes
  # the filament curve and meander rather than scatter independently.
  # jitter_lat_nm is the per-step SD (nm); larger values = more curling.
  if (jitter_lat_nm > 0) {
    jit1 <- cumsum(rnorm(n, 0, jitter_lat_nm))
    jit2 <- if (use_3d) cumsum(rnorm(n, 0, jitter_lat_nm)) else rep(0, n)
  } else {
    jit1 <- rep(0, n)
    jit2 <- rep(0, n)
  }

  # Helical displacement: cosine in p1, sine in p2 → right-handed helix
  if (helix_r > 0 && helix_pitch > 0) {
    phi  <- 2 * pi * s / helix_pitch
    hel1 <- helix_r * cos(phi)
    hel2 <- helix_r * sin(phi)
  } else {
    hel1 <- rep(0, n)
    hel2 <- rep(0, n)
  }

  d1 <- jit1 + hel1   # total displacement in p1 direction
  d2 <- jit2 + hel2   # total displacement in p2 direction (0 in 2D XY)

  df <- data.frame(
    x = s * ax + d1 * p1x + d2 * p2x,
    y = s * ay + d1 * p1y + d2 * p2y
  )
  if (use_3d) df$z <- s * az + d1 * p1z + d2 * p2z
  # Store filament axis tangent for direction-aware linker displacement.
  # All molecules share the same global axis direction (straight filament).
  df$ldir_x    <- ax
  df$ldir_y    <- ay
  df$ldir_z    <- if (use_3d) az else 0
  df$ldir_type <- "perp"
  df
}

# Hollow sphere - uniform surface sampling with adjustable shell thickness.
# In 2D collapses to a ring (circle cross-section).
gen_hollow_sphere <- function(n, radius, shell_thick, use_3d = FALSE) {
  if (!use_3d) {
    r     <- rnorm(n, mean = radius, sd = max(0.1, shell_thick / 4))
    theta <- runif(n, 0, 2 * pi)
    df <- data.frame(x = r * cos(theta), y = r * sin(theta))
  } else {
    phi   <- acos(2 * runif(n) - 1)
    theta <- runif(n, 0, 2 * pi)
    r     <- radius + runif(n, -shell_thick / 2, shell_thick / 2)
    df <- data.frame(x = r * sin(phi) * cos(theta),
                     y = r * sin(phi) * sin(theta),
                     z = r * cos(phi))
  }
  df$ldir_type <- "radial"
  df$struct_cx <- 0; df$struct_cy <- 0
  if (use_3d) df$struct_cz <- 0
  df
}

# Double ring - two parallel rings separated along Z (or offset along Y in 2D).
# Useful for nuclear pore complex (~35 nm separation), centrioles, etc.
gen_double_ring <- function(n, diameter, ring_width, separation_nm,
                             z_thickness = 0, use_3d = FALSE) {
  radius <- diameter / 2
  n1 <- round(n / 2); n2 <- n - n1
  make_one <- function(n_pts, z_pos) {
    theta <- runif(n_pts, 0, 2 * pi)
    r     <- rnorm(n_pts, mean = radius, sd = max(0.1, ring_width / 4))
    df    <- data.frame(x = r * cos(theta), y = r * sin(theta))
    if (use_3d)
      df$z <- rnorm(n_pts, z_pos, max(0.1, z_thickness / 4))
    df
  }
  df <- rbind(make_one(n1, -separation_nm / 2),
              make_one(n2,  separation_nm / 2))
  df$ldir_type <- "radial"
  df$struct_cx <- 0; df$struct_cy <- 0
  if (use_3d) df$struct_cz <- 0
  df
}

# Arc - a partial ring spanning arc_extent_deg degrees.
# arc_start_deg sets the angular position of the start of the arc.
gen_arc <- function(n, diameter, ring_width, arc_start_deg, arc_extent_deg,
                     z_thickness = 0, use_3d = FALSE) {
  radius     <- diameter / 2
  start_rad  <- arc_start_deg  * pi / 180
  extent_rad <- max(1, arc_extent_deg) * pi / 180
  theta <- runif(n, start_rad, start_rad + extent_rad)
  r     <- rnorm(n, mean = radius, sd = max(0.1, ring_width / 4))
  df    <- data.frame(x = r * cos(theta), y = r * sin(theta))
  if (use_3d)
    df$z <- rnorm(n, 0, max(0.1, z_thickness / 4))
  df$ldir_type <- "radial"
  df$struct_cx <- 0; df$struct_cy <- 0
  if (use_3d) df$struct_cz <- 0
  df
}

# Oligomeric ring stack - fixed, evenly-spaced subunit positions (no randomness).
# Models protein complexes: hexamers, heptamers,
# nuclear pore subcomplexes, etc.
# GT positions are rigidly determined: exactly n_subunits × n_stacks molecules.
# angle_offset_deg: each successive ring is rotated by this amount (e.g., for
# symmetry mismatch between stacked assemblies).
# In 3D, rings are separated along Z by separation_nm.
# In 2D, all rings project onto XY (stacks overlap); use angle_offset_deg to
# make the rotational offset visible.
gen_oligoring_stack <- function(n_subunits, diameter, n_stacks,
                                 separation_nm, angle_offset_deg = 0,
                                 base_angle_deg = 0,
                                 use_3d = FALSE) {
  n_subunits <- max(2L, as.integer(n_subunits))
  n_stacks   <- max(1L, as.integer(n_stacks))
  radius     <- diameter / 2
  off_rad    <- angle_offset_deg * pi / 180
  base_rad   <- base_angle_deg   * pi / 180

  layers <- lapply(seq_len(n_stacks), function(k) {
    z_pos  <- (k - (n_stacks + 1) / 2) * separation_nm
    rot    <- base_rad + (k - 1L) * off_rad
    angles <- rot + 2 * pi * (seq_len(n_subunits) - 1L) / n_subunits
    df     <- data.frame(x = radius * cos(angles),
                         y = radius * sin(angles))
    if (use_3d) df$z <- z_pos
    df
  })
  df <- do.call(rbind, layers)
  df$ldir_type <- "radial"
  df$struct_cx <- 0; df$struct_cy <- 0
  if (use_3d) df$struct_cz <- 0
  rownames(df) <- NULL
  df
}

# ── Image-guided helpers ─────────────────────────────────────────────────────

read_image_gray <- function(path) {
  img   <- magick::image_read(path)
  img_g <- magick::image_convert(img, colorspace = "Gray")
  arr   <- magick::image_data(img_g, "gray")
  # dim(arr) = [1, width, height]  →  [H x W] matrix
  W <- dim(arr)[2]; H <- dim(arr)[3]
  mat <- t(matrix(as.integer(arr[1,,]), nrow = W, ncol = H)) / 255
  rng <- range(mat)
  if (diff(rng) > .Machine$double.eps)
    mat <- (mat - rng[1]) / diff(rng)
  list(mat = mat, img = img_g)
}

img_to_b64 <- function(img_magick, max_px = 800L) {
  info <- magick::image_info(img_magick)
  if (max(info$width, info$height) > max_px)
    img_magick <- magick::image_scale(
      img_magick, paste0(max_px, "x", max_px))
  raw_bytes <- magick::image_write(img_magick, format = "png")
  paste0("data:image/png;base64,", base64enc::base64encode(raw_bytes))
}

# Sample proteins from confocal image intensity.
# XY is always drawn from pixel intensity (proportional to brightness).
# In 3D mode, Z spread scales with local pixel intensity: bright regions
# receive more axial depth (z_spread_nm), dim regions stay near Z = 0.
gen_from_image <- function(img_gray, n_proteins, px_size_nm,
                            threshold, min_trace, max_trace,
                            use_3d = FALSE,
                            z_spread_nm = 200, z_bg_nm = 20) {
  m     <- img_gray; m[m < threshold] <- 0
  probs <- as.vector(m)
  if (sum(probs) == 0)
    stop("All pixels are below the threshold - lower the threshold value.")
  probs <- probs / sum(probs)
  H <- nrow(m); W <- ncol(m)

  idx     <- sample(H * W, n_proteins, replace = TRUE, prob = probs)
  row_idx <- ((idx - 1L) %% H) + 1L
  col_idx <- ((idx - 1L) %/% H) + 1L
  x_nm    <- (col_idx - 0.5 + runif(n_proteins, -0.5, 0.5)) * px_size_nm -
             (W / 2) * px_size_nm
  y_nm    <- ((H - row_idx + 0.5) + runif(n_proteins, -0.5, 0.5)) * px_size_nm -
             (H / 2) * px_size_nm
  local_int  <- probs[idx] / max(probs, na.rm = TRUE)   # 0–1
  n_locs_exp <- min_trace + (max_trace - min_trace) * local_int

  if (!use_3d) {
    df <- data.frame(x = x_nm, y = y_nm, n_locs_expected = n_locs_exp)
  } else {
    # 3D: Z spread = pmax(z_bg_nm, z_spread_nm × local_intensity)
    z_sd  <- pmax(max(0.5, z_bg_nm), z_spread_nm * local_int)
    z_nm  <- rnorm(n_proteins, 0, z_sd)
    z_nm  <- pmax(-1000, pmin(1000, z_nm))
    df <- data.frame(x = x_nm, y = y_nm, z = z_nm, n_locs_expected = n_locs_exp)
  }
  df$ldir_type <- NA_character_
  df
}

