## ─────────────────────────────────────────────────────────────────────────────
## sim_labeling.R  -  simulator: labelling, blinking and localizations
## ─────────────────────────────────────────────────────────────────────────────
##
## PURPOSE
##   From molecule positions to a MINFLUX localization table: label presets and linkage, fluorophore positions per molecule, blinking into trace ids, and localizations scattered with the localization precision.
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

# ── Labeling presets ─────────────────────────────────────────────────────────
# n_fluor: mean fluorophores per labeled molecule (Poisson-drawn, min 1).
# link_nm: mean distance fluorophore→true molecule centre (nm).
# In MINFLUX each blinking fluorophore = one trace ID.
# n_fluor is the *theoretical* count; actual detected TIDs are ~60–80 % of
# this due to incomplete blinking/detection efficiency.
#
# Key assumptions:
#   IgG arm (Fab tip to Fc base) ≈ 11 nm; secondary adds another ≈ 11 nm.
#   Nanobody body ≈ 2.5 nm → effective linker ~3 nm.
#   NHS-ester DOL for commercial IgG kits: 2–4 (mean ~3).
#   Monoclonal: 1 antibody per epitope (bivalent, but typically 1 per target).
#   Polyclonal: 3–5 antibodies per antigen (multiple epitopes accessible).
#   Secondary binding stoichiometry: 2–4 secondaries per primary IgG.
# dye_spread_nm: Gaussian sigma of individual dye positions within the label
# complex (around the shared complex-centre displaced by link_nm).
# n_fluor > 1 on the SAME complex → dyes cluster within dye_spread_nm;
# each molecule's complex is displaced independently by link_nm.
# Single source of truth lives in helpers_shared.R so the SMLM tab's inverse
# GT-estimator inverts exactly the numbers this forward model generates with.
LABEL_PRESETS <- minflux_label_presets()

# Expand labeled molecules to individual fluorophore positions.
# eff_pct:       labeling efficiency (% of molecules that get labeled)
# n_fluor:       mean fluorophores per labeled molecule (Poisson-drawn, min 1)
# link_nm:       mean distance from label-complex centre to molecule centre (nm)
# dye_spread_nm: Gaussian SD of individual dye positions within the complex
#
# Two-level displacement model:
#   Level 1 (per molecule): one complex-centre offset of ~link_nm, shared by
#     all fluorophores on that molecule.
#   Level 2 (per fluorophore): each dye scatters by dye_spread_nm around the
#     complex centre. For single-dye labels dye_spread_nm = 0.
# Each returned row = one fluorophore = one future trace ID.
expand_to_fluorophores <- function(gt_df, eff_pct, n_fluor,
                                    link_nm, dye_spread_nm = 0,
                                    ldir_aware = TRUE,
                                    protein_r_nm = 0) {
  use_3d <- "z" %in% names(gt_df)

  # Step 1: all GT molecules can be labeled; orig_rows tracks indices into
  # the full gt_df so mol_id can look up res$gt in the plotting code.
  orig_rows <- seq_len(nrow(gt_df))
  n_mol     <- nrow(gt_df)
  if (n_mol == 0) {
    base <- gt_df[, intersect(c("x","y","z","n_locs_expected"), names(gt_df)),
                  drop = FALSE]
    base$mol_id <- integer(0)
    return(base)
  }

  # Step 2: draw fluorophore count per molecule (Poisson, min 1)
  if (n_fluor <= 1L) {
    n_per <- rep(1L, n_mol)
  } else {
    n_per <- pmax(1L, rpois(n_mol, n_fluor))
  }

  # Step 3: expand molecule rows to fluorophore rows.
  # mol_id stores the original GT row index for use in the plotting code.
  mol_idx <- rep(seq_len(n_mol), n_per)
  keep    <- intersect(c("x","y","z","n_locs_expected"), names(gt_df))
  f_df    <- gt_df[mol_idx, keep, drop = FALSE]
  f_df$mol_id <- orig_rows[mol_idx]
  n_f     <- nrow(f_df)

  # Step 3b: per-fluorophore detection filter.
  # Labeling efficiency applies to each individual fluorophore independently -
  # this is the biologically relevant quantity for MINFLUX: what fraction of
  # the theoretical TIDs are actually detected?  Molecules where every
  # fluorophore is missed simply produce no TIDs (equivalent to being
  # unlabeled), which emerges naturally without a separate molecule-level step.
  if (eff_pct < 100 && n_f > 0) {
    detected <- runif(n_f) < (eff_pct / 100)
    f_df  <- f_df[detected, , drop = FALSE]
    n_f   <- nrow(f_df)
    if (n_f == 0) return(f_df)
    mol_idx <- mol_idx[detected]
  }

  # Local helper: draw n_mol structured displacement vectors of given magnitude.
  # Respects ldir_type: "perp" (filament) → lateral only; "radial" (ring/sphere)
  # → outward from structure centre; NA → isotropic.
  # Returns a list with elements dx, dy, and (when use_3d) dz.
  struct_displace <- function(mag_vec) {
    ldt <- if (ldir_aware && "ldir_type" %in% names(gt_df))
             gt_df$ldir_type[1L] else NA_character_

    if (!is.na(ldt) && ldt == "perp" &&
        all(c("ldir_x","ldir_y") %in% names(gt_df))) {
      if (use_3d) {
        tang <- cbind(gt_df$ldir_x, gt_df$ldir_y, gt_df$ldir_z)
        u    <- matrix(rnorm(n_mol * 3L), ncol = 3L)
        u    <- u - rowSums(u * tang) * tang
        u_n  <- sqrt(rowSums(u^2)); u_n[u_n == 0] <- 1
        u    <- u / u_n
        list(dx = mag_vec * u[,1L], dy = mag_vec * u[,2L], dz = mag_vec * u[,3L])
      } else {
        sgn <- sample(c(-1L, 1L), n_mol, replace = TRUE)
        list(dx = mag_vec * sgn * (-gt_df$ldir_y),
             dy = mag_vec * sgn * ( gt_df$ldir_x))
      }
    } else if (!is.na(ldt) && ldt == "radial") {
      cx <- if ("struct_cx" %in% names(gt_df)) gt_df$struct_cx else rep(0, n_mol)
      cy <- if ("struct_cy" %in% names(gt_df)) gt_df$struct_cy else rep(0, n_mol)
      rx <- gt_df$x - cx; ry <- gt_df$y - cy
      if (use_3d) {
        cz <- if ("struct_cz" %in% names(gt_df)) gt_df$struct_cz else rep(0, n_mol)
        rz <- gt_df$z - cz
        rr <- sqrt(rx^2 + ry^2 + rz^2); rr[rr == 0] <- 1
        list(dx = mag_vec * rx/rr, dy = mag_vec * ry/rr, dz = mag_vec * rz/rr)
      } else {
        rr <- sqrt(rx^2 + ry^2); rr[rr == 0] <- 1
        list(dx = mag_vec * rx/rr, dy = mag_vec * ry/rr)
      }
    } else {
      if (use_3d) {
        dirs <- matrix(rnorm(n_mol * 3L), ncol = 3L)
        lens <- sqrt(rowSums(dirs^2)); lens[lens == 0] <- 1
        dirs <- dirs / lens
        list(dx = mag_vec * dirs[,1L], dy = mag_vec * dirs[,2L],
             dz = mag_vec * dirs[,3L])
      } else {
        theta <- runif(n_mol, 0, 2 * pi)
        list(dx = mag_vec * cos(theta), dy = mag_vec * sin(theta))
      }
    }
  }

  # Step 3c (Level 0): epitope offset from protein centroid.
  # The antibody binds a surface epitope, not the geometric centre of the protein.
  # Displacement magnitude is exactly the protein radius (r = 0.66 × MW^(1/3) nm);
  # direction follows the same structure-aware rules as the linker arm.
  if (protein_r_nm > 0 && n_mol > 0) {
    d0 <- struct_displace(rep(protein_r_nm, n_mol))
    f_df$x <- f_df$x + d0$dx[mol_idx]
    f_df$y <- f_df$y + d0$dy[mol_idx]
    if (use_3d) f_df$z <- f_df$z + d0$dz[mol_idx]
  }

  # Step 4 (Level 1): one complex-centre displacement per MOLECULE.
  # All fluorophores on the same molecule share this offset.
  # Direction is constrained by the local structure geometry when ldir_aware=TRUE.
  if (link_nm > 0 && n_mol > 0) {
    sigma_link <- link_nm * sqrt(pi / 2)   # half-normal → mean ≈ link_nm
    d1 <- struct_displace(abs(rnorm(n_mol, 0, sigma_link)))
    f_df$x <- f_df$x + d1$dx[mol_idx]
    f_df$y <- f_df$y + d1$dy[mol_idx]
    if (use_3d) f_df$z <- f_df$z + d1$dz[mol_idx]
  }

  # Step 5 (Level 2): per-fluorophore scatter within the complex body.
  # Skipped when dye_spread_nm = 0 (single-dye labels, custom dye_spread = 0).
  if (dye_spread_nm > 0 && n_f > 0) {
    f_df$x <- f_df$x + rnorm(n_f, 0, dye_spread_nm)
    f_df$y <- f_df$y + rnorm(n_f, 0, dye_spread_nm)
    if (use_3d) f_df$z <- f_df$z + rnorm(n_f, 0, dye_spread_nm)
  }

  rownames(f_df) <- NULL
  f_df
}

# ── Blinking / TID splitting ────────────────────────────────────────────────
# Models dye re-activation: each fluorophore may appear as K additional TIDs,
# K ~ Poisson(lambda). Extra rows carry the same mol_id and position; each TID
# generates its own localizations independently via gen_localizations, giving
# apparent over-counting proportional to (1 + lambda) × n_fluor_per_molecule.
split_for_blinking <- function(fluoro_df, lambda) {
  if (lambda <= 0 || nrow(fluoro_df) == 0L) return(fluoro_df)
  extra <- rpois(nrow(fluoro_df), lambda)
  if (sum(extra) == 0L) return(fluoro_df)
  rbind(fluoro_df, fluoro_df[rep(seq_len(nrow(fluoro_df)), extra), , drop = FALSE])
}

# ── Localization generator ──────────────────────────────────────────────────
# Uses n_locs_expected column in fluoro_df when present (image-guided mode);
# otherwise draws from Poisson(trace_length) for all fluorophores.
# sigma_loc_xy / sigma_loc_z are the per-localization precision in nm.
#
# use_geometric:  replace Poisson with Geometric(p = 1/lambda).  Same mean
#   but heavy right tail: constant per-step bleaching gives many short traces
#   and a few long ones, matching real MINFLUX trace-length distributions.
#
# scale_by_fluor: divide each fluorophore's lambda by the number of
#   fluorophores on its molecule (requires mol_id column in fluoro_df).
#   Keeps total locs per molecule ≈ trace_length regardless of labeling
#   strategy - models a fixed observation budget per molecule.

gen_localizations <- function(fluoro_df, trace_length,
                               sigma_loc_xy, sigma_loc_z = sigma_loc_xy,
                               use_geometric  = FALSE,
                               scale_by_fluor = FALSE) {
  n      <- nrow(fluoro_df)
  use_3d <- "z" %in% names(fluoro_df)
  if (n == 0) {
    cols  <- if (use_3d) c("x","y","z","tid","tim")
             else        c("x","y","tid","tim")
    empty <- as.data.frame(matrix(ncol = length(cols), nrow = 0))
    names(empty) <- cols
    return(empty)
  }

  lambda <- if ("n_locs_expected" %in% names(fluoro_df))
              fluoro_df$n_locs_expected
            else
              rep(max(1, trace_length), n)

  # Divide each fluorophore's budget by how many siblings it has on its molecule
  if (scale_by_fluor && "mol_id" %in% names(fluoro_df)) {
    mol_counts <- ave(rep(1, n), fluoro_df$mol_id, FUN = sum)
    lambda     <- lambda / pmax(1, mol_counts)
  }
  lambda <- pmax(1, lambda)

  if (use_geometric) {
    # rgeom(n, p) = failures before first success (min 0); +1 → min trace = 1
    # mean = 1/p = lambda, heavy right tail from constant bleach probability
    p     <- pmin(0.999, 1 / lambda)
    n_loc <- pmax(1L, rgeom(n, p) + 1L)
  } else {
    n_loc <- pmax(1L, rpois(n, lambda))
  }

  idx    <- rep(seq_len(n), n_loc)
  total  <- length(idx)
  rep_df <- fluoro_df[idx, intersect(c("x","y","z"), names(fluoro_df)),
                      drop = FALSE]
  rep_df$x <- rep_df$x + rnorm(total, 0, sigma_loc_xy)
  rep_df$y <- rep_df$y + rnorm(total, 0, sigma_loc_xy)
  if (use_3d) {
    rep_df$z <- rep_df$z + rnorm(total, 0, sigma_loc_z)
    rep_df$z <- pmax(-1000, pmin(1000, rep_df$z))
  }
  rep_df$tid <- as.character(idx)
  rep_df$tim <- seq(0, by = 10, length.out = total)
  rownames(rep_df) <- NULL
  rep_df
}

# ── Linkage presets ──────────────────────────────────────────────────────────
LINK_PRESETS <- c(
  "Direct tag (HaloTag / SNAP / click-EdU)" = 2,
  "Nanobody"                                 = 4,
  "Primary antibody (IgG)"                  = 10,
  "Primary + secondary antibody"            = 20,
  "Custom"                                   = NA_real_
)

