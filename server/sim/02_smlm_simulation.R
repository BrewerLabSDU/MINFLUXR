# ╔══════════════════════════════════════════════════════════════════════════╗
# ║  PHASE 1 - Synthetic localizations                                      ║
# ╚══════════════════════════════════════════════════════════════════════════╝

# ── Phase 1 simulation (button-gated) ───────────────────────────────────────
sim_result <- eventReactive(input$sim_run, {
  if (isTRUE(input$sim_fix_seed)) {
    seed <- suppressWarnings(as.integer(input$sim_seed))
    if (!is.finite(seed)) seed <- 42L
  } else {
    seed <- sample.int(.Machine$integer.max, 1L)
  }
  set.seed(seed)

  if (isTRUE(input$sim_fix_label_seed)) {
    label_seed <- suppressWarnings(as.integer(input$sim_label_seed))
    if (!is.finite(label_seed)) label_seed <- 1L
  } else {
    label_seed <- sample.int(.Machine$integer.max, 1L)
  }

  use_3d <- identical(input$sim_view_mode, "3d")
  n_prot <- max(1L, as.integer(input$sim_n_proteins))

  # Acquisition mode adds the FOV, multi-structure placement, timing model and
  # rate-based background; Normal mode is the simple structure + labeling view.
  acq_mode <- identical(input$sim_mode %||% "normal", "acquisition")
  # FOV active → background spreads across the whole FOV (see sim_use_fov()).
  use_fov  <- sim_use_fov()
  fov_x    <- max(1, as.numeric(input$sim_fov_x %||% 4000))
  fov_y    <- max(1, as.numeric(input$sim_fov_y %||% 4000))
  fov_z    <- max(1, as.numeric(input$sim_fov_z %||% 1000))

  img_b64  <- NULL
  img_w_nm <- NULL
  img_h_nm <- NULL

  if (identical(input$sim_shape, "image")) {
    shiny::validate(shiny::need(
      !is.null(input$sim_img_file), "Please upload a confocal image."))
    img_read <- read_image_gray(input$sim_img_file$datapath)
    img_gray <- img_read$mat
    img_b64  <- img_to_b64(img_read$img)
    px_nm    <- max(1, as.numeric(input$sim_img_px_size))
    img_w_nm <- ncol(img_gray) * px_nm
    img_h_nm <- nrow(img_gray) * px_nm
    gt <- gen_from_image(
      img_gray    = img_gray,
      n_proteins  = n_prot,
      px_size_nm  = px_nm,
      threshold   = suppressWarnings(as.numeric(input$sim_img_threshold)),
      min_trace   = max(1L, as.integer(input$sim_img_min_trace)),
      max_trace   = max(1L, as.integer(input$sim_img_max_trace)),
      use_3d      = use_3d,
      z_spread_nm = max(0, as.numeric(input$sim_img_z_spread %||% 200)),
      z_bg_nm     = max(0, as.numeric(input$sim_img_z_bg     %||% 20))
    )
  } else {
    gt <- switch(input$sim_shape,
      cluster = gen_cluster(
        n = n_prot,
        radius       = max(1, as.numeric(input$sim_cluster_radius)),
        distribution = input$sim_cluster_dist,
        use_3d       = use_3d
      ),
      ring = gen_ring(
        n = n_prot,
        diameter    = max(1, as.numeric(input$sim_ring_diam)),
        ring_width  = max(0, as.numeric(input$sim_ring_width)),
        z_thickness = if (use_3d)
                        max(0, as.numeric(input$sim_ring_z_thick))
                      else 0,
        use_3d = use_3d
      ),
      ellipse = gen_ellipse(
        n = n_prot,
        semi_a = max(1, as.numeric(input$sim_ellipse_a)),
        semi_b = max(1, as.numeric(input$sim_ellipse_b)),
        semi_c = if (use_3d) max(1, as.numeric(input$sim_ellipse_c)) else 1,
        filled = identical(input$sim_ellipse_fill, "filled"),
        use_3d = use_3d
      ),
      filament = gen_filament(
        n             = n_prot,
        spacing_nm    = max(0.1, as.numeric(input$sim_fil_spacing)),
        angle_deg     = suppressWarnings(as.numeric(input$sim_fil_angle)),
        jitter_lat_nm = max(0, as.numeric(input$sim_fil_jitter_lat)),
        jitter_ax_nm  = max(0, as.numeric(input$sim_fil_jitter_ax)),
        helix_r       = if (isTRUE(input$sim_fil_helix))
                          max(0, as.numeric(input$sim_fil_helix_r))
                        else 0,
        helix_pitch   = max(0.1, as.numeric(input$sim_fil_helix_pitch)),
        use_3d        = use_3d,
        z_angle_deg   = if (use_3d)
                          suppressWarnings(as.numeric(input$sim_fil_z_angle))
                        else 0
      ),
      hollow_sphere = gen_hollow_sphere(
        n           = n_prot,
        radius      = max(1, as.numeric(input$sim_hsph_radius)),
        shell_thick = max(0, as.numeric(input$sim_hsph_thick)),
        use_3d      = use_3d
      ),
      double_ring = gen_double_ring(
        n              = n_prot,
        diameter       = max(1, as.numeric(input$sim_dring_diam)),
        ring_width     = max(0, as.numeric(input$sim_dring_width)),
        separation_nm  = max(0, as.numeric(input$sim_dring_sep)),
        z_thickness    = if (use_3d)
                           max(0, as.numeric(input$sim_dring_z_thick))
                         else 0,
        use_3d         = use_3d
      ),
      arc = gen_arc(
        n              = n_prot,
        diameter       = max(1, as.numeric(input$sim_arc_diam)),
        ring_width     = max(0, as.numeric(input$sim_arc_width)),
        arc_start_deg  = suppressWarnings(as.numeric(input$sim_arc_start)),
        arc_extent_deg = max(1, as.numeric(input$sim_arc_extent)),
        z_thickness    = if (use_3d)
                           max(0, as.numeric(input$sim_arc_z_thick))
                         else 0,
        use_3d         = use_3d
      ),
      oligoring = gen_oligoring_stack(
        n_subunits       = max(2L, as.integer(input$sim_oligoring_nsub       %||% 6L)),
        diameter         = max(1,  as.numeric(input$sim_oligoring_diam       %||% 12)),
        n_stacks         = max(1L, as.integer(input$sim_oligoring_nstack     %||% 1L)),
        separation_nm    = max(0.5, as.numeric(input$sim_oligoring_sep       %||% 7)),
        angle_offset_deg = suppressWarnings(as.numeric(input$sim_oligoring_angle_off  %||% 0)),
        base_angle_deg   = suppressWarnings(as.numeric(input$sim_oligoring_base_angle %||% 0)),
        use_3d           = use_3d
      )
    )
    # Oligoring GT count is fixed by structure; update n_prot so downstream
    # labeling code (labeling efficiency, steric check) uses the correct N.
    if (identical(input$sim_shape, "oligoring"))
      n_prot <- nrow(gt)
  }

  # Clamp GT Z to ±1000 nm regardless of shape parameters
  if (use_3d && "z" %in% names(gt))
    gt$z <- pmax(-1000, pmin(1000, gt$z))

  # ── Random XY orientation ────────────────────────────────────────────────
  # Rotates the entire GT cloud by a uniformly sampled angle.  Skipped in
  # multi-copy mode (each copy is rotated independently there) and in image
  # mode (the image defines orientation implicitly).
  if (isTRUE(input$sim_rotate) && !identical(input$sim_shape, "image") &&
      !isTRUE(input$sim_multi_copy)) {
    rot_range <- max(0, as.numeric(input$sim_rotate_range %||% 180))
    angle_rad <- runif(1, -rot_range * pi / 180, rot_range * pi / 180)
    cos_a <- cos(angle_rad); sin_a <- sin(angle_rad)
    x_new <- gt$x * cos_a - gt$y * sin_a
    y_new <- gt$x * sin_a + gt$y * cos_a
    gt$x  <- x_new
    gt$y  <- y_new
    # Rotate stored filament tangent with the same angle
    if ("ldir_type" %in% names(gt) && identical(gt$ldir_type[1L], "perp")) {
      lx_new <- gt$ldir_x * cos_a - gt$ldir_y * sin_a
      ly_new <- gt$ldir_x * sin_a + gt$ldir_y * cos_a
      gt$ldir_x <- lx_new; gt$ldir_y <- ly_new
    }
  }

  # ── Multi-structure placement ────────────────────────────────────────────
  # Place N copies of the base structure at random positions in a rectangular
  # FOV (centred at the origin), each independently rotated.  Two drivers:
  #   • Normal mode  → optional "Multi-copy tiling" (square copies FOV).
  #   • Acquisition mode → "Number of structures" placed in the X×Y FOV.
  # Disabled for image mode because the image already fills the full FOV.
  if (acq_mode && !identical(input$sim_shape, "image")) {
    do_tiling <- TRUE
    n_copies  <- max(1L, as.integer(input$sim_acq_n_structures %||% 1L))
    tile_w    <- fov_x; tile_h <- fov_y
    rot_ea    <- TRUE
  } else if (!acq_mode && isTRUE(input$sim_multi_copy) &&
             !identical(input$sim_shape, "image")) {
    do_tiling <- TRUE
    n_copies  <- max(2L, as.integer(input$sim_copies_n   %||% 5L))
    tile_w    <- tile_h <- max(1, as.numeric(input$sim_copies_fov %||% 2000))
    rot_ea    <- isTRUE(input$sim_copies_rotate)
  } else {
    do_tiling <- FALSE
  }

  if (do_tiling && n_copies > 1L) {
    copies <- lapply(seq_len(n_copies), function(i) {
      cp <- gt
      if (rot_ea) {
        a     <- runif(1, 0, 2 * pi)
        ca    <- cos(a); sa <- sin(a)
        x_rot <- cp$x * ca - cp$y * sa
        y_rot <- cp$x * sa + cp$y * ca
        cp$x  <- x_rot; cp$y <- y_rot
        # Rotate filament tangent with each per-copy rotation
        if ("ldir_type" %in% names(cp) && identical(cp$ldir_type[1L], "perp")) {
          lx_r <- cp$ldir_x * ca - cp$ldir_y * sa
          ly_r <- cp$ldir_x * sa + cp$ldir_y * ca
          cp$ldir_x <- lx_r; cp$ldir_y <- ly_r
        }
      }
      # Keep the first structure near the centre (~0-100 nm), like a real
      # MINFLUX ROI is usually centred on a feature; the rest tile randomly.
      if (i == 1L) {
        off_x <- runif(1, -50, 50)
        off_y <- runif(1, -50, 50)
      } else {
        off_x <- runif(1, -tile_w / 2, tile_w / 2)
        off_y <- runif(1, -tile_h / 2, tile_h / 2)
      }
      cp$x <- cp$x + off_x
      cp$y <- cp$y + off_y
      # Track per-copy structure centre for radial displacement
      if ("ldir_type" %in% names(cp) && identical(cp$ldir_type[1L], "radial")) {
        if ("struct_cx" %in% names(cp)) cp$struct_cx <- cp$struct_cx + off_x
        if ("struct_cy" %in% names(cp)) cp$struct_cy <- cp$struct_cy + off_y
      }
      if (use_3d && "z" %in% names(cp))
        cp$z <- pmax(-1000, pmin(1000, cp$z))
      cp
    })
    gt <- do.call(rbind, copies)
    rownames(gt) <- NULL
  }

  # ── Labeling: preset or custom ───────────────────────────────────────────
  set.seed(label_seed)
  label_key <- input$sim_label_type %||% "Custom"
  if (identical(label_key, "Custom")) {
    n_fluor       <- max(1L, as.integer(input$sim_fluor_custom     %||% 1L))
    link_nm       <- max(0,  as.numeric(input$sim_link_custom_nm   %||% 5))
    dye_spread_nm <- max(0,  as.numeric(input$sim_dye_spread_custom %||% 0))
    footprint_nm  <- max(0,  as.numeric(input$sim_footprint_custom  %||% (link_nm * 1.5)))
    blink_lambda  <- max(0,  as.numeric(input$sim_blink_mean        %||% 0.3))
    lbl_note      <- sprintf("%d fluorophore(s)/molecule, %.0f nm linker, %.0f nm dye spread",
                             n_fluor, link_nm, dye_spread_nm)
  } else {
    p             <- LABEL_PRESETS[[label_key]]
    n_fluor       <- p$n_fluor
    link_nm       <- p$link_nm
    dye_spread_nm <- p$dye_spread_nm
    footprint_nm  <- p$footprint_nm
    blink_lambda  <- if (!is.null(p$blink_lambda)) p$blink_lambda else 0.3
    lbl_note      <- p$note
  }
  eff_pct        <- max(1, min(100, as.numeric(input$sim_label_eff)))
  protein_mw_kda <- max(0, as.numeric(input$sim_protein_mw %||% 0))
  protein_r_nm   <- if (protein_mw_kda > 0) 0.66 * protein_mw_kda^(1/3) else 0
  fluoro  <- expand_to_fluorophores(gt, eff_pct, n_fluor, link_nm, dye_spread_nm,
                                     ldir_aware   = isTRUE(input$sim_ldir_aware),
                                     protein_r_nm = protein_r_nm)
  # Unique physical fluorophores (one dye = one position), before any blinking
  # or rate-based over-counting.  The confocal view renders from THIS so a
  # re-detected molecule is not rendered brighter than it physically is.
  fluoro_pool <- fluoro

  # ── Blinking / TID splitting (Normal mode only) ──────────────────────────
  # In Acquisition mode over-counting comes from the detection rate × duration
  # below, so the λ-based blink split is skipped there.
  n_fluoro_base <- nrow(fluoro)   # unique fluorophores (detectable pool)
  if (!acq_mode && isTRUE(input$sim_blink_model)) {
    bl_lambda <- max(0, as.numeric(input$sim_blink_mean %||% blink_lambda))
    fluoro    <- split_for_blinking(fluoro, bl_lambda)
  }

  # ── Localization precision: manual or CRLB-derived ──────────────────────
  # CRLB mode: σ = √[(k·σ_CRLB)² + σ₀²], σ_CRLB = L/√(8N)·g(SBR) from the
  # physical final-iteration L (smlm_sigma_loc in packages.R).
  if (identical(input$sim_prec_mode %||% "manual", "crlb")) {
    sigma_loc_xy <- max(0.1, crlb_inputs()$total)
    sigma_loc_z  <- if (use_3d) {
      z_fac <- max(0.1, as.numeric(input$sim_crlb_z_factor %||% 0.7))
      max(0.1, sigma_loc_xy * z_fac)
    } else sigma_loc_xy
  } else {
    sigma_loc_xy <- max(0.1, as.numeric(input$sim_precision_xy))
    sigma_loc_z  <- if (use_3d)
                      max(0.1, as.numeric(input$sim_precision_z %||% sigma_loc_xy))
                    else sigma_loc_xy
  }
  trace_length <- max(1, as.numeric(input$sim_trace_length))

  # ── Acquisition window + over-counting (Acquisition mode only) ────────────
  # The user sets the duration (hours + minutes) AND the structure detection
  # rate (traces/hour).  Requested structure traces = rate × duration.  Each
  # fluorophore yields a finite number of detections before PERMANENT bleaching
  # (mean = sim_acq_blinks), so the total is capped at that photon budget: a
  # short run under-samples the pool (low % resolved), a longer run resolves
  # more and over-counts, and beyond the budget the structure plateaus (dead)
  # while background keeps accumulating.  Detections are drawn WITHOUT
  # replacement from the per-fluorophore capacity, so molecules with more
  # fluorophores / more blinks are found sooner.
  # (The resolution/saturation curves are computed live in sim_acq_live();
  # here we only need the over-counted detection list for the MINFLUX view.)
  fov_area_nm2 <- fov_x * fov_y
  acq_T_sec    <- 0
  struct_rate  <- max(0, as.numeric(input$sim_struct_rate %||% 20))   # traces/hour
  n_pool       <- nrow(fluoro)        # unique fluorophores (detectable pool)
  if (acq_mode) {
    hrs       <- max(0, as.numeric(input$sim_acq_hours   %||% 0))
    mins      <- max(0, as.numeric(input$sim_acq_minutes %||% 0))
    acq_T_sec <- max(1, (hrs * 60 + mins) * 60)
    blinks    <- max(1, as.numeric(input$sim_acq_blinks %||% 5))   # detections/fluorophore
    dec_min   <- max(0.1, as.numeric(input$sim_acq_decay_min  %||% 15))
    flr_frac  <- max(0, min(1, as.numeric(input$sim_acq_floor_frac %||% 0.25)))
    det_model <- input$sim_acq_det_model %||% "decay"
    n_req     <- as.integer(round(                                # absolute-rate count
                   acq_expected_count(struct_rate / 60, acq_T_sec,  # per hour → per min
                                      dec_min, flr_frac, det_model)))
    if (n_pool > 0 && n_req > 0) {
      k_i        <- 1L + rpois(n_pool, max(0, blinks - 1))         # per-fluorophore budget
      struct_cap <- sum(k_i)
      n_s        <- min(n_req, struct_cap, 200000L)
      det_idx    <- if (struct_cap <= 5e6)
                      sample(rep(seq_len(n_pool), k_i), n_s)        # exact capacity
                    else
                      sample(seq_len(n_pool), n_s, replace = TRUE, prob = k_i)
      fluoro     <- fluoro[det_idx, , drop = FALSE]
      rownames(fluoro) <- NULL
    } else {
      fluoro     <- fluoro[0, , drop = FALSE]
    }
  }

  locs <- gen_localizations(fluoro, trace_length, sigma_loc_xy, sigma_loc_z,
                             use_geometric  = isTRUE(input$sim_tl_geometric),
                             scale_by_fluor = isTRUE(input$sim_tl_scale_fluor))

  # ── Background / noise TIDs ──────────────────────────────────────────────
  # Acquisition mode: count = background detection rate (traces/hour) × duration,
  # mirroring the structure - longer runs accumulate more non-specific TIDs.
  # Normal mode: a plain fixed count.  Background is spread across the FOV when
  # set, else the multi-copy FOV, else the GT bounding box + 10 % padding.
  # Each background detection is a fresh random position. TID numbering
  # continues from the last real TID.
  n_bg   <- 0L
  bg_pos <- NULL
  bg_rate <- max(0, as.numeric(input$sim_bg_rate %||% 5))   # traces/hour
  if (isTRUE(input$sim_bg_noise)) {
    copies_mode <- !acq_mode && isTRUE(input$sim_multi_copy) &&
                   !identical(input$sim_shape, "image")
    if (use_fov) {
      bg_xc <- 0; bg_yc <- 0
      bg_hx <- fov_x / 2; bg_hy <- fov_y / 2
    } else if (copies_mode) {
      fov_half <- max(1, as.numeric(input$sim_copies_fov %||% 2000)) / 2
      bg_xc <- 0; bg_yc <- 0
      bg_hx <- fov_half; bg_hy <- fov_half
    } else {
      pad   <- 0.1
      x_rng <- range(gt$x, na.rm = TRUE); y_rng <- range(gt$y, na.rm = TRUE)
      x_pad <- diff(x_rng) * pad; y_pad <- diff(y_rng) * pad
      bg_xc <- mean(x_rng); bg_yc <- mean(y_rng)
      bg_hx <- diff(x_rng) / 2 + x_pad; bg_hy <- diff(y_rng) / 2 + y_pad
    }

    if (acq_mode) {
      n_bg <- as.integer(round(bg_rate * acq_T_sec / 3600))   # per hour × hours
    } else {
      n_bg <- as.integer(input$sim_bg_n %||% 20L)
    }
    n_bg <- max(0L, min(200000L, n_bg))

    if (n_bg > 0L) {
      bg_df <- data.frame(x = runif(n_bg, bg_xc - bg_hx, bg_xc + bg_hx),
                          y = runif(n_bg, bg_yc - bg_hy, bg_yc + bg_hy))
      if (use_3d && "z" %in% names(gt)) {
        z_rng  <- range(gt$z, na.rm = TRUE)
        # FOV mode → fill the FOV Z slab; otherwise track the GT extent, with
        # a precision-based fallback for degenerate (planar) structures.
        z_half <- if (use_fov)
                    fov_z / 2
                  else if (diff(z_rng) < 1) max(50, sigma_loc_z * 10)
                  else diff(z_rng) / 2 * 1.1
        z_mid  <- if (use_fov) 0 else mean(z_rng)
        bg_df$z <- runif(n_bg, z_mid - z_half, z_mid + z_half)
      }
      bg_df$mol_id <- NA_integer_
      # Save XY(Z) positions so the confocal view can render bg spots
      bg_pos <- bg_df[, intersect(c("x","y","z"), names(bg_df)), drop = FALSE]
      # Each background TID gets its own locs via gen_localizations
      bg_locs <- gen_localizations(bg_df, trace_length, sigma_loc_xy,
                                    sigma_loc_z,
                                    use_geometric  = isTRUE(input$sim_tl_geometric),
                                    scale_by_fluor = FALSE)
      # Offset TID numbers so they don't clash with real TIDs
      bg_locs$tid <- as.character(as.integer(bg_locs$tid) + nrow(fluoro))
      locs <- rbind(locs, bg_locs)
    }
  }

  # ── RMSE: per-TID displacement error (GT position vs trace median) ──────────
  # Only computed for real TIDs (1 … nrow(fluoro)); background TIDs excluded.
  # Pre-split locs by TID for O(1) lookup per trace - essential when multi-
  # copy produces thousands of TIDs.  Capped at 1 000 sampled TIDs so the
  # computation stays fast even with very large multi-copy scenes.
  rmse_vals <- NULL
  if ("mol_id" %in% names(fluoro) && nrow(gt) > 0 && nrow(locs) > 0) {
    n_real    <- nrow(fluoro)
    tid_samp  <- if (n_real > 1000L) sample(n_real, 1000L) else seq_len(n_real)
    real_rows <- which(suppressWarnings(as.integer(locs$tid)) <= n_real &
                       !is.na(suppressWarnings(as.integer(locs$tid))))
    locs_by_tid <- split(real_rows, locs$tid[real_rows])
    errs <- vapply(tid_samp, function(k) {
      rows <- locs_by_tid[[as.character(k)]]
      if (is.null(rows) || length(rows) == 0) return(NA_real_)
      lk <- locs[rows, , drop = FALSE]
      gi <- fluoro$mol_id[k]
      if (is.na(gi) || gi < 1L || gi > nrow(gt)) return(NA_real_)
      dx <- median(lk$x) - gt$x[gi]
      dy <- median(lk$y) - gt$y[gi]
      dz <- if (use_3d && "z" %in% names(lk) && "z" %in% names(gt))
              median(lk$z) - gt$z[gi]
            else 0
      sqrt(dx^2 + dy^2 + dz^2)
    }, numeric(1L))
    rmse_vals <- errs[!is.na(errs)]
  }

  # ── Steric exclusion check ───────────────────────────────────────────────
  # Compute mean nearest-neighbour distance among GT molecule positions and
  # compare it to the label complex footprint.  A footprint larger than the
  # protein spacing means antibodies physically cannot label every molecule.
  steric_warn <- NULL
  gt_nnd      <- NA_real_
  if (nrow(gt) >= 2L) {
    gt_pts <- as.matrix(gt[, intersect(c("x","y","z"), names(gt)), drop = FALSE])
    samp_i <- if (nrow(gt_pts) > 300L) sample(nrow(gt_pts), 300L) else seq_len(nrow(gt_pts))
    d_mat  <- as.matrix(dist(gt_pts[samp_i, , drop = FALSE]))
    diag(d_mat) <- Inf
    gt_nnd <- mean(apply(d_mat, 1, min), na.rm = TRUE)
  }
  # Effective exclusion zone = label footprint + protein diameter.
  # Proteins physically can't overlap, so even before the antibody is considered
  # the centre-to-centre spacing must be ≥ 2 × protein radius.
  eff_exclusion_nm <- max(footprint_nm, 2 * protein_r_nm, na.rm = TRUE)
  if (!is.na(eff_exclusion_nm) && !is.na(gt_nnd)) {
    if (gt_nnd < eff_exclusion_nm * 0.5) {
      steric_warn <- list(
        level = "danger",
        msg   = sprintf(
          paste0("Severe steric exclusion: combined exclusion zone (~%.0f nm: ",
                 "%.0f nm label footprint + %.0f nm protein diameter) exceeds ",
                 "mean protein spacing (~%.0f nm). Full stoichiometric labeling is ",
                 "physically impossible at this density."),
          eff_exclusion_nm, footprint_nm, 2 * protein_r_nm, gt_nnd)
      )
    } else if (gt_nnd < eff_exclusion_nm) {
      steric_warn <- list(
        level = "warning",
        msg   = sprintf(
          paste0("Steric exclusion likely: combined exclusion zone (~%.0f nm) ",
                 "approaches mean protein spacing (~%.0f nm). Some molecules may be ",
                 "inaccessible; effective labeling efficiency could be lower than expected."),
          eff_exclusion_nm, gt_nnd)
      )
    }
  }

  # ── Snapshot detection times (Acquisition mode) ──────────────────────────
  # Assign front-loaded times to THIS realization's TIDs so the MINFLUX Time
  # colour mode and CSV export have realistic timestamps.  The Acquisition-tab
  # saturation/rate curves are computed live in sim_acq_live() (decoupled so
  # they react to rate/duration/background without re-Generating).
  if (acq_mode && nrow(locs) > 0 && acq_T_sec > 0) {
    all_tids   <- unique(locs$tid)
    n_real     <- nrow(fluoro)
    decay_min  <- max(0.1, as.numeric(input$sim_acq_decay_min  %||% 15))
    floor_frac <- max(0, min(1, as.numeric(input$sim_acq_floor_frac %||% 0.25)))
    det_model  <- input$sim_acq_det_model %||% "decay"
    tid_times  <- setNames(rep(NA_real_, length(all_tids)), all_tids)
    real_ids   <- intersect(as.character(seq_len(n_real)), all_tids)
    if (length(real_ids) > 0)
      tid_times[real_ids] <- sim_sample_detection_times(
        length(real_ids), acq_T_sec, decay_min * 60, floor_frac, det_model)
    bg_ids <- setdiff(all_tids, real_ids)
    if (length(bg_ids) > 0)
      tid_times[bg_ids] <- runif(length(bg_ids), 0, acq_T_sec)
    within   <- ave(seq_len(nrow(locs)), locs$tid,
                    FUN = function(i) (seq_along(i) - 1L) * 0.01)
    locs$tim <- unname(tid_times[locs$tid]) + within
  }

  list(gt = gt, fluoro = fluoro, fluoro_pool = fluoro_pool,
       locs = locs,
       link_nm = link_nm, n_fluor = n_fluor, lbl_note = lbl_note,
       dye_spread_nm = dye_spread_nm, blink_lambda = blink_lambda,
       eff_pct = eff_pct,
       protein_r_nm = protein_r_nm,
       sigma_loc_xy = sigma_loc_xy, sigma_loc_z = sigma_loc_z,
       n_prot = n_prot, n_fluoro = nrow(fluoro),
       n_fluoro_base = n_fluoro_base,
       n_locs = nrow(locs), n_bg = n_bg, use_3d = use_3d,
       img_b64 = img_b64,
       img_w_nm = img_w_nm, img_h_nm = img_h_nm,
       rmse_vals = rmse_vals, bg_pos = bg_pos,
       steric_warn = steric_warn,
       footprint_nm = footprint_nm,
       gt_nnd = gt_nnd, acq_mode = acq_mode,
       blink_active = (!acq_mode) && isTRUE(input$sim_blink_model),
       blink_lambda_used = if ((!acq_mode) && isTRUE(input$sim_blink_model))
                             max(0, as.numeric(input$sim_blink_mean %||% blink_lambda))
                           else 0)
}, ignoreNULL = TRUE)

# ── Live acquisition timeline ────────────────────────────────────────────────
# Recomputes the saturation / detection-rate curves from the last-generated
# fluorophore pool whenever the acquisition parameters change - so the
# Acquisition tab reacts to rate / duration / background / bleaching live, with
# no need to re-Generate the structure.  A fixed seed keeps the curve stable
# for a given parameter set (it only moves when an input actually changes).
sim_acq_live <- reactive({
  if (!identical(input$sim_mode %||% "normal", "acquisition")) return(NULL)
  res <- sim_result()
  if (is.null(res) || is.null(res$fluoro_pool)) return(NULL)
  pool_mol <- res$fluoro_pool$mol_id
  if (length(pool_mol) == 0) return(NULL)
  hrs   <- max(0, as.numeric(input$sim_acq_hours   %||% 0))
  mins  <- max(0, as.numeric(input$sim_acq_minutes %||% 0))

  # Emitter-crowding weight per molecule: a molecule with many other labeled
  # molecules inside the donut radius is harder to single out, so its effective
  # discovery rate is scaled by 1/(1 + neighbours).  Aligned to the sorted
  # unique mol_ids (= the order compute_acq_timeline builds f_m in).
  crowd_m <- NULL
  if (isTRUE(input$sim_crowding)) {
    u  <- sort(unique(pool_mol))
    gp <- res$gt[u, , drop = FALSE]
    cols <- intersect(c("x", "y", "z"), names(gp))
    if (length(u) >= 2 && length(cols) >= 2) {
      R   <- max(1, as.numeric(input$sim_crowd_radius %||% 50))
      pts <- as.matrix(gp[, cols, drop = FALSE])
      nn  <- RANN::nn2(pts, pts, k = min(nrow(pts), 64L),
                       searchtype = "radius", radius = R)
      neigh   <- rowSums(nn$nn.idx > 0) - 1L
      crowd_m <- 1 / (1 + neigh)
    }
  }

  set.seed(1L)                          # stable curve per parameter set
  compute_acq_timeline(
    pool_mol    = pool_mol,
    n_mol_total = nrow(res$gt),
    struct_rate = max(0, as.numeric(input$sim_struct_rate %||% 20)),
    T_sec       = max(1, (hrs * 60 + mins) * 60),
    blinks      = max(1, as.numeric(input$sim_acq_blinks %||% 5)),
    decay_min   = max(0.1, as.numeric(input$sim_acq_decay_min %||% 15)),
    floor_frac  = max(0, min(1, as.numeric(input$sim_acq_floor_frac %||% 0.25))),
    bg_rate     = max(0, as.numeric(input$sim_bg_rate %||% 5)),
    resolve_q   = as.numeric(input$sim_acq_resolve_q %||% 0.9),
    crowd_m     = crowd_m,
    model       = input$sim_acq_det_model %||% "decay")
})
