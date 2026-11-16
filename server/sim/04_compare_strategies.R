## ─────────────────────────────────────────────────────────────────────────────
## server/sim/04_compare_strategies.R  -  compare labelling strategies
## ─────────────────────────────────────────────────────────────────────────────
## Part of sim_server(). Not a standalone script: its statements are evaluated in
## the server function's environment, in file-name order (server/parts.R).
##
## compare_result: the same ground truth labelled two ways.
##
## INTERFACE (generated: Rscript tools/server_map.R --headers; do not edit by hand)
##   DEFINES for other parts:
##     compare_result
##   READS from other parts:
##     02_smlm_simulation: sim_result
## END INTERFACE
## ─────────────────────────────────────────────────────────────────────────────

# ── Compare-strategies reactive ─────────────────────────────────────────────
# Applies two labeling strategies to the same GT from sim_result().
# Geometry is NOT re-run; seeds come from current label_seed + offsets.
compare_result <- eventReactive(input$sim_compare_run, {
  res <- sim_result()
  shiny::validate(shiny::need(
    !is.null(res),
    "Run the main simulation first to generate a ground truth."))

  gt   <- res$gt
  use_3d <- res$use_3d
  sigma_loc_xy <- res$sigma_loc_xy
  sigma_loc_z  <- res$sigma_loc_z
  trace_length <- max(1, as.numeric(input$sim_trace_length))

  # Base seed: use the fixed label seed if set, else a fresh one per compare run
  base_seed <- if (isTRUE(input$sim_fix_label_seed)) {
    suppressWarnings(as.integer(input$sim_label_seed))
  } else {
    sample.int(.Machine$integer.max, 1L)
  }
  if (!is.finite(base_seed)) base_seed <- 1L

  cmp_protein_mw_kda <- max(0, as.numeric(input$sim_protein_mw %||% 0))
  cmp_protein_r_nm   <- if (cmp_protein_mw_kda > 0) 0.66 * cmp_protein_mw_kda^(1/3) else 0

  run_strategy <- function(label_key, eff_pct, seed_val, suffix) {
    set.seed(seed_val)
    if (identical(label_key, "Custom")) {
      n_fluor       <- max(1L, as.integer(input[[paste0("sim_cmp_", suffix, "_fluor")]] %||% 2L))
      link_nm       <- max(0,  as.numeric(input[[paste0("sim_cmp_", suffix, "_link")]]  %||% 10))
      dye_spread_nm <- max(0,  as.numeric(input[[paste0("sim_cmp_", suffix, "_spread")]] %||% 0))
      footprint_nm  <- max(0,  as.numeric(input[[paste0("sim_cmp_", suffix, "_foot")]]  %||% 8))
      blink_lambda  <- max(0,  as.numeric(input$sim_blink_mean %||% 0.3))
      lbl_note      <- sprintf("Custom: %d fluor / %.0f nm linker",
                               n_fluor, link_nm)
    } else {
      p             <- LABEL_PRESETS[[label_key]]
      n_fluor       <- p$n_fluor
      link_nm       <- p$link_nm
      dye_spread_nm <- p$dye_spread_nm
      footprint_nm  <- p$footprint_nm
      blink_lambda  <- p$blink_lambda %||% 0.3
      lbl_note      <- p$note
    }

    fluoro <- expand_to_fluorophores(
      gt, eff_pct, n_fluor, link_nm, dye_spread_nm,
      ldir_aware   = isTRUE(input$sim_ldir_aware),
      protein_r_nm = cmp_protein_r_nm)

    if (isTRUE(input$sim_blink_model))
      fluoro <- split_for_blinking(
        fluoro, max(0, as.numeric(input$sim_blink_mean %||% blink_lambda)))

    locs <- gen_localizations(fluoro, trace_length, sigma_loc_xy, sigma_loc_z,
                               use_geometric  = isTRUE(input$sim_tl_geometric),
                               scale_by_fluor = isTRUE(input$sim_tl_scale_fluor))

    # RMSE (capped at 500 TIDs for speed)
    rmse_vals <- NULL
    if ("mol_id" %in% names(fluoro) && nrow(gt) > 0 && nrow(locs) > 0) {
      n_samp <- min(nrow(fluoro), 500L)
      real_rows <- which(suppressWarnings(as.integer(locs$tid)) <= nrow(fluoro) &
                         !is.na(suppressWarnings(as.integer(locs$tid))))
      locs_by_tid <- split(real_rows, locs$tid[real_rows])
      errs <- vapply(seq_len(n_samp), function(k) {
        rows <- locs_by_tid[[as.character(k)]]
        if (is.null(rows) || length(rows) == 0L) return(NA_real_)
        lk <- locs[rows, , drop = FALSE]
        gi <- fluoro$mol_id[k]
        if (is.na(gi) || gi < 1L || gi > nrow(gt)) return(NA_real_)
        dx <- median(lk$x) - gt$x[gi]
        dy <- median(lk$y) - gt$y[gi]
        sqrt(dx^2 + dy^2)
      }, numeric(1L))
      rmse_vals <- errs[!is.na(errs)]
    }

    list(fluoro = fluoro, locs = locs, link_nm = link_nm,
         n_fluor = n_fluor, n_fluoro = nrow(fluoro),
         lbl_note = lbl_note, footprint_nm = footprint_nm,
         protein_r_nm = cmp_protein_r_nm,
         rmse_vals = rmse_vals)
  }

  eff_a <- max(1, min(100, as.numeric(input$sim_cmp_eff_a %||% 80)))
  eff_b <- max(1, min(100, as.numeric(input$sim_cmp_eff_b %||% 80)))
  a <- run_strategy(input$sim_cmp_label_a %||%
                      "Direct tag (SNAP / Halo / Click, DOL=1)",
                    eff_a, base_seed, "a")
  b <- run_strategy(input$sim_cmp_label_b %||%
                      "Primary monoclonal IgG + secondary IgG",
                    eff_b, base_seed + 1L, "b")

  list(gt = gt, use_3d = use_3d, a = a, b = b,
       sigma_loc_xy = sigma_loc_xy,
       lbl_a = input$sim_cmp_label_a %||%
                 "Direct tag (SNAP / Halo / Click, DOL=1)",
       lbl_b = input$sim_cmp_label_b %||%
                 "Primary monoclonal IgG + secondary IgG")
}, ignoreNULL = TRUE)
