output$download_two_color_settings <- downloadHandler(
  filename = function() {
    base <- if (!is.null(input$active_file) && nzchar(input$active_file)) tools::file_path_sans_ext(input$active_file) else "smlm"
    paste0(base, "_channel_settings_", format(Sys.time(), "%Y%m%d-%H%M%S"), ".json")
  },
  content = function(file) {
    s <- list(
      # Master toggle - as with cluster2, without it the restored 2-color
      # parameters had no effect because 2-color stayed switched off.
      enable_two_color = isTRUE(input$enable_two_color),
      color_by_fluo = isTRUE(input$color_by_fluo),
      use_pooled_dcr = isTRUE(input$use_pooled_dcr),
      dcr_split_manual = isTRUE(input$dcr_split_manual),
      dcr_group1 = unname(input$dcr_group1),
      dcr_group2 = unname(input$dcr_group2),
      color_group1 = input$color_group1,
      color_group2 = input$color_group2,

      
      # Clustering parameters are NOT written here. One pipeline means one set
      # of de-blink / Median Clustering controls, and those are saved by the
      # main settings file. Writing a second copy here is what let the two
      # drift apart in the first place.
      show_tc1_alpha = isTRUE(input$show_tc1_alpha),
      tc1_alpha = input$tc1_alpha,
      tc1_alpha_opacity = input$tc1_alpha_opacity,
      tc1_alpha_min_pts = input$tc1_alpha_min_pts
      # Median Cluster wrap / isosurface settings are saved with Median
      # Clustering, not here.
    )
    
    jsonlite::write_json(s, path = file, auto_unbox = TRUE, pretty = TRUE)
  }
)
output$download_filtered_csv <- downloadHandler(
  filename = function() {
    paste0(
      tools::file_path_sans_ext(input$active_file %||% "smlm"),
      "_filtered_",
      format(Sys.time(), "%Y%m%d_%H%M%S"),
      ".csv"
    )
  },
  
  content = function(file) {
    df <- plottable_data()

    shiny::validate(
      shiny::need(!is.null(df) && nrow(df) > 0, "No data to export.")
    )

    # Core MINFLUX columns - always included when present
    core_cols <- c("x", "y", "z", "tid", "efo", "cfr", "eco", "fbg", "dcr")
    # Optional columns included only when they exist and carry real data
    # dcr_pooled is exported so the pooled value SURVIVES a round-trip: it is
    # derived from the per-iteration dcr1, which a CSV does not carry, so
    # without this it could never be recomputed after re-import.
    optional_cols <- c("fluo", "tim", "itr", "dwell", "dcr_pooled",
                       "pcid", "scid", "n_in_cluster2",
                       ".dataset_tag", "analysis_group",
                       "mplex_channel", "mplex_label", "mplex_color")

    # Multiplex: record per row the TOTAL alignment shift already in x/y/z
    # (shift loaded with the file + nudges since), so a re-import can show
    # it and the file documents how its channels were aligned.
    act_mp <- isolate(uploaded_df_rv())
    mp_sp  <- attr(act_mp, "mplex_spec")
    if ("mplex_channel" %in% names(df) && !is.null(mp_sp)) {
      now <- .mplex_store_off(attr(act_mp, "mplex_store_name"))
      tot <- t(vapply(mp_sp, function(c) {
        pr <- if (length(c$prior) == 3L) unname(c$prior) else c(0, 0, 0)
        pr + unname(.mplex_off_of(now, .mplex_key(c)))
      }, numeric(3)))
      m <- match(df$mplex_channel, vapply(mp_sp, function(c) as.numeric(c$channel), numeric(1)))
      df$mplex_dx_nm <- tot[m, 1]; df$mplex_dy_nm <- tot[m, 2]; df$mplex_dz_nm <- tot[m, 3]
      optional_cols <- c(optional_cols, "mplex_dx_nm", "mplex_dy_nm", "mplex_dz_nm")
    }

    keep <- c(
      intersect(core_cols, names(df)),
      intersect(optional_cols, names(df))
    )
    # Drop optional columns that are entirely NA (e.g. clustering not run)
    keep <- keep[sapply(keep, function(col) !all(is.na(df[[col]])))]

    df <- df[, keep, drop = FALSE]

    # Friendly export names
    names(df)[names(df) == ".dataset_tag"]  <- "dataset_tag"
    names(df)[names(df) == "n_in_cluster2"] <- "n_in_scid"

    # Axial provenance: z below is written at this factor. Re-importing divides
    # it out before applying the then-current calibration, so a round-trip is
    # exact instead of relying on the 0.7 assumption in
    # detect_csv_z_calibration(). Written even for 2D data - it costs one
    # column and removes the guess for every future reader.
    # The COMMITTED factor - that is what the exported z is actually at. Using
    # the pending input value would mislabel the file whenever the control has
    # been typed into but not applied.
    df$z_scale_applied <- z_scale_applied()

    write.csv(df, file, row.names = FALSE)
  }
)
observeEvent(input$download_filtered_csv, {
  showNotification("Exported filtered dataset.", type = "message", duration = 3)
})
observeEvent(input$load_settings, {
  req(input$load_settings$datapath)
  raw <- jsonlite::fromJSON(input$load_settings$datapath)
  # A multiplex settings file has its own loader; applying it here would read
  # its fields as filter values.
  if (identical(raw$settings_type, MPLEX_SETTINGS_TYPE)) {
    showNotification("This is a multiplex settings file. Load it under Multiplexed MINFLUX.",
                     type = "warning", duration = 8)
    return()
  }
  
  # Settings that change results: from the file's smlm_settings block, or from
  # the older layout for files saved before it existed. Some are restored below
  # by the filter / crop / QC code from their older places (it clamps ranges to
  # the loaded data), so the block's values are copied there first; every other
  # setting is applied by restore_settings_plan(). Settings a file does not hold
  # are left as they are.
  sf  <- smlm_settings_from_file(raw)
  raw <- smlm_settings_to_legacy(raw, sf$settings,
                                 intersect(sf$present, .SETTINGS_RESTORED_WITH_FILTERS))

  base_fil <- filters_defaults_from_data(uploaded_df())
  fil <- modifyList(base_fil, raw, keep.null = TRUE)
  ov  <- modifyList(overlay_defaults(), raw$overlay %||% list(), keep.null = TRUE)
  hm  <- modifyList(heatmap_defaults(),  raw$heatmap %||% list(), keep.null = TRUE)
  
  # 1) Push everything EXCEPT tid_length_range first
  fil_no_tid <- fil
  fil_no_tid$tid_length_range <- NULL
  
  session$onFlushed(function() {
    apply_filters(fil_no_tid)
    # Restore view_SD from the saved file (apply_filters no longer touches it).
    if (!is.null(fil$view_SD)) updateCheckboxInput(session, "view_SD", value = isTRUE(fil$view_SD))
    # Restore the min-localizations-for-spread setting.
    if (!is.null(fil$sd_min_locs)) {
      v <- suppressWarnings(as.integer(fil$sd_min_locs))
      if (length(v) == 1 && is.finite(v)) updateNumericInput(session, "sd_min_locs", value = v)
    }

    # Restore trimming BEFORE the ranges below: it changes which localizations
    # exist, so the nm bounds have to land on the same point set they were set
    # on. Files written before this setting existed carry no field and were
    # saved with no trimming.
    updateCheckboxInput(session, "trim_outliers", value = isTRUE(fil$trim_outliers))
    if (!is.null(fil$trim_k)) {
      v <- suppressWarnings(as.numeric(fil$trim_k))
      if (length(v) == 1 && is.finite(v)) updateNumericInput(session, "trim_k", value = v)
    }

    # ── Restore spatial crop / trimmer (changes the point set) ─────────────
    cr <- raw$crop %||% list()
    # 2D crop box: a server-side reactiveVal, so no DOM round-trip needed.
    cx <- suppressWarnings(as.numeric(unlist(cr$crop2d_x)))
    cy <- suppressWarnings(as.numeric(unlist(cr$crop2d_y)))
    if (length(cx) == 2 && length(cy) == 2 && all(is.finite(c(cx, cy))))
      smlm_2d_crop_rv(list(x = cx, y = cy)) else smlm_2d_crop_rv(NULL)
    # 3D box/sphere crop: static inputs restore directly.
    updateCheckboxInput(session, "enable_trimmer", value = isTRUE(cr$enable_trimmer))
    if (!is.null(cr$trim_mode)) updateRadioButtons(session, "trim_mode", selected = cr$trim_mode)
    for (nm in c("trim_cx", "trim_cy", "trim_cz", "trim_r")) {
      v <- suppressWarnings(as.numeric(cr[[nm]]))
      if (length(v) == 1 && is.finite(v)) updateNumericInput(session, nm, value = v)
    }
    # Re-activate the committed crop region so the reloaded file is cropped as saved.
    ta <- cr$trim_applied
    if (!is.null(ta) && !is.null(ta$mode)) {
      if (identical(ta$mode, "sphere")) {
        cc <- suppressWarnings(as.numeric(unlist(ta$c))); rr <- suppressWarnings(as.numeric(ta$r))
        trim_applied_rv(if (length(cc) == 3 && all(is.finite(cc)) && is.finite(rr))
          list(mode = "sphere", c = cc, r = rr) else NULL)
      } else {
        num2 <- function(v) { v <- suppressWarnings(as.numeric(unlist(v))); if (length(v) == 2 && all(is.finite(v))) v else NULL }
        trim_applied_rv(list(mode = "box", x = num2(ta$x), y = num2(ta$y), z = num2(ta$z)))
      }
    } else trim_applied_rv(NULL)
    # Box ranges live in a renderUI (trim_box_ui); set them after it rebuilds.
    if (isTRUE(cr$enable_trimmer) && identical(cr$trim_mode %||% "box", "box")) {
      bx <- suppressWarnings(as.numeric(unlist(cr$trim_x_range)))
      by <- suppressWarnings(as.numeric(unlist(cr$trim_y_range)))
      bz <- suppressWarnings(as.numeric(unlist(cr$trim_z_range)))
      later::later(function() {
        if (length(bx) == 2 && all(is.finite(bx))) shinyWidgets::updateNumericRangeInput(session, "trim_x_range", value = bx)
        if (length(by) == 2 && all(is.finite(by))) shinyWidgets::updateNumericRangeInput(session, "trim_y_range", value = by)
        if (length(bz) == 2 && all(is.finite(bz))) shinyWidgets::updateNumericRangeInput(session, "trim_z_range", value = bz)
        plot_trigger(Sys.time())
      }, delay = 0.4)
    }

    # ── Restore fiducial-scan exclusion + its thresholds ───────────────────
    fd <- raw$fiducial %||% list()
    updateCheckboxInput(session, "fid_exclude",   value = isTRUE(fd$fid_exclude))
    updateCheckboxInput(session, "fid_all_files", value = isTRUE(fd$fid_all_files))
    for (nm in c("fid_efo_z", "fid_dur_x", "fid_min_locs", "fid_len_z", "fid_recur_r", "fid_recur_min")) {
      v <- suppressWarnings(as.numeric(fd[[nm]]))
      if (length(v) == 1 && is.finite(v)) updateNumericInput(session, nm, value = v)
    }

    # ── Restore analysis options that change results ───────────────────────
    an <- raw$analysis %||% list()
    for (nm in c("use_pooled_dcr", "spatial_edge_correct", "spatial_use_support",
                 "dens_norm_enable", "dens_thin_enable"))
      if (!is.null(an[[nm]])) {
        # A loaded column is not the user switching it (see dcr_track).
        if (nm == "use_pooled_dcr") dcr_track$pooled <- isTRUE(an[[nm]])
        updateCheckboxInput(session, nm, value = isTRUE(an[[nm]]))
      }
    # Which points draw the observation window changes every area on the tab.
    if (!is.null(an$dens_support_from))
      updateRadioButtons(session, "dens_support_from", selected = an$dens_support_from)

    # ── Restore Local density & hotspots ──────────────────────────────
    hs <- raw$hotspot %||% list()
    for (nm in c("hot_k", "hot_enh_thresh", "hot_quadrat", "hot_minpts",
                 "hot_target_frac", "hot_nsim", "hot_box_top")) {
      v <- suppressWarnings(as.numeric(hs[[nm]]))
      if (length(v) == 1 && is.finite(v)) updateNumericInput(session, nm, value = v)
    }
    for (nm in c("hot_unit", "hot_palette", "hot_metric"))
      if (!is.null(hs[[nm]])) updateSelectInput(session, nm, selected = hs[[nm]])
    if (!is.null(hs$hot_curve_mode))
      updateRadioButtons(session, "hot_curve_mode", selected = hs$hot_curve_mode)
    for (nm in c("hot_use_3d", "hot_show_foci"))
      if (!is.null(hs[[nm]])) updateCheckboxInput(session, nm, value = isTRUE(hs[[nm]]))

    # Update sd_filt so the actual spread filter reflects the loaded settings.
    # NA in the upper bound means the filter was at full range when saved
    # (Inf serialises to null/NA in JSON); restore it as Inf.
    load_sd_rng <- function(v) {
      if (is.null(v)) return(NULL)
      r <- as.numeric(unlist(v))
      if (length(r) != 2) return(NULL)
      if (is.na(r[2]) || is.infinite(r[2])) r[2] <- Inf
      r
    }
    loaded_x <- load_sd_rng(fil$sd_x_range)
    loaded_y <- load_sd_rng(fil$sd_y_range)
    loaded_z <- load_sd_rng(fil$sd_z_range)
    if (!is.null(loaded_x)) sd_filt$x <- loaded_x
    if (!is.null(loaded_y)) sd_filt$y <- loaded_y
    if (!is.null(loaded_z)) sd_filt$z <- loaded_z

    # Set tid_length_range after the filter UI is in the DOM; coerce to numeric
    # in case jsonlite parsed the JSON array as a list.
    if (!is.null(fil$tid_length_range) && length(fil$tid_length_range) == 2) {
      df <- isolate(metadata_filtered_df())
      tl <- table(df$tid)
      if (length(tl)) {
        min_len <- min(tl); max_len <- max(tl)
        v <- as.numeric(unlist(fil$tid_length_range))
        v[1] <- max(min_len, min(v[1], max_len))
        v[2] <- max(min_len, min(v[2], max_len))
        shinyWidgets::updateNumericRangeInput(session, "tid_length_range", value = v)
      }
    }

    # ── Restore every other setting that changes results ───────────────────
    # Channels, de-blinking, median clustering, channel scope, density unit.
    # Before this, a reloaded settings file restored the filters but clustered
    # with whatever the session had.
    restore_settings_plan(smlm_settings_restore_plan(
      sf$settings, setdiff(sf$present, .SETTINGS_RESTORED_WITH_FILTERS)))

    # Axial calibration rescales every z distance; restore and commit it, as the
    # Median Cluster settings file does.
    zs <- suppressWarnings(as.numeric(raw$z_scale))
    if (length(zs) == 1L && is.finite(zs) && zs > 0) {
      updateNumericInput(session, "z_scale", value = zs)
      commit_z_scale(zs, reason = "Axial calibration from settings")
    }
    # Channel colours (display only, but part of the saved two-colour block).
    tc <- raw$two_color %||% list()
    for (nm in c("color_group1", "color_group2"))
      if (is.character(tc[[nm]]) && length(tc[[nm]]) == 1L)
        colourpicker::updateColourInput(session, nm, value = tc[[nm]])

    apply_overlay(ov)
    apply_heatmap(hm)
    # Re-run the plotted/exported pipeline so a restored 2D crop or trimmer
    # takes effect (plottable_data is gated on plot_trigger).
    plot_trigger(Sys.time())
  }, once = TRUE)
})

# Settings the loader above restores through the filter / crop / QC code, from
# their places in the older file layout. Every other setting in
# SMLM_SETTINGS_SPEC is restored by restore_settings_plan().
.SETTINGS_RESTORED_WITH_FILTERS <- c(
  "efo_range", "cfr_range", "dcr_range", "eco_range", "tim_range", "dwell_range",
  "fluo_range", "tid_length_range", "sd_min_locs", "trim_outliers", "trim_k",
  "view_SD", "target_sd_x", "target_sd_y", "fid_exclude", "use_pooled_dcr",
  "dens_norm_enable", "dens_thin_enable")

# Apply a smlm_settings_restore_plan(): one widget update per setting, in order.
restore_settings_plan <- function(plan) {
  vals <- attr(plan, "values")
  for (i in seq_len(nrow(plan))) {
    id <- plan$id[i]; v <- vals[[i]]
    switch(plan$widget[i],
      checkbox = updateCheckboxInput(session, id, value = isTRUE(v)),
      numeric  = if (length(v) == 1L && !is.na(v)) updateNumericInput(session, id, value = v),
      range    = if (length(v) == 2L) {
        # A loaded DCR range is not the user typing one (see dcr_track).
        if (id %in% c("dcr_group1", "dcr_group2")) dcr_note_written(v)
        shinyWidgets::updateNumericRangeInput(session, id, value = v)
      },
      radio    = updateRadioButtons(session, id, selected = v),
      select   = updateSelectInput(session, id, selected = v))
  }
  invisible(plan)
}
cluster2_file_reset_key <- reactiveVal(0)

output$cluster2_settings_file_ui <- renderUI({
  cluster2_file_reset_key()  
  fileInput("cluster2_settings_file", "Load Median Cluster settings (.json)", accept = c(".json"))
})

observeEvent(input$cluster2_settings_file, {
  req(input$cluster2_settings_file$datapath)
  
  s <- jsonlite::read_json(
    input$cluster2_settings_file$datapath,
    simplifyVector = TRUE
  )
  
  apply_cluster2_settings(session, s)
  
  showNotification(
    paste0(
      "Median Cluster settings loaded: ",
      input$cluster2_settings_file$name
    ),
    type = "message",
    duration = 3
  )
})
two_color_file_reset_key <- reactiveVal(0)

output$two_color_settings_file_ui <- renderUI({
  two_color_file_reset_key()
  fileInput("two_color_settings_file", "Load channel settings (.json)", accept = c(".json"))
})
observeEvent(input$two_color_settings_file, {
  req(input$two_color_settings_file$datapath)
  s <- jsonlite::read_json(input$two_color_settings_file$datapath, simplifyVector = TRUE)
  apply_two_color_settings(session, s)
  
  showNotification(
    paste0("2-color settings loaded: ", input$two_color_settings_file$name),
    type = "message", duration = 3
  )
})
observeEvent(input$minpts2 < 4L, {
  mp <- as.integer(input$minpts2 %||% 4L)
  if (!is.na(mp) && mp < 4L) {
    updateNumericInput(session, "minpts2", value = 4L)
    showNotification("Median Cluster minPts raised to 4 to keep geometry stable (convex hull / fits need ≥4 points).",
                     type = "warning", duration = 4)
  }
})

observeEvent(input$reset_filters, {
  req(input$file)
  was_view_SD <- isTRUE(isolate(input$view_SD))
  vals <- filters_defaults_from_data(uploaded_df())
  # Reset SD filter synchronously - Inf upper bound means "pass all".
  sd_filt$x <- c(0, Inf); sd_filt$y <- c(0, Inf); sd_filt$z <- c(0, Inf)
  # Trimming removes localizations, so "reset filters" has to put them back -
  # otherwise the spread the user returns to is not the raw spread.
  updateCheckboxInput(session, "trim_outliers", value = FALSE)
  shinyjs::reset("filters_form")
  session$onFlushed(function() {
    sd_skip_arming(TRUE)
    apply_filters(vals)
    if (was_view_SD) updateCheckboxInput(session, "view_SD", value = TRUE)
    # Update meta_full_rng to the full-range values NOW so the all.equal bypass
    # in metadata_filtered_df passes when the plot re-evaluates after the delay.
    meta_full_rng$efo   <- vals$efo_range
    meta_full_rng$cfr   <- vals$cfr_range
    meta_full_rng$tim   <- vals$tim_range
    meta_full_rng$dcr   <- vals$dcr_range
    meta_full_rng$eco   <- vals$eco_range
    meta_full_rng$dwell <- vals$dwell_range
    meta_full_rng$fluo  <- vals$fluo_range
    # Fire plot_trigger after a short delay to let the browser finish processing
    # all the widget updates sent by apply_filters above.  Without this delay,
    # input$tid_length_range / input$efo_range etc. still hold old filtered
    # values when the plot re-evaluates.
    later::later(function() {
      sd_skip_arming(FALSE)
      if (!is.null(isolate(input$file))) plot_trigger(Sys.time())
    }, delay = 0.5)
  }, once = TRUE)
  shinyjs::runjs("const el = document.getElementById('load_settings'); if (el) { el.value=''; $(el).trigger('change'); }")
  showNotification("Filters reset to data ranges.", type = "message", duration = 3)
})
# When the SD panel opens its renderUI creates fresh widgets, which fire their
# initial values as input changes. Suppress those for 800 ms so they are not
# mistaken for a user-driven narrowing of the filter.
observeEvent(input$view_SD, {
  if (isTRUE(input$view_SD)) {
    sd_skip_arming(TRUE)
    later::later(function() sd_skip_arming(FALSE), delay = 0.8)
  }
}, ignoreInit = TRUE)
# Commit the user's chosen range into sd_filt - the single source of truth.
# Uses Inf for any axis whose widget is at the full-range maximum so that
# spread_filtered_df always passes those traces without a fragile value comparison.
observeEvent(list(input$sd_x_range, input$sd_y_range, input$sd_z_range), {
  req(input$view_SD)
  if (isolate(sd_skip_arming())) return()
  sdr   <- isolate(sd_grey_tbl())
  get_mx <- function(col) { v <- sdr[[col]]; v <- v[is.finite(v)]; if (!length(v)) Inf else max(v) }
  # True when the widget value is at (or indistinguishably close to) the data maximum.
  # A tiny relative tolerance handles JSON float round-trip precision loss, which can
  # cause rx[2] to come back as mx * (1 - ~1e-15), making the strict >= check fail and
  # setting sd_filt to a finite value that then excludes the highest-spread traces.
  at_max <- function(val, mx) is.infinite(val) || (is.finite(mx) && val >= mx - 1e-9 * max(1, abs(mx)))
  # Clearing a box makes numericRangeInput report NA, not NULL, so `%||%` never
  # sees it and the NA propagates through at_max() into `if`, which stops the
  # session with "missing value where TRUE/FALSE is required". An empty box is
  # read as "no bound on that side": blank lower = 0, blank upper = no upper
  # limit, which is the same thing the widget means when it sits at the data
  # maximum. A reversed pair is swapped rather than refused, so dragging the
  # upper below the lower cannot wedge the panel either.
  san <- function(v) {
    v  <- suppressWarnings(as.numeric(v))
    lo <- if (length(v) >= 1L && is.finite(v[1])) max(0, v[1]) else 0
    hi <- if (length(v) >= 2L && !is.na(v[2]))    v[2]         else Inf
    if (is.finite(hi) && hi < lo) { tmp <- lo; lo <- hi; hi <- tmp }
    c(lo, hi)
  }
  mx_x  <- get_mx("sd_x"); mx_y <- get_mx("sd_y")
  rx <- san(input$sd_x_range)
  ry <- san(input$sd_y_range)
  rz <- san(input$sd_z_range)
  sd_filt$x <- if (rx[1] <= 0 && at_max(rx[2], mx_x)) c(0, Inf) else rx
  sd_filt$y <- if (ry[1] <= 0 && at_max(ry[2], mx_y)) c(0, Inf) else ry
  if ("sd_z" %in% names(sdr)) {
    mx_z <- get_mx("sd_z")
    sd_filt$z <- if (rz[1] <= 0 && at_max(rz[2], mx_z)) c(0, Inf) else rz
  }
}, ignoreInit = TRUE)
observeEvent(input$reset_cluster2_settings, {
  apply_cluster2_settings(session, cluster2_defaults())
  
  # Force fileInput to be recreated
  cluster2_file_reset_key(cluster2_file_reset_key() + 1)
  
  showNotification("Median Cluster settings reset (and file cleared).", type = "message", duration = 3)
})

# Uncheck "within cluster2" when cluster2 DBSCAN is disabled
observeEvent(input$enable_cluster2, {
  if (!isTRUE(input$enable_cluster2) && isTRUE(input$nn_within_cluster2)) {
    updateCheckboxInput(session, "nn_within_cluster2", value = FALSE)
  }
}, ignoreInit = TRUE)

observeEvent(input$reset_two_color_settings, {
  apply_two_color_settings(session, two_color_defaults())
  
  # Force file input to be recreated
  two_color_file_reset_key(two_color_file_reset_key() + 1)
  
  showNotification("2-color settings reset to defaults.", type = "message", duration = 3)
})


observeEvent(input$load_overlay_settings, {
  req(input$load_overlay_settings)
  ov <- jsonlite::fromJSON(input$load_overlay_settings$datapath)
  cur_name <- if (!is.null(input$bg_image)) input$bg_image$name else NULL
  apply_overlay_settings(ov, session, current_bg_name = cur_name)
  showNotification("Overlay settings applied.", type = "message", duration = 4)
})
observeEvent(input$reset_overlay, {
  ov <- overlay_defaults()
  apply_overlay_settings(ov, session, current_bg_name = NULL)
  shinyjs::reset("overlay_form")
  showNotification("Overlay reset to defaults.", type = "message", duration = 3)
})
