## ─────────────────────────────────────────────────────────────────────────────
## server/smlm/04_settings_collect_apply.R  -  collect the settings for saving; defaults and apply functions
## ─────────────────────────────────────────────────────────────────────────────
## Part of smlm_server(). Not a standalone script: its statements are evaluated in
## the server function's environment, in file-name order (server/parts.R).
##
## get_ui_settings (what the settings file holds), the defaults of each
## settings group, and the functions that push a loaded or default settings
## list back into the widgets.
##
## INTERFACE (generated: Rscript tools/server_map.R --headers; do not edit by hand)
##   DEFINES for other parts:
##     apply_cluster2_settings, apply_filters, apply_heatmap, apply_overlay,
##     apply_overlay_settings, apply_two_color_settings, cluster2_defaults,
##     filters_defaults_from_data, get_ui_settings, heatmap_defaults,
##     overlay_defaults, two_color_defaults
##   READS from other parts:
##     01_view_upload: smlm_2d_crop_rv, trim_applied_rv, uploaded_df_rv
##     03_dataset_zcal: commit_z_scale
##     05_filter_pipeline: sd_filt
##     11_dcr_split: dcr_note_written
## END INTERFACE
## ─────────────────────────────────────────────────────────────────────────────

# --------- 3) COLLECT CURRENT UI SETTINGS ----------
get_ui_settings <- reactive({
  req(uploaded_df_rv())
  list(
    # Every setting that changes results, complete by construction
    # (R/smlm_settings.R). The loader reads this block first; the fields below
    # keep the older layout readable and carry everything that is not a setting
    # (display state, crop region, overlay, QC thresholds).
    settings_version = SMLM_SETTINGS_FILE_VERSION,
    smlm_settings    = unclass(smlm_settings_from_input(input)),

    efo_range        = input$efo_range,
    cfr_range        = input$cfr_range,
    tim_range        = input$tim_range,
    dcr_range        = input$dcr_range,
    eco_range        = input$eco_range,
    tid_length_range = input$tid_length_range,
    dwell_range      = input$dwell_range,
    fluo_range       = input$fluo_range,
    # SD filter - save from sd_filt (source of truth)
    # which is unreliable when the SD panel is closed (renderUI removes the widgets).
    sd_x_range       = sd_filt$x,
    sd_y_range       = sd_filt$y,
    sd_z_range       = sd_filt$z,
    sd_min_locs      = input$sd_min_locs,
    # Outlier trimming REMOVES LOCALIZATIONS, so it changes every per-trace
    # statistic downstream. Omitting it would make a reloaded file silently
    # analyse a different point set than the one the ranges were set on.
    trim_outliers    = isTRUE(input$trim_outliers),
    trim_k           = input$trim_k,
    make_positive    = input$make_positive,

    # ── Histogram + spread-panel display state ─────────────────────────────
    # apply_filters() already restores hist_variable/hist_bins/sd_* , but they
    # were never written out, so every load reset them to the data defaults.
    show_metadata_hist      = isTRUE(input$show_metadata_hist),
    hist_variable           = input$hist_variable,
    hist_bins               = input$hist_bins,
    hist_show_raw           = isTRUE(input$hist_show_raw),
    view_SD                 = isTRUE(input$view_SD),
    sd_bins                 = input$sd_bins,
    sd_hist_shared          = isTRUE(input$sd_hist_shared),
    sd_hist_range           = input$sd_hist_range,
    sd_overlay_all          = isTRUE(input$sd_overlay_all),
    sd_show_mean_sd         = isTRUE(input$sd_show_mean_sd),
    sd_show_overall_mean_sd = isTRUE(input$sd_show_overall_mean_sd),
    sd_use_q95              = isTRUE(input$sd_use_q95),
    target_sd_x             = input$target_sd_x,
    target_sd_y             = input$target_sd_y,
    target_sd_z             = input$target_sd_z,

    z_scale             = input$z_scale,

    enable_clustering   = input$enable_clustering,
    deblink_method      = input$deblink_method,
    eps_centroid        = input$eps_centroid,
    vote_frac           = input$vote_frac,
    vote_rho            = input$vote_rho,
    vote_link           = input$vote_link,
    centroid_k          = input$centroid_k,
    # part of the merge metric: eps only means the same thing alongside it
    deblink_anisotropic = input$deblink_anisotropic,
    # Channel handling. cluster_group_scope decides WHICH channels a run
    # covers; deblink_shared_spread decides whether the merge distance is
    # measured per channel or once over all of them. Both change what a
    # de-blinked molecule IS, so neither may be left to a default on reload.
    cluster_group_scope   = input$cluster_group_scope,
    deblink_shared_spread = isTRUE(input$deblink_shared_spread),
    # minPts is derived from tid_length_range and exclude_unclustered no longer
    # exists, so neither is saved - both would be stale on reload.
    
    enable_cluster2  = input$enable_cluster2,
    eps2             = input$eps2,
    minpts2          = input$minpts2,
    
    two_color = list(
      enable        = input$enable_two_color,
      use_fluo      = input$color_by_fluo,
      dcr_group1    = input$dcr_group1,
      dcr_group2    = input$dcr_group2,
      color_group1  = input$color_group1,
      color_group2  = input$color_group2
    ),
    
    # background overlay settings 
    overlay = list(
      show_confocal = input$show_confocal,
      show_bg       = input$show_bg,
      px_size   = input$px_size,
      img_x0    = input$img_x0,
      img_y0    = input$img_y0,
      img_invert_y = input$img_invert_y,
      img_flip_x   = input$img_flip_x,
      img_opacity  = input$img_opacity,
      img_rotation = input$img_rotation,
      img_lut      = input$img_lut,
      img_lut_custom = input$img_lut_custom,
      img_brightness = input$img_brightness,
      img_contrast   = input$img_contrast
    ),
    # heatmap settings
    heatmap = list(
      show_heatmap      = input$show_heatmap,
      heatmap_n         = input$heatmap_n,
      heatmap_opacity   = input$heatmap_opacity,
      heatmap_palette   = input$heatmap_palette,
      heatmap_bw_factor = input$heatmap_bw_factor,
      heatmap_bw_abs    = input$heatmap_bw_abs,
      heatmap_floor     = input$heatmap_floor,
      heatmap_use_raw   = input$heatmap_use_raw,
      heatmap_clip      = input$heatmap_clip,
      heatmap_gamma     = input$heatmap_gamma,
      heatmap_log       = input$heatmap_log,
      heatmap_showlegend= input$heatmap_showlegend
    ),

    # ── Spatial crop / trimmer - CHANGES THE DATA (the point set) ──────────
    # Both the 2D "Crop to view" box and the 3D box/sphere crop
    # remove localizations, so a reload without them analyses a wider region
    # than the ranges were set on. Saved so a round-trip keeps the same crop.
    crop = {
      c2d <- smlm_2d_crop_rv()
      list(
        crop2d_x       = if (!is.null(c2d)) c2d$x else NULL,
        crop2d_y       = if (!is.null(c2d)) c2d$y else NULL,
        enable_trimmer = isTRUE(input$enable_trimmer),
        trim_mode      = input$trim_mode,
        trim_x_range   = input$trim_x_range,
        trim_y_range   = input$trim_y_range,
        trim_z_range   = input$trim_z_range,
        trim_cx = input$trim_cx, trim_cy = input$trim_cy,
        trim_cz = input$trim_cz, trim_r  = input$trim_r,
        # the committed region is what actually crops (inputs above are the editor)
        trim_applied   = trim_applied_rv()
      )
    },

    # ── Fiducial-scan exclusion - CHANGES THE DATA when on ─────────────────
    # The thresholds are saved too so re-running the scan reproduces the flags.
    fiducial = list(
      fid_exclude   = isTRUE(input$fid_exclude),
      fid_efo_z     = input$fid_efo_z,
      fid_dur_x     = input$fid_dur_x,
      fid_min_locs  = input$fid_min_locs,
      fid_len_z     = input$fid_len_z,
      fid_recur_r   = input$fid_recur_r,
      fid_recur_min = input$fid_recur_min,
      fid_all_files = isTRUE(input$fid_all_files)
    ),

    # ── Analysis options that change results ───────────────────────────────
    analysis = list(
      use_pooled_dcr       = isTRUE(input$use_pooled_dcr),
      spatial_edge_correct = isTRUE(input$spatial_edge_correct),
      spatial_use_support  = isTRUE(input$spatial_use_support),
      dens_support_from    = input$dens_support_from,
      dens_norm_enable     = isTRUE(input$dens_norm_enable),
      dens_thin_enable     = isTRUE(input$dens_thin_enable)
    ),

    # ── Local density & hotspots ────────────────────────────────────
    # Every one of these changes a reported number, not just a picture: k and the
    # quadrat size set the scale the density is measured at, the footprint pair
    # sets the area every statistic divides by, and the eps calibration decides
    # what counts as a focus. A reload that silently reverted any of them would
    # reproduce the same plot from different data.
    hotspot = list(
      hot_unit        = input$hot_unit,
      hot_use_3d      = isTRUE(input$hot_use_3d),
      hot_k           = input$hot_k,
      hot_enh_thresh  = input$hot_enh_thresh,
      hot_quadrat     = input$hot_quadrat,
      hot_minpts      = input$hot_minpts,
      hot_target_frac = input$hot_target_frac,
      hot_nsim        = input$hot_nsim,
      hot_palette     = input$hot_palette,
      hot_show_foci   = isTRUE(input$hot_show_foci),
      hot_box_top     = input$hot_box_top,
      hot_metric      = input$hot_metric,
      hot_curve_mode  = input$hot_curve_mode
    )
  )
})
# --------- 4) DEFAULTS saved + APPLY functions  ----------
filters_defaults_from_data <- function(df) {
  rng_of <- function(v) { v <- v[is.finite(v)]; if (!length(v)) c(NA,NA) else range(v, na.rm=TRUE) }
  tid_rng <- {
    if ("tid" %in% names(df)) {
      tl <- as.numeric(table(df$tid)); if (length(tl)) range(tl, na.rm=TRUE) else c(NA,NA)
    } else c(NA,NA)
  }
  # SD defaults (computed on uploaded data)
  # Lower bound is always 0 to match the widget initial value; the bypass check
  # in spread_filtered_df compares input values against meta_full_rng which is
  # also initialized with c(0, max), so they must be consistent.
  max_or1 <- function(v) { v <- v[is.finite(v)]; if (!length(v)) return(1); mx <- max(v, na.rm = TRUE); if (!is.finite(mx) || mx <= 0) 1 else mx }
  sda <- sd_by_tid(df, min_n = input$sd_min_locs)
  rx  <- if (nrow(sda)) c(0, max_or1(sda$sd_x)) else c(NA, NA)
  ry  <- if (nrow(sda)) c(0, max_or1(sda$sd_y)) else c(NA, NA)
  rz  <- if (nrow(sda) && "sd_z" %in% names(sda) && any(is.finite(sda$sd_z))) c(0, max_or1(sda$sd_z)) else c(NA, NA)
  
  list(
    hist_variable="efo", hist_bins=200,
    efo_range = if ("efo"%in%names(df)) rng_of(df$efo) else c(NA,NA),
    cfr_range = if ("cfr"%in%names(df)) rng_of(df$cfr) else c(NA,NA),
    tim_range = if ("tim"%in%names(df)) rng_of(df$tim) else c(NA,NA),
    dcr_range = if ("dcr"%in%names(df)) rng_of(df$dcr) else c(NA,NA),
    eco_range = if ("eco"%in%names(df)) rng_of(df$eco) else c(NA,NA),
    dwell_range = if ("dwell"%in%names(df)) rng_of(df$dwell) else c(NA,NA),
    fluo_range  = if ("fluo"%in%names(df)) rng_of(df$fluo) else c(NA,NA),
    tid_length_range = tid_rng,
    #SD summary controls
    target_sd_x=NA_real_, target_sd_y=NA_real_, target_sd_z=NA_real_,
    sd_bins = 100,
    sd_overlay_all = TRUE,
    sd_show_mean_sd = TRUE,
    sd_show_overall_mean_sd = FALSE,
    
    # --- SD filter ranges ---
    sd_x_range = rx, sd_y_range = ry, sd_z_range = rz
  )
}
overlay_defaults <- reactive({
  list(show_confocal=FALSE, show_bg=FALSE, px_size=10, img_x0=0, img_y0=0,
       img_invert_y=TRUE, img_flip_x=FALSE, img_opacity=0.6, img_rotation="0",
       img_lut="original", img_lut_custom="#00FF00", img_brightness=100,
       img_contrast=0, bg_image_name=NULL)
})
heatmap_defaults <- reactive({
  list(show_heatmap=FALSE, heatmap_n=200, heatmap_opacity=0.35, heatmap_palette="YlOrRd",
       heatmap_bw_factor=1, heatmap_bw_abs=NA_real_, heatmap_floor=0, heatmap_use_raw=FALSE,
       heatmap_clip=0.99, heatmap_gamma=1, heatmap_log=FALSE, heatmap_showlegend=FALSE)
})
cluster2_defaults <- function() {
  # Data-affecting values come from SMLM_DEFAULTS (R/smlm_settings.R), the same
  # values a fresh session shows, so Reset cannot disagree with a new session.
  c(SMLM_DEFAULTS[c("cluster2_method", "cluster2_xi", "cluster2_prefilter",
                    "cluster2_prefilter_r", "cluster2_prefilter_min", "cluster2_prune",
                    "cluster2_prune_tol", "cluster2_core_only", "cluster2_min_members",
                    "eps2", "minpts2")],
  list(
    
    # 3D wrappers
    show_cluster2_shape = FALSE,
    show_cluster2_isosurface = FALSE,
    cluster2_iso_selected_only = TRUE,
    cluster2_iso_id = 1,
    cluster2_iso_color = "#00CED1",
    cluster2_iso_unicolor = FALSE,
    cluster2_iso_by_channel = TRUE,
    cluster2_iso_grid_res = 22,
    cluster2_iso_pad = 0,
    cluster2_iso_manual_offset = FALSE,
    cluster2_iso_offset = 30,
    cluster2_iso_opacity = 0.25,
    cluster2_iso_manual_sigma = FALSE,
    cluster2_iso_smooth_rel = 1,
    cluster2_iso_sigma = 30,
    
    # 2D wrapper
    show_cluster2_shape_2D = FALSE,
    
    # metrics toggles
    show_cluster2_metrics = FALSE,
    
    # 2D metrics
    show_cluster2_pca_axes = FALSE,   # note: you reuse same id in 2D & 3D panels
    show_fit_distances_2d = FALSE,
    fit_color_2d = "#00BFFF",
    fit_line_width_2d = 3,
    
    # 3D metrics
    show_cluster2_labels = FALSE,
    show_cluster2_plane_angle = FALSE,
    show_cluster2_plane_distance = FALSE,
    show_bending_radius = FALSE,
    fit_method = "heuristic",
    fit_weighted = FALSE,
    mle_center_mode = "free",
    mle_center_lambda = 0.1,
    fit_type = "none",
    helix_turns = 3,
    show_fit_tube = FALSE,
    fit_tube_radius = 5,
    show_fit_distances = FALSE,
    fit_opacity = 0.4,
    
    # NN within cluster2
    show_cluster2_nn_lines = FALSE,
    nn_k_cluster2 = 3,
    distance_range_cluster2 = c(0, 1250)
  ))
}
two_color_defaults <- function() {
  list(
    color_by_fluo = SMLM_DEFAULTS$color_by_fluo,
    dcr_split_manual = SMLM_DEFAULTS$dcr_split_manual,
    dcr_group1 = SMLM_DEFAULTS$dcr_group1,
    dcr_group2 = SMLM_DEFAULTS$dcr_group2,
    
    color_group1 = SMLM_DISPLAY_DEFAULTS$color_group1,
    color_group2 = SMLM_DISPLAY_DEFAULTS$color_group2,
    
    
    # CLUSTERING PARAMETERS ARE NOT HERE ANY MORE. There is one clustering
    # pipeline, so there is one set of de-blink / median-clustering controls,
    # and they live in the main settings file. This file now carries only what
    # DEFINES the channels and how they are drawn. Old files that still hold
    # the retired *_tc1 / *_tc2 keys are migrated on load (see below).
    show_tc1_alpha = FALSE,
    tc1_alpha = 120,
    tc1_alpha_opacity = 0.15,
    tc1_alpha_min_pts = 5
    # The Median Cluster wrap and isosurface are no longer duplicated here;
    # they are saved with the Median Clustering controls. Old files carrying
    # show_tc2_alpha / show_tc2_isosurface are migrated on load.
  )
}

apply_two_color_settings <- function(session, s) {
  getv <- function(name, default = NULL) if (!is.null(s[[name]])) s[[name]] else default
  
  # checkboxes
  cb_ids <- c(
    # master toggle first, so the restored 2-color parameters take effect
    "enable_two_color",
    "color_by_fluo", "use_pooled_dcr", "dcr_split_manual",
    "show_tc1_alpha"
  )
  for (id in cb_ids) if (!is.null(s[[id]])) updateCheckboxInput(session, id, value = isTRUE(s[[id]]))

  # These widget changes come from a file, so they must not read as the user
  # typing a range or switching the DCR column (see dcr_track).
  if (!is.null(s[["use_pooled_dcr"]])) dcr_track$pooled <- isTRUE(s[["use_pooled_dcr"]])
  dcr_note_written(getv("dcr_group1"), getv("dcr_group2"))

  # numeric range inputs (shinyWidgets)
  if (!is.null(getv("dcr_group1")) && length(getv("dcr_group1")) == 2)
    shinyWidgets::updateNumericRangeInput(session, "dcr_group1", value = as.numeric(getv("dcr_group1")))
  if (!is.null(getv("dcr_group2")) && length(getv("dcr_group2")) == 2)
    shinyWidgets::updateNumericRangeInput(session, "dcr_group2", value = as.numeric(getv("dcr_group2")))
  
  # colors
  if (!is.null(getv("color_group1"))) colourpicker::updateColourInput(session, "color_group1", value = as.character(getv("color_group1")))
  if (!is.null(getv("color_group2"))) colourpicker::updateColourInput(session, "color_group2", value = as.character(getv("color_group2")))
  
  # NN
  # Retired 2-colour NN controls migrate onto the shared NN panel: there is one
  # NN section now, and enabling channels adapts it rather than adding a second.
  if (!is.null(getv("nn_k_two_color")) && is.null(s[["nn_k"]]))
    updateNumericInput(session, "nn_k", value = as.integer(getv("nn_k_two_color")))
  if (!is.null(getv("tc_label_font_size")) && is.null(s[["nn_label_font_size"]]))
    updateNumericInput(session, "nn_label_font_size", value = as.numeric(getv("tc_label_font_size")))
  if (!is.null(getv("tc_label_offset")) && is.null(s[["nn_label_offset"]]))
    updateNumericInput(session, "nn_label_offset", value = as.numeric(getv("tc_label_offset")))
  rng <- getv("distance_range_two_color")
  if (!is.null(rng) && length(rng) == 2 && is.null(s[["distance_range"]]))
    updateSliderInput(session, "distance_range", value = as.numeric(rng))
  # "between" / "within" map onto the shared three-valued channel scope.
  if (!is.null(getv("two_color_scope")) && is.null(s[["nn_mplex_scope"]]))
    updateRadioButtons(session, "nn_mplex_scope",
      selected = if (identical(as.character(getv("two_color_scope")), "between")) "across" else "within")
  if (!is.null(getv("draw_nn_two_color")) && isTRUE(getv("draw_nn_two_color")))
    updateCheckboxInput(session, "show_nn_lines", value = TRUE)
  
  # ── LEGACY MIGRATION: retired *_tc1 / *_tc2 keys ────────────────────────
  # Settings files written before the pipelines were unified carry a second,
  # 2-colour-only copy of every clustering parameter. There is one pipeline
  # now, so those values are applied to the single set of controls rather than
  # dropped - silently ignoring them would reload a saved analysis with default
  # radii, which is the failure mode this migration exists to prevent.
  #
  # Only fills a control the new file did not already set, so a settings file
  # containing both old and new keys keeps the new ones.
  tc_legacy <- list(
    eps_tc1                    = "eps_centroid",
    deblink_method_tc1         = "deblink_method",
    vote_frac_tc1              = "vote_frac",
    vote_rho_tc1               = "vote_rho",
    eps_tc2                    = "eps2",
    minpts_tc2                 = "minpts2",
    cluster2_method_tc2        = "cluster2_method",
    cluster2_xi_tc2            = "cluster2_xi",
    cluster2_prefilter_tc2     = "cluster2_prefilter",
    cluster2_prefilter_r_tc2   = "cluster2_prefilter_r",
    cluster2_prefilter_min_tc2 = "cluster2_prefilter_min",
    cluster2_prune_tc2         = "cluster2_prune",
    cluster2_prune_tol_tc2     = "cluster2_prune_tol",
    cluster2_core_only_tc2     = "cluster2_core_only",
    cluster2_min_members_tc2   = "cluster2_min_members",
    enable_two_color_cluster1  = "enable_clustering",
    enable_two_color_cluster2  = "enable_cluster2"
  )
  migrated <- character(0)
  for (old_id in names(tc_legacy)) {
    v <- getv(old_id)
    if (is.null(v)) next
    new_id <- tc_legacy[[old_id]]
    if (!is.null(s[[new_id]])) next          # new key wins
    if (is.logical(v) || old_id %in% c("enable_two_color_cluster1", "enable_two_color_cluster2",
                                       "cluster2_prefilter_tc2", "cluster2_prune_tc2",
                                       "cluster2_core_only_tc2")) {
      updateCheckboxInput(session, new_id, value = isTRUE(v))
    } else if (is.character(v)) {
      updateRadioButtons(session, new_id, selected = as.character(v))
    } else {
      updateNumericInput(session, new_id, value = as.numeric(v))
    }
    migrated <- c(migrated, old_id)
  }
  if (length(migrated))
    showNotification(sprintf(
      "Loaded an older 2-colour settings file: %d retired clustering setting(s) were applied to the unified de-blink / Median Clustering controls.",
      length(migrated)), type = "message", duration = 8)

  # 2-colour display settings (these are still 2-colour-specific)
  if (!is.null(getv("tc1_alpha"))) updateNumericInput(session, "tc1_alpha", value = as.numeric(getv("tc1_alpha")))
  if (!is.null(getv("tc1_alpha_opacity"))) updateSliderInput(session, "tc1_alpha_opacity", value = as.numeric(getv("tc1_alpha_opacity")))
  if (!is.null(getv("tc1_alpha_min_pts"))) updateNumericInput(session, "tc1_alpha_min_pts", value = as.integer(getv("tc1_alpha_min_pts")))
  
  # RETIRED: the two-colour Median Cluster wrap and isosurface. They now live
  # under Median Clustering and draw per channel. An older file that had them
  # switched on turns on the unified controls instead, carrying the isosurface
  # geometry across so the surface looks as it did.
  tc2_moved <- character(0)
  if (isTRUE(getv("show_tc2_alpha"))) {
    updateCheckboxInput(session, "show_cluster2_shape", value = TRUE)
    tc2_moved <- c(tc2_moved, "alpha-shape wrap")
  }
  if (isTRUE(getv("show_tc2_isosurface"))) {
    updateCheckboxInput(session, "show_cluster2_isosurface", value = TRUE)
    updateCheckboxInput(session, "cluster2_iso_by_channel", value = TRUE)
    for (pr in list(c("tc2_iso_grid_res",   "cluster2_iso_grid_res"),
                    c("tc2_iso_pad",        "cluster2_iso_pad"),
                    c("tc2_iso_smooth_rel", "cluster2_iso_smooth_rel"),
                    c("tc2_iso_offset",     "cluster2_iso_offset"),
                    c("tc2_iso_sigma",      "cluster2_iso_sigma"),
                    c("tc2_iso_opacity",    "cluster2_iso_opacity"))) {
      v <- suppressWarnings(as.numeric(getv(pr[1])))
      if (length(v) == 1L && is.finite(v)) updateSliderInput(session, pr[2], value = v)
    }
    for (pr in list(c("tc2_iso_manual_offset", "cluster2_iso_manual_offset"),
                    c("tc2_iso_manual_sigma",  "cluster2_iso_manual_sigma")))
      if (!is.null(getv(pr[1])))
        updateCheckboxInput(session, pr[2], value = isTRUE(getv(pr[1])))
    tc2_moved <- c(tc2_moved, "isosurface")
  }
  if (length(tc2_moved))
    showNotification(sprintf(
      "This settings file used the two-colour Median Cluster %s. That now lives under Median Clustering, drawn separately for each channel, and has been switched on there.",
      paste(tc2_moved, collapse = " and ")), type = "message", duration = 10)
}



apply_overlay_settings <- function(ov, session, current_bg_name=NULL) {  # (same as your 2D)
  if (!is.null(ov$show_confocal)) updateCheckboxInput(session,"show_confocal", value=ov$show_confocal)
  if (!is.null(ov$show_bg))       updateCheckboxInput(session,"show_bg",       value=ov$show_bg)
  if (!is.null(ov$px_size))      updateNumericInput(session,"px_size", value=ov$px_size)
  if (!is.null(ov$img_x0))       updateNumericInput(session,"img_x0",  value=ov$img_x0)
  if (!is.null(ov$img_y0))       updateNumericInput(session,"img_y0",  value=ov$img_y0)
  if (!is.null(ov$img_invert_y)) updateCheckboxInput(session,"img_invert_y", value=ov$img_invert_y)
  if (!is.null(ov$img_flip_x))   updateCheckboxInput(session,"img_flip_x",   value=ov$img_flip_x)
  if (!is.null(ov$img_opacity))  updateSliderInput(session,"img_opacity",    value=ov$img_opacity)
  if (!is.null(ov$img_rotation)) updateSelectInput(session,"img_rotation",   selected=as.character(ov$img_rotation))
  if (!is.null(ov$img_lut))      updateSelectInput(session,"img_lut",        selected=as.character(ov$img_lut))
  if (!is.null(ov$img_lut_custom)) colourpicker::updateColourInput(session,"img_lut_custom", value=as.character(ov$img_lut_custom))
  if (!is.null(ov$img_brightness)) updateSliderInput(session,"img_brightness", value=ov$img_brightness)
  if (!is.null(ov$img_contrast))   updateSliderInput(session,"img_contrast",   value=ov$img_contrast)
  if (!is.null(ov$bg_image_name)) {
    cur <- if (!is.null(input$bg_image)) input$bg_image$name else NULL
    if (is.null(cur) || !identical(cur, ov$bg_image_name)) {
      showNotification(paste0("Overlay loaded. Re-upload '", ov$bg_image_name, "' to reproduce alignment."),
                       type="message", duration=6)
    }
  }
}
apply_filters <- function(vals) {
  # jsonlite::fromJSON may parse JSON arrays as lists rather than numeric vectors;
  # coerce all range values defensively so updateNumericRangeInput never receives a list.
  rng <- function(v) if (is.null(v)) NULL else as.numeric(unlist(v))
  updateSelectInput(session,"hist_variable", selected=vals$hist_variable)
  updateSliderInput(session,"hist_bins", value=vals$hist_bins)
  shinyWidgets::updateNumericRangeInput(session,"efo_range",  value=rng(vals$efo_range))
  shinyWidgets::updateNumericRangeInput(session,"cfr_range",  value=rng(vals$cfr_range))
  shinyWidgets::updateNumericRangeInput(session,"tim_range",  value=rng(vals$tim_range))
  shinyWidgets::updateNumericRangeInput(session,"dcr_range",  value=rng(vals$dcr_range))
  shinyWidgets::updateNumericRangeInput(session,"eco_range",  value=rng(vals$eco_range))
  shinyWidgets::updateNumericRangeInput(session,"dwell_range",value=rng(vals$dwell_range))
  shinyWidgets::updateNumericRangeInput(session,"fluo_range", value=rng(vals$fluo_range))
  shinyWidgets::updateNumericRangeInput(session,"tid_length_range", value=rng(vals$tid_length_range))
  updateNumericInput(session,"target_sd_x", value=vals$target_sd_x)
  updateNumericInput(session,"target_sd_y", value=vals$target_sd_y)
  updateNumericInput(session,"target_sd_z", value=vals$target_sd_z)
  updateSliderInput(session, "sd_bins", value = vals$sd_bins)
  updateCheckboxInput(session, "sd_hist_shared", value = isTRUE(vals$sd_hist_shared))
  if (!is.null(vals$sd_hist_range))
    shinyWidgets::updateNumericRangeInput(session, "sd_hist_range",
                                          value = rng(vals$sd_hist_range))
  updateCheckboxInput(session, "sd_overlay_all", value = vals$sd_overlay_all)
  updateCheckboxInput(session, "sd_show_mean_sd", value = vals$sd_show_mean_sd %||% TRUE)
  updateCheckboxInput(session, "sd_show_overall_mean_sd", value = vals$sd_show_overall_mean_sd %||% FALSE)
  if (!is.null(vals$sd_x_range)) shinyWidgets::updateNumericRangeInput(session, "sd_x_range", value = rng(vals$sd_x_range))
  if (!is.null(vals$sd_y_range)) shinyWidgets::updateNumericRangeInput(session, "sd_y_range", value = rng(vals$sd_y_range))
  if (!is.null(vals$sd_z_range)) shinyWidgets::updateNumericRangeInput(session, "sd_z_range", value = rng(vals$sd_z_range))

  # ── Panel display state ───────────────────────────────────────────────────
  # Only restore when the key is actually present, so an OLD settings file (saved
  # before these were persisted) leaves the current UI alone.
  # (Settings files written up to v0.7.0 may carry remove_rogues / rogue_rule /
  # rogue_factor. Rogue removal was taken out of the app in 2026-11; those keys
  # are ignored, and the file loads as before in every other respect.)
  if (!is.null(vals$show_metadata_hist)) updateCheckboxInput(session, "show_metadata_hist", value = isTRUE(vals$show_metadata_hist))
  if (!is.null(vals$hist_show_raw))      updateCheckboxInput(session, "hist_show_raw", value = isTRUE(vals$hist_show_raw))
  if (!is.null(vals$sd_use_q95))         updateCheckboxInput(session, "sd_use_q95", value = isTRUE(vals$sd_use_q95))
}
apply_overlay  <- function(vals) { apply_overlay_settings(vals, session) }
apply_heatmap  <- function(vals) {
  updateCheckboxInput(session,"show_heatmap", value=vals$show_heatmap)
  updateNumericInput(session,"heatmap_n", value=vals$heatmap_n)
  updateSliderInput(session,"heatmap_opacity", value=vals$heatmap_opacity)
  updateSelectInput(session,"heatmap_palette", selected=vals$heatmap_palette)
  updateSliderInput(session,"heatmap_bw_factor", value=vals$heatmap_bw_factor)
  updateNumericInput(session,"heatmap_bw_abs", value=vals$heatmap_bw_abs)
  updateSliderInput(session,"heatmap_floor", value=vals$heatmap_floor)
  updateCheckboxInput(session,"heatmap_use_raw", value=vals$heatmap_use_raw)
  updateSliderInput(session,"heatmap_clip", value=vals$heatmap_clip)
  updateSliderInput(session,"heatmap_gamma", value=vals$heatmap_gamma)
  updateCheckboxInput(session,"heatmap_log", value=vals$heatmap_log)
  updateCheckboxInput(session,"heatmap_showlegend", value=vals$heatmap_showlegend)
}
apply_cluster2_settings <- function(session, s) {
  getv <- function(name, default = NULL) {
    if (!is.null(s[[name]])) s[[name]] else default
  }
  
  # method + numeric inputs
  if (!is.null(getv("cluster2_method"))) updateRadioButtons(session, "cluster2_method", selected = as.character(getv("cluster2_method")))
  if (!is.null(getv("cluster2_xi"))) updateNumericInput(session, "cluster2_xi", value = as.numeric(getv("cluster2_xi")))
  if (!is.null(getv("cluster2_prefilter"))) updateCheckboxInput(session, "cluster2_prefilter", value = isTRUE(getv("cluster2_prefilter")))
  if (!is.null(getv("cluster2_prefilter_r"))) updateNumericInput(session, "cluster2_prefilter_r", value = as.numeric(getv("cluster2_prefilter_r")))
  if (!is.null(getv("cluster2_prefilter_min"))) updateNumericInput(session, "cluster2_prefilter_min", value = as.numeric(getv("cluster2_prefilter_min")))
  if (!is.null(getv("cluster2_prune"))) updateCheckboxInput(session, "cluster2_prune", value = isTRUE(getv("cluster2_prune")))
  if (!is.null(getv("cluster2_prune_tol"))) updateNumericInput(session, "cluster2_prune_tol", value = as.numeric(getv("cluster2_prune_tol")))
  if (!is.null(getv("cluster2_core_only"))) updateCheckboxInput(session, "cluster2_core_only", value = isTRUE(getv("cluster2_core_only")))
  if (!is.null(getv("cluster2_min_members"))) updateNumericInput(session, "cluster2_min_members", value = as.numeric(getv("cluster2_min_members")))
  if (!is.null(getv("eps2")))    updateNumericInput(session, "eps2", value = as.numeric(getv("eps2")))
  if (!is.null(getv("minpts2"))) updateNumericInput(session, "minpts2", value = as.numeric(getv("minpts2")))

  # ── PRIMARY (de-blinking) settings ──
  # These were written to the settings file but never restored, so a reloaded
  # analysis silently ran with the default radius. minpts / exclude_unclustered
  # are deliberately absent: minPts now derives from tid_length_range, and
  # exclude_unclustered no longer exists (de-blinking never removes molecules).
  for (nmv in c("vote_frac", "vote_rho"))
    if (!is.null(getv(nmv))) updateNumericInput(session, nmv, value = as.numeric(getv(nmv)))
  if (!is.null(getv("vote_link")))
    updateRadioButtons(session, "vote_link", selected = as.character(getv("vote_link")))
  if (!is.null(getv("eps_centroid")))
    updateNumericInput(session, "eps_centroid", value = as.numeric(getv("eps_centroid")))
  if (!is.null(getv("centroid_k")))
    updateNumericInput(session, "centroid_k", value = as.numeric(getv("centroid_k")))
  # Restored before enable_clustering, so the criterion is in place when it fires.
  if (!is.null(getv("deblink_method")))
    updateRadioButtons(session, "deblink_method", selected = as.character(getv("deblink_method")))
  if (!is.null(getv("deblink_anisotropic")))
    updateCheckboxInput(session, "deblink_anisotropic", value = isTRUE(getv("deblink_anisotropic")))
  if (!is.null(getv("cluster_group_scope")))
    updateSelectInput(session, "cluster_group_scope", selected = as.character(getv("cluster_group_scope")))
  if (!is.null(getv("deblink_shared_spread")))
    updateCheckboxInput(session, "deblink_shared_spread", value = isTRUE(getv("deblink_shared_spread")))
  # Restored last so the parameters above are in place before it takes effect.
  if (!is.null(getv("enable_clustering")))
    updateCheckboxInput(session, "enable_clustering", value = isTRUE(getv("enable_clustering")))

  # Axial calibration: rescales every z distance, so it must survive a reload.
  if (!is.null(getv("z_scale"))) {
    zs <- as.numeric(getv("z_scale"))
    updateNumericInput(session, "z_scale", value = zs)
    # Loading settings is itself an explicit action, so commit rather than
    # leaving the restored calibration pending behind a Plot Data press.
    commit_z_scale(zs, reason = "Axial calibration from settings")
  }


  # checkboxes
  cb_ids <- c(
    # master toggle first, so the restored parameters actually take effect
    "enable_cluster2", "show_cluster2_knn",
    "show_cluster2_shape", "show_cluster2_isosurface", "cluster2_iso_selected_only",
    "cluster2_iso_by_channel",
    "show_cluster2_shape_2D", "show_cluster2_metrics",
    "show_cluster2_pca_axes",
    "show_cluster2_labels", "show_cluster2_plane_angle", "show_cluster2_plane_distance",
    "show_bending_radius", "fit_weighted",
    "show_fit_distances", "show_fit_tube",
    "show_cluster2_nn_lines",
    "compute_rmsd_metrics", "membrane_protein_options", "show_c2_nn_labels"
  )
  for (id in cb_ids) {
    if (!is.null(s[[id]])) updateCheckboxInput(session, id, value = isTRUE(s[[id]]))
  }
  
  # selects / radios
  if (!is.null(getv("fit_method")))     updateRadioButtons(session, "fit_method", selected = as.character(getv("fit_method")))
  if (!is.null(getv("mle_center_mode")))updateRadioButtons(session, "mle_center_mode", selected = as.character(getv("mle_center_mode")))
  if (!is.null(getv("fit_type")))       updateSelectInput(session, "fit_type", selected = as.character(getv("fit_type")))
  # sliders / numeric misc
  if (!is.null(getv("mle_center_lambda"))) updateSliderInput(session, "mle_center_lambda", value = as.numeric(getv("mle_center_lambda")))
  if (!is.null(getv("helix_turns")))       updateSliderInput(session, "helix_turns", value = as.numeric(getv("helix_turns")))
  if (!is.null(getv("fit_tube_radius")))   updateNumericInput(session, "fit_tube_radius", value = as.numeric(getv("fit_tube_radius")))
  if (!is.null(getv("fit_opacity")))       updateNumericInput(session, "fit_opacity", value = as.numeric(getv("fit_opacity")))
  
  # isosurface controls
  if (!is.null(getv("cluster2_iso_id")))       updateNumericInput(session, "cluster2_iso_id", value = as.integer(getv("cluster2_iso_id")))
  if (!is.null(getv("cluster2_iso_grid_res"))) updateSliderInput(session, "cluster2_iso_grid_res", value = as.numeric(getv("cluster2_iso_grid_res")))
  if (!is.null(getv("cluster2_iso_pad")))      updateSliderInput(session, "cluster2_iso_pad", value = as.numeric(getv("cluster2_iso_pad")))
  if (!is.null(getv("cluster2_iso_offset")))   updateSliderInput(session, "cluster2_iso_offset", value = as.numeric(getv("cluster2_iso_offset")))
  if (!is.null(getv("cluster2_iso_manual_offset"))) updateCheckboxInput(session, "cluster2_iso_manual_offset", value = isTRUE(getv("cluster2_iso_manual_offset")))
  if (!is.null(getv("cluster2_iso_opacity")))  updateSliderInput(session, "cluster2_iso_opacity", value = as.numeric(getv("cluster2_iso_opacity")))
  if (!is.null(getv("cluster2_iso_sigma")))    updateSliderInput(session, "cluster2_iso_sigma", value = as.numeric(getv("cluster2_iso_sigma")))
  if (!is.null(getv("cluster2_iso_manual_sigma"))) updateCheckboxInput(session, "cluster2_iso_manual_sigma", value = isTRUE(getv("cluster2_iso_manual_sigma")))
  if (!is.null(getv("cluster2_iso_smooth_rel"))) updateSliderInput(session, "cluster2_iso_smooth_rel", value = as.numeric(getv("cluster2_iso_smooth_rel")))
  if (!is.null(getv("cluster2_iso_color")))    colourpicker::updateColourInput(session, "cluster2_iso_color", value = as.character(getv("cluster2_iso_color")))
  if (!is.null(getv("cluster2_iso_unicolor"))) updateCheckboxInput(session, "cluster2_iso_unicolor", value = isTRUE(getv("cluster2_iso_unicolor")))
  if (!is.null(getv("cluster2_iso_by_channel"))) updateCheckboxInput(session, "cluster2_iso_by_channel", value = isTRUE(getv("cluster2_iso_by_channel")))
  
  # color
  if (!is.null(getv("fit_color_2d"))) colourpicker::updateColourInput(session, "fit_color_2d", value = as.character(getv("fit_color_2d")))
  
  # fit/metric options added after this block was last updated
  if (!is.null(getv("arc_span_deg")))      updateSliderInput(session, "arc_span_deg", value = as.numeric(getv("arc_span_deg")))
  if (!is.null(getv("max_dist_lines_nm"))) updateNumericInput(session, "max_dist_lines_nm", value = as.numeric(getv("max_dist_lines_nm")))

  # NN controls
  if (!is.null(getv("nn_k_cluster2"))) updateNumericInput(session, "nn_k_cluster2", value = as.integer(getv("nn_k_cluster2")))
  if (!is.null(getv("c2_label_font_size"))) updateNumericInput(session, "c2_label_font_size", value = as.numeric(getv("c2_label_font_size")))
  if (!is.null(getv("c2_label_offset")))    updateNumericInput(session, "c2_label_offset", value = as.numeric(getv("c2_label_offset")))
  rng <- getv("distance_range_cluster2")
  if (!is.null(rng) && length(rng) == 2) {
    updateSliderInput(session, "distance_range_cluster2", value = as.numeric(rng))
  }
}
