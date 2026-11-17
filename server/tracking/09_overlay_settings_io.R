## ─────────────────────────────────────────────────────────────────────────────
## server/tracking/09_overlay_settings_io.R  -  image overlay; settings files; filtered CSV
## ─────────────────────────────────────────────────────────────────────────────
## Part of tracking_server(). Not a standalone script: its statements are evaluated in
## the server function's environment, in file-name order (server/parts.R).
##
## The confocal overlay with auto-alignment, and saving/loading the filter,
## playback and overlay settings.
##
## INTERFACE (generated: Rscript tools/server_map.R --headers; do not edit by hand)
##   DEFINES for other parts:
##     trk_image_overlay
##   READS from other parts:
##     01_state_sigma: trk_raw_stats_rv, trk_rechunk_stats_rv,
##                     trk_spatial_filter
##     03_playback_state_zcal: trk_z_applied
##     04_datasets_upload: trk_view_base
##     05_filter_ui: reset_ranges_trigger
##     06_filter_pipeline: .FILT_IDS, trk_filt, trk_spatial_filtered
##     10_triggers: data_ver, fire_frames, fire_layout, layout_trigger
##     11_data_frames: trk
## END INTERFACE
## ─────────────────────────────────────────────────────────────────────────────

# ---- Background image overlay (tracking tab) --------------------------------
# All IDs are trk_-prefixed; logic mirrors the SMLM tab exactly.

trk_confocal_meta <- reactive({
  file <- input$trk_bg_image
  if (is.null(file) || !isTRUE(input$trk_show_bg)) return(NULL)
  img <- tryCatch(magick::image_read(file$datapath), error = function(e) NULL)
  if (is.null(img)) return(NULL)
  rot <- suppressWarnings(as.integer(input$trk_img_rotation %||% 0))
  if (is.finite(rot) && (rot %% 360) != 0) img <- magick::image_rotate(img, rot)
  if (isTRUE(input$trk_img_flip_x)) img <- magick::image_flop(img)
  img     <- apply_background_lut(img, input$trk_img_lut %||% "original",
                                  custom_color = input$trk_img_lut_custom,
                                  brightness   = input$trk_img_brightness %||% 100,
                                  contrast     = input$trk_img_contrast   %||% 0)
  info    <- magick::image_info(img)
  raw_png <- magick::image_write(img, format = "png")
  list(
    uri       = paste0("data:image/png;base64,", base64enc::base64encode(raw_png)),
    width_px  = info$width,
    height_px = info$height
  )
})

# ── trk_image_overlay() - place a confocal image behind the tracks ──────────
# GOAL: turn an uploaded confocal/widefield image into a plotly background
# layer positioned in DATA COORDINATES (nm), so tracks overlay the cell image.
#
# The image is in pixels; the tracks are in nm. Converting between them needs
# the pixel size and an origin offset, both user-supplied:
#   nm = pixel_index * trk_px_size + (trk_img_x0 / trk_img_y0)
# returns NULL (draw nothing) if there is no image or the pixel size is unset -
# every caller treats NULL as "no overlay".
trk_image_overlay <- reactive({
  meta <- trk_confocal_meta()                                  # image + its dimensions
  if (is.null(meta)) return(NULL)                              # no image uploaded
  ps <- suppressWarnings(as.numeric(input$trk_px_size))
  if (!is.finite(ps) || ps <= 0) return(NULL)
  x0 <- suppressWarnings(as.numeric(input$trk_img_x0)); if (!is.finite(x0)) x0 <- 0
  y0 <- suppressWarnings(as.numeric(input$trk_img_y0)); if (!is.finite(y0)) y0 <- 0
  sizex <- meta$width_px  * ps
  sizey <- meta$height_px * ps
  if (isTRUE(input$trk_img_invert_y)) {
    y_anchor <- "top";    y_val <- y0 + sizey
  } else {
    y_anchor <- "bottom"; y_val <- y0
  }
  list(
    source  = meta$uri,
    xref = "x", yref = "y",
    x = x0,     y = y_val,
    sizex = sizex, sizey = sizey,
    xanchor = "left", yanchor = y_anchor,
    sizing  = "stretch",
    opacity = if (is.null(input$trk_img_opacity)) 0.6 else input$trk_img_opacity,
    layer   = "below"
  )
})

# ---- Filter settings save/load ----
output$trk_download_filter_settings <- downloadHandler(
  filename = function() {
    base <- trk_view_base()
    paste0(base, "_filter_settings_", format(Sys.time(), "%Y%m%d-%H%M%S"), ".json")
  },
  content = function(file) {
    sf <- trk_spatial_filter()
    s <- list(
      # ── Every setting that changes a Tracking result (R/trk_settings.R) ─────
      # The loader restores this block when present. The keys below it stay for
      # files read by older versions of the app.
      settings_version = 2L,
      settings_type    = "minfluxr_tracking",
      trk_settings     = unclass(isolate(trk_settings_from_input(input))),
      # ── Metadata (per-localization) ranges ──────────────────────────────
      # The stored filters (null = no filter), not the widgets' clamped display.
      efo_rng   = trk_filt$efo %||% NA,   cfr_rng   = trk_filt$cfr %||% NA,
      dcr_rng   = trk_filt$dcr %||% NA,   tim_rng   = trk_filt$tim %||% NA,
      eco_rng   = trk_filt$eco %||% NA,   dwell_rng = trk_filt$dwell %||% NA,
      fluo_rng  = trk_filt$fluo %||% NA,
      trk_tid_len    = trk_filt$tid_len %||% NA,
      trk_dur_filter = trk_filt$dur %||% NA,

      # ── Track CLEANING - these change the tracks themselves ─────────────
      # Cleaning runs before every metadata filter (split → trim → gate), so
      # omitting it meant a reloaded settings file reproduced the ranges but
      # silently analysed uncleaned tracks: different track IDs, lengths and
      # durations, hence different MSD/D. These matter more than the ranges.
      trk_split_gaps_enable   = isTRUE(input$trk_split_gaps_enable),
      trk_split_adaptive      = isTRUE(input$trk_split_adaptive),
      trk_split_gap_cycles    = input$trk_split_gap_cycles,
      trk_split_min_len       = input$trk_split_min_len,

      trk_endtrim_enable      = isTRUE(input$trk_endtrim_enable),
      trk_endtrim_dark        = isTRUE(input$trk_endtrim_dark),
      trk_endtrim_dark_cycles = input$trk_endtrim_dark_cycles,
      trk_endtrim_efo         = isTRUE(input$trk_endtrim_efo),
      trk_endtrim_efo_frac    = input$trk_endtrim_efo_frac,

      trk_quality_enable      = isTRUE(input$trk_quality_enable),
      trk_quality_darktol     = input$trk_quality_darktol,
      trk_quality_min_integ   = input$trk_quality_min_integ,

      # ── Spatial crop (also a data filter) ───────────────────────────────
      trk_crop_pct   = input$trk_crop_pct,
      trk_spatial_x  = if (!is.null(sf)) unname(sf$x) else NULL,
      trk_spatial_y  = if (!is.null(sf)) unname(sf$y) else NULL,
      trk_spatial_pct_thr = if (!is.null(sf)) sf$pct_thr else NULL,

      # ── Histogram panel display state ───────────────────────────────────
      trk_show_metadata_hist = isTRUE(input$trk_show_metadata_hist),
      trk_hist_param         = input$trk_hist_param,
      trk_hist_bins          = input$trk_hist_bins,
      trk_hist_show_raw      = isTRUE(input$trk_hist_show_raw)
    )
    writeLines(jsonlite::toJSON(s, pretty = TRUE, auto_unbox = TRUE, na = "null"), file)
  }
)
observeEvent(input$trk_load_filter_settings, {
  req(input$trk_load_filter_settings$datapath)
  s <- jsonlite::fromJSON(input$trk_load_filter_settings$datapath)

  # Every restore is guarded on the key being PRESENT, so an older settings file
  # (saved before cleaning was persisted) leaves the current cleaning setup alone
  # instead of silently switching split/trim/gate off.
  num <- function(v) suppressWarnings(as.numeric(unlist(v)))
  upd_num <- function(id) {
    v <- s[[id]]
    if (!is.null(v)) { n <- num(v)[1]; if (is.finite(n)) updateNumericInput(session, id, value = n) }
  }
  upd_sld <- function(id) {
    v <- s[[id]]
    if (!is.null(v)) { n <- num(v)[1]; if (is.finite(n)) updateSliderInput(session, id, value = n) }
  }
  upd_chk <- function(id) if (!is.null(s[[id]])) updateCheckboxInput(session, id, value = isTRUE(s[[id]]))

  # Range filters go straight into the store (a missing key leaves a filter as it
  # is, null clears it); the widgets then re-render from it.
  for (k in names(.FILT_IDS)) {
    id <- .FILT_IDS[[k]]
    if (!id %in% names(s)) next
    v <- num(s[[id]])
    trk_filt[[k]] <- if (length(v) == 2L && all(is.finite(v))) sort(v) else NULL
  }
  reset_ranges_trigger(isolate(reset_ranges_trigger()) + 1L)

  # cleaning: enables first, then their parameters
  for (id in c("trk_split_gaps_enable","trk_split_adaptive",
               "trk_endtrim_enable","trk_endtrim_dark","trk_endtrim_efo",
               "trk_quality_enable",
               "trk_show_metadata_hist","trk_hist_show_raw"))
    upd_chk(id)
  for (id in c("trk_split_gap_cycles","trk_split_min_len",
               "trk_endtrim_dark_cycles","trk_endtrim_efo_frac","trk_crop_pct"))
    upd_num(id)
  for (id in c("trk_quality_darktol","trk_quality_min_integ","trk_hist_bins"))
    upd_sld(id)
  if (!is.null(s$trk_hist_param))
    updateSelectInput(session, "trk_hist_param", selected = as.character(s$trk_hist_param))

  # Spatial crop: restore the stored box directly (it is a reactiveVal, not an
  # input) and bump the data version so the filter is actually re-applied.
  if (!is.null(s$trk_spatial_x) && !is.null(s$trk_spatial_y) &&
      length(unlist(s$trk_spatial_x)) == 2 && length(unlist(s$trk_spatial_y)) == 2) {
    pct <- num(s$trk_spatial_pct_thr)[1]
    if (!is.finite(pct) || pct <= 0 || pct > 1) pct <- 0.5
    trk_spatial_filter(list(x = num(s$trk_spatial_x), y = num(s$trk_spatial_y), pct_thr = pct))
    data_ver(isolate(data_ver()) + 1L)
    fire_frames()
  }

  # ── The settings block (files from v0.7.1 on): every data-affecting setting ──
  # Restored last, so where a key appears both here and at the top level the
  # block wins. The calibration profile is not reselected: selecting it copies
  # the profile's prefactor, pattern, L and wander floor into those widgets,
  # which could overwrite the exact values the file holds (the profile may have
  # been re-measured since). Those values are restored directly instead.
  blk <- trk_settings_from_file(s)
  n_restored <- 0L
  if (!is.null(blk)) {
    plan <- trk_settings_restore_plan(blk$settings[setdiff(blk$present, "trk_cal_profile")])
    for (i in seq_len(nrow(plan))) {
      id <- plan$id[i]; v <- plan$value[[i]]
      switch(plan$widget[i],
        numeric  = updateNumericInput(session, id, value = v),
        checkbox = updateCheckboxInput(session, id, value = isTRUE(v)),
        select   = updateSelectInput(session, id, selected = as.character(v)),
        radio    = updateRadioButtons(session, id, selected = as.character(v)),
        slider   = updateSliderInput(session, id, value = v),
        text     = updateTextInput(session, id, value = as.character(v)))
    }
    n_restored <- nrow(plan)
  }
  prof <- if (!is.null(blk)) blk$settings$trk_cal_profile else NULL

  cleaning_on <- any(vapply(c("trk_split_gaps_enable","trk_endtrim_enable","trk_quality_enable"),
                            function(k) isTRUE(s[[k]]), logical(1)))
  showNotification(
    paste0(if (n_restored) sprintf("Settings loaded (%d analysis settings restored).", n_restored)
           else "Filter settings loaded (a file from before v0.7.1: filters and cleaning only).",
           if (!is.null(prof) && nzchar(prof))
             sprintf(" Saved with calibration profile '%s'; its values were restored, the profile itself is not reselected.", prof)
           else "",
           if (cleaning_on) " Cleaning settings restored - press Plot Data to re-apply them." else "",
           if (!is.null(s$trk_spatial_x)) " Spatial crop restored." else ""),
    type = "message", duration = 6)
})

# ---- Playback settings save/load ----
output$trk_download_playback_settings <- downloadHandler(
  filename = function() {
    base <- trk_view_base()
    paste0(base, "_playback_settings_", format(Sys.time(), "%Y%m%d-%H%M%S"), ".json")
  },
  content = function(file) {
    s <- list(
      trk_play_mode      = input$trk_play_mode,
      trk_play_total     = input$trk_play_total,
      trk_bin_dt         = input$trk_bin_dt,
      trk_play_speed     = input$trk_play_speed,
      trk_target_markers = input$trk_target_markers,
      trk_trail_len      = input$trk_trail_len,
      trk_per_tid_cap    = input$trk_per_tid_cap,
      trk_trail_nmax     = input$trk_trail_nmax,
      trk_max_pts_tid    = input$trk_max_pts_tid,
      # Panel-visibility and trail toggles. The advanced numerics above were
      # already saved but their gating checkboxes were not, so a reloaded file
      # restored values the user could not see and trails stayed off.
      trk_show_playback_params = isTRUE(input$trk_show_playback_params),
      trk_playback_ctrl        = isTRUE(input$trk_playback_ctrl),
      trk_show_moving_trails   = isTRUE(input$trk_show_moving_trails),
      trk_show_advanced        = isTRUE(input$trk_show_advanced),
      # View mode lives just above the playback panel and decides whether any of
      # the settings below are even used (playback vs static full tracks), so it
      # belongs with them rather than in the filter file.
      trk_view_mode_main       = input$trk_view_mode_main,
      trk_full_connect         = isTRUE(input$trk_full_connect)
    )
    writeLines(jsonlite::toJSON(s, pretty = TRUE, auto_unbox = TRUE, na = "null"), file)
  }
)
observeEvent(input$trk_load_playback_settings, {
  req(input$trk_load_playback_settings$datapath)
  s <- jsonlite::fromJSON(input$trk_load_playback_settings$datapath)
  upd_num <- function(id) {
    v <- s[[id]]
    if (!is.null(v) && is.finite(as.numeric(v)))
      updateNumericInput(session, id, value = as.numeric(v))
  }
  if (!is.null(s$trk_play_mode))
    updateSelectInput(session, "trk_play_mode", selected = s$trk_play_mode)
  if (!is.null(s$trk_play_speed) && is.finite(as.numeric(s$trk_play_speed)))
    updateSliderInput(session, "trk_play_speed", value = as.numeric(s$trk_play_speed))
  for (id in c("trk_play_total","trk_bin_dt","trk_target_markers",
               "trk_trail_len",
               "trk_per_tid_cap","trk_trail_nmax","trk_max_pts_tid"))
    upd_num(id)
  # Gating checkboxes - restore only when present, so an older playback file
  # does not collapse the panels the user currently has open.
  for (id in c("trk_show_playback_params","trk_playback_ctrl",
               "trk_show_moving_trails","trk_show_advanced","trk_full_connect"))
    if (!is.null(s[[id]])) updateCheckboxInput(session, id, value = isTRUE(s[[id]]))
  if (!is.null(s$trk_view_mode_main))
    updateRadioButtons(session, "trk_view_mode_main", selected = as.character(s$trk_view_mode_main))
  showNotification("Playback settings loaded.", type = "message", duration = 3)
})

# Live guidance for the drift-speed window W. W must sit well BELOW a typical track
# duration: once W >= a track's duration every point in it shares the same window
# (lo=1, hi=n), so the whole track renders as ONE flat colour and the within-track
# gradient the mode exists for disappears. Suggest ~1/10 of the median duration, which
# leaves ~10 independent window positions along a typical track.
output$trk_speed_win_note <- renderUI({
  if (!identical(input$trk_color_mode, "speed")) return(NULL)
  rs <- trk_rechunk_stats_rv(); ro <- trk_raw_stats_rv()
  dur <- if (!is.null(rs) && length(rs$tid_duration)) rs$tid_duration
         else if (!is.null(ro) && length(ro$tid_duration)) ro$tid_duration else NULL
  dur <- dur[is.finite(dur) & dur > 0]
  if (!length(dur)) return(NULL)
  med  <- stats::median(dur)
  W    <- suppressWarnings(as.numeric(input$trk_color_speed_win))
  sugg <- signif(med / 10, 2)
  flat <- if (is.finite(W)) mean(dur <= W) else 0
  warn <- is.finite(W) && W > med / 3
  div(style = sprintf("color:%s; font-size:0.8em; line-height:1.35; margin-top:-6px;",
                      if (warn) "#b8860b" else "#666"),
      if (warn) "⚠ " else "",
      sprintf("Median track duration %s ms. ", format(signif(med, 3), trim = TRUE)),
      if (flat > 0) sprintf(
        "At W = %s ms, %.0f%% of tracks are shorter than the window - those render as one flat colour. ",
        format(W, trim = TRUE), 100 * flat) else "",
      sprintf("Suggested W ≈ %s ms (~1/10 of a typical track).", format(sugg, trim = TRUE))
  )
})

# Warn when the trail is shorter than one bin: each frame advances by one bin width,
# so a trail shorter than a bin can't bridge consecutive frames - it skips ahead and
# looks jumpy. Smooth, continuous trails need trail length >= bin width.
output$trk_trail_warn <- renderUI({
  if (!isTRUE(input$trk_show_moving_trails)) return(NULL)
  tl <- suppressWarnings(as.numeric(input$trk_trail_len))
  bw <- suppressWarnings(as.numeric(input$trk_bin_dt))
  if (!is.finite(tl) || !is.finite(bw) || bw <= 0 || tl >= bw) return(NULL)
  div(style = "color:#b8860b; font-size:12px; line-height:1.3; margin:-4px 0 6px 0;",
      sprintf(paste0("⚠ Trail length (%g ms) is below the bin width (%g ms). The trail ",
                     "advances one bin per frame, so a trail shorter than a bin skips ahead ",
                     "each frame and looks jumpy. Set trail length ≥ bin width for smooth trails."),
              tl, bw))
})

output$trk_download_overlay_settings <- downloadHandler(
  filename = function() {
    base <- trk_view_base()
    paste0(base, "_image_alignment_", format(Sys.time(), "%Y%m%d-%H%M%S"), ".json")
  },
  content  = function(file) {
    ov <- list(
      show_confocal = input$trk_show_confocal,
      show_bg       = input$trk_show_bg,
      bg_image_name = if (!is.null(input$trk_bg_image)) input$trk_bg_image$name else NULL,
      px_size       = input$trk_px_size,
      img_x0        = input$trk_img_x0,
      img_y0        = input$trk_img_y0,
      img_invert_y  = input$trk_img_invert_y,
      img_flip_x    = input$trk_img_flip_x,
      img_opacity   = input$trk_img_opacity,
      img_rotation  = input$trk_img_rotation,
      img_lut       = input$trk_img_lut,
      img_lut_custom = input$trk_img_lut_custom,
      img_brightness = input$trk_img_brightness,
      img_contrast   = input$trk_img_contrast
    )
    writeLines(jsonlite::toJSON(ov, pretty = TRUE, auto_unbox = TRUE, na = "null"), file)
  }
)

observeEvent(input$trk_load_overlay_settings, {
  req(input$trk_load_overlay_settings)
  ov <- tryCatch(jsonlite::fromJSON(input$trk_load_overlay_settings$datapath),
                 error = function(e) NULL)
  if (is.null(ov)) return()
  if (!is.null(ov$show_confocal)) updateCheckboxInput(session, "trk_show_confocal", value = ov$show_confocal)
  if (!is.null(ov$show_bg))       updateCheckboxInput(session, "trk_show_bg",       value = ov$show_bg)
  if (!is.null(ov$px_size))       updateNumericInput(session,  "trk_px_size",        value = ov$px_size)
  if (!is.null(ov$img_x0))        updateNumericInput(session,  "trk_img_x0",         value = ov$img_x0)
  if (!is.null(ov$img_y0))        updateNumericInput(session,  "trk_img_y0",         value = ov$img_y0)
  if (!is.null(ov$img_invert_y))  updateCheckboxInput(session, "trk_img_invert_y",  value = ov$img_invert_y)
  if (!is.null(ov$img_flip_x))    updateCheckboxInput(session, "trk_img_flip_x",    value = ov$img_flip_x)
  if (!is.null(ov$img_opacity))   updateSliderInput(session,   "trk_img_opacity",   value = ov$img_opacity)
  if (!is.null(ov$img_rotation))  updateSelectInput(session,   "trk_img_rotation",  selected = as.character(ov$img_rotation))
  if (!is.null(ov$img_lut))       updateSelectInput(session,   "trk_img_lut",       selected = as.character(ov$img_lut))
  if (!is.null(ov$img_lut_custom)) colourpicker::updateColourInput(session, "trk_img_lut_custom", value = as.character(ov$img_lut_custom))
  if (!is.null(ov$img_brightness)) updateSliderInput(session,   "trk_img_brightness", value = ov$img_brightness)
  if (!is.null(ov$img_contrast))   updateSliderInput(session,   "trk_img_contrast",   value = ov$img_contrast)
  if (!is.null(ov$bg_image_name)) {
    cur <- if (!is.null(input$trk_bg_image)) input$trk_bg_image$name else NULL
    if (is.null(cur) || !identical(cur, ov$bg_image_name))
      showNotification(paste0("Re-upload '", ov$bg_image_name, "' to reproduce alignment."),
                       type = "message", duration = 6)
  }
})

output$trk_download_filtered_csv <- downloadHandler(
  filename = function() {
    base <- trk_view_base()
    paste0(base, "_filtered_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
  },
  content = function(file) {
    df <- trk_spatial_filtered()
    shiny::validate(shiny::need(nrow(df) > 0, "No filtered data to export."))
    drop_cols <- intersect(c("col", "t"), names(df))
    if (length(drop_cols)) df <- df[, !names(df) %in% drop_cols, drop = FALSE]
    # tim is stored internally in ms; convert back to seconds for export
    if ("tim" %in% names(df)) df$tim <- df$tim / 1000
    # Axial provenance: z below is at this factor. Re-importing (here or in the
    # SMLM tab) divides it out instead of falling back to the 0.7 guess in
    # detect_csv_z_calibration(), so a round-trip is exact.
    df$z_scale_applied <- trk_z_applied()
    write.csv(df, file, row.names = FALSE)
  }
)
observeEvent(input$trk_download_filtered_csv, {
  showNotification("Exported filtered tracking data.", type = "message", duration = 3)
})

observeEvent(input$trk_reset_overlay, {
  updateCheckboxInput(session, "trk_show_confocal", value = FALSE)
  updateCheckboxInput(session, "trk_show_bg",       value = FALSE)
  updateNumericInput(session,  "trk_px_size",        value = 10)
  updateNumericInput(session,  "trk_img_x0",         value = 0)
  updateNumericInput(session,  "trk_img_y0",         value = 0)
  updateCheckboxInput(session, "trk_img_invert_y",  value = TRUE)
  updateCheckboxInput(session, "trk_img_flip_x",    value = FALSE)
  updateSliderInput(session,   "trk_img_opacity",   value = 0.6)
  updateSelectInput(session,   "trk_img_rotation",  selected = "0")
  updateSelectInput(session,   "trk_img_lut",       selected = "original")
  colourpicker::updateColourInput(session, "trk_img_lut_custom", value = "#00FF00")
  updateSliderInput(session,   "trk_img_brightness", value = 100)
  updateSliderInput(session,   "trk_img_contrast",   value = 0)
  shinyjs::reset("trk_overlay_form")
})

# Overlay changes: push via proxy relayout (no DOM/frames rebuild needed).
observeEvent(
  list(input$trk_show_bg, input$trk_bg_image, input$trk_px_size,
       input$trk_img_x0, input$trk_img_y0, input$trk_img_invert_y,
       input$trk_img_flip_x, input$trk_img_opacity, input$trk_img_rotation,
       input$trk_img_lut, input$trk_img_lut_custom,
       input$trk_img_brightness, input$trk_img_contrast),
  {
    req(!is.na(layout_trigger()))   # plot must already exist
    ovl <- tryCatch(isolate(trk_image_overlay()), error = function(e) NULL)
    plotlyProxyInvoke(
      plotlyProxy("trk_track2D", session),
      "relayout",
      list(images = if (!is.null(ovl)) list(ovl) else list())
    )
  },
  ignoreInit = TRUE
)

# ---- Auto-align overlay helpers (shared image-loading logic) ----
.load_trk_img <- function() {
  f <- input$trk_bg_image
  if (is.null(f)) return(NULL)
  img <- tryCatch(magick::image_read(f$datapath), error = function(e) NULL)
  if (is.null(img)) return(NULL)
  rot <- suppressWarnings(as.integer(input$trk_img_rotation %||% 0))
  if (is.finite(rot) && (rot %% 360) != 0) img <- magick::image_rotate(img, rot)
  if (isTRUE(input$trk_img_flip_x)) img <- magick::image_flop(img)
  img
}

observeEvent(input$trk_auto_foci, {
  img <- .load_trk_img(); req(!is.null(img))
  df  <- isolate(trk()); req(nrow(df) > 0)
  ps  <- suppressWarnings(as.numeric(input$trk_px_size))
  req(is.finite(ps), ps > 0)
  res <- tryCatch(
    auto_align_overlay(img, df, ps, mode = "foci",
                       invert_y    = isTRUE(input$trk_img_invert_y),
                       n_peaks        = as.integer(input$trk_foci_n_peaks        %||% 5L),
                       bright_pct     = (input$trk_foci_bright_pct     %||% 95)  / 100,
                       border_pct     = (input$trk_foci_border_pct     %||% 10)  / 100,
                       eps_pct        = (input$trk_foci_eps_pct        %||% 3)   / 100,
                       aspect_max     =  input$trk_foci_aspect_max     %||% 4.0,
                       density_thr    = (input$trk_foci_density_thr    %||% 20)  / 100,
                       density_cell   =  input$trk_foci_density_cell   %||% 100,
                       density_smooth =  input$trk_foci_density_smooth %||% 100),
    error = function(e) {
      showNotification(paste("Foci alignment error:", conditionMessage(e)),
                       type = "error", duration = 8)
      NULL
    }
  )
  if (is.null(res)) {
    showNotification("Foci alignment: could not detect matching peaks. Check image and pixel size.",
                     type = "warning", duration = 6)
    return()
  }
  # Store diagnostic peaks (img_peaks are in nm assuming x0=0; add x0/y0 for plot coords)
  n_blobs    <- if (!is.null(res$img_peaks))  nrow(res$img_peaks)  else 0L
  n_clusters <- if (!is.null(res$data_peaks)) nrow(res$data_peaks) else 0L
  if (n_blobs > 0L) {
    rv$trk_diag <- list(
      data_clusters = res$data_peaks,
      img_blobs     = data.frame(x = res$x0 + res$img_peaks$x,
                                 y = res$y0 + res$img_peaks$y)
    )
  }
  updateCheckboxInput(session, "trk_show_bg", value = TRUE)
  updateNumericInput(session, "trk_img_x0", value = res$x0)
  updateNumericInput(session, "trk_img_y0", value = res$y0)
  if (!is.na(layout_trigger())) fire_layout("foci_align")
  showNotification(
    paste0("Foci alignment applied (", n_blobs, " confocal blobs, ",
           n_clusters, " MINFLUX clusters detected). ",
           "Enable 'Show detection diagnostics' to visualise."),
    type = "message", duration = 6)
}, ignoreInit = TRUE)

observeEvent(input$trk_show_diag, {
  if (!is.na(layout_trigger())) fire_layout("show_diag")
}, ignoreInit = TRUE)

observeEvent(input$trk_auto_boundary, {
  img <- .load_trk_img(); req(!is.null(img))
  df  <- isolate(trk()); req(nrow(df) > 0)
  ps  <- suppressWarnings(as.numeric(input$trk_px_size))
  req(is.finite(ps), ps > 0)
  res <- tryCatch(
    auto_align_overlay(img, df, ps, mode = "boundary",
                       invert_y = isTRUE(input$trk_img_invert_y)),
    error = function(e) {
      showNotification(paste("Boundary alignment error:", conditionMessage(e)),
                       type = "error", duration = 8)
      NULL
    }
  )
  if (is.null(res)) {
    showNotification("Boundary alignment: could not compute centroid. Check image and pixel size.",
                     type = "warning", duration = 6)
    return()
  }
  updateCheckboxInput(session, "trk_show_bg", value = TRUE)
  updateNumericInput(session, "trk_img_x0", value = res$x0)
  updateNumericInput(session, "trk_img_y0", value = res$y0)
  showNotification("Boundary alignment applied - nudge to fine-tune.",
                   type = "message", duration = 4)
}, ignoreInit = TRUE)

# ---- Nudge buttons for X and Y origin ----
local({
  deltas <- c(m1000 = -1000, m100 = -100, m10 = -10,
              p10   =    10, p100 =  100, p1000 = 1000)
  for (sfx in names(deltas)) {
    local({
      s <- sfx          # bind into this env so the reactive expression
      d <- deltas[[s]]  # resolves the correct string after loop ends
      observeEvent(input[[paste0("trk_nx_", s)]], {
        v <- suppressWarnings(as.numeric(input$trk_img_x0))
        updateNumericInput(session, "trk_img_x0",
                           value = (if (is.finite(v)) v else 0) + d)
      }, ignoreInit = TRUE)
      observeEvent(input[[paste0("trk_ny_", s)]], {
        v <- suppressWarnings(as.numeric(input$trk_img_y0))
        updateNumericInput(session, "trk_img_y0",
                           value = (if (is.finite(v)) v else 0) + d)
      }, ignoreInit = TRUE)
    })
  }
})


