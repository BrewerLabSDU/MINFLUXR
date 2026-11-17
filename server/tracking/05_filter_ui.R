## ─────────────────────────────────────────────────────────────────────────────
## server/tracking/05_filter_ui.R  -  range-filter widgets
## ─────────────────────────────────────────────────────────────────────────────
## Part of tracking_server(). Not a standalone script: its statements are evaluated in
## the server function's environment, in file-name order (server/parts.R).
##
## The metadata, length and histogram range controls and Reset Ranges.
##
## INTERFACE (generated: Rscript tools/server_map.R --headers; do not edit by hand)
##   DEFINES for other parts:
##     reset_ranges_trigger
##   READS from other parts:
##     01_state_sigma: trk_raw_stats_rv, trk_rechunk_bbox,
##                     trk_rechunk_stats_rv, trk_spatial_filter
##     04_datasets_upload: trk_has_data, trk_raw_df, trk_view_label
##     06_filter_pipeline: .FILT_IDS, .filt_value, trk_meta_filtered
##     07_cycle_time_sequence: tid_summary_rv
##     10_triggers: data_ver
## END INTERFACE
## ─────────────────────────────────────────────────────────────────────────────

# ---- Dynamic range UIs ----
# Bounds come from trk_raw_stats_rv() (frozen at upload, whole dataset) so
# Sticky bar label: filename + current chunk indicator
output$trk_bar_label <- renderUI({
  req(trk_has_data())
  nm <- trk_view_label()
  # Chunks scope the playback animation only, so the label mentions them only in
  # Playback mode - in Full tracks the whole dataset is always shown.
  chunk_info <- if (!identical(input$trk_view_mode_main %||% "full", "playback")) "" else {
    sel <- input$trk_chunk_sel
    if (identical(sel, "__all__"))        " - All chunks"
    else if (isTruthy(sel))              paste0(" - ", basename(sel))
    else                                 ""
  }
  tags$span(style = "font-size:13px; color:#555;", nm, chunk_info)
})

# sliders reflect the full data extent regardless of which chunk is selected.
# Re-renders on: upload (trk_raw_stats_rv() changes) OR Reset Ranges press
# (reset_ranges_trigger increments).  Never re-renders on chunk switches, so
# user-set filter values survive chunk navigation.
make_range_ui <- function(col, label) {
  renderUI({
    reset_ranges_trigger()        # re-render on Reset Ranges press
    rs_orig    <- trk_raw_stats_rv()
    rs_rechunk <- trk_rechunk_stats_rv()
    # Bounds from original upload so the slider never permanently shrinks
    orig_v <- if (!is.null(rs_orig) && col %in% names(rs_orig)) rs_orig[[col]] else {
      df <- trk_raw_df(); if (!col %in% names(df)) return(NULL); df[[col]]
    }
    orig_v <- orig_v[is.finite(orig_v)]; if (!length(orig_v)) return(NULL)
    bounds <- range(orig_v, na.rm = TRUE)
    # Value: the stored filter, clamped to this dataset's range (full range if none).
    val <- .filt_value(col, bounds)
    numericRangeInput(paste0(col, "_rng"), label, value = val, min = bounds[1], max = bounds[2])
  })
}
output$trk_efo_ui   <- make_range_ui("efo",   "EFO range")
output$trk_cfr_ui   <- make_range_ui("cfr",   "CFR range")
output$trk_dcr_ui   <- make_range_ui("dcr",   "DCR range")
output$trk_tim_ui   <- make_range_ui("tim",   "TIM range (ms)")
output$trk_eco_ui   <- make_range_ui("eco",   "ECO range")
output$trk_dwell_ui <- make_range_ui("dwell", "Dwell range")
output$trk_fluo_ui  <- make_range_ui("fluo",  "FLUO range")

# Stage 2 (length) auto-detection. The slider BOUNDS reflect the metadata-filtered
# (post-split) population, but update at PLOT cadence only - keyed on data_ver()
# (bumps on upload / Plot Data / Reset Ranges), NOT live while other sliders move.
# So the range always describes the currently PLOTTED data and stops jumping as you
# drag. The VALUE is PRESERVED (clamped into the new bounds) so a manual minimum
# (e.g. 50) is not wiped; it snaps to full only on new upload / Reset Ranges /
# post-rechunk, because those clear the stored filter (.filt_value then returns the
# full range).
output$trk_tid_len_ui <- renderUI({
  reset_ranges_trigger()
  data_ver()                                                           # update at PLOT cadence (upload / Plot Data), not live
  seg_df <- tryCatch(isolate(trk_meta_filtered()), error = function(e) NULL)  # isolate: no live re-fire while dragging
  lv <- if (!is.null(seg_df) && nrow(seg_df) > 0 && "tid" %in% names(seg_df))
          as.numeric(table(seg_df$tid)) else numeric(0)
  if (!length(lv)) {                          # pipeline not ready → rechunk stats, else upload
    rs_re <- trk_rechunk_stats_rv(); rs_or <- trk_raw_stats_rv()
    lv <- if (!is.null(rs_re) && length(rs_re$tid_length)) rs_re$tid_length
          else if (!is.null(rs_or) && length(rs_or$tid_length)) rs_or$tid_length else numeric(0)
  }
  if (!length(lv)) return(NULL)
  bounds <- range(lv, na.rm = TRUE)
  val <- .filt_value("tid_len", bounds)      # stored filter, clamped (full range if none)
  numericRangeInput("trk_tid_len", "Trace ID length (localizations)",
                    value = val, min = bounds[1], max = bounds[2])
})

# Histogram zoom range - always anchored to the ORIGINAL upload stats so the
# range never shrinks after Filter & Rechunk.  If it used tid_summary_rv() or
# trk_raw_df() (both reflect the current, possibly filtered, manifest), the zoom
# window would clip the grey bars to the filtered range, making them look
# identical to the blue bars.
output$trk_hist_range_ui <- renderUI({
  req(input$trk_hist_param)
  param <- input$trk_hist_param
  rs <- trk_raw_stats_rv()
  if (param == "tid_length") {
    vals <- if (!is.null(rs) && length(rs$tid_length)) rs$tid_length else {
      df <- tryCatch(trk_raw_df(), error = function(e) NULL)
      req(!is.null(df)); as.numeric(table(df$tid))
    }
  } else if (param == "tid_duration") {
    vals <- if (!is.null(rs) && length(rs$tid_duration)) rs$tid_duration else {
      ts <- tid_summary_rv()
      if (!is.null(ts) && nrow(ts) > 0) ts$duration_ms else {
        df <- tryCatch(trk_raw_df(), error = function(e) NULL)
        req(!is.null(df))
        df %>% dplyr::group_by(tid) %>%
          dplyr::summarise(d = max(tim) - min(tim), .groups = "drop") %>% dplyr::pull(d)
      }
    }
  } else {
    vals <- if (!is.null(rs) && param %in% names(rs)) rs[[param]] else {
      df <- tryCatch(trk_raw_df(), error = function(e) NULL)
      req(!is.null(df), param %in% names(df)); df[[param]]
    }
    vals <- vals[is.finite(vals)]
  }
  req(length(vals) > 0)
  rng <- range(vals, na.rm = TRUE)
  numericRangeInput("trk_hist_range", "Histogram zoom range:",
                    value = rng, min = rng[1], max = rng[2])
})

# Reset all filter inputs to the original upload's full range, and restore
# rv$manifest_dir to the original upload so trk_raw_df sees all TIDs again.
# All filter UIs (metadata, tid_len, tid_dur) are renderUI-wrapped - incrementing
# reset_ranges_trigger forces them all to re-render with value = full rng from
# trk_raw_stats_rv().  updateNumericRangeInput does not reliably update
# renderUI-wrapped widgets and is intentionally not used here.
reset_tracking_filters <- function() {
  if (!is.null(rv$original_manifest_dir))
    rv$manifest_dir <- rv$original_manifest_dir
  trk_rechunk_stats_rv(NULL)  # clear rechunk stats so UIs revert to original upload ranges
  trk_rechunk_bbox(NULL)      # clear locked extent so axes revert to per-TID bbox
  trk_spatial_filter(NULL)    # clear any active spatial crop
  for (k in names(.FILT_IDS)) trk_filt[[k]] <- NULL   # no range filters
  reset_ranges_trigger(isolate(reset_ranges_trigger()) + 1L)
}

observeEvent(input$trk_reset_filters, {
  reset_tracking_filters()
}, ignoreInit = TRUE)


# Incrementing this forces trk_dur_filter_ui to re-render with value = full rng.
# Used by reset_tracking_filters() as a reliable alternative to updateNumericRangeInput,
# which does not reliably update renderUI-wrapped widgets.
reset_ranges_trigger  <- reactiveVal(0L)

# On new data upload, snap the TID-level filters to the new full range.
# trk_raw_stats_rv() is frozen at upload and never updated by rechunking, so
# this observer fires only on genuine new uploads.
observeEvent(trk_raw_stats_rv(), {
  rs <- trk_raw_stats_rv()
  if (is.null(rs)) return()
}, ignoreInit = TRUE)
