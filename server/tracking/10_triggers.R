## ─────────────────────────────────────────────────────────────────────────────
## server/tracking/10_triggers.R  -  render triggers
## ─────────────────────────────────────────────────────────────────────────────
## Part of tracking_server(). Not a standalone script: its statements are evaluated in
## the server function's environment, in file-name order (server/parts.R).
##
## layout_trigger, plot_trigger and data_ver, fire_layout / fire_frames, and
## the Plot Data button.
##
## INTERFACE (generated: Rscript tools/server_map.R --headers; do not edit by hand)
##   DEFINES for other parts:
##     cluster_color_map_stored, conf_color_map_stored, data_ver, fire_frames,
##     fire_layout, has_z, layout_trigger, plot_trigger, trk_axis_x_span,
##     trk_data_key, trk_default_axis_state, trk_locked_ranges, trk_plot_ready,
##     use_3d
##   READS from other parts:
##     01_state_sigma: is_playing, rechunk_request, rv
##     03_playback_state_zcal: loop
##     06_filter_pipeline: bump_clean
##     11_data_frames: frames_store, trk
## END INTERFACE
## ─────────────────────────────────────────────────────────────────────────────

# ---- Two-level trigger system -----------------------------------------------
# layout_trigger: fires only when the Plotly CHART STRUCTURE changes (number of
#   traces, colors, static lines, heatmap layer, play mode).  renderPlotly
#   depends on this  -  reconstructing the full WebGL chart is expensive and
#   blanks the display, so it must happen as rarely as possible.
#
# plot_trigger: fires whenever FRAME DATA needs rebuilding (bin_dt, n_frames,
#   trail params, marker quotas, auto-tune).  frames_store/trk/seq_data depend
#   on this.  renderPlotly does NOT depend on plot_trigger; after a frames
#   rebuild the new frame 1 is pushed to the chart via draw_frame() proxy.
layout_trigger    <- reactiveVal(NA)
plot_trigger      <- reactiveVal(NA)
# TRUE once the track plot has actually rendered WITH event_register(), so the
# relayout observer never reads event_data() before the "track" source is
# registered (that prints a benign, deferred "source not registered" warning
# suppressWarnings can't catch). Stays FALSE until a plot with data renders -
# layout_trigger alone isn't enough (it can be set while frames_store is empty).
trk_plot_ready    <- reactiveVal(FALSE)
# Always NULL since the Rg/Rmax panel was retired (2026-08-03) - the observer that
# filled it moved to retired_tracking_server.R, deleted 2026-11 (in the git history). Kept so the playback colour override
# below needs no change, and so reinstating the retired panel just works.
conf_color_map_stored    <- reactiveVal(NULL)  # pre-computed tid->hex when conf coloring is on
cluster_color_map_stored <- reactiveVal(NULL)  # pre-computed tid->hex when MSD cluster coloring is on

# TRUE when the loaded data has a usable z column (enables 3D toggle). Requires real
# axial spread, so a placeholder all-zero z column does not offer a 3D mode that
# would silently pool a zero-variance axis (see .has_real_z).
has_z <- reactive({
  .has_real_z(trk())
})

# TRUE when the user has selected 3D mode and data supports it
use_3d <- reactive({
  isTRUE(has_z()) && identical(input$trk_view_mode %||% "2d", "3d")
})
trk_axis_x_span   <- reactiveVal(NULL)
trk_locked_ranges <- reactiveVal(NULL)   # list(x = c(x0,x1), y = c(y0,y1)) - held fixed during playback
# Stable key that only changes on new upload/rechunk - used as uirevision so
# display-only re-renders (show/hide traces, dark mode) do NOT reset user zoom.
trk_data_key      <- reactiveVal("init")
trk_default_axis_state <- function() {
  pre <- isolate(frames_store())
  if (!is.list(pre) || !length(pre$frames)) return(NULL)

  f1 <- pre$frames[[1]]
  bx <- f1$markers_x[is.finite(f1$markers_x)]
  by <- f1$markers_y[is.finite(f1$markers_y)]
  if (!length(bx) || !length(by)) {
    sp <- pre$sp
    bx <- unlist(lapply(sp, `[[`, "x"), use.names = FALSE)
    by <- unlist(lapply(sp, `[[`, "y"), use.names = FALSE)
    bx <- bx[is.finite(bx)]
    by <- by[is.finite(by)]
  }
  if (!length(bx) || !length(by)) return(NULL)

  xr   <- range(bx)
  yr   <- range(by)
  xpad <- max(1, diff(xr) * 0.04)
  ypad <- max(1, diff(yr) * 0.04)
  list(
    xr = xr,
    yr = yr,
    xpad = xpad,
    ypad = ypad,
    x_span = diff(xr) + 2 * xpad
  )
}

# data_ver: bumped whenever the filtered dataset itself changes (filter inputs OR new
# upload/layout).  trk() is eventReactive on this so prep_tracking + tid_color_map
# are skipped on rebuild_tick-only (playback param) plot_trigger fires.
data_ver <- reactiveVal(0L)

# src is a human-readable label of the trigger (kept for self-documenting call sites).
fire_layout <- function(src = "?") {
  t <- Sys.time()
  loop$last_fire_layout <- TRUE
  layout_trigger(t)
  plot_trigger(t)
  data_ver(isolate(data_ver()) + 1L)   # layout = new data shape; rerun trk()
}
fire_frames <- function(src = "?") {
  loop$last_fire_layout <- FALSE
  plot_trigger(Sys.time())
}

# ── Plot Data - the single apply action ─────────────────────────────────────
# Applies everything the user has changed: cleaning (split / end-trim / continuity
# gate), the metadata / length / duration / spatial filters, re-chunks the
# survivors at the current TIDs-per-chunk setting, and re-renders.
#
# There is no separate rechunk button. The rechunk path already ended in
# trk_data_key(), which fires the layout, so the two actions were doing the same
# work with the rechunk merely persisting the result - splitting them into two
# presses only raised the question of which order to use them in.
#
# NO fire_layout() on the rechunk branch: trk_data_key() fires it, and calling
# both would render the chart twice. The fallback branch covers a dataset with no
# manifest (chunking is skipped for some file types), where there is nothing to
# rechunk and a plain rebuild is the whole job.
observeEvent(input$trk_plot_data, {
  bump_clean()      # (re)apply split / end-trim / continuity gate with current settings
  manf <- if (!is.null(rv$original_manifest_dir))
            file.path(rv$original_manifest_dir, "manifest.csv") else NA_character_
  if (!is.na(manf) && file.exists(manf)) rechunk_request(isolate(rechunk_request()) + 1L)
  else fire_layout("plot_data_btn")
})
observeEvent(list(
  input$trk_play_mode,
  input$trk_view_mode_main, input$trk_full_connect, input$trk_full_raw,
  input$trk_show_full_tracks_faint, input$trk_show_moving_trails, input$trk_color_mode,
  input$trk_dark_show, input$trk_color_speed_win, input$trk_chunk_scope_full
  # trk_show_heatmap deliberately NOT here - the dwell-time heatmap is heavy, so it
  # applies on the next Plot Data (heatmap_store is plot_trigger-gated). Toggling the
  # checkbox alone no longer forces a full rebuild that showed nothing until Plot Data.
), {
  if (!is.na(layout_trigger()) && !isTRUE(is_playing())) fire_layout("display_opts")
}, ignoreInit = TRUE)

# Frame-data-only changes -> frames rebuild, chart stays
rebuild_tick <- reactiveVal(0L)
observeEvent(list(
  input$trk_bin_dt,
  input$trk_target_markers, input$trk_per_tid_cap,
  input$trk_trail_len, input$trk_trail_nmax,
  input$trk_max_pts_tid, input$trk_play_total
), {
  if (!isTRUE(is_playing()) && !is.na(plot_trigger()) && !isTRUE(loop$auto_tuning)) {
    rebuild_tick(isolate(rebuild_tick()) + 1L)
  }
}, ignoreInit = TRUE)

observeEvent(debounce(rebuild_tick(), 200), {
  fire_frames()
}, ignoreInit = TRUE)
