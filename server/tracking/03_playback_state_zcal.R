# Non-reactive loop state for sequential mode (safe for later())
is_seq   <- reactive(identical(input$trk_play_mode, "seq"))
is_simul <- reactive(identical(input$trk_play_mode, "simul"))

# Bin width (ms) is the master time-resolution control. The frame count is
# derived from it (window / bin width) rather than entered separately - the
# old "Frames per playthrough" input was always back-computed into bin_dt.
# Used by sequential playback (frames per TID) and its step slider; the
# simultaneous path builds its own bin edges directly from trk_bin_dt.
n_frames_eff <- reactive({
  tot <- suppressWarnings(as.numeric(input$trk_play_total))
  bw  <- suppressWarnings(as.numeric(input$trk_bin_dt))
  if (!is.finite(tot) || tot <= 0 || !is.finite(bw) || bw <= 0) return(120L)
  max(2L, as.integer(round(tot / bw)))
})
loop <- new.env(parent = emptyenv())
loop$playing     <- FALSE
loop$sp          <- NULL   # list of data.frames per tid (x, y, t)
loop$tids        <- NULL   # character vector of tids in play order
loop$cur_tid     <- 1L
loop$cur_frame     <- 1L
loop$ggen      <- 0L      # generation for global (simultaneous) loop
loop$gplaying  <- FALSE   # playing flag for global loop
loop$trail_len   <- 0.5
loop$trail_nmax  <- 50L
loop$pt_size     <- 6
loop$col_map     <- NULL
loop$markers_idx <- 3L
loop$trail_start <- 2L
loop$ms_per_unit <- 0
loop$speed       <- 1    # "Playback speed (×)" multiplier, captured at play start
loop$gen         <- 0L   # sequential playback generation
loop$auto_tuning      <- FALSE  # TRUE while apply_auto_tune() is batching input updates
loop$last_fire_layout <- TRUE   # TRUE = last trigger was fire_layout (full rebuild)


# ---- Optional chunk picker (manifest.csv in same folder as upload) ----
output$trk_chunk_ui <- renderUI({
  req(trk_has_data())
  
  # If we're still preparing chunks, show a placeholder message
  if (isTRUE(rv$chunking)) {
    return(tags$div(
      "Preparing chunks... This can take a bit for very large files.",
      style = "color:#555;"
    ))
  }
  
  # Otherwise, if a manifest exists, show the chunk selector
  req(!is.null(rv$manifest_dir))
  dirp <- rv$manifest_dir
  manf <- file.path(dirp, "manifest.csv")
  if (file.exists(manf)) {
    man <- readr::read_csv(manf, show_col_types = FALSE)
    chunk_choices <- setNames(man$file, paste0("Chunk ", man$chunk, " (", man$n_tids, " tids)"))
    all_choices <- c("All TIDs" = "__all__", chunk_choices)
    # Chunk 1 if the user was scoped to a chunk before this rebuild, else all TIDs.
    sel_now <- if (identical(trk_chunk_pref(), "first") && nrow(man) >= 1L)
                 man$file[1] else "__all__"
    return(tagList(
      selectInput("trk_chunk_sel", "Playback chunk", choices = all_choices, selected = sel_now),
      helpText(style = "font-size:0.8em; color:#666; margin-top:-6px;",
        "Filtering, gap-splitting, analysis and the full-track overlay always use ALL ",
        "TIDs. This selector only limits which TIDs are ANIMATED during playback, to ",
        "keep the animation smooth - \"All TIDs\" animates everything.")
    ))
  }

  NULL
})
output$trk_msd_tid_ui <- renderUI({
  req(isTRUE(input$trk_show_msd), isTRUE(input$trk_msd_show_single))
  df <- trk(); req(nrow(df) > 0)
  min_pts <- max(5L, as.integer(input$trk_msd_min_pts %||% 20L))
  # Only offer TIDs that have enough points to compute MSD
  counts  <- table(df$tid)
  tids    <- sort(names(counts[counts >= min_pts]))
  req(length(tids) > 0)
  selectizeInput(
    "trk_msd_tid", "Choose TID for overlay:",
    choices = tids, selected = tids[1],
    options = list(maxOptions = 5000)
  )
})
# 2D/3D toggle and the XY/XZ/YZ plane selector are static inputs (ui_tracking.R).
# They used to be renderUIs that VANISHED on planar data. Two reasons that is
# worse than greying them:
#   - a control that disappears reads as a bug, not as "not applicable";
#   - the plane selector keyed off trk(), which only exists after Plot Data, so
#     it was absent on upload and appeared later.
# The observer below disables the z-dependent choices instead, using the SAME
# planar test as the axial calibration panel.
observe({
  planar <- identical(trk_z_state(), "planar")

  shinyjs::toggleState(selector = "input[name='trk_view_mode'][value='3d']",  condition = !planar)
  shinyjs::toggleState(selector = "input[name='trk_view_plane'][value='XZ']", condition = !planar)
  shinyjs::toggleState(selector = "input[name='trk_view_plane'][value='YZ']", condition = !planar)

  # Disabling a radio does not deselect it - force the view back, or a planar
  # upload while in 3D/XZ leaves the plot in a mode whose control is greyed out.
  if (planar) {
    if (!identical(isolate(input$trk_view_mode) %||% "2d", "2d"))
      updateRadioButtons(session, "trk_view_mode", selected = "2d")
    if (!identical(isolate(input$trk_view_plane) %||% "XY", "XY"))
      updateRadioButtons(session, "trk_view_plane", selected = "XY")
  }
})

output$trk_z_axes_note <- renderUI({
  if (!identical(trk_z_state(), "planar")) return(NULL)
  helpText(tags$b("2D dataset: "), "3D mode and the XZ / YZ side views are unavailable.")
})

# Switching view mode requires a full plot rebuild
observeEvent(input$trk_view_mode, {
  # Compare against what was last RENDERED, so a programmatic reset to "2d"
  # (planar guard above) or a value that already matches the plot does not
  # trigger a redundant rebuild.
  if (identical(input$trk_view_mode %||% "2d", loop$rendered_view %||% "2d")) return()
  fire_layout("view_mode_2d3d")
}, ignoreInit = TRUE)
# A plane change is a coordinate remap baked in frames_store, so it needs a rebuild -
# same creation-fire guard as the 2D/3D toggle (compare to the plane last rendered).
observeEvent(input$trk_view_plane, {
  if (identical(input$trk_view_plane %||% "XY", loop$rendered_plane %||% "XY")) return()
  fire_layout("view_plane")
}, ignoreInit = TRUE)
# ── Axial calibration: pending vs applied ────────────────────────────────
# input$trk_z_scale is only a PENDING value; trk_z_applied() is what the data
# is actually scaled by. Committing on Plot Data matches how every other
# tracking setting behaves (see the note above observeEvent(input$trk_plot_data))
# and avoids re-running the whole pipeline on every keystroke - the reason the
# old control needed a debounce.
trk_z_applied <- reactiveVal(TRK_Z_SCALE_DEFAULT)

observeEvent(input$trk_plot_data, {
  val <- suppressWarnings(as.numeric(input$trk_z_scale))
  if (!length(val) || !is.finite(val) || val <= 0) return()
  if (isTRUE(all.equal(val, trk_z_applied()))) return()
  trk_z_applied(val)
  # Record it either way, but only announce it when there is axial data for it
  # to act on - a 2D dataset should never surface a z message.
  if (identical(trk_z_state(), "3d"))
    showNotification(sprintf("Axial calibration applied: %.3f.", val),
                     type = "message", duration = 5)
  # No fire_layout here on purpose: the generic Plot Data handler below already
  # calls it, and priority = 10 guarantees this runs first - so that single
  # rebuild picks up the new factor. Firing here too would rebuild twice.
}, ignoreInit = TRUE, priority = 10)

# Three states, because "nothing uploaded yet" must not look like "this data is
# 2D". Reads the RAW loaded data, not trk(), so the panel settles as soon as a
# file is uploaded rather than waiting for the first Plot Data.
trk_z_state <- reactive({
  df <- tryCatch(trk_raw_loaded(), error = function(e) NULL)
  if (is.null(df) || !nrow(df)) return("nodata")
  if (isTRUE(.has_real_z(df))) "3d" else "planar"
})

# Dim (not hide) for planar data - see the SMLM equivalent for the reasoning.
observe({
  ok <- !identical(trk_z_state(), "planar")
  shinyjs::toggleState("trk_z_scale", condition = ok)
  shinyjs::toggleClass(id = "trk_z_scale_group", class = "dimmed", condition = !ok)
})

output$trk_z_scale_status <- renderUI({
  if (identical(trk_z_state(), "planar"))
    return(helpText(tags$b("Not applicable: "),
                    "the uploaded dataset has no axial (z) data."))

  pend <- suppressWarnings(as.numeric(input$trk_z_scale))
  appl <- trk_z_applied()
  if (!length(pend) || !is.finite(pend) || pend <= 0)
    return(helpText(sprintf("Applied: %.3f. Enter a positive factor.", appl)))
  if (isTRUE(all.equal(pend, appl)))
    helpText(sprintf("Applied: %.3f.", appl))
  else
    helpText(tags$b(sprintf("Pending: %.3f", pend)),
             sprintf(" (applied: %.3f) - press Plot Data to apply.", appl))
})
