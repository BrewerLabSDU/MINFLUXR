uploaded_df <- reactive({
  df <- uploaded_df_rv()
  shiny::validate(shiny::need(!is.null(df) && nrow(df) > 0, "No data loaded yet."))
  df
})

# =========================================================================
# commit_z_scale() - APPLY THE AXIAL CALIBRATION TO DATA ALREADY IN MEMORY
# =========================================================================
# GOAL
#   Rescale every loaded dataset from the currently-applied axial factor to a
#   new one. Rescaling by the ratio is exact and avoids re-reading files;
#   z_scale_applied() tracks what is baked into the store so the ratio is
#   always right.
#
# WHY IT IS NOT WIRED TO input$z_scale DIRECTLY
#   Applying it replaces the active dataset, which fires the "new dataset"
#   reset below. Committing only on Plot Data (or on an upload that records
#   its own factor) keeps the change deliberate and lets the guard below tell
#   a rescale apart from a genuine dataset change.
#
# ARGUMENTS  new : the factor to apply;  reason : shown in the notification
# RETURNS    TRUE if anything was rescaled
# =========================================================================
z_rescale_guard <- reactiveVal(FALSE)

commit_z_scale <- function(new, reason = "Axial calibration") {
  new <- suppressWarnings(as.numeric(new))
  old <- z_scale_applied()
  if (!length(new) || !is.finite(new) || new <= 0) return(invisible(FALSE))
  if (!is.finite(old) || old <= 0) { z_scale_applied(new); return(invisible(FALSE)) }
  if (isTRUE(all.equal(new, old))) return(invisible(FALSE))

  ratio <- new / old
  store <- smlm_file_store()

  # Nothing axial loaded: record the factor so a later 3D upload and any export
  # label are right, but do no work and stay silent. Rescaling zeros and
  # announcing it is the noise this guard exists to prevent.
  if (!any(vapply(store, has_axial_data, logical(1)))) {
    z_scale_applied(new)
    return(invisible(FALSE))
  }

  if (length(store)) {
    # Tell the reset observer this is a rescale of the SAME data, not a new
    # dataset - otherwise every metadata filter the user set is wiped.
    z_rescale_guard(TRUE)
    # Multiplex z alignment offsets are in the same nm, so they rescale with
    # the data; otherwise the next nudge would undo a stale (unscaled) shift.
    .scale_z <- function(tbl) lapply(tbl, function(v) { v[["z"]] <- v[["z"]] * ratio; v })
    mplex_offsets_rv(.scale_z(isolate(mplex_offsets_rv())))
    store <- lapply(store, function(d) {
      if (is.data.frame(d) && "z" %in% names(d)) d$z <- d$z * ratio
      if (!is.null(attr(d, "mplex_off"))) attr(d, "mplex_off") <- .scale_z(attr(d, "mplex_off"))
      sp <- attr(d, "mplex_spec")                  # alignment recorded in a loaded CSV
      if (!is.null(sp)) attr(d, "mplex_spec") <- lapply(sp, function(c) {
        if (length(c$prior) == 3L) c$prior[["z"]] <- c$prior[["z"]] * ratio; c })
      d
    })
    smlm_file_store(store)
    if (nzchar(isolate(input$mplex_nudge_ch) %||% ""))
      .mplex_show_offset(isolate(input$mplex_nudge_ch))
    cur <- isolate(input$active_file)
    if (!is.null(cur) && nzchar(cur) && !is.null(store[[cur]])) uploaded_df_rv(store[[cur]])
    # The z spread filter is expressed in nm, and nm along z just changed
    # meaning - carrying the old bound over would filter on a stale scale.
    sd_filt$z <- c(0, Inf)
  }
  z_scale_applied(new)
  showNotification(sprintf("%s set to %.3f - z distances rescaled by %.4f. Axial spread filter cleared.",
                           reason, new, ratio), type = "message", duration = 6)
  invisible(TRUE)
}

# Commit on Plot Data, so typing a value does not rescale mid-keystroke.
observeEvent(input$plot_data, {
  commit_z_scale(input$z_scale)
}, ignoreInit = TRUE)

# Does the ACTIVE dataset actually carry axial data? Three states, because
# "nothing loaded yet" must not look like "this data is 2D": an empty app
# showing a greyed-out panel reads as broken rather than as not-applicable.
smlm_z_state <- reactive({
  df <- uploaded_df_rv()
  if (is.null(df) || !nrow(df)) return("nodata")
  if (!"z" %in% names(df) || !any(is.finite(df$z)) ||
      stats::sd(df$z, na.rm = TRUE) <= 1e-6) return("planar")
  "3d"
})

# Dim (not hide) the panel for planar data. Hiding it would raise "where did
# the setting go?"; dimmed-and-disabled shows it exists and says why it is off.
observe({
  ok <- !identical(smlm_z_state(), "planar")
  shinyjs::toggleState("z_scale", condition = ok)
  shinyjs::toggleClass(id = "z_scale_group", class = "dimmed", condition = !ok)
})

# Every axis-bearing control follows the SAME planar test as the calibration
# panel, so a dataset can never be simultaneously "2D" for one control and
# "3D" for another. 3D mode and the XZ/YZ side views on planar data do not
# error - they render a flat line - which is exactly the kind of quietly wrong
# plot worth making unreachable.
observe({
  planar <- identical(smlm_z_state(), "planar")

  shinyjs::toggleState(selector = "input[name='mode3d'][value='3D']",       condition = !planar)
  shinyjs::toggleState(selector = "input[name='view_plane_2d'][value='XZ']", condition = !planar)
  shinyjs::toggleState(selector = "input[name='view_plane_2d'][value='YZ']", condition = !planar)

  # Disabling a radio does not deselect it. Without this, loading planar data
  # while sitting in 3D/XZ would leave the view stuck in a mode whose control
  # is greyed out - unreachable and unfixable from the UI.
  if (planar) {
    if (!identical(isolate(input$mode3d) %||% "2D", "2D"))
      updateRadioButtons(session, "mode3d", selected = "2D")
    if (!identical(isolate(input$view_plane_2d) %||% "XY", "XY"))
      updateRadioButtons(session, "view_plane_2d", selected = "XY")
  }
})

output$z_axes_note <- renderUI({
  if (!identical(smlm_z_state(), "planar")) return(NULL)
  helpText(tags$b("2D dataset: "), "3D mode and the XZ / YZ side views are unavailable.")
})

# Shows whether the typed value is live yet.
output$z_scale_status <- renderUI({
  if (identical(smlm_z_state(), "planar"))
    return(helpText(tags$b("Not applicable: "),
                    "the active dataset has no axial (z) data."))

  pend <- suppressWarnings(as.numeric(input$z_scale))
  appl <- z_scale_applied()
  if (!length(pend) || !is.finite(pend) || pend <= 0)
    return(helpText(sprintf("Applied: %.3f. Enter a positive factor.", appl)))
  if (isTRUE(all.equal(pend, appl)))
    helpText(sprintf("Applied: %.3f.", appl))
  else
    helpText(tags$b(sprintf("Pending: %.3f", pend)),
             sprintf(" (applied: %.3f) - press Plot Data to apply.", appl))
})


# --------- 2) RESET FILTERS ON NEW DATASET ----------
observeEvent(uploaded_df(), {
  # An axial rescale replaces the active dataset object but it is the SAME
  # data - resetting every filter here would discard the user's whole setup
  # just because they changed the z factor. commit_z_scale() has already
  # cleared the one filter whose units actually changed (axial spread).
  if (isTRUE(z_rescale_guard())) { z_rescale_guard(FALSE); return() }

  # Respect the lock: when the user has checked "Lock filters", keep the
  # current filter values when switching between loaded files.
  # We still need to arm sd_skip_arming so the SD range observer doesn't fire
  # during the widget rebuild (sd_ranges_ui re-renders because sd_grey_tbl()
  # depends on uploaded_df()).  Without this, the brief NULL window while
  # widgets are being recreated causes the observer to write c(0,Inf) into
  # sd_filt, resetting the locked spread filter.
  if (isTRUE(isolate(input$lock_filters))) {
    sd_skip_arming(TRUE)
    later::later(function() sd_skip_arming(FALSE), delay = 0.8)
    return()
  }

  # Freeze meta_full_rng to current metadata widget values so the all.equal()
  # bypass in metadata_filtered_df passes during the browser round-trip window
  # between upload and the browser processing the new renderUI widgets.
  freeze_if_valid <- function(val) {
    if (!is.null(val) && length(val) == 2 && all(is.finite(as.numeric(val)))) val else NULL
  }
  cur_efo   <- freeze_if_valid(isolate(input$efo_range))
  cur_cfr   <- freeze_if_valid(isolate(input$cfr_range))
  cur_tim   <- freeze_if_valid(isolate(input$tim_range))
  cur_dcr   <- freeze_if_valid(isolate(input$dcr_range))
  cur_eco   <- freeze_if_valid(isolate(input$eco_range))
  cur_dwell <- freeze_if_valid(isolate(input$dwell_range))
  cur_fluo  <- freeze_if_valid(isolate(input$fluo_range))
  if (!is.null(cur_efo))   meta_full_rng$efo   <- cur_efo
  if (!is.null(cur_cfr))   meta_full_rng$cfr   <- cur_cfr
  if (!is.null(cur_tim))   meta_full_rng$tim   <- cur_tim
  if (!is.null(cur_dcr))   meta_full_rng$dcr   <- cur_dcr
  if (!is.null(cur_eco))   meta_full_rng$eco   <- cur_eco
  if (!is.null(cur_dwell)) meta_full_rng$dwell <- cur_dwell
  if (!is.null(cur_fluo))  meta_full_rng$fluo  <- cur_fluo

  # Reset SD filter synchronously - no JS round-trip needed.
  sd_filt$x <- c(0, Inf); sd_filt$y <- c(0, Inf); sd_filt$z <- c(0, Inf)
  # Signal tid_length_range_ui to reset to the new data's full range instead of
  # carrying forward the previous dataset's range.
  tid_reset_flag(TRUE)

  # Whenever a NEW dataset is loaded, reset filters to its full ranges
  shinyjs::reset("filters_form")

  vals <- filters_defaults_from_data(uploaded_df())

  session$onFlushed(function() {
    # Suppress SD range observer fires triggered by apply_filters updating the
    # SD widgets (only fires when the SD panel is already open).
    sd_skip_arming(TRUE)
    apply_filters(vals)
    meta_full_rng$efo   <- vals$efo_range
    meta_full_rng$cfr   <- vals$cfr_range
    meta_full_rng$tim   <- vals$tim_range
    meta_full_rng$dcr   <- vals$dcr_range
    meta_full_rng$eco   <- vals$eco_range
    meta_full_rng$dwell <- vals$dwell_range
    meta_full_rng$fluo  <- vals$fluo_range
    later::later(function() sd_skip_arming(FALSE), delay = 0.8)
  }, once = TRUE)

  # Also clear any loaded filter-settings JSON file input (optional)
  shinyjs::runjs("
    const el = document.getElementById('load_settings');
    if (el) { el.value=''; $(el).trigger('change'); }
  ")

}, ignoreInit = TRUE)
