smlm_axis_x_span  <- reactiveVal(NULL)
smlm_zoom_ranges  <- reactiveVal(NULL)   # list(x=c(x0,x1), y=c(y0,y1)) from 2D plot zoom
smlm_2d_crop_rv   <- reactiveVal(NULL)   # active 2D crop box: list(x=c(x0,x1), y=c(y0,y1))
# Committed 3D crop region. The box/sphere inputs are only a PREVIEW; the crop
# is applied when the user presses "Apply crop", which commits the region here.
# list(mode="box", x=c(lo,hi), y=..., z=...) or list(mode="sphere", c=c(cx,cy,cz), r=r).
trim_applied_rv   <- reactiveVal(NULL)

smlm_default_x_span <- function() {
  df <- isolate(plottable_data())
  if (!is.data.frame(df) || !nrow(df)) return(NULL)

  lim <- get_equal_limits_2d(df$x, df$y)
  ovl <- tryCatch(isolate(image_overlay()), error = function(e) NULL)
  lim <- expand_limits_with_image(lim, ovl, invert_y = isTRUE(isolate(input$img_invert_y)))

  span <- suppressWarnings(as.numeric(diff(lim$x)))
  if (length(span) == 1 && is.finite(span) && span > 0) span else NULL
}

smlm_push_scale_bar <- function() {
  proxy <- plotlyProxy("scatter2D", session)
  if (!isTRUE(input$show_scale_bar) || !identical(input$mode3d, "2D")) {
    plotlyProxyInvoke(proxy, "relayout", list(shapes = list(), annotations = list()))
    return()
  }

  span <- smlm_axis_x_span()
  if (length(span) != 1 || !is.finite(span) || span <= 0) {
    span <- smlm_default_x_span()
    if (length(span) == 1 && is.finite(span) && span > 0) smlm_axis_x_span(span)
  }

  sb_len   <- suppressWarnings(as.numeric(input$scale_bar_length))
  sb_thick <- suppressWarnings(as.numeric(input$scale_bar_thickness))
  sb_font  <- suppressWarnings(as.numeric(input$scale_font_size))
  sb <- make_scale_bar(sb_len, span, sb_thick, dark = isTRUE(input$dark_mode),
                       scale_font_size = if (is.finite(sb_font) && sb_font > 0) sb_font else 18)
  if (!is.null(sb)) {
    plotlyProxyInvoke(proxy, "relayout",
                      list(shapes = sb$shapes, annotations = sb$annotations))
  }
}

# ── 3D auto-rotate ──────────────────────────────────────────────────────────
observeEvent(
  list(input$auto_rotate_3d, input$rotate_speed_3d, input$rotate_elev_3d, input$rotate_dist_3d),
  ignoreNULL = FALSE, {
    if (isTRUE(input$auto_rotate_3d)) {
      speed <- input$rotate_speed_3d %||% 2
      elev  <- input$rotate_elev_3d  %||% 25
      dist  <- input$rotate_dist_3d  %||% 2.5
      shinyjs::runjs(sprintf("startRotate3D(%s, %s, %s);", speed, elev, dist))
    } else {
      shinyjs::runjs("stopRotate3D();")
    }
  }
)

# ── 3D recording ────────────────────────────────────────────────────────────
observeEvent(input$record_start_3d, {
  w       <- input$record_width_3d   %||% 800
  h       <- input$record_height_3d  %||% 600
  fps     <- input$record_fps_3d     %||% 20
  bitrate <- input$record_bitrate_3d %||% 10
  shinyjs::runjs(sprintf("startRecord3D(%d, %d, %d, %s);", as.integer(w), as.integer(h), as.integer(fps), bitrate))
})

observeEvent(input$record_stop_3d, {
  shinyjs::runjs("stopRecord3D();")
})

#1 ) UPLOAD AND PREPROCESS DATA ----------
uploaded_df_rv  <- reactiveVal(NULL)
smlm_file_store <- reactiveVal(list())   # fname → data frame

# ── Axial (z) calibration ────────────────────────────────────────────────
# MINFLUX z is multiplied by a scaling factor at import to correct the
# refractive-index mismatch between immersion medium and sample. pyMINFLUX
# ships 0.7 as a default; the physically correct value depends on the actual
# index pair and imaging depth, and 0.67 is also in use - so it is exposed
# rather than hard-coded.
#
# It multiplies EVERY axial distance downstream: z trace spread, cluster
# extent in z, 3D nearest-neighbour distances, Rg. Treat it as a calibration
# to be measured on a known standard, not a display preference.
SMLM_Z_SCALE_DEFAULT <- 0.7
# The factor currently baked into smlm_file_store() - the single source of
# truth. input$z_scale is only a PENDING value until Plot Data commits it, so
# nothing may read the input directly to decide what the data is scaled by.
z_scale_applied <- reactiveVal(SMLM_Z_SCALE_DEFAULT)

# Helper: load one file from a temp path, return df or NULL
.load_one_file <- function(path, fname) {
  ext <- get_ext(fname)
  tryCatch({
    withProgress(message = paste0("Loading ", fname, "…"), value = 0, {
      incProgress(0.03, detail = "Preparing")
      loaded <- switch(
        ext,
        "csv" = {
          incProgress(0.20, detail = "Reading CSV...")
          read.csv(path, check.names = FALSE)
        },

        "json" = {
          incProgress(0.10, detail = "Reading JSON...")
          read_minflux_json_as_df(
            path,
            keep_final_only = TRUE,
            keep_valid_only = TRUE,
            final_itr = NA_integer_,
            cfr_strategy = "last_nonzero",
            # Read RAW: .npy/.json come straight off the instrument with no
            # axial correction applied. The app-wide factor is applied once,
            # below, for every source alike.
            z_scale = 1.0,
            dwell_time_ms = 1.0,
            pool_dcr = isTRUE(input$pool_dcr_json)
          )
        },
        "npy" = {
          incProgress(0.10, detail = "Reading .npy (binary)")
          read_minflux_npy_as_df(
            path,
            keep_final_only = TRUE,
            keep_valid_only = TRUE,
            final_itr = NA_integer_,
            cfr_strategy = "last_nonzero",
            # Read RAW: .npy/.json come straight off the instrument with no
            # axial correction applied. The app-wide factor is applied once,
            # below, for every source alike.
            z_scale = 1.0,
            dwell_time_ms = 1.0,
            pool_dcr = isTRUE(input$pool_dcr_json)
          )
        },
        stop("Unsupported file type: ", ext)
      )
      if (!is.data.frame(loaded) || !nrow(loaded)) stop("No rows loaded from file.")
      incProgress(0.15, detail = "Canonicalizing columns")
      loaded <- canonicalize_minflux_df(loaded)
      if (!"z" %in% names(loaded)) loaded$z <- 0

      # ── Normalize the axial calibration across all three sources ──────────
      # Everything reaching this point is reduced to RAW (uncorrected) z, then
      # the single app-wide factor is applied. Without this, .npy/.json would
      # follow the calibration control while CSVs stayed frozen at whatever
      # their producer baked in - and a CSV would get 0.7 applied twice.
      #
      #   .npy / .json  raw already (read above with z_scale = 1)
      #   .csv          divide out the factor its producer applied
      if (identical(ext, "csv")) {
        loaded  <- strip_csv_z_calibration(loaded)     # -> raw
        src_cal <- attr(loaded, "z_cal")
      } else {
        src_cal <- list(factor = 1.0, source = "raw instrument data", certain = TRUE)
      }
      # Applied at the COMMITTED factor, not the pending input value, so the
      # store is never left at a factor the user has not confirmed.
      loaded$z <- loaded$z * z_scale_applied()         # -> app calibration
      attr(loaded, "z_cal") <- src_cal

      # Say nothing about the axial calibration for planar data - a 2D upload
      # should not produce messages about z. has_axial_data() tests VARIATION,
      # not presence: 2D exports still carry a z column, filled with zeros.
      if (identical(ext, "csv") && isTRUE(src_cal$has_z)) {
        showNotification(
          sprintf("%s: axial calibration %.3f detected (%s)%s.",
                  fname, src_cal$factor, src_cal$source,
                  if (isTRUE(src_cal$certain)) "" else " - assumed"),
          type = "message", duration = 8)
      }
      incProgress(0.07, detail = "Finalizing")
      loaded
    })
  }, error = function(e) {
    showNotification(paste0("Could not load ", fname, ": ", conditionMessage(e)),
                     type = "error", duration = 15)
    NULL
  })
}

observeEvent(input$file, {
  req(input$file)
  store     <- smlm_file_store()
  new_names <- character(0)

  adopt_factor <- NULL   # a recorded (not guessed) factor from an uploaded file

  for (i in seq_len(nrow(input$file))) {
    fname <- input$file$name[i]
    if (fname %in% names(store)) next   # already loaded - skip duplicate

    df <- .load_one_file(input$file$datapath[i], fname)
    if (is.null(df)) next

    # A file of ours records the factor it was written at. Adopt it so the
    # control reflects the data being loaded and z comes back exactly as
    # exported. Only `certain` factors are adopted - the 0.7 assumed for
    # pyMINFLUX and legacy exports is a guess and must not silently become
    # the session's calibration.
    # ...and only when the file actually has axial data: adopting a factor for
    # a planar dataset would change a number that affects nothing, while
    # announcing it.
    zc <- attr(df, "z_cal")
    if (!is.null(zc) && isTRUE(zc$certain) && isTRUE(zc$has_z) &&
        is.finite(zc$factor) && zc$factor > 0)
      adopt_factor <- zc$factor

    df$.dataset_tag <- paste0(sample(c(0:9, letters[1:6]), 4, replace = TRUE), collapse = "")
    # A CSV exported from a multiplex keeps its channels (see .mplex_adopt_csv).
    if (identical(get_ext(fname), "csv")) {
      df <- .mplex_adopt_csv(df, fname)
      if (!is.null(attr(df, "mplex_spec")))
        showNotification(sprintf("%s: multiplex with %d channels recognised.",
                                 fname, length(attr(df, "mplex_spec"))),
                         type = "message", duration = 6)
    }
    store[[fname]] <- df
    new_names <- c(new_names, fname)

    if (identical(get_ext(fname), "json")) {
      final_used <- if ("itr" %in% names(df)) max(df$itr, na.rm = TRUE) else NA_integer_
      showNotification(
        sprintf("%s - JSON import: using final localization (itr=%s); CFR carried over per TID.",
                fname, final_used),
        type = "message", duration = 8)
    }
  }

  if (!length(new_names)) return()
  smlm_file_store(store)

  # Directly update uploaded_df_rv here so plotting works even when
  # input$active_file doesn't change (e.g. re-uploading the same filename
  # after removal, where the client-side value may still equal new_names[1]
  # and the observeEvent(input$active_file) observer never fires).
  selected_df <- store[[new_names[1]]]
  if (!is.null(selected_df)) {
    uploaded_df_rv(selected_df)
    .mplex_sync_active(selected_df)
    hist_base <- c("efo", "tid", "cfr", "tim", "eco", "dcr", "dwell", "fluo")
    hist_choices <- if ("dcr_pooled" %in% names(selected_df))
      c(hist_base, "dcr_pooled") else hist_base
    updateSelectInput(session, "hist_variable", choices = hist_choices,
                      selected = isolate(input$hist_variable) %||% "efo")
  }

  updateSelectInput(session, "active_file",
                    choices  = names(store),
                    selected = new_names[1])

  # Adopt a recorded axial factor, once, at upload. commit_z_scale() rescales
  # everything already in the store to match, so the newly loaded file ends up
  # at exactly the z it was exported with.
  if (!is.null(adopt_factor)) {
    updateNumericInput(session, "z_scale", value = adopt_factor)
    commit_z_scale(adopt_factor, reason = "Axial calibration from file")
  }
}, ignoreInit = TRUE)

# Switch active dataset when the file selector changes
observeEvent(input$active_file, {
  req(input$active_file, nzchar(input$active_file))
  df <- smlm_file_store()[[input$active_file]]
  if (is.null(df)) return()

  # Snapshot sd_filt BEFORE the data change so we can restore it after the
  # reactive graph settles.  This is the authoritative lock mechanism for SD
  # ranges - more reliable than the timing-based sd_skip_arming approach.
  snap_sd <- if (isTRUE(isolate(input$lock_filters)))
    isolate(list(x = sd_filt$x, y = sd_filt$y, z = sd_filt$z))
  else NULL

  uploaded_df_rv(df)
  .mplex_sync_active(df)

  # Update histogram variable choices - show dcr_pooled only when present
  hist_base <- c("efo", "tid", "cfr", "tim", "eco", "dcr", "dwell", "fluo")
  hist_choices <- if ("dcr_pooled" %in% names(df))
    c(hist_base, "dcr_pooled") else hist_base
  updateSelectInput(session, "hist_variable", choices = hist_choices,
                    selected = isolate(input$hist_variable) %||% "efo")

  if (!is.null(snap_sd)) {
    session$onFlushed(function() {
      sd_skip_arming(TRUE)
      sd_filt$x <- snap_sd$x
      sd_filt$y <- snap_sd$y
      sd_filt$z <- snap_sd$z
      # Update the widgets to display the restored values
      shinyWidgets::updateNumericRangeInput(session, "sd_x_range",
        value = c(max(0, snap_sd$x[1]), if (is.infinite(snap_sd$x[2])) 1e6 else snap_sd$x[2]))
      shinyWidgets::updateNumericRangeInput(session, "sd_y_range",
        value = c(max(0, snap_sd$y[1]), if (is.infinite(snap_sd$y[2])) 1e6 else snap_sd$y[2]))
      shinyWidgets::updateNumericRangeInput(session, "sd_z_range",
        value = c(max(0, snap_sd$z[1]), if (is.infinite(snap_sd$z[2])) 1e6 else snap_sd$z[2]))
      later::later(function() sd_skip_arming(FALSE), delay = 0.5)
    }, once = TRUE)
  }
}, ignoreInit = TRUE)

observeEvent(input$flush_file, {
  fname <- isolate(input$active_file)
  if (!nzchar(fname %||% "")) return()
  store <- smlm_file_store()
  store[[fname]] <- NULL
  smlm_file_store(store)
  remaining <- names(store)
  if (length(remaining) == 0) {
    uploaded_df_rv(NULL)
    updateSelectInput(session, "active_file", choices = character(0))
  } else {
    updateSelectInput(session, "active_file",
                      choices  = remaining,
                      selected = remaining[1])
  }
})
