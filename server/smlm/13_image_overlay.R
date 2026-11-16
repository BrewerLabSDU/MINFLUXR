# Reads, optionally rotates, and prepares a data-URI + pixel size
confocal_meta <- reactive({
  file <- input$bg_image
  if (is.null(file) || !isTRUE(input$show_bg)) return(NULL)
  
  img <- tryCatch(magick::image_read(file$datapath), error = function(e) NULL)
  if (is.null(img)) return(NULL)
  
  rot <- suppressWarnings(as.integer(if (is.null(input$img_rotation)) 0 else input$img_rotation))
  if (is.finite(rot) && (rot %% 360) != 0) {
    img <- magick::image_rotate(img, rot)
  }
  # 2) Optional: mirror horizontally (left↔right)
  if (isTRUE(input$img_flip_x)) {
    img <- magick::image_flop(img)   # horizontal mirror
  }
  # 3) Optional: recolour through a single-colour LUT (Green, Magenta, ...)
  #    plus brightness / contrast adjustment.
  img <- apply_background_lut(img, input$img_lut %||% "original",
                              custom_color = input$img_lut_custom,
                              brightness   = input$img_brightness %||% 100,
                              contrast     = input$img_contrast   %||% 0)

  info <- magick::image_info(img)
  raw_png <- magick::image_write(img, format = "png")
  list(
    uri = paste0("data:image/png;base64,", base64enc::base64encode(raw_png)),
    width_px  = info$width,
    height_px = info$height
  )
})

image_overlay <- reactive({
  meta <- confocal_meta()
  if (is.null(meta)) return(NULL)
  
  ps <- suppressWarnings(as.numeric(input$px_size));  if (!is.finite(ps) || ps <= 0) return(NULL)
  x0 <- suppressWarnings(as.numeric(input$img_x0));   if (!is.finite(x0)) x0 <- 0
  y0 <- suppressWarnings(as.numeric(input$img_y0));   if (!is.finite(y0)) y0 <- 0
  
  sizex <- meta$width_px  * ps
  sizey <- meta$height_px * ps
  
  if (isTRUE(input$img_invert_y)) {
    y_anchor <- "top";    y_val <- y0 + sizey
  } else {
    y_anchor <- "bottom"; y_val <- y0
  }
  
  list(
    source  = meta$uri,
    xref    = "x",   yref = "y",
    x       = x0,    y    = y_val,
    sizex   = sizex, sizey = sizey,
    xanchor = "left", yanchor = y_anchor,
    sizing  = "stretch",
    opacity = if (is.null(input$img_opacity)) 0.6 else input$img_opacity,
    layer   = "below"
  )
})


# ---- Auto-align overlay helpers ----
.load_smlm_img <- function() {
  f <- input$bg_image
  if (is.null(f)) return(NULL)
  img <- tryCatch(magick::image_read(f$datapath), error = function(e) NULL)
  if (is.null(img)) return(NULL)
  rot <- suppressWarnings(as.integer(
    if (is.null(input$img_rotation)) 0 else input$img_rotation))
  if (is.finite(rot) && (rot %% 360) != 0) img <- magick::image_rotate(img, rot)
  if (isTRUE(input$img_flip_x)) img <- magick::image_flop(img)
  img
}

observeEvent(input$auto_foci, {
  img <- .load_smlm_img(); req(!is.null(img))
  df  <- isolate(plottable_data()); req(nrow(df) > 0)
  ps  <- suppressWarnings(as.numeric(input$px_size))
  req(is.finite(ps), ps > 0)
  res <- tryCatch(
    auto_align_overlay(img, df, ps, mode = "foci",
                       invert_y    = isTRUE(input$img_invert_y),
                       n_peaks        = as.integer(input$foci_n_peaks        %||% 5L),
                       bright_pct     = (input$foci_bright_pct     %||% 95)  / 100,
                       border_pct     = (input$foci_border_pct     %||% 10)  / 100,
                       eps_pct        = (input$foci_eps_pct        %||% 3)   / 100,
                       aspect_max     =  input$foci_aspect_max     %||% 4.0,
                       density_thr    = (input$foci_density_thr    %||% 20)  / 100,
                       density_cell   =  input$foci_density_cell   %||% 100,
                       density_smooth =  input$foci_density_smooth %||% 100),
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
  if (!is.null(res$img_peaks)) {
    smlm_diag(list(
      data_clusters = res$data_peaks,
      img_blobs     = data.frame(x = res$x0 + res$img_peaks$x,
                                 y = res$y0 + res$img_peaks$y)
    ))
  }
  updateCheckboxInput(session, "show_bg", value = TRUE)
  updateNumericInput(session, "img_x0", value = res$x0)
  updateNumericInput(session, "img_y0", value = res$y0)
  smlm_diag_mats(list(raw = res$img_mat_debug, smooth = res$img_mat_smooth))
  if (isTRUE(input$show_diag)) {
    pv <- res$img_peaks_debug
    mr <- res$mat_range_debug
    dd <- res$img_dim_debug
    dim_str <- if (!is.null(dd))
      paste0(" | src ", dd[1], "×", dd[2], " → proc ", dd[3], "×", dd[4],
             " px, raw bytes: ", dd[5], " (", dd[6], ")")
    else ""
    showNotification(
      paste0("Confocal detection - gray range: ", mr[1], "–", mr[2],
             " | peaks: ", paste(round(as.numeric(pv), 0), collapse = ", "),
             dim_str),
      type = "message", duration = 20)
  }
  showNotification("Foci alignment applied - nudge to fine-tune.",
                   type = "message", duration = 4)
}, ignoreInit = TRUE)

observeEvent(input$auto_boundary, {
  img <- .load_smlm_img(); req(!is.null(img))
  df  <- isolate(plottable_data()); req(nrow(df) > 0)
  ps  <- suppressWarnings(as.numeric(input$px_size))
  req(is.finite(ps), ps > 0)
  res <- tryCatch(
    auto_align_overlay(img, df, ps, mode = "boundary",
                       invert_y = isTRUE(input$img_invert_y)),
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
  updateCheckboxInput(session, "show_bg", value = TRUE)
  updateNumericInput(session, "img_x0", value = res$x0)
  updateNumericInput(session, "img_y0", value = res$y0)
  showNotification("Boundary alignment applied - nudge to fine-tune.",
                   type = "message", duration = 4)
}, ignoreInit = TRUE)

# ---- Nudge buttons for X and Y origin ----
local({
  deltas <- c(m1000 = -1000, m100 = -100, m10 = -10,
              p10   =    10, p100 =  100, p1000 = 1000)
  for (sfx in names(deltas)) {
    local({
      s <- sfx
      d <- deltas[[s]]
      observeEvent(input[[paste0("nx_", s)]], {
        v <- suppressWarnings(as.numeric(input$img_x0))
        updateNumericInput(session, "img_x0",
                           value = (if (is.finite(v)) v else 0) + d)
      }, ignoreInit = TRUE)
      observeEvent(input[[paste0("ny_", s)]], {
        v <- suppressWarnings(as.numeric(input$img_y0))
        updateNumericInput(session, "img_y0",
                           value = (if (is.finite(v)) v else 0) + d)
      }, ignoreInit = TRUE)
    })
  }
})

output$download_overlay_settings <- downloadHandler(
  filename = function() {
    base <- if (!is.null(input$active_file) && nzchar(input$active_file)) tools::file_path_sans_ext(input$active_file) else "smlm"
    paste0(base, "_image_alignment_", format(Sys.time(), "%Y%m%d-%H%M%S"), ".json")
  },
  content  = function(file) {
    ov <- list(
      show_confocal = input$show_confocal,
      show_bg       = input$show_bg,
      bg_image_name = if (!is.null(input$bg_image)) input$bg_image$name else NULL,
      px_size       = input$px_size,
      img_x0        = input$img_x0,
      img_y0        = input$img_y0,
      img_invert_y  = input$img_invert_y,
      img_flip_x    = input$img_flip_x,
      img_opacity   = input$img_opacity,
      img_rotation  = input$img_rotation,
      img_lut       = input$img_lut,
      img_lut_custom = input$img_lut_custom,
      img_brightness = input$img_brightness,
      img_contrast   = input$img_contrast
    )
    writeLines(jsonlite::toJSON(ov, pretty = TRUE, auto_unbox = TRUE, na = "null"),
               file, useBytes = TRUE)
  }
)
# --- Foci detection debug plot (shown when Show diagnostics is checked) ----
output$foci_debug_plot <- renderPlot({
  req(isTRUE(input$show_diag))
  mats <- smlm_diag_mats(); req(!is.null(mats))
  raw    <- mats$raw;    req(!is.null(raw))
  smooth <- mats$smooth; req(!is.null(smooth))
  par(mfrow = c(1, 2), mar = c(2, 2, 2, 1))
  image(t(raw)[, nrow(raw):1],    col = grDevices::gray.colors(256),
        main = "Raw gray (border masked)", axes = FALSE)
  image(t(smooth)[, nrow(smooth):1], col = grDevices::gray.colors(256),
        main = "After smoothing (detection input)", axes = FALSE)
}, bg = "#222")
