## ─────────────────────────────────────────────────────────────────────────────
## server/tracking/24_calibration.R  -  panel: calibration from an immobile sample
## ─────────────────────────────────────────────────────────────────────────────
## Part of tracking_server(). Not a standalone script: its statements are evaluated in
## the server function's environment, in file-name order (server/parts.R).
##
## Measuring the prefactor a and the wander floor per file (R/est_calib.R),
## and calibration profiles.
##
## INTERFACE (generated: Rscript tools/server_map.R --headers; do not edit by hand)
##   DEFINES for other parts:
##     trk_cal_profiles
##   READS from other parts:
##     04_datasets_upload: trk_store
##     06_filter_pipeline: trk_filter_params
##     07_cycle_time_sequence: trk_seq_params
##     15_msd: .msd_settings
##     25_compare_cells: .cmp_read_file
## END INTERFACE
## ─────────────────────────────────────────────────────────────────────────────

# ============================================================================
# INSTRUMENT CALIBRATION - a and D_w from an immobile sample, kept as profiles
# ============================================================================
# GOAL: measure the photon-budget prefactor a and the wander floor D_w for one
# acquisition sequence (calib_cell, R/est_calib.R), keep the result as a named
# profile that records the sequence it belongs to, and let profiles be saved,
# loaded and applied, globally or per file in Compare Cells. A profile is the
# unit a lab measures once per sequence and reuses; without one, a falls back to
# the ideal-estimator value of the pattern and the readout says "uncalibrated".
trk_cal_profiles <- reactiveVal(list())   # name -> profile (list)
trk_cal_result   <- reactiveVal(NULL)     # per-cell table of the last measurement

# Loaded files to choose the immobile ones from.
observe({
  fs  <- names(trk_store())
  sel <- isolate(input$trk_cal_files)
  updateSelectInput(session, "trk_cal_files", choices = fs, selected = intersect(sel, fs))
})

# A loaded sequence file fills in the descriptors the profile records.
observeEvent(input$trk_seq_json, {
  sp <- tryCatch(trk_seq_params(), error = function(e) NULL)
  if (is.null(sp) || !isTRUE(sp$ok)) return(invisible(NULL))
  if (is.finite(sp$patDwellTime_ms %||% NA)) updateNumericInput(session, "trk_cal_dwell", value = round(sp$patDwellTime_ms * 1000))
  if (is.finite(sp$phtLimit %||% NA))        updateNumericInput(session, "trk_cal_pl", value = sp$phtLimit)
  if (is.finite(sp$damping %||% NA))         updateNumericInput(session, "trk_cal_damping", value = sp$damping)
}, ignoreInit = TRUE, priority = -1)

observeEvent(input$trk_cal_run, {
  files <- input$trk_cal_files
  if (!length(files)) {
    showNotification("Choose the immobile-sample files first.", type = "warning", duration = 5)
    return(invisible(NULL))
  }
  L   <- suppressWarnings(as.numeric(input$trk_L_tcp %||% DMAX_L_TCP_DEFAULT))
  if (!is.finite(L) || L <= 0) L <- DMAX_L_TCP_DEFAULT
  pat  <- .tcp_pattern_name(input$trk_tcp_pattern %||% TCP_PATTERN_DEFAULT)
  if (is.na(pat) || !nzchar(pat)) pat <- TCP_PATTERN_DEFAULT
  filt <- trk_filter_params()
  ms   <- .msd_settings(data.frame())                  # the MSD panel's window and weighting
  mopt <- list(window = ms$window, track_weight = ms$track_weight, min_pts = ms$min_pts)
  rows <- withProgress(message = "Measuring each immobile cell...", value = 0, {
    lapply(files, function(f) {
      incProgress(1 / length(files), detail = f)
      d <- .cmp_read_file(f)
      if (is.null(d)) d <- data.frame(tid = character(0), x = numeric(0), y = numeric(0),
                                      tim = numeric(0), eco = numeric(0))
      tryCatch(calib_cell(d, f, L, filt, msd = mopt,
                          n_vert = if (identical(pat, "triangle")) 3L else if (identical(pat, "square")) 4L else 6L),
               error = function(e) {
        r <- calib_cell(d[0, , drop = FALSE], f, L); r$note <- conditionMessage(e); r })
    })
  })
  curves <- lapply(rows, attr, "curve")
  tab <- do.call(rbind, lapply(rows, function(r) { attr(r, "curve") <- NULL; r }))
  attr(tab, "L") <- L; attr(tab, "pattern") <- pat; attr(tab, "window") <- ms$window
  attr(tab, "curves") <- stats::setNames(curves, tab$file)
  trk_cal_result(tab)
})

output$trk_cal_summary <- renderUI({
  tab <- trk_cal_result(); req(!is.null(tab))
  sm  <- calib_summary(tab)
  AMB <- "#B8860B"
  ideal <- .tcp_geom(attr(tab, "pattern"))$a_ideal
  tagList(
    tags$p(style = "margin:4px 0;",
      tags$b(sprintf("a = %.2f", sm$a)),
      sprintf("  (quartiles %.2f–%.2f, %d cell%s)  ·  ", sm$a_q25, sm$a_q75, sm$n_cells,
              if (sm$n_cells == 1L) "" else "s"),
      tags$b(sprintf("D_w = %.0f nm\u00b2/ms", sm$D_w)),
      sprintf("  (%.0f\u2013%.0f) over cycles %d\u2013%d, the MSD fit window", sm$Dw_q25, sm$Dw_q75,
              as.integer(attr(tab, "window")[1]), as.integer(attr(tab, "window")[2]))),
    helpText(style = "font-size:0.82em; color:#666; margin-top:0;",
      sprintf(paste0("\u03c3 includes the motion-blur term (R = %.3f per step, from the photon collection ",
                     "window); a would be %.2f\u2013%.2f if R were 30%% lower or higher. Over the first steps ",
                     "the immobile MSD rises at %.0f nm\u00b2/ms."),
              sm$R %||% NA_real_, sm$a_R_hi %||% NA_real_, sm$a_R_lo %||% NA_real_, sm$D_w_short %||% NA_real_)),
    helpText(style = "font-size:0.82em; color:#666; margin-top:0;",
      sprintf("L %.1f nm, %s. The ideal-estimator value for this pattern is %.2f.",
              attr(tab, "L"), attr(tab, "pattern"), ideal)),
    if (sm$n_cells < 3L) helpText(style = sprintf("color:%s; font-size:0.82em;", AMB),
      "Fewer than 3 cells: the quartiles say little about how well a is known."),
    if (any(nzchar(tab$note))) helpText(style = sprintf("color:%s; font-size:0.82em;", AMB),
      sprintf("%d file(s) gave no value; see the note column.", sum(nzchar(tab$note)))))
})

# The immobile MSD: every cell, their mean, the localization-noise level 4 sigma^2
# and the MSD fit window, on log-log axes.
output$trk_cal_msd_plot <- renderPlotly({
  tab <- trk_cal_result(); req(!is.null(tab))
  curves <- attr(tab, "curves"); win <- attr(tab, "window")
  req(length(Filter(Negate(is.null), curves)) > 0)
  dark   <- isTRUE(isolate(input$trk_dark_mode))
  bg_col <- if (dark) "#000000" else "#FFFFFF"
  fg_col <- if (dark) "#FFFFFF" else "#333333"
  fig <- plot_ly()
  for (nm in names(curves)) {
    cv <- curves[[nm]]; if (is.null(cv)) next
    ok <- cv$keep & is.finite(cv$msd) & cv$msd > 0
    fig <- fig %>% add_trace(x = cv$lag[ok], y = cv$msd[ok], type = "scatter", mode = "lines",
                             line = list(color = "rgba(140,140,140,0.45)", width = 1),
                             name = nm, legendgroup = "cells", showlegend = FALSE,
                             hovertemplate = paste0(nm, "<br>cycle %{x}<br>MSD %{y:.0f} nm\u00b2<extra></extra>"))
  }
  # mean of the cells' curves, each cell counting once
  cl <- Filter(function(c) !is.null(c) && nrow(c) > 0, curves)
  Lm <- max(vapply(cl, nrow, integer(1)))
  mm <- matrix(vapply(cl, function(c) { v <- rep(NA_real_, Lm); ok <- c$keep & is.finite(c$msd)
                                        v[c$lag[ok]] <- c$msd[ok]; v }, numeric(Lm)), nrow = Lm)
  pl <- data.frame(lag = seq_len(Lm), msd = rowMeans(mm, na.rm = TRUE))
  if (nrow(pl)) {
    okp <- is.finite(pl$msd) & pl$msd > 0
    fig <- fig %>% add_trace(x = pl$lag[okp], y = pl$msd[okp], type = "scatter", mode = "lines",
                             line = list(color = "#2c6fbb", width = 3), name = "mean of cells",
                             hovertemplate = "cycle %{x}<br>MSD %{y:.0f} nm\u00b2<extra>mean</extra>")
  }
  s2 <- stats::median(tab$sigma_nm, na.rm = TRUE)
  xr <- range(unlist(lapply(curves, function(c) if (is.null(c)) NULL else c$lag)), na.rm = TRUE)
  if (is.finite(s2))
    fig <- fig %>% add_trace(x = xr, y = rep(4 * s2^2, 2), type = "scatter", mode = "lines",
                             line = list(color = "#c03030", width = 1.5, dash = "dash"),
                             name = "4\u03c3\u00b2 (localization noise)", hoverinfo = "skip")
  shapes <- list(list(type = "rect", xref = "x", yref = "paper", x0 = log10(win[1]), x1 = log10(win[2]),
                      y0 = 0, y1 = 1, fillcolor = "rgba(44,111,187,0.08)", line = list(width = 0),
                      layer = "below"))
  fig %>% layout(
    xaxis = list(title = "lag (cycles)", type = "log", color = fg_col),
    yaxis = list(title = "immobile MSD (nm\u00b2)", type = "log", color = fg_col),
    shapes = shapes, paper_bgcolor = bg_col, plot_bgcolor = bg_col, font = list(color = fg_col),
    legend = list(orientation = "h", y = -0.25),
    annotations = list(
      list(x = log10(sqrt(win[1] * win[2])), y = 1, xref = "x", yref = "paper", yanchor = "bottom",
           text = "MSD fit window", showarrow = FALSE, font = list(size = 11, color = "#2c6fbb"))),
    uirevision = "trk_cal_msd")
})

output$trk_cal_results <- DT::renderDT({
  tab <- trk_cal_result(); req(!is.null(tab))
  show <- tab[, c("file", "n_tracks", "t_c_ms", "eco_med", "R", "msd1", "sigma_nm", "a", "D_w", "D_w_short", "note")]
  num  <- vapply(show, is.numeric, logical(1))
  show[num] <- lapply(show[num], function(v) signif(v, 3))
  DT::datatable(show, rownames = FALSE, selection = "none",
                options = list(dom = "t", paging = FALSE, scrollX = TRUE))
})

observeEvent(input$trk_cal_save, {
  tab <- trk_cal_result()
  nm  <- trimws(input$trk_cal_name %||% "")
  if (is.null(tab)) { showNotification("Measure first.", type = "warning", duration = 5); return(invisible(NULL)) }
  if (!nzchar(nm))  { showNotification("Give the profile a name.", type = "warning", duration = 5); return(invisible(NULL)) }
  sm <- calib_summary(tab)
  if (!sm$n_cells) { showNotification("No cell gave a value of a; nothing to save.", type = "error", duration = 6); return(invisible(NULL)) }
  num <- function(x) { v <- suppressWarnings(as.numeric(x)); if (length(v) && is.finite(v[1])) v[1] else NA_real_ }
  profs <- trk_cal_profiles()
  profs[[nm]] <- list(name = nm, a = sm$a, a_q25 = sm$a_q25, a_q75 = sm$a_q75,
                      D_w = sm$D_w, Dw_q25 = sm$Dw_q25, Dw_q75 = sm$Dw_q75, n_cells = sm$n_cells,
                      D_w_short = sm$D_w_short, window = as.integer(attr(tab, "window")),
                      R = sm$R, a_blur_lo = sm$a_R_hi, a_blur_hi = sm$a_R_lo,
                      pattern = attr(tab, "pattern"), L_nm = attr(tab, "L"),
                      dwell_us = num(input$trk_cal_dwell), photon_limit = num(input$trk_cal_pl),
                      damping = num(input$trk_cal_damping), note = trimws(input$trk_cal_note %||% ""),
                      date = format(Sys.Date()), files = tab$file[is.finite(tab$a)])
  trk_cal_profiles(profs)
  showNotification(sprintf("Profile '%s' saved (a = %.2f, D_w = %.0f nm²/ms).", nm, sm$a, sm$D_w),
                   type = "message", duration = 6)
})

# Keep the Acquisition Parameters selector in step with the profiles.
observeEvent(trk_cal_profiles(), {
  nms <- names(trk_cal_profiles())
  sel <- isolate(input$trk_cal_profile %||% "")
  updateSelectInput(session, "trk_cal_profile",
                    choices = c("None (uncalibrated)" = "", stats::setNames(nms, nms)),
                    selected = if (sel %in% nms) sel else "")
}, ignoreInit = TRUE)

output$trk_cal_profiles_tab <- DT::renderDT({
  profs <- trk_cal_profiles(); req(length(profs) > 0)
  f <- function(x, d = 2) if (is.null(x) || !is.finite(x)) "" else formatC(x, digits = d, format = "f")
  tab <- do.call(rbind, lapply(profs, function(p) data.frame(
    profile = p$name, a = sprintf("%s (%s–%s)", f(p$a), f(p$a_q25), f(p$a_q75)),
    D_w = f(p$D_w, 0),
    window = if (length(p$window) == 2L) sprintf("%d\u2013%d", as.integer(p$window[1]), as.integer(p$window[2])) else "",
    R = f(p$R, 3), cells = p$n_cells %||% NA, pattern = p$pattern %||% "", L_nm = f(p$L_nm, 1),
    dwell_us = f(p$dwell_us, 0), photon_limit = f(p$photon_limit, 0), damping = f(p$damping, 0),
    note = p$note %||% "", date = p$date %||% "", stringsAsFactors = FALSE)))
  DT::datatable(tab, rownames = FALSE, selection = "none",
                options = list(dom = "t", paging = FALSE, scrollX = TRUE))
})

output$trk_cal_profiles_dl <- downloadHandler(
  filename = function() paste0("calibration_profiles_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".json"),
  content  = function(file) writeLines(jsonlite::toJSON(unname(trk_cal_profiles()), pretty = TRUE,
                                                        auto_unbox = TRUE, na = "null", digits = NA), file))

observeEvent(input$trk_cal_profiles_ul, {
  up <- tryCatch(jsonlite::fromJSON(input$trk_cal_profiles_ul$datapath, simplifyVector = FALSE),
                 error = function(e) NULL)
  if (!is.null(up) && !is.null(up$name)) up <- list(up)          # a single profile
  ok <- Filter(function(p) is.list(p) && is.character(p$name) && nzchar(p$name) &&
                 is.numeric(p$a) && is.finite(p$a) && p$a > 0, up %||% list())
  if (!length(ok)) {
    showNotification("No valid calibration profile in that file.", type = "error", duration = 6)
    return(invisible(NULL))
  }
  num <- function(x) if (is.null(x)) NA_real_ else suppressWarnings(as.numeric(x))
  profs <- trk_cal_profiles()
  for (p in ok) {
    for (k in c("a", "a_q25", "a_q75", "D_w", "Dw_q25", "Dw_q75", "D_w_short", "R", "a_blur_lo", "a_blur_hi", "L_nm",
                "dwell_us", "photon_limit", "damping"))
      p[[k]] <- num(p[[k]])
    p$window <- suppressWarnings(as.integer(unlist(p$window)))
    p$files <- unlist(p$files)
    profs[[p$name]] <- p
  }
  trk_cal_profiles(profs)
  showNotification(sprintf("%d profile(s) loaded.", length(ok)), type = "message", duration = 5)
})

# Applying a profile sets a, D_w, the pattern and L; "None" returns to uncalibrated.
observeEvent(input$trk_cal_profile, {
  nm <- input$trk_cal_profile %||% ""
  if (!nzchar(nm)) {
    updateNumericInput(session, "trk_crlb_pref", value = NA)
    updateNumericInput(session, "trk_dw_floor", value = 0)
    return(invisible(NULL))
  }
  p <- trk_cal_profiles()[[nm]]; if (is.null(p)) return(invisible(NULL))
  updateNumericInput(session, "trk_crlb_pref", value = round(p$a, 3))
  if (is.finite(p$D_w %||% NA)) updateNumericInput(session, "trk_dw_floor", value = round(max(0, p$D_w)))
  pat <- .tcp_pattern_name(p$pattern %||% "")
  if (!is.na(pat) && nzchar(pat)) updateSelectInput(session, "trk_tcp_pattern", selected = pat)
  if (is.finite(p$L_nm %||% NA)) updateNumericInput(session, "trk_L_tcp", value = p$L_nm)
}, ignoreInit = TRUE)
