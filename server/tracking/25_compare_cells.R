# ============================================================================
# COMPARE CELLS - one value per file (= one cell), conditions compared across cells
# ============================================================================
# GOAL: label every loaded file with its condition and the settings of the
# sequence it was acquired with, analyse each file on its own (percell_summary,
# R/est_percell.R: own cycle time, own gap split, own L / pattern / a / D_w), and
# compare conditions with the cell as the unit of replication. Pooling tracks
# across cells would treat correlated tracks as independent and let one busy
# cell dominate a condition.

# Every loaded file, whichever view is active.
cmp_file_names <- reactive(names(trk_store()))

# One stored file, prepared the way trk_raw_loaded() / trk_raw_df() prepare the
# active view: character track ids, times in ms, the axial factor applied.
.cmp_read_file <- function(name) {
  st <- isolate(trk_store())[[name]]
  d  <- if (is.null(st)) NULL else read_chunk_file(st$path)
  if (is.null(d) || !nrow(d)) return(NULL)
  d$tid <- as.character(d$tid)
  for (cc in intersect(c("x", "y", "tim"), names(d))) d[[cc]] <- as.numeric(d[[cc]])
  d$tim <- d$tim * 1000                                   # s -> ms
  zf <- isolate(trk_z_applied())
  if (is.finite(zf) && zf > 0 && zf != 1 && "z" %in% names(d) && any(is.finite(d$z))) d$z <- d$z * zf
  d
}

cmp_tab     <- reactiveVal(NULL)   # per-file settings: file, condition, L_nm, pattern, a, D_w
cmp_results <- reactiveVal(NULL)   # per-cell numbers from the last run, with the settings used

# New upload: one row per file. Sequence settings start BLANK, and a blank means
# "the current Acquisition Parameters value", so a dataset acquired with one
# sequence follows the sequence file whenever it is loaded, and only settings
# that differ between files need entering. Adding or removing a file keeps the
# other rows.
observeEvent(cmp_file_names(), {
  nm  <- cmp_file_names()
  old <- isolate(cmp_tab())
  if (!length(nm)) { cmp_tab(NULL); cmp_results(NULL); return(invisible(NULL)) }
  tb <- data.frame(file = nm, condition = "", calibration = "", L_nm = NA_real_, pattern = "",
                   a = NA_real_, D_w = NA_real_, stringsAsFactors = FALSE)
  if (!is.null(old)) {                                    # adding/removing keeps other rows
    m <- match(tb$file, old$file); k <- which(!is.na(m))
    if (length(k)) tb[k, -1] <- old[m[k], -1]
  }
  cmp_tab(tb)
  r <- isolate(cmp_results())                             # drop results of removed files
  if (!is.null(r)) cmp_results(if (any(r$file %in% nm)) r[r$file %in% nm, , drop = FALSE] else NULL)
}, ignoreNULL = FALSE)

# The table with every blank filled from the current Acquisition Parameters:
# the values a run actually uses.
cmp_resolved <- reactive({
  tb <- cmp_tab(); if (is.null(tb)) return(NULL)
  num <- function(x, d) { v <- suppressWarnings(as.numeric(x)); if (length(v) && is.finite(v[1])) v[1] else d }
  pat <- tolower(input$trk_tcp_pattern %||% TCP_PATTERN_DEFAULT)
  if (!pat %in% names(TCP_PATTERN_GEOM)) pat <- TCP_PATTERN_DEFAULT
  # A named calibration profile fills its row's blanks first (it belongs to one
  # sequence), then the Acquisition Parameters fill what is still blank.
  profs <- trk_cal_profiles()
  for (i in seq_len(nrow(tb))) {
    pr <- profs[[tb$calibration[i]]]
    if (is.null(pr)) next
    if (!is.finite(tb$L_nm[i]) && is.finite(pr$L_nm %||% NA)) tb$L_nm[i] <- pr$L_nm
    if (!nzchar(tb$pattern[i])) { pp <- .tcp_pattern_name(pr$pattern %||% ""); if (!is.na(pp)) tb$pattern[i] <- pp }
    if (!is.finite(tb$a[i])   && is.finite(pr$a %||% NA))   tb$a[i]   <- pr$a
    if (!is.finite(tb$D_w[i]) && is.finite(pr$D_w %||% NA)) tb$D_w[i] <- max(0, pr$D_w)
  }
  tb$L_nm[!is.finite(tb$L_nm)] <- num(input$trk_L_tcp, DMAX_L_TCP_DEFAULT)
  tb$pattern[!nzchar(tb$pattern)] <- pat
  a_in <- suppressWarnings(as.numeric(input$trk_crlb_pref))
  for (i in which(!is.finite(tb$a)))                       # entered a, else the pattern's ideal
    tb$a[i] <- if (length(a_in) && is.finite(a_in[1]) && a_in[1] > 0) a_in[1]
               else .tcp_geom(tb$pattern[i])$a_ideal
  tb$D_w[!is.finite(tb$D_w)]   <- num(input$trk_dw_floor, 0)
  tb
})


output$trk_cmp_table <- DT::renderDT({
  tb <- cmp_tab(); req(!is.null(tb))
  DT::datatable(tb, rownames = FALSE, selection = "none",
                editable = list(target = "cell", disable = list(columns = 0)),
                options = list(dom = "t", paging = FALSE, scrollY = "260px",
                               scrollCollapse = TRUE, ordering = FALSE))
})

# Cell edits from the table. DT reports 0-based columns without row names.
observeEvent(input$trk_cmp_table_cell_edit, {
  tb <- cmp_tab(); info <- input$trk_cmp_table_cell_edit
  req(!is.null(tb), NROW(info) > 0)
  for (k in seq_len(nrow(info))) {
    col <- names(tb)[info$col[k] + 1L]
    v   <- .cmp_coerce(col, info$value[k])
    if (is.null(v)) {
      showNotification(sprintf("'%s' is not a valid value for %s.", info$value[k], col),
                       type = "warning", duration = 5)
      next
    }
    tb[info$row[k], col] <- v
  }
  cmp_tab(tb)
})

# Bulk fill: one column, every file whose name contains the text.
observeEvent(input$trk_cmp_fill_apply, {
  tb <- cmp_tab(); req(!is.null(tb))
  col <- input$trk_cmp_fill_col %||% "condition"
  v   <- .cmp_coerce(col, input$trk_cmp_fill_value %||% "")
  if (is.null(v)) {
    ok <- if (col == "pattern") paste0(paste(names(TCP_PATTERN_GEOM), collapse = ", "),
                                       ", or 3 / 4 / 6 vertices") else "a positive number"
    showNotification(sprintf("%s needs %s.", col, ok), type = "warning", duration = 5)
    return(invisible(NULL))
  }
  pat <- tolower(trimws(input$trk_cmp_fill_match %||% ""))
  hit <- if (nzchar(pat)) grepl(pat, tolower(tb$file), fixed = TRUE) else rep(TRUE, nrow(tb))
  if (!any(hit)) {
    showNotification("No file name contains that text.", type = "warning", duration = 5)
    return(invisible(NULL))
  }
  tb[hit, col] <- v
  cmp_tab(tb)
  showNotification(sprintf("%s set for %d file(s).", col, sum(hit)), type = "message", duration = 4)
})

# Save / load the settings table, so a dataset's conditions are entered once and
# travel with the analysis.
output$trk_cmp_table_dl <- downloadHandler(
  filename = function() paste0("cell_conditions_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
  content  = function(file) utils::write.csv(cmp_tab(), file, row.names = FALSE, na = ""))

observeEvent(input$trk_cmp_table_ul, {
  tb <- cmp_tab(); req(!is.null(tb))
  up <- tryCatch(utils::read.csv(input$trk_cmp_table_ul$datapath, stringsAsFactors = FALSE,
                                 check.names = FALSE, colClasses = "character"),
                 error = function(e) NULL)
  if (is.null(up) || !"file" %in% names(up)) {
    showNotification("The table needs a 'file' column.", type = "error", duration = 6)
    return(invisible(NULL))
  }
  m <- match(tb$file, up$file); n_bad <- 0L
  for (col in intersect(c("condition", "calibration", "L_nm", "pattern", "a", "D_w"), names(up)))
    for (i in which(!is.na(m))) {
      v <- .cmp_coerce(col, up[[col]][m[i]])
      if (is.null(v)) n_bad <- n_bad + 1L else tb[i, col] <- v
    }
  cmp_tab(tb)
  showNotification(sprintf("Settings loaded for %d of %d file(s)%s.", sum(!is.na(m)), nrow(tb),
                           if (n_bad) sprintf("; %d invalid value(s) skipped", n_bad) else ""),
                   type = if (any(!is.na(m))) "message" else "warning", duration = 6)
})

# Run: analyse every loaded file on its own, straight from the store, so the
# comparison does not depend on which view is active.
observeEvent(input$trk_cmp_run, {
  tb <- cmp_resolved(); req(!is.null(tb), nrow(tb) > 0)
  filt <- trk_filter_params()                             # the Display tab's filters
  nu <- suppressWarnings(as.numeric(input$trk_beam_fwhm %||% TCP_BEAM_FWHM_NM))
  if (!is.finite(nu) || nu < 0) nu <- 0
  rows <- withProgress(message = "Analysing each cell...", value = 0, {
    lapply(seq_len(nrow(tb)), function(i) {
      incProgress(1 / nrow(tb), detail = tb$file[i])
      di <- .cmp_read_file(tb$file[i])
      if (is.null(di)) di <- data.frame(tid = character(0), x = numeric(0), y = numeric(0), tim = numeric(0))
      r  <- tryCatch(percell_summary(di, tb$file[i], L = tb$L_nm[i], pattern = tb$pattern[i],
                                     a = tb$a[i], D_w = tb$D_w[i], nu = nu, filters = filt,
                                     msd = .msd_settings(di)),
                     error = function(e) {
                       x <- percell_summary(di[0, , drop = FALSE], tb$file[i], tb$L_nm[i])
                       x$note <- conditionMessage(e); x })
      if (nzchar(tb$calibration[i]) && is.null(trk_cal_profiles()[[tb$calibration[i]]]))
        r$note <- trimws(paste(r$note, sprintf("calibration '%s' is not loaded", tb$calibration[i])))
      cbind(r, calibration = tb$calibration[i], L_nm = tb$L_nm[i], pattern = tb$pattern[i],
            a = tb$a[i], D_w = tb$D_w[i], stringsAsFactors = FALSE)
    })
  })
  cmp_results(do.call(rbind, rows))
})

# Results joined to the CURRENT conditions, so relabelling needs no re-run; flags
# results whose sequence settings were edited after the run.
cmp_view <- reactive({
  r <- cmp_results(); tb <- cmp_resolved()
  if (is.null(r) || is.null(tb)) return(NULL)
  m <- match(r$file, tb$file)
  r$condition <- ifelse(is.na(m), "", tb$condition[m])
  stale <- any(!is.na(m) & (abs(r$L_nm - tb$L_nm[m]) > 1e-9 | r$pattern != tb$pattern[m] |
                              abs(r$a - tb$a[m]) > 1e-9 | abs(r$D_w - tb$D_w[m]) > 1e-9))
  attr(r, "stale") <- isTRUE(stale)
  r
})

output$trk_cmp_plot <- renderPlotly({
  r <- cmp_view(); req(!is.null(r))
  met <- input$trk_cmp_metric %||% "D_um2_s"
  d   <- r[nzchar(trimws(r$condition)) & is.finite(r[[met]]), , drop = FALSE]
  shiny::validate(shiny::need(nrow(d) > 0,
    "Give the files a condition in the table above to compare them."))
  st <- percell_condition_stats(d, met)
  xi <- match(d$condition, st$condition)
  set.seed(1L); jit <- (stats::runif(nrow(d)) - 0.5) * 0.3     # horizontal spread only
  dark   <- isTRUE(isolate(input$trk_dark_mode))
  bg_col <- if (dark) "#000000" else "#FFFFFF"
  fg_col <- if (dark) "#FFFFFF" else "#333333"
  ci_up  <- ifelse(is.finite(st$ci_hi), st$ci_hi - st$mean, 0)
  ci_dn  <- ifelse(is.finite(st$ci_lo), st$mean - st$ci_lo, 0)
  plot_ly() %>%
    add_trace(x = xi + jit, y = d[[met]], type = "scatter", mode = "markers",
              marker = list(size = 9, color = "rgba(70,130,200,0.75)",
                            line = list(color = "rgba(40,90,160,0.9)", width = 1)),
              text = d$file, name = "one cell",
              hovertemplate = "%{text}<br>%{y:.4g}<extra></extra>") %>%
    add_trace(x = seq_len(nrow(st)), y = st$mean, type = "scatter", mode = "markers",
              marker = list(symbol = "line-ew-open", size = 34, color = fg_col,
                            line = list(width = 3, color = fg_col)),
              error_y = list(type = "data", symmetric = FALSE, array = ci_up, arrayminus = ci_dn,
                             color = fg_col, thickness = 1.5, width = 10),
              name = "mean ± 95% CI across cells",
              hovertemplate = "mean %{y:.4g}<extra></extra>") %>%
    layout(
      xaxis = list(title = "", tickmode = "array", tickvals = seq_len(nrow(st)),
                   ticktext = sprintf("%s<br>(%d cell%s)", st$condition, st$n_cells,
                                      ifelse(st$n_cells == 1L, "", "s")),
                   range = c(0.4, nrow(st) + 0.6), color = fg_col, zeroline = FALSE),
      yaxis = list(title = PERCELL_METRICS[[met]], color = fg_col,
                   rangemode = "tozero"),
      paper_bgcolor = bg_col, plot_bgcolor = bg_col, font = list(color = fg_col),
      legend = list(orientation = "h", y = -0.22), uirevision = "trk_cmp")
})

# Condition table for the chosen metric, plus the flags that change how the
# numbers should be read.
output$trk_cmp_stats <- renderUI({
  r <- cmp_view(); req(!is.null(r))
  met <- input$trk_cmp_metric %||% "D_um2_s"
  st  <- percell_condition_stats(r, met)
  fmt <- function(x) ifelse(is.finite(x), formatC(x, digits = 3, format = "g"), "-")
  AMB <- "#B8860B"
  asg <- nzchar(trimws(r$condition))
  n_dw     <- sum(asg & is.finite(r$Dw_frac) & r$Dw_frac > 0.5)
  tagList(
    if (isTRUE(attr(r, "stale"))) helpText(style = sprintf("color:%s; font-size:0.85em;", AMB),
      "Sequence settings were edited after the last run: press Analyse each cell again."),
    if (nrow(st)) tags$table(class = "table table-condensed",
      style = "font-size:0.85em; margin:6px 0; width:auto;",
      tags$thead(tags$tr(lapply(c("condition", "cells", "mean", "SD", "95% CI"), tags$th))),
      tags$tbody(lapply(seq_len(nrow(st)), function(i) tags$tr(
        tags$td(st$condition[i]), tags$td(st$n_cells[i]), tags$td(fmt(st$mean[i])),
        tags$td(fmt(st$sd[i])), tags$td(sprintf("%s to %s", fmt(st$ci_lo[i]), fmt(st$ci_hi[i]))))))),
    if (nrow(st) && any(st$n_cells < 3L)) helpText(style = "font-size:0.82em; color:#666;",
      "A condition needs at least 3 cells for a meaningful CI."),
    if (n_dw) helpText(style = sprintf("color:%s; font-size:0.82em;", AMB), sprintf(
      "In %d cell(s) the wander floor is more than half of the measured D.", n_dw)),
    if (sum(!asg)) helpText(style = "font-size:0.82em; color:#666;", sprintf(
      "%d file(s) without a condition are listed below and left out of the comparison.", sum(!asg))))
})

output$trk_cmp_results <- DT::renderDT({
  r <- cmp_view(); req(!is.null(r))
  cols <- c("condition", "file", "calibration", "n_tracks", "len_med", "t_c_ms", "frac_single", "eco_med",
            "sigma_nm", "a", "D_msd_um2_s", "msd_lags", "D_um2_s", "D_corr_um2_s", "Dw_frac", "note")
  show <- r[, cols, drop = FALSE]
  num  <- vapply(show, is.numeric, logical(1))
  show[num] <- lapply(show[num], function(v) signif(v, 3))
  DT::datatable(show, rownames = FALSE, selection = "none",
                options = list(dom = "t", paging = FALSE, scrollX = TRUE,
                               scrollY = "300px", scrollCollapse = TRUE))
})

output$trk_cmp_results_dl <- downloadHandler(
  filename = function() paste0("per_cell_results_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
  content  = function(file) {
    r <- cmp_view(); if (is.null(r)) r <- data.frame()
    utils::write.csv(r, file, row.names = FALSE)
  })
