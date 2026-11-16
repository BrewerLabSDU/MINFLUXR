## ─────────────────────────────────────────────────────────────────────────────
## server/smlm/11_dcr_split.R  -  two-colour split by DCR
## ─────────────────────────────────────────────────────────────────────────────
## Part of smlm_server(). Not a standalone script: its statements are evaluated in
## the server function's environment, in file-name order (server/parts.R).
##
## Automatic DCR threshold from a two-component mixture, the pooled-DCR
## option, and the per-trace assignment confidence.
##
## INTERFACE (generated: Rscript tools/server_map.R --headers; do not edit by hand)
##   DEFINES for other parts:
##     dcr_conf_tbl, dcr_note_written, hist_channels_df
##   READS from other parts:
##     05_filter_pipeline: filtered_df
## END INTERFACE
## ─────────────────────────────────────────────────────────────────────────────

# ------ 2-color: auto-detect DCR split via GMM
dcr_gmm_result <- reactiveVal(NULL)

# DCR RANGE BOXES: WHO WROTE THE NUMBERS
# The boxes are filled by the GMM and can be typed over, and typing must switch
# on "Use typed DCR ranges". Shiny reports both kinds of change the same way, so
# every number the server sends is remembered here, and a value that is NOT on
# the list came from the keyboard. Several writes are kept, because the browser
# can echo an older write after a newer one has been sent.
# A plain environment: this is bookkeeping and must not trigger reactivity.
dcr_track <- new.env(parent = emptyenv())
dcr_track$seen   <- c(SMLM_DEFAULTS$dcr_group1, SMLM_DEFAULTS$dcr_group2)
dcr_track$pooled <- isTRUE(SMLM_DEFAULTS$use_pooled_dcr)

# GOAL: remember numbers the server is about to put in the boxes.
dcr_note_written <- function(...) {
  dcr_track$seen <- utils::tail(c(dcr_track$seen, as.numeric(unlist(list(...)))), 40L)
  invisible(NULL)
}

# GOAL: fit the GMM on the column the split uses and put the result in the boxes.
# Runs automatically (notify = FALSE) and from the Auto-detect button (notify =
# TRUE, which also reports the fit). When the fit fails the boxes are left as
# they are, and resolve_channels() falls back to them.
dcr_fit_and_fill <- function(notify = FALSE) {
  df <- tryCatch(filtered_df(), error = function(e) NULL)
  .col <- if (!is.null(df) && "tid" %in% names(df))
    dcr_split_column(df, input$use_pooled_dcr)
  if (is.null(.col)) {
    if (notify) {
      showNotification("No DCR data available for auto-detection.", type = "warning", duration = 5)
      dcr_gmm_result(list(success = FALSE, reason = "No DCR column in data"))
    }
    return(invisible(NULL))
  }
  dcr_meds <- df %>%
    dplyr::group_by(tid) %>%
    dplyr::summarise(dcr_med = stats::median(.data[[.col]], na.rm = TRUE), .groups = "drop") %>%
    dplyr::pull(dcr_med)

  res <- auto_detect_dcr_split(dcr_meds, col_used = .col)
  dcr_gmm_result(res)

  if (!isTRUE(res$success)) {
    if (notify)
      showNotification(paste("GMM auto-detect failed:", res$reason), type = "error", duration = 7)
    return(invisible(NULL))
  }

  split <- res$split
  dcr_note_written(c(0, split), c(split + 1e-4, 1))   # before sending, see dcr_track
  shinyWidgets::updateNumericRangeInput(session, "dcr_group1", value = c(0, split))
  shinyWidgets::updateNumericRangeInput(session, "dcr_group2", value = c(split + 1e-4, 1))
  if (notify) showNotification(
    paste0("Split set at DCR = ", split,
           "  |  G1 ", res$mean1, " (SD ", res$sd1, ")",
           "  |  G2 ", res$mean2, " (SD ", res$sd2, ")",
           "  |  d' = ", res$dprime, " (", res$quality, "), ",
           sprintf("%.1f%% expected crosstalk", 100 * res$err_rate),
           if (!identical(.col, "dcr_pooled"))
             "  |  using RAW dcr - try pooled DCR for a cleaner split" else ""),
    type = if (identical(res$quality, "poor")) "warning" else "message", duration = 9
  )
  invisible(NULL)
}

# The button always means "use the GMM": it leaves typed mode and refits.
observeEvent(input$auto_detect_dcr, {
  if (isTRUE(input$dcr_split_manual)) updateCheckboxInput(session, "dcr_split_manual", value = FALSE)
  dcr_fit_and_fill(notify = TRUE)
})

# Keep the boxes on the GMM's numbers whenever the GMM is in charge: on enabling
# the split, on a filter or file change, and on switching raw/pooled DCR.
observe({
  if (!isTRUE(input$enable_two_color) || isTRUE(input$color_by_fluo) ||
      isTRUE(input$dcr_split_manual)) return()
  dcr_fit_and_fill(notify = FALSE)
})

# A number the server did not write was typed, so the typed ranges take over.
# Incomplete input (a box being cleared) is ignored until it is a valid range.
observeEvent(list(input$dcr_group1, input$dcr_group2), {
  if (isTRUE(input$dcr_split_manual)) return()
  vals <- c(input$dcr_group1, input$dcr_group2)
  if (length(vals) != 4L || !all(is.finite(vals))) return()
  from_server <- vapply(vals, function(v) any(abs(v - dcr_track$seen) < 1e-6), logical(1))
  if (!all(from_server)) updateCheckboxInput(session, "dcr_split_manual", value = TRUE)
})

# Raw and pooled DCR are different scales, so ranges typed for one say nothing
# about the other. Changing the column hands the split back to the GMM, which
# refills the boxes for the new column. Compared with the last known value
# because this widget appears late (renderUI) and a settings load also sets it.
observeEvent(input$use_pooled_dcr, {
  now <- isTRUE(input$use_pooled_dcr)
  if (identical(now, dcr_track$pooled)) return()
  dcr_track$pooled <- now
  if (isTRUE(input$dcr_split_manual)) {
    updateCheckboxInput(session, "dcr_split_manual", value = FALSE)
    showNotification("The DCR column changed, so the typed ranges were replaced by a new GMM split.",
                     type = "message", duration = 6)
  }
})

# Per-trace split confidence: which individual molecules sit on the edge of the
# DCR split. The GMM summary only gives an aggregate ambiguous %, which cannot
# tell you WHICH molecules are uncertain - and those are exactly the ones that
# manufacture false short-range cross-channel pairs.
dcr_conf_tbl <- reactive({
  gmm <- dcr_gmm_result()
  if (is.null(gmm) || !isTRUE(gmm$success)) return(NULL)
  df <- tryCatch(filtered_df(), error = function(e) NULL)
  if (is.null(df) || !"tid" %in% names(df)) return(NULL)
  .col <- if (identical(gmm$col_used, "dcr_pooled") && "dcr_pooled" %in% names(df))
    "dcr_pooled" else if ("dcr" %in% names(df)) "dcr" else NULL
  if (is.null(.col)) return(NULL)

  med <- df %>%
    dplyr::group_by(tid) %>%
    dplyr::summarise(dcr_med = stats::median(.data[[.col]], na.rm = TRUE),
                     n_locs  = dplyr::n(), .groups = "drop") %>%
    dplyr::filter(is.finite(dcr_med))
  if (!nrow(med)) return(NULL)

  cf <- dcr_split_confidence(med$dcr_med, gmm)
  med$group      <- cf$group
  med$posterior  <- round(cf$post, 4)
  med$log_odds   <- round(cf$log_odds, 3)
  med$on_edge    <- cf$edge
  med$dcr_column <- .col
  med[order(med$posterior), , drop = FALSE]   # least certain first
})

output$dcr_conf_plot <- renderPlotly({
  tb  <- dcr_conf_tbl()
  gmm <- dcr_gmm_result()
  shiny::validate(shiny::need(!is.null(tb) && nrow(tb) > 0,
    "Run 'Auto-detect split (GMM)' to score each trace's assignment confidence."))
  cols <- ifelse(tb$on_edge, "#e67e22", ifelse(tb$group == 1L, "#1f77b4", "#d62728"))
  plotly::plot_ly(
    x = tb$dcr_med, y = tb$posterior, type = "scatter", mode = "markers",
    marker = list(color = cols, size = 6, opacity = 0.75,
                  line = list(color = "white", width = 0.4)),
    text = paste0("tid ", tb$tid, "<br>DCR ", round(tb$dcr_med, 4),
                  "<br>group ", tb$group, "<br>posterior ", round(tb$posterior, 3),
                  ifelse(tb$on_edge, "<br><b>on edge</b>", "")),
    hoverinfo = "text") %>%
    plotly::layout(
      title = list(text = sprintf(
        "Per-trace split confidence - %d of %d traces on the edge (posterior < %.2f), orange",
        sum(tb$on_edge, na.rm = TRUE), nrow(tb), gmm$post_cut), font = list(size = 12)),
      margin = list(t = 46),
      xaxis = list(title = paste0("Trace median ", gmm$col_used)),
      yaxis = list(title = "Posterior of assigned group", range = c(0.45, 1.02)),
      shapes = list(
        list(type = "line", x0 = gmm$split, x1 = gmm$split, y0 = 0.45, y1 = 1.02,
             line = list(color = "#666666", width = 1.5, dash = "dash")),
        list(type = "line", x0 = min(tb$dcr_med), x1 = max(tb$dcr_med),
             y0 = gmm$post_cut, y1 = gmm$post_cut,
             line = list(color = "#e67e22", width = 1, dash = "dot"))),
      showlegend = FALSE)
})

output$download_dcr_conf <- downloadHandler(
  filename = function() paste0("dcr_split_confidence_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
  content = function(file) {
    tb <- dcr_conf_tbl()
    shiny::validate(shiny::need(!is.null(tb) && nrow(tb) > 0, "Run the GMM auto-detect first."))
    write.csv(tb, file, row.names = FALSE)
  }
)

# Split summary + an explicit split-QUALITY verdict. The crosstalk numbers are
# what a cross-channel co-localization result has to be read against, so they are
# shown next to the split rather than buried: an apparent co-localization signal
# at short r is only meaningful if it is large compared with the contamination.
output$dcr_gmm_info <- renderUI({
  res <- dcr_gmm_result()
  typed_note <- if (isTRUE(input$dcr_split_manual))
    tags$div(style = "font-size:0.85em; margin-top:4px; font-weight:600;",
             "Using the typed ranges. Any GMM numbers shown here come from the last automatic fit.")
  if (is.null(res)) return(typed_note)
  if (!isTRUE(res$success))
    return(tagList(typed_note, tags$div(style = "font-size:0.85em; color:#c0392b;",
                                        paste("GMM failed:", res$reason))))

  col_lbl  <- if (identical(res$col_used, "dcr_pooled")) "pooled DCR" else "raw DCR"
  pooled   <- identical(res$col_used, "dcr_pooled")
  has_pooled_col <- {
    d <- tryCatch(filtered_df(), error = function(e) NULL)
    !is.null(d) && "dcr_pooled" %in% names(d)
  }
  qcol <- switch(res$quality, clean = "#1a8a1a", workable = "#b26a00", "#c0392b")
  qtxt <- switch(res$quality,
    clean    = "clean split",
    workable = "workable split",
    poor     = "POOR split - channels overlap heavily",
    "separation unknown")

  tagList(
    typed_note,
    tags$div(style = "font-size:0.85em; margin-top:4px;",
      tags$b("Split: "), res$split, " on ", col_lbl,
      tags$br(),
      "G1 mean ", res$mean1, " (SD ", res$sd1, ", ", round(100 * res$pi1), "%)",
      " | G2 mean ", res$mean2, " (SD ", res$sd2, ", ", round(100 * res$pi2), "%)",
      " | n = ", res$n_traces, " traces"),
    tags$div(style = paste0("font-size:0.85em; margin-top:3px; color:", qcol, "; font-weight:600;"),
      "d' = ", res$dprime, " - ", qtxt),
    tags$div(style = "font-size:0.82em; color:#555; margin-top:2px;",
      tags$b("Expected crosstalk: "),
      sprintf("%.1f%% of traces mis-assigned overall; group 1 is ~%.1f%% contaminated by dye 2, group 2 ~%.1f%% by dye 1. %.1f%% of traces are ambiguous (max posterior < %.2f).",
              100 * res$err_rate, 100 * res$contam1, 100 * res$contam2,
              100 * res$ambiguous_frac, res$post_cut)),
    tags$div(style = "font-size:0.8em; color:#7a5c00; margin-top:3px;",
      HTML(paste0(
        "Crosstalk pushes cross-channel NN / L₁₂ / g₁₂ <b>up at short r</b>, toward false ",
        "co-localization. Treat a co-localization signal as real only if it is large ",
        "compared with these numbers.",
        # Only suggest pooling where pooling is actually possible: it happens in the
        # .npy/.json readers. A .csv is a pyMINFLUX export of the final iteration and
        # has already been through their pipeline, so there is nothing to pool here.
        if (!pooled && has_pooled_col)
          paste0(" <b>This file has a pooled DCR column you are not using</b> - tick ",
                 "<i>Use pooled DCR</i> and re-upload; pooled DCR usually separates the dyes better.")
        else if (!pooled)
          paste0(" Splitting on raw DCR. If this came from <b>.npy/.json</b>, re-uploading with ",
                 "<i>Pool DCR</i> enabled (photon-weighted across the iteration sequence) usually ",
                 "sharpens the split; a <b>.csv</b> is already a pyMINFLUX final-iteration export, ",
                 "so there is nothing further to pool.")
        else ""))),
    # Keep this proportionate: the GMM is a reasonable estimator and the split is
    # what the data allows. The point is only that split QUALITY is not evidence of
    # two dyes existing - the experiment design is.
    tags$div(style = "margin-top:5px; padding:5px 8px; border-left:3px solid #b26a00; background:rgba(178,106,0,0.06); font-size:0.8em; color:#7f8c8d; line-height:1.45;",
      HTML(paste0(
        "<b>d' measures separation, not whether two dyes exist.</b>",
        " A two-component GMM will always produce 2 groups, even if there is only one dye."))))
  })

# GOAL: the filtered data with each localization's channel, for the histogram.
# Calls resolve_channels(), the function the plot's clustering pipeline uses, so
# the histogram shows the same split. A reactive of its own so that moving the
# bin slider or a color picker does not refit the GMM inside it. Returns NULL
# when the split is off or cannot be made; the histogram then draws plain bars.
hist_channels_df <- reactive({
  if (!isTRUE(input$enable_two_color)) return(NULL)
  df <- filtered_df()
  if (is.null(df) || !nrow(df)) return(NULL)
  tryCatch(resolve_channels(df, as_smlm_settings(input)), error = function(e) NULL)
})
