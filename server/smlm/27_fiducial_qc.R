## ─────────────────────────────────────────────────────────────────────────────
## server/smlm/27_fiducial_qc.R  -  panel: fiducial scan and exclusion
## ─────────────────────────────────────────────────────────────────────────────
## Part of smlm_server(). Not a standalone script: its statements are evaluated in
## the server function's environment, in file-name order (server/parts.R).
##
## The fiducial scan, its plot and table, and qc_exclude_from, which drops
## flagged traces from both data chains.
##
## INTERFACE (generated: Rscript tools/server_map.R --headers; do not edit by hand)
##   DEFINES for other parts:
##     fiducial_results, qc_exclude_from, smlm_nn_unit
##   READS from other parts:
##     01_view_upload: smlm_file_store
##     23_prepared_cache: deblink_on
## END INTERFACE
## ─────────────────────────────────────────────────────────────────────────────

# ==========================================================================
# QC 3 - Fiducial / bright-contaminant flagging (report only)
# ==========================================================================
# GOAL: flag traces that behave like gold fiducials / bright aggregates rather
# than single blinking molecules. Learned from a labelled pair of real datasets
# (one with a fiducial section, one clean), fiducials separate on TWO robust
# signatures, NOT on spread (spread does not distinguish them):
#   • BRIGHT   - efo is an extreme upper outlier. A gold bead reflects ~100× a
#                fluorophore (efo ~1e7 vs ~1e5), so log10(efo) sits many robust
#                MADs above the bulk. This is the single strongest signal.
#   • PERSISTENT + RECURRENT - the bead is re-localized far more than a blinking
#                dye (long trace) AND many such traces pile up at the SAME spot
#                (fiducials recur in one place; real labels are spread out).
# Scores are per-file and ROBUST (median + k·MAD on log; MAD tolerates up to ~50%
# contamination, so a heavy fiducial section does not hide itself). On the clean
# control file this flags 0 traces; on the fiducial file it flags a coherent,
# spatially-clustered set. REPORT ONLY unless "Exclude flagged" is ticked.
fiducial_results <- reactiveVal(NULL)
observeEvent(input$run_fiducial, {
  store <- smlm_file_store()
  if (!length(store)) { showNotification("No files loaded.", type = "warning"); return() }
  if (!isTRUE(input$fid_all_files)) {
    af <- input$active_file; if (!is.null(af) && af %in% names(store)) store <- store[af]
  }
  k_efo     <- max(2, input$fid_efo_z    %||% 4)      # brightness outlier (robust z)
  k_dur     <- max(1.5, input$fid_dur_x  %||% 5)      # duration, as a multiple of the file's q99
  min_locs  <- max(1, input$fid_min_locs %||% 5)      # min localizations for the brightness rule
  k_len     <- max(2, input$fid_len_z    %||% 4)      # length outlier (robust z), reported only
  recur_r   <- max(1, input$fid_recur_r  %||% 100)    # recurrence radius (nm), reported only
  recur_min <- max(1, input$fid_recur_min %||% 4)     # co-located traces, reported only
  rows <- list()
  withProgress(message = "Scanning for fiducials / contaminants...", value = 0, {
    nf <- length(store)
    for (i in seq_along(store)) {
      fname <- names(store)[[i]]; incProgress(1 / nf, detail = fname)
      # A single odd file must not abort the whole scan - skip it and continue.
      rows[[fname]] <- tryCatch(
        fiducial_scan_file(store[[fname]], fname, k_efo = k_efo, k_dur = k_dur,
                           min_locs = min_locs, recur_r = recur_r),
        error = function(e) NULL)
    }
  })
  rows <- Filter(Negate(is.null), rows)
  res <- if (length(rows)) do.call(rbind, rows) else NULL
  if (is.null(res) || !nrow(res)) {
    # Not an error - just nothing to look at. Store an empty (0-row) result so
    # the outputs show a calm "nothing found" state rather than a red message.
    fiducial_results(res %||% data.frame())
    showNotification("Scan complete - no traces to assess in the selected files.",
                     type = "message", duration = 5); return()
  }
  fiducial_results(res)
  nfl <- sum(res$flagged)
  showNotification(
    if (nfl == 0)
      sprintf("Scan complete - no fiducial/contaminant candidates found in %d trace(s) across %d file(s).",
              nrow(res), length(unique(res$dataset)))
    else
      sprintf("Flagged %d of %d traces (%.1f%%) as fiducial/contaminant candidates across %d file(s).",
              nfl, nrow(res), 100 * mean(res$flagged), length(unique(res$dataset))),
    type = "message", duration = 5)
})

output$fiducial_plot <- renderPlotly({
  res <- fiducial_results()
  # Calm empty state (no red validate) before a run or when nothing was found.
  if (is.null(res) || !nrow(res) || !"flagged" %in% names(res))
    # Declare the trace type (scatter/markers) so plotly does not print a
    # "No trace type specified" console warning for this empty placeholder.
    return(plotly::plot_ly(type = "scatter", mode = "markers") %>% plotly::layout(
      xaxis = list(visible = FALSE), yaxis = list(visible = FALSE),
      annotations = list(list(text = "Press 'Run' to scan for fiducials / contaminants.",
                              showarrow = FALSE, font = list(color = "#888", size = 13)))))
  # The two discriminating axes: brightness (efo, log - gold reflects ~100×) and
  # persistence (localizations per trace, log). Flagged traces sit in the bright
  # and/or long extremes. Falls back to spread on the x-axis when efo is absent.
  has_e <- "efo_med" %in% names(res) && any(is.finite(res$efo_med))
  xv <- if (has_e) pmax(res$efo_med, 1) else res$spread
  xl <- if (has_e) "Trace median EFO (log)" else "Trace spread (nm)"
  ok <- !res$flagged
  ht <- function(m) paste0("TID ", res$tid[m], "<br>", res$dataset[m],
    "<br>", res$n_loc[m], " locs · spread ", round(res$spread[m], 1), " nm",
    if (has_e) paste0("<br>efo ", format(round(res$efo_med[m]), big.mark = ","),
                      " (z ", round(res$z_efo[m], 1), ")") else "",
    "<br>len z ", round(res$z_len[m], 1),
    if ("recur" %in% names(res)) paste0(" · ", res$recur[m], " neighbours") else "")
  # Add each layer ONLY when it has points - an add_markers() with empty x/y is
  # what makes plotly warn "No trace type / no positional attributes specified".
  p <- plotly::plot_ly()
  if (any(ok))
    p <- p %>% plotly::add_markers(x = xv[ok], y = res$n_loc[ok], name = "traces",
      marker = list(color = "rgba(120,120,120,0.5)", size = 5,
                    line = list(color = "white", width = 0.3)),
      text = ht(ok), hoverinfo = "text")
  if (any(!ok))
    p <- p %>% plotly::add_markers(x = xv[!ok], y = res$n_loc[!ok], name = "flagged",
      marker = list(color = "#d62728", size = 8, symbol = "x"),
      text = ht(!ok), hoverinfo = "text")
  # If nothing was flagged, say so on the plot (the points are still shown).
  anns <- if (!any(!ok)) list(list(text = "No fiducials / contaminants detected",
    xref = "paper", yref = "paper", x = 0.5, y = 1, showarrow = FALSE,
    font = list(color = "#2a6", size = 13))) else NULL
  p %>% plotly::layout(
    xaxis = list(title = xl, type = if (has_e) "log" else "linear"),
    yaxis = list(title = "Localizations per trace", type = "log"),
    legend = list(orientation = "h"), annotations = anns)
})

output$fiducial_table <- renderTable({
  res <- fiducial_results()
  if (is.null(res) || !nrow(res) || !"flagged" %in% names(res)) return(NULL)
  f <- res[isTRUE(res$flagged) | res$flagged %in% TRUE, , drop = FALSE]
  if (!nrow(f)) return(data.frame(Message = "No traces flagged with the current thresholds."))
  f <- f[order(-f$z_efo, -f$n_loc), , drop = FALSE]
  utils::head(data.frame(
    dataset = f$dataset, tid = f$tid, n_loc = f$n_loc,
    efo_med = round(f$efo_med, 0), z_efo = round(f$z_efo, 1),
    z_len = round(f$z_len, 1), neighbours = f$recur,
    spread_nm = round(f$spread, 1), cfr_med = round(f$cfr_med, 3),
    check.names = FALSE), 200L)
}, striped = TRUE, spacing = "xs", width = "100%")

output$download_fiducial <- downloadHandler(
  filename = function() paste0("fiducial_scan_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
  content = function(file) {
    res <- fiducial_results(); shiny::validate(shiny::need(!is.null(res), "No results."))
    write.csv(res, file, row.names = FALSE)
  })

# Optional exclusion: when "Exclude flagged traces" is on, drop the traces the
# Fiducial scan flagged (for that file) from a data frame before it is analysed.
# NON-DESTRUCTIVE and reversible - the loaded data is untouched; unticking the box
# (or re-running the scan) restores everything. tid excludes are matched PER FILE
# (tids are only unique within a file), so `fname` must be the file the df is from.
# Applied both in the active-file pipeline and in every batch analysis.
qc_exclude_from <- function(df, fname) {
  if (!isTRUE(input$fid_exclude) || is.null(df) || !nrow(df) || !"tid" %in% names(df))
    return(df)
  fr <- tryCatch(fiducial_results(), error = function(e) NULL)
  if (is.null(fr) || !all(c("flagged","tid","dataset") %in% names(fr))) return(df)
  ex <- fr$tid[fr$flagged %in% TRUE & fr$dataset == fname]
  if (length(ex)) df[!(df$tid %in% ex), , drop = FALSE] else df
}

# ── Shared within-cluster NN statistics ─────────────────────────────────────
# Which unit the within-cluster NN is measured between: de-blinked molecules
# (primary `cluster`) when de-blinking is enabled and present, else trace IDs.
smlm_nn_unit <- function(df) {
  if (deblink_on() && "cluster" %in% names(df) &&
      any(!is.na(df$cluster) & df$cluster > 0L)) "cluster" else "tid"
}
