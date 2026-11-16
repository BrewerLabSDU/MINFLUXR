## ─────────────────────────────────────────────────────────────────────────────
## server/smlm/30_radius_gyration.R  -  panel: radius of gyration
## ─────────────────────────────────────────────────────────────────────────────
## Part of smlm_server(). Not a standalone script: its statements are evaluated in
## the server function's environment, in file-name order (server/parts.R).
##
## Run handler, histogram, small-N note and download.
##
## INTERFACE (generated: Rscript tools/server_map.R --headers; do not edit by hand)
##   DEFINES for other parts:
##     build_rog_hist_fig, rog_results
##   READS from other parts:
##     01_view_upload: smlm_file_store
##     20_dataset_groups: apply_dataset_groups, series_palette
##     23_prepared_cache: .small_n_note, cluster2_label, deblink_label,
##                        deblink_on, prepared_df
## END INTERFACE
## ─────────────────────────────────────────────────────────────────────────────

## ── Radius of Gyration ────────────────────────────────────────────────────
rog_results <- reactiveVal(NULL)
rog_params  <- reactiveVal(NULL)   # compute params, frozen at Run

observeEvent(input$run_rog, {
  store <- smlm_file_store()
  if (!length(store)) {
    showNotification("No files loaded.", type = "warning"); return()
  }
  # mode drives BOTH what we group by and what a "point" is:
  #   "cluster2"  structure size - group by median cluster, one point per
  #               de-blinked MOLECULE when de-blinking is on, else one point per
  #               trace ID. (The auto fallback is the whole point of this mode.)
  #   "cluster"   label footprint - group by de-blinked molecule, one point per trace.
  mode        <- input$rog_cluster_col %||% "cluster2"
  cluster_col <- if (mode == "cluster") "cluster" else "cluster2"
  store <- smlm_file_store()
  if (!length(store)) {
    showNotification("No files loaded.", type = "warning"); return()
  }
  all_res <- list()
  withProgress(message = "Computing RoG...", value = 0, {
    n_files <- length(store)
    for (i in seq_along(store)) {
      fname <- names(store)[[i]]
      incProgress(1 / n_files, detail = fname)
      df <- prepared_df(fname)   # cached filter -> thin -> cluster chain
      if (is.null(df) || !nrow(df) || !cluster_col %in% names(df)) next

      # Split AFTER preparation, so the channels come from resolve_channels()
      # and the split covers a 2-colour DCR split as well as a multiplex build.
      # A cluster belongs to exactly one channel by construction, so this only
      # partitions the results - it never recomputes them differently.
      parts <- if (isTRUE(input$rog_by_channel)) channel_parts(df, fname)
               else list(list(df = df, tag = fname))

      for (pt in parts) {
        # Only the plain structure-size mode collapses to one point per de-blinked
        # molecule (when de-blinking is on); the other two use raw traces as points.
        member_col <- if (mode == "cluster2" && deblink_on() &&
                          "cluster" %in% names(pt$df) && any(pt$df$cluster > 0, na.rm = TRUE))
                        "cluster" else "tid"
        r <- tryCatch(compute_cluster_rog(pt$df, cluster_col = cluster_col,
                                          member_col = member_col),
                      error = function(e) NULL)
        if (is.null(r) || !nrow(r)) next
        r$dataset <- pt$tag
        all_res[[pt$tag]] <- r
      }
    }
  })
  res <- do.call(rbind, Filter(Negate(is.null), all_res))

  if (is.null(res) || !nrow(res)) {
    hint <- if (cluster_col == "cluster")
      sprintf("Enable %s and press Plot for each file.", deblink_label())
    else
      sprintf("Enable %s and press Plot for each file.", cluster2_label())
    showNotification(paste0("No RoG computed. ", hint), type = "warning")
    rog_results(NULL); return()
  }
  rog_results(res)
  # Compute parameters frozen at Run: reading them live relabelled an existing
  # plot without recomputing it, so the axis could describe a different unit
  # than the data on screen.
  rog_params(list(mode = mode))
  showNotification(
    sprintf("RoG computed for %d cluster(s) across %d file(s).",
            nrow(res), length(unique(res$dataset))),
    type = "message", duration = 3)
})

build_rog_hist_fig <- function() {
  res <- rog_results()
  shiny::validate(shiny::need(!is.null(res) && nrow(res) > 0, "Run RoG first."))
  res <- apply_dataset_groups(res)
  shiny::validate(shiny::need(nrow(res) > 0, "No files assigned to groups (group mode)."))

  density_mode <- identical(input$rog_yaxis, "density")
  nbins     <- max(5L, as.integer(input$rog_bins %||% 30L))
  aggregate <- isTRUE(input$rog_aggregate)
  med_v     <- stats::median(res$rog, na.rm = TRUE)

  fig    <- plotly::plot_ly()
  shapes <- list()

  if (aggregate) {
    col <- "rgba(31,119,180,0.65)"
    fig <- fig %>% plotly::add_histogram(
      x        = res$rog,
      name     = sprintf("All datasets (N clusters = %d)", nrow(res)),
      nbinsx   = nbins,
      histnorm = if (density_mode) "probability density" else "",
      marker   = list(color = col, line = list(color = "rgba(31,119,180,1)", width = 1)),
      opacity  = 0.75
    )
    shapes <- list(list(
      type = "line", x0 = med_v, x1 = med_v, y0 = 0, y1 = 1, yref = "paper",
      line = list(color = "firebrick", width = 2, dash = "dash")
    ))
  } else {
    datasets <- if ("dataset" %in% names(res)) unique(res$dataset) else NULL
    n_ds     <- length(datasets)
    pal <- if (n_ds >= 1) series_palette(datasets) else NULL

    if (!is.null(pal)) {
      for (ds in datasets) {
        sub     <- res[res$dataset == ds, , drop = FALSE]
        col_hex <- pal[ds]
        med_ds  <- stats::median(sub$rog, na.rm = TRUE)
        fig <- fig %>% plotly::add_histogram(
          x        = sub$rog,
          name     = ds,
          nbinsx   = nbins,
          histnorm = if (density_mode) "probability density" else "",
          marker   = list(color = paste0(col_hex, "99"),
                          line  = list(color = col_hex, width = 0.5)),
          opacity  = 0.75
        )
        shapes <- c(shapes, list(list(
          type = "line", x0 = med_ds, x1 = med_ds, y0 = 0, y1 = 1, yref = "paper",
          line = list(color = col_hex, width = 2, dash = "dash"), opacity = 0.9
        )))
      }
    } else {
      fig <- fig %>% plotly::add_histogram(
        x        = res$rog,
        nbinsx   = nbins,
        histnorm = if (density_mode) "probability density" else "",
        marker   = list(color = "rgba(31,119,180,0.65)",
                        line  = list(color = "rgba(31,119,180,1)", width = 1)),
        name     = "RoG"
      )
      shapes <- list(list(
        type = "line", x0 = med_v, x1 = med_v, y0 = 0, y1 = 1, yref = "paper",
        line = list(color = "firebrick", width = 2, dash = "dash")
      ))
    }
  }

  fig <- fig %>%
    plotly::layout(
      barmode = "overlay",
      xaxis = list(title = if (identical((rog_params()$mode) %||% (input$rog_cluster_col %||% "cluster2"), "cluster"))
                             "Label footprint Rg (nm)" else "Median-cluster Rg - structure size (nm)"),
      yaxis = list(title = if (density_mode) "Density" else "Count"),
      shapes = shapes,
      annotations = list(
        list(x = med_v, y = 1, yref = "paper",
             text = sprintf(" overall median: %.1f nm", med_v),
             showarrow = FALSE, xanchor = "left",
             font = list(color = "firebrick", size = 11))
      ),
      legend = list(orientation = "h", y = -0.18),
      paper_bgcolor = "white", plot_bgcolor = "white"
    )
  fig
}
output$rog_hist <- renderPlotly(build_rog_hist_fig())

output$rog_small_n <- renderUI({
  res <- rog_results()
  if (is.null(res) || !nrow(res)) return(NULL)
  nv <- if ("n_units" %in% names(res)) res$n_units else
        if ("n_traces" %in% names(res)) res$n_traces else NULL
  if (is.null(nv)) return(NULL)
  .small_n_note(nv, "structure")
})

output$download_rog <- downloadHandler(
  filename = function() paste0("rog_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
  content = function(file) {
    res <- rog_results()
    shiny::validate(shiny::need(!is.null(res) && nrow(res) > 0, "No results."))
    write.csv(res, file, row.names = FALSE)
  }
)
