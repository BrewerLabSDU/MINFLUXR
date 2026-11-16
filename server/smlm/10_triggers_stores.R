## ─────────────────────────────────────────────────────────────────────────────
## server/smlm/10_triggers_stores.R  -  plot trigger and small result stores
## ─────────────────────────────────────────────────────────────────────────────
## Part of smlm_server(). Not a standalone script: its statements are evaluated in
## the server function's environment, in file-name order (server/parts.R).
##
## plot_trigger (what makes plottable_data recompute), diagnostics and
## isosurface stores, the fit log, and the observers that tie display switches
## to a replot.
##
## INTERFACE (generated: Rscript tools/server_map.R --headers; do not edit by hand)
##   DEFINES for other parts:
##     cluster2_iso_diag, cluster2_table_store, fit_log_rv, iso_geom,
##     iso_mesh_cache, plot_trigger, smlm_diag, smlm_diag_mats
##   READS from other parts:
##     01_view_upload: uploaded_df_rv
##     17_plottable_data: plottable_data
## END INTERFACE
## ─────────────────────────────────────────────────────────────────────────────

# --------- 8) TRIGGERS & STORES ----------
plot_trigger <- reactiveVal(NULL)
smlm_diag      <- reactiveVal(NULL)   # list(data_clusters, img_blobs) in data coords
smlm_diag_mats <- reactiveVal(NULL)  # list(raw=mat_c, smooth=mat_cs) for debug plot
# One row per isosurface drawn in the last 3-D render, so the controls can
# report the geometry they picked. See iso_status_ui().
cluster2_iso_diag <- reactiveVal(NULL)

# ==========================================================================
# ISOSURFACE: COMMITTED GEOMETRY + MESH CACHE
# ==========================================================================
# Wrapping every Median Cluster is the most expensive thing the 3-D view does,
# and in 2-colour mode it runs over both channels. Two separate problems made
# that painful, and they need different fixes:
#
#   1. DRAGGING A SLIDER re-ran the whole wrap on every intermediate value. The
#      geometry sliders are therefore COMMITTED, exactly like the clustering
#      parameters: they are snapshotted when Plot Data is pressed and the view
#      ignores them until then.
#   2. ANY re-render rebuilt the meshes, because a renderPlotly body re-runs on
#      all of its reactive dependencies - so changing point size or colour mode,
#      which have nothing to do with the surfaces, still paid for them. The
#      cache keys each mesh on its geometry alone, so those re-renders reuse it.
#
# The split follows the control panel: Shape is committed, Appearance is live.
# Opacity and colour stay instant because they are cache hits, not rebuilds.
iso_mesh_cache <- new.env(parent = emptyenv())

iso_geom <- eventReactive(plot_trigger(), {
  # A new commit invalidates every cached mesh. Clearing here rather than on a
  # timer keeps the cache bounded to one settings epoch.
  rm(list = ls(iso_mesh_cache), envir = iso_mesh_cache)
  list(
    c2_grid_res     = input$cluster2_iso_grid_res,
    c2_pad          = input$cluster2_iso_pad,
    c2_smooth_rel   = input$cluster2_iso_smooth_rel,
    c2_manual_off   = isTRUE(input$cluster2_iso_manual_offset),
    c2_offset       = input$cluster2_iso_offset,
    c2_manual_sigma = isTRUE(input$cluster2_iso_manual_sigma),
    c2_sigma        = input$cluster2_iso_sigma
  )
})
output$cluster2_iso_status <- renderUI(iso_status_ui(cluster2_iso_diag()))
cluster2_table_store <- reactiveVal(NULL)
# Cache of clustered data per file - populated whenever a file is plotted.
# Enables batch analysis across all previously-plotted files.
smlm_plotted_store <- reactiveVal(list())
shinyjs::hide("scatter2D")
observeEvent(plot_trigger(), {
  shinyjs::show("scatter2D")
}, ignoreNULL = TRUE, once = TRUE)
observeEvent(input$plot_data, { plot_trigger(Sys.time()) })

output$plot_bar_file_label <- renderUI({
  nm <- input$active_file %||% ""
  if (!nzchar(nm)) return(NULL)
  tags$span(style = "font-size: 13px; color: #555;", nm)
})

# view_plane_2d is now a static input (ui_smlm.R). The old renderUI dropped
# XZ/YZ from the choices when `any(is.finite(df$z))` was FALSE - which never
# fired for a 2D MINFLUX file, because those carry z = 0 and 0 IS finite. So
# side views were offered on planar data and produced a flat line. The planar
# observer below disables them on z VARIATION instead.
observeEvent(input$enable_clustering, {
  req(input$file)
  plot_trigger(Sys.time())
}, ignoreInit = TRUE)
fit_log_rv <- reactiveVal("")

# Keep smlm_plotted_store in sync: every time a file is plotted, cache its
# clustered data so batch analysis can iterate all plotted files.
observeEvent(plottable_data(), {
  df  <- tryCatch(plottable_data(), error = function(e) NULL)
  nm  <- isolate(input$active_file) %||% isolate(input$file)[[1]]
  if (is.null(df) || !nrow(df) || is.null(nm)) return()
  store <- smlm_plotted_store()
  store[[nm]] <- df
  smlm_plotted_store(store)
}, ignoreNULL = TRUE)


output$fit_log <- renderText({
  fit_log_rv()
})

# clear the box whenever metrics are disabled
observeEvent(
  list(input$show_cluster2_metrics, input$enable_cluster2),
  {
    if (!isTRUE(input$enable_cluster2) || !isTRUE(input$show_cluster2_metrics)) {
      fit_log_rv("")
    }
  },
  ignoreInit = TRUE
)
observeEvent(
  list(input$remove_rogues, input$view_SD),
  {
    if (!isTRUE(input$remove_rogues) || !isTRUE(input$view_SD)) {
      # only clear if the current message is the rogue-removal one
      cur <- fit_log_rv()
      if (is.character(cur) && length(cur) == 1 && grepl("^Rogue removal:", cur)) {
        fit_log_rv("")
      }
    }
  },
  ignoreInit = TRUE
)
observeEvent(input$enable_two_color, {
  if (isTRUE(input$enable_two_color) && !is.null(input$color_mode) && input$color_mode != "default") {
    updateSelectInput(session, "color_mode", selected = "default")
    showNotification(
      "2-color mode uses Default coloring. If you switch to a numeric color mode (z/time/density), that will override the 2-color colors.",
      type = "message", duration = 6
    )
  }
}, ignoreInit = TRUE)

# Show "Use Pooled DCR" checkbox only when dcr_pooled is present in the data
output$use_pooled_dcr_ui <- renderUI({
  df <- uploaded_df_rv()
  if (!is.null(df) && "dcr_pooled" %in% names(df)) {
    checkboxInput("use_pooled_dcr",
                  "Use Pooled DCR for color splitting",
                  value = SMLM_DEFAULTS$use_pooled_dcr)
  }
})
