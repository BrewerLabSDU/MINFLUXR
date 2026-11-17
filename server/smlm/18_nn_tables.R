## ─────────────────────────────────────────────────────────────────────────────
## server/smlm/18_nn_tables.R  -  nearest-neighbour tables and their downloads
## ─────────────────────────────────────────────────────────────────────────────
## Part of smlm_server(). Not a standalone script: its statements are evaluated in
## the server function's environment, in file-name order (server/parts.R).
##
## nn_table, the NN lines between Median Clusters, and the table downloads of
## the main plot.
##
## INTERFACE (generated: Rscript tools/server_map.R --headers; do not edit by hand)
##   DEFINES for other parts:
##     secondary_nn_lines
##   READS from other parts:
##     03_dataset_zcal: uploaded_df
##     10_triggers_stores: cluster2_table_store, plot_trigger
##     12_metadata_histogram: nn_table_cluster2
##     16_crop_range_ui: data
##     17_plottable_data: plottable_data
## END INTERFACE
## ─────────────────────────────────────────────────────────────────────────────

output$hist_range_ui <- renderUI({
  req(input$file, input$hist_variable)
  
  # Use RAW data to define min/max so the slider doesn't shrink/grow with filters
  df <- uploaded_df()
  varname <- input$hist_variable
  
  # Always reset to the full data range on re-render (new upload or variable change).
  # Preserving a stale curr value across uploads causes the histogram to appear pre-filtered.
  if (varname == "tid") {
    tlens <- as.numeric(table(df$tid))
    if (!length(tlens)) return(NULL)
    min_v <- min(tlens, na.rm = TRUE)
    max_v <- max(tlens, na.rm = TRUE)
    shinyWidgets::numericRangeInput(
      "hist_range", "Histogram Range:",
      value = c(min_v, max_v), min = min_v, max = max_v
    )
  } else {
    vals <- df[[varname]]
    if (is.null(vals) || !is.numeric(vals) || !length(vals)) return(NULL)
    vals <- vals[is.finite(vals)]
    if (!length(vals)) return(NULL)
    min_v <- min(vals, na.rm = TRUE)
    max_v <- max(vals, na.rm = TRUE)
    shinyWidgets::numericRangeInput(
      "hist_range", "Histogram Range:",
      value = c(min_v, max_v), min = min_v, max = max_v
    )
  }
})
# Enhanced NN table with all required metadata
nn_table <- reactive({
  req(plot_trigger())
  df <- plottable_data()
  
  scope <- if (is.null(input$nn_group_scope)) "all" else input$nn_group_scope
  
  if (isTRUE(input$enable_clustering)) {
    # cluster medians
    df_geo <- df %>%
      dplyr::filter(!is.na(cluster) & cluster != 0) %>%
      dplyr::group_by(cluster) %>%
      dplyr::summarise(
        id             = dplyr::first(cluster),
        tid            = dplyr::first(tid),
        x              = stats::median(x, na.rm = TRUE),
        y              = stats::median(y, na.rm = TRUE),
        z              = stats::median(z, na.rm = TRUE),
        color = dplyr::first(color),     
        cluster2_color = dplyr::first(cluster2_color),
        cluster2       = dplyr::first(cluster2),
        analysis_group = dplyr::first(analysis_group),
        channel = if ("channel" %in% names(df)) dplyr::first(channel) else NA_integer_,
        channel_color = if ("channel_color" %in% names(df)) dplyr::first(channel_color) else NA_character_,
        channel_label = if ("channel_label" %in% names(df)) dplyr::first(channel_label) else NA_character_,
        mplex_channel = if ("mplex_channel" %in% names(df)) dplyr::first(mplex_channel) else NA_integer_,
        .groups = "drop"
      )
    id_label <- "cluster"
  } else {
    # tid medians
    df_geo <- df %>%
      dplyr::group_by(tid) %>%
      dplyr::summarise(
        id             = dplyr::first(tid),
        x              = stats::median(x, na.rm = TRUE),
        y              = stats::median(y, na.rm = TRUE),
        z              = stats::median(z, na.rm = TRUE),
        color = dplyr::first(color),
        cluster2_color = dplyr::first(cluster2_color),
        cluster        = dplyr::first(cluster),
        cluster2       = dplyr::first(cluster2),
        analysis_group = dplyr::first(analysis_group),
        channel = if ("channel" %in% names(df)) dplyr::first(channel) else NA_integer_,
        channel_color = if ("channel_color" %in% names(df)) dplyr::first(channel_color) else NA_character_,
        channel_label = if ("channel_label" %in% names(df)) dplyr::first(channel_label) else NA_character_,
        mplex_channel = if ("mplex_channel" %in% names(df)) dplyr::first(mplex_channel) else NA_integer_,
        .groups = "drop"
      )
    id_label <- "tid"
  }

  if (nrow(df_geo) < 2) return(df_geo[0, ])

  coords <- as.matrix(df_geo[, c("x","y","z")])
  pm <- .nn_pairs_scoped(df_geo, coords, smlm_nn_options(input))
  if (!nrow(pm)) return(df_geo[0, ])

  rows <- list()
  for (r in seq_len(nrow(pm))) {
    i <- pm[r, 1]; ni <- pm[r, 2]
    {
      # ---- HERE: apply within/between/all filter
      g1 <- df_geo$analysis_group[i]
      g2 <- df_geo$analysis_group[ni]
      if (!pair_passes_scope(g1, g2, scope)) next

      dist <- pm[r, 3]
      if (dist < input$distance_range[1] || dist > input$distance_range[2]) next

      # within-cluster2 filter
      if (isTRUE(input$nn_within_cluster2) && isTRUE(input$enable_cluster2)) {
        c2i <- df_geo$cluster2[i]; c2ni <- df_geo$cluster2[ni]
        if (is.na(c2i) || is.na(c2ni) || c2i == 0 || c2ni == 0 || c2i != c2ni) next
      }

      rows[[length(rows) + 1]] <- data.frame(
        source_id  = df_geo$id[i],
        target_id  = df_geo$id[ni],
        distance   = round(dist, 2),
        source_x   = df_geo$x[i],  source_y = df_geo$y[i],  source_z = df_geo$z[i],
        target_x   = df_geo$x[ni], target_y = df_geo$y[ni], target_z = df_geo$z[ni],
        source_cluster  = if ("cluster" %in% names(df_geo)) df_geo$cluster[i] else NA,
        target_cluster  = if ("cluster" %in% names(df_geo)) df_geo$cluster[ni] else NA,
        source_cluster2 = df_geo$cluster2[i],
        target_cluster2 = df_geo$cluster2[ni],
        stringsAsFactors = FALSE
      )
    }
  }
  out <- if (length(rows)) dplyr::bind_rows(rows) else df_geo[0, ]
  names(out)[names(out) == "source_id"] <- paste0("source_", id_label)
  names(out)[names(out) == "target_id"] <- paste0("target_", id_label)
  out
})

output$download_nn_table <- downloadHandler(
  filename = function() paste0("nn_distances_", Sys.Date(), ".xlsx"),
  content  = function(file) {
    out <- list(distances = nn_table())
    # Second sheet when channels exist: one row per drawn pair, carrying which
    # channel each end belongs to. This replaces the separate "2-color NN table"
    # export, which was the same thing computed by a second code path.
    ed <- tryCatch(build_nn_segments(plottable_data(), smlm_nn_options(input), is3d = TRUE)$edges,
                   error = function(e) NULL)
    if (is.data.frame(ed) && nrow(ed) &&
        length(unique(stats::na.omit(c(ed$ch1, ed$ch2)))) > 1L)
      out$channel_pairs <- ed[, c("x1","y1","z1","x2","y2","z2",
                                  "dist","ch1","ch2","lab1","lab2")]
    openxlsx::write.xlsx(out, file, row.names = FALSE)
  },
  contentType = "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
)

# Secondary cluster median-to-median NN lines (green)
secondary_nn_lines <- reactive({
  df <- plottable_data()
  req(isTRUE(input$enable_cluster2))
  
  # majority group per cluster2 (for within/between filtering)
  mode_int <- function(x) {
    x <- x[is.finite(x)]
    if (!length(x)) return(NA_integer_)
    as.integer(names(which.max(table(x))))
  }
  
  df_geo <- df %>%
    dplyr::filter(!is.na(cluster2) & cluster2 != 0) %>%
    dplyr::group_by(cluster2) %>%
    dplyr::summarise(
      x = stats::median(x, na.rm = TRUE),
      y = stats::median(y, na.rm = TRUE),
      z = stats::median(z, na.rm = TRUE),
      color = dplyr::first(color),    
      analysis_group = mode_int(analysis_group),
      .groups = "drop"
    )
  if (nrow(df_geo) < 2) return(list())
  
  scope  <- if (is.null(input$nn_group_scope)) "all" else input$nn_group_scope
  coords <- as.matrix(df_geo[, c("x","y","z")])
  tmp    <- safe_nn2(coords, k = input$nn_k_cluster2 + 1)
  if (is.null(tmp$nn)) return(list())
  
  nn    <- tmp$nn
  usedK <- max(0, tmp$k_eff - 1)
  
  lines <- list()
  seen  <- new.env(parent = emptyenv())
  for (i in seq_len(nrow(coords))) {
    for (j in 1:usedK) {
      ni <- nn$nn.idx[i, j + 1]
      if (is.na(ni)) next
      
      # de-duplicate undirected pairs
      key <- paste0(min(i, ni), "_", max(i, ni))
      if (!is.null(seen[[key]])) next
      seen[[key]] <- TRUE
      
      # apply within/between/all filter
      g1 <- df_geo$analysis_group[i]
      g2 <- df_geo$analysis_group[ni]
      if (!pair_passes_scope(g1, g2, scope)) next
      
      dist <- as.numeric(nn$nn.dists[i, j + 1])
      if (dist < input$distance_range_cluster2[1] || dist > input$distance_range_cluster2[2]) next

      nn_off <- suppressWarnings(as.numeric(input$c2_label_offset %||% 5))
      if (!is.finite(nn_off) || nn_off < 0) nn_off <- 5
      midpoint <- perp_offset_xy(coords[i, ], coords[ni, ], nn_off)
      lines[[length(lines) + 1]] <- list(
        x = c(coords[i, 1], coords[ni, 1]),
        y = c(coords[i, 2], coords[ni, 2]),
        z = c(coords[i, 3], coords[ni, 3]),
        text = round(dist, 2),
        midpoint = midpoint
      )
    }
  }
  lines
})


# 3) Plane distance table ------------ 
output$download_plane_distance_table <- downloadHandler(
  filename = function() paste0("cluster2_xy_plane_distances_", Sys.Date(), ".xlsx"),
  content  = function(file) {
    df <- plottable_data()
    tbl <- compute_cluster2_plane_distances(df)
    openxlsx::write.xlsx(
      x         = tbl,
      file      = file,
      sheetName = "Plane distances",
      rowNames  = FALSE
    )
  },
  contentType = "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
)

# 4) Cluster2 metrics ---------------- 
output$download_cluster2_table <- downloadHandler(
  filename = function() paste0("cluster2_fit_metrics_", Sys.Date(), ".xlsx"),
  content  = function(file) {
    tbl <- cluster2_table_store()
    if (is.null(tbl) || nrow(tbl) == 0) stop("No Median Cluster metrics to download.")
    openxlsx::write.xlsx(
      x         = tbl,
      file      = file,
      sheetName = "Median Cluster metrics",
      rowNames  = FALSE
    )
  },
  contentType = "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
)

# 5) Cluster2 NN table --------------- 
output$download_nn_table2 <- downloadHandler(
  filename = function() paste0("nn_distances_cluster2_", Sys.Date(), ".xlsx"),
  content  = function(file) {
    openxlsx::write.xlsx(
      x         = nn_table_cluster2(),
      file      = file,
      sheetName = "NN Median Cluster",
      rowNames  = FALSE
    )
  },
  contentType = "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
)

df_full_raw <- reactive({
  data()
})
