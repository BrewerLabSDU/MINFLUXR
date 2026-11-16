## ─────────────────────────────────────────────────────────────────────────────
## server/smlm/06_settings_download.R  -  download the main and the Median Clustering settings files
## ─────────────────────────────────────────────────────────────────────────────
## Part of smlm_server(). Not a standalone script: its statements are evaluated in
## the server function's environment, in file-name order (server/parts.R).
##
## The two download handlers that write settings JSON.
##
## INTERFACE (generated: Rscript tools/server_map.R --headers; do not edit by hand)
##   DEFINES for other parts:
##     (nothing)
##   READS from other parts:
##     04_settings_collect_apply: get_ui_settings
## END INTERFACE
## ─────────────────────────────────────────────────────────────────────────────

# ===========================  END SMLM FILTER PIPELINE  =====================

# --------- 7) SETTINGS LOAD / RESET ----------
output$download_ui_settings <- downloadHandler(
  filename = function() {
    base <- if (!is.null(input$active_file) && nzchar(input$active_file)) tools::file_path_sans_ext(input$active_file) else "smlm"
    paste0(base, "_filter_settings_", format(Sys.time(), "%Y%m%d-%H%M%S"), ".json")
  },
  content  = function(file) {
    settings <- get_ui_settings()
    writeLines(jsonlite::toJSON(settings, pretty=TRUE, auto_unbox=TRUE, na="null"), file, useBytes=TRUE)
  }
)
output$download_cluster2_settings <- downloadHandler(
  filename = function() {
    base <- if (!is.null(input$active_file) && nzchar(input$active_file)) tools::file_path_sans_ext(input$active_file) else "smlm"
    paste0(base, "_cluster2_settings_", format(Sys.time(), "%Y%m%d-%H%M%S"), ".json")
  },
  content = function(file) {
    s <- list(
      # Master toggle - without it, loading a Median Cluster settings file left
      # clustering switched OFF and none of the restored parameters took effect.
      enable_cluster2 = isTRUE(input$enable_cluster2),
      show_cluster2_knn = isTRUE(input$show_cluster2_knn),

      cluster2_method = input$cluster2_method %||% SMLM_DEFAULTS$cluster2_method,
      cluster2_xi = input$cluster2_xi,
      cluster2_prefilter = isTRUE(input$cluster2_prefilter),
      cluster2_prefilter_r = input$cluster2_prefilter_r,
      cluster2_prefilter_min = input$cluster2_prefilter_min,
      cluster2_prune = isTRUE(input$cluster2_prune),
      cluster2_prune_tol = input$cluster2_prune_tol,
      cluster2_core_only = isTRUE(input$cluster2_core_only),
      cluster2_min_members = input$cluster2_min_members,
      eps2 = input$eps2,
      minpts2 = input$minpts2,

      show_cluster2_shape = isTRUE(input$show_cluster2_shape),
      show_cluster2_isosurface = isTRUE(input$show_cluster2_isosurface),
      cluster2_iso_selected_only = isTRUE(input$cluster2_iso_selected_only),
      cluster2_iso_id = input$cluster2_iso_id,
      cluster2_iso_color = input$cluster2_iso_color,
      cluster2_iso_unicolor = isTRUE(input$cluster2_iso_unicolor),
      cluster2_iso_by_channel = isTRUE(input$cluster2_iso_by_channel),
      cluster2_iso_grid_res = input$cluster2_iso_grid_res,
      cluster2_iso_pad = input$cluster2_iso_pad,
      cluster2_iso_manual_offset = isTRUE(input$cluster2_iso_manual_offset),
      cluster2_iso_offset = input$cluster2_iso_offset,
      cluster2_iso_opacity = input$cluster2_iso_opacity,
      cluster2_iso_sigma = input$cluster2_iso_sigma,
      cluster2_iso_manual_sigma = isTRUE(input$cluster2_iso_manual_sigma),
      cluster2_iso_smooth_rel = input$cluster2_iso_smooth_rel,
      
      show_cluster2_shape_2D = isTRUE(input$show_cluster2_shape_2D),
      
      show_cluster2_metrics = isTRUE(input$show_cluster2_metrics),
      
      show_cluster2_pca_axes = isTRUE(input$show_cluster2_pca_axes),
      
      show_cluster2_labels = isTRUE(input$show_cluster2_labels),
      show_cluster2_plane_angle = isTRUE(input$show_cluster2_plane_angle),
      show_cluster2_plane_distance = isTRUE(input$show_cluster2_plane_distance),
      show_bending_radius = isTRUE(input$show_bending_radius),
      
      fit_method = input$fit_method,
      fit_weighted = isTRUE(input$fit_weighted),
      mle_center_mode = input$mle_center_mode,
      mle_center_lambda = input$mle_center_lambda,
      
      fit_type = input$fit_type,
      helix_turns = input$helix_turns,
      show_fit_tube = isTRUE(input$show_fit_tube),
      fit_tube_radius = input$fit_tube_radius,
      show_fit_distances = isTRUE(input$show_fit_distances),

      fit_opacity = input$fit_opacity,

      # Fit/metric options that were added after this save list was last updated
      compute_rmsd_metrics     = isTRUE(input$compute_rmsd_metrics),
      membrane_protein_options = isTRUE(input$membrane_protein_options),
      arc_span_deg             = input$arc_span_deg,
      max_dist_lines_nm        = input$max_dist_lines_nm,

      show_cluster2_nn_lines = isTRUE(input$show_cluster2_nn_lines),
      nn_k_cluster2 = input$nn_k_cluster2,
      # NN label styling
      show_c2_nn_labels  = isTRUE(input$show_c2_nn_labels),
      c2_label_font_size = input$c2_label_font_size,
      c2_label_offset    = input$c2_label_offset,
      distance_range_cluster2 = unname(input$distance_range_cluster2) # avoid jsonlite named-vector warning
    )
    
    jsonlite::write_json(s, path = file, auto_unbox = TRUE, pretty = TRUE)
  }
)
