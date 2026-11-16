## ─────────────────────────────────────────────────────────────────────────────
## server/smlm/19_scatter3d.R  -  main 3D plot
## ─────────────────────────────────────────────────────────────────────────────
## Part of smlm_server(). Not a standalone script: its statements are evaluated in
## the server function's environment, in file-name order (server/parts.R).
##
## output$scatter3D: one render function. It also fits the shapes and
## isosurfaces it draws and fills the Median Cluster table store.
##
## INTERFACE (generated: Rscript tools/server_map.R --headers; do not edit by hand)
##   DEFINES for other parts:
##     (nothing)
##   READS from other parts:
##     03_dataset_zcal: uploaded_df
##     10_triggers_stores: cluster2_iso_diag, cluster2_table_store,
##                         fit_log_rv, iso_geom, iso_mesh_cache, plot_trigger
##     14_gt_estimator: gt_estimate
##     16_crop_range_ui: add_3d_crop_preview, data
##     17_plottable_data: plottable_data
##     18_nn_tables: secondary_nn_lines
## END INTERFACE
## ─────────────────────────────────────────────────────────────────────────────

output$scatter3D <- renderPlotly({
  req(input$file)
  req(plot_trigger(), input$mode3d == "3D")
  
  df <- plottable_data()
  shiny::validate(shiny::need(nrow(df) > 0, "No points to display - please check filters or upload a valid file."))
  
  df_raw <- uploaded_df()
  
  # Match coordinate shift if make_positive was applied
  if (isolate(input$make_positive)) {
    # Use same shift as used in plottable_data()
    min_x <- min(df_raw$x, na.rm = TRUE)
    min_y <- min(df_raw$y, na.rm = TRUE)
    min_z <- min(df_raw$z, na.rm = TRUE)
    
    df_raw <- df_raw %>%
      mutate(
        x = x - min_x,
        y = y - min_y,
        z = z - min_z
      )
  }
  
  # Use consistent axis limits for visual stability
  limits <- get_equal_limits(df$x, df$y, df$z)
  
  
  df$color <- as.character(df$color)
  
  # Axis theming in 3D
  bg_col   <- if (isTRUE(input$dark_mode)) "#000000" else "#FFFFFF"
  fg_col   <- if (isTRUE(input$dark_mode)) "#FFFFFF" else "#63666A"
  show_g   <- !isTRUE(input$hide_grid)
  grid_col <- if (isTRUE(input$dark_mode)) "rgba(255,255,255,0.18)" else "rgba(0,0,0,0.18)"
  zero_col <- if (isTRUE(input$dark_mode)) "rgba(255,255,255,0.25)" else "rgba(0,0,0,0.25)"
  
  # Build the point table we actually plot
  disp <- smlm_display_from_input(input)
  pts <- build_pts_3d(df, disp, deblinked = isTRUE(input$enable_clustering))
  
  # Compute colors for that table (default uses your existing palette in pts$color)
  pts$color <- build_colors_for_pts(pts, disp, default_palette_col = "color")
  
  cmap <- color_mapping_for_pts(pts, disp)
  
  # Plot the points ONCE
  fig <- plotly::plot_ly()
  
  id_col <- if (isTRUE(input$enable_clustering)) "cluster" else "tid"
  
  if (identical(id_col, "cluster")) {
    if (!("cluster" %in% names(pts)) || all(is.na(pts$cluster) | pts$cluster == 0)) {
      id_col <- "tid"
    }
  }
  
  # Key vector for points (use character keys ALWAYS)
  gid_key_pts <- as.character(pts[[id_col]])
  marker_hex  <- as.character(cmap$hex)
  # drop cluster 0
  keep <- !is.na(gid_key_pts) & nzchar(gid_key_pts)
  if (identical(id_col, "cluster")) keep <- keep & gid_key_pts != "0"
  
  gid_key_pts <- gid_key_pts[keep]
  marker_hex  <- marker_hex[keep]
  
  gid_hex_map <- tapply(marker_hex, gid_key_pts, function(x) x[1], simplify = TRUE)
  gid_hex_map <- unlist(gid_hex_map, use.names = TRUE)
  
  group_ids <- names(gid_hex_map)
  if (!length(group_ids)) return(fig)
  
  

  
  mplex_draw <- identical(input$color_mode %||% "default", "multiplex") &&
                "mplex_label" %in% names(pts) && any(!is.na(pts$mplex_label))
  if (mplex_draw) {
    # one trace per channel so each gets a named, colored legend entry
    ptso  <- pts[order(pts$mplex_channel), , drop = FALSE]
    chans <- unique(ptso$mplex_channel[!is.na(ptso$mplex_channel)])
    for (chv in chans) {
      sub <- ptso[!is.na(ptso$mplex_channel) & ptso$mplex_channel == chv, , drop = FALSE]
      if (!nrow(sub)) next
      fig <- fig %>% plotly::add_trace(
        data = sub, x = ~x, y = ~y, z = ~z,
        type = "scatter3d", mode = "markers",
        name = sub$mplex_label[1],
        marker = list(size = input$median_point_size %||% 4, color = sub$mplex_color[1]),
        text = ~hover_text,
        hovertemplate = "%{text}<extra></extra>",
        hoverinfo = "text",
        showlegend = TRUE, inherit = FALSE
      )
    }
  } else if (isTRUE(cmap$showscale) && identical(cmap$kind, "scale")) {
    fig <- fig %>% plotly::add_trace(
      data = pts,
      x = ~x, y = ~y, z = ~z,
      type = "scatter3d", mode = "markers",
      name = if (isTRUE(input$show_medians)) "Localization Medians" else "Localizations",
      marker = list(
        size = input$median_point_size %||% 4,
        color = cmap$colorvar,
        colorscale = cmap$colorscale,
        cmin = cmap$cmin, cmax = cmap$cmax,
        showscale = TRUE,
        colorbar = list(
          title     = cmap$title,
          titlefont = list(color = fg_col),
          tickfont  = list(color = fg_col)
        )
      ),
      text = ~hover_text,
      hovertemplate = "%{text}<extra></extra>",
      hoverinfo = "text",
      showlegend = TRUE, inherit = FALSE
    )
  } else {
    # hex colors (default / unicolor / degenerate numeric)
    hex <- if (!is.null(cmap$hex)) cmap$hex else rep("#bdbdbd", nrow(pts))
    fig <- fig %>% plotly::add_trace(
      data = pts,
      x = ~x, y = ~y, z = ~z,
      type = "scatter3d", mode = "markers",
      marker = list(size = input$median_point_size %||% 4, color = hex),
      text = ~hover_text,
      hovertemplate = "%{text}<extra></extra>",
      hoverinfo = "text",
      name = if (isTRUE(input$show_medians)) "Localization Medians" else "Localizations",
      showlegend = TRUE, inherit = FALSE
    )
  }
  fig <- fig %>% plotly::layout(
    title = "3D Scatter Plot",
    scene = list(
      xaxis = list(title = "X", range = limits$x),
      yaxis = list(title = "Y", range = limits$y),
      zaxis = list(title = "Z", range = limits$z),
      aspectmode = "manual",
      aspectratio = list(x = 1, y = 1, z = 1)
    ),
    showlegend = TRUE,
    legend = list(
      font = list(color = fg_col),
      orientation = "h",
      x = 0, y = 0.97,
      xanchor = "left", yanchor = "bottom"
    )
  )

  # Low-opacity preview of the 3D crop region (box / sphere) from the live inputs.
  fig <- add_3d_crop_preview(fig)

  cluster2_medians <- NULL
  if (isTRUE(input$enable_cluster2) && "cluster2" %in% names(df)) {
    # Member unit follows de-blinking, matching the two-colour isosurface and
    # every other per-structure statistic (cluster2_member_col / smlm_nn_unit).
    # This used to be hardcoded to `tid`, so with de-blinking on the envelope
    # was built from blink events rather than molecules.
    iso_unit <- cluster2_member_col(df)
    cluster2_medians <- df %>%
      dplyr::filter(!is.na(cluster2), cluster2 != 0) %>%
      { if (identical(iso_unit, "cluster"))
          dplyr::filter(., !is.na(cluster), cluster > 0L) else . } %>%
      # Grouped by (unit, channel), not the unit alone. A de-blinked molecule
      # never spans channels, but a trace ID can in multiplex mode (that split
      # is per localization), and one median over both would sit between them.
      dplyr::group_by(.data[[iso_unit]],
                      .ch = if ("channel" %in% names(.)) channel else NA_integer_) %>%
      dplyr::summarise(
        x = stats::median(x, na.rm = TRUE),
        y = stats::median(y, na.rm = TRUE),
        z = stats::median(z, na.rm = TRUE),
        cluster2 = dplyr::first(cluster2),
        # carried so the wrap can be keyed and coloured by channel (see
        # wrap_groups()); constant within a group by construction
        channel_label = if ("channel_label" %in% names(df))
                          dplyr::first(channel_label) else NA_character_,
        channel_color = if ("channel_color" %in% names(df))
                          dplyr::first(channel_color) else NA_character_,
        .groups = "drop"
      ) %>%
      dplyr::rename(channel = .ch) %>%
      dplyr::filter(is.finite(x), is.finite(y), is.finite(z), is.finite(cluster2), cluster2 != 0)
  }
  
  
  # PRIMARY ALPHA SHAPE RENDERING
  if (isTRUE(input$show_convex_hull) && nrow(df) > 3) {
    
    id_col <- if (isTRUE(input$enable_clustering)) "cluster" else "tid"
    if (!id_col %in% names(pts)) return(fig)   # since we now hull-wrap pts
    
    # group_ids + gid_hex_map must already exist (built from pts + cmap above)
    if (!exists("group_ids") || !length(group_ids) ||
        !exists("gid_hex_map") || !length(gid_hex_map)) {
      return(fig)
    }
    
    # No alpha radius is asked for: alpha_surface_3d() derives one per group
    # from that group's own point spacing. See its header in R/smlm_shapes.R for
    # why a single global radius was both meaningless (a 30 nm molecule and a
    # 500 nm structure need different radii) and inert (the old code drew the
    # whole tetrahedralization, which is identical at every alpha).
    for (gid_key in as.character(group_ids)) {

      subset_df <- pts[as.character(pts[[id_col]]) == gid_key, , drop = FALSE]
      if (nrow(subset_df) < 4) next

      srf <- alpha_surface_3d(as.matrix(subset_df[, c("x","y","z"), drop = FALSE]))
      if (is.null(srf)) next

      # name-based lookup (NOW guaranteed)
      hull_col <- unname(gid_hex_map[[as.character(gid_key)]])

      if (is.null(hull_col) || !nzchar(hull_col)) hull_col <- "#bdbdbd"

      fig <- fig %>% plotly::add_trace(
        type = "mesh3d",
        # x/y/z come from srf, NOT from subset_df: alpha_surface_3d may have
        # subsampled, and i/j/k index the points it returned.
        x = srf$x, y = srf$y, z = srf$z,
        i = srf$i, j = srf$j, k = srf$k,
        opacity = input$hull_opacity %||% 0.2,

        # Force solid mesh color (Plotly sometimes ignores `color` for mesh3d)
        color       = I(hull_col),
        facecolor   = rep(hull_col, srf$n_tri),         # per-triangle color
        vertexcolor = rep(hull_col, srf$n_pts),         # per-vertex color

        flatshading = TRUE,
        lighting = list(ambient = 1, diffuse = 0, specular = 0, roughness = 1, fresnel = 0),
        name = paste("Hull", gid_key),
        hoverinfo = "skip",
        showlegend = FALSE,
        inherit = FALSE
      )
    }
  }
  
  
  
  
  # --- 4B) Primary NN lines on plot (uses same k/range as the download table) ---
  # General NN
  # --- General NN (trace medians OR cluster medians) ---
  # 3D NN lines. This used to be a second, inline copy of build_nn_segments() -
  # the same geometry, pair selection, dedupe and range filtering written out
  # again - which is how the 2D and 3D views could disagree about which pairs
  # were neighbours. One builder now serves both.
  if (isTRUE(input$show_nn_lines) && nrow(df) > 1) {
    nn <- build_nn_segments(df, smlm_nn_options(input), is3d = TRUE)
    if (!is.null(nn)) {
      ln_col <- if (isTRUE(input$dark_mode)) "white" else "black"
      ed <- nn$edges
      multi_ch <- !is.null(ed) && nrow(ed) &&
        length(unique(stats::na.omit(c(ed$ch1, ed$ch2)))) > 1L
      if (multi_ch) {
        # Line colour says which channels the pair spans: within-channel takes
        # that channel's colour, cross-channel is magenta.
        for (cc in unique(ed$color)) {
          sel <- which(ed$color == cc)
          fig <- fig %>% plotly::add_trace(
            x = as.numeric(t(cbind(ed$x1[sel], ed$x2[sel], NA))),
            y = as.numeric(t(cbind(ed$y1[sel], ed$y2[sel], NA))),
            z = as.numeric(t(cbind(ed$z1[sel], ed$z2[sel], NA))),
            type = "scatter3d", mode = "lines",
            line = list(color = cc, width = 3),
            hoverinfo = "none", showlegend = FALSE, inherit = FALSE
          )
        }
      } else {
        fig <- fig %>% plotly::add_trace(
          x = nn$seg_x, y = nn$seg_y, z = nn$seg_z,
          type = "scatter3d", mode = "lines",
          line = list(color = ln_col, width = 2),
          hoverinfo = "none", showlegend = FALSE, inherit = FALSE
        )
      }
      if (length(nn$labs) && isTRUE(input$show_nn_labels)) {
        lab_col <- if (isTRUE(input$dark_mode)) "orange" else "red"
        nn_lfs <- suppressWarnings(as.numeric(input$nn_label_font_size %||% 12))
        if (!is.finite(nn_lfs) || nn_lfs < 1) nn_lfs <- 12
        fig <- fig %>% plotly::add_trace(
          x = nn$lab_x, y = nn$lab_y, z = nn$lab_z,
          type = "scatter3d", mode = "text",
          text = nn$labs, textfont = list(color = lab_col, size = nn_lfs),
          hoverinfo = "none", showlegend = FALSE, inherit = FALSE
        )
      }
    }
  }
  # ── 2-color alpha shapes (de-blinked molecules) ─────────────────────────
  # Only the per-molecule wrap lives here. The Median Cluster wrap and
  # isosurface are drawn by the Median Clustering controls, per channel (see
  # wrap_groups()); a second copy here was redundant and had drifted.
  if (isTRUE(input$enable_two_color) && isTRUE(input$show_tc1_alpha) &&
      all(c("x","y","z") %in% names(df))) {
    
    add_alpha_for <- function(cluster_col, alpha, opacity, min_pts) {
      ids <- sort(unique(df[[cluster_col]]))
      ids <- ids[is.finite(ids) & ids != 0]
      if (!length(ids)) return(fig)
      
      for (cid in ids) {
        sub <- df[df[[cluster_col]] == cid, , drop = FALSE]
        if (nrow(sub) < min_pts) next
        
        pts <- as.matrix(sub[, c("x","y","z")])
        pts <- pts[stats::complete.cases(pts), , drop = FALSE]
        if (nrow(pts) < min_pts) next
        
        msh <- make_alpha_mesh3d(pts, alpha = alpha)
        if (is.null(msh)) next
        
        fig <- fig %>% plotly::add_trace(
          type = "mesh3d",
          x = msh$x, y = msh$y, z = msh$z,
          i = msh$i, j = msh$j, k = msh$k,
          opacity = opacity,
          hoverinfo = "skip",
          showscale = FALSE,
          inherit = FALSE
        )
      }
      fig
    }
    
    if (isTRUE(input$show_tc1_alpha) && "cluster" %in% names(df)) {
      fig <- add_alpha_for(
        cluster_col = "cluster",
        alpha   = input$tc1_alpha %||% 120,
        opacity = input$tc1_alpha_opacity %||% 0.15,
        min_pts = input$tc1_alpha_min_pts %||% 30
      )
    }
    
  }
  if (!isTRUE(input$show_cluster2_isosurface)) cluster2_iso_diag(NULL)

  # ── cluster2 isosurface (3D) ────────────────────────────────────────────
  if (isTRUE(input$show_cluster2_isosurface) &&
      isTRUE(input$enable_cluster2) &&
      !is.null(cluster2_medians) &&
      nrow(cluster2_medians) > 0 &&
      isTRUE(input$mode3d == "3D")) {
    
    # Geometry from the committed snapshot (Plot Data); opacity stays live.
    ig       <- iso_geom()
    grid_res <- ig$c2_grid_res %||% 22
    pad      <- ig$c2_pad %||% 0
    op       <- input$cluster2_iso_opacity %||% 0.25
    # NULL = auto, i.e. 3 x localization precision floored at half the member
    # gap. Anchoring to the microscope rather than to the labelling is what
    # keeps two clusters comparable; see cluster_density_field().
    iso_off  <- if (!isTRUE(ig$c2_manual_off)) NULL else {
      o <- suppressWarnings(as.numeric(ig$c2_offset %||% 30))
      if (!is.finite(o) || o <= 0) NULL else o
    }
    sigma    <- if (!isTRUE(ig$c2_manual_sigma)) NULL else {
      s <- suppressWarnings(as.numeric(ig$c2_sigma %||% 30))
      if (!is.finite(s) || s <= 0) NULL else s
    }

    # See the tc2 block: one dataset-wide precision, so every Median Cluster in
    # the view is wrapped on the same physical ruler.
    # Per-axis: MINFLUX z is usually TIGHTER than x/y (measured z:lateral 0.60
    # on the a real imaging file files), so an isotropic wrap over-wraps axially.
    iso_precision <- smlm_loc_precision_nm(df, per_axis = TRUE)
    c2_diag <- new.env(parent = emptyenv()); c2_diag$rows <- list()

    # One wrap per Median Cluster PER CHANNEL, so a surface can never enclose
    # molecules of two colours. See wrap_groups().
    groups <- wrap_groups(cluster2_medians, "cluster2")
    # Decided on the WHOLE dataset, before any single cluster is picked out:
    # a selected channel-B surface must still say "(B)" and still be red.
    # Counted over all points, the same test that colors the points. Counting
    # the clusters instead gets it wrong when only one channel formed any.
    multi_ch <- n_frame_channels(plottable_data()) > 1L

    if (isTRUE(input$cluster2_iso_selected_only)) {
      sel <- .int_input(input$cluster2_iso_id, 1L, 1L)
      groups <- groups[groups$id == sel, , drop = FALSE]
    }

    if (nrow(groups)) {

    # build per-cluster color map using the same palette as cluster2 points
    qual_pals_iso <- RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == "qual", ]
    base_colors_iso <- unlist(mapply(RColorBrewer::brewer.pal, qual_pals_iso$maxcolors, rownames(qual_pals_iso)))
    pal_iso <- rep(base_colors_iso, length.out = max(groups$id) + 1)

    for (gi in seq_len(nrow(groups))) {
      g   <- groups[gi, , drop = FALSE]
      cid <- g$id
      sub <- wrap_group_rows(cluster2_medians, "cluster2", g)
      if (nrow(sub) < 4) next

      pts <- as.matrix(sub[, c("x","y","z")])
      pts <- pts[stats::complete.cases(pts), , drop = FALSE]
      if (nrow(pts) < 4) next

      # Gaussian KDE over the member centroids, contoured at the level that
      # sits `offset` nm outside them. See cluster_density_field() for why the
      # field is a density rather than a nearest-distance transform, and why
      # the offset is anchored to localization precision rather than spacing.
      # Colour precedence: an explicit single colour the user picked, then the
      # channel the members belong to, then the cluster palette. Channel colour
      # is what makes this usable in multiplex mode, where the cluster palette
      # says nothing about which target a structure is.
      # By-channel colouring only means something with more than one channel;
      # with one it would paint every surface the same colour.
      iso_col <- if (isTRUE(input$cluster2_iso_unicolor)) {
        input$cluster2_iso_color %||% "#00CED1"
      } else if (multi_ch && isTRUE(input$cluster2_iso_by_channel)) {
        cc <- iso_channel_color(sub)
        if (is.na(cc)) pal_iso[(cid %% length(pal_iso)) + 1] else cc
      } else {
        pal_iso[(cid %% length(pal_iso)) + 1]
      }

      fig <- add_cluster_isosurface(
        fig, pts, color = iso_col, opacity = op,
        name = if (multi_ch && nzchar(g$label %||% ""))
                 paste0("Median Cluster ", cid, " (", g$label, ")")
               else paste0("Median Cluster ", cid),
        grid_res = grid_res, pad = pad, sigma = sigma, offset = iso_off,
        precision = iso_precision, diag_env = c2_diag,
        sigma_scale = ig$c2_smooth_rel %||% 1, cache = iso_mesh_cache
      )
    }
    } # end else(!length(ids))
    cluster2_iso_diag(c2_diag$rows)
  }


  if (isTRUE(input$show_xy_plane)) {
    df_raw <- isolate(data())
    df_full <- plottable_data()  # or df_full <- data()
    limits <- get_equal_limits(df_full$x, df_full$y, df_full$z)
    df_for_plane <- if (isTRUE(input$plane_from_unclustered)) {
      df_full %>% filter(is.na(cluster2) | cluster2 == 0)
    } else {
      df_full
    }
    
    plane <- generate_xy_plane(df_for_plane, limits, z_offset = input$plane_offset)
    fig <- fig %>% add_surface(
      x = plane$x,
      y = plane$y,
      z = plane$z,
      opacity = 0.3,
      showscale = FALSE,
      colorscale = list(list(0, "gray"), list(1, "gray")),
      name = "XY Tilted Plane"
    )
  }
  if (isTRUE(input$enable_cluster2) && isTRUE(input$show_cluster2_metrics)) {
    df_full <- plottable_data()
    limits <- get_equal_limits(df_full$x, df_full$y, df_full$z)
    df_for_plane <- if (isTRUE(input$plane_from_unclustered)) {
      df_full %>% filter(is.na(cluster2) | cluster2 == 0)
    } else {
      df_full
    }
    
    plane <- generate_xy_plane(df_for_plane, limits, z_offset = input$plane_offset)
    
    membrane_mode <- isTRUE(input$membrane_protein_options)
    
    result <- add_cluster2_lines_and_angles(
      fig,
      df_full,
      plane = plane,
      use_pca_normal = input$show_xy_plane,
      cluster2_table_store = cluster2_table_store,
      enable_clustering = input$enable_clustering,
      helix_turns = input$helix_turns,
      dark_mode = input$dark_mode,
      bending_radius      = membrane_mode && isTRUE(input$show_bending_radius),
      show_plane_angles   = membrane_mode && isTRUE(input$show_cluster2_plane_angle),
      show_plane_distance = membrane_mode && isTRUE(input$show_cluster2_plane_distance),
      show_pca_axes = input$show_cluster2_pca_axes,
      show_info_labels = input$show_cluster2_labels,
      selected_fit_type   = input$fit_type,
      selected_fit_method = input$fit_method,
      compute_rmsd_metrics = isTRUE(input$compute_rmsd_metrics),
      fit_weighted        = isTRUE(input$fit_weighted),
      mle_center_mode     = input$mle_center_mode,
      mle_center_lambda   = input$mle_center_lambda,
      arc_span_hint       = input$arc_span_deg %||% 0
    )
    
    fig <- result$fig
    fit_log_rv(result$log)
  }
  # ---- GLOBAL SURFACE (from representative points) ----
  if (isTRUE(input$view_surface)) {
    df_surface <- plottable_data()
    
    # 1) optionally use only unclustered (must happen on raw points BEFORE aggregation)
    if (isTRUE(input$plane_from_unclustered)) {
      df_surface <- df_surface %>% dplyr::filter(is.na(cluster2) | cluster2 == 0)
      # If you meant "unclustered in cluster1", use:
      # df_surface <- df_surface %>% dplyr::filter(is.na(cluster) | cluster == 0)
    }
    
    # 2) exclude TIDs (still on raw points, before aggregation)
    ex_txt <- input$surface_exclude_tids %||% ""
    ex_ids <- unique(unlist(strsplit(ex_txt, "[,;\\s]+")))
    ex_ids <- ex_ids[nzchar(ex_ids)]
    
    if (length(ex_ids) && "tid" %in% names(df_surface)) {
      df_surface <- df_surface %>%
        dplyr::mutate(tid_chr = as.character(tid)) %>%
        dplyr::filter(!(.data$tid_chr %in% ex_ids))
    }

    # 2b) filter by DCR/fluo group (group 1 or 2 only)
    surface_group_sel <- input$surface_color_group %||% "all"
    if (surface_group_sel != "all") {
      # If plottable_data's analysis_group is all NA (e.g. grouping_mode="none")
      # but two-color is active, derive group tags via assign_two_color_group.
      has_valid_group <- "analysis_group" %in% names(df_surface) &&
                         any(!is.na(df_surface$analysis_group))
      if (!has_valid_group && isTRUE(input$enable_two_color)) {
        gmap <- tryCatch(assign_two_color_group(df_surface, as_smlm_settings(input)), error = function(e) NULL)
        if (is.data.frame(gmap) && "analysis_group" %in% names(gmap)) {
          df_surface <- df_surface %>%
            dplyr::select(-dplyr::any_of("analysis_group")) %>%
            dplyr::left_join(gmap %>% dplyr::select(tid, analysis_group), by = "tid")
        }
      }
      g_int <- suppressWarnings(as.integer(surface_group_sel))
      if ("analysis_group" %in% names(df_surface) && is.finite(g_int)) {
        df_surface <- df_surface %>%
          dplyr::filter(.data$analysis_group == g_int)
      }
    }

    # 3) build representative point set:
    #    - if de-blinking enabled -> de-blinked molecule medians
    #    - else -> tid medians
    if (isTRUE(input$enable_clustering) && "cluster" %in% names(df_surface)) {
      df_surface_pts <- df_surface %>%
        dplyr::filter(!is.na(cluster) & cluster != 0) %>%
        dplyr::group_by(cluster) %>%
        dplyr::summarise(
          x = stats::median(x, na.rm = TRUE),
          y = stats::median(y, na.rm = TRUE),
          z = stats::median(z, na.rm = TRUE),
          .groups = "drop"
        )
    } else {
      df_surface_pts <- df_surface %>%
        dplyr::group_by(tid) %>%
        dplyr::summarise(
          x = stats::median(x, na.rm = TRUE),
          y = stats::median(y, na.rm = TRUE),
          z = stats::median(z, na.rm = TRUE),
          .groups = "drop"
        )
    }
    
    df_surface_pts <- df_surface_pts %>%
      dplyr::mutate(across(c(x,y,z), as.numeric)) %>%
      dplyr::filter(is.finite(x), is.finite(y), is.finite(z))
    
    method <- input$surface_method %||% "blnk"
    if (method != "blnk" && nrow(df_surface_pts) >= 10) {
      
      if (method == "loess") {
        if (nrow(df_surface_pts) > 100000) df_surface_pts <- dplyr::slice_sample(df_surface_pts, n = 100000)
        
        surf <- tryCatch(
          create_global_surface(
            df_surface_pts,
            grid_n = input$surface_grid_n %||% 120,
            span   = input$surface_span   %||% 0.25
          ),
          error = function(e) NULL
        )
        
        if (!is.null(surf)) {
          fig <- fig %>% plotly::add_surface(
            x = surf$x, y = surf$y, z = surf$z,
            opacity = input$surface_opacity %||% 0.35,
            showscale = TRUE,
            colorscale = plotly_colorscale(
              input$surface_palette %||% "Viridis",
              reverse = isTRUE(input$surface_palette_reverse)
            ),
            colorbar = list(
              title     = "Surface z (nm)",
              titlefont = list(color = fg_col),
              tickfont  = list(color = fg_col)
            ),
            name = "Surface (LOESS)",
            hoverinfo = "skip",
            hovertemplate = ""
          )
        }
        
      } else if (method == "poly") {
        poly <- fit_poly_surface(
          df_surface_pts,
          degree = as.integer(input$poly_degree %||% 2L),
          grid_n  = as.integer(input$surface_grid_n %||% 120L)
        )
        if (!is.null(poly)) {
          fig <- fig %>% plotly::add_surface(
            x = poly$x, y = poly$y, z = poly$z,
            opacity = input$surface_opacity %||% 0.35,
            showscale = TRUE,
            colorscale = plotly_colorscale(
              input$surface_palette %||% "Viridis",
              reverse = isTRUE(input$surface_palette_reverse)
            ),
            colorbar = list(
              title     = "Surface z (nm)",
              titlefont = list(color = fg_col),
              tickfont  = list(color = fg_col)
            ),
            name = "Surface (Polynomial)",
            hoverinfo = "skip",
            hovertemplate = ""
          )
        }
        
      } else if (method == "tps") {
        if (nrow(df_surface_pts) > 50000) df_surface_pts <- dplyr::slice_sample(df_surface_pts, n = 50000)
        
        grid_n <- as.integer(input$surface_grid_n %||% 120L)
        gx <- seq(min(df_surface_pts$x), max(df_surface_pts$x), length.out = grid_n)
        gy <- seq(min(df_surface_pts$y), max(df_surface_pts$y), length.out = grid_n)
        
        gr <- expand.grid(x = gx, y = gy)
        
        fit2s <- fields::Tps(
          x = as.matrix(df_surface_pts[, c("x","y")]),
          Y = df_surface_pts$z,
          lambda = input$tps_lambda %||% 0.01
        )
        
        zhat <- predict(fit2s, as.matrix(gr[, c("x","y")]))
        

        zmat <- matrix(zhat, nrow = length(gx), ncol = length(gy), byrow = FALSE)
        zmat <- t(zmat)  # now rows = y, cols = x
        
        surf <- list(x = gx, y = gy, z = zmat)
        
        if (!is.null(surf)) {
          fig <- fig %>% plotly::add_surface(
            x = surf$x, y = surf$y, z = surf$z,
            opacity = input$surface_opacity %||% 0.35,
            showscale = TRUE,
            colorscale = plotly_colorscale(
              input$surface_palette %||% "Viridis",
              reverse = isTRUE(input$surface_palette_reverse)
            ),
            colorbar = list(
              title     = "Surface z (nm)",
              titlefont = list(color = fg_col),
              tickfont  = list(color = fg_col)
            ),
            name = "Surface (TPS)",
            hoverinfo = "skip",
            hovertemplate = ""
          )
        }
      }
      
    }
  }
  
  
  if (isTRUE(input$show_cluster2_nn_lines) && isTRUE(input$enable_cluster2) && !is.null(df$cluster2)) {
    lines <- secondary_nn_lines()
    if (length(lines)) {
      seg_x <- c(); seg_y <- c(); seg_z <- c()
      lab_x <- c(); lab_y <- c(); lab_z <- c(); labs <- c()
      for (ln in lines) {
        seg_x <- c(seg_x, ln$x, NA_real_)
        seg_y <- c(seg_y, ln$y, NA_real_)
        seg_z <- c(seg_z, ln$z, NA_real_)
        lab_x <- c(lab_x, ln$midpoint[1]); lab_y <- c(lab_y, ln$midpoint[2]); lab_z <- c(lab_z, ln$midpoint[3]); labs <- c(labs, ln$text)
      }
      fig <- fig %>% plotly::add_trace(
        x = seg_x, y = seg_y, z = seg_z,
        type = "scatter3d", mode = "lines",
        line = list(color = "green", width = 2),
        hoverinfo = "none", showlegend = FALSE, inherit = FALSE
      )
      if (length(labs) && isTRUE(input$show_c2_nn_labels)) {
        nn_lfs <- suppressWarnings(as.numeric(input$c2_label_font_size %||% 12))
        if (!is.finite(nn_lfs) || nn_lfs < 1) nn_lfs <- 12
        fig <- fig %>% plotly::add_trace(
          x = lab_x, y = lab_y, z = lab_z,
          type = "scatter3d", mode = "text",
          text = labs, textfont = list(color = "green", size = nn_lfs),
          hoverinfo = "none", showlegend = FALSE, inherit = FALSE
        )
      }
    }
  }
  # ==== Cluster2 fit drawing (Alpha or Smooth) =====================
  if (isTRUE(input$enable_cluster2) && !is.null(input$fit_type) && input$fit_type != "none") {
    df_full <- plottable_data()
    df_cluster2 <- df_full %>% dplyr::filter(!is.na(cluster2) & cluster2 != 0)
    clusters <- df_cluster2 %>% dplyr::group_by(cluster2) %>% dplyr::group_split()
    if (!length(clusters)) {
      fit_log_rv("No Median Cluster structures available for fitting. Check Median Clustering DBSCAN settings (eps2/minpts2).")
    }
    cap_radius_for_cluster <- function(radius, X, cap_mult = 2.5) {
      if (!is.finite(radius) || radius <= 0) return(NA_real_)
      ctr <- colMeans(X, na.rm = TRUE)
      d <- sqrt(rowSums((X - matrix(ctr, nrow(X), 3, byrow = TRUE))^2))
      d <- d[is.finite(d)]
      if (!length(d)) return(radius)
      cap <- as.numeric(stats::quantile(d, 0.95, na.rm = TRUE)) * cap_mult
      if (!is.finite(cap) || cap <= 0) cap <- max(d, na.rm = TRUE)
      min(radius, cap)
    }
    cap_axes_for_cluster <- function(axes, X, cap_mult = 2.5) {
      ax <- as.numeric(axes)
      if (length(ax) != 3) return(rep(NA_real_, 3))
      ctr <- colMeans(X, na.rm = TRUE)
      d <- sqrt(rowSums((X - matrix(ctr, nrow(X), 3, byrow = TRUE))^2))
      d <- d[is.finite(d)]
      if (!length(d)) return(ax)
      cap <- as.numeric(stats::quantile(d, 0.95, na.rm = TRUE)) * cap_mult
      if (!is.finite(cap) || cap <= 0) cap <- max(d, na.rm = TRUE)
      pmax(1e-6, pmin(ax, cap))
    }
    # GOAL: triangulate the surface of a FIT (sphere / ellipsoid), given the
    # dense uniform point sampling of it that the fit functions return.
    # Those surfaces are convex by construction, so the convex hull IS the
    # exact triangulation - no alpha shape is needed or wanted here. (This
    # used to build an ashape3d first and select `triang[, 4] >= 2`; column 4
    # is `on.ch`, a 0/1 flag, so that test never passed and every call fell
    # through to this hull anyway, after paying for a 4000-point
    # tetrahedralization. Data clouds, which are NOT convex, are wrapped by
    # alpha_surface_3d() instead.)
    surface_triangles <- function(pts) {
      tryCatch({
        hh <- geometry::convhulln(pts, options = "Qt")
        as.matrix(hh[, 1:3, drop = FALSE])
      }, error = function(e) NULL)
    }
    build_tid_fit_medians <- function(d) {
      d %>%
        dplyr::group_by(.data$tid) %>%
        dplyr::summarise(
          x = stats::median(.data$x, na.rm = TRUE),
          y = stats::median(.data$y, na.rm = TRUE),
          z = stats::median(.data$z, na.rm = TRUE),
          spread = stats::mad(
            sqrt(
              (.data$x - stats::median(.data$x, na.rm = TRUE))^2 +
              (.data$y - stats::median(.data$y, na.rm = TRUE))^2 +
              (.data$z - stats::median(.data$z, na.rm = TRUE))^2
            ),
            na.rm = TRUE
          ),
          n_locs = dplyr::n(),
          .groups = "drop"
        ) %>%
        dplyr::filter(is.finite(.data$x), is.finite(.data$y), is.finite(.data$z))
    }
    
    # Cache MLE fit objects keyed by cluster_id so the distance-lines loop
    # below can reuse them without refitting the same data a second time.
    mle_fit_cache <- list()

    for (cluster_data in clusters) {
      cluster_id <- unique(cluster_data$cluster2)
      if (is.na(cluster_id) || nrow(cluster_data) < 4) next

      # --- choose color for this cluster’s fits ---
      fit_col <- {
        co <- unique(na.omit(as.character(cluster_data$cluster2_color)))
        if (length(co) && nzchar(co[1])) co[1] else "#00BFFF"   # default fallback
      }
      
      # medians used by MLE path weights, unchanged
      if (isTRUE(input$enable_clustering) && "cluster" %in% names(cluster_data)) {
        med <- cluster_data %>%
          dplyr::filter(!is.na(.data$cluster) & .data$cluster != 0) %>%
          dplyr::group_by(.data$cluster) %>%
          dplyr::summarise(
            mx     = stats::median(.data$x, na.rm = TRUE),
            my     = stats::median(.data$y, na.rm = TRUE),
            mz     = stats::median(.data$z, na.rm = TRUE),
            spread = if (dplyr::n() > 1)
              sqrt(mean((.data$x - mx)^2 + (.data$y - my)^2 + (.data$z - mz)^2, na.rm = TRUE)) else 0,
            n_locs = dplyr::n(),
            .groups = "drop"
          ) %>%
          dplyr::rename(x = mx, y = my, z = mz)
      } else {
        med <- cluster_data %>%
          dplyr::group_by(.data$tid) %>%
          dplyr::summarise(
            mx     = stats::median(.data$x, na.rm = TRUE),
            my     = stats::median(.data$y, na.rm = TRUE),
            mz     = stats::median(.data$z, na.rm = TRUE),
            spread = if (dplyr::n() > 1)
              sqrt(mean((.data$x - mx)^2 + (.data$y - my)^2 + (.data$z - mz)^2, na.rm = TRUE)) else 0,
            n_locs = dplyr::n(),
            .groups = "drop"
          ) %>%
          dplyr::rename(x = mx, y = my, z = mz)
      }

      min_fit_n <- if (identical(input$fit_type, "helix")) 5L else 4L
      if (nrow(med) < min_fit_n) {
        med <- build_tid_fit_medians(cluster_data)
      }
      if (nrow(med) < min_fit_n) next
      X <- as.matrix(med[, c("x","y","z")])
      w <- if (isTRUE(input$fit_weighted)) {
        s   <- med$spread
        n_l <- if ("n_locs" %in% names(med)) med$n_locs else rep(1L, nrow(med))
        s_ref <- stats::median(s[n_l > 1 & s > 0], na.rm = TRUE)
        if (!is.finite(s_ref) || s_ref <= 0) s_ref <- max(s[s > 0], na.rm = TRUE)
        if (!is.finite(s_ref) || s_ref <= 0) s_ref <- 1
        s_eff <- ifelse(n_l > 1 & s > 0, s, s_ref)
        n_l / s_eff^2
      } else NULL
      
      if (input$fit_method == "heuristic") {
        # --- Heuristic path ---
        if (input$fit_type == "sphere") {
          fit <- generate_sphere_from_cluster2(
            cluster_data,
            use_cluster1 = isTRUE(input$enable_clustering)
          )
          if (!is.null(fit)) {
            fit$radius <- cap_radius_for_cluster(fit$radius, X)
            if (!is.finite(fit$radius) || fit$radius <= 0) next
            pts <- fit$virtual_points
            surf_tris <- surface_triangles(pts)
            if (!is.null(surf_tris) && nrow(surf_tris) > 0) {
              fig <- fig %>% add_trace(
                type = "mesh3d",
                x = pts[,1], y = pts[,2], z = pts[,3],
                i = surf_tris[,1]-1, j = surf_tris[,2]-1, k = surf_tris[,3]-1,
                opacity = input$fit_opacity,
                color = I(fit_col),
                facecolor   = rep(fit_col, nrow(surf_tris)),
                vertexcolor = rep(fit_col, nrow(pts)),
                flatshading = TRUE,
                lighting = list(ambient = 1, diffuse = 0, specular = 0, roughness = 1, fresnel = 0),
                name = paste("Spherical Fit", cluster_id),
                hoverinfo = "skip",
                hovertemplate = "",
                showlegend = TRUE
              )
            }
          }

        } else if (input$fit_type == "ellipsoid") {
          fit <- generate_ellipsoid_from_cluster2(
            cluster_data,
            use_cluster1 = isTRUE(input$enable_clustering)
          )
          if (!is.null(fit)) {
            fit$axes <- cap_axes_for_cluster(fit$axes, X)
            if (!all(is.finite(fit$axes)) || any(fit$axes <= 0)) next
            pts <- fit$virtual_points
            surf_tris <- surface_triangles(pts)
            if (!is.null(surf_tris) && nrow(surf_tris) > 0) {
              fig <- fig %>% add_trace(
                type = "mesh3d",
                x = pts[,1], y = pts[,2], z = pts[,3],
                i = surf_tris[,1]-1, j = surf_tris[,2]-1, k = surf_tris[,3]-1,
                opacity = input$fit_opacity,
                color = I(fit_col),
                facecolor   = rep(fit_col, nrow(surf_tris)),
                vertexcolor = rep(fit_col, nrow(pts)),
                flatshading = TRUE,
                lighting = list(ambient = 1, diffuse = 0, specular = 0, roughness = 1, fresnel = 0),
                name = paste("Ellipsoid Fit", cluster_id),
                hoverinfo = "skip",
                hovertemplate = "",
                showlegend = TRUE
              )
            }
          }
          
        } else if (input$fit_type == "helix") {
          fit <- generate_helix_from_cluster2(
            cluster_data,
            n_points = 3000,
            turns = input$helix_turns,
            use_cluster1 = isTRUE(input$enable_clustering)
          )
          if (!is.null(fit)) {
            if (isTRUE(input$show_fit_tube)) {
              geo_h  <- fit_helix_geometry(X, turns = input$helix_turns)
              tube_r <- max(1, input$fit_tube_radius %||% 5)
              tmax_h <- 2 * pi * max(1, input$helix_turns)
              msh <- tryCatch(make_helix_tube_mesh(
                center = geo_h$proj_center, axis = geo_h$axis,
                radius = geo_h$radius, pitch = geo_h$pitch, phi0 = geo_h$phi0,
                tmin = 0, tmax = tmax_h, tube_r = tube_r
              ), error = function(e) NULL)
              if (!is.null(msh)) {
                fig <- fig %>% add_trace(
                  type = "mesh3d", x = msh$x, y = msh$y, z = msh$z,
                  i = msh$i, j = msh$j, k = msh$k,
                  opacity = input$fit_opacity, color = I(fit_col),
                  facecolor = rep(fit_col, length(msh$i)),
                  vertexcolor = rep(fit_col, length(msh$x)),
                  flatshading = TRUE,
                  lighting = list(ambient = 1, diffuse = 0, specular = 0, roughness = 1, fresnel = 0),
                  name = paste("Helix Fit", cluster_id),
                  showlegend = TRUE, inherit = FALSE,
                  hoverinfo = "skip", hovertemplate = ""
                )
              }
            } else {
              hp <- fit$points
              fig <- fig %>% add_trace(
                x = hp$x, y = hp$y, z = hp$z,
                type = "scatter3d", mode = "lines",
                line = list(color = fit_col, width = 4),
                name = paste("Helix Fit", cluster_id),
                showlegend = TRUE, inherit = FALSE,
                hoverinfo = "skip", hovertemplate = ""
              )
            }
          }

        } else if (input$fit_type == "ring") {
          fit <- generate_ring_from_cluster2(
            cluster_data,
            n_points = 500,
            use_cluster1 = isTRUE(input$enable_clustering)
          )
          if (!is.null(fit)) {
            if (isTRUE(input$show_fit_tube)) {
              tube_r <- max(1, input$fit_tube_radius %||% 5)
              msh <- tryCatch(make_ring_tube_mesh(
                fit$center, fit$normal, fit$radius, tube_r
              ), error = function(e) NULL)
              if (!is.null(msh)) {
                fig <- fig %>% add_trace(
                  type = "mesh3d", x = msh$x, y = msh$y, z = msh$z,
                  i = msh$i, j = msh$j, k = msh$k,
                  opacity = input$fit_opacity, color = I(fit_col),
                  facecolor = rep(fit_col, length(msh$i)),
                  vertexcolor = rep(fit_col, length(msh$x)),
                  flatshading = TRUE,
                  lighting = list(ambient = 1, diffuse = 0, specular = 0, roughness = 1, fresnel = 0),
                  name = paste("Ring Fit", cluster_id),
                  showlegend = TRUE, inherit = FALSE,
                  hoverinfo = "skip", hovertemplate = ""
                )
              }
            } else {
              rp <- fit$points
              fig <- fig %>% add_trace(
                x = rp$x, y = rp$y, z = rp$z,
                type = "scatter3d", mode = "lines",
                line = list(color = fit_col, width = 4),
                name = paste("Ring Fit", cluster_id),
                showlegend = TRUE, inherit = FALSE,
                hoverinfo = "skip", hovertemplate = ""
              )
            }
          }

        } else if (input$fit_type == "arc") {
          fit <- generate_arc_from_cluster2(
            cluster_data,
            n_points = 300,
            arc_span_hint = input$arc_span_deg %||% 0,
            use_cluster1 = isTRUE(input$enable_clustering)
          )
          if (!is.null(fit)) {
            if (isTRUE(input$show_fit_tube)) {
              tube_r <- max(1, input$fit_tube_radius %||% 5)
              msh <- tryCatch(make_arc_tube_mesh(
                fit$center, fit$normal, fit$radius,
                fit$t_center, fit$arc_half, tube_r
              ), error = function(e) NULL)
              if (!is.null(msh)) {
                fig <- fig %>% add_trace(
                  type = "mesh3d", x = msh$x, y = msh$y, z = msh$z,
                  i = msh$i, j = msh$j, k = msh$k,
                  opacity = input$fit_opacity, color = I(fit_col),
                  facecolor = rep(fit_col, length(msh$i)),
                  vertexcolor = rep(fit_col, length(msh$x)),
                  flatshading = TRUE,
                  lighting = list(ambient = 1, diffuse = 0, specular = 0, roughness = 1, fresnel = 0),
                  name = paste("Arc Fit", cluster_id),
                  showlegend = TRUE, inherit = FALSE,
                  hoverinfo = "skip", hovertemplate = ""
                )
              }
            } else {
              rp <- fit$points
              fig <- fig %>% add_trace(
                x = rp$x, y = rp$y, z = rp$z,
                type = "scatter3d", mode = "lines",
                line = list(color = fit_col, width = 4),
                name = paste("Arc Fit", cluster_id),
                showlegend = TRUE, inherit = FALSE,
                hoverinfo = "skip", hovertemplate = ""
              )
            }
          }
        }

      } else {
        # --------- MLE path with center handling ---------
        center_args <- build_center_args(X, input$mle_center_mode, input$mle_center_lambda)
        
        if (input$fit_type == "sphere") {
          f <- do.call(mle_fit_sphere, c(list(X = X, w = w), center_args))
          mle_fit_cache[[as.character(cluster_id)]] <- f   # cache for distance-lines loop
          if (!is.null(f)) {
            f$radius <- cap_radius_for_cluster(f$radius, X)
            if (!is.finite(f$radius) || f$radius <= 0) next
            pts <- mle_sample_sphere(f$center, f$radius, n = 4000)
            surf_tris <- surface_triangles(pts)
            if (!is.null(surf_tris) && nrow(surf_tris) > 0) {
              fig <- fig %>% add_trace(
                type = "mesh3d",
                x = pts[,1], y = pts[,2], z = pts[,3],
                i = surf_tris[,1]-1, j = surf_tris[,2]-1, k = surf_tris[,3]-1,
                opacity = input$fit_opacity,
                color   = I(fit_col),
                facecolor   = rep(fit_col, nrow(surf_tris)),
                vertexcolor = rep(fit_col, nrow(pts)),
                flatshading = TRUE,
                lighting = list(ambient = 1, diffuse = 0, specular = 0, roughness = 1, fresnel = 0),
                name = paste("Sphere (MLE)", cluster_id),
                hoverinfo = "skip",
                hovertemplate = "",
                showlegend = TRUE
              )
            }
          }

        } else if (input$fit_type == "ellipsoid") {
          f <- do.call(mle_fit_ellipsoid, c(list(X = X, w = w), center_args))
          mle_fit_cache[[as.character(cluster_id)]] <- f   # cache for distance-lines loop
          if (!is.null(f)) {
            f$axes <- cap_axes_for_cluster(f$axes, X)
            if (!all(is.finite(f$axes)) || any(f$axes <= 0)) next
            pts <- mle_sample_ellipsoid(f$center, f$axes, f$R, n = 4000)
            surf_tris <- surface_triangles(pts)
            if (!is.null(surf_tris) && nrow(surf_tris) > 0) {
              fig <- fig %>% add_trace(
                type = "mesh3d",
                x = pts[,1], y = pts[,2], z = pts[,3],
                i = surf_tris[,1]-1, j = surf_tris[,2]-1, k = surf_tris[,3]-1,
                opacity = input$fit_opacity,
                color   = I(fit_col),
                facecolor   = rep(fit_col, nrow(surf_tris)),
                vertexcolor = rep(fit_col, nrow(pts)),
                flatshading = TRUE,
                lighting = list(ambient = 1, diffuse = 0, specular = 0, roughness = 1, fresnel = 0),
                name = paste("Ellipsoid (MLE)", cluster_id),
                showlegend = TRUE,
                hoverinfo = "skip",
                hovertemplate = ""
              )
            }
          }
          
        } else if (input$fit_type == "helix") {
          f <- do.call(mle_fit_helix, c(list(X = X, w = w, turns_hint = input$helix_turns), center_args))
          mle_fit_cache[[as.character(cluster_id)]] <- f
          if (!is.null(f)) {
            if (isTRUE(input$show_fit_tube)) {
              tube_r <- max(1, input$fit_tube_radius %||% 5)
              pad    <- 0.25
              tmin_h <- min(f$t) - 2*pi*pad
              tmax_h <- max(f$t) + 2*pi*pad
              msh <- tryCatch(make_helix_tube_mesh(
                center = f$center, axis = f$axis,
                radius = f$radius, pitch = f$pitch, phi0 = f$phi0,
                tmin = tmin_h, tmax = tmax_h, tube_r = tube_r
              ), error = function(e) NULL)
              if (!is.null(msh)) {
                fig <- fig %>% add_trace(
                  type = "mesh3d", x = msh$x, y = msh$y, z = msh$z,
                  i = msh$i, j = msh$j, k = msh$k,
                  opacity = input$fit_opacity, color = I(fit_col),
                  facecolor = rep(fit_col, length(msh$i)),
                  vertexcolor = rep(fit_col, length(msh$x)),
                  flatshading = TRUE,
                  lighting = list(ambient = 1, diffuse = 0, specular = 0, roughness = 1, fresnel = 0),
                  name = paste("Helix (MLE)", cluster_id), showlegend = TRUE, inherit = FALSE,
                  hoverinfo = "skip", hovertemplate = ""
                )
              }
            } else {
              hp <- mle_sample_helix(f, n = 3000)
              fig <- fig %>% add_trace(
                x = hp$V1, y = hp$V2, z = hp$V3,
                type = "scatter3d", mode = "lines",
                line = list(color = fit_col, width = 4),
                name = paste("Helix (MLE)", cluster_id), showlegend = TRUE, inherit = FALSE,
                hoverinfo = "skip", hovertemplate = ""
              )
            }
          }

        } else if (input$fit_type == "ring") {
          f <- do.call(mle_fit_ring, c(list(X = X, w = w), center_args))
          mle_fit_cache[[as.character(cluster_id)]] <- f
          if (!is.null(f)) {
            if (isTRUE(input$show_fit_tube)) {
              tube_r <- max(1, input$fit_tube_radius %||% 5)
              msh <- tryCatch(make_ring_tube_mesh(
                f$center, f$normal, f$radius, tube_r
              ), error = function(e) NULL)
              if (!is.null(msh)) {
                fig <- fig %>% add_trace(
                  type = "mesh3d", x = msh$x, y = msh$y, z = msh$z,
                  i = msh$i, j = msh$j, k = msh$k,
                  opacity = input$fit_opacity, color = I(fit_col),
                  facecolor = rep(fit_col, length(msh$i)),
                  vertexcolor = rep(fit_col, length(msh$x)),
                  flatshading = TRUE,
                  lighting = list(ambient = 1, diffuse = 0, specular = 0, roughness = 1, fresnel = 0),
                  name = paste("Ring (MLE)", cluster_id), showlegend = TRUE, inherit = FALSE,
                  hoverinfo = "skip", hovertemplate = ""
                )
              }
            } else {
              rp <- mle_sample_ring(f, n = 500)
              fig <- fig %>% add_trace(
                x = rp$x, y = rp$y, z = rp$z,
                type = "scatter3d", mode = "lines",
                line = list(color = fit_col, width = 4),
                name = paste("Ring (MLE)", cluster_id), showlegend = TRUE, inherit = FALSE,
                hoverinfo = "skip", hovertemplate = ""
              )
            }
          }

        } else if (input$fit_type == "arc") {
          f <- do.call(mle_fit_arc, c(list(X = X, w = w, arc_span_hint = input$arc_span_deg %||% 0), center_args))
          mle_fit_cache[[as.character(cluster_id)]] <- f
          if (!is.null(f)) {
            if (isTRUE(input$show_fit_tube)) {
              tube_r <- max(1, input$fit_tube_radius %||% 5)
              msh <- tryCatch(make_arc_tube_mesh(
                f$center, f$normal, f$radius,
                f$t_center, f$arc_half, tube_r
              ), error = function(e) NULL)
              if (!is.null(msh)) {
                fig <- fig %>% add_trace(
                  type = "mesh3d", x = msh$x, y = msh$y, z = msh$z,
                  i = msh$i, j = msh$j, k = msh$k,
                  opacity = input$fit_opacity, color = I(fit_col),
                  facecolor = rep(fit_col, length(msh$i)),
                  vertexcolor = rep(fit_col, length(msh$x)),
                  flatshading = TRUE,
                  lighting = list(ambient = 1, diffuse = 0, specular = 0, roughness = 1, fresnel = 0),
                  name = paste("Arc (MLE)", cluster_id), showlegend = TRUE, inherit = FALSE,
                  hoverinfo = "skip", hovertemplate = ""
                )
              }
            } else {
              rp <- mle_sample_arc(f, n = 300)
              fig <- fig %>% add_trace(
                x = rp$x, y = rp$y, z = rp$z,
                type = "scatter3d", mode = "lines",
                line = list(color = fit_col, width = 4),
                name = paste("Arc (MLE)", cluster_id), showlegend = TRUE, inherit = FALSE,
                hoverinfo = "skip", hovertemplate = ""
              )
            }
          }
        }
      }
    }
  }
  
  if (!is.null(input$fit_type) &&
      input$fit_type != "none" &&
      isTRUE(input$show_fit_distances)) {
    
    df_full <- plottable_data()
    df_cluster2 <- df_full %>% dplyr::filter(!is.na(cluster2) & cluster2 != 0)
    clusters <- df_cluster2 %>% dplyr::group_by(cluster2) %>% dplyr::group_split()
    
    for (cluster_data in clusters) {
      # medians for distance computation
      if (isTRUE(input$enable_clustering)) {
        # Use de-blinked molecule medians
        med <- cluster_data %>%
          dplyr::filter(!is.na(cluster) & cluster != 0) %>%
          dplyr::group_by(cluster) %>%
          dplyr::summarise(
            mx     = median(x), my = median(y), mz = median(z),
            spread = if (dplyr::n() > 1)
              sqrt(mean((x - mx)^2 + (y - my)^2 + (z - mz)^2)) else 0,
            n_locs = dplyr::n(),
            .groups = "drop"
          ) %>%
          dplyr::rename(x = mx, y = my, z = mz)
      } else {
        # Default: TID medians
        med <- cluster_data %>%
          dplyr::group_by(tid) %>%
          dplyr::summarise(
            mx     = median(x), my = median(y), mz = median(z),
            spread = if (dplyr::n() > 1)
              sqrt(mean((x - mx)^2 + (y - my)^2 + (z - mz)^2)) else 0,
            n_locs = dplyr::n(),
            .groups = "drop"
          ) %>%
          dplyr::rename(x = mx, y = my, z = mz)
      }
      if (nrow(med) < 4) next
      X <- as.matrix(med[, c("x","y","z")])
      w <- if (isTRUE(input$fit_weighted)) {
        s   <- med$spread
        n_l <- if ("n_locs" %in% names(med)) med$n_locs else rep(1L, nrow(med))
        s_ref <- stats::median(s[n_l > 1 & s > 0], na.rm = TRUE)
        if (!is.finite(s_ref) || s_ref <= 0) s_ref <- max(s[s > 0], na.rm = TRUE)
        if (!is.finite(s_ref) || s_ref <= 0) s_ref <- 1
        s_eff <- ifelse(n_l > 1 & s > 0, s, s_ref)
        n_l / s_eff^2
      } else NULL
      
      if (input$fit_method == "heuristic") {
        if (input$fit_type == "sphere") {
          g <- fit_sphere_geometry(X)
          r <- cap_radius_for_cluster(g$radius, X)
          if (is.finite(r) && r > 0) {
            pr <- mle_project_sphere(X, g$center, r)
            med$nearest_x <- pr$nearest[,1]; med$nearest_y <- pr$nearest[,2]; med$nearest_z <- pr$nearest[,3]
            med$distance_to_fit <- pr$dist
            fig <- add_distance_lines_to_plot(fig, med, max_dist = if (isTRUE(is.finite(input$max_dist_lines_nm) && input$max_dist_lines_nm > 0)) input$max_dist_lines_nm else Inf)
          }
        } else if (input$fit_type == "ellipsoid") {
          g <- fit_ellipsoid_geometry(X)
          if (!is.null(g)) {
            g$axes <- cap_axes_for_cluster(g$axes, X)
            if (all(is.finite(g$axes)) && all(g$axes > 0)) {
              pr <- mle_project_ellipsoid(X, g$center, g$axes, g$R)
              med$nearest_x <- pr$nearest[,1]; med$nearest_y <- pr$nearest[,2]; med$nearest_z <- pr$nearest[,3]
              med$distance_to_fit <- pr$dist
              fig <- add_distance_lines_to_plot(fig, med, max_dist = if (isTRUE(is.finite(input$max_dist_lines_nm) && input$max_dist_lines_nm > 0)) input$max_dist_lines_nm else Inf)
            }
          }
        } else if (input$fit_type == "helix") {
          g <- fit_helix_geometry(X, turns = input$helix_turns, n_points = 3000)
          t_h <- .project_t_helix(X, g$proj_center, g$axis, g$radius, g$pitch, g$phi0,
                                   check_neighbors = TRUE)
          H_h <- t(vapply(t_h, .helix_point, numeric(3),
                          c = g$proj_center, u = g$axis, r = g$radius, p = g$pitch, phi0 = g$phi0))
          med$nearest_x <- H_h[,1]; med$nearest_y <- H_h[,2]; med$nearest_z <- H_h[,3]
          med$distance_to_fit <- sqrt(rowSums((X - H_h)^2))
          fig <- add_distance_lines_to_plot(fig, med, max_dist = if (isTRUE(is.finite(input$max_dist_lines_nm) && input$max_dist_lines_nm > 0)) input$max_dist_lines_nm else Inf)

        } else if (input$fit_type == "ring") {
          g  <- fit_ring_geometry(X)
          pr <- mle_project_ring(X, g$center, g$normal, g$radius)
          med$nearest_x <- pr$nearest[,1]; med$nearest_y <- pr$nearest[,2]; med$nearest_z <- pr$nearest[,3]
          med$distance_to_fit <- pr$dist
          fig <- add_distance_lines_to_plot(fig, med, max_dist = if (isTRUE(is.finite(input$max_dist_lines_nm) && input$max_dist_lines_nm > 0)) input$max_dist_lines_nm else Inf)

        } else if (input$fit_type == "arc") {
          g  <- fit_arc_geometry(X, arc_span_hint = input$arc_span_deg %||% 0)
          pr <- mle_project_arc(X, g$center, g$normal, g$radius, g$t_center, g$arc_half)
          med$nearest_x <- pr$nearest[,1]; med$nearest_y <- pr$nearest[,2]; med$nearest_z <- pr$nearest[,3]
          med$distance_to_fit <- pr$dist
          fig <- add_distance_lines_to_plot(fig, med, max_dist = if (isTRUE(is.finite(input$max_dist_lines_nm) && input$max_dist_lines_nm > 0)) input$max_dist_lines_nm else Inf)
        }

      } else {
        # --- MLE: reuse cached fit from the drawing loop above (same data, same settings).
        # Falls back to refitting only if the cache is missing (e.g. show_fit_distances
        # was enabled on a render where the shape loop was skipped).
        cl_key      <- as.character(unique(cluster_data$cluster2))
        center_args <- build_center_args(X, input$mle_center_mode, input$mle_center_lambda)

        if (input$fit_type == "sphere") {
          f <- mle_fit_cache[[cl_key]] %||%
                 do.call(mle_fit_sphere, c(list(X = X, w = w), center_args))
          if (!is.null(f)) {
            f$radius <- cap_radius_for_cluster(f$radius, X)
            if (!is.finite(f$radius) || f$radius <= 0) next
            pr <- mle_project_sphere(X, f$center, f$radius)
            med$nearest_x <- pr$nearest[,1]; med$nearest_y <- pr$nearest[,2]; med$nearest_z <- pr$nearest[,3]
            med$distance_to_fit <- pr$dist
            fig <- add_distance_lines_to_plot(fig, med, max_dist = if (isTRUE(is.finite(input$max_dist_lines_nm) && input$max_dist_lines_nm > 0)) input$max_dist_lines_nm else Inf)
          }

        } else if (input$fit_type == "ellipsoid") {
          f <- mle_fit_cache[[cl_key]] %||%
                 do.call(mle_fit_ellipsoid, c(list(X = X, w = w), center_args))
          if (!is.null(f)) {
            f$axes <- cap_axes_for_cluster(f$axes, X)
            if (!all(is.finite(f$axes)) || any(f$axes <= 0)) next
            pr <- mle_project_ellipsoid(X, f$center, f$axes, f$R)
            med$nearest_x <- pr$nearest[,1]; med$nearest_y <- pr$nearest[,2]; med$nearest_z <- pr$nearest[,3]
            med$distance_to_fit <- pr$dist
            fig <- add_distance_lines_to_plot(fig, med, max_dist = if (isTRUE(is.finite(input$max_dist_lines_nm) && input$max_dist_lines_nm > 0)) input$max_dist_lines_nm else Inf)
          }

        } else if (input$fit_type == "helix") {
          f <- mle_fit_cache[[cl_key]] %||%
                 do.call(mle_fit_helix, c(list(X = X, w = w, turns_hint = input$helix_turns), center_args))
          if (!is.null(f)) {
            t_mle <- .project_t_helix(X, f$center, f$axis, f$radius, f$pitch, f$phi0,
                                      t0 = f$t, check_neighbors = TRUE)
            H_mle <- t(vapply(t_mle, .helix_point, numeric(3),
                              c = f$center, u = f$axis, r = f$radius, p = f$pitch, phi0 = f$phi0))
            med$nearest_x <- H_mle[,1]; med$nearest_y <- H_mle[,2]; med$nearest_z <- H_mle[,3]
            med$distance_to_fit <- sqrt(rowSums((X - H_mle)^2))
            fig <- add_distance_lines_to_plot(fig, med, max_dist = if (isTRUE(is.finite(input$max_dist_lines_nm) && input$max_dist_lines_nm > 0)) input$max_dist_lines_nm else Inf)
          }

        } else if (input$fit_type == "ring") {
          f <- mle_fit_cache[[cl_key]] %||%
                 do.call(mle_fit_ring, c(list(X = X, w = w), center_args))
          if (!is.null(f)) {
            pr <- mle_project_ring(X, f$center, f$normal, f$radius)
            med$nearest_x <- pr$nearest[,1]; med$nearest_y <- pr$nearest[,2]; med$nearest_z <- pr$nearest[,3]
            med$distance_to_fit <- pr$dist
            fig <- add_distance_lines_to_plot(fig, med, max_dist = if (isTRUE(is.finite(input$max_dist_lines_nm) && input$max_dist_lines_nm > 0)) input$max_dist_lines_nm else Inf)
          }

        } else if (input$fit_type == "arc") {
          f <- mle_fit_cache[[cl_key]] %||%
                 do.call(mle_fit_arc, c(list(X = X, w = w, arc_span_hint = input$arc_span_deg %||% 0), center_args))
          if (!is.null(f)) {
            pr <- mle_project_arc(X, f$center, f$normal, f$radius, f$t_center, f$arc_half)
            med$nearest_x <- pr$nearest[,1]; med$nearest_y <- pr$nearest[,2]; med$nearest_z <- pr$nearest[,3]
            med$distance_to_fit <- pr$dist
            fig <- add_distance_lines_to_plot(fig, med, max_dist = if (isTRUE(is.finite(input$max_dist_lines_nm) && input$max_dist_lines_nm > 0)) input$max_dist_lines_nm else Inf)
          }
        }
      }
    }
  }

  if (input$dark_mode) {
    fig <- fig %>%
      layout(
        paper_bgcolor = 'black',
        plot_bgcolor = 'black',
        scene = list(
          bgcolor = "black",
          xaxis = list(
            color = "white", gridcolor = "white", zerolinecolor = "white",
            titlefont = list(color = "white"), tickfont = list(color = "white")
          ),
          yaxis = list(
            color = "white", gridcolor = "white", zerolinecolor = "white",
            titlefont = list(color = "white"), tickfont = list(color = "white")
          ),
          zaxis = list(
            color = "white", gridcolor = "white", zerolinecolor = "white",
            titlefont = list(color = "white"), tickfont = list(color = "white")
          )
        )
      )
  }
  grid_color <- if (input$hide_grid) "rgba(0,0,0,0)" else if (input$dark_mode) "white" else "black"
  
  fig <- fig %>%
    layout(
      paper_bgcolor = if (input$dark_mode) "black" else "white",
      plot_bgcolor  = if (input$dark_mode) "black" else "white",
      scene = list(
        bgcolor = if (input$dark_mode) "black" else "white",
        xaxis = list(
          showgrid = !input$hide_grid,
          zeroline = !input$hide_grid,
          showticklabels = !input$hide_grid,
          title = if (input$hide_grid) "" else "X",
          color = if (input$dark_mode) "white" else "black"
        ),
        yaxis = list(
          showgrid = !input$hide_grid,
          zeroline = !input$hide_grid,
          showticklabels = !input$hide_grid,
          title = if (input$hide_grid) "" else "Y",
          color = if (input$dark_mode) "white" else "black"
        ),
        zaxis = list(
          showgrid = !input$hide_grid,
          zeroline = !input$hide_grid,
          showticklabels = !input$hide_grid,
          title = if (input$hide_grid) "" else "Z",
          color = if (input$dark_mode) "white" else "black"
        )
      )
    )
  
  if (isTRUE(input$show_cluster2_shape) && isTRUE(input$enable_cluster2) && !is.null(df$cluster2)) {
    # One wrap per Median Cluster PER CHANNEL, the same rule the isosurface
    # follows, so an alpha shape can never join two colours (wrap_groups()).
    a_groups   <- wrap_groups(df, "cluster2")
    # All points, the same test that colors the points (see n_frame_channels).
    a_multi_ch <- n_frame_channels(plottable_data()) > 1L
    for (gi in seq_len(nrow(a_groups))) {
      g   <- a_groups[gi, , drop = FALSE]
      cl2 <- g$id
      group_df <- wrap_group_rows(df, "cluster2", g)
      if (nrow(group_df) >= 4) {
        # Radius derived per cluster from its own point spacing - see
        # alpha_surface_3d() in R/smlm_shapes.R.
        srf <- alpha_surface_3d(as.matrix(group_df[, c("x", "y", "z")]))
        if (!is.null(srf)) {

          # pick a single solid color: first non-NA cluster2 color (fallback to light gray)
          col_hex <- group_df$cluster2_color[which(!is.na(group_df$cluster2_color) & nzchar(group_df$cluster2_color))][1]
          if (is.na(col_hex) || !nzchar(col_hex)) col_hex <- "#bdbdbd"
          # With several channels, the shared "Color wraps by channel" option
          # takes the channel's colour instead, as the isosurface does.
          if (a_multi_ch && isTRUE(input$cluster2_iso_by_channel) &&
              !is.na(g$color) && nzchar(g$color)) col_hex <- g$color

          fig <- fig %>% plotly::add_trace(
            type = "mesh3d",
            x = srf$x, y = srf$y, z = srf$z,   # srf may be subsampled; i/j/k index these
            i = srf$i, j = srf$j, k = srf$k,
            opacity = input$hull_opacity,
            # Colour as a flat intensity + one-colour scale, the encoding
            # iso_mesh_trace() uses. `color = I(hex)` on a mesh3d is DROPPED by
            # plotly_build() (with or without I()), so this wrap used to render
            # in plotly's default mesh colour whatever col_hex said.
            intensity = rep(0, length(srf$x)),
            colorscale = list(list(0, col_hex), list(1, col_hex)),
            showscale = FALSE,
            name = if (a_multi_ch && nzchar(g$label %||% ""))
                     paste0("Median Cluster ", cl2, " (", g$label, ")")
                   else paste("Median Cluster", cl2),
            hoverinfo = "skip",
            hovertemplate = "",
            showlegend = FALSE,
            inherit = FALSE
          )
        }
      }
    }
  }

  # ---- Estimated ground-truth candidates (hollow spheres, green->red) ----
  ge_3d <- gt_estimate()
  if (!is.null(ge_3d) && nrow(ge_3d$points)) {
    gpts <- ge_3d$points
    lo3  <- ge_3d$meta$link_nm %||% 0
    hi3  <- max(ge_3d$meta$footprint_nm %||% 0, 2 * lo3,
                stats::quantile(gpts$pos_sd_nm, 0.9, na.rm = TRUE))
    col3 <- gt_uncertainty_colors(gpts$pos_sd_nm, lo3, hi3)
    sz3  <- suppressWarnings(as.numeric(input$gt_circle_size))
    sz3  <- if (!is.finite(sz3) || sz3 <= 0) 7 else max(3, sz3 * 0.5)
    struct3 <- ifelse(gpts$cluster2 > 0, paste0("Structure ", gpts$cluster2), "ungrouped")
    txt3 <- sprintf("<b>%s</b><br>GT candidate · ± %.0f nm%s<br>%d contributing unit(s) · %d TID / %d loc",
                    struct3, gpts$pos_sd_nm, ifelse(gpts$fused, " (fused: multiple GTs)", ifelse(gpts$count_fused %||% FALSE, " (more TIDs than one molecule)", "")),
                    gpts$n_units, gpts$n_tid, gpts$n_loc)

    # Contribution fan-out lines in 3D.
    if (isTRUE(input$gt_show_lines %||% TRUE) && !is.null(ge_3d$links) && nrow(ge_3d$links)) {
      lk <- ge_3d$links
      lx <- as.vector(rbind(lk$ux, lk$gx, NA_real_))
      ly <- as.vector(rbind(lk$uy, lk$gy, NA_real_))
      lz <- as.vector(rbind(lk$uz, lk$gz, NA_real_))
      line_col <- if (isTRUE(input$dark_mode)) "rgba(255,255,255,0.35)" else "rgba(0,0,0,0.35)"
      fig <- fig %>% plotly::add_trace(
        x = lx, y = ly, z = lz, type = "scatter3d", mode = "lines",
        line = list(color = line_col, width = 2), name = "GT contributions",
        showlegend = TRUE, inherit = FALSE, hoverinfo = "skip"
      )
    }
    fig <- fig %>% plotly::add_trace(
      x = gpts$x, y = gpts$y, z = gpts$z,
      type = "scatter3d", mode = "markers",
      marker = list(symbol = "circle-open", size = sz3,
                    color = col3, line = list(color = col3, width = 2)),
      name = "Estimated GT", showlegend = TRUE, inherit = FALSE,
      text = txt3, hoverinfo = "text"
    )
  }

  fig %>% plotly::event_register('plotly_relayout')
})
