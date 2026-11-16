## ─────────────────────────────────────────────────────────────────────────────
## server/smlm/17_plottable_data.R  -  plottable_data: THE analysis-ready data of the active file
## ─────────────────────────────────────────────────────────────────────────────
## Part of smlm_server(). Not a standalone script: its statements are evaluated in
## the server function's environment, in file-name order (server/parts.R).
##
## Crops, rogue removal, grouping and clustering (apply_clustering_pipeline),
## then colours and hover text. Recomputes on plot_trigger only.
##
## INTERFACE (generated: Rscript tools/server_map.R --headers; do not edit by hand)
##   DEFINES for other parts:
##     plottable_data
##   READS from other parts:
##     01_view_upload: smlm_2d_crop_rv, uploaded_df_rv, z_scale_applied
##     02_multiplex: .mplex_clu_memo, .mplex_off_of, .mplex_residual_shift,
##                   .mplex_store_off
##     05_filter_pipeline: filtered_df
##     10_triggers_stores: fit_log_rv, plot_trigger
##     16_crop_range_ui: apply_trimmer
## END INTERFACE
## ─────────────────────────────────────────────────────────────────────────────

# make_two_color_tid_map() was a third copy of the DCR/fluo trace split, kept
# alongside assign_two_color_group() and assign_analysis_group(). All three are
# now resolve_channels() in R/smlm_channels.R, which is the only place a channel
# is defined.


# two_color_nn_edges() lived here: a second NN builder that existed only to
# colour lines by channel pair. build_nn_segments() returns that colouring now,
# so the shared NN panel does the whole job and this one is gone.




# ==========================================================================
# plottable_data() - THE ANALYSIS-READY DATASET
# ==========================================================================
#
# GOAL
#   The single hand-off point between the filter pipeline and everything that
#   analyses or draws the data. Nearly every analysis below begins with
#   `df <- plottable_data()`. If a filter change is not reaching an analysis,
#   this is the first place to look.
#
#   It is the SMLM counterpart of trk() in the tracking tab, but it does more:
#   as well as delivering the filtered data it performs the CLUSTERING, so
#   every consumer sees the same cluster / cluster2 assignment. That is
#   deliberate - if each analysis clustered independently they could disagree.
#
# WHAT IT ADDS ON TOP OF filtered_df()
#   1. apply_trimmer() ..... optional outlier/edge trimming
#   2. 2D crop ............. applies the zoom rectangle SERVER-SIDE, so the
#                            crop affects the analyses too and behaves
#                            identically in 2D and 3D (rather than relying on
#                            a browser round-trip that only changes the view)
#   3. CLUSTERING .......... de-blinking and/or median clustering
#                            (structure) via apply_clustering_pipeline().
#                            Multiplex datasets are clustered PER CHANNEL.
#   4. colour columns ...... per-point colour for the chosen colour mode
#
# WHY eventReactive(plot_trigger())
#   Clustering a large localization table is expensive. Recomputing on every
#   keystroke would make the app unusable, so the analyses refresh when Plot is
#   pressed. Consequence worth knowing: after changing a clustering parameter
#   you must press Plot before the numbers update.
#
# RETURNS the filtered data.frame plus (when enabled):
#   cluster ..... de-blinked molecule ID (no 0: de-blinking assigns every trace)
#   cluster2 .... secondary/structure ID (0 = noise)
#   colour columns, and the multiplex channel/label columns when active
# ==========================================================================
plottable_data <- eventReactive(plot_trigger(), {
  run_set <- smlm_settings_from_input(input)   # one settings snapshot for this run
  df <- filtered_df()
  # Multiplex channel alignment: shift not yet in the active copy (see the
  # MULTIPLEX CHANNEL ALIGNMENT block). Before the crops, so a crop box drawn
  # on the aligned view selects the aligned points.
  df <- .mplex_residual_shift(df)
  df <- apply_trimmer(df, input)
  # 2D zoom crop - applied server-side so it works in both 2D and 3D without
  # relying on the browser round-trip from update* calls.
  crop2d <- isolate(smlm_2d_crop_rv())
  if (!is.null(crop2d) && nrow(df) > 0)
    df <- df[df$x >= crop2d$x[1] & df$x <= crop2d$x[2] &
             df$y >= crop2d$y[1] & df$y <= crop2d$y[2], , drop = FALSE]
  if (nrow(df) == 0) return(df)
  # Remove rogue traces
  if (isTRUE(input$remove_rogues)) {
    rr <- remove_rogue_traces(
      df,
      rule   = input$rogue_rule %||% "factor_median",
      factor = input$rogue_factor %||% 2.5,
      min_n  = input$sd_min_locs %||% SMLM_DEFAULTS$sd_min_locs,
    )
    df <- rr$df
    #log text
    fit_log_rv(
      paste0(
        "Rogue removal: removed ", rr$removed_n, " trace(s). Rule: ", rr$note
      )
    )
  }
  
  if (nrow(df) == 0) return(df)
  
  df <- assign_analysis_group(df, run_set)
  
  # Safety check
  if (nrow(df) == 0) return(df)
  
  # --- Precompute per-trace median DCR (used by unmixer + optional color scale) ---
  if ("dcr" %in% names(df)) {
    tid_dcr <- df %>%
      dplyr::group_by(tid) %>%
      dplyr::summarise(tid_median_dcr = median(dcr, na.rm = TRUE), .groups = "drop")
    df <- df %>% dplyr::left_join(tid_dcr, by = "tid")
  }
  
  # 2) Base color assignment (use a single name everywhere)
  qual_col_pals <- RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == "qual", ]
  base_colors   <- unlist(mapply(RColorBrewer::brewer.pal,
                                 qual_col_pals$maxcolors,
                                 rownames(qual_col_pals)))
  
  # ── CLUSTERING: one call, shared with every batch tool ──────────────────
  # This block used to inline both pipelines - the single-colour one and a
  # second copy for 2-colour - which is why the old comments kept insisting
  # the two "must agree". They cannot disagree any more: the plot and the
  # batch tools now run the SAME function, and channels (none / 2-colour split
  # / multiplex) come from resolve_channels() inside it.
  # Multiplex: clustering runs per channel, and a channel nudge only
  # translates a channel, so it cannot change any cluster. When the same rows
  # meet the same settings, reuse the last result and move its coordinates by
  # the change in offset. Anything else (a crop that now selects other points,
  # a changed setting, another dataset) misses the key and re-clusters.
  mp_key <- NULL
  # The pipeline reads `input` lazily (as_smlm_settings), so this reactive depends
  # only on the settings the current branch uses. Only the multiplex memo needs
  # every setting at once, for its key; there the one object serves key AND run.
  mp_set  <- input
  mp_act  <- isolate(uploaded_df_rv())
  mp_spec <- attr(mp_act, "mplex_spec")
  mp_off  <- .mplex_store_off(attr(mp_act, "mplex_store_name"))
  if ("mplex_channel" %in% names(df) && !is.null(mp_spec)) {
    mp_set <- smlm_settings_from_input(input)
    mp_key <- digest::digest(list(
      df$tid, df$tim, df$mplex_channel,                      # which rows
      if (".dataset_tag" %in% names(df)) df$.dataset_tag[1], # which build
      isolate(z_scale_applied()),
      smlm_settings_key(mp_set),                             # every clustering setting
      input$remove_rogues, input$rogue_rule, input$rogue_factor
    ), algo = "xxhash64")
  }
  memo <- .mplex_clu_memo
  if (!is.null(mp_key) && identical(memo$key, mp_key)) {
    cl <- memo$df
    for (c in mp_spec) {                                     # move to the current offsets
      k <- .mplex_key(c)
      delta <- .mplex_off_of(mp_off, k) - .mplex_off_of(memo$off, k)
      if (all(abs(delta) < 1e-9)) next
      rows <- which(cl$mplex_channel == c$channel)
      cl$x[rows] <- cl$x[rows] + delta[["x"]]
      cl$y[rows] <- cl$y[rows] + delta[["y"]]
      if ("z" %in% names(cl)) cl$z[rows] <- cl$z[rows] + delta[["z"]]
    }
    df <- cl
  } else {
    df <- apply_clustering_pipeline(df, mp_set)
  }
  if (!is.null(mp_key)) {
    memo$key <- mp_key; memo$df <- df; memo$off <- mp_off
  } else {
    memo$key <- NULL; memo$df <- NULL                        # free the memory
  }
  # Colour is presentation, so it is applied here rather than inside the cached
  # chain - see attach_channel_colors(). This is what keeps the colour pickers
  # out of the prepared-data cache key.
  df <- attach_channel_colors(df, smlm_channel_colors(input))
  n_ch <- n_frame_channels(df)

  # ── Colours ─────────────────────────────────────────────────────────────
  # Colouring is a display concern, so it stays here rather than moving into
  # the pipeline.
  if (isTRUE(input$enable_clustering)) {
    pal <- rep(base_colors, length.out = max(1, max(df$cluster, na.rm = TRUE) + 1))
    df$color <- ifelse(df$cluster == 0, "#bdbdbd",
                       pal[(df$cluster %% length(pal)) + 1])
  } else {
    extended_palette_tid <- rep(base_colors, length.out = nrow(df))
    df$color <- extended_palette_tid[
      (as.numeric(as.factor(df$tid)) %% length(extended_palette_tid)) + 1]
  }

  if (isTRUE(input$enable_cluster2)) {
    # Members per Median Cluster, counted in whatever unit fed the clustering.
    id_col <- if (isTRUE(input$enable_clustering)) "cluster" else "tid"
    c2_counts <- df %>%
      dplyr::filter(!is.na(cluster2) & cluster2 != 0) %>%
      dplyr::group_by(cluster2) %>%
      dplyr::summarise(n_in_cluster2 = dplyr::n_distinct(.data[[id_col]]), .groups = "drop")
    df <- df %>%
      dplyr::select(-dplyr::any_of("n_in_cluster2")) %>%
      dplyr::left_join(c2_counts, by = "cluster2")
    df$n_in_cluster2 <- as.integer(df$n_in_cluster2)

    max_c2 <- suppressWarnings(max(df$cluster2, na.rm = TRUE))
    if (is.finite(max_c2) && max_c2 > 0) {
      pal2 <- rep(base_colors, length.out = max_c2 + 1)
      df$cluster2_color <- ifelse(is.na(df$cluster2) | df$cluster2 == 0, NA_character_,
                                  pal2[(df$cluster2 %% length(pal2)) + 1])
    } else {
      df$cluster2_color <- NA_character_
    }
  } else {
    df$cluster2_color <- NA_character_
    df$n_in_cluster2  <- NA_integer_
  }
  df$curvature <- NA_real_

  # Channel colours override the per-cluster palette, because when channels
  # exist the thing the user is looking for is which channel a point is in.
  # Numeric colour modes (z / time / density) still override this downstream.
  if (n_ch > 1L && "channel_color" %in% names(df) && any(!is.na(df$channel_color)))
    df$color <- ifelse(is.na(df$channel_color), df$color, df$channel_color)

  if (!"two_color_group" %in% names(df)) df$two_color_group <- NA_integer_
  # cluster_tc1 / cluster_tc2 are gone: they were aliases of cluster / cluster2
  # kept while the two pipelines were being merged, and every reader now uses
  # the real columns.
  
  
  # If requested, (re)attach TID-median DCR for hover text
  if (isTRUE(input$color_tid_by_median_dcr) && "dcr" %in% names(df)) {
    tid_dcr <- df %>%
      dplyr::group_by(tid) %>%
      dplyr::summarise(tid_median_dcr = median(dcr, na.rm = TRUE), .groups = "drop")
    df <- df %>% dplyr::left_join(tid_dcr, by = "tid")
  }
  # median DCR to the summarized df (trace or cluster) so hover is stable
  if (input$show_medians) {
    if (input$enable_clustering) {
      # cluster-level medians for DCR
      cl_dcr <- df %>%
        dplyr::filter(!is.na(cluster) & cluster != 0) %>%
        dplyr::group_by(cluster) %>%
        dplyr::summarise(cluster_median_dcr = median(dcr, na.rm = TRUE), .groups = "drop")
      df <- df %>% dplyr::left_join(cl_dcr, by = "cluster")
    } else {
      # trace-level medians for DCR
      tid_dcr <- df %>%
        dplyr::group_by(tid) %>%
        dplyr::summarise(tid_median_dcr = median(dcr, na.rm = TRUE), .groups = "drop")
      df <- df %>% dplyr::left_join(tid_dcr, by = "tid")
    }
  }
  # 7) Hover text 
  show_dcr_hover  <- isTRUE(input$enable_two_color) && ("dcr" %in% names(df))
  show_fluo_hover <- isTRUE(input$color_by_fluo) && ("fluo" %in% names(df))

  # A-#/B-# labels come straight off cluster/cluster2 now - there is one
  # clustering, so the "2-colour cluster" was never a different object.
  ch_lbl  <- if ("channel_label" %in% names(df)) df$channel_label else NULL
  tc1_lbl <- format_tc_id(df$cluster,  ch_lbl)
  tc2_lbl <- format_tc_id(df$cluster2, ch_lbl)
  
  df <- df %>% dplyr::mutate(
    hover_text = paste0(
      if ("tid" %in% names(.)) paste0("Trace ID: ", tid) else "",
      "<br>X: ", round(x, 2),
      "<br>Y: ", round(y, 2),
      if ("z" %in% names(.)) paste0("<br>Z: ", round(z, 2)) else "",
      if (!all(is.na(cluster)))  paste0("<br>De-blinked ID: ", cluster) else "",
      if (!all(is.na(cluster2))) paste0(
        "<br>Median Cluster: ", cluster2,
        ifelse(!is.na(n_in_cluster2), paste0(" (n=", n_in_cluster2, ")"), "")
      ) else "",

      # show group label only in the mode that defines it
      if (show_fluo_hover && "fluo_med" %in% names(.))
        paste0("<br>FLUO group: ", suppressWarnings(as.integer(round(fluo_med))))
      else "",
      
      # DCR hover text only when two-color is enabled
      if (show_dcr_hover && "tid_median_dcr" %in% names(.))
        paste0("<br>TID median DCR: ", round(tid_median_dcr, 2))
      else "",
      
      if (show_dcr_hover && "cluster_median_dcr" %in% names(.))
        paste0("<br>Cluster median DCR: ", round(cluster_median_dcr, 2))
      else "",
      
      if (show_dcr_hover && "dcr" %in% names(.))
        paste0("<br>Raw DCR: ", round(dcr, 2))
      else "",
      ifelse(
        isTRUE(input$enable_two_color) & !is.na(tc1_lbl),
        paste0("<br>Tc-De-blinked ID: ", tc1_lbl),
        ""),
      ifelse(
        isTRUE(input$enable_two_color) & !is.na(tc2_lbl),
        paste0("<br>Tc-Median Cluster: ", tc2_lbl),
        ""
      )
      )
  )

  # Multiplex: append the channel label to hover when present
  if ("mplex_label" %in% names(df))
    df$hover_text <- paste0(df$hover_text, "<br>Channel: ", df$mplex_label)

  # Derive globally-unique pcid / scid using the per-file random tag
  .tag <- if (".dataset_tag" %in% names(df)) df$.dataset_tag[[1L]] else "0000"
  df$pcid <- ifelse(
    !is.na(df$cluster) & df$cluster > 0L,
    paste0(.tag, "-", sprintf("%04d", df$cluster)),
    NA_character_
  )
  df$scid <- ifelse(
    !is.na(df$cluster2) & df$cluster2 > 0L,
    paste0(.tag, "-", sprintf("%04d", df$cluster2)),
    NA_character_
  )

  df
})
