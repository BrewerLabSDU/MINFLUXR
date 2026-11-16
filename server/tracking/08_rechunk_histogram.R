## ─────────────────────────────────────────────────────────────────────────────
## server/tracking/08_rechunk_histogram.R  -  filter and rechunk on Plot Data; metadata histogram
## ─────────────────────────────────────────────────────────────────────────────
## Part of tracking_server(). Not a standalone script: its statements are evaluated in
## the server function's environment, in file-name order (server/parts.R).
##
## The observer that runs the whole filter pipeline on the original upload and
## rewrites the playback chunks, and the metadata histogram.
##
## INTERFACE (generated: Rscript tools/server_map.R --headers; do not edit by hand)
##   DEFINES for other parts:
##     (nothing)
##   READS from other parts:
##     01_state_sigma: .tids_per_chunk, rechunk_request, trk_chunk_pref,
##                     trk_raw_stats_rv, trk_rechunk_bbox,
##                     trk_rechunk_stats_rv, trk_spatial_filter
##     05_filter_ui: reset_dur_filter_flag, trk_len_reset_flag
##     06_filter_pipeline: trk_dur_filtered_d, trk_filter_params
##     07_cycle_time_sequence: trk_cluster_cycle, trk_cycle_time
##     10_triggers: trk_data_key
## END INTERFACE
## ─────────────────────────────────────────────────────────────────────────────

# ---- Filter + rechunk (driven by Plot Data) ---------------------------------
# Always reads from original_manifest_dir (frozen at upload) so repeated
# rechunks are idempotent - Reset Ranges + Rechunk restores the full dataset.
# Writes filtered result to a new dir that becomes the new manifest_dir.
observeEvent(rechunk_request(), {
  req(!is.null(rv$original_manifest_dir))
  manf <- file.path(rv$original_manifest_dir, "manifest.csv")
  req(file.exists(manf))
  # Read BEFORE anything is rebuilt: was the user scoped to a chunk, or on all TIDs?
  sel_before <- isolate(input$trk_chunk_sel)

  withProgress(message = "Applying filters: loading chunks...", value = 0, {
    man <- readr::read_csv(manf, show_col_types = FALSE)
    n_chunks <- nrow(man)

    # The SAME pipeline the live view runs (filter_pipeline, R/est_clean.R), on all
    # of the original data at once and with the stored filters, so what Plot Data
    # bakes into the new chunks cannot drift from what the analyses read. (This
    # used to be a second, hand-written copy of every stage, which had drifted.)
    incProgress(0.1, detail = "Reading chunks")
    df_all <- dplyr::bind_rows(Filter(Negate(is.null), lapply(man$file, read_chunk_file)))
    if (!nrow(df_all)) {
      showNotification("No data to filter  -  rechunk aborted.", type = "error", duration = 8)
      return(invisible(NULL))
    }
    # Same post-processing as trk_raw_df
    df_all$tid <- as.character(df_all$tid)
    df_all$x   <- as.numeric(df_all$x)
    df_all$y   <- as.numeric(df_all$y)
    df_all$tim <- as.numeric(df_all$tim) * 1000           # seconds -> ms
    fp <- trk_filter_params()
    fp$t_c <- isolate(trk_cycle_time())
    cc_r   <- tryCatch(isolate(trk_cluster_cycle()), error = function(e) NULL)
    fp$sw  <- if (!is.null(cc_r) && is.finite(cc_r$sw)) cc_r$sw else 0
    incProgress(0.3, detail = "Applying filters")
    df_all <- filter_pipeline(df_all, fp)
    if (nrow(df_all) == 0) {
      showNotification("No data survived the filters  -  rechunk aborted.", type = "error", duration = 8)
      return(invisible(NULL))
    }

    # Spatial filter - apply if active, using the same pct_thr threshold
    sf <- isolate(trk_spatial_filter())
    if (!is.null(sf) && "x" %in% names(df_all) && "y" %in% names(df_all)) {
      in_box   <- df_all$x >= sf$x[1] & df_all$x <= sf$x[2] &
                  df_all$y >= sf$y[1] & df_all$y <= sf$y[2]
      all_tids <- unique(df_all$tid)
      keep_sp  <- vapply(all_tids, function(tid) {
        idx <- df_all$tid == tid
        mean(in_box[idx], na.rm = TRUE) >= sf$pct_thr
      }, logical(1))
      df_all <- df_all[df_all$tid %in% all_tids[keep_sp], , drop = FALSE]
      if (nrow(df_all) == 0) {
        showNotification("No TIDs survived the spatial filter  -  rechunk aborted.", type = "error", duration = 8)
        return(invisible(NULL))
      }
    }

    # Convert tim back to seconds (trk_raw_df will multiply by 1000 on load)
    df_all$tim <- df_all$tim / 1000

    # Compute new tid_sum and sort longest-first
    tid_sum_new <- df_all %>%
      dplyr::group_by(tid) %>%
      dplyr::summarise(
        duration_ms = (max(tim) - min(tim)) * 1000,
        n_pts = dplyr::n(), .groups = "drop"
      ) %>%
      dplyr::arrange(dplyr::desc(duration_ms))
    df_all <- df_all[order(match(df_all$tid, tid_sum_new$tid)), , drop = FALSE]

    incProgress(0.6, detail = "Writing new chunks...")
    out_dir  <- file.path(tempdir(), paste0("rechunked_", as.integer(Sys.time())))
    user_tpc <- .tids_per_chunk()

    manf_new <- split_by_tid_to_dir(df_all, out_dir, tids_per_chunk = user_tpc, fmt = CHUNK_FORMAT)
    if (!file.exists(manf_new)) {
      showNotification("Rechunking failed (split_by_tid_to_dir returned no manifest).", type = "error", duration = 8)
      return(invisible(NULL))
    }

    readr::write_csv(tid_sum_new, file.path(out_dir, "tid_summary.csv"))
    # Carry the scoping, not the index: on a chunk before -> chunk 1 of the new
    # chunking; on "All TIDs" before -> still all. Set before the manifest swap so
    # the selector rebuild below already sees it.
    trk_chunk_pref(if (isTruthy(sel_before) && !identical(sel_before, "__all__"))
                     "first" else "__all__")
    old_rechunk         <- rv$last_rechunk_dir   # superseded by out_dir; removed below
    rv$last_rechunk_dir <- out_dir
    rv$manifest_dir     <- out_dir
    trk_data_key(as.character(as.numeric(Sys.time())))

    # Compute stats across ALL remaining TIDs so filter UIs auto-update to the
    # new data extent (not just the first chunk). tim is currently in seconds;
    # multiply by 1000 to convert to ms format used in the filter UIs and stats displays.
    loc_vec_r <- function(col) if (col %in% names(df_all)) df_all[[col]][is.finite(df_all[[col]])] else numeric(0)
    trk_rechunk_stats_rv(list(
      tid_length   = as.numeric(table(df_all$tid)),
      tid_duration = tid_sum_new$duration_ms,
      efo          = loc_vec_r("efo"),
      cfr          = loc_vec_r("cfr"),
      dcr          = loc_vec_r("dcr"),
      tim          = loc_vec_r("tim") * 1000,
      eco          = loc_vec_r("eco"),
      dwell        = loc_vec_r("dwell"),
      fluo         = loc_vec_r("fluo")
    ))
    # Lock plot axes to the full spatial extent of all remaining TIDs so the
    # view stays stable as individual TIDs are played through chunks.
    rx <- df_all$x[is.finite(df_all$x)]; ry <- df_all$y[is.finite(df_all$y)]
    if (length(rx) && length(ry))
      trk_rechunk_bbox(list(x = range(rx), y = range(ry)))

    # Auto-fit the total playback window to the longest remaining trace so the whole
    # longest track plays without being cut off. tid_sum_new is sorted longest-first,
    # so [1] is the max duration (ms). The user can still lower it afterwards.
    max_dur <- suppressWarnings(tid_sum_new$duration_ms[1])
    if (length(max_dur) == 1 && is.finite(max_dur) && max_dur > 0)
      updateNumericInput(session, "trk_play_total", value = signif(max_dur, 4))

    # The split (if any) is now BAKED into the saved segment TIDs. Live splitting
    # can stay ON: trk_gap_split strips the "_s<N>" suffix before splitting, so the
    # reloaded segments are recombined and re-split to the SAME segments rather than
    # split a second time. (That double split is what used to shrink the duration
    # upper bound on every rechunk.) Pressing Filter & Rechunk again is idempotent.
    # Flags stay TRUE so the length/duration sliders snap to the new full range
    # on their next (data-driven) re-render - no forced re-render here, which
    # would fire during the async checkbox round-trip while the split is still on.
    trk_len_reset_flag(TRUE)
    reset_dur_filter_flag(TRUE)

    n_tids <- length(unique(df_all$tid))
    n_new  <- nrow(readr::read_csv(manf_new, show_col_types = FALSE))
    showNotification(
      sprintf("%s TIDs kept, in %d chunk(s) @ %d TIDs/chunk. Showing %s.",
              format(n_tids, big.mark = ","), n_new, user_tpc,
              if (identical(isolate(trk_chunk_pref()), "first") && n_new >= 1)
                "chunk 1" else "all TIDs"),
      type = "message", duration = 8
    )

    # Each press writes a fresh chunk set to tempdir(). Now that this runs on
    # every Plot Data, the superseded one is removed rather than left to
    # accumulate copies of the dataset. Done last, so nothing is deleted unless
    # the new set was written and swapped in successfully. The frozen upload set
    # is never a candidate.
    if (!is.null(old_rechunk) && !identical(old_rechunk, out_dir) &&
        !identical(old_rechunk, rv$original_manifest_dir) && dir.exists(old_rechunk))
      unlink(old_rechunk, recursive = TRUE)
  })
}, ignoreInit = TRUE)

# SMLM-style metadata histogram: raw (grey) + filtered (blue) overlay
output$trk_meta_hist <- renderPlotly({
  param <- input$trk_hist_param
  req(!is.null(param), param != "none")
  n_bins   <- max(5L, as.integer(input$trk_hist_bins %||% 60L))
  show_raw <- isTRUE(input$trk_hist_show_raw)
  zoom     <- input$trk_hist_range  # may be NULL before UI renders

  # Grey bars: frozen at upload time, unaffected by filtering or rechunking
  raw_stats <- trk_raw_stats_rv()
  raw_vec   <- if (!is.null(raw_stats) && param %in% names(raw_stats))
                 raw_stats[[param]] else numeric(0)

  # Blue bars: fully-filtered data - debounced so the histogram only recomputes
  # after the user stops moving sliders, not on every tick.
  filt_df  <- tryCatch(trk_dur_filtered_d(), error = function(e) NULL)
  filt_vec <- if (param == "tid_length") {
    if (!is.null(filt_df) && nrow(filt_df) > 0 && "tid" %in% names(filt_df))
      as.numeric(table(filt_df$tid))
    else numeric(0)
  } else if (param == "tid_duration") {
    if (!is.null(filt_df) && nrow(filt_df) > 0 && all(c("tid","tim") %in% names(filt_df)))
      filt_df %>%
        dplyr::group_by(tid) %>%
        dplyr::summarise(dur = max(tim) - min(tim), .groups = "drop") %>%
        dplyr::pull(dur)
    else numeric(0)
  } else {
    if (!is.null(filt_df) && nrow(filt_df) > 0 && param %in% names(filt_df))
      filt_df[[param]][is.finite(filt_df[[param]])]
    else numeric(0)
  }

  # Subset both vectors to the zoom range (same approach as SMLM tab)
  if (!is.null(zoom) && length(zoom) == 2 && all(is.finite(zoom)) && zoom[1] < zoom[2]) {
    raw_vec  <- raw_vec[raw_vec  >= zoom[1] & raw_vec  <= zoom[2]]
    filt_vec <- filt_vec[filt_vec >= zoom[1] & filt_vec <= zoom[2]]
  }
  shiny::validate(shiny::need(
    length(filt_vec) > 0 || (show_raw && length(raw_vec) > 0),
    "No data to plot for the current filters/zoom range."
  ))

  xlab <- switch(param,
    tid_length   = "Trace ID length (localizations)",
    tid_duration = "Trace ID duration (ms)",
    toupper(param)
  )
  ylab <- if (param %in% c("tid_length", "tid_duration")) "Number of traces" else "Count"

  fig <- plot_ly()

  if (show_raw && length(raw_vec) > 0) {
    fig <- fig %>% add_histogram(
      x        = raw_vec,
      nbinsx   = n_bins,
      bingroup = "x",
      name     = "Raw data",
      opacity  = 0.35,
      marker   = list(color = "lightgray")
    )
  }
  if (length(filt_vec) > 0) {
    fig <- fig %>% add_histogram(
      x        = filt_vec,
      nbinsx   = n_bins,
      bingroup = "x",
      name     = "Filtered data",
      opacity  = 0.85,
      marker   = list(color = "steelblue")
    )
  }
  fig %>% layout(
    barmode = "overlay",
    xaxis   = list(title = xlab),
    yaxis   = list(title = ylab),
    legend  = list(orientation = "h", x = 0, y = 1.1)
  )
})
