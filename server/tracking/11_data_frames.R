## ─────────────────────────────────────────────────────────────────────────────
## server/tracking/11_data_frames.R  -  trk(): THE analysis-ready dataset; playback frames
## ─────────────────────────────────────────────────────────────────────────────
## Part of tracking_server(). Not a standalone script: its statements are evaluated in
## the server function's environment, in file-name order (server/parts.R).
##
## trk() (what every analysis reads), the chunk scope of playback, seq_data
## and frames_store.
##
## INTERFACE (generated: Rscript tools/server_map.R --headers; do not edit by hand)
##   DEFINES for other parts:
##     frames_store, seq_data, trk
##   READS from other parts:
##     01_state_sigma: tick
##     03_playback_state_zcal: is_simul
##     06_filter_pipeline: trk_gap_split, trk_spatial_filtered
##     07_cycle_time_sequence: trk_cycle_time
##     10_triggers: cluster_color_map_stored, conf_color_map_stored,
##                  data_ver, layout_trigger, plot_trigger, use_3d
## END INTERFACE
## ─────────────────────────────────────────────────────────────────────────────

# Filter changes do NOT live-update the track plot (the old filter_tick auto-replot
# was removed). The metadata / length / duration filters - like the split / trim /
# quality-gate cleaning stages - apply to the plotted data ONLY when the user presses
# Plot Data (fire_layout → bumps data_ver → reruns trk()) or Filter & Rechunk. Live,
# filter changes still drive the metadata HISTOGRAM (reads the debounced filtered
# stages directly) and the length/duration range read-outs (re-keyed on data_ver, so
# they describe the plotted data) - but the WebGL track render stays put until an
# explicit Plot Data / Rechunk. See observeEvent(input$trk_plot_data) above.


# ---- Prepare data once per filter-data change ----
# trk() is keyed on data_ver, NOT plot_trigger.  This means it skips
# prep_tracking + tid_color_map for rebuild_tick-only (playback param) triggers
# where the filtered dataframe is identical to the previous run.
# frames_store() still runs on every plot_trigger; it reads isolate(trk()) and
# gets the cached result when data_ver did not change.


# ── trk() - THE analysis-ready dataset ──────────────────────────────────────
# GOAL: the single hand-off point between the filter pipeline and everything
# that analyses or draws the data. Almost every bundle below starts with
# `df <- trk()`. If a filter change is not showing up in an analysis, this is
# the first place to look.
#
# It does three things:
#   1. takes the fully-filtered data (last pipeline stage, trk_spatial_filtered)
#   2. canonicalises it via prep_tracking() - sorted by (tid, time), non-finite
#      rows dropped, per-track `tim0` clock added
#   3. assigns each track a stable colour, with split segments of the same
#      parent track sharing a hue family
#
# WHY eventReactive(data_ver()): recomputing on every upstream twitch would be
# far too expensive. `data_ver()` is a version counter bumped only when the
# data genuinely changes, so trk() recomputes once per real change. Everything
# inside is read with isolate() so only data_ver() triggers it.
trk <- eventReactive(data_ver(), {
  df <- isolate(trk_spatial_filtered()); req(nrow(df) > 0)   # last pipeline stage
  df$tid <- as.character(df$tid)                             # TIDs are character everywhere
  # normalize_per_tid: in simultaneous playback every track must start at t=0
  # so they animate together; in sequential mode real timestamps are kept.
  df <- prep_tracking(df, normalize_per_tid = isolate(isTRUE(is_simul())))
  cmap <- tid_color_map(df$tid)                              # stable colour per TID

  # For split-track segments (TIDs matching "<parent>_s<N>"), derive a darker
  # shade of the first segment's colour so sibling segments share a hue family.
  # _s1 keeps the base colour; _s2 is 80% brightness, _s3 is 65%, etc.
  split_ids <- grep("_s\\d+$", names(cmap), value = TRUE, perl = TRUE)
  if (length(split_ids)) {
    parents <- unique(sub("_s\\d+$", "", split_ids, perl = TRUE))
    for (par in parents) {
      segs <- sort(as.integer(sub(paste0("^", par, "_s"), "", grep(
        paste0("^", par, "_s\\d+$"), names(cmap), value = TRUE, perl = TRUE
      ), perl = TRUE)))
      if (!length(segs)) next
      base_nm  <- paste0(par, "_s", segs[1])
      base_hex <- cmap[[base_nm]]
      for (sn in segs[-1]) {
        nm_sn <- paste0(par, "_s", sn)
        if (!nm_sn %in% names(cmap)) next
        rgb_m <- grDevices::col2rgb(base_hex) / 255
        hsv_m <- grDevices::rgb2hsv(rgb_m[1], rgb_m[2], rgb_m[3])
        v_new <- max(0.55, hsv_m["v", 1] * (1 - 0.09 * (sn - 1L)))
        cmap[[nm_sn]] <- grDevices::hsv(hsv_m["h", 1], hsv_m["s", 1], v_new)
      }
    }
  }

  df$col <- unname(cmap[as.character(df$tid)])
  attr(df, "col_map") <- cmap

  # Per-localization dark-cycle count + colour (df is sorted by tid,tim here).
  # dt uses RAW tim (ms); normalization only offsets per-TID so diffs are real.
  t_c <- tryCatch(isolate(trk_cycle_time()), error = function(e) NA_real_)
  n2   <- nrow(df)
  same <- if (n2 > 1) c(FALSE, df$tid[-1] == df$tid[-n2]) else FALSE
  dt   <- if (n2 > 1) c(NA_real_, diff(df$tim)) else NA_real_
  dt[!same] <- NA_real_
  if (is.finite(t_c) && t_c > 0) {
    ndk <- pmax(0L, as.integer(round(dt / t_c)) - 1L)   # k-1 dark cycles; k=round(dt/t_c)
  } else ndk <- rep(0L, n2)
  ndk[is.na(ndk)] <- 0L                                 # first loc of each TID → integrated
  df$n_dark   <- ndk
  df$col_dark <- dark_cycle_cols(ndk)
  # (Z is already corrected at the pipeline root in trk_raw_df() by the physical
  # Z factor, so it flows in here through the filter chain - not re-applied.)
  df
}, ignoreNULL = TRUE)

# Sequential data (always order by raw 'tim'; decimate for speed)
# ── Chunk-scoped PLAYBACK ───────────────────────────────────────────────────
# Chunking exists only to keep the animation smooth: the filter pipeline,
# analysis, and the full-track static overlay all use ALL TIDs, but the moving
# markers/trails are limited to the selected chunk's TIDs. Returns the TIDs exactly
# as the chunk file stores them - original ("42") for upload chunks, segment
# ("42_s1") once a rechunk has baked the split in; scope_playback() below accepts
# either. NULL = animate everything ("All TIDs" selected).
trk_playback_chunk_tids <- reactive({
  sel <- input$trk_chunk_sel
  if (!isTruthy(sel) || identical(sel, "__all__") || !file.exists(sel)) return(NULL)
  ch <- read_chunk_file(sel, columns = "tid")
  if (is.null(ch) || !"tid" %in% names(ch) || !nrow(ch)) return(NULL)
  unique(as.character(ch$tid))
})


# ── seq_data() - prepare tracks for PLAYBACK ────────────────────────────────
# GOAL: turn the analysis dataset into a per-track structure the animation can
# step through cheaply. Playback must stay smooth, so this does the expensive
# work once (scope to chunk, thin, split by track) instead of per frame.
#
# Two deliberate cost controls:
#   * only the SELECTED playback chunk is animated (analyses still use all TIDs)
#   * each track is decimated to at most `trk_max_pts_tid` points
# Both are DISPLAY-only; no estimator ever sees this thinned data.
#
# Recomputes only when plot_trigger() fires, not on every input change.
seq_data <- eventReactive(plot_trigger(), {
  df <- isolate(trk()); req(nrow(df) > 0)
  df <- scope_playback(df, isolate(trk_playback_chunk_tids()))   # animate selected chunk only
  # decimate on raw 'tim' to keep timing shape but lighten the load
  df_d <- decimate_per_tid(df, max_pts = isolate(input$trk_max_pts_tid), tcol = "tim")
  df_d <- df_d[order(df_d$tid, df_d$t), ]  # 't' is raw tim here
  sp   <- split(df_d, df_d$tid)
  tids <- names(sp)
  # play TIDs in order of earliest time
  ord  <- order(sapply(sp, function(tb) min(tb$t, na.rm = TRUE)))
  list(sp = sp[tids[ord]], tids = tids[ord])
}, ignoreNULL = TRUE)
frames_store <- eventReactive(plot_trigger(), {
  df <- isolate(trk()); req(nrow(df) > 0)

  # Drift-speed colouring is computed HERE, not in trk(): trk() is keyed on data_ver
  # only, so it would not refresh when the colour mode or window changes. frames_store
  # re-runs on plot_trigger, which fire_layout("display_opts") bumps. Only computed
  # when the speed mode is active, since it is an O(N) pass.
  speed_max <- NA_real_
  if (identical(isolate(input$trk_color_mode %||% "tid"), "speed")) {
    w_ms <- suppressWarnings(as.numeric(isolate(input$trk_color_speed_win %||% 200)))
    if (!is.finite(w_ms) || w_ms <= 0) w_ms <- 200
    df$v_drift <- trk_point_drift_speed(df, win_ms = w_ms)
    # Robust upper limit so a few outlier windows don't flatten the whole scale.
    speed_max <- suppressWarnings(as.numeric(stats::quantile(df$v_drift, 0.95, na.rm = TRUE)))
    if (!is.finite(speed_max) || speed_max <= 0)
      speed_max <- suppressWarnings(max(df$v_drift, na.rm = TRUE))
    if (!is.finite(speed_max) || speed_max <= 0) speed_max <- 1
    df$col_speed <- speed_cols(df$v_drift, speed_max)
  }

  # ── 2D viewing plane (XY / XZ / YZ) ────────────────────────────────────────
  # Bake the plane choice in as a COORDINATE REMAP, so the whole 2D pipeline
  # (full-tracks + playback frames + axes) plots the chosen plane without any
  # per-branch z handling - the same trick the SMLM Display tab uses. Only in 2D
  # mode and only when z is real; 3D mode always keeps native XY(Z). Placed AFTER
  # the drift-speed pass so v_drift stays a true in-plane physical speed.
  plane <- if (!isTRUE(isolate(use_3d()))) (isolate(input$trk_view_plane) %||% "XY") else "XY"
  # Fall back to XY unless z is present AND actually varies (defensive: mirrors the
  # plane-selector gate, so a stale XZ/YZ selection on planar data can't flat-line it).
  z_fin <- if ("z" %in% names(df)) df$z[is.finite(df$z)] else numeric(0)
  if (plane %in% c("XZ", "YZ") && !(length(z_fin) > 0 && diff(range(z_fin)) > 1e-6)) plane <- "XY"
  remap_plane <- function(d) {
    if (is.null(d) || !("z" %in% names(d))) return(d)
    if (plane == "XZ")      d$y <- d$z              # X horizontal, Z vertical
    else if (plane == "YZ") { d$x <- d$y; d$y <- d$z }  # Y horizontal, Z vertical
    d
  }
  df <- remap_plane(df)

  # Optional faint background of ALL raw (pre-filter) tracks - shown in BOTH full
  # and playback modes so you can see what the filters removed. Derived from the
  # POST-SPLIT raw object (trk_gap_split): when splitting is on the faint tracks
  # are split too (so they line up with the analyzed segments); when it's off,
  # gap_split returns the raw data unchanged. Computed once here (decimated), only
  # when the toggle is on.
  sp_raw <- if (isTRUE(isolate(input$trk_full_raw))) {
    mpb <- suppressWarnings(as.integer(isolate(input$trk_max_pts_tid)))
    if (!is.finite(mpb) || mpb < 50) mpb <- 400L
    rdf <- tryCatch(isolate(trk_gap_split()), error = function(e) NULL)
    if (!is.null(rdf) && nrow(rdf) && all(c("tid","tim","x","y") %in% names(rdf))) {
      rd <- decimate_per_tid(rdf, max_pts = mpb, tcol = "tim")
      rd <- remap_plane(rd)   # keep the faint raw tracks in the same plane as the data
      split(rd, rd$tid)
    } else NULL
  } else NULL

  # ── FULL-TRACKS MODE: skip ALL playback frame/bin/trail machinery ──────────
  # The static render (view_mode_main == "full") reads only sp_static / col_map /
  # tids - never the per-frame data. Building frames here is what made length,
  # duration and split edits feel slow in the default view, so short-circuit:
  # compute just the full-track overlay + colours and return one empty placeholder
  # frame (keeps downstream req(length(frames) > 0) happy). Frames are built only
  # when the user switches to Playback, which fires fire_layout() → full rebuild.
  if (identical(isolate(input$trk_view_mode_main), "full")) {
    use3d_f <- isolate(isTRUE(use_3d()))
    max_pts <- suppressWarnings(as.integer(isolate(input$trk_max_pts_tid)))
    if (!is.finite(max_pts) || max_pts < 50) max_pts <- 400L

    # Progress feedback: building the per-track overlay (decimate + split into one
    # trace per track) is the visible cost of each Plot Data / upload-settle render.
    # Full tracks normally draws every filtered TID - the chunk is a playback
    # cost control. When the user asks for it, scope this view to the selected
    # chunk as well, using the same TIDs the animation would show. Analysis and
    # filtering are untouched either way: they always read all TIDs.
    df_full <- if (isTRUE(isolate(input$trk_chunk_scope_full)))
                 scope_playback(df, isolate(trk_playback_chunk_tids())) else df

    full_build <- withProgress(message = "Building track view…", value = 0.2, {
      df_static <- decimate_per_tid(df_full, max_pts = max_pts, tcol = "tim",
                                    extra_cols = c("col_dark", "n_dark",
                                                   "col_speed", "v_drift"))
      setProgress(0.5, detail = sprintf("%s tracks", format(length(unique(df_static$tid)), big.mark = ",")))
      sp_static <- split(df_static, df_static$tid)
      # (Z is already corrected in trk() by the physical Z factor - no display-only
      # stretch here, otherwise it would be applied twice.)
      tids_all     <- names(sp_static)
      base_col_map <- if (!is.null(attr(df, "col_map"))) attr(df, "col_map") else tid_color_map(tids_all)
      cluster_cm   <- isolate(cluster_color_map_stored())
      conf_cm      <- isolate(conf_color_map_stored())
      override_cm  <- if (!is.null(cluster_cm)) cluster_cm else conf_cm
      if (!is.null(override_cm))
        for (nm in names(override_cm)) if (nm %in% names(base_col_map)) base_col_map[[nm]] <- override_cm[[nm]]
      setProgress(1)
      list(sp_static = sp_static, tids_all = tids_all,
           base_col_map = base_col_map, sp_raw = sp_raw)
    })
    sp_static    <- full_build$sp_static
    tids_all     <- full_build$tids_all
    base_col_map <- full_build$base_col_map
    sp_raw       <- full_build$sp_raw

    return(list(
      static_lines = NULL,
      frames       = list(f1 = list(
        trails_per_tid = list(), trails_x = numeric(0), trails_y = numeric(0),
        markers_x = numeric(0), markers_y = numeric(0), markers_z = NULL,
        marker_cols = character(0), step_dt = NA_real_, bin_width = NA_real_
      )),
      t_grid     = NA_real_,
      col_map    = base_col_map,
      tids       = tids_all,
      sp         = sp_static,    # axis fallback / faint overlay read this
      sp_static  = sp_static,
      sp_raw     = sp_raw,       # optional faint background of ALL raw tracks
      trail_len  = 0, trail_nmax = 1L, bin_width = NA_real_,
      using_bins = FALSE, step_dt = NA_real_, has_3d = use3d_f,
      plane      = plane,        # 2D viewing plane the coordinates were remapped to
      speed_max  = speed_max     # drift-speed colour-scale limit (NA unless speed mode)
    ))
  }
  # ── PLAYBACK MODE: full windowing + per-frame binning below ────────────────

  # df         = ALL TIDs → used for the static full-track overlay (sp_static)
  # df_play    = selected chunk only → used for the animated markers/trails
  df_play <- scope_playback(df, isolate(trk_playback_chunk_tids()))
  use3d <- isolate(isTRUE(use_3d()))

  # Playback time basis
  tcol_play <- if (isTRUE(isolate(is_simul()))) "tim0" else "tim"

  # Fixed total playback window (in the chosen basis)
  total <- suppressWarnings(as.numeric(isolate(input$trk_play_total)))
  if (!is.finite(total) || total <= 0) total <- 1000
  
  # Determine [tmin, tmax] in the playback basis (over the animated chunk)
  rng_play_all <- range(df_play[[tcol_play]], na.rm = TRUE)
  if (isTRUE(isolate(is_simul()))) {
    tmin <- 0
    tmax <- min(tmin + total, max(rng_play_all, na.rm = TRUE))
  } else {
    tmin <- min(rng_play_all, na.rm = TRUE)
    tmax <- min(tmin + total, max(rng_play_all, na.rm = TRUE))
  }
  if (!is.finite(tmin) || !is.finite(tmax) || tmax <= tmin) {
    tmin <- rng_play_all[1]; tmax <- rng_play_all[2]
  }
  
  # Filter the ANIMATED data (chunk-scoped) to the fixed window
  dwin <- df_play[df_play[[tcol_play]] >= tmin & df_play[[tcol_play]] <= tmax, , drop = FALSE]
  req(nrow(dwin) > 0)
  
  # Light pre-decimation for TRAILS only
  max_pts_tid  <- suppressWarnings(as.integer(isolate(input$trk_max_pts_tid)))
  if (!is.finite(max_pts_tid) || max_pts_tid < 50) max_pts_tid <- 400
  # Carry the per-point colour columns so trails can be coloured by the active mode.
  dwin_tr_dec  <- decimate_per_tid(dwin, max_pts = max_pts_tid, tcol = tcol_play,
                                   extra_cols = c("col_dark", "col_speed"))
  dwin_tr_dec  <- dwin_tr_dec[order(dwin_tr_dec$tid, dwin_tr_dec$t), ]
  
  # Split full-res for markers, and decimated for trails
  sp_full <- split(dwin[order(dwin$tid, dwin[[tcol_play]]), ], dwin$tid)  # FULL (markers)
  sp_tr   <- split(dwin_tr_dec, dwin_tr_dec$tid)                          # DECIMATED (trails)

  # Full-range (unwindowed) spatial data for the static full-tracks trace.
  # Using df (not dwin) so tracks outside the playback window still appear.
  df_static    <- decimate_per_tid(df, max_pts = max_pts_tid, tcol = "tim")
  sp_static    <- split(df_static, df_static$tid)

  # (Z is already corrected in trk() by the physical Z factor, so the per-TID
  # spatial lists inherit it - no separate display stretch is applied here.)

  # Colors: start with random per-TID colors, then overlay confinement or cluster map.
  # Priority: cluster > confinement > base random.
  tids         <- names(sp_full)
  base_col_map <- if (!is.null(attr(df, "col_map"))) attr(df, "col_map") else tid_color_map(tids)
  cluster_cm  <- isolate(cluster_color_map_stored())
  conf_cm     <- isolate(conf_color_map_stored())
  override_cm <- if (!is.null(cluster_cm)) cluster_cm else conf_cm
  col_map <- if (!is.null(override_cm)) {
    cm <- base_col_map
    for (nm in names(override_cm)) if (nm %in% names(cm)) cm[[nm]] <- override_cm[[nm]]
    cm
  } else base_col_map

  # Per-point colour column for the active mode (NULL = one colour per TID).
  pcol <- switch(isolate(input$trk_color_mode %||% "tid"),
                 dark = "col_dark", speed = "col_speed", NULL)
  color_by_dark <- !is.null(pcol)   # any per-point mode: markers carry the colour
  # Display filter by integration status: "all" | "integrated" (n_dark==0) | "dark" (n_dark>=1)
  dark_show <- isolate(input$trk_dark_show %||% "all")

  # Frame grid from bin width OR #frames over [tmin, tmax]
  bin_dt <- suppressWarnings(as.numeric(isolate(input$trk_bin_dt)))
  if (!is.finite(bin_dt) || bin_dt <= 0) bin_dt <- (tmax - tmin) / 120  # fallback
  edges <- seq(tmin, tmax, by = bin_dt)
  if (length(edges) < 2) edges <- c(tmin, tmax)
  if (tail(edges, 1) < (tmax - 1e-12)) edges <- c(edges, tmax)
  
  n_frames <- length(edges) - 1L
  bin_width <- if (length(edges) >= 2) (edges[2] - edges[1]) else NA_real_
  req(n_frames >= 1L)
  
  # Quota params
  target_markers <- suppressWarnings(as.integer(isolate(input$trk_target_markers)))
  if (!is.finite(target_markers) || target_markers < 0) target_markers <- 0L
  per_tid_cap    <- suppressWarnings(as.integer(isolate(input$trk_per_tid_cap)))
  if (!is.finite(per_tid_cap)    || per_tid_cap    < 0) per_tid_cap    <- 0L
  
  # Trail params
  trail_len  <- suppressWarnings(as.numeric(isolate(input$trk_trail_len)));  if (!is.finite(trail_len)  || trail_len  < 0) trail_len  <- 0
  trail_nmax <- suppressWarnings(as.integer(isolate(input$trk_trail_nmax))); if (!is.finite(trail_nmax) || trail_nmax < 1) trail_nmax <- 1L
  
  # ----- Precompute frame membership for markers (O(N_total) once) -----
  # For each TID, split its row-indices by bin index so the per-frame loop
  # does a cheap list-lookup instead of scanning all points every frame.
  # findInterval returns 1..n_frames for points inside [tmin,tmax]; 0 or
  # n_frames+1 for out-of-window points (excluded via the 'valid' mask).
  frame_rows_per_tid <- lapply(tids, function(tid) {
    tb  <- sp_full[[tid]]
    bi  <- findInterval(tb[[tcol_play]], edges)
    valid <- bi >= 1L & bi <= n_frames
    # Integration-status display filter (markers only; trails keep the full path).
    if (dark_show != "all" && "n_dark" %in% names(tb)) {
      valid <- valid & (if (dark_show == "integrated") tb$n_dark == 0L else tb$n_dark >= 1L)
    }
    if (!any(valid)) return(list())
    split(which(valid), bi[valid])   # named list: "k" -> integer row-indices
  })

  # ----- Precompute sorted trail-time vectors for O(log n) window lookups -----
  sp_tr_t <- lapply(tids, function(tid) {
    tb <- sp_tr[[tid]]
    if (is.null(tb) || !nrow(tb)) numeric(0) else tb$t
  })

  # Build frames (combined trail + markers)
  frames <- vector("list", n_frames)
  names(frames) <- paste0("f", seq_len(n_frames))

  for (k in seq_len(n_frames)) {
    e0    <- edges[k]
    k_str <- as.character(k)

    # --- Markers (FULL), with quotas per bin ---
    # Fast O(1) lookup from precomputed membership instead of O(n_pts) which()
    per_tid_idx <- lapply(frame_rows_per_tid, function(fr) fr[[k_str]] %||% integer(0))
    bin_counts  <- lengths(per_tid_idx)

    total_in_bin <- sum(bin_counts)
    if (target_markers > 0 && total_in_bin > target_markers) {
      # proportional allocation + largest remainder
      prop   <- bin_counts / total_in_bin
      raw_q  <- prop * target_markers
      q      <- pmin(bin_counts, floor(raw_q))
      if (per_tid_cap > 0) q <- pmin(q, per_tid_cap)
      left   <- max(0L, as.integer(round(target_markers - sum(q))))
      if (left > 0L) {
        gain <- raw_q - q
        add_order <- order(gain, decreasing = TRUE)
        j <- 1L
        while (left > 0L && length(add_order) > 0L) {
          i <- add_order[j]
          if (q[i] < bin_counts[i] && (per_tid_cap == 0L || q[i] < per_tid_cap)) {
            q[i] <- q[i] + 1L; left <- left - 1L
          }
          j <- if (j < length(add_order)) j + 1L else 1L
        }
      }
      quotas <- q
    } else {
      quotas <- bin_counts
      if (per_tid_cap > 0) quotas <- pmin(quotas, per_tid_cap)
    }

    mx <- my <- numeric(0); mz <- numeric(0); mcols <- character(0)

    for (i in seq_along(tids)) {
      w <- per_tid_idx[[i]]
      if (!length(w)) next
      kq <- quotas[i]
      if (kq >= length(w)) {
        keep <- w
      } else if (kq > 0) {
        keep <- w[unique(round(seq(1, length(w), length.out = kq)))]
      } else {
        keep <- integer(0)
      }
      if (length(keep)) {
        tb <- sp_full[[tids[i]]]
        mx    <- c(mx, tb$x[keep])
        my    <- c(my, tb$y[keep])
        if (use3d) mz <- c(mz, tb$z[keep])
        if (!is.null(pcol) && pcol %in% names(tb)) {
          mcols <- c(mcols, as.character(tb[[pcol]][keep]))
        } else {
          mcols <- c(mcols, rep(unname(as.character(col_map[[tids[i]]])), length(keep)))
        }
      }
    }

    # --- Per-TID trail data: binary search via findInterval (O(log n) per TID) ---
    trails_per_tid <- vector("list", length(tids))
    tx <- ty <- numeric(0)  # combined, used only for dedup check
    if (trail_len > 0) {
      for (i in seq_along(tids)) {
        tt <- sp_tr_t[[i]]
        if (!length(tt)) {
          trails_per_tid[[i]] <- list(x = numeric(0), y = numeric(0))
          next
        }
        # Binary search: points with t in (e0 - trail_len, e0]
        lo_idx <- findInterval(e0 - trail_len, tt) + 1L
        hi_idx <- findInterval(e0, tt, rightmost.closed = TRUE)
        if (hi_idx >= lo_idx) {
          wtr <- seq.int(lo_idx, hi_idx)
          if (length(wtr) > trail_nmax) wtr <- tail(wtr, trail_nmax)
          tb_tr <- sp_tr[[tids[i]]]
          tr_entry <- list(x = tb_tr$x[wtr], y = tb_tr$y[wtr])
          if (use3d && "z" %in% names(tb_tr)) tr_entry$z <- tb_tr$z[wtr]
          # Per-point colour mode: colour the whole trail by its LEADING (most recent)
          # point. plotly 2D lines take a single colour per trace - per-vertex line
          # colour is not supported - but a trail only spans trail_len, so its leading
          # value is "the value right now" and the trail recolours as playback advances.
          if (!is.null(pcol) && pcol %in% names(tb_tr))
            tr_entry$col <- as.character(tb_tr[[pcol]][wtr[length(wtr)]])
          trails_per_tid[[i]] <- tr_entry
          tx <- c(tx, tb_tr$x[wtr], NA_real_)
          ty <- c(ty, tb_tr$y[wtr], NA_real_)
        } else {
          trails_per_tid[[i]] <- list(x = numeric(0), y = numeric(0))
        }
      }
    } else {
      for (i in seq_along(tids))
        trails_per_tid[[i]] <- list(x = numeric(0), y = numeric(0))
    }

    frames[[k]] <- list(
      trails_per_tid = trails_per_tid,
      trails_x       = tx,   # combined, used only for dedup
      trails_y       = ty,
      markers_x      = mx,
      markers_y      = my,
      markers_z      = if (use3d) mz else NULL,
      marker_cols    = mcols,
      step_dt        = if (length(edges) >= 2) (edges[2] - edges[1]) else NA_real_,
      bin_width      = bin_width
    )
  }
  
  # --- Collapse identical consecutive frames (same markers + same trails) ---
  same_vec <- function(a, b, tol = 1e-12) {
    length(a) == length(b) && (length(a) == 0 || max(abs(a - b)) <= tol)
  }
  keep <- logical(length(frames)); prev <- NULL
  for (k in seq_along(frames)) {
    fr <- frames[[k]]
    if (is.null(prev)) {
      keep[k] <- TRUE
    } else {
      dup <- same_vec(fr$markers_x, prev$markers_x) &&
        same_vec(fr$markers_y, prev$markers_y) &&
        same_vec(fr$trails_x,  prev$trails_x)  &&
        same_vec(fr$trails_y,  prev$trails_y)
      keep[k] <- !dup
    }
    prev <- fr
  }
  t_grid_full <- edges[-length(edges)]  # start-of-bin time for each frame
  frames  <- frames[keep]
  t_grid_full <- t_grid_full[keep]     # match filtered frames
  
  list(
    static_lines = NULL,
    frames      = frames,
    t_grid      = t_grid_full,
    col_map     = col_map,
    tids        = tids,
    sp          = sp_tr,      # trail animation uses windowed decimated data
    sp_static   = sp_static,  # full-tracks overlay uses unwindowed data
    sp_raw      = sp_raw,     # optional faint background of ALL raw tracks
    trail_len   = trail_len,
    trail_nmax  = trail_nmax,
    bin_width   = bin_width,
    using_bins  = TRUE,
    step_dt     = if (length(edges) >= 2) (edges[2] - edges[1]) else NA_real_,
    has_3d      = use3d,     # TRUE when z column is present and 3D mode is active
    plane       = plane,     # 2D viewing plane the coordinates were remapped to
    color_by_dark = color_by_dark,  # per-point colour mode: markers carry it (overlaid even in trails mode)
    speed_max   = speed_max         # robust upper limit for the drift-speed colour scale
  )
}, ignoreNULL = TRUE)

# --- show frame dt and time basis ----------------------------------------
output$trk_frame_meta <- renderUI({
  req(!is.na(plot_trigger()))  # plot_trigger: updates after auto-tune frame rebuilds too
  pre <- frames_store(); req(is.list(pre), length(pre$frames) > 0)
  tcol <- if (isTRUE(is_simul())) "tim0" else "tim"
  dt_show <- pre$step_dt
  if (!is.finite(dt_show) || dt_show <= 0) dt_show <- pre$bin_width
  # Live frame + marker count (moved off the plot title: updating a title needs a
  # relayout, which re-solves scaleanchor and rescales the FOV every frame).
  k_now <- suppressWarnings(as.integer(tick()))
  cur   <- if (is.finite(k_now) && k_now >= 1L && k_now <= length(pre$frames))
             pre$frames[[k_now]] else NULL
  htmltools::div(
    style = "color:#555;margin-top:6px;margin-bottom:6px;",
    htmltools::tags$b("Time basis: "), tcol, htmltools::HTML("&nbsp;&middot;&nbsp;"),
    htmltools::tags$b("Frames: "), length(pre$frames), htmltools::HTML("&nbsp;&middot;&nbsp;"),
    htmltools::tags$b("Frame dt ~ "), fmt_num(dt_show), " ms",
    # Own line + nowrap so changing digit counts can't reflow the sidebar mid-playback.
    if (!is.null(cur)) htmltools::div(
      style = "white-space:nowrap;",
      htmltools::tags$b("Frame: "), k_now, " / ", length(pre$frames),
      htmltools::HTML("&nbsp;&middot;&nbsp;"),
      htmltools::tags$b("Markers: "), length(cur$markers_x)
    )
  )
})

outputOptions(output, "trk_frame_meta", suspendWhenHidden = FALSE)

# Wrap plotlyOutput in a renderUI so that every fire_layout() destroys and
# recreates the DOM element. This forces Plotly to call newPlot() (not react()),
# which properly initialises the WebGL context on chunk switches and structural
# changes. Without this, Plotly's in-place react() path blanks scattergl charts.
output$trk_track2D_wrapper <- renderUI({
  req(!is.na(layout_trigger()))
  plotlyOutput("trk_track2D", height = "1000px")
})
