## ─────────────────────────────────────────────────────────────────────────────
## trk_playback.R  -  tracking: playback timing, colour scales and the dwell heatmap
## ─────────────────────────────────────────────────────────────────────────────
##
## PURPOSE
##   Display helpers of the track viewer: the timer interval of the play loop,
##   dark-cycle and speed colour scales, the colour-bar trace, default playback
##   settings for a dataset, frame scoping, and the dwell-time heatmap.
##
##   Plain functions and constants with no Shiny state, moved unchanged out of
##   the tracking_server() closure (server_tracking.R). The server gathers the
##   inputs and calls these.
##
## UNITS: nm and ms; D in nm^2/ms unless a name says um2_s.
## ─────────────────────────────────────────────────────────────────────────────

# Base animation cadence (~16 fps at 1x speed). Replaces the old "Tick interval"
# input: it drives the simultaneous (prebinned) loop and floors the sequential
# loop. The "Playback speed (×)" slider scales the per-frame delay around this.
TICK_INTERVAL_SEC <- 0.06

# Hard lower bound on per-frame delay (~50 fps). Stops the loop from outpacing
# the websocket (which is what made playback "refuse to stop"), no matter how
# high the speed multiplier or how small "ms per time unit" gets.
MIN_INTERVAL_SEC  <- 0.02

# Effective per-frame delay (seconds) for the current speed multiplier.
play_interval_for <- function(speed) {
  s <- suppressWarnings(as.numeric(speed))
  if (!is.finite(s) || s <= 0) s <- 1
  max(MIN_INTERVAL_SEC, TICK_INTERVAL_SEC / s)
}

# ── Dark-cycle colouring ─────────────────────────────────────────────────────
# A step spanning k = round(dt / cycle_time) acquisition cycles has k-1 DARK cycles
# (k=1 = integrated, i.e. a localization every cycle). Each localization inherits the
# dark-cycle count of the step LEADING INTO it. Colour scheme: integrated = blue,
# 1+ dark cycles = light→deep red gradient (capped at DARK_NMAX; higher = deeper red).
DARK_NMAX <- 6L

.dark_hex <- function(nmax = DARK_NMAX) {
  reds <- grDevices::colorRampPalette(c("#fca9a4", "#e34a33", "#7f000d"))(nmax)
  c("#2c7fb8", reds)                       # length nmax+1; index = n_dark + 1
}

dark_cycle_cols <- function(n_dark, nmax = DARK_NMAX) {
  pal <- .dark_hex(nmax)
  i <- pmin(pmax(as.integer(n_dark), 0L), nmax) + 1L
  i[!is.finite(i)] <- 1L
  pal[i]
}

dark_cycle_colorscale <- function(nmax = DARK_NMAX) {
  pal <- .dark_hex(nmax)                   # colours for values 0..nmax
  lapply(0:nmax, function(k) list(k / nmax, pal[k + 1]))
}

# ── Speed colouring (drift velocity) ─────────────────────────────────────────
# Per-step speed |Δr|/Δt is NOT usable here: Brownian motion has no instantaneous
# velocity, so r/Δt diverges as Δt→0 and mostly measures the sampling rate (see the
# header note, and why the Speed panel's path speed can't separate diffusion from
# transport). We therefore colour by the same COHERENT-transport measure the Speed &
# Transport panel uses: v_drift = |net displacement| / span over a CENTRED window W.
# Diffusion's net displacement grows only as sqrt(t) so its apparent v ~ sqrt(4D/W)
# stays low, while sustained directed motion reads high - which is what makes a
# confined stretch read blue and a directed excursion read red.
# Vectorised per track via findInterval on the window edges → O(N) overall.
trk_point_drift_speed <- function(df, win_ms = 200, min_pts = 3L) {
  n_all <- nrow(df)
  v <- rep(NA_real_, n_all)
  if (!is.finite(win_ms) || win_ms <= 0 || n_all < 2L) return(v)
  has_z <- "z" %in% names(df) && any(is.finite(df$z))
  hw <- win_ms / 2
  for (ii in split(seq_len(n_all), df$tid)) {
    n <- length(ii); if (n < 2L) next
    t <- df$tim[ii]; xx <- df$x[ii]; yy <- df$y[ii]
    zz <- if (has_z) df$z[ii] else NULL
    lo <- pmin(pmax(findInterval(t - hw, t) + 1L, 1L), n)   # first pt >= t-W/2
    hi <- pmin(pmax(findInterval(t + hw, t),      1L), n)   # last  pt <= t+W/2
    span <- t[hi] - t[lo]
    dr <- sqrt((xx[hi] - xx[lo])^2 + (yy[hi] - yy[lo])^2 +
               (if (has_z) (zz[hi] - zz[lo])^2 else 0))
    okw <- is.finite(span) & span > 0 & is.finite(dr) & (hi - lo + 1L) >= min_pts
    vv <- rep(NA_real_, n); vv[okw] <- dr[okw] / span[okw]   # nm/ms = µm/s
    v[ii] <- vv
  }
  v
}

# Blue (slow / confined) → red (fast / directed). NA (window can't be formed, e.g.
# at track ends) stays GREY rather than being painted a misleading "slow" blue.
.speed_hex <- function(n = 64L)
  grDevices::colorRampPalette(c("#2c7fb8", "#41b6c4", "#fecc5c", "#fd8d3c", "#e31a1c"))(n)

speed_cols <- function(v, vmax) {
  pal <- .speed_hex()
  out <- rep("#9e9e9e", length(v))
  ok  <- is.finite(v) & isTRUE(is.finite(vmax) && vmax > 0)
  if (any(ok)) {
    i <- 1L + as.integer(round(pmin(pmax(v[ok] / vmax, 0), 1) * (length(pal) - 1L)))
    out[ok] <- pal[i]
  }
  out
}

speed_colorscale <- function() {
  pal <- .speed_hex(9L)
  lapply(seq_along(pal), function(k) list((k - 1) / (length(pal) - 1), pal[k]))
}

# Colorbar-only trace. Per-point colours are baked as hex strings (needed by the
# restyle-driven playback path), and hex arrays carry no colourbar - so we append one
# invisible 2-point trace whose sole job is to render the scale. Added LAST so it can
# never shift trail_start / markers_idx in trace_map.
# ax/ay/az anchor the trace on REAL coordinates: plotly (scattergl especially) drops a
# trace whose points are all NA, and a dropped trace draws no colourbar. The markers
# are made invisible with opacity 0 + size 0 instead.
add_colorbar_trace <- function(fig, mode, fg_col, vmax = NA_real_, use3d = FALSE,
                               ax = 0, ay = 0, az = 0) {
  if (!mode %in% c("dark", "speed")) return(fig)
  if (!is.finite(ax) || !is.finite(ay)) { ax <- 0; ay <- 0 }
  if (!is.finite(az)) az <- 0
  if (identical(mode, "dark")) {
    cs <- dark_cycle_colorscale(); cmin <- 0; cmax <- DARK_NMAX
    cb <- list(title = list(text = "dark<br>cycles", font = list(color = fg_col)),
               tickmode = "array", tickvals = 0:DARK_NMAX,
               ticktext = c(as.character(0:(DARK_NMAX - 1L)), paste0(DARK_NMAX, "+")),
               tickfont = list(color = fg_col), len = 0.5, thickness = 14)
  } else {
    if (!is.finite(vmax) || vmax <= 0) return(fig)
    cs <- speed_colorscale(); cmin <- 0; cmax <- vmax
    cb <- list(title = list(text = "drift speed<br>(µm/s)", font = list(color = fg_col)),
               tickfont = list(color = fg_col), len = 0.5, thickness = 14)
  }
  args <- list(x = c(ax, ax), y = c(ay, ay),
               type = if (use3d) "scatter3d" else "scattergl", mode = "markers",
               marker = list(size = 0, opacity = 0, color = c(cmin, cmax),
                             colorscale = cs, cmin = cmin, cmax = cmax,
                             showscale = TRUE, colorbar = cb),
               hoverinfo = "skip", showlegend = FALSE, name = "scale")
  if (use3d) args$z <- c(az, az)
  do.call(add_trace, c(list(fig), args))
}

# Subset a display data.frame to the chunk's TIDs.
#
# A chunk file holds one of TWO TID conventions, and the match has to accept both:
#   * UPLOAD chunks hold the original TIDs ("42"). Gap-splitting happens live
#     afterwards, so the display frame carries segments ("42_s1") that must still
#     be matched to their parent chunk - hence the base-TID comparison.
#   * RE-CHUNKED sets have the split BAKED IN (the rechunk writes "42_s1" as the
#     TID), so the display frame and the chunk file agree exactly and the base
#     comparison alone never matches - "42" is not in {"42_s1", ...}. Every row
#     then failed, the guard below returned the whole frame, and selecting a
#     chunk appeared to do nothing.
# Testing both forms covers each case without over-matching: a chunk holding
# "42_s1" still does not pull in "42_s2", which may belong to a different chunk.
#
# Falls back to the full frame if the scope would leave nothing to animate.
scope_playback <- function(d, chunk_tids) {
  if (is.null(chunk_tids) || !nrow(d) || !"tid" %in% names(d)) return(d)
  tid_d <- as.character(d$tid)
  base  <- sub("_s\\d+$", "", tid_d, perl = TRUE)
  keep  <- tid_d %in% chunk_tids | base %in% chunk_tids
  if (!any(keep)) return(d)            # chunk fully filtered out - don't blank playback
  d[keep, , drop = FALSE]
}

fmt_num <- function(x) ifelse(is.finite(x), format(signif(x, 4), trim = TRUE), "NA")

hex_to_rgba <- function(hex, alpha = 0.5) {
  rgb_vals <- grDevices::col2rgb(hex)
  sprintf("rgba(%d,%d,%d,%.2f)", rgb_vals[1L], rgb_vals[2L], rgb_vals[3L], alpha)
}

smart_defaults <- function(df, mode = c("simul","seq"), target_frames = 200, target_fps = 15) {
  mode <- match.arg(mode)
  tcol <- if (mode == "simul") "tim0" else "tim"
  
  # robust time span
  rng <- range(df[[tcol]], na.rm = TRUE)
  total_units <- if (mode == "simul") max(df$tim0, na.rm = TRUE) else diff(rng)
  if (!is.finite(total_units) || total_units <= 0) total_units <- 1
  
  # round up to a nice number for the playback window
  nice_total <- max(pretty(c(0, total_units), n = 5))
  
  # frame count / width
  n_frames <- target_frames
  bin_dt   <- nice_total / n_frames
  
  # trail ~5% of window, at least one frame wide
  trail_len <- max(nice_total * 0.05, bin_dt)
  
  # marker quota ~ avg points per frame (*1.5 safety), clipped
  n_pts <- nrow(df)
  target_markers <- max(200, min(8000, round((n_pts / n_frames) * 1.5)))
  
  # sequential timing: choose ms_per_unit to get ~target_fps
  # use median per-TID duration for stability
  durs <- df |>
    dplyr::group_by(tid) |>
    dplyr::summarise(dur = max(.data[[tcol]]) - min(.data[[tcol]]), .groups = "drop") |>
    dplyr::pull(dur)
  med_dur <- suppressWarnings(stats::median(durs[is.finite(durs) & durs > 0]))
  if (!is.finite(med_dur) || med_dur <= 0) med_dur <- total_units
  step_units   <- med_dur / n_frames           # units per frame (seq)
  frame_ms     <- 1000 / target_fps            # ~15 fps
  ms_per_unit  <- max(1, round(frame_ms / max(step_units, .Machine$double.eps)))
  
  list(
    nice_total      = nice_total,
    n_frames        = n_frames,
    bin_dt          = bin_dt,
    trail_len       = trail_len,
    target_markers  = target_markers,
    ms_per_unit     = ms_per_unit
  )
}  

#--------HEATMAP-------
hm_colorscale <- function(name = "fire") {
  switch(name,
    magenta = list(
      c(0.00, "#000000"), c(0.40, "#3b0f3f"),
      c(0.70, "#a100a1"), c(1.00, "#ffe6ff")
    ),
    fire = list(
      c(0.00, "#000000"), c(0.25, "#4d0000"),
      c(0.50, "#b30000"), c(0.75, "#ff8c00"),
      c(1.00, "#ffff00")
    ),
    # all NUMERIC_PALETTES entries fall through to plotly_colorscale
    plotly_colorscale(name)
  )
}

apply_heatmap_contrast <- function(zmat, clip = c(1, 99), gamma = 1) {
  zvec <- as.numeric(zmat)
  zvec <- zvec[is.finite(zvec)]
  if (!length(zvec)) return(zmat)
  
  lo <- stats::quantile(zvec, clip[1] / 100, names = FALSE, type = 7)
  hi <- stats::quantile(zvec, clip[2] / 100, names = FALSE, type = 7)
  if (!is.finite(lo) || !is.finite(hi) || lo >= hi) {
    lo <- min(zvec, na.rm = TRUE)
    hi <- max(zvec, na.rm = TRUE) + 1e-9
  }
  
  zclip <- pmin(pmax(zmat, lo), hi)
  z01   <- (zclip - lo) / (hi - lo)
  
  if (is.finite(gamma) && gamma > 0 && abs(gamma - 1) > 1e-9) {
    # gamma < 1 brightens faint values; gamma > 1 compresses lows
    z01 <- z01 ^ gamma
  }
  z01
}

build_dwell_heatmap <- function(df, nbins = 80) {
  # df needs columns x, y, tid, tim
  df <- df[is.finite(df$x) & is.finite(df$y) & is.finite(df$tim), ]
  if (!nrow(df)) return(NULL)
  
  df <- df[order(df$tid, df$tim), ]
  # dwell weight = time until next sample within same tid
  df <- df %>%
    group_by(tid) %>%
    mutate(dt = lead(tim) - tim) %>%
    ungroup()
  df$dt[!is.finite(df$dt) | df$dt < 0] <- 0
  
  xr <- range(df$x, na.rm = TRUE); if (xr[1] == xr[2]) xr <- xr + c(-0.5, 0.5)
  yr <- range(df$y, na.rm = TRUE); if (yr[1] == yr[2]) yr <- yr + c(-0.5, 0.5)
  
  xb <- seq(xr[1], xr[2], length.out = nbins + 1)
  yb <- seq(yr[1], yr[2], length.out = nbins + 1)
  ix <- findInterval(df$x, xb, all.inside = TRUE)
  iy <- findInterval(df$y, yb, all.inside = TRUE)
  
  # full grid with zeros
  df$ix <- factor(pmax(1, pmin(nbins, ix)), levels = 1:nbins)
  df$iy <- factor(pmax(1, pmin(nbins, iy)), levels = 1:nbins)
  
  z <- xtabs(df$dt ~ df$iy + df$ix)     # rows=y, cols=x
  z <- as.matrix(z)
  
  xcent <- 0.5 * (xb[-1] + xb[-length(xb)])
  ycent <- 0.5 * (yb[-1] + yb[-length(yb)])
  list(x = xcent, y = ycent, z = z)
}

auto_hm_bins <- function(npts) {
  if (npts < 2e5) 300 else if (npts < 1e6) 200 else 120
}

