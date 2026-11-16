# ---- Playback engine ----
# --- Play / Pause ---
observeEvent(input$trk_btn_play, {
  if (isTRUE(is_playing())) return(invisible(NULL))
  if (is_seq()) {
    # validate data BEFORE mutating loop state  -  req() would leave state inconsistent
    sd  <- isolate(seq_data())
    pre <- isolate(frames_store())
    tm  <- isolate(trace_map())
    if (!is.list(sd) || !length(sd$tids)) return(invisible(NULL))
    if (!is.list(pre)) return(invisible(NULL))
    if (is.null(tm))   return(invisible(NULL))
    # start/restart sequential
    loop$gen       <- loop$gen + 1L
    loop$playing   <- TRUE
    
    loop$sp            <- sd$sp
    loop$tids          <- sd$tids
    loop$trail_len     <- pre$trail_len
    loop$trail_nmax    <- pre$trail_nmax
    loop$pt_size       <- input$trk_pt_size
    loop$col_map       <- pre$col_map
    loop$markers_idx   <- tm$markers_idx
    loop$trail_start   <- tm$trail_start
    loop$frames_per_tid <- max(2L, n_frames_eff())
    loop$speed         <- as.numeric(input$trk_play_speed)
    loop$interval_sec  <- play_interval_for(loop$speed)
    loop$ms_per_unit   <- 0   # sequential mode uses a fixed frame interval scaled by the Playback speed slider
    loop$times <- lapply(loop$tids, function(nm) {
      tb <- loop$sp[[nm]]; r <- range(tb$t, na.rm = TRUE)
      if (!is.finite(r[1]) || !is.finite(r[2]) || r[1] == r[2]) r <- c(0, 1)
      seq(r[1], r[2], length.out = loop$frames_per_tid)
    })
    loop$cur_tid   <- 1L
    loop$cur_frame <- 1L
    is_playing(TRUE)
    
    run_sequential_loop_plain(loop$gen)
    
  } else {
    # start/restart global
    loop$ggen <- loop$ggen + 1L
    loop$gplaying <- TRUE 
    is_playing(TRUE)
    start_global_loop(loop$ggen)
  }
}, ignoreInit = TRUE)
start_global_loop <- function(id) {
  pre <- isolate(frames_store())
  if (!is.list(pre) || !length(pre$frames)) return(invisible(NULL))
  # capture maxk once; the per-frame delay is read live from loop$speed so the
  # "Playback speed (×)" slider takes effect mid-playback (see loop_once below).
  maxk     <- length(pre$frames)
  step     <- 1L
  loop$speed    <- as.numeric(isolate(input$trk_play_speed))
  loop$gplaying <- TRUE
  draw_frame(1L, force = TRUE)
  
  k_local <- 0L   # private counter, don't read tick()
  
  loop_once <- function() {
    # Wrap entire body so an unexpected error never silently kills the loop
    tryCatch({
      if (!isTRUE(loop$gplaying) || id != loop$ggen) return(invisible(NULL))
      
      k_local <<- k_local + step
      if (k_local < 1L) k_local <<- 1L
      if (k_local > maxk) {
        # Playback finished naturally  -  reset state inside domain so reactive
        # writes flush correctly from a later() callback context
        loop$gplaying <- FALSE
        loop$ggen <- loop$ggen + 1L
        shiny::withReactiveDomain(session$domain, {
          isolate({ is_playing(FALSE); slider_guard(0L) })
        })
        return(invisible(NULL))
      }

      # run Shiny ops inside domain, and reactive reads inside isolate()
      shiny::withReactiveDomain(session$domain, {
        isolate({
          tick(k_local)
          draw_frame(k_local)
        })
      })
      
      if (isTRUE(loop$gplaying) && id == loop$ggen) {
        later::later(loop_once, play_interval_for(loop$speed))
      }
    }, error = function(e) {
      # On error: log to console and keep the loop alive if still playing
      message("draw_frame error (frame ", k_local, "): ", conditionMessage(e))
      if (isTRUE(loop$gplaying) && id == loop$ggen) {
        later::later(loop_once, play_interval_for(loop$speed))
      }
    })
  }
  
  # kick off immediately
  later::later(loop_once, 0)
}
observeEvent(input$trk_chunk_sel, {
  # Skip the creation fire after upload/rechunk (value == the just-rendered chunk);
  # data_key already rebuilt. Only a genuine chunk switch gets through.
  if (identical(input$trk_chunk_sel %||% "__all__", loop$rendered_chunk %||% "__all__")) return()
  if (isTruthy(input$trk_chunk_sel)) fire_layout("chunk_sel")
}, ignoreInit = TRUE)
# Fires whenever new data is loaded (upload or rechunk), regardless of whether
# trk_chunk_sel changed. Ensures the chart always rebuilds with the new dataset.
# This is the SINGLE authoritative first-plot trigger on load; the chunk_sel /
# view_mode / z_scale observers compare against loop$rendered_* (set by the render)
# so their input-CREATION after load does not fire extra redundant rebuilds.
observeEvent(trk_data_key(), {
  if (!identical(trk_data_key(), "init")) fire_layout("data_key")
}, ignoreInit = TRUE)
observeEvent(input$trk_btn_pause, {
  loop$gplaying <- FALSE   # stop global loop (belt-and-suspenders with ggen check)
  loop$playing  <- FALSE   # stop sequential loop
  loop$ggen <- loop$ggen + 1L
  loop$gen  <- loop$gen  + 1L
  is_playing(FALSE)
  slider_guard(0L)
}, ignoreInit = TRUE)

observeEvent(input$trk_btn_stop, {
  is_playing(FALSE)
  loop$playing  <- FALSE
  loop$gplaying <- FALSE
  loop$ggen     <- loop$ggen + 1L
  loop$gen      <- loop$gen + 1L
  slider_guard(0L)
  
  loop$cur_tid   <- 1L
  loop$cur_frame <- 1L
  
  if (isTRUE(isolate(is_seq()))) {
    sd <- isolate(seq_data())
    if (!is.null(sd) && length(sd$tids) > 0) draw_tid_point(sd$tids[1], 1L)
  } else {
    pre <- isolate(frames_store())
    if (!is.null(pre) && length(pre$frames) > 0) {
      tick(1L)
      draw_frame(1L, force = TRUE)   # <-- force a redraw while "not playing"
    }
  }
})
apply_auto_tune <- function() {
  df <- isolate(trk()); req(nrow(df) > 0)
  mode <- if (isTRUE(isolate(is_simul()))) "simul" else "seq"
  s <- smart_defaults(df, mode = mode, target_frames = 200, target_fps = 15)

  # Set the plain (non-reactive) flag BEFORE queuing updateNumericInput calls.
  # Browser round-trips from those calls arrive ~100-400ms later; rebuild_tick
  # must stay suppressed until AFTER fire_frames() has consumed the new values.
  # fire_frames() is delayed 500 ms so the new input values are in place when
  # frames_store() re-evaluates.  auto_tuning resets 1.5 s after that.
  loop$auto_tuning <- TRUE
  later::later(function() { loop$auto_tuning <- FALSE }, delay = 1.5)

  updateNumericInput(session, "trk_play_total",      value = s$nice_total)
  updateNumericInput(session, "trk_trail_len",       value = signif(s$trail_len, 3))
  updateNumericInput(session, "trk_target_markers",  value = s$target_markers)
  updateNumericInput(session, "trk_bin_dt",          value = s$nice_total / s$n_frames)

  # Delay fire_frames until after input round-trips arrive (~200-400 ms),
  # so frames_store() reads the newly-tuned bin_dt / n_frames / total values.
  later::later(function() { fire_frames() }, delay = 0.5)
}

observeEvent(input$trk_auto_tune, {
  apply_auto_tune()
}, ignoreInit = TRUE)
# --- Global (prebinned) mode stays on invalidateLater (fast) ---
run_sequential_loop_plain <- function(id) {
  tryCatch({
    # cancel if paused or superseded
    if (id != loop$gen || !isTRUE(loop$playing)) return(invisible(NULL))
    
    tids <- loop$tids; sp <- loop$sp
    if (is.null(tids) || is.null(sp) || !length(tids)) return(invisible(NULL))
    
    ti <- loop$cur_tid
    fi <- loop$cur_frame
    if (ti < 1L || ti > length(tids)) ti <- 1L
    if (fi < 1L || fi > loop$frames_per_tid) fi <- 1L
    nm <- tids[ti]; tb <- sp[[nm]]
    
    if (is.null(tb) || !nrow(tb)) {
      if (ti >= length(tids)) {
        loop$playing <- FALSE
        loop$gen <- loop$gen + 1L
        shiny::withReactiveDomain(session$domain, {
          isolate({ is_playing(FALSE); slider_guard(0L) })
        })
        return(invisible(NULL))
      }
      loop$cur_tid <- ti + 1L
      loop$cur_frame <- 1L
      later::later(function() { run_sequential_loop_plain(id) }, 0.05)
      return(invisible(NULL))
    }
    
    # current frame time + where we're going next
    tnow    <- loop$times[[ti]][fi]
    next_fi <- if (fi < loop$frames_per_tid) fi + 1L else fi
    next_ti <- ti
    if (fi >= loop$frames_per_tid && ti < length(tids)) {
      next_ti <- ti + 1L
      next_fi <- 1L
    }
    
    # draw this frame
    w <- which(tb$t >= (tnow - loop$trail_len) & tb$t <= tnow)
    if (length(w) > loop$trail_nmax) w <- tail(w, loop$trail_nmax)
    xx <- if (length(w)) c(tb$x[w], NA_real_) else numeric(0)
    yy <- if (length(w)) c(tb$y[w], NA_real_) else numeric(0)
    
    idx <- nearest_idx(tb$t, tnow)
    mx  <- if (is.na(idx)) numeric(0) else tb$x[idx]
    my  <- if (is.na(idx)) numeric(0) else tb$y[idx]
    
    col   <- as.character(loop$col_map[[nm]])
    proxy <- plotlyProxy("trk_track2D", session)
    if (!is.na(loop$trail_start))
      plotlyProxyInvoke(proxy, "restyle", list(x = list(xx), y = list(yy), `line.color` = list(col)), loop$trail_start)
    plotlyProxyInvoke(proxy, "restyle", list(x = list(mx), y = list(my), `marker.size` = loop$pt_size, `marker.color` = list(col)), loop$markers_idx)
    # Skip relayout during playback to avoid scaleanchor jitter
    
    # advance state
    if (fi < loop$frames_per_tid) {
      loop$cur_frame <- fi + 1L
    } else if (ti < length(tids)) {
      loop$cur_tid   <- ti + 1L
      loop$cur_frame <- 1L
    } else {
      loop$playing <- FALSE
      loop$gen <- loop$gen + 1L
      shiny::withReactiveDomain(session$domain, {
        isolate({ is_playing(FALSE); slider_guard(0L) })
      })
      return(invisible(NULL))
    }
    
    # compute delay (ms_per_unit -> variable; else fixed play_interval)
    # The "Playback speed (×)" multiplier scales every branch: higher speed =>
    # shorter delay. Everything is floored at MIN_INTERVAL_SEC so the loop can
    # never outpace the websocket (the old 0.01s firehose is what made playback
    # "refuse to stop" - it starved the Stop click + reactive flush).
    speed     <- if (is.finite(loop$speed) && loop$speed > 0) loop$speed else 1
    floor_sec <- play_interval_for(speed)
    if (is.finite(loop$ms_per_unit) && loop$ms_per_unit > 0) {
      tnext    <- loop$times[[ti]][min(next_fi, loop$frames_per_tid)]
      dt_units <- if (next_ti != ti) 0 else (tnext - tnow)
      delay_sec <- max(floor_sec, abs(dt_units) * (loop$ms_per_unit / 1000) / speed)
      if (next_ti != ti) delay_sec <- max(delay_sec, 0.15 / speed)
    } else {
      delay_sec <- floor_sec
    }
    
    if (id == loop$gen && isTRUE(loop$playing)) {
      later::later(function() { run_sequential_loop_plain(id) }, delay_sec)
    }
  }, error = function(e) {
    message("Sequential loop error (tid=", loop$cur_tid, " frame=", loop$cur_frame, "): ", conditionMessage(e))
    if (id == loop$gen && isTRUE(loop$playing)) {
      later::later(function() { run_sequential_loop_plain(id) }, loop$interval_sec)
    }
  })
}

# Live speed updates: both loops read loop$speed on each scheduled tick, so
# dragging the slider during playback changes the rate without restarting.
observeEvent(input$trk_play_speed, {
  loop$speed <- as.numeric(input$trk_play_speed)
  if (isTRUE(loop$playing)) loop$interval_sec <- play_interval_for(loop$speed)
}, ignoreInit = TRUE)

observeEvent(input$trk_play_mode, {
  is_playing(FALSE); loop$playing <- FALSE; loop$gplaying <- FALSE
  loop$ggen <- loop$ggen + 1L
  loop$gen    <- loop$gen + 1L
  slider_guard(0L)
})

# ---- Slider (seek without auto-pause) ----
output$trk_frame_slider_ui <- renderUI({
  req(!is.na(layout_trigger()))
  if (is_seq()) {
    sd <- seq_data(); req(is.list(sd), length(sd$tids) > 0)
    ti <- current_tid_idx(); if (ti < 1L || ti > length(sd$tids)) ti <- 1L
    nm <- sd$tids[ti]
    sliderInput(
      "trk_frame_idx", paste0("Step (tid = ", nm, ")"),
      min = 1, max = max(2L, n_frames_eff()), value = 1, step = 1, width = "100%"
    )
  } else {
    pre <- frames_store(); req(length(pre$frames) > 0)
    sliderInput("trk_frame_idx", "Frame",
                min = 1, max = length(pre$frames), value = 1, step = 1, width = "100%")
  }
})


observeEvent(input$trk_frame_idx, {
  req(!is.na(layout_trigger()))
  # Consume programmatic echoes FIRST, before the is_playing() gate. draw_frame bumps
  # slider_guard for every updateSliderInput it pushes; if the is_playing() gate ran
  # first those credits were never spent during playback, so the counter accumulated
  # and was then wiped by slider_guard(0L) at stop/end - after which late-arriving
  # echoes from the last frames sailed through and re-drove draw_frame.
  guard <- isolate(slider_guard())
  if (guard > 0L) {
    slider_guard(guard - 1L)
    return(invisible(NULL))
  }
  if (isTRUE(is_playing())) return(invisible(NULL))

  if (is_seq()) {
    # unchanged: sequential seeks by time within the current TID
    sd <- seq_data(); req(is.list(sd), length(sd$tids) > 0)
    ti <- current_tid_idx(); if (ti < 1L || ti > length(sd$tids)) ti <- 1L
    nm <- sd$tids[ti]; tb <- sd$sp[[nm]]; req(nrow(tb) > 0)
    
    k <- max(1L, min(as.integer(input$trk_frame_idx), max(2L, n_frames_eff())))
    tnow <- if (!is.null(loop$times) && length(loop$times) >= ti) {
      loop$times[[ti]][k]
    } else {
      r <- range(tb$t, na.rm = TRUE)
      seq(r[1], r[2], length.out = max(2L, n_frames_eff()))[k]
    }
    i <- nearest_idx(tb$t, tnow)
    draw_tid_point(nm, i)
    
  } else {
    # simultaneous seeks into the prebuilt frames
    pre <- frames_store(); req(length(pre$frames) > 0)
    k <- max(1L, min(as.integer(input$trk_frame_idx), length(pre$frames)))
    tick(k)                     # optional: keep any tick-based UI in sync
    draw_frame(k, force = TRUE) # <- force draw while paused
  }
})

# ---- Reset state + push frame 1 on data rebuild --------------------------------
# For fire_layout: renderPlotly handles the full rebuild including frame 1 data.
#   The plot_trigger observer just resets play state and absorbs the slider
#   round-trip from trk_frame_slider_ui re-rendering.
# For fire_frames (Plot Data, rebuild_tick, auto-tune): renderPlotly does NOT
#   re-fire (frames_store is now isolated there). After frames_store re-evaluates,
#   we push frame 1 data + new axis ranges via proxy so the chart updates in-place
#   with no WebGL reinit and no blank.
observeEvent(plot_trigger(), {
  is_playing(FALSE); tick(0L)
  loop$playing  <- FALSE
  loop$gplaying <- FALSE
  loop$ggen     <- loop$ggen + 1L
  loop$gen      <- loop$gen  + 1L
  was_layout <- isTRUE(loop$last_fire_layout)
  slider_guard(if (was_layout) 3L else 0L)

  if (!was_layout) {
    # fire_frames path: push new frame 1 + axes via proxy after frames_store settles.
    # Skip entirely in full-tracks mode - the static plot needs no proxy update.
    if (identical(isolate(input$trk_view_mode_main), "full")) return()
    # Capture generation so a stale callback from a prior session cannot jump the
    # display to frame 1 after a new upload/playback has already begun.
    snap_ggen <- loop$ggen
    snap_gen  <- loop$gen
    later::later(function() {
      shiny::withReactiveDomain(session$domain, {
        tryCatch({
          # Abort if a new upload or play session has superseded this callback.
          if (loop$ggen != snap_ggen || loop$gen != snap_gen) return()
          pre <- isolate(frames_store())
          if (!is.list(pre) || !length(pre$frames)) return()
          tm <- isolate(trace_map())
          if (is.null(tm)) return()

          # push frame 1 marker data (colors already baked into frames)
          draw_frame(1L, force = TRUE)

          # If per-TID trail traces exist, update their line colors too
          # (trails are set at renderPlotly time; conf coloring needs a proxy push).
          # Skipped in per-point colour modes (dark cycles / speed): there draw_frame
          # above already coloured each trail from its leading point, and pushing the
          # per-TID colours here would clobber it back to one colour per track.
          if (identical(tm$trails_mode, "per_tid") && tm$n_tid > 0 &&
              !isTRUE(pre$color_by_dark)) {
            line_cols <- vector("list", tm$n_tid)
            for (i in seq_len(tm$n_tid)) {
              tid_i <- if (i <= length(pre$tids)) pre$tids[i] else ""
              line_cols[[i]] <- unname(pre$col_map[[as.character(tid_i)]])
            }
            plotlyProxyInvoke(plotlyProxy("trk_track2D", session), "restyle",
              list(`line.color` = line_cols),
              seq(tm$trail_start, tm$trail_start + tm$n_tid - 1L))
          }

          # push updated axis ranges so filtered/new data is fully visible
          f1     <- pre$frames[[1]]
          bx     <- f1$markers_x[is.finite(f1$markers_x)]
          by     <- f1$markers_y[is.finite(f1$markers_y)]
          if (!length(bx) || !length(by)) {
            sp <- pre$sp
            bx <- unlist(lapply(sp, `[[`, "x"), use.names = FALSE); bx <- bx[is.finite(bx)]
            by <- unlist(lapply(sp, `[[`, "y"), use.names = FALSE); by <- by[is.finite(by)]
          }
          if (length(bx) > 0 && length(by) > 0) {
            xr <- range(bx); yr <- range(by)
            xpad <- max(1, diff(xr) * 0.04); ypad <- max(1, diff(yr) * 0.04)
            trk_axis_x_span(diff(xr) + 2 * xpad)
            trk_locked_ranges(list(x = c(xr[1]-xpad, xr[2]+xpad),
                                  y = c(yr[1]-ypad, yr[2]+ypad)))
            plotlyProxyInvoke(plotlyProxy("trk_track2D", session), "relayout", list(
              `xaxis.range[0]` = xr[1]-xpad, `xaxis.range[1]` = xr[2]+xpad,
              `xaxis.autorange` = FALSE,
              `yaxis.range[0]` = yr[1]-ypad, `yaxis.range[1]` = yr[2]+ypad,
              `yaxis.autorange` = FALSE
            ))
          }
        }, error = function(e) NULL)
      })
    }, delay = 0.2)
  } else {
    # fire_layout path: renderPlotly rebuilds the chart fresh. Push a grid/colour
    # relayout shortly after so the axis styling is consistent regardless of what
    # Plotly does internally during the DOM recreate + newPlot sequence.
    later::later(function() {
      shiny::withReactiveDomain(session$domain, {
        tryCatch(trk_axis_relayout(), error = function(e) NULL)
      })
    }, delay = 0.35)
  }
}, ignoreNULL = TRUE, ignoreInit = TRUE)
observeEvent(input$trk_play_mode, { loop$playing <- FALSE })   # switching modes stops loop

# ---- Live marker size update ----
observeEvent(input$trk_pt_size, {
  req(!is.na(layout_trigger()))
  tm <- isolate(trace_map()); req(!is.null(tm))
  plotlyProxyInvoke(plotlyProxy("trk_track2D", session), "restyle",
                    list(`marker.size` = input$trk_pt_size), tm$markers_idx)
})

trk_axis_relayout <- function() {
  req(!is.na(layout_trigger()))
  hide  <- isTRUE(input$trk_hide_grid)
  dark  <- isTRUE(input$trk_dark_mode)
  is3d  <- isTRUE(use_3d())
  fg    <- if (dark) "#FFFFFF" else "#63666A"
  gc    <- if (dark) "rgba(255,255,255,0.18)" else "rgba(0,0,0,0.18)"
  bg    <- if (dark) "#000000" else "#FFFFFF"
  # Keep the 2D axis titles matching the active viewing plane so a dark-mode /
  # grid relayout does not silently revert an XZ/YZ side view to "X/Y".
  plane2 <- if (!is3d) (input$trk_view_plane %||% "XY") else "XY"
  plab   <- switch(plane2,
    XZ = list(x = "X (nm)", y = "Z (nm)"),
    YZ = list(x = "Y (nm)", y = "Z (nm)"),
    list(x = "X (nm)", y = "Y (nm)"))

  if (is3d) {
    # Recompute equal ranges so the relayout doesn't reset autorange.
    pre  <- isolate(frames_store())
    sp   <- if (is.list(pre)) pre$sp else list()
    all_x <- unlist(lapply(sp, `[[`, "x"), use.names = FALSE)
    all_y <- unlist(lapply(sp, `[[`, "y"), use.names = FALSE)
    all_z <- unlist(lapply(sp, `[[`, "z"), use.names = FALSE)
    lim3  <- get_equal_limits(all_x, all_y, all_z)

    ax3_style <- list(
      showgrid = !hide, zeroline = !hide, showticklabels = !hide,
      gridcolor = gc, color = fg, backgroundcolor = bg
    )
    make_ax3 <- function(title_txt, lim) {
      c(list(title = title_txt, range = lim, autorange = FALSE), ax3_style)
    }
    plotlyProxyInvoke(
      plotlyProxy("trk_track2D", session),
      "relayout",
      list(
        paper_bgcolor = bg,
        font = list(color = fg),
        scene = list(
          bgcolor = bg,
          xaxis = make_ax3(if (hide) "" else "X (nm)", lim3$x),
          yaxis = make_ax3(if (hide) "" else "Y (nm)", lim3$y),
          zaxis = make_ax3(if (hide) "" else "Z (nm)", lim3$z),
          aspectmode  = "manual",
          aspectratio = list(x = 1, y = 1, z = 1)
        )
      )
    )
  } else {
    ovl <- if (plane2 == "XY") tryCatch(isolate(trk_image_overlay()), error = function(e) NULL) else NULL
    # Include locked ranges so the xaxis/yaxis object-form relayout doesn't
    # discard range + autorange that were set by the initial layout or proxy.
    lr         <- isolate(trk_locked_ranges())
    x_range_kv <- if (!is.null(lr)) list(range = lr$x, autorange = FALSE) else list()
    y_range_kv <- if (!is.null(lr)) list(range = lr$y, autorange = FALSE) else list()
    plotlyProxyInvoke(
      plotlyProxy("trk_track2D", session),
      "relayout",
      list(
        paper_bgcolor = bg, plot_bgcolor = bg,
        font = list(color = fg),
        xaxis = c(list(
          title = if (hide) "" else plab$x,
          showgrid = !hide, zeroline = !hide, showticklabels = !hide,
          gridcolor = gc, color = fg
        ), x_range_kv),
        yaxis = c(list(
          title = if (hide) "" else plab$y,
          showgrid = !hide, zeroline = !hide, showticklabels = !hide,
          gridcolor = gc, color = fg,
          scaleanchor = "x", scaleratio = 1
        ), y_range_kv),
        images = if (!is.null(ovl)) list(ovl) else list()
      )
    )
  }

  # Restyle faint static trace color when dark mode toggles
  tm <- isolate(trace_map())
  if (!is.null(tm) && !is.null(tm$static_idx) && !is.na(tm$static_idx)) {
    faint_color <- if (dark) {
      if (is3d) "rgba(255,255,255,0.55)" else "rgba(255,255,255,0.25)"
    } else {
      if (is3d) "rgba(80,80,80,0.55)" else "rgba(80,80,80,0.15)"
    }
    plotlyProxyInvoke(
      plotlyProxy("trk_track2D", session),
      "restyle",
      list(`line.color` = faint_color),
      list(tm$static_idx)
    )
  }
}

observeEvent(input$trk_hide_grid,  { trk_axis_relayout() }, ignoreInit = TRUE)
observeEvent(input$trk_dark_mode,  { trk_axis_relayout() }, ignoreInit = TRUE)

# Scale bar: push via proxy whenever any scale-bar input changes.
trk_push_scale_bar <- function() {
  req(!is.na(layout_trigger()))
  proxy <- plotlyProxy("trk_track2D", session)
  if (!isTRUE(input$trk_show_scale_bar)) {
    plotlyProxyInvoke(proxy, "relayout",
                      list(shapes = list(), annotations = list()))
    return()
  }
  sb_len   <- suppressWarnings(as.numeric(input$trk_scale_bar_length))
  sb_thick <- suppressWarnings(as.numeric(input$trk_scale_bar_thickness))
  sb_font  <- suppressWarnings(as.numeric(input$trk_scale_font_size))
  sb <- make_scale_bar(sb_len, trk_axis_x_span() %||% 1000, sb_thick, dark = isTRUE(input$trk_dark_mode),
                       scale_font_size = if (is.finite(sb_font) && sb_font > 0) sb_font else 18)
  if (!is.null(sb))
    plotlyProxyInvoke(proxy, "relayout",
                      list(shapes = sb$shapes, annotations = sb$annotations))
}

observe({
  # Only read the relayout event once the track plot has rendered + registered
  # AND the Tracking tab is on screen, so plotly's deferred "source not
  # registered" warning is never scheduled. trk_plot_ready stays FALSE until a
  # plot with data renders (layout_trigger alone can be set with no frames yet).
  if (!isTRUE(trk_plot_ready()) ||
      !identical(input$main_navbar, "tab_tracking")) return()
  req(!is.na(layout_trigger()))
  ed <- suppressWarnings(plotly::event_data("plotly_relayout", source = "track"))
  if (!is.list(ed)) return()

  x0 <- suppressWarnings(as.numeric(ed[["xaxis.range[0]"]]))
  x1 <- suppressWarnings(as.numeric(ed[["xaxis.range[1]"]]))
  if (length(x0) == 1 && length(x1) == 1 && is.finite(x0) && is.finite(x1) && x1 > x0) {
    trk_axis_x_span(x1 - x0)
    y0 <- suppressWarnings(as.numeric(ed[["yaxis.range[0]"]]))
    y1 <- suppressWarnings(as.numeric(ed[["yaxis.range[1]"]]))
    if (length(y0) == 1 && length(y1) == 1 && is.finite(y0) && is.finite(y1) && y1 > y0)
      trk_zoom_ranges(list(x = c(x0, x1), y = c(y0, y1)))
  } else if (isTRUE(ed[["xaxis.autorange"]])) {
    axis_state <- trk_default_axis_state()
    if (!is.null(axis_state)) trk_axis_x_span(axis_state$x_span)
    trk_zoom_ranges(NULL)
  } else {
    return()
  }

  if (isTRUE(input$trk_show_scale_bar)) trk_push_scale_bar()
})

# ---- Spatial crop: "Crop to view" / "Clear" --------------------------------
observeEvent(input$trk_crop_to_view, {
  # The spatial crop filters the real x,y. In an XZ/YZ side view the plotted
  # vertical axis is Z, so the zoom box's "y" is a z-range - cropping on it would
  # silently filter the wrong axis. Require the XY plane for cropping.
  if (!identical(input$trk_view_plane %||% "XY", "XY")) {
    showNotification("Crop to view only works in the XY plane. Switch the viewing plane to XY first.",
                     type = "warning", duration = 6)
    return()
  }
  zr  <- trk_zoom_ranges()
  pct <- suppressWarnings(as.numeric(input$trk_crop_pct))
  if (is.null(zr)) {
    showNotification("Zoom into a region first, then press Crop to view.", type = "warning", duration = 5)
    return()
  }
  if (!is.finite(pct) || pct <= 0 || pct > 100) pct <- 50
  trk_spatial_filter(list(x = zr$x, y = zr$y, pct_thr = pct / 100))
  data_ver(isolate(data_ver()) + 1L)
  fire_frames()
  showNotification(
    sprintf("Spatial crop active: x [%.0f, %.0f] nm, y [%.0f, %.0f] nm (≥%.0f%% per TID).",
            zr$x[1], zr$x[2], zr$y[1], zr$y[2], pct),
    type = "message", duration = 5
  )
}, ignoreInit = TRUE)

observeEvent(input$trk_clear_spatial, {
  trk_spatial_filter(NULL)
  data_ver(isolate(data_ver()) + 1L)
  fire_frames()
  showNotification("Spatial crop cleared.", type = "message", duration = 3)
}, ignoreInit = TRUE)

# Info UI: show active box or a hint
output$trk_spatial_info <- renderUI({
  sf <- trk_spatial_filter()
  if (is.null(sf)) {
    tags$p(style = "font-size:0.82em; color:#888; margin:2px 0;",
           "No spatial filter active. Zoom into a region and press Crop to view.")
  } else {
    tags$p(style = "font-size:0.82em; color:#2070a0; margin:2px 0;",
           sprintf("Active: x [%.0f, %.0f] nm  |  y [%.0f, %.0f] nm  |  ≥%.0f%% per TID",
                   sf$x[1], sf$x[2], sf$y[1], sf$y[2], sf$pct_thr * 100))
  }
})

observeEvent(list(input$trk_show_scale_bar, input$trk_scale_bar_length,
                  input$trk_scale_bar_thickness, input$trk_scale_font_size,
                  input$trk_dark_mode),
  { trk_push_scale_bar() }, ignoreInit = TRUE)
