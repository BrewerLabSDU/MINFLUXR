# --------- 7) Trace spread inspection  ----------
sd_grey_tbl <- reactive({
  sd_by_tid(raw_df(), min_n = input$sd_min_locs)
})

sd_black_tbl <- reactive({
  sd_by_tid(filtered_df(), min_n = input$sd_min_locs)
})

output$sd_summary_ui <- renderUI({
  req(input$view_SD)
  df  <- filtered_df()    
  med <- compute_sd_medians(df, min_n = input$sd_min_locs)
  
  # One axis per line, in the pyMINFLUX reporting form: median +/- MAD with the
  # interquartile range after it. The MAD is the width of the trace-to-trace
  # distribution, NOT an error on the median - the IQR is shown next to it so
  # the asymmetry of that distribution stays visible rather than being hidden
  # behind a symmetric +/-.
  fmt_axis <- function(axis, med_v, mad_v, q1, q3) {
    if (!is.finite(med_v %||% NA_real_)) return(NULL)
    txt <- sprintf("Median trace spread %s: %.2f ± %.2f nm", axis, med_v, mad_v %||% NA_real_)
    if (is.finite(q1 %||% NA_real_) && is.finite(q3 %||% NA_real_))
      txt <- sprintf("%s  (IQR %.2f–%.2f)", txt, q1, q3)
    helpText(txt)
  }

  els <- list(
    fmt_axis("X", med$x, med$mad_x, med$q1_x, med$q3_x),
    fmt_axis("Y", med$y, med$mad_y, med$q1_y, med$q3_y)
  )
  if (!is.null(med$z)) {
    els <- append(els, list(fmt_axis("Z", med$z, med$mad_z, med$q1_z, med$q3_z)))
  }
  els <- Filter(Negate(is.null), els)
  if (length(els) && is.finite(med$n_traces %||% NA_real_))
    els <- append(els, list(helpText(sprintf("Across %d traces with ≥ %d localizations.",
                                             as.integer(med$n_traces),
                                             as.integer(input$sd_min_locs %||% SMLM_DEFAULTS$sd_min_locs)))))
  do.call(tagList, els)
  })
output$sd_ranges_ui <- renderUI({
  req(input$view_SD)
  sdr <- sd_grey_tbl()

  get_max <- function(v) {
    v <- v[is.finite(v)]
    if (!length(v)) return(1)
    mx <- max(v, na.rm = TRUE)
    if (!is.finite(mx) || mx <= 0) 1 else mx
  }
  # Initialise widget from sd_filt (the committed filter state) so that a
  # previously loaded or user-set range is shown when the panel is reopened.
  # Falls back to [0, mx] when the filter is at full range (Inf upper bound).
  # When "lock filters" is active, preserve the locked upper bound exactly and
  # extend the widget max to fit it - prevents clamping from silently resetting
  # sd_filt via the sd_x_range observer when the new file has a smaller max SD.
  locked <- isTRUE(input$lock_filters)
  init_val <- function(filt_rng, mx) {
    lo <- max(0, filt_rng[1])
    hi <- if (is.infinite(filt_rng[2])) mx
          else if (locked) filt_rng[2]
          else min(mx, filt_rng[2])
    c(lo, hi)
  }
  eff_max <- function(mx, hi) if (locked && is.finite(hi) && hi > mx) hi else mx

  mx_x <- get_max(sdr$sd_x)
  mx_y <- get_max(sdr$sd_y)
  # Use isolate() so that changes to sd_filt do NOT cause sd_ranges_ui to
  # re-render.  Without isolate, setting sd_filt inside the SD range observer
  # would invalidate sd_ranges_ui → new widgets fire → observer fires again →
  # loop.  With isolate, sd_filt is read once when the panel opens and the
  # output is not tracked as a reactive dependency of sd_filt.
  cur <- isolate(list(x = sd_filt$x, y = sd_filt$y, z = sd_filt$z))
  vx  <- init_val(cur$x, mx_x)
  vy  <- init_val(cur$y, mx_y)

  tagList(
    shinyWidgets::numericRangeInput("sd_x_range", "Allowable trace spread in X",
      value = vx, min = 0, max = eff_max(mx_x, vx[2])),
    shinyWidgets::numericRangeInput("sd_y_range", "Allowable trace spread in Y",
      value = vy, min = 0, max = eff_max(mx_y, vy[2])),
    if ("sd_z" %in% names(sdr) && any(is.finite(sdr$sd_z))) {
      mx_z <- get_max(sdr$sd_z)
      vz   <- init_val(cur$z, mx_z)
      shinyWidgets::numericRangeInput("sd_z_range", "Allowable trace spread in Z",
        value = vz, min = 0, max = eff_max(mx_z, vz[2]))
    },
    helpText("An empty box means no bound on that side: blank lower reads as 0, ",
             "blank upper as no upper limit.")
  )
})
observeEvent(input$view_SD, {
  req(input$file)
  df  <- filtered_df()
  med <- compute_sd_medians(df, min_n = input$sd_min_locs)
  isolate({
    updateNumericInput(session, "target_sd_x", value = round(med$x %||% NA_real_, 2))
    updateNumericInput(session, "target_sd_y", value = round(med$y %||% NA_real_, 2))
    if (!is.null(med$z)) {
      updateNumericInput(session, "target_sd_z", value = round(med$z %||% NA_real_, 2))
    }
  })
})
# If user toggles q95 checkbox, immediately refresh the preset ranges
observeEvent(input$sd_use_q95, {
  req(input$view_SD)
  sdr <- sd_grey_tbl()
  req(is.data.frame(sdr), nrow(sdr) > 0)
  
  set_preset <- function(col, inputId) {
    if (!col %in% names(sdr)) return()
    v <- as.numeric(sdr[[col]])
    v <- v[is.finite(v)]
    if (!length(v)) return()
    
    mx  <- max(v, na.rm = TRUE); if (!is.finite(mx) || mx <= 0) mx <- 1
    top <- if (isTRUE(input$sd_use_q95)) stats::quantile(v, 0.95, na.rm = TRUE, names = FALSE) else mx
    if (!is.finite(top) || top <= 0) top <- mx
    
    shinyWidgets::updateNumericRangeInput(session, inputId, value = c(0, round(top, 3)))
  }
  
  set_preset("sd_x", "sd_x_range")
  set_preset("sd_y", "sd_y_range")
  if ("sd_z" %in% names(sdr) && any(is.finite(sdr$sd_z))) set_preset("sd_z", "sd_z_range")
}, ignoreInit = TRUE)
has_sd_z <- reactive({
  sda <- sd_all_tbl()
  is.data.frame(sda) && "sd_z" %in% names(sda) && any(is.finite(sda$sd_z))
})
sd_all_tbl <- reactive({
  df0 <- uploaded_df()
  shiny::validate(shiny::need(!is.null(df0) && nrow(df0) > 0, "No data"))
  sd_by_tid(df0, min_n = input$sd_min_locs)            # RAW
})

sd_filtered_tbl <- reactive({
  dfF <- filtered_df()
  shiny::validate(shiny::need(!is.null(dfF) && nrow(dfF) > 0, "No data after current filters"))
  sd_by_tid(dfF, min_n = input$sd_min_locs)
})

# Effective spread display ranges: committed filter bounds with Inf replaced
# by the actual data maximum.  Mirrors exactly what the spread sliders show,
# even when the "Inspect trace spread" panel has never been opened.
smlm_sd_display_ranges <- reactive({
  sda    <- tryCatch(sd_all_tbl(), error = function(e) NULL)
  mx_fn  <- function(v) { v <- v[is.finite(v)]; if (length(v)) max(v) else NA_real_ }
  eff    <- function(filt_val, mx) if (is.infinite(filt_val)) mx else filt_val
  mx_x   <- if (is.data.frame(sda)) mx_fn(sda$sd_x) else NA_real_
  mx_y   <- if (is.data.frame(sda)) mx_fn(sda$sd_y) else NA_real_
  mx_z   <- if (is.data.frame(sda) && "sd_z" %in% names(sda)) mx_fn(sda$sd_z) else NA_real_
  list(
    x = c(sd_filt$x[1], eff(sd_filt$x[2], mx_x)),
    y = c(sd_filt$y[1], eff(sd_filt$y[2], mx_y)),
    z = c(sd_filt$z[1], eff(sd_filt$z[2], mx_z))
  )
})
observe({
  for (id in c("sd_hist_x", "sd_hist_y", "sd_hist_z")) {
    outputOptions(output, id, suspendWhenHidden = FALSE)
  }
})


# Dynamic container that shows the right set of histograms (X,Y and Z if present)
output$sd_hist_ui <- renderUI({
  req(input$view_SD)
  # Decide if Z is present from uploaded data
  df_u <- uploaded_df()
  shiny::validate(shiny::need(!is.null(df_u) && nrow(df_u) > 0, "No data uploaded"))
  has_z <- "z" %in% names(df_u) && any(is.finite(df_u$z))
  
  els <- list(
    plotlyOutput("sd_hist_x", height = 220),
    plotlyOutput("sd_hist_y", height = 220)
  )
  if (isTRUE(has_sd_z())) {
    els <- append(els, list(plotlyOutput("sd_hist_z", height = 220)))
  }
  do.call(tagList, els)
})

# ==========================================================================
# Shared histogram ranges across sigma X / Y / Z - OPT IN
# ==========================================================================
# GOAL: put the three spread histograms on ONE set of bin edges and ONE count
# axis, so the panels can be read against each other and dropped into a figure
# as a set.
#
# WHY IT IS NOT THE DEFAULT. Per-axis ranges are the better view for reading a
# single axis: MINFLUX z precision is genuinely different from xy, so a shared
# range spends most of the x panel's width on empty space. Sharing is a
# presentation choice, so it is opt-in and changes nothing else.
#
# The bin EDGES are shared too, not just the axis limits - histograms on
# different edges are not comparable however the axes are labelled, because
# the bar heights answer different questions.
#
# RETURNS NULL when off, so each renderer falls back to its own range.
# ==========================================================================

# Median +/- normalized MAD, ONE LINE PER CHANNEL, each in its channel's colour.
#
# The bars were split by channel before this was, which made the panel say two
# different things at once: coloured distributions with a single pooled median
# drawn across them. That pooled median is not a property of either channel -
# it lands between them, weighted by whichever channel happens to have more
# traces, so on an unbalanced sample (which is normal: one target is usually
# denser than the other) it sits close to the bigger channel and describes the
# smaller one badly. It is also the number the spread FILTER is set from, and
# that filter is shared across channels, so reading it per channel is what
# tells you whether one threshold cuts them evenly.
#
# Falls back to the pooled line when there is only one channel.
.sd_add_channel_medians <- function(fig, cur, colname, y_top, show_iqr) {
  cm <- tryCatch(sd_channel_map(), error = function(e) NULL)
  if (is.null(cm) || !"tid" %in% names(cur)) return(NULL)
  j <- dplyr::left_join(cur, cm, by = "tid")
  if (!"channel" %in% names(j) || all(is.na(j$channel))) return(NULL)
  chs <- sort(unique(j$channel[!is.na(j$channel)]))
  if (length(chs) < 2L) return(NULL)

  # Stagger the label heights so two nearby medians do not overprint. Channels
  # that differ by only a few percent are exactly the case worth reading.
  n <- length(chs)
  for (i in seq_along(chs)) {
    sub <- j[which(j$channel == chs[i]), , drop = FALSE]
    v   <- sub[[colname]]
    if (!length(v[is.finite(v)])) next
    col <- sub$channel_color[1]; if (is.na(col)) col <- "#666666"
    lab <- sigma_label_med_pm_sd(v)
    if (!is.null(lab)) lab <- paste0(sub$channel_label[1], "  ", lab)
    yl  <- y_top * (1 - 0.10 * (i - 1) / max(1, n - 1))
    fig <- add_sigma_vline(fig, stats::median(v, na.rm = TRUE), yl, col, lab)
    if (isTRUE(show_iqr)) fig <- add_sigma_iqr_lines(fig, v, yl, col)
  }
  fig
}

# ── TRACE SPREAD BY CHANNEL ───────────────────────────────────────────────
# Trace spread is the ONE metadata distribution that is split by channel, and
# deliberately the only one: EFO, CFR, DCR, ECO and dwell are properties of the
# acquisition, shared by every channel in it, so splitting them would add noise
# without adding information. Spread is different - it IS the localization
# precision, precision follows photon yield, and photon yield is a property of
# the dye. It is also what sets the de-blink merge distance, so seeing the two
# distributions apart is what tells you whether one shared distance is honest.
sd_channel_map <- reactive({
  df <- tryCatch(filtered_df(), error = function(e) NULL)
  if (is.null(df) || !nrow(df)) return(NULL)
  d <- tryCatch(attach_channel_colors(resolve_channels(df, as_smlm_settings(input)), smlm_channel_colors(input)),
                error = function(e) NULL)
  if (is.null(d) || !nrow(d) || !"channel" %in% names(d)) return(NULL)
  if (length(unique(d$channel[!is.na(d$channel)])) < 2L) return(NULL)
  d %>%
    dplyr::group_by(tid) %>%
    dplyr::summarise(channel = dplyr::first(channel),
                     channel_label = dplyr::first(channel_label),
                     channel_color = dplyr::first(channel_color),
                     .groups = "drop")
})

# One bar series per channel, in that channel's own colour, replacing the single
# black "Filtered" series. Bars are stacked, so the silhouette still reads as the
# same total distribution the black bars showed.
.sd_add_channel_bars <- function(fig, cur, colname, sh, brk, A, in_rng) {
  cm <- tryCatch(sd_channel_map(), error = function(e) NULL)
  if (is.null(cm) || !"tid" %in% names(cur)) return(NULL)
  j <- dplyr::left_join(cur, cm, by = "tid")
  if (!"channel" %in% names(j) || all(is.na(j$channel))) return(NULL)
  chs <- sort(unique(j$channel[!is.na(j$channel)]))
  if (length(chs) < 2L) return(NULL)

  # STACKED, BUT DRAWN IN OVERLAY MODE - and the distinction is the whole fix.
  #
  # plotly's barmode is a LAYOUT property, so it applies to EVERY bar trace on
  # the panel. Asking for "stack" to stack the channels also stacked the grey
  # raw-data overlay, which put the grey in FRONT of the colours. The panel
  # needs both behaviours at once: grey behind, channels on top of each other.
  #
  # So the stack is built by hand. Draw the cumulative total in the last
  # channel's colour, then the running total one channel shorter over it, down
  # to the first channel. Each bar starts at zero and covers the taller one
  # behind it, so the visible bands are exactly the per-channel counts - a real
  # stack - while barmode stays "overlay" and the grey drawn earlier stays put.
  # Bars are opaque: a translucent stack would blend with the grey behind it and
  # neither colour would be the channel's own any more.
  per_ch <- lapply(chs, function(ch) {
    sub <- j[which(j$channel == ch), , drop = FALSE]
    cnt <- sd_counts(sd_hist_clip(sub[[colname]], sh), brk)
    cnt$counts[!in_rng] <- 0
    cc <- sub$channel_color[1]
    list(counts = cnt$counts,
         color  = if (is.na(cc)) "#666666" else cc,
         label  = sub$channel_label[1])
  })

  for (i in rev(seq_along(per_ch))) {
    cum <- Reduce(`+`, lapply(per_ch[seq_len(i)], `[[`, "counts"))
    fig <- fig %>% plotly::add_bars(
      x = A$centers, y = cum, width = A$width,
      marker = list(color = per_ch[[i]]$color, line = list(width = 0)),
      name = paste0(per_ch[[i]]$label, " (filtered)"), hoverinfo = "skip")
  }
  fig
}

sd_hist_common <- reactive({
  if (!isTRUE(input$sd_hist_shared)) return(NULL)
  all <- tryCatch(sd_all_tbl(),      error = function(e) NULL)
  cur <- tryCatch(sd_filtered_tbl(), error = function(e) NULL)
  if (!is.data.frame(all) || !nrow(all)) return(NULL)

  cols <- c("sd_x", "sd_y")
  if (isTRUE(has_sd_z()) && "sd_z" %in% names(all)) cols <- c(cols, "sd_z")
  cols <- intersect(cols, names(all))
  if (!length(cols)) return(NULL)

  bins <- as.integer(input$sd_bins %||% 60)

  # A user range wins when it is usable; otherwise fall back to the data's own
  # shared span, so clearing the boxes cannot leave the panels blank.
  ur     <- suppressWarnings(as.numeric(input$sd_hist_range))
  manual <- length(ur) == 2 && all(is.finite(ur)) && ur[2] > ur[1]

  brk <- if (manual) seq(ur[1], ur[2], length.out = max(2L, bins + 1L))
         else        sd_breaks(unlist(all[cols], use.names = FALSE), bins)
  if (is.null(brk) || length(brk) < 2) return(NULL)
  sh0 <- list(manual = manual, xrange = range(brk, na.rm = TRUE))

  # Counts have to be recomputed on the SHARED edges to find the common
  # ceiling - the per-axis maxima are on different edges and would not
  # transfer. Raw counts are included whether or not the grey overlay is on,
  # so toggling the overlay cannot rescale the shared axis under the user.
  ymax <- 0
  n_out <- 0L; n_tot <- 0L
  for (cc in cols) {
    v <- all[[cc]]; v <- v[is.finite(v)]
    n_tot <- n_tot + length(v)
    n_out <- n_out + sum(v < sh0$xrange[1] | v > sh0$xrange[2])
    ymax <- max(ymax, sd_counts(sd_hist_clip(all[[cc]], sh0), brk)$counts, na.rm = TRUE)
    if (is.data.frame(cur) && cc %in% names(cur))
      ymax <- max(ymax, sd_counts(sd_hist_clip(cur[[cc]], sh0), brk)$counts, na.rm = TRUE)
  }
  list(brk    = brk,
       manual = manual,
       xrange = sh0$xrange,
       yrange = c(0, if (is.finite(ymax) && ymax > 0) ymax * 1.05 else 1),
       n_out  = n_out, n_tot = n_tot)
})

# A figure must not hide data silently - say how much the chosen range leaves out.
output$sd_hist_range_note <- renderUI({
  sh <- sd_hist_common()
  if (is.null(sh) || !isTRUE(sh$manual)) return(NULL)
  if (sh$n_out <= 0)
    return(helpText("All per-trace spreads fall inside this range."))
  helpText(HTML(sprintf(
    "<b>%s of %s per-trace spread values (%.1f%%) fall outside this range</b> and are not drawn. They are still in the data - this is a display range only, and it does not filter anything.",
    format(sh$n_out, big.mark = ","), format(sh$n_tot, big.mark = ","),
    100 * sh$n_out / max(1L, sh$n_tot))))
})

output$sd_hist_x <- renderPlotly({
  req(input$view_SD)
  all <- sd_all_tbl(); cur <- sd_filtered_tbl()
  shiny::validate(shiny::need(nrow(all) > 0, "No SD data"))

  bins <- as.integer(input$sd_bins %||% 60)
  sh   <- sd_hist_common()
  brk  <- if (!is.null(sh)) sh$brk else sd_breaks(all$sd_x, bins)
  A    <- sd_counts(sd_hist_clip(all$sd_x, sh), brk)
  B    <- sd_counts(sd_hist_clip(cur$sd_x, sh), brk)
  
  # active filter range
  r <- input$sd_x_range %||% c(-Inf, Inf)
  in_rng <- A$centers >= r[1] & A$centers <= r[2]
  
  # clamp filtered counts to 0 outside range (so black bars disappear there)
  B$counts[!in_rng] <- 0
  
  fig <- plotly::plot_ly()
  
  if (isTRUE(input$sd_overlay_all)) {
    fig <- fig %>% plotly::add_bars(
      x = A$centers, y = A$counts, width = A$width,
      # Drawn FIRST so it sits BEHIND the filtered bars, and lightly, so what
      # the filter removed reads as data disappearing rather than as a series
      # in its own right. Same idiom as the metadata histogram.
      marker = list(color = "rgba(160,160,160,0.35)"),
      name = "Raw data", hoverinfo = "skip"
    )
  }
  
  # Black foreground: filtered (B already clamped)
  # One series per channel when the data has channels, one black series when it
  # does not. Same totals either way.
  figc <- .sd_add_channel_bars(fig, cur, "sd_x", sh, brk, A, in_rng)
  fig <- if (!is.null(figc)) figc else fig %>% plotly::add_bars(
    x = A$centers, y = B$counts, width = A$width,
    marker = list(color = "rgba(0,0,0,0.85)"),
    name = "Filtered", hoverinfo = "skip"
  )
  ymax  <- max(c(A$counts, B$counts), na.rm = TRUE)
  # Under a shared axis the median label has to sit inside the SHARED
  # ceiling, or it is drawn off the top of the panel.
  y_top <- if (!is.null(sh)) sh$yrange[2]
           else if (is.finite(ymax) && ymax > 0) ymax * 1.05 else 1
  
  # Filtered median ± nMAD (BLUE)
  if (isTRUE(input$sd_show_mean_sd)) {
    # One median per channel when the data has channels; the pooled median
    # otherwise. A pooled median drawn over channel-coloured bars would be a
    # number belonging to neither of them.
    figm <- .sd_add_channel_medians(fig, cur, "sd_x", y_top,
                                    show_iqr = TRUE)
    if (!is.null(figm)) {
      fig <- figm
    } else {
      muF  <- stats::median(cur$sd_x, na.rm = TRUE)
      labF <- sigma_label_med_pm_sd(cur$sd_x)
      fig  <- add_sigma_vline(fig, muF, y_top, "steelblue", labF)
    # Quartiles flank the median only for the filtered set - drawing them for
    # the raw overlay too would put four thin lines on one panel.
      fig  <- add_sigma_iqr_lines(fig, cur$sd_x, y_top, "steelblue")
    }
  }

  # Raw median ± nMAD (TEAL, dashed)
  if (isTRUE(input$sd_show_overall_mean_sd)) {
    muA  <- stats::median(all$sd_x, na.rm = TRUE)
    labA <- sigma_label_med_pm_sd(all$sd_x)
    fig  <- add_sigma_vline(fig, muA, y_top, "cadetblue", labA, dash = "dash")
  }

  fig %>% plotly::layout(
    barmode = "overlay",
    title = "Trace spread \u03C3 X (nm)",
    xaxis = c(list(title = "\u03C3 X (nm)"), if (!is.null(sh)) list(range = sh$xrange)),
    yaxis = c(list(title = "Count"),         if (!is.null(sh)) list(range = sh$yrange))
  )
})


output$sd_hist_y <- renderPlotly({
  req(input$view_SD)
  all <- sd_all_tbl(); cur <- sd_filtered_tbl()
  shiny::validate(shiny::need(nrow(all) > 0, "No SD data"))
  
  bins <- as.integer(input$sd_bins %||% 60)
  sh   <- sd_hist_common()
  brk  <- if (!is.null(sh)) sh$brk else sd_breaks(all$sd_y, bins)
  A    <- sd_counts(sd_hist_clip(all$sd_y, sh), brk)
  B    <- sd_counts(sd_hist_clip(cur$sd_y, sh), brk)
  
  r <- input$sd_y_range %||% c(-Inf, Inf)
  in_rng <- A$centers >= r[1] & A$centers <= r[2]
  B$counts[!in_rng] <- 0
  
  fig <- plotly::plot_ly()
  
  if (isTRUE(input$sd_overlay_all)) {
    fig <- fig %>% plotly::add_bars(
      x = A$centers, y = A$counts, width = A$width,
      # Drawn FIRST so it sits BEHIND the filtered bars, and lightly, so what
      # the filter removed reads as data disappearing rather than as a series
      # in its own right. Same idiom as the metadata histogram.
      marker = list(color = "rgba(160,160,160,0.35)"),
      name = "Raw data", hoverinfo = "skip"
    )
  }
  
  # Black foreground: filtered (B already clamped)
  # One series per channel when the data has channels, one black series when it
  # does not. Same totals either way.
  figc <- .sd_add_channel_bars(fig, cur, "sd_y", sh, brk, A, in_rng)
  fig <- if (!is.null(figc)) figc else fig %>% plotly::add_bars(
    x = A$centers, y = B$counts, width = A$width,
    marker = list(color = "rgba(0,0,0,0.85)"),
    name = "Filtered", hoverinfo = "skip"
  )
  
  ymax  <- max(c(A$counts, B$counts), na.rm = TRUE)
  # Under a shared axis the median label has to sit inside the SHARED
  # ceiling, or it is drawn off the top of the panel.
  y_top <- if (!is.null(sh)) sh$yrange[2]
           else if (is.finite(ymax) && ymax > 0) ymax * 1.05 else 1
  
  # Filtered median ± nMAD (BLUE)
  if (isTRUE(input$sd_show_mean_sd)) {
    # One median per channel when the data has channels; the pooled median
    # otherwise. A pooled median drawn over channel-coloured bars would be a
    # number belonging to neither of them.
    figm <- .sd_add_channel_medians(fig, cur, "sd_y", y_top,
                                    show_iqr = TRUE)
    if (!is.null(figm)) {
      fig <- figm
    } else {
      muF  <- stats::median(cur$sd_y, na.rm = TRUE)
      labF <- sigma_label_med_pm_sd(cur$sd_y)
      fig  <- add_sigma_vline(fig, muF, y_top, "steelblue", labF)
      fig  <- add_sigma_iqr_lines(fig, cur$sd_y, y_top, "steelblue")
    }
  }

  # Raw median ± nMAD (TEAL, dashed)
  if (isTRUE(input$sd_show_overall_mean_sd)) {
    muA  <- stats::median(all$sd_y, na.rm = TRUE)
    labA <- sigma_label_med_pm_sd(all$sd_y)
    fig  <- add_sigma_vline(fig, muA, y_top, "cadetblue", labA, dash = "dash")
  }

  fig %>% plotly::layout(
    barmode = "overlay",
    title = "Trace spread \u03C3 Y (nm)",
    xaxis = c(list(title = "\u03C3 Y (nm)"), if (!is.null(sh)) list(range = sh$xrange)),
    yaxis = c(list(title = "Count"),         if (!is.null(sh)) list(range = sh$yrange))
  )
})


output$sd_hist_z <- renderPlotly({
  req(input$view_SD, has_sd_z())
  all <- sd_all_tbl(); cur <- sd_filtered_tbl()
  shiny::validate(shiny::need("sd_z" %in% names(all) && any(is.finite(all$sd_z)), "No Z in dataset"))
  
  bins <- as.integer(input$sd_bins %||% 60)
  sh   <- sd_hist_common()
  brk  <- if (!is.null(sh)) sh$brk else sd_breaks(all$sd_z, bins)
  A    <- sd_counts(sd_hist_clip(all$sd_z, sh), brk)
  B    <- sd_counts(sd_hist_clip(cur$sd_z, sh), brk)
  
  r <- input$sd_z_range %||% c(-Inf, Inf)
  in_rng <- A$centers >= r[1] & A$centers <= r[2]
  B$counts[!in_rng] <- 0
  
  fig <- plotly::plot_ly()
  
  if (isTRUE(input$sd_overlay_all)) {
    fig <- fig %>% plotly::add_bars(
      x = A$centers, y = A$counts, width = A$width,
      # Drawn FIRST so it sits BEHIND the filtered bars, and lightly, so what
      # the filter removed reads as data disappearing rather than as a series
      # in its own right. Same idiom as the metadata histogram.
      marker = list(color = "rgba(160,160,160,0.35)"),
      name = "Raw data", hoverinfo = "skip"
    )
  }
  
  # Black foreground: filtered (B already clamped)
  # One series per channel when the data has channels, one black series when it
  # does not. Same totals either way.
  figc <- .sd_add_channel_bars(fig, cur, "sd_z", sh, brk, A, in_rng)
  fig <- if (!is.null(figc)) figc else fig %>% plotly::add_bars(
    x = A$centers, y = B$counts, width = A$width,
    marker = list(color = "rgba(0,0,0,0.85)"),
    name = "Filtered", hoverinfo = "skip"
  )
  
  ymax  <- max(c(A$counts, B$counts), na.rm = TRUE)
  # Under a shared axis the median label has to sit inside the SHARED
  # ceiling, or it is drawn off the top of the panel.
  y_top <- if (!is.null(sh)) sh$yrange[2]
           else if (is.finite(ymax) && ymax > 0) ymax * 1.05 else 1
  
  # Filtered median ± nMAD (BLUE)
  if (isTRUE(input$sd_show_mean_sd)) {
    # One median per channel when the data has channels; the pooled median
    # otherwise. A pooled median drawn over channel-coloured bars would be a
    # number belonging to neither of them.
    figm <- .sd_add_channel_medians(fig, cur, "sd_z", y_top,
                                    show_iqr = TRUE)
    if (!is.null(figm)) {
      fig <- figm
    } else {
      muF  <- stats::median(cur$sd_z, na.rm = TRUE)
      labF <- sigma_label_med_pm_sd(cur$sd_z)
      fig  <- add_sigma_vline(fig, muF, y_top, "steelblue", labF)
      fig  <- add_sigma_iqr_lines(fig, cur$sd_z, y_top, "steelblue")
    }
  }

  # Raw median ± nMAD (TEAL, dashed)
  if (isTRUE(input$sd_show_overall_mean_sd)) {
    muA  <- stats::median(all$sd_z, na.rm = TRUE)
    labA <- sigma_label_med_pm_sd(all$sd_z)
    fig  <- add_sigma_vline(fig, muA, y_top, "cadetblue", labA, dash = "dash")
  }

  fig %>% plotly::layout(
    barmode = "overlay",
    title = "Trace spread \u03C3 Z (nm)",
    xaxis = c(list(title = "\u03C3 Z (nm)"), if (!is.null(sh)) list(range = sh$xrange)),
    yaxis = c(list(title = "Count"),         if (!is.null(sh)) list(range = sh$yrange))
  )
})
