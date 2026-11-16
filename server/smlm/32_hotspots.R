## ─────────────────────────────────────────────────────────────────────────────
## server/smlm/32_hotspots.R  -  panel: hotspots (local density enhancement)
## ─────────────────────────────────────────────────────────────────────────────
## Part of smlm_server(). Not a standalone script: its statements are evaluated in
## the server function's environment, in file-name order (server/parts.R).
##
## Run handler, status, map, ECDF, per-file metric and density plots, table
## and downloads.
##
## INTERFACE (generated: Rscript tools/server_map.R --headers; do not edit by hand)
##   DEFINES for other parts:
##     build_hotspot_density_fig, build_hotspot_ecdf_fig,
##     build_hotspot_map_fig, build_hotspot_metric_fig, spatial_outline_pts
##   READS from other parts:
##     01_view_upload: smlm_file_store
##     20_dataset_groups: file_group_map, series_palette
##     22_spatial_window: spatial_gate_close, spatial_report_windows,
##                        spatial_support_of
##     23_prepared_cache: .stale_badge, prepared_df, resolve_unit
## END INTERFACE
## ─────────────────────────────────────────────────────────────────────────────

## ══════════════════════════════════════════════════════════════════════════
## LOCAL DENSITY & HOTSPOTS
## ══════════════════════════════════════════════════════════════════════════
## Computes, per file: the cell footprint, every unit's local density as a
## MULTIPLE OF THAT CELL'S OWN MEAN, a CSR-calibrated DBSCAN eps, and the foci
## that eps finds. See the block banner "8c. LOCAL DENSITY & HOTSPOTS" in
## R/smlm_hotspots.R for why these statistics and not nearest-neighbour distance.
##
## THE FOOTPRINT IS NOT OPTIONAL HERE, which is the one place this panel departs
## from the rest of the tab. Elsewhere `spatial_use_support` is an opt-in that
## improves an estimate; here every statistic is "per unit of this cell's own
## area", so a bounding box does not lose power, it silently changes the number:
## a cell imaged half in frame reads at half its true density and every
## enhancement doubles. The panel therefore resolves its own footprint with
## `force = TRUE` and says so in the readout.

hotspot_rows   <- reactiveVal(NULL)   # one summary row per file
hotspot_detail <- reactiveVal(NULL)   # per-file pts / support / enh / foci, for the map
hotspot_params <- reactiveVal(NULL)   # resolved params, frozen at Run
hotspot_unit   <- reactiveVal(NULL)   # the unit the stored results were built on


## GOAL: pick the point set that draws the observation window, per the shared
## `dens_support_from` setting. `locs` is the localization table, `units` the
## one-point-per-trace/molecule centroids. See spatial_window_of() for the
## measured trade-off; the analyses themselves always run on the units.
spatial_outline_pts <- function(locs_df, unit_pts) {
  if (identical(input$dens_support_from %||% "units", "locs") &&
      is.data.frame(locs_df) && nrow(locs_df)) return(locs_df)
  if (is.data.frame(unit_pts) && nrow(unit_pts)) unit_pts else locs_df
}



observeEvent(input$run_hotspot, {
  unit    <- resolve_unit(input$hot_unit, "cluster", "tid", "Local density & hotspots")
  ulab    <- hotspot_unit_label(unit)
  use_3d  <- isTRUE(input$hot_use_3d)
  k       <- max(1L, as.integer(input$hot_k %||% 5L))
  quad    <- max(10, as.numeric(input$hot_quadrat %||% 250))
  thresh  <- max(1, as.numeric(input$hot_enh_thresh %||% 5))
  minpts  <- max(2L, as.integer(input$hot_minpts %||% 5L))
  target  <- min(0.5, max(0.001, as.numeric(input$hot_target_frac %||% 0.05)))
  nsim    <- max(3L, as.integer(input$hot_nsim %||% 15L))
  # The window is NOT this panel's to define. Pitch, closing radius and estimator
  # mode all come from "Density & area normalization", so this panel measures the
  # same cell area that Ripley's L and g(r) measure. See spatial_window_of().
  h_nm    <- max(10, as.numeric(input$dens_grid_h %||% 100))
  needs_cluster <- unit == "cluster"

  store <- smlm_file_store()
  shiny::validate(shiny::need(length(store) > 0, "No files loaded."))

  # ---- 1. prepare every file, then take one point per unit -------------------
  items <- lapply(names(store), function(nm) {
    df <- prepared_df(nm, needs_cluster = needs_cluster)
    if (is.null(df) || !nrow(df)) return(NULL)
    pts <- unit_centroids(df, unit = unit, use_3d = use_3d)
    if (is.null(pts) || nrow(pts) < max(3L, minpts)) return(NULL)
    list(tag = nm, df = df, pts = pts)
  })
  items <- Filter(Negate(is.null), items)
  shiny::validate(shiny::need(length(items) > 0,
    "No file had enough units. Check the filters, or lower minPts."))

  # ---- 2. ONE gap-fill radius for the whole set -----------------------------
  # Shared, for the same reason the curve panels share theirs: the footprint sets
  # the area, the area sets the density, and a per-file radius would put each
  # cell on its own density scale - which is precisely the comparison this panel
  # exists to make. Resolved on the UNIT positions (not the raw localizations),
  # because the units are what the density is measured on.
  # THE FOOTPRINT MUST NOT BE BUILT FROM THE UNITS, so the sweep runs on every
  # LOCALIZATION, exactly as the curve panels and the density panel do. The
  # support is a property of where signal was detected; estimating it from unit
  # centroids makes the area depend on the unit count, which is the very thing
  # being normalised out. Measured (see compute_unit_density): two synthetic
  # cells with identical geometry and a 4x labelling difference gave
  # centroid-based supports of 5.01e6 against 3.82e6 nm^2, turning a true 4.0x
  # density ratio into 3.05x.
  cc <- spatial_gate_close(lapply(items, function(it)
          list(df = spatial_outline_pts(it$df, it$pts), tag = it$tag)),
          use_3d, force = TRUE, h_nm = h_nm)
  close_nm <- if (!is.null(cc$close_nm) && is.finite(cc$close_nm)) cc$close_nm else
    max(0, as.numeric(input$dens_close_nm %||% 200))

  # ---- 3. per file: footprint -> enhancement -> calibrated eps -> foci -------
  out <- withProgress(message = sprintf("Local density: %d file(s)...", length(items)),
                      value = 0, {
    lapply(seq_along(items), function(i) {
      it <- items[[i]]
      incProgress(1 / length(items), detail = it$tag)
      # The SAME call the curve panels make, on the SAME point set (all
      # localizations - see the note on the radius above), so the window is
      # logged once and the "these files were not measured the same way" check
      # covers this panel too. The unit centroids are what the STATISTICS run
      # on; they are not what defines the cell.
      sup <- tryCatch(spatial_support_of(spatial_outline_pts(it$df, it$pts),
                                         use_3d, tag = it$tag,
                                         close_nm = close_nm, force = TRUE),
                      error = function(e) NULL)
      hotspot_for_points(it$pts, support = sup, k = k, quadrat_nm = quad,
                         enh_thresh = thresh, minPts = minpts, target_frac = target,
                         nsim = nsim, use_3d = use_3d, tag = it$tag)
    })
  })
  out <- Filter(Negate(is.null), out)
  shiny::validate(shiny::need(length(out) > 0, "Local density returned no data."))
  names(out) <- vapply(out, function(o) o$tag, character(1))

  rows <- do.call(rbind, lapply(out, `[[`, "row"))
  rows$close_nm <- close_nm
  hotspot_rows(rows)
  hotspot_unit(unit)
  hotspot_detail(out)
  hotspot_params(list(unit = input$hot_unit, k = input$hot_k,
                      `quadrat (nm)` = input$hot_quadrat,
                      `enhancement threshold` = input$hot_enh_thresh,
                      minPts = input$hot_minpts,
                      `false-cluster target` = input$hot_target_frac,
                      `CSR draws` = input$hot_nsim,
                      `3-D` = use_3d))

  # The map file picker only ever offers files that actually produced a result.
  spatial_report_windows(vapply(items, function(it) as.character(it$tag),
                                character(1)))

  updateSelectInput(session, "hot_map_file", choices = names(out),
                    selected = if (!is.null(input$hot_map_file) &&
                                   input$hot_map_file %in% names(out))
                      input$hot_map_file else names(out)[1])

  # A bounding-box window makes every number here a per-FRAME density rather than
  # a per-CELL one, which is a different quantity, not a noisier version of the
  # same one. It happens two ways: the estimator mode is set to bbox or hull (no
  # mask), or the footprint guard rejected the window for that file.
  fb <- rows$dataset[grepl("bbox", rows$window_mode %||% "", fixed = TRUE)]
  if (length(fb) == nrow(rows))
    showNotification(HTML(sprintf(
      "<b>Every window is a bounding box</b>, so densities are per imaged frame, not per cell. Set the estimator to <b>grid</b> or <b>kde</b> under Density &amp; area normalization. (%s)",
      rows$window_mode[1])), type = "warning", duration = NULL)
  else if (length(fb))
    showNotification(HTML(sprintf(
      "<b>No footprint for %s.</b> Those files fell back to the bounding box and are not on the same scale as the rest. Raise the gap-fill radius or check the support map.",
      paste(utils::head(fb, 4), collapse = ", "))), type = "warning", duration = NULL)
})

output$hotspot_stale <- renderUI(.stale_badge(
  list(unit = input$hot_unit, k = input$hot_k, `quadrat (nm)` = input$hot_quadrat,
       `enhancement threshold` = input$hot_enh_thresh, minPts = input$hot_minpts,
       `false-cluster target` = input$hot_target_frac, `CSR draws` = input$hot_nsim,
       `3-D` = isTRUE(input$hot_use_3d)),
  hotspot_params()))
# hot_curve_mode and hot_metric are deliberately NOT in the stale badge: they
# re-draw stored results and never change what was computed.

# ── Readout: what the run resolved, in one line ────────────────────────────
output$hotspot_status <- renderUI({
  rows <- hotspot_rows()
  if (is.null(rows) || !nrow(rows))
    return(helpText("Press Run to measure local density and detect hotspots."))
  det <- hotspot_detail()
  floors <- sum(vapply(det, function(o) isTRUE(o$eps_cal$hit_floor), logical(1)))
  tags$div(style = paste0("margin:4px 0 8px; padding:6px 11px; border-left:4px solid #cfe0ee;",
                          "background:#f6f9fc; font-size:0.86em; color:#445;"),
    sprintf("%d file(s), measured on %s. Window: %s, gap fill %.0f nm, pitch %.0f nm (from Density & area normalization). Footprint %.1f-%.1f um2, density %.1f-%.1f per um2.",
            nrow(rows), hotspot_unit_label(hotspot_unit()),
            input$dens_support_mode %||% "grid",
            rows$close_nm[1], max(10, as.numeric(input$dens_grid_h %||% 100)),
            min(rows$area_um2), max(rows$area_um2),
            min(rows$density), max(rows$density)),
    # eps is reported separately because a file whose calibration failed leaves
    # NA here, and min()/max() over an all-NA column returns Inf with a warning.
    if (any(is.finite(rows$eps_nm)))
      sprintf(" Calibrated eps %.0f-%.0f nm.",
              min(rows$eps_nm, na.rm = TRUE), max(rows$eps_nm, na.rm = TRUE))
    else " eps was not calibrated on any file.",
    if (floors > 0) tags$div(style = "color:#a04000; margin-top:4px;",
      sprintf("%d file(s) hit the bottom of the eps grid - a random pattern still clustered at the smallest eps tried, so their focus numbers are not calibrated.", floors)))
})

# ── The map: footprint + units coloured by local-density enhancement ───────
build_hotspot_map_fig <- function() {
  det <- hotspot_detail()
  shiny::validate(shiny::need(!is.null(det) && length(det) > 0, "Run this panel first."))
  fname <- input$hot_map_file %||% names(det)[1]
  shiny::validate(shiny::need(fname %in% names(det), "Pick a dataset."))
  o <- det[[fname]]

  dark   <- isTRUE(input$dark_mode)
  bg     <- if (dark) "#1e1e1e" else "white"
  fgcol  <- if (dark) "#e6e6e6" else "#333a45"
  maskc  <- if (dark) "rgba(255,255,255,0.09)" else "rgba(120,132,148,0.13)"

  pts <- o$pts[o$le$keep, , drop = FALSE]
  enh <- o$enh
  # Log colour scale: enhancement spans two orders of magnitude in a clustered
  # cell, so a linear scale renders every background unit the same colour and
  # only the single densest spot is visible.
  lo  <- log10(max(min(enh, na.rm = TRUE), 0.1))
  hi  <- log10(max(enh, na.rm = TRUE))
  if (!is.finite(lo) || !is.finite(hi) || hi <= lo) { lo <- -1; hi <- 1 }

  fig <- plotly::plot_ly()
  map_shapes <- list(); map_annots <- list()
  # Footprint underneath, as a flat mask - this is the "crop to the cell" the
  # density numbers are actually divided by, so it must be visible, not implied.
  if (!is.null(o$support) && !is.null(o$support$mask)) {
    m2   <- mask_project_2d(o$support$mask)
    cell <- o$support$h
    d    <- dim(m2)
    xc <- o$support$origin[1] + (seq_len(d[1]) - 0.5) * cell
    yc <- o$support$origin[2] + (seq_len(d[2]) - 0.5) * cell
    zin <- ifelse(m2, 1, NA_real_)
    fig <- fig %>% plotly::add_heatmap(
      x = xc, y = yc, z = t(zin), showscale = FALSE, hoverinfo = "skip",
      colorscale = list(list(0, maskc), list(1, maskc)))
  }

  fig <- fig %>% plotly::add_markers(
    x = pts$x, y = pts$y, type = "scatter", mode = "markers",
    marker = list(size = 5, color = log10(pmax(enh, 0.1)),
                  colorscale = plotly_colorscale(input$hot_palette %||% "Turbo"),
                  cmin = lo, cmax = hi, line = list(width = 0),
                  colorbar = list(title = list(text = "local density\n/ cell mean"),
                                  tickvals = log10(c(0.1, 0.3, 1, 3, 10, 30, 100)),
                                  ticktext = c("0.1x", "0.3x", "1x", "3x", "10x", "30x", "100x"),
                                  len = 0.8)),
    text = sprintf("%.1fx cell mean", enh), hoverinfo = "text",
    name = hotspot_unit_label(hotspot_unit()))

  # ── Box the detected foci, at their real size ─────────────────────────
  # A marker was the wrong mark for this. `marker$size` is in PIXELS, so a ring
  # sat at a fixed screen size whatever the focus measured and whatever the zoom:
  # on foci spanning 100-490 nm it read as a dot in the middle of the group
  # rather than as its outline. These are rectangle SHAPES in data coordinates,
  # so each box is the size of the focus it encloses and stays correct on zoom.
  #
  # Every focus is boxed, not only the biggest, and each carries its count. That
  # is deliberate: "largest" ranks by member count while the eye reads DENSITY,
  # and the two disagree in half the two-condition files - a big loose focus outranks a
  # small tight one. Boxing all of them with their counts removes the question of
  # which one the panel picked.
  if (isTRUE(input$hot_show_foci) && !is.null(o$foci) && nrow(o$foci$table)) {
    tb  <- o$foci$table
    # Boxing every focus is unreadable on a busy cell - one two-condition file yields 23.
    # The table is ordered largest first, so taking the head keeps the ones worth
    # looking at. 0 means "all", for when the spatial arrangement of the small
    # ones is the point.
    ntop <- as.integer(input$hot_box_top %||% 5L)
    if (is.finite(ntop) && ntop > 0L) tb <- utils::head(tb, ntop)
    fp  <- o$pts[o$foci$keep, , drop = FALSE]
    cl  <- o$foci$cluster
    ring <- if (dark) "#ffd166" else "#c0392b"
    shp <- list(); ann <- list()
    for (ii in seq_len(nrow(tb))) {
      m <- as.matrix(fp[cl == tb$focus[ii], c("x", "y"), drop = FALSE])
      if (!nrow(m)) next
      rx <- range(m[, 1]); ry <- range(m[, 2])
      # A little air, floored so a very tight focus still gets a visible box.
      pad <- max(0.10 * max(diff(rx), diff(ry)), 20)
      # The biggest focus gets the heavier line, so the table's headline number
      # is findable on the map without counting.
      lw <- if (ii == 1L) 2.2 else 1.2
      shp[[length(shp) + 1L]] <- list(
        type = "rect", xref = "x", yref = "y",
        x0 = rx[1] - pad, x1 = rx[2] + pad,
        y0 = ry[1] - pad, y1 = ry[2] + pad,
        line = list(color = ring, width = lw), fillcolor = "rgba(0,0,0,0)",
        layer = "above")
      ann[[length(ann) + 1L]] <- list(
        x = mean(rx), y = ry[2] + pad, xref = "x", yref = "y",
        text = sprintf("%d", tb$n[ii]), showarrow = FALSE,
        yanchor = "bottom", font = list(size = 11, color = ring))
    }
    map_shapes <- shp
    map_annots <- ann
    # An invisible trace over the focus centroids carries the hover detail; the
    # shapes themselves cannot hold a tooltip.
    fig <- fig %>% plotly::add_markers(
      x = tb$x, y = tb$y, type = "scatter", mode = "markers",
      marker = list(size = 12, color = "rgba(0,0,0,0.01)",
                    line = list(width = 0)),
      text = sprintf("focus %d: %d %s, Rg %.0f nm, aniso %.2f",
                     tb$focus, tb$n, hotspot_unit_label(hotspot_unit()),
                     tb$rg, tb$aniso),
      hoverinfo = "text", name = "foci", showlegend = FALSE)
  }

  plotly::layout(fig,
    paper_bgcolor = bg, plot_bgcolor = bg,
    font = list(color = fgcol),
    showlegend = FALSE,
    title = list(text = sprintf("%s - %d %s, %.1f um2, %.1f/um2",
                                fname, o$row$n_units, hotspot_unit_label(hotspot_unit()),
                                o$row$area_um2, o$row$density),
                 font = list(size = 13)),
    shapes = map_shapes, annotations = map_annots,
    xaxis = list(title = "x (nm)", zeroline = FALSE),
    # scaleanchor keeps the cell shape honest; a stretched footprint reads as a
    # different cell.
    yaxis = list(title = "y (nm)", zeroline = FALSE,
                 scaleanchor = "x", scaleratio = 1))
}
output$hotspot_map <- renderPlotly(build_hotspot_map_fig())

# ── The enhancement distribution, one curve per file ──────────────────────
## Cumulative rather than a histogram: the conditions are indistinguishable over
## the bulk of the distribution and separate only in the upper tail, so a
## cumulative curve shows where they separate without a bin-width choice.
##
## TWO THINGS ABOUT THE x AXIS, both learned from the real data.
##
## 1. THE LOW TAIL IS REAL BUT UNINFORMATIVE, AND IT DOMINATES A LOG AXIS.
##    An isolated unit in a sparse corner has a large k-th neighbour distance, so
##    its enhancement runs down to ~0.05x. Left alone the axis spans nearly four
##    decades and about 40% of the width carries the bottom 1% of the units,
##    which is where the plot "wants to go far below anything reasonable". The
##    lower limit is therefore the smallest 1st percentile across the drawn
##    files, floored at 0.1x. Nothing is dropped from the curves - only the view
##    is clipped - and the floor is stated on the axis title.
##
## 2. THE SEPARATION IS ENTIRELY IN THE UPPER TAIL, which a cumulative curve
##    compresses into its last few percent. The "upper tail" mode plots the
##    SURVIVAL function 1 - F on log-log instead, where "1 unit in 20 is at 10x
##    or denser" is a readable distance between curves rather than two lines
##    converging on 1.
build_hotspot_ecdf_fig <- function() {
  det <- hotspot_detail()
  shiny::validate(shiny::need(!is.null(det) && length(det) > 0, "Run this panel first."))
  dark <- isTRUE(input$dark_mode)
  bg   <- if (dark) "#1e1e1e" else "white"
  fg   <- if (dark) "#e6e6e6" else "#333a45"
  tail_mode <- !identical(input$hot_curve_mode %||% "tail", "cdf")
  ulab <- hotspot_unit_label(hotspot_unit())

  grouped <- (input$analysis_by %||% "file") == "group"
  gm  <- file_group_map()
  # Resolve every file's label FIRST, "Ungrouped" included, and build the palette
  # from that complete set. Building it from the grouped files only left an
  # ungrouped file asking for a colour that had never been created.
  labs <- vapply(names(det), function(nm)
    if (grouped) .name_or(gm, nm, "Ungrouped") else nm, character(1))
  pal <- series_palette(unique(labs))
  thresh <- max(1, as.numeric(input$hot_enh_thresh %||% 5))

  # One shared x window for every curve, so the files are visually comparable.
  p1s <- vapply(det, function(o) {
    e <- o$enh[is.finite(o$enh) & o$enh > 0]
    if (length(e) < 3L) NA_real_ else unname(stats::quantile(e, 0.01))
  }, numeric(1))
  x_lo <- max(0.1, min(p1s, na.rm = TRUE))
  x_hi <- max(vapply(det, function(o) max(o$enh[is.finite(o$enh)], na.rm = TRUE),
                     numeric(1)), na.rm = TRUE)
  if (!is.finite(x_lo) || x_lo <= 0) x_lo <- 0.1
  if (!is.finite(x_hi) || x_hi <= x_lo) x_hi <- x_lo * 100

  fig <- plotly::plot_ly()
  seen <- character(0)
  for (nm in names(det)) {
    e <- sort(det[[nm]]$enh)
    e <- e[is.finite(e) & e > 0]
    if (length(e) < 3L) next
    lab <- unname(labs[nm])
    col <- .name_or(pal, lab, "#888888")
    F_e <- seq_along(e) / length(e)
    if (tail_mode) {
      # Survival, and drop the final point: 1 - F is exactly 0 there, which a log
      # axis cannot draw and which would clip the curve short.
      y <- 1 - F_e
      k <- y > 0
      xx <- e[k]; yy <- y[k]
      htm <- paste0(nm, "<br>%{y:.3f} of ", ulab,
                    " at %{x:.2f}x or denser<extra></extra>")
    } else {
      xx <- e; yy <- F_e
      htm <- paste0(nm, "<br>%{y:.2f} of ", ulab,
                    " below %{x:.2f}x<extra></extra>")
    }
    if (!length(xx)) next
    fig <- fig %>% plotly::add_lines(
      x = xx, y = yy, name = lab,
      legendgroup = lab, showlegend = !(lab %in% seen),
      line = list(color = col, width = 2), hovertemplate = htm)
    seen <- c(seen, lab)
  }

  # 1x (this cell's own average) and the reporting threshold, as references.
  vline <- function(at, col, dash) list(
    type = "line", xref = "x", yref = "paper",
    x0 = log10(at), x1 = log10(at), y0 = 0, y1 = 1,
    line = list(dash = dash, width = 1, color = col))

  # EXPLICIT TICKS. Plotly's own log ticks come out as 5e-02 / 5e-01 / 5e+00,
  # which reads as scientific notation rather than as "half the cell average" /
  # "five times the cell average" - the quantity this axis carries. A 1-2-5
  # ladder labelled with an x is what the number means.
  ladder <- c(0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000)
  tk <- ladder[ladder >= x_lo * 0.999 & ladder <= x_hi * 1.15]
  if (length(tk) < 3L) tk <- ladder
  xax <- list(type = "log", range = c(log10(x_lo), log10(x_hi * 1.15)),
              tickvals = tk,
              ticktext = ifelse(tk < 1, sprintf("%.1fx", tk), sprintf("%gx", tk)),
              title = sprintf("local density / cell mean (view from %.1fx)", x_lo))

  # In tail mode the y axis is also a fraction, so give it the same treatment.
  yfr <- c(0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1)
  yax <- if (tail_mode)
    list(type = "log", tickvals = yfr, ticktext = sprintf("%g", yfr),
         title = sprintf("fraction of %s at least this dense", ulab))
  else
    list(title = sprintf("cumulative fraction of %s", ulab), range = c(0, 1))

  plotly::layout(fig,
    paper_bgcolor = bg, plot_bgcolor = bg, font = list(color = fg),
    title = list(text = sprintf("Local-density enhancement (%s)", ulab),
                 font = list(size = 13)),
    legend = list(orientation = "h", y = -0.22),
    xaxis = xax, yaxis = yax,
    shapes = list(vline(1, "#888888", "dot"), vline(thresh, "#c0392b", "dash")))
}
output$hotspot_ecdf <- renderPlotly(build_hotspot_ecdf_fig())

# ── One metric per file, as a dot plot by group ────────────────────────────
.HOTSPOT_METRICS <- c(
  "Density unevenness (x random)"      = "cv_ratio",
  "CV of local density (raw)"          = "cv_local",
  "Largest focus (x random 95th pct)"  = "max_focus_ratio",
  "Largest focus (n)"                  = "max_focus_n",
  "Enhancement, 90th percentile"       = "enh_q90",
  "Enhancement, 99th percentile"       = "enh_q99",
  "Fraction above the threshold"       = "frac_above",
  "Index of dispersion (quadrats)"     = "vmr",
  "Number of foci"                     = "n_foci",
  "Fraction in foci"                   = "frac_in_foci",
  "Largest focus Rg (nm)"              = "max_focus_rg",
  "Calibrated eps (nm)"                = "eps_nm",
  "Median NN distance (nm)"            = "nn_med",
  "Clark-Evans R"                      = "clark_evans",
  "Density (per um2)"                  = "density")

build_hotspot_metric_fig <- function() {
  rows <- hotspot_rows()
  shiny::validate(shiny::need(!is.null(rows) && nrow(rows) > 0, "Run this panel first."))
  mkey <- input$hot_metric %||% "cv_ratio"
  shiny::validate(shiny::need(mkey %in% names(rows), "Unknown metric."))
  lab  <- names(.HOTSPOT_METRICS)[match(mkey, .HOTSPOT_METRICS)]
  if (!length(lab) || is.na(lab)) lab <- mkey
  dark <- isTRUE(input$dark_mode)
  bg   <- if (dark) "#1e1e1e" else "white"
  fg   <- if (dark) "#e6e6e6" else "#333a45"

  gm  <- file_group_map()
  grp <- unname(gm[rows$dataset]); grp[is.na(grp)] <- "Ungrouped"
  lv  <- unique(grp)
  pal <- series_palette(lv)
  v   <- rows[[mkey]]

  fig <- plotly::plot_ly()
  for (g in lv) {
    s <- grp == g
    # Jitter only on x, so the value stays exactly readable off the y axis.
    fig <- fig %>% plotly::add_markers(
      x = match(g, lv) + stats::runif(sum(s), -0.09, 0.09), y = v[s],
      name = g, marker = list(size = 11, color = .name_or(pal, g, "#888"),
                              line = list(width = 1, color = bg)),
      text = rows$dataset[s],
      hovertemplate = paste0("%{text}<br>", lab, ": %{y:.3g}<extra></extra>"))
    # Median bar: one line, no box - with 4-6 cells a boxplot draws quartiles
    # from too few points to mean anything.
    md <- stats::median(v[s], na.rm = TRUE)
    if (is.finite(md))
      fig <- fig %>% plotly::add_segments(
        x = match(g, lv) - 0.22, xend = match(g, lv) + 0.22, y = md, yend = md,
        line = list(color = fg, width = 3), showlegend = FALSE, hoverinfo = "skip")
  }
  # Three of these metrics have a known value for an unclustered pattern, so
  # draw it. Without the line the reader has no way to tell "clustered" from
  # "large number".
  ref <- switch(mkey, cv_ratio = 1, vmr = 1, clark_evans = 1, max_focus_ratio = 1,
                cv_local = csr_cv_local_density(input$hot_k %||% 5), NULL)
  if (!is.null(ref) && is.finite(ref))
    fig <- fig %>% plotly::add_segments(
      x = 0.4, xend = length(lv) + 0.6, y = ref, yend = ref,
      line = list(color = "#888888", width = 1, dash = "dash"),
      showlegend = FALSE, hoverinfo = "skip") %>%
      plotly::add_annotations(x = length(lv) + 0.55, y = ref, text = "random",
        showarrow = FALSE, xanchor = "right", yanchor = "bottom",
        font = list(size = 10, color = "#888888"))
  plotly::layout(fig,
    paper_bgcolor = bg, plot_bgcolor = bg, font = list(color = fg),
    title = list(text = lab, font = list(size = 13)), showlegend = FALSE,
    xaxis = list(title = "", tickvals = seq_along(lv), ticktext = lv,
                 range = c(0.4, length(lv) + 0.6)),
    yaxis = list(title = lab))
}
output$hotspot_metric <- renderPlotly(build_hotspot_metric_fig())

# ── The metric-versus-density plot ─────────────────────────────────────────
## THE POINT OF THIS PLOT: with transient transfection every cell carries a
## different plasmid load, so a metric that rises with density is reporting
## expression, not organisation. Plot it and look. If the groups separate
## VERTICALLY at matched density the effect is real; if both groups lie on one
## line the metric is a density readout. Median NN does exactly that, which is
## why it is in the metric list - it is the worked example of the failure.
build_hotspot_density_fig <- function() {
  rows <- hotspot_rows()
  shiny::validate(shiny::need(!is.null(rows) && nrow(rows) > 0, "Run this panel first."))
  mkey <- input$hot_metric %||% "cv_ratio"
  # A metric restored from an older settings file may name a column this build no
  # longer has; without the check the plot dies on rows[[mkey]] instead of saying so.
  shiny::validate(shiny::need(mkey %in% names(rows), "Unknown metric."))
  lab  <- names(.HOTSPOT_METRICS)[match(mkey, .HOTSPOT_METRICS)]
  if (!length(lab) || is.na(lab)) lab <- mkey
  dark <- isTRUE(input$dark_mode)
  bg   <- if (dark) "#1e1e1e" else "white"
  fg   <- if (dark) "#e6e6e6" else "#333a45"

  gm  <- file_group_map()
  grp <- unname(gm[rows$dataset]); grp[is.na(grp)] <- "Ungrouped"
  lv  <- unique(grp); pal <- series_palette(lv)
  v   <- rows[[mkey]]

  fig <- plotly::plot_ly()
  for (g in lv) {
    s <- grp == g
    fig <- fig %>% plotly::add_markers(
      x = rows$density[s], y = v[s], name = g,
      marker = list(size = 11, color = .name_or(pal, g, "#888"),
                    line = list(width = 1, color = bg)),
      text = rows$dataset[s],
      hovertemplate = "%{text}<br>%{x:.1f}/um2<br>%{y:.3g}<extra></extra>")
    # Trend line per group, drawn only with 3+ cells - with two points a line is
    # not a trend, it is a line through two points.
    if (sum(s & is.finite(v)) >= 3L) {
      ld <- log10(rows$density[s]); yy <- v[s]
      ok <- is.finite(ld) & is.finite(yy)
      if (sum(ok) >= 3L && stats::sd(ld[ok]) > 0) {
        cf <- stats::coef(stats::lm(yy[ok] ~ ld[ok]))
        xr <- range(rows$density[s][ok])
        xs <- seq(xr[1], xr[2], length.out = 20)
        fig <- fig %>% plotly::add_lines(
          x = xs, y = cf[1] + cf[2] * log10(xs), showlegend = FALSE,
          line = list(color = .name_or(pal, g, "#888"), width = 1, dash = "dash"),
          hoverinfo = "skip")
      }
    }
  }
  ref <- switch(mkey, cv_ratio = 1, vmr = 1, clark_evans = 1, max_focus_ratio = 1,
                cv_local = csr_cv_local_density(input$hot_k %||% 5), NULL)
  if (!is.null(ref) && is.finite(ref)) {
    xr <- range(rows$density[is.finite(rows$density)])
    fig <- fig %>% plotly::add_lines(
      x = xr, y = c(ref, ref), showlegend = FALSE, hoverinfo = "skip",
      line = list(color = "#888888", width = 1, dash = "dash"))
  }
  plotly::layout(fig,
    paper_bgcolor = bg, plot_bgcolor = bg, font = list(color = fg),
    title = list(text = sprintf("%s vs labelling density", lab), font = list(size = 13)),
    xaxis = list(title = sprintf("%s per um2 of footprint",
                                 hotspot_unit_label(hotspot_unit())), type = "log"),
    yaxis = list(title = lab))
}
output$hotspot_vs_density <- renderPlotly(build_hotspot_density_fig())

# ── The per-file table ─────────────────────────────────────────────────────
output$hotspot_table <- renderTable({
  rows <- hotspot_rows()
  shiny::validate(shiny::need(!is.null(rows) && nrow(rows) > 0, "Run this panel first."))
  gm <- file_group_map()
  data.frame(
    Dataset   = rows$dataset,
    Group     = ifelse(is.na(gm[rows$dataset]), "-", gm[rows$dataset]),
    N         = rows$n_units,
    `Area um2` = round(rows$area_um2, 2),
    `Density`  = round(rows$density, 1),
    `Unevenness` = round(rows$cv_ratio, 2),
    `CV local` = round(rows$cv_local, 2),
    `Enh q90`  = round(rows$enh_q90, 1),
    `Frac above` = sprintf("%.1f%%", 100 * rows$frac_above),
    `Dispersion` = round(rows$vmr, 2),
    `eps nm`   = round(rows$eps_nm),
    `Foci`     = rows$n_foci,
    `Largest`  = rows$max_focus_n,
    `Largest / random` = round(rows$max_focus_ratio %||% NA_real_, 2),
    `Rg nm`    = round(rows$max_focus_rg),
    `Aniso`    = round(rows$max_focus_aniso, 2),
    `NN med`   = round(rows$nn_med),
    `Clark-Evans` = round(rows$clark_evans, 2),
    check.names = FALSE, stringsAsFactors = FALSE)
}, striped = TRUE, hover = TRUE, spacing = "xs", width = "100%")

output$download_hotspot <- downloadHandler(
  filename = function() paste0("local_density_hotspots_",
                               format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
  content = function(file) {
    rows <- hotspot_rows()
    shiny::validate(shiny::need(!is.null(rows) && nrow(rows) > 0, "No results."))
    gm <- file_group_map()
    rows$group <- unname(gm[rows$dataset])
    write.csv(rows, file, row.names = FALSE)
  })

# Per-focus detail: one row per detected hotspot across every file, which is what
# you need to plot a focus-size distribution rather than one number per cell.
output$download_hotspot_foci <- downloadHandler(
  filename = function() paste0("hotspot_foci_",
                               format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
  content = function(file) {
    det <- hotspot_detail()
    shiny::validate(shiny::need(!is.null(det) && length(det) > 0, "No results."))
    gm <- file_group_map()
    out <- do.call(rbind, lapply(names(det), function(nm) {
      tb <- det[[nm]]$foci$table
      if (is.null(tb) || !nrow(tb)) return(NULL)
      tb$dataset <- nm
      tb$group   <- unname(gm[nm]) %||% NA_character_
      tb$eps_nm  <- det[[nm]]$eps_cal$eps %||% NA_real_
      tb
    }))
    shiny::validate(shiny::need(!is.null(out) && nrow(out) > 0, "No foci detected."))
    write.csv(out, file, row.names = FALSE)
  })
