# ==========================================================================
# QC 1 - Localization precision cloud (within-trace scatter about the centroid)
# ==========================================================================
# GOAL: the direct measure of localization PRECISION. Pool every localization
# minus its OWN trace centroid, recentre at (0,0), and look at the scatter.
# Its RMS radius is σ_loc. The SHAPE flags problems: elongation = drift within a
# trace or direction-dependent estimator precision; heavy tails = mislinks /
# double events; skew = a systematic bias. Uses raw uploaded data (before
# metadata filters), so it reports the precision as acquired.
#
# NOT a PSF width. MINFLUX fits no PSF - the doughnut is the EXCITATION
# pattern and position comes from photon-count ratios about its central zero,
# so this scatter is not a spot width and does not scale like one.
prec_cloud_results <- reactiveVal(NULL)
observeEvent(input$run_prec_cloud, {
  run_set <- smlm_settings_from_input(input)   # one settings snapshot for this run
  store <- smlm_file_store()
  if (!length(store)) { showNotification("No files loaded.", type = "warning"); return() }
  # Precision is ALWAYS assessed on the single ACTIVE dataset - never pooled across
  # files (an aggregate σ would blur per-acquisition precision). Restrict to it here.
  af <- input$active_file
  if (is.null(af) || !af %in% names(store)) af <- names(store)[[1]]
  store <- store[af]
  min_loc <- max(3L, as.integer(input$prec_cloud_min_loc %||% 5L))
  use_z   <- isTRUE(input$prec_cloud_use_z)
  # Filtered on the SAME localizations the Display tab plots (metadata, trace
  # length and spread filters, fiducial exclusion). If filtering fails the cloud
  # falls back to the unfiltered file, while the drift estimate skips it.
  filt <- lapply(names(store), function(fname)
    tryCatch(apply_meta_filters_batch(qc_exclude_from(store[[fname]], fname), run_set),
             error = function(e) NULL))
  names(filt) <- names(store)
  dfs <- lapply(names(store), function(fname) filt[[fname]] %||% store[[fname]])
  names(dfs) <- names(store)
  res <- withProgress(message = "Computing precision cloud...", value = 0,
    precision_cloud(dfs, min_loc = min_loc, use_z = use_z, drift_dfs = filt,
                    progress = function(frac, detail) incProgress(frac, detail = detail)))
  if (is.null(res)) {
    showNotification("No traces with enough localizations (lower 'Min localizations / trace').",
                     type = "warning"); prec_cloud_results(NULL); return()
  }
  prec_cloud_results(res)
  showNotification(sprintf("Precision cloud: %s localizations; σ_xy = %.2f nm.",
                           format(res$n_locs, big.mark = ","), res$sigma_xy),
                   type = "message", duration = 4)
})

output$prec_cloud_plot <- renderPlotly({
  r <- prec_cloud_results()
  shiny::validate(shiny::need(!is.null(r), "Press 'Run' to compute the precision cloud."))

  # Each axis pair gets its own SQUARE panel: density heatmap + 1σ/2σ ellipses, a
  # faint zero crosshair through (0,0), no background grid, and a shared ±L nm frame.
  # The panels are built manually on explicit independent axes (x/y, x2/y2, x3/y3)
  # rather than plotly::subplot() - subplot does not reliably re-map per-panel
  # scaleanchor, which cross-linked the zoom across panels. Here each y_i is scale-
  # anchored to its own x_i, so aspect is equal (1 nm x == 1 nm y) AND zooming one
  # panel leaves the others untouched. Ellipse semi-axes are the per-axis robust σ.
  # Colour encodes the PROBABILITY MASS, line style encodes data-vs-model:
  # solid = where that fraction of the localizations actually is, dotted =
  # where it would be if the cloud were Gaussian. The two coincide on clean
  # data and separate where the model fails, which is the entire diagnostic.
  c_50 <- "#c03030"; c_90 <- "#2c7fb8"   # 50% red, 90% blue
  c_rb <- "#1a9850"                      # per-trace robust fit = green
  # TWO KINDS OF RING, EACH LABELLED AS WHAT IT IS
  #   SOLID  = where the localizations actually are: the contours enclosing
  #            50% and 90% of them. Round masses because a reader can act on
  #            them.
  #   DOTTED = the fitted Gaussian at 1 and 2 sigma. Semi-axes come from the
  #            eigen-decomposition of the covariance, so each axis gets its own
  #            sigma and the ring is a rotated ELLIPSE, not a circle - which is
  #            what makes its tilt readable as drift.
  #
  #   The two families are NOT the same probability mass: a 2D Gaussian holds
  #   1 - exp(-k^2/2) inside k sigma, i.e. 39.3% at 1 sigma and 86.5% at 2. So
  #   on a perfectly Gaussian cloud the solid 50% ring sits ~18% outside the
  #   dotted 1 sigma ring, and the solid 90% ~7% outside the dotted 2 sigma.
  #   Judge SHAPE and TILT by these rings, not a millimetre-exact overlap; the
  #   tail comparison is carried numerically by the kurtosis in the stats line.
  P_IN <- c(0.50, 0.90)
  # Legend names, also used as legendgroup keys - keep them short and literal.
  NM_50  <- "50 % of locs";  NM_1S <- "1σ (Gaussian fit)"
  NM_90  <- "90 % of locs";  NM_2S <- "2σ (Gaussian fit)"
  NM_TR  <- "1σ per trace"
  has_z <- isTRUE(r$has_z) && any(is.finite(r$dz))
  # Sp = pooled covariance (the population the heatmap shows)
  # Sr = robust per-trace covariance (the headline σ)
  panels <- if (has_z) list(
    list(a = r$dx, b = r$dy, Sp = r$cov_pooled$xy, Sr = r$cov_trace$xy, la = "Δx (nm)", lb = "Δy (nm)"),
    list(a = r$dx, b = r$dz, Sp = r$cov_pooled$xz, Sr = r$cov_trace$xz, la = "Δx (nm)", lb = "Δz (nm)"),
    list(a = r$dy, b = r$dz, Sp = r$cov_pooled$yz, Sr = r$cov_trace$yz, la = "Δy (nm)", lb = "Δz (nm)")
  ) else list(
    list(a = r$dx, b = r$dy, Sp = r$cov_pooled$xy, Sr = r$cov_trace$xy,
         la = "Δx from trace centroid (nm)", lb = "Δy from trace centroid (nm)"))
  # ── Frame half-width, derived from the data ───────────────────────────────
  # This used to be a hard-coded 30 nm, which only frames the plot correctly
  # when σ happens to be ~10 nm. On any dataset with a wider cloud the heatmap
  # filled the whole panel and the 1σ/2σ rings ran clean off it, leaving stray
  # arcs cutting across the corners with no visible centre - the plot looked
  # broken and over-zoomed because it WAS over-zoomed.
  # The frame follows the CORE - the robust, one-vote-per-trace σ - at a plain
  # 3σ. It must not be driven by the POOLED covariance: σ varies between
  # traces (photon count), so the pooled ellipse is pulled outward by a
  # minority of wide traces and sizing the frame to it leaves the visible
  # cloud small in the middle with only the zero crosshair reaching the edges.
  # A clean Gaussian therefore lands at exactly ±3σ and fills the panel; the
  # 50%/90% contours sit at 1.18σ/2.15σ, comfortably inside.
  # The second term is the safety valve: the radius that actually contains the
  # outer mass, so a genuinely wide or heavy-tailed cloud still gets its outer
  # contour inside the frame instead of having it
  # silently dropped. It is a quantile, not a moment, so outliers cannot run
  # away with it. ONE L is shared by all panels so they stay comparable.
  # Each panel is framed on ITS OWN pair, not on a single shared L. In 3D
  # σ_z is routinely 2-3x σ_xy, and one shared frame is then sized by z -
  # which leaves the Δx-Δy cloud as a small dot in the middle of an empty
  # panel. The axis numbers, which are printed per panel, carry the
  # cross-panel comparison instead.
  # NOTE no rounding to "nice" numbers: plotly chooses tick POSITIONS inside
  # whatever range it is given, so rounding the range buys nothing and only
  # overshoots (36 -> 40 pushes the cloud back toward the middle third).
  ring_reach <- function(S, k) {          # semi-MAJOR bounds the ring on both axes
    e <- cov_ellipse_params(S, k = k)
    if (is.null(e)) NA_real_ else e$ra
  }
  panel_lim <- function(P) {
    need <- c(3.0 * ring_reach(P$Sr, 1),
              1.15 * stats::quantile(sqrt(P$a^2 + P$b^2), P_IN[2],
                                     na.rm = TRUE, names = FALSE))
    need <- need[is.finite(need)]
    core <- max(if (length(need)) signif(max(need), 3) else 30, 5)  # never over-zoom
    # The 2σ model ring may reach past a frame sized to the CORE, and when it
    # does it was being dropped and the panel silently lost a legend entry.
    # Let the frame grow to admit it, but CAP that growth: sizing the frame to
    # the pooled fit without a cap is what once collapsed the cloud to a dot
    # and washed the heatmap out, because a minority of wide traces can
    # inflate the pooled covariance without limit. Beyond the cap the ring is
    # still dropped - a ring that far out cannot be shown without wrecking the
    # view of the data it is meant to describe.
    r2 <- ring_reach(P$Sp, 2)
    if (is.finite(r2) && r2 > core) core <- min(1.05 * r2, 1.4 * core)
    core
  }
  Lv <- vapply(panels, panel_lim, numeric(1))

  np  <- length(panels)
  gap <- if (np > 1) 0.06 else 0          # horizontal gap between panels (paper frac)
  w   <- (1 - gap * (np - 1)) / np        # each panel's width (paper frac)

  fig <- plotly::plot_ly()
  shapes <- list(); axis_args <- list()
  for (i in seq_len(np)) {
    P <- panels[[i]]
    L <- Lv[i]                            # this panel's own half-range
    xid <- if (i == 1) "x" else paste0("x", i)
    yid <- if (i == 1) "y" else paste0("y", i)
    keep <- is.finite(P$a) & is.finite(P$b)
    # ── Density: kernel estimate, not a binned histogram ───────────────────
    # A 2D histogram estimates every bin independently, so each carries its
    # own full Poisson noise: at 60x60 bins and ~40k points that is ~9% at the
    # peak but ~26% at the 2σ radius and ~90% at 3σ. The outer contours - the
    # ones the eye reads as "the shape of the cloud" - were therefore mostly
    # shot noise, which is exactly the ragged, under-smoothed look. Re-binning
    # cannot fix that (coarser blurs the core, finer worsens the tails); a
    # kernel estimate can, because each point informs a neighbourhood.
    # zsmooth = "best" then interpolates the grid so no cell edges show.
    s_rob <- if (!is.null(P$Sr)) sqrt(mean(diag(P$Sr))) else NA_real_
    dh <- density_hdr_2d(P$a[keep], P$b[keep], lim = L, sigma = s_rob, probs = P_IN)
    if (!is.null(dh)) {
      fig <- fig %>% plotly::add_heatmap(
        x = dh$x, y = dh$y, z = t(dh$z),      # kde2d is z[x, y]; plotly wants z[y, x]
        xaxis = xid, yaxis = yid, colorscale = "Hot", reversescale = TRUE,
        showscale = FALSE, zsmooth = "best", hoverinfo = "skip", showlegend = FALSE)
    }
    # Ellipses are drawn from explicit outline coordinates: plotly's "circle"
    # shape is defined by a bounding box and can only ever be axis-aligned,
    # which made a diagonally elongated cloud look circular.
    # Rings are TRACES, not layout shapes, so the legend can switch them off.
    # legendgroup ties the three panels together (one click hides a ring
    # everywhere) and only panel 1 contributes the legend entry.
    # add_PATHS, never add_lines: add_lines() SORTS the points by x before
    # connecting them. These are closed curves, so sorting turns every ring
    # into a zigzag that fills its own interior - the rings render as solid
    # blobs and the dashed ones as scattered stipple. add_paths() keeps the
    # order the points were traced in.
    ring <- function(x, y, colr, dsh, nm) {
      fig <<- plotly::add_paths(fig, x = x, y = y, xaxis = xid, yaxis = yid,
        name = nm, legendgroup = nm, showlegend = (i == 1L),
        hoverinfo = "skip", inherit = FALSE,
        line = c(list(color = colr, width = 2), if (nzchar(dsh)) list(dash = dsh)))
    }
    # A ring whose major axis leaves the frame is DROPPED, not drawn. Drawn, it
    # appears only as arcs cutting the corners with no visible centre, which
    # reads as a broken plot; and the frame must stay sized to the core, so
    # growing it to fit is not the alternative. The stats line under the plot
    # already names this case ("pooled σ >> per-trace σ").
    ell <- function(S, k, colr, dsh, nm) {
      e <- cov_ellipse_params(S, k = k)
      if (is.null(e) || !is.finite(e$ra) || e$ra > L) return(invisible(NULL))
      xy <- ellipse_xy(e$ra, e$rb, e$theta)
      if (is.null(xy)) return(invisible(NULL))
      ring(xy$x, xy$y, colr, dsh, nm)
    }
    axl <- function(x0, x1, y0, y1) list(type = "line", xref = xid, yref = yid,
      x0 = x0, x1 = x1, y0 = y0, y1 = y1, line = list(color = "rgba(0,0,0,0.4)", width = 1))
    # Empirical highest-density contours, drawn SOLID at the same two
    # probability masses the dotted Gaussian rings enclose. Same mass, same
    # colour, different line style: wherever a solid line sits outside its
    # dotted partner, the cloud has more localizations out there than a
    # Gaussian of the fitted width would put there. A count-level contour
    # could not be compared to anything.
    # Match each empirical contour to its mass by VALUE, not by a hard-coded
    # midpoint - with the inner mass at 0.50 a "prob < 0.5" test would colour
    # it as the outer one. A mass can yield several disjoint contours; they
    # share one legendgroup so they toggle together, and only the first
    # carries the legend entry.
    seen <- character(0)
    for (q in dh$paths %||% list()) {
      inner <- q$prob <= P_IN[1] + 1e-9
      nm <- if (inner) NM_50 else NM_90
      show <- (i == 1L) && !(nm %in% seen); seen <- c(seen, nm)
      fig <- plotly::add_paths(fig, x = q$x, y = q$y, xaxis = xid, yaxis = yid,
        name = nm, legendgroup = nm, showlegend = show, hoverinfo = "skip",
        inherit = FALSE,
        line = list(color = if (inner) c_50 else c_90, width = 2))
    }
    ell(P$Sp, 1, c_50, "dot",  NM_1S)   # 1σ ellipse, pooled fit
    ell(P$Sp, 2, c_90, "dot",  NM_2S)   # 2σ ellipse, pooled fit
    # Same 1σ, but fitted with one vote per trace instead of one per
    # localization: where the green sits well inside the dotted red, the
    # pooled width is being carried by a minority of wide or long traces.
    ell(P$Sr, 1, c_rb, "dash", NM_TR)
    shapes <- c(shapes, list(axl(0, 0, -L, L), axl(-L, L, 0, 0)))
    dom <- c((i - 1) * (w + gap), (i - 1) * (w + gap) + w)
    # Anchor the X axis to Y and let X give up DOMAIN (width) to satisfy the
    # 1:1 scale. Anchoring the other way round let plotly satisfy the ratio by
    # widening X's RANGE instead, which produced panels spanning x -20..40
    # against y -20..20: asymmetric, off-centre, and 1.5x wider than tall. An
    # ellipse read against axes like that looks rotated even when it is not.
    # Both axes now show exactly -L..L, so the data window really is square
    # and a tilt on screen is a tilt in the data.
    axis_args[[if (i == 1) "xaxis" else paste0("xaxis", i)]] <- list(
      title = P$la, range = c(-L, L), zeroline = FALSE, showgrid = FALSE,
      domain = dom, anchor = yid, scaleanchor = yid, scaleratio = 1,
      constrain = "domain")
    # Y carries the range and the full height; X is the axis that yields (above).
    # Only ONE of the pair may declare scaleanchor - declaring it on both makes
    # the constraint circular and plotly resolves it by widening a range.
    axis_args[[if (i == 1) "yaxis" else paste0("yaxis", i)]] <- list(
      title = P$lb, range = c(-L, L), zeroline = FALSE, showgrid = FALSE,
      domain = c(0, 1), anchor = xid)
  }
  # The rings are real traces now, so the legend entries are the rings
  # themselves: clicking one hides it across every panel (they share a
  # legendgroup), and double-click isolates it. Solid is where the
  # localizations are, dotted is where a Gaussian of the fitted width would
  # put the same fraction, dashed gives every trace one vote instead of every
  # localization.
  p <- do.call(plotly::layout, c(list(fig, shapes = shapes,
    legend = list(orientation = "h", x = 0.5, xanchor = "center",
                  y = -0.24, yanchor = "top",
                  itemclick = "toggle", itemdoubleclick = "toggleothers"),
    margin = list(t = 10, b = 100)), axis_args))
  # Autoscale-to-data interacts with scaleanchor to blow the ranges up (±100) and
  # desync the panel sizes. Make double-click / reset return to the fixed ±L frame
  # instead, and drop the autoscale button so the only "fit" action is that reset -
  # the panels then stay identical squares. Box-zoom (per panel) still works.
  plotly::config(p, doubleClick = "reset", displaylogo = FALSE,
                 modeBarButtonsToRemove = list("autoScale2d"))
})

output$prec_cloud_stats <- renderUI({
  r <- prec_cloud_results(); if (is.null(r)) return(NULL)
  warn <- character(0)
  # Keyed off the EIGENVALUE ratio, not max(σx,σy)/min(σx,σy): the latter scores
  # ~1.0 on a 45°-tilted cloud and would call drifted data isotropic. The angle
  # is reported because drift has a direction and that is the useful clue.
  if (is.finite(r$aniso_ev) && r$aniso_ev > 1.3)
    warn <- c(warn, sprintf("elongated %.2f:1 along %.0f° (drift within traces, or direction-dependent estimator precision)",
                            r$aniso_ev, r$tilt_deg))
  # Tails are judged against what the per-trace σ SPREAD already predicts, not
  # against 3. A pooled cloud of Gaussians with different widths is
  # leptokurtic by construction, and MINFLUX σ varies with photon count, so
  # "kurtosis > 5" fires on clean data (7-15 on artifact-free simulations at a
  # realistic σ spread). Only the excess over the predicted value is evidence.
  # Both an absolute and a relative margin must be cleared. The relative one
  # carries the very-short-trace case: at 5 localizations per trace the
  # per-trace variances have 4 degrees of freedom, so the prediction is itself
  # noisy and a bare "excess > 2" still misfires on clean simulated data.
  kobs <- suppressWarnings(max(r$kurt_x, r$kurt_y, na.rm = TRUE))
  if (is.finite(kobs) && is.finite(r$kurt_exp) &&
      kobs - r$kurt_exp > 2 && kobs > 1.5 * r$kurt_exp)
    warn <- c(warn, sprintf("heavy tails beyond what the σ spread explains (kurtosis %.1f vs %.1f expected - possible mislinks / double events)",
                            kobs, r$kurt_exp))
  # If the length-weighted pooled σ is much larger than the per-trace median, a
  # few long/wide traces (fiducials!) dominate - flag it and point at the scan.
  if (is.finite(r$sigma_xy_med) && r$sigma_xy > 1.5 * r$sigma_xy_med)
    warn <- c(warn, sprintf("pooled σ (%.1f) ≫ per-trace σ (%.1f): the pooled number is carried by a minority of wide or long traces%s",
                            r$sigma_xy, r$sigma_xy_med,
                            # Both a fiducial AND ordinary photon-count spread
                            # produce this gap; naming only fiducials sent the
                            # user hunting for contamination that need not exist.
                            if (is.finite(r$sigma_spread) && r$sigma_spread > 0.5)
                              " - σ also varies widely between traces (photon-count spread), so check that before assuming contamination"
                            else " - run the Fiducial scan"))
  htmltools::div(style = "font-size:0.9em; margin-top:4px;",
    htmltools::tags$b(sprintf("σ_xy = %.2f nm (median per trace, robust)", r$sigma_xy_med)),
    if (r$has_z && is.finite(r$sigma_z_med)) sprintf(" · σ_z %.2f nm", r$sigma_z_med) else "",
    sprintf("  |  elongation %.2f:1 at %.0f° (ρ %.2f)  |  pooled σ_xy %.2f (σ_x %.2f · σ_y %.2f · axis ratio %.2f)  |  σ varies %.0f%% across traces (IQR/median)  |  skew x/y %.2f/%.2f  ·  kurtosis x/y %.1f/%.1f (%.1f expected from that σ spread)  |  %s locs, %s traces",
            r$aniso_ev, r$tilt_deg, r$rho_xy,
            r$sigma_xy, r$sigma_x, r$sigma_y, r$anisotropy,
            100 * (r$sigma_spread %||% NA_real_), r$skew_x, r$skew_y,
            r$kurt_x, r$kurt_y, r$kurt_exp %||% NA_real_,
            format(r$n_locs, big.mark = ","),
            format(r$n_traces, big.mark = ",")),
    if (length(warn)) htmltools::div(style = "color:#b26a00; margin-top:2px;",
      "⚠ ", paste(warn, collapse = "; "), ".") else
      htmltools::div(style = "color:#2a6; margin-top:2px;", "Cloud is compact and near-isotropic."),
    # ── Sample stability, one line per file ─────────────────────────────────
    # Drift accumulates with time and has a direction; jitter does not. They
    # are reported separately because they say different things about a
    # sample, and because jitter - not drift - is what sets how far apart two
    # blinks of one molecule land.
    if (length(r$drift %||% list())) {
      htmltools::tagList(
        htmltools::div(style = "margin-top:6px; font-weight:600;", "Sample stability"),
        lapply(names(r$drift), function(fn) {
          d <- r$drift[[fn]]
          txt <- sprintf("%s: drift %.1f nm over %.0f min (95%% CI %.1f-%.1f), %.1f nm across a typical %.0f s gap between blinks; jitter %.1f nm; %.0f%% of pair displacement is directional. From %s pairs%s.",
                         fn, d$drift_run, d$span / 60, d$drift_run_lo, d$drift_run_hi,
                         d$drift_typ, d$dt_typ, d$jitter, 100 * (d$frac_directed %||% NA_real_),
                         format(d$n_pairs, big.mark = ","),
                         if (is.finite(d$n_units)) sprintf(" in %s multi-blink molecules",
                                                           format(d$n_units, big.mark = ",")) else "")
          # Highlighted only when the drift is big enough to MATTER: a few tens
          # of nm across a whole run is excellent stability, and colouring it
          # as a warning would train the user to ignore the line. The bar is
          # the drift across a typical gap between blinks exceeding the
          # localization spread itself.
          htmltools::div(style = if (isTRUE(d$significant) &&
                                     isTRUE(d$drift_typ > (r$sigma_xy_med %||% Inf)))
                           "color:#b26a00;" else "color:#555;", txt)
        }),
        htmltools::div(style = "color:#777; font-size:0.85em;",
          "Only molecules whose blinks were merged can be seen here, so drift large enough to break a merge is understated. This is a sample readout and is not used to set the de-blink threshold.")
      )
    } else NULL)
})


# ==========================================================================
# SAMPLE DRIFT PLOT - displacement from each molecule's FIRST blink
# ==========================================================================
# GOAL: show, rather than assert, whether the sample moved. Every molecule
# that blinked more than once is placed at the origin at its first blink; each
# later blink is one point, positioned relative to that origin and coloured by
# the gap since it.
#
# HOW TO READ IT: a cloud centred on the origin is a still sample - the
# scatter is jitter, which does not accumulate. Points marching off in one
# direction, further the longer the gap, are drift. The two are hard to tell
# apart from a single number, which is the whole reason for the plot.
#
# The binned means carry the reading. Averaging kills the isotropic part, so
# under pure jitter every bin sits on the origin whatever its time gap, while
# under drift the bins step outward in order. That is a far stronger cue than
# the raw scatter, which drift barely changes when jitter dominates.
output$prec_drift_plot <- renderPlotly({
  r <- prec_cloud_results()
  validate(need(!is.null(r), "Run the precision cloud first."))
  dl <- r$drift %||% list()
  validate(need(length(dl) > 0, "Sample drift could not be measured (needs traces with times, and molecules that blinked more than once)."))
  st <- do.call(rbind, lapply(names(dl), function(fn) {
    s <- dl[[fn]]$steps; if (is.null(s) || !nrow(s)) return(NULL)
    s$file <- fn
    # Molecule ids are per-file integers, so two files would otherwise share
    # ids and be drawn as one molecule.
    if (!is.null(s$mol)) s$mol <- paste(fn, s$mol, sep = " / ")
    s
  }))
  validate(need(!is.null(st) && nrow(st) >= 5, "Too few multi-blink molecules to plot."))

  ax  <- intersect(c("dx", "dy", "dz"), names(st))
  lab <- c(dx = "Δx", dy = "Δy", dz = "Δz")
  pairs_ax <- if (length(ax) >= 3) list(c("dx","dy"), c("dx","dz"), c("dy","dz")) else list(ax[1:2])

  # One frame for all panels so the eye compares them directly; driven by the
  # data rather than the fit, and by a quantile so a single wild pair cannot
  # zoom the plot out to nothing.
  L <- max(stats::quantile(abs(unlist(st[ax])), 0.99, na.rm = TRUE, names = FALSE) * 1.15, 5)

  # Only draw a fitted drift arrow when there is exactly ONE file: drift has a
  # direction, and two acquisitions do not share one, so a single arrow over
  # pooled files would be meaningless.
  one <- if (length(dl) == 1L) dl[[1]] else NULL
  cmode  <- input$drift_plot_color %||% "time"
  # Colouring by molecule needs a molecule to name. Spatial pairing (no
  # de-blinked IDs available) leaves `mol` NA, so fall back to time.
  by_mol <- identical(cmode, "mol") && "mol" %in% names(st) && any(!is.na(st$mol))
  # Acquisition time needs t_abs; older cached results will not carry it.
  by_abs <- identical(cmode, "abs") && "t_abs" %in% names(st) && any(is.finite(st$t_abs))
  cvals  <- (if (by_abs) st$t_abs else st$dt) / 60
  cpal   <- if (by_abs) "Turbo" else "Viridis"
  ctitle <- if (by_abs) "run time (min)" else "gap (min)"

  # ── Geometric median displacement per time-gap bin ───────────────────────
  # Estimated ONCE in full dimensionality and then projected into each panel,
  # so all three panels show the same estimated point viewed three ways. A
  # per-panel 2D fit would give the XY marker a different underlying point
  # from the XZ marker, which cannot be read as one displacement.
  #
  # GEOMETRIC median, not the mean and not a per-axis median. The mean has no
  # breakdown point and these bins really are contaminated - a merged
  # "molecule" occasionally holds two different fluorophores, dropping a
  # tens-of-nm displacement into an otherwise ~10 nm cloud. A per-axis median
  # is robust but NOT rotation-equivariant, so its marker would move if x and
  # y were relabelled, which is fatal for a plot whose content is a direction.
  # The geometric median is both robust and equivariant.
  binc <- NULL
  nb <- min(5L, max(2L, floor(nrow(st) / 8)))
  br <- unique(stats::quantile(st$dt, seq(0, 1, length.out = nb + 1), na.rm = TRUE))
  if (length(br) > 2) {
    bin  <- cut(st$dt, br, include.lowest = TRUE)
    idx  <- split(seq_len(nrow(st)), bin)
    idx  <- idx[lengths(idx) > 0]
    if (length(idx) >= 2) {
      cen <- t(vapply(idx, function(ix)
        geometric_median(as.matrix(st[ix, ax, drop = FALSE])), numeric(length(ax))))
      colnames(cen) <- ax
      binc <- list(cen = cen,
                   t = vapply(idx, function(ix) stats::median(st$dt[ix], na.rm = TRUE), numeric(1)))
    }
  }

  np  <- length(pairs_ax)
  gap <- if (np > 1) 0.07 else 0
  w   <- (1 - gap * (np - 1)) / np
  fig <- plotly::plot_ly(); shapes <- list(); axis_args <- list()

  for (i in seq_len(np)) {
    a <- pairs_ax[[i]][1]; b <- pairs_ax[[i]][2]
    xid <- if (i == 1) "x" else paste0("x", i)
    yid <- if (i == 1) "y" else paste0("y", i)

    if (by_mol) {
      # One colour per molecule, cycled: with hundreds of molecules a colour
      # cannot IDENTIFY one, but it does group a molecule's blinks so the eye
      # can follow them. Each molecule is joined from the origin in time order,
      # which turns every multi-blink fluorophore into its own small track -
      # under drift those tracks all lean the same way.
      pal <- grDevices::hcl.colors(12, "Dark 3")
      mm  <- split(seq_len(nrow(st)), st$mol)
      # Cap the number drawn. Each molecule is its own plotly trace, times the
      # number of panels, so an unbounded count stalls the browser rather than
      # informing anyone. The molecules with the MOST blinks are kept, since a
      # three- or four-blink track shows a direction and a single step cannot.
      cap <- 200L
      n_mol_all <- length(mm)
      if (n_mol_all > cap) mm <- mm[order(-lengths(mm))][seq_len(cap)]
      for (mi in seq_along(mm)) {
        ix <- mm[[mi]]
        fig <- fig %>% plotly::add_trace(
          x = c(0, st[[a]][ix]), y = c(0, st[[b]][ix]), xaxis = xid, yaxis = yid,
          type = "scatter", mode = "lines+markers",
          line = list(color = pal[(mi - 1L) %% length(pal) + 1L], width = 1),
          marker = list(size = 6, opacity = 0.75,
                        color = pal[(mi - 1L) %% length(pal) + 1L]),
          name = if (n_mol_all > cap)
                   sprintf("molecule track (%s of %s shown)", cap,
                           format(n_mol_all, big.mark = ","))
                 else "molecule track",
          legendgroup = "mol",
          showlegend = (i == 1 && mi == 1), inherit = FALSE,
          text = paste("de-blinked ID", st$mol[ix][1]),
          hovertemplate = paste0("%{text}<br>", lab[[a]], " %{x:.1f} nm<br>",
                                 lab[[b]], " %{y:.1f} nm<extra></extra>"))
      }
    } else {
      # Two continuous colourings, both drawn from the app's own palette set:
      #   time - how long AFTER its molecule's first blink (Viridis)
      #   abs  - WHEN in the acquisition it happened (Turbo, blue -> red)
      # They answer different questions. The gap shows accumulation; the
      # acquisition time shows whether early and late parts of the run pull in
      # different directions, which a gap colouring cannot separate.
      fig <- fig %>% plotly::add_markers(
        x = st[[a]], y = st[[b]], xaxis = xid, yaxis = yid,
        marker = list(size = 6, opacity = 0.55, color = cvals,
                      colorscale = plotly_colorscale(cpal),
                      showscale = (i == 1),
                      colorbar = list(title = ctitle, len = 0.9, thickness = 12)),
        name = "blink displacement", legendgroup = "pts", showlegend = FALSE,
        inherit = FALSE, text = sprintf("%.1f min", cvals),
        hovertemplate = paste0(lab[[a]], " %{x:.1f} nm<br>", lab[[b]],
                               " %{y:.1f} nm<br>", ctitle,
                               " %{text}<extra></extra>"))
    }

    # The per-bin centres, projected into this panel (computed once above).
    if (!is.null(binc)) {
      mx <- binc$cen[, a]; my <- binc$cen[, b]; mt <- binc$t
      keep <- is.finite(mx) & is.finite(my)
      # add_trace, not add_paths: add_paths forces mode = "lines" and plotly
      # then warns and re-adds markers on every call. Neither sorts by x, so
      # the point order is preserved either way (add_LINES is the one to
      # avoid - it sorts, and would scramble this trajectory).
      if (sum(keep) >= 2)
        fig <- fig %>% plotly::add_trace(
          x = c(0, mx[keep]), y = c(0, my[keep]), xaxis = xid, yaxis = yid,
          type = "scatter", mode = "lines+markers",
          line = list(color = "#c03030", width = 2),
          marker = list(size = 10, color = "#c03030", symbol = "circle"),
          name = "geometric median per time bin",
          legendgroup = "bin", showlegend = (i == 1), inherit = FALSE,
          text = c("start", sprintf("%.0f min", mt[keep] / 60)),
          hovertemplate = "%{text}<br>%{x:.1f}, %{y:.1f} nm<extra></extra>")
    }

    # Fitted drift, extended over the observed range of gaps.
    if (!is.null(one) && is.finite(one$rate) && one$rate > 0) {
      v <- one$velocity; tmax <- max(st$dt, na.rm = TRUE)
      va <- v[[sub("^d", "", a)]]; vb <- v[[sub("^d", "", b)]]
      if (is.finite(va) && is.finite(vb))
        fig <- fig %>% plotly::add_paths(
          x = c(0, va * tmax), y = c(0, vb * tmax), xaxis = xid, yaxis = yid,
          line = list(color = "#2c7fb8", width = 2, dash = "dash"),
          name = "fitted drift", legendgroup = "fit", showlegend = (i == 1),
          hoverinfo = "skip", inherit = FALSE)
    }

    axl <- function(x0, x1, y0, y1) list(type = "line", xref = xid, yref = yid,
      x0 = x0, x1 = x1, y0 = y0, y1 = y1, line = list(color = "rgba(0,0,0,0.35)", width = 1))
    shapes <- c(shapes, list(axl(0, 0, -L, L), axl(-L, L, 0, 0)))

    dom <- c((i - 1) * (w + gap), (i - 1) * (w + gap) + w)
    # X yields domain to satisfy the 1:1 scale; only one axis of the pair may
    # declare scaleanchor or plotly widens a range to resolve the loop.
    axis_args[[if (i == 1) "xaxis" else paste0("xaxis", i)]] <- list(
      title = paste(lab[[a]], "from first blink (nm)"), range = c(-L, L),
      zeroline = FALSE, showgrid = FALSE, domain = dom, anchor = yid,
      scaleanchor = yid, scaleratio = 1, constrain = "domain")
    axis_args[[if (i == 1) "yaxis" else paste0("yaxis", i)]] <- list(
      title = paste(lab[[b]], "(nm)"), range = c(-L, L),
      zeroline = FALSE, showgrid = FALSE, domain = c(0, 1), anchor = xid)
  }

  do.call(plotly::layout, c(list(fig, shapes = shapes,
    showlegend = TRUE,
    legend = list(orientation = "h", x = 0.5, xanchor = "center", y = -0.22,
                  yanchor = "top"),
    margin = list(t = 10, b = 90)), axis_args)) %>%
    plotly::config(displaylogo = FALSE)
})

output$download_prec_cloud <- downloadHandler(
  filename = function() paste0("precision_cloud_perfile_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
  content = function(file) {
    r <- prec_cloud_results(); shiny::validate(shiny::need(!is.null(r), "No results."))
    write.csv(r$per_file, file, row.names = FALSE)
  })
