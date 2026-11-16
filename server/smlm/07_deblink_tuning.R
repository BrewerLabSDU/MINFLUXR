## ─────────────────────────────────────────────────────────────────────────────
## server/smlm/07_deblink_tuning.R  -  tuning aids for de-blinking and Median Clustering
## ─────────────────────────────────────────────────────────────────────────────
## Part of smlm_server(). Not a standalone script: its statements are evaluated in
## the server function's environment, in file-name order (server/parts.R).
##
## kNN distance plot for the Median Clustering eps, the centroid-merge
## threshold suggestion, per-channel trace spread, and the exploratory
## molecule-count-versus-radius curve.
##
## INTERFACE (generated: Rscript tools/server_map.R --headers; do not edit by hand)
##   DEFINES for other parts:
##     (nothing)
##   READS from other parts:
##     05_filter_pipeline: filtered_df
##     11_dcr_split: dcr_conf_tbl
##     17_plottable_data: plottable_data
## END INTERFACE
## ─────────────────────────────────────────────────────────────────────────────

# ── kNN-distance diagnostic for tuning cluster2 eps ──
# Plots, for the same centroid set cluster2 operates on (trace medians, or
# de-blinked molecule medians when de-blinking is on), the sorted distance
# to the (MinPts-1)-th neighbour. The "elbow" separates dense structures from
# scattered background = a sensible eps. Red dashed line = current eps2.
output$cluster2_knn_plot <- renderPlot({
  req(isTRUE(input$enable_cluster2))
  df <- plottable_data()
  req(!is.null(df), nrow(df) > 0)

  use_deblinked <- isTRUE(input$enable_clustering) && "cluster" %in% names(df)
  meds <- if (use_deblinked) {
    df %>%
      dplyr::filter(!is.na(.data$cluster) & .data$cluster != 0) %>%
      dplyr::group_by(.data$cluster) %>%
      dplyr::summarise(mx = stats::median(x, na.rm = TRUE),
                       my = stats::median(y, na.rm = TRUE),
                       mz = stats::median(z, na.rm = TRUE), .groups = "drop")
  } else {
    df %>%
      dplyr::group_by(.data$tid) %>%
      dplyr::summarise(mx = stats::median(x, na.rm = TRUE),
                       my = stats::median(y, na.rm = TRUE),
                       mz = stats::median(z, na.rm = TRUE), .groups = "drop")
  }

  # Must match what median clustering actually measures (always 3D), or the elbow
  # you read eps2 off is on the wrong scale - measured 1.46-1.89x too small on real
  # 3D MINFLUX data when this followed the 2D display toggle.
  cols <- c("mx","my","mz")
  pts <- as.matrix(meds[, cols, drop = FALSE])
  pts <- pts[stats::complete.cases(pts), , drop = FALSE]
  validate(need(nrow(pts) >= 3, "Not enough centroids for a kNN-distance plot."))

  k <- max(1L, as.integer(input$minpts2 %||% SMLM_DEFAULTS$minpts2) - 1L)
  k <- min(k, nrow(pts) - 1L)
  dbscan::kNNdistplot(pts, k = k)
  graphics::title(main = sprintf("kNN-distance (k = %d, n = %d centroids)", k, nrow(pts)))
  if (isTRUE(input$cluster2_method == "dbscan")) {
    e <- suppressWarnings(as.numeric(input$eps2))
    if (length(e) == 1 && is.finite(e)) graphics::abline(h = e, col = "red", lty = 2, lwd = 2)
  }
})


# ── Centroid merge: threshold suggestion + diagnostic ─────────────────────
# Same shape as the raw-point pair, but measured on CENTROID distances, which
# is what that criterion actually compares.
centroid_suggestion <- reactiveVal(NULL)
# Last measured trace spread, in whatever metric is active. Cached so that
# editing k can convert to nanometres without re-measuring on every keystroke.
centroid_sigma_rv <- reactiveVal(NA_real_)

compute_centroid_suggestion <- function() {
  df <- tryCatch(filtered_df(), error = function(e) NULL)
  if (is.null(df) || !nrow(df)) return(NULL)
  w <- if (isTRUE(input$deblink_anisotropic))
    deblink_axis_weights(df, min_n = input$sd_min_locs %||% SMLM_DEFAULTS$sd_min_locs)$w
    else NULL
  r <- suggest_centroid_threshold(
    df,
    min_n = input$sd_min_locs %||% SMLM_DEFAULTS$sd_min_locs,
    use_z = TRUE,          # decided from the data, never from the plot mode
    w     = w,
    k     = { kk <- suppressWarnings(as.numeric(input$centroid_k))
              if (!length(kk) || !is.finite(kk) || kk <= 0) 2.8 else kk })
  if (!is.null(r) && is.finite(r$sigma)) centroid_sigma_rv(r$sigma)
  r
}

# ── k drives the distance ────────────────────────────────────────────────
# The merge itself only ever reads the DISTANCE. k is a convenience for
# expressing that distance in units of the measured trace spread, so it has to
# actually move the distance - otherwise editing k changes nothing and the
# control is decorative. Debounced so typing does not re-measure per keystroke.
centroid_k_d <- debounce(reactive(input$centroid_k), 600)
observeEvent(centroid_k_d(), {
  req(isTRUE(input$enable_clustering),
      identical(input$deblink_method %||% SMLM_DEFAULTS$deblink_method, "centroid"))
  kk <- suppressWarnings(as.numeric(centroid_k_d()))
  if (!length(kk) || !is.finite(kk) || kk <= 0) return()
  sg <- centroid_sigma_rv()
  if (!is.finite(sg)) {                     # never measured yet - measure now
    r <- tryCatch(compute_centroid_suggestion(), error = function(e) NULL)
    if (is.null(r)) return()
    sg <- r$sigma
  }
  newd <- round(kk * sg, 1)
  cur  <- suppressWarnings(as.numeric(input$eps_centroid))
  if (length(cur) == 1 && is.finite(cur) && isTRUE(all.equal(cur, newd))) return()
  updateNumericInput(session, "eps_centroid", value = newd)
}, ignoreInit = TRUE)

observeEvent(input$suggest_eps_centroid, {
  r <- compute_centroid_suggestion()
  if (is.null(r)) {
    centroid_suggestion(NULL)
    showNotification("Not enough traces with enough localizations to measure the spread.",
                     type = "warning", duration = 6); return()
  }
  df <- tryCatch(filtered_df(), error = function(e) NULL)
  w  <- if (isTRUE(input$deblink_anisotropic))
    deblink_axis_weights(df, min_n = input$sd_min_locs %||% SMLM_DEFAULTS$sd_min_locs)$w
    else NULL
  # Raw per-axis spread, UNWHITENED. The whitened numbers are identical on
  # every axis by construction, so quoting them three times says nothing.
  r$sigma_ax_raw <- tryCatch(
    suggest_centroid_threshold(df, min_n = input$sd_min_locs %||% SMLM_DEFAULTS$sd_min_locs,
                               use_z = TRUE, w = NULL, k = r$k)$sigma_ax,
    error = function(e) NULL)
  r$axis_w <- w

  # DELIBERATELY NO DRIFT TERM HERE. Drift is measured, but in the precision
  # QC panel, not in this suggestion - see sample_drift_estimate(). Feeding it
  # back into the merge radius would be circular: the only molecules whose
  # drift can be observed are the ones already merged, i.e. the ones that
  # drifted LESS than the very threshold being set. The threshold is the trace
  # spread times k, and nothing else.
  updateNumericInput(session, "centroid_k", value = r$k)
  centroid_suggestion(r)
  updateNumericInput(session, "eps_centroid", value = r$suggested)
  showNotification(sprintf("Centroid threshold %.1f nm (k = %.2f x σ %.2f nm).",
                           r$suggested, r$k, r$sigma),
                   type = "message", duration = 7)
})

output$centroid_suggest_summary <- renderUI({
  r <- centroid_suggestion(); if (is.null(r)) return(NULL)
  axr <- r$sigma_ax_raw %||% r$sigma_ax
  # Derived from the LIVE input, not from the value the suggestor proposed:
  # if the user has since typed a distance by hand, this must say so.
  cur_thr <- suppressWarnings(as.numeric(input$eps_centroid))
  if (!length(cur_thr) || !is.finite(cur_thr)) cur_thr <- r$suggested
  shape <- if (isTRUE(input$deblink_anisotropic) && !is.null(r$axis_w))
    # Whitened per-axis spreads are equal by construction, so the useful
    # number is the real-space merge volume, which is an ellipsoid.
    sprintf("ellipsoid %s nm", paste(sprintf("%.1f", cur_thr / r$axis_w[seq_len(r$dims)]),
                                     collapse = " / "))
    else sprintf("sphere r = %.1f nm", cur_thr)
  tagList(
    helpText(sprintf("Trace spread %s nm  |  threshold %.1f nm = %.2f x spread (%s)  |  %s traces, %.0f%% merge.",
                     paste(sprintf("%.1f", axr), collapse = " / "),
                     cur_thr, if (is.finite(r$sigma) && r$sigma > 0) cur_thr / r$sigma else NA_real_,
                     shape, format(r$n_traces, big.mark = ","),
                     100 * (r$frac_merged %||% NA))),
    helpText("Sample drift is reported in the Localization Precision Cloud QC panel; it is not folded into this threshold.")
  )
})

centroid_gap_data <- reactive({
  req(isTRUE(input$view_centroid_gap), isTRUE(input$enable_clustering))
  compute_centroid_suggestion()
})


output$centroid_gap_plot <- renderPlot({
  r <- centroid_gap_data()
  validate(need(!is.null(r), "Not enough traces to measure the spread."))
  validate(need(sum(is.finite(r$within) & r$within > 0) >= 2,
                "Not enough within-trace distances to plot."))
  .draw_centroid_gap(r, suppressWarnings(as.numeric(input$eps_centroid)))
})

# ── PER-CHANNEL TRACE SPREAD ──────────────────────────────────────────────
# The 2-colour panel used to run this same estimator POOLED over both channels
# and set a separate eps_tc1 from it. Its own comment admitted the problem:
# "two dyes can genuinely differ in precision, and if they do this under-serves
# the tighter channel". They do differ - photon yield differs, so localization
# precision differs, so trace spread differs - and the merge distance is k x
# that spread. So it is measured per channel now, and the pooled number is a
# labelled fallback rather than the default.
#
# k stays SHARED. It says how many spreads apart two blinks can be and still be
# one molecule, which is a statement about blinking, not about the dye.
channel_spread_table <- reactive({
  df <- tryCatch(filtered_df(), error = function(e) NULL)
  if (is.null(df) || !nrow(df)) return(NULL)
  df <- tryCatch(attach_channel_colors(resolve_channels(df, as_smlm_settings(input)), smlm_channel_colors(input)),
                 error = function(e) NULL)
  if (is.null(df) || !nrow(df) || !"channel" %in% names(df)) return(NULL)
  thr  <- tryCatch(deblink_threshold_for_channel(df, as_smlm_settings(input)), error = function(e) NULL)
  if (is.null(thr)) return(NULL)
  ct   <- attr(df, "channel_table")
  meas <- attr(thr, "measured"); sig <- attr(thr, "sigma")
  ct$sigma    <- as.numeric(sig[as.character(ct$channel)])
  ct$distance <- as.numeric(unlist(thr[as.character(ct$channel)]))
  ct$measured <- as.logical(meas[as.character(ct$channel)])
  ct
})

# ── WHAT THE CHANNEL RULE COSTS ───────────────────────────────────────────
# "Never merge across a channel" is the right rule, and it is not free. Both
# costs below bias the de-blinked molecule count DOWNWARD, and neither was
# visible anywhere before this readout.
#
# 1. AMBIGUOUS DCR. A trace is assigned to a channel by its MEDIAN DCR, so two
#    blinks of one molecule that land on opposite sides of the split can never
#    be merged - the molecule is counted twice, once per channel. The size of
#    the effect is already computed as the GMM posterior behind the split
#    confidence plot; it is just never stated where the merge is configured.
#
#    Note this is the SAME root cause as the crosstalk warning on the split
#    itself, but a different consequence. Crosstalk pushes cross-channel
#    statistics UP at short r (false co-localization); here it pushes molecule
#    counts DOWN. One mis-assignment, two directions.
#
# 2. STRADDLING TRACES (multiplex only). A multiplex channel is assigned per
#    LOCALIZATION, so a trace blinking across a TIM breakpoint sits in two
#    channels and is split into two molecules. Merging it would invent a
#    molecule that was never simultaneously present, so the split is correct -
#    but it is a real count, not a rounding error, and it belongs on screen.
output$channel_integrity_summary <- renderUI({
  items <- list()

  # 1. ambiguous DCR band
  if (isTRUE(input$enable_two_color)) {
    tb <- tryCatch(dcr_conf_tbl(), error = function(e) NULL)
    if (is.data.frame(tb) && nrow(tb)) {
      n_edge <- sum(tb$on_edge, na.rm = TRUE)
      if (n_edge > 0)
        items <- c(items, list(tags$div(
          sprintf("%s of %s traces sit in the ambiguous DCR band (%.1f%%). A molecule whose blinks straddle the split is counted once in each channel, so the molecule count is biased low by about that fraction.",
                  format(n_edge, big.mark = ","), format(nrow(tb), big.mark = ","),
                  100 * n_edge / nrow(tb)))))
    }
  }

  # 2. traces straddling a multiplex boundary
  df <- tryCatch(filtered_df(), error = function(e) NULL)
  if (is.data.frame(df) && nrow(df) && "mplex_channel" %in% names(df) &&
      any(!is.na(df$mplex_channel))) {
    per_tid <- tapply(df$mplex_channel, df$tid, function(v) length(unique(v[!is.na(v)])))
    n_str <- sum(per_tid > 1L, na.rm = TRUE)
    if (n_str > 0)
      items <- c(items, list(tags$div(style = "margin-top:4px;",
        sprintf("%s of %s traces span more than one multiplex channel and are split into one molecule per channel. Merging them would claim a molecule was present in two exchange rounds at once.",
                format(n_str, big.mark = ","), format(length(per_tid), big.mark = ",")))))
  }

  if (!length(items)) return(NULL)
  tags$div(style = "margin-top:8px; padding:6px 9px; border-left:3px solid #e0c98a; background:#fdfaf2; font-size:0.82em; color:#5a4a20;",
           tags$b("What the channel rule costs here"), tags$br(), items)
})

output$channel_spread_summary <- renderUI({
  ct <- tryCatch(channel_spread_table(), error = function(e) NULL)
  if (is.null(ct) || nrow(ct) < 2L) return(NULL)

  rows <- lapply(seq_len(nrow(ct)), function(i) {
    col <- ct$color[i]; if (is.na(col)) col <- "#888888"
    tags$div(style = "display:flex; align-items:center; gap:8px; font-size:0.86em; margin:2px 0;",
      tags$span(style = sprintf("width:10px;height:10px;border-radius:2px;background:%s;display:inline-block;", col)),
      tags$b(as.character(ct$label[i])),
      tags$span(if (isTRUE(ct$measured[i]))
                  sprintf("spread %.1f nm  ->  merge at %.1f nm", ct$sigma[i], ct$distance[i])
                else sprintf("spread not measurable  ->  shared %.1f nm", ct$distance[i])),
      tags$span(style = "color:#7b8794;", sprintf("(%s locs)", format(ct$n[i], big.mark = ","))))
  })

  # The spread RATIO is the number that decides whether one shared distance is
  # defensible, so state it rather than making the user divide two numbers.
  sg <- ct$sigma[is.finite(ct$sigma)]
  note <- if (length(sg) >= 2) {
    rr <- max(sg) / min(sg)
    if (isTRUE(input$deblink_shared_spread))
      helpText(style = "color:#7a5c00;", sprintf(
        "Using ONE shared distance while the channels differ in spread by %.2fx. The broader channel is under-merged and the tighter one over-merged.", rr))
    else if (rr >= 1.15)
      helpText(sprintf("Channels differ in spread by %.2fx, so they are merged at different distances.", rr))
    else
      helpText(sprintf("Channels agree in spread to within %.2fx - a shared distance would give nearly the same result here.", rr))
  } else NULL

  tagList(tags$div(style = "margin:4px 0 2px; font-weight:600; font-size:0.86em;",
                   "Merge distance per channel"), rows, note)
})

# ── EXPLORATORY: molecule count vs radius ─────────────────────────────────
# Gated on its own checkbox: the sweep runs one DBSCAN per radius, so it must
# not fire on every filter change.
deblink_knee_data <- reactive({
  req(isTRUE(input$view_deblink_knee), isTRUE(input$enable_clustering))
  df <- tryCatch(filtered_df(), error = function(e) NULL)
  req(!is.null(df), nrow(df) > 0)
  # Same metric the clustering uses, or the radii on the x-axis mean nothing.
  w <- if (isTRUE(input$deblink_anisotropic))
    deblink_axis_weights(df, min_n = input$sd_min_locs %||% SMLM_DEFAULTS$sd_min_locs)$w
    else NULL
  # Sweep whatever criterion is ACTIVE - a curve over raw-point eps says
  # nothing about a centroid threshold, and vice versa.
  # Sweep range follows the DATA, not a hard-coded 40 nm: on the centroid
  # criterion the interesting structure runs well past that.
  # Centroid criterion only - see deblink_count_curve().
  req(identical(input$deblink_method %||% SMLM_DEFAULTS$deblink_method, "centroid"))
  cs  <- tryCatch(compute_centroid_suggestion(), error = function(e) NULL)
  top <- if (!is.null(cs))
    max(40, ceiling(stats::quantile(cs$between, 0.75, names = FALSE)), 4 * cs$suggested)
    else 40
  deblink_count_curve(df, minpts = deblink_minpts(input$tid_length_range),
                      eps_seq = seq(1, min(top, 400), length.out = 60),
                      w = w, method = "centroid")
})

output$deblink_knee_plot <- renderPlot({
  r <- deblink_knee_data()
  validate(need(!is.null(r), "Not enough data to sweep the radius."))
  cv <- r$curve
  graphics::plot(cv$eps, cv$n_molecules, type = "l", lwd = 2.5, col = "#1f78b4",
                 xlab = if (identical(input$deblink_method %||% SMLM_DEFAULTS$deblink_method, "vote"))
                          "vote fraction f" else "centroid distance threshold (nm)",
                 ylab = "de-blinked molecules",
                 main = "Molecule count vs merge threshold")
  graphics::points(cv$eps, cv$n_molecules, pch = 16, cex = 0.5, col = "#1f78b4")
  graphics::abline(h = r$n_traces, col = "#969696", lty = 3, lwd = 2)   # no merging at all
  if (is.finite(r$knee)) graphics::abline(v = r$knee, col = "#e08214", lwd = 2, lty = 4)
  e_cur <- suppressWarnings(as.numeric(
    if (identical(input$deblink_method %||% SMLM_DEFAULTS$deblink_method, "vote")) input$vote_frac
    else input$eps_centroid))
  if (length(e_cur) == 1 && is.finite(e_cur))
    graphics::abline(v = e_cur, col = "#d7191c", lwd = 2.5, lty = 2)
  graphics::legend("topright", bty = "n", cex = 0.9,
    legend = c("molecules", "traces (no merging)", "knee", "current threshold"),
    col = c("#1f78b4", "#969696", "#e08214", "#d7191c"),
    lty = c(1, 3, 4, 2), lwd = c(2.5, 2, 2, 2.5))
})

output$deblink_knee_summary <- renderUI({
  r <- tryCatch(deblink_knee_data(), error = function(e) NULL)
  if (is.null(r)) return(NULL)
  cv <- r$curve
  at <- function(e) { i <- which.min(abs(cv$eps - e)); cv$n_molecules[i] }
  e_cur <- suppressWarnings(as.numeric(
    if (identical(input$deblink_method %||% SMLM_DEFAULTS$deblink_method, "vote")) input$vote_frac
    else input$eps_centroid))
  tagList(
    helpText(sprintf("%s traces → %s molecules at the knee (%.0f nm)%s.",
                     format(r$n_traces, big.mark = ","),
                     if (is.finite(r$knee)) format(at(r$knee), big.mark = ",") else "-",
                     r$knee,
                     if (length(e_cur) == 1 && is.finite(e_cur))
                       sprintf("; %s at your current %.0f nm",
                               format(at(e_cur), big.mark = ","), e_cur) else "")),
    if (r$n_locs < 60000L) NULL else
      helpText("Curve computed on a 60,000-localization subsample."),
    # Measured: the knee detector is calibrated for the raw-point curve, whose
    # steep initial drop comes from traces FRAGMENTING at small eps. The
    # centroid curve has no such regime - it simply starts at "no merging" and
    # declines - so the curvature maximum lands mid-slope and reads far too
    # low (synthetic: knee 5 nm -> 422 molecules, against 360 truth and 360 at
    # the σ-based suggestion). Do not tune the centroid threshold from it.
    if (identical(input$deblink_method %||% SMLM_DEFAULTS$deblink_method, "centroid"))
      helpText(tags$b("Unreliable for this criterion: "),
               "the knee is calibrated for the raw-point curve shape and reads too low ",
               "on centroid distances. Use 'Suggest threshold from trace spread' instead; ",
               "the curve itself is still worth reading.")
    else NULL
  )
})
