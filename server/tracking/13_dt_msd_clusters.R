# --- dt extraction per TID (based on raw 'tim') --------------------------
dt_vectors <- reactive({
  df <- trk(); req(nrow(df) > 0)
  df <- df[is.finite(df$tim), ]
  df <- df[order(df$tid, df$tim), ]
  sp  <- split(df$tim, df$tid)
  lapply(sp, function(tt) {
    if (length(tt) < 2) return(numeric(0))
    dt <- diff(tt)
    dt[is.finite(dt) & dt >= 0]
  })
})
# --- dt histogram ---------------------------------------------------------

# --- Δt histogram in analysis tab (richer version with stats) ---------------
output$trk_dt_hist_jdd <- renderPlotly({
  req(isTRUE(input$trk_show_tsa), tsa_ran())
  dts  <- dt_vectors()
  vec  <- unlist(dts, use.names = FALSE)
  vec  <- vec[is.finite(vec) & vec > 0]
  shiny::validate(shiny::need(length(vec) > 0,
    "No Δt values - need at least 2 time points per TID."))
  n_bins  <- 80L    # only used by the no-cycle-time fallback binning below
  dark    <- isTRUE(isolate(input$trk_dark_mode))
  bg_col  <- if (dark) "#000000" else "#FFFFFF"
  fg_col  <- if (dark) "#FFFFFF" else "#333333"

  dt_med  <- median(vec, na.rm = TRUE)
  # Bin-FREE typical Δt (half-sample mode) - independent of the bin slider.
  dt_hsm  <- half_sample_mode(vec)
  rng     <- range(vec)                                 # x-extent, for the cycle markers

  # The app's best cycle-time estimate anchors the bins and the cycle markers. Read
  # reactively so both follow a cycle change.
  t_c <- tryCatch(trk_cycle_time(), error = function(e) NA_real_)
  # SW offset: peaks sit at SW + k·t_c, so the bin grid has to be shifted by it.
  # Prefer the fitted value; fall back to (typical Δt − t_c), which is the same
  # quantity read straight off the 1-cycle peak.
  cc_sw <- {
    cc <- tryCatch(trk_cluster_cycle(), error = function(e) NULL)
    s  <- if (!is.null(cc) && is.finite(cc$sw)) cc$sw else NA_real_
    if (!is.finite(s) && is.finite(dt_hsm) && is.finite(t_c)) max(0, dt_hsm - t_c) else s
  }
  if (!is.finite(cc_sw)) cc_sw <- 0

  # ── Binning: one bar per cycle count ──────────────────────────────────────
  # Δt is not continuous: it takes values SW + k·t_c for integer k. Fixed-width bins
  # alias against that comb - at 80 bins the grid runs ~1.63 bins per cycle, so each
  # peak lands at a different offset inside its bin and the bars drift visibly away
  # from the cycle markers, with empty bins interleaved. Instead the bin EDGES sit on
  # the classification boundaries SW + (k ± ½)·t_c, so every bar is exactly one cycle
  # count centred on its peak, and the picture is identical to the
  # round((Δt − SW)/t_c) the analysis performs. Falls back to fixed-width bins only
  # when no cycle time is available.
  use_cycle_bins <- is.finite(t_c) && t_c > 0
  if (use_cycle_bins) {
    kv    <- round((vec - cc_sw) / t_c)
    kv    <- kv[is.finite(kv) & kv >= 1]
    kmax  <- max(1L, min(as.integer(max(kv)), 400L))
    kcnt  <- tabulate(pmin(kv, kmax), nbins = kmax)
    bar_x <- cc_sw + seq_len(kmax) * t_c      # bar centre = the peak itself
    bar_y <- kcnt
    bar_w <- t_c * 0.9
    keepb <- bar_y > 0
    bar_x <- bar_x[keepb]; bar_y <- bar_y[keepb]
  } else {
    h_tmp <- hist(vec, breaks = n_bins, plot = FALSE)
    bar_x <- h_tmp$mids; bar_y <- h_tmp$counts
    bar_w <- diff(h_tmp$breaks)[1]
  }

  # Counts span ~100× between the 1-cycle peak and the dark tail, so the tail is
  # invisible on a linear axis; log y makes the whole decay readable.
  log_y <- isTRUE(input$trk_dt_log_y %||% FALSE)
  ymax  <- if (log_y) 10^(log10(max(bar_y, na.rm = TRUE)) * 1.06)
           else max(bar_y, na.rm = TRUE) * 1.08
  ybase <- if (log_y) 0.7 else 0        # log axes cannot start at 0

  fig <- plot_ly() %>%
    add_bars(x = bar_x, y = bar_y, width = bar_w,
             marker = list(color = "rgba(100,150,220,0.7)",
                           line  = list(color = "rgba(70,110,180,0.9)", width = 0.5)),
             hovertemplate = if (use_cycle_bins)
               "Δt %{x:.4g} ms<br>%{y:,} steps<extra></extra>"
             else "Δt %{x:.4g} ms<br>%{y:,}<extra></extra>",
             showlegend = FALSE)
  fig <- fig %>%
    add_trace(x = c(dt_med, dt_med), y = c(ybase, ymax),
              type = "scatter", mode = "lines",
              line   = list(color = "#E03020", width = 2, dash = "dash"),
              name = sprintf("Median %.3g ms", dt_med)) %>%
    add_trace(x = c(dt_hsm, dt_hsm), y = c(ybase, ymax),
              type = "scatter", mode = "lines",
              line   = list(color = "#20A040", width = 2, dash = "dashdot"),
              name = sprintf("Typical Δt (bin-free) %.3g ms", dt_hsm))

  # ── Cycle-multiple markers ────────────────────────────────────────────────
  # A k-cycle interval sits at Δt = SW + k·t_c, i.e. modal + (k−1)·t_c (the peaks
  # are NOT at k·t_c - that omits the once-per-step SW overhead and lands every
  # marker short of its bell). So anchor marker 1 on the observed modal Δt (the
  # tallest = 1-cycle peak, = the green Typical-Δt line) and step by the SELECTED
  # cycle time t_c. Consequence: if the selected cycle time matches the data, all
  # markers land on the bells; if a manual/hardware override is wrong, the higher
  # markers visibly drift off - a built-in check that the chosen cycle time is right.
  # Layout shapes (not traces) so they stay out of the legend.
  shp <- list(); ann <- list()
  if (isTRUE(input$trk_dt_mark_cycles %||% TRUE) &&
      is.finite(t_c) && t_c > 0 && is.finite(dt_hsm)) {
    n_cap <- suppressWarnings(as.integer(input$trk_cycle_dark_max %||% 10L))
    if (!is.finite(n_cap) || n_cap < 1L) n_cap <- 10L
    Nmax  <- max(1L, min(n_cap, 12L, as.integer(floor((rng[2] - dt_hsm) / t_c)) + 1L))
    grid_c <- if (dark) "rgba(200,200,200,0.35)" else "rgba(90,90,90,0.35)"
    for (nn in seq_len(Nmax)) {
      xc <- dt_hsm + (nn - 1L) * t_c        # SW + nn·t_c, anchored on the modal peak
      shp[[length(shp) + 1L]] <- list(type = "line", x0 = xc, x1 = xc,
        y0 = ybase, y1 = ymax,
        line = list(color = grid_c, width = 1, dash = "dot"), layer = "below")
      ann[[length(ann) + 1L]] <- list(x = xc, y = ymax, text = paste0(nn, "×"),
        showarrow = FALSE, yanchor = "bottom",
        font = list(size = 9, color = fg_col))
    }
  }

  fig %>% layout(
    xaxis  = list(title = "Δt between consecutive localizations (ms)",
                  color = fg_col),
    yaxis  = if (log_y)
               list(title = "Count (log)", color = fg_col, type = "log",
                    range = c(log10(ybase), log10(ymax)))
             else list(title = "Count", color = fg_col, range = c(0, ymax)),
    paper_bgcolor = bg_col, plot_bgcolor = bg_col,
    font   = list(color = fg_col),
    # Vertical legend to the RIGHT of the plot (x > 1 = outside the axes) so it no
    # longer overlaps the bars / cycle labels along the top.
    legend = list(orientation = "v", x = 1.02, xanchor = "left",
                  y = 1, yanchor = "top", font = list(color = fg_col)),
    margin = list(r = 10),
    bargap = 0.05,
    shapes = shp,
    annotations = ann,
    uirevision = "dt_jdd"
  )
})

output$trk_dt_stats_ui <- renderUI({
  req(isTRUE(input$trk_show_tsa), tsa_ran())
  dts <- dt_vectors()
  vec <- unlist(dts, use.names = FALSE)
  vec <- vec[is.finite(vec) & vec > 0]
  req(length(vec) > 0)
  dt_hsm  <- half_sample_mode(vec)      # bin-free typical Δt
  dt_med  <- median(vec)
  dt_p10  <- quantile(vec, 0.10)
  dt_p25  <- quantile(vec, 0.25)
  dt_p75  <- quantile(vec, 0.75)
  hw      <- trk_hw_cycle_time()
  # Peak-based cycle time: the Δt peaks sit at k·cycle + SW, so a least-squares fit
  # of the peak centres against the integer cycle index k gives slope = cycle time,
  # intercept = SW. This is the "fit the 1×,2×,3×… peaks and read the cycle time off
  # the spacing" idea, and is the same quantity Bela/Vogler recover from the EOD/SW
  # time-signature fit. Computed here regardless of the acquisition-tab toggle so it
  # is always visible next to the density peaks it comes from.
  rc <- tryCatch(trk_regression_cycle(), error = function(e) NULL)
  cc <- tryCatch(trk_cluster_cycle(),    error = function(e) NULL)
  tags$div(style = "font-size:0.84em; color:#555; margin:4px 0;",
    tags$div(
      tags$b("Typical Δt (bin-free, half-sample mode): "),
      tags$span(style = "color:#20803a;", sprintf("%.4g ms", dt_hsm)),
      sprintf("   ·   median %.4g   ·   p10 %.4g   ·   p25 %.4g   ·   p75 %.4g   ·   n = %s",
              dt_med, dt_p10, dt_p25, dt_p75, format(length(vec), big.mark = ","))
    ),
    if (!is.null(rc) && is.finite(rc$cycle) && rc$cycle > 0)
      tags$div(style = "color:#3a5a8a; margin-top:2px;",
        tags$b("Peak-fit cycle time: "),
        tags$span(style = "font-weight:bold;", sprintf("%.4g ms", rc$cycle)),
        sprintf("   ·   SW offset %.3g ms   ·   R² = %.4f   ·   fit on %d cycle peaks",
                rc$sw, rc$r2 %||% NA_real_, rc$n),
        if (!is.null(cc) && is.finite(cc$frac_one_cycle))
          sprintf("   ·   %.0f%% are single (integrated) cycles", 100 * cc$frac_one_cycle) else "",
        # The cycle time is smaller than every gap you can actually observe, which
        # reads as a discrepancy unless the decomposition is spelled out.
        if (is.finite(rc$cycle) && is.finite(rc$sw))
          tags$div(style = "color:#666; font-weight:normal; margin-top:1px;",
            sprintf("Cycle time %.4g + SW %.3g = %.4g ms, the typical \u0394t between consecutive localizations. The cycle time itself is the repeat unit, not an observable gap.",
                    rc$cycle, rc$sw, rc$cycle + rc$sw))
        else NULL)
    else if (!is.null(cc) && is.finite(cc$cycle) && cc$cycle > 0)
      tags$div(style = "color:#3a5a8a; margin-top:2px;",
        tags$b("Peak-spacing cycle time: "),
        tags$span(style = "font-weight:bold;", sprintf("%.4g ms", cc$cycle)),
        sprintf("   ·   from %d Δt peaks (need ≥2 for the regression fit)", length(cc$centers %||% numeric(0)))),
    if (is.finite(hw) && hw > 0)
      tags$div(style = "color:#888; margin-top:2px;",
        sprintf("Hardware cycle = %.4g ms. A typical \u0394t near twice this means dark cycles dominate.", hw)),
    tags$div(style = "color:#888; margin-top:2px; font-size:0.95em;",
      "This estimate is independent of the histogram bin slider (which is display-only).",
      " The peak-fit cycle time regresses the Δt density peaks (1×, 2×, 3×… the cycle) against",
      " their cycle index - enable ", tags$em("Fit cycle time from data (EOD/SW regression)"),
      " in Acquisition Parameters to use it everywhere.")
  )
})


#---- MSD export ----
msd_export <- reactive({
  ub <- msd_unified_bundle(); req(ub)
  linear_tau <- isTRUE(input$trk_uni_linear_tau)
  tau_all <- if (linear_tau) ub$tau_lin else ub$tau
  ok  <- is.finite(tau_all) & is.finite(ub$msd) & ub$keep_lag
  fr  <- tryCatch(msd_unified_fit(ub), error = function(e) NULL)
  td  <- tryCatch(trk_unified_tdiff(), error = function(e) NULL)

  data <- data.frame(
    cycle_lag = ub$lag[ok],
    tau_ms    = tau_all[ok],
    MSD_um2   = ub$msd[ok] / 1e6,
    n_pairs   = if (!is.null(ub$npairs)) ub$npairs[ok] else NA_integer_,
    n_tracks  = ub$n_trk_lag[ok],
    stringsAsFactors = FALSE
  )

  params <- list(
    mode               = "MSD (cumulative cycle-index binning)",
    time_unit          = "ms",
    space_unit         = "µm",
    tau_axis           = if (linear_tau) "linearized lag × cycle" else "real mean Δt per bin",
    ensemble_weighting = if (isTRUE(ub$track_weighted)) "track-weighted (eMSD)" else "pair-weighted",
    fit_method         = if (isTRUE(input$trk_uni_iterative)) "iterative OLSF" else "anchored OLSF + Michalet",
    fit_window_cycles  = sprintf("%s-%s", input$trk_uni_min_lag %||% 5, input$trk_uni_max_lag %||% 50),
    cycle_time_ms      = ub$t_cycle,
    msd_ndim           = ub$ndim,
    eDIFF_um2_s        = if (!is.null(fr)) fr$D / 1e3 else NA_real_,
    tDIFF_um2_s        = if (!is.null(td)) td$tDIFF else NA_real_,
    n_tracks           = ub$n_tracks,
    msd_formula        = sprintf("MSD = %.4gD*tau + %.4g*sigma^2 (d_eff = %.4g)",
                                 ub$coef_eff %||% (2L * ub$ndim),
                                 ub$coef_eff %||% (2L * ub$ndim),
                                 ub$d_eff %||% ub$ndim)
  )

  list(data = data, params = params)
})
# Per-track D and alpha fits - computed when histogram OR clustering is open
# Per-track fits feeding the D histogram, the D/CI clustering and the per-track
# CSV export - all derived from the SAME Unified MSD calculation (per-track OLSF
# D + confinement index), via the cached trk_unified_tdiff reactive.
msd_pertk_fits <- reactive({
  req(isTRUE(input$trk_show_msd))
  req(isTRUE(input$trk_msd_show_hist) || isTRUE(input$trk_msd_cluster))
  td <- tryCatch(trk_unified_tdiff(), error = function(e) NULL)
  shiny::validate(shiny::need(!is.null(td) && !is.null(td$pertrack) && nrow(td$pertrack) > 0,
    "Enable the MSD and click 'Plot MSD' to compute per-track diffusion."))
  td$pertrack
})

# --- MSD cluster analysis ---------------------------------------------------
# k-means on (log10 D, log10 CI); clusters sorted ascending by mean log10(D) so that
# cluster 1 = most confined / slowest, cluster k = most mobile / fastest.
MSD_CLUSTER_PAL <- c("#3060C0","#C03030","#2CA02C","#FF7F0E","#9467BD","#8C564B")

msd_cluster_rv <- reactive({
  req(isTRUE(input$trk_show_msd), isTRUE(input$trk_msd_cluster))
  fits <- msd_pertk_fits()
  req(!is.null(fits) && nrow(fits) >= 3L)

  k <- max(2L, min(6L, as.integer(input$trk_msd_n_clusters %||% 3L)))
  req(nrow(fits) >= k)

  # Feature matrix: alpha (primary) + log10(CI) (secondary).
  # After z-scoring, alpha is doubled so it is the dominant axis.
  # With k=2 the split is primarily along alpha (confined vs free/directed);
  # CI adds refinement for k >= 3.
  req(!is.null(fits$conf_index))
  ci_safe   <- pmax(fits$conf_index, 1e-4)
  logD_safe <- log10(pmax(fits$D_um2_s, 1e-9))
  feat <- scale(cbind(logD_safe, log10(ci_safe)))
  feat[!is.finite(feat)] <- 0
  feat[, 1] <- feat[, 1] * 2  # log10(D) is the primary separator

  set.seed(42L)
  km <- tryCatch(
    kmeans(feat, centers = k, nstart = 25L, iter.max = 200L),
    error = function(e) NULL
  )
  req(!is.null(km))

  # Sort clusters by mean log10(D) ascending: cluster 1 = most confined, cluster k = most mobile
  mean_logD <- tapply(logD_safe, km$cluster, mean, na.rm = TRUE)
  rank_ord  <- order(mean_logD)
  new_lbl  <- integer(nrow(fits))
  for (r in seq_len(k)) new_lbl[km$cluster == rank_ord[r]] <- r

  pal     <- MSD_CLUSTER_PAL[seq_len(k)]
  col_map <- setNames(pal[new_lbl], as.character(fits$tid))

  fits$cluster    <- new_lbl
  lbl_names <- c("confined","intermediate","mobile","fast","faster","fastest")
  fits$cluster_lbl <- paste0("Cluster ", new_lbl, " (", lbl_names[new_lbl], ")")

  list(col_map = col_map, palette = pal, n_clusters = k,
       fits = fits, cluster_ids = setNames(new_lbl, as.character(fits$tid)))
})

# ── Per-cluster motion model ────────────────────────────────────────────────
# The k-means split above groups tracks by D and confinement index, but it does not
# say WHAT KIND of motion each group shows. Here each cluster's own ensemble MSD
# (mean of its members' TAMSDs) is fitted with all four models and the best one is
# selected by AICc, so a cluster can be reported as free / confined / directed /
# anomalous rather than just "slow" or "fast".
#
# AICc (not raw residuals) because the models differ in parameter count: linear has
# 2, the others 3. n here is the number of MSD lag points in the fit window, which is
# small, so the small-sample correction matters. Lags are correlated, so this is a
# ranking heuristic, not a likelihood test - reported as such.
msd_cluster_models <- reactive({
  cl <- tryCatch(msd_cluster_rv(), error = function(e) NULL)
  ub <- tryCatch(msd_unified_bundle(), error = function(e) NULL)
  if (is.null(cl) || is.null(ub) || is.null(ub$track_msd)) return(NULL)
  min_lag <- max(1L, as.integer(input$trk_uni_min_lag %||% 4L))
  max_lag <- max(min_lag + 2L, as.integer(input$trk_uni_max_lag %||% 50L))
  tau_all <- if (isTRUE(input$trk_uni_linear_tau)) ub$tau_lin else ub$tau
  win <- ub$lag >= min_lag & ub$lag <= max_lag & ub$keep_lag & is.finite(tau_all)
  if (sum(win) < 5L) return(NULL)

  ids <- cl$cluster_ids
  rows <- lapply(seq_len(cl$n_clusters), function(cc) {
    tids_c <- names(ids)[ids == cc]
    cols   <- which(colnames(ub$track_msd) %in% tids_c)
    if (!length(cols)) return(NULL)
    # Cluster ensemble MSD = mean of member TAMSDs (same convention as the eMSD).
    y <- rowMeans(ub$track_msd[, cols, drop = FALSE], na.rm = TRUE)[win]
    s <- tau_all[win]
    okp <- is.finite(s) & is.finite(y) & y > 0
    if (sum(okp) < 5L) return(NULL)
    s <- s[okp]; y <- y[okp]; np <- ub$npairs[win][okp]

    aicc_of <- function(fr, k_par) {
      if (is.null(fr) || is.null(fr$tau_fit)) return(NA_real_)
      yhat <- stats::approx(fr$tau_fit, fr$msd_fit, xout = s, rule = 2)$y
      rss  <- sum((y - yhat)^2); n <- length(y)
      if (!is.finite(rss) || rss <= 0 || n <= k_par + 2L) return(NA_real_)
      n * log(rss / n) + 2 * k_par + (2 * k_par * (k_par + 1)) / (n - k_par - 1)
    }
    # The free model is fitted here by plain OLS over the SAME window as the other
    # three, NOT by the OLSF used for the headline D. OLSF deliberately fits only the
    # first p (Michalet-optimal) lags, so its residuals over the full window are not
    # comparable with full-window non-linear fits - scoring it that way made every
    # free cluster look like a 3-parameter model (free won 0/10 on simulated free
    # curves; 9/10 once fitted on the same window). The headline D is unaffected;
    # this OLS is used only to rank the models against each other.
    cv <- 2L * ub$ndim
    free_fr <- tryCatch({
      lf <- stats::lm(y ~ s)
      list(D = unname(stats::coef(lf)[2]) / cv,
           sigma_loc = sqrt(max(0, unname(stats::coef(lf)[1])) / cv),
           tau_fit = s, msd_fit = unname(stats::fitted(lf)))
    }, error = function(e) NULL)
    cand <- list(
      free      = list(fr = free_fr, k = 2L),
      confined  = list(fr = tryCatch(.msd_model_fit(s, y, np, ub$ndim, "confined",
                                      FALSE, ub$t_cycle, coef_override = ub$coef_eff, R = trk_blur()$R_msd),
                                     error = function(e) NULL), k = 3L),
      directed  = list(fr = tryCatch(.msd_model_fit(s, y, np, ub$ndim, "directed",
                                      FALSE, ub$t_cycle, coef_override = ub$coef_eff, R = trk_blur()$R_msd),
                                     error = function(e) NULL), k = 3L),
      anomalous = list(fr = tryCatch(.msd_model_fit(s, y, np, ub$ndim, "anomalous",
                                      FALSE, ub$t_cycle, coef_override = ub$coef_eff, R = trk_blur()$R_msd),
                                     error = function(e) NULL), k = 3L))
    aic <- vapply(names(cand), function(nm) aicc_of(cand[[nm]]$fr, cand[[nm]]$k), numeric(1))
    if (!any(is.finite(aic))) return(NULL)
    auto_best <- names(cand)[which.min(replace(aic, !is.finite(aic), Inf))]
    # Per-cluster override: the whole point is to be able to say "Cluster 1 is free,
    # Cluster 2 is confined" and see both fits, so an explicit choice wins over AICc.
    # "auto" (default) keeps the AICc pick.
    pick <- input[[paste0("trk_cmodel_", cc)]] %||% "auto"
    if (!identical(pick, "auto") && !is.null(cand[[pick]]$fr)) best <- pick else best <- auto_best
    bfr  <- cand[[best]]$fr
    forced <- !identical(pick, "auto") && !identical(best, auto_best)
    # A confined fit that never reached its plateau is really just a line.
    if (identical(best, "confined") && isTRUE(bfr$degenerate)) best <- "free (no plateau)"
    detail <- switch(sub(" .*", "", best),
      confined  = sprintf("L_c = %.0f nm", bfr$L_conf),
      directed  = sprintf("v = %.3g µm/s", bfr$v_um_s),
      anomalous = sprintf("α = %.2f", bfr$alpha),
      "")
    d_aic <- sort(aic[is.finite(aic)])
    list(
      row = data.frame(
        cluster = cc, n_tracks = length(tids_c),
        model = best, auto_model = auto_best, forced = forced,
        D_um2_s = if (is.finite(bfr$D)) bfr$D / 1e3 else NA_real_,
        detail = detail,
        # Gap to the runner-up: small = the models are not really distinguishable here.
        dAICc = if (length(d_aic) >= 2L) unname(d_aic[2] - d_aic[1]) else NA_real_,
        stringsAsFactors = FALSE),
      # Observed cluster eMSD + the fitted curve, for the overlay plot (µm², ms).
      curve = list(cluster = cc, tau = s, msd = y,
                   fit_tau = bfr$tau_fit, fit_msd = bfr$msd_fit))
  })
  rows <- Filter(Negate(is.null), rows)
  if (!length(rows)) return(NULL)
  list(table  = do.call(rbind, lapply(rows, `[[`, "row")),
       curves = lapply(rows, `[[`, "curve"))
})

# One model selector per cluster. Auto = the AICc pick; anything else forces that
# model on that cluster, so different populations can be modelled differently.
output$trk_msd_cluster_model_ui <- renderUI({
  req(isTRUE(input$trk_show_msd), isTRUE(input$trk_msd_cluster))
  cl <- tryCatch(msd_cluster_rv(), error = function(e) NULL)
  if (is.null(cl)) return(NULL)
  ch <- c("Auto (best AICc)" = "auto", "Free / Brownian" = "free",
          "Confined" = "confined", "Directed" = "directed", "Anomalous" = "anomalous")
  tags$div(style = "margin-top:4px;",
    tags$div(style = "font-size:0.85em; color:#555; margin-bottom:2px;",
             "Model per cluster"),
    lapply(seq_len(cl$n_clusters), function(cc) {
      div(style = "display:inline-block; width:49%; padding-right:4px; vertical-align:top;",
          selectInput(paste0("trk_cmodel_", cc),
                      label = HTML(sprintf(
                        "<span style='color:%s;font-weight:bold;'>Cluster %d</span>",
                        cl$palette[cc], cc)),
                      choices = ch, selected = "auto", width = "100%"))
    }))
})

output$trk_msd_cluster_models <- renderUI({
  req(isTRUE(input$trk_show_msd), isTRUE(input$trk_msd_cluster))
  res <- tryCatch(msd_cluster_models(), error = function(e) NULL)
  if (is.null(res)) return(NULL)
  tb  <- res$table
  pal <- tryCatch(msd_cluster_rv()$palette, error = function(e) NULL)
  tags$div(style = "font-size:0.85em; margin-top:6px;",
    tags$b("Motion model per cluster"),
    tags$table(style = "width:100%; border-collapse:collapse; margin-top:3px;",
      tags$tr(style = "color:#666; text-align:left;",
        tags$th("Cluster"), tags$th("Tracks"), tags$th("Model"),
        tags$th("D (µm²/s)"), tags$th(""), tags$th("ΔAICc")),
      lapply(seq_len(nrow(tb)), function(i) {
        r <- tb[i, ]
        tags$tr(style = "border-top:1px solid #eee;",
          tags$td(tags$span(style = sprintf("color:%s; font-weight:bold;",
                                            if (!is.null(pal)) pal[r$cluster] else "#333"),
                            sprintf("%d", r$cluster))),
          tags$td(r$n_tracks),
          tags$td(tags$b(r$model),
                  # Make it obvious when a forced choice differs from the AICc pick.
                  if (isTRUE(r$forced))
                    tags$span(style = "color:#b26a00;",
                              sprintf(" (set; AICc: %s)", r$auto_model)) else NULL),
          tags$td(if (is.finite(r$D_um2_s)) sprintf("%.3g", r$D_um2_s) else "-"),
          tags$td(r$detail),
          tags$td(if (is.finite(r$dAICc)) sprintf("%.1f", r$dAICc) else "-"))
      })),
    tags$div(style = "color:#888; margin-top:3px;",
      "Each cluster's ensemble MSD is fitted with every model; Auto picks the lowest AICc. ",
      tags$b("\u0394AICc"), " below about 2: the two best models are not distinguishable. It ranks models; it is not a test.")
  )
})

# Cluster eMSD curves with each cluster's fitted model overlaid.
output$trk_msd_cluster_curves <- renderPlotly({
  req(isTRUE(input$trk_show_msd), isTRUE(input$trk_msd_cluster))
  res <- tryCatch(msd_cluster_models(), error = function(e) NULL)
  shiny::validate(shiny::need(!is.null(res), "Click 'Plot MSD' to compute cluster models."))
  cl   <- msd_cluster_rv()
  dark <- isTRUE(isolate(input$trk_dark_mode))
  bg <- if (dark) "#000000" else "#FFFFFF"; fg <- if (dark) "#FFFFFF" else "#333333"
  tb <- res$table
  fig <- plot_ly()
  for (cv in res$curves) {
    col <- cl$palette[cv$cluster]
    lbl <- tb$model[match(cv$cluster, tb$cluster)]
    fig <- fig %>%
      add_trace(x = cv$tau, y = cv$msd / 1e6, type = "scatter", mode = "markers",
                marker = list(size = 5, color = col),
                name = sprintf("Cluster %d eMSD", cv$cluster))
    if (!is.null(cv$fit_tau))
      fig <- fig %>%
        add_trace(x = cv$fit_tau, y = pmax(0, cv$fit_msd) / 1e6,
                  type = "scatter", mode = "lines",
                  line = list(color = col, width = 2, dash = "dash"),
                  name = sprintf("Cluster %d: %s", cv$cluster, lbl))
  }
  fig %>% layout(
    xaxis = list(title = "τ (ms)", color = fg),
    yaxis = list(title = "MSD (µm²)", color = fg, rangemode = "tozero"),
    paper_bgcolor = bg, plot_bgcolor = bg, font = list(color = fg),
    legend = list(font = list(color = fg, size = 10), orientation = "h",
                  x = 0, xanchor = "left", y = -0.34, yanchor = "top"),
    margin = list(l = 55, r = 20, t = 10, b = 110),
    uirevision = "msd_cluster_curves")
})

# Observer: keep cluster_color_map_stored in sync
observe({
  if (!isTRUE(input$trk_msd_cluster_color_tracks)) {
    cluster_color_map_stored(NULL)
    return()
  }
  cl <- tryCatch(msd_cluster_rv(), error = function(e) NULL)
  cluster_color_map_stored(if (!is.null(cl)) cl$col_map else NULL)
})

# Rebuild frames whenever stored cluster map changes
observeEvent(cluster_color_map_stored(), {
  fire_frames()
}, ignoreInit = TRUE, ignoreNULL = FALSE)

# Also rebuild when checkbox is toggled or k changes
observeEvent(input$trk_msd_cluster_color_tracks, {
  fire_frames()
}, ignoreInit = TRUE)
observeEvent(input$trk_msd_n_clusters, {
  if (isTRUE(input$trk_msd_cluster_color_tracks)) fire_frames()
}, ignoreInit = TRUE)
