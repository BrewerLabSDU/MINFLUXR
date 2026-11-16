## ── Median Cluster (cluster2) Crowding vs Detection-Time Shift ─────────────────────────────
# One row per Cluster2 structure (built from its member centroids - de-blinked
# molecules when de-blinking is on, else trace IDs):
#   nn_median  = median k-NN distance among members (smaller = denser) ± MAD
#   time_shift = median(member tim_med) - median(unclustered background tim_med),
#                computed per file (positive = structure detected later)
# Tests whether denser (more crowded) structures are detected later.
crowd_results    <- reactiveVal(NULL)
crowd_null_stats <- reactiveVal(NULL)   # list(rho_obs, p, B, outside) from the permutation null

observeEvent(input$run_crowd, {
  store <- smlm_file_store()
  if (!length(store)) { showNotification("No files loaded.", type = "warning"); return() }
  if (!cluster2_on()) {
    showNotification(sprintf("Enable %s in the Display tab first.", cluster2_label()),
                     type = "warning"); return()
  }
  if (!isTRUE(input$crowd_all_files)) {
    af <- input$active_file
    if (!is.null(af) && af %in% names(store)) store <- store[af]
  }

  min_n  <- max(3L, as.integer(input$crowd_min_members %||% 4L))
  kk     <- max(1L, as.integer(input$crowd_k %||% 1L))
  use_3d <- isTRUE(input$crowd_use_3d)

  # ── no-crowding null (label permutation) ──────────────────────────────────
  # The panel's y-axis has a trend built in even with NO crowding effect: denser
  # clusters hold more traces, and the MEDIAN of a larger sample drawn from a
  # decaying detection-rate distribution lands later purely by sampling. To make
  # the plot readable we need the "zero crowding" expectation drawn on it.
  #
  # WHY A PERMUTATION AND NOT THE SIMULATOR: shuffling the detection times among
  # the traces of a file keeps the real cluster membership, the real sizes, the
  # real geometry (so every point's x stays exactly where it is) and the real
  # detection-time distribution of that acquisition - it destroys ONLY the
  # association being tested. A simulator null would instead require assuming a
  # bleaching/detection model, which is itself calibrated from data of this kind;
  # that would make the test partly circular. The permutation is exact,
  # assumption-free and cheaper.
  do_null <- isTRUE(input$crowd_null)
  B       <- if (do_null) max(49L, min(2000L, as.integer(input$crowd_nsim %||% 200L))) else 0L

  rows      <- list()
  null_rows <- list()   # per file: structures x B matrix of null time shifts
  # Why-it-found-nothing counters. An empty result has several distinct causes,
  # so each one is counted and named in the warning.
  diag <- c(no_tim = 0L, filtered_out = 0L, no_c2 = 0L, no_z = 0L,
            no_clusters = 0L, no_background = 0L, too_small = 0L)
  withProgress(message = "Computing crowding vs detection-time shift...", value = 0, {
    n_files <- length(store)
    for (i in seq_along(store)) {
      fname <- names(store)[[i]]
      incProgress(1 / n_files, detail = fname)
      df <- qc_exclude_from(store[[fname]], fname)
      if (!is.data.frame(df) || !nrow(df) || !"tim" %in% names(df)) { diag["no_tim"] <- diag["no_tim"] + 1L; next }
      df  <- prepared_df(fname)   # cached filter -> thin -> cluster chain
      got <- crowding_for_file(df, fname, id_col = smlm_nn_unit(df), min_n = min_n,
                               k = kk, use_3d = use_3d, B = B)
      diag[names(got$diag)] <- diag[names(got$diag)] + got$diag
      if (!is.null(got$rows)) rows[[length(rows) + 1L]] <- got$rows
      if (!is.null(got$null)) null_rows[[length(null_rows) + 1L]] <- got$null
    }
  })

  res <- if (length(rows)) do.call(rbind, rows) else NULL
  if (is.null(res) || !nrow(res)) {
    why <- c(
      if (diag["no_tim"])        sprintf("%d file(s) have no `tim` column", diag["no_tim"]),
      if (diag["filtered_out"])  sprintf("%d file(s) left empty by the Display-tab metadata filters", diag["filtered_out"]),
      if (diag["no_c2"])         sprintf("%d file(s) produced no cluster2 column (is Median Clustering enabled?)", diag["no_c2"]),
      if (diag["no_z"])          sprintf("%d file(s) have no z (uncheck 'Use 3-D coordinates')", diag["no_z"]),
      if (diag["no_clusters"])   sprintf("%d file(s) yielded no Median Cluster structures at all", diag["no_clusters"]),
      if (diag["no_background"]) sprintf("%d file(s) had no UNCLUSTERED traces to use as the time reference", diag["no_background"]),
      if (diag["too_small"])     sprintf("%d structure(s) had fewer than %d traces", diag["too_small"], min_n))
    showNotification(
      HTML(paste0("<b>No structures to plot.</b><br>",
                  if (length(why)) paste(why, collapse = "<br>") else
                    "No qualifying Median Cluster structures found.")),
      type = "warning", duration = 12)
    crowd_results(NULL); crowd_null_stats(NULL); return()
  }
  # ── attach the null band + the global trend test ──────────────────────────
  tst <- crowding_null_test(res, if (do_null && length(null_rows)) do.call(rbind, null_rows) else NULL)
  res <- tst$res
  crowd_null_stats(tst$stats)

  crowd_results(res)
  ns <- crowd_null_stats()
  showNotification(
    sprintf("Crowding vs detection-time shift: %d structures across %d file(s).%s",
            nrow(res), length(unique(res$dataset)),
            if (!is.null(ns) && is.finite(ns$p))
              sprintf(" Trend vs no-crowding null: rho = %.3f, p = %.3f (%d permutations).",
                      ns$rho_obs, ns$p, ns$B) else ""),
    type = "message", duration = 6)
})

build_crowd_plot_fig <- function() {
  res <- crowd_results()
  shiny::validate(shiny::need(!is.null(res) && nrow(res) > 0,
                              "Press 'Run' to compute crowding vs detection-time shift."))
  n_all <- nrow(res)
  res <- apply_dataset_groups(res)
  shiny::validate(shiny::need(nrow(res) > 0, "No files assigned to groups (group mode)."))
  # Group mode can drop files; the stored rho/p were computed over the full run,
  # so only quote them when every structure is still on screen.
  stats_match <- nrow(res) == n_all
  fig <- plotly::plot_ly()
  datasets <- unique(res$dataset)
  n_ds     <- length(datasets)
  pal <- series_palette(datasets)

  # ── no-crowding null band (drawn first, so points sit on top) ─────────────
  # Each cluster has its own null interval (its size is fixed), so the band is
  # per-cluster quantiles ordered along x and lightly smoothed for readability.
  if ("null_lo" %in% names(res) && any(is.finite(res$null_lo))) {
    nb  <- res[is.finite(res$null_lo) & is.finite(res$nn_median), , drop = FALSE]
    nb  <- nb[order(nb$nn_median), , drop = FALSE]
    smooth_band <- function(y) {
      if (nrow(nb) < 8) return(y)
      tr <- tryCatch(stats::loess(y ~ nb$nn_median, span = 0.9, degree = 1), error = function(e) NULL)
      if (is.null(tr)) return(y)
      p <- stats::predict(tr); ifelse(is.finite(p), p, y)
    }
    lo <- smooth_band(nb$null_lo); hi <- smooth_band(nb$null_hi); md <- smooth_band(nb$null_med)
    fig <- fig %>%
      plotly::add_lines(x = nb$nn_median, y = hi, name = "No-crowding null (95%)",
                        line = list(color = "rgba(120,120,120,0)", width = 0),
                        hoverinfo = "skip", showlegend = FALSE) %>%
      plotly::add_lines(x = nb$nn_median, y = lo, name = "No-crowding null (95%)",
                        fill = "tonexty", fillcolor = "rgba(120,120,120,0.18)",
                        line = list(color = "rgba(120,120,120,0)", width = 0),
                        hoverinfo = "skip") %>%
      plotly::add_lines(x = nb$nn_median, y = md, name = "Null median",
                        line = list(color = "#777777", width = 1.5, dash = "dot"),
                        hoverinfo = "skip")
  }

  for (ds in datasets) {
    sub <- res[res$dataset == ds, , drop = FALSE]
    ex  <- if ("nn_mad" %in% names(sub)) list(type = "data", array = sub$nn_mad,
                                              color = pal[ds], thickness = 1, width = 0,
                                              opacity = 0.4) else NULL
    fig <- fig %>% plotly::add_markers(
      x = sub$nn_median, y = sub$time_shift,
      name = if (n_ds > 1) ds else "Structures",
      error_x = ex,   # ± MAD of the within-cluster NN spacing
      marker = list(color = pal[ds], size = 7, opacity = 0.7,
                    line = list(color = "white", width = 0.5)),
      text = paste0("Median cluster ", sub$cluster2, "<br>n = ", sub$n_tids,
                    " ", if ("unit" %in% names(sub))
                           ifelse(sub$unit == "molecule", "de-blinked molecules", "trace IDs") else "",
                    "<br>NN median = ", round(sub$nn_median, 1), " ± ",
                    round(if ("nn_mad" %in% names(sub)) sub$nn_mad else 0, 1), " nm",
                    "<br>time shift = ", round(sub$time_shift, 2)),
      hoverinfo = "text"
    )
  }

  # optional loess trend - per group when grouped, else over the pooled cloud
  add_crowd_loess <- function(fig, d, col, nm) {
    if (!(isTRUE(input$crowd_trend) && nrow(d) >= 8)) return(fig)
    ord <- order(d$nn_median)
    tr  <- tryCatch(stats::loess(time_shift ~ nn_median, data = d,
                     span = 0.9, degree = 1, family = "symmetric"), error = function(e) NULL)
    if (is.null(tr)) return(fig)
    xx <- d$nn_median[ord]; yy <- stats::predict(tr, newdata = data.frame(nn_median = xx)); ok <- is.finite(yy)
    if (sum(ok) >= 2) fig <- fig %>% plotly::add_lines(x = xx[ok], y = yy[ok], name = nm,
                               line = list(color = col, width = 2))
    fig
  }
  if ((input$analysis_by %||% "file") == "group" && length(datasets) > 1) {
    for (ds in datasets) fig <- add_crowd_loess(fig, res[res$dataset == ds, , drop = FALSE], unname(pal[ds]), paste(ds, "loess"))
  } else {
    fig <- add_crowd_loess(fig, res, "#D7301F", "Trend (loess)")
  }

  # Verdict line: the permutation test on the TREND (rho), plus how many points
  # left their own null interval. Stated explicitly so the slope is never read
  # against zero - the null band, not zero, is the correct reference here.
  ns  <- crowd_null_stats()
  ttl <- if (!is.null(ns) && is.finite(ns$p) && stats_match) {
    list(text = sprintf(
      "Crowding vs detection-time shift  |  trend vs no-crowding null: Spearman rho = %.3f, p = %.3f (%d permutations)  |  %.0f%% of structures outside their null band (5%% expected)",
      ns$rho_obs, ns$p, ns$B, 100 * ns$outside), font = list(size = 12))
  } else if (!is.null(ns) && is.finite(ns$p)) {
    list(text = "Crowding vs detection-time shift  |  null band shown per structure; trend test not shown (group filter changed the set)",
         font = list(size = 12))
  } else NULL

  fig %>% plotly::layout(
    title  = ttl,
    margin = list(l = 70, t = if (is.null(ttl)) 40 else 60),
    xaxis  = list(title = "Median NN distance among members (nm)  -  smaller = more crowded"),
    yaxis  = list(title = "Detection-time shift (tim)",
                  zeroline = TRUE, zerolinecolor = "#888888", zerolinewidth = 1),
    legend = list(orientation = "h")
  )
}
output$crowd_plot <- renderPlotly(build_crowd_plot_fig())

output$download_crowd <- downloadHandler(
  filename = function()
    paste0("cluster2_crowding_vs_detection_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
  content = function(file) {
    res <- crowd_results()
    shiny::validate(shiny::need(!is.null(res) && nrow(res) > 0, "No results."))
    write.csv(res, file, row.names = FALSE)
  }
)

## ── Within-Cluster2 NN vs Cluster Size ────────────────────────────────────
# One row per Cluster2: size = # member units (trace IDs, or de-blinked molecules
# if de-blinking is on); nn_median = median 1-NN distance among the
# member centroids. Tests whether spacing is constant with size or "balloons".
size_nn_results <- reactiveVal(NULL)
size_nn_params  <- reactiveVal(NULL)   # compute params, frozen at Run

observeEvent(input$run_size_nn, {
  store <- smlm_file_store()
  if (!length(store)) { showNotification("No files loaded.", type = "warning"); return() }
  if (!cluster2_on()) {
    showNotification(sprintf("Enable %s in the Display tab first.", cluster2_label()),
                     type = "warning"); return()
  }
  if (!isTRUE(input$size_nn_all_files)) {
    af <- input$active_file
    if (!is.null(af) && af %in% names(store)) store <- store[af]
  }
  min_n       <- max(2L, as.integer(input$size_nn_min_members %||% 3L))
  kk          <- max(1L, as.integer(input$size_nn_k %||% 1L))
  use_3d      <- isTRUE(input$size_nn_use_3d)

  rows <- list()
  # Same why-it-found-nothing counters as the crowding panel: an empty result has
  # several distinct causes and the old message named only the size threshold.
  diag <- c(empty = 0L, filtered_out = 0L, no_c2 = 0L, no_z = 0L,
            no_clusters = 0L, too_small = 0L, degenerate = 0L)
  withProgress(message = "Computing within-cluster NN vs size...", value = 0, {
    nf <- length(store)
    for (i in seq_along(store)) {
      fname <- names(store)[[i]]; incProgress(1 / nf, detail = fname)
      df <- qc_exclude_from(store[[fname]], fname)
      if (!is.data.frame(df) || !nrow(df)) { diag["empty"] <- diag["empty"] + 1L; next }
      df  <- prepared_df(fname)   # cached filter -> thin -> cluster chain
      got <- size_nn_for_file(df, fname, id_col = smlm_nn_unit(df), min_n = min_n,
                              k = kk, use_3d = use_3d)
      diag[names(got$diag)] <- diag[names(got$diag)] + got$diag
      if (!is.null(got$rows)) rows[[length(rows) + 1L]] <- got$rows
    }
  })

  res <- if (length(rows)) do.call(rbind, rows) else NULL
  if (is.null(res) || !nrow(res)) {
    why <- c(
      if (diag["empty"])        sprintf("%d file(s) were empty", diag["empty"]),
      if (diag["filtered_out"]) sprintf("%d file(s) left empty by the Display-tab metadata filters", diag["filtered_out"]),
      if (diag["no_c2"])        sprintf("%d file(s) produced no cluster2 column (is Median Clustering enabled?)", diag["no_c2"]),
      if (diag["no_z"])         sprintf("%d file(s) have no z (uncheck 'Use 3-D coordinates')", diag["no_z"]),
      if (diag["no_clusters"])  sprintf("%d file(s) yielded no Median Cluster structures at all", diag["no_clusters"]),
      if (diag["too_small"])    sprintf("%d structure(s) had fewer than %d members", diag["too_small"], min_n),
      if (diag["degenerate"])   sprintf("%d structure(s) had fewer than 2 usable coordinates", diag["degenerate"]))
    showNotification(
      HTML(paste0("<b>No structures to plot.</b><br>",
                  if (length(why)) paste(why, collapse = "<br>") else
                    "No qualifying Median Cluster structures found.")),
      type = "warning", duration = 12)
    size_nn_results(NULL); return()
  }
  size_nn_results(res)
  # Compute parameters frozen at Run: reading them live relabelled an existing
  # plot without recomputing it, so the axis could describe a different unit
  # than the data on screen.
  size_nn_params(list(k = kk, xmode = input$size_nn_xmode %||% "count"))
  rho <- suppressWarnings(stats::cor(res$size, res$nn_median, method = "spearman"))
  showNotification(
    sprintf("Within-cluster NN vs size: %d structures across %d file(s); Spearman rho = %.2f.",
            nrow(res), length(unique(res$dataset)), rho),
    type = "message", duration = 4)
})

build_size_nn_plot_fig <- function() {
  res <- size_nn_results()
  shiny::validate(shiny::need(!is.null(res) && nrow(res) > 0,
                              "Press 'Run' to compute within-cluster NN vs size."))
  res <- apply_dataset_groups(res)
  shiny::validate(shiny::need(nrow(res) > 0, "No files assigned to groups (group mode)."))
  # res$unit holds "molecule" / "tid" (written by smlm_nn_unit above). The old
  # test compared against "primary_cluster", a name from when de-blinking was
  # raw-point DBSCAN, so it never matched and the label always read "trace IDs"
  # even when the panel was plotting de-blinked molecules.
  unit_lbl <- if (all(res$unit == "molecule")) "de-blinked molecules" else "trace IDs"
  xmode    <- (size_nn_params()$xmode) %||% (input$size_nn_xmode %||% "count")

  if (identical(xmode, "hull")) {
    res <- res[is.finite(res$hull_size), , drop = FALSE]
    shiny::validate(shiny::need(nrow(res) > 0,
      "No clusters with a valid convex hull (need >= 3 members in 2D / >= 4 in 3D)."))
    res$xv  <- res$hull_size
    is3d    <- all(res$dim == 3)
    xlab    <- if (is3d) "Convex-hull volume (nm^3)" else "Convex-hull area (nm^2)"
    xname   <- if (is3d) "hull vol" else "hull area"
    slope_u <- if (is3d) "nm per nm^3" else "nm per nm^2"
  } else {
    res$xv  <- res$size
    xlab    <- sprintf("Cluster size (# %s)", unit_lbl)
    xname   <- "size"
    slope_u <- "nm/member"
  }

  rho <- suppressWarnings(stats::cor(res$xv, res$nn_median, method = "spearman"))
  sl  <- tryCatch(unname(stats::coef(stats::lm(nn_median ~ xv, data = res))[2]),
                  error = function(e) NA_real_)

  fig <- plotly::plot_ly()
  datasets <- unique(res$dataset); n_ds <- length(datasets)
  pal <- series_palette(datasets)

  for (ds in datasets) {
    sub <- res[res$dataset == ds, , drop = FALSE]
    ey  <- if ("nn_mad" %in% names(sub)) list(type = "data", array = sub$nn_mad,
                                              color = pal[ds], thickness = 1, width = 0,
                                              opacity = 0.4) else NULL
    fig <- fig %>% plotly::add_markers(x = sub$xv, y = sub$nn_median,
      name = if (n_ds > 1) ds else "Structures",
      error_y = ey,   # ± MAD of the within-cluster NN spacing
      marker = list(color = pal[ds], size = 7, opacity = 0.7, line = list(color = "white", width = 0.5)),
      text = paste0("Median cluster ", sub$cluster2, "<br>", xname, " ", signif(sub$xv, 4),
                    "<br>median NN ", round(sub$nn_median, 1), " ± ",
                    round(if ("nn_mad" %in% names(sub)) sub$nn_mad else 0, 1), " nm<br>", sub$dataset),
      hoverinfo = "text")
  }

  grouped <- (input$analysis_by %||% "file") == "group" && length(datasets) > 1

  # dashed linear fit (always) + optional loess, over a data subset
  add_trends <- function(fig, d, col, prefix) {
    if (nrow(d) < 3) return(fig)
    sl_d  <- tryCatch(unname(stats::coef(stats::lm(nn_median ~ xv, data = d))[2]),
                      error = function(e) NA_real_)
    rho_d <- suppressWarnings(stats::cor(d$xv, d$nn_median, method = "spearman"))
    xr    <- range(d$xv, na.rm = TRUE)
    yr    <- tryCatch(stats::predict(stats::lm(nn_median ~ xv, data = d),
                      newdata = data.frame(xv = xr)), error = function(e) NULL)
    if (!is.null(yr)) {
      nm <- if (nzchar(prefix))
              sprintf("%s fit (slope %s, rho %.2f)", prefix,
                      if (is.finite(sl_d)) signif(sl_d, 3) else "NA", rho_d)
            else "Linear fit"
      fig <- fig %>% plotly::add_lines(x = xr, y = yr, name = nm,
               line = list(color = col, width = 2, dash = "dash"))
    }
    if (isTRUE(input$size_nn_trend) && nrow(d) >= 8) {
      ord <- order(d$xv)
      tr  <- tryCatch(stats::loess(nn_median ~ xv, data = d, span = 0.9, degree = 1, family = "symmetric"),
                      error = function(e) NULL)
      if (!is.null(tr)) {
        xx <- d$xv[ord]; yy <- stats::predict(tr, newdata = data.frame(xv = xx)); ok <- is.finite(yy)
        if (sum(ok) >= 2)
          fig <- fig %>% plotly::add_lines(x = xx[ok], y = yy[ok],
                   name = if (nzchar(prefix)) paste(prefix, "loess") else "Trend (loess)",
                   line = list(color = col, width = 2))
      }
    }
    fig
  }

  if (grouped) {
    for (ds in datasets) fig <- add_trends(fig, res[res$dataset == ds, , drop = FALSE], unname(pal[ds]), ds)
    ttl <- "Within-cluster NN vs size - per-group fits (slope & rho in legend)"
  } else {
    fig <- add_trends(fig, res, "#555555", "")
    ttl <- sprintf("Spearman rho = %.2f   |   linear slope = %s %s",
                   rho, if (is.finite(sl)) signif(sl, 3) else "NA", slope_u)
  }
  kk   <- max(1L, as.integer((size_nn_params()$k) %||% (input$size_nn_k %||% 1L)))
  ylab <- if (kk == 1L) "Within-cluster median NN distance (nm)"
          else sprintf("Within-cluster median spacing (mean of %d nearest, nm)", kk)
  fig %>% plotly::layout(
    title  = list(text = ttl, font = list(size = 13)),
    xaxis  = list(title = xlab,
                  type = if (isTRUE(input$size_nn_logx)) "log" else "linear"),
    yaxis  = list(title = ylab),
    legend = list(orientation = "h"))
}
output$size_nn_plot <- renderPlotly(build_size_nn_plot_fig())

output$size_nn_small_n <- renderUI({
  res <- size_nn_results()
  if (is.null(res) || !nrow(res)) return(NULL)
  # `size` is this panel's member count (see the run_size_nn result frame).
  nv <- if ("size" %in% names(res)) res$size else
        if ("n_tids" %in% names(res)) res$n_tids else NULL
  if (is.null(nv)) return(NULL)
  .small_n_note(nv, "structure")
})

output$download_size_nn <- downloadHandler(
  filename = function() paste0("within_cluster2_nn_vs_size_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
  content = function(file) {
    res <- size_nn_results()
    shiny::validate(shiny::need(!is.null(res) && nrow(res) > 0, "No results."))
    write.csv(res, file, row.names = FALSE)
  }
)
