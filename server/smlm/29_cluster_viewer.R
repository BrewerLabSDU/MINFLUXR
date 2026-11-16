## ─────────────────────────────────────────────────────────────────────────────
## server/smlm/29_cluster_viewer.R  -  panel: cluster viewer and shape fit (exploratory)
## ─────────────────────────────────────────────────────────────────────────────
## Part of smlm_server(). Not a standalone script: its statements are evaluated in
## the server function's environment, in file-name order (server/parts.R).
##
## Snapshot of the clusters to browse, the current cluster, the shape fits
## with model selection, and the viewer plot.
##
## INTERFACE (generated: Rscript tools/server_map.R --headers; do not edit by hand)
##   DEFINES for other parts:
##     (nothing)
##   READS from other parts:
##     01_view_upload: smlm_file_store
##     23_prepared_cache: cluster2_label, deblink_on, prepared_df
## END INTERFACE
## ─────────────────────────────────────────────────────────────────────────────

## ── Cluster Viewer (2D / 3D + shape fit) ──────────────────────────────────
# Inspect one Median Cluster structure at a time. A Run button snapshots the
# clusters; Prev/Next and view options then work on that snapshot, so it does
# not refresh on filter/plot changes until Run is pressed again.
# MLE best-fit via mle_fit_*; overlay from the fitted parameters.
#
# SPANS EVERY LOADED FILE, not just the plotted one. Median Clustering is a
# global setting that applies to all loaded datasets, so the viewer was the odd
# one out in reading only plottable_data(): with eight files loaded and the
# active one a non-clustering control, it reported "no Median Clusters" and
# refused to open - while seven other files were full of them.
#
# Cluster IDs are only unique WITHIN a file, so structures are keyed by
# file + id (`.vkey`). Without that, cluster 3 of file A and cluster 3 of file B
# would be the same entry and the viewer would show their points merged.
viewer_snapshot <- reactiveVal(NULL)   # list(df, ids, keys) captured on Run
viewer_idx      <- reactiveVal(1L)

observeEvent(input$run_viewer, {
  all_files <- !isTRUE(input$viewer_active_only)
  store <- smlm_file_store()
  if (!length(store)) { showNotification("No files loaded.", type = "warning")
                        viewer_snapshot(NULL); return() }

  nms <- names(store)
  if (!all_files) {
    af <- input$active_file
    nms <- if (!is.null(af) && af %in% nms) af else nms[1]
  }

  parts <- list(); skipped <- character(0)
  withProgress(message = "Loading Median Clusters...", value = 0, {
    for (i in seq_along(nms)) {
      nm <- nms[[i]]
      incProgress(1 / length(nms), detail = nm)
      d <- tryCatch(prepared_df(nm, needs_cluster = TRUE), error = function(e) NULL)
      if (is.null(d) || !nrow(d) || !"cluster2" %in% names(d)) { skipped <- c(skipped, nm); next }
      d <- d[!is.na(d$cluster2) & d$cluster2 > 0L, , drop = FALSE]
      if (!nrow(d)) { skipped <- c(skipped, nm); next }
      d$.vfile <- nm
      d$.vkey  <- paste0(nm, "#", d$cluster2)
      parts[[nm]] <- d
    }
  })

  if (!length(parts)) {
    showNotification(sprintf("No Median Clusters in any loaded file. Enable %s and check the filters.",
                             cluster2_label()), type = "warning")
    viewer_snapshot(NULL); return()
  }
  df <- dplyr::bind_rows(parts)
  keys <- df %>%
    dplyr::distinct(.vkey, .vfile, cluster2) %>%
    dplyr::arrange(.vfile, cluster2)

  viewer_snapshot(list(df = df, ids = keys$.vkey, keys = keys))
  viewer_idx(1L)
  msg <- sprintf("Loaded %d Median Cluster structure(s) from %d file(s).",
                 nrow(keys), length(parts))
  if (length(skipped))
    msg <- paste0(msg, sprintf(" %d file(s) had none.", length(skipped)))
  showNotification(msg, type = "message", duration = 4)
})

viewer_ids <- reactive({ s <- viewer_snapshot(); if (is.null(s)) integer(0) else s$ids })

observeEvent(input$viewer_prev, {
  n <- length(viewer_ids()); if (!n) return()
  viewer_idx(((viewer_idx() - 2L) %% n) + 1L)   # wrap-around
})
observeEvent(input$viewer_next, {
  n <- length(viewer_ids()); if (!n) return()
  viewer_idx((viewer_idx() %% n) + 1L)
})

viewer_current <- reactive({
  s   <- viewer_snapshot()
  shiny::validate(shiny::need(!is.null(s) && length(s$ids) > 0,
    "Press 'Run / Load clusters' to load Median Cluster structures from the loaded files."))
  ids <- s$ids; n <- length(ids)
  idx <- min(max(viewer_idx(), 1L), n)
  key <- ids[idx]
  df  <- s$df
  # Subset on the file+id key, never on cluster2 alone - see the note on Run.
  sub <- df[df$.vkey == key, , drop = FALSE]
  id  <- sub$cluster2[1]
  src <- sub$.vfile[1]
  # THE UNIT FOLLOWS DE-BLINKING, like every other panel.
  #
  # This used to group by `tid` unconditionally, so the viewer showed - and the
  # shape fit consumed - trace-ID medians even with de-blinking switched on.
  # That is the same defect the tc2 isosurface had: one fluorophore reappears as
  # several trace IDs, so a molecule that blinked six times contributed six
  # points to the surface and one that blinked once contributed one. The fitted
  # shape then tracks BLINKING as much as geometry, and the long-lived molecules
  # pull hardest. De-blinking exists precisely to collapse that, and the viewer
  # was ignoring it.
  unit_col <- if (deblink_on() && "cluster" %in% names(sub) &&
                  any(!is.na(sub$cluster) & sub$cluster > 0L)) "cluster" else "tid"
  csub <- if (identical(unit_col, "cluster"))
    sub[!is.na(sub$cluster) & sub$cluster > 0L, , drop = FALSE] else sub

  # Photon budget summed over the UNIT (all localizations of the trace, or of the
  # whole de-blinked molecule). MINFLUX localization variance scales as
  # 1/photons and the unit's centroid averages its localizations, so the
  # centroid's variance goes as 1/sum(eco) - the inverse-variance weight for
  # that centroid IS its total photon count.
  cent <- csub %>% dplyr::group_by(.data[[unit_col]]) %>%
    dplyr::summarise(
      x = stats::median(x, na.rm = TRUE),
      y = stats::median(y, na.rm = TRUE),
      z = if ("z" %in% names(csub)) stats::median(z, na.rm = TRUE) else 0,
      w_eco = if ("eco" %in% names(csub)) sum(eco, na.rm = TRUE) else NA_real_,
      .groups = "drop")
  w <- cent$w_eco
  if (!length(w) || all(!is.finite(w)) || sum(w, na.rm = TRUE) <= 0) {
    w <- NULL                                   # no usable eco -> unweighted
  } else {
    w[!is.finite(w) | w <= 0] <- stats::median(w[is.finite(w) & w > 0])
  }
  zc    <- if ("z" %in% names(sub)) sub$z else rep(0, nrow(sub))
  cloud <- cbind(x = sub$x, y = sub$y, z = zc)
  cloud <- cloud[stats::complete.cases(cloud), , drop = FALSE]
  list(idx = idx, n = n, id = id, file = src, sub = sub,
       unit = unit_col,
       mat = as.matrix(cent[, c("x", "y", "z")]),
       w_eco = w,
       cloud = cloud)
})

# Full MLE fit (optim orthogonal-distance) of all 5 shapes to the trace-median
# centroids - same pathway the Display tab uses when fit_method = "mle", not
# the heuristic moment-based geometry. Returns RMSDs, best (min-RMSD) shape and
# the fitted parameters (used to draw the overlay).
# Heavy part: the 5 MLE fits + RMSD + AIC/AICc/BIC. Independent of the ranking
# criterion, so switching AICc<->BIC<->RMSD does NOT re-run the optimisers.
viewer_fit_core <- reactive({
  if (!isTRUE(input$viewer_enable_fit)) return(NULL)   # inspect-only unless fitting enabled
  cur <- viewer_current()
  X   <- cur$mat
  if (nrow(X) < 4) return(NULL)
  turns <- input$helix_turns %||% 3
  cargs <- tryCatch(build_center_args(
    X, mode = input$mle_center_mode %||% "free", lambda = input$mle_center_lambda %||% 0.1),
    error = function(e) list())

  # Inverse-variance (photon) weighting. The default iid-isotropic assumption is
  # wrong in detail for MINFLUX: precision varies with photon count, so a
  # 50-photon trace should not pull the surface as hard as a 5000-photon one.
  # w_i = total eco of trace i (see viewer_current). n_eff is Kish's effective
  # sample size - it exposes the cost of weighting, i.e. whether a couple of very
  # bright traces are now carrying the fit.
  W <- if (isTRUE(input$viewer_weight_eco)) cur$w_eco else NULL
  n_eff <- if (!is.null(W)) sum(W)^2 / sum(W^2) else nrow(X)

  fitS <- tryCatch(do.call(mle_fit_sphere,    c(list(X = X, w = W), cargs)), error = function(e) NULL)
  fitE <- tryCatch(do.call(mle_fit_ellipsoid, c(list(X = X, w = W), cargs)), error = function(e) NULL)
  fitH <- tryCatch(do.call(mle_fit_helix,     c(list(X = X, w = W, turns_hint = turns), cargs)), error = function(e) NULL)
  fitR <- tryCatch(do.call(mle_fit_ring,      c(list(X = X, w = W), cargs)), error = function(e) NULL)
  fitA <- tryCatch(do.call(mle_fit_arc,       c(list(X = X, w = W), cargs)), error = function(e) NULL)

  # Weighted RMSD: sqrt(sum(w*d^2)/sum(w)). With W = NULL this is the plain RMSD.
  rmsd_of <- function(d) {
    ok <- is.finite(d)
    if (!any(ok)) return(NA_real_)
    ww <- if (is.null(W)) rep(1, length(d)) else W
    sqrt(sum(ww[ok] * d[ok]^2) / sum(ww[ok]))
  }
  rS <- if (!is.null(fitS)) rmsd_of(as.numeric(mle_project_sphere(X, fitS$center, fitS$radius)$dist)) else NA_real_
  rE <- if (!is.null(fitE)) rmsd_of(as.numeric(mle_project_ellipsoid(X, fitE$center, fitE$axes, fitE$R)$dist)) else NA_real_
  rR <- if (!is.null(fitR)) rmsd_of(as.numeric(mle_project_ring(X, fitR$center, fitR$normal, fitR$radius)$dist)) else NA_real_
  rA <- if (!is.null(fitA)) rmsd_of(as.numeric(mle_project_arc(X, fitA$center, fitA$normal, fitA$radius, fitA$t_center, fitA$arc_half)$dist)) else NA_real_
  rH <- NA_real_
  if (!is.null(fitH) && !is.null(fitH$t)) {
    t_h <- .project_t_helix(X, fitH$center, fitH$axis, fitH$radius, fitH$pitch, fitH$phi0,
                            t0 = fitH$t, check_neighbors = TRUE)
    H_h <- t(vapply(t_h, .helix_point, numeric(3),
                    c = fitH$center, u = fitH$axis, r = fitH$radius, p = fitH$pitch, phi0 = fitH$phi0))
    rH <- rmsd_of(sqrt(rowSums((X - H_h)^2)))
  }
  rmsd <- c(sphere = rS, ellipsoid = rE, helix = rH, ring = rR, arc = rA)

  # ── model selection: AIC / AICc / BIC under isotropic Gaussian residuals ─────
  # sigma^2_hat = RMSD^2, so -2logL = n*log(2*pi*sigma2) + n. k = geometric DOF + 1
  # (the variance). AICc adds the small-sample correction and -> Inf when a model
  # has too few points to estimate (n <= k+1), which refuses over-parameterised fits.
  # When photon-weighted, point i has variance sigma^2/w_i and -2logL gains a
  # -sum(log w_i) term - identical for every shape (same weights), so it cancels in
  # every AIC/AICc/BIC *difference*. The criteria stay comparable across shapes;
  # only their absolute values shift, so never compare them across weighting modes.
  n   <- nrow(X)
  pk  <- c(sphere = 4, ellipsoid = 9, helix = 8, ring = 6, arc = 8)  # geometric DOF
  shp <- names(rmsd)
  M <- t(vapply(shp, function(s) {
    r <- rmsd[[s]]
    if (!is.finite(r)) return(c(rmsd = r, aic = NA_real_, aicc = NA_real_, bic = NA_real_, params = pk[[s]]))
    sigma2 <- max(r, 1e-6)^2
    k    <- pk[[s]] + 1L                       # + variance parameter
    n2ll <- n * log(2 * pi) + n * log(sigma2) + n
    aic  <- n2ll + 2 * k
    aicc <- if (n - k - 1L > 0) aic + (2 * k * (k + 1)) / (n - k - 1L) else Inf
    bic  <- n2ll + k * log(n)
    c(rmsd = r, aic = aic, aicc = aicc, bic = bic, params = pk[[s]])
  }, numeric(5)))
  rownames(M) <- shp

  list(rmsd = rmsd, M = M, n = n,
       weighted = !is.null(W), n_eff = n_eff,
       fits = list(sphere = fitS, ellipsoid = fitE, helix = fitH, ring = fitR, arc = fitA))
})

# Cheap part: pick the best shape for the chosen criterion (just a which.min).
viewer_fit <- reactive({
  core <- viewer_fit_core()
  if (is.null(core)) return(NULL)
  M <- core$M
  crit_col <- switch(input$viewer_rank %||% "aicc",
                     aicc = "aicc", bic = "bic", rmsd = "rmsd", "aicc")
  pick_best <- function(col) {
    v <- M[, col]; v[!is.finite(v)] <- Inf
    if (all(!is.finite(v))) return(NA_character_)
    rownames(M)[which.min(v)]
  }
  best <- pick_best(crit_col)
  if (is.na(best) && crit_col == "aicc") best <- pick_best("aic")   # tiny clusters
  if (is.na(best)) best <- pick_best("rmsd")
  c(core, list(best = best, criterion = crit_col))
})



output$viewer_header <- renderText({
  cur  <- viewer_current()
  vf   <- viewer_fit()
  best <- if (!is.null(vf)) vf$best else NA_character_
  rb   <- if (!is.null(vf) && !is.na(best)) as.numeric(vf$rmsd[[best]]) else NA_real_
  crit <- toupper(if (!is.null(vf)) vf$criterion else "aicc")
  wtxt <- if (!is.null(vf) && isTRUE(vf$weighted))
    sprintf("   |   photon-weighted (n_eff %.1f of %d)", vf$n_eff, vf$n) else ""
  unit_lbl <- if (identical(cur$unit, "cluster")) "de-blinked molecules" else "trace IDs"
  sprintf("Cluster %d / %d   |   %s   |   id %d   |   %d %s%s%s",
          cur$idx, cur$n, cur$file %||% "?", cur$id, nrow(cur$mat), unit_lbl,
          if (is.finite(rb)) sprintf("   |   best MLE fit (%s): %s (%sRMSD %.1f nm)", crit, best,
                                     if (!is.null(vf) && isTRUE(vf$weighted)) "weighted " else "", rb) else "",
          wtxt)
})

output$viewer_rmsd_table <- renderTable({
  vf <- viewer_fit()
  shiny::validate(shiny::need(!is.null(vf), "Need >= 4 traces for shape fitting."))
  M      <- vf$M
  shapes <- rownames(M)
  fnum   <- function(v) ifelse(is.finite(v), round(v, 1), NA_real_)
  data.frame(
    Shape        = shapes,
    params       = as.integer(M[, "params"]),
    `RMSD (nm)`  = round(M[, "rmsd"], 2),
    AICc         = fnum(M[, "aicc"]),
    BIC          = fnum(M[, "bic"]),
    Best         = ifelse(shapes == vf$best, "<-- best", ""),
    check.names  = FALSE
  )
}, striped = TRUE, spacing = "xs")

output$viewer_plot <- renderPlotly({
  cur   <- viewer_current()
  mat   <- cur$mat
  cloud <- cur$cloud
  # Series name follows the unit, so the legend cannot claim "trace medians"
  # while showing de-blinked molecule centroids.
  unit_nm <- if (identical(cur$unit, "cluster")) "De-blinked molecules" else "Trace medians"
  mode  <- input$viewer_mode %||% "3d"

  vf    <- viewer_fit()
  sel   <- input$viewer_fit_shape %||% "best"
  shape <- if (sel == "none") NA_character_
           else if (sel == "best") { if (!is.null(vf)) vf$best else NA_character_ }
           else sel
  f          <- if (!is.na(shape) && !is.null(vf)) vf$fits[[shape]] else NULL
  is_surface <- !is.na(shape) && shape %in% c("sphere", "ellipsoid")
  is_curve   <- !is.na(shape) && shape %in% c("ring", "arc", "helix")
  fit_col    <- "#D7301F"

  if (identical(mode, "pca")) {
    shiny::validate(shiny::need(nrow(mat) >= 3, "Need >= 3 units for a PCA projection."))
    # PCA basis from the UNIT centroids so unit length (localization count)
    # does not weight the axes; the full cloud is then projected onto that basis.
    ctr  <- colMeans(mat)
    pca  <- stats::prcomp(sweep(mat, 2, ctr), center = FALSE, scale. = FALSE)
    proj <- function(P) sweep(P, 2, ctr) %*% pca$rotation
    cl   <- proj(cloud); md <- proj(mat)
    fig  <- plotly::plot_ly() %>%
      plotly::add_markers(x = cl[, 1], y = cl[, 2],
        marker = list(color = "#2C7FB8", size = 3, opacity = 0.30),
        name = "Localizations") %>%
      plotly::add_markers(x = md[, 1], y = md[, 2],
        marker = list(color = "#08306B", size = 7, line = list(color = "white", width = 0.5)),
        name = unit_nm)
    if (!is.null(f) && is_curve) {
      cp <- .viewer_curve_points(f, shape)
      if (!is.null(cp)) { ps <- proj(cp)
        fig <- fig %>% plotly::add_trace(x = ps[, 1], y = ps[, 2],
          type = "scatter", mode = "lines",
          line = list(color = fit_col, width = 3), name = paste("Fit:", shape)) }
    } else if (!is.null(f) && is_surface) {
      mh <- .viewer_surface_mesh(f, shape)
      if (!is.null(mh)) { pv <- proj(mh$verts); h <- grDevices::chull(pv)
        poly <- pv[c(h, h[1]), , drop = FALSE]                # closed silhouette
        fig <- fig %>% plotly::add_trace(x = poly[, 1], y = poly[, 2],
          type = "scatter", mode = "lines", fill = "toself",
          fillcolor = "rgba(215,48,31,0.12)",
          line = list(color = fit_col, width = 3), name = paste("Fit:", shape)) }
    }
    fig %>% plotly::layout(
      xaxis  = list(title = "PC1 (nm)", scaleanchor = "y", scaleratio = 1),
      yaxis  = list(title = "PC2 (nm)"),
      legend = list(orientation = "h"))
  } else {
    fig <- plotly::plot_ly() %>%
      plotly::add_trace(x = cloud[, 1], y = cloud[, 2], z = cloud[, 3],
        type = "scatter3d", mode = "markers",
        marker = list(color = "#2C7FB8", size = 2, opacity = 0.30), name = "Localizations") %>%
      plotly::add_trace(x = mat[, 1], y = mat[, 2], z = mat[, 3],
        type = "scatter3d", mode = "markers",
        marker = list(color = "#08306B", size = 4), name = unit_nm)
    if (!is.null(f) && is_curve) {
      cp <- .viewer_curve_points(f, shape)
      if (!is.null(cp))
        fig <- fig %>% plotly::add_trace(x = cp[, 1], y = cp[, 2], z = cp[, 3],
          type = "scatter3d", mode = "lines",
          line = list(color = fit_col, width = 5), name = paste("Fit:", shape))
    } else if (!is.null(f) && is_surface) {
      mh <- .viewer_surface_mesh(f, shape)
      if (!is.null(mh))
        fig <- fig %>% plotly::add_trace(type = "mesh3d",
          x = mh$verts[, 1], y = mh$verts[, 2], z = mh$verts[, 3],
          i = mh$i, j = mh$j, k = mh$k,
          facecolor = rep(fit_col, length(mh$i)), opacity = 0.20,
          flatshading = FALSE, hoverinfo = "skip",
          name = paste("Fit:", shape), showlegend = TRUE)
    }
    fig %>% plotly::layout(
      scene = list(xaxis = list(title = "x (nm)"),
                   yaxis = list(title = "y (nm)"),
                   zaxis = list(title = "z (nm)"),
                   aspectmode = "data"),
      legend = list(orientation = "h"))
  }
})
