# --- Main 2D plot -----------------------------------------------------------

# Shared "exploratory + references" footer for the GT / shape read-outs.
gt_exploratory_note <- function() {
  div(style = "margin-top:8px;padding:6px 8px;border-left:3px solid #e67e22;background:rgba(230,126,34,0.06);font-size:11.5px;color:#7f8c8d;line-height:1.4;",
      HTML(paste0(
        "<b>Be aware:</b> This is a model-based estimate, not a direct measurement; validate against your own controls. ",
        "Methods: gyration-tensor shape anisotropy (Rudnick &amp; Gaspari, <i>J. Phys. A</i> <b>19</b>, L191, 1986); ",
        "linear/planar/spherical anisotropy coefficients (Westin et&nbsp;al., <i>MICCAI/Med. Image Anal.</i>, 1997/2002); ",
        "minimum-spanning-tree backbone (Prim, 1957). ",
        "The labelling inverse is a heuristic of this app, checked against the simulator rather than a known structure. Cite it as exploratory.")))
}

# ── GT counting-standard calibration ───────────────────────────────────────
# Populate the standard-dataset selector from loaded files.
observeEvent(smlm_file_store(), {
  fs  <- names(smlm_file_store())
  sel <- isolate(input$gt_calib_file)
  updateSelectInput(session, "gt_calib_file", choices = fs,
    selected = if (!is.null(sel) && sel %in% fs) sel else (isolate(input$active_file) %||% fs[1]))
}, ignoreNULL = FALSE)

# Measure observed TIDs-per-molecule from the chosen standard: apply the SAME
# metadata filters + clustering settings the sample uses, then measure per
# Median Cluster (each = one standard structure of known stoichiometry).
gt_calib <- reactiveVal(NULL)
observeEvent(input$gt_calib_run, {
  run_set <- smlm_settings_from_input(input)   # one settings snapshot for this run
  store <- smlm_file_store(); fn <- input$gt_calib_file
  if (is.null(fn) || !nzchar(fn %||% "") || !fn %in% names(store)) {
    showNotification("Pick a loaded standard dataset first.", type = "warning"); return()
  }
  k <- max(1L, as.integer(input$gt_calib_k %||% 1L))
  res <- tryCatch({
    df <- apply_clustering_pipeline(apply_meta_filters_batch(store[[fn]], run_set), run_set)
    measure_label_multiplicity(df, k_per_structure = k)
  }, error = function(e) { message("[GT-calib] ", conditionMessage(e)); NULL })
  if (is.null(res) || !isTRUE(res$ok)) {
    gt_calib(NULL)
    showNotification("No Median Clusters in the standard. Enable de-blinking and Median Clustering, then re-run.",
                     type = "warning", duration = 6); return()
  }
  res$file <- fn; gt_calib(res)
  showNotification(sprintf("Standard '%s': %.1f TIDs/molecule over %d structure(s).",
                           fn, res$m_tid, res$n_struct), type = "message", duration = 5)
})

output$gt_calib_readout <- renderUI({
  r <- gt_calib()
  if (is.null(r) || !isTRUE(r$ok))
    return(helpText(style = "color:#888;", "Measure a standard to calibrate the count."))
  active <- isTRUE(input$gt_use_calib)
  div(style = sprintf("font-size:12px; margin:4px 0; padding:5px 8px; border-left:3px solid %s; background:%s;",
                      if (active) "#1a8a1a" else "#b0b8c0", if (active) "rgba(26,138,26,0.06)" else "#f6f7f8"),
    HTML(sprintf("Measured <b>%.1f</b> TIDs / molecule (MAD %.1f, IQR %.0f&ndash;%.0f) from <b>%d</b> structure%s in <b>%s</b>%s. %s",
      r$m_tid, r$m_tid_mad, r$m_tid_iqr[1], r$m_tid_iqr[2], r$n_struct,
      if (r$n_struct == 1) "" else "s", r$file,
      if (r$k > 1) sprintf(" (÷ %d molecules each)", r$k) else "",
      if (active) "<b style=\"color:#1a8a1a\">In use for the count.</b>" else "Tick the box below to use it.")))
})

# ── Inverse GT-position estimate (labeling-aware) ───────────────────────────
# Shared by the 2D/3D overlays and the readout panel. Returns NULL when the
# feature is off or there is nothing to estimate.
gt_estimate <- reactive({
  if (!isTRUE(input$show_gt_estimate)) return(NULL)
  df <- plottable_data()
  if (is.null(df) || !nrow(df)) return(NULL)

  presets     <- minflux_label_presets()
  preset_name <- input$gt_label_preset %||% "Primary monoclonal IgG + secondary IgG"
  preset      <- presets[[preset_name]]
  if (is.null(preset)) return(NULL)

  # Custom: build the preset from the numeric inputs.
  if (identical(preset_name, "Custom")) {
    preset <- list(
      n_fluor       = max(1, suppressWarnings(as.numeric(input$gt_fluor_custom   %||% 1))),
      link_nm       = max(0, suppressWarnings(as.numeric(input$gt_link_custom     %||% 5))),
      dye_spread_nm = max(0, suppressWarnings(as.numeric(input$gt_spread_custom   %||% 0))),
      footprint_nm  = max(0, suppressWarnings(as.numeric(input$gt_footprint_custom %||% 8))),
      blink_lambda  = max(0, suppressWarnings(as.numeric(input$gt_blink_custom    %||% 0.3))),
      note          = "Custom labeling parameters"
    )
  }

  eff <- suppressWarnings(as.numeric(input$gt_label_eff)); if (!is.finite(eff)) eff <- 70
  lvl <- input$gt_model_level %||% "auto"
  # Merge radius: blank / NA / negative -> auto (NULL); 0 -> no merge; >0 -> override.
  cr_in <- suppressWarnings(as.numeric(input$gt_complex_r))
  cr    <- if (length(cr_in) && is.finite(cr_in) && cr_in >= 0) cr_in else NULL
  # Inward displacement magnitude: blank -> auto (link_nm).
  dn_in <- suppressWarnings(as.numeric(input$gt_displace_nm))
  dn    <- if (length(dn_in) && is.finite(dn_in) && dn_in >= 0) dn_in else NULL
  tryCatch(
    estimate_gt_from_labeling(
      df, preset, eff_pct = eff, level = lvl, complex_r = cr,
      scope_cluster2  = isTRUE(input$gt_scope_cluster2 %||% TRUE),
      displace_inward = isTRUE(input$gt_displace_inward %||% TRUE),
      displace_nm     = dn,
      guided_merge    = isTRUE(input$gt_guided_merge),
      support_strength = { v <- suppressWarnings(as.numeric(input$gt_support_strength)); if (is.finite(v)) v else 1 },
      m_tid_override  = if (isTRUE(input$gt_use_calib)) {
        g <- gt_calib(); if (!is.null(g) && isTRUE(g$ok)) g$m_tid else NULL
      } else NULL
    ),
    error = function(e) { message("[GT-est] ", conditionMessage(e)); NULL }
  )
})

# Add the hollow-circle GT candidate overlay to a 2D figure (view-plane aware).
add_gt_overlay_2d <- function(fig, vp) {
  ge <- gt_estimate()
  if (is.null(ge) || !nrow(ge$points)) return(fig)
  pts <- ge$points
  remap_x <- function(x, y, z) switch(vp, XZ = x, YZ = y, x)
  remap_y <- function(x, y, z) switch(vp, XZ = z, YZ = z, y)
  px <- remap_x(pts$x, pts$y, pts$z); py <- remap_y(pts$x, pts$y, pts$z)
  lo  <- ge$meta$link_nm %||% 0
  hi  <- max(ge$meta$footprint_nm %||% 0, 2 * lo, stats::quantile(pts$pos_sd_nm, 0.9, na.rm = TRUE))
  col <- gt_uncertainty_colors(pts$pos_sd_nm, lo, hi)
  sz  <- suppressWarnings(as.numeric(input$gt_circle_size)); if (!is.finite(sz) || sz <= 0) sz <- 15
  struct_lbl <- ifelse(pts$cluster2 > 0, paste0("Structure ", pts$cluster2), "ungrouped")
  txt <- sprintf("<b>%s</b><br>GT candidate · ± %.0f nm%s<br>%d contributing unit(s) · %d TID / %d loc",
                 struct_lbl, pts$pos_sd_nm, ifelse(pts$fused, " (fused: multiple GTs)", ifelse(pts$count_fused %||% FALSE, " (more TIDs than one molecule)", "")),
                 pts$n_units, pts$n_tid, pts$n_loc)

  # Contribution fan-out lines: each member unit centroid -> its GT centroid.
  if (isTRUE(input$gt_show_lines %||% TRUE) && !is.null(ge$links) && nrow(ge$links)) {
    lk <- ge$links
    lx <- as.vector(rbind(remap_x(lk$ux, lk$uy, lk$uz),
                          remap_x(lk$gx, lk$gy, lk$gz), NA_real_))
    ly <- as.vector(rbind(remap_y(lk$ux, lk$uy, lk$uz),
                          remap_y(lk$gx, lk$gy, lk$gz), NA_real_))
    line_col <- if (isTRUE(input$dark_mode)) "rgba(255,255,255,0.35)" else "rgba(0,0,0,0.35)"
    fig <- fig %>% plotly::add_trace(
      x = lx, y = ly, type = "scatter", mode = "lines",
      line = list(color = line_col, width = 1), name = "GT contributions",
      showlegend = TRUE, inherit = FALSE, hoverinfo = "skip"
    )
  }
  fig %>% plotly::add_trace(
    x = px, y = py, type = "scatter", mode = "markers",
    marker = list(symbol = "circle-open", size = sz,
                  color = col, line = list(color = col, width = 3)),
    name = "Estimated GT", showlegend = TRUE, inherit = FALSE,
    text = txt, hoverinfo = "text"
  )
}

# Text/HTML readout of the GT estimate: counts, uncertainty floor, crowding.
output$gt_estimate_readout <- renderUI({
  if (!isTRUE(input$show_gt_estimate))
    return(helpText("Enable the estimate to see candidate GT positions and counts."))
  ge <- gt_estimate()
  if (is.null(ge) || !nrow(ge$summary))
    return(helpText("No estimable structures (need clustered data; try enabling de-blinking and Median clustering)."))
  s   <- ge$summary
  m   <- ge$meta
  tot_units  <- sum(s$n_units)
  tot_gt     <- sum(s$n_gt_est)
  tot_stoich <- sum(s$n_stoich)
  tot_lo     <- sum(s$n_gt_lo);  tot_hi <- sum(s$n_gt_hi)
  any_crowd  <- any(s$crowded)
  any_mm     <- any(s$stoich_mismatch)
  lvl_lbl    <- if (identical(m$level, "primary")) "de-blinked molecules" else "raw trace IDs"

  rows <- lapply(seq_len(nrow(s)), function(i) {
    r <- s[i, ]
    lab <- if (r$cluster2 == 0) "All (ungrouped)" else paste("Structure", r$cluster2)
    warn <- if (isTRUE(r$crowded))
      span(style = "color:#e67e22;font-weight:600;", " ⚠ crowded / under-resolved") else ""
    mm <- if (isTRUE(r$stoich_mismatch))
      span(style = "color:#c0392b;font-weight:600;",
           sprintf(" ⚠ %s: %.1f TIDs per GT, %.1f expected",
                   if (identical(r$count_verdict, "over-split")) "over-split" else "under-resolved",
                   r$tids_per_gt, m$tids_per_mol %||% m$m_tid)) else ""
    tags$li(HTML(sprintf(
      "<b>%s</b>: <b>%.0f</b> GT (spatial) vs <b>%.0f</b> (stoichiometry) · CI %.0f–%.0f · %d TID / %d unit%s%s%s%s",
      lab, r$n_gt_est, r$n_stoich, r$n_gt_lo, r$n_gt_hi, r$n_tid_total, r$n_units,
      ifelse(r$n_units == 1, "", "s"),
      ifelse(r$n_fused > 0, sprintf(", %d fused", r$n_fused), ""),
      as.character(warn), as.character(mm))))
  })

  tagList(
    div(style = "font-size:15px;margin-bottom:6px;",
        HTML(sprintf("Estimated total: <b>%.0f</b> GT molecules (spatial) · <b>%.0f</b> (stoichiometry) · range %.0f–%.0f",
                tot_gt, tot_stoich, tot_lo, tot_hi))),
    tags$ul(style = "margin:4px 0 8px 0;padding-left:18px;", rows),
    if (any_mm) div(style = "font-size:12px;color:#c0392b;margin-bottom:6px;",
        HTML("The trace total does not fit the spatial count in ≥1 structure. Under-resolved: candidates hold several molecules, so the spatial count is a lower bound. Over-split: re-blinks or dyes of one molecule were kept apart; raise the merge radius or check de-blinking. Either can also mean the preset or efficiency is off.")),
    div(style = "font-size:12px;color:#7f8c8d;line-height:1.45;",
        HTML(sprintf(
          "Modeled from %s. Units pooled within <b>%.0f nm</b> (label-complex radius)%s into one molecule%s. Position floor <b>±%.0f nm</b> (linker, irreducible) · <b>%.1f TIDs/molecule (%s)</b> · detection η=%.0f%% · σ<sub>loc</sub>≈%.1f nm.<br>Marker colour = position uncertainty, combining the linker floor, localization sampling (σ<sub>loc</sub>/√n<sub>loc</sub>), a <b>support penalty</b> (fewer TIDs → less certain), inward-displacement direction error, and a fusion penalty.<br><i>Counts are conditional on the labeling calibration; positions cannot be resolved below the linker floor. The two counts bracket the honest range.%s</i>",
          lvl_lbl, m$complex_r_nm %||% 0,
          if (isTRUE(m$guided_merge)) " (stoichiometry-guided)" else "",
          if ((m$displace_nm %||% 0) > 0) sprintf(", then displaced <b>%.0f nm</b> inward toward the structure", m$displace_nm) else "",
          m$link_nm, m$tids_per_mol %||% m$m_tid,
          if (identical(m$m_tid_source, "measured")) "measured from a standard" else "from preset",
          100 * m$eta, m$sigma_loc_nm,
          ifelse(any_crowd, " Crowded structures undercount - treat those counts as lower bounds.", ""))),
    gt_exploratory_note())
  )
})

# ── Phase B: per-structure shape / polymer morphology ───────────────────────
gt_morphology <- reactive({
  if (!isTRUE(input$gt_show_shape)) return(NULL)
  df <- plottable_data()
  if (is.null(df) || !nrow(df)) return(NULL)
  sp <- suppressWarnings(as.numeric(input$gt_subunit_spacing))
  if (!is.finite(sp) || sp <= 0) sp <- NA_real_
  # use the estimated GT candidates as a clean backbone when available
  bpts <- tryCatch({ ge <- gt_estimate(); if (!is.null(ge)) ge$points else NULL },
                   error = function(e) NULL)
  tryCatch(
    analyze_cluster2_morphology(df, backbone_pts = bpts, subunit_spacing_nm = sp),
    error = function(e) { message("[GT-morph] ", conditionMessage(e)); NULL }
  )
})

output$gt_morphology_readout <- renderUI({
  if (!isTRUE(input$gt_show_shape))
    return(helpText("Enable to classify each cluster's shape (filament / ring / punctum / sheet / diffuse) and report neutral size descriptors."))
  mo <- gt_morphology()
  if (is.null(mo) || !nrow(mo))
    return(helpText("No analysable structures (need Median Clusters with ≥4 traces)."))
  icon_for <- function(s) switch(s, filament = "――", ring = "○",
                                 punctum = "●", sheet = "▧", "∴")
  rows <- lapply(seq_len(nrow(mo)), function(i) {
    r <- mo[i, ]
    spc <- if (is.finite(r$spacing_nn_nm)) sprintf(", spacing %.0f nm", r$spacing_nn_nm) else ""
    detail <- if (r$shape == "filament")
      sprintf("contour <b>%.0f nm</b> (straight %.0f, width %.0f)%s, %.1f points/µm%s",
              r$contour_nm, r$extent_nm, r$width_nm, spc, r$lin_density_um,
              if (is.finite(r$n_subunit_est)) sprintf(" → <b>~%.0f subunits</b>", r$n_subunit_est) else "")
    else if (r$shape == "ring")
      sprintf("radius <b>%.0f nm</b> (circumference %.0f nm), R<sub>g</sub> %.0f%s", r$ring_r_nm, r$contour_nm, r$rg_nm, spc)
    else
      sprintf("extent <b>%.0f nm</b>, width %.0f, R<sub>g</sub> %.0f%s", r$extent_nm, r$width_nm, r$rg_nm, spc)
    tags$li(HTML(sprintf("<b>Structure %d</b> - <b>%s</b> %s: %s",
                         r$cluster2, r$shape, icon_for(r$shape), detail)))
  })
  tagList(
    tags$ul(style = "margin:4px 0 8px 0;padding-left:18px;", rows),
    div(style = "font-size:12px;color:#7f8c8d;line-height:1.45;",
        HTML("Shape is <b>inferred</b> from the gyration-tensor eigenvalues and a hollow-centre test for rings. Contour length follows the minimum spanning tree, so it tracks curves; spacing is the median nearest-neighbour distance."),
        gt_exploratory_note())
  )
})
