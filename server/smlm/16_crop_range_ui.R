#--- 3D-Only server parts:


df_trim_ref <- reactive({
  df <- filtered_df()
  df
})

output$trim_box_ui <- renderUI({
  df <- df_trim_ref(); shiny::validate(shiny::need(nrow(df) > 0, "No data"))
  rx <- range(df$x, na.rm = TRUE); ry <- range(df$y, na.rm = TRUE); rz <- range(df$z, na.rm = TRUE)
  tagList(
    shinyWidgets::numericRangeInput("trim_x_range", "Keep X range", value = rx, min = rx[1], max = rx[2]),
    shinyWidgets::numericRangeInput("trim_y_range", "Keep Y range", value = ry, min = ry[1], max = ry[2]),
    shinyWidgets::numericRangeInput("trim_z_range", "Keep Z range", value = rz, min = rz[1], max = rz[2])
  )
})

output$trim_sphere_ui <- renderUI({
  df <- df_trim_ref(); shiny::validate(shiny::need(nrow(df) > 0, "No data"))
  c0 <- colMeans(df[,c("x","y","z")], na.rm = TRUE)
  r0 <- stats::quantile(sqrt(rowSums((as.matrix(df[,c("x","y","z")]) - 
                                        matrix(c0, nrow(df), 3, TRUE))^2)), probs = 0.95, na.rm = TRUE)
  tagList(
    numericInput("trim_cx", "Center X", value = round(c0[1], 2)),
    numericInput("trim_cy", "Center Y", value = round(c0[2], 2)),
    numericInput("trim_cz", "Center Z", value = round(c0[3], 2)),
    numericInput("trim_r",  "Radius",   value = round(r0, 2), min = 0)
  )
})

# When "Lock filters" is on, preserve the current widget value across file switches.
# The widget is still rebuilt (so bounds reflect the new file's data range), but
# the displayed value comes from the last-set filter rather than the new data's range.
.locked_val <- function(cur_input, rng) {
  if (isTRUE(input$lock_filters)) isolate(cur_input) %||% rng else rng
}
output$efo_range_ui <- renderUI({
  df  <- uploaded_df()
  rng <- range(df$efo, na.rm = TRUE)
  shinyWidgets::numericRangeInput("efo_range", "EFO Range",
    value = .locked_val(input$efo_range, rng), min = rng[1], max = rng[2])
})
output$cfr_range_ui <- renderUI({
  df  <- uploaded_df()
  rng <- range(df$cfr, na.rm = TRUE)
  shinyWidgets::numericRangeInput("cfr_range", "CFR Range",
    value = .locked_val(input$cfr_range, rng), min = rng[1], max = rng[2])
})
output$tim_range_ui <- renderUI({
  df  <- uploaded_df()
  rng <- range(df$tim, na.rm = TRUE)
  shinyWidgets::numericRangeInput("tim_range", "TIM Range",
    value = .locked_val(input$tim_range, rng), min = rng[1], max = rng[2])
})
output$dcr_range_ui <- renderUI({
  df  <- uploaded_df()
  rng <- range(df$dcr, na.rm = TRUE)
  shinyWidgets::numericRangeInput("dcr_range", "DCR Range",
    value = .locked_val(input$dcr_range, rng), min = rng[1], max = rng[2])
})
output$dwell_range_ui <- renderUI({
  df  <- uploaded_df()
  if (!"dwell" %in% names(df)) return(NULL)
  rng <- range(df$dwell, na.rm = TRUE)
  shinyWidgets::numericRangeInput("dwell_range", "DWELL Range",
    value = .locked_val(input$dwell_range, rng), min = rng[1], max = rng[2])
})
output$fluo_range_ui <- renderUI({
  df  <- uploaded_df()
  if (!"fluo" %in% names(df)) return(NULL)
  rng <- range(df$fluo, na.rm = TRUE)
  shinyWidgets::numericRangeInput("fluo_range", "FLUO Range",
    value = .locked_val(input$fluo_range, rng), min = rng[1], max = rng[2])
})
output$eco_range_ui <- renderUI({
  df  <- uploaded_df()
  rng <- range(df$eco, na.rm = TRUE)
  shinyWidgets::numericRangeInput("eco_range", "ECO Range",
    value = .locked_val(input$eco_range, rng), min = rng[1], max = rng[2])
})

data <- reactive({
  req(input$file)
  df <- filtered_df()
  df
})
cross_dcr_nn_table <- reactive({
  # Only compute when the toggle is on
  if (!isTRUE(input$nn_cross_dcr)) return(data.frame())
  
  df <- plottable_data()
  if (!nrow(df)) return(data.frame())
  
  # We need a DCR column to define groups 1/2
  if (!("dcr" %in% names(df))) return(data.frame())
  
  # DCR ranges must exist and be finite
  g1 <- input$dcr_group1; g2 <- input$dcr_group2
  if (is.null(g1) || is.null(g2) || length(g1) != 2 || length(g2) != 2 ||
      !all(is.finite(c(g1, g2)))) return(data.frame())
  
  # Per-trace medians (geometry + dcr); use dcr_pooled when requested
  .col <- if (isTRUE(input$use_pooled_dcr) && "dcr_pooled" %in% names(df)) "dcr_pooled" else "dcr"
  med <- df %>%
    dplyr::group_by(tid) %>%
    dplyr::summarise(
      x = stats::median(x, na.rm = TRUE),
      y = stats::median(y, na.rm = TRUE),
      z = stats::median(z, na.rm = TRUE),
      dcr_med = stats::median(.data[[.col]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::filter(is.finite(x), is.finite(y), is.finite(z), is.finite(dcr_med))
  
  if (!nrow(med)) return(data.frame())
  
  # Label groups strictly by DCR (this feature is explicitly “cross DCR”)
  med$grp <- dplyr::case_when(
    med$dcr_med >= g1[1] & med$dcr_med <  g1[2] ~ 1L,
    med$dcr_med >= g2[1] & med$dcr_med <= g2[2] ~ 2L,
    TRUE ~ NA_integer_
  )
  med <- med %>% dplyr::filter(grp %in% c(1L, 2L))
  if (!nrow(med) || length(unique(med$grp)) < 2) return(data.frame())
  
  # Split groups
  A <- med %>% dplyr::filter(grp == 1L)
  B <- med %>% dplyr::filter(grp == 2L)
  if (!nrow(A) || !nrow(B)) return(data.frame())
  
  # For each A, nearest neighbor in B (3D)
  coordsA <- as.matrix(A[, c("x","y","z")])
  coordsB <- as.matrix(B[, c("x","y","z")])
  nn <- tryCatch(RANN::nn2(data = coordsB, query = coordsA, k = 1),
                 error = function(e) NULL)
  if (is.null(nn)) return(data.frame())
  
  rng <- input$distance_range %||% c(0, Inf)
  
  out <- lapply(seq_len(nrow(A)), function(i) {
    j <- nn$nn.idx[i, 1]
    d <- nn$nn.dists[i, 1]
    if (!is.finite(d) || d < rng[1] || d > rng[2]) return(NULL)
    data.frame(
      from_x = A$x[i], from_y = A$y[i], from_z = A$z[i],
      to_x   = B$x[j], to_y   = B$y[j], to_z   = B$z[j],
      distance = d
    )
  })
  
  out <- do.call(rbind, out)
  if (is.null(out)) data.frame() else out
})

  
# Applies the COMMITTED 3D crop region (from trim_applied_rv, set by "Apply
# crop"), gated by the master "Crop to region (3D)" checkbox. The live box/
# sphere inputs drive only the preview, never the data, so tweaking coordinates
# never crops until the user applies.
apply_trimmer <- function(df, input) {
  if (!isTRUE(input$enable_trimmer) || nrow(df) == 0) return(df)
  ar <- trim_applied_rv()
  if (is.null(ar)) return(df)
  keep <- rep(TRUE, nrow(df))
  if (identical(ar$mode, "sphere")) {
    cc <- ar$c; r <- ar$r
    if (length(cc) == 3 && all(is.finite(c(cc, r))) && r >= 0)
      keep <- (df$x - cc[1])^2 + (df$y - cc[2])^2 + (df$z - cc[3])^2 <= r^2
  } else {
    if (!is.null(ar$x)) keep <- keep & df$x >= ar$x[1] & df$x <= ar$x[2]
    if (!is.null(ar$y)) keep <- keep & df$y >= ar$y[1] & df$y <= ar$y[2]
    if (!is.null(ar$z)) keep <- keep & df$z >= ar$z[1] & df$z <= ar$z[2]
  }
  df[keep, , drop = FALSE]
}

# Commit the current box/sphere inputs as the active crop region.
trim_region_from_inputs <- function() {
  mode <- input$trim_mode %||% "box"
  num2 <- function(v) { v <- suppressWarnings(as.numeric(unlist(v))); if (length(v) == 2 && all(is.finite(v))) v else NULL }
  if (identical(mode, "sphere")) {
    cc <- suppressWarnings(as.numeric(c(input$trim_cx, input$trim_cy, input$trim_cz)))
    r  <- suppressWarnings(as.numeric(input$trim_r))
    if (length(cc) == 3 && all(is.finite(cc)) && is.finite(r) && r >= 0)
      list(mode = "sphere", c = cc, r = r) else NULL
  } else {
    list(mode = "box", x = num2(input$trim_x_range),
         y = num2(input$trim_y_range), z = num2(input$trim_z_range))
  }
}
observeEvent(input$trim_apply, {
  reg <- trim_region_from_inputs()
  if (is.null(reg)) { showNotification("Set the crop coordinates first.", type = "warning"); return() }
  trim_applied_rv(reg); plot_trigger(Sys.time())
  showNotification(if (identical(reg$mode, "sphere"))
    sprintf("3D crop applied: sphere r %.0f nm at (%.0f, %.0f, %.0f).", reg$r, reg$c[1], reg$c[2], reg$c[3])
    else "3D crop applied (box).", type = "message", duration = 4)
}, ignoreInit = TRUE)
observeEvent(input$trim_clear, {
  trim_applied_rv(NULL); plot_trigger(Sys.time())
  showNotification("3D crop cleared.", type = "message", duration = 3)
}, ignoreInit = TRUE)

# Redraw the preview shape a short moment after the coordinates change (only
# while previewing in 3D), so the low-opacity box/sphere tracks the inputs
# without re-rendering on every keystroke.
trim_preview_tick <- shiny::debounce(reactive({
  list(input$trim_mode, input$trim_x_range, input$trim_y_range, input$trim_z_range,
       input$trim_cx, input$trim_cy, input$trim_cz, input$trim_r, input$trim_preview)
}), 400)
observeEvent(trim_preview_tick(), {
  if (isTRUE(input$enable_trimmer) && isTRUE(input$trim_preview) && identical(input$mode3d, "3D"))
    plot_trigger(Sys.time())
}, ignoreInit = TRUE)

# Add the low-opacity crop-region preview to a 3D figure from the LIVE inputs.
add_3d_crop_preview <- function(fig) {
  if (!isTRUE(input$enable_trimmer) || !isTRUE(input$trim_preview)) return(fig)
  reg <- trim_region_from_inputs(); if (is.null(reg)) return(fig)
  col <- "#2c7fb8"
  if (identical(reg$mode, "sphere")) {
    cc <- reg$c; r <- reg$r
    nu <- 20L; nv <- 14L
    th <- seq(0, 2 * pi, length.out = nu); ph <- seq(0, pi, length.out = nv)
    g  <- expand.grid(iu = seq_len(nu), iv = seq_len(nv))
    vx <- cc[1] + r * sin(ph[g$iv]) * cos(th[g$iu])
    vy <- cc[2] + r * sin(ph[g$iv]) * sin(th[g$iu])
    vz <- cc[3] + r * cos(ph[g$iv])
    idx <- function(iu, iv) (iv - 1L) * nu + iu
    ii <- jj <- kk <- integer(0)
    for (iv in seq_len(nv - 1L)) for (iu in seq_len(nu - 1L)) {
      a <- idx(iu, iv); b <- idx(iu + 1L, iv); c2 <- idx(iu, iv + 1L); d <- idx(iu + 1L, iv + 1L)
      ii <- c(ii, a - 1L, b - 1L); jj <- c(jj, b - 1L, d - 1L); kk <- c(kk, c2 - 1L, c2 - 1L)
    }
    return(fig %>% plotly::add_trace(type = "mesh3d", x = vx, y = vy, z = vz,
      i = ii, j = jj, k = kk, opacity = 0.15, color = I(col), facecolor = rep(col, length(ii)),
      flatshading = TRUE, hoverinfo = "skip", showlegend = FALSE, showscale = FALSE, name = "crop preview"))
  }
  x <- reg$x; y <- reg$y; z <- reg$z
  if (is.null(x) || is.null(y) || is.null(z)) return(fig)
  vx <- c(x[1], x[2], x[2], x[1], x[1], x[2], x[2], x[1])
  vy <- c(y[1], y[1], y[2], y[2], y[1], y[1], y[2], y[2])
  vz <- c(z[1], z[1], z[1], z[1], z[2], z[2], z[2], z[2])
  ii <- c(0,0,4,4,0,0,3,3,0,0,1,1); jj <- c(1,2,5,6,1,5,2,6,3,7,2,6); kk <- c(2,3,6,7,5,4,6,7,7,4,6,5)
  fig %>% plotly::add_trace(type = "mesh3d", x = vx, y = vy, z = vz,
    i = ii, j = jj, k = kk, opacity = 0.15, color = I(col), facecolor = rep(col, length(ii)),
    flatshading = TRUE, hoverinfo = "skip", showlegend = FALSE, showscale = FALSE, name = "crop preview")
}
