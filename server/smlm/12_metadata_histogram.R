# ------ HISTOGRAM for filters
output$metadata_histogram <- renderPlotly({
  req(input$show_metadata_hist, input$hist_variable, uploaded_df_rv())
  
  varname <- input$hist_variable
  
  get_vals <- function(df, var) {
    if (var == "tid") {
      tlens <- as.numeric(table(df$tid))
      tlens[is.finite(tlens)]
    } else {
      v <- df[[var]]
      v[is.finite(v)]
    }
  }
  
  df_full <- uploaded_df()
  
  df_keep <- filtered_df()
  
  vals_full <- get_vals(df_full, varname)
  vals_keep <- get_vals(df_keep, varname)
  
  if (!is.null(input$hist_range)) {
    lo <- input$hist_range[1]; hi <- input$hist_range[2]
    vals_full <- vals_full[vals_full >= lo & vals_full <= hi]
    vals_keep <- vals_keep[vals_keep >= lo & vals_keep <= hi]
  }
  
  if (length(vals_full) == 0 && length(vals_keep) == 0) return(NULL)
  
  bg_col   <- if (isTRUE(input$dark_mode)) "#1e1e1e" else "white"
  font_col <- if (isTRUE(input$dark_mode)) "white" else "black"
  
  fig <- plot_ly()

  # With the 2-color split on, DCR and pooled DCR are stacked by channel in the
  # channel colors. NULL (any other variable, split off or failed) falls through
  # to the plain histogram below.
  split_bins <- if (varname %in% c("dcr", "dcr_pooled"))
    channel_hist_layers(
      attach_channel_colors(hist_channels_df(), smlm_channel_colors(input)),
      varname, vals_keep,
      vals_raw = if (isTRUE(input$hist_show_raw)) vals_full else NULL,
      n_bins = input$hist_bins
    )

  if (!is.null(split_bins)) {
    # One bar series. Empty bins are left out to keep hover and payload clean;
    # marker line width 0 keeps narrow bins from turning into outlines.
    add_split_bars <- function(fig, y, base, name, color, opacity) {
      keep <- y > 0
      if (!any(keep)) return(fig)
      add_bars(
        fig, x = split_bins$mids[keep], y = y[keep], base = base[keep],
        width = split_bins$width, name = name, opacity = opacity,
        marker = list(color = color, line = list(width = 0)),
        hovertext = sprintf("%s: %d<br>%s %.3f to %.3f", name, as.integer(y[keep]), varname,
                            split_bins$breaks[-length(split_bins$breaks)][keep],
                            split_bins$breaks[-1L][keep]),
        hoverinfo = "text"
      )
    }
    if (!is.null(split_bins$raw))
      fig <- add_split_bars(fig, split_bins$raw, numeric(length(split_bins$raw)),
                            "Raw data", "lightgray", 0.35)
    for (lay in split_bins$layers)
      fig <- add_split_bars(fig, lay$count, lay$base, lay$name, lay$color, 0.85)
  } else {
    # Optional raw overlay (grey)
    if (isTRUE(input$hist_show_raw)) {
      fig <- fig %>%
        add_histogram(
          x = vals_full, nbinsx = input$hist_bins, name = "Raw data",
          opacity = 0.35, marker = list(color = "lightgray"), bingroup = "x"
        )
    }

    # Always show the "filtered" bars (blue)
    fig <- fig %>%
      add_histogram(
        x = vals_keep, nbinsx = input$hist_bins, name = "Filtered data",
        opacity = 0.85, marker = list(color = "steelblue"), bingroup = "x"
      )
  }

  fig %>%
    layout(
      barmode = "overlay",
      title = paste("Histogram of", if (varname == "tid") "Trace lengths" else varname),
      xaxis = list(title = if (varname == "tid") "Trace length (# per TID)" else varname),
      yaxis = list(title = "Count"),
      paper_bgcolor = bg_col, plot_bgcolor = bg_col, font = list(color = font_col),
      legend = list(orientation = "h", x = 0, y = 1.1)
    )
})
#2-color NN helpers










#--- 2D-Only server parts:

# --- 2D helpers --------------------------------------------------------------

# Build a dataURI for the uploaded background image and optionally rotate (0/90/180/270)



# ---- NN Table for Cluster2 (2D) ----
nn_table_cluster2 <- reactive({
  df <- plottable_data()
  req(isTRUE(input$enable_cluster2))
  
  df_geo <- df %>%
    dplyr::filter(!is.na(cluster2) & cluster2 != 0) %>%
    dplyr::group_by(cluster2) %>%
    dplyr::summarise(
      x = stats::median(x, na.rm = TRUE),
      y = stats::median(y, na.rm = TRUE),
      color = dplyr::first(color),    
      .groups = "drop"
    )
  if (nrow(df_geo) < 2) return(df_geo[0, ])
  
  coords <- as.matrix(df_geo[, c("x","y")])
  k_eff  <- min(input$nn_k_cluster2 + 1, nrow(coords))
  nn     <- RANN::nn2(coords, k = k_eff)
  
  rows <- list()
  for (i in 1:nrow(coords)) {
    for (j in 2:k_eff) {
      ni   <- nn$nn.idx[i, j]
      dist <- as.numeric(nn$nn.dists[i, j])
      if (dist >= input$distance_range_cluster2[1] && dist <= input$distance_range_cluster2[2]) {
        rows[[length(rows) + 1]] <- data.frame(
          source_cluster2 = df_geo$cluster2[i],
          target_cluster2 = df_geo$cluster2[ni],
          distance = round(dist, 2),
          source_x = df_geo$x[i], source_y = df_geo$y[i],
          target_x = df_geo$x[ni], target_y = df_geo$y[ni]
        )
      }
    }
  }
  if (length(rows)) dplyr::bind_rows(rows) else df_geo[0, ]
})


# ---- Reset filters to full ranges ----
observeEvent(uploaded_df(), {
  # When filters are locked, sd_ranges_ui already initialises the widgets from
  # sd_filt (via isolate(sd_filt) in the renderUI).  Calling updateNumericRangeInput
  # here would override those correct locked values with the new file's max range.
  if (isTRUE(isolate(input$lock_filters))) return()

  df0 <- uploaded_df()
  req(df0)

  sda <- sd_by_tid(df0, min_n = input$sd_min_locs)
  if (!is.data.frame(sda) || !nrow(sda)) return()

  max_or1 <- function(v) {
    v <- v[is.finite(v)]
    if (!length(v)) return(1)
    mx <- max(v, na.rm = TRUE)
    if (!is.finite(mx) || mx <= 0) 1 else mx
  }

  shinyWidgets::updateNumericRangeInput(session, "sd_x_range", value = c(0, max_or1(sda$sd_x)))
  shinyWidgets::updateNumericRangeInput(session, "sd_y_range", value = c(0, max_or1(sda$sd_y)))

  if ("sd_z" %in% names(sda) && any(is.finite(sda$sd_z))) {
    shinyWidgets::updateNumericRangeInput(session, "sd_z_range", value = c(0, max_or1(sda$sd_z)))
  } else {
    # if z doesn't exist, set a wide range so it never filters anything
    shinyWidgets::updateNumericRangeInput(session, "sd_z_range", value = c(0, Inf))
  }
}, ignoreInit = TRUE)


raw_df <- reactive({
  df <- uploaded_df()
  shiny::validate(shiny::need(!is.null(df) && nrow(df) > 0, "No data uploaded"))
  df
})

observeEvent(input$hist_variable, {
  df <- filtered_df()
  if (input$hist_variable == "tid") {
    tlens <- as.numeric(table(df$tid))
    if (!length(tlens)) return()
    rng <- range(tlens, na.rm = TRUE)
    if (all(is.finite(rng))) {
      shinyWidgets::updateNumericRangeInput(session, "hist_range", value = rng)
    }
  } else {
    vals <- df[[input$hist_variable]]
    if (is.null(vals) || !is.numeric(vals) || !length(vals)) return()
    rng <- range(vals, na.rm = TRUE)
    shinyWidgets::updateNumericRangeInput(session, "hist_range", value = rng)
  }
})
