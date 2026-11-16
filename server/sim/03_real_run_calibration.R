## ─────────────────────────────────────────────────────────────────────────────
## server/sim/03_real_run_calibration.R  -  calibrate the acquisition model from a real run
## ─────────────────────────────────────────────────────────────────────────────
## Part of sim_server(). Not a standalone script: its statements are evaluated in
## the server function's environment, in file-name order (server/parts.R).
##
## Reading a real file, the quality gate, and pushing the measured rates into
## the acquisition controls.
##
## INTERFACE (generated: Rscript tools/server_map.R --headers; do not edit by hand)
##   DEFINES for other parts:
##     sim_real_rv
##   READS from other parts:
##     (nothing)
## END INTERFACE
## ─────────────────────────────────────────────────────────────────────────────

# ── Calibrate from a real MINFLUX run ────────────────────────────────────────
# Load a .npy / .json / .csv run, measure its detection-rate profile, duration
# and FOV, pre-fill the Acquisition inputs, and keep the measured rate curve to
# overlay on the Acquisition plot.
sim_real_rv  <- reactiveVal(NULL)   # analysis of the GATED run (drives everything)
sim_real_raw <- reactiveVal(NULL)   # the run exactly as uploaded, never gated


# The gated run.  Rebuilt whenever the file, a gate range or the axial factor
# changes; every calibrated number downstream is measured on THIS subset.
# The stored z is RAW, so the factor is applied here — one place, so sigma_z
# and fov_z can never disagree with the number shown in the input.
sim_real_gated <- reactive({
  raw <- sim_real_raw()
  if (is.null(raw)) return(NULL)
  g <- gate_real_run(raw$df,
                     efo_rng = input$sim_real_efo_range,
                     cfr_rng = input$sim_real_cfr_range,
                     len_rng = input$sim_real_len_range)
  if (is.null(g)) return(NULL)
  # Until the widget exists, fall back to the factor detected for this file,
  # so the very first calibration already uses the right one.
  zc <- suppressWarnings(as.numeric(input$sim_real_z_cal %||% raw$cal$factor))
  if (!is.finite(zc) || zc <= 0) zc <- 1
  if ("z" %in% names(g$df)) g$df$z <- g$df$z * zc
  g$z_cal <- zc
  g$z_src <- raw$cal
  g$has_z <- isTRUE(raw$has_z)
  g
})

# Axial calibration for the loaded run.  Rebuilt per upload (so a new file can
# never inherit the previous file's factor) and pre-filled from what the file
# itself says: .npy/.json are raw and need the correction; a CSV is returned
# to raw first, then offered the factor its producer had applied.
# NOTE: this is the refractive-index axial correction, NOT sim_crlb_z_factor
# (the sigma_z/sigma_xy precision ratio). Both default near 0.7 and they are
# unrelated quantities.
output$sim_real_zcal_ui <- renderUI({
  raw <- sim_real_raw()
  if (is.null(raw)) return(NULL)
  cal <- raw$cal
  # A 2D run has no axial structure to calibrate — say nothing about z.
  if (!isTRUE(raw$has_z))
    return(tags$div(style = "font-size:0.78em; color:#999; margin-bottom:6px;",
                    "2D run — no axial calibration needed."))
  tagList(
    numericInput("sim_real_z_cal", "Axial calibration factor (z x)",
                 value = cal$factor, min = 0.05, max = 5, step = 0.05),
    tags$div(style = "font-size:0.78em; color:#666; margin-top:-6px;",
      sprintf("Detected: %s%s. Applied to z before σ_z and the FOV depth are measured.",
              cal$source, if (isTRUE(cal$certain)) "" else " - assumed"))
  )
})

# Gate widgets, built from the uploaded run's own ranges so the default value
# spans the whole data and is therefore a no-op.  Columns the file does not
# carry (a plain CSV often has no efo / cfr) say so instead of silently
# passing everything.
output$sim_real_gate_ui <- renderUI({
  raw <- sim_real_raw()
  if (is.null(raw)) return(NULL)
  df     <- raw$df
  tl_rng <- range(as.integer(table(df$tid)))

  # No `step` here on purpose: numericRangeInput errors on step = NULL, and
  # EFO/CFR have no natural increment - the browser's default suits both.
  num_rng <- function(id, label, nm) {
    if (!.col_ok(df, nm))
      return(tags$div(style = "font-size:0.78em; color:#999; margin-bottom:6px;",
                      sprintf("%s - not in this file, gate unavailable.",
                              toupper(nm))))
    v <- suppressWarnings(as.numeric(df[[nm]]))
    r <- range(v[is.finite(v)])
    shinyWidgets::numericRangeInput(id, label, value = r,
                                    min = r[1], max = r[2])
  }

  tags$details(
    class = "sim-group", open = NA,
    tags$summary(tags$b("Calibration gate")),
    helpText("Restricts which localizations and traces the precision and",
             "trace-length statistics are measured from. EFO and CFR gate",
             "individual localizations; the length gate then drops whole",
             "traces by how many localizations survived."),
    num_rng("sim_real_efo_range", "EFO range (Hz)",  "efo"),
    num_rng("sim_real_cfr_range", "CFR range",       "cfr"),
    shinyWidgets::numericRangeInput(
      "sim_real_len_range", "Localizations per trace",
      value = tl_rng, min = tl_rng[1], max = tl_rng[2], step = 1),
    uiOutput("sim_real_gate_info")
  )
})

# Before / after counts for both gate stages.  The gate is only useful if its
# severity is visible, so this always reports rather than only on trouble.
output$sim_real_gate_info <- renderUI({
  g <- sim_real_gated()
  if (is.null(g)) return(NULL)
  pct <- function(a, b) if (b > 0) sprintf("%.0f%%", 100 * a / b) else "-"

  too_few <- g$n_loc_out < 10L || g$n_tr_sigma < 3L
  tags$div(
    class = if (too_few) "alert alert-warning" else "well well-sm",
    style = "padding:6px 10px; font-size:0.78em; margin:6px 0;",
    tags$b("Gate"), tags$br(),
    sprintf("Localizations: %d → %d  (%s kept)",
            g$n_loc_in, g$n_loc_out, pct(g$n_loc_out, g$n_loc_in)), tags$br(),
    sprintf("Traces: %d → %d  (%s kept); %d with ≥ 3 locs set σ",
            g$n_tr_in, g$n_tr_out, pct(g$n_tr_out, g$n_tr_in), g$n_tr_sigma),
    if (any(g$n_na > 0)) tagList(tags$br(),
      sprintf("%s value(s) missing and left ungated.",
              paste(sprintf("%s: %d", toupper(names(g$n_na)[g$n_na > 0]),
                            g$n_na[g$n_na > 0]), collapse = ", "))),
    if (too_few) tagList(tags$br(), tags$b(
      "Too little left to calibrate - widen the gate.")),
    # Gating is not neutral: every one of these selections keeps the brighter,
    # better-centred, longer-lived molecules, which are also the most precise.
    if (length(g$narrowed) > 0 && !too_few) tagList(tags$br(), tags$span(
      style = "color:#854d0e;",
      sprintf("%s narrowed - this keeps the brighter / better-centred / longer traces, so the calibrated σ is tighter than the run as a whole.",
              paste(g$narrowed, collapse = ", "))))
  )
})

# Split the measured total rate into structure + constant background given the
# fraction of traces that are on-structure, and push the fitted MODEL plus its
# struct_rate / bg_rate / level. Background is the constant offset; the
# structure is the time-varying part on top, whose shape (decay / constant /
# rise) analyze_real_run already chose. struct_rate is set to the structure's
# CHARACTERISTIC rate (peak for decay, plateau for rise), and level to the
# matching [0,1] ratio, so the simulated curve reproduces the measured profile.
apply_real_rates <- function(a, frac) {
  frac  <- max(0, min(1, frac))
  model <- a$det_model %||% "decay"
  # Background: constant offset, capped below the run minimum so the structure
  # component stays positive at both ends.
  bg <- max(0, min((1 - frac) * a$rate_mean, 0.9 * min(a$rate_peak, a$rate_floor)))
  s_early <- max(0.1, a$rate_peak  - bg)      # structure rate near the start
  s_late  <- max(0,   a$rate_floor - bg)      # structure rate near the end
  if (model == "constant") {
    struct_rate <- max(0.1, a$rate_mean - bg); level <- 1
  } else if (model == "rise") {               # plateau = late, starts low
    struct_rate <- max(0.1, s_late);  level <- max(0, min(1, s_early / struct_rate))
  } else {                                    # decay: peak = early, → floor
    struct_rate <- max(0.1, s_early); level <- max(0, min(1, s_late / struct_rate))
  }
  updateRadioButtons(session, "sim_acq_det_model", selected = model)
  updateNumericInput(session, "sim_struct_rate", value = round(struct_rate))
  updateNumericInput(session, "sim_bg_rate",     value = round(bg))
  updateSliderInput (session, "sim_acq_floor_frac", value = round(level, 2))
}

observeEvent(input$sim_real_file, {
  f <- input$sim_real_file
  req(f)
  ext <- get_ext(f$name)
  # Read to RAW, uncorrected z from every source, exactly like the SMLM tab —
  # then ONE axial factor is applied in sim_real_gated(), under user control.
  # Reading with the readers' 0.7 default instead would bake the correction in
  # for .npy/.json while leaving a .csv at whatever its producer applied, so
  # the same run would give a different sigma_z depending on its format.
  got <- tryCatch(
    withProgress(message = "Reading run…", value = 0.2, {
      switch(ext,
        npy  = list(df = read_minflux_npy_as_df(
                           f$datapath, keep_final_only = TRUE,
                           keep_valid_only = TRUE, z_scale = 1.0),
                    cal = list(factor = PYMINFLUX_Z_SCALE_DEFAULT,
                               source = "raw instrument data (.npy)",
                               certain = TRUE)),
        json = list(df = read_minflux_json_as_df(
                           f$datapath, keep_final_only = TRUE,
                           keep_valid_only = TRUE, z_scale = 1.0),
                    cal = list(factor = PYMINFLUX_Z_SCALE_DEFAULT,
                               source = "raw instrument data (.json)",
                               certain = TRUE)),
        # A CSV arrives with SOME factor already in its z. strip_csv_z_
        # calibration() divides it back out and reports what it found —
        # our own export records the factor explicitly, pyMINFLUX does not.
        csv  = {
          raw <- strip_csv_z_calibration(read.csv(f$datapath, check.names = FALSE))
          list(df = raw, cal = attr(raw, "z_cal"))
        },
        stop("Unsupported file type (use .npy, .json or .csv)."))
    }),
    error = function(e) {
      showNotification(paste("Load failed:", conditionMessage(e)),
                       type = "error", duration = 8); NULL
    })
  if (is.null(got)) return()
  # Same column canonicalization the SMLM tab applies, so alternative headers
  # (trace_id, x_nm, photons, …) are accepted here too.
  df <- canonicalize_minflux_df(got$df)
  if (!all(c("tid", "tim", "x", "y") %in% names(df))) {
    showNotification("Could not read run (need tid, tim, x, y).",
                     type = "error", duration = 6); return()
  }
  # Store the run ungated and at RAW z, together with the detected axial
  # provenance.  Building the widgets from it makes sim_real_gated() fire, and
  # the calibration observer does the rest — so upload, a gate change and an
  # axial-factor change all follow exactly the same path.
  sim_real_raw(list(df = df, cal = got$cal, ext = ext,
                    has_z = has_axial_data(df)))
  showNotification("Run loaded - calibrating from the current gate.",
                   type = "message", duration = 4)
})

# Push every measured value into the Acquisition / precision inputs.  Called on
# upload AND on every gate change, so the inputs always reflect the subset the
# gate currently selects.
push_real_calibration <- function(a) {
  # Pre-fill the Acquisition inputs from the measured values.
  hrs  <- floor(a$duration_hr)
  mins <- min(59L, as.integer(round((a$duration_hr - hrs) * 60)))
  updateRadioButtons(session, "sim_mode", selected = "acquisition")
  updateNumericInput(session, "sim_acq_hours",     value = hrs)
  updateNumericInput(session, "sim_acq_minutes",   value = mins)
  updateNumericInput(session, "sim_acq_decay_min", value = round(a$tau_min, 1))
  updateNumericInput(session, "sim_fov_x", value = ceiling(a$fov_x * 1.1 / 100) * 100)
  updateNumericInput(session, "sim_fov_y", value = ceiling(a$fov_y * 1.1 / 100) * 100)
  if (a$fov_z > 1)
    updateNumericInput(session, "sim_fov_z",
                       value = ceiling(a$fov_z * 1.1 / 100) * 100)
  # Structure/background split from the spatial classification (overridable).
  updateSliderInput(session, "sim_real_struct_frac",
                    value = round(a$struct_frac, 2))
  apply_real_rates(a, a$struct_frac)

  # Measured localization precision → Manual σ; photons, k and σ₀ → CRLB.
  if (is.finite(a$sigma_x) && is.finite(a$sigma_y))
    updateNumericInput(session, "sim_precision_xy",
                       value = round(mean(c(a$sigma_x, a$sigma_y)), 1))
  if (is.finite(a$sigma_z))
    updateNumericInput(session, "sim_precision_z", value = round(a$sigma_z, 1))
  if (is.finite(a$photons))
    updateNumericInput(session, "sim_crlb_photons", value = round(a$photons))
  # L stays the user's instrument setting. k is expressed relative to the
  # ideal bound at the CURRENT L and background, so it reproduces the data.
  # bound1 = σ_CRLB at N = 1, i.e. L·g/√8.
  L_now  <- suppressWarnings(as.numeric(isolate(input$sim_crlb_L)))
  bg_now <- suppressWarnings(as.numeric(isolate(input$sim_crlb_bg)))
  if (!is.finite(L_now))  L_now  <- 40
  if (!is.finite(bg_now)) bg_now <- 0.05
  bound1 <- smlm_sigma_loc(L_now, 1, bg_now)$crlb
  pf     <- a$prec_fit
  if (!is.null(pf) && isTRUE(pf$ok)) {
    updateNumericInput(session, "sim_crlb_k",     value = round(pf$c_phot / bound1, 2))
    updateNumericInput(session, "sim_crlb_floor", value = round(pf$floor, 1))
  } else if (is.finite(a$photons) && a$photons > 0 &&
             is.finite(a$sigma_x) && is.finite(a$sigma_y)) {
    # Photon counts too uniform to separate the terms: assume the ideal
    # bound (k = 1) and give the rest of the measured variance to the floor.
    sig_med <- mean(c(a$sigma_x, a$sigma_y))
    crlb_m  <- bound1 / sqrt(a$photons)
    updateNumericInput(session, "sim_crlb_k",     value = 1)
    updateNumericInput(session, "sim_crlb_floor",
                       value = round(sqrt(max(0, sig_med^2 - crlb_m^2)), 1))
  }
}

# Re-calibrate whenever the run or the gate changes.  Everything downstream
# reads sim_real_rv(), so the whole tab follows the gate from this one place.
# A gate that leaves too little to measure keeps the previous calibration
# rather than silently resetting the inputs - the gate info box says so.
observeEvent(sim_real_gated(), {
  g <- sim_real_gated()
  req(g)
  a <- analyze_real_run(g$df)
  if (is.null(a)) return()
  sim_real_rv(a)
  push_real_calibration(a)
}, ignoreNULL = TRUE)

# Manual override of the structure/background split (re-derives the rates).
observeEvent(input$sim_real_struct_frac, {
  a <- sim_real_rv()
  if (!is.null(a)) apply_real_rates(a, input$sim_real_struct_frac)
}, ignoreInit = TRUE)

output$sim_real_status <- renderUI({
  a <- sim_real_rv()
  if (is.null(a))
    return(tags$small(style = "color:#888;",
      "Load a .npy / .json MINFLUX run to calibrate the model from real data."))
  # Counts here are those of the GATED subset, not the file - the gate box
  # directly above reports the before/after so the two can never be confused.
  # When the gate is starved, these numbers are the last calibration that
  # succeeded, so they must be labelled as stale rather than as current.
  g     <- sim_real_gated()
  stale <- !is.null(g) && (g$n_loc_out < 10L || g$n_tr_sigma < 3L)
  suffix <- if (stale) " (last valid gate)"
            else if (!is.null(g) && length(g$narrowed) > 0) " (gated)"
            else ""
  tags$div(class = if (stale) "alert alert-warning" else "alert alert-success",
           style = "padding:6px 10px; font-size:0.82em; margin:4px 0;",
    tags$b(sprintf("%d traces over %.1f h%s", a$n_traces, a$duration_hr, suffix)),
    tags$br(),
    sprintf("Total rate: peak %.0f, mean %.0f, floor %.0f /hr  (τ≈%.0f min)",
            a$rate_peak, a$rate_mean, a$rate_floor, a$tau_min),
    tags$br(),
    sprintf("Spatial split: %.0f%% on-structure / %.0f%% background",
            100 * a$struct_frac, 100 * (1 - a$struct_frac)),
    tags$br(),
    sprintf("FOV: %.0f × %.0f%s nm", a$fov_x, a$fov_y,
            if (a$fov_z > 1) sprintf(" × %.0f", a$fov_z) else ""),
    # ── Point-cloud profile ──────────────────────────────────────────────
    # Median ± nMAD per axis: the SPREAD is what makes the cloud MINFLUX-like,
    # so it is shown next to the median rather than collapsed into it.
    if (is.finite(a$sigma_x)) tagList(
      tags$hr(style = "margin:5px 0;"),
      tags$b("Point cloud"),
      sprintf(" (%d traces, ≥3 locs)", a$n_sigma_traces), tags$br(),
      sprintf("σ ± nMAD: x %.1f ± %.1f, y %.1f ± %.1f%s nm",
              a$sigma_x, a$sigma_x_mad, a$sigma_y, a$sigma_y_mad,
              if (is.finite(a$sigma_z))
                sprintf(", z %.1f ± %.1f", a$sigma_z, a$sigma_z_mad) else ""),
      tags$br(),
      sprintf("Axis ratio: y/x %.2f%s", a$ratio_yx,
              if (is.finite(a$ratio_zx)) sprintf(", z/x %.2f", a$ratio_zx) else ""),
      tags$br(),
      sprintf("Locs per trace: median %.0f (IQR %.0f–%.0f)",
              a$len_med, a$len_q[1], a$len_q[2]),
      if (is.finite(a$photons)) tagList(tags$br(),
        sprintf("%.0f photons/loc (n=%d)", a$photons, a$eco_n)),
      # Precision model fit: photon term at the median N, and the floor.
      {
        pf <- a$prec_fit
        if (!is.null(pf) && isTRUE(pf$ok)) tagList(tags$br(), sprintf(
          "σ² vs 1/N over %d photon bins: photon term %.1f nm at median N, floor σ₀ %.1f%s nm",
          pf$n_bins, pf$c_phot / sqrt(a$photons), pf$floor,
          if (is.finite(pf$se_floor)) sprintf(" ± %.1f", pf$se_floor) else ""))
        else if (is.finite(a$photons)) tagList(tags$br(), tags$span(
          style = "color:#8a6d3b;",
          "Photon counts too uniform to fit k and σ₀ separately; k set to 1 and the rest assigned to σ₀."))
      },
      tags$br(),
      tags$span(style = "color:#666;", sprintf(
        "σ is small-n debiased (×%.3f at the median trace length). Within-trace SD also contains real motion and drift, so it bounds precision rather than measuring it.",
        1 / a$debias_med))))
})
