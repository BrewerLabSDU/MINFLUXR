# ---------------------------------------------------------------------------
# ---- Change-point segmentation + GMM clustering ---------------------------
# Splits each track into diffusion-state segments via PELT on log(r²/Δt),
# then fits a GMM (mclust, BIC-selected G) to the per-segment D distribution.
# ---------------------------------------------------------------------------
# ==========================================================================
# cpd_segments_rv() + cpd_gmm_rv() - CHANGE-POINT SEGMENTATION + STATE COUNT
# ==========================================================================
#
# GOAL
#   Find the MOMENTS a molecule changes its mobility, without assuming in
#   advance how many states exist. This is the model-free counterpart to the
#   HMM: the HMM fits K states globally, whereas this first detects abrupt
#   changes within each track and only then asks how many groups the
#   resulting segments fall into.
#
# THE MATH
#   (a) THE SIGNAL. From <r²> = 2dDΔt (header §2b), the per-step quantity
#                        r² / Δt   ∝   local D
#       is an (extremely noisy) instantaneous diffusivity. Dividing by Δt is
#       what makes irregular MINFLUX sampling comparable - this is the
#       SQUARED-displacement form of the sqrt(Δt) rule.
#
#       We segment log(r²/Δt), not r²/Δt, because:
#         - single-step D is roughly exponential and strongly right-skewed;
#           the log makes it ~Gaussian, which a mean-shift test assumes;
#         - a change in D is MULTIPLICATIVE (e.g. "10x slower"), and in log
#           space that becomes an ADDITIVE mean shift - what PELT detects.
#
#   (b) PELT (Killick, Fearnhead & Eckley 2012) finds the change-point set
#       minimising  Σ_segments cost(segment) + β·(n change-points),
#       exactly and in O(n). The BIC penalty β stops it inventing a change for
#       every noisy step; `minseglen` sets the shortest allowed segment.
#
#   (c) PER-SEGMENT D uses the mean of ratios, NOT the ratio of means:
#            D = mean(r²/Δt) / (2d)          (correct)
#        not D = mean(r²) / (2d·mean(Δt))    (biased when r² and Δt co-vary -
#                                             Jensen's inequality)
#       Then subtract the noise floor at the D level:
#            D_true = max(0, D − σ²·mean(1/Δt))
#       so an immobile segment reads ~0 rather than a small positive D.
#
#   (d) HOW MANY STATES: pool the per-segment D values and fit a Gaussian
#       mixture in log10(D), component count chosen by BIC (`cpd_gmm_rv`).
#       Exploratory: segments of different lengths from one population can
#       return two components (est_cpd.R header).
#
# OBJECTS IT CALLS (inputs)
#   trk_spatial_filtered() (last pipeline stage), has_z(), trk_sigma_loc_eff()
#   inputs: trk_cpd_compute (Run trigger), trk_cpd_min_seg, trk_cpd_max_g,
#           trk_cpd_sub_noise
#
# STEP BY STEP
#   1. Per track build r², Δt, and log(r²/Δt).
#   2. Run PELT -> change-point indices -> segment boundaries.
#   3. Per segment: D via mean-of-ratios, minus the noise floor.
#   4. Pool all segments; fit the BIC-selected mixture in log10 D.
#
# RETURNS
#   cpd_segments_rv : data.frame(tid, seg_idx, n_states, n_steps,
#                                D_nm2_ms, D_um2_s) - one row per segment
#   cpd_gmm_rv      : fitted mclust object; its components ARE the inferred
#                     diffusive states
#
# CAVEAT: single-step D is very noisy, so PELT can over-segment (`min_seg` is
# the guard). Cross-check the state count against the HMM - same question,
# different route, they should broadly agree.
# ==========================================================================
CPD_PAL <- c("#3060C0","#C03030","#2CA02C","#FF7F0E","#9467BD","#8C564B")

cpd_segments_rv <- eventReactive(input$trk_cpd_compute, {
  req(isTRUE(input$trk_show_cpd))
  df <- tryCatch(trk_spatial_filtered(), error = function(e) NULL)
  req(!is.null(df), nrow(df) > 0,
      all(c("x","y","tid","tim") %in% names(df)))

  min_seg    <- max(3L, as.integer(input$trk_cpd_min_seg %||% 8L))
  use_z_cpd  <- isTRUE(has_z())
  ndim       <- if (use_z_cpd) 3L else 2L
  coef_D     <- 2L * ndim   # 4 (2D) or 6 (3D): MSD = coef_D * D * dt

  # Optional σ_loc noise-floor subtraction on the per-segment D. A pure-noise step
  # has E[|Δr|²] = coef_D·σ² (each of the two localizations contributes σ² per axis),
  # so the apparent D of an immobile segment is σ²·⟨1/Δt⟩. Subtract it at the D level
  # (not per step - that would rectify and bias upward), clamped ≥ 0. Same floor the
  # D-map uses, so immobile segments no longer read as a small non-zero D. No-op when
  # σ_loc is unset or the toggle is off.
  sig_cpd    <- tryCatch(trk_sigma_loc_eff(), error = function(e) NA_real_)
  sub_noise_cpd <- isTRUE(input$trk_cpd_sub_noise) && is.finite(sig_cpd) && sig_cpd > 0
  noise_D_off   <- if (sub_noise_cpd) sig_cpd^2 else 0   # × ⟨1/Δt⟩ per segment below

  t_c_cpd <- tryCatch(trk_cycle_time(), error = function(e) NA_real_)
  # Motion blur needs the cycle time. Without it the fit would silently run with no
  # blur correction, which reads D about a third low at MINFLUX cycle times, so the
  # panel asks for the calibration instead of guessing.
  shiny::validate(shiny::need(is.finite(t_c_cpd) && t_c_cpd > 0,
    "Cycle time unknown - run Data Quality & Calibration first (the blur correction needs it)."))
  # Segmentation core lives in est_cpd.R (pure, tested).
  result <- withProgress(message = "Segmenting tracks…", value = 0,
    cpd_segment_tracks(df, min_seg = min_seg, use_z = use_z_cpd,
                       noise_sig2 = noise_D_off,
                       R_abs = trk_blur()$R_abs,
                       progress = function(f) incProgress(f)))
  shiny::validate(shiny::need(
    !is.null(result) && nrow(result) >= 3L,
    "No segments found - try lowering Min segment length."
  ))
  result
}, ignoreNULL = TRUE)

# ── cpd_gmm_rv() - how many diffusive states are there? ─────────────────────
# GOAL: take the per-segment diffusion coefficients produced by the
# change-point step and ask whether they fall into distinct populations.
# A Gaussian mixture is fitted and the number of components chosen by BIC, so
# the state count is selected by the data rather than asserted by the user.
#
# Fitted in log10(D) space: diffusion coefficients span orders of magnitude and
# are right-skewed, so log space is where the populations look Gaussian.
cpd_gmm_rv <- reactive({
  req(isTRUE(input$trk_show_cpd))
  segs <- tryCatch(cpd_segments_rv(), error = function(e) NULL)   # per-segment D values
  req(!is.null(segs), nrow(segs) >= 3L)
  cpd_gmm(segs, max_G = input$trk_cpd_max_g %||% 4L)              # est_cpd.R
})

output$trk_cpd_hist <- renderPlotly({
  req(isTRUE(input$trk_show_cpd))
  segs <- tryCatch(cpd_segments_rv(), error = function(e) {
    shiny::validate(shiny::need(FALSE, conditionMessage(e))); NULL
  })
  req(!is.null(segs))
  gmm <- tryCatch(cpd_gmm_rv(), error = function(e) NULL)

  dark   <- isTRUE(isolate(input$trk_dark_mode))
  bg_col <- if (dark) "#000000" else "#FFFFFF"
  fg_col <- if (dark) "#FFFFFF" else "#333333"

  log10_D_fin <- log10(segs$D_um2_s)
  log10_D_fin <- log10_D_fin[is.finite(log10_D_fin)]
  n_bins <- 40L

  # add_histogram, not inline plot_ly(type="histogram", nbinsx=): the GMM overlay
  # scatter lines added below would otherwise inherit nbinsx from the plot_ly level
  # and warn "'scatter' objects don't have these attributes: 'nbinsx'".
  fig <- plot_ly() %>% add_histogram(
    x    = log10_D_fin,
    nbinsx = n_bins,
    marker = list(color = "rgba(100,120,210,0.55)",
                  line  = list(color = "rgba(80,100,190,0.8)", width = 0.5)),
    name = "Segments"
  )

  if (!is.null(gmm) && gmm$G >= 1L) {
    x_seq  <- seq(min(log10_D_fin) - 0.5, max(log10_D_fin) + 0.5, length.out = 400)
    props  <- gmm$parameters$pro
    means  <- gmm$parameters$mean
    sds    <- sqrt(gmm$parameters$variance$sigmasq)
    bw     <- diff(range(log10_D_fin)) / n_bins
    scale  <- length(log10_D_fin) * bw     # histogram → density scale factor

    total_y <- numeric(length(x_seq))
    for (g in seq_len(gmm$G)) {
      y_g      <- props[g] * dnorm(x_seq, means[g], sds[g]) * scale
      total_y  <- total_y + y_g
      lbl      <- sprintf("Pop. %d  D =%.3g µm²/s  (%.1f%%)",
                          g, 10^means[g], props[g] * 100)
      fig <- fig %>% add_trace(
        x = x_seq, y = y_g, type = "scatter", mode = "lines",
        line = list(color = CPD_PAL[g], width = 2.5),
        name = lbl
      )
    }
    fig <- fig %>% add_trace(
      x = x_seq, y = total_y, type = "scatter", mode = "lines",
      line = list(color = "rgba(200,30,30,0.75)", width = 1.5, dash = "dash"),
      name = "GMM total"
    )
  }

  fig %>% layout(
    barmode       = "overlay",
    paper_bgcolor = bg_col, plot_bgcolor = bg_col,
    font          = list(color = fg_col),
    xaxis = list(title = "log₁₀(D)  (µm²/s)", color = fg_col,
                 gridcolor = if (dark) "rgba(255,255,255,0.15)" else "rgba(0,0,0,0.12)"),
    yaxis = list(title = "Segments", color = fg_col,
                 gridcolor = if (dark) "rgba(255,255,255,0.15)" else "rgba(0,0,0,0.12)"),
    legend = list(orientation = "h", x = 0, y = 1.12, font = list(size = 11))
  )
})

output$trk_cpd_summary <- renderUI({
  req(isTRUE(input$trk_show_cpd))
  segs <- tryCatch(cpd_segments_rv(), error = function(e) NULL)
  req(!is.null(segs))
  gmm  <- tryCatch(cpd_gmm_rv(),      error = function(e) NULL)

  n_tracks   <- length(unique(segs$tid))
  n_segs     <- nrow(segs)
  n_multi    <- sum(tapply(segs$n_states, segs$tid, max) > 1L)

  gmm_line <- if (!is.null(gmm)) {
    G <- gmm$G
    paste0("GMM selected ", G, " component",
           if (G > 1L) "s" else "", " by BIC:  ",
           paste0(sapply(seq_len(G), function(g)
             sprintf("Pop.%d D =%.3g µm²/s (%.1f%%)",
                     g, 10^gmm$parameters$mean[g],
                     gmm$parameters$pro[g] * 100)),
             collapse = "  |  "))
  } else "GMM not yet computed."

  tags$div(class = "well well-sm", style = "font-size:0.84em; margin-top:6px;",
    tags$b("Segmentation:"), " ",
    n_tracks, " tracks  •  ",
    n_segs, " segments  •  ",
    n_multi, " tracks with ≥2 states detected", tags$br(),
    tags$b("GMM: "), gmm_line
  )
})

output$trk_dl_cpd_csv <- downloadHandler(
  filename = function() paste0("cpd_segments_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
  content  = function(file) {
    segs <- tryCatch(cpd_segments_rv(), error = function(e) NULL)
    gmm  <- tryCatch(cpd_gmm_rv(),      error = function(e) NULL)
    if (!is.null(segs) && !is.null(gmm)) {
      log10_D   <- log10(segs$D_um2_s)
      segs$gmm_component <- ifelse(
        is.finite(log10_D),
        mclust::predict.Mclust(gmm, log10_D)$classification,
        NA_integer_
      )
    }
    write.csv(segs, file, row.names = FALSE)
  }
)
