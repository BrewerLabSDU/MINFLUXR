## ─────────────────────────────────────────────────────────────────────────────
## server/tracking/17_jdd.R  -  panel: jump-distance distribution
## ─────────────────────────────────────────────────────────────────────────────
## Part of tracking_server(). Not a standalone script: its statements are evaluated in
## the server function's environment, in file-name order (server/parts.R).
##
## jdd_bundle (R/est_jdd.R), plot, estimates and download.
##
## INTERFACE (generated: Rscript tools/server_map.R --headers; do not edit by hand)
##   DEFINES for other parts:
##     jdd_bundle
##   READS from other parts:
##     01_state_sigma: .sigma_stamp, trk_sigma_loc_eff
##     07_cycle_time_sequence: trk_blur, trk_cycle_time
##     11_data_frames: trk
## END INTERFACE
## ─────────────────────────────────────────────────────────────────────────────

# ============================================================================
# jdd_bundle() - JUMP-DISTANCE DISTRIBUTION (how many diffusing populations?)
# ============================================================================
#
# GOAL
#   The MSD gives ONE average D for the whole dataset. But a sample often holds
#   several populations at once (e.g. a bound slow fraction + a free fast one).
#   This asks a different question: what does the DISTRIBUTION of single-step
#   jump lengths look like, and how many populations does it take to explain it?
#
# THE MATH
#   For 2D free diffusion each coordinate step is Gaussian,
#        Δx, Δy ~ N(0, s²)  with  s² = 2·D·Δt
#   so the step LENGTH r = sqrt(Δx² + Δy²) follows a RAYLEIGH distribution:
#
#        p(r) = (r / s²) · exp( −r² / (2 s²) ),        s² = 2·D·Δt
#
#   Its mean square is <r²> = 2s² = 4DΔt = 2dDΔt (d = 2), consistent with the
#   MSD relation in the header. Inverting for D:
#
#        D = s² / (2·Δt)        and after sqrt(Δt) normalisation:  D = s'² / d
#
#   N co-existing populations -> a WEIGHTED MIXTURE of N chi-distributed components,
#        p(r) = Σ_k w_k · C_d · r^(d−1) / s_k^d · exp( −r² / (2 s_k²) ),   Σ w_k = 1
#   fitted here by maximum likelihood. Each fitted s_k gives that population's
#   D_k, and w_k its abundance. d = spatial dimensions: d=2 is the classic Rayleigh
#   (Schütz 1997), d=3 the Maxwell form used when z is included. s_k is the PER-AXIS
#   spread either way, so D_k = s_k²/2 regardless of d.
#
# THE MINFLUX ADAPTATION (two routes, pick one)
#   Jumps must be pooled at a COMMON Δt or the mixture is meaningless.
#   (a) sqrt(Δt) normalisation (default): divide every jump by sqrt of its own
#       Δt, which removes the Δt dependence for free diffusion (see header §2b).
#   (b) "Integrated steps only": keep single-cycle steps, so every jump already
#       shares one true Δt and NO normalisation assumption is needed. Cleaner -
#       prefer it when the motion may not be free.
#
# OBJECTS IT CALLS (inputs)
#   trk(), trk_cycle_time(), trk_sigma_loc_eff()
#   inputs: trk_jdd_run (Run trigger), trk_jdd_min_pts, trk_jdd_n_comp,
#           trk_jdd_norm_dt, trk_jdd_integrated_only, trk_jdd_min_dt/max_dt,
#           trk_jdd_sigma_loc, trk_jdd_step_mode (+ bin_width / subsamp_n)
#
# STEP BY STEP
#   1. Split into tracks; optionally coarsen each track first
#      (consecutive steps | bin-then-median | every N-th localization).
#   2. Per track, form step lengths r and their Δt.
#   3. Keep steps inside [min_dt, max_dt]; if "integrated only", additionally
#      require round(Δt/t_c) == 1.
#   4. If normalising, r <- r / sqrt(Δt).
#   5. Histogram the pooled jumps (upper tail truncated at the 99.5th pct so a
#      few outliers do not stretch the axis).
#   6. Fit the N-component chi-d mixture by MLE (Rayleigh in 2D, Maxwell in 3D).
#      Parameters are optimised as log(s_k) and softmax weights, so widths stay
#      positive and weights sum to 1 without constrained optimisation.
#   7. Convert widths to D, subtracting the localization-noise contribution
#      (`loc_bias`), and clamp at 0 - a component at the floor is noise, not motion.
#   8. Compute r_noise = sqrt(2d)·σ_loc, the RMS length of a pure-noise step:
#      anything left of it on the plot is noise-dominated.
#
# RETURNS a list:
#   mids, dens ...... histogram (x, density) for plotting
#   sigs[k], ws[k] .. fitted per-axis widths and mixture weights
#   D_vals[k] ....... per-component D (noise-corrected, clamped >= 0)
#   n_comp, norm_dt, integrated_only, t_cycle, dt0, use_z, d_fit (2 or 3)
#   n_jumps, sigma_loc, r_noise, min_dt_used, max_dt_used, step_mode
#
# CAVEAT: the mixture assumes every population is FREELY diffusing. Confined or
# anomalous motion does not follow this law, and extra components may then be
# fitting curve shape rather than real populations. Keep N small.
# ============================================================================
jdd_bundle <- eventReactive(input$trk_jdd_run, {
  req(isTRUE(input$trk_show_jdd))
  df <- trk(); req(nrow(df) > 0)
  min_pts   <- max(2L, as.integer(input$trk_jdd_min_pts))
  norm_dt   <- isTRUE(input$trk_jdd_norm_dt)
  t_c_jdd   <- tryCatch(trk_cycle_time(), error = function(e) NA_real_)
  # Integrated-only: keep ONLY single-cycle steps (round(Δt/t_c)==1). All jumps are
  # then at the same true Δt, so the fit needs no √Δt-pooling assumption at all - the
  # rigorous MINFLUX-native mode (at the cost of dropping dark-gap-adjacent steps).
  integ_only <- isTRUE(input$trk_jdd_integrated_only) && is.finite(t_c_jdd) && t_c_jdd > 0
  n_comp    <- max(1L, min(4L, as.integer(input$trk_jdd_n_comp)))
  min_dt    <- suppressWarnings(as.numeric(input$trk_jdd_min_dt))
  max_dt    <- suppressWarnings(as.numeric(input$trk_jdd_max_dt))
  sigma_loc <- suppressWarnings(as.numeric(input$trk_jdd_sigma_loc))
  step_mode <- input$trk_jdd_step_mode %||% TRK_DEFAULTS$trk_jdd_step_mode
  bin_width <- suppressWarnings(as.numeric(input$trk_jdd_bin_width))
  subsamp_n <- suppressWarnings(as.integer(input$trk_jdd_subsamp_n))
  if (!is.finite(min_dt)    || min_dt    <= 0) min_dt    <- 0
  if (!is.finite(max_dt)    || max_dt    <= 0) max_dt    <- Inf
  sig_info <- if (is.finite(sigma_loc) && sigma_loc > 0) .sigma_stamp(sigma_loc, "value entered in this panel")
              else NULL
  if (!is.finite(sigma_loc) || sigma_loc <= 0) {   # 0 = auto from app-wide σ_loc
    s_app <- tryCatch(trk_sigma_loc_eff(), error = function(e) NA_real_)
    sigma_loc <- if (is.finite(s_app) && s_app > 0) s_app else 0
    sig_info <- .sigma_stamp(if (sigma_loc > 0) sigma_loc else NA_real_)
  }
  if (!is.finite(bin_width) || bin_width <= 0) bin_width <- 1
  if (!is.finite(subsamp_n) || subsamp_n <  2) subsamp_n <- 4L

  # Real axial data only. A placeholder all-zero z column would make this fit the 3D
  # Maxwell law to what are really 2D distances (dz ≡ 0), biasing D by −33%.
  # Motion blur needs the cycle time. Without it the fit would silently run with no
  # blur correction, which reads D about a third low at MINFLUX cycle times, so the
  # panel asks for the calibration instead of guessing.
  shiny::validate(shiny::need(is.finite(t_c_jdd) && t_c_jdd > 0,
    "Cycle time unknown - run Data Quality & Calibration first (the blur correction needs it)."))
  use_z_jdd <- .has_real_z(df)
  # Collection and fit live in est_jdd.R (pure, tested).
  withProgress(message = "Computing JDD…", value = 0, {
  cj <- jdd_collect_jumps(df, min_pts = min_pts, step_mode = step_mode, bin_width = bin_width,
                          subsamp_n = subsamp_n, min_dt = min_dt, max_dt = max_dt,
                          integ_only = integ_only, t_c = t_c_jdd, norm_dt = norm_dt,
                          use_z = use_z_jdd,
                          progress = function(f) incProgress(f, detail = "collecting steps"))
  shiny::validate(shiny::need(length(cj$jumps) >= 10,
    "Too few jumps - lower 'Min points per TID' or relax filters."))
  jf <- jdd_fit(cj$jumps, cj$dt_all, n_comp = n_comp, use_z = use_z_jdd, norm_dt = norm_dt,
                sigma_loc = sigma_loc,
                R_abs = trk_blur()$R_abs,
                progress = function(f) incProgress(f, detail = "fitting jump-distance mixture"))

  structure(list(mids = jf$mids, dens = jf$dens, sigs = jf$sigs, ws = jf$ws, D_vals = jf$D_vals,
       n_comp = n_comp, norm_dt = norm_dt, n_jumps = length(cj$jumps),
       use_z = use_z_jdd, d_fit = jf$d_fit, dt0 = jf$dt0, inv_dt = jf$inv_dt, sigma_loc = sigma_loc,
       r_noise = jf$r_noise, R_abs = trk_blur()$R_abs,
       step_mode = step_mode, integrated_only = integ_only, t_cycle = t_c_jdd,
       min_dt_used = if (min_dt > 0) min_dt else NA_real_,
       max_dt_used = if (is.finite(max_dt)) max_dt else NA_real_),
    sigma_info = sig_info)
  })  # withProgress
})

output$trk_jdd_plot <- renderPlotly({
  bd <- jdd_bundle()
  shiny::validate(shiny::need(!is.null(bd), "No jump data."))
  dark   <- isTRUE(isolate(input$trk_dark_mode))
  bg_col <- if (dark) "#000000" else "#FFFFFF"
  fg_col <- if (dark) "#FFFFFF" else "#333333"

  r_seq <- seq(0, max(bd$mids) * 1.15, length.out = 500)

  # Same chi-d PDF the MLE fitted (Rayleigh in 2D, Maxwell in 3D) - the drawn curve
  # must be the fitted curve, or the overlay silently disagrees with the table.
  .d_ov <- bd$d_fit %||% 2L
  .C_ov <- 1 / (2^(.d_ov / 2 - 1) * gamma(.d_ov / 2))
  rayleigh_comp <- function(r, sig, w)
    w * .C_ov * r^(.d_ov - 1) / sig^.d_ov * exp(-r^2 / (2 * sig^2))

  total_dens <- Reduce(`+`, lapply(seq_len(bd$n_comp), function(k)
    rayleigh_comp(r_seq, bd$sigs[k], bd$ws[k])
  ))

  comp_cols <- c("#E07040", "#4080D0", "#50B060", "#B040B0")

  fig <- plot_ly() %>%
    add_bars(
      x = bd$mids, y = bd$dens,
      name = "Jumps",
      marker = list(color = "rgba(120,160,220,0.6)",
                    line  = list(color = "rgba(80,120,180,0.8)", width = 0.5))
    ) %>%
    add_lines(
      x = r_seq, y = total_dens,
      name = "Fit (total)",
      line = list(color = "#E03020", width = 2.5)
    )

  if (bd$n_comp > 1) {
    for (k in seq_len(bd$n_comp)) {
      fig <- fig %>% add_lines(
        x = r_seq,
        y = rayleigh_comp(r_seq, bd$sigs[k], bd$ws[k]),
        name = sprintf("Comp %d  D=%.2f", k, bd$D_vals[k]),
        line = list(color = comp_cols[k], width = 1.5, dash = "dash")
      )
    }
  }

  if (is.finite(bd$r_noise %||% NA_real_) && bd$r_noise > 0 && bd$r_noise < max(r_seq)) {
    ymax <- max(c(bd$dens, total_dens), na.rm = TRUE)
    fig <- fig %>% add_lines(
      x = c(bd$r_noise, bd$r_noise), y = c(0, ymax),
      name = "loc-noise floor",
      line = list(color = "#888888", width = 1.5, dash = "dot"),
      hovertemplate = sprintf("loc-noise floor \u2248 %.1f (left = noise-dominated)<extra></extra>", bd$r_noise)
    )
  }

  x_lab <- if (bd$norm_dt) "Step size / \u221a\u0394t  (nm / \u221ams)"
           else            "Step size  (nm)"

  fig %>% layout(
    xaxis = list(title = x_lab,   color = fg_col, zeroline = FALSE),
    yaxis = list(title = "Density", color = fg_col, zeroline = FALSE),
    paper_bgcolor = bg_col, plot_bgcolor = bg_col,
    font    = list(color = fg_col),
    legend  = list(font = list(color = fg_col)),
    bargap  = 0.05,
    uirevision = "jdd"
  )
})

output$trk_jdd_estimates <- renderUI({
  bd <- jdd_bundle()
  if (is.null(bd)) return(NULL)

  show_D    <- isTRUE(bd$norm_dt)
  # \u03C3 is the per-axis spread in both 2D and 3D (chi-d PDF), so the /2 is dimension-
  # independent; only the fitted PDF shape changes (Rayleigh 2D / Maxwell 3D).
  d_formula <- sprintf("D = (\u03C3\u00B2 \u2212 2\u03C3\u1D38\u1D52\u1D9C\u00B2\u27E81/\u0394t\u27E9) / (2(1 \u2212 2R_abs\u27E81/\u0394t\u27E9))   [%s]",
                       if (isTRUE(bd$use_z)) "3D, Maxwell" else "2D, Rayleigh")
  dt0_str   <- if (is.finite(bd$dt0)) sprintf("%.2g ms", bd$dt0) else "n/a"
  max_dt_str <- if (is.finite(bd$max_dt_used)) sprintf("%.4g ms", bd$max_dt_used) else "none"

  sigma_note <- if (show_D && bd$sigma_loc > 0 && is.finite(bd$dt0))
    tags$p(sprintf("\u03C3\u1D38\u1D52\u1D9C = %.1f nm  |  median \u0394t = %s  |  noise subtracted: %.2g nm\u00B2/ms  |  blur factor %.2f",
                   bd$sigma_loc, dt0_str, 2 * bd$sigma_loc^2 * (bd$inv_dt %||% (1 / bd$dt0)),
                   1 - 2 * (bd$R_abs %||% 0) * (bd$inv_dt %||% (1 / bd$dt0))),
           style = "font-size:0.82em; color:#888; margin:2px 0;")
  else NULL

  col_headers <- if (show_D) c("", "Weight", "Width (\u03C3)", "Diffusivity (D, nm\u00B2/ms)")
                 else        c("", "Weight", "Width (\u03C3)")

  rows <- lapply(seq_len(bd$n_comp), function(k) {
    cells <- list(
      tags$td(sprintf("Component %d", k),          style = "padding:2px 8px;"),
      tags$td(sprintf("%.1f%%", bd$ws[k] * 100),   style = "padding:2px 8px;"),
      tags$td(sprintf("\u03C3 = %.1f nm/\u221Ams", bd$sigs[k]), style = "padding:2px 8px;")
    )
    if (show_D) {
      d_val <- bd$D_vals[k]
      d_str <- if (d_val <= 0 || d_val < 1e-6)
        tags$span("D ≈ 0", tags$span(style="color:#c05020; font-size:0.88em;",
          " (within loc. noise - step sizes explained by σₗₒₓ alone)"))
      else
        sprintf("D = %.3g", d_val)
      cells <- c(cells, list(tags$td(d_str, style = "padding:2px 8px;")))
    }
    do.call(tags$tr, cells)
  })

  no_d_warn <- if (!show_D)
    tags$p("\u26A0 D values are not shown: without \u221A\u0394t normalisation \u03C3 has units nm, not nm/\u221Ams, so \u03C3\u00B2/2 is not a diffusion coefficient.",
           style = "font-size:0.82em; color:#c05020; margin:4px 0;")
  else NULL

  min_dt_str <- if (is.finite(bd$min_dt_used) && bd$min_dt_used > 0)
                   sprintf("%.4g ms", bd$min_dt_used) else "none"
  mode_str <- switch(bd$step_mode %||% "consec",
    binmed  = "bin-then-median",
    subsamp = "subsampled",
    "consecutive steps"
  )
  tags$div(
    tags$p(sprintf("n = %d jumps  |  %s  |  %s  |  \u0394t window: [%s, %s]",
                   bd$n_jumps, d_formula, mode_str, min_dt_str, max_dt_str),
           style = "font-size:0.85em; color:#888; margin:4px 0;"),
    sigma_note,
    no_d_warn,
    tags$table(
      style = "font-size:0.88em; border-collapse:collapse;",
      tags$thead(tags$tr(
        lapply(col_headers, function(h)
          tags$th(h, style = "padding:2px 8px; border-bottom:1px solid #ccc;"))
      )),
      tags$tbody(rows)
    )
  )
})

output$trk_dl_jdd_csv <- downloadHandler(
  filename = function() paste0("jdd_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
  content  = function(file) {
    bd <- jdd_bundle()
    shiny::validate(shiny::need(!is.null(bd), "No JDD computed."))
    # chi-d PDF as fitted (Rayleigh 2D / Maxwell 3D) - keep the export in step
    # with the plot and the fitted parameters.
    .d_ov <- bd$d_fit %||% 2L
    .C_ov <- 1 / (2^(.d_ov / 2 - 1) * gamma(.d_ov / 2))
    rayleigh_comp <- function(r, sig, w)
      w * .C_ov * r^(.d_ov - 1) / sig^.d_ov * exp(-r^2 / (2 * sig^2))
    out <- data.frame(mid_nm = bd$mids, density = bd$dens)
    out$fit_total <- Reduce(`+`, lapply(seq_len(bd$n_comp), function(k)
      rayleigh_comp(bd$mids, bd$sigs[k], bd$ws[k])
    ))
    for (k in seq_len(bd$n_comp))
      out[[paste0("comp_", k)]] <- rayleigh_comp(bd$mids, bd$sigs[k], bd$ws[k])
    write.csv(out, file, row.names = FALSE)
  }
)
