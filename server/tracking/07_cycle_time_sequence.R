# ===========================  END FILTER PIPELINE  =========================

# ===========================================================================
#  CYCLE-TIME & TSA SUPPORT  (read by Stage 2 above and the TSA section)
#  Cycle-time estimate (cluster / hardware / manual), Time-Signature Analysis
#  readouts, sequence-JSON import, and the adaptive gap-split threshold.
# ===========================================================================

# Half-sample mode (Robertson–Cryer / Bickel): a bin-free Δt mode
# estimator.  Keeps the densest half-sample until ≤3 points remain.
# half_sample_mode() lives in R/est_clean.R.

# Typical Δt (ms) of the RAW data - the empirical acquisition interval, bin-free
# via the half-sample mode. Computed on trk_raw_df (upstream of ALL cleaning) so
# gap-split / end-trim / quality can use the cycle time with no circular
# dependency, and so metadata-removed localizations don't distort the Δt.
trk_pipeline_modal_dt <- reactive({
  df <- .trk_timing_df()
  if (is.null(df) || !nrow(df) || !all(c("tid", "tim") %in% names(df))) return(NA_real_)
  df  <- df[order(df$tid, df$tim), ]
  dts <- unlist(tapply(df$tim, df$tid, function(t) diff(sort(t))), use.names = FALSE)
  dts <- dts[is.finite(dts) & dts > 0]
  if (!length(dts)) return(NA_real_)
  half_sample_mode(dts)
})

# ---- Data-driven cycle estimate (a light TSA) ------------------------------
# Groups consecutive Δt into integer-cycle CLUSTERS and takes the median
# spacing between adjacent cluster centres as the cycle time.  Because the
# cluster *spacing* is the cycle time regardless of which cluster is most
# populated, this recovers the true (atomic) cycle even when the data are poor
# and the typical interval is itself 2–3 cycles (dense regions / frequent
# dropouts) - exactly where the half-sample mode would mis-report the cycle.
# Returns list(cycle, sw, n_clusters, centers, modal, frac_one_cycle).
# Reads trk_raw_df (upstream of all cleaning → no circular dependency; raw Δt is
# the true acquisition timing).
# The acquisition timing comes from the ORIGINAL upload of the active dataset, so
# Plot Data (which rechunks the filtered survivors) cannot shift the cycle time.
.trk_timing_df <- function() {
  df <- tryCatch(trk_raw_original(), error = function(e) NULL)
  if (is.null(df) || !nrow(df)) df <- tryCatch(trk_raw_df(), error = function(e) NULL)
  df
}

trk_cluster_cycle <- reactive({
  df <- .trk_timing_df()
  if (is.null(df) || !nrow(df) || !all(c("tid", "tim") %in% names(df))) return(NULL)
  cycle_cluster(df$tid, df$tim)                     # R/est_clean.R
})

# Regression-based cycle/SW estimate - Bela's EOD/SW fit, reduced to a line fit.
# The cluster centres lie at ~ k·cycle + SW for integer cycle counts k, so a
# least-squares fit of centre ~ k gives slope = cycle time and intercept = SW.
# This is the same quantity Bela obtains by fitting EOD_OVERHEAD / SW_OVERHEAD to
# the time signal (EOD is just a reparametrisation of the slope through the
# hardware vertices). Using ALL centres at once is more precise than the median
# spacing and is robust to an occasional missing cycle. Returns list(cycle, sw,
# n, r2, eod) - eod filled only when the hardware vertices/dwell are known.
trk_regression_cycle <- reactive({
  cc <- tryCatch(trk_cluster_cycle(), error = function(e) NULL)
  if (is.null(cc) || length(cc$centers) < 2L ||
      !is.finite(cc$cycle) || cc$cycle <= 0) return(NULL)
  centers <- sort(as.numeric(cc$centers))
  idx <- round(centers / cc$cycle)                   # integer cycle index per centre
  keep <- is.finite(idx) & idx >= 1L & !duplicated(idx)
  centers <- centers[keep]; idx <- idx[keep]
  if (length(unique(idx)) < 2L) return(NULL)
  fit <- tryCatch(stats::lm(centers ~ idx), error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  cf <- stats::coef(fit)
  slope <- unname(cf[2]); intercept <- unname(cf[1])
  if (!is.finite(slope) || slope <= 0) return(NULL)
  r2 <- tryCatch(summary(fit)$r.squared, error = function(e) NA_real_)
  # Optional EOD per move, if the hardware geometry is supplied: cycle = dwell +
  # EOD × vertices  →  EOD = (cycle − dwell) / vertices  (µs).
  eod <- NA_real_
  if (isTRUE(input$trk_hw_enable)) {
    dwell_us <- suppressWarnings(as.numeric(input$trk_hw_dwell_us))
    pat <- input$trk_hw_pattern %||% "6"
    nv  <- if (identical(pat, "other")) suppressWarnings(as.integer(input$trk_hw_vertices))
           else suppressWarnings(as.integer(pat))
    if (is.finite(dwell_us) && is.finite(nv) && nv >= 1L)
      eod <- (slope * 1000 - dwell_us) / nv
  }
  list(cycle = slope, sw = max(0, intercept), n = length(idx), r2 = r2, eod = eod)
})

# Cycle time (ms) from the hardware model - TSA formula derived from Vogler:
#   cycle = T_DWELL + EOD_OVERHEAD × NUM_VERTICES

trk_hw_cycle_time <- reactive({
  if (!isTRUE(input$trk_hw_enable)) return(NA_real_)
  dwell_us <- suppressWarnings(as.numeric(input$trk_hw_dwell_us))
  eod_us   <- suppressWarnings(as.numeric(input$trk_hw_eod_us))
  pat      <- input$trk_hw_pattern %||% "6"
  nv <- if (identical(pat, "other"))
          suppressWarnings(as.integer(input$trk_hw_vertices))
        else suppressWarnings(as.integer(pat))
  if (!is.finite(dwell_us) || dwell_us <= 0 || !is.finite(nv) || nv < 1) return(NA_real_)
  if (!is.finite(eod_us) || eod_us < 0) eod_us <- 0
  (dwell_us + eod_us * nv) / 1000   # µs → ms
})

# ---- SINGLE SOURCE OF TRUTH for the acquisition cycle time (ms) -------------
# Priority: (1) hardware model if enabled & valid, (2) manual override in the
# Advanced MSD panel if set > 0, (3) data-driven CLUSTER cycle (robust atomic
# cycle), (4) half-sample-mode modal Δt as a last resort.  Everything that
# needs a cycle time reads this so they can never disagree.
trk_cycle_time <- reactive({
  hw <- trk_hw_cycle_time()
  if (is.finite(hw) && hw > 0) return(hw)
  man <- suppressWarnings(as.numeric(input$trk_msd_cycle_ms))
  if (is.finite(man) && man > 0) return(man)
  # Optional: Bela-style regression cycle time (toggle), preferred over the
  # cluster median when enabled and valid.
  if (isTRUE(input$trk_cycle_use_regression)) {
    rc <- tryCatch(trk_regression_cycle(), error = function(e) NULL)
    if (!is.null(rc) && is.finite(rc$cycle) && rc$cycle > 0) return(rc$cycle)
  }
  cc <- tryCatch(trk_cluster_cycle(), error = function(e) NULL)
  if (!is.null(cc) && is.finite(cc$cycle) && cc$cycle > 0) return(cc$cycle)
  trk_pipeline_modal_dt()
})

# Standard-MSD fit window bounds, entered in CYCLES, converted to ms via the
# cycle time (Vogler's fit_limits, default 5–50).  0 means "off".
msd_min_lag_ms <- reactive({
  cyc <- suppressWarnings(as.numeric(input$trk_msd_min_lag_cyc))
  tc  <- trk_cycle_time()
  if (length(cyc) == 0 || !is.finite(cyc) || cyc <= 0 || !is.finite(tc) || tc <= 0) return(0)
  cyc * tc
})
# Upper end of the same fit window, also entered in CYCLES and converted to ms.
# Same "0 means off" convention as msd_min_lag_ms above.
msd_fit_max_ms <- reactive({
  cyc <- suppressWarnings(as.numeric(input$trk_msd_fit_max_cyc))   # user value, in cycles
  tc  <- trk_cycle_time()                                          # ms per cycle
  # Any missing/nonsensical input (or unknown cycle time) disables the bound.
  if (length(cyc) == 0 || !is.finite(cyc) || cyc <= 0 || !is.finite(tc) || tc <= 0) return(0)
  cyc * tc                                                         # cycles -> milliseconds
})

# Status / cross-check line shown in the Acquisition Parameters panel
output$trk_cycle_time_status <- renderUI({
  req(trk_has_data())
  mdt <- trk_pipeline_modal_dt()
  hw  <- trk_hw_cycle_time()
  tc  <- trk_cycle_time()
  rc0 <- tryCatch(trk_regression_cycle(), error = function(e) NULL)
  reg_active <- isTRUE(input$trk_cycle_use_regression) &&
                !is.null(rc0) && is.finite(rc0$cycle) && rc0$cycle > 0
  src <- if (is.finite(hw) && hw > 0) "hardware model"
         else {
           man <- suppressWarnings(as.numeric(input$trk_msd_cycle_ms))
           if (is.finite(man) && man > 0) "manual override"
           else if (reg_active) "regression fit (EOD/SW)"
           else "cluster cycle (data)"
         }
  if (!is.finite(tc) || tc <= 0)
    return(tags$p(style = "font-size:0.82em; color:#888; margin:2px 0;",
                  "Cycle time: not yet determined (load data)."))
  rows <- list(
    tags$p(style = "font-size:0.84em; color:#2070a0; margin:2px 0;",
           sprintf("Cycle time in use: %.4g ms  (source: %s)", tc, src))
  )
  if (is.finite(hw) && hw > 0 && is.finite(mdt) && mdt > 0) {
    # The atomic cycle (model) is the MINIMUM interval, so it should sit
    # slightly BELOW the typical/modal Δt (which includes small per-cycle
    # variability).  ratio = model/modal ≈ 0.8–1.0 is the healthy regime.
    ratio <- hw / mdt
    verdict <- if (ratio >= 0.75 && ratio <= 1.05)
                 "consistent ✓ (atomic cycle just below typical Δt, as expected)"
               else if (ratio > 1.05)
                 "⚠ model exceeds typical Δt - dwell/overhead likely too high"
               else
                 "⚠ model well below typical Δt - check dwell, or many dark cycles inflate Δt"
    col <- if (ratio >= 0.75 && ratio <= 1.05) "#2a8a3e"
           else if (ratio > 1.05) "#c0392b" else "#b07000"
    msg <- sprintf("Atomic cycle (model) %.4g ms vs typical Δt %.4g ms (ratio %.2f) - %s",
                   hw, mdt, ratio, verdict)
    rows <- c(rows, list(tags$p(style = sprintf("font-size:0.82em; color:%s; margin:2px 0;", col), msg)))
  } else if (is.finite(mdt) && mdt > 0) {
    rows <- c(rows, list(tags$p(style = "font-size:0.82em; color:#888; margin:2px 0;",
           sprintf("Empirical modal Δt = %.4g ms.", mdt))))
  }
  # Data-driven cluster cycle and the RAW single-cycle fraction (cluster_cycle
  # reads raw data). Descriptive only, no verdict.
  cc <- tryCatch(trk_cluster_cycle(), error = function(e) NULL)
  if (!is.null(cc) && is.finite(cc$cycle) && cc$cycle > 0 && is.finite(cc$modal)) {
    fo <- cc$frac_one_cycle
    rows <- c(rows, list(
      tags$p(style = "font-size:0.82em; color:#555; margin:2px 0;",
             sprintf("Cluster (atomic) cycle = %.4g ms  |  typical interval (HSM) = %.4g ms  |  %d cycle-clusters detected",
                     cc$cycle, cc$modal, cc$n_clusters)),
      if (is.finite(fo))
        tags$p(style = "font-size:0.82em; color:#555; margin:2px 0;",
               sprintf("Single-cycle (integrated) steps, raw: %.0f%%", fo * 100))))
  }
  # Regression (EOD/SW) cycle - Bela-style line fit, shown for side-by-side
  # comparison with the cluster median (and used as the source when toggled on).
  if (!is.null(rc0) && is.finite(rc0$cycle) && rc0$cycle > 0) {
    eod_txt <- if (is.finite(rc0$eod)) sprintf("  |  EOD ≈ %.2f µs/move", rc0$eod) else ""
    r2_txt  <- if (is.finite(rc0$r2)) sprintf("  |  R² = %.4f", rc0$r2) else ""
    rows <- c(rows, list(
      tags$p(style = sprintf("font-size:0.82em; color:%s; margin:2px 0;",
                             if (reg_active) "#2070a0" else "#777"),
             sprintf("%s Regression fit: cycle = %.4g ms  |  SW = %.4g ms%s  (%d cluster centres)%s",
                     if (reg_active) "▶" else "○",
                     rc0$cycle, rc0$sw, eod_txt, rc0$n, r2_txt))))
  }
  do.call(tags$div, rows)
})

# ms-equivalent hint shown under the cycle-based gap threshold
output$trk_split_gap_ms_hint <- renderUI({
  req(isTRUE(input$trk_split_gaps_enable))
  cyc <- suppressWarnings(as.numeric(input$trk_split_gap_cycles))
  mdt <- trk_cycle_time()
  if (!is.finite(cyc) || cyc <= 0 || !is.finite(mdt) || mdt <= 0) return(NULL)
  tags$p(style = "font-size:0.81em; color:#2070a0; margin:2px 0;",
         sprintf("≈ %.4g ms  (at cycle time = %.4g ms)", cyc * mdt, mdt))
})

# Adaptive (Bela) threshold readout: the data-driven cycle threshold plus how
# many tracks/localizations survive the split + min-length filter.
output$trk_split_adaptive_hint <- renderUI({
  req(isTRUE(input$trk_split_gaps_enable), isTRUE(input$trk_split_adaptive))
  thr <- trk_adaptive_gap_cycles()
  tot <- trk_tot_cycles()
  if (!is.finite(thr) || length(tot) < 2)
    return(tags$p(style = "color:#888; font-size:0.81em;", "Need data to compute the adaptive threshold."))
  med <- as.integer(stats::median(tot)); sdv <- stats::sd(tot)
  cyc <- trk_cycle_time()
  n_split <- sum(tot >= thr)
  df_in    <- trk_raw_df()   # split is the first stage now → its input is raw
  df_split <- trk_gap_split()
  df_final <- tryCatch(trk_spatial_filtered(), error = function(e) df_split)
  n_tid_in    <- if (nrow(df_in)    > 0) length(unique(df_in$tid))    else 0L
  n_tid_split <- if (nrow(df_split) > 0) length(unique(df_split$tid)) else 0L
  n_tid_final <- if (nrow(df_final) > 0) length(unique(df_final$tid)) else 0L
  tags$div(style = "font-size:0.81em; color:#2a6; background:#f0f7f0; border-radius:4px; padding:5px 7px; margin:3px 0;",
    tags$b(sprintf("Auto threshold = %d cycles", thr)),
    sprintf("  (median %d + SD %.2f + 1)", med, sdv), tags$br(),
    sprintf("≈ %.4g ms · %s of %s steps exceed it",
            thr * cyc, format(n_split, big.mark = ","), format(length(tot), big.mark = ",")), tags$br(),
    sprintf("Tracks: %s in → %s after split+min-length → %s displayed",
            format(n_tid_in, big.mark = ","), format(n_tid_split, big.mark = ","),
            format(n_tid_final, big.mark = ",")),
    if (n_tid_final < n_tid_split)
      tags$span(style = "color:#888;", tags$br(),
        "(length / duration / spatial filters remove the rest - their ranges now expand to cover segments)"))
})

# End-trim funnel hint (localizations removed from track ends).
output$trk_endtrim_hint <- renderUI({
  req(isTRUE(input$trk_endtrim_enable))
  din  <- tryCatch(trk_gap_split_d(),   error = function(e) NULL)
  dout <- tryCatch(trk_end_trimmed_d(), error = function(e) NULL)
  if (is.null(din) || is.null(dout) || !nrow(din))
    return(tags$p(style = "color:#888; font-size:0.81em;", "Need data to compute end-trim."))
  n_in <- nrow(din); n_out <- nrow(dout)
  tags$div(style = "font-size:0.81em; color:#2a6; background:#f0f7f0; border-radius:4px; padding:5px 7px; margin:3px 0;",
    tags$b(sprintf("Localizations: %s → %s kept", format(n_in, big.mark = ","), format(n_out, big.mark = ","))),
    sprintf("  (%s trimmed from track ends, %.1f%%)",
            format(n_in - n_out, big.mark = ","), 100 * (n_in - n_out) / max(1L, n_in)))
})

# Track-quality gate funnel hint.
output$trk_quality_hint <- renderUI({
  req(isTRUE(input$trk_quality_enable))
  din  <- tryCatch(trk_gap_split_d(),        error = function(e) NULL)
  dout <- tryCatch(trk_quality_filtered_d(), error = function(e) NULL)
  if (is.null(din) || is.null(dout) || !nrow(din))
    return(tags$p(style = "color:#888; font-size:0.81em;", "Need data to compute track quality."))
  n_in  <- length(unique(din$tid)); n_out <- length(unique(dout$tid))
  t_c   <- tryCatch(trk_cycle_time(), error = function(e) NA_real_)
  thr   <- max(0, min(1, suppressWarnings(as.numeric(input$trk_quality_min_integ %||% 0.6))))
  H     <- suppressWarnings(as.numeric(input$trk_quality_darktol %||% 3))
  tags$div(style = "font-size:0.81em; color:#2a6; background:#f0f7f0; border-radius:4px; padding:5px 7px; margin:3px 0;",
    tags$b(sprintf("Tracks: %s → %s kept (%.0f%%)",
                   format(n_in, big.mark = ","), format(n_out, big.mark = ","),
                   100 * n_out / max(1L, n_in))),
    sprintf("  at continuity ≥ %.2f (tolerance %.1f)", thr, H),
    if (is.finite(t_c)) sprintf("  ·  cycle = %.4g ms", t_c) else "")
})


# .crlb_doughnut() (photon-budget σ, pattern-aware) and .tcp_geom() live in est_cve.R.

# TCP diameter of one iteration, in nm. Sequence files store the pattern size as
# the dimensionless `patGeoFactor`; L = patGeoFactor x TCP_L_PER_GEOFACTOR_NM.
# A 3D iteration may carry a [lateral, axial] pair (e.g. zline gives [0.8, 4.0]);
# the LATERAL element is the one that governs escape from the pattern, and it is
# also the L the D_max form is written for, so the first element is used.
.itr_L_nm <- function(it) {
  pg <- it[["patGeoFactor"]]
  if (is.null(pg)) return(NA_real_)
  pg <- suppressWarnings(as.numeric(unlist(pg)))
  pg <- pg[is.finite(pg) & pg > 0]
  if (!length(pg)) return(NA_real_)
  base <- suppressWarnings(as.numeric(input$trk_geofactor_base %||% TCP_L_PER_GEOFACTOR_NM))
  if (!is.finite(base) || base <= 0) base <- TCP_L_PER_GEOFACTOR_NM
  pg[1] * base
}

# ── trk_seq_params() - read the Imspector imaging sequence JSON ─────────────
# GOAL: pull the ACQUISITION SETTINGS straight from the microscope's sequence
# file, so the app does not have to infer everything from the data alone.
# Used as a cross-check on the data-driven cycle time, and to pre-fill the
# gap-split threshold (from maxOffTime) and the hardware pattern.
#
# The file's `Itr` array is one entry per MINFLUX iteration; the last one is
# the tracking iteration we care about.
# returns: list(ok = FALSE, msg = ...) on any problem, else the parsed fields.
#          Callers must always check $ok first.
trk_seq_params <- reactive({
  req(input$trk_seq_json)
  fp <- input$trk_seq_json$datapath
  seq <- tryCatch(jsonlite::fromJSON(fp, simplifyVector = FALSE),
                  error = function(e) NULL)
  if (is.null(seq) || is.null(seq$Itr) || length(seq$Itr) == 0)
    return(list(ok = FALSE, msg = "Not a recognised MINFLUX sequence (no 'Itr' array)."))
  itrs <- seq$Itr
  getp  <- function(it, k) { v <- it[[k]]; if (is.null(v)) NA else v }
  getm  <- function(it, k) { m <- it[["Mode"]]; if (is.null(m)) return(NA); v <- m[[k]]; if (is.null(v)) NA else v }
  pats  <- vapply(itrs, function(it) { v <- getm(it, "pattern"); if (is.na(v)) NA_character_ else as.character(v) }, character(1))
  dims  <- suppressWarnings(vapply(itrs, function(it) as.integer(getm(it, "dim")), integer(1)))
  # maxOffTime may sit on the final iteration or at the top level, and is often
  # the string "unspecified" on stock tracking sequences. Prefer the iteration
  # value, fall back to the top-level one, else leave NA for the 3 ms default.
  moff  <- suppressWarnings(as.numeric(getp(itrs[[length(itrs)]], "maxOffTime")))   # seconds
  if (!is.finite(moff)) moff <- suppressWarnings(as.numeric(getp(seq, "maxOffTime")))
  prep  <- suppressWarnings(as.numeric(getp(itrs[[length(itrs)]], "patRepeat")))
  pht   <- suppressWarnings(as.numeric(getp(itrs[[length(itrs)]], "phtLimit")))
  pdw   <- suppressWarnings(as.numeric(getp(itrs[[length(itrs)]], "patDwellTime"))) # seconds (timeout)
  # Timing inputs for the blur model (est_blur.R): the CFR check adds a centre dwell
  # of ctrDwellFactor x dwell when ccrLimit > 0, and `damping` > 0 smooths the
  # exported positions across steps, which no blur model describes.
  ccr   <- suppressWarnings(as.numeric(getp(itrs[[length(itrs)]], "ccrLimit")))
  cdf   <- suppressWarnings(as.numeric(getp(seq, "ctrDwellFactor")))
  dmp   <- suppressWarnings(as.numeric(getp(seq, "damping")))
  final_pat <- tail(pats[!is.na(pats)], 1L)
  if (length(final_pat) == 0) final_pat <- NA_character_
  dim_final <- suppressWarnings(as.integer(tail(dims[!is.na(dims)], 1L)))
  # TCP diameter per iteration, and the final one - the iteration a tracking
  # sequence reiterates once locked on, which is the L the D_max ceiling needs.
  # The whole ladder is returned too, so a sequence that loops an earlier
  # iteration can be read off and entered by hand.
  L_by_itr <- vapply(itrs, .itr_L_nm, numeric(1))
  L_final  <- if (length(L_by_itr)) L_by_itr[length(L_by_itr)] else NA_real_
  # `headstart` is the loop-back offset: -1 (every tracking sequence seen) means
  # the last iteration is the one repeated. Reported, not acted on.
  hs <- suppressWarnings(as.integer(getp(seq, "headstart")))
  list(
    ok = TRUE,
    n_itr = length(itrs),
    patterns = pats,
    final_pattern = final_pat,
    vertices = .pattern_vertices(final_pat),
    dim = if (length(dim_final) && is.finite(dim_final)) dim_final else NA_integer_,
    maxOffTime_ms = if (is.finite(moff)) moff * 1000 else NA_real_,
    patRepeat = prep,
    phtLimit = pht,
    patDwellTime_ms = if (is.finite(pdw)) pdw * 1000 else NA_real_,
    cfr_on = is.finite(ccr) && ccr > 0,
    ctrDwellFactor = if (is.finite(cdf)) cdf else NA_real_,
    damping = if (is.finite(dmp)) dmp else NA_real_,
    L_by_itr = L_by_itr,
    L_tcp_nm = L_final,
    headstart = if (is.finite(hs)) hs else NA_integer_
  )
})

# ── trk_blur() - ONE MOTION-BLUR COEFFICIENT FOR EVERY TRACKING ESTIMATOR ───
# GOAL: the blur R of a localization, from the photon-collection timeline of the
# acquisition (est_blur.R), converted once into the forms the estimators use:
# R_abs (ms, JDD / D-map / change-point / HMM), R_msd (relative to the cycle, MSD
# fits and the confinement null) and R_step (the CVE).
# Inputs, best first: the sequence JSON (pattern, dwell, pattern repeats, CFR
# centre dwell); otherwise the data (collection window eco/efo, the hardware
# pattern setting). The step dt is the mean single-cycle interval of the data.
# eta (passes per localization) = eco / (efo * dwell) when the dwell is known.
trk_blur <- reactive({
  df  <- trk(); t_c <- tryCatch(trk_cycle_time(), error = function(e) NA_real_)
  out <- list(R_step = BLUR_R_STEP_FALLBACK, R_abs = NA_real_, R_msd = MSD_R_FALLBACK,
              source = "fallback (standard 2D sequence)", damping = NA_real_, detail = NULL)
  if (is.null(df) || !nrow(df) || !is.finite(t_c) || t_c <= 0) return(out)
  d  <- df[order(df$tid, df$tim), c("tid", "tim", intersect(c("eco", "efo"), names(df)))]
  m  <- nrow(d); same <- c(FALSE, d$tid[-1] == d$tid[-m])
  dt <- c(NA_real_, diff(d$tim)); dt[!same] <- NA_real_
  one <- is.finite(dt) & round(dt / t_c) == 1
  dt_int <- if (any(one)) mean(dt[one]) else t_c
  med <- function(v) { v <- v[is.finite(v) & v > 0]; if (length(v)) stats::median(v) else NA_real_ }
  t_coll <- if (all(c("eco", "efo") %in% names(d))) med(d$eco / d$efo) * 1000 else NA_real_   # ms
  sp <- tryCatch(if (is.null(input$trk_seq_json)) NULL else trk_seq_params(), error = function(e) NULL)
  has_seq <- !is.null(sp) && isTRUE(sp$ok) && is.finite(sp$patDwellTime_ms)
  b <- NULL; note <- NULL
  if (has_seq) {
    eta <- if (is.finite(t_coll)) max(1, round(t_coll / sp$patDwellTime_ms)) else 1
    b <- blur_R_step(dt_int, sp$patDwellTime_ms, pattern = sp$final_pattern,
                     n_vert = sp$vertices, pat_repeat = sp$patRepeat %||% 1, eta = eta,
                     cfr = isTRUE(sp$cfr_on), ctr_factor = sp$ctrDwellFactor %||% 0)
    src <- sprintf("sequence file (%s, %.0f µs dwell%s%s)", sp$final_pattern %||% "?",
                   1000 * sp$patDwellTime_ms, if (isTRUE(sp$cfr_on)) ", CFR" else "",
                   if (eta > 1) sprintf(", %d passes", eta) else "")
    # A sequence whose collection alone outlasts the measured step does not describe
    # this data (e.g. an imaging sequence loaded with tracking data); one whose final
    # pattern is axial has no lateral R. Both fall back to the data, and say so.
    if (is.null(b) || !is.finite(b$R)) {
      note <- sprintf("the sequence's final pattern (%s) gives no lateral blur", sp$final_pattern %||% "?")
      b <- NULL
    } else if (isTRUE(b$overrun)) {
      note <- sprintf("the sequence needs %.0f µs per localization but the data step is %.0f µs, so it does not match this data",
                      1000 * b$t_loc_ms, 1000 * dt_int)
      b <- NULL
    }
  }
  if (is.null(b) && is.finite(t_coll)) {
    # Lateral polygon only: an axial pattern loaded from a sequence can set the
    # vertex control to 2, which says nothing about the lateral collection.
    nv <- suppressWarnings(as.integer(input$trk_hw_vertices %||% 6L))
    if (!is.finite(nv) || nv < 3L) nv <- 6L
    b <- blur_R_step(dt_int, t_coll, n_vert = nv)
    src <- sprintf("data (collection window eco/efo %.0f µs, %d vertices)", 1000 * t_coll, nv)
  }
  if (!is.null(b) && is.finite(b$R)) {
    o <- blur_offsets(b$R, dt_int, t_c)
    out <- c(o, list(source = src, detail = b))
  } else {
    out$R_abs <- blur_offsets(out$R_step, dt_int, t_c)$R_abs
  }
  out$note    <- note
  out$damping <- if (!is.null(sp) && isTRUE(sp$ok)) sp$damping else NA_real_
  out$dt_int  <- dt_int
  out
})

# Apply the parsed sequence: set the gap-split threshold from maxOffTime, enable
# gap splitting, and (when the lateral pattern is unambiguous) pre-fill the
# hardware scan pattern. The cycle time stays data-driven.
observeEvent(input$trk_seq_json, {
  sp <- trk_seq_params()
  if (is.null(sp) || !isTRUE(sp$ok)) return(invisible(NULL))
  # gap-split threshold = round(maxOffTime / cycle time)  [in cycles]
  cyc <- isolate(trk_cycle_time())
  if (is.finite(sp$maxOffTime_ms) && is.finite(cyc) && cyc > 0) {
    n_cyc <- max(2, round(sp$maxOffTime_ms / cyc))
    updateNumericInput(session, "trk_split_gap_cycles", value = n_cyc)
    updateCheckboxInput(session, "trk_split_gaps_enable", value = TRUE)
  }
  # pre-fill the hardware pattern if it maps to a clean lateral polygon
  nv <- sp$vertices
  if (is.finite(nv) && nv %in% c(3L, 4L, 6L)) {
    updateSelectInput(session, "trk_hw_pattern", selected = as.character(nv))
  } else if (is.finite(nv)) {
    updateSelectInput(session, "trk_hw_pattern", selected = "other")
    updateNumericInput(session, "trk_hw_vertices", value = nv)
  }
  # D_max ceiling: take L from the sequence rather than a default. The ceiling
  # scales as L^2, so this is the input that most needs to come from the file.
  if (is.finite(sp$L_tcp_nm) && sp$L_tcp_nm > 0)
    updateNumericInput(session, "trk_L_tcp", value = round(sp$L_tcp_nm, 1))
  # TCP geometry for the photon-budget sigma (off-centre coefficient, axes
  # spanned). Only the patterns with a modelled geometry are pre-selected; any
  # other final pattern leaves the current choice alone.
  fp <- tolower(sp$final_pattern %||% "")
  if (fp %in% names(TCP_PATTERN_GEOM))
    updateSelectInput(session, "trk_tcp_pattern", selected = fp)
}, ignoreInit = TRUE)

output$trk_seq_json_status <- renderUI({
  if (is.null(input$trk_seq_json)) return(NULL)
  sp <- trk_seq_params()
  if (is.null(sp) || !isTRUE(sp$ok))
    return(tags$p(style = "color:#b00; font-size:0.85em;", sp$msg %||% "Could not parse file."))
  cyc <- trk_cycle_time()
  gap_cyc <- if (is.finite(sp$maxOffTime_ms) && is.finite(cyc) && cyc > 0)
               max(2, round(sp$maxOffTime_ms / cyc)) else NA
  upat <- unique(sp$patterns[!is.na(sp$patterns)])
  tags$div(style = "font-size:0.84em; background:#eef5ee; border:1px solid #cdd; border-radius:4px; padding:6px 8px; margin-top:-2px;",
    tags$b("Sequence loaded "), sprintf("(%d iterations)", sp$n_itr), tags$br(),
    sprintf("Patterns: %s", paste(upat, collapse = ", ")), tags$br(),
    sprintf("Tracking iteration: %s  ·  dim = %s  ·  L_TCP = %s",
            if (length(sp$final_pattern) && !is.na(sp$final_pattern)) sp$final_pattern else "?",
            if (is.finite(sp$dim)) sprintf("%dD", sp$dim) else "?",
            if (is.finite(sp$L_tcp_nm)) sprintf("%.1f nm", sp$L_tcp_nm) else "?"), tags$br(),
    if (any(is.finite(sp$L_by_itr)))
      tags$span(sprintf("Iteration L: %s nm",
                        paste(ifelse(is.finite(sp$L_by_itr), sprintf("%.0f", sp$L_by_itr), "?"),
                              collapse = " · ")), tags$br()) else NULL,
    if (is.finite(sp$maxOffTime_ms))
      sprintf("maxOffTime = %.3g ms  →  gap-split @ %s cycles",
              sp$maxOffTime_ms, ifelse(is.na(gap_cyc), "?", as.character(gap_cyc)))
    else "maxOffTime: not found",
    if (is.finite(sp$phtLimit)) tags$span(tags$br(), sprintf("Photon limit = %g  ·  patRepeat = %g",
                                                             sp$phtLimit, sp$patRepeat)) )
})

# ── TSA: timing summary ─────────────────────────────────────────────────────
# The robust cluster cycle against the modal Δt, and the single-cycle fraction
# before and after the filters. Descriptive only: the single-cycle fraction
# mostly reflects the sequence (photon limit against photons per pass) and the
# dye's brightness, so it carries no verdict.
output$trk_tsa_quality <- renderUI({
  req(isTRUE(input$trk_show_tsa), tsa_ran())
  cc <- tryCatch(trk_cluster_cycle(), error = function(e) NULL)
  if (is.null(cc) || !is.finite(cc$cycle) || cc$cycle <= 0)
    return(tags$p(style = "color:#888;", "Upload tracking data to run the time-signature analysis."))
  modal <- cc$modal
  ratio <- if (is.finite(modal) && cc$cycle > 0) modal / cc$cycle else NA_real_
  frac1 <- cc$frac_one_cycle                       # RAW (metadata-filtered)
  # Filtered fraction from the same data the Integrated/Dark charts below use
  # (dt_vectors -> trk(), after gate / end-trim / length / duration / spatial filters).
  ccb       <- tryCatch(cycle_class_bundle(), error = function(e) NULL)
  frac_filt <- if (!is.null(ccb) && is.finite(ccb$frac_integ)) ccb$frac_integ else NA_real_
  cyc_in_use <- trk_cycle_time()
  tags$div(style = "border-left:5px solid #2c6fbb; padding:6px 10px; background:#fafafa;",
    tags$div(style = "font-size:1.05em; font-weight:bold; color:#2c6fbb;", "Timing summary"),
    tags$table(style = "font-size:0.86em; margin-top:4px;",
      tags$tr(tags$td(style="padding-right:12px;", "Atomic cycle (cluster)"),
              tags$td(tags$b(sprintf("%.4g ms", cc$cycle)))),
      tags$tr(tags$td("Modal Δt (raw)"),
              tags$td(sprintf("%.4g ms", modal))),
      tags$tr(tags$td("Modal / cluster ratio"),
              tags$td(tags$b(sprintf("%.2f×", ratio)))),
      if (is.finite(frac_filt)) tags$tr(tags$td("Single-cycle fraction (filtered)"),
              tags$td(tags$b(sprintf("%.0f%%", 100 * frac_filt)))),
      if (is.finite(frac1)) tags$tr(tags$td("Single-cycle fraction (raw)"),
              tags$td(sprintf("%.0f%%", 100 * frac1))),
      tags$tr(tags$td("Cycle time used downstream"),
              tags$td(tags$b(sprintf("%.4g ms", cyc_in_use))))))
})

# Per-interval TOTAL acquisition-cycle counts across all meta-filtered tracks:
# TOT_i = round((Δt_i − SW_overhead) / cycle_time).  This is Bela's "TOT" signal
# (integrated + dark cycles spanned by each consecutive step).
trk_tot_cycles <- reactive({
  df  <- .trk_timing_df()   # split operates on raw → estimate its threshold on raw
  if (is.null(df)) return(integer(0))
  cyc <- trk_cycle_time()
  if (nrow(df) == 0 || !is.finite(cyc) || cyc <= 0) return(integer(0))
  cc  <- tryCatch(trk_cluster_cycle(), error = function(e) NULL)
  sw  <- if (!is.null(cc) && is.finite(cc$sw)) cc$sw else 0
  # Bela computes suggest_split_thresholds AFTER filter_by_length(50), so in
  # adaptive mode only tracks meeting the minimum length count.
  min_len <- if (isTRUE(input$trk_split_adaptive))
               max(1L, as.integer(input$trk_split_min_len %||% 50L)) else NULL
  total_cycles(df$tid, df$tim, cyc, sw, min_len)    # R/est_clean.R
})

# Bela's adaptive split threshold (in cycles): med_th = floor(median + std + 1)
# of the TOT signal, capped by the hardware "particle lost" ceiling
# high_th = round(maxOffTime / cycle).  maxOffTime comes from the uploaded
# sequence JSON if present, else Bela's 3 ms default.
trk_adaptive_gap_cycles <- reactive({
  moff_ms <- tryCatch({ sp <- trk_seq_params(); if (!is.null(sp) && isTRUE(sp$ok) &&
                          is.finite(sp$maxOffTime_ms)) sp$maxOffTime_ms else 3.0 },
                      error = function(e) 3.0)
  adaptive_gap_cycles(trk_tot_cycles(), trk_cycle_time(), moff_ms)   # R/est_clean.R
})

# (The filter pipeline stages 1.5–4 are defined together near Stage 1, under
# the "TRACKING FILTER PIPELINE" banner. trk_tot_cycles / trk_adaptive_gap_cycles
# above are the cycle-time support those stages read.)

output$trk_chunk_stat <- renderUI({
  req(input$trk_chunk_sel)
  if (input$trk_chunk_sel == "__all__") {
    return(tags$div(style = "color:#555;", "All TIDs selected  -  loading full file."))
  }
  df <- read_chunk_file(input$trk_chunk_sel)
  if (is.null(df)) return(NULL)
  tids <- length(unique(as.character(df$tid)))
  tags$div(style="color:#555;",
           paste0("Selected chunk: ", format(nrow(df), big.mark=","), " rows * ",
                  tids, " TIDs"))
})

# ---- Per-TID summary (from tid_summary.csv saved at chunking time) ----------
tid_summary_rv <- reactive({
  req(!is.null(rv$manifest_dir))
  f <- file.path(rv$manifest_dir, "tid_summary.csv")
  if (!file.exists(f)) return(NULL)
  readr::read_csv(f, show_col_types = FALSE)
})

# Lifetime summary stat shown above the metadata filter section
output$trk_lifetime_stat <- renderUI({
  ts <- tid_summary_rv()
  if (is.null(ts) || !nrow(ts)) return(NULL)
  avg <- mean(ts$duration_ms, na.rm = TRUE)
  med <- stats::median(ts$duration_ms, na.rm = TRUE)
  n   <- nrow(ts)
  tags$div(
    style = "color:#555; margin-bottom:4px; font-size:0.9em;",
    tags$b("Trace lifetimes -"),
    sprintf("%d traces \u00B7 mean %.0f ms \u00B7 median %.0f ms", n, avg, med)
  )
})

# Stage 3 (duration) auto-detection. BOUNDS reflect the length-filtered population
# but update at PLOT cadence only - keyed on data_ver() (upload / Plot Data / Reset
# Ranges), NOT live. So the range describes the currently PLOTTED data and stops
# jumping as sliders move. The VALUE is PRESERVED (clamped into the new bounds) so a
# manual minimum (e.g. 50 ms) sticks; it snaps to full only on new upload / Reset
# Ranges / post-rechunk, signalled by reset_dur_filter_flag.
output$trk_dur_filter_ui <- renderUI({
  reset_ranges_trigger()
  data_ver()                                                           # update at PLOT cadence (upload / Plot Data), not live
  df <- tryCatch(isolate(trk_len_filtered()), error = function(e) NULL)
  dv <- if (!is.null(df) && nrow(df) > 0 && all(c("tid","tim") %in% names(df)))
          (df %>% dplyr::group_by(tid) %>%
             dplyr::summarise(dur = max(tim) - min(tim), .groups = "drop"))$dur
        else numeric(0)
  if (!length(dv)) {                          # pipeline not ready → rechunk stats, else upload
    rs_re <- trk_rechunk_stats_rv(); rs_or <- trk_raw_stats_rv()
    dv <- if (!is.null(rs_re) && length(rs_re$tid_duration)) rs_re$tid_duration
          else if (!is.null(rs_or) && length(rs_or$tid_duration)) rs_or$tid_duration else numeric(0)
  }
  if (!length(dv)) return(NULL)
  bounds <- range(dv, na.rm = TRUE)
  val <- .filt_value("dur", bounds)          # stored filter, clamped (full range if none)
  shiny::isolate(reset_dur_filter_flag(FALSE))
  numericRangeInput("trk_dur_filter", "Trace duration filter (ms)",
                    value = val, min = bounds[1], max = bounds[2])
})
