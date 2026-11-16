## ─────────────────────────────────────────────────────────────────────────────
## server/tracking/02_data_quality_cve.R  -  panel: data quality and calibration (CVE)
## ─────────────────────────────────────────────────────────────────────────────
## Part of tracking_server(). Not a standalone script: its statements are evaluated in
## the server function's environment, in file-name order (server/parts.R).
##
## trk_cve_bundle (CVE sigma and D, per-track sigma, linking QC, photon-budget
## sigma, headline D), publishing sigma, and the panel's readouts.
##
## INTERFACE (generated: Rscript tools/server_map.R --headers; do not edit by hand)
##   DEFINES for other parts:
##     trk_cve_bundle
##   READS from other parts:
##     01_state_sigma: .continuity_scores, trk_raw_stats_rv,
##                     trk_sigma_loc_cve, trk_sigma_loc_photon
##     04_datasets_upload: trk_raw_original
##     06_filter_pipeline: .a_in_use, trk_dur_filtered_d, trk_end_trimmed_d,
##                         trk_gap_split_d, trk_len_filtered_d,
##                         trk_meta_filtered_d, trk_quality_filtered_d,
##                         trk_spatial_filtered
##     07_cycle_time_sequence: trk_blur, trk_cluster_cycle, trk_cycle_time,
##                             trk_seq_params
##     11_data_frames: trk
##     15_msd: .msd_settings
##     24_calibration: trk_cal_profiles
## END INTERFACE
## ─────────────────────────────────────────────────────────────────────────────

# ============================================================================
# trk_cve_bundle() - LOCALIZATION PRECISION (σ_loc) + DATA-QUALITY CALIBRATION
# ============================================================================
#
# GOAL
#   Measure how precisely the microscope located each molecule (σ_loc, nm) and
#   score how clean the acquisition timing is. σ_loc becomes the app-wide
#   "noise floor" that the D-map, JDD, confinement, HMM, directionality and
#   speed tools all subtract. Run this FIRST - everything else calibrates off it.
#
# THE IDEA (CVE - Vestergaard, Blainey & Flyvbjerg 2014)
#   Localization error appears in TWO consecutive steps with OPPOSITE sign: the
#   noise at point i is the END of step i-1 and the START of step i. That makes
#   neighbouring steps anti-correlated by exactly −σ², while real diffusion
#   increments are independent and contribute nothing to that correlation. So
#
#       σ_loc² = −< Δx_i · Δx_{i+1} >          (per axis, averaged over x,y[,z])
#
#   That is the R = 0 (instantaneous exposure) case. A MINFLUX localization
#   integrates photons across the TCP cycle, so the code below uses the blur-aware
#   form σ² = R·<Δx²> + (2R − 1)·cov₁ (see "Motion blur" in the bundle).
#   Bonus: D falls out too, noise-corrected, with no MSD fit at all:
#       D = <Δ²>/(2Δt) + cov₁/Δt
#
# THE MINFLUX ADAPTATION
#   The algebra above assumes a CONSTANT Δt. MINFLUX Δt is irregular. So the
#   covariance is taken ONLY over pairs where BOTH steps are single-cycle
#   ("integrated", n_cyc == 1) and share a midpoint inside one track. That
#   reconstructs constant-Δt sub-runs, at the cost of discarding steps next to
#   a dark gap.
#
# OBJECTS IT CALLS (inputs)
#   trk_gap_split() ......... data after gap-splitting (cleaning stage 1)
#   trk_spatial_filtered() .. the fully filtered data, for the "after filters" score
#   trk_raw_original() ...... the frozen upload, for the "raw" reference score
#   trk_cycle_time() ........ base cycle time, to classify integrated vs dark
#   .continuity_scores(), .boot_ci()
#   inputs: trk_cve_run (Run trigger), trk_sig_min_integ, trk_sig_skip,
#           trk_quality_darktol
#
# STEP BY STEP
#   1. Build per-step Δx, Δy[, Δz] and Δt within each track (never across tracks).
#   2. Classify each step: n_cyc = round(Δt / t_c); n_cyc == 1 is "integrated".
#   3. Data-quality scores: fraction of integrated steps, raw vs after-filters,
#      plus a dark-tolerance-weighted continuity score.
#   4. POOLED σ_loc: average Δx_i·Δx_{i+1} over consecutive integrated pairs,
#      negate, square-root. Also the bias-corrected D_cve.
#   5. PER-TRACK σ_loc: repeat per track (tracks with enough integrated steps),
#      skipping the first few localizations - MINFLUX needs a moment to lock on
#      and those opening points read worse. Gives median, IQR, CV, bootstrap CI,
#      and a σ-vs-EFO correlation (precision SHOULD improve with photons; a
#      negative correlation validates the whole estimate).
#   6. PHOTON WEIGHTING (optional): fit σ ∝ EFO^p across tracks, assign every
#      localization a precision from its own EFO, and inverse-variance-weight
#      the pooled σ and D so bright, well-localized steps count for more.
#   7. QC flags: implausible jumps (|Δr| > 5× expected for that gap - likely
#      mis-links) and suspicious tracks (degenerate or outlier σ).
#
# RETURNS a list. The fields consumers actually use:
#   sigma ......... pooled σ_loc (nm)
#   sig_pt_med .... per-track median σ_loc    D_cve ..... noise-corrected D
#   degenerate .... TRUE if σ² came out <= 0 (drift dominates -> unresolvable)
#   frac_integ, cont_raw, cont_filt ... data-quality scores
#   pertrack_sigma .... per-track table (tid, sigma, n_integ, phot_med, var_step)
#   r1_norm / noise_frac ... how much of the step variance is localization noise
#   n_implaus, n_suspect, ...
#
# WHEN IT FAILS: strong directed motion makes consecutive steps POSITIVELY
# correlated, so σ² goes negative and σ_loc is genuinely unresolvable. That is
# reported as `degenerate` - enter σ_loc manually in that case.
# ============================================================================
trk_cve_bundle <- eventReactive(input$trk_cve_run, {
  # The plotted data of the active dataset (filtered, applied with Plot Data), as
  # every other analysis reads it, so all panels describe the same tracks.
  df <- tryCatch(trk(), error = function(e) NULL)
  req(!is.null(df), nrow(df) > 0, all(c("x","y","tid","tim") %in% names(df)))
  t_c   <- tryCatch(trk_cycle_time(), error = function(e) NA_real_)
  # Real axial data only - a placeholder all-zero z column would dilute both σ_loc
  # and D_cve (see .has_real_z).
  has_z <- .has_real_z(df)
  # Step arrays from est_cve.R, shared with cve_core() and the diagnostics below.
  st <- cve_steps(df, t_c, has_z)
  df <- st$df; n <- st$n; same <- st$same; dt <- st$dt
  dx <- st$dx; dy <- st$dy; dz <- st$dz; n_cyc <- st$n_cyc
  integ_step <- st$integ_step; n_steps <- st$n_steps; frac_integ <- st$frac_integ

  # Same integrated fraction, but on the FULLY-FILTERED data (post quality-gate /
  # end-trim / length / duration / spatial) - the data the analyses actually use.
  # Lets the panel show raw → after-filters and whether the filters helped.
  frac_integ_of <- function(d) {
    if (is.null(d) || !nrow(d) || !all(c("tid","tim") %in% names(d))) return(NA_real_)
    d <- d[order(d$tid, d$tim), ]; m <- nrow(d)
    sm  <- c(FALSE, d$tid[-1] == d$tid[-m])
    dtt <- c(NA_real_, diff(d$tim)); dtt[!sm] <- NA_real_
    nc  <- if (is.finite(t_c) && t_c > 0) round(dtt / t_c) else rep(1L, m)
    isi <- sm & is.finite(nc) & (nc == 1L); isi[is.na(isi)] <- FALSE
    ns  <- sum(sm, na.rm = TRUE)
    if (ns > 0) sum(isi) / ns else NA_real_
  }
  frac_integ_filt <- frac_integ_of(tryCatch(trk_spatial_filtered(), error = function(e) NULL))
  # RAW reference = TRUE frozen upload data (never filtered / rechunked), NOT df
  # (= trk_gap_split, which reads the current - possibly filtered - manifest). This
  # keeps "raw" at the uploaded value (e.g. 52%) so the raw→analyzed jump is real.
  frac_integ_raw  <- frac_integ_of(tryCatch(trk_raw_original(),      error = function(e) NULL))

  # Dataset-level cycle-CONTINUITY (tolerance-weighted; same weighting as the
  # continuity gate) on TRUE RAW (frozen upload) vs ANALYZED (filtered). This reflects
  # the dark-cycle tolerance, so a dataset with many gentle 1–2-cycle gaps scores
  # higher than its bare integrated fraction.
  cont_of <- function(d, H) {
    if (is.null(d) || !nrow(d) || !all(c("tid","tim") %in% names(d))) return(NA_real_)
    d <- d[order(d$tid, d$tim), ]; m <- nrow(d)
    sm  <- c(FALSE, d$tid[-1] == d$tid[-m])
    dtt <- c(NA_real_, diff(d$tim)); dtt[!sm] <- NA_real_
    kk  <- pmax(1L, as.integer(round(dtt / t_c)))
    ww  <- if (is.finite(H) && H > 0) 2^(-(kk - 1) / H) else as.numeric(kk == 1L)
    if (sum(sm, na.rm = TRUE) > 0) mean(ww[sm], na.rm = TRUE) else NA_real_
  }
  H_disp    <- suppressWarnings(as.numeric(input$trk_quality_darktol %||% 3))
  cont_raw  <- cont_of(tryCatch(trk_raw_original(), error = function(e) NULL), H_disp)  # TRUE raw (frozen at upload)
  cont_filt <- cont_of(tryCatch(trk_spatial_filtered(), error = function(e) NULL), H_disp)  # analyzed (live filtered)

  # Collection window from the photon budget: t_coll = ECO / EFO.
  .med_of <- function(cn) if (cn %in% names(df)) {
                            v <- df[[cn]][is.finite(df[[cn]]) & df[[cn]] > 0]
                            if (length(v)) stats::median(v) else NA_real_
                          } else NA_real_
  eco_m <- .med_of("eco"); efo_m <- .med_of("efo")
  t_coll <- if (is.finite(eco_m) && is.finite(efo_m) && efo_m > 0) eco_m / efo_m * 1000 else NA_real_
  n_vert <- {
    sp <- tryCatch(if (is.null(input$trk_seq_json)) NULL else trk_seq_params(), error = function(e) NULL)
    v  <- if (!is.null(sp) && isTRUE(sp$ok)) sp$vertices else NA_integer_
    if (!is.finite(v)) suppressWarnings(as.integer(input$trk_hw_vertices %||% 6L)) else v
  }
  if (!is.finite(n_vert) || n_vert < 1L) n_vert <- 6L
  # CVE core (covariance, blur bounds, σ, D) lives in est_cve.R.
  cc <- cve_core(dx, dy, dz, dt, integ_step, has_z, t_c, t_coll, n_vert,
                 R_step = trk_blur()$R_step)
  shiny::validate(shiny::need(!is.null(cc),
    "Too few consecutive integrated steps for a CVE estimate - relax filters or check the cycle time."))
  pair <- cc$pair; prod_x <- cc$prod_x; prod_y <- cc$prod_y; prod_z <- cc$prod_z
  vals <- rep(NA_real_, cc$n_pairs)                       # only its length is read below
  cov1 <- cc$cov1; dt_int <- cc$dt_int; R_lo <- cc$R_lo; R_hi <- cc$R_hi; R_blur <- cc$R_blur
  msq <- cc$msq; sigma2 <- cc$sigma2; sigma <- cc$sigma
  sigma_xy <- cc$sigma_xy; sigma_z <- cc$sigma_z; sig_aniso <- cc$sig_aniso
  D_cve <- cc$D_cve; r1_norm <- cc$r1_norm; noise_frac <- cc$noise_frac

  # ── Per-track σ_loc heterogeneity ──────────────────────────────────────────
  # One CVE σ per track (tracks with ≥ n_min integrated steps), skipping the
  # first `skip` localizations - MINFLUX iterates onto the molecule, so the
  # opening localizations carry lock-on transient and read worse. The spread of
  # this distribution says whether a single global σ_loc is a good approximation.
  n_min_pt <- max(10L, as.integer(input$trk_sig_min_integ %||% 50L))
  skip_pt  <- max(0L,  as.integer(input$trk_sig_skip %||% 5L))
  pos  <- sequence(rle(df$tid)$lengths)                     # within-track index (df sorted)
  post <- pos > skip_pt                                     # localization kept (past lock-on)
  step_post  <- same & post & c(FALSE, post[-n])            # step between two kept locs
  integ_post <- step_post & is.finite(n_cyc) & (n_cyc == 1L); integ_post[is.na(integ_post)] <- FALSE
  pair_pt <- integ_post & c(integ_post[-1], FALSE)          # consecutive integrated kept pair
  tid_pt  <- df$tid[pair_pt]
  pv <- c(prod_x[pair_pt], prod_y[pair_pt]); tv <- c(tid_pt, tid_pt)
  if (has_z) { pv <- c(pv, prod_z[pair_pt]); tv <- c(tv, tid_pt) }
  okpv <- is.finite(pv); pv <- pv[okpv]; tv <- tv[okpv]
  # ECO = photon counts (correct for σ ∝ 1/√N); EFO = emission frequency (fallback).
  phot_col <- if ("eco" %in% names(df) && any(is.finite(df$eco) & df$eco > 0)) "eco"
              else if ("efo" %in% names(df) && any(is.finite(df$efo) & df$efo > 0)) "efo"
              else NA_character_
  # Per-track σ carries the same blur correction. Note what it can and cannot
  # say: at R ≈ 1/6 the R*⟨Δx²⟩ term is ~90% of σ², so per-track differences here
  # track the per-track STEP SIZE far more than the photon budget.
  # Kept for the outlier/degenerate QC below, not as a precision measurement -
  # the per-localization photon estimate is what the panel plots.
  sq_x_pt <- dx^2; sq_y_pt <- dy^2; sq_z_pt <- if (has_z) dz^2 else rep(NA_real_, n)
  msq_by_tid <- {
    v <- c(sq_x_pt[integ_post], sq_y_pt[integ_post])
    t2 <- c(df$tid[integ_post], df$tid[integ_post])
    if (has_z) { v <- c(v, sq_z_pt[integ_post]); t2 <- c(t2, df$tid[integ_post]) }
    ok <- is.finite(v); tapply(v[ok], t2[ok], mean)
  }
  pertrack_sigma <- if (length(pv)) {
    cov_by <- tapply(pv, tv, mean)                    # per-track cov1
    s2_by  <- vapply(names(cov_by), function(nm) {
                mq <- msq_by_tid[[nm]]
                if (is.null(mq) || !is.finite(mq)) -1
                else R_blur * mq + (2 * R_blur - 1) * cov_by[[nm]]
              }, numeric(1))
    n_by   <- tapply(integ_post, df$tid, sum)[names(s2_by)]
    # Photon predictor for the precision model. σ_loc ≈ L_TCP/(2√N) depends on the
    # photon COUNT N, so ECO (effective counts) is the physically correct variable.
    # EFO is a FREQUENCY (emission rate); it only tracks N when the dwell time is
    # constant, which MINFLUX tracking violates (dwell varies, especially around dark
    # cycles), so it is a fallback only. Using ECO also gives a testable prediction:
    # the log-log slope should be ≈ -0.5.
    phot_by <- if (!is.na(phot_col)) {
                 kl <- post & is.finite(df[[phot_col]]) & df[[phot_col]] > 0
                 tapply(df[[phot_col]][kl], df$tid[kl], stats::median)[names(s2_by)]
               } else stats::setNames(rep(NA_real_, length(s2_by)), names(s2_by))
    # Per-track step variance, reported with the per-track table.
    v_by <- tapply(dx[integ_post]^2 + dy[integ_post]^2,
                   df$tid[integ_post], mean, na.rm = TRUE)[names(s2_by)]
    data.frame(tid = names(s2_by), sigma = sqrt(pmax(0, as.numeric(s2_by))),
               n_integ = as.integer(n_by), phot_med = as.numeric(phot_by),
               var_step = as.numeric(v_by),
               degenerate = as.numeric(s2_by) <= 0, stringsAsFactors = FALSE)
  } else NULL
  if (!is.null(pertrack_sigma))
    pertrack_sigma <- pertrack_sigma[is.finite(pertrack_sigma$n_integ) &
                                     pertrack_sigma$n_integ >= n_min_pt, , drop = FALSE]
  sig_pt_ok  <- if (!is.null(pertrack_sigma) && nrow(pertrack_sigma))
                  pertrack_sigma$sigma[!pertrack_sigma$degenerate] else numeric(0)
  sig_pt_med <- if (length(sig_pt_ok) >= 5L) stats::median(sig_pt_ok) else NA_real_
  sig_pt_iqr <- if (length(sig_pt_ok) >= 5L)
                  as.numeric(diff(stats::quantile(sig_pt_ok, c(.25, .75)))) else NA_real_
  sig_pt_cv  <- if (is.finite(sig_pt_med) && sig_pt_med > 0 && length(sig_pt_ok) >= 5L)
                  stats::sd(sig_pt_ok) / mean(sig_pt_ok) else NA_real_
  # 95% CI on the per-track-median σ_loc (bootstrap over tracks).
  sig_ci <- if (length(sig_pt_ok) >= 5L) .boot_ci(sig_pt_ok, stats::median, B = 1000L) else c(lo = NA_real_, hi = NA_real_)

  # ── Linking / sanity QC ─────────────────────────────────────────────────────
  # (1) Implausible jumps: a step whose |Δr| far exceeds what free diffusion +
  #     localization noise allow for its OWN gap (E[|Δr|²] = 2d·D·Δt + 2d·σ²).
  #     A big jump after one cycle is suspicious; after a long dark gap it isn't
  #     (the Δt term grows), so this naturally accounts for dark cycles. Flags
  #     mis-links / re-acquisition errors. (2) Suspicious tracks: per-track σ_loc
  #     that is degenerate (σ²<0) or a far outlier - candidate bad/mis-linked tracks.
  sig_use <- if (is.finite(sig_pt_med)) sig_pt_med else sigma
  d_dim   <- if (has_z) 3 else 2
  n_implaus <- NA_integer_; frac_implaus <- NA_real_
  if (is.finite(D_cve) && D_cve > 0 && is.finite(sig_use)) {
    dr2  <- dx^2 + dy^2 + if (has_z) dz^2 else 0            # NA where !same
    exp2 <- 2 * d_dim * D_cve * dt + 2 * d_dim * sig_use^2  # expected E[|Δr|²] per step
    ratio2 <- dr2 / exp2
    implaus <- same & is.finite(ratio2) & (ratio2 > 25)     # |Δr| > 5× expected RMS
    n_implaus    <- sum(implaus, na.rm = TRUE)
    frac_implaus <- if (n_steps > 0) n_implaus / n_steps else NA_real_
  }
  n_suspect <- if (!is.null(pertrack_sigma) && nrow(pertrack_sigma) >= 5L) {
    sg  <- pertrack_sigma$sigma
    ok  <- !pertrack_sigma$degenerate & is.finite(sg)
    med0 <- stats::median(sg[ok]); iqr0 <- stats::IQR(sg[ok])
    sum(pertrack_sigma$degenerate | (is.finite(sg) & sg > med0 + 3 * iqr0), na.rm = TRUE)
  } else NA_integer_
  n_tracks_pt <- if (!is.null(pertrack_sigma)) nrow(pertrack_sigma) else NA_integer_

  # ── Photon channel: DIAGNOSTIC ONLY (photon weighting has been REMOVED) ─────
  # There used to be an inverse-variance "photon-weighted" CVE here (fit σ ∝ photons^p,
  # then weight steps by the predicted σ). It is gone: photon counts differ between
  # tracks for several reasons besides precision (mobility, bleaching, emitter
  # brightness, background), so weighting by them would weight the tracks by
  # whatever drives those differences. The CVE weights every step equally.
  #   noise_frac - how much of the step variance is localization noise (−2·r₁).

  # ── σ_loc from the photon budget - the value the app actually uses ─────────
  # σ_i = L_TCP / (a*sqrt(N_i)) per localization, N = eco, a the measured effective
  # prefactor, plus the pattern's off-centre term (.crlb_doughnut, est_cve.R). It is
  # per-localization rather than pooled, and does not inherit the step-size
  # inflation that population heterogeneity and drift push into the CVE's ⟨Δx²⟩
  # term. The CVE σ above is kept as a diagnostic only.

  L_tcp <- suppressWarnings(as.numeric(input$trk_L_tcp %||% DMAX_L_TCP_DEFAULT))
  if (!is.finite(L_tcp) || L_tcp <= 0) L_tcp <- DMAX_L_TCP_DEFAULT
  ecov  <- if ("eco" %in% names(df)) df$eco[is.finite(df$eco) & df$eco > 0] else numeric(0)
  # The emitter drifts off the TCP centre while photons are collected; D from the
  # CVE above and the single-cycle interval set how far. The pattern sets how
  # fast precision falls off with that offset (TCP_PATTERN_GEOM).
  nu_beam <- suppressWarnings(as.numeric(input$trk_beam_fwhm %||% TCP_BEAM_FWHM_NM))
  if (!is.finite(nu_beam) || nu_beam <= 0) nu_beam <- TCP_BEAM_FWHM_NM
  tcp_pat <- input$trk_tcp_pattern %||% TCP_PATTERN_DEFAULT
  D_off   <- if (is.finite(D_cve) && D_cve > 0) D_cve else 0
  t_off   <- if (is.finite(dt_int) && dt_int > 0) dt_int else 0
  r_off   <- sqrt(.tcp_geom(tcp_pat)$dim * D_off * t_off)   # rms offset over the exposure, nm
  a_pref <- .a_in_use(tcp_pat)
  sig_phot_i   <- if (length(ecov)) .crlb_doughnut(L_tcp, nu_beam, ecov, D = D_off, t_loc = t_off,
                                                   pref = a_pref, pattern = tcp_pat) else numeric(0)
  sig_phot_med <- if (length(sig_phot_i)) stats::median(sig_phot_i) else NA_real_
  sig_phot_iqr <- if (length(sig_phot_i) >= 5L)
                    as.numeric(diff(stats::quantile(sig_phot_i, c(.25, .75)))) else NA_real_
  eco_med      <- if (length(ecov)) stats::median(ecov) else NA_real_
  # Headline D: the upper lag regime of the ensemble MSD, with the MSD panel's fit
  # settings. The CVE
  # above uses single-cycle steps, the lower regime, and is kept as a check.
  ms  <- .msd_settings(df)
  cch <- tryCatch(trk_cluster_cycle(), error = function(e) NULL)
  hd  <- tryCatch(msd_headline(df, t_c, if (!is.null(cch) && is.finite(cch$sw)) cch$sw else 0,
                               ms$window, ms$track_weight, ms$iterative, ms$linear_tau, ms$R,
                               if (ms$has_z) ms$coef[["d3"]] else ms$coef[["d2"]], ms$min_pts, ms$has_z),
                  error = function(e) NULL)

  list(sigma = sigma, sigma2 = sigma2, cov1 = cov1, n_pairs = length(vals),
       R_blur = R_blur, msq = msq, n_vert = n_vert, t_coll = t_coll,
       sig_phot = sig_phot_i, sig_phot_med = sig_phot_med, sig_phot_iqr = sig_phot_iqr,
       eco_med = eco_med, L_tcp = L_tcp, n_phot = length(ecov),
       nu_beam = nu_beam, r_off = r_off, tcp_pattern = tcp_pat,
       sig_simple = if (length(ecov)) stats::median(L_tcp / (a_pref * sqrt(ecov))) else NA_real_,
       a_pref = a_pref,
       sigma_xy = sigma_xy, sigma_z = sigma_z, sig_aniso = sig_aniso,
       D_cve = D_cve, dt_int = dt_int, t_cycle = t_c, has_z = has_z,
       frac_integ = frac_integ_raw, frac_integ_filt = frac_integ_filt,
       cont_raw = cont_raw, cont_filt = cont_filt, cont_H = H_disp,
       n_steps = n_steps, degenerate = sigma2 <= 0,
       pertrack_sigma = pertrack_sigma, sig_pt_med = sig_pt_med,
       sig_pt_iqr = sig_pt_iqr, sig_pt_cv = sig_pt_cv, sig_pt_lo = sig_ci[["lo"]],
       sig_pt_hi = sig_ci[["hi"]], n_pt = length(sig_pt_ok),
       n_implaus = n_implaus, frac_implaus = frac_implaus,
       n_suspect = n_suspect, n_tracks_pt = n_tracks_pt,
       phot_col = phot_col,   # photon column used for the diagnostics
       r1_norm = r1_norm, noise_frac = noise_frac,
       D_msd = if (!is.null(hd)) hd$D else NA_real_, msd_window = ms$window,
       msd_fit_lags = if (!is.null(hd)) hd$fit_lags else NULL,
       msd_fit_tau  = if (!is.null(hd)) hd$fit_tau else NULL)
})

# Publish the calibration result as the app-wide σ_loc.
#
# The PHOTON estimate is the one used. The CVE σ is not: even with the blur
# correction applied, ~90% of its σ² comes from the R*⟨Δx²⟩ term, so it is
# essentially sqrt(R) x rms step size - a function of the assumed R and of the
# step distribution, which population heterogeneity and drift both inflate. The
# data's own contribution, cov1, is a few-percent residual. It stays on screen as
# a cross-check on L, not as a precision measurement.
observeEvent(trk_cve_bundle(), {
  b <- trk_cve_bundle(); if (is.null(b)) return()
  if (is.finite(b$sig_phot_med) && b$sig_phot_med > 0) trk_sigma_loc_photon(b$sig_phot_med)
  val <- if (is.finite(b$sigma) && b$sigma > 0 && !isTRUE(b$degenerate)) b$sigma else NA_real_
  if (is.finite(val)) trk_sigma_loc_cve(val)
}, ignoreInit = TRUE)

# GOAL: report the CVE calibration as numbers plus short verdicts.
# The reasoning that used to be printed here on every run (why photons are not
# used as weights, what sigma_loc feeds downstream, why only integrated pairs
# qualify) now lives in the panel's collapsible Method box and in README
# "CVE (Covariance-Based Estimator)". On-screen text stays numeric.
output$trk_cve_status <- renderUI({
  b <- tryCatch(trk_cve_bundle(), error = function(e) NULL)
  if (is.null(b))
    return(htmltools::div(style = "color:#888; font-size:0.9em;",
      "Press ", htmltools::tags$b("Estimate ", sigloc(), " & quality"), " to calibrate."))

  AMB <- "#B8860B"; RED <- "#c03030"; DIM <- "#888"
  row  <- function(...) htmltools::div(style = "margin:1px 0;", ...)
  dim_ <- function(...) htmltools::span(style = sprintf("color:%s;", DIM), ...)
  warn <- function(...) htmltools::span(style = sprintf("color:%s;", AMB), ...)
  num  <- function(x, fmt, na = "-") if (is.finite(x)) sprintf(fmt, x) else na

  # ── sigma_loc ──────────────────────────────────────────────────────────────
  # The PHOTON estimate is the value the app uses; the CVE follows as a check.
  sigma_line <- if (is.finite(b$sig_phot_med))
    row(htmltools::tags$b(sigloc(), sprintf(" %.2f nm", b$sig_phot_med)),
        dim_(sprintf("  photon budget  \u00b7  %s  \u00b7  L = %.1f nm  \u00b7  a = %.2f%s  \u00b7  offset %.1f nm  \u00b7  median ECO %.0f  \u00b7  %s locs",
                     b$tcp_pattern %||% TCP_PATTERN_DEFAULT, b$L_tcp, b$a_pref %||% NA_real_,
                     if (is.finite(b$nu_beam %||% NA_real_) && b$nu_beam > 0)
                       sprintf("  \u00b7  \u03bd = %.0f nm", b$nu_beam) else "",
                     b$r_off %||% 0,
                     b$eco_med %||% NA_real_, format(b$n_phot %||% 0L, big.mark = ","))),
        if (is.finite(b$sig_simple %||% NA_real_))
          dim_(sprintf("  \u00b7  %.2f before the motion term", b$sig_simple)) else NULL)
  else
    row(warn(htmltools::HTML("&#9888; "), "no ECO column, so ", sigloc(),
             " cannot come from the photon budget; the CVE value below is used instead."))


  # Motion blur used by every estimator (trk_blur(), est_blur.R), and where it came
  # from. Damping in the sequence smooths positions ACROSS steps, which no blur model
  # describes, so it is flagged rather than corrected.
  bl <- tryCatch(trk_blur(), error = function(e) NULL)
  blur_line <- if (!is.null(bl)) row(
    dim_(sprintf("motion blur R = %.3f per step, from %s", bl$R_step, bl$source)),
    if (!is.null(bl$note)) warn(sprintf("  (sequence file not used: %s)", bl$note)) else NULL,
    if (is.finite(bl$damping %||% NA_real_) && bl$damping > 0)
      warn(htmltools::HTML(sprintf(
        "  &#9888; sequence damping = %g: positions are smoothed across steps, so no blur model applies and D from single steps is approximate",
        bl$damping))) else NULL) else NULL
  # CVE sigma, a sanity read on the covariance only: most of it is the assumed blur
  # term, and where blur outweighs the noise it has no solution.
  # The covariance-based sigma is a diagnostic only (nothing downstream uses it),
  # so an unresolvable value is reported plainly rather than as a data problem.
  cve_line <- if (isTRUE(b$degenerate))
    row(dim_("CVE σ: no solution (blur outweighs the noise in the step covariance); diagnostic only"))
  else if (is.finite(b$sigma))
    row(dim_(sprintf("CVE σ %.2f nm (blur-aware, R = %.3f), diagnostic only",
                     b$sigma, b$R_blur %||% NA_real_)))
  else NULL

  # ── lateral vs axial (3D only) - the pooled sigma averages x, y and z, so when
  #    the axial precision differs from the lateral one it misstates the floor the
  #    x-y-only tools (D-map, 2D JDD/MSD) subtract. The axial spread is usually the
  #    TIGHTER one in MINFLUX, so the warning is two-sided: an axial-tighter file
  #    (ratio < 1/1.3) under-subtracts in those tools just as an axial-wider one
  #    over-subtracts.
  aniso_line <- if (isTRUE(b$has_z) && is.finite(b$sigma_xy) && is.finite(b$sigma_z))
    row(dim_(sprintf("lateral %.2f  \u00b7  axial %.2f nm", b$sigma_xy, b$sigma_z)),
        if (is.finite(b$sig_aniso) && (b$sig_aniso > 1.3 || b$sig_aniso < 1 / 1.3))
          warn(htmltools::HTML(sprintf(
            "  &#9888; axial is %.1f\u00d7 %s. Enter %.2f nm (lateral) for the D-map and 2D JDD.",
            if (b$sig_aniso > 1) b$sig_aniso else 1 / b$sig_aniso,
            if (b$sig_aniso > 1) "wider" else "tighter",
            b$sigma_xy))) else NULL)
  else NULL


  # ── data quality: same integrated-fraction metric as TSA / cycle-time status ──
  fi_raw  <- b$frac_integ %||% NA_real_;      fi_filt <- b$frac_integ_filt %||% NA_real_
  c_raw   <- b$cont_raw   %||% NA_real_;      c_filt  <- b$cont_filt      %||% NA_real_
  pair_txt <- function(a, z, fmt) if (is.finite(a) && is.finite(z)) sprintf(paste0(fmt, "\u2192", fmt), a, z)
                                  else if (is.finite(a)) sprintf(fmt, a) else "-"
  # D in the MOTION's degrees of freedom: a per-axis CVE spreads 2-DOF motion over
  # 3 recorded axes, so surface motion tracked in 3D reads 2/3 of the truth.
  d_obs <- if (isTRUE(b$has_z)) 3L else 2L
  d_mot <- { mo <- input$trk_motion %||% "volume"
             if (identical(mo, "filament")) 1L else if (identical(mo, "surface")) min(2L, d_obs) else d_obs }
  d_use <- b$D_cve * d_obs / d_mot
  # Descriptive only: the single-cycle fraction mostly reflects the sequence (photon
  # limit against photons per pass) and the dye's brightness, so it carries no verdict.
  qual_line <- row(
    dim_(sprintf("single-cycle steps %s  \u00b7  continuity %s%s",
                 pair_txt(100 * fi_raw, 100 * fi_filt, "%.0f%%"),
                 pair_txt(c_raw, c_filt, "%.2f"),
                 if (is.finite(c_raw) && is.finite(c_filt) && (c_filt - c_raw) > 0.01) " \u25b2" else "")))
  # Where a came from: a calibration profile, the ideal-estimator fallback
  # (uncalibrated), or a value typed by hand.
  prof_line <- if (is.finite(b$sig_phot_med %||% NA_real_)) {
    a_used <- b$a_pref %||% NA_real_
    ideal  <- .tcp_geom(b$tcp_pattern %||% TCP_PATTERN_DEFAULT)$a_ideal
    prof   <- trk_cal_profiles()[[input$trk_cal_profile %||% ""]]
    if (!is.null(prof) && isTRUE(abs(a_used - prof$a) < 1e-3))
      row(dim_(sprintf("calibrated: %s  \u00b7  a = %.2f (quartiles %.2f\u2013%.2f over %d cells)",
                       prof$name, prof$a, prof$a_q25 %||% NA_real_, prof$a_q75 %||% NA_real_,
                       as.integer(prof$n_cells %||% NA))))
    else if (isTRUE(abs(a_used - ideal) < 1e-9))
      row(dim_(sprintf(paste0("uncalibrated: a = %.2f is the ideal-estimator value for this pattern; ",
                              "measure it for this sequence under Instrument Calibration"), a_used)))
    else
      row(dim_(sprintf("a = %.2f entered by hand  \u00b7  at the ideal-estimator value %.2f, %s would be %.2f nm",
                       a_used, ideal, "\u03c3", b$sig_phot_med * a_used / ideal)))
  } else NULL

  # The OLSF fits the first few lags of the window, so the label names the lags
  # (and lag times) actually fitted, then the window they were taken from.
  fit_txt <- if (length(b$msd_fit_lags) == 2L && length(b$msd_fit_tau) == 2L)
    sprintf("lags %d\u2013%d (%.2g\u2013%.2g ms) of window %d\u2013%d",
            as.integer(b$msd_fit_lags[1]), as.integer(b$msd_fit_lags[2]),
            b$msd_fit_tau[1], b$msd_fit_tau[2],
            as.integer(b$msd_window[1]), as.integer(b$msd_window[2]))
  else sprintf("window %d\u2013%d", as.integer(b$msd_window[1]), as.integer(b$msd_window[2]))
  headline_line <- if (is.finite(b$D_msd %||% NA_real_))
    row(htmltools::tags$b(sprintf("D %.3g \u00b5m\u00b2/s", b$D_msd / 1e3)),
        dim_(sprintf("  MSD, %s  \u00b7  lag-1 check (CVE) %.3g%s",
                     fit_txt, d_use / 1e3,
                     if (is.finite(d_use) && b$D_msd > 0) sprintf(" (\u00d7%.2f)", d_use / b$D_msd) else "")))
  else row(warn("D: too few lags with data in the MSD fit window; only the lag-1 CVE value is available"))

  calib_line <- row(dim_(sprintf(
    "cycle %s  \u00b7  CVE D %s \u00b5m\u00b2/s%s  \u00b7  %s steps",
    num(b$t_cycle, "%.4g ms", "unknown"), num(d_use / 1e3, "%.3g"),
    if (d_mot != d_obs) sprintf(" (%d/%d for %s)", d_obs, d_mot, input$trk_motion %||% "volume") else "",
    format(b$n_steps, big.mark = ","))))

  # Tracker wander floor: the apparent diffusion the loop produces on an immobile
  # target. Reported beside D, never subtracted silently - D_w scales with the
  # photon budget, so a single constant would over-correct dim tracks.
  dw_line <- {
    dw <- suppressWarnings(as.numeric(input$trk_dw_floor %||% 0))
    d_head <- if (is.finite(b$D_msd %||% NA_real_)) b$D_msd else d_use   # the headline D
    prof <- trk_cal_profiles()[[input$trk_cal_profile %||% ""]]
    win_note <- if (!is.null(prof) && length(prof$window) == 2L &&
                    !isTRUE(all(as.integer(prof$window) == as.integer(b$msd_window))))
      sprintf(" (floor measured over cycles %d\u2013%d, D over %d\u2013%d)", as.integer(prof$window[1]),
              as.integer(prof$window[2]), as.integer(b$msd_window[1]), as.integer(b$msd_window[2])) else ""
    if (!is.finite(dw) || dw <= 0 || !is.finite(d_head)) NULL else {
      frac <- dw / d_head
      row(dim_(sprintf("wander floor %.3g µm²/s  →  D − Dw = %.3g%s",
                       dw / 1e3, max(0, d_head - dw) / 1e3, win_note)),
          if (frac > 0.5) warn(htmltools::HTML(sprintf(
            "  &#9888; floor is %.0f%% of D (not resolved above the tracker)", 100 * frac)))
          else if (frac > 0.2) warn(sprintf("  (%.0f%% of D)", 100 * frac))
          else dim_(sprintf("  (%.0f%% of D)", 100 * frac)))
    }
  }

  # ── linking QC: implausible jumps and per-track sigma outliers ──────────────
  qc_line <- if (is.finite(b$frac_implaus %||% NA_real_) || is.finite(b$n_suspect %||% NA_real_)) {
    bad <- (is.finite(b$frac_implaus) && b$frac_implaus > 0.01) ||
           (is.finite(b$n_suspect) && is.finite(b$n_tracks_pt) && b$n_tracks_pt > 0 &&
            b$n_suspect / b$n_tracks_pt > 0.05)
    txt <- paste(c(
      if (is.finite(b$n_implaus)) sprintf("%s implausible jumps (%.2f%%)", format(b$n_implaus, big.mark = ","), 100 * b$frac_implaus),
      if (is.finite(b$n_suspect)) sprintf("%d/%d suspicious tracks", b$n_suspect, b$n_tracks_pt)), collapse = "  \u00b7  ")
    row(htmltools::span(style = sprintf("color:%s;", if (bad) AMB else DIM),
      htmltools::tags$b("Linking QC "), txt, if (bad) ". Check for mis-links at high density or across long gaps." else ""))
  } else NULL

  htmltools::div(style = "font-size:0.92em; margin:4px 0; line-height:1.45;",
    headline_line, sigma_line, prof_line, blur_line, cve_line, aniso_line, qual_line, calib_line,
    dw_line, qc_line)
})

# Per-track σ_loc heterogeneity histogram (σ_loc is NOT a single instrument
# constant - this shows the spread; a tight peak justifies one global value).
# Distribution of the sigma_loc the app actually uses: one value per localization,
# sigma_i = L/(a*sqrt(N_i)). Plotting this rather than the per-track CVE sigma keeps
# the curve and the headline number the same quantity - the dashed median IS the
# figure reported above, so the two can never disagree.
output$trk_sigma_pertrack_plot <- renderPlotly({
  b <- tryCatch(trk_cve_bundle(), error = function(e) NULL)
  sg <- if (!is.null(b)) b$sig_phot else NULL
  sg <- sg[is.finite(sg) & sg > 0]
  shiny::validate(shiny::need(length(sg) >= 20,
    "Needs an ECO column: sigma_loc is computed per localization from the photon count."))
  dark <- isTRUE(isolate(input$trk_dark_mode))
  bg <- if (dark) "#000000" else "#FFFFFF"; fg <- if (dark) "#FFFFFF" else "#333333"
  med <- stats::median(sg); q <- stats::quantile(sg, c(.25, .75))
  plot_ly() %>%
    add_histogram(x = sg, nbinsx = 60,
                  marker = list(color = "rgba(120,160,220,0.7)",
                                line = list(color = "rgba(80,120,180,0.9)", width = 0.4)),
                  hovertemplate = "%{x:.2f} nm<br>%{y} localizations<extra></extra>",
                  name = "sigma") %>%
    layout(
      xaxis = list(title = "σ<sub>loc</sub> per localization (nm)", color = fg, zeroline = FALSE),
      yaxis = list(title = "localizations", color = fg),
      shapes = list(
        list(type = "rect", x0 = q[1], x1 = q[2], y0 = 0, y1 = 1, yref = "paper",
             fillcolor = "rgba(120,160,220,0.13)", line = list(width = 0), layer = "below"),
        list(type = "line", x0 = med, x1 = med, y0 = 0, y1 = 1, yref = "paper",
             line = list(color = "#c03030", width = 2, dash = "dash"))),
      annotations = list(list(
        x = med, y = 1, yref = "paper", yanchor = "bottom", showarrow = FALSE,
        text = sprintf("median %.2f nm  ·  IQR %.2f–%.2f", med, q[1], q[2]),
        font = list(color = "#c03030", size = 11))),
      paper_bgcolor = bg, plot_bgcolor = bg, font = list(color = fg),
      showlegend = FALSE, bargap = 0.04, uirevision = "sigma_pl")
})

# Transparency funnel: how many localizations / tracks survive each pipeline
# stage. Analyzing well-tracked molecules is a legitimate, DISCLOSED choice -
# this box is that disclosure. Reads debounced stages so it doesn't thrash.
output$trk_quality_funnel <- renderUI({
  cnt <- function(d) if (is.null(d) || !nrow(d) || !"tid" %in% names(d)) c(0, 0)
                     else c(nrow(d), length(unique(d$tid)))
  rs  <- trk_raw_stats_rv()
  raw <- if (!is.null(rs) && length(rs$tid_length)) c(sum(rs$tid_length), length(rs$tid_length)) else c(NA, NA)
  stg <- list(
    list("Gap-split",         cnt(tryCatch(trk_gap_split_d(),          error = function(e) NULL)), isTRUE(input$trk_split_gaps_enable)),
    list("End-trim",          cnt(tryCatch(trk_end_trimmed_d(),        error = function(e) NULL)), isTRUE(input$trk_endtrim_enable)),
    list("Quality gate",      cnt(tryCatch(trk_quality_filtered_d(),   error = function(e) NULL)), isTRUE(input$trk_quality_enable)),
    list("Metadata filter",   cnt(tryCatch(trk_meta_filtered_d(),      error = function(e) NULL)), TRUE),
    list("Length filter",     cnt(tryCatch(trk_len_filtered_d(),       error = function(e) NULL)), TRUE),
    list("Duration filter",   cnt(tryCatch(trk_dur_filtered_d(),       error = function(e) NULL)), TRUE),
    list("Analyzed (spatial)",cnt(tryCatch(trk_spatial_filtered(),     error = function(e) NULL)), TRUE)
  )
  base <- if (is.finite(raw[1]) && raw[1] > 0) raw[1] else stg[[1]][[2]][1]
  fmt  <- function(v) if (is.finite(v)) format(round(v), big.mark = ",") else "-"
  row  <- function(label, loc, tid, active, pct_base) {
    col <- if (isTRUE(active)) "#333" else "#aaa"
    htmltools::tags$tr(
      htmltools::tags$td(label, style = sprintf("padding:1px 10px 1px 0; color:%s;", col)),
      htmltools::tags$td(fmt(loc), style = sprintf("padding:1px 10px; text-align:right; color:%s;", col)),
      htmltools::tags$td(if (is.finite(pct_base) && pct_base > 0 && is.finite(loc))
                           sprintf("%.0f%%", 100 * loc / pct_base) else "",
                         style = sprintf("padding:1px 10px; text-align:right; color:%s;", col)),
      htmltools::tags$td(fmt(tid), style = sprintf("padding:1px 10px; text-align:right; color:%s;", col)),
      htmltools::tags$td(if (!isTRUE(active)) "(off)" else "",
                         style = "padding:1px 6px; color:#aaa; font-size:0.9em;"))
  }
  rows <- list(if (is.finite(raw[1])) row("Uploaded (raw)", raw[1], raw[2], TRUE, base) else NULL)
  for (s in stg) rows <- c(rows, list(row(s[[1]], s[[2]][1], s[[2]][2], s[[3]], base)))
  htmltools::div(
    htmltools::tags$table(style = "font-size:0.83em; border-collapse:collapse;",
      htmltools::tags$thead(htmltools::tags$tr(lapply(
        c("Stage", "localizations", "% kept", "tracks", ""),
        function(h) htmltools::tags$th(h, style = "padding:1px 10px; text-align:right; border-bottom:1px solid #ccc; color:#666;")))),
      htmltools::tags$tbody(rows)),
    htmltools::tags$p(style = "font-size:0.8em; color:#888; margin:5px 0 0 0;",
      "Analyses run on the last row. Filtering for well-tracked molecules is a disclosed choice - ",
      "report this funnel in your methods. Track count can rise at Gap-split (segments) then fall."))
})

# GOAL: show whether track continuity follows mobility in the analysed data, i.e.
# whether keeping "well-tracked" molecules means keeping slow ones.
# A molecule that moves fast relative to the TCP is lost for a cycle more often, so
# its track collects dark cycles; a continuity threshold can then remove fast
# molecules under the name of quality. The table ranks the analysed tracks by the
# same continuity score the gate uses (same half-life) and gives each quarter its
# own pooled CVE D. It is computed after "Estimate σ_loc & quality", on the data
# the analyses see, so it describes the population actually being reported.
# The ±20% band for the verdict is a reading aid, not a significance test.
output$trk_continuity_mobility <- renderUI({
  b <- tryCatch(trk_cve_bundle(), error = function(e) NULL)
  if (is.null(b))
    return(htmltools::tags$p(style = "font-size:0.83em; color:#888;",
      "Run Estimate σ_loc & quality first."))
  df  <- tryCatch(trk(), error = function(e) NULL)
  t_c <- tryCatch(trk_cycle_time(), error = function(e) NA_real_)
  H   <- suppressWarnings(as.numeric(input$trk_quality_darktol %||% 3))
  sc  <- .continuity_scores(df, t_c, H)
  tab <- continuity_mobility_table(df, sc, t_c, has_z = .has_real_z(df))
  if (is.null(tab))
    return(htmltools::tags$p(style = "font-size:0.83em; color:#888;",
      "Too few tracks to split into quarters (at least 40 needed)."))

  # Two readings, kept apart. r_mol: the gradient across molecules (best quarter vs
  # the other three), which is what a per-track statistic sees. r_pool: what keeping
  # only the best quarter would do to the POOLED D, which weights tracks by their
  # steps; when that quarter holds most localizations, r_mol can be large while the
  # reported number barely moves.
  r_mol  <- attr(tab, "ratio_best_rest")
  r_pool <- attr(tab, "ratio_best_all")
  share4 <- tab$loc_share[4]
  coupled <- is.finite(r_mol) && abs(log(r_mol)) > log(1.2)
  verdict <- if (!is.finite(r_mol)) NULL
    else if (!coupled)
      htmltools::tags$p(style = "font-size:0.85em; color:#1a8a1a; margin:4px 0;",
        sprintf("D is within 20%% across continuity (most continuous quarter %.2f\u00d7 the rest): the continuity gate is not selecting on speed.", r_mol))
    else
      htmltools::tags$div(style = "font-size:0.85em; margin:4px 0;",
        htmltools::tags$p(style = "color:#B8860B; margin:0 0 3px 0;",
          htmltools::HTML(sprintf(
            "&#9888; <b>Across molecules</b>, the most continuous quarter is %.0f%% %s than the rest (D %.3g vs %.3g \u00b5m\u00b2/s): per-track statistics depend on how strictly tracks are gated.",
            100 * abs(r_mol - 1), if (r_mol < 1) "slower" else "faster",
            tab$D_um2_s[4], attr(tab, "D_rest_um2_s")))),
        htmltools::tags$p(style = "color:#555; margin:0;",
          htmltools::HTML(sprintf(
            "<b>On the pooled D</b>, that quarter holds %.0f%% of the localizations, so keeping only it moves D from %.3g to %.3g \u00b5m\u00b2/s (%+.0f%%).",
            100 * share4, attr(tab, "D_all_um2_s"), tab$D_um2_s[4], 100 * (r_pool - 1)))))

  fmt <- function(v, d) ifelse(is.finite(v), formatC(v, format = "f", digits = d), "-")
  rows <- lapply(seq_len(nrow(tab)), function(i) htmltools::tags$tr(
    lapply(c(tab$quartile[i], tab$n_tracks[i],
             sprintf("%s-%s", fmt(tab$score_lo[i], 2), fmt(tab$score_hi[i], 2)),
             paste0(fmt(100 * tab$loc_share[i], 1), "%"),
             paste0(fmt(100 * tab$integ_frac[i], 0), "%"),
             fmt(tab$D_um2_s[i], 3), fmt(tab$sigma_nm[i], 2)),
           function(v) htmltools::tags$td(v, style = "padding:1px 8px; text-align:right;"))))
  htmltools::tagList(
    verdict,
    htmltools::tags$table(style = "font-size:0.83em; border-collapse:collapse;",
      htmltools::tags$thead(htmltools::tags$tr(lapply(
        c("Continuity", "tracks", "score", "of locs", "integrated steps", "D (µm²/s)", "σ (nm)"),
        function(h) htmltools::tags$th(h, style = "padding:1px 8px; text-align:right; border-bottom:1px solid #ccc; color:#666;")))),
      htmltools::tags$tbody(rows)),
    htmltools::tags$p(style = "font-size:0.8em; color:#888; margin:5px 0 0 0;",
      sprintf("Four equal groups by continuity score (half-life %.3g cycles); D is the pooled CVE of each.", H)))
})
