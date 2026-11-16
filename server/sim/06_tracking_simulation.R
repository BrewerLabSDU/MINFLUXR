# ══════════════════════════════════════════════════════════════════════════
# Diffusion / Tracking ground-truth simulator (Phase 1)
# See docs/simulation_module.md. Free surface Brownian on a sphere or rod,
# sampled the MINFLUX way (cycle clock + dark cycles + heteroscedastic σ_loc +
# photobleaching) → MINFLUX-format tracks + a ground-truth D. The CVE round-trip
# in the summary demonstrates the validation loop in-place.
# ══════════════════════════════════════════════════════════════════════════
# --- Calibrate the sampling model from a real MINFLUX file -----------------
# Reads a real .csv/.npy/.json and auto-fills cycle time, photons, dark-cycle
# probability and mean track length so the simulator matches a real experiment.
sim_trk_calib_rv <- reactiveVal(NULL)
observeEvent(input$sim_trk_calib, {
  fi <- input$sim_trk_calib_file
  shiny::validate(shiny::need(!is.null(fi), "Choose a MINFLUX file first."))
  p <- fi$datapath; nm <- tolower(fi$name)
  df <- withProgress(message = "Reading MINFLUX file…", value = 0, {
    tryCatch({
      if (grepl("\\.npy$",  nm)) read_minflux_npy_as_df(p, keep_final_only = TRUE, keep_valid_only = TRUE)
      else if (grepl("\\.json$", nm)) read_minflux_json_as_df(p, keep_final_only = TRUE, keep_valid_only = TRUE)
      else if (grepl("\\.csv$",  nm)) {
        x <- as.data.frame(data.table::fread(p, showProgress = FALSE))
        if (all(c("tid","tim") %in% names(x))) x else canonicalize_minflux_df(x)
      } else NULL
    }, error = function(e) NULL)
  })
  shiny::validate(shiny::need(!is.null(df) && nrow(df) > 0 && all(c("tid","tim") %in% names(df)),
    "Could not read tid / tim from that file."))
  df <- df[is.finite(df$tim), ]
  # Calibrate ONLY on analysable traces (≥ MIN_LEN localizations, no upper bound), so
  # the short-junk traces don't drag the statistics down. The representative length is
  # the MEDIAN of the kept traces.
  MIN_LEN  <- max(2L, as.integer(input$sim_trk_calib_minlen %||% 100L))
  len_all  <- table(df$tid)
  keep_tid <- names(len_all)[len_all >= MIN_LEN]
  shiny::validate(shiny::need(length(keep_tid) >= 5, sprintf("Fewer than 5 traces have ≥%d localizations.", MIN_LEN)))
  n_total  <- length(len_all)
  df <- df[as.character(df$tid) %in% keep_tid, ]        # keep only the long traces
  len     <- as.numeric(table(df$tid)); med_len <- stats::median(len)
  photons <- if ("eco" %in% names(df)) stats::median(df$eco[is.finite(df$eco) & df$eco > 0]) else NA_real_
  o    <- order(df$tid, df$tim); df <- df[o, ]
  tid <- df$tid; tim <- df$tim; nrw <- length(tid)
  same <- c(FALSE, tid[-1] == tid[-nrw]); dt <- c(NA_real_, diff(tim)) * 1000; dt[!same] <- NA_real_  # ms
  dtv  <- dt[is.finite(dt) & dt > 0]
  shiny::validate(shiny::need(length(dtv) > 20, "Too few within-track steps to calibrate."))
  # ── Cycle time and SW from the Δt quantisation ────────────────────────────
  # Real Δt obey Δt = k·t_c + SW exactly, so the MODAL Δt is t_c + SW, not t_c.
  # Take the mode as a first guess, assign each Δt its cycle count k, then read
  # t_c off the SPACING between the levels and SW off the intercept. Taking the
  # mode for the cycle time inflates it by SW (~5% on the tracking test files).
  hb    <- hist(dtv[dtv < stats::quantile(dtv, 0.98)], breaks = 400, plot = FALSE)
  modal <- hb$mids[which.max(hb$counts)]                         # modal Δt = t_c + SW
  k0    <- round(dtv / modal)
  ok    <- k0 >= 1 & k0 <= 10 & abs(dtv - k0 * modal) < 0.3 * modal   # on-level Δt only
  cyc <- modal; sw <- 0
  if (sum(ok) > 20) {
    cf <- stats::coef(stats::lm(dtv[ok] ~ k0[ok]))
    if (is.finite(cf[2]) && cf[2] > 0 && cf[1] >= 0 && cf[1] < cf[2]) {
      cyc <- unname(cf[2]); sw <- unname(cf[1])
    }
  }
  ncyc <- pmax(1L, as.integer(round((dtv - sw) / cyc)))
  p_dark <- sum(ncyc - 1L) / sum(ncyc); integ <- mean(ncyc == 1L)
  # ── D and σ_loc from the single-cycle steps (CVE, blur-aware) ─────────────
  # σ² = R⟨Δx²⟩ + (2R−1)·cov₁ per axis with R = 1/6 (full-cycle integration);
  # D itself is blur-independent. Lateral axes only: real 2D tracking files
  # carry a placeholder z.
  dstep <- same & is.finite(dt) & abs(round((dt - sw) / cyc) - 1) < 1e-6
  dstep[is.na(dstep)] <- FALSE
  dx <- c(NA_real_, diff(df$x)); dy <- c(NA_real_, diff(df$y)); dx[!same] <- NA_real_; dy[!same] <- NA_real_
  pair <- dstep & c(dstep[-1], FALSE)
  cov1 <- mean(c(dx[pair] * c(dx[-1], NA_real_)[pair],
                 dy[pair] * c(dy[-1], NA_real_)[pair]), na.rm = TRUE)
  msq  <- mean(c(dx[dstep]^2, dy[dstep]^2), na.rm = TRUE)
  D_meas <- if (is.finite(msq) && is.finite(cov1)) (msq / (2 * (cyc + sw)) + cov1 / (cyc + sw)) / 1000 else NA_real_
  sig_meas <- if (is.finite(msq) && is.finite(cov1)) sqrt(max(0, msq / 6 - (2/3) * cov1)) else NA_real_
  # -- Photon budget: collection window, photon limit, repeats ---------------
  # tau = eco/efo is the collection window, written in the file per localization and
  # quantised in units of one pattern repeat. It gives the blur coefficient and the
  # repeat structure without any assumption.
  pht_lim <- NA_real_; expf <- NA_real_; nrep <- NA_real_
  if (all(c("eco","efo") %in% names(df))) {
    good <- is.finite(df$eco) & is.finite(df$efo) & df$efo > 0 & df$eco > 0
    if (sum(good) > 100) {
      tau_v   <- df$eco[good] / df$efo[good] * 1000                  # ms
      tau0    <- unname(stats::quantile(tau_v, 0.05))                # one pattern repeat
      pht_lim <- unname(stats::quantile(df$eco[good], 0.001))        # hard floor = phtLimit
      expf    <- tau0 / cyc
      nrep    <- unname(stats::quantile(tau_v / tau0, 0.999))
    }
  }
  # -- Brightness: between molecules, within a track, correlation, fade ------
  b_btw <- NA_real_; b_wth <- NA_real_; b_phi <- NA_real_; b_dec <- NA_real_
  if ("efo" %in% names(df)) {
    dd  <- df[is.finite(df$efo) & df$efo > 0, ]
    spl <- split(log(dd$efo), dd$tid)
    spl <- spl[vapply(spl, length, 1L) >= 20L]
    if (length(spl) >= 5L) {
      b_btw <- stats::sd(vapply(spl, stats::median, 0))              # molecule to molecule
      b_wth <- stats::median(vapply(spl, stats::sd, 0))              # inside one track
      b_phi <- stats::median(vapply(spl, function(v)
                 suppressWarnings(stats::cor(v[-length(v)], v[-1])), 0), na.rm = TRUE)
      # Fade: log drop from the first tenth of a track to the last, POOLED over
      # localizations rather than averaged over tracks. The per-track median measures
      # something else - short tracks mostly end by being lost rather than by
      # bleaching, so they barely fade - and calibrating on it under-states the fade
      # by a factor of about seven on this data.
      dec <- unlist(lapply(spl, function(v) {
               f <- (seq_along(v) - 1) / max(1, length(v) - 1); pmin(9L, floor(f * 10)) }))
      vv  <- unlist(spl, use.names = FALSE)
      b_dec <- stats::median(vv[dec == 0L]) - stats::median(vv[dec == 9L])
    }
  }
  # -- Track-length heterogeneity --------------------------------------------
  # The generator compounds a geometric lifetime with a lognormal mean; mean/median of
  # the resulting lengths is a monotone function of that spread, inverted here.
  ratio  <- mean(len) / stats::median(len)
  lifecv <- stats::approx(c(1.45, 1.68, 2.14, 2.82, 3.60, 5.01), c(0, 0.5, 1, 1.5, 2, 3),
                          xout = ratio, rule = 2)$y
  # -- D-lifetime coupling ----------------------------------------------------
  # Per-track D against track length. Both are noisy, so the correlation the simulator
  # reproduces is weaker than the one measured; the value filled in is de-attenuated.
  rho_meas <- NA_real_; rho_set <- NA_real_
  if (all(c("x","y") %in% names(df))) {
    pt <- data.frame(tid = df$tid, dx = dx, dy = dy,
                     st = dstep, pr = dstep & c(dstep[-1], FALSE))
    sp <- split(pt, pt$tid)
    dd <- vapply(sp, function(z) {
            ii <- which(z$st); jj <- which(z$pr)
            if (length(ii) < 30 || length(jj) < 10) return(NA_real_)
            (mean(c(z$dx[ii]^2, z$dy[ii]^2), na.rm = TRUE) / (2*(cyc+sw)) +
             mean(c(z$dx[jj]*z$dx[jj+1], z$dy[jj]*z$dy[jj+1]), na.rm = TRUE) / (cyc+sw)) / 1000
          }, 0)
    nn <- vapply(sp, nrow, 1L)
    gd <- is.finite(dd) & dd > 0
    if (sum(gd) >= 10) {
      rho_meas <- suppressWarnings(stats::cor(log(nn[gd]), log(dd[gd]), method = "spearman"))
      rho_set  <- max(-0.95, min(0.95, 2 * rho_meas))
    }
  }
  setn <- function(id, v, dig = 2) if (is.finite(v)) updateNumericInput(session, id, value = round(v, dig))
  if (is.finite(cyc))     updateNumericInput(session, "sim_trk_cycle",   value = signif(cyc, 3))
  if (is.finite(sw))      updateNumericInput(session, "sim_trk_sw",      value = signif(sw, 2))
  if (is.finite(p_dark))  updateNumericInput(session, "sim_trk_pdark",   value = round(p_dark, 2))
  if (is.finite(photons)) updateNumericInput(session, "sim_trk_photons", value = round(photons))
  if (is.finite(med_len)) updateNumericInput(session, "sim_trk_meanlen", value = round(med_len))
  if (is.finite(D_meas) && D_meas > 0)
    updateNumericInput(session, "sim_trk_D", value = signif(D_meas, 2))
  if (length(keep_tid) >= 1) updateNumericInput(session, "sim_trk_ntracks", value = length(keep_tid))
  setn("sim_trk_phtlimit", pht_lim, 0); setn("sim_trk_expfrac", expf)
  setn("sim_trk_taurep",   nrep, 0)
  setn("sim_trk_phtbtw",   b_btw); setn("sim_trk_phtcv", b_wth)
  setn("sim_trk_phtphi",   max(0, b_phi)); setn("sim_trk_phtdecay", max(0, b_dec))
  setn("sim_trk_lifecv",   lifecv); setn("sim_trk_Dlife", rho_set)
  sim_trk_calib_rv(list(file = fi$name, n_tracks = length(keep_tid), n_total = n_total,
    min_len = MIN_LEN, n_loc = nrow(df), med_len = med_len, mean_len = mean(len),
    photons = photons, cyc = cyc, sw = sw, modal = modal, p_dark = p_dark, integ = integ,
    D = D_meas, sigma = sig_meas, med_dt = stats::median(dtv),
    pht_lim = pht_lim, expf = expf, nrep = nrep, b_btw = b_btw, b_wth = b_wth,
    b_phi = b_phi, b_dec = b_dec, lifecv = lifecv, rho_meas = rho_meas, rho_set = rho_set))
  showNotification(sprintf("Calibrated on %d traces (≥%d locs) of %d total.",
                           length(keep_tid), MIN_LEN, n_total), type = "message", duration = 6)
})

output$sim_trk_calib_report <- renderUI({
  cb <- sim_trk_calib_rv(); if (is.null(cb)) return(NULL)
  rowf <- function(k, v) htmltools::tags$tr(
    htmltools::tags$td(style = "padding:1px 10px 1px 0; color:#666;", k),
    htmltools::tags$td(htmltools::tags$b(v)))
  htmltools::div(style = "font-size:0.8em; color:#333; margin-top:4px;",
    htmltools::tags$div(style = "color:#228833;", sprintf("✓ Calibrated from %s", cb$file)),
    htmltools::tags$table(style = "border-collapse:collapse;",
      rowf(sprintf("Traces used (≥%d locs)", cb$min_len %||% 50L),
           sprintf("%s of %s  ·  %s locs", format(cb$n_tracks, big.mark = ","),
                   format(cb$n_total %||% cb$n_tracks, big.mark = ","), format(cb$n_loc, big.mark = ","))),
      rowf("Cycle time", sprintf("%.4g ms  (SW %.3g ms · modal Δt %.4g)", cb$cyc, cb$sw %||% 0, cb$modal %||% NA)),
      rowf("Photons/loc (median eco)", if (is.finite(cb$photons)) sprintf("%.0f", cb$photons) else "-"),
      rowf("Dark-cycle prob.", sprintf("%.2f  (%.0f%% integrated)", cb$p_dark, 100 * cb$integ)),
      rowf("Track length (median)", sprintf("%.0f · mean %.0f", cb$med_len, cb$mean_len)),
      rowf("D (CVE, lateral)", if (is.finite(cb$D %||% NA)) sprintf("%.3g µm²/s", cb$D) else "-"),
      rowf("σ_loc (CVE, measured R)", if (is.finite(cb$sigma %||% NA)) sprintf("%.2f nm", cb$sigma) else "-"),
      rowf("Collection window", if (is.finite(cb$expf %||% NA))
             sprintf("%.0f%% of the cycle · up to %.0f repeats · photon limit %.0f",
                     100*cb$expf, cb$nrep, cb$pht_lim) else "-"),
      rowf("Brightness", if (is.finite(cb$b_btw %||% NA))
             sprintf("spread %.2f between / %.2f within · correlation %.2f · fades to %.0f%%",
                     cb$b_btw, cb$b_wth, cb$b_phi, 100*exp(-cb$b_dec)) else "-"),
      rowf("Length spread", sprintf("CV %.1f  (mean/median = %.2f)", cb$lifecv, cb$mean_len/cb$med_len)),
      rowf("D vs track length", if (is.finite(cb$rho_meas %||% NA))
             sprintf("%+.2f measured → %+.2f set (noise attenuates it)", cb$rho_meas, cb$rho_set) else "-")),
    htmltools::tags$p(style = "color:#888; margin-top:3px;",
      "Every acquisition field is filled from this file. Set the sample geometry, the motion model ",
      "and the TCP diameter, then Simulate."))
})

sim_trk_result <- eventReactive(input$sim_trk_run, {
  set.seed(as.integer(input$sim_trk_seed %||% 1L))
  geom     <- input$sim_trk_geom %||% "rod"
  mode     <- input$sim_trk_mode %||% "surface"      # surface | volume (intracellular)
  mthick   <- max(0, as.numeric(input$sim_trk_mthick %||% 6))              # membrane thickness (nm)
  # Membrane tilt: the sheet is primarily in XY but inclined by mtilt°, so the plane
  # is spanned by e1,e2 (normal n̂ near ẑ). The path is generated flat (local u,v,w)
  # then rotated by R_mem = [e1 e2 n̂] → z varies across the plane. 0° = flat.
  mtilt <- max(0, min(89, as.numeric(input$sim_trk_mtilt %||% 15))) * pi / 180
  maz   <- runif(1, 0, 2 * pi)
  n_hat <- c(sin(mtilt) * cos(maz), sin(mtilt) * sin(maz), cos(mtilt))
  e1    <- { cx <- c(-n_hat[2], n_hat[1], 0); m <- sqrt(sum(cx^2)); if (m < 1e-9) c(1,0,0) else cx / m }  # ẑ×n̂ (in XY)
  e2    <- c(n_hat[2]*e1[3] - n_hat[3]*e1[2], n_hat[3]*e1[1] - n_hat[1]*e1[3], n_hat[1]*e1[2] - n_hat[2]*e1[1])  # n̂×e1
  R_mem <- cbind(e1, e2, n_hat)                                            # local (u,v,w) → world
  # motion DOF (used by the CVE check): free=3, membrane=2 (lateral), else surface=2/volume=3
  d_mot    <- if (identical(geom, "free")) 3L else if (identical(geom, "membrane")) 2L
              else if (identical(mode, "volume")) 3L else 2L
  D_um2_s  <- max(0, as.numeric(input$sim_trk_D %||% 0.05))
  D_nm2ms  <- D_um2_s * 1000                        # µm²/s → nm²/ms
  n_trk    <- max(1L, as.integer(input$sim_trk_ntracks %||% 200L))
  mean_len <- max(5L, as.integer(input$sim_trk_meanlen %||% 80L))
  cyc      <- max(0.01, as.numeric(input$sim_trk_cycle %||% 0.192))   # ms, spacing between cycles
  # Fixed per-localization overhead (settling / window). Real files obey
  # Delta-t = k*cycle + SW exactly, so the modal Delta-t is cycle + SW, not cycle.
  sw       <- max(0, as.numeric(input$sim_trk_sw %||% 0))              # ms
  p_dark   <- min(0.95, max(0, as.numeric(input$sim_trk_pdark %||% 0.28)))
  photons  <- max(1, as.numeric(input$sim_trk_photons %||% 38))
  L_tcp    <- max(1, as.numeric(input$sim_trk_Ltcp %||% 100.8))
  blur     <- isTRUE(input$sim_trk_blur)                               # integrate over the cycle
  # Photon collection window as a fraction of the cycle. tau together with the step
  # duration sets the blur coefficient R = tau/(6*dt): a full-cycle window gives 1/6,
  # a short burst tends to 0.
  expfrac  <- min(1, max(0.05, as.numeric(input$sim_trk_expfrac %||% 0.55)))
  # A localization followed by a dark gap integrates ACROSS that gap (repeated pattern
  # passes) instead of skipping it, up to this many repeats. Measured on real files as
  # tau = eco/efo: 1 repeat before a single-cycle step, 2 before a two-cycle step,
  # saturating near 3. 1 disables it.
  tau_rep  <- max(1L, as.integer(input$sim_trk_taurep %||% 1L))
  # Lognormal spread of D between molecules (0 = every track identical). The input D
  # stays the MEDIAN of the population.
  D_cv     <- max(0, as.numeric(input$sim_trk_Dcv %||% 0))
  D_sdlog  <- if (D_cv > 0) sqrt(log(1 + D_cv^2)) else 0
  # Re-acquisition shrinkage: the first localization after a dark gap is pulled this
  # fraction of the way back toward the previous reported position.
  shrink   <- min(1, max(0, as.numeric(input$sim_trk_shrink %||% 0)))
  # Photon model: the cycle ends once phtLimit photons are in, so eco has a hard
  # floor there; brightness differences between molecules give the long right tail.
  # Rate is lognormal (median = the photons input) and the count Poisson around it,
  # so sdlog = 0 reduces to plain Poisson.
  pht_min  <- max(1L, as.integer(input$sim_trk_phtlimit %||% 1L))
  pht_sd   <- max(0, as.numeric(input$sim_trk_phtcv %||% 0))       # within-track sdlog
  pht_btw  <- max(0, as.numeric(input$sim_trk_phtbtw %||% 0))      # between-molecule sdlog
  pht_phi  <- min(0.95, max(0, as.numeric(input$sim_trk_phtphi %||% 0)))  # AR(1) on log rate
  # Emission fades over a molecule's lifetime: real tracks end at ~44% of their
  # starting rate, and their gap fraction climbs from 18% to 43% as they fade.
  pht_dec  <- max(0, as.numeric(input$sim_trk_phtdecay %||% 0))
  # Photon budget: the tracker repeats the pattern until phtLimit photons are in, and
  # each extra repeat costs a cycle. Dark cycles then FOLLOW from brightness instead of
  # being drawn from p_dark, and eco gets its hard floor at the limit for free.
  budget   <- isTRUE(input$sim_trk_budget)
  # Spread of the bleaching lifetime between molecules. 0 = one geometric lifetime
  # for the whole population; above 0 the per-molecule mean lifetime is lognormal,
  # which makes the recorded-length distribution heavy-tailed like a real file.
  life_cv  <- max(0, as.numeric(input$sim_trk_lifecv %||% 0))
  life_sd  <- if (life_cv > 0) sqrt(log(1 + life_cv^2)) else 0
  # Correlation between a molecule's D and its lifetime. Real files anticorrelate at
  # about −0.4: the slow molecules are the ones tracked longest, which is exactly why
  # a step-weighted pooled D lands below the per-track median.
  D_life   <- min(0.95, max(-0.95, as.numeric(input$sim_trk_Dlife %||% 0)))
  proj2d   <- isTRUE(input$sim_trk_project2d)
  sd_step  <- sqrt(2 * D_nm2ms * cyc)               # per-axis step SD per cycle (nm)
  R    <- max(1, as.numeric(input$sim_trk_R %||% 500))
  rr   <- max(1, as.numeric(input$sim_trk_r %||% 400))
  Lcyl <- max(1, as.numeric(input$sim_trk_Lcyl %||% 2000))
  circ <- 2 * pi * rr
  ea   <- max(1, as.numeric(input$sim_trk_ea %||% 500))
  eb   <- max(1, as.numeric(input$sim_trk_eb %||% 500))
  ec   <- max(1, as.numeric(input$sim_trk_ec %||% 900))
  abc  <- if (identical(geom, "sphere")) c(R, R, R) else c(ea, eb, ec)   # ellipsoid semi-axes
  # Long-axis orientation for the rod: E. coli on a coverslip lie flat → axis in the
  # image plane. Generated in canonical coords (axis ∥ X) then permuted to the choice.
  rod_axis <- input$sim_trk_rod_axis %||% "x"
  rod_perm <- switch(rod_axis, x = c(1L,2L,3L), y = c(2L,1L,3L), z = c(2L,3L,1L))
  # ── FIELD OF VIEW ──────────────────────────────────────────────────────────
  # More than one cell: scatter them over the field at random positions and random
  # in-plane orientations, and give each track a cell. Placement is rejection-based
  # so cells do not sit on top of each other. A single cell stays at the origin, so
  # the estimator-validation view is unchanged.
  n_cells  <- max(1L, as.integer(input$sim_trk_ncells %||% 1L))
  fov_nm   <- max(1, as.numeric(input$sim_trk_fov %||% 10)) * 1000
  cell_len <- if (identical(geom, "rod")) Lcyl + 2 * rr else 2 * max(abc)
  cell_pos <- matrix(0, n_cells, 2); cell_ang <- numeric(n_cells)
  if (n_cells > 1L) {
    for (ci in seq_len(n_cells)) {
      for (try in 1:200) {
        cand <- c(runif(1, -fov_nm/2, fov_nm/2), runif(1, -fov_nm/2, fov_nm/2))
        if (ci == 1L) break
        dmin <- min(sqrt(rowSums((cell_pos[seq_len(ci-1), , drop = FALSE] -
                                  rep(cand, each = ci-1))^2)))
        if (dmin > cell_len) break
      }
      cell_pos[ci, ] <- cand; cell_ang[ci] <- runif(1, 0, 2*pi)
    }
  }

  # ── unified geometry operators (rod = spherocylinder body + hemispherical caps;
  #    else ellipsoid, sphere = a=b=c). project = nearest surface point; reflect =
  #    back inside the volume; normal = unit outward normal (for the tangent step). ─
  is_rod <- identical(geom, "rod")
  geo_normal <- function(p) {
    if (is_rod) { a <- p[1]
      n <- if (a >= 0 && a <= Lcyl) c(0, p[2], p[3]) else { cx <- if (a < 0) 0 else Lcyl; p - c(cx,0,0) }
      m <- sqrt(sum(n^2)); if (m == 0) c(0,1,0) else n / m
    } else { n <- p / abc^2; n / sqrt(sum(n^2)) }
  }
  geo_project <- function(p) {
    if (is_rod) { a <- p[1]
      if (a >= 0 && a <= Lcyl) { rad <- c(p[2],p[3]); rho <- sqrt(sum(rad^2)); rad <- if (rho == 0) c(rr,0) else rad*rr/rho; c(a, rad) }
      else { cx <- if (a < 0) 0 else Lcyl; d <- p - c(cx,0,0); dd <- sqrt(sum(d^2)); d <- if (dd == 0) c(if (a<0) -rr else rr,0,0) else d*rr/dd; c(cx,0,0)+d }
    } else { u <- p/abc; abc * u / sqrt(sum(u^2)) }
  }
  geo_reflect <- function(p) {
    if (is_rod) { a <- p[1]
      if (a >= 0 && a <= Lcyl) { rad <- c(p[2],p[3]); rho <- sqrt(sum(rad^2)); if (rho > rr) rad <- rad*(2*rr - rho)/rho; c(a, rad) }
      else { cx <- if (a < 0) 0 else Lcyl; d <- p - c(cx,0,0); dd <- sqrt(sum(d^2)); if (dd > rr) d <- d*(2*rr - dd)/dd; c(cx,0,0)+d }
    } else { u <- p/abc; nu <- sqrt(sum(u^2)); if (nu > 1) u <- u*(2-nu)/nu; abc*u }
  }
  geo_rand_surf <- function() if (is_rod) { s <- runif(1,0,Lcyl); phi <- runif(1,0,2*pi); c(s, rr*cos(phi), rr*sin(phi)) } else geo_project(rnorm(3) * abc)
  geo_rand_in   <- function() if (is_rod) {
      repeat { p <- c(runif(1,-rr,Lcyl+rr), runif(1,-rr,rr), runif(1,-rr,rr)); a <- p[1]
        if ((if (a >= 0 && a <= Lcyl) (p[2]^2+p[3]^2) <= rr^2 else if (a < 0) sum(p^2) <= rr^2 else sum((p-c(Lcyl,0,0))^2) <= rr^2)) return(p) }
    } else { d <- rnorm(3); d <- d/sqrt(sum(d^2)); abc * d * runif(1)^(1/3) }

  # ── motion modifier (Phase 3): free | confined (OU tether) | directed (drift) |
  #    twostate (mobile↔immobile CTMC). Applied to the proposed 3D point, then the
  #    geometry constraint re-projects/reflects it. ──────────────────────────────
  motion  <- input$sim_trk_motion  %||% "free"
  # Anomalous diffusion: MSD grows as t^alpha via correlated fBm increments. D still
  # sets the single-cycle MSD, so only the growth with lag changes.
  alpha   <- min(1.95, max(0.05, as.numeric(input$sim_trk_alpha %||% 1)))
  is_anom <- identical(motion, "anomalous")
  # Hop diffusion: free inside a corral, reflecting off its barrier, crossing into the
  # neighbour with probability p_perm per encounter. Free at short lags, slower at long
  # ones, so alpha < 1 in between without any long-range memory.
  L_hop   <- max(1, as.numeric(input$sim_trk_Lhop %||% 60))
  p_perm  <- min(1, max(0, as.numeric(input$sim_trk_pperm %||% 0.05)))
  L_c     <- max(1, as.numeric(input$sim_trk_Lc %||% 150))         # confinement radius (nm)
  v_nmms  <- max(0, as.numeric(input$sim_trk_vel %||% 0))          # directed speed (µm/s = nm/ms)
  kswitch <- min(1, max(0, as.numeric(input$sim_trk_kswitch %||% 0.03)))  # per-cycle prob, mobile -> slow
  # Slow-state D as a fraction of D, and the per-cycle prob of leaving the slow state.
  # Separate values make asymmetric dwells and a finite D contrast possible, which the
  # HMM recovery sweep needs (defaults reproduce the old symmetric mobile/immobile case).
  slow_frac <- min(1, max(0, as.numeric(input$sim_trk_slow_frac %||% 0.01)))
  kback     <- min(1, max(0, as.numeric(input$sim_trk_kswitch_back %||% kswitch)))

  drift0  <- if (is_rod) c(1,0,0)                                          # rod drift ∥ long axis
             else if (identical(geom, "membrane")) { a <- runif(1, 0, 2*pi); c(cos(a), sin(a), 0) }  # in-plane
             else { v <- rnorm(3); v/sqrt(sum(v^2)) }

  is_free <- identical(geom, "free"); is_mem <- identical(geom, "membrane")
  # dur[i] = real elapsed time of the step leaving cycle i (ms). The molecule keeps
  # diffusing through the SW overhead, so SW must enter the step duration and not
  # only the clock; otherwise D is under-recovered by cyc/(cyc+SW) (~5%).
  gen_path <- function(n_cyc, D_trk) {
    surface <- identical(mode, "surface")
    tau_c   <- if (D_trk > 0) L_c^2 / D_trk else Inf                # OU relaxation time (ms)
    p <- if (is_free) c(0,0,0)
         else if (is_mem) c(0, 0, if (mthick > 0) runif(1, -mthick/2, mthick/2) else 0)
         else if (surface) geo_rand_surf() else geo_rand_in()
    home <- p; st <- 1L
    pos <- matrix(NA_real_, n_cyc, 3); pos[1,] <- p; states <- integer(n_cyc); states[1] <- 1L
    # Anomalous motion draws its increments up-front: fBm increments are CORRELATED,
    # so they cannot be generated one step at a time like Brownian ones.
    inc <- if (is_anom)
             cbind(.sim_fgn(n_cyc, alpha/2), .sim_fgn(n_cyc, alpha/2), .sim_fgn(n_cyc, alpha/2))
           else NULL
    # Recorded/dark and the clock are decided inside the loop: dur[i] (the elapsed
    # time of the step leaving cycle i) needs to know whether cycle i was recorded,
    # and a displacement-dependent dark probability needs the step that arrived at
    # cycle i. Deciding them in the loop keeps both exact.
    # ── BRIGHTNESS (photons/ms) over the molecule's life ────────────────────
    # log-rate = median + between-molecule offset − fade over the lifetime + an
    # AR(1) wobble, so brightness is correlated from one localization to the next.
    tau0  <- expfrac * cyc                                  # one pattern repeat (ms)
    lrate <- log(max(1e-6, photons / tau0))
    if (pht_btw > 0) lrate <- lrate + stats::rnorm(1, 0, pht_btw)
    frac  <- if (n_cyc > 1L) (seq_len(n_cyc) - 1L) / (n_cyc - 1L) else 0
    lr    <- lrate - pht_dec * frac
    if (pht_sd > 0) {
      e <- numeric(n_cyc); z <- stats::rnorm(n_cyc)
      e[1] <- pht_sd * z[1]
      if (n_cyc > 1L) for (i in 2:n_cyc)
        e[i] <- pht_phi * e[i-1] + sqrt(1 - pht_phi^2) * pht_sd * z[i]
      lr <- lr + e
    }
    rate <- exp(lr)                                          # photons per ms

    integ <- logical(n_cyc); integ[1] <- TRUE
    dur   <- rep(cyc, n_cyc); dur[1] <- cyc + sw
    eco_i <- rep(NA_real_, n_cyc); rep_i <- rep(NA_integer_, n_cyc)
    if (budget) {
      # Repeat the pattern until the photon limit is reached, up to tau_rep repeats.
      # A window that never reaches the limit is abandoned: those cycles are dark and
      # the gap simply grows. This is what makes eco floor at the limit exactly.
      integ[] <- FALSE
      i <- 1L
      while (i <= n_cyc) {
        pend <- 0; j <- 0L
        while (j < tau_rep && i + j <= n_cyc) {
          j <- j + 1L
          pend <- pend + stats::rpois(1, rate[i + j - 1L] * tau0)
          if (pend >= pht_min) break
        }
        if (pend >= pht_min) { integ[i] <- TRUE; eco_i[i] <- pend; rep_i[i] <- j }
        i <- i + j
      }
      dur <- cyc + sw * integ
    }
    for (i in 2:n_cyc) {
      Dstep <- if (identical(motion,"twostate") && st == 0L) D_trk * slow_frac else D_trk
      dti   <- dur[i - 1L]                                                # this step's duration (ms)
      sdi   <- sqrt(2 * Dstep * dti)
      g0 <- if (is_anom) inc[i, ] * sdi else rnorm(3, 0, sdi)             # correlated when anomalous
      if (is_free)      { p2 <- p + g0 }                                  # unbounded 3D
      else if (is_mem)  { p2 <- p + c(g0[1:2], 0) }                        # 2D lateral (z handled below)
      else if (surface) { nn <- geo_normal(p); p2 <- p + (g0 - sum(g0*nn)*nn) }  # tangent step
      else              { p2 <- p + g0 }                                  # 3D volume
      if (identical(motion,"directed")) p2 <- p2 + drift0 * v_nmms * dti
      if (identical(motion,"confined")) p2 <- home + (p2 - home) * exp(-dti / tau_c)  # OU pull to home
      if (identical(motion,"hop")) {                       # corral with a permeable barrier
        dv <- p2 - home; nd <- sqrt(sum(dv^2))
        if (nd > L_hop) {
          if (runif(1) < p_perm) home <- home + (dv / nd) * (2 * L_hop)   # cross into the neighbour
          else                   p2   <- home + (dv / nd) * (2 * L_hop - nd)  # reflect off the barrier
        }
      }
      p  <- if (is_free) p2
            else if (is_mem) c(p2[1], p2[2],                              # keep z within the membrane thickness
                               if (mthick > 0) max(-mthick/2, min(mthick/2, p[3] + rnorm(1, 0, sdi * 0.2))) else 0)
            else if (surface) geo_project(p2) else geo_reflect(p2)
      if (identical(motion,"twostate") && runif(1) < (if (st == 1L) kswitch else kback)) st <- 1L - st
      pos[i,] <- p; states[i] <- st
      # Was this cycle recorded? Under the photon budget that was already decided from
      # brightness above; otherwise fall back to a flat per-cycle probability.
      if (!budget) {
        integ[i] <- runif(1) > p_dark
        dur[i]   <- cyc + sw * integ[i]
      }
    }
    # ── MOTION BLUR ──────────────────────────────────────────────────────────
    # A MINFLUX localization is the time-average of the path over its cycle, not
    # the position at the start of it. For Brownian motion the average over one
    # cycle, conditioned on both end points, is the midpoint plus bridge noise of
    # variance D·cyc/6 per axis. That is the R = 1/6 blur factor the noise-corrected
    # CVE / MSD estimators are written for. Computed here, in the canonical frame,
    # so the bridge noise can be kept in the surface tangent plane; the last cycle
    # has no successor and stays unblurred.
    bp <- NULL
    if (blur) {
      t_cyc <- c(0, cumsum(dur)[-n_cyc])
      # Each recorded cycle is averaged over its own window, ending at the NEXT
      # recorded cycle. A gap makes both the window (up to tau_rep repeats) and the
      # step longer, which is the structure real files show through tau = eco/efo.
      rec  <- which(integ)
      nxt  <- seq_len(n_cyc); Tstep <- rep(cyc, n_cyc); taui <- rep(expfrac * cyc, n_cyc)
      if (length(rec) > 1L) {
        nxt[rec[-length(rec)]]   <- rec[-1]
        Tstep[rec[-length(rec)]] <- diff(t_cyc[rec])
        reps_used <- if (budget) rep_i[rec[-length(rec)]] else pmin(diff(rec), tau_rep)
        taui[rec[-length(rec)]]  <- reps_used * expfrac * cyc
      }
      mid <- pos + (taui / (2 * Tstep)) * (pos[nxt, , drop = FALSE] - pos)
      Dst <- if (identical(motion, "twostate")) ifelse(states == 0L, D_trk * slow_frac, D_trk)
             else rep(D_trk, n_cyc)
      # Time-average over a window tau <= T, conditioned on both end points:
      #   mean = x_i + (tau/(2T))*(x_j - x_i),  var = (2*D*tau/3)*(1 - 3*tau/(4T))
      # which reduces to the midpoint plus D*T/6 when tau = T (R = 1/6).
      br  <- matrix(rnorm(n_cyc * 3), ncol = 3) *
             sqrt(pmax(0, (2 * Dst * taui / 3) * (1 - 3 * taui / (4 * Tstep))))
      if (is_mem) br[, 3] <- 0                                    # in-plane only
      else if (surface && !is_free) {                             # tangent plane only
        nrm <- t(vapply(seq_len(n_cyc), function(i) geo_normal(pos[i, ]), numeric(3)))
        br  <- br - (rowSums(br * nrm)) * nrm
      }
      bp <- mid + br
      bp[n_cyc, ] <- pos[n_cyc, ]                  # no successor: unblurred
    }
    if (is_rod) { pos <- pos[, rod_perm, drop = FALSE]
                  if (!is.null(bp)) bp <- bp[, rod_perm, drop = FALSE] }
    if (is_mem) { pos <- pos %*% t(R_mem)            # tilt the flat local sheet into 3D
                  if (!is.null(bp)) bp <- bp %*% t(R_mem) }
    attr(pos, "states")  <- states
    attr(pos, "blurpos") <- bp
    attr(pos, "integ")   <- integ
    attr(pos, "eco")     <- eco_i
    attr(pos, "reps")    <- rep_i
    # Start time of each cycle: cycles last `cyc`, plus SW once per RECORDED cycle, so
    # two recorded cycles k apart are exactly k·cyc + SW apart.
    attr(pos, "t_cyc")   <- c(0, cumsum(dur)[-n_cyc])
    pos
  }
  # Survivorship: a molecule that moves > R_reacq during a dark gap is lost - the
  # track ends there. Reproduces the dark-cycle bias (long jumps culled → apparent
  # confinement). Off unless enabled.
  surv    <- isTRUE(input$sim_trk_surv)
  R_reacq <- max(1, as.numeric(input$sim_trk_rreacq %||% 300))

  rows <- vector("list", n_trk); gt <- vector("list", n_trk); k <- 0L
  withProgress(message = "Simulating tracks…", value = 0, {
  for (ti in seq_len(n_trk)) {
    if (ti %% 20L == 0L) incProgress(20 / n_trk)
    # mean_len = target MEDIAN recorded localizations/track. Recorded count is ≈
    # geometric, whose median ≈ ln(2)·mean, and only (1−p_dark) of cycles are recorded,
    # so the bleach lifetime in CYCLES has mean = median / (ln2 · (1−p_dark)) to land
    # the recorded MEDIAN on mean_len despite dark-cycling.
    # D and lifetime are drawn as a correlated pair, so a sample where the slow
    # molecules are the ones tracked longest can be reproduced.
    z1 <- stats::rnorm(1)
    z2 <- D_life * z1 + sqrt(1 - D_life^2) * stats::rnorm(1)
    life  <- mean_len / (log(2) * max(0.05, 1 - p_dark))
    if (life_sd > 0) life <- life * exp(life_sd * z2)
    n_cyc <- max(3L, stats::rgeom(1, 1 / max(2, life)) + 1L)
    # Each molecule may carry its own D (lognormal, median = the D input).
    D_trk <- if (D_sdlog > 0) D_nm2ms * exp(D_sdlog * z1) else D_nm2ms
    pos   <- gen_path(n_cyc, D_trk)                   # ground-truth path at cycle resolution
    integ <- attr(pos, "integ"); t_cyc <- attr(pos, "t_cyc")
    idx   <- which(integ)
    if (surv && length(idx) > 1L) {                   # survivorship: cut at first over-range dark-gap jump
      jmp  <- sqrt(rowSums((pos[idx[-1], , drop = FALSE] - pos[idx[-length(idx)], , drop = FALSE])^2))
      lost <- which(jmp > R_reacq)
      if (length(lost)) idx <- idx[seq_len(lost[1])]  # keep up to (and incl.) the pre-jump loc
    }
    if (length(idx) < 5L) next
    # Photons per localization: the budget already counted them (floored at the limit
    # by construction); otherwise draw them from the brightness model directly.
    N <- if (budget) attr(pos, "eco")[idx] else {
           r <- if (pht_sd > 0) photons * exp(stats::rnorm(length(idx), 0, pht_sd)) else photons
           pmax(pht_min, stats::rpois(length(idx), r))
         }
    tau_i <- if (budget) attr(pos, "reps")[idx] * expfrac * cyc else rep(expfrac * cyc, length(idx))
    # MINFLUX CRLB for a hexagonal TCP of diameter L with N detected photons.
    # Cross-checked on real tracking data: L = 100.8 nm (patGeoFactor 0.28) and the
    # median eco of ~38 photons give 4.1 nm, exactly the blur-aware CVE sigma of the
    # tracking test files.
    sig   <- L_tcp / (TCP_CRLB_PREFACTOR * sqrt(N))  # per-localization σ_loc (nm)
    bpos  <- attr(pos, "blurpos")                    # blurred path, or NULL when off
    true_pos <- if (!is.null(bpos)) bpos[idx, , drop = FALSE] else pos[idx, , drop = FALSE]
    obs   <- true_pos + matrix(rnorm(length(idx) * 3), ncol = 3) * sig
    # ── RE-ACQUISITION SHRINKAGE ───────────────────────────────────────────────
    # After a dark gap the molecule sits off the pattern centre, where the estimator
    # is biased back toward it, so the first localization after a gap is reported
    # closer to the previous one than it really is. Sequential: the centre is the
    # PREVIOUS REPORTED position, which is all the instrument has.
    if (shrink > 0 && length(idx) > 1L) {
      gap <- c(FALSE, diff(idx) > 1L)                # first loc after a dark gap
      for (j in which(gap)) obs[j, ] <- obs[j - 1L, ] + (1 - shrink) * (obs[j, ] - obs[j - 1L, ])
    }
    # Place the track on its cell: rotate in-plane, then offset to the cell centre.
    if (n_cells > 1L) {
      ci <- 1L + (ti - 1L) %% n_cells
      ca <- cos(cell_ang[ci]); sa <- sin(cell_ang[ci])
      ctr <- if (identical(geom, "rod")) c(Lcyl/2, 0, 0)[rod_perm] else c(0,0,0)
      ox <- obs[, 1] - ctr[1]; oy <- obs[, 2] - ctr[2]
      obs[, 1] <- ox * ca - oy * sa + cell_pos[ci, 1]
      obs[, 2] <- ox * sa + oy * ca + cell_pos[ci, 2]
      obs[, 3] <- obs[, 3] - ctr[3]
    }
    k <- k + 1L
    rows[[k]] <- data.frame(
      tid = k, tim = t_cyc[idx] / 1000,              # SECONDS (importer ×1000 → ms)
      x = obs[, 1], y = obs[, 2], z = obs[, 3],
      efo = N / (tau_i / 1000), eco = N, stringsAsFactors = FALSE)
    # Ground-truth state for HMM validation, numbered like the HMM (1 = slow, 2 = mobile).
    if (identical(motion, "twostate")) rows[[k]]$state_true <- attr(pos, "states")[idx] + 1L
    gt[[k]] <- data.frame(tid = k, D_true_um2_s = D_trk / 1000, geometry = geom, mode = mode,
                          motion = motion, L_c_nm = if (identical(motion,"confined")) L_c else NA_real_,
                          v_um_s = if (identical(motion,"directed")) v_nmms else NA_real_,
                          n_locs = length(idx), sigma_loc_med = stats::median(sig),
                          stringsAsFactors = FALSE)
  }
  })
  if (k == 0L) return(NULL)
  df <- dplyr::bind_rows(rows[seq_len(k)])
  if (proj2d) df$z <- NULL
  list(df = df, gt = dplyr::bind_rows(gt[seq_len(k)]),
       D_um2_s = D_um2_s, cyc = cyc, sw = sw, blur = blur, p_dark = p_dark, proj2d = proj2d,
       expfrac = expfrac, D_cv = D_cv, shrink = shrink,
       life_cv = life_cv, pht_min = pht_min, pht_sd = pht_sd, n_cells = n_cells, fov_nm = fov_nm,
       tau_rep = tau_rep, budget = budget, alpha = if (is_anom) alpha else NA_real_,
       pht_dec = pht_dec, L_hop = L_hop, p_perm = p_perm,
       D_life = D_life,
       mode = mode, d_mot = d_mot, geom = geom,
       motion = motion, L_c = L_c, v_nmms = v_nmms, surv = surv, R_reacq = R_reacq)
})

# CVE round-trip: recover D from consecutive single-cycle (Δt = cyc) steps, using the
# TOTAL displacement with the MOTION dimensionality. D = (⟨|Δr|²⟩ + 2⟨Δr_i·Δr_{i+1}⟩)
# / (2·d·Δt). SURFACE diffusion has d=2 DOF but spreads over 3 Cartesian axes (a
# per-axis CVE would under-estimate D by ⅓); VOLUME diffusion has d=3. When projected
# to 2D only n_cart=2 axes are seen, so d = min(motion-DOF, n_cart). n_cart Cartesian
# axes carry the noise (cov = −n_cart·σ²).
sim_trk_cve <- reactive({
  res <- sim_trk_result(); if (is.null(res)) return(NULL)
  df <- res$df; cyc <- res$cyc
  dt_step <- cyc + (res$sw %||% 0)                    # true single-cycle Δt = cycle + SW
  axes  <- c("x", "y", if (!is.null(df$z)) "z")
  ncart <- length(axes)
  d_mot <- min(res$d_mot %||% 2L, ncart)              # 2 surface / 3 volume, capped by observed axes
  r2_sum <- 0; r2_n <- 0L; cov_sum <- 0; cov_n <- 0L
  for (tb in split(df, df$tid)) {
    tb  <- tb[order(tb$tim), ]
    dt  <- c(NA_real_, diff(tb$tim) * 1000)           # ms
    one <- is.finite(dt) & abs(dt - dt_step) < cyc * 0.25  # single-cycle step (ending at this row)
    # per-step displacement vector (rows = steps ending at row i)
    dmat <- sapply(axes, function(ax) c(NA_real_, diff(tb[[ax]])))
    if (is.null(dim(dmat))) dmat <- matrix(dmat, ncol = ncart)
    r2  <- rowSums(dmat^2)                             # |Δr|² total
    s1  <- one & is.finite(r2)
    r2_sum <- r2_sum + sum(r2[s1]); r2_n <- r2_n + sum(s1)
    nxt  <- c(s1[-1], FALSE); pair <- s1 & nxt         # consecutive single-cycle pair
    dot  <- rowSums(dmat * rbind(dmat[-1, , drop = FALSE], NA_real_))  # Δr_i · Δr_{i+1}
    cov_sum <- cov_sum + sum(dot[pair], na.rm = TRUE); cov_n <- cov_n + sum(pair)
  }
  if (r2_n < 1L) return(NULL)
  r2m  <- r2_sum / r2_n                                # ⟨|Δr|²⟩ (nm²)
  cov1 <- if (cov_n > 0L) cov_sum / cov_n else 0       # ⟨Δr_i·Δr_{i+1}⟩ = −n_cart·σ²
  # Vestergaard's CVE: D is blur-independent, σ is not. With motion blur on, the
  # localization is a time-average over the collection window τ, so R = τ/(6·Δt) and
  # σ² = R⟨Δx²⟩ + (2R−1)⟨Δx_i Δx_{i+1}⟩ per axis; σ² = −cov is the R = 0 (snapshot)
  # special case. Reading it with the wrong R is what makes a real 4 nm σ_loc look
  # like 1.4 nm.
  Rblur <- if (isTRUE(res$blur)) (res$expfrac %||% 1) * cyc / (6 * dt_step) else 0
  # (single-cycle steps only, so the window is one repeat regardless of tau_rep)
  sig2  <- Rblur * (r2m / ncart) + (2 * Rblur - 1) * (cov1 / ncart)
  list(D_um2_s  = ((r2m + 2 * cov1) / (2 * d_mot * dt_step)) / 1000,   # CVE D (µm²/s)
       sigma_nm = sqrt(max(0, sig2)), R_blur = Rblur, n_pairs = cov_n)
})

# ── sim_trk_kcurve() - displacement vs dark-gap length ──────────────────────
# GOAL: the diagnostic the new sampling knobs exist to drive. Mean per-axis ⟨Δx²⟩
# against the number of cycles k a step spans, next to the free-diffusion
# expectation. Free Brownian motion grows linearly in k; real tracking data often
# flattens after k = 2, and re-acquisition shrinkage / speed-dependent dark cycles
# are the candidate causes this table is meant to discriminate.
# returns: data.frame(k, n, msq, expect, pct) or NULL.
sim_trk_kcurve <- reactive({
  res <- sim_trk_result(); if (is.null(res)) return(NULL)
  df  <- res$df; cyc <- res$cyc; sw <- res$sw %||% 0
  axes  <- c("x", "y", if (!is.null(df$z)) "z")
  ncart <- length(axes)
  d_mot <- min(res$d_mot %||% 2L, ncart)
  kk <- list(); dd <- list()
  for (tb in split(df, df$tid)) {
    tb <- tb[order(tb$tim), ]; if (nrow(tb) < 2L) next
    dt <- diff(tb$tim) * 1000
    kk[[length(kk) + 1L]] <- rep(pmax(1L, round((dt - sw) / cyc)), ncart)
    dd[[length(dd) + 1L]] <- unlist(lapply(axes, function(a) diff(tb[[a]])))
  }
  if (!length(kk)) return(NULL)
  k <- unlist(kk); dx <- unlist(dd)
  ok <- is.finite(k) & is.finite(dx) & k >= 1L & k <= 6L
  if (!any(ok)) return(NULL)
  agg <- stats::aggregate(list(msq = dx[ok]^2), by = list(k = k[ok]),
                          FUN = function(v) c(m = mean(v), n = length(v)))
  out <- data.frame(k = agg$k, msq = agg$msq[, "m"], n = agg$msq[, "n"])
  # Free-diffusion expectation per axis. Surface motion has d_mot DOF spread over
  # ncart Cartesian axes, so each axis carries d_mot/ncart of the motion. The blur
  # term −(2/3)·D·τ is subtracted because the simulated positions are time-averages.
  w_gt  <- pmax(0, res$gt$n_locs - 1)
  D_ref <- (if (sum(w_gt) > 0) stats::weighted.mean(res$gt$D_true_um2_s, w_gt) else res$D_um2_s) * 1000
  sig2  <- stats::median(res$gt$sigma_loc_med, na.rm = TRUE)^2
  tau   <- if (isTRUE(res$blur)) (res$expfrac %||% 1) * cyc else 0
  fr    <- d_mot / ncart
  out$expect <- 2 * D_ref * fr * (out$k * cyc + sw) + 2 * sig2 - (2/3) * D_ref * fr * tau
  out$pct    <- 100 * out$msq / out$expect
  out[order(out$k), ]
})

output$sim_trk_plot <- renderPlotly({
  res <- sim_trk_result()
  shiny::validate(shiny::need(!is.null(res), "Set parameters and click 'Simulate tracks'."))
  df <- res$df; is3d <- !is.null(df$z)
  show <- head(unique(df$tid), 60L); d <- df[df$tid %in% show, ]
  if (is3d)
    plotly::plot_ly(d, x = ~x, y = ~y, z = ~z, split = ~factor(tid),
                    type = "scatter3d", mode = "lines",
                    line = list(width = 2), showlegend = FALSE) %>%
      plotly::layout(title = "Simulated MINFLUX tracks (first 60)",
        scene = list(xaxis = list(title = "x (nm)"), yaxis = list(title = "y (nm)"),
                     zaxis = list(title = "z (nm)"),
                     aspectmode = "data"))   # equal nm-per-unit on X/Y/Z (no axis stretching)
  else
    plotly::plot_ly(d, x = ~x, y = ~y, split = ~factor(tid),
                    type = "scatter", mode = "lines",
                    line = list(width = 1), showlegend = FALSE) %>%
      plotly::layout(title = "Simulated MINFLUX tracks - 2D projection (first 60)",
        xaxis = list(title = "x (nm)"),
        yaxis = list(title = "y (nm)", scaleanchor = "x"))
})

output$sim_trk_summary <- renderUI({
  res <- sim_trk_result(); if (is.null(res)) return(NULL)
  df <- res$df; cve <- sim_trk_cve()
  sig_true <- stats::median(res$gt$sigma_loc_med, na.rm = TRUE)   # ground-truth σ_loc (nm)
  dt_all <- unlist(lapply(split(df$tim, df$tid), function(t) diff(sort(t)) * 1000), use.names = FALSE)
  frac_int <- if (length(dt_all))
                mean(abs(dt_all - (res$cyc + (res$sw %||% 0))) < res$cyc * 0.25) else NA_real_
  # With a heterogeneous population the pooled CVE cannot recover the median D: it
  # recovers the STEP-WEIGHTED MEAN of the drawn per-track D's, so that is what the
  # round-trip is scored against (the two coincide when the spread is 0).
  w_gt  <- pmax(0, res$gt$n_locs - 1)
  D_ref <- if (sum(w_gt) > 0) stats::weighted.mean(res$gt$D_true_um2_s, w_gt) else res$D_um2_s
  recov <- if (!is.null(cve)) {
    rr  <- cve$D_um2_s / D_ref
    col <- if (is.finite(rr) && abs(rr - 1) < 0.25) "#228833" else "#BB7700"
    htmltools::tags$span(style = sprintf("color:%s; font-weight:bold;", col),
      sprintf("%.4g µm²/s  (%.0f%% of truth)", cve$D_um2_s, 100 * rr))
  } else "-"
  motion_lbl <- switch(res$motion %||% "free",
    anomalous = sprintf("anomalous fBm (α = %.2f)", res$alpha),
    hop = sprintf("hop diffusion (corral %.0f nm, permeability %.2f)", res$L_hop, res$p_perm),
    confined = sprintf("confined (OU tether, L_c = %.0f nm)", res$L_c),
    directed = sprintf("directed (v = %.3g µm/s)", res$v_nmms),
    twostate = "two-state (mobile↔slow)",
    "free")
  mode_lbl <- if ((res$geom %||% "") %in% c("free","membrane")) ""
              else if (identical(res$mode, "volume")) " · volume (intracellular)" else " · surface"
  geom_lbl <- sprintf("%s%s · %s%s", res$geom %||% "?", mode_lbl,
                      motion_lbl, if (isTRUE(res$surv)) sprintf(" · survivorship R=%.0fnm", res$R_reacq) else "")
  htmltools::div(style = "font-size:0.9em; color:#333; margin-top:8px;",
    htmltools::tags$table(style = "border-collapse:collapse;",
      htmltools::tags$tr(htmltools::tags$td(style="padding:2px 14px 2px 0;color:#666;", "Model"),
        htmltools::tags$td(htmltools::tags$b(geom_lbl))),
      htmltools::tags$tr(htmltools::tags$td(style="padding:2px 14px 2px 0;color:#666;", "Tracks / localizations"),
        htmltools::tags$td(sprintf("%d / %s", length(unique(df$tid)), format(nrow(df), big.mark=",")))),
      htmltools::tags$tr(htmltools::tags$td(style="padding:2px 14px 2px 0;color:#666;", "Ground-truth D"),
        htmltools::tags$td(htmltools::tags$b(
          if ((res$D_cv %||% 0) > 0)
            sprintf("%.4g µm²/s median · CV %.2f · step-weighted mean %.4g · per-track 10–90%% %.3g–%.3g",
                    res$D_um2_s, res$D_cv, D_ref,
                    stats::quantile(res$gt$D_true_um2_s, 0.1), stats::quantile(res$gt$D_true_um2_s, 0.9))
          else sprintf("%.4g µm²/s", res$D_um2_s)))),
      htmltools::tags$tr(htmltools::tags$td(style="padding:2px 14px 2px 0;color:#666;", "Recovered D (CVE round-trip)"),
        htmltools::tags$td(recov)),
      htmltools::tags$tr(htmltools::tags$td(style="padding:2px 14px 2px 0;color:#666;", "Geometric bias factor"),
        htmltools::tags$td(if (!is.null(cve) && is.finite(cve$D_um2_s) && res$D_um2_s > 0) {
          bf <- cve$D_um2_s / D_ref
          htmltools::tags$span(style = if (abs(bf - 1) > 0.05) "color:#b26a00;font-weight:bold;" else "",
            sprintf("%.2f×  →  divide apparent D by %.2f to correct", bf, bf))
        } else "-")),
      htmltools::tags$tr(htmltools::tags$td(style="padding:2px 14px 2px 0;color:#666;", "σ_loc  truth / recovered"),
        htmltools::tags$td(sprintf("%.2f nm / %s", sig_true,
          if (!is.null(cve)) sprintf("%.2f nm  (CVE, R = %.2g)", cve$sigma_nm, cve$R_blur %||% 0) else "-"))),
      htmltools::tags$tr(htmltools::tags$td(style="padding:2px 14px 2px 0;color:#666;", "Sampling clock"),
        htmltools::tags$td(sprintf("Δt = k × %.4g ms + %.3g ms (SW) · blur %s%s%s",
          res$cyc, res$sw %||% 0,
          if (isTRUE(res$blur))
            sprintf("on (τ = %.0f%% of cycle, R = %.3g)", 100 * (res$expfrac %||% 1),
                    (res$expfrac %||% 1) * res$cyc / (6 * (res$cyc + (res$sw %||% 0))))
          else "off (snapshot)",
          if ((res$shrink %||% 0) > 0) sprintf(" · gap shrinkage %.0f%%", 100 * res$shrink) else "", ""))),
      htmltools::tags$tr(htmltools::tags$td(style="padding:2px 14px 2px 0;color:#666;", "Photons / dark cycles"),
        htmltools::tags$td(sprintf("eco median %.0f (limit %d) · %s",
          stats::median(df$eco, na.rm = TRUE), res$pht_min %||% 1L,
          if (isTRUE(res$budget))
            sprintf("dark %.0f%%, emergent from the photon budget", 100 * (1 - frac_int))
          else sprintf("dark %.0f%% as set", 100 * res$p_dark)))),
      htmltools::tags$tr(htmltools::tags$td(style="padding:2px 14px 2px 0;color:#666;", "Integrated (single-cycle) fraction"),
        htmltools::tags$td(if (is.finite(frac_int)) sprintf("%.0f%%", 100 * frac_int) else "-")),
      htmltools::tags$tr(htmltools::tags$td(style="padding:2px 14px 2px 0;color:#666;", "Median Δt / median track length"),
        htmltools::tags$td(sprintf("%.3g ms / %.0f locs",
          if (length(dt_all)) stats::median(dt_all) else NA_real_, stats::median(as.numeric(table(df$tid))))))),
    htmltools::tags$p(style = "color:#888; font-size:0.9em; margin-top:6px;",
      if (identical(res$motion, "confined")) {
        tau_c   <- if (res$D_um2_s > 0) (res$L_c^2 / (res$D_um2_s * 1000)) else Inf   # ms
        med_dur <- stats::median(unlist(lapply(split(df$tim, df$tid), function(t) (max(t) - min(t)) * 1000)))
        verdict <- if (is.finite(tau_c) && tau_c < 0.5 * med_dur) "\u2713 resolvable" else "\u26a0 not resolvable: looks free"
        sprintf("Confined (OU, L_c = %.0f nm): the CVE gives the local D; confinement shows as an MSD plateau \u2248 2\u00b7d\u00b7L_c\u00b2. \u03c4_c = L_c\u00b2/D \u2248 %.3g ms against a median track of %.3g ms: %s.",
                res$L_c, tau_c, med_dur, verdict)
      }
      else if (identical(res$motion, "directed"))
        "Directed: the CVE gives the diffusive part; the drift shows as an upward-curving MSD and a sustained drift velocity."
      else if (identical(res$motion, "twostate"))
        "Two-state: mixed mobile/immobile steps - the CVE gives a population-average D; load in the Tracking tab and fit the HMM to recover the two states."
      else
        "Free diffusion: the CVE round-trip should land near 100 % of truth - that is the validation loop.",
      " Download the tracks (CSV) and load them in the Tracking tab to test the full MSD / confinement pipeline against this known ground truth."),
    # Displacement-vs-gap-length table: the k-curve diagnostic of the plan doc.
    local({
      kc <- tryCatch(sim_trk_kcurve(), error = function(e) NULL)
      if (is.null(kc) || !nrow(kc)) return(NULL)
      cell <- function(x, b = FALSE) htmltools::tags$td(
        style = "padding:1px 10px 1px 0; text-align:right;", if (b) htmltools::tags$b(x) else x)
      htmltools::div(style = "margin-top:8px;",
        htmltools::tags$div(style = "color:#666; font-size:0.86em;",
          "Per-axis ⟨Δx²⟩ (nm²) by dark-gap length k"),
        htmltools::tags$table(style = "border-collapse:collapse; font-size:0.84em;",
          htmltools::tags$tr(cell("k"), lapply(kc$k, function(v) cell(v, TRUE))),
          htmltools::tags$tr(cell("steps"), lapply(kc$n, function(v) cell(format(v, big.mark = ",")))),
          htmltools::tags$tr(cell("⟨Δx²⟩"), lapply(kc$msq, function(v) cell(sprintf("%.0f", v)))),
          htmltools::tags$tr(cell("free expectation"), lapply(kc$expect, function(v) cell(sprintf("%.0f", v)))),
          htmltools::tags$tr(cell("% of free"), lapply(kc$pct, function(v)
            htmltools::tags$td(style = sprintf("padding:1px 10px 1px 0; text-align:right; %s",
              if (is.finite(v) && abs(v - 100) > 15) "color:#b26a00;font-weight:bold;" else ""),
              sprintf("%.0f%%", v))))),
        htmltools::tags$p(style = "color:#888; font-size:0.84em; margin-top:3px;",
          "Free diffusion holds this row near 100% at every k. Real tracking files fall well below ",
          "that by k=6, which is what the gap-shrinkage and speed-dependent-dark controls are for."))
    }),
    # Projection artifact: 2D imaging of a CURVED-surface protein under-measures the
    # azimuthal motion (its out-of-plane component is dropped), so the recovered D is a
    # fixed fraction of truth (~75% for a side-on cylinder) - NOT a bug, and the same
    # bias applies to real 2D data. The factor is itself a useful geometric calibration.
    if ((res$geom %||% "") %in% c("sphere","rod","ellipsoid") && identical(res$mode, "surface") && isTRUE(res$proj2d))
      htmltools::tags$p(style = "color:#b26a00; font-size:0.9em; margin-top:4px;",
        htmltools::HTML("&#9888; "),
        htmltools::tags$b("2D projection of a curved surface:"),
        " recovered D is expected BELOW truth (≈75% for a side-on cylinder) because the ",
        "out-of-plane part of the azimuthal motion is lost - this is a real geometric bias that ",
        "also affects your real 2D data. Use this ratio as the correction factor, or acquire in 3D. ",
        "It applies to free AND confined alike, so it is not the confinement signal."))
})

output$sim_trk_dl_csv <- downloadHandler(
  filename = function() paste0("sim_minflux_tracks_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
  content  = function(file) {
    res <- sim_trk_result()
    shiny::validate(shiny::need(!is.null(res), "Simulate tracks first."))
    write.csv(res$df, file, row.names = FALSE)
  })
output$sim_trk_dl_gt <- downloadHandler(
  filename = function() paste0("sim_ground_truth_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
  content  = function(file) {
    res <- sim_trk_result()
    shiny::validate(shiny::need(!is.null(res), "Simulate tracks first."))
    write.csv(res$gt, file, row.names = FALSE)
  })
