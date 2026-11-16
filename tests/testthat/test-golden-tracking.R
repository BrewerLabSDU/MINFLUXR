# Golden-output regression test for the Tracking tab.
#
# CLAIM: for a fixed simulated dataset and the UI's default settings, every
# analysis panel of the Tracking tab returns the numbers stored in
# tests/testthat/_golden/tracking_server.rds. The dataset goes through the tab's
# own upload path, Plot Data and each panel's Run button in a test session, so
# the whole chain is covered: reader, chunking, cleaning, filters, cycle time,
# blur, noise floor and the estimator.
#
# This guards refactoring. A golden may only change together with a deliberate
# change of method; regenerate with MINFLUXR_UPDATE_GOLDEN=true and say why in
# the commit.

test_that("Tracking tab: every analysis panel reproduces its golden numbers", {
  skip_if_no_app()
  load_app()
  files <- .golden_track_files()
  res <- new.env()
  shiny::testServer(tracking_server, {
    trk_session_start(session, files[1])
    press <- function(...) { set_inputs(session, list(...)); session$flushReact() }
    res$data        <- .grab(as.data.frame(trk())[, c("tid", "x", "y", "tim")])
    res$cycle_time  <- .grab(list(t_c = trk_cycle_time(), cluster = trk_cluster_cycle(),
                                   tot = trk_tot_cycles(), blur = trk_blur()))
    # Data quality first: it publishes the photon sigma the later panels subtract.
    press(trk_cve_run = 1)
    res$cve         <- .grab(trk_cve_bundle())
    session$flushReact()
    # The FIRST run must publish sigma: the panels below subtract it, and a user
    # presses the button once (fixed 2026-11; before, only the second press did).
    res$sigma_loc_first_run <- .grab(list(eff = trk_sigma_loc_eff(), src = trk_sigma_loc_src()))
    expect_true(is.finite(trk_sigma_loc_eff()))
    expect_identical(trk_sigma_loc_src(), "photon budget")
    press(trk_cve_run = 2)
    res$sigma_loc   <- .grab(list(eff = trk_sigma_loc_eff(), src = trk_sigma_loc_src()))
    press(trk_show_msd = TRUE, trk_msd_unified = TRUE, trk_uni_bootstrap = TRUE, trk_msd_cluster = TRUE)
    set.seed(1); press(trk_msd_compute = 1)
    res$msd_bundle  <- .grab(msd_unified_bundle())
    res$msd_fit     <- .grab(msd_unified_fit(msd_unified_bundle()))
    res$msd_model   <- .grab(msd_unified_modelfit(msd_unified_bundle()))
    res$msd_tdiff   <- .grab(trk_unified_tdiff())
    res$msd_boot    <- .grab(trk_unified_boot())
    res$msd_cluster <- .grab(msd_cluster_rv())
    res$msd_cluster_models <- .grab(msd_cluster_models())
    res$msd_cluster_diff   <- .grab(trk_unified_cluster_diff())
    res$msd_bindiag <- .grab(trk_msd_bindiag())
    res$msd_export  <- .grab(msd_export())
    press(trk_msd_geom = TRUE)
    res$msd_geom_null <- .grab(msd_geom_null())
    press(trk_show_tsa = TRUE, trk_tsa_run = 1)
    res$tsa         <- .grab(cycle_class_bundle())
    res$dt_vectors  <- .grab(dt_vectors())
    press(trk_show_jdd = TRUE, trk_jdd_run = 1)
    res$jdd         <- .grab(jdd_bundle())
    press(trk_show_cpd = TRUE, trk_cpd_compute = 1)
    res$cpd         <- .grab(cpd_segments_rv())
    res$cpd_gmm     <- .grab(cpd_gmm_rv())
    press(trk_show_conf = TRUE); set.seed(1); press(trk_conf_ci_run = 1)
    res$conf_ci     <- .grab(conf_ci_bundle())
    press(trk_show_disp = TRUE); set.seed(1); press(trk_disp_run = 1)
    res$disp        <- .grab(trk_disp_bundle())
    res$disp_null   <- .grab(trk_disp_geom_null())
    press(trk_show_dir = TRUE, trk_dir_run = 1)
    res$dir         <- .grab(dir_bundle())
    press(trk_show_speed = TRUE, trk_speed_run = 1)
    res$transport   <- .grab(trk_export_summary())
    press(trk_show_hmm = TRUE); set.seed(1); press(trk_hmm_run = 1)
    res$hmm         <- .grab({ h <- hmm_bundle(); h[setdiff(names(h), "fit")] })
    press(trk_show_seg_msd = TRUE, trk_seg_msd_run = 1)
    res$seg_msd     <- .grab(seg_msd_rv())
    press(trk_show_dmap = TRUE, trk_dmap_run = 1)
    res$dmap        <- .grab(dmap_bundle())
    res$dmax        <- .grab(dmax_estimate(trk(), 100.8, 0.05))
  })
  res <- as.list(res)
  expect_identical(.golden_errors(res), character(0))
  expect_golden("tracking_server", res)
})

test_that("Tracking tab: cleaning and filters reproduce their golden dataset", {
  skip_if_no_app()
  load_app()
  files <- .golden_track_files()
  res <- new.env()
  shiny::testServer(tracking_server, {
    trk_session_start(session, files[1])
    stage_counts <- function() vapply(
      list(raw = trk_raw_df(), split = trk_gap_split(), trim = trk_end_trimmed(),
           gate = trk_quality_filtered(), meta = trk_meta_filtered(), len = trk_len_filtered(),
           dur = trk_dur_filtered(), spatial = trk_spatial_filtered()),
      function(d) c(rows = nrow(d), tids = length(unique(d$tid))), numeric(2))
    res$stages_default <- .grab(stage_counts())
    set_inputs(session, list(trk_split_gaps_enable = TRUE, trk_split_adaptive = FALSE, trk_split_gap_cycles = 3,
                             trk_split_min_len = 20, trk_endtrim_enable = TRUE, trk_endtrim_dark = TRUE,
                             trk_quality_enable = TRUE, trk_quality_min_integ = 0.5))
    trk_filt$eco <- c(30, 400); trk_filt$tid_len <- c(25, 1000)
    trk_spatial_filter(list(x = c(0, 3000), y = c(0, 3000), pct_thr = 0.5))
    session$setInputs(trk_plot_data = 2); session$flushReact()
    res$stages_cleaned <- .grab(stage_counts())
    res$data_cleaned   <- .grab(as.data.frame(trk())[, c("tid", "x", "y", "tim")])
    res$filter_params  <- .grab(trk_filter_params())
  })
  res <- as.list(res)
  expect_identical(.golden_errors(res), character(0))
  expect_golden("tracking_pipeline", res)
})

test_that("Tracking tab: calibration and Compare Cells reproduce their golden tables", {
  skip_if_no_app()
  load_app()
  files <- .golden_track_files()
  res <- new.env()
  shiny::testServer(tracking_server, {
    trk_session_start(session, files)
    res$store <- .grab(names(trk_store()))
    session$setInputs(trk_cal_files = names(trk_store())); session$flushReact()
    session$setInputs(trk_cal_run = 1); session$flushReact()
    res$calibration <- .grab({ t <- trk_cal_result(); attr(t, "curves") <- NULL; as.data.frame(t) })
    res$calibration_summary <- .grab(calib_summary(trk_cal_result()))
    session$setInputs(trk_cmp_run = 1); session$flushReact()
    res$compare_cells <- .grab(cmp_results())
  })
  res <- as.list(res)
  expect_identical(.golden_errors(res), character(0))
  expect_golden("tracking_cells", res)
})

# CLAIM: the estimators give the same numbers when called directly, without a
# session. This pins R/est_*.R on its own, so a server change and an estimator
# change cannot hide each other.
test_that("Tracking estimators (R/est_*.R) reproduce their golden numbers", {
  skip_if_no_app()
  load_app()
  df <- golden_tracks(); df$tim <- df$tim * 1000            # ms, as the tab holds it
  t_c <- 0.192; sw <- 0.011; R <- 0.106
  res <- list()
  res$cycle   <- .grab(cycle_cluster(df$tid, df$tim))
  res$cve     <- .grab({ c <- cve_estimate(df, t_c); c[setdiff(names(c), c("pair", "prod_x", "prod_y", "prod_z"))] })
  res$msd     <- .grab(msd_headline(df, t_c, sw, R = R))
  res$msd_ens <- .grab({ s <- msd_track_sums(df, t_c, sw, 60L, 50L); msd_ensemble(s, t_c, sw) })
  res$jdd     <- .grab({ cj <- jdd_collect_jumps(df, min_pts = 10, integ_only = TRUE, t_c = t_c)
                         jdd_fit(cj$jumps, cj$dt_all, n_comp = 2, sigma_loc = 4, R_abs = R * t_c) })
  res$dmap    <- .grab(dmap_compute(df, 8L, 20L, t_c, TRUE, 4, R_abs = R * t_c))
  res$cpd     <- .grab(cpd_segment_tracks(df, min_seg = 8L, noise_sig2 = 16, R_abs = R * t_c))
  res$clean   <- .grab({ d <- split_gaps(df, 3 * t_c, 20L); d <- trim_track_ends(d, t_c); d <- gate_continuity(d, t_c, 0.5)
                         list(n = nrow(d), tids = length(unique(d$tid)), scores = continuity_scores(d, t_c)) })
  res$hmm     <- .grab({ pr <- hmm_ct_prepare(df, ndim = 2L, sigma = 4, min_pts = 10L)
                         f <- hmm_ct_fit(pr, K = 2L, R_abs = R * t_c, t_c = t_c)
                         f[intersect(names(f), c("D", "Q", "loglik", "pi", "dwell_ms", "bic", "aic"))] })
  res$percell <- .grab(percell_summary(df, "cell_one", L = 100.8))
  expect_identical(.golden_errors(res), character(0))
  expect_golden("tracking_estimators", res)
})
