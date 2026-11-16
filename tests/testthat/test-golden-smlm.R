# Golden-output regression test for the SMLM tab.
#
# CLAIM: for two fixed simulated two-channel files and the UI's default settings
# (with de-blinking, the two-colour split and Median Clustering switched on),
# both data chains and every analysis panel of the SMLM tab return the numbers
# stored in tests/testthat/_golden/smlm_*.rds. The files are placed in the tab's
# file store, Plot is pressed and each panel's Run button is pressed in a test
# session.
#
# This guards refactoring. A golden may only change together with a deliberate
# change of method; regenerate with MINFLUXR_UPDATE_GOLDEN=true and say why in
# the commit.

# Columns of the analysis table that carry results (colours and hover text are
# presentation and are left out).
.smlm_cols <- function(df) {
  keep <- intersect(c("tid", "x", "y", "z", "tim", "cluster", "cluster2", "two_color_group",
                      "channel", "channel_label", "n_in_cluster2", "analysis_group"), names(df))
  as.data.frame(df)[, keep, drop = FALSE]
}

test_that("SMLM tab: the active-file chain reproduces its golden numbers", {
  skip_if_no_app()
  load_app()
  store <- golden_smlm_store()
  res <- new.env()
  shiny::testServer(smlm_server, {
    smlm_session_start(session, store,
                       list(show_gt_estimate = TRUE, view_deblink_knee = TRUE, view_SD = TRUE))
    res$filtered        <- .grab(.smlm_cols(filtered_df()))
    res$plottable       <- .grab(.smlm_cols(plottable_data()))
    res$nn_table        <- .grab(nn_table())
    res$nn_cluster2     <- .grab(nn_table_cluster2())
    res$nn_lines        <- .grab(secondary_nn_lines())
    res$cross_dcr_nn    <- .grab(cross_dcr_nn_table())
    res$sd_all          <- .grab(sd_all_tbl())
    res$sd_filtered     <- .grab(sd_filtered_tbl())
    res$channel_spread  <- .grab(channel_spread_table())
    res$deblink_knee    <- .grab(deblink_knee_data())
    res$gt_estimate     <- .grab(gt_estimate())
    res$gt_morphology   <- .grab(gt_morphology())
    res$settings        <- .grab({ s <- get_ui_settings(); s[["smlm_settings"]] })
    set.seed(1); session$setInputs(auto_detect_dcr = 1); session$flushReact()
    res$dcr_gmm         <- .grab(dcr_gmm_result())
    res$dcr_conf        <- .grab(dcr_conf_tbl())
    set.seed(1); session$setInputs(suggest_eps_centroid = 1); session$flushReact()
    res$centroid_suggestion <- .grab(centroid_suggestion())
    set_inputs(session, list(viewer_enable_fit = TRUE)); set.seed(1)
    session$setInputs(run_viewer = 1); session$flushReact()
    res$viewer_ids      <- .grab(viewer_ids())
    res$viewer_current  <- .grab({ v <- viewer_current(); v[intersect(names(v), c("mat", "id", "n"))] })
    res$viewer_fit      <- .grab({ f <- viewer_fit(); if (is.null(f)) NULL else f[intersect(names(f), c("table", "best", "rmsd"))] })
  })
  res <- as.list(res)
  expect_identical(.golden_errors(res), character(0))
  expect_golden("smlm_active_chain", res)
})

test_that("SMLM tab: every batch analysis panel reproduces its golden numbers", {
  skip_if_no_app()
  load_app()
  store <- golden_smlm_store()
  res <- new.env()
  shiny::testServer(smlm_server, {
    defaults <- ui_input_defaults(smlm_ui)
    all_files <- grep("_all_files$", names(defaults), value = TRUE)     # every panel on both files
    extra <- c(stats::setNames(rep(list(TRUE), length(all_files)), all_files),
               list(ripley_envelope = TRUE, ripley_nsim = 19, pcf_envelope = TRUE, pcf_nsim = 19,
                    pcf_fit_xi = TRUE))
    smlm_session_start(session, store, extra)
    press <- function(id) { set.seed(1); do.call(session$setInputs, stats::setNames(list(1), id)); session$flushReact() }
    res$prepared_A   <- .grab(.smlm_cols(prepared_df(names(store)[1])))
    res$prepared_B   <- .grab(.smlm_cols(prepared_df(names(store)[2])))
    res$prepared_noclust <- .grab(.smlm_cols(prepared_df(names(store)[1], needs_cluster = FALSE)))
    res$groups       <- .grab(file_group_map())
    press("run_batch_nn");   res$batch_nn    <- .grab(nn_batch_results())
    res$batch_nn_filtered <- .grab(nn_batch_filtered())
    press("run_tid_count");  res$tid_count   <- .grab(tid_count_results())
    press("run_det_time");   res$det_time    <- .grab(det_time_results())
    press("run_tid_time");   res$tid_time    <- .grab(tid_time_results())
    press("run_prec_cloud"); res$prec_cloud  <- .grab(prec_cloud_results())
    press("run_crowd");      res$crowd       <- .grab(crowd_results())
    res$crowd_null   <- .grab(crowd_null_stats())
    press("run_size_nn");    res$size_nn     <- .grab(size_nn_results())
    press("run_rog");        res$rog         <- .grab(rog_results())
    press("run_ripley");     res$ripley      <- .grab(ripley_results())
    res$ripley_env   <- .grab(ripley_env_results())
    press("run_pcf");        res$pcf         <- .grab(pcf_results())
    res$pcf_env      <- .grab(pcf_env_results())
    res$pcf_xi       <- .grab(pcf_xi_results())
    res$spatial_log  <- .grab(spatial_window_log())
    press("run_hotspot");    res$hotspot     <- .grab(hotspot_rows())
    res$hotspot_detail <- .grab(hotspot_detail())
    press("run_cross_nn");   res$cross_nn    <- .grab(cross_nn_results())
    press("run_ripley12");   res$ripley12    <- .grab(ripley12_results())
    res$ripley12_env <- .grab(ripley12_env_results())
    press("run_pcf12");      res$pcf12       <- .grab(pcf12_results())
    res$pcf12_env    <- .grab(pcf12_env_results())
    press("run_kde_coloc");  res$kde_coloc   <- .grab(kde_coloc_results())
    # The two panels that change the point set of all the others come last.
    press("run_fiducial");   res$fiducial    <- .grab(fiducial_results())
    press("run_density");    res$density     <- .grab(density_profile_rv())
    res$density_sweep <- .grab(density_sweep_rv())
    set_inputs(session, list(dens_norm_enable = TRUE, dens_thin_enable = TRUE)); session$flushReact()
    res$density_target <- .grab(list(target = density_target(), p = density_p_map()))
    set.seed(1)
    res$prepared_thinned <- .grab(.smlm_cols(prepared_df(names(store)[1])))
  })
  res <- as.list(res)
  expect_identical(.golden_errors(res), character(0))
  expect_golden("smlm_batch_panels", res)
})

test_that("SMLM tab: filters, trimming and crops reproduce their golden datasets", {
  skip_if_no_app()
  load_app()
  store <- golden_smlm_store()
  res <- new.env()
  shiny::testServer(smlm_server, {
    smlm_session_start(session, store)
    res$default_active <- .grab(.smlm_cols(plottable_data()))
    rng <- range(store[[1]]$efo)
    set_inputs(session, list(efo_range = c(rng[1] * 1.3, rng[2] * 0.8), trim_outliers = TRUE,
                             tid_length_range = c(8, 40), view_SD = TRUE))
    sd_filt$x <- c(0, 6); sd_filt$y <- c(0, 6)
    smlm_2d_crop_rv(list(x = c(100, 500), y = c(100, 600)))
    session$flushReact()
    session$setInputs(plot_data = 2); session$flushReact()
    res$stage_rows <- .grab(c(meta = nrow(metadata_filtered_df()), trim = nrow(trimmed_df()),
                              spread = nrow(spread_filtered_df()), len = nrow(filtered_df())))
    res$filtered_active  <- .grab(.smlm_cols(plottable_data()))
    res$filtered_batch_A <- .grab(.smlm_cols(prepared_df(names(store)[1])))
    res$filtered_batch_B <- .grab(.smlm_cols(prepared_df(names(store)[2])))
    res$trim_flags       <- .grab(sum(trim_flags()))
    # The crop is the ONLY thing the main plot applies that the all-files chain
    # does not: inside the crop box the two hold the same localizations.
    cr <- smlm_2d_crop_rv()
    b  <- prepared_df(names(store)[1], needs_cluster = FALSE)
    b  <- b[b$x >= cr$x[1] & b$x <= cr$x[2] & b$y >= cr$y[1] & b$y <= cr$y[2], ]
    expect_setequal(paste(b$tid, signif(b$x, 12)), paste(plottable_data()$tid, signif(plottable_data()$x, 12)))
  })
  res <- as.list(res)
  expect_identical(.golden_errors(res), character(0))
  expect_golden("smlm_pipeline", res)
})

# CLAIM: the filters set in the Display tab mean the same on the main plot and in
# every all-files panel. For the active file the two chains keep exactly the same
# localizations (crops aside: a crop is set per dataset on the main plot),
# and the same thresholds are applied to every other file.
test_that("SMLM tab: the plot's filter chain and the all-files chain agree", {
  skip_if_no_app()
  load_app()
  store <- golden_smlm_store()
  key <- function(d) paste(d$tid, signif(d$x, 12), signif(d$y, 12))
  shiny::testServer(smlm_server, {
    smlm_session_start(session, store)
    A <- names(store)[1]; B <- names(store)[2]
    same_points <- function() expect_setequal(key(filtered_df()), key(prepared_df(A, needs_cluster = FALSE)))
    same_points()                                                  # nothing set
    expect_equal(nrow(filtered_df()), nrow(store[[A]]))

    # 1. every stage switched on
    rng <- range(store[[A]]$efo)
    set_inputs(session, list(efo_range = c(rng[1] * 1.3, rng[2] * 0.8), trim_outliers = TRUE,
                             tid_length_range = c(8, 40)))
    sd_filt$x <- c(0, 5); sd_filt$y <- c(1, 5.5)                   # upper bounds and a lower bound
    session$flushReact()
    same_points()
    n_all <- nrow(filtered_df())
    expect_lt(n_all, nrow(store[[A]]))

    # the spread ranges reach the all-files chain: opening them lets traces back in
    sd_filt$x <- c(0, Inf); sd_filt$y <- c(0, Inf); session$flushReact()
    same_points()
    expect_gt(nrow(prepared_df(A, needs_cluster = FALSE)), n_all)

    # 2. file B gets the same thresholds, as numbers
    b <- prepared_df(B, needs_cluster = FALSE)
    expect_true(all(b$efo >= rng[1] * 1.3 & b$efo <= rng[2] * 0.8))
    expect_true(all(table(b$tid) >= 8 & table(b$tid) <= 40))

    # 3. untouched sliders filter no file. The sliders sit at the ACTIVE file's
    # full range; file B reaches beyond it and must not be clipped to it.
    set_inputs(session, list(trim_outliers = FALSE, tid_length_range = c(1, 1000),
                             efo_range = rng))
    meta_full_rng$efo <- rng
    session$flushReact()
    expect_true(any(store[[B]]$efo > rng[2] | store[[B]]$efo < rng[1]))
    expect_equal(nrow(prepared_df(B, needs_cluster = FALSE)), nrow(store[[B]]))
    expect_equal(nrow(filtered_df()), nrow(store[[A]]))
    # ... and once the user narrows it, it applies to B as well
    session$setInputs(efo_range = c(rng[1], rng[2] * 0.9)); session$flushReact()
    expect_lt(nrow(prepared_df(B, needs_cluster = FALSE)), nrow(store[[B]]))
    same_points()
  })
})

# CLAIM: the stage functions do what their headers say, on a table small enough
# to check by hand.
test_that("SMLM filter stages: untouched ranges, open bounds, undefined spread", {
  df <- data.frame(tid = rep(c("a", "b", "c"), c(6, 6, 2)),
                   x = c(0, 1, 0, 1, 0, 1,  0, 30, 0, 30, 0, 30,  5, 5),
                   y = 0, efo = c(1:6, 11:16, 21, 22))
  # a range equal to the stated full range is untouched, whatever the file holds
  expect_equal(nrow(smlm_filter_meta(df, list(efo = c(1, 16)), list(efo = c(1, 16)))), 14)
  # with no stated full range the file's own range decides
  expect_equal(nrow(smlm_filter_meta(df, list(efo = c(1, 22)))), 14)
  expect_equal(nrow(smlm_filter_meta(df, list(efo = c(1, 16)))), 12)
  expect_equal(nrow(smlm_filter_meta(df, list(efo = NULL, cfr = c(0, 1)))), 14)   # NULL / absent column
  # length: both ends, and an open upper end
  expect_equal(unique(smlm_filter_length(df, c(3, 10))$tid), c("a", "b"))
  expect_equal(unique(smlm_filter_length(df, c(3, Inf))$tid), c("a", "b"))
  expect_equal(nrow(smlm_filter_length(df, NULL)), 14)
  # spread: trace b (sd about 16 nm in x) fails an upper bound of 5; trace c has
  # too few localizations for a spread and is not judged
  out <- smlm_filter_spread(df, c(0, 5), c(0, Inf), c(0, Inf), min_n = 5)
  expect_equal(sort(unique(out$tid)), c("a", "c"))
  expect_equal(nrow(smlm_filter_spread(df, c(0, Inf), c(0, Inf), c(0, Inf), min_n = 5)), 14)
  # a lower bound alone is a filter too: only the scattered trace passes it
  expect_equal(sort(unique(smlm_filter_spread(df, c(5, Inf), c(0, Inf), c(0, Inf), 5)$tid)), c("b", "c"))
  # the composed chain with spread ranges given ignores the legacy target_sd settings
  s <- smlm_settings(view_SD = TRUE, target_sd_x = 0.1, target_sd_y = 0.1)
  expect_equal(nrow(apply_meta_filters_batch(df, s, sd_ranges = list(x = c(0, Inf), y = c(0, Inf), z = c(0, Inf)))), 14)
})
