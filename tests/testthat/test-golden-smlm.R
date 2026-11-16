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
                             tid_length_range = c(8, 40), remove_rogues = TRUE,
                             view_SD = TRUE, target_sd_x = 6, target_sd_y = 6))
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
  })
  res <- as.list(res)
  expect_identical(.golden_errors(res), character(0))
  expect_golden("smlm_pipeline", res)
})
