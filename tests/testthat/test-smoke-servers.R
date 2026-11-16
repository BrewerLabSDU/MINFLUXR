# Smoke tests: load data, trigger every panel, render every output.
#
# CLAIM: with data loaded and each panel's Run button pressed, every output the
# three server functions define (plots, tables, HTML readouts, download handlers)
# has the render status stored in tests/testthat/_golden/smoke_*.rds:
#   ok      rendered a value
#   silent  stopped by req()/validate() (the panel is switched off, or needs an
#           input a test session cannot provide)
#   error   raised an error
# A refactor that makes a panel stop rendering, or start failing, changes its
# status and fails this test. The numbers behind the panels are pinned by
# test-golden-*.R; this file covers the render and export code around them.
#
# Runs in shiny::testServer (no browser). Output ids are read from the server
# code by tools/server_map.R, so a new output is covered without editing this file.

# GOAL: render one output and classify what happened.
.render_status <- function(output, id) {
  tryCatch({ v <- output[[id]]; "ok" },
           shiny.silent.error = function(e) "silent",
           validation = function(e) "silent",
           error = function(e) paste("error:", substr(conditionMessage(e), 1, 80)))
}

# GOAL: status of every output of one server, as a named character vector.
.render_all <- function(output, ids) {
  set.seed(1)
  vapply(sort(ids), function(id) .render_status(output, id), "")
}

# GOAL: the output ids one server function defines, from its source.
.server_outputs <- function(tab) {
  old <- setwd(.repo_root); on.exit(setwd(old), add = TRUE)
  env <- new.env(); sys.source(file.path("tools", "server_map.R"), envir = env)
  ids <- unique(unlist(strsplit(env$server_statement_map(tab)$outputs, ";", fixed = TRUE)))
  ids[nzchar(ids) & ids != "<dynamic>"]
}

# GOAL: compare a status vector with its golden copy, naming the outputs that moved.
.expect_status_golden <- function(name, status) {
  path <- file.path(.golden_dir(), paste0(name, ".rds"))
  if (identical(Sys.getenv("MINFLUXR_UPDATE_GOLDEN"), "true")) {
    dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
    saveRDS(status, path, version = 2)
    testthat::skip(paste("golden written:", basename(path)))
  }
  ref <- readRDS(path)
  ids <- union(names(ref), names(status))
  a <- ref[ids]; b <- status[ids]; a[is.na(a)] <- "(absent)"; b[is.na(b)] <- "(absent)"
  moved <- ids[a != b]
  testthat::expect_identical(
    if (length(moved)) paste0(moved, ": ", a[match(moved, ids)], " -> ", b[match(moved, ids)]) else character(0),
    character(0))
}

test_that("Tracking tab: every output keeps its render status", {
  skip_if_no_app()
  load_app()
  files  <- .golden_track_files()
  ids    <- .server_outputs("tracking")
  status <- NULL
  shiny::testServer(tracking_server, {
    trk_session_start(session, files)
    d <- ui_input_defaults(tracking_ui)
    shows <- grep("^trk_show_", names(d), value = TRUE)                  # open every panel
    set_inputs(session, stats::setNames(rep(list(TRUE), length(shows)), shows))
    set_inputs(session, list(trk_msd_unified = TRUE, trk_msd_cluster = TRUE, trk_hist_param = "efo",
                             trk_cal_files = names(trk_store())))
    session$flushReact()
    for (b in c("trk_cve_run", "trk_cve_run", "trk_msd_compute", "trk_tsa_run", "trk_jdd_run", "trk_cpd_compute",
                "trk_conf_ci_run", "trk_disp_run", "trk_dir_run", "trk_speed_run", "trk_hmm_run",
                "trk_seg_msd_run", "trk_dmap_run", "trk_cal_run", "trk_cmp_run")) {
      set.seed(1)
      do.call(session$setInputs, stats::setNames(list((input[[b]] %||% 0) + 1), b)); session$flushReact()
    }
    status <<- .render_all(output, ids)
  })
  expect_gt(sum(status == "ok"), 60)
  .expect_status_golden("smoke_tracking", status)
})

test_that("SMLM tab: every output keeps its render status", {
  skip_if_no_app()
  load_app()
  store  <- golden_smlm_store()
  ids    <- .server_outputs("smlm")
  status <- NULL
  shiny::testServer(smlm_server, {
    d <- ui_input_defaults(smlm_ui)
    all_files <- grep("_all_files$", names(d), value = TRUE)
    smlm_session_start(session, store, c(
      stats::setNames(rep(list(TRUE), length(all_files)), all_files),
      list(show_gt_estimate = TRUE, view_deblink_knee = TRUE, view_SD = TRUE, show_metadata_hist = TRUE,
           hist_variable = "efo", pcf_fit_xi = TRUE, viewer_enable_fit = TRUE, ripley_envelope = TRUE,
           ripley_nsim = 19, pcf_envelope = TRUE, pcf_nsim = 19)))
    for (b in c("run_batch_nn", "run_tid_count", "run_det_time", "run_tid_time", "run_prec_cloud", "run_crowd",
                "run_size_nn", "run_rog", "run_ripley", "run_pcf", "run_hotspot", "run_cross_nn", "run_ripley12",
                "run_pcf12", "run_kde_coloc", "run_viewer", "run_fiducial", "run_density", "auto_detect_dcr")) {
      set.seed(1)
      do.call(session$setInputs, stats::setNames(list(1), b)); session$flushReact()
    }
    status <<- .render_all(output, ids)
  })
  expect_gt(sum(status == "ok"), 80)
  .expect_status_golden("smoke_smlm", status)
})

test_that("Simulator tab: every output keeps its render status", {
  skip_if_no_app()
  load_app()
  ids    <- .server_outputs("sim")
  status <- NULL
  shiny::testServer(sim_server, {
    sim_session_start(session)
    for (b in c("sim_run", "sim_compare_run", "sim_trk_run")) {
      do.call(session$setInputs, stats::setNames(list(1), b)); session$flushReact()
    }
    status <<- .render_all(output, ids)
  })
  expect_gt(sum(status == "ok"), 12)
  .expect_status_golden("smoke_sim", status)
})
