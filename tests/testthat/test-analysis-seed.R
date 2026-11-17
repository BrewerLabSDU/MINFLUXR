# CLAIM: the analyses that draw random numbers are reproducible. Pressing Run
# twice gives the same numbers, whatever ran in between; a different seed gives
# different numbers; and the seed travels in the settings files (it is in the
# tabs' settings specs, which the file tests cover).

test_that("tracking: random analyses repeat exactly and follow the seed", {
  skip_if_no_app()
  load_app()
  files <- .golden_track_files()
  res <- new.env()
  shiny::testServer(tracking_server, {
    trk_session_start(session, files[1])
    session$setInputs(trk_show_msd = TRUE, trk_msd_unified = TRUE, trk_msd_compute = 1)   # the null builds on it
    session$flushReact()
    run <- function(n) {
      session$setInputs(trk_cve_run = n); session$flushReact()
      session$setInputs(trk_show_conf = TRUE, trk_conf_ci_nsim = 40, trk_conf_ci_run = n); session$flushReact()
      list(ci = trk_cve_bundle()[c("sig_pt_lo", "sig_pt_hi")],
           null = conf_ci_bundle()[, c("ci_lo05", "ci_hi95")])
    }
    res$a <- run(1)
    stats::runif(1000)                                   # anything that moves the session's generator
    res$b <- run(2)
    session$setInputs(trk_seed = 2); session$flushReact()
    res$c <- run(3)
  })
  expect_identical(res$a, res$b)
  expect_false(identical(res$a$null, res$c$null))
})

test_that("SMLM: the Ripley envelope repeats exactly and follows the seed", {
  skip_if_no_app()
  load_app()
  store <- golden_smlm_store()
  res <- new.env()
  shiny::testServer(smlm_server, {
    smlm_session_start(session, store, list(ripley_all_files = TRUE, ripley_envelope = TRUE,
                                            ripley_nsim = 19, ripley_unit = "tid"))
    for (k in 1:3) {
      if (k == 2) stats::runif(1000)
      if (k == 3) { session$setInputs(smlm_seed = 9); session$flushReact() }
      session$setInputs(run_ripley = k); session$flushReact()
      res[[letters[k]]] <- ripley_env_results()[, c("lo", "hi")]
    }
  })
  expect_identical(res$a, res$b)
  expect_false(identical(res$a, res$c))
})
