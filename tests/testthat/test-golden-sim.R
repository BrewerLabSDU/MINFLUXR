# Golden-output regression test for the Simulator tab.
#
# CLAIM: with a fixed seed and the UI's default settings, the SMLM simulator, the
# strategy comparison and the tracking simulator return the data stored in
# tests/testthat/_golden/sim_*.rds. The simulators draw many random numbers in a
# fixed order, so this also pins the ORDER of the draws: moving generator code
# out of the server must not reorder them.
#
# Regenerate with MINFLUXR_UPDATE_GOLDEN=true only for a deliberate change.

# GOAL: the data-carrying parts of a simulator result (data frames, numbers and
# short vectors), leaving out images and plot objects.
.sim_parts <- function(r) {
  if (!is.list(r) || is.data.frame(r)) return(r)
  r[vapply(r, function(v) is.data.frame(v) || (is.atomic(v) && length(v) <= 5e4) ||
                          (is.list(v) && !inherits(v, c("plotly", "htmlwidget", "shiny.tag"))), logical(1))]
}

test_that("Simulator: SMLM ground truth and localizations reproduce their golden data", {
  skip_if_no_app()
  load_app()
  res <- new.env()
  shiny::testServer(sim_server, {
    sim_session_start(session)
    session$setInputs(sim_run = 1); session$flushReact()
    res$smlm_default <- .grab(.sim_parts(sim_result()))
    res$acq_live     <- .grab(.sim_parts(sim_acq_live()))
    res$crlb         <- .grab(crlb_inputs())
    # every structure generator once, so each has a pinned output
    shapes <- c("cluster", "ring", "arc", "double_ring", "ellipse", "hollow_sphere", "oligoring")
    for (i in seq_along(shapes)) {
      session$setInputs(sim_shape = shapes[i]); session$flushReact()
      session$setInputs(sim_run = i + 1L); session$flushReact()
      res[[paste0("smlm_", shapes[i])]] <- .grab({
        r <- sim_result(); .sim_parts(r)[intersect(names(r), c("locs", "gt", "df", "mol"))] })
    }
  })
  res <- as.list(res)
  expect_identical(.golden_errors(res), character(0))
  expect_golden("sim_smlm", res)
})

test_that("Simulator: tracking simulation reproduces its golden data", {
  skip_if_no_app()
  load_app()
  res <- new.env()
  shiny::testServer(sim_server, {
    sim_session_start(session)
    session$setInputs(sim_trk_run = 1); session$flushReact()
    res$trk_default <- .grab(.sim_parts(sim_trk_result()))
    res$trk_cve     <- .grab(sim_trk_cve())
    res$trk_kcurve  <- .grab(sim_trk_kcurve())
  })
  res <- as.list(res)
  expect_identical(.golden_errors(res), character(0))
  expect_golden("sim_tracking", res)
})

test_that("Simulator: strategy comparison reproduces its golden data", {
  skip_if_no_app()
  load_app()
  res <- new.env()
  shiny::testServer(sim_server, {
    sim_session_start(session)
    session$setInputs(sim_run = 1); session$flushReact()          # the comparison needs a ground truth
    session$setInputs(sim_compare_run = 1); session$flushReact()
    res$compare <- .grab(.sim_parts(compare_result()))
  })
  res <- as.list(res)
  expect_identical(.golden_errors(res), character(0))
  expect_golden("sim_compare", res)
})
