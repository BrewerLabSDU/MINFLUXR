# The observation window of the spatial statistics, end to end through the SMLM
# server: Ripley's L with the bounding box versus the footprint ("support")
# window, on a sample that fills only part of its box.
#
# CLAIM: when the labelling does not fill the frame, the footprint window is used
# by BOTH the data curve and the CSR envelope. The data curve of a random pattern
# then reads as random; with the bounding box it reads as strongly clustered,
# because the empty part of the box is counted as area.

# GOAL: CSR molecules in an L-shaped region covering 66% of a 6 x 6 um box, each
# seen as one trace (no re-blinks), so the expected L(r) - r is 0 at every r.
.l_shape_store <- function(seed = 9) {
  set.seed(seed)
  xy <- cbind(stats::runif(5000, 0, 6000), stats::runif(5000, 0, 6000))
  xy <- xy[xy[, 1] < 2500 | xy[, 2] < 2500, ][1:900, ]
  xy <- rbind(xy, c(0, 0), c(6000, 6000))               # pin the box corners
  df <- sim_smlm_traces(xy, blink_mean = 0, n_loc = 8L, sigma = 3, wobble = 0)
  df$mol <- NULL; df$tid <- as.character(df$tid)
  n <- nrow(df)
  df$tim <- stats::runif(n, 0, 1000); df$efo <- 6e4; df$cfr <- 0.3; df$eco <- 150; df$dcr <- 0.4
  df$.dataset_tag <- "lshp"
  list(`lshape.csv` = df)
}

test_that("the footprint window reaches the Ripley data curve and its envelope", {
  skip_if_no_app()
  load_app()
  store <- .l_shape_store()
  res <- new.env()
  shiny::testServer(smlm_server, {
    smlm_session_start(session, store, list(enable_two_color = FALSE, enable_cluster2 = FALSE,
                                            enable_clustering = FALSE))
    for (sup in c(FALSE, TRUE)) {
      set_inputs(session, list(spatial_use_support = sup, spatial_edge_correct = TRUE,
                               ripley_all_files = TRUE, ripley_envelope = TRUE, ripley_nsim = 19,
                               ripley_unit = "tid", ripley_r_max = 800, ripley_r_steps = 8))
      # The footprint's gap-fill sweep smooths a short curve with loess, which warns
      # on a sample this small; the warnings are not what is tested here.
      set.seed(1); suppressWarnings({ session$setInputs(run_ripley = 1 + sup); session$flushReact() })
      res[[if (sup) "fp" else "bb"]] <- list(curve = ripley_results(), env = ripley_env_results(),
                                             win = spatial_window_log())
    }
  })
  bb <- res$bb; fp <- res$fp
  expect_identical(bb$win$window, "bounding box")
  expect_identical(fp$win$window, "footprint")
  expect_lt(abs(fp$win$frac_of_bbox - 0.66), 0.05)          # the window found the sample's area

  last <- function(x) x$curve$Lr_minus_r[nrow(x$curve)]     # at r = 800 nm
  expect_gt(last(bb), 100)                                  # box: random points read as clustered
  expect_lt(abs(last(fp)), 15)                              # footprint: they read as random
  # the data curve sits inside its own envelope with the footprint, not with the box
  expect_true(last(fp) > fp$env$lo[nrow(fp$env)] - 5 && last(fp) < fp$env$hi[nrow(fp$env)] + 5)
  expect_gt(last(bb), bb$env$hi[nrow(bb$env)] + 50)
  # the envelope is simulated in the window the data uses, so it differs too
  expect_false(isTRUE(all.equal(bb$env$hi, fp$env$hi)))
})
