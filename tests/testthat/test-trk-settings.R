# The tracking settings spec (R/trk_settings.R).
source(file.path(.repo_root, "R", "trk_settings.R"))

# CLAIM: the spec and the UI agree. Every setting has a widget in ui_tracking.R
# of the kind the spec says, starting at the spec's default. A failure means one
# of the two was changed without the other.
test_that("every tracking setting has a widget that starts at the spec default", {
  skip_if_no_app()
  load_app()
  d <- ui_input_defaults(tracking_ui)
  ui_txt <- paste(readLines(file.path(.repo_root, "ui_tracking.R"), warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  ctor <- c(numeric = "numericInput", checkbox = "checkboxInput", select = "selectInput",
            radio = "radioButtons", slider = "sliderInput", text = "textInput")
  for (id in names(TRK_SETTINGS_SPEC)) {
    sp <- TRK_SETTINGS_SPEC[[id]]
    expect_true(grepl(sprintf("%s\\(\\s*\"%s\"", ctor[[sp$widget]], id), ui_txt), info = paste(id, "widget"))
    if (identical(id, "trk_cal_profile")) next                 # choices are filled by the server
    expect_identical(.trk_coerce(id, d[[id]]), sp$default, info = id)
  }
})

test_that("settings coerce, reject unknown names and hash stably", {
  s <- trk_settings(trk_uni_min_lag = "7", trk_cpd_sub_noise = 0, trk_tcp_pattern = "pentagon")
  expect_identical(s$trk_uni_min_lag, 7L)
  expect_false(s$trk_cpd_sub_noise)
  expect_identical(s$trk_tcp_pattern, "hexagon")              # not a pattern: the default
  expect_error(trk_settings(trk_uni_min_lagg = 7), "unknown setting")
  expect_true(is.na(trk_settings(trk_crlb_pref = NA)$trk_crlb_pref))
  expect_identical(trk_settings_key(trk_settings()), trk_settings_key(trk_settings(trk_uni_min_lag = 5)))
  expect_false(identical(trk_settings_key(trk_settings()), trk_settings_key(trk_settings(trk_seed = 2))))
  expect_identical(trk_settings_from_input(list())$trk_msd_min_pts, 50L)
})

test_that("a file block restores only what it recorded", {
  raw <- list(settings_version = 2L, trk_settings = list(trk_uni_min_lag = 8, trk_L_tcp = 151.2,
                                                         trk_crlb_pref = NULL, trk_bogus = 1))
  # JSON null (NULL here) for a setting where blank is a value restores it as blank
  f <- trk_settings_from_file(raw)
  expect_setequal(f$present, c("trk_uni_min_lag", "trk_L_tcp", "trk_crlb_pref"))
  expect_true(is.na(f$settings$trk_crlb_pref))
  expect_identical(f$settings$trk_uni_min_lag, 8L)
  expect_null(trk_settings_from_file(list(efo_rng = c(1, 2))))   # an older file: no block
  plan <- trk_settings_restore_plan(trk_settings(trk_uni_iterative = FALSE)[c("trk_uni_min_lag", "trk_uni_iterative")])
  expect_identical(plan$widget[1], "checkbox")                    # switches first
})

# CLAIM: a seeded run is reproducible and leaves the session's generator alone.
test_that("with_analysis_seed reproduces draws and restores the generator", {
  set.seed(99); before <- .Random.seed
  a <- with_analysis_seed(3, stats::runif(5))
  expect_identical(.Random.seed, before)                          # untouched
  b <- with_analysis_seed(3, stats::runif(5))
  expect_identical(a, b)
  expect_false(identical(a, with_analysis_seed(4, stats::runif(5))))
  x1 <- stats::runif(1); set.seed(99); x2 <- stats::runif(1)
  expect_identical(x1, x2)                                        # the session stream continued as if nothing ran
  # NA = no seeding: the call consumes the session stream as before
  set.seed(5); with_analysis_seed(NA, stats::runif(1)); y1 <- stats::runif(1)
  set.seed(5); stats::runif(1); y2 <- stats::runif(1)
  expect_identical(y1, y2)
  expect_error(with_analysis_seed(1, stop("boom")), "boom")
  expect_identical(.Random.seed, { set.seed(5); stats::runif(2); .Random.seed })
})
