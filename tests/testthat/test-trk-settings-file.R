# The tracking settings file, end to end: save in one session, load in a fresh one.
#
# CLAIM: every setting in TRK_SETTINGS_SPEC survives the round trip. The file is
# written with every setting changed from its default; loading it in a new
# session issues one widget update per setting with the saved value (the
# calibration profile excepted: its values are restored directly). A file from
# before the block existed still loads, restoring only its filters.

# GOAL: a value different from the default for every setting.
.changed_settings <- function() {
  out <- list()
  for (id in names(TRK_SETTINGS_SPEC)) {
    sp <- TRK_SETTINGS_SPEC[[id]]; d <- sp$default
    out[[id]] <- switch(sp$type,
      flag   = !isTRUE(d),
      int    = if (is.na(d)) 7L else d + 1L,
      num    = if (is.na(d)) 3.21 else d + 0.25,
      choice = setdiff(sp$choices, d)[1],
      text   = if (identical(id, "trk_hw_pattern")) "3" else "profile_x")
  }
  out
}

# GOAL: record update*() calls made by the server, keyed by input id. The parts
# run in an environment whose parent is the global one, so functions put there
# are found before shiny's.
.spy_updates <- function() {
  calls <- new.env()
  rec <- function(kind, arg) function(session, inputId, ...) {
    a <- list(...); assign(inputId, list(kind = kind, value = a[[arg]]), envir = calls)
  }
  fns <- list(updateNumericInput = rec("numeric", "value"), updateCheckboxInput = rec("checkbox", "value"),
              updateSelectInput = rec("select", "selected"), updateRadioButtons = rec("radio", "selected"),
              updateSliderInput = rec("slider", "value"), updateTextInput = rec("text", "value"))
  for (n in names(fns)) assign(n, fns[[n]], envir = globalenv())
  list(calls = calls, undo = function() rm(list = names(fns), envir = globalenv()))
}

test_that("a saved tracking settings file restores every setting in a new session", {
  skip_if_no_app()
  load_app()
  changed <- .changed_settings()
  got <- new.env()
  shiny::testServer(tracking_server, {
    session$flushReact()
    set_inputs(session, ui_input_defaults(tracking_ui))
    set_inputs(session, changed)
    # the download lives in a folder the session removes when it ends: copy it now
    path <- tempfile(fileext = ".json")
    file.copy(output$trk_download_filter_settings, path)
    got$path <- path
  })
  path <- got$path
  expect_true(file.exists(path))
  raw <- jsonlite::fromJSON(path)
  expect_identical(raw$settings_version, 2L)
  f <- trk_settings_from_file(raw)
  expect_setequal(f$present, names(TRK_SETTINGS_SPEC))
  for (id in names(changed)) expect_identical(f$settings[[id]], .trk_coerce(id, changed[[id]]), info = id)

  spy <- .spy_updates(); on.exit(spy$undo(), add = TRUE)
  shiny::testServer(tracking_server, {
    session$flushReact()
    set_inputs(session, ui_input_defaults(tracking_ui))
    session$setInputs(trk_load_filter_settings = data.frame(name = basename(path), datapath = path,
                                                            stringsAsFactors = FALSE))
    session$flushReact()
  })
  for (id in setdiff(names(changed), "trk_cal_profile")) {
    got <- spy$calls[[id]]
    expect_false(is.null(got), info = paste(id, "not restored"))
    if (is.null(got)) next
    expect_identical(got$kind, TRK_SETTINGS_SPEC[[id]]$widget, info = id)
    expect_identical(.trk_coerce(id, got$value), .trk_coerce(id, changed[[id]]), info = id)
  }
  expect_null(spy$calls[["trk_cal_profile"]])                   # not reselected (see the loader)
})

test_that("a blank prefactor in the file clears a typed one", {
  raw <- jsonlite::fromJSON(jsonlite::toJSON(list(trk_settings = unclass(trk_settings())),
                                             auto_unbox = TRUE, na = "null"))
  f <- trk_settings_from_file(raw)
  expect_true("trk_crlb_pref" %in% f$present)
  expect_true(is.na(f$settings$trk_crlb_pref))
  expect_identical(f$settings$trk_seed, 1L)
})

test_that("a settings file from before the block loads its filters only", {
  skip_if_no_app()
  load_app()
  old <- tempfile(fileext = ".json")
  writeLines(jsonlite::toJSON(list(efo_rng = c(1000, 50000), trk_split_gaps_enable = TRUE),
                              auto_unbox = TRUE), old)
  spy <- .spy_updates(); on.exit(spy$undo(), add = TRUE)
  shiny::testServer(tracking_server, {
    session$flushReact()
    set_inputs(session, ui_input_defaults(tracking_ui))
    session$setInputs(trk_load_filter_settings = data.frame(name = "old.json", datapath = old,
                                                            stringsAsFactors = FALSE))
    session$flushReact()
    expect_identical(trk_filt$efo, c(1000, 50000))
  })
  expect_true(isTRUE(spy$calls$trk_split_gaps_enable$value))    # the old top-level key still works
  expect_null(spy$calls$trk_uni_min_lag)                        # nothing else touched
})
