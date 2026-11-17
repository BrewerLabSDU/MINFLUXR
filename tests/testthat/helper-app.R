# App-level test support: load the whole app once, start a test session with the
# inputs the UI really starts with, and compare numerical results with stored
# golden fingerprints. Used by test-golden-*.R and test-smoke-*.R.
#
# These tests load every package and source all server files, which takes a few
# seconds. They run by default because they are the safety net for refactoring;
# set MINFLUXR_SKIP_APP_TESTS=true to leave them out.

# GOAL: skip a test when the app-level tests are switched off.
skip_if_no_app <- function()
  testthat::skip_if(identical(Sys.getenv("MINFLUXR_SKIP_APP_TESTS"), "true"),
                    "app-level test (MINFLUXR_SKIP_APP_TESTS=true)")

# GOAL: source the app into the global environment, as Shiny does, once per test
# run. global.R uses relative paths, so the working directory is the repo root
# while it is sourced.
load_app <- function() {
  if (isTRUE(getOption("minfluxr.app_loaded"))) return(invisible(TRUE))
  old <- setwd(.repo_root); on.exit(setwd(old), add = TRUE)
  suppressMessages(suppressWarnings(source("global.R")))
  options(minfluxr.app_loaded = TRUE)
  invisible(TRUE)
}

# ── ui_input_defaults() ──────────────────────────────────────────────────────
# GOAL: the value every static input widget of a UI starts with, read from the
# UI's tag tree. shiny::testServer() starts with all inputs NULL, so server code
# falls back to its `%||%` defaults, which are not always the widget's value. A
# session seeded with these values computes what the app computes.
# Covers numeric, text, colour, checkbox, select, radio, checkbox group, slider
# and numeric range inputs. Action buttons and file inputs are left out (unset
# is their state at start). Inputs created by renderUI are not in the tree.
# RETURNS a named list, id -> value.
ui_input_defaults <- function(ui) {
  out <- list()
  cls <- function(tag) paste(unlist(tag$attribs[names(tag$attribs) == "class"]), collapse = " ")
  att <- function(tag, nm) { v <- tag$attribs[[nm]]; if (is.null(v)) NULL else as.character(v)[1] }
  num <- function(v) if (is.null(v) || !nzchar(v)) NA_real_ else suppressWarnings(as.numeric(v))
  # all descendant tags of one kind (stops at nothing: widgets are small)
  find_tags <- function(x, name) {
    res <- list()
    walk <- function(node) {
      if (inherits(node, "shiny.tag")) {
        if (identical(node$name, name)) res[[length(res) + 1L]] <<- node
        for (ch in node$children) walk(ch)
      } else if (is.list(node)) for (ch in node) walk(ch)
    }
    walk(x); res
  }
  walk <- function(node) {
    if (inherits(node, "shiny.tag")) {
      id <- att(node, "id"); cl <- cls(node)
      handled <- FALSE
      if (!is.null(id)) {
        if (grepl("shiny-input-radiogroup", cl, fixed = TRUE)) {
          ins <- find_tags(node$children, "input")
          sel <- Filter(function(i) !is.null(i$attribs$checked), ins)
          out[[id]] <<- if (length(sel)) att(sel[[1]], "value") else NULL
          handled <- TRUE
        } else if (grepl("shiny-input-checkboxgroup", cl, fixed = TRUE)) {
          ins <- find_tags(node$children, "input")
          sel <- Filter(function(i) !is.null(i$attribs$checked), ins)
          out[[id]] <<- if (length(sel)) vapply(sel, att, "", nm = "value") else NULL
          handled <- TRUE
        } else if (grepl("shiny-numeric-range-input", cl, fixed = TRUE)) {
          ins <- find_tags(node$children, "input")
          out[[id]] <<- vapply(ins, function(i) num(att(i, "value")), numeric(1))
          handled <- TRUE
        } else if (identical(node$name, "select")) {
          ch <- ui_select_choices(node)
          multiple <- !is.null(node$attribs$multiple)
          out[[id]] <<- if (any(ch$selected)) ch$value[ch$selected]
                        else if (!multiple && length(ch$value)) ch$value[1] else NULL
          handled <- TRUE
        } else if (identical(node$name, "textarea")) {
          out[[id]] <<- paste(unlist(node$children), collapse = ""); handled <- TRUE
        } else if (identical(node$name, "input")) {
          type <- att(node, "type") %||% "text"
          if (grepl("js-range-slider", cl, fixed = TRUE)) {
            from <- num(att(node, "data-from")); to <- att(node, "data-to")
            out[[id]] <<- if (is.null(to)) from else c(from, num(to))
          } else if (identical(type, "checkbox")) {
            out[[id]] <<- !is.null(node$attribs$checked)
          } else if (identical(type, "number")) {
            out[[id]] <<- num(att(node, "value"))
          } else if (type %in% c("text", "password")) {
            out[[id]] <<- att(node, "value") %||% ""
          }
          handled <- TRUE
        }
      }
      if (!handled) for (ch in node$children) walk(ch)
    } else if (is.list(node)) for (ch in node) walk(ch)
  }
  walk(ui)
  out
}

# GOAL: the choices of one <select> tag. Shiny writes the <option> elements as
# one HTML string, so they are read from the rendered markup.
# RETURNS list(value = character, selected = logical).
ui_select_choices <- function(select_tag) {
  html <- paste(as.character(htmltools::tagList(select_tag$children)), collapse = "")
  m <- regmatches(html, gregexpr("<option[^>]*>", html))[[1]]
  val <- sub('.*value="([^"]*)".*', "\\1", m)
  val[!grepl('value="', m, fixed = TRUE)] <- ""
  for (ent in list(c("&amp;", "&"), c("&lt;", "<"), c("&gt;", ">"), c("&quot;", "\""), c("&#39;", "'")))
    val <- gsub(ent[1], ent[2], val, fixed = TRUE)
  list(value = val, selected = grepl("\\bselected\\b", m))
}

# GOAL: set a list of inputs on a testServer session in one flush.
set_inputs <- function(session, vals) do.call(session$setInputs, vals)

# ── Golden fingerprints ──────────────────────────────────────────────────────
# GOAL: reduce any result (number, vector, matrix, data.frame, nested list) to a
# small list of numbers that changes when the result changes: length, count of
# finite values, sum, sum of squares, five quantiles and the first values of
# each numeric vector; the sorted value table of each character or logical
# vector. Stored instead of the result itself so the golden files stay small
# and a mismatch names the field that moved.
golden_fingerprint <- function(x, depth = 0L) {
  if (is.null(x)) return(list(type = "NULL"))
  if (is.function(x) || is.environment(x)) return(list(type = class(x)[1]))
  if (inherits(x, "Mclust"))
    x <- x[intersect(c("G", "modelName", "bic", "loglik", "parameters", "classification"), names(x))]
  if (is.data.frame(x)) {
    x <- as.data.frame(x)
    return(list(type = "data.frame", nrow = nrow(x), names = names(x),
                cols = lapply(x, golden_fingerprint, depth = depth + 1L)))
  }
  if (is.list(x)) {
    if (depth > 6L) return(list(type = "list", n = length(x)))
    nm <- names(x) %||% rep("", length(x))
    return(list(type = "list", n = length(x), names = nm,
                items = lapply(x, golden_fingerprint, depth = depth + 1L)))
  }
  if (is.factor(x)) x <- as.character(x)
  if (inherits(x, c("POSIXt", "Date"))) x <- as.numeric(x)
  if (is.numeric(x) || is.complex(x)) {
    if (is.complex(x)) x <- Mod(x)
    v <- as.numeric(x); f <- v[is.finite(v)]
    return(list(type = "numeric", dim = dim(x), n = length(v), n_finite = length(f),
                n_na = sum(is.na(v)), sum = sum(f), sumsq = sum(f^2),
                q = if (length(f)) unname(stats::quantile(f, c(0, .25, .5, .75, 1), names = FALSE, type = 7))
                    else numeric(0),
                head = utils::head(v, 8)))
  }
  if (is.character(x) || is.logical(x)) {
    if (!length(x)) return(list(type = class(x)[1], n = 0L))
    tb <- table(x, useNA = "ifany")
    tb <- tb[order(names(tb))]
    if (length(tb) > 40L)                 # e.g. track ids: keep the count and a digest
      return(list(type = class(x)[1], n = length(x), n_unique = length(tb),
                  digest = digest::digest(list(names(tb), as.integer(tb)), algo = "xxhash64")))
    return(list(type = class(x)[1], n = length(x), values = names(tb), counts = as.integer(tb)))
  }
  list(type = class(x)[1])
}

.golden_dir <- function() file.path(testthat::test_path(), "_golden")

# GOAL: compare a named list of results with the stored golden fingerprints of
# the same name. Numbers must agree to `tolerance` (relative). With
# MINFLUXR_UPDATE_GOLDEN=true the stored file is (re)written instead and the
# test is skipped, so goldens are only ever changed on purpose.
# name   : file stem under tests/testthat/_golden/
# values : named list, one entry per result
expect_golden <- function(name, values, tolerance = 1e-8) {
  stopifnot(is.list(values), !is.null(names(values)), all(nzchar(names(values))))
  fp   <- lapply(values, golden_fingerprint)
  path <- file.path(.golden_dir(), paste0(name, ".rds"))
  if (identical(Sys.getenv("MINFLUXR_UPDATE_GOLDEN"), "true")) {
    dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
    saveRDS(fp, path, version = 2)
    testthat::skip(paste("golden written:", basename(path)))
  }
  if (!file.exists(path))
    testthat::fail(paste0("no golden file for '", name, "'; create it with MINFLUXR_UPDATE_GOLDEN=true"))
  else {
    ref <- readRDS(path)
    testthat::expect_setequal(names(fp), names(ref))
    for (k in intersect(names(ref), names(fp)))
      testthat::expect_equal(fp[[k]], ref[[k]], tolerance = tolerance, info = paste0(name, " :: ", k))
  }
  invisible(fp)
}

# GOAL: evaluate one result with a fixed RNG state; an error becomes a string so
# one broken panel does not hide the others.
.grab <- function(expr) {
  set.seed(1)
  tryCatch(force(expr), error = function(e) paste("ERROR:", conditionMessage(e)))
}

# GOAL: the entries of a result collector that hold an error message.
.golden_errors <- function(res) {
  bad <- Filter(function(v) is.character(v) && length(v) == 1L && startsWith(v, "ERROR:"), res)
  if (!length(bad)) return(character(0))
  paste(names(bad), unlist(bad))
}

# ── Shared simulated datasets ────────────────────────────────────────────────
# GOAL: a small two-population tracking dataset with dark cycles, photon noise
# and blur, written in the app's CSV layout (tim in SECONDS on disk). Fixed seed.
golden_tracks <- function(seed = 11L, n_tracks = 60L, n_locs = 220L) {
  set.seed(seed)
  eco <- pmax(20, stats::rlnorm(5e4, log(60), 0.5))
  sig <- 100.8 / (4 * sqrt(eco))
  Q   <- matrix(c(-0.004, 0.004, 0.006, -0.006), 2, byrow = TRUE)
  df  <- sim_ct_state_tracks(n_tracks, n_locs, D = c(40, 400), Q = Q, t_c = 0.192, sw = 0.011,
                             sigma = sig, R_blur = 0.106, seed = seed)
  n <- nrow(df)
  # spread the tracks over a field so spatial panels (D-map, crop) have extent
  off <- cbind(stats::runif(n_tracks, 0, 4000), stats::runif(n_tracks, 0, 4000))
  k   <- match(df$tid, unique(df$tid))
  df$x <- df$x + off[k, 1]; df$y <- df$y + off[k, 2]
  df$eco <- round(sample(eco, n, TRUE)); df$efo <- df$eco / 0.12e-3 * stats::runif(n, 0.8, 1.2)
  df$cfr <- stats::runif(n, 0, 0.8); df$dcr <- stats::runif(n, 0.2, 0.8)
  df$tim <- (df$tim + k * 50) / 1000                       # seconds on disk, tracks staggered
  df$z <- 0
  df[, c("tid", "tim", "x", "y", "z", "efo", "cfr", "eco", "dcr")]
}

# GOAL: two small two-channel SMLM files with known molecules on a jittered grid
# (channel by DCR), plus one clustered file, as a named store. Fixed seed.
golden_smlm_store <- function(seed = 5L) {
  set.seed(seed)
  one <- function(n_side, spacing, tag) {
    mol <- grid_molecules(n_side, spacing, jitter = spacing / 4)
    df  <- sim_smlm_traces(mol, blink_mean = 1.2, n_loc = 12L, sigma = 4,
                           dcr = rep(c(0.3, 0.62), length.out = nrow(mol)))
    df$mol <- NULL
    n <- nrow(df)
    df$tim <- sort(stats::runif(n, 0, 1800))[rank(df$tid, ties.method = "first")]
    df$efo <- stats::rlnorm(n, log(60000), 0.3); df$cfr <- stats::runif(n, 0, 0.8)
    df$eco <- round(stats::rlnorm(n, log(150), 0.3)); df$itr <- 9L
    df$tid <- as.character(df$tid)
    df$.dataset_tag <- tag
    df
  }
  list(`cellA.csv` = one(9, 70, "aaaa"), `cellB.csv` = one(8, 90, "bbbb"))
}

# ── Test sessions ────────────────────────────────────────────────────────────
# Each helper is called INSIDE shiny::testServer(<server>, { ... }): that
# expression sees the server closure's objects and the test's own variables.

# GOAL: bring a Tracking test session to the state "files uploaded, Plot Data
# pressed", starting from the UI's default inputs. Called inside
# shiny::testServer(tracking_server, { ... }), whose expression sees both the
# server closure's objects and the test's own variables.
trk_session_start <- function(session, files, extra_inputs = list()) {
  session$flushReact()                         # the start-up flush a real session has
  set_inputs(session, ui_input_defaults(tracking_ui))
  if (length(extra_inputs)) set_inputs(session, extra_inputs)
  session$setInputs(trk_file = data.frame(
    name = basename(files), size = file.size(files), type = "text/csv",
    datapath = files, stringsAsFactors = FALSE))
  session$flushReact()
  session$setInputs(trk_plot_data = 1)
  session$flushReact()
}

# GOAL: write the golden tracking dataset as two CSV files (the second is a
# subset, to give Compare Cells and the calibration panel two "cells").
.golden_track_files <- function() {
  load_app()                                   # the simulator lives in R/est_hmm.R
  df <- golden_tracks()
  dir <- file.path(tempdir(), "minfluxr_golden_trk"); dir.create(dir, showWarnings = FALSE)
  f1 <- file.path(dir, "cell_one.csv"); f2 <- file.path(dir, "cell_two.csv")
  utils::write.csv(df, f1, row.names = FALSE)
  keep <- utils::head(unique(df$tid), 25)
  utils::write.csv(df[df$tid %in% keep, ], f2, row.names = FALSE)
  c(f1, f2)
}

# GOAL: bring an SMLM test session to "files loaded, first one active, Plot
# pressed" from the UI's default inputs. Called inside
# shiny::testServer(smlm_server, { ... }). The file input cannot be driven in a
# test session, so the store is filled directly, as the upload observer would.
smlm_session_start <- function(session, store, extra_inputs = list()) {
  env <- session$env                           # the server function's environment
  session$flushReact()
  set_inputs(session, ui_input_defaults(smlm_ui))
  env$smlm_file_store(store)
  session$setInputs(active_file = names(store)[1])
  session$flushReact()
  if (is.null(tryCatch(env$uploaded_df(), error = function(e) NULL))) env$uploaded_df_rv(store[[1]])
  base <- list(enable_clustering = TRUE, enable_two_color = TRUE,
               enable_cluster2 = TRUE, eps2 = 150, minpts2 = 4)
  base[names(extra_inputs)] <- extra_inputs              # a caller's value wins
  set_inputs(session, base)
  session$flushReact()
  session$setInputs(plot_data = 1)
  session$flushReact()
}

# GOAL: start a Simulator test session from the UI's default inputs, with both
# seeds fixed. Called inside shiny::testServer(sim_server, { ... }).
sim_session_start <- function(session, extra_inputs = list()) {
  session$flushReact()
  set_inputs(session, ui_input_defaults(sim_ui))
  set_inputs(session, c(list(sim_fix_seed = TRUE, sim_seed = 42, sim_fix_label_seed = TRUE,
                             sim_label_seed = 7, sim_trk_seed = 3), extra_inputs))
  session$flushReact()
}
