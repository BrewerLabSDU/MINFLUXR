# Structure checks for the part files of the three server functions (server/).
# They parse the code and run nothing, so they are fast.

# GOAL: the mapper's functions, loaded with the repo root as working directory.
.server_map_env <- function() {
  env <- new.env()
  sys.source(file.path(.repo_root, "tools", "server_map.R"), envir = env)
  env
}

# CLAIM: no part has grown back into a monolith. 1,600 lines leaves room for the
# one statement that is that long on its own (the 3D plot render).
test_that("every server part stays under 1,600 lines", {
  old <- setwd(.repo_root); on.exit(setwd(old), add = TRUE)
  map <- .server_map_env()
  for (tab in names(map$SERVER_CLOSURES)) {
    n <- vapply(map$server_part_interfaces(tab), `[[`, 1L, "n_lines")
    expect_true(all(n <= 1600L), info = paste(tab, paste(names(n)[n > 1600L], collapse = ", ")))
    expect_gt(length(n), 1L)
  }
})

# CLAIM: no closure-level name is defined twice within a tab. All parts of a tab
# share one environment, so a second definition silently replaces the first and
# the first never runs (this happened in server_smlm.R before).
test_that("no server-level name is defined twice within a tab", {
  old <- setwd(.repo_root); on.exit(setwd(old), add = TRUE)
  map <- .server_map_env()
  for (tab in names(map$SERVER_CLOSURES))
    expect_identical(map$closure_duplicates(map$server_statement_map(tab)), character(0), info = tab)
})

# CLAIM: the DEFINES / READS block at the top of each part matches the code. It
# is the part's interface; if this fails after an edit, refresh the blocks with
#   Rscript tools/server_map.R --headers
# and look at what changed: a new READS line is a new dependency between panels.
test_that("each part's interface header matches its code", {
  old <- setwd(.repo_root); on.exit(setwd(old), add = TRUE)
  map <- .server_map_env()
  for (tab in names(map$SERVER_CLOSURES)) {
    ifs <- map$server_part_interfaces(tab)
    for (f in names(ifs)) {
      x <- readLines(file.path(map$SERVER_CLOSURES[[tab]]$parts, f), warn = FALSE, encoding = "UTF-8")
      a <- match(map$.PART_IFACE_OPEN, x); b <- match(map$.PART_IFACE_CLOSE, x)
      expect_false(is.na(a) || is.na(b), info = paste(tab, f, "has no interface block"))
      if (!is.na(a) && !is.na(b))
        expect_identical(x[a:b], map$.part_interface_lines(ifs[[f]]), info = paste(tab, f))
    }
  }
})

# CLAIM: a part reads only `input`, `output`, `session`, names defined in a part
# of the same tab, and global functions. In particular the tabs do not reach
# into each other.
test_that("the three tabs share no server-level name", {
  old <- setwd(.repo_root); on.exit(setwd(old), add = TRUE)
  map <- .server_map_env()
  nm <- lapply(names(map$SERVER_CLOSURES), function(tab) attr(map$server_statement_map(tab), "defined"))
  # plot_trigger exists in both the SMLM and the tracking tab, as separate objects
  # with different meanings (docs/architecture_review.md, 2.3 item 6).
  expect_identical(intersect(nm[[1]], nm[[2]]), "plot_trigger")
  expect_identical(intersect(nm[[1]], nm[[3]]), character(0))
  expect_identical(intersect(nm[[2]], nm[[3]]), character(0))
})
