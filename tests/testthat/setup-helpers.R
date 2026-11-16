# Source the pure helper files once for every test file. Located relative to this
# file so the tests run from the repo root or from tests/.
.repo_root <- normalizePath(file.path(testthat::test_path(), "..", ".."))
suppressMessages({
  source(file.path(.repo_root, "helpers_shared.R"))
  source(file.path(.repo_root, "helpers_tracking.R"))
  source(file.path(.repo_root, "helpers_smlm.R"))
  source(file.path(.repo_root, "R", "smlm_settings.R"))
})
# Slow Monte-Carlo recovery tests run only when asked for.
skip_if_not_slow <- function()
  testthat::skip_if_not(identical(Sys.getenv("MINFLUXR_SLOW_TESTS"), "true"),
                        "slow test (set MINFLUXR_SLOW_TESTS=true)")
# The source files of the SMLM server: the assembly file and its part files
# (server/smlm/), relative to the repo root. Tests that guard the server's code
# by reading it must read all of them.
smlm_server_files <- function()
  c("server_smlm.R", file.path("server", "smlm",
      sort(list.files(file.path(.repo_root, "server", "smlm"), "^[0-9]{2}_.*[.]R$"))))
