## ─────────────────────────────────────────────────────────────────────────────
## tests/run_tests.R  -  entry point for the MINFLUXR validation tests
## ─────────────────────────────────────────────────────────────────────────────
##
## PURPOSE
##   1. Recovery tests for the estimators and SMLM computations in R/ and the
##      helper files. Each test states a claim the methods doc makes and checks
##      it against simulated ground truth, so the claim cannot silently stop
##      being true.
##   2. Golden-output tests (test-golden-*.R): for fixed simulated data, every
##      analysis panel of the three tabs must return the numbers stored in
##      tests/testthat/_golden/. They guard refactoring: a change of structure
##      must not change a result.
##   3. Smoke tests (test-smoke-servers.R): every output of the three server
##      functions keeps its render status after all panels have been triggered.
##   Groups 2 and 3 drive the server functions in shiny::testServer (no browser).
##
## RUN
##   From the repo root:  Rscript tests/run_tests.R
##   MINFLUXR_SLOW_TESTS=true      also run the slow Monte-Carlo tests
##   MINFLUXR_SKIP_APP_TESTS=true  leave out groups 2 and 3 (they load the whole app)
##   MINFLUXR_UPDATE_GOLDEN=true   rewrite the golden files; only for a deliberate
##                                 change of method, and say why in the commit
## ─────────────────────────────────────────────────────────────────────────────
library(testthat)
root <- if (file.exists("global.R")) "." else ".."          # repo root or tests/
test_dir(file.path(root, "tests", "testthat"), reporter = "summary", stop_on_failure = TRUE)
