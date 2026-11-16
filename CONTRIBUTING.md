# Contributing to MINFLUXR

Thank you for helping improve MINFLUXR. Bug reports, method suggestions and code are all welcome.
Everyone taking part is expected to follow the [Code of Conduct](CODE_OF_CONDUCT.md).

## Reporting a problem

Open an issue at https://github.com/ViktorMebus/MINFLUXR/issues with:

1. The MINFLUXR version, your R version (`R.version.string`) and operating system.
2. The tab and panel, and what you clicked or changed.
3. What you expected and what happened, including any error text from the R console.
4. If possible, a file that reproduces it. Please do not post unpublished data. A file exported from
   the Simulator tab is usually enough.

## Suggesting an analysis

Describe the question it answers, the method or paper it comes from, and which MINFLUX property it
needs to handle (trace structure, quantised time steps, per-axis precision).

## Contributing code

1. Fork the repository and create a branch from `master`.
2. Put the computation in a plain function in `R/` (`est_*.R` and `trk_*.R` for tracking,
   `smlm_*.R` for SMLM, `sim_*.R` for the simulator). The server code lives in part files under
   `server/<tab>/`, one panel each; keep those for wiring inputs to outputs. `server/README.md`
   explains how the parts fit together.
3. Open each function with a comment saying what it does, its inputs and its output, and comment the
   steps inside it.
4. Add a test in `tests/` that checks the function against simulated data with a known answer. For a
   new panel, also add its result to the golden test of its tab (`tests/testthat/test-golden-*.R`).
5. A new analysis panel needs a method box saying what the analysis is for, what to be aware of and
   how its output can be used, and an entry in `docs/user_guide.md`.
6. Before opening a pull request, run:

       Rscript tools/server_map.R --headers
       Rscript tests/run_tests.R
       Rscript tools/check_dependencies.R

   The first refreshes the interface block at the top of each server part. If the tests report a
   changed golden value, that is a changed result: regenerate it with `MINFLUXR_UPDATE_GOLDEN=true`
   only when the change is intended, and say why in the commit.
7. Add a line under "Unreleased" in `CHANGELOG.md`.

By contributing you agree that your contribution is licensed under the GNU General Public License
v3.0 or later, the licence of this project.
