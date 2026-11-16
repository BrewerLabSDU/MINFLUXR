## ─────────────────────────────────────────────────────────────────────────────
## server/parts.R  -  load and run the part files of a server function
## ─────────────────────────────────────────────────────────────────────────────
##
## PURPOSE
##   smlm_server(), tracking_server() and sim_server() used to be single
##   functions of 4,000 to 14,000 lines. Their bodies now live in ordered part
##   files (server/smlm/, server/tracking/, server/sim/), one panel or one core
##   concern per file. These two functions put them back together.
##
## HOW IT WORKS
##   A part file is NOT a standalone script and defines no function of its own.
##   It holds a run of the statements that used to sit in the server function's
##   body. At app start each part is parsed once (load_server_parts). When a
##   session starts, the server function evaluates every statement of every part,
##   in file-name order, in its own environment (run_server_parts). The result is
##   the same as if the statements were written in the function body:
##     * every reactive, observer and function of a tab shares one environment,
##       so a part can use names defined in any other part of the same tab;
##     * statements run in the same order as before, so observers are created in
##       the same order;
##     * input and output ids are not namespaced.
##   The header of each part lists what it defines for other parts and what it
##   reads from them; tests/testthat/test-server-parts.R checks those lists
##   against the code. See server/README.md.
## ─────────────────────────────────────────────────────────────────────────────

# GOAL: parse the part files of one server function, once, at app start.
# dir : folder holding the parts, relative to the app folder. Parts are the files
#       named NN_<name>.R; they are loaded in file-name order.
# RETURNS a named list (file name -> parsed expressions). Source references are
#   kept, so an error message points into the part file.
load_server_parts <- function(dir) {
  files <- sort(list.files(dir, pattern = "^[0-9]{2}_.*[.]R$", full.names = TRUE))
  if (!length(files)) stop("no server part files found in ", dir)
  stats::setNames(lapply(files, parse, keep.source = TRUE, encoding = "UTF-8"), basename(files))
}

# GOAL: evaluate every statement of every part in `env`, in order.
# parts : result of load_server_parts()
# env   : the calling server function's environment, i.e. environment()
run_server_parts <- function(parts, env) {
  for (part in parts) for (stmt in part) eval(stmt, env)
  invisible(NULL)
}
