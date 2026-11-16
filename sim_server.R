## ─────────────────────────────────────────────────────────────────────────────
## sim_server.R  -  Simulation tab server logic
## Synthetic MINFLUX localizations (shapes OR image-guided) and synthetic tracks
## ─────────────────────────────────────────────────────────────────────────────
##
## This file only assembles sim_server(). Its body is in server/sim/, evaluated in
## this order in the server function's environment (server/parts.R,
## server/README.md):
##
##   01_view_scalebars.R                view state and scale bars
##   02_smlm_simulation.R               SMLM simulation
##   03_real_run_calibration.R          calibrate the acquisition model from a real run
##   04_compare_strategies.R            compare labelling strategies
##   05_smlm_outputs.R                  SMLM simulator outputs
##   06_tracking_simulation.R           tracking simulation
##
## The generators themselves hold no session state and live in R/sim_*.R.
## ─────────────────────────────────────────────────────────────────────────────

SIM_SERVER_PARTS <- load_server_parts(file.path("server", "sim"))

sim_server <- function(input, output, session) {
  run_server_parts(SIM_SERVER_PARTS, environment())
}
