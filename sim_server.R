## ─────────────────────────────────────────────────────────────────────────────
## sim_server.R  -  Simulation tab server logic
## Synthetic MINFLUX localizations (shapes OR image-guided)
## ─────────────────────────────────────────────────────────────────────────────

SIM_SERVER_PARTS <- load_server_parts(file.path("server", "sim"))

sim_server <- function(input, output, session) {
  run_server_parts(SIM_SERVER_PARTS, environment())
}
