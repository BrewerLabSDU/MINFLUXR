## ─────────────────────────────────────────────────────────────────────────────
## trk_geometry.R  -  tracking: cell-geometry labels
## ─────────────────────────────────────────────────────────────────────────────
##
## PURPOSE
##   The display label of the geometry selector, the dimensions line of a
##   geometry, and the colour of a displacement class.
##
##   Plain functions and constants with no Shiny state, moved unchanged out of
##   the tracking_server() closure (server_tracking.R). The server gathers the
##   inputs and calls these.
##
## UNITS: nm and ms; D in nm^2/ms unless a name says um2_s.
## ─────────────────────────────────────────────────────────────────────────────

# Human-readable name of the shared geometry, for panels that only echo it.
.geom_label <- function(g) switch(g %||% "none",
  none = "None - unbounded", sphere = "Sphere", rod = "Rod / spherocylinder",
  ellipsoid = "Ellipsoid", membrane = "Flat membrane", filament = "Filament / linear track", g)

# Labels for the current geometry's size parameters (used in the summary).
.disp_geom_dims <- function(geom, p) switch(geom,
  sphere    = sprintf("R = %.0f nm", p$R),
  rod       = sprintf("r = %.0f nm, L = %.0f nm", p$r, p$Lcyl),
  ellipsoid = sprintf("semi-axes %.0f / %.0f / %.0f nm", p$a, p$b, p$c),
  membrane  = sprintf("tilt = %.0f°", p$mtilt),
  free      = "unbounded", "")

# Blue = sub-diffusive (confined/geometric), grey = free-consistent, red = super-
# diffusive (directed) - matching the Speed colouring (blue slow, red fast).
.disp_class_col <- function(cls) c(`sub-diffusive` = "#2c7fb8",
                                   `free-consistent` = "#9e9e9e",
                                   `super-diffusive` = "#e31a1c",
                                   `confined` = "#2c7fb8",
                                   `free-on-geometry` = "#9e9e9e",
                                   `directed` = "#e31a1c")[cls]

