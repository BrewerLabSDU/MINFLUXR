## ─────────────────────────────────────────────────────────────────────────────
## trk_io.R  -  tracking: loading and table coercion
## ─────────────────────────────────────────────────────────────────────────────
##
## PURPOSE
##   Reading one tracking file into the canonical table, the labels of the
##   dataset selector, the vertex count of a TCP pattern name, and type coercion
##   for the Compare Cells table.
##
##   Plain functions and constants with no Shiny state, moved unchanged out of
##   the tracking_server() closure (server_tracking.R). The server gathers the
##   inputs and calls these.
##
## UNITS: nm and ms; D in nm^2/ms unless a name says um2_s.
## ─────────────────────────────────────────────────────────────────────────────

# Loader for one file → data.frame (or NULL).  Returns attr "json_itr".
load_one_trk <- function(p, nm) {
  d <- NULL; jitr <- NA_integer_
  if (grepl("\\.csv(\\.gz)?$", nm)) {
    d <- data.table::fread(p, showProgress = FALSE, nThread = max(1, parallel::detectCores()-1)) |>
         as.data.frame()
    # A CSV arrives with some axial factor already applied; return z to raw so
    # trk_raw_df() applies the user's calibration once (it was applied twice).
    d <- strip_csv_z_calibration(d)
    d <- canonicalize_minflux_df(d)
  } else if (grepl("\\.(qs|fst)$", nm)) {
    # read_chunk_file reports a missing 'qs' rather than throwing, so a machine
    # without it degrades to "this file was skipped" instead of a stack trace.
    d <- read_chunk_file(p)
    if (is.null(d)) return(NULL)
  } else if (grepl("\\.rds$", nm)) {
    d <- readRDS(p)
    if (!is.data.frame(d)) return(NULL)
    d <- canonicalize_minflux_df(d)
  } else if (grepl("\\.json$", nm)) {
    # z_scale = 1: chunks store RAW z. The axial calibration is applied once,
    # in trk_raw_df(), so every format lands on the same footing.
    d <- read_minflux_json_as_df(p, keep_final_only = TRUE, keep_valid_only = TRUE,
                                 z_scale = 1.0)
    jitr <- if ("itr" %in% names(d) && nrow(d) > 0) max(d$itr, na.rm = TRUE) else NA_integer_
  } else if (grepl("\\.npy$", nm)) {
    d <- read_minflux_npy_as_df(p, keep_final_only = TRUE, keep_valid_only = TRUE,
                                z_scale = 1.0)
    jitr <- if ("itr" %in% names(d) && nrow(d) > 0) max(d$itr, na.rm = TRUE) else NA_integer_
  } else {
    return(NULL)
  }
  attr(d, "json_itr") <- jitr
  d
}

# Choices for the Active dataset selector: every file, plus the merged view.
.trk_view_choices <- function(store) {
  nm <- names(store)
  ch <- stats::setNames(nm, nm)
  if (length(nm) >= 2L)
    ch <- c(ch, stats::setNames("__merged__", sprintf("All %d loaded files, merged", length(nm))))
  ch
}

# ── TSA: MINFLUX sequence (.json) import ────────────────────────────────────
# Parses an Imspector imaging-sequence JSON. The acquisition cycle time is still
# recovered from the DATA (the robust cluster method) - the sequence is used to
# cross-check the scan geometry and, most usefully, to set the gap-split
# threshold from the hardware's own maxOffTime ("particle lost" timeout).
.pattern_vertices <- function(p) {
  if (is.null(p) || is.na(p)) return(NA_integer_)
  switch(tolower(p),
         hexagon = 6L, square = 4L, triangle = 3L, octahedron = 6L,
         zline = 2L, zline2 = 2L, line = 2L,
         NA_integer_)
}

# One value typed for its column, or NULL when it does not fit (the caller then
# leaves the table unchanged and says so).
# An empty value clears a setting back to blank (= Acquisition Parameters).
.cmp_coerce <- function(col, val) {
  if (is.null(val) || (length(val) == 1L && is.na(val))) val <- ""
  if (col %in% c("L_nm", "a", "D_w")) {
    if (!nzchar(trimws(as.character(val)))) return(NA_real_)
    v <- suppressWarnings(as.numeric(val))
    if (!length(v) || !is.finite(v) || v < 0 || (col != "D_w" && v == 0)) return(NULL)
    v
  } else if (col == "pattern") {
    v <- .tcp_pattern_name(val)                      # names, short forms, or 3 / 4 / 6
    if (is.na(v)) return(NULL)
    v
  } else trimws(as.character(val))
}

