## ─────────────────────────────────────────────────────────────────────────────
## server/tracking/06_filter_pipeline.R  -  THE FILTER PIPELINE
## ─────────────────────────────────────────────────────────────────────────────
## Part of tracking_server(). Not a standalone script: its statements are evaluated in
## the server function's environment, in file-name order (server/parts.R).
##
## The stored filter settings (trk_filt) and the stages raw -> split ->
## end-trim -> continuity gate -> metadata -> length -> duration -> spatial
## crop.
##
## INTERFACE (generated: Rscript tools/server_map.R --headers; do not edit by hand)
##   DEFINES for other parts:
##     .a_in_use, .FILT_IDS, .filt_value, bump_clean, trk_dur_filtered_d,
##     trk_end_trimmed_d, trk_filt, trk_filter_params, trk_gap_split,
##     trk_gap_split_d, trk_len_filtered, trk_len_filtered_d,
##     trk_meta_filtered, trk_meta_filtered_d, trk_quality_filtered_d,
##     trk_spatial_filtered
##   READS from other parts:
##     01_state_sigma: trk_spatial_filter
##     04_datasets_upload: trk_raw_df
##     07_cycle_time_sequence: trk_adaptive_gap_cycles, trk_cycle_time,
##                             trk_seq_params
## END INTERFACE
## ─────────────────────────────────────────────────────────────────────────────

# ---- Three-stage filtering pipeline ----------------------------------------
# Pipeline: Stage 1 (metadata) → Stage 2 (trace ID length) → Stage 3 (duration).
# Each stage uses metrics computed from the PREVIOUS stage's output so that
# upstream filters propagate correctly (e.g. metadata filtering shrinks trace
# localization counts, which then drives the length filter).
# All filter inputs use NULL pass-through (input not yet rendered → all data
# passes) and a stale-range guard (filter range outside data range → all data
# passes).  TID values are coerced to character to guard against type mismatches
# when readr parses numeric TIDs from tid_summary.csv.

# Stage 1: per-localization metadata filters (efo, cfr, dcr, tim, eco, dwell, fluo)
# Stage 4 (metadata filter): now runs AFTER track cleaning (split → end-trim →
# quality gate), so efo/cfr/dcr/… filters apply to the cleaned tracks. Track
# cleaning is done on the RAW timing (dark-cycle detection is not distorted by
# metadata-removed localizations), which is the correct order.
trk_meta_filtered <- reactive({
  df <- trk_quality_filtered()
  shiny::validate(shiny::need(all(c("x","y","tid","tim") %in% names(df)),
                              "Data must contain x, y, tid, tim."))
  filter_ranges(df, trk_filt_ranges())                    # R/est_clean.R
})

# ===========================================================================
#  FILTER SETTINGS - stored as the user set them, independent of any dataset
# ===========================================================================
# A range filter ("EFO between 50k and 200k") applies to every loaded dataset
# alike, so it is kept here as the user entered it. The pipeline reads these
# values; the widgets only DISPLAY them, clamped to the active dataset's range.
# Storing the clamped widget values instead would silently change the filter
# every time a dataset with a narrower range was shown. Reset Ranges clears them.
.FILT_IDS <- c(efo = "efo_rng", cfr = "cfr_rng", dcr = "dcr_rng", tim = "tim_rng",
               eco = "eco_rng", dwell = "dwell_rng", fluo = "fluo_rng",
               tid_len = "trk_tid_len", dur = "trk_dur_filter")
trk_filt    <- reactiveValues()                    # key -> c(lo, hi); NULL = no filter
.filt_shown <- new.env(parent = emptyenv())        # key -> list(val, bounds) last displayed

# The value a widget shows: the stored filter clamped into the dataset's bounds,
# or the full bounds when there is no filter. Remembered, so the echo of our own
# render is not mistaken for a user edit.
.filt_value <- function(key, bounds) {
  sp  <- isolate(trk_filt[[key]])
  val <- if (is.null(sp)) bounds
         else c(min(max(sp[1], bounds[1]), bounds[2]), max(min(sp[2], bounds[2]), bounds[1]))
  assign(key, list(val = val, bounds = bounds), envir = .filt_shown)
  val
}

# A user edit becomes the stored filter; covering the full range clears it.
for (.k in names(.FILT_IDS)) local({
  key <- .k; id <- .FILT_IDS[[.k]]
  observeEvent(input[[id]], {
    v <- suppressWarnings(as.numeric(input[[id]]))
    if (length(v) != 2L || any(!is.finite(v))) return(invisible(NULL))
    sh <- if (exists(key, envir = .filt_shown, inherits = FALSE)) get(key, envir = .filt_shown) else NULL
    if (!is.null(sh) && isTRUE(all.equal(v, sh$val))) return(invisible(NULL))   # our own render
    b <- sh$bounds
    trk_filt[[key]] <- if (!is.null(b) && v[1] <= b[1] && v[2] >= b[2]) NULL else sort(v)
    if (!is.null(sh)) assign(key, list(val = v, bounds = b), envir = .filt_shown)
  }, ignoreInit = TRUE)
})
rm(.k)

# The prefactor a in use: the value entered (or set by a calibration profile),
# else the ideal-estimator value of the TCP pattern, which is uncalibrated.
.a_in_use <- function(pattern = input$trk_tcp_pattern) {
  v <- suppressWarnings(as.numeric(input$trk_crlb_pref))
  if (length(v) && is.finite(v[1]) && v[1] > 0) v[1]
  else .tcp_geom(pattern %||% TCP_PATTERN_DEFAULT)$a_ideal
}

# The per-localization range filters, as filter_ranges() takes them.
trk_filt_ranges <- reactive(stats::setNames(
  lapply(c("efo", "cfr", "dcr", "tim", "eco", "dwell", "fluo"), function(k) trk_filt[[k]]),
  c("efo", "cfr", "dcr", "tim", "eco", "dwell", "fluo")))

# Every filter setting, for running the whole pipeline on any dataset
# (filter_pipeline() in R/est_clean.R): Plot Data's rechunk and Compare Cells.
# t_c / sw belong to the dataset and are filled in by the caller.
trk_filter_params <- function() isolate(list(
  split   = list(enable = isTRUE(input$trk_split_gaps_enable), adaptive = isTRUE(input$trk_split_adaptive),
                 gap_cycles = input$trk_split_gap_cycles, min_len = input$trk_split_min_len %||% TRK_DEFAULTS$trk_split_min_len,
                 maxoff_ms = tryCatch({ sp <- trk_seq_params()
                                        if (!is.null(sp) && isTRUE(sp$ok) && is.finite(sp$maxOffTime_ms))
                                          sp$maxOffTime_ms else 3.0 }, error = function(e) 3.0)),
  trim    = list(enable = isTRUE(input$trk_endtrim_enable), dark = isTRUE(input$trk_endtrim_dark),
                 dark_cycles = input$trk_endtrim_dark_cycles %||% TRK_DEFAULTS$trk_endtrim_dark_cycles, efo = isTRUE(input$trk_endtrim_efo),
                 efo_frac = suppressWarnings(as.numeric(input$trk_endtrim_efo_frac %||% TRK_DEFAULTS$trk_endtrim_efo_frac))),
  gate    = list(enable = isTRUE(input$trk_quality_enable),
                 min_score = suppressWarnings(as.numeric(input$trk_quality_min_integ %||% TRK_DEFAULTS$trk_quality_min_integ)),
                 half_life = suppressWarnings(as.numeric(input$trk_quality_darktol %||% TRK_DEFAULTS$trk_quality_darktol))),
  ranges  = trk_filt_ranges(), tid_len = trk_filt$tid_len, dur = trk_filt$dur))

# ===========================================================================
#  TRACKING FILTER PIPELINE  (read top-to-bottom; each stage feeds the next)
#    raw → split → end-trim → continuity gate → metadata → length → duration
#  The cycle-time + TSA machinery lives further down, "CYCLE-TIME & TSA SUPPORT".
# ===========================================================================

# Cleaning-apply trigger. Split / end-trim / continuity-gate are HEAVY (per-track
# passes over ALL localizations), so they must not recompute on every slider
# touch. They read their parameters via isolate() and recompute only when this
# bumps - on Plot Data, a new upload, or Filter & Rechunk. (trk_raw_df is still a
# live dependency, so a new manifest also refreshes them.) Metadata / length /
# duration filters stay LIVE (they read the cached cleaning output).
clean_trigger <- reactiveVal(0L)
bump_clean    <- function() clean_trigger(isolate(clean_trigger()) + 1L)

# Stage 1 (gap splitting) - the FIRST cleaning stage, on the RAW data, so dark-
# cycle detection sees the true acquisition timing (metadata filtering, which
# can remove localizations and widen gaps, now happens LATER). Splits tracks
# wherever an interval is longer than the threshold (in acquisition cycles). New
# TIDs "<orig>_s1", "<orig>_s2", … In ADAPTIVE (Vogler/Bela) mode the threshold
# is data-driven (median+SD+1 of the total-cycle counts, capped by maxOffTime)
# and segments shorter than the min-length are dropped.
trk_gap_split <- reactive({
  clean_trigger()                        # apply only on Plot Data / upload / rechunk
  df <- trk_raw_df()
  if (!isTRUE(isolate(input$trk_split_gaps_enable)) || nrow(df) == 0) return(df)
  adaptive   <- isTRUE(isolate(input$trk_split_adaptive))
  gap_cycles <- if (adaptive) isolate(trk_adaptive_gap_cycles())
                else suppressWarnings(as.numeric(isolate(input$trk_split_gap_cycles)))
  modal_dt   <- isolate(trk_cycle_time())
  if (!is.finite(gap_cycles) || gap_cycles <= 0 ||
      !is.finite(modal_dt) || modal_dt <= 0) return(df)
  gap_thresh <- gap_cycles * modal_dt   # ms
  # Min-length applies in BOTH modes now (was adaptive-only): drop short fragments
  # so a fixed 7-cycle split doesn't manufacture 2–3-loc stubs and flood the plot.
  min_len <- max(1L, as.integer(isolate(input$trk_split_min_len %||% TRK_DEFAULTS$trk_split_min_len)))
  # The split itself (idempotent on an existing "_s<N>" suffix, radix-sorted,
  # short fragments dropped) is split_gaps() in R/est_clean.R, shared with the
  # per-cell comparison.
  split_gaps(df, gap_thresh, min_len)
})

# Stage 2.45: end-trim (optional; pass-through when disabled).
# Photobleaching degrades the ENDS of a track: fewer photons (lower EFO), worse
# precision, and blinky dark-gap-separated localizations. These inflate Rmax /
# MSD tails and produce spurious "exports at track ends". Peel leading & trailing
# localizations that are reached across a dark gap (≥ the dark-cycle threshold) or
# are dim (EFO < frac × the track's median), from each end inward, until a solid,
# well-connected localization is reached. Idempotent - a second pass stops at the
# same boundary. Runs per segment (after gap-split), so it also cleans mid-track
# segment ends, not just the bleaching tail.
trk_end_trimmed <- reactive({
  df <- trk_gap_split()                  # inherits the clean_trigger recompute
  if (!isTRUE(isolate(input$trk_endtrim_enable)) || nrow(df) == 0 || !"tim" %in% names(df))
    return(df)
  # Length floor: gap splitting promises segments of at least `min_len`
  # localizations, so trimming stops there rather than breaking that promise.
  floor_len <- if (isTRUE(isolate(input$trk_split_gaps_enable)))
                 max(1L, as.integer(isolate(input$trk_split_min_len) %||% TRK_DEFAULTS$trk_split_min_len)) else 1L
  trim_track_ends(df, tryCatch(isolate(trk_cycle_time()), error = function(e) NA_real_),
                  use_dark    = isTRUE(isolate(input$trk_endtrim_dark)),
                  dark_cycles = isolate(input$trk_endtrim_dark_cycles %||% TRK_DEFAULTS$trk_endtrim_dark_cycles),
                  use_efo     = isTRUE(isolate(input$trk_endtrim_efo)),
                  efo_frac    = suppressWarnings(as.numeric(isolate(input$trk_endtrim_efo_frac %||% TRK_DEFAULTS$trk_endtrim_efo_frac))),
                  floor_len   = floor_len)                 # R/est_clean.R
})

# Stage 2.5: track-CONTINUITY gate (optional; pass-through when disabled).
# MINFLUX tracks vary in quality: a well-tracked molecule is localized nearly
# every cycle, whereas a dark-cycle-heavy track has long gaps → larger positional
# uncertainty, mis-linking risk, a bleaching-degraded tail. Rather than a hard
# integrated/dark cut, each step is weighted by its dark-cycle count with a
# tunable tolerance (half-life): a few 1–2-cycle gaps barely hurt, long dark gaps
# do. Keeps tracks whose mean continuity score ≥ the threshold.
trk_quality_filtered <- reactive({
  df <- trk_end_trimmed()                # inherits the clean_trigger recompute
  if (!isTRUE(isolate(input$trk_quality_enable)) || nrow(df) == 0 || !"tim" %in% names(df))
    return(df)
  gate_continuity(df, tryCatch(isolate(trk_cycle_time()), error = function(e) NA_real_),
                  min_score = suppressWarnings(as.numeric(isolate(input$trk_quality_min_integ) %||% TRK_DEFAULTS$trk_quality_min_integ)),
                  half_life = suppressWarnings(as.numeric(isolate(input$trk_quality_darktol) %||% TRK_DEFAULTS$trk_quality_darktol)))
})

# Stage 5: trace ID length - localization count per TID, measured AFTER cleaning
# (split/trim/quality) and the metadata filter.
trk_len_filtered <- reactive({
  filter_tid_length(trk_meta_filtered(), trk_filt$tid_len)   # R/est_clean.R
})

# Stage 4: trace duration - max(tim)-min(tim) per TID measured AFTER stages 1,2 and 3.
trk_dur_filtered <- reactive({
  filter_duration(trk_len_filtered(), trk_filt$dur)          # R/est_clean.R
})

# Stage 5: spatial crop - keeps TIDs where >= pct_thr fraction of their
# localizations fall inside the user-defined view box. Passes through unchanged
# when no spatial filter is active. This is the LAST stage; trk() and the plot
# read from here.
trk_spatial_filtered <- reactive({
  df <- trk_dur_filtered()
  sf <- trk_spatial_filter()
  if (is.null(sf) || nrow(df) == 0) return(df)
  in_box <- df$x >= sf$x[1] & df$x <= sf$x[2] &
            df$y >= sf$y[1] & df$y <= sf$y[2]
  tids   <- unique(df$tid)
  keep   <- vapply(tids, function(tid) {
    idx <- df$tid == tid
    mean(in_box[idx], na.rm = TRUE) >= sf$pct_thr
  }, logical(1))
  df[df$tid %in% tids[keep], , drop = FALSE]
})

# Debounced views of the filter chain - UIs and the histogram use these so they
# don't recompute on every slider tick; they settle with fire_frames (400 ms).
# The length slider auto-detects from the split output (gap_split_d); the
# duration slider auto-detects from the length-filtered output (len_filtered_d) -
# so each slider's range reflects the LIVE previous stage (incl. split segments).
trk_gap_split_d        <- debounce(trk_gap_split,        350)
trk_end_trimmed_d      <- debounce(trk_end_trimmed,      350)
trk_quality_filtered_d <- debounce(trk_quality_filtered, 350)
trk_meta_filtered_d    <- debounce(trk_meta_filtered,    350)
trk_len_filtered_d     <- debounce(trk_len_filtered,     350)
trk_dur_filtered_d     <- debounce(trk_dur_filtered,     350)
