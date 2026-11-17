## ─────────────────────────────────────────────────────────────────────────────
## server/tracking/01_state_sigma.R  -  session state and the app-wide sigma_loc
## ─────────────────────────────────────────────────────────────────────────────
## Part of tracking_server(). Not a standalone script: its statements are evaluated in
## the server function's environment, in file-name order (server/parts.R).
##
## The reactive values and flags of the tab, the run gates of TSA and speed,
## and trk_sigma_loc_eff, the noise floor every short-lag estimator subtracts.
##
## INTERFACE (generated: Rscript tools/server_map.R --headers; do not edit by hand)
##   DEFINES for other parts:
##     .continuity_scores, .sigma_stamp, .tids_per_chunk, current_tid_idx,
##     is_playing, rechunk_request, rv, slider_guard, speed_ran, tick,
##     trace_map, trk_chunk_pref, trk_raw_stats_rv, trk_rechunk_bbox,
##     trk_rechunk_stats_rv, trk_sigma_loc_cve, trk_sigma_loc_eff,
##     trk_sigma_loc_photon, trk_sigma_loc_rv, trk_sigma_loc_src,
##     trk_spatial_filter, trk_zoom_ranges, tsa_ran
##   READS from other parts:
##     06_filter_pipeline: .a_in_use, bump_clean
## END INTERFACE
## ─────────────────────────────────────────────────────────────────────────────

# --- Reactives / state ---
rv <- reactiveValues(
  prev_tid_idx          = NULL,
  play_gen              = 0L,   # global playback generation (simultaneous mode)
  manifest_dir          = NULL, # current working chunks (may be filtered)
  original_manifest_dir = NULL, # frozen at upload - the rechunk always reads from here
  last_rechunk_dir      = NULL, # previous rechunk output, deleted once superseded
  bigfile_notice        = FALSE,
  chunking              = FALSE,
  trk_diag              = NULL  # list(data_clusters, img_blobs) in data coords
)
# ── Chunk size actually used ────────────────────────────────────────────────
# GOAL: turn the "TIDs per chunk" input into the value the chunker gets, for
# both the upload path and Filter & Rechunk, so the two cannot disagree.
#
# A DELIBERATELY SMALL SETTING IS HONOURED. Small chunks are the whole point of
# the control: playback animates one chunk at a time, so a modest machine that
# cannot render 100 tracks needs 30. Only a blank or non-numeric entry falls
# back to the default. The floor is TPC_MIN rather than the default, so a value
# below it is clamped up by a little instead of being thrown back to 100.
#
# WHY THERE IS A FLOOR AT ALL: each chunk becomes one file on disk, so a value
# near 1 would emit thousands of files for a large dataset and make every
# manifest read slow.
TPC_MIN <- 5L
.tids_per_chunk <- function() {
  v <- suppressWarnings(as.integer(isolate(input$trk_tids_per_chunk)))
  if (!length(v) || !is.finite(v[1]) || v[1] < 1L) return(TIDS_PER_CHUNK_DEFAULT)   # unset / blank / garbage
  max(TPC_MIN, v[1])
}

# ── Which chunk the selector lands on when it is rebuilt ────────────────────
# "__all__" = animate every TID; "first" = land on chunk 1.
#
# Rebuilding the selector (upload, Filter & Rechunk, Reset ranges) has to pick a
# value, and the previous selection cannot simply be carried over: re-chunking at
# a different size renumbers the chunks, so "chunk 5" may not exist afterwards and
# would not hold the same tracks if it did. What IS meaningful is whether the user
# was scoped to a chunk at all, so that is what is preserved - someone who was on a
# chunk stays on a chunk (the first one, which always exists), and someone on
# "All TIDs" stays on all. Landing a machine that could not render everything back
# on "All TIDs" is what this avoids.
trk_chunk_pref <- reactiveVal("__all__")

# Plot Data asks for a rechunk through this counter rather than calling the
# routine directly, so the 260-line filter/rechunk observer stays exactly where
# it is and keeps its own progress bar and error handling.
rechunk_request <- reactiveVal(0L)

is_playing         <- reactiveVal(FALSE)
tick               <- reactiveVal(0L)
slider_guard       <- reactiveVal(0L)
trace_map          <- reactiveVal(NULL)
current_tid_idx    <- reactiveVal(1L)
# Frozen at upload time; never overwritten by chunk-switching or Filter & Rechunk.
# Provides stable grey-bar reference for the metadata histogram.
trk_raw_stats_rv      <- reactiveVal(NULL)
# Updated after each Filter & Rechunk; holds stats across ALL remaining TIDs
# (not just the current chunk).  NULL means no rechunk has happened yet.
# Filter UIs use this as their value when non-NULL; Reset Ranges clears it.
trk_rechunk_stats_rv  <- reactiveVal(NULL)
# Full spatial extent of all remaining TIDs after rechunk.
# list(x = c(xmin, xmax), y = c(ymin, ymax)) - used to lock plot axes.
trk_rechunk_bbox      <- reactiveVal(NULL)
# Current plot axis ranges captured from plotly_relayout - used by "Crop to view".
# list(x = c(x0, x1), y = c(y0, y1)) or NULL before first zoom interaction.
trk_zoom_ranges    <- reactiveVal(NULL)
# Active spatial box filter - NULL means no spatial filter is applied.
# list(x = c(x0, x1), y = c(y0, y1), pct_thr = numeric[0-1])
trk_spatial_filter <- reactiveVal(NULL)

# "TSA has been run at least once" gate. The Integrated/Dark bundle is
# eventReactive-gated on trk_tsa_run, but the sibling display outputs (Δt
# histogram, stats, quality verdict, per-track progression) read dt_vectors()
# directly - this flag keeps them dormant until Run is pressed. Cleared when TSA
# is disabled; stays set across filter changes (re-press Run to refresh, same as
# the other run-gated analyses).
tsa_ran <- reactiveVal(FALSE)
# Run TSA (and Estimate σ_loc below) also APPLY the current cleaning settings, so
# the raw→analyzed diagnostics reflect split/trim/continuity without a separate
# Plot Data press. bump_clean() is a no-op when cleaning is unchanged.
observeEvent(input$trk_tsa_run, { bump_clean(); tsa_ran(TRUE) })
observeEvent(input$trk_cve_run, { bump_clean() }, ignoreInit = TRUE)
observeEvent(input$trk_show_tsa, {
  if (!isTRUE(input$trk_show_tsa)) tsa_ran(FALSE)
}, ignoreInit = TRUE)

# Same run-gate for the Instantaneous Speed section: both the TID selector
# (which ranks all tracks by export strength - heavy) and the plot stay dormant
# until Run is pressed, then remain interactive on TID selection.
speed_ran <- reactiveVal(FALSE)
observeEvent(input$trk_speed_run, { speed_ran(TRUE) })
observeEvent(input$trk_show_speed, {
  if (!isTRUE(input$trk_show_speed)) speed_ran(FALSE)
}, ignoreInit = TRUE)

# ── App-wide localization precision σ_loc (the short-lag noise floor) ────────
# Two independent sources feed the app-wide value, in priority order:
#   1. manual override (input$trk_sigma_loc_nm > 0)
#   2. the photon-budget estimate, L_TCP/(4*sqrt(eco))  (trk_sigma_loc_photon)
#   3. the blur-aware CVE estimate                      (trk_sigma_loc_cve)
#   4. the unified MSD fit's intercept                  (trk_sigma_loc_rv)
# Every short-lag tool (D-map, JDD, confinement, HMM, directionality, speed)
# reads trk_sigma_loc_eff() so the noise floor is consistent everywhere.
#
# WHY THE PHOTON MODEL OUTRANKS THE CVE. The CVE σ is dominated by its
# R*⟨Δx²⟩ blur term (~90% of σ² at R ≈ 1/6), so it mostly reports sqrt(R) x rms
# step - sensitive to the assumed R, and inflated by population heterogeneity and
# drift. L/(4*sqrt(N)) is a per-localization instrument model that assumes nothing
# about the motion. The CVE is kept as an independent cross-check on L.
trk_sigma_loc_rv     <- reactiveVal(NA_real_)   # nm, from the last unified MSD fit
trk_sigma_loc_cve    <- reactiveVal(NA_real_)   # nm, blur-aware CVE (cross-check)
trk_sigma_loc_photon <- reactiveVal(NA_real_)   # nm, median L/(4*sqrt(eco))
trk_sigma_loc_eff <- reactive({
  manual <- suppressWarnings(as.numeric(input$trk_sigma_loc_nm %||% 0))
  if (is.finite(manual) && manual > 0) return(manual)
  ph <- trk_sigma_loc_photon()
  if (is.finite(ph) && ph > 0) return(ph)
  cve <- trk_sigma_loc_cve()
  if (is.finite(cve) && cve > 0) return(cve)
  s <- trk_sigma_loc_rv()
  if (is.finite(s) && s > 0) s else NA_real_
})
# Which source is currently active (for status readouts).
trk_sigma_loc_src <- reactive({
  manual <- suppressWarnings(as.numeric(input$trk_sigma_loc_nm %||% 0))
  if (is.finite(manual) && manual > 0) "manual override"
  else if (is.finite(trk_sigma_loc_photon())) "photon budget"
  else if (is.finite(trk_sigma_loc_cve())) "CVE (cross-check)"
  else if (is.finite(trk_sigma_loc_rv()))  "MSD fit"
  else "not set"
})

# GOAL: the record of the sigma_loc one panel used, taken inside its Run handler
# so the panel can later say what it subtracted (sigma_used_text(), R/trk_stats.R;
# rendered by server/tracking/26_sigma_used.R). Attached to the panel's result as
# attribute "sigma_info", which leaves the result itself unchanged.
#   sigma  : the value the panel used (nm); NA when none
#   source : where it came from; NULL = the app-wide source now in force
#   role   : "subtracted" or "noise floor"
#   used   : FALSE when the panel's own switch turned the correction off
.sigma_stamp <- function(sigma, source = NULL, role = "subtracted", used = TRUE, note = NULL) {
  isolate({
    pat <- input$trk_tcp_pattern %||% TCP_PATTERN_DEFAULT
    typed <- suppressWarnings(as.numeric(input$trk_crlb_pref))
    list(sigma = if (length(sigma)) as.numeric(sigma[1]) else NA_real_,
         source = source %||% trk_sigma_loc_src(), role = role, used = used, note = note,
         a = .a_in_use(pat), a_calibrated = length(typed) == 1L && is.finite(typed) && typed > 0,
         profile = input$trk_cal_profile %||% "")
  })
}


# Per-track cycle-CONTINUITY score. A step spanning k acquisition cycles (k=1 =
# integrated, k≥2 = k−1 dark cycles) is weighted w(k) = 2^(−(k−1)/H), H = the
# dark-cycle "half-life" (tolerance): H small = harsh (only integrated counts),
# H large = lenient (a few dark cycles barely penalised, so 1–2-cycle steps are
# forgiven while 10–20-cycle steps are heavily down-weighted). The per-track score
# is the mean weight over its steps, in [0,1] (1 = every step integrated). Returns
# a named vector (tid → score). Half-life 0/NA → binary integrated fraction.
.continuity_scores <- continuity_scores                  # R/est_clean.R
