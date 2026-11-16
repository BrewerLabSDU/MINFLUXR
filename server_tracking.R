## """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
## server_tracking.R   -   Tracking analysis server logic
## All input/output IDs prefixed with trk_ to avoid collisions with SMLM IDs.
## Called from server.R inside the combined server function.
## """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
##
## GOAL OF THIS FILE
##   Everything behind the "MINFLUXR - Tracking" tab: load a MINFLUX
##   single-particle-tracking dataset, clean it, animate it, and run the
##   diffusion analyses (MSD, CVE, JDD, HMM, confinement, D-map, ...).
##   It is one large Shiny server function, `tracking_server()`.
##
## HOW TO READ A SHINY SERVER FILE
##   Reactives are LAZY: the order they appear in the file is NOT the order they
##   run. What matters is which reactive calls which. So to follow the data,
##   follow the CALLS. The dependency chain is drawn below.
##
## ── 1. THE DATA FLOW (the single most important thing in this file) ──────────
##
##   upload (.json / .csv / .npy)
##        |  split into per-TID chunk files on disk + manifest.csv
##        |  (helpers_tracking.R::split_by_tid_to_dir - keeps RAM bounded)
##        v
##   trk_raw_df ................ loads EVERY chunk of the current manifest
##        v
##   ---- THE FILTER PIPELINE (each stage feeds the next; ~line 1171) ----------
##   trk_gap_split ............. Stage 1  cut tracks at dark-cycle gaps
##   trk_end_trimmed ........... Stage 2  peel dim / gap-reached track ends
##   trk_quality_filtered ...... Stage 3  drop dark-cycle-heavy tracks
##   trk_meta_filtered ......... Stage 4  efo / cfr / dcr / eco / dwell ranges
##   trk_len_filtered .......... Stage 5  localizations per track
##   trk_dur_filtered .......... Stage 6  track duration
##   trk_spatial_filtered ...... Stage 7  spatial crop
##        v
##   trk() ..................... THE analysis-ready dataset. Almost every
##                               analysis below starts by calling trk().
##
##   Cleaning stages are HEAVY, so they read their parameters via isolate() and
##   only recompute when the "apply cleaning" trigger fires 
##
## ── 2. THE OTHER FOUNDATION: TIME ────────────────────────────────────────────
##   MINFLUX has no fixed clock. Almost every estimator here needs the base
##   acquisition CYCLE TIME to convert irregular timestamps into integer cycle
##   counts. That is what makes the analyses comparable.
##     trk_cycle_time .......... the resolved cycle time (hardware > manual >
##                               data-driven cluster > regression)
##     trk_cluster_cycle ....... data-driven estimate (robust to dark cycles)
##     n_cyc = round(dt / t_c) . cycles spanned by a step; 1 = "integrated"
##   A step with n_cyc > 1 skipped (n_cyc - 1) "dark" cycles.
##
## ── 2b. THE MATH IN ONE PLACE (every estimator below rests on these) ─────────
##
##   NOTATION   r = |Δr| step length (nm) | Δt = elapsed time (ms)
##              d = dimensions (2 or 3)   | D = diffusion coefficient (nm²/ms)
##              σ = σ_loc localization precision (nm)
##
##   (1) FREE DIFFUSION.  Squared displacement grows LINEARLY with time:
##
##           <r²> = 2·d·D·Δt
##
##       so the displacement itself grows as the SQUARE ROOT of time:
##           r ~ sqrt(2·d·D·Δt)
##
##   (2) WHY WE DIVIDE BY sqrt(Δt)
##       MINFLUX steps span different Δt, so before pooling them we must remove
##       the Δt dependence. From (1) that dependence is sqrt(Δt):
##
##           r / sqrt(Δt)  ->  scale sqrt(2dD)   CONSTANT in Δt   (what we use)
##           r / Δt        ->  a VELOCITY, scales as sqrt(2dD/Δt)
##                             -> DIVERGES as Δt -> 0             (never use)
##
##       Brownian motion is continuous but nowhere differentiable: it has no
##       instantaneous velocity. A "speed" r/Δt therefore measures your sampling
##       rate as much as the molecule. (This is exactly why the Speed panel's
##       path speed cannot separate diffusion from transport.)
##
##       Equivalent forms - the same operation, two representations:
##           (r / sqrt(Δt))² = r² / Δt  ∝ D
##       So: JDD + HMM normalise LENGTHS by sqrt(Δt) (Rayleigh-distributed),
##           D-map + change-point normalise SQUARED lengths by Δt (local D).
##
##       CAVEAT: this invariance assumes FREE diffusion. Under anomalous motion
##       <r²> ~ Δt^α the correct exponent is Δt^(α/2), so sqrt(Δt) over-corrects
##       for α<1. And localization noise does NOT scale with sqrt(Δt) (it is
##       constant), leaving a Δt-dependent noise term - hence the σ_loc
##       corrections, and why "integrated steps only" (one fixed Δt, no
##       normalisation needed) is the cleanest mode where it is available.
##
##   (3) LOCALIZATION NOISE adds a constant offset:
##
##           <r²>_measured = 2·d·D·Δt  +  2·d·σ²
##
##       -> MSD fit: SLOPE gives D, INTERCEPT gives σ.
##       -> An immobile molecule still shows apparent D = σ²/Δt. That is the
##          "noise floor" the D-map and change-point subtract.
##       -> A pure-noise step has RMS length sqrt(2d)·σ.
##
##   (4) CYCLE COUNTS. MINFLUX has no fixed clock, so time is measured in whole
##       acquisition cycles: n_cyc = round(Δt / t_c). n_cyc = 1 is "integrated"
##       (no gap); n_cyc > 1 skipped (n_cyc − 1) dark cycles.
##
##   UNITS: nm and ms internally. D in nm²/ms; divide by 1e3 for µm²/s.
##
## ── 3. MAP OF THE FILE (grep the name; line numbers drift, names do not) ─────
##   Setup / state ............. reactiveValues + reactiveVal block, top
##   sigma_loc (noise floor) ... trk_sigma_loc_eff, trk_sigma_loc_src
##   Calibration (CVE) ......... trk_cve_bundle          ~154
##   Upload & chunking ......... trk_raw_df, trk_raw_original ~870
##   FILTER PIPELINE ........... banner at                ~1171
##   Cycle time & TSA .......... banner at                ~1382
##   Playback / animation ...... seq_data, frames_store   ~2925
##   MSD (main engine) ......... msd_unified_bundle       ~6202
##     per-track D / ergodicity  trk_unified_tdiff        ~6410
##   Change-point + GMM ........ cpd_segments_rv          ~4616
##   Jump-distance (JDD) ....... jdd_bundle               ~7078
##   Confinement index (MSD) ... conf_ci_bundle           ~8060
##     (the Rg/Rmax shape test was retired 2026-08-03 -> retired_tracking_server.R)
##   Directionality ............ dir_bundle               ~7627
##   Speed & transport ......... trk_export_summary       ~7978
##   Diffusion-state HMM ....... hmm_bundle               ~8298
##   Per-state MSD ............. seg_msd_rv               ~8826
##   Spatial D-map ............. dmap_bundle              ~8982
##
## ── 4. CONVENTIONS & GOTCHAS (things that have bitten us) ────────────────────
##   * `tim` is MILLISECONDS in here (seconds on disk, x1000 right after load).
##     Coordinates are nanometres.
##   * "TID" = trace ID = one tracked molecule. Always character.
##   * `*_bundle` / `*_rv` reactives are eventReactive: they only recompute when
##     their Run button is pressed. That is deliberate - they are expensive.
##   * Chunking is PLAYBACK-only. The analyses always run on all TIDs; the
##     chunk selector just limits what gets animated so playback stays smooth.
##   * jsonlite is not attached - always write jsonlite:: (its `validate` would
##     shadow shiny::validate).
##   * Some sprintf() format strings contain literal \u escapes (e.g. °)
##     rather than the character itself. Match them exactly when editing.
##
##   Method provenance and the MINFLUX-appropriateness of every estimator are
##   documented in docs/MINFLUX_analysis_methods.md.
## """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

TRACKING_SERVER_PARTS <- load_server_parts(file.path("server", "tracking"))

tracking_server <- function(input, output, session) {
  run_server_parts(TRACKING_SERVER_PARTS, environment())
} # end tracking_server


