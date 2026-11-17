## ─────────────────────────────────────────────────────────────────────────────
## R/trk_settings.R  -  the settings that decide what the Tracking tab computes
## ─────────────────────────────────────────────────────────────────────────────
##
## PURPOSE
##   One definition of every input that changes a Tracking result: track
##   cleaning, the acquisition parameters (cycle time, TCP pattern and diameter,
##   prefactor, beam, z factor), the manual sigma_loc, the cell geometry, the
##   compute settings of every analysis panel, and the seed of the random
##   analyses. For each it records the default, the type, the widget and the
##   pipeline stage. It is the tracking counterpart of R/smlm_settings.R.
##
##   What it is used for:
##     * the settings file (Download / Load filter settings) carries every one of
##       these in a `trk_settings` block, so a reloaded file reproduces the
##       analysis, not only the filters;
##     * the server's fallbacks for an input that has not rendered read
##       TRK_DEFAULTS, so they can no longer disagree with each other;
##     * tests/testthat/test-trk-settings.R checks that every widget in
##       ui_tracking.R starts at the default recorded here.
##
##   The range filters (EFO, length, duration, ...) are not inputs: they live in
##   the server's trk_filt store and are saved beside this block. Display-only
##   inputs (colours, trails, point size, heatmap look, which panels are open)
##   are not here: they change the picture, never a number.
##
## FUNCTIONS
##   trk_settings(...) ............. defaults, overridden by `...`, validated
##   trk_settings_from_input(input)  the same, read from a Shiny `input`
##   trk_settings_from_file(raw) ... the block of a saved settings file, or NULL
##   trk_settings_key(s) ........... hash of the settings, for cache keys
##   trk_settings_restore_plan(s) .. which widget update restores each setting
##   TRK_SETTINGS_SPEC / TRK_DEFAULTS
## ─────────────────────────────────────────────────────────────────────────────

# GOAL: one row per setting.
#   type   flag | num | int | choice | text     (how a value is coerced)
#   widget numeric | checkbox | select | radio | slider | text (how it is restored)
#   na_ok  NA is a valid value (an unset prefactor)
#   stage  where the setting acts; documentation and grouping in the file
.ts <- function(default, type, widget, stage, min = -Inf, max = Inf, choices = NULL, na_ok = FALSE)
  list(default = default, type = type, widget = widget, stage = stage, min = min, max = max,
       choices = choices, na_ok = na_ok)

TRK_SETTINGS_SPEC <- list(
  # ── Track cleaning (applied on Plot Data, before every range filter) ─────
  trk_split_gaps_enable   = .ts(FALSE, "flag", "checkbox", "cleaning"),
  trk_split_adaptive      = .ts(FALSE, "flag", "checkbox", "cleaning"),
  trk_split_gap_cycles    = .ts(7,     "num",  "numeric",  "cleaning", min = 0),
  trk_split_min_len       = .ts(50L,   "int",  "numeric",  "cleaning", min = 1),
  trk_endtrim_enable      = .ts(FALSE, "flag", "checkbox", "cleaning"),
  trk_endtrim_dark        = .ts(TRUE,  "flag", "checkbox", "cleaning"),
  trk_endtrim_dark_cycles = .ts(3L,    "int",  "numeric",  "cleaning", min = 1),
  trk_endtrim_efo         = .ts(TRUE,  "flag", "checkbox", "cleaning"),
  trk_endtrim_efo_frac    = .ts(0.5,   "num",  "numeric",  "cleaning", min = 0, max = 1),
  trk_quality_enable      = .ts(FALSE, "flag", "checkbox", "cleaning"),
  trk_quality_darktol     = .ts(3,     "num",  "slider",   "cleaning", min = 0),
  trk_quality_min_integ   = .ts(0.6,   "num",  "slider",   "cleaning", min = 0, max = 1),
  trk_crop_pct            = .ts(50,    "num",  "numeric",  "cleaning", min = 0, max = 100),

  # ── Acquisition: cycle time, sequence, TCP pattern, photon model ─────────
  trk_msd_cycle_ms         = .ts(0,     "num",  "numeric",  "acquisition", min = 0),
  trk_cycle_use_regression = .ts(FALSE, "flag", "checkbox", "acquisition"),
  trk_hw_enable            = .ts(FALSE, "flag", "checkbox", "acquisition"),
  trk_hw_dwell_us          = .ts(100,   "num",  "numeric",  "acquisition", min = 0),
  trk_hw_eod_us            = .ts(5.3,   "num",  "numeric",  "acquisition", min = 0),
  trk_hw_pattern           = .ts("6",   "text", "select",   "acquisition"),
  trk_hw_vertices          = .ts(6L,    "int",  "numeric",  "acquisition", min = 1),
  trk_L_tcp                = .ts(100,   "num",  "numeric",  "acquisition", min = 0),
  trk_tcp_pattern          = .ts("hexagon", "choice", "select", "acquisition",
                                 choices = c("hexagon", "square", "triangle", "octahedron")),
  trk_crlb_pref            = .ts(NA_real_, "num", "numeric", "acquisition", min = 0, na_ok = TRUE),
  trk_beam_fwhm            = .ts(0,     "num",  "numeric",  "acquisition", min = 0),
  trk_geofactor_base       = .ts(360,   "num",  "numeric",  "acquisition", min = 0),
  trk_dw_floor             = .ts(0,     "num",  "numeric",  "acquisition", min = 0),
  trk_z_scale              = .ts(0.7,   "num",  "numeric",  "acquisition", min = 0),
  trk_cal_profile          = .ts("",    "text", "select",   "acquisition"),

  # ── The app-wide noise floor and the cell geometry ───────────────────────
  trk_sigma_loc_nm = .ts(0,        "num",    "numeric", "sigma", min = 0),
  trk_motion       = .ts("volume", "choice", "radio",   "geometry", choices = c("volume", "surface", "filament")),
  trk_geom_type    = .ts("none",   "choice", "select",  "geometry",
                         choices = c("none", "sphere", "rod", "ellipsoid", "membrane", "filament")),
  trk_disp_R       = .ts(500,  "num", "numeric", "geometry", min = 0),
  trk_disp_r       = .ts(400,  "num", "numeric", "geometry", min = 0),
  trk_disp_L       = .ts(2000, "num", "numeric", "geometry", min = 0),
  trk_disp_a       = .ts(500,  "num", "numeric", "geometry", min = 0),
  trk_disp_b       = .ts(500,  "num", "numeric", "geometry", min = 0),
  trk_disp_c       = .ts(900,  "num", "numeric", "geometry", min = 0),
  trk_disp_mtilt   = .ts(0,    "num", "numeric", "geometry"),

  # ── Data quality & calibration (CVE) ─────────────────────────────────────
  trk_sig_min_integ = .ts(50L, "int", "numeric", "cve", min = 10),
  trk_sig_skip      = .ts(5L,  "int", "numeric", "cve", min = 0),

  # ── MSD ──────────────────────────────────────────────────────────────────
  trk_msd_min_pts      = .ts(50L,     "int",    "numeric",  "msd", min = 5),
  trk_msd_ndim         = .ts("auto",  "choice", "radio",    "msd", choices = c("auto", "2", "3")),
  trk_uni_min_lag      = .ts(5L,      "int",    "numeric",  "msd", min = 1),
  trk_uni_max_lag      = .ts(50L,     "int",    "numeric",  "msd", min = 1),
  trk_uni_iterative    = .ts(TRUE,    "flag",   "checkbox", "msd"),
  trk_uni_linear_tau   = .ts(TRUE,    "flag",   "checkbox", "msd"),
  trk_uni_track_weight = .ts(TRUE,    "flag",   "checkbox", "msd"),
  trk_uni_drop_negd    = .ts(FALSE,   "flag",   "checkbox", "msd"),
  trk_uni_model        = .ts("linear", "choice", "select",  "msd",
                             choices = c("linear", "confined", "directed", "anomalous")),
  trk_uni_bootstrap    = .ts(FALSE,   "flag",   "checkbox", "msd"),
  trk_msd_geom         = .ts(FALSE,   "flag",   "checkbox", "msd"),
  trk_msd_cluster      = .ts(FALSE,   "flag",   "checkbox", "msd"),
  trk_msd_n_clusters   = .ts(3L,      "int",    "numeric",  "msd", min = 1),
  trk_dmax_t_hw        = .ts(0.05,    "num",    "numeric",  "msd", min = 0),

  # ── Time-signature analysis ──────────────────────────────────────────────
  trk_cycle_dark_max = .ts(10L, "int", "numeric", "tsa", min = 4),

  # ── Jump-distance distribution ───────────────────────────────────────────
  trk_jdd_min_pts         = .ts(5L,       "int",    "numeric",  "jdd", min = 2),
  trk_jdd_n_comp          = .ts(2L,       "int",    "numeric",  "jdd", min = 1, max = 4),
  trk_jdd_step_mode       = .ts("consec", "choice", "radio",    "jdd", choices = c("consec", "binmed", "subsamp")),
  trk_jdd_bin_width       = .ts(1,        "num",    "numeric",  "jdd", min = 0),
  trk_jdd_subsamp_n       = .ts(4L,       "int",    "numeric",  "jdd", min = 2),
  trk_jdd_norm_dt         = .ts(TRUE,     "flag",   "checkbox", "jdd"),
  trk_jdd_integrated_only = .ts(FALSE,    "flag",   "checkbox", "jdd"),
  trk_jdd_min_dt          = .ts(0,        "num",    "numeric",  "jdd", min = 0),
  trk_jdd_max_dt          = .ts(0,        "num",    "numeric",  "jdd", min = 0),
  trk_jdd_sigma_loc       = .ts(0,        "num",    "numeric",  "jdd", min = 0),

  # ── Change point ─────────────────────────────────────────────────────────
  trk_cpd_min_seg   = .ts(8L,   "int",  "numeric",  "cpd", min = 3),
  trk_cpd_max_g     = .ts(4L,   "int",  "numeric",  "cpd", min = 1, max = 6),
  trk_cpd_sub_noise = .ts(TRUE, "flag", "checkbox", "cpd"),

  # ── Confinement index and displacement versus diffusion ──────────────────
  trk_conf_ci_nsim = .ts(200L,  "int",    "numeric",  "confinement", min = 10),
  trk_disp_dims    = .ts("auto", "choice", "radio",   "displacement", choices = c("auto", "2", "3")),
  trk_disp_dman    = .ts(0.1,   "num",    "numeric",  "displacement", min = 0),
  trk_disp_dsrc    = .ts("cve", "choice", "select",   "displacement", choices = c("cve", "msd", "manual")),
  trk_disp_geom    = .ts(FALSE, "flag",   "checkbox", "displacement"),
  trk_disp_min_pts = .ts(10L,   "int",    "numeric",  "displacement", min = 2),

  # ── Directionality ───────────────────────────────────────────────────────
  trk_dir_min_pts  = .ts(6L, "int", "numeric", "directionality", min = 4),
  trk_dir_max_dark = .ts(2L, "int", "numeric", "directionality", min = 0),

  # ── Speed and transport ──────────────────────────────────────────────────
  trk_speed_min_pts      = .ts(20L, "int", "numeric", "speed", min = 2),
  trk_speed_smooth       = .ts(5L,  "int", "numeric", "speed", min = 1),
  trk_speed_sigma_nm     = .ts(0,   "num", "numeric", "speed", min = 0),
  trk_speed_drift_win    = .ts(10,  "num", "numeric", "speed", min = 0),
  trk_speed_drift_minpts = .ts(5L,  "int", "numeric", "speed", min = 2),
  trk_speed_alpha_win    = .ts(21L, "int", "numeric", "speed", min = 3),
  trk_speed_alpha_maxlag = .ts(6L,  "int", "numeric", "speed", min = 2),
  trk_speed_fold         = .ts(3,   "num", "numeric", "speed", min = 0),
  trk_speed_min_run      = .ts(3L,  "int", "numeric", "speed", min = 1),
  trk_speed_edge_guard   = .ts(1,   "num", "numeric", "speed", min = 0),
  trk_speed_onset_thr    = .ts(0,   "num", "numeric", "speed", min = 0),

  # ── Diffusion-state HMM and per-state MSD ────────────────────────────────
  trk_hmm_n_states    = .ts(2L,      "int",  "numeric",  "hmm", min = 1),
  trk_hmm_min_pts     = .ts(10L,     "int",  "numeric",  "hmm", min = 2),
  trk_hmm_boot        = .ts(0L,      "int",  "numeric",  "hmm", min = 0),
  trk_hmm_max_steps   = .ts(150000L, "int",  "numeric",  "hmm", min = 0),
  trk_hmm_bic         = .ts(FALSE,   "flag", "checkbox", "hmm"),
  trk_seg_msd_min_pts = .ts(15L,     "int",  "numeric",  "hmm", min = 2),

  # ── Spatial D-map ────────────────────────────────────────────────────────
  trk_dmap_bins            = .ts(10L,  "int",  "numeric",  "dmap", min = 5),
  trk_dmap_min_pts         = .ts(5L,   "int",  "numeric",  "dmap", min = 3),
  trk_dmap_integrated_only = .ts(TRUE, "flag", "checkbox", "dmap"),
  trk_dmap_sub_noise       = .ts(TRUE, "flag", "checkbox", "dmap"),

  # ── Reproducibility ──────────────────────────────────────────────────────
  # Seed of every random draw in the tab's analyses (bootstrap intervals,
  # simulated nulls, subsampling). Each Run starts from this seed, so a result
  # does not depend on what ran before it. NA = draw a fresh seed each run.
  trk_seed = .ts(1L, "int", "numeric", "reproducibility", na_ok = TRUE)
)
rm(.ts)

# GOAL: coerce one value to its setting's type. Returns the default for a missing
# or unusable value. Numbers are kept as entered (no clamping): the app passes on
# what the user typed, as it did before this file existed.
.trk_coerce <- function(id, v) {
  sp <- TRK_SETTINGS_SPEC[[id]]; def <- sp$default
  if (is.null(v) || !length(v)) return(def)
  switch(sp$type,
    flag = isTRUE(as.logical(v)[1]),
    num = , int = {
      x <- suppressWarnings(as.numeric(v)[1])
      if (is.na(x)) return(if (sp$na_ok) (if (sp$type == "int") NA_integer_ else NA_real_) else def)
      if (sp$type == "int") as.integer(round(x)) else x
    },
    choice = { x <- as.character(v)[1]; if (x %in% sp$choices) x else def },
    text = as.character(v)[1])
}

# Defaults pass through the same coercion as every value read later (5L, not 5),
# so a settings hash does not change when a widget reports its untouched default.
for (.id in names(TRK_SETTINGS_SPEC))
  TRK_SETTINGS_SPEC[[.id]]["default"] <- list(.trk_coerce(.id, TRK_SETTINGS_SPEC[[.id]]$default))
rm(.id)
TRK_DEFAULTS <- lapply(TRK_SETTINGS_SPEC, `[[`, "default")

# GOAL: the settings for one tracking run from code: the defaults, with `...`
# replacing some. An unknown name is an error (a typo would otherwise be ignored).
trk_settings <- function(...) {
  over <- list(...)
  bad <- setdiff(names(over), names(TRK_SETTINGS_SPEC))
  if (length(bad)) stop("trk_settings: unknown setting(s): ", paste(bad, collapse = ", "), call. = FALSE)
  out <- TRK_DEFAULTS
  for (id in names(over)) out[id] <- list(.trk_coerce(id, over[[id]]))
  structure(out, class = c("trk_settings", "list"))
}

# GOAL: the same, read from a Shiny `input` (or any list with these names). A NULL
# field (widget not rendered, or a file that predates it) takes the default.
# Reads every field, so call it inside an observer, handler or isolate().
trk_settings_from_input <- function(input) {
  out <- TRK_DEFAULTS
  for (id in names(TRK_SETTINGS_SPEC)) out[id] <- list(.trk_coerce(id, input[[id]]))
  structure(out, class = c("trk_settings", "list"))
}

# GOAL: a stable hash of the settings, for cache keys and "settings changed" checks.
trk_settings_key <- function(s) {
  if (!inherits(s, "trk_settings")) s <- trk_settings_from_input(s)
  digest::digest(unclass(s)[names(TRK_SETTINGS_SPEC)], algo = "xxhash64")
}

# GOAL: the trk_settings block of a parsed settings file, coerced, or NULL when
# the file has none (written before v0.8). Settings missing from the block are
# left out of the result rather than set to their default, so loading an older
# block changes only what it recorded. Unknown names are ignored.
# raw : the list from jsonlite::fromJSON() of the whole file
# RETURNS list(settings = named list, present = character, version)
trk_settings_from_file <- function(raw) {
  blk <- raw$trk_settings
  if (!is.list(blk) || !length(blk)) return(NULL)
  ids <- intersect(names(blk), names(TRK_SETTINGS_SPEC))
  vals <- lapply(ids, function(id) {
    v <- blk[[id]]
    if (is.null(v) || (length(v) == 1L && is.na(v) && !TRK_SETTINGS_SPEC[[id]]$na_ok)) return(NULL)
    .trk_coerce(id, v)
  })
  names(vals) <- ids
  vals <- vals[!vapply(vals, is.null, logical(1))]
  list(settings = vals, present = names(vals), version = raw$settings_version %||% NA)
}

# GOAL: how to put each setting back into its widget, in a safe order: switches
# first (they show or hide the controls that follow), then the values.
# RETURNS data.frame(id, widget, value) with value a list column.
trk_settings_restore_plan <- function(settings) {
  ids <- names(settings)
  w   <- vapply(ids, function(id) TRK_SETTINGS_SPEC[[id]]$widget, "")
  ord <- order(w != "checkbox", match(ids, names(TRK_SETTINGS_SPEC)))
  data.frame(id = ids[ord], widget = unname(w[ord]), value = I(unname(settings[ids[ord]])),
             stringsAsFactors = FALSE)
}

# GOAL: evaluate `expr` with R's random generator started from `seed`, then put
# the session's generator back as it was. Every random analysis runs inside this,
# so its result is a function of its inputs and the seed alone, and running it
# does not change the draws of anything that runs later. seed NA/NULL = do not
# seed (the generator just continues, as before seeds existed).
with_analysis_seed <- function(seed, expr) {
  seed <- suppressWarnings(as.integer(seed)[1])
  if (!length(seed) || is.na(seed)) return(expr)
  had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  old <- if (had) get(".Random.seed", envir = globalenv(), inherits = FALSE) else NULL
  on.exit(if (had) assign(".Random.seed", old, envir = globalenv())
          else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE))
            rm(".Random.seed", envir = globalenv()), add = TRUE)
  set.seed(seed)
  expr
}
