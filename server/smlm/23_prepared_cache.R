## ─────────────────────────────────────────────────────────────────────────────
## server/smlm/23_prepared_cache.R  -  prepared_df: THE analysis-ready data of any file (batch chain)
## ─────────────────────────────────────────────────────────────────────────────
## Part of smlm_server(). Not a standalone script: its statements are evaluated in
## the server function's environment, in file-name order (server/parts.R).
##
## The cached filter, thinning and clustering chain every batch panel reads,
## its cache key, the stale-results badges, and resolving a panel's unit
## against the de-blinking switch.
##
## INTERFACE (generated: Rscript tools/server_map.R --headers; do not edit by hand)
##   DEFINES for other parts:
##     .small_n_note, .stale_badge, batch_filters, cluster2_label, cluster2_on,
##     deblink_label, deblink_on, prepared_df, resolve_unit
##   READS from other parts:
##     01_view_upload: smlm_file_store, z_scale_applied
##     05_filter_pipeline: meta_full_rng, sd_filt
##     21_density_normalisation: apply_density_thinning, density_p_map
##     24_batch_nn: nn_batch_params_live
##     27_fiducial_qc: fiducial_results, qc_exclude_from
##     31_ripley_pcf: pcf_params, ripley_params
## END INTERFACE
## ─────────────────────────────────────────────────────────────────────────────

## ── Small-N reliability advisory for the per-structure geometry panels ──────
## Every shape metric is estimated from a cluster's own members and is biased at
## low member counts, in DIFFERENT directions, so a bare member count does not tell
## a reader how much to trust a value. Measured against known ground truth:
##   eccentricity  a TRUE SPHERE reads median ecc 0.81 at n = 10, 0.53 at n = 40,
##                 0.34 at n = 200, 0.18 at n = 2000 - PCA orders axes by sample
##                 variance, so the ordering itself manufactures elongation
##   alpha volume  29.9% of truth at n = 10, 54.6% at n = 20, 67.6% at n = 31,
##                 86.2% at n = 80, 96.0% at n = 300 - unusable below ~20
##   sphere radius 147.4 nm at n = 40 against a truth of 150; 149.9 at n = 1000
##   radius of gyration  Bessel-corrected here, but ~12% low at n = 2 from Jensen
## Advisory only: it never blanks a value, because which metric matters decides
## what counts as "enough" and hiding the numbers would hide the size dependence.
.small_n_note <- function(nvec, what = "structures") {
  n <- nvec[is.finite(nvec)]
  if (!length(n)) return(NULL)
  br <- c(sum(n < 20), sum(n >= 20 & n < 40), sum(n >= 40 & n < 200), sum(n >= 200))
  pct <- round(100 * br / length(n))
  sev <- br[1] / length(n)
  col <- if (sev > 0.5) "#b26a00" else if (sev > 0.2) "#7a5c00" else "#777"
  tags$div(style = paste0("font-size:0.82em; color:", col,
                          "; margin-top:4px; line-height:1.45;"),
    tags$b(sprintf("Members per %s: ", what)),
    sprintf("%d%% under 20, %d%% 20-39, %d%% 40-199, %d%% 200+ (median %.0f).",
            pct[1], pct[2], pct[3], pct[4], stats::median(n)),
    if (sev > 0.2) tags$span(
      tags$br(),
      "Shape measurements get less reliable the fewer members a structure has, and they do not ",
      "all fail in the same direction. Enclosed volume comes out too small, and is not usable ",
      "below about 20 members. Eccentricity comes out too large: even a perfect sphere reads as ",
      "elongated when it is built from few points, because the axes are ordered by the scatter of ",
      "those points. Compare an eccentricity against the sphere floor for its own member count, ",
      "given in the exported table, rather than against zero.")
    else NULL)
}

## ── Per-file prepared-data cache ────────────────────────────────────────────
## GOAL: run the filter + de-blink + median-clustering chain ONCE per file per
## settings combination, instead of once per file per PANEL.
##
## WHY: every batch analysis independently repeated
##   qc_exclude_from -> apply_meta_filters_batch -> apply_density_thinning
##   -> apply_clustering_pipeline
## on every loaded file. Measured on 4 real acquisitions (61,737 localizations),
## the panels that cluster took 6.1-11.3 s each while the three that do NOT cluster
## (trace detection 1.1 s, fiducial 0.8 s, precision cloud 4.0 s) were fast - so
## roughly 5 s of every clustering panel was re-deriving identical clusters.
##
## KEY DESIGN, and the thing to preserve if you touch this: prepared_df() builds
## ONE settings object per call (smlm_settings_from_input, R/smlm_settings.R),
## hashes it into the key, and passes THAT SAME OBJECT down the chain. The key
## therefore describes exactly the values the chain ran with. Under-keying (stale
## clusters served silently) cannot come from the settings: a chain function that
## asks for a name outside SMLM_SETTINGS_SPEC gets an error, not a value the key
## never saw. This replaces a hand-kept list of 45 input names and a runtime
## regex scan that checked it.
##
## What the settings do not cover is derived state, added to the key here: the
## data version, the density-thinning keep-probabilities and the fiducial
## exclusion list.
prep_cache_rv <- reactiveVal(list())   # key -> prepared data.frame

# ── The Display tab's filters, applied to ANY stored file ─────────────────────
# GOAL: every all-files panel filters a file exactly as the main plot filters the
# active one. Two pieces of filter state live outside the settings object and are
# handed to the shared chain here:
#   meta_full_rng  the active dataset's full metadata ranges; a slider still at
#                  that range is untouched and filters no file
#   sd_filt        the committed trace-spread ranges (the spread sliders)
# Use batch_filters() wherever a stored file is filtered; do not call
# apply_meta_filters_batch() directly from the server, or the panel will ignore
# the spread filter and clip other files to the active file's metadata range.
smlm_filter_state <- function()
  list(full_rng  = shiny::reactiveValuesToList(meta_full_rng),
       sd_ranges = list(x = sd_filt$x, y = sd_filt$y, z = sd_filt$z))
batch_filters <- function(df, settings) {
  st <- smlm_filter_state()
  apply_meta_filters_batch(df, settings, full_rng = st$full_rng, sd_ranges = st$sd_ranges)
}

# Everything the prepared data depends on, collapsed to one comparable string.
# `s` is the settings object the chain will run with; build it once and pass it
# to both this and the chain.
.prep_key <- function(fname, needs_cluster, two_color, s = smlm_settings_from_input(input)) {
  # density_p_map() is NULL until the density panel has been run; signif() on NULL
  # errors, which silently made every prepared_df() call fail and every panel
  # return no data.
  pm <- tryCatch(density_p_map(), error = function(e) NULL)
  pm_txt <- if (is.numeric(pm) && length(pm))
    paste(names(pm), signif(unname(pm), 8), collapse = ",") else ""

  fr <- tryCatch(fiducial_results(), error = function(e) NULL)
  ex <- if (!is.null(fr) && all(c("flagged", "tid", "dataset") %in% names(fr)))
    sort(fr$tid[fr$flagged %in% TRUE & fr$dataset == fname]) else NULL
  # DATA VERSION. The file name alone does not identify the data: a multiplex
  # is rebuilt under the same name, the axial calibration rescales z in place,
  # and a channel alignment nudge shifts coordinates in place.
  sdf <- smlm_file_store()[[fname]]
  ver <- paste(if (is.data.frame(sdf) && ".dataset_tag" %in% names(sdf)) sdf$.dataset_tag[1],
               z_scale_applied(),
               paste(utils::capture.output(utils::str(attr(sdf, "mplex_off"))), collapse = ""))
  # The filter state that is not in the settings object (see smlm_filter_state).
  flt <- digest::digest(smlm_filter_state(), algo = "xxhash64")
  paste(fname, needs_cluster, two_color, ver, smlm_settings_key(s), pm_txt,
        paste(ex, collapse = ","), flt, sep = "##")
}

## Prepared data for one file: cached filter -> thin -> cluster chain.
## `needs_cluster = FALSE` skips clustering entirely (the panels that only need
## trace positions), and is cached separately so it cannot be served a clustered
## frame or vice versa.
prepared_df <- function(fname, needs_cluster = TRUE, two_color = NULL) {
  store <- smlm_file_store()
  if (!fname %in% names(store)) return(NULL)
  s <- smlm_settings_from_input(input)                 # one object: key AND chain
  if (is.null(two_color)) two_color <- isTRUE(s$enable_two_color)
  key <- .prep_key(fname, needs_cluster, two_color, s)
  cache <- prep_cache_rv()
  if (!is.null(cache[[key]])) return(cache[[key]])

  df <- qc_exclude_from(store[[fname]], fname)
  if (!is.data.frame(df) || !nrow(df)) return(NULL)
  df <- tryCatch(batch_filters(df, s), error = function(e) NULL)
  if (is.null(df) || !nrow(df)) return(NULL)
  df <- apply_density_thinning(df, fname, s)
  if (is.null(df) || !nrow(df)) return(NULL)
  if (isTRUE(needs_cluster)) {
    # ONE pipeline. `two_color` survives only as part of the cache key, because
    # enabling the 2-colour split repartitions the data and so must invalidate
    # a cached frame - it no longer selects between two implementations.
    df <- tryCatch(apply_clustering_pipeline(df, s), error = function(e) NULL)
    if (is.null(df) || !nrow(df)) return(NULL)
  }
  # Bounded: one entry per (file x settings). Settings churn would otherwise grow
  # this without limit over a long session.
  if (length(cache) > 60L) cache <- list()
  cache[[key]] <- df
  prep_cache_rv(cache)
  df
}

## ── "Results are stale" indicator ───────────────────────────────────────────
## Compute controls deliberately take effect only on Run, which is what stops an
## accidental scroll through r max from recomputing anything. The cost is that a
## changed control looks ignored, so each panel exposes a small badge comparing the
## live compute settings against the ones its results were computed with.
## Display-only controls (loess smoothing, bin counts, log axis, trend lines,
## per-file overlays) are NOT listed here - those stay live on purpose.
.stale_badge <- function(live, used) {
  if (is.null(used)) return(NULL)
  live <- live[!vapply(live, is.null, logical(1))]
  used <- used[!vapply(used, is.null, logical(1))]
  keys <- intersect(names(live), names(used))
  diff <- keys[!vapply(keys, function(k)
    isTRUE(all.equal(live[[k]], used[[k]])), logical(1))]
  if (!length(diff)) return(NULL)
  tags$div(style = paste0("margin:4px 0; padding:4px 9px; border-left:4px solid #e67e22;",
                          "background:#fff8f0; font-size:0.82em; color:#7a5c00;"),
    tags$b("Settings changed since this was computed"), " - press Run to apply: ",
    paste(diff, collapse = ", "), ".")
}

output$ripley_stale <- renderUI(.stale_badge(
  list(unit = input$ripley_unit, `r max` = input$ripley_r_max,
       `r steps` = input$ripley_r_steps, `3-D` = isTRUE(input$ripley_use_3d),
       `all files` = isTRUE(input$ripley_all_files),
       envelope = isTRUE(input$ripley_envelope), simulations = input$ripley_nsim,
       `edge correction` = isTRUE(input$spatial_edge_correct),
       `support window` = isTRUE(input$spatial_use_support)),
  ripley_params()))

output$pcf_stale <- renderUI(.stale_badge(
  list(unit = input$pcf_unit, `r max` = input$pcf_r_max,
       `r steps` = input$pcf_r_steps, `3-D` = isTRUE(input$pcf_use_3d),
       `all files` = isTRUE(input$pcf_all_files),
       envelope = isTRUE(input$pcf_envelope), simulations = input$pcf_nsim,
       `edge correction` = isTRUE(input$spatial_edge_correct),
       `support window` = isTRUE(input$spatial_use_support)),
  pcf_params()))

output$nn_batch_stale <- renderUI(.stale_badge(
  list(unit = input$nn_batch_unit, k = input$nn_batch_k), nn_batch_params_live()))

## ── De-blinking is a GLOBAL setting: resolve every panel's unit against it ──
## The per-panel "Unit" selectors used to be independent of `enable_clustering`,
## which produced two silent failures:
##   1. de-blinking ON but a panel left on trace IDs -> that panel quietly ignored
##      de-blinking, so different panels reported on different units;
##   2. a molecule unit chosen while de-blinking is OFF -> apply_clustering_pipeline
##      never creates a `cluster` column, every file yields no points, and the
##      panel came back empty with no explanation. (Verified on two datasets:
##      Ripley, g(r) and Batch NN all returned nothing.)
## `resolve_unit()` maps the "auto" choice - now the default everywhere - onto the
## molecule unit when de-blinking is on and trace IDs otherwise, so switching
## de-blinking on once applies it to every analysis and every dataset. An explicit
## molecule unit with de-blinking off is downgraded to trace IDs and warned about
## rather than silently returning nothing.
## ONE SWITCH PAIR. There used to be two: the single-colour pair
## (enable_clustering / enable_cluster2) and a 2-colour pair
## (enable_two_color_cluster1 / enable_two_color_cluster2), with prepared_df()
## running whichever pipeline matched the mode. That meant a panel gating on
## `enable_cluster2` refused to run for a 2-colour user who had 2-colour Median
## Clustering enabled, and told them to tick a checkbox that would have had no
## effect. There is one pipeline now, so there is one pair of switches and
## nothing left to resolve between them.
two_color_on <- reactive(isTRUE(input$enable_two_color))
deblink_on   <- reactive(isTRUE(input$enable_clustering))
cluster2_on  <- reactive(isTRUE(input$enable_cluster2))

deblink_label  <- reactive("de-blinking")
cluster2_label <- reactive("Median Clustering")

resolve_unit <- function(sel, mol_value, tid_value = "tid", panel = "This panel") {
  sel <- sel %||% "auto"
  if (identical(sel, "auto")) return(if (deblink_on()) mol_value else tid_value)
  if (identical(sel, mol_value) && !deblink_on()) {
    showNotification(sprintf(
      "%s: de-blinked molecules requested but de-blinking is off - using trace IDs. Enable %s in the Display tab to group blinks into molecules.",
      panel, deblink_label()), type = "warning", duration = 7)
    return(tid_value)
  }
  sel
}
