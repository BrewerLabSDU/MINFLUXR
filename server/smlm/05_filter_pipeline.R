## ─────────────────────────────────────────────────────────────────────────────
## server/smlm/05_filter_pipeline.R  -  THE FILTER PIPELINE of the active file
## ─────────────────────────────────────────────────────────────────────────────
## Part of smlm_server(). Not a standalone script: its statements are evaluated in
## the server function's environment, in file-name order (server/parts.R).
##
## Stage 1 metadata ranges, stage 1b outlier localizations, stage 2 trace
## spread, stage 3 trace length. Output: filtered_df, which plottable_data
## turns into the analysis data.
##
## INTERFACE (generated: Rscript tools/server_map.R --headers; do not edit by hand)
##   DEFINES for other parts:
##     filtered_df, metadata_filtered_df, sd_filt, sd_skip_arming,
##     tid_reset_flag
##   READS from other parts:
##     03_dataset_zcal: uploaded_df
##     27_fiducial_qc: qc_exclude_from
## END INTERFACE
## ─────────────────────────────────────────────────────────────────────────────

# --------- 5) RANGE UI ( ----------
# Store committed default ranges for metadata + spread filters (used to avoid "filtering on load").
# Updated in two steps per upload: (1) frozen to current inputs at the top of the upload observer
# so the all.equal bypass passes during the browser round-trip window, then (2) updated to the
# new data's ranges inside onFlushed after apply_filters runs.
# Bypass check for metadata filters: NULL = no bypass; finite c(lo,hi) = committed
# full-range value that the all.equal check in metadata_filtered_df compares against.
meta_full_rng <- reactiveValues(
  efo = NULL, cfr = NULL, tim = NULL, dcr = NULL, eco = NULL, dwell = NULL, fluo = NULL
)
# Single source of truth for the SD (spread) filter.
# c(0, Inf) = "no filter" - passes all traces regardless of spread.
# Set to Inf on upload/reset; set to a finite upper bound only when the user
# explicitly narrows the slider. Using Inf avoids all floating-point edge cases
# that arise when comparing widget values against data maxima.
sd_filt <- reactiveValues(x = c(0, Inf), y = c(0, Inf), z = c(0, Inf))
# When TRUE, suppress the SD range observer (widget initialisation fires that
# should not count as user-driven narrowing).
sd_skip_arming <- reactiveVal(FALSE)
# When TRUE, tid_length_range_ui resets to [min, max] on next render instead
# of preserving the current value.  Set on new upload; cleared after one render.
tid_reset_flag <- reactiveVal(FALSE)


# ==========================================================================
# FILTER PIPELINE - Stage 1 of 3: per-LOCALIZATION quality
# ==========================================================================
# GOAL: drop individual localizations that fail the microscope's own quality
# metrics, before any trace-level judgement is made.
#
# These are the per-point columns MINFLUX reports (see the canonical schema in
# helpers_shared.R):
#     efo  photon rate      -> low efo = few photons = poor precision
#     cfr  centre-freq ratio-> molecule poorly centred in the doughnut
#     eco  photon count     dcr  channel ratio (the 2-colour axis)
#     tim  acquisition time -> lets you cut the drifty tail of a long run
#
# KEY DESIGN CHOICE - `rng_filter` only applies a filter when the user's range
# actually DIFFERS from the dataset's full range. A slider sitting at its
# default must be a true no-op: comparing floating-point widget values against
# data maxima otherwise drops boundary points for no reason.
#
# INPUT : uploaded_df()   OUTPUT: same columns, fewer rows
# ==========================================================================
metadata_filtered_df <- reactive({
  # Drop Fiducial-scan-flagged traces first (only when the exclusion is on;
  # reversible - see qc_exclude_from). Active file = input$active_file.
  df <- qc_exclude_from(uploaded_df(), input$active_file)

  # Apply range filter ONLY if user range differs from the dataset's full range
  rng_filter <- function(dat, nm, rng_val, full_rng) {
    if (!nm %in% names(dat)) return(dat)
    if (is.null(rng_val) || length(rng_val) != 2 || any(!is.finite(as.numeric(rng_val)))) return(dat)
    
    v <- suppressWarnings(as.numeric(dat[[nm]]))
    if (!length(v)) return(dat)
    
    # Prefer stored full range; fall back to computed full finite range if needed
    if (is.null(full_rng) || length(full_rng) != 2 || any(!is.finite(as.numeric(full_rng)))) {
      vv <- v[is.finite(v)]
      if (!length(vv)) return(dat)
      full_rng <- range(vv, na.rm = TRUE)
    }
    
    # If unchanged from defaults, treat as "no filter"
    if (isTRUE(all.equal(as.numeric(rng_val), as.numeric(full_rng), tolerance = 0))) return(dat)
    
    # When filtering is active: keep NA rows (so they don't vanish unexpectedly)
    dplyr::filter(dat, is.na(.data[[nm]]) | dplyr::between(.data[[nm]], rng_val[1], rng_val[2]))
  }
  
  df %>%
    rng_filter("efo",   input$efo_range,   meta_full_rng$efo) %>%
    rng_filter("cfr",   input$cfr_range,   meta_full_rng$cfr) %>%
    rng_filter("tim",   input$tim_range,   meta_full_rng$tim) %>%
    rng_filter("dcr",   input$dcr_range,   meta_full_rng$dcr) %>%
    rng_filter("eco",   input$eco_range,   meta_full_rng$eco) %>%
    rng_filter("dwell", input$dwell_range, meta_full_rng$dwell) %>%
    rng_filter("fluo",  input$fluo_range,  meta_full_rng$fluo)
})

# ===========================================================================
#  SMLM FILTER PIPELINE  (read top-to-bottom; each stage feeds the next)
#    Stage 1  metadata  metadata_filtered_df   (defined just above)
#    Stage 2  spread    spread_filtered_df
#    Stage 3  length    filtered_df            -> feeds plotting / clustering
#  The Trace-ID-Length slider UI (tid_length_range_ui) is rendered right after
#  Stage 3 because it derives its range from the spread-filtered data.
# ===========================================================================

# ==========================================================================
# FILTER PIPELINE - Stage 2 of 3: per-TRACE spread (the SD filter)
# ==========================================================================
# GOAL: remove whole traces whose localizations are too SCATTERED to be one
# stationary molecule.
#
# THE REASONING. In SMLM a `tid` should be one fluorophore sitting still while
# blinking, so its localizations should scatter only by the localization
# precision - a few nm. A trace with a much larger spread is not a molecule:
# it is usually stage drift during a long blink, or two nearby molecules whose
# localizations were merged into one trace ID. Either way its median position
# is meaningless, and the median is exactly what the clustering and every
# spatial statistic downstream will use.
#
# So the filter is per-axis SD of each trace's localizations, computed AFTER
# Stage 1 (a trace should be judged on its surviving good points, not on
# points already rejected as low quality).
#
# `sd_filt` holds the committed bounds. An upper bound of Inf means "no
# filter" - an explicit sentinel, chosen so no floating-point comparison
# against a data maximum is needed to express "keep everything".
#
# INPUT : metadata_filtered_df()   OUTPUT: whole traces removed
# ==========================================================================
# ==========================================================================
# FILTER PIPELINE - Stage 1b of 3: outlier LOCALIZATIONS inside a trace
# ==========================================================================
# GOAL: drop the individual localizations that do not belong to their trace,
# BEFORE any per-trace statistic is computed from them.
#
# WHY IT SITS HERE. Every later stage measures a trace: its spread (Stage 2),
# its length (Stage 3), its median position (clustering). All of those are
# computed from the surviving points, so a stray localization has to be gone
# before them or it corrupts each one in turn. Placing it after the metadata
# filters is deliberate for the same reason those come first - a point already
# rejected on quality should not get a vote on whether its neighbours are
# outliers.
#
# This is the ONLY place localizations are removed for spatial reasons.
# Trace-level judgement stays in Stage 2 and the rogue rules; the two are kept
# apart so that trimming can never quietly rescue a trace that is not one
# molecule - the trim FRACTION is reported precisely so that case stays visible.
#
# OFF by default: opening the spread panel shows the untrimmed data first, and
# resetting the filters returns to it.
#
# INPUT : metadata_filtered_df()   OUTPUT: outlier localizations removed
# ==========================================================================
trim_flags <- reactive({
  df <- metadata_filtered_df()
  if (!isTRUE(input$trim_outliers) || !is.data.frame(df) || !nrow(df))
    return(NULL)
  tryCatch(
    trace_point_outliers(df, k = input$trim_k %||% SMLM_DEFAULTS$trim_k,
                         keep_n = input$sd_min_locs %||% SMLM_DEFAULTS$sd_min_locs),
    error = function(e) NULL)
})

trimmed_df <- reactive({
  df <- metadata_filtered_df()
  fl <- trim_flags()
  if (is.null(fl) || !is.data.frame(df) || length(fl) != nrow(df)) return(df)
  df[!fl, , drop = FALSE]
})

# Measured, not theoretical - the Gaussian rejection rate for a given k is
# wrong by about a factor of ten on real MINFLUX residuals.
output$trim_summary_ui <- renderUI({
  fl <- trim_flags()
  df <- metadata_filtered_df()
  if (is.null(fl) || !is.data.frame(df) || length(fl) != nrow(df))
    return(helpText("No data to trim yet."))
  s <- trim_summary(df, fl)

  # Trimming shortens traces, so a trace can drop below the Trace-ID-Length
  # minimum and be removed by Stage 3 even though it passed before. That is
  # the correct reading - a trace with 4 valid localizations IS a 4-point
  # trace - but it must not happen silently, so it is counted here rather
  # than prevented. Raising the trim guard to the length minimum instead
  # would make WHICH POINTS ARE OUTLIERS depend on how long traces are
  # required to be, and would keep known-bad localizations purely to satisfy
  # a count - corrupting the very median and spread the length filter exists
  # to protect.
  lo <- suppressWarnings(as.numeric((input$tid_length_range %||% c(0, Inf))[1]))
  demoted <- NA_integer_
  if (is.finite(lo) && lo > 0) {
    nb <- table(as.character(df$tid))
    na <- table(as.character(df$tid[!fl]))
    na <- na[match(names(nb), names(na))]; na[is.na(na)] <- 0
    demoted <- sum(nb >= lo & na < lo)
  }

  helpText(HTML(paste0(sprintf(
    "Rejected <b>%s</b> of %s localizations (<b>%.2f%%</b>) in %s of %s traces (%.0f%%).<br/>%s",
    format(s$n_flag, big.mark = ","), format(s$n_pts, big.mark = ","), s$pct,
    format(s$n_tr_touched, big.mark = ","), format(s$n_tr, big.mark = ","), s$pct_tr,
    if (s$n_tr_heavy > 0)
      sprintf(paste0("<b>%s trace(s) lost ≥ 30%% of their localizations</b> - that is ",
                     "not a clean trace with a few strays, and is worth inspecting before ",
                     "trusting its position."), format(s$n_tr_heavy, big.mark = ","))
    else "No trace lost 30% or more of its localizations."),
    if (!is.na(demoted) && demoted > 0)
      sprintf(paste0("<br/><b>%s trace(s) fell below the trace-length minimum of %g</b> ",
                     "because of trimming and are now removed by that filter."),
              format(demoted, big.mark = ","), lo)
    else ""
  )))
})

# Stage 2: spread filter - per-trace SD applied after metadata row filters
spread_filtered_df <- reactive({
  df <- trimmed_df()
  if (!is.data.frame(df) || !nrow(df) || !"tid" %in% names(df)) return(df)

  # sd_filt holds the committed filter bounds set by the user.
  # c(0, Inf) upper bound = "no filter" - all traces pass regardless of spread.
  # This avoids every floating-point comparison edge case that arose from
  # matching widget values against data maxima computed from different subsets.
  rx <- c(sd_filt$x[1], sd_filt$x[2])
  ry <- c(sd_filt$y[1], sd_filt$y[2])
  rz <- c(sd_filt$z[1], sd_filt$z[2])

  # "No filter on this axis" needs BOTH ends open, not just the upper one. An
  # upper bound of Inf with a lower bound still set means "at least this much,
  # no ceiling" - treating Inf alone as no-filter would silently discard the
  # lower bound the user typed.
  open <- function(r) is.infinite(r[2]) && r[1] <= 0
  if (open(rx) && open(ry) && open(rz)) return(df)

  # per-tid SDs computed from CURRENT df (after metadata filters)
  sda <- sd_by_tid(df, min_n = input$sd_min_locs)
  if (!is.data.frame(sda) || !nrow(sda)) return(df[0, , drop = FALSE])

  if (!"sd_z" %in% names(sda)) sda$sd_z <- NA_real_

  # A trace whose spread is undefined on an axis (too few localizations) is
  # never judged on that axis - unchanged behaviour.
  in_rng <- function(v, r) is.na(v) | (v >= r[1] & (is.infinite(r[2]) | v <= r[2]))
  pass_x <- open(rx) | in_rng(sda$sd_x, rx)
  pass_y <- open(ry) | in_rng(sda$sd_y, ry)
  pass_z <- open(rz) | in_rng(sda$sd_z, rz)

  keep_tids <- sda$tid[ pass_x & pass_y & pass_z ]

  df %>% dplyr::filter(tid %in% keep_tids)
})
# ==========================================================================
# FILTER PIPELINE - Stage 3 of 3: per-TRACE localization count
# ==========================================================================
# GOAL: keep only traces with a sensible number of localizations. This is the
# last stage; its output is what the analyses consume.
#
# BOTH ENDS OF THE RANGE ARE MEANINGFUL:
#   too FEW  localizations -> the trace's median position is poorly determined,
#                             and it may be a spurious detection rather than a
#                             real fluorophore.
#   too MANY localizations -> often not a single molecule at all: a persistent
#                             bright contaminant, or several molecules sharing
#                             one trace ID. Such a trace would also dominate
#                             any per-localization statistic by sheer count.
#
# Counting happens AFTER Stages 1-2 for the same reason as before: the count
# that matters is of localizations that actually passed quality control.
#
# INPUT : spread_filtered_df()   OUTPUT: the analysis-ready trace set
# ==========================================================================
# Stage 3: trace ID length filter - applied after spread filter
filtered_df <- reactive({
  df <- spread_filtered_df()

  # tid length filter
  if ("tid" %in% names(df) && !is.null(input$tid_length_range)) {
    tl <- table(df$tid)
    keep <- names(tl)[tl >= input$tid_length_range[1] & tl <= input$tid_length_range[2]]
    df <- df %>% dplyr::filter(tid %in% keep)
  }
  # Make positive
  if (isTRUE(input$make_positive)) {
    min_x <- suppressWarnings(min(df$x, na.rm = TRUE))
    min_y <- suppressWarnings(min(df$y, na.rm = TRUE))
    if (is.finite(min_x) && min_x < 0) df$x <- df$x - min_x
    if (is.finite(min_y) && min_y < 0) df$y <- df$y - min_y
    if ("z" %in% names(df)) {
      min_z <- suppressWarnings(min(df$z, na.rm = TRUE))
      if (is.finite(min_z) && min_z < 0) df$z <- df$z - min_z
    }
  }

  df
})

# Trace-ID-Length slider UI - range derived from the spread-filtered data
# (Stage 2 output), so it sits at the end of the pipeline. Stage 3 (filtered_df)
# reads input$tid_length_range from here.
output$tid_length_range_ui <- renderUI({
  df <- spread_filtered_df()
  req("tid" %in% names(df))

  tl <- table(df$tid)
  shiny::validate(shiny::need(length(tl) > 0, "No trace lengths available (after metadata + spread filters)."))

  min_len <- min(tl)
  max_len <- max(tl)

  # Clamp the user's current value against the FULL uploaded data's max trace length,
  # not the spread-filtered max.  Without this, applying a spread filter (which removes
  # long traces) would silently lower the upper bound here.  On reset the clamped value
  # would be re-used, permanently excluding any trace whose length exceeded the old
  # spread-filtered maximum - even after the spread filter is cleared.
  df_full     <- uploaded_df()
  tl_full     <- if (!is.null(df_full) && "tid" %in% names(df_full)) table(df_full$tid) else tl
  max_len_full <- if (length(tl_full) > 0) max(tl_full) else max_len

  # On a fresh upload tid_reset_flag is TRUE → always reset to the new data's
  # full range.  On any other re-render (spread filter change, settings load) →
  # preserve the user's current value so they don't lose their TID filter.
  # The flag is read with isolate() so it does NOT create a reactive dependency
  # here; it is cleared immediately after use so the next render (e.g. triggered
  # by a spread filter tweak) falls through to the preserve-curr branch.
  if (isTRUE(isolate(tid_reset_flag()))) {
    isolate(tid_reset_flag(FALSE))
    val <- c(min_len, max_len)
  } else {
    curr <- isolate(input$tid_length_range)
    if (!is.null(curr) && length(curr) == 2 && all(is.finite(as.numeric(curr)))) {
      v1 <- max(min_len, min(as.numeric(curr[1]), max_len_full))
      v2 <- max(min_len, min(as.numeric(curr[2]), max_len_full))
      val <- c(min(v1, v2), max(v1, v2))
    } else {
      val <- c(min_len, max_len)
    }
  }

  shinyWidgets::numericRangeInput(
    "tid_length_range", "Trace ID Length Range",
    value = val, min = min_len, max = max_len
  )
})
