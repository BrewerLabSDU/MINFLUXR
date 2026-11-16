## ── Support-aware observation window for the spatial statistics ────────────
## One switch, because it is one methodological choice: whether Ripley's L, g(r),
## L12, g12 and the co-localization null measure intensity over the SAMPLE SUPPORT
## or over the axis-aligned bounding box. See the SUPPORT-AWARE OBSERVATION WINDOW
## banner in R/smlm_support.R for the measured effect (a CSR pattern in an L-shaped
## cell read as clustered by +68.9 nm of L-r with the box, -13.5 nm with the
## support; co-localization false positives 17-37% -> 0%).
##
## The support settings are the Labelling Density panel's; if that panel has never
## been opened, these defaults apply so the switch works on its own.
## TWO SEPARATE CONCERNS, deliberately not bundled:
##
## `spatial_edge_correct` - always safe, and the bigger win. The bounding box is
##   a perfectly valid window; it is just missing its boundary correction. For a
##   cell that FILLS the scan field the box IS the cell, so this alone is the whole
##   fix (measured: clustered data in a 5 um cell, bbox peak L-r +196 vs the true
##   window's +202).
##
## `spatial_use_support` - only for a cell that genuinely does NOT fill the frame
##   (concave, or partly out of field). Measured on CSR in an L-shaped cell, the box
##   invents clustering (peak L-r +116 where the truth is 0) and the footprint
##   window removes it (+0.7).
##
## ⚠ The window must be a FOOTPRINT, never the density HDR. A density-derived
## support is "where the sample is dense", which on clustered data is the CLUSTERS,
## and using that as the window conditions on the very thing being tested: measured
## 29% of the cell, crushing peak L-r from +202 to +49 and dragging the tail
## negative - the "straight line into the negative" symptom. So this deliberately
## does NOT follow the Labelling Density panel's estimator; it always uses the
## occupancy grid at its closing plateau, which is a footprint by construction.
## ── ONE CLOSING RADIUS FOR THE WHOLE COMPARISON SET ────────────────────────
## The gap-fill radius used to be chosen per file, from that file's own closing
## sweep, by plateau_close_nm() = the FIRST sweep step growing under 5%. Both
## halves of that are wrong for a between-file comparison:
##
##   PER FILE. Two files of the same biology get different windows, so their
##   areas mean different things and K (which is proportional to the area) is not
##   comparable. Measured on one cell with 26 identical structures, varying ONLY
##   the localizations per structure: the auto radius came out 400, 1600 or
##   2400 nm and the window 12.6-36.8 um^2, a 2.9x spread on identical geometry.
##
##   FIRST-DIP. A patchy sample's closing curve is not monotone-then-flat, it is
##   bumpy: the area flattens once the STRUCTURES fill in, long before the gaps
##   BETWEEN them bridge. Whether a file dips under the tolerance early is luck
##   of its density, and the file that does gets a structure-hugging window and
##   lands next to the 50% rejection threshold (measured: largest_frac 54%).
##
## So the radius is now resolved ONCE for the whole set: the smallest swept
## radius at which EVERY file resolves to one dominant piece, still capped by
## frac_of_bbox so it cannot swallow the frame. Rejection then stops being
## something to fall back from - it is impossible by construction - and if no
## radius satisfies every file the comparison is REFUSED rather than run with
## mixed windows, which is not a comparison at all.
SPATIAL_MAX_BBOX_FRAC <- 0.97   # same cap plateau_close_nm() uses

## Cost note: this runs ONE closing sweep per file - the same work the per-file
## path already did - and builds NO masks and runs NO connected-components. An
## earlier version searched a radius grid with a mask + component test per file
## per radius, which ran before any progress bar and looked like a hang.
##
## Why max-of-plateaus is the right shared radius: closing is MONOTONE, so a
## larger radius can never fragment a window that was already connected. Taking
## the largest per-file plateau therefore leaves every file at least as connected
## as it would have been on its own, while giving them all one radius.
## `force` and `h_nm` exist for the Local Density & Hotspots panel, which needs a
## footprint unconditionally (its statistics are all "per unit of THIS cell's own
## area", so a bounding box does not merely lose power, it changes the number) and
## renders at a finer grid pitch than the curve panels. Everything else calls this
## with the defaults and behaves exactly as before.
spatial_resolve_close <- function(df_list, use_3d = FALSE, force = FALSE,
                                  h_nm = NULL) {
  if (!isTRUE(input$spatial_use_support) && !isTRUE(force))
    return(list(ok = TRUE, close_nm = NULL))
  h    <- max(10, as.numeric(h_nm %||% input$dens_grid_h %||% 100))
  smode <- input$dens_support_mode %||% "grid"
  keep <- if (use_3d) c("x", "y", "z") else c("x", "y")

  prep <- function(df) {
    if (is.null(df) || !nrow(df) || !all(keep %in% names(df))) return(NULL)
    sp <- df[, keep, drop = FALSE]
    if (!use_3d) sp$z <- 0
    sp <- sp[stats::complete.cases(sp), , drop = FALSE]
    if (nrow(sp) < 3L) NULL else sp
  }
  items <- Filter(function(it) !is.null(it$sp),
                  lapply(df_list, function(it)
                    list(sp = prep(it$df), tag = as.character(it$tag %||% ""))))
  if (!length(items)) return(list(ok = TRUE, close_nm = NULL))

  # Radius read off the sweep by the BRIDGING JUMP, not by plateau_close_nm().
  # plateau_close_nm() takes the first step growing under its tolerance, which on
  # a patchy sample is the moment the STRUCTURES fill in - long before the gaps
  # BETWEEN them close. It is a knife-edge: over 20 independent 6-file sets from
  # one generator it picked 400 nm in 8 of them and 2400 nm in the other 12, and
  # every 400 nm set left all six files fragmented. The largest relative jump in
  # the sweep IS the gaps bridging, so take the first settled radius after it -
  # stable at 2400 nm in 20/20, with nothing fragmented.
  bridge_close <- function(sw) {
    if (is.null(sw) || nrow(sw) < 3L) return(NA_real_)
    s  <- sw[order(sw$close_nm), , drop = FALSE]
    ok <- is.finite(s$frac_of_bbox) & s$frac_of_bbox <= SPATIAL_MAX_BBOX_FRAC
    s  <- s[ok, , drop = FALSE]
    if (nrow(s) < 3L) return(NA_real_)
    j <- which.max(s$rel_step[-1]) + 1L               # the bridging jump
    k <- which(seq_len(nrow(s)) > j & is.finite(s$rel_step) & s$rel_step < 0.10)
    if (length(k)) s$close_nm[k[1]] else s$close_nm[nrow(s)]
  }
  per <- vapply(items, function(it) tryCatch({
    p <- bridge_close(support_closing_sweep(it$sp, use_3d = use_3d, h_nm = h))
    if (is.finite(p)) p else NA_real_
  }, error = function(e) NA_real_), numeric(1))

  if (!any(is.finite(per)))
    return(list(ok = TRUE, close_nm = max(0, as.numeric(input$dens_close_nm %||% 200))))

  shared <- max(per, na.rm = TRUE)
  list(ok = TRUE, close_nm = shared,
       per_file = per,
       tags = vapply(items, function(it) it$tag, character(1)),
       # the "your files disagree about their own area scale" signal
       spread = if (sum(is.finite(per)) > 1)
                  max(per, na.rm = TRUE) / max(min(per[is.finite(per)]), 1e-9) else 1)
}

## ONE WINDOW PER FILE, FOR EVERY PANEL THAT DIVIDES BY AN AREA.
## The estimator settings all come from "Density & area normalization" at the top
## of the tab, so Ripley's L, g(r) and Local Density & Hotspots measure the same
## cell the same way. Three panels resolving three different areas for one cell
## cannot be reconciled afterwards: K(r) scales with the window area and every
## hotspot statistic is per unit of it, so the numbers would silently disagree.
##
## `force` is for the hotspot panel, which needs a footprint UNCONDITIONALLY.
## Elsewhere the footprint is an opt-in that improves an estimate; there it is
## the definition of the quantity, since a bounding box makes a half-framed cell
## read at half its density and doubles every enhancement.
## WHICH POINTS DRAW THE OUTLINE is `dens_support_from`, shared with the density
## panel so all four windows match. Neither choice is free of a bias:
##
##   "units"  (default) - trace / molecule centroids. What the analyses actually
##            count, so the outline is the region those units occupy. Its bias is
##            that a sparsely labelled cell yields a raggeder, smaller footprint,
##            so the area carries some of the unit count (measured on synthetic
##            pairs: 5.01e6 vs 3.82e6 nm2 for identical geometry at a 4x labelling
##            difference, turning a 4.0x density ratio into 3.05x).
##
##   "locs" - every localization. A denser sample of the same region, so a
##            smoother outline, but trace LENGTH is photophysics (how long one
##            fluorophore stayed on), and it enters the closing sweep. Measured on
##            the two-condition set: the localization outline inflated the area enough to
##            loosen the CSR-calibrated eps, and the largest-focus size went from
##            cleanly separating the conditions (6-12 against 19-51 units) to
##            overlapping (10-24 against 19-52), by sweeping diffuse background
##            into the control's foci.
##
## Default is "units": hotspots are hotspots OF CENTROIDS, and the failure above
## is a false negative on the headline metric, which is worse than the density
## ratio compression.
spatial_window_of <- function(df, use_3d = FALSE, close_nm = NULL, force = FALSE) {
  if (!isTRUE(input$spatial_edge_correct) && !isTRUE(input$spatial_use_support) &&
      !isTRUE(force))
    return(NULL)
  keep <- if (use_3d) c("x", "y", "z") else c("x", "y")
  if (is.null(df) || !nrow(df) || !all(keep %in% names(df))) return(NULL)
  sp <- df[, keep, drop = FALSE]
  if (!use_3d) sp$z <- 0
  sp <- sp[stats::complete.cases(sp), , drop = FALSE]
  if (nrow(sp) < 3L) return(NULL)
  h <- max(10, as.numeric(input$dens_grid_h %||% 100))
  smode <- input$dens_support_mode %||% "grid"

  if (!isTRUE(input$spatial_use_support) && !isTRUE(force))
    return(bbox_support(sp, use_3d = use_3d, h_nm = h))

  # Shared radius when the caller resolved one for the set; the per-file plateau
  # remains only as the single-file / un-resolved fallback.
  cn <- if (!is.null(close_nm) && is.finite(close_nm)) close_nm else tryCatch({
    sw <- support_closing_sweep(sp, use_3d = use_3d, h_nm = h)
    p  <- plateau_close_nm(sw)
    if (is.finite(p)) p else max(0, as.numeric(input$dens_close_nm %||% 200))
  }, error = function(e) max(0, as.numeric(input$dens_close_nm %||% 200)))
  # The estimator mode is the density block's, not a hard-coded "grid". `hull`
  # and `bbox` return no mask, and the randomisation-based callers (the CSR
  # envelopes, the calibrated eps) cannot work without one, so those two fall
  # back to the bounding box here rather than half-working downstream.
  fp <- tryCatch(estimate_support(sp, use_3d = use_3d, mode = smode, h_nm = h,
                                 close_nm = cn, return_mask = TRUE,
                                 bw_nm = max(10, as.numeric(input$dens_kde_bw %||% 150)),
                                 frac  = min(1, max(0.5, as.numeric(input$dens_kde_frac %||% 0.99)))),
                 error = function(e) NULL)
  if (is.null(fp) || is.null(fp$mask)) return(bbox_support(sp, use_3d = use_3d, h_nm = h))

  # GUARD: reject a window that is many pieces of COMPARABLE size (structures),
  # not one that merely has a few detached slivers (a real cell edge). Measured on
  # four real MINFLUX acquisitions the footprint came out in 2/3/1/1 pieces with
  # >=93% of the window in the largest, so a count-based rule rejected good data.
  # The occupancy footprint stops inside the real cell edge at sparse labelling, so
  # a null simulated inside it is denser than the data and random points read as
  # regular (measured: L-r = -6.9 nm at r = 200 on CSR, 33% below the pointwise
  # band). Grow it back to the region the points actually fill. Not applied to the
  # hotspot panel (`force`), whose eps calibration randomises inside this same mask
  # and is therefore already self-consistent.
  if (!isTRUE(force)) {
    fp <- tryCatch(debias_support_window(fp, nrow(sp)), error = function(e) fp)
  }
  cmp <- mask_components(if (use_3d) fp$mask else mask_project_2d(fp$mask), fp$h)
  fp$n_components  <- cmp$n
  fp$largest_frac  <- cmp$largest_frac
  if (is.finite(cmp$largest_frac) && cmp$largest_frac < 0.5) {
    showNotification(sprintf(
      "Footprint rejected: %d pieces, %.0f%% of the area in the largest, so it traces the structures rather than the cell. Using the bounding box instead; check the support map.",
      cmp$n, 100 * cmp$largest_frac), type = "warning", duration = 12)
    bb <- bbox_support(sp, use_3d = use_3d, h_nm = h)
    bb$rejected_footprint <- TRUE
    return(bb)
  }
  fp
}
# Back-compatible name used by the analysis panels. Also records what each file's
# window came out as, so the readout below can show it rather than the user having
# to infer it from a warning that may or may not have fired.
spatial_window_log <- reactiveVal(NULL)
# The radius the current run resolved for the whole set, shown once in the readout.
spatial_common_close <- reactiveVal(NULL)
spatial_support_of <- function(df, use_3d = FALSE, tag = NULL, close_nm = NULL,
                               force = FALSE) {
  w <- spatial_window_of(df, use_3d, close_nm = close_nm, force = force)
  if (!is.null(w) && !is.null(tag)) {
    # Edge correction on a bounding box is only right when the sample fills the box.
    # When the footprint window is off, measure (without using it) how much of the
    # box a footprint would cover, so the readout can prompt for the footprint
    # window. Uses the footprint at the LARGEST swept gap fill (24 grid steps), the
    # most generous outline: clustered data filling its frame reads ~80% there
    # (closing cannot fill indentations at the box edge), so only a box that is
    # mostly empty at the micrometre scale falls below the prompt threshold.
    # Plan view only; ~0.1 s per file.
    fp_frac <- NA_real_
    if (identical(w$mode, "bbox") && !isTRUE(input$spatial_use_support) && !isTRUE(force) &&
        isTRUE(input$spatial_edge_correct) && !is.null(df) && all(c("x", "y") %in% names(df))) {
      fp_frac <- tryCatch({
        sp <- df[stats::complete.cases(df[, c("x", "y")]), c("x", "y"), drop = FALSE]
        sp$z <- 0                                             # plan view
        h  <- max(10, as.numeric(input$dens_grid_h %||% 100))
        sw <- support_closing_sweep(sp, use_3d = FALSE, h_nm = h)
        if (is.null(sw) || !nrow(sw)) NA_real_
        else sw$frac_of_bbox[which.max(sw$close_nm)]          # most generous footprint
      }, error = function(e) NA_real_)
    }
    row <- data.frame(
      dataset      = tag,
      window       = if (identical(w$mode, "bbox"))
                       (if (isTRUE(w$rejected_footprint)) "bounding box (footprint rejected)"
                        else "bounding box") else "footprint",
      area_nm2     = w$support %||% NA_real_,
      frac_of_bbox = w$frac_of_bbox %||% 1,
      gap_fill_nm  = w$close_nm %||% NA_real_,
      pieces       = w$n_components %||% NA_integer_,
      largest_frac = w$largest_frac %||% NA_real_,
      footprint_frac = fp_frac,        # diagnostic only, filled when the box is in use
      stringsAsFactors = FALSE)
    cur <- spatial_window_log()
    cur <- if (is.null(cur)) row else rbind(cur[cur$dataset != tag, , drop = FALSE], row)
    spatial_window_log(cur)
  }
  w
}

## Resolves the shared radius and stores it for the readout. It never blocks and
## it does NOT warn about the areas differing between files - different cells
## genuinely have different areas, and that is what dividing by the area is for.
## The thing that actually breaks a comparison is files being measured by
## DIFFERENT METHODS, which is only knowable once the windows are built; see
## spatial_report_windows(), called after the run.
##
## Shown under its own progress message: this is a couple of seconds of sweeping
## before any curve exists, and in silence it reads as a hang.
spatial_gate_close <- function(df_list, use_3d = FALSE, force = FALSE,
                               h_nm = NULL) {
  cc <- if (isTRUE(input$spatial_use_support) || isTRUE(force))
    withProgress(message = "Measuring sample area...", value = 0.5,
                 spatial_resolve_close(df_list, use_3d, force = force, h_nm = h_nm))
    else spatial_resolve_close(df_list, use_3d, force = force, h_nm = h_nm)
  spatial_common_close(cc)
  cc
}

## Called AFTER a run, with the tags it produced. Warns only when the set was not
## measured the same way - i.e. some files got the footprint and others fell back
## to the bounding box. That mix is the real problem: the bounding box counts
## empty field as sample, K is proportional to the area, so the fallback file's
## curve is inflated relative to the rest and the sets are not on one scale.
## Deliberately posted after the run so it cannot sit on top of the progress bar.
spatial_report_windows <- function(tags) {
  if (!isTRUE(input$spatial_use_support)) return(invisible(NULL))
  lg <- spatial_window_log()
  if (is.null(lg) || !nrow(lg)) return(invisible(NULL))
  lg <- lg[lg$dataset %in% tags, , drop = FALSE]
  if (nrow(lg) < 2L) return(invisible(NULL))
  fell_back <- grepl("bounding box", lg$window, fixed = TRUE)
  if (!any(fell_back) || all(fell_back)) return(invisible(NULL))
  showNotification(HTML(sprintf(
    "<b>These files were not measured the same way.</b> %s fell back to the bounding box while the rest used the footprint, so their curves are not on the same scale. Check the support map, or turn the footprint window off for all files.",
    paste(utils::head(lg$dataset[fell_back], 4), collapse = ", "))),
    type = "warning", duration = NULL)
  invisible(NULL)
}

output$spatial_window_readout <- renderUI({
  if (!isTRUE(input$spatial_edge_correct) && !isTRUE(input$spatial_use_support))
    return(helpText("Off - spatial statistics use the bounding box with no boundary correction."))
  lg <- spatial_window_log()
  if (is.null(lg) || !nrow(lg))
    return(helpText("Run Ripley's L or g(r) to see the window each file resolved to."))
  rows <- lapply(seq_len(nrow(lg)), function(i) {
    warn <- grepl("rejected", lg$window[i], fixed = TRUE)
    tags$tr(
      tags$td(style = "padding:1px 12px 1px 0;", lg$dataset[i]),
      tags$td(style = paste0("padding:1px 12px 1px 0;",
                             if (warn) "color:#b26a00;font-weight:600;" else ""), lg$window[i]),
      tags$td(style = "padding:1px 12px 1px 0;", sprintf("%.3g nm2", lg$area_nm2[i])),
      tags$td(style = "padding:1px 12px 1px 0;", sprintf("%.0f%% of box", 100 * lg$frac_of_bbox[i])),
      tags$td(style = "padding:1px 12px 1px 0; color:#777;",
              if (is.na(lg$pieces[i])) "-" else
                sprintf("%d piece(s), %.0f%% in largest", lg$pieces[i], 100 * lg$largest_frac[i])))
  })
  cc <- spatial_common_close()
  small <- which(is.finite(lg$frac_of_bbox) & lg$frac_of_bbox < 0.40 &
                   !grepl("bounding box", lg$window, fixed = TRUE))
  # Edge-corrected bounding box on a sample that leaves much of the box empty: the
  # correction then amplifies the framing error, so prompt for the footprint window.
  # 0.70 is a prompt threshold (a sample filling 70% of its box inflates K by
  # ~1.4x). Synthetic check: CSR filling the frame 100%, clustered filling the frame
  # 79%, L-shaped cell 75%, so it fires on half-framed cells and diagonal rods.
  empty_box <- if ("footprint_frac" %in% names(lg))
                 which(is.finite(lg$footprint_frac) & lg$footprint_frac < 0.70 &
                       grepl("bounding box", lg$window, fixed = TRUE))
               else integer(0)
  tags$div(style = "font-size:0.85em; margin-top:6px;",
    tags$div(style = "font-weight:600; margin-bottom:3px;", "Window used per file:"),
    tags$table(style = "border-collapse:collapse;", rows),
    if (!is.null(cc) && isTRUE(cc$ok) && !is.null(cc$close_nm) && is.finite(cc$close_nm))
      tags$p(style = "margin:5px 0 0; color:#2c3e50;",
        tags$b(sprintf("Gap fill: %.0f nm, shared by every file in this run.", cc$close_nm)),
        " One radius means every file's area is measured the same way. The areas themselves are",
        " expected to differ - that is what dividing by the area is for."),
    if (length(unique(lg$window)) > 1L)
      tags$p(style = "margin:5px 0 0; color:#b26a00;",
        tags$b("Mixed methods above. "),
        "Some files used the footprint and others the bounding box. The bounding box counts empty",
        " field as sample and K scales with the window area, so those curves are inflated relative",
        " to the rest and the two groups are not on one scale."),
    if (length(empty_box))
      tags$p(style = "margin:5px 0 0; color:#b26a00;",
        tags$b("Turn on the footprint window. "),
        sprintf("The sample covers about %s of the bounding box in %s. ",
                paste(sprintf("%.0f%%", 100 * lg$footprint_frac[empty_box]), collapse = ", "),
                paste(lg$dataset[empty_box], collapse = ", ")),
        "The box counts that empty field as sample, which reads as clustering, and edge",
        " correction on such a box makes the error larger."),
    if (length(small))
      tags$p(style = "margin:5px 0 0; color:#b26a00;",
        tags$b("Check the support map. "),
        sprintf("%s came out under 40%% of the bounding box. ",
                paste(lg$dataset[small], collapse = ", ")),
        "That is expected if the cell genuinely occupies a small part of the field, but a window",
        " that small can instead be tracing the structures - which understates clustering,",
        " because a window drawn round the clusters conditions on the thing being measured."),
    tags$p(style = "margin:5px 0 0; color:#7a5c00;",
      "Each file's window is built from its own localizations, so the area estimation method is",
      " identical across files."))
})
