# ── Zoom-aware scale bar (mirrors the SMLM tab) ──────────────────────────────
# The bar lives in paper coordinates with width = length_nm / visible_x_span,
# so it must be re-pushed whenever the user zooms (the visible span changes).
sim_axis_x_span <- reactiveVal(NULL)

# TRUE once each plot has actually rendered WITH event_register(), so the
# relayout observers never call event_data() before the source is registered.
# That call would otherwise print a benign "source not registered" warning,
# deferred via session$onFlushed() so suppressWarnings() can't catch it. They
# stay FALSE until the user clicks Generate (the render short-circuits on a
# validate() until then, so the source is never registered before that).
sim_plot_ready <- reactiveVal(FALSE)
sim_diff_ready <- reactiveVal(FALSE)

# Single source of truth: is the FOV / ROI window active?  TRUE in Acquisition
# mode, or in Normal mode when the user pins an ROI ("Pin axes to a FOV").
# Used everywhere the FOV matters - axis pinning, background spread, and the
# confocal window - so they can never disagree.
sim_use_fov <- function() {
  identical(input$sim_mode %||% "normal", "acquisition") ||
    isTRUE(input$sim_view_fov)
}

# Full-view x-span: the FOV X when the FOV is active, else the data extent.
sim_full_x_span <- function() {
  if (isolate(sim_use_fov()))
    return(max(1, as.numeric(isolate(input$sim_fov_x) %||% 4000)))
  res <- tryCatch(isolate(sim_result()), error = function(e) NULL)
  if (is.null(res)) return(NULL)
  ref  <- if (nrow(res$locs) > 0) res$locs else res$gt
  span <- suppressWarnings(diff(range(ref$x, na.rm = TRUE)))
  if (length(span) == 1 && is.finite(span) && span > 0) span else NULL
}

sim_push_scale_bar <- function() {
  proxy <- plotlyProxy("sim_plot", session)
  if (!isTRUE(input$sim_show_scale_bar) ||
      !identical(input$sim_view_mode, "2d")) {
    plotlyProxyInvoke(proxy, "relayout",
                      list(shapes = list(), annotations = list()))
    return()
  }
  span <- sim_axis_x_span()
  if (length(span) != 1 || !is.finite(span) || span <= 0) {
    span <- sim_full_x_span()
    if (length(span) == 1 && is.finite(span) && span > 0) sim_axis_x_span(span)
  }
  if (length(span) != 1 || !is.finite(span) || span <= 0) return()
  sb_font <- suppressWarnings(as.numeric(input$sim_scale_font_size))
  sb <- make_scale_bar(
    suppressWarnings(as.numeric(input$sim_scale_bar_length)),
    span,
    suppressWarnings(as.numeric(input$sim_scale_bar_thickness)),
    dark = isTRUE(input$sim_dark_mode),
    scale_font_size = if (is.finite(sb_font) && sb_font > 0) sb_font else 18)
  if (!is.null(sb))
    plotlyProxyInvoke(proxy, "relayout",
                      list(shapes = sb$shapes, annotations = sb$annotations))
}

# Same zoom-aware scale bar for the Confocal view (always white on black).
sim_diff_x_span <- reactiveVal(NULL)   # zoomed visible x-span
sim_diff_full   <- reactiveVal(NULL)   # full rendered-image x-span

sim_push_diff_scale_bar <- function() {
  proxy <- plotlyProxy("sim_diff_plot", session)
  if (!isTRUE(input$sim_show_scale_bar)) {
    plotlyProxyInvoke(proxy, "relayout",
                      list(shapes = list(), annotations = list()))
    return()
  }
  span <- sim_diff_x_span()
  if (length(span) != 1 || !is.finite(span) || span <= 0) span <- sim_diff_full()
  if (length(span) != 1 || !is.finite(span) || span <= 0) return()
  sb_font <- suppressWarnings(as.numeric(input$sim_scale_font_size))
  sb <- make_scale_bar(
    suppressWarnings(as.numeric(input$sim_scale_bar_length)),
    span,
    suppressWarnings(as.numeric(input$sim_scale_bar_thickness)),
    dark = TRUE,
    scale_font_size = if (is.finite(sb_font) && sb_font > 0) sb_font else 18)
  if (!is.null(sb))
    plotlyProxyInvoke(proxy, "relayout",
                      list(shapes = sb$shapes, annotations = sb$annotations))
}
