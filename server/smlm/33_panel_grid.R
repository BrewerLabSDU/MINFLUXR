## ══════════════════════════════════════════════════════════════════════════
## MULTI-PANEL VIEW  (3 x 2 slots, one analysis per slot)
## ══════════════════════════════════════════════════════════════════════════
## GOAL
##   Put any six of the Analysis-tab plots side by side, so a set of files can be
##   read across several statistics at once instead of by scrolling between
##   panels. Each slot has its own picker; the slot draws whatever that picker
##   names.
##
## HOW IT AVOIDS BEING A SECOND IMPLEMENTATION
##   It does not recompute anything and it owns no analysis of its own. Every
##   entry in the registry below is the SAME builder function the full-size panel
##   renders (`build_<name>_fig()`), so a slot shows exactly what its panel shows,
##   including its "Run this first" message when the panel has not been run. The
##   builders were factored out of the existing renderPlotly() calls for this;
##   their bodies were not otherwise touched.
##
##   CONSEQUENCE, and the thing to remember when using it: a slot is NOT
##   independent of its panel. It reads that panel's controls and its stored
##   results, so changing a setting in the panel changes the slot, and a panel
##   that has not been run shows a placeholder rather than an empty plot. There is
##   deliberately no per-slot copy of the settings - two copies of one control is
##   how the same statistic ends up displayed twice with different parameters.

# id -> label shown in the picker, and the builder that draws it. The `section`
# entry only groups the dropdown, so a long list stays navigable.
.PANEL_SOURCES <- list(
  none              = list(label = "(empty)",                          section = NULL,             fig = NULL),
  hotspot_map       = list(label = "Hotspot map (footprint + enhancement)", section = "Local density", fig = function() build_hotspot_map_fig()),
  hotspot_ecdf      = list(label = "Local-density enhancement curve",   section = "Local density",  fig = function() build_hotspot_ecdf_fig()),
  hotspot_metric    = list(label = "Hotspot metric per cell",           section = "Local density",  fig = function() build_hotspot_metric_fig()),
  hotspot_density   = list(label = "Hotspot metric vs density",         section = "Local density",  fig = function() build_hotspot_density_fig()),
  ripley_plot       = list(label = "Ripley's L(r) - r",                 section = "Spatial",        fig = function() build_ripley_plot_fig()),
  pcf_plot          = list(label = "Pair correlation g(r)",             section = "Spatial",        fig = function() build_pcf_plot_fig()),
  nn_batch_hist     = list(label = "Batch NN distance histogram",       section = "Spatial",        fig = function() build_nn_batch_hist_fig()),
  crowd_plot        = list(label = "Crowding",                          section = "Spatial",        fig = function() build_crowd_plot_fig()),
  density_plot      = list(label = "Labelling density profile",         section = "Density",        fig = function() build_density_plot_fig()),
  density_map       = list(label = "Support map",                       section = "Density",        fig = function() build_density_map_fig()),
  density_sweep_plot = list(label = "Support closing sweep",            section = "Density",        fig = function() build_density_sweep_plot_fig()),
  rog_hist          = list(label = "Radius of gyration",                section = "Median cluster", fig = function() build_rog_hist_fig()),
  size_nn_plot      = list(label = "Cluster size vs NN",                section = "Median cluster", fig = function() build_size_nn_plot_fig()),
  tid_count_hist    = list(label = "Units per median cluster",          section = "Median cluster", fig = function() build_tid_count_hist_fig()),
  det_time_plot     = list(label = "Detections over time",              section = "Time",           fig = function() build_det_time_plot_fig()),
  tid_time_plot     = list(label = "Trace IDs over time",               section = "Time",           fig = function() build_tid_time_plot_fig()),
  cross_nn_hist     = list(label = "Cross-channel NN",                  section = "Two colour",     fig = function() build_cross_nn_hist_fig()),
  coloc_curve       = list(label = "Co-localization curve",             section = "Two colour",     fig = function() build_coloc_curve_fig()),
  ripley12_plot     = list(label = "Bivariate Ripley's L12",            section = "Two colour",     fig = function() build_ripley12_plot_fig()),
  pcf12_plot        = list(label = "Bivariate g12(r)",                  section = "Two colour",     fig = function() build_pcf12_plot_fig()),
  kde_coloc_plot    = list(label = "KDE co-localization",               section = "Two colour",     fig = function() build_kde_coloc_plot_fig())
)

# Choices for the slot pickers, grouped by section. "(empty)" stays ungrouped and
# first, so a slot can always be cleared.
panel_source_choices <- function() {
  secs <- unique(Filter(Negate(is.null),
                        lapply(.PANEL_SOURCES, function(s) s$section)))
  out <- list(`(empty)` = "none")
  for (sc in secs) {
    ids <- names(.PANEL_SOURCES)[vapply(.PANEL_SOURCES,
            function(s) identical(s$section, sc), logical(1))]
    grp <- stats::setNames(as.list(ids),
                           vapply(ids, function(i) .PANEL_SOURCES[[i]]$label, character(1)))
    out[[sc]] <- grp
  }
  out
}

output$panel_slot_pickers <- renderUI({
  ch <- panel_source_choices()
  # Defaults spell out the workflow this panel was built for: see the pattern,
  # then its distribution, then the two statistics that read its tail.
  defaults <- c("hotspot_map", "hotspot_ecdf", "pcf_plot",
                "hotspot_metric", "hotspot_density", "ripley_plot")
  tagList(lapply(1:6, function(i)
    column(2, selectInput(paste0("panel_slot_", i), paste("Panel", i),
                          choices = ch, selected = defaults[i], width = "100%"))))
})

## One renderer per slot. Written as an explicit loop with `local()` so each
## iteration captures its own `i`; without it all six outputs would close over
## the final value of the loop variable and every slot would draw panel 6.
for (.slot in 1:6) local({
  i  <- .slot
  id <- paste0("panel_slot_", i)
  output[[paste0(id, "_plot")]] <- renderPlotly({
    sel <- input[[id]] %||% "none"
    src <- .PANEL_SOURCES[[sel]]
    shiny::validate(shiny::need(!is.null(src) && !is.null(src$fig),
                                sprintf("Panel %d: pick an analysis.", i)))
    src$fig()
  })
})
rm(.slot)

## The 3 x 2 grid itself. Built server-side rather than written out in the UI so
## the row height stays ONE control instead of six copies of the same number.
output$panel_grid <- renderUI({
  hgt <- paste0(max(200, min(900, as.numeric(input$panel_row_height %||% 340))), "px")
  row_of <- function(idx) fluidRow(lapply(idx, function(i)
    column(4, plotlyOutput(paste0("panel_slot_", i, "_plot"), height = hgt))))
  tagList(row_of(1:3), tags$div(style = "height:10px;"), row_of(4:6))
})
