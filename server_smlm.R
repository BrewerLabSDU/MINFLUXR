## ─────────────────────────────────────────────────────────────────────────────
## server_smlm.R  -  SMLM analysis server logic
## Called from server.R inside the combined server function.
## ─────────────────────────────────────────────────────────────────────────────
##
## GOAL OF THIS FILE
##   Everything behind the SMLM (imaging) tab: load localization data, filter it,
##   group it into molecules and structures, and run the spatial analyses
##   (nearest neighbour, Ripley, clustering, shape fitting, two-colour
##   co-localization, ground-truth estimation). One Shiny server function,
##   assembled from the part files in server/smlm/ (section 3).
##
##   The Tracking tab asks how molecules MOVE. This tab asks where they ARE and
##   how they are ARRANGED.
##
## HOW TO READ A SHINY SERVER FILE
##   Reactives are LAZY - file order is NOT execution order. Follow the CALLS.
##   The dependency chain is drawn below.
##
## ── 1. THE DATA FLOW ─────────────────────────────────────────────────────────
##
##   upload (.json / .csv / .npy)
##        |   readers + canonicalisation live in helpers_shared.R
##        v
##   uploaded_df_rv -> uploaded_df() ..... the raw, unfiltered dataset
##        v
##   ---- THE FILTER PIPELINE (three stages, each feeds the next) --------------
##   metadata_filtered_df .... Stage 1  per-LOCALIZATION quality ranges
##                                      (efo / cfr / dcr / eco / tim ...).
##                                      Drops individual points.
##   spread_filtered_df ...... Stage 2  per-TRACE spread (SD) filter.
##                                      Drops whole traces that are too diffuse
##                                      - i.e. a "trace" that is really drift or
##                                      two molecules merged.
##   filtered_df ............. Stage 3  per-TRACE localization count.
##                                      Drops traces that are too short/long.
##        v
##   plottable_data() ........ THE analysis-ready dataset. Adds the trimmer, the
##                             2D crop, and - importantly - runs the CLUSTERING.
##                             Almost every analysis starts by calling this.
##
##   ORDER MATTERS AND IS DELIBERATE: rows are filtered first (Stage 1), then
##   whole traces are judged on what survived (Stages 2-3). Filtering traces
##   before cleaning their points would judge them on data you already rejected.
##
## ── 2. THE GROUPING HIERARCHY (the key concept of this tab) ──────────────────
##   localizations -> tid -> cluster (primary) -> cluster2 (secondary)
##
##     tid ....... one continuous blink of one fluorophore
##     cluster ... DE-BLINKING. Merges repeat appearances of one molecule into
##                 one object, at a threshold of a few times the measured trace
##                 spread (~10-30 nm). Complete-linkage on trace centroids, or
##                 a majority vote of the localizations - NOT DBSCAN. It merges
##                 and never removes, so this column has NO id 0.
##     cluster2 .. MEDIAN CLUSTERING over the medians above, LARGE eps
##                 (~100-300 nm), by DBSCAN / HDBSCAN / OPTICS. Job: STRUCTURE -
##                 group molecules into the biological object they form.
##
##   WHY IT EXISTS: a fluorophore blinks many times, and an antibody carries
##   several dyes 10-40 nm from its target. Statistics run on RAW localizations
##   therefore measure the LABEL, not the biology. In `cluster2`, id 0 means
##   noise/unassigned and must be excluded before counting; `cluster` has no 0
##   because de-blinking assigns every trace (see deblink_assign_singletons).
##   (Implementation: apply_clustering_pipeline() in R/smlm_clustering.R.)
##
## ── 3. WHERE THE CODE IS ─────────────────────────────────────────────────────
##   This file only assembles smlm_server(). Its body is in server/smlm/, one
##   panel or core concern per part file, evaluated in this order in the server
##   function's environment (server/parts.R, server/README.md). Each part's
##   header lists what it defines for the others and what it reads from them.
##
##   01_view_upload.R                   view state, 3D rotate/record, file upload
##   02_multiplex.R                     multiplexed MINFLUX: build, align, save
##   03_dataset_zcal.R                  active dataset, axial calibration, filter reset on new data
##   04_settings_collect_apply.R        collect the settings for saving; defaults and apply functions
##   05_filter_pipeline.R               THE FILTER PIPELINE of the active file
##   06_settings_download.R             download the main and the Median Clustering settings files
##   07_deblink_tuning.R                tuning aids for de-blinking and Median Clustering
##   08_settings_load_reset.R           load and reset settings; filtered CSV export
##   09_trace_spread.R                  trace spread inspection panel
##   10_triggers_stores.R               plot trigger and small result stores
##   11_dcr_split.R                     two-colour split by DCR
##   12_metadata_histogram.R            metadata histogram and Median Cluster NN table
##   13_image_overlay.R                 confocal image overlay
##   14_gt_estimator.R                  ground-truth estimate (exploratory)
##   15_scatter2d.R                     main 2D plot
##   16_crop_range_ui.R                 3D crop controls and the metadata range sliders
##   17_plottable_data.R                plottable_data: THE analysis-ready data of the active file
##   18_nn_tables.R                     nearest-neighbour tables and their downloads
##   19_scatter3d.R                     main 3D plot
##   20_dataset_groups.R                dataset groups and the per-file summary
##   21_density_normalisation.R         labelling-density profile and thinning
##   22_spatial_window.R                observation window of the spatial statistics
##   23_prepared_cache.R                prepared_df: THE analysis-ready data of any file (batch chain)
##   24_batch_nn.R                      panel: batch nearest-neighbour histogram
##   25_counts_detection_time.R         panels: traces per cluster, detection over time
##   26_precision_cloud.R               panel: precision cloud and sample stability
##   27_fiducial_qc.R                   panel: fiducial scan and exclusion
##   28_crowding_size_nn.R              panels: crowding versus detection time, NN versus cluster size
##   29_cluster_viewer.R                panel: cluster viewer and shape fit (exploratory)
##   30_radius_gyration.R               panel: radius of gyration
##   31_ripley_pcf.R                    panels: Ripley L and pair correlation g(r)
##   32_hotspots.R                      panel: hotspots (local density enhancement)
##   33_panel_grid.R                    figure grid
##   34_two_colour.R                    panels: two-colour statistics
##
##   The statistics themselves (Ripley L/L12, PCF, KDE co-localization, RoG,
##   batch NN, clustering) live in R/smlm_*.R; the parts wire them to inputs
##   and outputs. To find a reactive by name: Rscript tools/server_map.R out/
##
## ── 4. CONVENTIONS & GOTCHAS ─────────────────────────────────────────────────
##   * Units: nanometres. `tim` stays in SECONDS here (the Tracking tab converts
##     to ms; this tab does not).
##   * plottable_data is an eventReactive on plot_trigger(): analyses refresh
##     when Plot is pressed, not on every keystroke.
##   * Multiplex: when active, clustering is forced to run PER CHANNEL so
##     molecules from different exchange rounds are never merged.
##   * Retired code was in deadcode_smlm.R, deleted 2026-11 (in the git history). Several reactives were
##     previously defined TWICE at this scope - the later definition silently
##     shadowed the earlier, so the first never ran. Before adding a reactive,
##     check the name is not already used further down the file.
##
##   Method provenance and MINFLUX-appropriateness of each statistic:
##   docs/MINFLUX_analysis_methods.md (Part II §6).
## ─────────────────────────────────────────────────────────────────────────────

SMLM_SERVER_PARTS <- load_server_parts(file.path("server", "smlm"))

smlm_server <- function(input, output, session) {
  run_server_parts(SMLM_SERVER_PARTS, environment())

  list(
    sd_filt           = sd_filt,
    sd_display_ranges = smlm_sd_display_ranges,
    file_store        = smlm_file_store,
    plottable_data    = plottable_data,
    uploaded_df       = uploaded_df,
    meta_full_rng     = meta_full_rng,
    sd_filtered_tbl   = sd_filtered_tbl,
    sd_all_tbl        = sd_all_tbl
  )
} # end smlm_server





