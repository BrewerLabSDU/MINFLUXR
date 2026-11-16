## ─────────────────────────────────────────────────────────────────────────────
## server_smlm.R  -  SMLM analysis server logic
## Called from server.R inside the combined server function.
## ─────────────────────────────────────────────────────────────────────────────
##
## GOAL OF THIS FILE
##   Everything behind the SMLM (imaging) tab: load localization data, filter it,
##   group it into molecules and structures, and run the spatial analyses
##   (nearest neighbour, Ripley, clustering, shape fitting, two-colour
##   co-localization, ground-truth estimation). One large Shiny server function.
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
## ── 3. MAP OF THE ANALYSES (grep the name; line numbers drift) ───────────────
##   Upload & settings ......... uploaded_df, get_ui_settings, raw_df, df_full_raw
##   Filter pipeline ........... metadata_filtered_df / spread_filtered_df /
##                               filtered_df   (+ the sd_* trace-spread helpers)
##   Analysis-ready data ....... plottable_data
##   Nearest neighbour ......... nn_table, nn_table_cluster2, secondary_nn_lines,
##                               nn_batch_filtered (batch across files/channels)
##   Two-colour / multiplex .... mplex_split_active, mplex_cross_pair,
##                               cross_dcr_nn_table
##   Ground-truth estimate ..... gt_estimate, gt_morphology   [EXPLORATORY]
##   Cluster viewer / shapes ... viewer_ids, viewer_current, viewer_fit_core,
##                               viewer_fit                    [EXPLORATORY]
##   Dataset groups ............ file_group_map, group_color_map
##   Image overlay ............. confocal_meta, image_overlay
##
##   The heavy spatial statistics themselves (Ripley L/L12, PCF, KDE
##   co-localization, RoG, batch NN) live in R/smlm_*.R - this file wires
##   them to inputs and outputs.
##
## ── 4. CONVENTIONS & GOTCHAS ─────────────────────────────────────────────────
##   * Units: nanometres. `tim` stays in SECONDS here (the Tracking tab converts
##     to ms; this tab does not).
##   * plottable_data is an eventReactive on plot_trigger(): analyses refresh
##     when Plot is pressed, not on every keystroke.
##   * Multiplex: when active, clustering is forced to run PER CHANNEL so
##     molecules from different exchange rounds are never merged.
##   * Retired code is in deadcode_smlm.R (NOT sourced). Several reactives were
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





