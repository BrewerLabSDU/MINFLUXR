# MINFLUXR architecture review

Status: sections 1 to 5 were written before any change. Sections 6 to 8 describe
the refactor on branch `refactor/modules` and what it left.

All `file:line` references are to the tree at commit `7e8b6f9` (v0.7.0 plus one
commit), before any refactoring. After the split, find the same code by name:
`Rscript tools/server_map.R <out_dir>` writes a CSV listing every reactive,
observer, output and function of each server closure with its current location.

How this review was made. The three server closures are 30,000 lines, so they were
not only read top to bottom. `tools/server_map.R` parses each closure and lists
every top-level statement with its line span, the closure names it reads and the
input ids it reads. The dependency tables below come from that parse. The data
chains, every Run handler of the tracking estimators, the settings code and the
`R/est_*.R` files were read in full. The plot builders (`build_*_fig`,
`output$scatter3D`, `output$trk_track2D`) were read for structure only.

Limit of the dependency parse: it matches names, so a local variable with the same
name as a closure object (`has_z`, `data`, `loop`, `tick`) counts as a read of that
object. Edge counts are therefore upper bounds. Chains quoted in section 2 were
checked by hand.

---

## 1. Data flow

### 1.1 Load (shared by the three tabs)

| Step | Function | File |
|---|---|---|
| extension | `get_ext()` | helpers_shared.R:109 |
| `.npy` binary | `read_npy()` then `read_minflux_npy_as_df()` | helpers_shared.R:286, 736 |
| `.json` | `read_minflux_json_as_df()` | helpers_shared.R:367 |
| iteration collapse, CFR/ECO carry | `.minflux_finalize_df()`, `collapse_minflux_iterations()` | helpers_shared.R:587, 855 |
| `.csv` | `read.csv` / `read_minflux_csv()` | helpers_shared.R:190 |
| column aliases, missing columns | `canonicalize_minflux_df()` | helpers_shared.R:1034 |
| axial calibration of a CSV | `detect_csv_z_calibration()`, `strip_csv_z_calibration()` | helpers_shared.R:943, 1019 |

The canonical table is `x, y, z` in nm, `tid` character, `tim`, `itr, efo, cfr, eco,
dcr` (documented at the top of helpers_shared.R). `tim` stays in seconds in the SMLM
tab and is multiplied by 1000 right after load in the tracking tab.

### 1.2 SMLM tab (`smlm_server`, server_smlm.R)

```
input$file ──> .load_one_file() [194] ──> smlm_file_store (reactiveVal, name -> df) [175]
                                              │
                    input$active_file [364] ──┤ picks one entry
                                              v
                                    uploaded_df_rv [174] ──> uploaded_df() [1241]
                                              │
              ACTIVE-FILE CHAIN               │                BATCH CHAIN (all files)
  metadata_filtered_df [2160]                 │     prepared_df(fname) [9079], cached in
    qc_exclude_from (fiducials) [10672]       │     prep_cache_rv [9048], key .prep_key [9053]
    7 range filters, no-op at full range      │       qc_exclude_from
  trim_flags / trimmed_df [2254, 2264]        │       apply_meta_filters_batch (R/smlm_filters.R:854)
  spread_filtered_df [2316]  (sd_filt x/y/z)  │         7 range filters, always applied
  filtered_df [2373]  (length, make_positive) │         trim outliers, LENGTH, then SPREAD (target_sd_x/y)
  plottable_data [5856]  eventReactive(plot_trigger)   apply_density_thinning [8186]
    .mplex_residual_shift [906]                        apply_clustering_pipeline (R/smlm_clustering.R)
    apply_trimmer (3D crop) [5715]
    2D crop (smlm_2d_crop_rv)
    remove_rogue_traces
    assign_analysis_group
    apply_clustering_pipeline  <── the same function ──┘
    colours, hover text, pcid/scid
```

Consumers of the **active-file chain** (`plottable_data()`): `output$scatter2D`
[5061], `output$scatter3D` [6391], `nn_table` [6142], `secondary_nn_lines` [6256],
`nn_table_cluster2` [4488], `cross_dcr_nn_table` [5642], `gt_estimate` [4865],
`gt_morphology` [5008], Cluster Viewer (`viewer_snapshot` set at [11114]),
`centroid_suggestion`/`deblink_knee_data` [2593, 2820], `channel_spread_table`
[2714], the trace-spread panel (`sd_all_tbl`, `sd_filtered_tbl` [3522, 3528]),
`output$metadata_histogram` [4373], `output$download_filtered_csv` [2934], KDE
co-localization (active file only) [13790].

Consumers of the **batch chain** (`prepared_df()`), each a Run button that writes a
`*_results` reactiveVal read by a `build_*_fig()` and a download handler:

| Panel | Run observer | Computation (R/) | Result store |
|---|---|---|---|
| Density profile / thinning | `run_density` [8074] | `density_profile_files` (smlm_support.R) | `density_profile_rv` |
| Batch NN | `run_batch_nn` [9209] | `batch_nn_for_file` (smlm_nn.R) | `nn_batch_results` |
| Traces per cluster | `run_tid_count` [9498] | inline dplyr | `tid_count_results` |
| Median Cluster detection over time | `run_det_time` [9682] | inline | `det_time_results` |
| Trace detection over time | `run_tid_time` [9846] | inline | `tid_time_results` |
| Precision cloud / drift | `run_prec_cloud` [9992] | `precision_cloud` (smlm_support.R) | `prec_cloud_results` |
| Fiducial scan | `run_fiducial` [10555] | `fiducial_scan_file` | `fiducial_results` |
| Crowding vs time shift | `run_crowd` [10699] | `crowding_for_file`, `crowding_null_test` | `crowd_results`, `crowd_null_stats` |
| Within-cluster NN vs size | `run_size_nn` [10908] | `size_nn_for_file` | `size_nn_results` |
| Radius of gyration | `run_rog` [11450] | `compute_rog*` (smlm_structures.R) | `rog_results` |
| Ripley L | `run_ripley` [11632] | `spatial_curve_for_points` (smlm_spatial.R) | `ripley_results`, `ripley_env_results` |
| Pair correlation g(r), xi | `run_pcf` [12017], `pcf_xi_results` [12135] | same | `pcf_results`, `pcf_env_results` |
| Hotspots | `run_hotspot` [12569] | `hotspot_for_points` (smlm_hotspots.R) | `hotspot_rows`, `hotspot_detail` |
| Cross-channel NN | `run_cross_nn` [13264] | smlm_nn.R | `cross_nn_results` |
| L12, g12 | `run_ripley12` [13479], `run_pcf12` [13649] | `cross_curve_for_points` | `ripley12_results`, `pcf12_results` |

Settings: `SMLM_SETTINGS_SPEC` (R/smlm_settings.R) names the 45 data-affecting
inputs. `get_ui_settings()` [1463] writes the JSON (`settings_version` 2 with an
`smlm_settings` block, plus `crop`, `fiducial`, `analysis` sub-lists).
`observeEvent(input$load_settings)` [3007] restores it through
`smlm_settings_from_file()`, `smlm_settings_to_legacy()`, `apply_filters()` [1947]
and `restore_settings_plan()` [3201]. Cluster2, two-colour, overlay and multiplex
each have a separate file and loader [2457/3225, 2900/3250, 4762/3377, 1065/1101].

Exports: filtered CSV [2934], NN tables [6237, 6358], plane distances [6326],
cluster2 metrics [6342], per-file summary [8012], one CSV per analysis panel.

### 1.3 Tracking tab (`tracking_server`, server_tracking.R)

```
input$trk_file [1375] ──> load_one_trk [1210] ──> trk_store (name -> fst on disk) [1187]
input$trk_active / "__merged__" ──> .trk_activate [1338] ──> .trk_ingest [1258]
     split_by_tid_to_dir (helpers_tracking.R:410): chunk files + manifest.csv
     rv$original_manifest_dir frozen here
                    v
trk_raw_loaded [1520] (reads the ORIGINAL manifest) ──> trk_raw_df [1574] (z factor)
                    v                  clean_trigger [1867] gates the three cleaning stages
trk_gap_split [1877]        split_gaps            R/est_clean.R
trk_end_trimmed [1906]      trim_track_ends
trk_quality_filtered [1929] gate_continuity
trk_meta_filtered [1778]    filter_ranges   <── trk_filt (reactiveValues, the stored filters) [1796]
trk_len_filtered [1940]     filter_tid_length
trk_dur_filtered [1945]     filter_duration
trk_spatial_filtered [1953] crop: keep tracks with >= pct_thr of points in the box
                    v
trk() [3639]  eventReactive(data_ver): prep_tracking + colours. THE analysis dataset.
```

The same stages as one function, `filter_pipeline()` (R/est_clean.R:274), serve the
rechunk on Plot Data [2621] and Compare Cells [11450].

Time base: `trk_cluster_cycle` [2022] (`cycle_cluster`), `trk_regression_cycle`
[2036], `trk_hw_cycle_time` [2068] and the manual input are resolved by
`trk_cycle_time` [2086]. `trk_seq_params` [2305] parses the sequence JSON, `trk_blur`
[2372] turns it into the blur coefficients every estimator uses.

Noise floor: `trk_sigma_loc_eff` [258] picks manual > photon budget > CVE > MSD
intercept, from three reactiveVals written by observers [654, 7160].

| Panel | Trigger | Bundle | Computation |
|---|---|---|---|
| Data quality / CVE | `trk_cve_run` | `trk_cve_bundle` [405] | `cve_steps`, `cve_core`, `.crlb_doughnut`, `msd_headline` (R/) + 170 lines in the closure |
| MSD | `trk_msd_compute` | `msd_unified_bundle` [6428], `msd_unified_fit` [6593], `trk_unified_tdiff` [6662] | `msd_track_sums`, `msd_fit_olsf*`, `.msd_model_fit` (R/est_msd.R); ensemble averaging in the closure |
| MSD clusters | reactive | `msd_cluster_rv` [5375], `msd_cluster_models` [5429] | closure |
| Change point | `trk_cpd_compute` | `cpd_segments_rv` [5702], `cpd_gmm_rv` [5750] | R/est_cpd.R |
| TSA / cycle classes | `trk_tsa_run` | `cycle_class_bundle` [7428], `dt_vectors` [5097] | closure, `compute_cycle_stats` [7512] |
| JDD | `trk_jdd_run` | `jdd_bundle` [7764] | R/est_jdd.R |
| Confinement index | `trk_conf_ci_run` | `conf_ci_bundle` [8017] | closure + `sim_free_ci_null` (helpers_tracking.R:501) |
| Displacement vs diffusion | run | `trk_disp_bundle` [8227], `trk_disp_geom_null` [8379] | closure + `sim_free_geom_rnet` |
| Directionality | run | `dir_bundle` [8731] | closure |
| Speed and transport | `trk_speed_run` | `trk_export_summary` [9137] | closure fns `trk_speed_one`, `trk_drift_velocity_one`, `trk_local_alpha_one` [8934-9053] |
| HMM | `trk_hmm_run` | `hmm_bundle` [9632] | R/est_hmm.R |
| Per-state MSD | `trk_seg_msd_run` | `seg_msd_rv` [10037] | closure |
| D-map | `trk_dmap_run` | `dmap_bundle` [10255] | R/est_dmap.R |
| Calibration | `trk_cal_run` [11071] | `trk_cal_result` | R/est_calib.R |
| Compare Cells | `trk_cmp_run` [11450] | `cmp_results` | R/est_percell.R |

Playback is a separate system: `seq_data` [3743], `frames_store` [3755],
`output$trk_track2D` [4192], `draw_frame` [4817], the play loop [10403-10834] and the
plain environment `loop` [995].

Settings: three JSON files, filters [2909/2961], playback [3024/3056], overlay
[3125/3151]. Calibration profiles and the Compare Cells table are CSV.

### 1.4 Simulator (`sim_server`, sim_server.R)

`sim_result` [1124] (eventReactive, 500 lines) builds the SMLM ground truth and
localizations from the closure functions `gen_*` [498-720], `expand_to_fluorophores`
[817], `gen_localizations` [972] and the acquisition model [111-240].
`sim_trk_result` [3616] (394 lines) builds tracking data. `analyze_real_run` [244]
and `sim_trk_calib` [3433] calibrate from a real file. Exports are CSV in the
app's own format [3181, 3379, 4239, 4246]. The simulator shares no reactive with the
other two tabs: data moves between tabs only through a downloaded CSV.

---

## 2. Reactive dependency structure

### 2.1 Size and shape

| Closure | Lines | Statements | reactive / eventReactive | reactiveVal | observers | outputs | closure functions |
|---|---|---|---|---|---|---|---|
| `smlm_server` | 13,896 | 470 | 56 / 2 | 67 | 94 | 139 | 78 (3,475 lines) |
| `tracking_server` | 11,432 | 397 | 52 / 17 | 39 | 70 | 110 | 65 (1,272 lines) |
| `sim_server` | 4,254 | 84 | 5 / 3 | 8 | 7 | 24 | 35 (1,126 lines) |

Lines inside `render*` and observers: SMLM 4,799 + 2,417, tracking 4,302 + 1,195,
simulator 1,348 + 261. `output$scatter3D` alone is 1,493 lines [6391],
`output$trk_track2D` 603 [4192], `output$scatter2D` 412 [5061], `sim_result` 500.

### 2.2 Shared across panels (number of statements that read each)

SMLM: `smlm_file_store` 37, `uploaded_df_rv` 22, `plot_trigger` 21, `uploaded_df`
21, `filtered_df` 20, `plottable_data` 20, `file_group_map` 16, `series_palette` 13,
`sd_filt` 11, `mplex_offsets_rv` 10, `prepared_df` 9, `qc_exclude_from` 9,
`density_profile_rv` 9, `z_scale_applied` 8.

Tracking: `loop` 37, `trk_cycle_time` 35, `trk` 26, `layout_trigger` 15, `trk_blur`
15, `msd_unified_bundle` 13, `rv` 13, `frames_store` 12, `trk_cluster_cycle` 12,
`trk_sigma_loc_eff` 10, `trk_filt` 9.

Simulator: `sim_result` 10, `sim_trk_result` 6. Nothing else is read by more than 4.

So each tab has a small core that everything reads (the data chain, the time base,
the triggers) and a ring of panels that read the core and little else. That is the
shape the target structure in section 4 uses.

### 2.3 Long or surprising invalidation chains

1. **SMLM, 9 deep**: `uploaded_df -> metadata_filtered_df -> trim_flags ->
   trimmed_df -> spread_filtered_df -> filtered_df -> plottable_data ->
   gt_estimate -> gt_morphology`. The first six are live reactives, so every slider
   move recomputes the filter chain. Only `plottable_data` waits for
   `plot_trigger`.
2. **Tracking, 11 deep before `trk()`**: `trk_cluster_cycle ->
   trk_regression_cycle -> trk_cycle_time -> trk_tot_cycles ->
   trk_adaptive_gap_cycles -> trk_gap_split -> ... -> trk_spatial_filtered`. The
   cycle time is read with `isolate()` inside the cleaning stages [1884, 1914,
   1933], so a changed cycle time does not re-run them until `clean_trigger` bumps.
   The dependency exists in the data but not in the graph.
3. **`plot_trigger` (SMLM) has 12 writers** [3007, 3269, 4041, 4054, 5514, 5746,
   ...] including `.mplex_build` and the settings loader. Any of them re-runs
   `plottable_data` and with it clustering and both scatter plots.
4. **`uploaded_df_rv` has 6 writers** [284, 364, 404, 627, 744, 1267]. Replacing it
   rebuilds every filter `renderUI` and resets the filters [1387, 4528]. The
   multiplex alignment code works around this with `.mplex_residual_shift` [906]
   rather than writing the reactiveVal.
5. **Hidden edges through observers.** The noise floor travels
   `trk_cve_bundle -> observer [654] -> trk_sigma_loc_photon / _cve ->
   trk_sigma_loc_eff -> JDD, D-map, CPD, HMM, confinement, speed`, and
   `msd_unified_bundle -> observer [7160] -> trk_sigma_loc_rv -> same`. Running
   the CVE or the MSD therefore changes what six other panels will compute on
   their next Run, with no reactive edge between them and no stale badge.
   SMLM has the same pattern: `run_density -> density_profile_rv -> density_p_map
   -> .prep_key -> every batch panel`, and `run_fiducial -> fiducial_results ->
   qc_exclude_from -> both data chains`.
6. **Two triggers with one name.** `plot_trigger` is a reactiveVal in both closures
   [smlm 3991, tracking 3375] with different meanings (SMLM: recompute the analysis
   data; tracking: rebuild playback frames). They are separate objects, but the
   name invites mistakes.
7. **`trk()` is `eventReactive(data_ver())`** and `data_ver` is bumped only inside
   `fire_layout()` [3439]. Every analysis panel therefore depends on the playback
   layout trigger.
8. **`frames_store`, `seq_data`, `heatmap_store`** read `trk()` through
   `plot_trigger`; `output$trk_track2D` writes four reactiveVals while rendering
   (`trace_map`, `trk_axis_x_span`, `trk_locked_ranges`, `trk_plot_ready`). Renders
   with side effects are the reason the upload cascade needed
   `loop$rendered_view` guards [10502].

### 2.4 State that is not reactive

- `loop` [tracking 995]: an environment with 18 fields, read by 37 statements. It
  holds playback state and render bookkeeping.
- `.filt_shown` [tracking 1797]: environment remembering what each filter widget
  last displayed, used to tell a user edit from the echo of a render.
- `.mplex_clu_memo` [smlm 928], `iso_mesh_cache` [smlm 4016], `dcr_track`
  [smlm 4129], `.smlm2d_view` [smlm 5056]: environments used as caches or flags.
- `prep_cache_rv` [smlm 9048] is a reactiveVal used as a cache: `prepared_df()`
  writes it during a read.

---

## 3. Problems

### 3.1 Duplicated logic

**P1. Two SMLM filter chains that do not apply the same filters.** The active-file
chain (server_smlm.R:2160-2395, 5856-5886) and the batch chain
(R/smlm_filters.R:854-901) differ in five ways:

| | Active file | Batch (`prepared_df`) |
|---|---|---|
| range filters | skipped when the slider equals the data range [2181] | always applied [860] |
| spread filter | `sd_filt` x, y, z ranges with lower bounds [2324-2348] | `target_sd_x`, `target_sd_y` upper limits, only when `view_SD` is on, no z [886-898] |
| order | spread, then length [2316, 2373] | length, then spread [879, 886] |
| crops (3D trimmer, 2D crop) | applied [5863-5869] | not applied |
| rogue removal, `make_positive` | applied [5872, 2383] | not applied |

The main plot and its NN tables can therefore describe a different point set from
Ripley, g(r), Batch NN and the other batch panels run on the same file with the
same settings. See question Q1 in section 5.

**P2. `run_tests.R` header is stale.** tests/run_tests.R:9-11 says CVE, MSD OLSF,
JDD, HMM, change-point and D-map are still nested in the server closure. They are
in R/est_cve.R, est_msd.R, est_jdd.R, est_hmm.R, est_cpd.R and est_dmap.R, the
server calls them (server_tracking.R:415, 470, 5731, 6455, 6604, 7800, 9709, 10269),
and tests/testthat/test-estimators-tracking.R tests them. No copy of the core
estimators remains in the server.

**P3. What still duplicates R/ in the tracking closure.**
- `frac_integ_of()` server_tracking.R:423-432 repeats `cve_steps()`
  (R/est_cve.R:33-48).
- `cont_of()` server_tracking.R:443-451 repeats `continuity_scores()`
  (R/est_clean.R:197-210) as a dataset mean.
- `msd_unified_bundle` server_tracking.R:6459-6512 repeats `msd_ensemble()`
  (R/est_msd.R:381-393): TAMSD, pair or track weighting, `keep_lag` at 10% of the
  tracks. `msd_headline()` uses the R/ version, the MSD panel uses the closure
  version. The two agree today; nothing keeps them in step.
- The fit-window parsing `max(1L, as.integer(input$trk_uni_min_lag %||% ..))` and
  the choice between `msd_fit_olsf_iterative` and `msd_fit_olsf` appear six times:
  5343, 5433, 6583, 6594-6613, 6665-6683, 6717-6734.
- `msd_fit_result()` server_tracking.R:5873-5941 is a second anomalous/linear MSD
  fitter (nls, no blur term) beside `.msd_model_fit()` (R/est_msd.R:250).
- `dmax_estimate()` 5977-6002 recomputes median dt and `t_coll = eco/efo`, as do
  `trk_cve_bundle` [457-462] and `trk_blur` [2372].

**P4. SMLM and tracking each carry their own copy of the same feature.**
- Image overlay: `confocal_meta`/`image_overlay`/`.load_smlm_img`/auto-align/nudge
  (server_smlm.R:4584-4788) and `trk_confocal_meta`/`trk_image_overlay`/
  `.load_trk_img`/auto-align/nudge (server_tracking.R:2852-2907, 3237-3347).
- Scale bar: `smlm_push_scale_bar` [121], `trk_push_scale_bar` [10942],
  `sim_push_scale_bar` and `sim_push_diff_scale_bar` [sim 42, 72].
- Axial calibration: `commit_z_scale`/`smlm_z_state`/`z_scale_status`
  (server_smlm.R:1267-1383) and `trk_z_applied`/`trk_z_state`/`trk_z_scale_status`
  (server_tracking.R:1129-1176).
- Range-filter UIs: seven `*_range_ui` (server_smlm.R:5592-5636) and
  `make_range_ui` (server_tracking.R:1630-1652).
- Metadata histogram: server_smlm.R:4373 and server_tracking.R:2770.
- `SMLM_Z_SCALE_DEFAULT <- 0.7` (server_smlm.R:187, in the closure) and
  `TRK_Z_SCALE_DEFAULT <- 0.7` (packages.R:169, global).

**P5. Repeated Run-handler boilerplate (SMLM).** Each of the 16 batch observers
repeats: resolve files (active or all), loop with `withProgress`, `prepared_df`,
`apply_dataset_groups`, bind rows, store params for the stale badge. About 40 of
each handler's 60 to 130 lines are this frame.

### 3.2 Computation inside observers and render functions

- `trk_cve_bundle` server_tracking.R:481-576: per-track sigma, bootstrap CI,
  linking QC, suspect-track count.
- `conf_ci_bundle` 8017-8116, `trk_disp_bundle` 8227-8344, `trk_disp_geom_null`
  8379-8458, `dir_bundle` 8731-8827, `trk_export_summary` 9137-9268, `seg_msd_rv`
  10037-10110, `msd_cluster_rv`/`msd_cluster_models` 5375-5518,
  `trk_unified_tdiff` 6662-6702, `cycle_class_bundle` 7428-7453: the whole
  estimator is in the reactive.
- `trk_seq_params` 2305-2361 parses the sequence JSON in a reactive.
- `output$trk_dt_hist_jdd` 5111-5238 and `output$trk_dt_stats_ui` 5240-5296 each
  compute the half-sample mode and cycle statistics while rendering.
- `output$trk_msd_unified_estimates` 7202-7371 and `output$trk_cve_status` 666-829
  compute derived numbers (ratios, verdict thresholds) while building HTML.
- SMLM: `run_tid_count` 9498-9567, `run_det_time` 9682-9735, `run_tid_time`
  9846-9902 compute with inline dplyr; `output$prec_cloud_plot` 10025-10246 and
  `output$prec_drift_plot` 10344-10528 bin and fit inside `renderPlotly`;
  `pcf_xi_results` 12135-12178; `plottable_data` 5895-6106 mixes clustering,
  colours and hover text; `output$scatter3D` 6391-7883 fits shapes, builds
  isosurfaces and computes the curvature and plane-distance tables it then stores
  in `cluster2_table_store` for download.
- Simulator: all of it. `sim_result` 1124-1623 and `sim_trk_result` 3616-4009 are
  the simulators; the 28 pure generator functions [111-1112, 3402] are defined
  inside the closure, so the tests cannot call the simulator and use their own
  (tests/testthat/helper-smlm-sim.R, `sim_ct_state_tracks` in R/est_hmm.R).

### 3.3 Hidden global state

- `set.seed()` on the global RNG inside the simulator: sim_server.R:1131, 1346,
  1659, 2106, 3617. After a simulation, every later random draw in the session
  (bootstrap CIs, envelopes, null simulations) is a deterministic function of the
  simulator seed.
- No analysis that draws random numbers takes a seed: `.boot_ci`
  (server_tracking.R:301), `sim_free_ci_null`, `sim_free_geom_msd`,
  `sim_free_geom_rnet` (helpers_tracking.R:501, 601, 243), the spatial envelopes,
  the crowding permutation null, density thinning (`apply_density_thinning`
  server_smlm.R:8186), `.dataset_tag` (server_smlm.R:311). Exported numbers that
  depend on them are not reproducible from a settings file.
- `.plugin_tabs`, `.plugin_servers` mutated with `<<-` (global.R:39-45).
- `options(repos=...)`, `options(shiny.maxRequestSize=...)` and package
  installation at source time (packages.R:4-18).
- helpers_smlm.R sources R/smlm_*.R into its own environment and locates itself
  through `sys.frame()`.
- The environments listed in 2.4.

### 3.4 Dead code

- Never sourced: `deadcode_smlm.R` (297 lines), `deadcode_tracking.R` (1,202),
  `retired_tracking_server.R` (328), `retired_ui.R` (81).
- Outputs with no UI binding (never run): `download_curvature_table`
  (server_smlm.R:6372), `trk_msd_tid_ui` (server_tracking.R:1065), `t_window_ui`
  (3349), `trk_msd_min_lag_cycle_hint` (6240), `trk_msd_cycle_status` (6284).
- reactiveVals that are written but never consumed: `smlm_plotted_store`
  (server_smlm.R:4036, kept in sync by the observer at 4062 and read nowhere
  else), `slider_guard` (server_tracking.R:195, incremented or reset at 4904,
  5003, 10474, 10522, 10531 and never tested), `reset_dur_filter_flag` (1743,
  set at 1731, 1760, 2746 and cleared at 2612, never tested),
  `trk_len_reset_flag` (1747, same pattern), `current_pt_idx` (208, no use at
  all).
- reactiveVals that are read but never written: `current_tid_idx`
  (server_tracking.R:197, read at 10690 and 10721, always 1) and
  `conf_color_map_stored` (3385, documented as always NULL).
  `cluster2_table_store` (server_smlm.R:4033) is written only from inside the
  `scatter3D` render [6868], so its download is empty until the 3D plot has drawn.
- Reactives with no reader: `msd_min_lag_ms` (server_tracking.R:2104),
  `msd_fit_max_ms` (2112).
- `deblink_label`, `cluster2_label` (server_smlm.R:9177-9178) are reactives that
  return a constant.
- `fire_layout(src)` / `fire_frames(src)`: the argument is unused (3439, 3446).

### 3.5 Bugs found while reading (not fixed)

- **B1. server_tracking.R:3497** `observeEvent(debounce(rebuild_tick(), 200), ...)`
  passes the value of `rebuild_tick()` to `debounce()`, which expects a reactive.
  Each change builds a new debounce object, the handler fires at once, and one
  orphan observer is left behind per change. The intended 200 ms debounce of
  playback-parameter edits does not happen.
- **B2. Fallback defaults that disagree with each other or with the UI.** They
  only matter when the input is NULL (before the widget exists, and in
  `testServer`), but they make a test session compute something the app never does.
  - `trk_uni_min_lag`: UI 5 (ui_tracking.R:926); fallback 5 at server_tracking.R:
    5343, 6583, 10089 and 4 at 5433, 6552, 6594, 6623, 6665, 6717.
  - `trk_msd_min_pts`: UI 50; fallback 20 at 1068, 50 at 6588.
  - `trk_cycle_dark_max`: UI 10; fallback 10 at 5204, 20 at 7446.
  - `trk_uni_track_weight`, `trk_uni_iterative`, `trk_uni_linear_tau`: NULL means
    FALSE in `msd_unified_bundle`/`msd_unified_fit` [6484, 6598, 6602] and TRUE in
    `.msd_settings` [6585-6587].
  - `dens_support_mode`: "grid" at server_smlm.R:8081, 8465, 8691 and "kde" at 8226.
- **B3. `n_cyc` is defined three ways.** `round(dt / t_c)` in `cve_steps`
  (R/est_cve.R:41) and `jdd_collect_jumps` (R/est_jdd.R:85); `pmax(1,
  round(dt / t_c))` in `dmap_compute` (R/est_dmap.R:68) and `continuity_scores`
  (R/est_clean.R:204); `pmax(1, round((dt - sw) / t_c))` in `msd_track_sums`
  (R/est_msd.R:352). They agree while SW is well under half a cycle and no step is
  shorter than half a cycle.

- **B4. The first "Estimate σ_loc" run does not publish σ.** The observer that
  copies the CVE result into the app-wide noise floor (server_tracking.R:654-659)
  is `observeEvent(trk_cve_bundle(), ..., ignoreInit = TRUE)`. At session start
  `trk_cve_bundle()` stops with `req()` (button not pressed), which does not count
  as the initial run, so the first *successful* run is the one `ignoreInit`
  swallows. After one press `trk_sigma_loc_eff()` is still NA and
  `trk_sigma_loc_src()` says "not set"; JDD, D-map, change-point, HMM and the
  confinement null then run with no noise subtraction. A second press publishes
  it. Seen in `testServer` (pinned as `sigma_loc_first_run` in the tracking
  golden) and it follows from how `bindEvent` handles `ignoreInit`; please
  confirm once in the running app.
- **B5. The change-point CSV download fails when any segment has D = 0.**
  server_tracking.R:5857-5864 calls `mclust::predict.Mclust(gmm, log10_D)` on all
  segments inside `ifelse()`, including those with `log10(0) = -Inf`, and stops
  with "NA/NaN/Inf in foreign function call". D = 0 segments exist whenever noise
  subtraction clamps (see S2), so the export fails exactly for data with an
  immobile fraction. Pinned as an error status in the smoke test.

### 3.6 Inconsistencies between the SMLM and tracking tabs

| Aspect | SMLM | Tracking |
|---|---|---|
| `tim` unit | seconds | milliseconds |
| loaded data | in memory, `smlm_file_store` | fst chunks on disk, `trk_store` |
| several files | active file plus "all files" batch panels | active file or merged view; Compare Cells per file |
| source of truth for range filters | the widget value, with a no-op rule against `meta_full_rng` | `trk_filt` store, widgets only display it |
| apply model | filters live; clustering on Plot | cleaning and filters on Plot Data |
| crop | drops points | keeps or drops whole tracks by fraction inside |
| settings file | one versioned JSON with a typed spec, plus three side files | three unversioned JSON files, no spec |
| defaults | `SMLM_DEFAULTS` in one place | literals in UI and `%||%` fallbacks |
| stale-results badge | Ripley, g(r), Batch NN, hotspots | none |
| results of a Run | reactiveVal written by an observer | eventReactive |
| z factor default | closure constant | global constant |
| input id prefix | none | `trk_` (except `efo_rng` ... `fluo_rng`) |

---

## 4. Target structure

### 4.1 What a "module" can be here

Three constraints decide the design.

1. **Input ids are a public interface.** The settings JSON is keyed by input id
   (`SMLM_SETTINGS_SPEC`), 268 `conditionalPanel` conditions and the custom
   JavaScript refer to ids by name, and the plugin API passes the root `input`.
   `moduleServer()` with `NS()` renames every id inside the module. A v0.7.0
   settings file would then need an id translation layer, and every JS condition
   would need rewriting.
2. **The UI cannot be tested automatically on this machine.** shinytest2 / headless
   Edge is not usable here, so a change that renames ids can only be verified by
   hand. `testServer` verifies server logic, not that a widget still reaches it.
3. **Reactives refer to each other by lexical scope**, often to names defined later
   in the file.

So the target is modules in the sense of *one file, one responsibility, a declared
interface, callable in `testServer`*, without namespacing the ids:

```
R/                      pure computation, no shiny, tested against ground truth
  est_*.R               tracking estimators            (exists)
  smlm_*.R              SMLM computations              (exists)
  sim_*.R               simulator generators           (to extract from sim_server.R)
  trk_*.R               tracking panel computations    (to extract)
server/
  smlm/NN_<panel>.R     one panel or one core concern each, <= ~1,500 lines
  tracking/NN_<panel>.R
  sim/NN_<panel>.R
server_smlm.R           assembles the SMLM parts into smlm_server()
server_tracking.R       same for tracking
sim_server.R            same for the simulator
ui/ ...                 same split for the UI files (later)
```

Each part file starts with a header naming what it defines for others and what it
reads from other parts; a test regenerates that list from the code and fails when
the header is out of date, so the interface cannot drift silently.

Stage A (this refactor): parts are evaluated in the server function's own
environment, in the original order. Nothing about scoping, observer creation
order or ids changes, which is what makes the step provably behaviour-preserving
(the parsed statements of the parts are compared with the parsed statements of
the original closure).

Stage B (later, panel by panel): a leaf panel becomes
`smlm_ripley_server(input, output, session, core)` where `core` is a named list
of the reactives it declared as reads. The declared-reads header of Stage A is
exactly the argument list of Stage B.

Stage C (optional, needs a browser test rig): `moduleServer` + `NS` for panels
whose ids are not in the settings spec.

### 4.2 Migration order, lowest risk first

| # | Step | Risk | Check |
|---|---|---|---|
| 1 | Golden-output and smoke tests | none (tests only) | they pass on the untouched code |
| 2 | Move the pure closure functions of the simulator to R/sim_*.R | very low: no reactive touched | parse-level identity of the moved functions, goldens |
| 3 | Move the pure closure functions of tracking to R/ | very low | same |
| 4 | Split the three closures into ordered part files | low: same statements, same order, same environment | statement-by-statement identity, goldens, smoke |
| 5 | Extract estimator bodies from tracking bundles into R/ functions (CVE diagnostics, MSD ensemble, confinement, displacement, directionality, transport, per-state MSD, TSA) | medium: code moves across a function boundary | goldens for each bundle before and after |
| 6 | Extract the simulator bodies (`sim_result`, `sim_trk_result`) into R/ | medium: RNG call order must be preserved exactly | simulate with fixed seeds before and after, `identical()` |
| 7 | Extract SMLM inline computations (trace counts, detection over time, precision-cloud binning, xi) | medium | goldens |
| 8 | One Run-handler frame for the SMLM batch panels (P5) | medium | goldens |
| 9 | Stage B for leaf panels | medium | smoke tests per panel |
| 10 | Unify the two SMLM filter chains (P1) | **changes results**: needs Q1 answered | not started without a decision |
| 11 | Shared overlay / scale bar / z-calibration code for both tabs (P4) | medium-high: UI behaviour | manual |
| 12 | Split `scatter3D`, `scatter2D`, `trk_track2D` into figure builders | high: 2,500 lines of plot code with side effects, no automatic check | manual |
| 13 | Tracking settings spec like `SMLM_SETTINGS_SPEC`, versioned file | medium, file format must stay readable | round-trip tests |
| 14 | Stage C | high | needs a browser test rig |

---

## 5. Questions and scientific issues for Viktor

Nothing in this section was changed in the code.

### Questions (two code paths, intended one unclear)

**Q1. Which SMLM filter chain is the intended one?** (P1.) Concretely: (a) should
the batch panels apply the spread filter the user set with the x/y/z range
sliders, or the `target_sd_x/y` numbers? (b) spread before length, or length
before spread? A trace that loses points to the spread stage cannot exist, but
the order matters once trimming is on, because the batch chain counts length
before the spread test and the active chain after. (c) should a 2D/3D crop and
rogue removal restrict the batch panels for the active file? (d) should a range
slider left at the data range be a no-op in the batch chain as well? The slider
value is rounded by the widget, so "always apply" can drop the localizations at
the extreme of each metadata column.
*Test that shows it:* load one file, set a z spread range, enable rogue removal,
press Plot; compare `sort(unique(plottable_data()$tid))` with
`sort(unique(prepared_df(active)$tid))`. They should be equal if the chains are
meant to agree.

**Q2. MSD fit window lower bound: 4 or 5?** (B2.) The UI default is 5. Which value
did earlier reported numbers use? Only matters where the input was NULL.

**Q2b. Should the first σ_loc run publish?** (B4.) If yes, removing
`ignoreInit = TRUE` at server_tracking.R:659 fixes it; it changes what the
noise-subtracting panels compute after a single press, so it is yours to decide.

**Q3. Is the 200 ms debounce at server_tracking.R:3497 wanted?** (B1.) Fixing it
changes how fast playback parameters apply, not any number.

### Scientific issues

**S1. Change-point segmentation is not noise-aware, and the test that guards it has
no dark cycles.** `cpd_segment_tracks()` (R/est_cpd.R:64-80) runs PELT on
`log(r²/Δt)`. With localization noise, `E[r²/Δt] = 2dD + 2dσ²/Δt`: the noise term
depends on Δt. A stretch of a track with many dark cycles (large Δt) therefore has
a lower mean of `log(r²/Δt)` than a stretch of single-cycle steps at the same D,
and PELT can cut there. The effect is largest when `σ²/t_c` is comparable to D,
which is the slow-molecule regime where a state change is most interesting, and
dark cycles cluster towards the end of a track (bleaching), so the artefact would
look like "slows down before disappearing". The noise is subtracted only after
segmentation (line 97). The guarding test (test-estimators-tracking.R:86-94) uses
`gap_p = 1`, that is, no dark cycles at all.
*Proposed test:* `sim_ct_state_tracks()` with one state, D = 20 nm²/ms, σ = 6 nm,
`gap_p = c(0.5, 0.25, 0.15, 0.1)` and a second variant where gaps are concentrated
in the last third of each track. Expect the fraction of tracks cut to stay under
the 15% of the existing test, and no correlation between segment D (before noise
subtraction) and segment `mean(1/Δt)`. I expect both to fail.

**S2. Immobile segments drop out of the state count.** A segment's D is clamped at
0 after noise subtraction (R/est_cpd.R:97-98), and `cpd_gmm()` then takes
`log10(D)` and discards non-finite values (lines 118-119). For a truly immobile
segment the noise-subtracted D is negative about half the time, so about half of
the immobile population is removed before the mixture is fitted, and the other
half enters with a D that is a positive noise fluctuation. The immobile weight is
under-counted and its fitted D has no meaning. The D-map has the same clamp per
cell (R/est_dmap.R:102), which biases the mean D of slow regions upward.
*Proposed test:* 50/50 mixture of D = 0 and D = 200 nm²/ms, σ = 5 nm; compare the
fraction of segments entering the GMM and the fitted weights with 0.5.

**S3. The noise-subtracting estimators inherit the uncalibrated σ.** Already on your
list (2026-11-10): with no calibration profile, `a = sqrt(8)` gives a photon σ
about 1.8 times the CVE σ on GP67H-013 and JDD / D-map / CPD read D at a third of
the CVE value. Structurally (2.3 item 5) nothing tells the user that these panels
depend on a number another panel wrote. A warning needs no method change.
*Proposed test:* in `testServer`, run JDD before and after the CVE run and assert
the result changes only together with a visible "σ source" line.

**S4. Range filters at the data extreme.** Related to Q1(d): `hard_filter` in the
batch chain compares against the widget's rounded range. If the widget rounds the
maximum down, the batch panels lose the localizations with the largest EFO / CFR /
time that the main plot keeps. Small, but it is a silent difference between two
views of "the same" data.
*Proposed test:* a file whose `efo` maximum is not representable at the slider's
step; compare row counts of the two chains with all sliders untouched.

**S5. Results that depend on unseeded random numbers.** (3.3.) Bootstrap CIs,
envelopes, the crowding null, the geometry nulls and density thinning give
different numbers on each run, and thinning changes the point set of every batch
panel. A saved settings file therefore does not reproduce an exported table that
used any of them. Not a bias, but a reproducibility gap for a paper.
*Proposed test:* run the Ripley envelope twice in one session with the same
settings and compare `global_p`.

---

## 6. What was done

Branch `refactor/modules`, on top of `7e8b6f9`. No analysis method, default or
number was changed; the golden files written from the untouched code were never
regenerated.

| Commit | Step of 4.2 | What |
|---|---|---|
| `0879bd5` | - | This review and `tools/server_map.R` |
| `af5dfb9` | 1 | Golden-output and smoke tests |
| `f157e3e` | 2 | 30 simulator functions and presets moved from `sim_server()` to `R/sim_*.R` |
| `273f24a` | 3 | 39 stateless tracking functions moved from `tracking_server()` to `R/trk_*.R`, `R/est_transport.R` |
| `463a365` | 4 | The three server functions split into 65 ordered part files under `server/` |
| `3807136` | 4 | Part headers with generated DEFINES / READS blocks, `server/README.md`, structure tests |
| `173dee2` | 5 (part) | Six tracking computations moved from reactives into `R/`, with ground-truth tests |

### 6.1 Safety net (step 1)

- `tests/testthat/helper-app.R`: loads the app once; `ui_input_defaults()` reads the
  start value of every static widget from the UI tag tree, so a `testServer`
  session starts with the inputs the app starts with instead of the `%||%`
  fallbacks (B2); `expect_golden()` compares a fingerprint of each result (counts,
  sums, sums of squares, quantiles, first values per numeric vector; value tables
  for text) at relative tolerance 1e-8.
- `test-golden-tracking.R`: upload through the tab's own path, Plot Data, every
  Run button; 30 results of the main session, the cleaning/filter stages with
  non-default settings, calibration and Compare Cells, and the `R/est_*.R`
  estimators called directly.
- `test-golden-smlm.R`: the active-file chain (19 results), every batch panel on
  two files (33 results, envelopes on), and both chains with filters, trimming,
  rogue removal and a crop.
- `test-golden-sim.R`: SMLM simulation for the default and each structure
  generator, the strategy comparison and the tracking simulation, at fixed seeds.
  These pin the order of the random draws.
- `test-smoke-servers.R`: after all panels have been triggered, every output of
  each server (110 tracking, 139 SMLM, 24 simulator; ids read from the code) is
  rendered and its status (ok / silent / error) compared with the stored one.
- Mutation checks: changing a quantile probability from 0.25 to 0.2501 in the MSD
  curve, and one argument in a server call, each fail the tracking golden.

shinytest2 was not used: a headless browser is not usable on this machine. The
checks are `testServer` plus one HTTP fetch of the running app's UI (200, all
three tabs present).

### 6.2 Structure (steps 2 to 4)

| | Before | After |
|---|---|---|
| `server_smlm.R` | 13,996 | 147 (header and assembly) + 34 parts, 14,571 lines, largest 1,516 |
| `server_tracking.R` | 11,573 | 161 + 25 parts, 11,044 lines, largest 1,303 |
| `sim_server.R` | 4,254 | 24 + 6 parts, 3,220 lines, largest 1,220 |
| `R/` | 15,539 in 23 files | 17,884 in 35 files |
| `tests/testthat/` | 2,795 in 36 files | 3,756 in 43 files, plus 13 golden files (125 kB) |

The part totals include about 23 header lines per part. The one part over 1,500
lines is `server/smlm/19_scatter3d.R`, which is a single statement (the 3D render).

How the moves were verified, beyond the tests:
- Functions moved to `R/`: each was checked with `codetools::findGlobals` to read
  no `input`, `output`, `session` or closure object other than functions moved
  with it; the moved text parses to the same expression as the original
  statement; the closure's remaining statements are the original ones minus the
  moved ones, in order; no top-level name is defined in two files.
- The split: the parsed statements of all parts in file order, plus the
  statements left in the server function, are identical to the parsed statements
  of the original closure; and every non-blank line of the original body
  (comments included) is present in the parts, in order.

Parts are "Stage A" of section 4.1: evaluated in the server function's
environment in the original order. `server/parts.R` documents the mechanism and
`server/README.md` the working rules. Each part's header carries a generated
interface block; `test-server-parts.R` fails when a block is out of date, a part
exceeds 1,600 lines, or a server-level name is defined twice in a tab.

Three existing guard tests read the text of `server_smlm.R` (`test-smlm-cache-key.R`,
`test-smlm-settings-adapter.R`, `test-smlm-settings.R`). After the split two of
them failed and the third passed without looking at any server code. All three
now read the assembly file plus the parts (`smlm_server_files()` in
`setup-helpers.R`); their claims are unchanged. `tools/check_dependencies.R`
scans `server/` and `R/` as well. A guard that greps one file by name is a
pattern to watch when files move.

### 6.3 Computation moved out of reactives (step 5, partly)

| Function | File | Was in |
|---|---|---|
| `dir_compute()` | `R/est_directionality.R` | `dir_bundle` |
| `cycle_class_counts()` | `R/trk_stats.R` | `cycle_class_bundle` |
| `msd_unified_curves()` | `R/est_msd.R` | `msd_unified_bundle` |
| `msd_pertrack_fits()` | `R/est_msd.R` | `trk_unified_tdiff` |
| `cve_integrated_fraction()`, `cve_continuity_mean()` | `R/est_cve.R` | local functions of `trk_cve_bundle` |

`test-tracking-extracted.R` tests them against known answers. One of its tests
pins `msd_unified_curves()` to `msd_ensemble()`, so the duplication noted in P3
can no longer drift without a failure.

One difference in reactive wiring, with no effect on results: `trk_unified_tdiff`
now reads `trk_blur()` once per run instead of once per fitted track.

### 6.4 Settings files

No statement of the settings save/load code was changed (it was moved as text
into `server/smlm/04`, `06`, `08` and `server/tracking/09`). The existing
round-trip test `test-smlm-settings-file.R` (slow suite) passes. There is no
automated test that loads a settings file written by the released v0.7.0 build;
loading one once in the refactored app would close that.

## 7. Deliberately left

- **True Shiny modules (`moduleServer` + `NS`).** Not done, for the reasons in 4.1.
  The done-criterion "tabs/panels converted into Shiny modules" is met only in
  the sense of one file per responsibility with a checked interface, not in the
  sense of namespaced modules.
- **Parts are contiguous runs of the original file.** Statements were not
  reordered, so a few concerns are still spread over several parts (SMLM
  settings over `04`, `06`, `08`; the two filter-reset observers in `03` and
  `12`). Regrouping changes observer creation order and was left for a step
  that can be checked by hand in the browser.
- **The rest of step 5.** Still computing inside reactives: the per-track sigma
  and linking QC of `trk_cve_bundle`, `conf_ci_bundle`, `trk_disp_bundle`,
  `trk_disp_geom_null`, `trk_export_summary`, `seg_msd_rv`, `msd_cluster_rv`,
  `msd_cluster_models`, `trk_seq_params`, and the six copies of the MSD
  fit-window parsing.
- **Steps 6 to 14 of 4.2**: the simulator bodies (`sim_result`, `sim_trk_result`),
  SMLM inline computations, a shared Run-handler frame, the duplicated overlay /
  scale-bar / z-calibration code of the two tabs, splitting the three large plot
  renders, a settings spec for tracking.
- **UI files** (`ui_smlm.R` 2,628 lines, `ui_tracking.R` 1,632, `sim_ui.R` 1,132)
  are untouched: no automatic check of a UI change is available here.
- **Three small SMLM helpers** (`gt_exploratory_note`, `.small_n_note`,
  `.stale_badge`) and the multiplex parsers are stateless but were left in the
  server, in line with the earlier decision to keep UI-note builders there.
- **Dead code** (3.4) was not removed. Nothing listed there was deleted or
  changed; the never-sourced files are still in the repo root.
- **Bugs B1 to B5, questions Q1 to Q3 and issues S1 to S5** are recorded, not
  fixed. B4, B5 and the fallback values of B2 are pinned by the golden and smoke
  tests as current behaviour, so fixing one will (correctly) require
  regenerating the affected golden entry.
- The scan for global helper functions that are never called did not complete
  and is not reported.
- **File references in the other docs.** `docs/MINFLUX_analysis_methods.md`,
  `docs/simulator.md` and `docs/simulation_module.md` still name
  `server_smlm.R`, `server_tracking.R` and `sim_server.R` (some with line
  numbers that were already out of date) as the place of the code. The function
  and reactive names they give are unchanged and can be found with
  `tools/server_map.R`; the docs themselves were not edited.

## 8. Remaining risks

- **What the tests cannot see.** `testServer` does not run JavaScript, so
  `conditionalPanel` conditions, the plotly proxy calls of the play loop, relayout
  events, file inputs and `renderUI`-created widgets are not exercised. The split
  did not change any of that code, but it has not been clicked through. A manual
  pass in the browser (upload a file in each tab, Plot, play, one Run button
  per panel, save and load settings) is the missing check before merging.
- **Start-up cost.** The parts are parsed once at app start and evaluated per
  session with `eval()`. Session start was not timed.
- **Coverage gaps of the goldens.**
  - The simulated tracks all have the same length, so the rule that drops lags
    reached by under 10% of the tracks is never active (found by a mutation that
    did not fail). Tracks of mixed length would cover it.
  - 2D data only in tracking (no real z), no multiplex dataset, no image
    overlay, no sequence JSON, no calibration profile applied, no merged view.
  - `scatter3D` is not rendered in the smoke test (it requires `input$file`),
    so the shape fits and isosurfaces it computes are uncovered.
  - Render status only for plots and HTML: a figure that draws the wrong thing
    still counts as "ok".
- **Goldens are tied to this machine's numerics** at 1e-8. Another BLAS or R
  version may need a looser tolerance for the mixture and HMM fits.
- **Test run time** rose from about 35 s to about 3.5 min for the default suite
  (the app-level tests load every package and run every panel).
  `MINFLUXR_SKIP_APP_TESTS=true` gives the old behaviour.
- **Tool heuristics.** The interface blocks come from a static scope analysis.
  A name used only through `get()`, `session$env` or a string is not seen. None
  was found in the server code. The block is a checked description, not an
  enforcement: all parts of a tab can still reach each other's objects.
- **Tests that reach into the closure** (`prepared_df`, `trk_cve_bundle`, ...
  by name inside `testServer`) keep working because the parts share the server
  function's environment. They will need touching when a part becomes a
  function with arguments.

## 9. Fixes made after the review (2026-11-16)

Decided with Viktor after reading sections 3 and 5. Each is one commit; the golden
entries that moved are named in the commit message, and nothing else moved.

| Item | Commit | What changed |
|---|---|---|
| B4 | `7243aac` | sigma_loc is published on the first run. `ignoreInit = TRUE` removed from the CVE observer and from the MSD observer, which had the same defect. |
| B5 | `e12b950` | The change-point CSV export classifies only segments with a finite log10 D; the others get NA. |
| P1 / Q1 / S4 | `b4eea25` | One SMLM filter chain. See below. |
| docs | `409ad95` | The method and simulator docs point at `server/<tab>/` and `R/`. |

### 9.1 The SMLM filter chain (P1)

Intended order, as stated by Viktor: metadata ranges, outlier trimming, trace
length, trace spread; and the filters set in the Display tab apply to every
loaded file.

What the review missed and the fix found: `target_sd_x` / `target_sd_y`, which the
batch chain used for its spread filter, have **no widget** in the UI. In the app
the spread filter therefore never reached any all-files panel. And because the
batch chain always applied the metadata ranges, sliders left untouched clipped
every other file to the active file's range.

Now (R/smlm_filters.R, "THE SMLM FILTER STAGES"):
- `smlm_filter_meta()`, `smlm_filter_trim()`, `smlm_filter_length()`,
  `smlm_filter_spread()` are the only implementation of each stage. The active
  file's reactives call them one by one (the panels need the intermediate
  results); `apply_meta_filters_batch()` composes them.
- The server filters a stored file through `batch_filters()`
  (server/smlm/23_prepared_cache.R), which passes the spread ranges (`sd_filt`)
  and the active dataset's full metadata ranges (`meta_full_rng`). Both are in
  the cache key of `prepared_df()`.
- A metadata slider equal to the active dataset's full range is untouched and
  filters no file. Once narrowed it applies to all files as numbers.
- Length before or after spread gives the same result: both drop whole traces
  judged on the same trimmed points. The active chain's reactives were left in
  their existing order (spread, then length) because the length slider takes its
  range from the spread-filtered data; the composed chain uses the stated order.
- `target_sd_x` / `target_sd_y` stay in the settings spec so older settings files
  and scripts read as before. They are used only when `apply_meta_filters_batch()`
  is called without spread ranges, which the app no longer does.

Decided afterwards:
- **Rogue removal was removed from the app.** It predates the spread filter and
  the outlier trimming, which cover its purpose. The checkbox, its two controls,
  the step in `plottable_data`, `remove_rogue_traces()` and the settings entries
  are gone. Settings files that carry `remove_rogues`, `rogue_rule`,
  `rogue_factor` still load; those keys are ignored.
- **Crops stay per dataset.** A crop (3D box/sphere, 2D crop to view) is a choice
  made on the main plot for one dataset and is not applied to the all-files
  panels. A golden test asserts the crop is the only difference: inside the crop
  box the two chains hold the same localizations.

Remaining differences, by design:
- `make_positive` shifts coordinates for display in the main plot only.
- Density thinning applies to the all-files chain only (it exists to make files
  comparable).

### 9.2 Change-point segmentation and dark cycles (S1, S2)

Tested instead of argued: `tests/validation/cpd_dark_cycle_clustering.R`.
One free population, fraction of tracks cut:

| D (nm²/ms) | σ (nm) | no dark cycles | spread evenly | clustered late |
|---|---|---|---|---|
| 20 | 3 | 0.125 | 0.090 | 0.290 |
| 20 | 6 | 0.105 | 0.130 | 0.590 |
| 215 | 3 | 0.065 | 0.070 | 0.065 |
| 215 | 6 | 0.060 | 0.045 | 0.120 |

So S1 as written was too broad. Evenly spread dark cycles do no harm. A false
state change appears only when dark cycles cluster in one part of the track AND
localization noise dominates the step (slow molecules). At D around 0.2 µm²/s
the effect is small. Nothing was changed in the estimator.

S2 (segments with D clamped to 0 leave the mixture) matters only for a sample
with an immobile fraction: in a 50/50 immobile/mobile simulation 55% of the
immobile segments were dropped and the mixture put the immobile share at 0.35
instead of 0.53. Not changed.
