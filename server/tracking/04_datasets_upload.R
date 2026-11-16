# ============================================================================
# DATASETS - uploads add files to a store; one VIEW of it is active
# ============================================================================
# Every uploaded file is read once and kept on disk (trk_store). The active
# view is either one file or all loaded files merged (track ids then prefixed
# F1_, F2_, ... so they cannot collide). Activating a view runs the ingest
# below - chunking, frozen raw statistics, the Filter & Rechunk reset point -
# exactly as an upload used to, so everything downstream sees a fresh dataset.
# Removing a file drops it from the store and activates another one.
trk_store     <- reactiveVal(list())          # name -> list(path, n, json_itr)
trk_view      <- reactiveVal(character(0))    # names of the files in the active view
trk_store_dir <- file.path(tempdir(), paste0("trk_store_", session$token))
.trk_last_active <- ""                        # last ingested view (guards the selector echo)
.trk_ingest_n    <- 0L                        # makes every ingest's chunk dir unique

# TRUE once there is a dataset to show (any view, or an adopted chunk manifest).
trk_has_data <- reactive(length(trk_view()) > 0L || !is.null(rv$manifest_dir))

# Name of the active view for labels and download file names.
trk_view_label <- reactive({
  v <- trk_view()
  if (!length(v)) "tracking" else if (length(v) == 1L) v
  else sprintf("%d files merged (%s, …)", length(v), v[1])
})
trk_view_base <- reactive({
  v <- trk_view()
  if (!length(v)) return("tracking")
  b <- tools::file_path_sans_ext(v[1])
  if (length(v) > 1L) paste0(b, "_merged", length(v)) else b
})


# ── .trk_ingest() ──────────────────────────────────────────────────────────
# GOAL: make one data.frame the active dataset: per-track summary, frozen raw
# statistics (histogram backgrounds, filter reset), chunks on disk, manifest.
# This is the body the upload observer used to run on its merged data.
.trk_ingest <- function(df, label) {
  rv$bigfile_notice <- FALSE
  rv$chunking <- TRUE
  on.exit({ rv$chunking <- FALSE }, add = TRUE)
  .trk_ingest_n <<- .trk_ingest_n + 1L
  safe_base <- gsub("[^A-Za-z0-9_-]+", "_", label)
  out_dir <- file.path(tempdir(), sprintf("chunks_%s_%d_%d", safe_base, as.integer(Sys.time()), .trk_ingest_n))

  # Compute per-TID summary (duration in ms, n_pts) and sort longest-first
  # tim * 1000 matches the ms conversion trk_raw_df() applies after chunk loading
  tid_sum <- df %>%
    dplyr::group_by(tid) %>%
    dplyr::summarise(
      duration_ms = (max(tim, na.rm = TRUE) - min(tim, na.rm = TRUE)) * 1000,
      n_pts        = dplyr::n(),
      .groups      = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(duration_ms))

  # Reorder df rows so chunks end up sorted by duration (longest first)
  df <- df[order(match(df$tid, tid_sum$tid)), , drop = FALSE]

  # Freeze raw stats for histogram grey bars and filter reset.
  # Done once per dataset - never overwritten by chunk switching or rechunking.
  # tim stored in ms (raw file seconds; trk_raw_df multiplies by 1000).
  loc_vec <- function(col) if (col %in% names(df)) df[[col]][is.finite(df[[col]])] else numeric(0)
  trk_raw_stats_rv(list(
    tid_length   = as.numeric(table(df$tid)),
    tid_duration = tid_sum$duration_ms,
    efo          = loc_vec("efo"),
    cfr          = loc_vec("cfr"),
    dcr          = loc_vec("dcr"),
    tim          = loc_vec("tim") * 1000,
    eco          = loc_vec("eco"),
    dwell        = loc_vec("dwell"),
    fluo         = loc_vec("fluo")
  ))
  trk_rechunk_stats_rv(NULL)  # clear any leftover rechunk stats from prior dataset
  trk_rechunk_bbox(NULL)      # clear locked spatial extent from prior dataset
  trk_spatial_filter(NULL)    # a crop drawn on another dataset does not apply here
  trk_chunk_pref("__all__")   # a new dataset does not inherit a chunk scoping

  user_tpc <- .tids_per_chunk()

  manf <- split_by_tid_to_dir(
    df, out_dir,
    tids_per_chunk = user_tpc,
    fmt = CHUNK_FORMAT
  )
  if (file.exists(manf)) {
    rv$manifest_dir          <- dirname(manf)
    rv$original_manifest_dir <- rv$manifest_dir
    trk_data_key(as.character(as.numeric(Sys.time())))
    # Save per-TID summary alongside the manifest
    readr::write_csv(tid_sum, file.path(rv$manifest_dir, "tid_summary.csv"))
    rv$chunking <- FALSE
  } else {
    # Fallback: compress to a single chunk, in whatever format the app writes.
    # This used to call qs::qsave() unconditionally, which crashed the upload on
    # any machine without qs even though CHUNK_FORMAT had already fallen back.
    out_one <- file.path(out_dir, paste0(safe_base, ".", CHUNK_FORMAT))
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    if (identical(CHUNK_FORMAT, "qs") && requireNamespace("qs", quietly = TRUE))
      qs::qsave(df, out_one, preset = "high")
    else
      fst::write_fst(as.data.frame(df), out_one, compress = 50)
    man <- data.frame(chunk = 1L, n_tids = length(unique(df$tid)), file = out_one)
    readr::write_csv(man, file.path(out_dir, "manifest.csv"))
    readr::write_csv(tid_sum, file.path(out_dir, "tid_summary.csv"))
    rv$manifest_dir          <- out_dir
    rv$original_manifest_dir <- out_dir
    trk_data_key(as.character(as.numeric(Sys.time())))
    rv$chunking <- FALSE
  }
  invisible(NULL)
}

# ── .trk_activate() ────────────────────────────────────────────────────────
# GOAL: build the data of one view from the store and ingest it.
# sel : a file name, or "__merged__" for all loaded files
.trk_activate <- function(sel) {
  store <- isolate(trk_store())
  files <- if (identical(sel, "__merged__")) names(store) else intersect(sel, names(store))
  if (!length(files)) return(invisible(NULL))
  dfs <- lapply(seq_along(files), function(i) {
    d <- read_chunk_file(store[[files[i]]]$path)
    if (is.null(d)) return(NULL)
    if (length(files) > 1L) d$tid <- paste0("F", i, "_", as.character(d$tid))
    d
  })
  df <- dplyr::bind_rows(Filter(Negate(is.null), dfs))
  if (!nrow(df)) {
    showNotification("The selected dataset could not be read.", type = "error", duration = 8)
    return(invisible(NULL))
  }
  .trk_last_active <<- sel
  trk_view(files)
  withProgress(message = sprintf("Preparing %s...", if (length(files) > 1L)
                                   sprintf("%d merged files", length(files)) else files),
               value = 0.5, .trk_ingest(df, isolate(trk_view_base())))
  if (length(files) > 1L)
    showNotification(sprintf("Merged %d files into one dataset (%s localizations).",
                             length(files), format(nrow(df), big.mark = ",")),
                     type = "message", duration = 6)
  invisible(NULL)
}

# Point the selector at a view and ingest it directly. Setting the selector to
# the value it already holds does not fire its observer (e.g. re-uploading a
# file of the same name), so the ingest cannot be left to the observer.
.trk_set_active <- function(sel) {
  updateSelectInput(session, "trk_active", choices = .trk_view_choices(isolate(trk_store())),
                    selected = sel)
  .trk_activate(sel)
}

# ---- Upload: read every file into the store -------------------------------
observeEvent(input$trk_file, {
  req(input$trk_file)
  n_files   <- nrow(input$trk_file)
  paths     <- input$trk_file$datapath
  names_o   <- as.character(input$trk_file$name)          # as uploaded, for display
  names_v   <- tolower(names_o)                           # for extension checks
  size_bytes <- sum(suppressWarnings(as.numeric(input$trk_file$size)), na.rm = TRUE)
  if (!is.finite(size_bytes) || size_bytes <= 0) {
    info <- tryCatch(file.info(paths), error = function(e) NULL)
    if (!is.null(info)) size_bytes <- sum(info$size, na.rm = TRUE)
  }
  size_mb <- round(as.numeric(size_bytes) / (1024^2), 1)

  # If a single file is uploaded next to an existing manifest, adopt it and bail.
  # (Skipped for multi-file uploads.) The adopted chunks bypass the store.
  if (n_files == 1L) {
    dirp <- dirname(paths[1])
    if (file.exists(file.path(dirp, "manifest.csv"))) {
      rv$manifest_dir          <- dirp
      rv$original_manifest_dir <- dirp
      trk_data_key(as.character(as.numeric(Sys.time())))
      rv$chunking <- FALSE
      if (!rv$bigfile_notice) {
        showNotification("Found existing chunk manifest; using chunked loading.", type = "message", duration = 6)
        rv$bigfile_notice <- TRUE
      }
      return(invisible(NULL))
    }
  }

  dir.create(trk_store_dir, showWarnings = FALSE, recursive = TRUE)
  store <- trk_store(); added <- character(0); json_itr <- NA_integer_
  keep_all <- c("tid","x","y","z","tim","itr","efo","cfr","dcr","eco","dwell","fluo")
  withProgress(message = sprintf("Loading %d file(s) (%.1f MB)...", n_files, size_mb), value = 0, {
    for (fi in seq_len(n_files)) {
      incProgress(1 / n_files, detail = names_o[fi])
      di <- tryCatch(load_one_trk(paths[fi], names_v[fi]), error = function(e) NULL)
      if (is.null(di) || !is.data.frame(di) || !nrow(di)) {
        showNotification(sprintf("%s could not be read and was skipped.", names_o[fi]),
                         type = "warning", duration = 8)
        next
      }
      if (!all(c("tid","x","y","tim") %in% names(di))) {
        showNotification(sprintf("%s is missing required columns (tid, x, y, tim) and was skipped.",
                                 names_o[fi]), type = "error", duration = 8)
        next
      }
      ji <- attr(di, "json_itr"); if (!is.na(ji)) json_itr <- ji
      di$tid <- as.character(di$tid)
      # Keep only the columns tracking needs  -  reduces stored file size significantly.
      # No automatic CFR pre-filter: CFR is not a strong artefact indicator for
      # tracking, and the metadata filter offers a CFR range when it is wanted.
      di <- di[, intersect(keep_all, names(di)), drop = FALSE]
      old <- store[[names_o[fi]]]
      if (!is.null(old)) unlink(old$path)                    # same name: replaced
      f <- file.path(trk_store_dir, sprintf("f%03d_%s.fst", length(store) + length(added) + 1L,
                                            gsub("[^A-Za-z0-9_-]+", "_", tools::file_path_sans_ext(names_o[fi]))))
      fst::write_fst(di, f, compress = 50)
      store[[names_o[fi]]] <- list(path = f, n = nrow(di), json_itr = ji)
      added <- c(added, names_o[fi])
    }
  })
  if (!length(added)) return(invisible(NULL))
  trk_store(store)

  # Show JSON iteration notification (mirrors SMLM tab behaviour)
  if (!is.na(json_itr))
    showNotification(sprintf(".json or .npy import: using final localization (itr = %d) for x/y.",
                             as.integer(json_itr)), type = "message", duration = 8)
  if (length(store) > 1L)
    showNotification(sprintf("%d datasets loaded. Showing %s; switch or merge them under Active dataset.",
                             length(store), added[1]), type = "message", duration = 7)
  .trk_set_active(added[1])                                 # one file at a time by default
}, ignoreInit = TRUE)

# User picks another view.
observeEvent(input$trk_active, {
  sel <- input$trk_active
  req(nzchar(sel %||% ""))
  if (identical(sel, .trk_last_active)) return(invisible(NULL))   # already active
  .trk_activate(sel)
}, ignoreInit = TRUE)

# Remove the selected file from the store and show another one.
observeEvent(input$trk_remove, {
  sel <- input$trk_active %||% ""
  store <- trk_store()
  if (identical(sel, "__merged__")) {
    showNotification("Select a single dataset to remove it.", type = "warning", duration = 5)
    return(invisible(NULL))
  }
  if (!nzchar(sel) || is.null(store[[sel]])) return(invisible(NULL))
  unlink(store[[sel]]$path)
  store[[sel]] <- NULL
  trk_store(store)
  if (length(store)) {
    .trk_set_active(names(store)[1])
  } else {
    .trk_last_active <<- ""
    trk_view(character(0))
    rv$manifest_dir <- NULL; rv$original_manifest_dir <- NULL
    trk_raw_stats_rv(NULL)
    updateSelectInput(session, "trk_active", choices = character(0))
    trk_data_key(as.character(as.numeric(Sys.time())))
  }
  showNotification(sprintf("%s removed.", sel), type = "message", duration = 4)
})

# ── Adopt an axial factor recorded in an uploaded CSV ─────────────────────
# Reads a small head of each file rather than the whole thing, so this stays
# cheap even for a multi-GB export. A sample rather than a single row because
# z VARIATION is what decides whether the calibration applies at all, and one
# row cannot show variation - a planar file must be passed over in silence.
#
# Done here rather than inside trk_raw_loaded() because that is a reactive:
# writing to an input from inside one is a side effect in the wrong place.
#
# Only `certain` factors are adopted - i.e. a real z_scale_applied column. The
# 0.7 assumed for pyMINFLUX and legacy exports is a guess, and a guess must not
# silently become the session's calibration.
observeEvent(input$trk_file, {
  req(input$trk_file)
  adopt <- NULL
  for (i in seq_len(nrow(input$trk_file))) {
    if (!grepl("\\.csv(\\.gz)?$", tolower(input$trk_file$name[i]))) next
    head_df <- tryCatch(
      as.data.frame(data.table::fread(input$trk_file$datapath[i], nrows = 2000L,
                                      showProgress = FALSE)),
      error = function(e) NULL)
    if (is.null(head_df) || !nrow(head_df)) next
    if (!has_axial_data(head_df, min_sd = 1)) next   # planar: nothing to adopt
    cal <- detect_csv_z_calibration(head_df)
    if (isTRUE(cal$certain) && is.finite(cal$factor) && cal$factor > 0) adopt <- cal$factor
  }
  if (is.null(adopt) || isTRUE(all.equal(adopt, trk_z_applied()))) return()

  updateNumericInput(session, "trk_z_scale", value = adopt)
  # Committed straight away: an upload is itself an explicit action, and the
  # pipeline is being rebuilt anyway, so there is nothing to defer to Plot Data.
  trk_z_applied(adopt)
  showNotification(sprintf("Axial calibration %.3f adopted from the uploaded file.", adopt),
                   type = "message", duration = 6)
}, ignoreInit = TRUE, priority = 20)   # before the loaders read the data

# ---- Load CSV/fast formats OR a pre-split chunk ----
trk_raw_loaded <- reactive({
  req(trk_has_data())
  required <- c("tid","x","y","tim")
  
  # Source for the WHOLE filter pipeline: every chunk of the ORIGINAL upload of the
  # active dataset (rv$original_manifest_dir), never a single chunk. Plot Data
  # rechunks the filtered survivors into rv$manifest_dir for the playback chunk
  # selector; the pipeline must not read that set, or it runs a second time on
  # already-filtered data (removed localizations leave holes that the gap split
  # then cuts at), and results would depend on whether Plot Data was pressed.
  df <- NULL
  src_dir <- rv$original_manifest_dir %||% rv$manifest_dir
  if (!is.null(src_dir)) {
    manf <- file.path(src_dir, "manifest.csv")
    if (file.exists(manf)) {
      man <- tryCatch(readr::read_csv(manf, show_col_types = FALSE), error = function(e) NULL)
      if (!is.null(man) && nrow(man) > 0 && "file" %in% names(man)) {
        dfs <- lapply(man$file, read_chunk_file)
        df  <- dplyr::bind_rows(Filter(Negate(is.null), dfs))
      }
    }
  }
  # Every dataset is ingested into chunks (.trk_ingest), so there is always a
  # manifest once anything is loaded.
  shiny::validate(shiny::need(!is.null(df) && nrow(df) > 0, "No dataset loaded."))
  
  # --- common post-processing (applies to both chunk & single-file) ---
  miss <- setdiff(required, names(df))
  shiny::validate(shiny::need(length(miss) == 0, paste("Missing columns:", paste(miss, collapse = ", "))))
  
  df$tid <- as.character(df$tid)
  df$x   <- as.numeric(df$x)
  df$y   <- as.numeric(df$y)
  df$tim <- as.numeric(df$tim)
  
  # tim is always in seconds  -  convert to milliseconds
  df$tim <- df$tim * 1000

  df
})

# ── Physical Z correction factor (applied at the pipeline ROOT) ──────────────
# A refractive-index / mounting-media mismatch scales the apparent axial (z)
# coordinate by a constant. Applying it in a thin layer directly over the chunk
# loader means EVERY downstream reactive - track cleaning, gap-splitting, the CVE
# σ_loc estimator, trk(), and the 2D/3D/side-view display - reads corrected z from
# this single point, so the correction reaches all z-dependent parameters at once.
#
# Kept SEPARATE from trk_raw_loaded() so changing the factor re-runs only this
# cheap rescale, never the expensive chunk reload. The z dependency is intentionally
# NOT isolated: it must invalidate the whole filter pipeline so the cached stages
# (which trk() reads via isolate()) recompute with the new z rather than serving a
# stale cached value. z is not a filter/histogram variable, so rescaling it here has
# no effect on the metadata panels or the raw-stats ranges.
trk_raw_df <- reactive({
  df <- trk_raw_loaded()
  # trk_raw_loaded() is RAW, uncorrected z for every format (readers use
  # z_scale = 1; CSVs are stripped on load), so this is the ONE place the
  # axial calibration enters the tracking pipeline.
  zf <- trk_z_applied()
  if (is.finite(zf) && zf > 0 && zf != 1.0 && "z" %in% names(df) && any(is.finite(df$z)))
    df$z <- df$z * zf
  df
})

# TRUE, frozen raw - the ORIGINAL uploaded data, NEVER touched by the filters or by
# Filter & Rechunk. trk_raw_df() reads the CURRENT manifest (which becomes the
# FILTERED set after a rechunk), so it can't be used as the "raw" reference. This
# reactive always reads rv$original_manifest_dir (frozen at upload, restored by
# Reset Ranges), so "raw" stays the uploaded data (e.g. 52% integrated) regardless
# of any filtering/rechunk. Used for the RAW side of the σ_loc estimator and the
# TSA before/after. Only recomputes on a genuinely new upload.
trk_raw_original <- reactive({
  dir <- rv$original_manifest_dir
  # Fallback: unchunked / manifest-less loads → best-effort current raw
  if (is.null(dir) || !file.exists(file.path(dir, "manifest.csv")))
    return(tryCatch(trk_raw_df(), error = function(e) NULL))
  man <- tryCatch(readr::read_csv(file.path(dir, "manifest.csv"), show_col_types = FALSE),
                  error = function(e) NULL)
  if (is.null(man) || !nrow(man) || !"file" %in% names(man)) return(NULL)
  df <- dplyr::bind_rows(Filter(Negate(is.null), lapply(man$file, read_chunk_file)))
  if (is.null(df) || !nrow(df) || !all(c("tid","tim") %in% names(df))) return(NULL)
  df$tid <- as.character(df$tid)
  for (cc in intersect(c("x","y","tim"), names(df))) df[[cc]] <- as.numeric(df[[cc]])
  df$tim <- df$tim * 1000                  # s → ms (match trk_raw_df post-processing)
  df
})
