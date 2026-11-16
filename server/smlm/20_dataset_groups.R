## ── Dataset Groups (shared across Analysis-tab analyses) ──────────────────
# Keep the 6 group file-selectors populated with the loaded file names.
observeEvent(smlm_file_store(), {
  fs <- names(smlm_file_store())
  for (i in 1:6) {
    sel <- isolate(input[[paste0("grp_files_", i)]])
    updateSelectInput(session, paste0("grp_files_", i),
                      choices = fs, selected = intersect(sel, fs))
  }
}, ignoreNULL = FALSE)


# file -> group-name map (only the first n_groups slots, only assigned files)
file_group_map <- reactive({
  ng <- .int_input(input$n_groups, 6L, 1L, 6L)
  m <- character(0)
  for (i in seq_len(ng)) {
    files <- input[[paste0("grp_files_", i)]]
    if (!length(files)) next
    nm <- input[[paste0("grp_name_", i)]]
    if (is.null(nm) || !nzchar(trimws(nm))) nm <- paste("Group", i)
    m[files] <- trimws(nm)
  }
  m
})

# group-name -> colour (from the per-group colour pickers)
.GROUP_DEFAULT_COLORS <- c("#1f77b4", "#d62728", "#2ca02c", "#9467bd", "#ff7f0e", "#17becf")
group_color_map <- reactive({
  ng <- .int_input(input$n_groups, 6L, 1L, 6L)
  cols <- character(0)
  for (i in seq_len(ng)) {
    nm <- input[[paste0("grp_name_", i)]]
    if (is.null(nm) || !nzchar(trimws(nm))) nm <- paste("Group", i)
    cols[trimws(nm)] <- input[[paste0("grp_color_", i)]] %||% .GROUP_DEFAULT_COLORS[i]
  }
  cols
})

# Colour palette for a set of series labels. In group mode, use the group
# colour pickers (so every analysis is consistent); otherwise the auto palette.
series_palette <- function(labels) {
  labels <- unique(as.character(labels))
  if ((input$analysis_by %||% "file") == "group") {
    gm  <- group_color_map()
    out <- unname(gm[labels])
    miss <- is.na(out)
    if (any(miss)) {
      fb <- if (sum(miss) <= 8) RColorBrewer::brewer.pal(max(3L, sum(miss)), "Set2")[seq_len(sum(miss))]
            else viridisLite::viridis(sum(miss))
      out[miss] <- fb
    }
    return(stats::setNames(as.character(out), labels))
  }
  # File mode: CATEGORICAL colours. Files are nominal, so a sequential ramp
  # (the old Set2-then-viridis) was the wrong encoding - it made every file a
  # shade of blue-green with neighbours hard to tell apart.
  stats::setNames(qual_palette(length(labels)), labels)
}

# Per-FILE colours that respect the grouping: in group mode each file becomes a
# shade of its group's colour (so the hue reads as the condition and the shade as
# the replicate); in file mode each file gets its own categorical colour.
# Returns a named vector indexed by FILE name.
file_palette <- function(files) {
  files <- unique(as.character(files))
  if (!length(files)) return(character(0))
  if ((input$analysis_by %||% "file") != "group")
    return(stats::setNames(qual_palette(length(files)), files))
  gm <- file_group_map()
  gc <- group_color_map()
  out <- stats::setNames(rep(NA_character_, length(files)), files)
  for (g in unique(unname(gm[intersect(files, names(gm))]))) {
    fs   <- intersect(files, names(gm)[gm == g])
    base <- unname(gc[g])
    if (is.na(base)) base <- QUAL_COLORS[1]
    out[fs] <- shade_colors(base, length(fs))
  }
  miss <- is.na(out)
  if (any(miss)) out[miss] <- .adjust_lum(qual_palette(sum(miss)), 0.25)
  out
}

# ── CHANNEL FACET for univariate analyses ─────────────────────────────────
#
# Splits an already-prepared per-file frame into one subset per channel, so a
# distribution can be read per channel instead of pooled. Returns a list of
# (label, df); with the facet off, or on single-channel data, that list is one
# entry with an empty label and the caller behaves exactly as before.
#
# WHERE THIS IS SAFE, AND WHERE IT IS NOT.
#   Safe for DISTRIBUTIONS that do not divide by an area: NN distances, radius
#   of gyration, molecules per Median Cluster. Subsetting changes what is in
#   the sample but not what the statistic means.
#
#   NOT safe, as written, for Ripley's L or g(r). Those divide by the observed
#   area, and the observation window is currently derived from the points that
#   are passed in - so a channel with fewer localizations would get a smaller
#   footprint, a higher apparent density, and a DEFLATED curve. The window has
#   to be built once from the pooled file and then shared by both channels,
#   which is a change to spatial_gate_close()'s plumbing rather than a call to
#   this function. Until that exists the spatial statistics stay pooled: a
#   pooled curve is honest, a per-channel curve with a per-channel window is
#   not.


# Relabel a results data frame's `dataset` column to its group. In group mode,
# files not assigned to any group are dropped. Analyses call this at read time,
# so toggling files<->groups updates plots without re-running.
apply_dataset_groups <- function(res, col = "dataset") {
  if (is.null(res) || !nrow(res) || !col %in% names(res)) return(res)
  if ((input$analysis_by %||% "file") != "group") return(res)
  m <- file_group_map()
  if (!length(m)) return(res[0, , drop = FALSE])
  res <- res[res[[col]] %in% names(m), , drop = FALSE]
  if (nrow(res)) res[[col]] <- unname(m[res[[col]]])
  res
}

## ── Per-file summary export (unit of replication = the file / cell) ────────
## GOAL: one CSV with per-file summaries of every distribution analysis that has
## been run, tagged with each file's group, so between-group statistics use one
## value per cell. Group-mode histograms pool units across files, which is fine
## for looking but pseudo-replicates any test run on them; see per_file_summary()
## in R/smlm_summary.R. Reads the stored results at download time, so it reflects
## the last Run of each panel. The result stores are defined further down this
## function; they exist by the time the handler can fire.
output$download_per_file_summary <- downloadHandler(
  filename = function() paste0("per_file_summary_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
  content  = function(file) {
    rv <- function(f) tryCatch(f(), error = function(e) NULL)   # a store may be unset
    cross <- rv(cross_nn_results)
    if (is.data.frame(cross) && nrow(cross) && "direction" %in% names(cross))
      cross$dataset <- paste(cross$dataset, cross$direction, sep = " · ")  # one row per direction
    results <- list(
      batch_nn         = list(df = rv(nn_batch_results),  metrics = "distance",   unit = "NN distance (nm)"),
      units_per_median_cluster = list(df = rv(tid_count_results), metrics = "n_tids", unit = "Median Cluster"),
      radius_of_gyration = list(df = rv(rog_results),     metrics = c("rog", "n_units"), unit = "cluster"),
      crowding         = list(df = rv(crowd_results),     metrics = c("nn_median", "time_shift", "n_tids"), unit = "Median Cluster"),
      within_cluster_nn_vs_size = list(df = rv(size_nn_results), metrics = c("nn_median", "size"), unit = "Median Cluster"),
      cross_channel_nn = list(df = cross,                  metrics = "distance",   unit = "cross-channel NN distance (nm)"))
    out <- per_file_summary(results, group_map = file_group_map(),
                            files = names(smlm_file_store() %||% list()))
    if (!nrow(out))
      out <- data.frame(note = "No distribution analysis has been run yet. Run a panel, then download again.")
    utils::write.csv(out, file, row.names = FALSE)
  })
