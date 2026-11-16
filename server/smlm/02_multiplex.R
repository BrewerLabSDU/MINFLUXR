# ── MULTIPLEXED MINFLUX ─────────────────────────────────────────────────────
# Build an N-channel combined dataset (overlay loaded files, or split the active
# file by TIM breakpoints) carrying mplex_channel / mplex_label / mplex_color.
# Registered as one synthetic entry in smlm_file_store so the plot and every
# Analysis tool operate on it unchanged. See plan: graceful-marinating-pudding.
# Plain name: it becomes part of downloaded file names (settings, CSV), where a
# symbol such as a star is rejected or mangled by some systems.
mplex_name   <- "Multiplex"
MPLEX_TID_OFF <- 1e10            # per-channel tid offset (> u4 tid max ~4.3e9)
mplex_spec_rv <- reactiveVal(NULL)  # list(channel,label,color,n_loc,n_trace)

MPLEX_MAX_BREAKS <- 6L              # up to 6 splits → 7 time-window channels
MPLEX_MAX_WIN    <- 7L              # up to 7 explicit windows (= 7 channels)
.mplex_parse_breaks <- function(txt) {
  if (is.null(txt) || !nzchar(txt)) return(numeric(0))
  v <- suppressWarnings(as.numeric(strsplit(txt, "[,;[:space:]]+")[[1]]))
  v <- sort(unique(v[is.finite(v)]))
  utils::head(v, MPLEX_MAX_BREAKS)
}
# Parse "start-end, start-end, …" into a sorted N×2 matrix of included windows.
.mplex_parse_windows <- function(txt) {
  if (is.null(txt) || !nzchar(txt)) return(matrix(numeric(0), ncol = 2))
  parts <- strsplit(txt, "[,;]")[[1]]
  rows <- lapply(parts, function(p) {
    v <- suppressWarnings(as.numeric(strsplit(trimws(p), "[-–:[:space:]]+")[[1]]))
    v <- v[is.finite(v)]
    if (length(v) >= 2) sort(v[1:2]) else NULL
  })
  rows <- Filter(Negate(is.null), rows)
  if (!length(rows)) return(matrix(numeric(0), ncol = 2))
  m <- do.call(rbind, rows)
  m <- m[order(m[, 1]), , drop = FALSE]
  utils::head(m, MPLEX_MAX_WIN)
}
# Assign a tim vector to channel indices per the active TIM method.
# Returns list(ch = integer channel per tim (NA = removed/wash), nch, method,
# breaks/wins). Used by the colors UI, the preview histogram and the builder.
# The method / breakpoints / windows default to the controls; the builder
# passes them explicitly when rebuilding from a settings file.
.mplex_tim_assign <- function(tim, method = input$mplex_tim_method,
                              breaks_txt = input$mplex_tim_breaks,
                              windows_txt = input$mplex_tim_windows) {
  method <- method %||% "breaks"
  if (identical(method, "windows")) {
    wins <- .mplex_parse_windows(windows_txt)
    ch <- rep(NA_integer_, length(tim))
    for (w in seq_len(nrow(wins))) {
      inw <- is.na(ch) & tim >= wins[w, 1] & tim <= wins[w, 2]
      ch[inw] <- w
    }
    list(ch = ch, nch = nrow(wins), method = "windows", wins = wins)
  } else {
    br <- .mplex_parse_breaks(breaks_txt)
    list(ch = findInterval(tim, br) + 1L, nch = length(br) + 1L,
         method = "breaks", breaks = br)
  }
}

# Keep the overlay file picker in sync with loaded files (excluding the
# synthetic multiplex entry itself).
observeEvent(smlm_file_store(), {
  nm <- setdiff(names(smlm_file_store()), mplex_name)
  updateSelectInput(session, "mplex_files", choices = nm,
                    selected = intersect(isolate(input$mplex_files), nm))
}, ignoreNULL = FALSE, ignoreInit = TRUE)

output$mplex_tim_info <- renderUI({
  df <- uploaded_df_rv()
  if (is.null(df) || !"tim" %in% names(df) || !nrow(df))
    return(helpText("Load a file with a tim column to split by time."))
  rng <- range(df$tim, na.rm = TRUE)
  helpText(sprintf("Active file tim range: %.2f – %.2f s (span %.2f s).",
                   rng[1], rng[2], diff(rng)))
})

output$mplex_files_colors_ui <- renderUI({
  fs <- input$mplex_files
  if (is.null(fs) || !length(fs))
    return(helpText("Pick one or more files above; each becomes a channel."))
  defs <- mplex_default_colors(length(fs))
  lapply(seq_along(fs), function(j) {
    lid <- paste0("mplex_flab_", j); cid <- paste0("mplex_fcol_", j)
    cur_l <- isolate(input[[lid]]); cur_c <- isolate(input[[cid]])
    fluidRow(
      column(7, textInput(lid, NULL,
        value = if (!is.null(cur_l) && nzchar(cur_l)) cur_l else fs[j])),
      column(5, colourpicker::colourInput(cid, NULL,
        value = if (!is.null(cur_c) && nzchar(cur_c)) cur_c else defs[j]))
    )
  })
})

output$mplex_tim_colors_ui <- renderUI({
  method <- input$mplex_tim_method %||% "breaks"
  df  <- isolate(uploaded_df_rv())
  has_tim <- !is.null(df) && "tim" %in% names(df) && nrow(df)
  if (identical(method, "windows")) {
    wins <- .mplex_parse_windows(input$mplex_tim_windows)
    nch  <- nrow(wins)
    deflab_fn <- function(j) sprintf("Ch%d (%.0f–%.0f s)", j, wins[j, 1], wins[j, 2])
  } else {
    br  <- .mplex_parse_breaks(input$mplex_tim_breaks)
    nch <- length(br) + 1L
    edges <- c(if (has_tim) min(df$tim, na.rm = TRUE) else 0, br,
               if (has_tim) max(df$tim, na.rm = TRUE) else NA_real_)
    deflab_fn <- function(j) sprintf("Ch%d (%.0f–%.0f s)", j, edges[j], edges[j + 1])
  }
  if (nch < 1L) return(helpText("Define at least one window / breakpoint above."))
  defs <- mplex_default_colors(nch)
  lapply(seq_len(nch), function(j) {
    lid <- paste0("mplex_tlab_", j); cid <- paste0("mplex_tcol_", j)
    cur_l <- isolate(input[[lid]]); cur_c <- isolate(input[[cid]])
    fluidRow(
      column(7, textInput(lid, NULL,
        value = if (!is.null(cur_l) && nzchar(cur_l)) cur_l else deflab_fn(j))),
      column(5, colourpicker::colourInput(cid, NULL,
        value = if (!is.null(cur_c) && nzchar(cur_c)) cur_c else defs[j]))
    )
  })
})

# Traces-over-time preview: one point per trace (median tim) binned over time,
# bars colored by the time-window channel they fall in, dashed lines at splits.
# Helps the user place breakpoints. Live on the active file + breakpoints input.
output$mplex_tim_hist <- plotly::renderPlotly({
  df <- uploaded_df_rv()
  shiny::validate(shiny::need(!is.null(df) && "tim" %in% names(df) && nrow(df),
                              "Load a file with a tim column."))
  tt <- df %>% dplyr::group_by(tid) %>%
    dplyr::summarise(tim = stats::median(tim, na.rm = TRUE), .groups = "drop")
  tt <- tt$tim[is.finite(tt$tim)]
  shiny::validate(shiny::need(length(tt) > 1 && diff(range(tt)) > 0,
                              "Not enough time spread to plot."))
  sp   <- .mplex_tim_assign(tt)
  nch  <- max(1L, sp$nch)
  cols <- mplex_default_colors(nch)
  for (j in seq_len(nch)) {
    v <- input[[paste0("mplex_tcol_", j)]]
    if (!is.null(v) && nzchar(v)) cols[j] <- v
  }
  rng   <- range(tt)
  nbins <- max(2L, as.integer(input$mplex_tim_bins %||% 60L))
  edges <- seq(rng[1], rng[2], length.out = nbins + 1L)
  h     <- graphics::hist(tt, breaks = edges, plot = FALSE)
  binch <- .mplex_tim_assign(h$mids)$ch         # channel per bin (NA = wash)
  barcols <- ifelse(is.na(binch), "#cccccc", cols[pmin(binch, nch)])
  # boundary lines: window edges (windows mode) or breakpoints (breaks mode)
  bounds <- if (identical(sp$method, "windows")) as.numeric(t(sp$wins)) else sp$breaks
  shapes <- lapply(bounds, function(b) list(
    type = "line", x0 = b, x1 = b, y0 = 0, y1 = 1, yref = "paper",
    line = list(color = "#333333", width = 2, dash = "dash")))
  plotly::plot_ly(
    x = h$mids, y = h$counts, type = "bar",
    marker = list(color = barcols),
    hovertemplate = "t = %{x:.0f} s<br>%{y} traces<extra></extra>"
  ) %>% plotly::layout(
    bargap = 0,
    xaxis = list(title = "time (s)"),
    yaxis = list(title = "# trace IDs"),
    showlegend = FALSE,
    margin = list(t = 10, b = 36, l = 48, r = 8),
    shapes = shapes
  )
})

# =========================================================================
# .mplex_build(cfg) - BUILD THE COMBINED MULTIPLEX DATASET
# =========================================================================
# GOAL
#   One builder for every way a multiplex is made: the Build button (cfg read
#   from the controls) and "Load multiplex settings" (cfg read from a file).
#   Taking the recipe as an argument, rather than reading the controls, is
#   what lets a settings file rebuild without waiting for the per-channel
#   label / colour widgets to be rendered first.
#
# ARGS  cfg : list(mode = "files" | "tim",
#                  files  = file names, one channel each      [files mode]
#                  tim_file, tim_method, tim_breaks, tim_windows [tim mode]
#                  labels, colors = one per channel (optional; defaults used))
# RETURNS TRUE on success; FALSE after telling the user what is missing.
# SIDE EFFECTS registers the result as the "Multiplex" store entry, makes it
#   active, remembers cfg (for Save) and re-applies remembered alignment offsets.
# =========================================================================
mplex_cfg_rv <- reactiveVal(NULL)   # recipe of the last successful build

# The controls, as a cfg. The TIM source is the active file, unless the active
# file IS the multiplex (a rebuild): then it is the file the last build split.
.mplex_cfg_from_inputs <- function() {
  mode <- input$mplex_mode %||% "files"
  if (identical(mode, "files")) {
    fs <- input$mplex_files %||% character(0)
    return(list(mode = "files", files = fs,
      labels = vapply(seq_along(fs), function(j) input[[paste0("mplex_flab_", j)]] %||% "", ""),
      colors = vapply(seq_along(fs), function(j) input[[paste0("mplex_fcol_", j)]] %||% "", "")))
  }
  act <- input$active_file %||% ""
  if (identical(act, mplex_name)) act <- isolate(mplex_cfg_rv())$tim_file %||% ""
  sp  <- .mplex_tim_assign(numeric(0))                  # only for the channel count
  list(mode = "tim", tim_file = act,
       tim_method  = input$mplex_tim_method %||% "breaks",
       tim_breaks  = input$mplex_tim_breaks %||% "",
       tim_windows = input$mplex_tim_windows %||% "",
       labels = vapply(seq_len(sp$nch), function(j) input[[paste0("mplex_tlab_", j)]] %||% "", ""),
       colors = vapply(seq_len(sp$nch), function(j) input[[paste0("mplex_tcol_", j)]] %||% "", ""))
}

.mplex_build <- function(cfg) {
  store <- smlm_file_store()
  spec  <- list()
  pick  <- function(v, j, default) {                     # label/colour j, or the default
    x <- if (length(v) >= j) v[[j]] else NULL
    if (is.null(x) || is.na(x) || !nzchar(x)) default else x
  }

  if (identical(cfg$mode, "files")) {
    fs <- cfg$files
    if (!length(fs)) {
      showNotification("Pick at least one file to overlay.", type = "warning"); return(FALSE)
    }
    miss <- setdiff(fs, names(store))
    if (length(miss)) {
      showNotification(paste0("Load these files first: ", paste(miss, collapse = ", ")),
                       type = "error", duration = 12); return(FALSE)
    }
    defs  <- mplex_default_colors(length(fs))
    parts <- list()
    for (j in seq_along(fs)) {
      d <- store[[fs[j]]]
      if (is.null(d) || !nrow(d)) next
      lab <- pick(cfg$labels, j, fs[j]); col <- pick(cfg$colors, j, defs[j])
      d <- d %>% dplyr::select(-dplyr::any_of(
        c("mplex_channel", "mplex_label", "mplex_color", ".dataset_tag")))
      attributes(d)[c("mplex_off", "mplex_spec", "mplex_store_name")] <- NULL
      d$tid <- suppressWarnings(as.numeric(d$tid)) + j * MPLEX_TID_OFF
      d$mplex_channel <- j; d$mplex_label <- lab; d$mplex_color <- col
      parts[[length(parts) + 1L]] <- d
      spec[[length(spec) + 1L]] <- list(channel = j, label = lab, color = col,
        n_loc = nrow(d), n_trace = dplyr::n_distinct(d$tid),
        src = paste0("file:", fs[j]))   # alignment offsets are keyed on the source
    }
    if (!length(parts)) {
      showNotification("Selected files contained no data.", type = "error"); return(FALSE)
    }
    combined <- dplyr::bind_rows(parts)

  } else {
    d <- store[[cfg$tim_file %||% ""]]
    if (is.null(d) || !nrow(d)) {
      showNotification(sprintf("Load the file to split by time first: %s", cfg$tim_file %||% "(none)"),
                       type = "warning", duration = 12); return(FALSE)
    }
    if (!"tim" %in% names(d)) {
      showNotification("That file has no tim column.", type = "error"); return(FALSE)
    }
    d <- d %>% dplyr::select(-dplyr::any_of(c("mplex_channel", "mplex_label", "mplex_color")))
    attributes(d)[c("mplex_off", "mplex_spec", "mplex_store_name")] <- NULL
    sp  <- .mplex_tim_assign(d$tim, cfg$tim_method, cfg$tim_breaks, cfg$tim_windows)
    nch <- sp$nch
    if (nch < 1L || (identical(sp$method, "breaks") && nch < 2L)) {
      showNotification(
        if (identical(sp$method, "windows"))
          "Enter at least one time window (start-end)."
        else "Enter at least one TIM breakpoint (seconds).",
        type = "warning"); return(FALSE)
    }
    ch   <- sp$ch
    keep <- !is.na(ch)            # drop wash localizations (NA = outside windows)
    if (!any(keep)) {
      showNotification("No localizations fall inside the selected windows.",
                       type = "error"); return(FALSE)
    }
    d  <- d[keep, , drop = FALSE]
    ch <- ch[keep]
    defs <- mplex_default_colors(nch)
    labs <- vapply(seq_len(nch), function(j) pick(cfg$labels, j, paste0("Ch", j)), character(1))
    cols <- vapply(seq_len(nch), function(j) pick(cfg$colors, j, defs[j]), character(1))
    d$mplex_channel <- ch
    d$mplex_label   <- labs[ch]
    d$mplex_color   <- cols[ch]
    d$tid <- suppressWarnings(as.numeric(d$tid)) + ch * MPLEX_TID_OFF
    combined <- d
    spec <- lapply(seq_len(nch), function(j) {
      sub <- d[d$mplex_channel == j, , drop = FALSE]
      list(channel = j, label = labs[j], color = cols[j],
           n_loc = nrow(sub), n_trace = dplyr::n_distinct(sub$tid),
           src = paste0("tim:", j))       # alignment offsets are keyed on the source
    })
    cfg$labels <- labs; cfg$colors <- cols
  }
  if (identical(cfg$mode, "files")) {             # record what was actually used
    cfg$labels <- vapply(spec, `[[`, "", "label"); cfg$colors <- vapply(spec, `[[`, "", "color")
  }

  if (!"z" %in% names(combined)) combined$z <- 0
  combined$.dataset_tag <- paste0(
    sample(c(0:9, letters[1:6]), 4, replace = TRUE), collapse = "")
  # The dataset describes itself: its channels, its store name and (below)
  # the alignment baked into it. Everything downstream reads these.
  attr(combined, "mplex_spec")       <- spec
  attr(combined, "mplex_store_name") <- mplex_name
  attr(combined, "mplex_off")        <- list()      # freshly built: nothing shifted yet

  st <- smlm_file_store(); st[[mplex_name]] <- combined; smlm_file_store(st)
  mplex_cfg_rv(cfg)
  mplex_target_rv(mplex_name)
  mplex_spec_rv(spec)
  # Re-apply the offsets remembered for any source still in the multiplex, so
  # rebuilding (e.g. after renaming a channel) keeps the alignment.
  combined <- .mplex_apply_offsets(update_active = FALSE) %||% combined
  uploaded_df_rv(combined)
  updateSelectInput(session, "active_file", choices = names(st), selected = mplex_name)
  updateSelectInput(session, "color_mode", selected = "multiplex")
  plot_trigger(Sys.time())
  showNotification(sprintf("Built multiplex: %d channel(s), %d localizations.",
                           length(spec), nrow(combined)),
                   type = "message", duration = 6)
  TRUE
}

observeEvent(input$mplex_build, {
  .mplex_build(.mplex_cfg_from_inputs())
}, ignoreInit = TRUE)

observeEvent(input$mplex_clear, {
  st <- smlm_file_store()
  if (!is.null(st[[mplex_name]])) { st[[mplex_name]] <- NULL; smlm_file_store(st) }
  if (identical(mplex_target_rv(), mplex_name)) {  # the wanted offsets are kept for a rebuild
    mplex_target_rv(NULL); mplex_spec_rv(NULL)
  }
  mplex_cfg_rv(NULL)
  rem <- names(smlm_file_store())
  if (length(rem)) {
    updateSelectInput(session, "active_file", choices = rem, selected = rem[1])
    uploaded_df_rv(smlm_file_store()[[rem[1]]])
    .mplex_sync_active(smlm_file_store()[[rem[1]]])
  } else {
    updateSelectInput(session, "active_file", choices = character(0))
    uploaded_df_rv(NULL)
  }
  if (is.null(mplex_spec_rv()) && identical(isolate(input$color_mode), "multiplex"))
    updateSelectInput(session, "color_mode", selected = "default")
})

output$mplex_summary <- renderText({
  sp <- mplex_spec_rv()
  if (is.null(sp) || !length(sp)) return("No multiplex built yet.")
  src <- mplex_target_rv()
  head <- if (!is.null(src) && !identical(src, mplex_name)) paste0("From file: ", src, "\n") else ""
  paste0(head, paste(vapply(sp, function(c) {
    line <- sprintf("Ch%d  %s:  %s loc,  %s traces",
                    c$channel, c$label,
                    format(c$n_loc, big.mark = ","), format(c$n_trace, big.mark = ","))
    pr <- c$prior
    if (length(pr) == 3L && any(abs(pr) > 1e-9))   # alignment already in the file
      line <- paste0(line, sprintf("  (aligned in file: %+.1f, %+.1f, %+.1f nm)",
                                   pr[[1]], pr[[2]], pr[[3]]))
    line
  }, character(1)), collapse = "\n"))
})

# =========================================================================
# MULTIPLEX CHANNEL ALIGNMENT - nudge one channel in x / y (/ z)
# =========================================================================
# GOAL
#   Exchange rounds of a multiplex acquisition can land slightly offset from
#   each other (the sample settles while buffers are exchanged). This block
#   lets the user translate individual channels of the multiplex so the
#   rounds line up again.
#
# WHICH DATASET IS "THE MULTIPLEX"
#   Either the built "Multiplex" entry, or a loaded CSV that was exported
#   from a multiplex (it carries mplex_channel / mplex_label). mplex_target_rv
#   names the store entry; mplex_spec_rv holds its channels. Each such store
#   entry describes itself through attributes:
#     attr "mplex_spec"        its channel list
#     attr "mplex_store_name"  its name in smlm_file_store
#     attr "mplex_off"         the offsets BAKED INTO its coordinates
#   .mplex_sync_active() points target + spec at the active dataset whenever
#   the active dataset is a multiplex.
#
# DESIGN
#   * The loaded source files are never modified. "Clear" + "Build" restores
#     the original coordinates exactly; "Reset" subtracts the shift again
#     (float residue ~1e-13 nm).
#   * Two offset tables, both named lists  source-key -> c(x, y, z) nm:
#       mplex_offsets_rv            what the user WANTS (drives the controls)
#       attr(store entry, "mplex_off")  what is BAKED INTO that entry
#     Applying adds (wanted - baked) to that channel's rows. Working in deltas
#     means repeated nudges never re-read or copy the sources.
#   * Offsets are keyed on the SOURCE ("file:<name>", "tim:<j>",
#     "csv:<file>:<label>"), not on the channel number, so rebuilding the
#     multiplex with a different file selection cannot move the wrong round.
#     Clearing the multiplex keeps the wanted offsets; the next Build
#     re-applies them.
#   * The STORE entry is the truth (batch tools and file switching read it).
#   * The ACTIVE copy (uploaded_df_rv) is deliberately NOT replaced on a nudge.
#     Replacing it re-runs the whole filter pipeline and rebuilds every filter
#     widget at the new data's full range, which resets the user's filters.
#     Instead plottable_data() adds the RESIDUAL shift (store offsets minus the
#     offsets baked into the active copy) right after filtering. The filters
#     are translation-free (metadata, per-trace spread, trace length), so
#     applying the shift after them gives the same point set.
#   * Clustering is per channel on a multiplex, and a per-channel translation
#     keeps every within-channel distance. plottable_data() therefore reuses
#     the last clustering when only the shift changed (see .mplex_clu_memo).
#   * z offsets are only offered for axial data (smlm_z_state() == "3d"), the
#     same planar test every other axis-bearing control uses.
# =========================================================================
mplex_offsets_rv <- reactiveVal(list())   # wanted offsets, by source key
mplex_target_rv  <- reactiveVal(NULL)     # store name of the current multiplex

.mplex_zero <- c(x = 0, y = 0, z = 0)
# Offset for one source key from a table, zero when absent.
.mplex_off_of <- function(tbl, key) {
  v <- tbl[[key]]
  if (is.null(v) || length(v) != 3L || any(!is.finite(v))) .mplex_zero else v
}
# Offsets baked into a store entry (list() when none / not a multiplex).
.mplex_store_off <- function(name) {
  if (is.null(name)) return(list())
  attr(isolate(smlm_file_store())[[name]], "mplex_off") %||% list()
}

# GOAL: point the alignment controls at the active dataset when it is a
# multiplex (built or loaded from CSV). Non-multiplex data leaves the current
# multiplex selection alone, so switching to a single file and back is free.
.mplex_sync_active <- function(df) {
  if (!is.data.frame(df) || !"mplex_label" %in% names(df)) return(invisible(FALSE))
  sp <- attr(df, "mplex_spec"); nm <- attr(df, "mplex_store_name")
  if (is.null(sp) || is.null(nm)) return(invisible(FALSE))
  mplex_target_rv(nm)
  mplex_spec_rv(sp)
  if (!isTRUE(isolate(input$mplex_enable))) updateCheckboxInput(session, "mplex_enable", value = TRUE)
  updateSelectInput(session, "color_mode", selected = "multiplex")
  invisible(TRUE)
}


# GOAL: bring the current multiplex in the STORE to the WANTED offsets.
# ARGS  update_active : replot the active view (FALSE from the builder, which
#                       does that itself).
# RETURNS the updated data.frame (offsets in attr "mplex_off"), or NULL when
#         there is no multiplex.
.mplex_apply_offsets <- function(update_active = TRUE) {
  nm <- isolate(mplex_target_rv())
  st <- isolate(smlm_file_store())
  d  <- if (is.null(nm)) NULL else st[[nm]]
  sp <- attr(d, "mplex_spec")
  if (is.null(sp) || is.null(d) || !nrow(d)) return(NULL)
  want    <- isolate(mplex_offsets_rv())
  baked   <- attr(d, "mplex_off") %||% list()
  has_z   <- "z" %in% names(d) && has_axial_data(d)   # never shift a planar z = 0 column
  changed <- FALSE
  for (c in sp) {
    key <- .mplex_key(c)
    w <- .mplex_off_of(want, key)
    a <- .mplex_off_of(baked, key)
    if (!has_z) w["z"] <- a["z"]                       # z untouched on planar data
    delta <- w - a
    if (all(abs(delta) < 1e-9)) next                   # this channel is already there
    rows <- which(d$mplex_channel == c$channel)        # this channel's localizations
    d$x[rows] <- d$x[rows] + delta[["x"]]
    d$y[rows] <- d$y[rows] + delta[["y"]]
    if (has_z) d$z[rows] <- d$z[rows] + delta[["z"]]
    baked[[key]] <- w                                  # record what is now baked in
    changed <- TRUE
  }
  if (!changed) return(d)
  attr(d, "mplex_off") <- baked                        # the data says how it was shifted
  st[[nm]] <- d
  smlm_file_store(st)
  # The active copy stays as it is (see DESIGN); a replot lets plottable_data()
  # pick up the residual shift. A multiplex sitting in the store behind another
  # active file just keeps its new coordinates.
  if (update_active && identical(isolate(input$active_file), nm))
    plot_trigger(Sys.time())
  d
}

# GOAL: add the shift the active copy does not carry yet.
#   The active copy (uploaded_df_rv) was taken from the store at some point and
#   carries the offsets baked in at that moment in attr "mplex_off". Its store
#   entry may since have moved on. The difference, per channel, is added here.
# ARGS  df : any frame derived from the active copy (keeps mplex_channel).
# RETURNS df with the residual shift applied; unchanged for non-multiplex data.
.mplex_residual_shift <- function(df) {
  if (!is.data.frame(df) || !nrow(df) || !"mplex_channel" %in% names(df)) return(df)
  act <- isolate(uploaded_df_rv())
  sp  <- attr(act, "mplex_spec")
  if (is.null(act) || is.null(sp)) return(df)
  base  <- attr(act, "mplex_off") %||% list()          # offsets inside the active copy
  now   <- .mplex_store_off(attr(act, "mplex_store_name"))
  has_z <- "z" %in% names(df) && has_axial_data(act)
  for (c in sp) {
    key   <- .mplex_key(c)
    delta <- .mplex_off_of(now, key) - .mplex_off_of(base, key)
    if (all(abs(delta) < 1e-9)) next
    rows <- which(df$mplex_channel == c$channel)
    df$x[rows] <- df$x[rows] + delta[["x"]]
    df$y[rows] <- df$y[rows] + delta[["y"]]
    if (has_z) df$z[rows] <- df$z[rows] + delta[["z"]]
  }
  df
}

# Last multiplex clustering, reused while only the channel shift changes.
# Plain environment: a cache, and must not trigger anything when written.
.mplex_clu_memo <- new.env()

# Apply shortly after the last change, so a run of quick clicks costs one
# re-plot rather than one per click.
mplex_offsets_debounced <- debounce(reactive(mplex_offsets_rv()), 300)
observeEvent(mplex_offsets_debounced(), {
  .mplex_apply_offsets()
}, ignoreInit = TRUE)

# TRUE when the active dataset has axial data; drives the z row of the pad.
output$mplex_align_3d <- reactive(identical(smlm_z_state(), "3d"))
outputOptions(output, "mplex_align_3d", suspendWhenHidden = FALSE)

# Channel picker: label shown, source key as the value. Channel 2 is selected
# by default so channel 1 acts as the fixed reference.
observeEvent(mplex_spec_rv(), {
  sp <- mplex_spec_rv()
  if (is.null(sp) || !length(sp)) {
    updateSelectInput(session, "mplex_nudge_ch", choices = character(0)); return()
  }
  keys <- vapply(sp, .mplex_key, character(1))
  labs <- vapply(sp, function(c) sprintf("Ch%d  %s", c$channel, c$label), character(1))
  cur  <- isolate(input$mplex_nudge_ch)
  sel  <- if (!is.null(cur) && cur %in% keys) cur else keys[min(2L, length(keys))]
  updateSelectInput(session, "mplex_nudge_ch", choices = stats::setNames(keys, labs),
                    selected = sel)
  .mplex_show_offset(sel)                              # the boxes follow the new channel list
}, ignoreNULL = FALSE)

# Show the selected channel's wanted offset in the three number boxes.
.mplex_show_offset <- function(key) {
  v <- .mplex_off_of(isolate(mplex_offsets_rv()), key)
  updateNumericInput(session, "mplex_nudge_dx", value = round(v[["x"]], 2))
  updateNumericInput(session, "mplex_nudge_dy", value = round(v[["y"]], 2))
  updateNumericInput(session, "mplex_nudge_dz", value = round(v[["z"]], 2))
}
observeEvent(input$mplex_nudge_ch, {
  req(nzchar(input$mplex_nudge_ch))
  .mplex_show_offset(input$mplex_nudge_ch)
})

# GOAL: add one step to the selected channel's wanted offset (arrow buttons).
.mplex_nudge <- function(axis, sign) {
  key <- input$mplex_nudge_ch
  if (is.null(key) || !nzchar(key)) return()
  step <- suppressWarnings(as.numeric(input$mplex_nudge_step))
  if (!length(step) || !is.finite(step) || step <= 0) return()
  tbl <- mplex_offsets_rv()
  v <- .mplex_off_of(tbl, key)
  v[[axis]] <- v[[axis]] + sign * step
  tbl[[key]] <- v
  mplex_offsets_rv(tbl)
  .mplex_show_offset(key)                              # keep the boxes in step
}
# +y is "up" and +x is "right", matching the 2D XY plot axes.
observeEvent(input$mplex_nudge_left,  .mplex_nudge("x", -1))
observeEvent(input$mplex_nudge_right, .mplex_nudge("x", +1))
observeEvent(input$mplex_nudge_up,    .mplex_nudge("y", +1))
observeEvent(input$mplex_nudge_down,  .mplex_nudge("y", -1))
observeEvent(input$mplex_nudge_zup,   .mplex_nudge("z", +1))
observeEvent(input$mplex_nudge_zdown, .mplex_nudge("z", -1))

# Typed offsets: debounced so a half-typed number is not applied. When the
# boxes were only updated to mirror a button or a channel switch, the values
# already equal the table and this is a no-op.
mplex_typed <- debounce(reactive(list(input$mplex_nudge_dx,
                                      input$mplex_nudge_dy,
                                      input$mplex_nudge_dz)), 600)
observeEvent(mplex_typed(), {
  key <- isolate(input$mplex_nudge_ch)
  if (is.null(key) || !nzchar(key)) return()
  v <- vapply(mplex_typed(), function(e)                # NULL / "" / NA -> NA
    suppressWarnings(as.numeric(e %||% NA)[1]), numeric(1))
  tbl <- isolate(mplex_offsets_rv())
  cur <- .mplex_off_of(tbl, key)
  v[!is.finite(v)] <- cur[!is.finite(v)]              # an emptied box keeps its value
  names(v) <- names(.mplex_zero)
  if (all(abs(v - round(cur, 2)) < 1e-9)) return()     # boxes only mirror the table
  tbl[[key]] <- v
  mplex_offsets_rv(tbl)
}, ignoreInit = TRUE)

observeEvent(input$mplex_nudge_reset, {
  key <- input$mplex_nudge_ch
  if (is.null(key) || !nzchar(key)) return()
  tbl <- mplex_offsets_rv(); tbl[[key]] <- .mplex_zero
  mplex_offsets_rv(tbl)
  .mplex_show_offset(key)
})
observeEvent(input$mplex_nudge_reset_all, {
  mplex_offsets_rv(lapply(mplex_offsets_rv(), function(v) .mplex_zero))
  if (nzchar(input$mplex_nudge_ch %||% "")) .mplex_show_offset(input$mplex_nudge_ch)
})

# Readout of every channel's offset as currently baked into the data.
output$mplex_offsets_txt <- renderText({
  sp <- mplex_spec_rv()
  plot_trigger()                                       # refresh after each apply
  if (is.null(sp) || !length(sp)) return("")
  baked  <- .mplex_store_off(mplex_target_rv())
  show_z <- identical(smlm_z_state(), "3d")
  paste(vapply(sp, function(c) {
    v <- .mplex_off_of(baked, .mplex_key(c))
    if (show_z)
      sprintf("Ch%d  %s:  dx %+.1f  dy %+.1f  dz %+.1f nm", c$channel, c$label,
              v[["x"]], v[["y"]], v[["z"]])
    else
      sprintf("Ch%d  %s:  dx %+.1f  dy %+.1f nm", c$channel, c$label, v[["x"]], v[["y"]])
  }, character(1)), collapse = "\n")
})

# =========================================================================
# SAVE / LOAD MULTIPLEX SETTINGS
# =========================================================================
# GOAL
#   Store the recipe of the built multiplex (which files or time windows,
#   channel names and colours) together with the alignment offsets, so the
#   same aligned multiplex can be rebuilt in a later session in one step.
#
# FILE FORMAT  JSON with settings_type = "minfluxr_multiplex". The marker lets
#   the loaders tell file kinds apart: this loader rejects any other JSON, and
#   the filter-settings loader points the user here when given one of these.
#   Offsets are stored per channel with their source key, in nm at the saved
#   axial calibration (z_scale_applied); dz is rescaled on load if the
#   calibration has changed since.
#
# WHAT IT NEEDS ON LOAD  the same source files, loaded under the same names
#   (files mode) or the file that was split (time mode). Anything missing is
#   named and nothing is built.
# =========================================================================
MPLEX_SETTINGS_TYPE <- "minfluxr_multiplex"

# Save is offered only for a BUILT multiplex; a multiplex loaded from CSV is
# its own record and has no recipe to rebuild.
output$mplex_has_recipe <- reactive(!is.null(mplex_cfg_rv()))
outputOptions(output, "mplex_has_recipe", suspendWhenHidden = FALSE)

output$download_mplex_settings <- downloadHandler(
  filename = function() paste0("multiplex_settings_", format(Sys.time(), "%Y%m%d-%H%M%S"), ".json"),
  content  = function(file) {
    cfg <- mplex_cfg_rv()
    sp  <- mplex_spec_rv()
    want <- mplex_offsets_rv()
    chans <- lapply(sp, function(c) {
      v <- .mplex_off_of(want, .mplex_key(c))
      list(channel = c$channel, label = c$label, color = c$color, source = .mplex_key(c),
           dx_nm = unname(v[["x"]]), dy_nm = unname(v[["y"]]), dz_nm = unname(v[["z"]]))
    })
    s <- list(
      settings_type   = MPLEX_SETTINGS_TYPE,
      version         = 1L,
      saved           = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      mode            = cfg$mode,
      files           = if (identical(cfg$mode, "files")) I(cfg$files) else NULL,
      tim_file        = if (identical(cfg$mode, "tim")) cfg$tim_file else NULL,
      tim_method      = if (identical(cfg$mode, "tim")) cfg$tim_method else NULL,
      tim_breaks      = if (identical(cfg$mode, "tim")) cfg$tim_breaks else NULL,
      tim_windows     = if (identical(cfg$mode, "tim")) cfg$tim_windows else NULL,
      z_scale_applied = z_scale_applied(),
      channels        = chans
    )
    s <- s[!vapply(s, is.null, logical(1))]
    jsonlite::write_json(s, path = file, auto_unbox = TRUE, pretty = TRUE, digits = NA)
  }
)

mplex_settings_reset_key <- reactiveVal(0)          # re-render clears the file box
output$mplex_settings_file_ui <- renderUI({
  mplex_settings_reset_key()
  fileInput("mplex_settings_file", "Load multiplex settings (.json)", accept = ".json")
})

# GOAL: rebuild a saved multiplex, with its names, colours and alignment.
observeEvent(input$mplex_settings_file, {
  req(input$mplex_settings_file$datapath)
  on.exit(mplex_settings_reset_key(mplex_settings_reset_key() + 1), add = TRUE)
  s <- tryCatch(jsonlite::read_json(input$mplex_settings_file$datapath, simplifyVector = FALSE),
                error = function(e) NULL)
  if (!is.list(s) || !identical(s$settings_type, MPLEX_SETTINGS_TYPE)) {
    showNotification("This is not a multiplex settings file.", type = "error", duration = 8)
    return()
  }
  chans <- s$channels %||% list()
  cfg <- list(mode = s$mode %||% "files",
              labels = vapply(chans, function(c) c$label %||% "", ""),
              colors = vapply(chans, function(c) c$color %||% "", ""))
  if (identical(cfg$mode, "files")) {
    cfg$files <- as.character(unlist(s$files))
  } else {
    cfg$tim_file    <- s$tim_file %||% ""
    cfg$tim_method  <- s$tim_method %||% "breaks"
    cfg$tim_breaks  <- s$tim_breaks %||% ""
    cfg$tim_windows <- s$tim_windows %||% ""
  }

  # Offsets first, so the build applies them. dz follows the calibration.
  sz <- suppressWarnings(as.numeric(s$z_scale_applied %||% NA)[1])
  zr <- if (is.finite(sz) && sz > 0) z_scale_applied() / sz else 1
  tbl <- mplex_offsets_rv()
  for (c in chans) {
    if (is.null(c$source)) next
    tbl[[c$source]] <- c(x = as.numeric(c$dx_nm %||% 0), y = as.numeric(c$dy_nm %||% 0),
                         z = as.numeric(c$dz_nm %||% 0) * zr)
  }
  mplex_offsets_rv(tbl)

  if (!isTRUE(.mplex_build(cfg))) return()

  # Mirror the recipe in the controls, so a later Build reproduces it.
  updateCheckboxInput(session, "mplex_enable", value = TRUE)
  updateRadioButtons(session, "mplex_mode", selected = cfg$mode)
  if (identical(cfg$mode, "files")) {
    updateSelectInput(session, "mplex_files", selected = cfg$files)
  } else {
    updateRadioButtons(session, "mplex_tim_method", selected = cfg$tim_method)
    updateTextInput(session, "mplex_tim_breaks",  value = cfg$tim_breaks)
    updateTextInput(session, "mplex_tim_windows", value = cfg$tim_windows)
  }
  # Per-channel label / colour boxes are rendered from the controls above, so
  # they exist only after the next flush.
  pre <- if (identical(cfg$mode, "files")) c("mplex_flab_", "mplex_fcol_") else c("mplex_tlab_", "mplex_tcol_")
  labs <- cfg$labels; cols <- cfg$colors
  later::later(function() {
    for (j in seq_along(labs)) {
      if (nzchar(labs[j])) updateTextInput(session, paste0(pre[1], j), value = labs[j])
      if (nzchar(cols[j])) colourpicker::updateColourInput(session, paste0(pre[2], j), value = cols[j])
    }
  }, delay = 0.6)
  showNotification(paste0("Multiplex settings loaded: ", input$mplex_settings_file$name),
                   type = "message", duration = 5)
})

# Split the active multiplex dataset into a named list (label → channel df) so
# the per-file Analysis loops produce one result per channel. NULL when the
# active dataset is not a multiplex. Used by the *_by_channel toggles.
mplex_split_active <- reactive({
  df <- uploaded_df_rv()                           # keeps the reactive dependency
  plot_trigger()                                   # refresh after each nudge is applied
  df <- .mplex_residual_shift(df)                  # include any alignment nudge
  if (is.null(df) || !"mplex_label" %in% names(df) || !nrow(df)) return(NULL)
  ord  <- order(df$mplex_channel)
  labs <- unique(df$mplex_label[ord])
  stats::setNames(
    lapply(labs, function(l) df[df$mplex_label == l, , drop = FALSE]), labs)
})

# ── Split one prepared file into its channels, for the spatial statistics ──


# TRUE when the active dataset is a multiplex (drives the *_by_channel UIs).
output$mplex_active <- reactive({
  df <- uploaded_df_rv()
  !is.null(df) && "mplex_label" %in% names(df)
})
outputOptions(output, "mplex_active", suspendWhenHidden = FALSE)

# How many channels the loaded data currently resolves to. Drives every
# conditionalPanel that only makes sense once there is more than one - a scope
# control on single-channel data is a control with one legal answer.
# Deliberately counts SOURCES rather than calling resolve_channels(): this runs
# before Plot and must not depend on the clustering chain.
output$n_channels <- reactive({
  df <- uploaded_df_rv()
  if (is.null(df)) return(1L)
  n_m <- if ("mplex_channel" %in% names(df))
    length(unique(df$mplex_channel[!is.na(df$mplex_channel)])) else 1L
  n_c <- if (isTRUE(input$enable_two_color)) 2L else 1L
  max(1L, n_m) * n_c
})
outputOptions(output, "n_channels", suspendWhenHidden = FALSE)

# Scope picker, built from the channels that actually exist. There is no
# "pooled" entry: pooling across channels is exactly what the channel column
# exists to prevent, so offering it would be offering a wrong answer.
output$channel_scope_ui <- renderUI({
  df <- tryCatch(filtered_df(), error = function(e) NULL)
  ct <- NULL
  if (!is.null(df) && nrow(df))
    ct <- tryCatch(attr(attach_channel_colors(resolve_channels(df, as_smlm_settings(input)), smlm_channel_colors(input)),
                        "channel_table"), error = function(e) NULL)
  if (is.null(ct) || nrow(ct) < 2L) return(NULL)
  ch <- stats::setNames(paste0("only:", ct$channel),
                        paste0("Only ", ct$label))
  selectInput("cluster_group_scope", "Channel scope",
              choices  = c("Within each channel (recommended)" = "within", ch),
              selected = input$cluster_group_scope %||% SMLM_DEFAULTS$cluster_group_scope)
})

# Channel pair (A, B) for the pairwise cross-channel analyses, or NULL. When
# set, the cross tools map channel A→group 1 / B→group 2 and reuse the entire
# existing 2-color machinery (which honours a pre-assigned analysis_group).
mplex_cross_pair <- reactive({
  df <- uploaded_df_rv()
  if (is.null(df) || !"mplex_label" %in% names(df)) return(NULL)
  a <- input$mplex_cross_a; b <- input$mplex_cross_b
  if (is.null(a) || is.null(b) || !nzchar(a) || !nzchar(b) || identical(a, b))
    return(NULL)
  c(a, b)
})

# Keep the channel selectors (cross-channel A/B + Batch-NN X/Y) populated.
observeEvent(mplex_spec_rv(), {
  sp <- mplex_spec_rv()
  labs <- if (is.null(sp)) character(0) else vapply(sp, function(c) c$label, character(1))
  sel1 <- labs[1] %||% NULL
  sel2 <- if (length(labs) >= 2) labs[2] else NULL
  updateSelectInput(session, "mplex_cross_a", choices = labs, selected = sel1)
  updateSelectInput(session, "mplex_cross_b", choices = labs, selected = sel2)
  updateSelectInput(session, "nn_batch_chx",  choices = labs, selected = sel1)
  updateSelectInput(session, "nn_batch_chy",  choices = labs, selected = sel2)
}, ignoreNULL = FALSE)
