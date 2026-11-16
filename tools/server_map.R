## ─────────────────────────────────────────────────────────────────────────────
## tools/server_map.R  -  structural map of a Shiny server closure
## ─────────────────────────────────────────────────────────────────────────────
##
## GOAL
##   The server files are single closures of many thousand lines. File order is
##   not execution order, so reading top to bottom does not show what depends on
##   what. This script parses a server function and lists every top-level
##   statement in its body (reactive, reactiveVal, function, output, observer)
##   with its line span, the closure-level names it reads, the input ids it
##   reads, and the section banner it sits under. It is the source of the
##   dependency tables in docs/architecture_review.md and the tool used to plan
##   and check the split into panel files.
##
## USE
##   Rscript tools/server_map.R                 # maps the three server closures
##   Rscript tools/server_map.R out_dir         # writes <out_dir>/<name>_map.csv
##   From R:  source("tools/server_map.R"); m <- map_server_closure(file, fn)
##
## It only parses the code (no package is loaded, nothing is run).
## ─────────────────────────────────────────────────────────────────────────────

# GOAL: the name of the function a call invokes, as a string ("" when the call
# head is itself a call, e.g. f(x)(y)). pkg::fn is reduced to fn.
.call_head <- function(e) {
  if (!is.call(e)) return("")
  h <- e[[1]]
  if (is.symbol(h)) return(as.character(h))
  if (is.call(h) && is.symbol(h[[1]]) && as.character(h[[1]]) %in% c("::", ":::") && is.symbol(h[[3]]))
    return(as.character(h[[3]]))
  ""
}

# GOAL: the positions of a call (or pairlist) that can be walked: everything but
# NULL and the empty argument of calls like x[, 1], which errors when touched.
.walkable <- function(x) {
  lst <- as.list(x)
  which(!vapply(seq_along(lst), function(i)
    is.null(lst[[i]]) || identical(lst[[i]], quote(expr = )), logical(1)))
}

# GOAL: every input id an expression reads, from both spellings: input$id and
# input[["id"]]. Ids built at run time (input[[paste0(...)]]) come back as "<dynamic>".
.input_ids <- function(e) {
  out <- character(0)
  walk <- function(x) {
    if (is.call(x)) {
      h <- .call_head(x)
      if (h %in% c("$", "[[") && length(x) >= 3 && identical(x[[2]], as.name("input"))) {
        id <- x[[3]]
        out <<- c(out, if (is.symbol(id) || is.character(id)) as.character(id) else "<dynamic>")
      }
      for (i in .walkable(x)) walk(x[[i]])
    } else if (is.pairlist(x) || is.expression(x)) {
      for (i in .walkable(x)) walk(x[[i]])
    }
  }
  walk(e)
  unique(out)
}

# GOAL: every output id an expression assigns (output$id <- ..., output[["id"]] <- ...),
# at any depth, so outputs created inside lapply()/local() are found too.
.output_ids <- function(e) {
  out <- character(0)
  walk <- function(x) {
    if (is.call(x)) {
      h <- .call_head(x)
      if (h %in% c("<-", "=", "<<-") && is.call(x[[2]]) &&
          .call_head(x[[2]]) %in% c("$", "[[") && identical(x[[2]][[2]], as.name("output"))) {
        id <- x[[2]][[3]]
        out <<- c(out, if (is.symbol(id) || is.character(id)) as.character(id) else "<dynamic>")
      }
      for (i in .walkable(x)) walk(x[[i]])
    }
  }
  walk(e)
  unique(out)
}

# GOAL: classify one top-level statement of a server body.
# RETURNS list(kind, name, trigger). kind is one of: reactive, eventReactive,
# reactiveVal, reactiveValues, debounce, function, value, output, observeEvent,
# observe, call.
.classify <- function(e) {
  h <- .call_head(e)
  if (h %in% c("<-", "=", "<<-")) {
    lhs <- e[[2]]; rhs <- e[[3]]
    if (is.call(lhs) && .call_head(lhs) %in% c("$", "[[") && identical(lhs[[2]], as.name("output")))
      return(list(kind = "output", name = paste0("output$", paste(as.character(lhs[[3]]), collapse = "")),
                  trigger = .call_head(rhs)))
    nm <- paste(deparse(lhs, width.cutoff = 200L), collapse = "")
    rh <- .call_head(rhs)
    # x <- reactive(...) %>% debounce(500) and x <- debounce(reactive(...), 500)
    if (rh %in% c("%>%", "|>") ) rh <- .call_head(rhs[[3]])
    kind <- if (rh %in% c("reactive", "eventReactive", "reactiveVal", "reactiveValues")) rh
            else if (rh %in% c("debounce", "throttle")) "debounce"
            else if (rh == "function") "function"
            else if (rh %in% c("bindEvent", "bindCache")) "reactive"
            else "value"
    trig <- if (kind == "eventReactive" && length(rhs) >= 2)
      paste(deparse(rhs[[2]], width.cutoff = 200L), collapse = "") else ""
    return(list(kind = kind, name = nm, trigger = trig))
  }
  if (h %in% c("observeEvent", "observe")) {
    trig <- if (h == "observeEvent" && length(e) >= 2)
      paste(deparse(e[[2]], width.cutoff = 200L), collapse = "") else ""
    return(list(kind = h, name = "", trigger = trig))
  }
  list(kind = "call", name = h, trigger = "")
}

# GOAL: the structural map of one server closure.
# ARGS  file    : path to the R file
#       fn_name : name of the server function defined at the top level of `file`
#                 (NULL = the file is a fragment whose top-level statements ARE the
#                 closure body, as in the server/ panel files)
#       known   : extra closure-level names to count as dependencies (used when a
#                 closure is spread over several fragment files)
# RETURNS data.frame, one row per top-level statement: line span, kind, name,
#   trigger, section banner, deps (closure names read), inputs (input ids read),
#   outputs (output ids assigned).
map_server_closure <- function(file, fn_name = NULL, known = character(0)) {
  src   <- readLines(file, warn = FALSE, encoding = "UTF-8")
  exprs <- parse(file, keep.source = TRUE, encoding = "UTF-8")
  if (is.null(fn_name)) {
    stmts <- as.list(exprs)
    refs  <- attr(exprs, "srcref")
  } else {
    hit <- which(vapply(exprs, function(e) .call_head(e) %in% c("<-", "=") &&
                          identical(e[[2]], as.name(fn_name)), logical(1)))
    if (length(hit) != 1L) stop("expected exactly one top-level definition of ", fn_name)
    body  <- exprs[[hit]][[3]][[3]]                 # function(...) { body }
    stmts <- as.list(body)[-1]
    refs  <- attr(body, "srcref")[-1]
  }
  cls   <- lapply(stmts, .classify)
  names_defined <- unique(c(known, vapply(cls, function(k) if (k$kind != "output") k$name else "", "")))
  names_defined <- names_defined[nzchar(names_defined) & grepl("^[.A-Za-z][.A-Za-z0-9_]*$", names_defined)]

  # Section banners: comment lines that carry a heading mark.
  banner_idx <- grep("^\\s*#+\\s*(──|==+ |--+ |[0-9]+ ?\\)|══)\\s*\\S", src)
  banner_txt <- trimws(gsub("[#─═=\\-]+", " ", src[banner_idx]))
  keep <- nzchar(banner_txt); banner_idx <- banner_idx[keep]; banner_txt <- banner_txt[keep]

  rows <- lapply(seq_along(stmts), function(i) {
    e  <- stmts[[i]]; k <- cls[[i]]
    l1 <- refs[[i]][1]; l2 <- refs[[i]][3]
    syms <- unique(all.names(e, functions = TRUE, unique = TRUE))
    deps <- setdiff(intersect(syms, names_defined), k$name)
    b    <- banner_idx[banner_idx <= l1]
    data.frame(
      file = basename(file), line_start = l1, line_end = l2, n_lines = l2 - l1 + 1L,
      kind = k$kind, name = k$name, trigger = substr(k$trigger, 1, 120),
      section = if (length(b)) substr(banner_txt[which(banner_idx == max(b))], 1, 70) else "",
      deps = paste(sort(deps), collapse = ";"),
      inputs = paste(sort(.input_ids(e)), collapse = ";"),
      outputs = paste(sort(.output_ids(e)), collapse = ";"),
      stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  attr(out, "defined") <- names_defined
  out
}

# GOAL: for each closure-level name, how many statements read it (fan-in from
# the reader's side = how widely shared the object is).
closure_fan_in <- function(map) {
  deps <- strsplit(map$deps, ";", fixed = TRUE)
  tab  <- sort(table(unlist(deps)), decreasing = TRUE)
  data.frame(name = names(tab), n_readers = as.integer(tab), stringsAsFactors = FALSE)
}

# GOAL: names assigned more than once at closure level. The later definition
# silently replaces the earlier one, which then never runs (this has happened).
closure_duplicates <- function(map) {
  nm <- map$name[nzchar(map$name) & map$kind != "call"]
  unique(nm[duplicated(nm)])
}

if (sys.nframe() == 0L) {
  args    <- commandArgs(trailingOnly = TRUE)
  out_dir <- if (length(args)) args[1] else NULL
  targets <- list(c("server_smlm.R", "smlm_server"), c("server_tracking.R", "tracking_server"),
                  c("sim_server.R", "sim_server"))
  for (t in targets) {
    if (!file.exists(t[1])) next
    m <- map_server_closure(t[1], t[2])
    cat(sprintf("\n%s :: %s  -  %d top-level statements\n", t[1], t[2], nrow(m)))
    print(table(m$kind))
    d <- closure_duplicates(m)
    if (length(d)) cat("DUPLICATE closure-level names:", paste(d, collapse = ", "), "\n")
    if (!is.null(out_dir)) {
      dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
      utils::write.csv(m, file.path(out_dir, paste0(t[2], "_map.csv")), row.names = FALSE)
      utils::write.csv(closure_fan_in(m), file.path(out_dir, paste0(t[2], "_fan_in.csv")), row.names = FALSE)
    }
  }
}
