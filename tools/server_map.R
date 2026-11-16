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
##   Rscript tools/server_map.R --headers       # refresh the DEFINES/READS block of every part
##   From R:  source("tools/server_map.R"); m <- server_statement_map("smlm")
##   Run from the repo root.
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
  names_defined <- unique(c(known, vapply(cls, function(k) if (k$kind %in% c("output", "call")) "" else k$name, "")))
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

# The three server closures: the file that defines the function, and the folder
# that holds its part files once the closure is split (see server/README.md).
SERVER_CLOSURES <- list(
  smlm     = list(file = "server_smlm.R",     fn = "smlm_server",     parts = file.path("server", "smlm")),
  tracking = list(file = "server_tracking.R", fn = "tracking_server", parts = file.path("server", "tracking")),
  sim      = list(file = "sim_server.R",      fn = "sim_server",      parts = file.path("server", "sim"))
)

# GOAL: the statement map of one tab's server closure, wherever its code lives:
# in the single server file, or spread over ordered part files. With part files,
# a name defined in any part counts as a dependency in every part, as it does at
# run time (all parts share one environment).
# tab : "smlm", "tracking" or "sim". Run from the repo root.
server_statement_map <- function(tab) {
  s <- SERVER_CLOSURES[[tab]]
  if (is.null(s)) stop("unknown tab: ", tab)
  parts <- if (dir.exists(s$parts)) sort(list.files(s$parts, "[.]R$", full.names = TRUE)) else character(0)
  if (!length(parts)) return(map_server_closure(s$file, s$fn))
  first <- lapply(parts, map_server_closure, fn_name = NULL)
  known <- unique(unlist(lapply(first, attr, "defined")))
  out   <- do.call(rbind, lapply(parts, map_server_closure, fn_name = NULL, known = known))
  attr(out, "defined") <- known
  out
}

## ── Part interfaces ──────────────────────────────────────────────────────────
## Once a closure is split into part files, each part's header states what it
## DEFINES for other parts and what it READS from them. The functions below
## compute those lists from the code, write them into the headers, and let a
## test confirm the headers are current.

# GOAL: the closure-level names one statement reads. Scope-aware (codetools): a
# local variable that happens to share a closure object's name is not a read.
# The statement is wrapped as the body of a function so its own assignment
# target counts as local.
.stmt_reads <- function(e, names) {
  g <- tryCatch({
    f <- eval(call("function", NULL, e), baseenv())
    codetools::findGlobals(f, merge = TRUE)
  }, error = function(err) all.names(e, unique = TRUE))
  intersect(g, names)
}

# GOAL: for every part file of a tab, the names it defines that another part
# reads, and the names it reads from other parts (with the part they come from).
# RETURNS a named list per part file: list(defines = chr, reads = named list
#   source part -> chr, n_lines, n_statements).
server_part_interfaces <- function(tab) {
  s <- SERVER_CLOSURES[[tab]]
  parts <- sort(list.files(s$parts, "^[0-9]{2}_.*[.]R$", full.names = TRUE))
  parsed <- lapply(parts, parse, keep.source = FALSE, encoding = "UTF-8")
  defs <- lapply(parsed, function(p) {
    nm <- vapply(as.list(p), function(e) { k <- .classify(e); if (k$kind %in% c("output", "call")) "" else k$name }, "")
    unique(nm[grepl("^[.A-Za-z][.A-Za-z0-9_]*$", nm)])
  })
  owner <- stats::setNames(rep(basename(parts), lengths(defs)), unlist(defs))
  all_names <- names(owner)
  reads <- lapply(parsed, function(p) unique(unlist(lapply(as.list(p), .stmt_reads, names = all_names))))
  out <- lapply(seq_along(parts), function(i) {
    me <- basename(parts[i])
    ext <- sort(setdiff(reads[[i]], defs[[i]]))
    used_elsewhere <- unique(unlist(reads[-i]))
    list(defines = sort(intersect(defs[[i]], used_elsewhere)),
         reads = split(ext, owner[ext]),
         n_lines = length(readLines(parts[i], warn = FALSE)), n_statements = length(parsed[[i]]))
  })
  stats::setNames(out, basename(parts))
}

.PART_IFACE_OPEN  <- "## INTERFACE (generated: Rscript tools/server_map.R --headers; do not edit by hand)"
.PART_IFACE_CLOSE <- "## END INTERFACE"

# GOAL: the generated interface block of one part, as comment lines.
.part_interface_lines <- function(iface) {
  wrap <- function(prefix, x) {
    if (!length(x)) return(paste0(prefix, "(nothing)"))
    w <- strwrap(paste(x, collapse = ", "), width = 78 - nchar(prefix))
    paste0(c(prefix, rep(paste0("##", strrep(" ", nchar(prefix) - 2L)), length(w) - 1L)), w)
  }
  rd <- if (length(iface$reads)) unlist(lapply(names(iface$reads), function(src)
          wrap(sprintf("##     %s: ", sub("[.]R$", "", src)), iface$reads[[src]]))) else "##     (nothing)"
  c(.PART_IFACE_OPEN,
    "##   DEFINES for other parts:", sub("^##   ", "##     ", wrap("##   ", iface$defines)),
    "##   READS from other parts:", rd,
    .PART_IFACE_CLOSE)
}

# GOAL: write (or refresh) the interface block in every part header of a tab.
# A part with no block yet gets a full header built from `titles` (file name ->
# c(title, purpose)); an existing header keeps its hand-written lines.
write_part_headers <- function(tab, titles = list()) {
  s <- SERVER_CLOSURES[[tab]]
  ifs <- server_part_interfaces(tab)
  rule <- paste0("## ", strrep("─", 77))
  for (f in names(ifs)) {
    path <- file.path(s$parts, f)
    x <- readLines(path, warn = FALSE, encoding = "UTF-8")
    blk <- .part_interface_lines(ifs[[f]])
    a <- match(.PART_IFACE_OPEN, x); b <- match(.PART_IFACE_CLOSE, x)
    if (!is.na(a) && !is.na(b)) {
      x <- c(x[seq_len(a - 1L)], blk, x[-seq_len(b)])
    } else {
      t <- titles[[f]] %||% c(sub("[.]R$", "", f), "")
      x <- c(rule, sprintf("## %s/%s  -  %s", gsub("\\\\", "/", s$parts), f, t[1]), rule,
             sprintf("## Part of %s(). Not a standalone script: its statements are evaluated in", s$fn),
             "## the server function's environment, in file-name order (server/parts.R).",
             "##", paste0("## ", strwrap(t[2], 76)), "##", blk, rule, "", x)
    }
    con <- file(path, open = "wb"); writeLines(enc2utf8(x), con, sep = "\n", useBytes = TRUE); close(con)
  }
  invisible(ifs)
}

`%||%` <- function(a, b) if (is.null(a)) b else a

if (sys.nframe() == 0L) {
  args    <- commandArgs(trailingOnly = TRUE)
  if ("--headers" %in% args) {                 # refresh the interface block of every part
    for (tab in names(SERVER_CLOSURES)) if (dir.exists(SERVER_CLOSURES[[tab]]$parts)) write_part_headers(tab)
    cat("part headers refreshed
"); quit(save = "no")
  }
  out_dir <- if (length(args)) args[1] else NULL
  for (tab in names(SERVER_CLOSURES)) {
    s <- SERVER_CLOSURES[[tab]]
    if (!file.exists(s$file)) next
    m <- server_statement_map(tab)
    cat(sprintf("
%s  -  %d top-level statements in %d file(s)
", s$fn, nrow(m), length(unique(m$file))))
    print(table(m$kind))
    d <- closure_duplicates(m)
    if (length(d)) cat("DUPLICATE closure-level names:", paste(d, collapse = ", "), "
")
    if (!is.null(out_dir)) {
      dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
      utils::write.csv(m, file.path(out_dir, paste0(s$fn, "_map.csv")), row.names = FALSE)
      utils::write.csv(closure_fan_in(m), file.path(out_dir, paste0(s$fn, "_fan_in.csv")), row.names = FALSE)
    }
  }
}
