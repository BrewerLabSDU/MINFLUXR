# server/

The bodies of the three server functions, split into part files.

| Folder | Server function | Assembled in |
|---|---|---|
| `smlm/` | `smlm_server()` | `server_smlm.R` |
| `tracking/` | `tracking_server()` | `server_tracking.R` |
| `sim/` | `sim_server()` | `sim_server.R` |

`parts.R` holds the two functions that load and run the parts.

## What a part file is

A part is a run of statements that used to sit in the server function's body:
reactives, observers, outputs and the functions that need them. It is not a
standalone script and it does not define a function. When a session starts, the
server function evaluates every part, in file-name order, in its own
environment. So:

- All parts of one tab share one environment. A part can use any name another
  part of the same tab defines, including one defined in a later file, as long as
  it is used inside a reactive, observer or function (those run later).
- The number in the file name is the order of evaluation. It matters for the few
  statements that run at once (a `debounce()` of an earlier reactive, an
  `outputOptions()` call, plain values) and for the order in which observers are
  created. Do not renumber parts without checking that.
- Input and output ids are global, as before. Nothing is namespaced, so settings
  files, `conditionalPanel` conditions and plugins see the same ids.

Computation that needs no session state belongs in `R/`, not here.

## The header of a part

Each part starts with its title, a sentence on what it is for, and a generated
block:

```
## INTERFACE (generated: Rscript tools/server_map.R --headers; do not edit by hand)
##   DEFINES for other parts:   names this part provides to the rest of the tab
##   READS from other parts:    names it uses, grouped by the part that defines them
## END INTERFACE
```

That block is the part's interface. `tests/testthat/test-server-parts.R` fails
when it no longer matches the code. After editing a part:

```
Rscript tools/server_map.R --headers
```

and read the diff. A new line under READS is a new dependency between panels.

## Finding things

`Rscript tools/server_map.R out/` writes one CSV per tab listing every reactive,
observer, output and function with its file, line span, the names it reads and
the input ids it reads.

## Adding a panel

1. Put the computation in `R/` as a plain function and test it there.
2. Add a part file with the next free number (or append to the part of a related
   panel). Give it the header lines; leave the interface block to the tool.
3. Check the name of every new reactive is unused in the tab
   (`test-server-parts.R` checks this too: a second definition silently replaces
   the first).
4. Add the panel to the golden test of its tab.

## Why not Shiny modules with `NS()`

See `docs/architecture_review.md`, section 4.1. In short: the ids are a public
interface here (settings files, JavaScript conditions, plugins), and a renamed id
cannot be verified without a browser test. The interface block gives each part a
declared, checked dependency list without renaming anything; turning a part into
a function that takes those dependencies as arguments is the next step and can
be done one panel at a time.
