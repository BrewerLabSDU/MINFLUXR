## ─────────────────────────────────────────────────────────────────────────────
## server/tracking/26_sigma_used.R  -  which sigma_loc each panel used
## ─────────────────────────────────────────────────────────────────────────────
## Part of tracking_server(). Not a standalone script: its statements are evaluated in
## the server function's environment, in file-name order (server/parts.R).
##
## One line under each panel whose result depends on the localization precision,
## saying which sigma_loc it used when it ran and where that value came from. The
## record is taken in each panel's Run handler (.sigma_stamp, 01_state_sigma.R)
## and worded by sigma_used_text() (R/trk_stats.R).
##
## INTERFACE (generated: Rscript tools/server_map.R --headers; do not edit by hand)
##   DEFINES for other parts:
##     (nothing)
##   READS from other parts:
##     (nothing)
## END INTERFACE
## ─────────────────────────────────────────────────────────────────────────────

# Panel -> the result whose "sigma_info" it shows. The output id is
# <panel>_sigma_used; ui_tracking.R places it under the panel's summary.
.SIGMA_USED_PANELS <- c(trk_jdd = "jdd_bundle", trk_dmap = "dmap_bundle", trk_cpd = "cpd_segments_rv",
                        trk_hmm = "hmm_bundle", trk_conf_ci = "conf_ci_bundle", trk_disp = "trk_disp_bundle",
                        trk_dir = "dir_bundle", trk_speed = "trk_export_summary")
for (.p in names(.SIGMA_USED_PANELS)) local({
  res_name <- .SIGMA_USED_PANELS[[.p]]
  output[[paste0(.p, "_sigma_used")]] <- renderUI({
    r <- tryCatch(get(res_name)(), error = function(e) NULL)
    if (is.null(r)) return(NULL)
    info <- attr(r, "sigma_info")
    txt <- sigma_used_text(info)
    if (is.null(txt)) return(NULL)
    warn <- !isTRUE(info$used) || !is.finite(info$sigma) || info$sigma <= 0 ||
            (identical(info$source, "photon budget") && !isTRUE(info$a_calibrated))
    tags$div(style = paste0("font-size:0.82em; margin:4px 0 6px; color:", if (warn) "#8a6d00;" else "#666;"), txt)
  })
})
rm(.p)
