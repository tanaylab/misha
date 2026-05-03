# R/genome-build-installers.R
# Per-set interval installers (genes/rmsk/cgi/cytoband) and the shared .save_intervals() helper.

# Save a data.frame as a misha intervals set under the active groot.
# - Filters to ALLGENOME contigs (counts dropped rows in a single message).
# - Errors on existing set unless overwrite = TRUE.
# - Skips empty frames with a message (no error).
.save_intervals <- function(name, df, overwrite = FALSE, verbose = TRUE) {
    if (!nrow(df)) {
        if (verbose) message(sprintf("  %s: 0 rows after filtering, skipping", name))
        return(invisible(NULL))
    }
    if (gintervals.exists(name)) {
        if (!overwrite) {
            stop(
                sprintf(
                    "Intervals set '%s' already exists. Use overwrite = TRUE to replace.",
                    name
                ),
                call. = FALSE
            )
        }
        gintervals.rm(name, force = TRUE)
    }
    chroms <- as.character(get("ALLGENOME", envir = .misha)[[1]]$chrom)
    n_in <- nrow(df)
    df <- df[df$chrom %in% chroms, , drop = FALSE]
    n_dropped <- n_in - nrow(df)
    if (n_dropped > 0L && verbose) {
        message(sprintf(
            "  %s: dropped %d/%d rows on contigs not in ALLGENOME",
            name, n_dropped, n_in
        ))
    }
    if (!nrow(df)) {
        if (verbose) message(sprintf("  %s: 0 rows remained after filtering, skipping", name))
        return(invisible(NULL))
    }
    df$chrom <- factor(df$chrom, levels = chroms)
    gintervals.save(name, df)
    if (verbose) message(sprintf("  %s: saved %d intervals", name, nrow(df)))
    invisible(NULL)
}

# DEPRECATED shim: kept until existing backends are migrated to the new installers.
# Remove when backends are rewired in Phase D/E (Task 22).
.save_post_build_intervals <- function(intervals_set_name, df, verbose = TRUE) {
    .save_intervals(intervals_set_name, df, overwrite = FALSE, verbose = verbose)
}
