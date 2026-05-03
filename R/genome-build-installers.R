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

# Save the combined `<prefix>rmsk` and a `<prefix>rmsk_<class>` per unique class.
# Class normalization: lowercase; trailing '?' -> '_qmark'.
# df must already be chrom-translated and contain columns:
#   chrom start end strand name class family
.install_rmsk_set <- function(df, prefix = "", overwrite = FALSE, verbose = TRUE) {
    cols <- c("chrom", "start", "end", "strand", "name", "class", "family")
    if (!all(cols %in% names(df))) {
        stop(sprintf(
            ".install_rmsk_set: df missing columns: %s",
            paste(setdiff(cols, names(df)), collapse = ", ")
        ), call. = FALSE)
    }
    df <- df[, cols, drop = FALSE]

    # Bump big-set option for the duration; restore on exit.
    old_big <- getOption("gbig.intervals.size")
    options(gbig.intervals.size = max(if (is.null(old_big)) 0 else old_big, 1e9))
    on.exit(options(gbig.intervals.size = old_big), add = TRUE)

    .save_intervals(paste0(prefix, "rmsk"), df,
        overwrite = overwrite, verbose = verbose
    )

    classes <- unique(df$class)
    classes <- classes[!is.na(classes) & nzchar(classes)]
    for (cls in classes) {
        sub <- df[!is.na(df$class) & df$class == cls, , drop = FALSE]
        suffix <- gsub("\\?$", "_qmark", tolower(cls))
        .save_intervals(sprintf("%srmsk_%s", prefix, suffix), sub,
            overwrite = overwrite, verbose = verbose
        )
    }
    invisible(NULL)
}
