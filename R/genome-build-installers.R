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

.install_cgi_set <- function(df, prefix = "", overwrite = FALSE, verbose = TRUE) {
    .save_intervals(paste0(prefix, "cgi"), df, overwrite = overwrite, verbose = verbose)
}

.install_cytoband_set <- function(df, prefix = "", overwrite = FALSE, verbose = TRUE) {
    .save_intervals(paste0(prefix, "cytoband"), df, overwrite = overwrite, verbose = verbose)
}

# Install genes-derived interval sets onto the active groot.
# `asset` is one entry from a fetcher's output:
#   list(file = "<path>", format = "gtf"|"gff3"|"genepred",
#        gtf_source = "ncbiRefSeq" or NULL,
#        translate = function(rows, chrom_col) -> rows
#                    (chrom translator closure; identity if no chromAlias))
# `gene_sets` is a named char vec mapping role -> on-disk name; NA skips a role.
# Roles must include: tss, exons, utr3, utr5.
.install_genes_set <- function(asset, prefix = "", gene_sets = NULL,
                               overwrite = FALSE, verbose = TRUE) {
    stopifnot(!is.null(gene_sets), all(c("tss", "exons", "utr3", "utr5") %in% names(gene_sets)))

    # 1. Convert to genePred if needed.
    workdir <- tempfile("misha_install_genes_")
    dir.create(workdir, recursive = TRUE)
    on.exit(unlink(workdir, recursive = TRUE), add = TRUE)

    src <- asset$file
    if (asset$format == "gtf") {
        converter <- .gtf_to_genepred_resolve_or_install()
        gp_raw <- file.path(workdir, "raw.genePred")
        ret <- system2(converter,
            args = c("-genePredExt", shQuote(src), shQuote(gp_raw)),
            stdout = if (verbose) "" else FALSE,
            stderr = if (verbose) "" else FALSE
        )
        if (ret != 0L) stop(sprintf("gtfToGenePred failed (exit %d) on %s", ret, src), call. = FALSE)
        src <- gp_raw
    } else if (asset$format == "gff3") {
        converter <- .gff3_to_genepred_resolve_or_install()
        gp_raw <- file.path(workdir, "raw.genePred")
        ret <- system2(converter,
            args = c(shQuote(src), shQuote(gp_raw)),
            stdout = if (verbose) "" else FALSE,
            stderr = if (verbose) "" else FALSE
        )
        if (ret != 0L) stop(sprintf("gff3ToGenePred failed (exit %d) on %s", ret, src), call. = FALSE)
        src <- gp_raw
    }
    # else format == "genepred": use src directly.

    # 2. Trim 15-col extended -> 12-col classic if needed (existing helper).
    gp_trim <- file.path(workdir, "trim.genePred")
    .normalize_ucsc_genepred(src, gp_trim)

    # 3. Translate chrom column (col 2) via asset$translate, if provided.
    if (!is.null(asset$translate)) {
        rows <- utils::read.table(gp_trim,
            sep = "\t", header = FALSE,
            stringsAsFactors = FALSE, quote = "", comment.char = ""
        )
        rows <- asset$translate(rows, 2L)
        gp_xlat <- file.path(workdir, "xlat.genePred")
        utils::write.table(rows, gp_xlat,
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
        )
        gp_trim <- gp_xlat
    }

    # 4. Import into the groot. Returns a named list of data frames:
    #    tss, exons, utr3, utr5 (NULL if the set is empty).
    if (verbose) message("  Importing genes via gintervals.import_genes() ...")
    result <- gintervals.import_genes(gp_trim)

    # 5. Save each role under the (prefixed, renamed) name. NA skips the role.
    roles <- c("tss", "exons", "utr3", "utr5")
    for (role in roles) {
        target_name <- gene_sets[[role]]
        if (is.na(target_name)) next
        df <- result[[role]]
        if (is.null(df) || !nrow(df)) {
            if (verbose) message(sprintf("  %s: empty, skipping", role))
            next
        }
        full_name <- paste0(prefix, target_name)
        # For dotted names (hierarchy), gintervals.save needs the directory to exist.
        parts <- strsplit(full_name, ".", fixed = TRUE)[[1]]
        if (length(parts) > 1L) {
            rel_dir <- do.call(file.path, as.list(parts[-length(parts)]))
            abs_dir <- file.path(get("GROOT", envir = .misha), "tracks", rel_dir)
            if (!dir.exists(abs_dir)) dir.create(abs_dir, recursive = TRUE, mode = "0755")
        }
        if (gintervals.exists(full_name)) {
            if (!overwrite) {
                stop(sprintf("Intervals set '%s' already exists. Use overwrite = TRUE to replace.", full_name),
                    call. = FALSE
                )
            }
            gintervals.rm(full_name, force = TRUE)
        }
        gintervals.save(full_name, df)
        if (verbose) message(sprintf("  %s: saved %d intervals", full_name, nrow(df)))
    }

    # 6. Reload so new sets are visible.
    gdb.reload()
    invisible(NULL)
}
