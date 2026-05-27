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
    # Dotted set names map to filesystem subdirectories under tracks/.
    # Create them if needed so gintervals.save() doesn't error on a fresh groot.
    parts <- strsplit(name, ".", fixed = TRUE)[[1L]]
    if (length(parts) > 1L) {
        rel_dir <- do.call(file.path, as.list(parts[-length(parts)]))
        abs_dir <- file.path(get("GROOT", envir = .misha), "tracks", rel_dir)
        if (!dir.exists(abs_dir)) {
            dir.create(abs_dir, recursive = TRUE, mode = "0755")
        }
    }
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

    # Split once instead of N full-table scans (was ~20 passes over ~5M rows
    # for hg38). split() drops NA groups by default; explicitly drop empties.
    valid <- !is.na(df$class) & nzchar(df$class)
    if (any(valid)) {
        groups <- split(df[valid, , drop = FALSE], df$class[valid])
        for (cls in names(groups)) {
            suffix <- gsub("\\?$", "_qmark", tolower(cls))
            .save_intervals(sprintf("%srmsk_%s", prefix, suffix), groups[[cls]],
                overwrite = overwrite, verbose = verbose
            )
        }
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
        # -warnAndContinue + -maxConvertErrors=-1 keep the conversion alive when
        # NCBI's RefSeq GFF has a handful of records the tool can't convert
        # (typical case: immunoglobulin V/D/J segments whose CDS coords don't
        # fit inside parent exons). Without these flags a 5-record problem
        # aborts the entire ~70k-record run.
        ret <- system2(converter,
            args = c("-warnAndContinue", "-maxConvertErrors=-1", shQuote(src), shQuote(gp_raw)),
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
        # data.table is in Suggests; use fread/fwrite when present (~5x faster
        # on ~50k-row genePred), fall back to base read.table/write.table.
        use_dt <- requireNamespace("data.table", quietly = TRUE)
        if (use_dt) {
            # Force cols 9/10 (exonStarts/exonEnds) to character so a
            # genePred whose every row is single-exon (e.g. "100,",
            # "500,") doesn't get type-inferred as numeric and lose the
            # trailing commas the C++ importer requires.
            rows <- data.table::fread(gp_trim,
                sep = "\t", header = FALSE, quote = "",
                showProgress = FALSE, data.table = FALSE,
                colClasses = list(character = c(9, 10))
            )
        } else {
            rows <- utils::read.table(gp_trim,
                sep = "\t", header = FALSE,
                stringsAsFactors = FALSE, quote = "", comment.char = ""
            )
        }
        rows <- asset$translate(rows, 2L)
        # Drop rows whose chrom didn't translate to a groot name (NA from
        # rev_idx misses, "" from rows whose canonical column was unset by
        # the 3-pass resolution). gintervals.import_genes (read_genes_file
        # in src/IntervalsImport.cpp) rejects empty CHROM as
        # "invalid file format" on line 1; NA strings it skips silently,
        # but filtering both is cleaner and avoids the asymmetry.
        keep <- !is.na(rows[[2L]]) & nzchar(rows[[2L]])
        if (any(!keep)) {
            if (verbose) {
                message(sprintf(
                    "  %d genePred rows dropped: their chrom didn't translate to a groot contig (keeping %d).",
                    sum(!keep), sum(keep)
                ))
            }
            rows <- rows[keep, , drop = FALSE]
        }
        gp_xlat <- file.path(workdir, "xlat.genePred")
        if (use_dt) {
            data.table::fwrite(rows, gp_xlat,
                sep = "\t", quote = FALSE, col.names = FALSE,
                na = "NA", showProgress = FALSE
            )
        } else {
            utils::write.table(rows, gp_xlat,
                sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
            )
        }
        gp_trim <- gp_xlat
    }

    # 4. Build a (name -> gene symbol) sidecar from the trimmed genePred so the
    #    installed sets carry gene names. Col 1 is the transcript/RNA accession
    #    (the join key the annots mechanism keys on), col 12 is name2 (the gene
    #    symbol). Deduplicate by col 1: gintervals.import_genes' annots reader
    #    rejects a repeated id. Empty name2 stays an empty field (na = ""), so
    #    sources without symbols still build, just with a blank geneName.
    use_dt <- requireNamespace("data.table", quietly = TRUE)
    if (use_dt) {
        name_map <- data.table::fread(gp_trim,
            sep = "\t", header = FALSE, quote = "", showProgress = FALSE,
            data.table = FALSE, select = c(1L, 12L), colClasses = "character",
            na.strings = NULL
        )
    } else {
        name_map <- utils::read.table(gp_trim,
            sep = "\t", header = FALSE, quote = "", comment.char = "",
            colClasses = "character", na.strings = character(0)
        )[, c(1L, 12L), drop = FALSE]
    }
    name_map <- name_map[!duplicated(name_map[[1L]]), , drop = FALSE]
    annots_file <- file.path(workdir, "gene_names.annots")
    utils::write.table(name_map, annots_file,
        sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE, na = ""
    )

    # 5. Import into the groot. Returns a named list of data frames:
    #    tss, exons, utr3, utr5 (NULL if the set is empty), each with the
    #    name / geneName annotation columns attached.
    if (verbose) message("  Importing genes via gintervals.import_genes() ...")
    result <- gintervals.import_genes(gp_trim, annots_file,
        annots.names = c("name", "geneName")
    )

    # 6. Save each role under the (prefixed, renamed) name. NA skips the role.
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

    # 7. Reload so new sets are visible.
    gdb.reload()
    invisible(NULL)
}
