# R/genome-build-seq.R
# Per-source FASTA-to-groot builders. Dispatched from gdb.build_genome().

# Dispatch to per-source seq builder. Returns a list with files_record and an
# optional `chrom_alias_df` (used downstream to seed chrom_aliases.tsv).
# `target_chroms` (optional) lets ucsc-hub auto-pick the chromAlias column
# whose values best cover those names; ignored by other backends.
# `target_lengths` + `match_by_length=TRUE` together opt into the force-align
# path in ucsc-hub: every target chrom is mapped to an alias row by name or
# unique-length pairing, and the FASTA is renamed to target_chroms outright.
.build_seq <- function(recipe, path, target_chroms = NULL, target_lengths = NULL,
                       format = NULL, prefetched_alias = NULL,
                       match_by_length = TRUE, verbose = TRUE) {
    fn <- switch(recipe$source,
        ucsc = .build_seq_ucsc,
        `ucsc-hub` = .build_seq_ucsc_hub,
        ncbi = .build_seq_ncbi,
        manual = .build_seq_manual,
        s3 = .build_seq_s3,
        local = .build_seq_local,
        stop(sprintf("No seq builder for source '%s'", recipe$source), call. = FALSE)
    )
    if (recipe$source == "ucsc-hub") {
        fn(recipe, path,
            target_chroms = target_chroms,
            target_lengths = target_lengths,
            format = format,
            prefetched_alias = prefetched_alias,
            match_by_length = match_by_length,
            verbose = verbose
        )
    } else {
        # target_chroms / target_lengths / match_by_length are ucsc-hub-only;
        # gdb.build_genome already rejects them for other sources before this
        # dispatch is reached.
        fn(recipe, path,
            target_chroms = target_chroms, format = format,
            verbose = verbose
        )
    }
}

.build_seq_ucsc_hub <- function(recipe, path, target_chroms = NULL,
                                target_lengths = NULL, format,
                                prefetched_alias = NULL,
                                match_by_length = TRUE, verbose) {
    accession <- recipe$accession
    base <- .hub_url_for(accession)
    files <- .hub_list_dir(base, verbose = verbose)
    if (is.null(files)) {
        stop(sprintf("UCSC mammal hub directory not found: %s", base), call. = FALSE)
    }
    fa_file <- files[grepl(sprintf("^%s\\.(fa|fasta)\\.gz$", accession), files)]
    if (!length(fa_file)) {
        # Fall back to .2bit: not currently supported; require .fa.gz.
        stop(sprintf(
            "No %s.fa.gz at %s. (.2bit not yet supported.)",
            accession, base
        ), call. = FALSE)
    }
    workdir <- tempfile("misha_hub_seq_")
    dir.create(workdir, recursive = TRUE)
    on.exit(unlink(workdir, recursive = TRUE), add = TRUE)
    local_fa <- file.path(workdir, fa_file[[1L]])
    .download_to(paste0(base, fa_file[[1L]]), local_fa, verbose = verbose)

    alias_df <- NULL
    if (!is.null(prefetched_alias)) {
        alias_df <- prefetched_alias$df
    } else {
        alias_file <- files[grepl("\\.chromAlias\\.txt$", files)]
        if (length(alias_file)) {
            local_alias <- file.path(workdir, basename(alias_file[[1L]]))
            .download_to(paste0(base, alias_file[[1L]]), local_alias, verbose = verbose)
            alias_df <- .parse_ucsc_chromalias(local_alias)
        }
    }

    # Force-align path: caller supplied both target_chroms and target_lengths
    # and opted into match_by_length. Every target chrom is placed on an alias
    # row by name match across columns or unique-on-both-sides length pairing
    # (errors otherwise), and FASTA headers are rewritten to target_chroms
    # outright. Rows in the alias that aren't in target_chroms keep their
    # FASTA header (alias_df[[src_col]] value), so non-HAL contigs remain
    # queryable under their original accession.
    force_align <- !is.null(alias_df) && !is.null(target_chroms) &&
        !is.null(target_lengths) && isTRUE(match_by_length)
    if (force_align) {
        row_lengths <- if (!is.null(prefetched_alias)) {
            prefetched_alias$row_lengths
        } else {
            NULL
        }
        if (is.null(row_lengths)) {
            sizes_file <- files[grepl("\\.chrom\\.sizes\\.txt$", files)]
            if (!length(sizes_file)) {
                stop("target_lengths requires chrom.sizes.txt at the hub; not found.",
                    call. = FALSE
                )
            }
            local_sizes <- file.path(workdir, basename(sizes_file[[1L]]))
            .download_to(paste0(base, sizes_file[[1L]]), local_sizes, verbose = verbose)
            sizes_df <- utils::read.table(local_sizes,
                sep = "\t", header = FALSE,
                col.names = c("name", "length"),
                stringsAsFactors = FALSE, comment.char = "", quote = ""
            )
            row_lengths <- .alias_row_lengths_from_sizes(alias_df, sizes_df)
            if (is.null(row_lengths)) {
                stop("Could not align chrom.sizes.txt to a chromAlias column for length-based FASTA rename.",
                    call. = FALSE
                )
            }
        }
        target_per_row <- .assign_target_chroms_per_row(
            alias_df, target_chroms, target_lengths, row_lengths
        )
        hdrs <- .read_fasta_headers(local_fa, n = min(20L, nrow(alias_df)))
        src_col <- .detect_alias_column(alias_df, hdrs)
        if (is.na(src_col)) {
            stop("Could not auto-detect FASTA's source column in chromAlias for force-align rename.",
                call. = FALSE
            )
        }
        src_col_chr <- as.character(src_col)
        # Map src_col value -> target_chrom (empty means "no rename", per the
        # existing .rename_fasta_headers contract that already handles the
        # chrom_naming-with-empty-cells case).
        map <- setNames(target_per_row, alias_df[[src_col_chr]])
        renamed <- file.path(workdir, "renamed.fa")
        .rename_fasta_headers(local_fa, renamed, map, verbose = verbose)
        local_fa <- renamed
        if (verbose) {
            n_placed <- sum(nzchar(target_per_row))
            message(sprintf(
                "  Force-aligned %d of %d FASTA contigs to target_chroms (name + length pairing).",
                n_placed, nrow(alias_df)
            ))
        }
    } else {
        # If `target_chroms` is supplied, auto-pick the chromAlias column whose
        # values cover those names best -- this lets callers say "I want the groot
        # to use whatever names HAL/halStats reports for this species" without
        # having to figure out which column that corresponds to. Wins over any
        # caller-supplied chrom_naming.
        chrom_naming <- recipe$chrom_naming %||% "ucsc"
        if (!is.null(alias_df) && !is.null(target_chroms)) {
            detected <- .detect_alias_column(alias_df, target_chroms,
                min_coverage = 0.5
            )
            if (is.na(detected)) {
                scores <- attr(detected, "scores")
                warning(sprintf(
                    "target_chroms covered no chromAlias column above 50%% (best=%s); falling back to chrom_naming='%s'.",
                    paste(sprintf("%s=%d", names(scores), scores), collapse = ","),
                    chrom_naming
                ), call. = FALSE)
            } else {
                chrom_naming <- as.character(detected)
                if (verbose) {
                    message(sprintf(
                        "  Auto-detected chrom_naming = '%s' (%.2f%% of target_chroms in this column).",
                        chrom_naming, 100 * attr(detected, "overlap")
                    ))
                }
            }
        }
        if (!is.null(alias_df) && chrom_naming != "accession") {
            # The hub FASTA is keyed by the column UCSC chose for headers
            # (typically the 'genbank' column, the assembly's contig accessions).
            # We auto-detect by sniffing first FASTA header against alias columns.
            hdrs <- .read_fasta_headers(local_fa, n = min(20L, nrow(alias_df)))
            src_col <- .detect_alias_column(alias_df, hdrs)
            if (is.na(src_col)) {
                warning("Could not auto-detect FASTA's chrom column in chromAlias; ",
                    "FASTA headers kept as-is.",
                    call. = FALSE
                )
            } else {
                target_col <- .resolve_hub_target_col(chrom_naming, src_col, names(alias_df))
                if (is.na(target_col)) {
                    warning(attr(target_col, "reason"), " FASTA kept as-is.",
                        call. = FALSE
                    )
                } else if (target_col != src_col) {
                    map <- setNames(alias_df[[target_col]], alias_df[[src_col]])
                    renamed <- file.path(workdir, "renamed.fa")
                    .rename_fasta_headers(local_fa, renamed, map, verbose = verbose)
                    local_fa <- renamed
                }
            }
        }
    }

    gdb.create(
        groot = path, fasta = local_fa, genes.file = NULL,
        annots.file = NULL, format = format, verbose = verbose
    )
    list(
        files_record = list(fasta = list(url = paste0(base, fa_file[[1L]]))),
        chrom_alias_df = alias_df
    )
}

# Read the first n FASTA headers (without the '>') from a possibly-gzipped FASTA.
.read_fasta_headers <- function(file, n = 20L) {
    con <- if (grepl("\\.gz$", file)) gzfile(file, "rt") else file(file, "rt")
    on.exit(close(con), add = TRUE)
    out <- character(0)
    chunk <- 10000L
    repeat {
        lines <- readLines(con, n = chunk, warn = FALSE)
        if (!length(lines)) break
        hdrs <- lines[startsWith(lines, ">")]
        ids <- sub("^>([^[:space:]]+).*$", "\\1", hdrs, perl = TRUE)
        out <- c(out, ids)
        if (length(out) >= n) {
            return(utils::head(out, n))
        }
    }
    out
}

.build_seq_ucsc <- function(recipe, path, target_chroms = NULL, format, verbose) {
    assembly <- recipe$assembly
    base <- sprintf("%s/%s", .UCSC_GOLDENPATH, assembly)
    workdir <- tempfile("misha_ucsc_seq_")
    dir.create(workdir, recursive = TRUE)
    on.exit(unlink(workdir, recursive = TRUE), add = TRUE)
    local_fa <- file.path(workdir, sprintf("%s.fa.gz", assembly))
    .download_to(sprintf("%s/bigZips/%s.fa.gz", base, assembly), local_fa, verbose = verbose)
    gdb.create(
        groot = path, fasta = local_fa, genes.file = NULL,
        annots.file = NULL, format = format, verbose = verbose
    )
    list(
        files_record = list(fasta = list(url = sprintf("%s/bigZips/%s.fa.gz", base, assembly))),
        chrom_alias_df = NULL
    )
}

.build_seq_ncbi <- function(recipe, path, target_chroms = NULL, format, verbose) {
    accession <- recipe$accession
    chrom_naming <- recipe$chrom_naming %||% .NCBI_DEFAULT_CHROM_NAMING
    workdir <- tempfile("misha_ncbi_seq_")
    dir.create(workdir, recursive = TRUE)
    on.exit(unlink(workdir, recursive = TRUE), add = TRUE)
    zip_path <- file.path(workdir, "datasets.zip")
    .download_to(
        .ncbi_datasets_zip_url(accession,
            include = c("GENOME_FASTA", "SEQUENCE_REPORT")
        ),
        zip_path,
        verbose = verbose
    )
    extract_dir <- file.path(workdir, "extract")
    dir.create(extract_dir)
    utils::unzip(zip_path, exdir = extract_dir)

    fasta_files <- list.files(extract_dir,
        pattern = "\\.(fna|fasta|fa)(\\.gz)?$",
        recursive = TRUE, full.names = TRUE
    )
    if (!length(fasta_files)) {
        stop(sprintf("No FASTA in NCBI Datasets payload for %s", accession), call. = FALSE)
    }
    fasta_file <- fasta_files[[1L]]
    if (grepl("\\.gz$", fasta_file)) fasta_file <- .gunzip_to_file(fasta_file)

    seqrep_files <- list.files(extract_dir,
        pattern = "sequence_report\\.jsonl$",
        recursive = TRUE, full.names = TRUE
    )
    seqrep <- if (length(seqrep_files)) .parse_ncbi_sequence_report(seqrep_files[[1L]]) else NULL
    rename_map <- NULL
    if (!is.null(seqrep)) {
        rename_map <- .build_ncbi_rename_map(seqrep, chrom_naming)
    } else if (chrom_naming != "accession") {
        warning(sprintf(
            "No sequence_report in NCBI payload for %s; falling back to chrom_naming='accession'.",
            accession
        ), call. = FALSE)
        chrom_naming <- "accession"
    }

    if (!is.null(rename_map) && chrom_naming != "accession") {
        renamed <- file.path(workdir, "renamed.fna")
        .rename_fasta_headers(fasta_file, renamed, rename_map, verbose = verbose)
        fasta_file <- renamed
    }

    gdb.create(
        groot = path, fasta = fasta_file, genes.file = NULL,
        annots.file = NULL, format = format, verbose = verbose
    )

    # Persist chrom_aliases.tsv in the long format the rest of the build/
    # install pipeline uses. After the FASTA rename above, the groot chroms
    # equal rename_map[refseqAccession]; inject that as a synthetic
    # ".canonical" column so .merge_chrom_aliases_tsv treats every other
    # column as an alias of the groot name.
    if (!is.null(seqrep) && !is.null(rename_map)) {
        alias_df <- .ncbi_seqrep_to_alias_df(seqrep)
        alias_df[[".canonical"]] <- unname(rename_map[seqrep$refseqAccession])
        .merge_chrom_aliases_tsv(path, alias_df, ".canonical")
    }

    list(
        files_record = list(fasta = list(name = basename(fasta_file))),
        chrom_alias_df = NULL # written above via .merge_chrom_aliases_tsv
    )
}

.build_seq_manual <- function(recipe, path, target_chroms = NULL, format, verbose) {
    gdb.create(
        groot = path, fasta = recipe$fasta, genes.file = NULL,
        annots.file = NULL, format = format, verbose = verbose
    )
    list(
        files_record = list(fasta = list(url = recipe$fasta)),
        chrom_alias_df = NULL
    )
}

.build_seq_s3 <- function(recipe, path, target_chroms = NULL, format, verbose) {
    parent_dir <- dirname(path)
    if (!dir.exists(parent_dir)) dir.create(parent_dir, recursive = TRUE)
    if (basename(path) != recipe$assembly) {
        warning(sprintf(
            "S3 backend extracts to %s/%s; got path=%s.",
            parent_dir, recipe$assembly, path
        ), call. = FALSE)
    }
    gdb.create_genome(recipe$assembly, path = parent_dir)
    list(files_record = list(s3 = recipe$assembly), chrom_alias_df = NULL)
}

.build_seq_local <- function(recipe, path, target_chroms = NULL, format, verbose) {
    if (!dir.exists(recipe$path)) {
        stop(sprintf("Local groot does not exist: %s", recipe$path), call. = FALSE)
    }
    if (!file.exists(file.path(recipe$path, "chrom_sizes.txt"))) {
        stop(sprintf("Path %s does not look like a misha groot", recipe$path), call. = FALSE)
    }
    gdb.init(recipe$path, rescan = TRUE)
    list(files_record = list(local = recipe$path), chrom_alias_df = NULL)
}
