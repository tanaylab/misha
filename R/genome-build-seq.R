# R/genome-build-seq.R
# Per-source FASTA-to-groot builders. Dispatched from gdb.build_genome().

# Dispatch to per-source seq builder. Returns a list with files_record and an
# optional `chrom_alias_df` (used downstream to seed chrom_aliases.tsv).
.build_seq <- function(recipe, path, format = NULL, verbose = TRUE) {
    fn <- switch(recipe$source,
        ucsc = .build_seq_ucsc,
        `ucsc-hub` = .build_seq_ucsc_hub,
        ncbi = .build_seq_ncbi,
        manual = .build_seq_manual,
        s3 = .build_seq_s3,
        local = .build_seq_local,
        stop(sprintf("No seq builder for source '%s'", recipe$source), call. = FALSE)
    )
    fn(recipe, path, format = format, verbose = verbose)
}

.build_seq_ucsc_hub <- function(recipe, path, format, verbose) {
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
    alias_file <- files[grepl("\\.chromAlias\\.txt$", files)]
    if (length(alias_file)) {
        local_alias <- file.path(workdir, basename(alias_file[[1L]]))
        .download_to(paste0(base, alias_file[[1L]]), local_alias, verbose = verbose)
        alias_df <- .parse_ucsc_chromalias(local_alias)
    }

    # If chrom_naming requested and chromAlias has the column, rename FASTA
    # headers; otherwise keep as-shipped (RefSeq accessions).
    chrom_naming <- recipe$chrom_naming %||% "ucsc"
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
            target_col <- if (chrom_naming == "ucsc") {
                "ucsc"
            } else if (chrom_naming == "sequence_name") {
                "assembly"
            } else {
                src_col
            }
            if (!target_col %in% names(alias_df)) {
                warning(sprintf(
                    "chromAlias has no '%s' column; FASTA kept as-is.",
                    target_col
                ), call. = FALSE)
            } else if (target_col != src_col) {
                map <- setNames(alias_df[[target_col]], alias_df[[src_col]])
                renamed <- file.path(workdir, "renamed.fa")
                .rename_fasta_headers(local_fa, renamed, map, verbose = verbose)
                local_fa <- renamed
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
