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

.build_seq_ucsc <- function(recipe, path, format, verbose) {
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

.build_seq_ncbi <- function(recipe, path, format, verbose) {
    accession <- recipe$accession
    chrom_naming <- recipe$chrom_naming %||% .NCBI_DEFAULT_CHROM_NAMING
    workdir <- tempfile("misha_ncbi_seq_")
    dir.create(workdir, recursive = TRUE)
    on.exit(unlink(workdir, recursive = TRUE), add = TRUE)
    zip_path <- file.path(workdir, "datasets.zip")
    .download_to(.ncbi_datasets_zip_url(accession), zip_path, verbose = verbose)
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

    if (!is.null(seqrep) && !is.null(rename_map)) {
        .write_chrom_aliases_tsv(path, seqrep, rename_map)
    }

    list(
        files_record = list(fasta = list(name = basename(fasta_file))),
        chrom_alias_df = NULL # NCBI uses sequence_report -> chrom_aliases.tsv directly
    )
}

.build_seq_manual <- function(recipe, path, format, verbose) {
    gdb.create(
        groot = path, fasta = recipe$fasta, genes.file = NULL,
        annots.file = NULL, format = format, verbose = verbose
    )
    list(
        files_record = list(fasta = list(url = recipe$fasta)),
        chrom_alias_df = NULL
    )
}

.build_seq_s3 <- function(recipe, path, format, verbose) {
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

.build_seq_local <- function(recipe, path, format, verbose) {
    if (!dir.exists(recipe$path)) {
        stop(sprintf("Local groot does not exist: %s", recipe$path), call. = FALSE)
    }
    if (!file.exists(file.path(recipe$path, "chrom_sizes.txt"))) {
        stop(sprintf("Path %s does not look like a misha groot", recipe$path), call. = FALSE)
    }
    gdb.init(recipe$path, rescan = TRUE)
    list(files_record = list(local = recipe$path), chrom_alias_df = NULL)
}
