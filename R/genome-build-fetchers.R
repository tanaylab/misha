# R/genome-build-fetchers.R
# Per-source asset fetchers. Return a uniform list consumed by per-set installers.

# Internal: parse an Apache directory-index HTML body via regex.
# Exposed for testing without network.
.hub_list_dir_parse <- function(body) {
    hits <- regmatches(body, gregexpr('href="([^"?/][^"]*)"', body, perl = TRUE))[[1L]]
    hits <- sub('^href="', "", hits)
    hits <- sub('"$', "", hits)
    hits[!hits %in% c("/", "../") & !startsWith(hits, "?")]
}

# Fetch a UCSC hub-style directory listing. Returns a character vector of
# filenames (or NULL if URL is unreachable or 404). Uses curl (already in
# Imports) — no xml2/rvest.
.hub_list_dir <- function(url, verbose = TRUE) {
    if (verbose) message(sprintf("Listing %s ...", url))
    h <- curl::new_handle()
    resp <- tryCatch(curl::curl_fetch_memory(url, handle = h), error = function(e) NULL)
    if (is.null(resp) || resp$status_code >= 400) {
        return(NULL)
    }
    body <- rawToChar(resp$content)
    .hub_list_dir_parse(body)
}

# URL for a UCSC mammal-hub directory by GenBank/RefSeq accession.
.hub_url_for <- function(accession) {
    m <- regmatches(accession, regexec(
        "^(GC[FA])_([0-9]{3})([0-9]{3})([0-9]{3})\\.[0-9]+$",
        accession
    ))[[1L]]
    if (length(m) != 5L) {
        stop(sprintf("Invalid accession '%s'; expected GC[FA]_<9 digits>.<n>", accession),
            call. = FALSE
        )
    }
    sprintf(
        "https://hgdownload.soe.ucsc.edu/hubs/%s/%s/%s/%s/%s/",
        m[[2L]], m[[3L]], m[[4L]], m[[5L]], accession
    )
}

# From a list of filenames in a hub's genes/ subdir, pick the first matching
# entry of `priority`. Returns list(file = "<path>", source = "<which>") or NULL.
# Matches both bare (`.ensGene.gtf.gz`) and version-suffixed
# (`.ensGene.2020_05.gtf.gz`) filenames; UCSC publishes both shapes.
.pick_gtf <- function(files, priority) {
    base <- sub("\\.gtf(\\.gz)?$", "", files)
    for (src in priority) {
        pat <- sprintf("\\.%s(\\.[^.]+)*$", src)
        hit <- files[grepl(pat, base)]
        if (length(hit)) {
            return(list(file = hit[[1L]], source = src))
        }
    }
    NULL
}

# Fetch every asset needed for the requested sets from a UCSC mammal hub.
# Downloads files into a workdir (caller manages cleanup via on.exit).
# Returns a list with named entries per asset, each:
#   list(file = <local path>, format = <tag>, ...)
# plus chrom_alias = list(file=..., df=<parsed data.frame>).
.hub_fetch_assets <- function(recipe, sets, workdir, gtf_priority,
                              prefetched_alias = NULL, verbose = TRUE) {
    base <- .hub_url_for(recipe$accession)
    files <- .hub_list_dir(base, verbose = verbose)
    if (is.null(files)) {
        stop(sprintf(
            "UCSC mammal hub directory not found: %s\nTry source: ncbi for this accession.",
            base
        ), call. = FALSE)
    }
    out <- list()

    if (!is.null(prefetched_alias)) {
        out$chrom_alias <- list(
            file = prefetched_alias$alias_file,
            df = prefetched_alias$df,
            row_lengths = prefetched_alias$row_lengths
        )
    } else {
        # chromAlias is always fetched (we need it even for non-translation cases).
        alias_file <- files[grepl("\\.chromAlias\\.txt$", files)]
        if (!length(alias_file)) {
            out$chrom_alias <- NULL
            if (verbose) message("  No chromAlias.txt at hub; alias-driven translation disabled.")
        } else {
            local <- file.path(workdir, basename(alias_file[[1L]]))
            .download_to(paste0(base, alias_file[[1L]]), local, verbose = verbose)
            alias_df <- .parse_ucsc_chromalias(local)
            out$chrom_alias <- list(file = local, df = alias_df)
            # chrom.sizes.txt (2-col TSV: name<TAB>length, keyed on the FASTA's
            # source column -- typically refseq for GCF, genbank for GCA). Used by
            # gdb.install_intervals when match_by_length=TRUE to fill alias rows
            # whose chosen canonical column is empty (e.g. MT row with no genbank).
            sizes_file <- files[grepl("\\.chrom\\.sizes\\.txt$", files)]
            if (length(sizes_file)) {
                local_sizes <- file.path(workdir, basename(sizes_file[[1L]]))
                .download_to(paste0(base, sizes_file[[1L]]), local_sizes, verbose = verbose)
                sizes <- utils::read.table(local_sizes,
                    sep = "\t", header = FALSE,
                    col.names = c("name", "length"),
                    stringsAsFactors = FALSE, comment.char = "",
                    quote = ""
                )
                out$chrom_alias$row_lengths <- .alias_row_lengths_from_sizes(alias_df, sizes)
            }
        }
    }

    # Pick one filename out of `files` matching `pattern`, download it from
    # `base`, and return list(file=<local>, format=<format>) -- or warn-and-NULL
    # when nothing matches. `label_pattern` is the human-readable fragment for
    # the warning ("repeatMasker.out.gz", not the regex).
    grab_optional <- function(set_label, pattern, label_pattern, format) {
        hit <- files[grepl(pattern, files)]
        if (!length(hit)) {
            warning(sprintf(
                "'%s' requested but no %s at %s; skipping.",
                set_label, label_pattern, base
            ), call. = FALSE)
            return(NULL)
        }
        local <- file.path(workdir, basename(hit[[1L]]))
        .download_to(paste0(base, hit[[1L]]), local, verbose = verbose)
        list(file = local, format = format)
    }

    if ("rmsk" %in% sets) {
        out$rmsk <- grab_optional(
            "rmsk", "\\.repeatMasker\\.out\\.gz$", "repeatMasker.out.gz", "rmsk-out"
        )
    }

    if ("genes" %in% sets) {
        genes_files <- .hub_list_dir(paste0(base, "genes/"), verbose = verbose)
        if (is.null(genes_files)) {
            warning(sprintf("'genes' requested but no genes/ subdir at %s; skipping.", base),
                call. = FALSE
            )
        } else {
            picked <- .pick_gtf(genes_files, priority = gtf_priority)
            if (is.null(picked)) {
                warning(sprintf(
                    "'genes' requested but none of %s found at %s/genes/; skipping.",
                    paste(gtf_priority, collapse = ", "), base
                ), call. = FALSE)
            } else {
                local <- file.path(workdir, basename(picked$file))
                .download_to(paste0(base, "genes/", picked$file), local, verbose = verbose)
                out$genes <- list(file = local, format = "gtf", gtf_source = picked$source)
            }
        }
    }

    if ("cgi" %in% sets) {
        out$cgi <- grab_optional(
            "cgi", "\\.cpgIslandExt\\.txt\\.gz$", "cpgIslandExt.txt.gz", "ucsc-cpg-11col"
        )
    }

    if ("cytoband" %in% sets) {
        warning("'cytoband' is not available from UCSC mammal hubs; skipping.", call. = FALSE)
    }

    out
}

# Build the same URLs the existing UCSC golden-path backend uses, download
# only the annotation files needed for `sets`, and return uniform asset shape.
# (FASTA download is the seq-builder's job, not the fetcher's.)
.ucsc_fetch_assets <- function(recipe, sets, workdir, verbose = TRUE) {
    assembly <- recipe$assembly
    base <- sprintf("%s/%s", .UCSC_GOLDENPATH, assembly)
    out <- list()

    # UCSC golden path doesn't ship chromAlias in the same shape as hubs.
    # For now, omit. Existing seq-build path uses .compute_chrom_aliases.
    out$chrom_alias <- NULL

    if ("genes" %in% sets) {
        url <- sprintf("%s/database/ncbiRefSeq.txt.gz", base)
        local <- file.path(workdir, "ncbiRefSeq.txt.gz")
        .download_to(url, local, verbose = verbose)
        # Trim extended-genePred 16->12 cols.
        gp <- file.path(workdir, "ncbiRefSeq.12col.txt")
        .normalize_ucsc_genepred(local, gp)
        out$genes <- list(file = gp, format = "genepred", gtf_source = "ncbiRefSeq")

        annots_url <- sprintf("%s/database/ncbiRefSeqLink.txt.gz", base)
        raw <- file.path(workdir, "ncbiRefSeqLink.raw.txt.gz")
        .download_to(annots_url, raw, verbose = verbose)
        healed <- file.path(workdir, "ncbiRefSeqLink.healed.txt")
        .heal_ucsc_tsv_escapes(raw, healed)
        norm <- file.path(workdir, "ncbiRefSeqLink.normalized.txt")
        .normalize_ucsc_tsv(healed, norm, n_cols = length(.UCSC_NCBI_REFSEQ_LINK_COLS))
        out$genes_annots <- list(
            file = norm, format = "ucsc-refseq-link",
            names = .UCSC_NCBI_REFSEQ_LINK_COLS
        )
    }
    download_table <- function(table_name, format) {
        url <- sprintf("%s/database/%s.txt.gz", base, table_name)
        local <- file.path(workdir, paste0(table_name, ".txt.gz"))
        .download_to(url, local, verbose = verbose)
        list(file = local, format = format)
    }
    if ("rmsk" %in% sets) out$rmsk <- download_table("rmsk", "rmsk-ucsc-17col")
    if ("cgi" %in% sets) out$cgi <- download_table("cpgIslandExt", "ucsc-cpg-11col")
    if ("cytoband" %in% sets) out$cytoband <- download_table("cytoBandIdeo", "ucsc-cytoband-5col")
    out
}

# NCBI Datasets fetcher. The existing backend already extracts the zip and finds
# the GFF; we reuse that logic. For seq-only builds, only chromAlias-equivalent
# (sequence_report.jsonl) is fetched.
.ncbi_fetch_assets <- function(recipe, sets, workdir, verbose = TRUE) {
    accession <- recipe$accession
    # Pre-flight: hit /dataset_report (a few KB) to learn whether 'genes' can
    # actually be installed before downloading 800+ MB. About 80% of the
    # Phylo447/Zoonomia community-submitted assemblies (Sanger TOL et al.) on
    # NCBI ship without any annotation at all; pulling the FASTA only to find
    # there's no GFF is pure waste. The check is best-effort: any failure
    # (network, parse) falls through to the original code path.
    if ("genes" %in% sets) {
        report <- tryCatch(.ncbi_dataset_report(accession),
            error = function(e) NULL
        )
        if (!is.null(report)) {
            info <- .ncbi_parse_annotation_info(report)
            hint <- if (!info$has_annotation) {
                .ncbi_suggest_annotated_alternative(info$organism_tax_id, accession)
            } else {
                ""
            }
            res <- .ncbi_resolve_sets_with_preflight(sets, info,
                accession = accession, hint = hint
            )
            for (w in res$warnings) warning(w, call. = FALSE)
            sets <- res$sets
        }
    }
    if (!length(sets)) {
        return(list(chrom_alias = NULL))
    }
    zip_path <- file.path(workdir, "datasets.zip")
    .download_to(.ncbi_datasets_zip_url(accession), zip_path, verbose = verbose)
    extract_dir <- file.path(workdir, "extract")
    dir.create(extract_dir, recursive = TRUE)
    utils::unzip(zip_path, exdir = extract_dir)

    out <- list(chrom_alias = NULL)

    seqrep <- list.files(extract_dir,
        pattern = "sequence_report\\.jsonl$",
        recursive = TRUE, full.names = TRUE
    )
    if (length(seqrep)) {
        out$ncbi_sequence_report <- list(
            file = seqrep[[1L]],
            df = .parse_ncbi_sequence_report(seqrep[[1L]])
        )
    }

    if ("genes" %in% sets) {
        gff <- list.files(extract_dir,
            pattern = "\\.gff(\\.gz)?$",
            recursive = TRUE, full.names = TRUE
        )
        if (length(gff)) {
            f <- gff[[1L]]
            if (grepl("\\.gz$", f)) f <- .gunzip_to_file(f)
            out$genes <- list(file = f, format = "gff3", gtf_source = "RefSeq")
        } else {
            warning(sprintf(
                "'genes' requested but no GFF in NCBI payload for %s; skipping.",
                accession
            ), call. = FALSE)
        }
    }
    if ("rmsk" %in% sets) {
        warning("'rmsk' from NCBI is not implemented (v1); skipping. Use ucsc-hub or manual.",
            call. = FALSE
        )
    }
    if ("cgi" %in% sets) {
        warning("'cgi' is not available from NCBI Datasets; skipping.", call. = FALSE)
    }
    if ("cytoband" %in% sets) {
        warning("'cytoband' is not available from NCBI Datasets; skipping.", call. = FALSE)
    }
    out
}

# Download `url` to <workdir>/<stem>_input; gunzip in place when the URL ends
# in .gz. Returns the final local path (post-gunzip if applicable).
.fetch_with_optional_gunzip <- function(url, workdir, stem, verbose = TRUE) {
    local <- file.path(workdir, paste0(stem, "_input"))
    .download_to(url, local, verbose = verbose)
    if (grepl("\\.gz$", url)) {
        local <- .gunzip_to_file(local, paste0(local, ".unzipped"))
    }
    local
}

# Manual fetcher: takes the recipe URLs verbatim, downloads them, returns
# uniform asset shape. Used by the manual source backend.
# Note: %||% is defined in R/db-core.R (shared utility).
.manual_fetch_assets <- function(recipe, sets, workdir, verbose = TRUE) {
    out <- list(chrom_alias = NULL)
    # Trivial-format sets: download (gunzipping when needed), tag with format.
    simple_specs <- list(
        rmsk     = "rmsk-ucsc-17col",
        cgi      = "ucsc-cpg-11col",
        cytoband = "ucsc-cytoband-5col"
    )
    for (key in names(simple_specs)) {
        url <- recipe[[key]]
        if (!key %in% sets || is.null(url)) next
        local <- .fetch_with_optional_gunzip(url, workdir, key, verbose = verbose)
        out[[key]] <- list(file = local, format = simple_specs[[key]])
    }
    if ("genes" %in% sets && !is.null(recipe$genes)) {
        local <- .fetch_with_optional_gunzip(recipe$genes, workdir, "genes", verbose = verbose)
        fmt <- recipe$genes_format %||% "genepred"
        if (!fmt %in% c("genepred", "gtf", "gff3")) {
            stop(sprintf(
                "Unsupported genes_format '%s'. Supported: genepred, gtf, gff3.",
                fmt
            ), call. = FALSE)
        }
        out$genes <- list(file = local, format = fmt, gtf_source = "manual")
    }
    out
}

# Pre-flight chromAlias coverage gate for ucsc-hub builds. Fetches only the
# small alias and chrom.sizes files (no FASTA) and runs .coverage_gate so
# gdb.build_genome can fail before pulling multi-GB sequence data. Returns
# the prefetched alias bundle so callers can pass it to .build_seq_ucsc_hub
# and gdb.install_intervals to avoid re-download.
.hub_preflight_coverage <- function(accession, target_chroms = NULL,
                                    chrom_naming = NULL, min_coverage,
                                    workdir, verbose = TRUE) {
    base <- .hub_url_for(accession)
    files <- .hub_list_dir(base, verbose = verbose)
    if (is.null(files)) {
        stop(sprintf("UCSC mammal hub directory not found: %s", base),
            call. = FALSE
        )
    }
    alias_file <- files[grepl("\\.chromAlias\\.txt$", files)]
    if (!length(alias_file)) {
        # No chromAlias means no gate; nothing to pre-flight, the build will
        # proceed exactly as before.
        return(NULL)
    }
    local_alias <- file.path(workdir, basename(alias_file[[1L]]))
    .download_to(paste0(base, alias_file[[1L]]), local_alias, verbose = verbose)
    alias_df <- .parse_ucsc_chromalias(local_alias)

    sizes_file <- files[grepl("\\.chrom\\.sizes\\.txt$", files)]
    row_lengths <- NULL
    local_sizes <- NULL
    if (length(sizes_file)) {
        local_sizes <- file.path(workdir, basename(sizes_file[[1L]]))
        .download_to(paste0(base, sizes_file[[1L]]), local_sizes, verbose = verbose)
        sizes <- utils::read.table(local_sizes,
            sep = "\t", header = FALSE,
            col.names = c("name", "length"),
            stringsAsFactors = FALSE, comment.char = "", quote = ""
        )
        row_lengths <- .alias_row_lengths_from_sizes(alias_df, sizes)
    }

    # Determine the canonical column the build will end up with.
    canonical_col <- NA_character_
    if (!is.null(target_chroms)) {
        detected <- .detect_alias_column(alias_df, target_chroms,
            min_coverage = 0.5
        )
        if (!is.na(detected)) {
            canonical_col <- as.character(detected)
        }
    }
    if (is.na(canonical_col)) {
        cn <- chrom_naming %||% "ucsc"
        # Friendly aliases mirror .resolve_hub_target_col() in
        # R/genome-build-chromalias.R: "sequence_name" maps to the alias's
        # "assembly" column; "ucsc" stays as-is; anything else is treated as
        # a literal column name.
        target <- switch(cn,
            ucsc          = "ucsc",
            sequence_name = "assembly",
            cn
        )
        if (cn == "accession" || !target %in% names(alias_df)) {
            # Either the build will keep the FASTA's source column (which we
            # can't sniff without the FASTA), or chrom_naming maps to a column
            # the alias doesn't have. Pick the column whose values cover the
            # largest share of all alias values across the frame -- a rough
            # "most-populated column" proxy. The post-build gate inside
            # gdb.install_intervals will make the precise call once the FASTA
            # is on disk; this fallback exists only to fail obviously-broken
            # cases (e.g. an alias with no usable column at all) before the
            # download.
            best <- .detect_alias_column(alias_df,
                target_chroms = unlist(alias_df, use.names = FALSE),
                min_coverage = 0
            )
            canonical_col <- as.character(best)
        } else {
            canonical_col <- target
        }
    }

    # Build a per-row target vector: the canonical col's value for non-empty
    # rows, NA for rows with no canonical value. NA targets fail %in% against
    # every column (including the canonical col's empty cell), so those rows
    # contribute their length to the denominator without contributing to any
    # column's coverage score. Net effect: a 16 kb mitochondrion missing from
    # the canonical col costs 16 kb of denom, exactly the bp-weighted gap we
    # want the gate to detect.
    #
    # Note: this is deliberately stricter than the post-build gate inside
    # gdb.install_intervals, which can rescue empty canonical cells via
    # length-matching against the actual built groot (the match_by_length
    # path). Pre-flight has no FASTA yet, so length-matching isn't available.
    # Trade-off: rare false-positive fail-early on builds that would have
    # succeeded post-rename, in exchange for never paying the multi-GB FASTA
    # download for a build that ultimately couldn't pass the strict gate at
    # the user's chosen min_coverage.
    canonical_vals <- alias_df[[canonical_col]]
    empties <- is.na(canonical_vals) | !nzchar(canonical_vals)
    groot_chroms <- canonical_vals
    groot_chroms[empties] <- NA_character_
    groot_lengths <- row_lengths

    .coverage_gate(alias_df, groot_chroms, groot_lengths,
        min_coverage = min_coverage,
        label = sprintf("groot (%s)", canonical_col)
    )

    list(
        df = alias_df,
        row_lengths = row_lengths,
        alias_file = local_alias,
        sizes_file = local_sizes,
        canonical_col = canonical_col
    )
}
