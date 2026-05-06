# R/genome-build-chromalias.R
# Source-agnostic chromosome-alias parsing, column detection, and translation.

# Resolve which chromAlias column a ucsc-hub build should use as the canonical
# chrom-name source, given a user-facing `chrom_naming` label, the FASTA's own
# source column (auto-detected), and the columns available in chromAlias.
#
# Friendly aliases:
#   "ucsc"          -> column "ucsc"
#   "sequence_name" -> column "assembly"
#   "accession"     -> keep src_col (no FASTA rename)
# Anything else is treated as a literal chromAlias column name (e.g. "genbank",
# "refseq", "ncbi"). HAL/Cactus pipelines typically key on GenBank accessions,
# so chrom_naming = "genbank" matches what halStats reports.
#
# Returns the target column name, or NA_character_ with a `reason` attribute
# describing why no column was resolvable.
.resolve_hub_target_col <- function(chrom_naming, src_col, alias_cols) {
    if (chrom_naming == "accession") {
        return(src_col)
    }
    target <- switch(chrom_naming,
        ucsc          = "ucsc",
        sequence_name = "assembly",
        chrom_naming
    )
    if (!target %in% alias_cols) {
        return(structure(NA_character_,
            reason = sprintf(
                "chromAlias has no '%s' column (available: %s).",
                target, paste(alias_cols, collapse = ", ")
            )
        ))
    }
    target
}

# Parse a UCSC chromAlias.txt(.gz). Header line begins with '# '; remaining lines
# are tab-separated, one row per contig, one column per naming scheme.
# Returns a data.frame with one named character column per scheme.
.parse_ucsc_chromalias <- function(path) {
    if (!file.exists(path)) {
        stop(sprintf("chromAlias file does not exist: %s", path), call. = FALSE)
    }
    con <- if (grepl("\\.gz$", path)) gzfile(path, "rt") else file(path, "rt")
    on.exit(close(con), add = TRUE)
    # First non-empty line is the header (starts with '#').
    header_line <- ""
    repeat {
        l <- readLines(con, n = 1L, warn = FALSE)
        if (!length(l)) {
            stop(sprintf("chromAlias file %s is empty", path), call. = FALSE)
        }
        if (nzchar(l)) {
            header_line <- l
            break
        }
    }
    schemes <- strsplit(sub("^#\\s*", "", header_line), "\t", fixed = TRUE)[[1L]]
    schemes <- trimws(schemes)
    body <- readLines(con, warn = FALSE)
    body <- body[nzchar(body)]
    fields <- strsplit(body, "\t", fixed = TRUE)
    n_cols <- length(schemes)
    short <- vapply(fields, length, integer(1)) < n_cols
    if (any(short)) {
        # Pad short rows with empty strings so all rows have n_cols entries.
        fields[short] <- lapply(
            fields[short],
            function(x) c(x, rep("", n_cols - length(x)))
        )
    }
    mat <- do.call(rbind, lapply(fields, `[`, seq_len(n_cols)))
    df <- as.data.frame(mat, stringsAsFactors = FALSE)
    names(df) <- schemes
    df
}

# Map each alias_df row to a length, given a sizes data.frame (cols: name,
# length) keyed on whichever column the FASTA used for headers (typically
# refseq for GCF or genbank for GCA hubs). Returns a numeric vector aligned
# with alias_df rows; NA where the row's source-column value isn't in sizes.
.alias_row_lengths_from_sizes <- function(alias_df, sizes_df) {
    # Pick the alias column whose values match sizes_df$name best. Use a
    # generous threshold -- the chrom.sizes file should match one column
    # closely but doesn't have to be perfect.
    src <- .detect_alias_column(alias_df, sizes_df$name, min_coverage = 0.5)
    if (is.na(src)) {
        return(NULL)
    }
    name_to_length <- setNames(as.numeric(sizes_df$length), sizes_df$name)
    unname(name_to_length[alias_df[[as.character(src)]]])
}

# Returns the column whose values cover at least `min_coverage` of target_chroms.
# Among columns meeting the threshold, the one with highest coverage wins; ties
# broken by column order. Default `min_coverage = 1.0` keeps strict semantics.
#
# Coverage is bp-weighted when `chrom_lengths` is supplied (a numeric aligned
# with `target_chroms`) -- a long-tail of small unmapped contigs (e.g. a 16 kb
# mitochondrion missing from UCSC's genbank column out of a 3 Gb genome) then
# costs the score ~0.0005% instead of ~1/N. Without lengths, falls back to
# count-based: fraction of distinct chrom names covered.
#
# Returns the column name with attributes(scores=, overlap=, bp_weighted=) on
# success, or NA_character_ with the same attributes when no column meets the
# threshold (caller produces the diagnostic).
.detect_alias_column <- function(alias_df, target_chroms, min_coverage = 1.0,
                                 chrom_lengths = NULL) {
    if (!length(target_chroms)) {
        stop(".detect_alias_column called with empty target chrom set", call. = FALSE)
    }
    bp_weighted <- !is.null(chrom_lengths)
    if (bp_weighted) {
        if (length(chrom_lengths) != length(target_chroms)) {
            stop(sprintf(
                ".detect_alias_column: chrom_lengths (%d) must align with target_chroms (%d)",
                length(chrom_lengths), length(target_chroms)
            ), call. = FALSE)
        }
        weights <- as.numeric(chrom_lengths)
        total <- sum(weights)
        scores <- vapply(alias_df, function(col) {
            sum(weights[target_chroms %in% col])
        }, numeric(1))
    } else {
        unique_targets <- unique(target_chroms)
        scores <- vapply(
            alias_df,
            function(col) sum(unique_targets %in% col),
            integer(1)
        )
        total <- length(unique_targets)
    }
    coverages <- scores / total
    valid <- which(coverages >= min_coverage)
    if (!length(valid)) {
        return(structure(NA_character_,
            scores = scores,
            overlap = max(coverages),
            bp_weighted = bp_weighted
        ))
    }
    best <- valid[which.max(coverages[valid])]
    structure(names(alias_df)[best],
        scores = scores,
        overlap = unname(coverages[best]),
        bp_weighted = bp_weighted
    )
}

# For each alias row whose `canonical` value is empty, look up the unique
# groot chrom whose length matches the row's length. A match counts only
# when the length appears exactly once on BOTH the groot side and the alias
# side -- ambiguous lengths (typical for tiny scaffolds in Cactus assemblies)
# are left empty rather than guessed.
#
# canonical, alias_row_lengths must be aligned with alias_df rows.
# groot_chroms, groot_lengths must be aligned with each other.
.length_match_fill <- function(canonical, alias_row_lengths,
                               groot_chroms, groot_lengths) {
    if (length(canonical) != length(alias_row_lengths)) {
        stop(".length_match_fill: canonical/alias_row_lengths length mismatch",
            call. = FALSE
        )
    }
    needs <- is.na(canonical) | !nzchar(canonical)
    if (!any(needs)) {
        return(canonical)
    }
    # Lengths that appear exactly once on each side. Stringify because lengths
    # are integers/numerics that don't survive named-list indexing well.
    g_counts <- table(groot_lengths)
    g_unique <- names(g_counts)[g_counts == 1L]
    a_counts <- table(alias_row_lengths)
    a_unique <- names(a_counts)[a_counts == 1L]
    fill_lengths <- intersect(g_unique, a_unique)
    if (!length(fill_lengths)) {
        return(canonical)
    }
    # Build "length -> groot chrom" lookup over the unique-on-both side.
    keep <- as.character(groot_lengths) %in% fill_lengths
    groot_lookup <- setNames(groot_chroms[keep], as.character(groot_lengths[keep]))
    fill_idx <- which(needs & as.character(alias_row_lengths) %in% fill_lengths)
    canonical[fill_idx] <- groot_lookup[as.character(alias_row_lengths[fill_idx])]
    canonical
}

# Build a reverse index for cross-column per-row alias translation: maps every
# non-empty cell of `alias_df` (excluding `canonical_col`) directly to the
# row's canonical name. First occurrence of any duplicate value wins.
.build_alias_rev_index <- function(alias_df, canonical_col) {
    cols <- setdiff(names(alias_df), canonical_col)
    canonical <- as.character(alias_df[[canonical_col]])
    rev_idx <- new.env(hash = TRUE, parent = emptyenv())
    for (col in cols) {
        vals <- alias_df[[col]]
        for (i in seq_along(vals)) {
            v <- vals[i]
            if (!is.na(v) && nzchar(v) && !exists(v, envir = rev_idx, inherits = FALSE)) {
                assign(v, canonical[i], envir = rev_idx)
            }
        }
    }
    rev_idx
}

# Per-row asset translation. For each value in `rows[[chrom_col]]`, look up
# the canonical name via `rev_idx` (built once per gdb.install_intervals call
# by .build_alias_rev_index). Unmatched chroms become NA -- caller decides
# whether to drop or warn.
.translate_chroms_per_row <- function(rows, chrom_col, rev_idx) {
    asset_names <- as.character(rows[[chrom_col]])
    out <- vapply(asset_names, function(nm) {
        if (is.na(nm) || !nzchar(nm) || !exists(nm, envir = rev_idx, inherits = FALSE)) {
            return(NA_character_)
        }
        get(nm, envir = rev_idx)
    }, character(1), USE.NAMES = FALSE)
    rows[[chrom_col]] <- out
    rows
}

# Translate `rows[[chrom_col]]` from source_col naming scheme to groot_col naming
# scheme via alias_df. If source_col == groot_col, no-op. Caller MUST ensure
# every rows[[chrom_col]] value is present in alias_df[[source_col]] (use
# .detect_alias_column with 100% threshold first); this function does not handle
# missing keys.
.translate_chroms <- function(rows, chrom_col, alias_df, source_col, groot_col) {
    if (identical(source_col, groot_col)) {
        return(rows)
    }
    map <- setNames(alias_df[[groot_col]], alias_df[[source_col]])
    rows[[chrom_col]] <- map[rows[[chrom_col]]]
    rows
}

# Extend (or create) <groot>/chrom_aliases.tsv with rows mapping every alias
# column of alias_df to the canonical name in groot_col. Existing rows are
# preserved; new rows are added; conflicts (existing alias mapping to a
# different canonical) are warned and skipped.
.merge_chrom_aliases_tsv <- function(groot, alias_df, groot_col) {
    out_path <- file.path(groot, "chrom_aliases.tsv")
    existing <- if (file.exists(out_path)) {
        utils::read.table(out_path,
            sep = "\t", header = TRUE,
            stringsAsFactors = FALSE, quote = "", comment.char = ""
        )
    } else {
        data.frame(
            canonical = character(0), alias = character(0),
            source = character(0), stringsAsFactors = FALSE
        )
    }
    canonical <- alias_df[[groot_col]]
    new_rows <- list()
    conflicts <- 0L
    for (col in setdiff(names(alias_df), groot_col)) {
        vals <- alias_df[[col]]
        keep <- nzchar(vals) & !is.na(vals)
        if (!any(keep)) next
        df <- data.frame(
            canonical = canonical[keep],
            alias = vals[keep],
            source = col,
            stringsAsFactors = FALSE
        )
        # Detect conflicts against existing entries.
        idx <- match(df$alias, existing$alias)
        is_existing <- !is.na(idx)
        if (any(is_existing)) {
            conflict_mask <- existing$canonical[idx[is_existing]] != df$canonical[is_existing]
            conflicts <- conflicts + sum(conflict_mask, na.rm = TRUE)
            df <- df[!is_existing | !conflict_mask, , drop = FALSE]
            df <- df[!df$alias %in% existing$alias, , drop = FALSE] # dedupe non-conflicts too
        }
        new_rows[[length(new_rows) + 1L]] <- df
    }
    if (conflicts > 0L) {
        warning(sprintf(
            "chrom_aliases.tsv: %d alias->canonical conflict(s) preserved as existing.",
            conflicts
        ), call. = FALSE)
    }
    out <- rbind(existing, do.call(rbind, new_rows))
    utils::write.table(out, out_path, sep = "\t", quote = FALSE, row.names = FALSE)
    invisible(out_path)
}
