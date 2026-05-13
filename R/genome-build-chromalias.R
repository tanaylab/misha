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

# Parse an NCBI assembly_report.txt (the same file UCSC mammal hubs mirror
# at `<acc>_assembly_report.txt`). Format: a run of `# Key: Value` metadata
# lines, then a `# <col>\t<col>\t...` header line, then tab-separated data
# rows. The header is the LAST `#`-prefixed line. Standard columns:
# Sequence-Name, Sequence-Role, Assigned-Molecule, Assigned-Molecule-Location/Type,
# GenBank-Accn, Relationship, RefSeq-Accn, Assembly-Unit, Sequence-Length,
# UCSC-style-name.
#
# Returns a data.frame with one character column per header field, with
# NCBI's "na" placeholder (used in GenBank-Accn / RefSeq-Accn / UCSC-style-
# name when no counterpart exists) collapsed to "". Returns NULL when the
# file has no data rows so callers can no-op cleanly.
.parse_ucsc_assembly_report <- function(path) {
    if (!file.exists(path)) {
        stop(sprintf("Assembly report not found: %s", path), call. = FALSE)
    }
    lines <- readLines(path, warn = FALSE)
    lines <- lines[nzchar(lines)]
    is_comment <- startsWith(lines, "#")
    if (!any(is_comment)) {
        return(NULL)
    }
    header_idx <- max(which(is_comment))
    if (header_idx == length(lines)) {
        return(NULL)
    }
    cols <- strsplit(sub("^#\\s*", "", lines[[header_idx]]), "\t", fixed = TRUE)[[1L]]
    cols <- trimws(cols)
    data_lines <- lines[(header_idx + 1L):length(lines)]
    data_lines <- data_lines[!startsWith(data_lines, "#")]
    if (!length(data_lines)) {
        return(NULL)
    }
    fields <- strsplit(data_lines, "\t", fixed = TRUE)
    n_cols <- length(cols)
    short <- vapply(fields, length, integer(1)) < n_cols
    if (any(short)) {
        fields[short] <- lapply(
            fields[short],
            function(x) c(x, rep("", n_cols - length(x)))
        )
    }
    mat <- do.call(rbind, lapply(fields, `[`, seq_len(n_cols)))
    df <- as.data.frame(mat, stringsAsFactors = FALSE)
    names(df) <- cols
    # The mirrored copies at the UCSC hub use CRLF line endings; readLines
    # strips '\n' but keeps the trailing '\r' on the last field. Strip
    # whitespace on every column so downstream %in% checks aren't tripped
    # by invisible '\r'.
    for (col in names(df)) {
        df[[col]] <- trimws(df[[col]])
    }
    # NCBI's "not applicable" placeholder -- normalize to empty so downstream
    # %in% checks don't false-match it as a chrom name.
    for (col in c("GenBank-Accn", "RefSeq-Accn", "UCSC-style-name")) {
        if (col %in% names(df)) {
            df[[col]][df[[col]] == "na"] <- ""
        }
    }
    df
}

# Normalize a free-form report column name to a chromAlias-style identifier:
# lowercase, runs of non-alphanumerics replaced with a single `_`, trailing
# `_` stripped. "UCSC-style-name" -> "ucsc_style_name".
.normalize_report_colname <- function(x) {
    s <- gsub("[^A-Za-z0-9]+", "_", tolower(x))
    sub("_+$", "", s)
}

# Merge an NCBI assembly_report into alias_df, keeping every report column
# as a SEPARATE alias_df column (under a normalized name) so .detect_alias_
# column scores each naming scheme independently. This is deliberate: even
# columns that look equivalent (chromAlias.ucsc vs report.UCSC-style-name)
# routinely differ -- chromAlias uses RefSeq-derived names for unplaced
# scaffolds (chr1_NW_xxx_random) while the report uses GenBank-derived
# names (chr1_AABRxxx_random). The HAL groot may follow either convention;
# letting both compete on coverage gives the gate the right answer.
#
# Picks a join key automatically (alias_df$refseq <-> report$RefSeq-Accn
# first, falling back to genbank). Appends report rows missing from
# alias_df (e.g. unplaced scaffolds the chromAlias drops); for those rows
# only the join key column on the chromAlias side gets a value (the link),
# other chromAlias columns stay empty so a single column never carries
# mixed naming conventions.
#
# No-op when alias_df is NULL, report is empty, or no join key reaches 50%.
.merge_assembly_report_into_alias <- function(alias_df, report_df) {
    if (is.null(alias_df) || is.null(report_df) || nrow(report_df) == 0L) {
        return(alias_df)
    }
    rep_norm <- report_df
    names(rep_norm) <- vapply(names(rep_norm), .normalize_report_colname, character(1))
    # Pick join key: prefer refseq, fall back to genbank. Use whichever has
    # at least 50% match rate of non-empty alias values into the report.
    join <- NULL
    for (pair in list(c("refseq", "refseq_accn"), c("genbank", "genbank_accn"))) {
        a_col <- pair[[1L]]
        r_col <- pair[[2L]]
        if (!(a_col %in% names(alias_df)) || !(r_col %in% names(rep_norm))) next
        a_vals <- as.character(alias_df[[a_col]])
        a_vals[is.na(a_vals)] <- ""
        non_empty <- nzchar(a_vals)
        if (!any(non_empty)) next
        n_match <- sum(a_vals[non_empty] %in% rep_norm[[r_col]])
        if (n_match / sum(non_empty) >= 0.5) {
            join <- list(alias = a_col, report = r_col)
            break
        }
    }
    if (is.null(join)) {
        return(alias_df)
    }
    idx <- match(alias_df[[join$alias]], rep_norm[[join$report]])
    # Add every report column except the join key (whose data already lives
    # in alias_df[[join$alias]]). On rare name collisions with an existing
    # alias column, prefix with "report_" to keep both naming schemes
    # available.
    new_cols <- setdiff(names(rep_norm), join$report)
    out_name <- function(col) {
        if (col %in% names(alias_df)) paste0("report_", col) else col
    }
    out_names <- vapply(new_cols, out_name, character(1))
    for (i in seq_along(new_cols)) {
        col <- new_cols[[i]]
        vals <- character(nrow(alias_df))
        hit <- !is.na(idx)
        vals[hit] <- as.character(rep_norm[[col]])[idx[hit]]
        alias_df[[out_names[[i]]]] <- vals
    }
    # Append report rows not represented in alias_df. A row is "represented"
    # if any of its identifier values (RefSeq-Accn, GenBank-Accn,
    # UCSC-style-name, Sequence-Name) appears anywhere in alias_df's
    # standard chromAlias columns. This is stricter than "join key matches"
    # because the join key alone misses rows like NCBI's MT contig where
    # RefSeq-Accn is "na" but GenBank-Accn is real (e.g.
    # GCF_000002285.5 dog: MT has GenBank CM023446.1, RefSeq "na").
    report_id_cols <- intersect(
        c("refseq_accn", "genbank_accn", "ucsc_style_name", "sequence_name"),
        names(rep_norm)
    )
    alias_id_cols <- intersect(
        c("refseq", "genbank", "ucsc", "assembly", "ncbi"),
        names(alias_df)
    )
    alias_value_set <- unique(unlist(
        lapply(alias_id_cols, function(c) {
            v <- as.character(alias_df[[c]])
            v[nzchar(v)]
        }),
        use.names = FALSE
    ))
    is_represented <- rep(FALSE, nrow(rep_norm))
    has_any_id <- rep(FALSE, nrow(rep_norm))
    for (col in report_id_cols) {
        v <- as.character(rep_norm[[col]])
        nz <- nzchar(v)
        has_any_id <- has_any_id | nz
        is_represented <- is_represented | (nz & v %in% alias_value_set)
    }
    extra_mask <- !is_represented & has_any_id
    if (any(extra_mask)) {
        rep_extra <- rep_norm[extra_mask, , drop = FALSE]
        extra_df <- as.data.frame(
            matrix("", nrow = nrow(rep_extra), ncol = ncol(alias_df)),
            stringsAsFactors = FALSE
        )
        names(extra_df) <- names(alias_df)
        # Fill the chromAlias-style columns from their report equivalents.
        # Notably, NOT ucsc <- UCSC-style-name: those two columns can use
        # different naming conventions for unplaced scaffolds (chromAlias
        # uses RefSeq-derived "_NW_xxx_random", report uses GenBank-derived
        # "_AABRxxx_random"), and a single column shouldn't carry a mix.
        # The report's UCSC-style-name lives in the separate
        # ucsc_style_name column below.
        safe_syn <- list(
            refseq   = "refseq_accn",
            genbank  = "genbank_accn",
            assembly = "sequence_name"
        )
        for (a_col in names(safe_syn)) {
            r_col <- safe_syn[[a_col]]
            if (a_col %in% names(extra_df) && r_col %in% names(rep_extra)) {
                extra_df[[a_col]] <- as.character(rep_extra[[r_col]])
            }
        }
        for (i in seq_along(new_cols)) {
            out <- out_names[[i]]
            if (out %in% names(extra_df)) {
                extra_df[[out]] <- as.character(rep_extra[[new_cols[[i]]]])
            }
        }
        alias_df <- rbind(alias_df, extra_df)
    }
    alias_df
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

# Assign each alias_df row to a target_chrom name, or to "" if the row can't
# be placed. Pass 1 is name match across every chromAlias column (first
# column whose value is in target_chroms wins for that row). Pass 2 fills
# any remaining row whose length uniquely pairs with a target_chrom on both
# sides via .length_match_fill. Errors if any target_chrom is left unplaced
# -- the contract is "every target chrom maps to exactly one alias row,
# else the requested canonical naming can't be realized."
#
# Used in three places:
#   * .build_seq_ucsc_hub force-align -- rename hub FASTA headers to
#     target_chroms regardless of whether target_chroms appears in any
#     chromAlias column.
#   * .hub_preflight_coverage -- dry-run as a fail-early check before the
#     multi-GB FASTA download.
#   * gdb.install_intervals (called from gdb.build_genome with
#     target_lengths) -- inject the per-row assignment as a virtual
#     ".target_chroms" column so the canonical naming pipeline naturally
#     emits target_chroms.
.assign_target_chroms_per_row <- function(alias_df, target_chroms,
                                          target_lengths, alias_row_lengths) {
    if (!length(target_chroms)) {
        stop(".assign_target_chroms_per_row called with empty target chrom set",
            call. = FALSE
        )
    }
    if (length(target_lengths) != length(target_chroms)) {
        stop(sprintf(
            ".assign_target_chroms_per_row: target_lengths (%d) must align with target_chroms (%d)",
            length(target_lengths), length(target_chroms)
        ), call. = FALSE)
    }
    if (length(alias_row_lengths) != nrow(alias_df)) {
        stop(sprintf(
            ".assign_target_chroms_per_row: alias_row_lengths (%d) must align with alias_df rows (%d)",
            length(alias_row_lengths), nrow(alias_df)
        ), call. = FALSE)
    }
    canonical <- character(nrow(alias_df))
    target_set <- unique(target_chroms)
    # Pass 1: name match across every chromAlias column.
    for (col in names(alias_df)) {
        vals <- as.character(alias_df[[col]])
        vals[is.na(vals)] <- ""
        hits <- !nzchar(canonical) & vals %in% target_set
        canonical[hits] <- vals[hits]
    }
    # Pass 2: unique-on-both-sides length pairing for any remaining empty
    # rows. .length_match_fill returns NA when no unique pair exists; collapse
    # back to "" so the missing-target check below sees a uniform empty marker.
    canonical <- .length_match_fill(
        canonical, alias_row_lengths,
        target_chroms, target_lengths
    )
    canonical[is.na(canonical)] <- ""
    missing <- setdiff(target_chroms, canonical)
    if (length(missing)) {
        stop(sprintf(
            "Cannot align %d of %d target_chroms to chromAlias rows: no name match in any column and no unique-on-both-sides length pair.\nFirst 5 unplaced: %s\nHint: verify target_lengths match the assembly (halStats --sequenceStats), or relax target_chroms to the subset you can guarantee.",
            length(missing), length(target_chroms),
            paste(utils::head(missing, 5L), collapse = ", ")
        ), call. = FALSE)
    }
    canonical
}

# Third-pass refinement of the canonical column: for any row whose canonical
# value isn't in the groot, look across the row's OTHER alias columns for a
# value that IS in the groot, and use that. Catches the case where the
# chosen canonical column's value for a row is missing or follows a
# convention the groot doesn't use, but some other alias column (e.g.
# GenBank-Accn) carries the exact groot chrom name. Never reuses a groot
# chrom already in canonical (no canonical-name collisions).
#
# Returns the canonical vector with overrides applied; caller can compare to
# the input to count how many rows were touched.
.name_match_override <- function(canonical, alias_df, canonical_col, groot_chroms) {
    not_in_groot <- !(canonical %in% groot_chroms)
    if (!any(not_in_groot)) {
        return(canonical)
    }
    for (col in setdiff(names(alias_df), canonical_col)) {
        if (!any(not_in_groot)) break
        vals <- as.character(alias_df[[col]])
        vals[is.na(vals)] <- ""
        can_replace <- not_in_groot & nzchar(vals) & vals %in% groot_chroms &
            !(vals %in% canonical)
        if (any(can_replace)) {
            canonical[can_replace] <- vals[can_replace]
            not_in_groot[can_replace] <- FALSE
        }
    }
    canonical
}

# Format a per-contig diagnostic for groot chroms still unmapped after all
# canonical-resolution passes. Categorizes each unmapped chrom as:
#   * "in alias column(s) X, Y but not canonical" -- a bug or naming edge
#     case (post-pass-3 this is rare).
#   * "length-ambiguous (N alias rows x M groot chroms at length L)" --
#     more than one alias row competes for this groot chrom's length, so
#     no safe pairing is possible.
#   * "no alias row at this length" -- this contig isn't represented in
#     the assembly's chromAlias or assembly_report at all (chromAlias was
#     built from a different assembly version).
# Includes a copy-pasteable line the user can append to chrom_aliases.tsv.
.diagnose_unmapped_chroms <- function(unmapped, alias_df, canonical_col,
                                      groot_chroms, groot_lengths,
                                      alias_row_lengths, groot_path,
                                      max_show = 5L) {
    lines <- character()
    for (chrom in utils::head(unmapped, max_show)) {
        L <- groot_lengths[match(chrom, groot_chroms)]
        L_str <- format(L, big.mark = ",")
        found_cols <- character()
        for (col in names(alias_df)) {
            if (identical(col, canonical_col)) next
            if (chrom %in% alias_df[[col]]) found_cols <- c(found_cols, col)
        }
        if (length(found_cols)) {
            lines <- c(lines, sprintf(
                "    %s (%s bp): present in alias column(s) %s but not in the canonical column.",
                chrom, L_str, paste(found_cols, collapse = ", ")
            ))
        } else if (!is.null(alias_row_lengths)) {
            n_alias_same <- sum(alias_row_lengths == L, na.rm = TRUE)
            n_groot_same <- sum(groot_lengths == L)
            if (n_alias_same == 0L) {
                lines <- c(lines, sprintf(
                    "    %s (%s bp): no alias row at this length -- contig absent from chromAlias and assembly_report. chromAlias was likely built from a different assembly version.",
                    chrom, L_str
                ))
            } else if (n_alias_same > 1L || n_groot_same > 1L) {
                lines <- c(lines, sprintf(
                    "    %s (%s bp): length-ambiguous (%d alias rows x %d groot chroms at this length); length pairing refuses to guess.",
                    chrom, L_str, n_alias_same, n_groot_same
                ))
            } else {
                lines <- c(lines, sprintf(
                    "    %s (%s bp): no name match, no unique length pair.",
                    chrom, L_str
                ))
            }
        } else {
            lines <- c(lines, sprintf(
                "    %s (%s bp): no name match (no chrom.sizes available for length pairing).",
                chrom, L_str
            ))
        }
    }
    sprintf(
        "First %d unmapped:\n%s\nManual fix: append to %s/chrom_aliases.tsv (tab-separated 'canonical<TAB>alias<TAB>source'). For a self-mapping when you only know the groot name, use:\n    %s\t%s\tmanual\nMisha will pick up new rows on the next gdb.init(rescan = TRUE).",
        min(length(unmapped), max_show),
        paste(lines, collapse = "\n"),
        groot_path,
        unmapped[1L], unmapped[1L]
    )
}

# Companion to .length_match_fill for the case where canonical is non-empty
# but the value doesn't actually appear in the groot (e.g. canonical says
# "chrM" but the groot has the MT contig under its GenBank accession
# "AY172581.1"). For each such misaligned row, if the row's length pairs
# uniquely on both sides with a groot chrom that isn't already placed
# anywhere in canonical, REPLACE the canonical value with that groot chrom.
# Strictly conservative: ambiguous lengths are left alone, and groot chroms
# already represented in canonical are never reused (avoids two rows
# colliding on the same canonical name).
#
# canonical, alias_row_lengths must be aligned with alias_df rows.
# groot_chroms, groot_lengths must be aligned with each other.
.length_match_override <- function(canonical, alias_row_lengths,
                                   groot_chroms, groot_lengths) {
    if (length(canonical) != length(alias_row_lengths)) {
        stop(".length_match_override: canonical/alias_row_lengths length mismatch",
            call. = FALSE
        )
    }
    misaligned <- nzchar(canonical) & !is.na(canonical) &
        !(canonical %in% groot_chroms)
    if (!any(misaligned)) {
        return(canonical)
    }
    g_counts <- table(groot_lengths)
    g_unique <- names(g_counts)[g_counts == 1L]
    a_counts <- table(alias_row_lengths)
    a_unique <- names(a_counts)[a_counts == 1L]
    pair_lengths <- intersect(g_unique, a_unique)
    if (!length(pair_lengths)) {
        return(canonical)
    }
    keep <- as.character(groot_lengths) %in% pair_lengths
    groot_lookup <- setNames(groot_chroms[keep], as.character(groot_lengths[keep]))
    # Drop groot chroms already in canonical -- we never reuse one. (Empty
    # canonical cells were filled by .length_match_fill upstream, so any
    # name already there is a "real" placement we shouldn't disturb.)
    groot_lookup <- groot_lookup[!groot_lookup %in% canonical]
    if (!length(groot_lookup)) {
        return(canonical)
    }
    replace_idx <- which(misaligned &
        as.character(alias_row_lengths) %in% names(groot_lookup))
    canonical[replace_idx] <- groot_lookup[as.character(alias_row_lengths[replace_idx])]
    canonical
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
# Returned as a named character vector for O(1) bulk lookup via match()/`[`;
# previously returned an environment, which forced per-element exists()/get().
.build_alias_rev_index <- function(alias_df, canonical_col) {
    cols <- setdiff(names(alias_df), canonical_col)
    canonical <- as.character(alias_df[[canonical_col]])
    if (!length(cols)) {
        return(setNames(character(0), character(0)))
    }
    keys <- unlist(
        lapply(cols, function(col) as.character(alias_df[[col]])),
        use.names = FALSE
    )
    vals <- rep(canonical, times = length(cols))
    keep <- !is.na(keys) & nzchar(keys)
    keys <- keys[keep]
    vals <- vals[keep]
    first <- !duplicated(keys)
    setNames(vals[first], keys[first])
}

# Per-row asset translation. For each value in `rows[[chrom_col]]`, look up
# the canonical name via `rev_idx` (built once per gdb.install_intervals call
# by .build_alias_rev_index). Unmatched chroms become NA -- caller decides
# whether to drop or warn.
.translate_chroms_per_row <- function(rows, chrom_col, rev_idx) {
    asset_names <- as.character(rows[[chrom_col]])
    # Vector indexing handles NA / "" / unmatched names by returning NA, which
    # is exactly the contract callers expect.
    rows[[chrom_col]] <- unname(rev_idx[asset_names])
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
        raw <- utils::read.table(out_path,
            sep = "\t", header = TRUE,
            stringsAsFactors = FALSE, quote = "", comment.char = ""
        )
        # Wide-format (NCBI builds prior to v5.6.30 wrote
        # canonical/refseqAccession/genbankAccession/sequenceName/chrName/role/length).
        # Pivot to long so the merge can dedupe / append uniformly.
        if (!all(c("alias", "source") %in% names(raw)) && "canonical" %in% names(raw)) {
            alias_cols <- setdiff(names(raw), c("canonical", "role", "length"))
            long_parts <- lapply(alias_cols, function(col) {
                vals <- as.character(raw[[col]])
                keep <- nzchar(vals) & !is.na(vals) & vals != raw$canonical
                if (!any(keep)) {
                    return(NULL)
                }
                data.frame(
                    canonical = raw$canonical[keep],
                    alias = vals[keep],
                    source = col,
                    stringsAsFactors = FALSE
                )
            })
            long_parts <- Filter(Negate(is.null), long_parts)
            raw <- if (length(long_parts)) {
                unique(do.call(rbind, long_parts))
            } else {
                data.frame(
                    canonical = character(0), alias = character(0),
                    source = character(0), stringsAsFactors = FALSE
                )
            }
        }
        raw
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
        # Count conflicts (existing alias mapping to a different canonical) for
        # the user-facing warning, then drop every alias already in `existing`.
        # Non-conflict duplicates are dropped too: the existing row wins.
        idx <- match(df$alias, existing$alias)
        is_existing <- !is.na(idx)
        if (any(is_existing)) {
            conflict_mask <- existing$canonical[idx[is_existing]] != df$canonical[is_existing]
            conflicts <- conflicts + sum(conflict_mask, na.rm = TRUE)
            df <- df[!is_existing, , drop = FALSE]
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

# Run the chromAlias coverage gate. Returns the winning column name (with
# attributes from .detect_alias_column) on success; stops with a diagnostic
# matching the existing wording on failure. NULL alias_df is a no-op so
# callers don't need to check.
.coverage_gate <- function(alias_df, target_chroms, target_lengths,
                           min_coverage, label = "groot") {
    if (is.null(alias_df)) {
        return(invisible(NULL))
    }
    col <- .detect_alias_column(alias_df, target_chroms,
        min_coverage = min_coverage, chrom_lengths = target_lengths
    )
    if (!is.na(col)) {
        return(col)
    }
    scores <- attr(col, "scores")
    bp_weighted <- isTRUE(attr(col, "bp_weighted"))
    # Match .detect_alias_column's count-branch denom (unique names); bp branch
    # weights are kept per-row to mirror its bp denom.
    denom <- if (bp_weighted) sum(target_lengths) else length(unique(target_chroms))
    unit <- if (bp_weighted) "bp coverage" else "name coverage"
    stop(sprintf(
        "chromAlias has no column with %.0f%% %s of %s.\nPer-column coverage: %s\nFirst 5 unmapped chroms: %s\nHints:\n  * Lower `min_coverage` to relax (e.g. min_coverage = 0.99).\n  * The registry may resolve this name to a different assembly version than your groot was built from. Try an explicit alternate accession via `source = list(source = \"ucsc-hub\", accession = \"GCF_<old_or_new_version>\")` -- e.g. CanFam3.1 (GCF_000002285.3) vs mCanLor1.2 (GCF_000002285.5) cover the same species but different assemblies.\n  * If the groot was built from de novo contigs (e.g. SPAdes/MEGAHIT k-mer-named) no public chromAlias will match; alias-driven annotation isn't applicable.",
        100 * min_coverage, unit, label,
        paste(sprintf("%s=%.4f%%", names(scores), 100 * scores / denom),
            collapse = ", "
        ),
        paste(utils::head(setdiff(target_chroms, unlist(alias_df, use.names = FALSE)), 5L),
            collapse = ", "
        )
    ), call. = FALSE)
}
