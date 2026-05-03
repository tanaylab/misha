# R/genome-build-chromalias.R
# Source-agnostic chromosome-alias parsing, column detection, and translation.

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

# Returns the (unique) column whose values contain *every* target_chrom.
# Ties broken by column order. NA with attributes(scores=, overlap=) when no
# column achieves 100% coverage — caller produces the diagnostic.
.detect_alias_column <- function(alias_df, target_chroms) {
    if (!length(target_chroms)) {
        stop(".detect_alias_column called with empty target chrom set", call. = FALSE)
    }
    scores <- vapply(
        alias_df,
        function(col) length(intersect(col, target_chroms)),
        integer(1)
    )
    full <- names(alias_df)[scores == length(target_chroms)]
    if (!length(full)) {
        return(structure(NA_character_,
            scores = scores,
            overlap = max(scores) / length(target_chroms)
        ))
    }
    full[[1L]]
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
.merge_chrom_aliases_tsv <- function(groot, alias_df, groot_col, source_label = "build") {
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
