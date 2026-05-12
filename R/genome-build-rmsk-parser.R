# R/genome-build-rmsk-parser.R
# RepeatMasker .out format parser (16-column whitespace, 3 header lines).

# Parses a RepeatMasker .out (or .out.gz) file. Returns a data.frame with
# chrom, start (0-based), end, strand (1L/-1L), name, class, family. Rows
# with < 11 whitespace-delimited fields are dropped with an accumulated
# warning at the end of the file.
.parse_rm_out <- function(file, verbose = TRUE) {
    con <- if (grepl("\\.gz$", file)) gzfile(file, "rt") else file(file, "rt")
    on.exit(close(con), add = TRUE)
    readLines(con, n = 3L) # consume the 3 header lines
    chunks <- list()
    chunk_size <- 50000L
    n_dropped <- 0L
    n_keep_cols <- 11L # only the first 11 fields are used downstream
    repeat {
        lines <- readLines(con, n = chunk_size, warn = FALSE)
        if (!length(lines)) break
        f <- strsplit(trimws(lines), "\\s+", perl = TRUE)
        nf <- lengths(f)
        keep <- nf >= n_keep_cols
        n_dropped <- n_dropped + sum(!keep)
        f <- f[keep]
        if (!length(f)) next
        # Truncate every row to the first 11 fields so they share a length,
        # then matrixify in one pass. ~10x faster than 7 separate vapply calls.
        f <- lapply(f, `[`, seq_len(n_keep_cols))
        chunks[[length(chunks) + 1L]] <- matrix(
            unlist(f, use.names = FALSE),
            ncol = n_keep_cols, byrow = TRUE
        )
    }
    if (n_dropped > 0L && verbose) {
        warning(
            sprintf(
                ".parse_rm_out: dropped %d malformed row(s) with < %d fields",
                n_dropped, n_keep_cols
            ),
            call. = FALSE
        )
    }
    if (!length(chunks)) {
        return(data.frame(
            chrom = character(0), start = integer(0), end = integer(0),
            strand = integer(0), name = character(0),
            class = character(0), family = character(0),
            stringsAsFactors = FALSE
        ))
    }
    mat <- do.call(rbind, chunks)
    cf <- mat[, 11L]
    # Vectorized "class/family" split via regex; avoids per-row strsplit/[.
    slash <- regexpr("/", cf, fixed = TRUE)
    has_slash <- slash > 0L
    class <- ifelse(has_slash, substr(cf, 1L, slash - 1L), cf)
    family <- ifelse(
        has_slash,
        substring(cf, slash + 1L),
        NA_character_
    )
    data.frame(
        chrom = mat[, 5L],
        start = as.integer(mat[, 6L]) - 1L,
        end = as.integer(mat[, 7L]),
        strand = ifelse(mat[, 9L] == "+", 1L, -1L),
        name = mat[, 10L],
        class = class,
        family = family,
        stringsAsFactors = FALSE
    )
}
