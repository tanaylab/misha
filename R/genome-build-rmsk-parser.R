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
    rows <- list()
    chunk_size <- 50000L
    n_dropped <- 0L
    repeat {
        lines <- readLines(con, n = chunk_size, warn = FALSE)
        if (!length(lines)) break
        f <- strsplit(trimws(lines), "\\s+", perl = TRUE)
        nf <- vapply(f, length, integer(1))
        keep <- nf >= 11L
        n_dropped <- n_dropped + sum(!keep)
        f <- f[keep]
        if (!length(f)) next
        rows[[length(rows) + 1L]] <- data.frame(
            chrom = vapply(f, `[`, character(1), 5L),
            start = as.integer(vapply(f, `[`, character(1), 6L)) - 1L,
            end = as.integer(vapply(f, `[`, character(1), 7L)),
            strand = ifelse(vapply(f, `[`, character(1), 9L) == "+", 1L, -1L),
            name = vapply(f, `[`, character(1), 10L),
            cf = vapply(f, `[`, character(1), 11L),
            stringsAsFactors = FALSE
        )
    }
    if (n_dropped > 0L && verbose) {
        warning(
            sprintf(
                ".parse_rm_out: dropped %d malformed row(s) with < 11 fields",
                n_dropped
            ),
            call. = FALSE
        )
    }
    if (!length(rows)) {
        return(data.frame(
            chrom = character(0), start = integer(0), end = integer(0),
            strand = integer(0), name = character(0),
            class = character(0), family = character(0),
            stringsAsFactors = FALSE
        ))
    }
    df <- do.call(rbind, rows)
    cf_split <- strsplit(df$cf, "/", fixed = TRUE)
    df$class <- vapply(cf_split, `[`, character(1), 1L)
    df$family <- vapply(
        cf_split,
        function(x) if (length(x) >= 2L) x[[2L]] else NA_character_,
        character(1)
    )
    df$cf <- NULL
    df
}
