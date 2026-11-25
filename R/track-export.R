#' Export a track to bedGraph format
#'
#' Exports a track or track expression to a UCSC bedGraph file.
#'
#' This function evaluates a track expression over the specified genomic
#' intervals and writes the result in standard bedGraph format (4-column,
#' tab-separated: chrom, start, end, value). NaN values are omitted from
#' the output.
#'
#' The function supports physical tracks, virtual tracks, and arbitrary track
#' expressions (e.g. \code{"dense_track * 2"}). 2D tracks are not supported.
#'
#' If the output file path ends in \code{.gz}, the output is gzip-compressed.
#'
#' @param track track name or track expression (character string)
#' @param file output file path. If it ends in \code{.gz}, output is
#'   gzip-compressed.
#' @param intervals genomic intervals to export. If \code{NULL} (default),
#'   the entire genome (\code{.misha$ALLGENOME}) is used.
#' @param iterator iterator bin size. If \code{NULL} (default), the iterator
#'   is determined automatically from the track expression.
#' @param name track name for the bedGraph header line. If \code{NULL}
#'   (default), uses the \code{track} parameter value.
#' @return \code{NULL} (invisible). Called for its side effect of writing a
#'   file.
#' @seealso \code{\link{gextract}}, \code{\link{gtrack.export_bigwig}},
#'   \code{\link{gtrack.info}}
#' @keywords ~track ~export ~bedgraph
#' @examples
#' \dontrun{
#' gdb.init_examples()
#'
#' # Export a dense track
#' gtrack.export_bedgraph("dense_track", "/tmp/dense.bedgraph")
#'
#' # Export with specific intervals
#' intervs <- gintervals(1, 0, 1000)
#' gtrack.export_bedgraph("dense_track", "/tmp/dense_chr1.bedgraph",
#'     intervals = intervs
#' )
#'
#' # Export a track expression
#' gtrack.export_bedgraph("dense_track * 2", "/tmp/scaled.bedgraph",
#'     iterator = 100
#' )
#'
#' # Export compressed
#' gtrack.export_bedgraph("dense_track", "/tmp/dense.bedgraph.gz")
#' }
#'
#' @export gtrack.export_bedgraph
gtrack.export_bedgraph <- function(track, file, intervals = NULL, iterator = NULL, name = NULL) {
    if (missing(track) || missing(file)) {
        stop("Usage: gtrack.export_bedgraph(track, file, intervals = NULL, iterator = NULL, name = NULL)",
            call. = FALSE
        )
    }
    .gtrack_export_bedgraph(track, file, intervals, iterator, name)
}

.gtrack_export_bedgraph <- function(track, file, intervals, iterator, name) {
    .gcheckroot()

    if (!is.character(track) || length(track) != 1) {
        stop("'track' must be a single character string (track name or expression)",
            call. = FALSE
        )
    }

    if (!is.character(file) || length(file) != 1) {
        stop("'file' must be a single character string (output file path)",
            call. = FALSE
        )
    }

    # Check for 2D tracks: if the track name matches a known track, verify it's not 2D
    tryCatch(
        {
            info <- gtrack.info(track)
            if (!is.null(info$dimensions) && info$dimensions == 2) {
                stop("2D tracks are not supported by bedGraph export", call. = FALSE)
            }
        },
        error = function(e) {
            # If gtrack.info fails, the track might be a track expression
            # or virtual track, which is fine. Only re-throw if it's our 2D error.
            if (grepl("2D tracks are not supported", e$message)) {
                stop(e$message, call. = FALSE)
            }
        }
    )

    # Check that we can write to the output path
    output_dir <- dirname(file)
    if (!dir.exists(output_dir)) {
        stop(sprintf("Cannot write to '%s': directory does not exist", file),
            call. = FALSE
        )
    }

    if (is.null(name)) {
        name <- track
    }

    # Extract the data. Overlapping intervals would extract the same bins twice, and bedGraph does
    # not allow overlapping records, so they are merged first.
    extract_canonic <- function(intervals) {
        gextract(track, intervals = gintervals.canonic(intervals, unify_touching_intervals = FALSE), iterator = iterator)
    }
    if (is.null(intervals)) {
        data <- gextract(track, intervals = .misha$ALLGENOME, iterator = iterator)
    } else if (is.character(intervals) && length(intervals) == 1 && .gintervals.is_bigset(intervals, FALSE)) {
        # gintervals.canonic cannot read a big intervals set: one chromosome at a time
        data <- do.call(rbind, lapply(as.character(gintervals.chrom_sizes(intervals)$chrom), function(chrom) {
            extract_canonic(gintervals.load(intervals, chrom = chrom))
        }))
    } else {
        data <- extract_canonic(intervals)
    }

    # The value column name is the track expression itself
    value_col <- track

    # gextract returns NULL when the intervals hold no values (e.g. a sparse track)
    if (is.null(data)) {
        warning("The track has no values in the intervals; output file will contain only the header",
            call. = FALSE
        )
    } else if (all(is.nan(data[[value_col]]))) {
        warning("All values are NaN; output file will contain only the header",
            call. = FALSE
        )
    }

    # Rows in genome order (chromosome as in the genome database), then start. NULL when they
    # already are (gextract's output usually is). NaN rows are dropped while writing.
    ord <- NULL
    if (!is.null(data) && nrow(data) > 1) {
        genome_chroms <- as.character(gintervals.all()$chrom)
        chroms <- data$chrom
        chrom_order <- if (is.factor(chroms)) {
            match(levels(chroms), genome_chroms)[as.integer(chroms)]
        } else {
            match(chroms, genome_chroms)
        }
        ord <- order(chrom_order, data$start)
        if (!is.unsorted(ord)) {
            ord <- NULL
        }
    }

    # Write output
    use_gz <- grepl("\\.gz$", file)
    con <- if (use_gz) gzfile(file, "wt") else base::file(file, "wt")
    on.exit(close(con))

    # Write header
    writeLines(sprintf("track type=bedGraph name=\"%s\"", name), con)

    # Write data lines (tab-separated: chrom, start, end, value).
    if (!is.null(data)) {
        .gtrack_write_bedgraph_lines(con, data$chrom, data$start, data$end, data[[value_col]], ord)
    }

    invisible(NULL)
}

# The arguments of as.character() of a double that C++ needs to format as R does: the scipen
# and OutDec options
.gbedgraph_format_options <- function() {
    scipen <- suppressWarnings(as.integer(getOption("scipen"))[1])
    if (is.na(scipen)) {
        scipen <- 0L
    }
    list(scipen = scipen, dec = as.character(getOption("OutDec", "."))[1])
}

# The bedGraph lines of R: paste(chrom, sprintf("%.0f", start), sprintf("%.0f", end), value).
# Coordinates are doubles, formatted as plain integers so large round values (e.g. 1e8) are not
# written in scientific notation, which UCSC's bedGraphToBigWig parses as 1 -> out-of-order
# intervals / corrupt bigwig.
.gbedgraph_lines_r <- function(chrom, start, end, value) {
    paste(chrom, sprintf("%.0f", start), sprintf("%.0f", end), value, sep = "\t")
}

# TRUE if C++ formats doubles as this R session does. It follows R 4.3 and 4.4 built with long
# double, so it is checked against R itself on values that take each branch of R's formatting.
.gbedgraph_cpp_matches_r <- function(fmt) {
    if (!isTRUE(capabilities("long.double"))) {
        return(FALSE)
    }
    probe <- c(
        0, -0, 1, -1.5, 0.1, 1 / 3, -2 / 3, 0.1 + 0.2, 1e5, 1e-5, 1e-4, 1e15, 1e16, 123456.7,
        1234567.1, 100000, 99999.99999999999, 9.9999999999999995e-5, 0.000123456789,
        123456789012345678, 1e-300, 5e-324, .Machine$double.xmax, 3.4e38, pi * 1e10,
        exp(1) * 1e-10, 9.995, 0.15, 2^60, -1e-310, NA, Inf, -Inf
    )
    n <- length(probe)
    got <- .gcall(
        "gtrack_bedgraph_lines", rep("c", n), rep(0, n), rep(1, n), as.numeric(probe), NULL,
        1, 1e6, fmt$scipen, fmt$dec, .misha_env()
    )[[1]]
    identical(got, paste0(.gbedgraph_lines_r(rep("c", n), rep(0, n), rep(1, n), probe), "\n", collapse = ""))
}

# Writes the bedGraph lines of the rows ord (NULL: all of them, in order) whose value is not NaN
# to con. C++ builds them as R would, 2^20 rows at a time, unless it does not format doubles as
# this R session does; then R builds them.
.gtrack_write_bedgraph_lines <- function(con, chrom, start, end, value, ord) {
    fmt <- .gbedgraph_format_options()
    if (.gbedgraph_cpp_matches_r(fmt)) {
        from <- 1
        num_rows <- if (is.null(ord)) length(value) else length(ord)
        while (from <= num_rows) {
            res <- .gcall(
                "gtrack_bedgraph_lines", chrom, start, end, as.numeric(value), ord,
                from, 2^20, fmt$scipen, fmt$dec, .misha_env()
            )
            writeLines(res[[1]], con, sep = "")
            from <- res[[2]]
        }
        return(invisible(NULL))
    }

    rows <- if (is.null(ord)) seq_along(value) else ord
    rows <- rows[!is.nan(value[rows])]
    if (length(rows)) {
        writeLines(.gbedgraph_lines_r(chrom[rows], start[rows], end[rows], value[rows]), con)
    }
    invisible(NULL)
}


#' Export a track to BigWig format
#'
#' Exports a track or track expression to a bigWig file.
#'
#' The track expression is evaluated over the intervals as by \code{\link{gextract}},
#' one batch of chromosomes at a time, and the values are written to the file as they are
#' extracted, so the whole genome's values are never held in memory. Overlapping intervals are
#' merged first, so each base is written once. NaN values are not written. Values are stored as
#' 32-bit floats, as in every bigWig file.
#'
#' The file is written by misha itself; no external tool is needed. Runs of bins of equal size
#' (as from a dense track or a numeric iterator) are stored as fixedStep sections, other values as
#' bedGraph sections. The file has zoom levels (summaries at coarser resolutions), which genome
#' browsers and UCSC's \code{bigWigSummary} use for zoomed-out views. Chromosomes are written in
#' the order of their names (as the bigWig format requires), and every chromosome of the genome
#' is listed in the file, also those without values.
#'
#' The file is written under a temporary name in the same directory and renamed to \code{file}
#' once it is complete, so an export that fails or is interrupted leaves a file already at
#' \code{file} unchanged.
#'
#' @param track track name or track expression (character string)
#' @param file output file path (typically ending in \code{.bw} or
#'   \code{.bigwig}).
#' @param intervals genomic intervals to export. If \code{NULL} (default),
#'   the entire genome (\code{.misha$ALLGENOME}) is used.
#' @param iterator iterator bin size. If \code{NULL} (default), the iterator
#'   is determined automatically from the track expression.
#' @return \code{NULL} (invisible). Called for its side effect of writing a
#'   file.
#' @seealso \code{\link{gextract}}, \code{\link{gtrack.export_bedgraph}},
#'   \code{\link{gtrack.import}}, \code{\link{gtrack.info}}
#' @keywords ~track ~export ~bigwig
#' @examples
#' \dontshow{
#' options(gmax.processes = 2)
#' }
#'
#' gdb.init_examples()
#'
#' f <- tempfile(fileext = ".bw")
#' gtrack.export_bigwig("dense_track", f)
#'
#' # A region, and a track expression in 100 bp bins
#' gtrack.export_bigwig("dense_track * 2", f,
#'     intervals = gintervals(1, 0, 1e5), iterator = 100
#' )
#'
#' # Read it back
#' head(gintervals.import_bigwig(f, gintervals(1, 0, 1000)))
#' unlink(f)
#'
#' @export gtrack.export_bigwig
gtrack.export_bigwig <- function(track, file, intervals = NULL, iterator = NULL) {
    if (missing(track) || missing(file)) {
        stop("Usage: gtrack.export_bigwig(track, file, intervals = NULL, iterator = NULL)",
            call. = FALSE
        )
    }
    .gcheckroot()

    if (!is.character(track) || length(track) != 1) {
        stop("'track' must be a single character string (track name or expression)",
            call. = FALSE
        )
    }

    if (!is.character(file) || length(file) != 1) {
        stop("'file' must be a single character string (output file path)",
            call. = FALSE
        )
    }

    # Check output directory exists
    output_dir <- dirname(file)
    if (!dir.exists(output_dir)) {
        stop(sprintf("Cannot write to '%s': directory does not exist", file),
            call. = FALSE
        )
    }

    # Check for 2D tracks
    tryCatch(
        {
            info <- gtrack.info(track)
            if (!is.null(info$dimensions) && info$dimensions == 2) {
                stop("2D tracks are not supported by BigWig export", call. = FALSE)
            }
        },
        error = function(e) {
            if (grepl("2D tracks are not supported", e$message)) {
                stop(e$message, call. = FALSE)
            }
        }
    )

    # The chromosomes in the order of the file: by name, byte by byte, as the bigWig format's
    # chromosome tree is keyed ("1", "10", "2"), which is not the genome order of a database whose
    # chromosomes are 1, 2, ..., 10. Chromosome id i of the file is file_chroms[i + 1].
    genome <- gintervals.all()
    genome_chroms <- as.character(genome$chrom)
    ord <- order(genome_chroms, method = "radix")
    file_chroms <- genome_chroms[ord]

    groups <- .gbigwig_batches(intervals, iterator, genome[ord, c("chrom", "start", "end")], file_chroms)

    # The file is written under a temporary name in the same directory and renamed over 'file'
    # only once it is complete, so an export that fails (or is interrupted) at any point leaves a
    # file already at the path as it was, and no partial file behind. A symbolic link is written
    # through, to the file it points to, as when the file was written in place.
    dest <- if (nzchar(Sys.readlink(file))) normalizePath(file, mustWork = FALSE) else file
    if (dir.exists(dest)) {
        stop(sprintf("Cannot write to '%s': it is a directory", file), call. = FALSE)
    }
    if (file.access(dirname(dest), 2) != 0 || (file.exists(dest) && file.access(dest, 2) != 0)) {
        stop(sprintf("Cannot write to '%s': permission denied", file), call. = FALSE)
    }
    tmp <- tempfile(pattern = paste0(".", basename(dest), "."), tmpdir = dirname(dest), fileext = ".tmp")
    success <- FALSE
    on.exit(if (!success) unlink(tmp), add = TRUE)

    state <- .gcall("gtrack_bigwig_create", tmp, file_chroms, genome$end[ord], .misha_env())
    for (group in groups) {
        for (intervs in group()) {
            data <- gextract(track, intervals = intervs, iterator = iterator, colnames = "value")
            if (is.null(data) || !nrow(data)) {
                next
            }
            chroms <- data$chrom
            ids <- if (is.factor(chroms)) {
                match(levels(chroms), file_chroms)[as.integer(chroms)]
            } else {
                match(as.character(chroms), file_chroms)
            }
            rows <- order(ids, data$start)
            if (!is.unsorted(rows)) {
                rows <- seq_along(ids)
            }
            state <- .gcall(
                "gtrack_bigwig_append", tmp, state, ids[rows] - 1L, data$start[rows], data$end[rows],
                as.numeric(data$value[rows]), .misha_env()
            )
        }
    }
    num_records <- .gcall("gtrack_bigwig_finish", tmp, state, .misha_env())
    if (file.exists(dest)) {
        Sys.chmod(tmp, file.mode(dest), use_umask = FALSE)
    }
    if (!suppressWarnings(file.rename(tmp, dest))) {
        stop(sprintf("Cannot write to '%s': renaming the finished file into place failed", file), call. = FALSE)
    }
    success <- TRUE

    if (num_records == 0) {
        warning("The track has no values in the intervals; the bigWig file holds no data", call. = FALSE)
    }

    invisible(NULL)
}

# The scope of a bigWig export, in the order of the file (by the rank of the chromosome in
# file_chroms, then start), as a list of functions that each return a list of batches of intervals
# (data frames): one function for an intervals set read one chromosome at a time, one in all
# otherwise. Overlapping intervals are merged. A batch holds whole chromosomes of about 1e8 bases in
# all (an interval is not cut); with a numeric iterator, about 2^22 bins, the intervals being cut at
# multiples of 2^22 bins (bins do not cross them, so the values are those of the whole interval).
.gbigwig_batches <- function(intervals, iterator, genome, file_chroms) {
    batch_bases <- 1e8
    cut <- is.numeric(iterator) && length(iterator) == 1 && is.finite(iterator) && iterator >= 1 && iterator == round(iterator)
    if (cut) {
        batch_bases <- iterator * 2^22
    }

    # the batches of intervals, cut and grouped (NULL or a data frame of any order)
    split_batches <- function(intervs) {
        if (is.null(intervs) || !nrow(intervs)) {
            return(list())
        }
        intervs <- intervs[, c("chrom", "start", "end")]
        intervs <- intervs[order(match(as.character(intervs$chrom), file_chroms), intervs$start), ]
        first <- floor(intervs$start / batch_bases)
        last <- floor((intervs$end - 1) / batch_bases)
        n <- last - first + 1
        if (cut && any(n > 1)) {
            idx <- rep(seq_len(nrow(intervs)), n)
            k <- sequence(n) - 1
            intervs <- data.frame(
                chrom = intervs$chrom[idx],
                start = pmax(intervs$start[idx], (first[idx] + k) * batch_bases),
                end = pmin(intervs$end[idx], (first[idx] + k + 1) * batch_bases)
            )
        }
        widths <- intervs$end - intervs$start
        batch <- floor((cumsum(widths) - widths) / batch_bases)
        lapply(unname(split(seq_len(nrow(intervs)), batch)), function(rows) intervs[rows, ])
    }

    if (is.null(intervals)) {
        return(list(function() split_batches(genome)))
    }

    if (is.character(intervals) && length(intervals) == 1 && .gintervals.is_bigset(intervals, FALSE)) {
        # gintervals.canonic cannot read a big intervals set: one chromosome at a time
        chroms <- as.character(gintervals.chrom_sizes(intervals)$chrom)
        chroms <- chroms[order(match(chroms, file_chroms))]
        return(lapply(chroms, function(chrom) {
            force(chrom)
            function() {
                split_batches(gintervals.canonic(gintervals.load(intervals, chrom = chrom), unify_touching_intervals = FALSE))
            }
        }))
    }

    canonic <- gintervals.canonic(intervals, unify_touching_intervals = FALSE)
    list(function() split_batches(canonic))
}
