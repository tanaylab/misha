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

# gtrack.export_bedgraph with the chromosomes written in the order of chroms (names), genome order if NULL
.gtrack_export_bedgraph <- function(track, file, intervals, iterator, name, chroms = NULL) {
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

    # Extract the data. Overlapping intervals would extract the same bins twice, and bedGraph (and
    # bedGraphToBigWig) does not allow overlapping records, so they are merged first.
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

    # Rows in genome order (chromosome as in the genome database, or as in chroms), then start.
    # NULL when they already are (gextract's output usually is). NaN rows are dropped while writing.
    ord <- NULL
    if (!is.null(data) && nrow(data) > 1) {
        genome_chroms <- if (is.null(chroms)) as.character(gintervals.all()$chrom) else chroms
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
# written in scientific notation, which bedGraphToBigWig parses as 1 -> out-of-order intervals /
# corrupt bigwig.
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
# to con. C++ builds them as R would, in chunks of about 64 MB, unless it does not format
# doubles as this R session does; then R builds them.
.gtrack_write_bedgraph_lines <- function(con, chrom, start, end, value, ord) {
    fmt <- .gbedgraph_format_options()
    if (.gbedgraph_cpp_matches_r(fmt)) {
        from <- 1
        num_rows <- if (is.null(ord)) length(value) else length(ord)
        while (from <= num_rows) {
            res <- .gcall(
                "gtrack_bedgraph_lines", chrom, start, end, as.numeric(value), ord,
                from, 64 * 1024^2, fmt$scipen, fmt$dec, .misha_env()
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
#' Exports a track or track expression to BigWig format by first creating a
#' temporary bedGraph file and then converting it using \code{bedGraphToBigWig}
#' (or \code{wigToBigWig} as a fallback).
#'
#' This function requires the UCSC \code{bedGraphToBigWig} utility to be
#' installed and available on the system PATH, or bundled with the misha
#' package. If not found, the function will raise an error with installation
#' instructions.
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
#'   \code{\link{gtrack.info}}
#' @keywords ~track ~export ~bigwig
#' @examples
#' \dontrun{
#' gdb.init_examples()
#'
#' # Export to BigWig (requires bedGraphToBigWig)
#' gtrack.export_bigwig("dense_track", "/tmp/dense.bw")
#'
#' # With specific region
#' gtrack.export_bigwig("dense_track", "/tmp/dense_chr1.bw",
#'     intervals = gintervals(1, 0, 1e6)
#' )
#' }
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

    # Locate bedGraphToBigWig converter
    converter <- NULL
    converter_name <- "bedGraphToBigWig"

    # 1. Check if bundled with the package
    bundled <- system.file("bedGraphToBigWig", package = "misha")
    if (nzchar(bundled) && file.exists(bundled)) {
        converter <- bundled
    }

    # 2. Check PATH
    if (is.null(converter)) {
        path_result <- Sys.which("bedGraphToBigWig")
        if (nzchar(path_result)) {
            converter <- path_result
        }
    }

    # 3. Fallback to wigToBigWig
    if (is.null(converter)) {
        bundled_wig <- system.file("wigToBigWig", package = "misha")
        if (nzchar(bundled_wig) && file.exists(bundled_wig)) {
            converter <- bundled_wig
            converter_name <- "wigToBigWig"
        } else {
            path_wig <- Sys.which("wigToBigWig")
            if (nzchar(path_wig)) {
                converter <- path_wig
                converter_name <- "wigToBigWig"
            }
        }
    }

    if (is.null(converter)) {
        stop(
            "bedGraphToBigWig or wigToBigWig not found on PATH.\n",
            "Install with: conda install -c bioconda ucsc-bedgraphtobigwig\n",
            "Or download from https://hgdownload.cse.ucsc.edu/admin/exe/ - note that ",
            "UCSC's prebuilt binaries require glibc 2.34 or newer, so on older ",
            "distributions the conda build is the one that works.",
            call. = FALSE
        )
    }

    # Create temporary files
    tmp_bedgraph <- tempfile(fileext = ".bedgraph")
    tmp_chromsizes <- tempfile(fileext = ".chrom.sizes")
    on.exit(
        {
            unlink(tmp_bedgraph)
            unlink(tmp_chromsizes)
        },
        add = TRUE
    )

    # Write bedGraph (uncompressed temp file)
    # bedGraphToBigWig wants the chromosomes sorted by name in the C locale ("10" before "9"),
    # which is not the genome order of a database whose chromosomes are 1, 2, ..., 10
    .gtrack_export_bedgraph(track, tmp_bedgraph, intervals, iterator, NULL,
        chroms = sort(as.character(gintervals.all()$chrom), method = "radix")
    )

    # Write chrom.sizes from the genome database
    # Format sizes as plain integers (no scientific notation) since
    # bedGraphToBigWig requires integer values
    genome_intervals <- gintervals.all()
    chrom_lines <- paste0(
        genome_intervals$chrom, "\t",
        sprintf("%.0f", genome_intervals$end)
    )
    writeLines(chrom_lines, tmp_chromsizes)

    # Run conversion
    result <- system2(converter,
        args = c(tmp_bedgraph, tmp_chromsizes, file),
        stdout = TRUE, stderr = TRUE
    )
    exit_code <- attr(result, "status")

    if (!is.null(exit_code) && exit_code != 0) {
        stop(sprintf(
            "%s failed (exit code %d): %s",
            converter_name, exit_code,
            paste(result, collapse = "\n")
        ), call. = FALSE)
    }

    invisible(NULL)
}
