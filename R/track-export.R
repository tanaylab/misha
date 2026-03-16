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

    # Set up intervals
    if (is.null(intervals)) {
        intervals <- .misha$ALLGENOME
    }

    # Extract data using gextract
    data <- gextract(track, intervals = intervals, iterator = iterator)

    # The value column name is the track expression itself
    value_col <- track

    # Remove rows with NaN values
    non_nan <- !is.nan(data[[value_col]])
    data <- data[non_nan, , drop = FALSE]

    if (nrow(data) == 0) {
        warning("All values are NaN; output file will contain only the header",
            call. = FALSE
        )
    }

    # Sort by chromosome (genome order) then by start
    # Get chromosome order from the genome database
    genome_chroms <- gintervals.all()$chrom
    chrom_order <- match(data$chrom, genome_chroms)
    ord <- order(chrom_order, data$start)
    data <- data[ord, , drop = FALSE]

    # Write output
    use_gz <- grepl("\\.gz$", file)
    con <- if (use_gz) gzfile(file, "wt") else base::file(file, "wt")
    on.exit(close(con))

    # Write header
    writeLines(sprintf("track type=bedGraph name=\"%s\"", name), con)

    # Write data lines (tab-separated: chrom, start, end, value)
    if (nrow(data) > 0) {
        lines <- paste(data$chrom, data$start, data$end, data[[value_col]], sep = "\t")
        writeLines(lines, con)
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
            "bedGraphToBigWig or wigToBigWig not found. ",
            "Install from UCSC tools: https://hgdownload.cse.ucsc.edu/admin/exe/",
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
    gtrack.export_bedgraph(track, tmp_bedgraph,
        intervals = intervals,
        iterator = iterator
    )

    # Write chrom.sizes from the genome database
    # Format sizes as plain integers (no scientific notation) since
    # bedGraphToBigWig requires integer values
    genome_intervals <- gintervals.all()
    chrom_lines <- paste0(
        genome_intervals$chrom, "\t",
        formatC(genome_intervals$end, format = "d")
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
