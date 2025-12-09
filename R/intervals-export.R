#' Export Intervals to BED Format
#'
#' @description
#' Exports a misha intervals set to BED format. Automatically detects whether to
#' export as BED3 (basic) or BED6 (with strand) based on available columns.
#'
#' @param intervals Either a data.frame with genomic intervals (must have chrom, start, end columns)
#'   or a character string naming a saved intervals set (bigset) to load and export.
#' @param file Output BED file path. Should end in .bed.
#' @param track_name Optional track name for BED header line. If provided, a track definition
#'   line will be added at the start of the file.
#' @param description Optional description for BED header line. Only used if \code{track_name}
#'   is also provided.
#'
#' @return Invisibly returns the file path (useful for piping).
#'
#' @details
#' This function exports genomic intervals (without track values) to BED format. The output
#' format is automatically determined:
#' \itemize{
#'   \item BED3: If only chrom, start, end are available (3 columns)
#'   \item BED6: If strand column is present (6 columns: chrom, start, end, name, score, strand)
#' }
#'
#' For BED6 format, the name field is set to "." and score to 0 as placeholders, since
#' intervals sets don't have these values.
#'
#' **Coordinate System:** Both misha and BED format use 0-based half-open intervals [start, end),
#' so coordinates are written directly without conversion.
#'
#' If \code{intervals} is a character string, the function will attempt to load it as a
#' saved intervals set using \code{\link{gintervals.load}}.
#'
#' @examples
#' \dontrun{
#' # Export intervals data.frame to BED3
#' my_intervals <- gintervals(c(1, 1, 2), c(1000, 5000, 2000), c(2000, 6000, 3000))
#' gintervals.export_bed(my_intervals, "output.bed")
#'
#' # Export with track header
#' gintervals.export_bed(my_intervals, "output.bed",
#'     track_name = "MyRegions",
#'     description = "Interesting genomic regions"
#' )
#'
#' # Export intervals with strand information (BED6)
#' stranded_intervals <- gintervals(c(1, 1, 2), c(1000, 5000, 2000), c(2000, 6000, 3000))
#' stranded_intervals$strand <- c(1, -1, 1)
#' gintervals.export_bed(stranded_intervals, "stranded.bed")
#'
#' # Export saved intervals set by name
#' gintervals.save("my_regions", my_intervals)
#' gintervals.export_bed("my_regions", "regions.bed")
#' }
#'
#' @seealso \code{\link{gtrack.import}} for importing BED files as tracks,
#'   \code{\link{gintervals.load}} for loading saved intervals sets,
#'   \code{\link{gintervals.save}} for saving intervals sets
#'
#' @export
gintervals.export_bed <- function(intervals, file, track_name = NULL, description = NULL) {
    if (!is.character(file) || length(file) != 1) {
        stop("'file' must be a single character string", call. = FALSE)
    }

    file_dir <- dirname(file)
    if (!dir.exists(file_dir)) {
        stop("Directory does not exist: ", file_dir, call. = FALSE)
    }

    if (is.character(intervals) && length(intervals) == 1) {
        if (gintervals.exists(intervals)) {
            intervals <- gintervals.load(intervals)
        } else {
            stop("Intervals set '", intervals, "' does not exist", call. = FALSE)
        }
    }

    if (!is.data.frame(intervals)) {
        stop("'intervals' must be a data.frame or a saved intervals set name", call. = FALSE)
    }

    if ("chrom1" %in% names(intervals) || "chrom2" %in% names(intervals)) {
        stop("Cannot export 2D intervals to BED format", call. = FALSE)
    }

    required_cols <- c("chrom", "start", "end")
    if (!all(required_cols %in% names(intervals))) {
        missing <- setdiff(required_cols, names(intervals))
        stop("Intervals must have chrom, start, and end columns. Missing: ",
            paste(missing, collapse = ", "),
            call. = FALSE
        )
    }

    if (nrow(intervals) == 0) {
        stop("Cannot export empty intervals set", call. = FALSE)
    }

    # BED uses 0-based coordinates (same as misha), no conversion needed
    bed_data <- data.frame(
        chrom = as.character(intervals$chrom),
        start = intervals$start,
        end = intervals$end,
        stringsAsFactors = FALSE
    )

    if ("strand" %in% names(intervals)) {
        # BED6 requires name and score columns before strand
        bed_data$name <- "."
        bed_data$score <- 0
        bed_data$strand <- ifelse(intervals$strand == 1, "+", ifelse(intervals$strand == -1, "-", "."))
    }

    if (!is.null(track_name)) {
        if (!is.character(track_name) || length(track_name) != 1) {
            stop("'track_name' must be a single character string", call. = FALSE)
        }

        header <- sprintf('track name="%s"', track_name)

        if (!is.null(description)) {
            if (!is.character(description) || length(description) != 1) {
                stop("'description' must be a single character string", call. = FALSE)
            }
            header <- sprintf('%s description="%s"', header, description)
        }

        writeLines(header, file)
        append <- TRUE
    } else {
        append <- FALSE
    }

    if (isTRUE(getOption("misha.use_readr", TRUE)) && .misha_readr_available()) {
        .misha_write_with_readr(
            bed_data,
            file,
            append = append,
            col_names = FALSE,
            na = "."
        )
    } else {
        write.table(bed_data,
            file,
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE,
            append = append
        )
    }

    invisible(file)
}
