#' Export Track to BigWig Format
#'
#' @description
#' Exports a misha track to BigWig format. Works with dense tracks, sparse tracks,
#' virtual tracks, and track expressions.
#'
#' @param track Track name or track expression to export. Can be a regular track name,
#'   a virtual track (e.g., "gc", "pwm.max"), or a track expression (e.g., "track1 + track2").
#' @param file Output BigWig file path. Should end in .bw or .bigWig.
#' @param intervals Genomic intervals to export. Default is ALLGENOME (entire genome).
#' @param iterator Iterator for computing track values. Can be:
#'   \itemize{
#'     \item Numeric: bin size in base pairs (e.g., \code{iterator = 1000} for 1kb bins)
#'     \item Character: track name to use as iterator intervals
#'     \item data.frame: explicit intervals to use as iterator
#'     \item NULL (default): use track's native resolution
#'   }
#' @param chrom_sizes Named numeric vector of chromosome sizes. If NULL (default),
#'   chromosome sizes are automatically extracted from ALLGENOME.
#' @param na_value How to handle NaN/NA values (BigWig format doesn't support them):
#'   \itemize{
#'     \item NULL (default): Remove intervals with NaN/NA values
#'     \item numeric: Replace NaN/NA values with this value (e.g., \code{na_value = 0})
#'   }
#' @param fixed_summaries Logical. If TRUE, compute summaries at fixed Ensembl-style zoom
#'   levels (30X, 65X, 130X, 260X, 450X, 648X, 950X, 1296X, 4800X, 19200X). If FALSE (default),
#'   zoom levels are dynamically determined based on data characteristics. Fixed summaries
#'   may be preferable for consistent visualization across different genome browsers.
#'
#' @return Invisibly returns the file path (useful for piping).
#'
#' @details
#' BigWig format does not support missing values (NaN/NA). By default (\code{na_value=NULL}),
#' intervals with NaN values are removed from the export. Alternatively, set \code{na_value}
#' to a numeric value (e.g., 0) to replace NaN with that value instead of removing the intervals.
#'
#' This function requires the \pkg{rtracklayer} package from Bioconductor. Install it with:
#' \code{BiocManager::install("rtracklayer")}
#'
#' The function uses \code{\link{gextract}} internally, which provides a unified interface
#' for extracting data from all track types (dense, sparse, virtual, expressions). The
#' \code{iterator} parameter is particularly useful for:
#' \itemize{
#'   \item Reducing output file size by binning high-resolution tracks
#'   \item Computing track expressions in fixed-size windows
#'   \item Exporting virtual tracks (which require an iterator to be evaluated)
#' }
#'
#' @examples
#' \dontrun{
#' # Export entire track at native resolution
#' gtrack.export_bw("my.track", "output.bw")
#'
#' # Export track binned at 1kb resolution
#' gtrack.export_bw("my.track", "output_1kb.bw", iterator = 1000)
#'
#' # Export virtual track (GC content in 500bp windows)
#' gtrack.export_bw("gc", "gc_content.bw", iterator = 500)
#'
#' # Export track expression
#' gtrack.export_bw("track1 + track2", "combined.bw", iterator = 1000)
#'
#' # Export specific chromosomes
#' gtrack.export_bw("my.track", "chr1.bw", intervals = gintervals(1))
#'
#' # Replace NaN values with 0 instead of removing them
#' gtrack.export_bw("coverage.track", "coverage.bw", na_value = 0, iterator = 100)
#'
#' # Use fixed Ensembl-style zoom levels for consistent browser visualization
#' gtrack.export_bw("my.track", "output.bw", fixed_summaries = TRUE)
#' }
#'
#' @seealso \code{\link{gtrack.import}} for importing tracks including BigWig files,
#'   \code{\link{gextract}} for extracting track data
#'
#' @export
gtrack.export_bw <- function(track, file, intervals = get("ALLGENOME", envir = .misha), iterator = NULL,
                             chrom_sizes = NULL, na_value = NULL, fixed_summaries = FALSE) {
    required_pkgs <- c("rtracklayer", "GenomicRanges", "IRanges", "GenomeInfoDb")
    missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]

    if (length(missing_pkgs) > 0) {
        stop("The following Bioconductor packages are required for BigWig export:\n",
            paste(missing_pkgs, collapse = ", "), "\n",
            "Install them with: BiocManager::install(c('",
            paste(missing_pkgs, collapse = "', '"), "'))",
            call. = FALSE
        )
    }

    if (!is.null(na_value) && (!is.numeric(na_value) || length(na_value) != 1)) {
        stop("Parameter 'na_value' must be NULL or a single numeric value", call. = FALSE)
    }

    if (!is.logical(fixed_summaries) || length(fixed_summaries) != 1) {
        stop("Parameter 'fixed_summaries' must be TRUE or FALSE", call. = FALSE)
    }

    file_dir <- dirname(file)
    if (!dir.exists(file_dir)) {
        stop("Directory does not exist: ", file_dir, call. = FALSE)
    }

    if (!is.null(chrom_sizes)) {
        if (!is.numeric(chrom_sizes) || is.null(names(chrom_sizes))) {
            stop("Parameter 'chrom_sizes' must be a named numeric vector", call. = FALSE)
        }
    } else {
        allgenome <- get("ALLGENOME", envir = .misha)
        chrom_sizes <- setNames(allgenome$end, allgenome$chrom)
    }

    # Extract track data using gextract - handles all track types (dense, sparse, virtual, expressions) and iterator
    data <- tryCatch(
        {
            gextract(track, intervals = intervals, iterator = iterator)
        },
        error = function(e) {
            stop("Failed to extract track data: ", e$message, call. = FALSE)
        }
    )

    if (is.null(data) || nrow(data) == 0) {
        stop("No data extracted from track. Check track name and intervals.", call. = FALSE)
    }

    if ("chrom1" %in% names(data) || "chrom2" %in% names(data)) {
        stop("Cannot export 2D tracks to BigWig format (use gextract with 1D projection)", call. = FALSE)
    }

    if (!all(c("chrom", "start", "end") %in% names(data))) {
        stop("Extracted data is missing required columns (chrom, start, end)", call. = FALSE)
    }

    standard_cols <- c("chrom", "start", "end", "intervalID")
    value_cols <- setdiff(names(data), standard_cols)

    if (length(value_cols) == 0) {
        stop("No value column found in extracted data", call. = FALSE)
    }

    value_col <- value_cols[1]

    values <- data[[value_col]]
    na_mask <- is.na(values)
    n_na <- sum(na_mask)

    if (n_na > 0) {
        if (is.null(na_value)) {
            pct_na <- 100 * n_na / length(values)
            if (pct_na > 10) {
                warning(sprintf(
                    "Removing %d intervals with NaN values (%.1f%% of data). Consider setting na_value=0 to replace NaN instead.",
                    n_na, pct_na
                ), call. = FALSE)
            }
            data <- data[!na_mask, ]
        } else {
            data[[value_col]][na_mask] <- na_value
        }
        values <- data[[value_col]]
    }

    if (nrow(data) == 0) {
        if (n_na > 0 && is.null(na_value)) {
            stop("All values are NaN. Set na_value=<number> to replace NaN or check your track/expression.", call. = FALSE)
        } else {
            stop("No valid (non-NaN) values to export. Consider setting na_value=0 to replace NaN instead of removing.", call. = FALSE)
        }
    }

    # Coordinate conversion: misha uses 0-based half-open [start, end), GRanges uses 1-based closed [start, end]
    # Conversion: GRanges_start = misha_start + 1, GRanges_end = misha_end
    granges <- tryCatch(
        {
            data_chroms <- unique(as.character(data$chrom))

            match_idx <- match(data_chroms, names(chrom_sizes))
            seq_lengths <- ifelse(!is.na(match_idx), chrom_sizes[match_idx], NA_real_)
            names(seq_lengths) <- data_chroms

            if (any(is.na(seq_lengths))) {
                max_ends <- aggregate(data$end, by = list(chrom = as.character(data$chrom)), FUN = max)
                max_ends_vec <- setNames(max_ends$x, max_ends$chrom)
                seq_lengths[is.na(seq_lengths)] <- max_ends_vec[names(seq_lengths)[is.na(seq_lengths)]]
            }

            seqinfo <- GenomeInfoDb::Seqinfo(
                seqnames = data_chroms,
                seqlengths = as.integer(seq_lengths)
            )

            gr <- GenomicRanges::GRanges(
                seqnames = as.character(data$chrom),
                ranges = IRanges::IRanges(
                    start = data$start + 1, # Convert from 0-based to 1-based
                    end = data$end
                ),
                score = values,
                seqinfo = seqinfo
            )

            gr
        },
        error = function(e) {
            stop("Failed to create GRanges object: ", e$message, call. = FALSE)
        }
    )

    tryCatch(
        {
            rtracklayer::export.bw(granges, file, fixedSummaries = fixed_summaries)
        },
        error = function(e) {
            stop("Failed to write BigWig file: ", e$message, call. = FALSE)
        }
    )

    if (!file.exists(file)) {
        stop("BigWig file was not created: ", file, call. = FALSE)
    }

    invisible(file)
}

#' Export Track to BED Format
#'
#' @description
#' Exports a misha track to BED format. Works with dense tracks, sparse tracks,
#' virtual tracks, and track expressions. Output format (BED3/BED5/BED6) is
#' automatically detected based on available data.
#'
#' @param track Track name or track expression to export. Can be a regular track name,
#'   a virtual track (e.g., "gc", "kmer.count"), or a track expression (e.g., "track1 + track2").
#' @param file Output BED file path. Should end in .bed.
#' @param intervals Genomic intervals to export. Default is ALLGENOME (entire genome).
#' @param iterator Iterator for computing track values. Can be:
#'   \itemize{
#'     \item Numeric: bin size in base pairs (e.g., \code{iterator = 1000} for 1kb bins)
#'     \item Character: track name to use as iterator intervals
#'     \item data.frame: explicit intervals to use as iterator
#'     \item NULL (default): use track's native resolution
#'   }
#' @param na_value How to handle NaN/NA values in the score column:
#'   \itemize{
#'     \item NULL (default): Write "." for NaN/NA values (standard BED convention)
#'     \item numeric: Replace NaN/NA values with this value (e.g., \code{na_value = 0})
#'   }
#' @param track_name Optional track name for BED header line. If provided, a track definition
#'   line will be added at the start of the file.
#' @param description Optional description for BED header line. Only used if \code{track_name}
#'   is also provided.
#'
#' @return Invisibly returns the file path (useful for piping).
#'
#' @details
#' The output format is automatically determined based on the extracted data:
#' \itemize{
#'   \item BED3: Only chrom, start, end (no track values extracted)
#'   \item BED5: Track values present → chrom, start, end, name, score
#'   \item BED6: Track values and strand present → chrom, start, end, name, score, strand
#' }
#'
#' For BED5/BED6 formats, the name field is auto-generated as "interval_N".
#'
#' **Coordinate System:** Both misha and BED format use 0-based half-open intervals [start, end),
#' so coordinates are written directly without conversion. This is different from BigWig export
#' which requires start+1 conversion.
#'
#' The function uses \code{\link{gextract}} internally, which provides a unified interface
#' for extracting data from all track types. The \code{iterator} parameter is particularly
#' useful for:
#' \itemize{
#'   \item Binning high-resolution tracks into fixed-size windows
#'   \item Computing track expressions in specific genomic regions
#'   \item Evaluating virtual tracks (which require an iterator)
#' }
#'
#' @examples
#' \dontrun{
#' # Export sparse track at native resolution
#' gtrack.export_bed("my.track", "output.bed")
#'
#' # Export track binned at 1kb resolution
#' gtrack.export_bed("my.track", "output_1kb.bed", iterator = 1000)
#'
#' # Export virtual track (GC content in 500bp windows)
#' gtrack.export_bed("gc", "gc_content.bed", iterator = 500)
#'
#' # Export track expression
#' gtrack.export_bed("track1 + track2", "combined.bed", iterator = 1000)
#'
#' # Export specific chromosomes with track header
#' gtrack.export_bed("my.track", "chr1.bed",
#'     intervals = gintervals(1),
#'     track_name = "MyTrack",
#'     description = "Chromosome 1 data"
#' )
#'
#' # Replace NaN values with 0
#' gtrack.export_bed("coverage.track", "coverage.bed", na_value = 0, iterator = 100)
#' }
#'
#' @seealso \code{\link{gtrack.export_bw}} for BigWig export,
#'   \code{\link{gtrack.import}} for importing tracks including BED files,
#'   \code{\link{gextract}} for extracting track data,
#'   \code{\link{gintervals.export_bed}} for exporting intervals without track values
#'
#' @export
gtrack.export_bed <- function(track, file, intervals = get("ALLGENOME", envir = .misha),
                              iterator = NULL, na_value = NULL,
                              track_name = NULL, description = NULL) {
    if (!is.character(file) || length(file) != 1) {
        stop("'file' must be a single character string", call. = FALSE)
    }

    file_dir <- dirname(file)
    if (!dir.exists(file_dir)) {
        stop("Directory does not exist: ", file_dir, call. = FALSE)
    }

    if (!is.null(na_value)) {
        if (!is.numeric(na_value) || length(na_value) != 1 || is.na(na_value)) {
            stop("'na_value' must be a single numeric value", call. = FALSE)
        }
    }

    if (!is.null(track_name) && (!is.character(track_name) || length(track_name) != 1)) {
        stop("'track_name' must be a single character string", call. = FALSE)
    }

    if (!is.null(description) && (!is.character(description) || length(description) != 1)) {
        stop("'description' must be a single character string", call. = FALSE)
    }

    data <- tryCatch(
        {
            gextract(track, intervals = intervals, iterator = iterator)
        },
        error = function(e) {
            stop("Failed to extract track data: ", e$message, call. = FALSE)
        }
    )

    if (is.null(data) || nrow(data) == 0) {
        stop("No data extracted from track", call. = FALSE)
    }

    if ("chrom1" %in% names(data) || "chrom2" %in% names(data)) {
        stop("Cannot export 2D tracks to BED format", call. = FALSE)
    }

    standard_cols <- c("chrom", "start", "end", "intervalID")
    value_cols <- setdiff(names(data), standard_cols)

    # BED uses 0-based coordinates (same as misha), no conversion needed
    bed_data <- data.frame(
        chrom = as.character(data$chrom),
        start = data$start,
        end = data$end,
        stringsAsFactors = FALSE
    )

    if (length(value_cols) > 0) {
        values <- data[[value_cols[1]]]

        if (!is.null(na_value)) {
            values[is.na(values)] <- na_value
        }

        bed_data$name <- paste0("interval_", seq_len(nrow(data)))
        bed_data$score <- values

        if ("strand" %in% names(data)) {
            bed_data$strand <- ifelse(data$strand == 1, "+",
                ifelse(data$strand == -1, "-", ".")
            )
        }
    }

    if (!is.null(track_name)) {
        header <- sprintf('track name="%s"', track_name)

        if (!is.null(description)) {
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
        # Use na = "." to write NaN/NA as "." (BED standard)
        write.table(bed_data,
            file,
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE,
            append = append,
            na = "."
        )
    }

    invisible(file)
}
