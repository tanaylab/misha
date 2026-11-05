# Interval annotation functions

#' Annotates 1D intervals using nearest neighbors
#'
#' Annotates one-dimensional intervals by finding nearest neighbors in another
#' set of intervals and adding selected columns from the neighbors to the
#' original intervals.
#'
#' The function wraps and extends \code{gintervals.neighbors} to provide
#' convenient column selection/renaming, optional distance inclusion, distance
#' thresholding with custom NA values, multiple neighbors per interval, and
#' deterministic tie-breaking. Currently supports 1D intervals only.
#'
#' @param intervals Intervals to annotate (1D).
#' @param annotation_intervals Source intervals containing annotation data (1D).
#' @param annotation_columns Character vector of column names to copy from
#'   \code{annotation_intervals}. If \code{NULL} (default), all non-basic
#'   columns are used, i.e. everything beyond the coordinate/strand columns
#'   among: chrom, start, end, chrom1, start1, end1, chrom2, start2, end2, strand.
#' @param column_names Optional custom names for the annotation columns. If
#'   provided, must have the same length as \code{annotation_columns}. Defaults
#'   to using the original names.
#' @param dist_column Name of the distance column to include. Use \code{NULL} to
#'   omit the distance column. Defaults to "dist".
#' @param max_dist Maximum absolute distance. When finite, neighbors with
#'   \code{|dist| > max_dist} result in annotation columns being set to
#'   \code{na_value} for those rows, while the row itself is retained.
#' @param na_value Value(s) to use for annotations when beyond \code{max_dist}
#'   or when no neighbor is found. Can be a single scalar recycled for all
#'   columns, or a named list/vector supplying per-column values matching
#'   \code{column_names}.
#' @param maxneighbors Maximum number of neighbors per interval (duplicates
#'   intervals as needed). Defaults to 1.
#' @param tie_method Tie-breaking when distances are equal: one of
#'   "first" (arbitrary but stable), "min.start" (smaller neighbor start first),
#'   or "min.end" (smaller neighbor end first). Applies when
#'   \code{maxneighbors > 1}.
#' @param overwrite When \code{FALSE} (default), errors if selected annotation
#'   columns would overwrite existing columns in \code{intervals}. When
#'   \code{TRUE}, conflicting base columns are replaced by the annotation
#'   columns.
#' @param keep_order If \code{TRUE} (default), preserves the original order of
#'   \code{intervals} rows in the output.
#' @param intervals.set.out intervals set name where the function result is
#' optionally outputted
#' @param ... Additional arguments forwarded to \code{gintervals.neighbors}
#'   (e.g., \code{mindist}, \code{maxdist}).
#'
#' @return A data frame containing the original intervals plus the requested
#'   annotation columns (and optional distance column). If
#'   \code{maxneighbors > 1}, rows may be duplicated per input interval to
#'   accommodate multiple neighbors.
#'
#' @details
#' - When \code{annotation_columns = NULL}, all non-basic columns present in
#'   \code{annotation_intervals} are included.
#' - Setting \code{dist_column = NULL} omits the distance column.
#' - If no neighbor is found for an interval, annotation columns are filled with
#'   \code{na_value} and the distance (when present) is \code{NA_real_}.
#' - Column name collisions are handled as follows: when \code{overwrite=FALSE}
#'   a clear error is emitted; when \code{overwrite=TRUE}, base columns with the
#'   same names are replaced by annotation columns.
#'
#' @examples
#' # Prepare toy data
#' intervs <- gintervals(1, c(1000, 5000), c(1100, 5050))
#' ann <- gintervals(1, c(900, 5400), c(950, 5500))
#' ann$remark <- c("a", "b")
#' ann$score <- c(10, 20)
#'
#' # Basic usage with default columns (all non-basic columns)
#' gintervals.annotate(intervs, ann)
#'
#' # Select specific columns, with custom names and distance column name
#' gintervals.annotate(
#'     intervs, ann,
#'     annotation_columns = c("remark"),
#'     column_names = c("ann_remark"),
#'     dist_column = "ann_dist"
#' )
#'
#' # Distance threshold with scalar NA replacement
#' gintervals.annotate(
#'     intervs, ann,
#'     annotation_columns = c("remark"),
#'     max_dist = 200,
#'     na_value = "no_ann"
#' )
#'
#' # Multiple neighbors with deterministic tie-breaking
#' nbrs <- gintervals.annotate(
#'     gintervals(1, 1000, 1100),
#'     {
#'         x <- gintervals(1, c(800, 1200), c(900, 1300))
#'         x$label <- c("left", "right")
#'         x
#'     },
#'     annotation_columns = "label",
#'     maxneighbors = 2,
#'     tie_method = "min.start"
#' )
#' nbrs
#'
#' # Overwrite existing columns in the base intervals
#' intervs2 <- intervs
#' intervs2$remark <- c("orig1", "orig2")
#' gintervals.annotate(intervs2, ann, annotation_columns = "remark", overwrite = TRUE)
#' @export
gintervals.annotate <- function(intervals,
                                annotation_intervals,
                                annotation_columns = NULL,
                                column_names = NULL,
                                dist_column = "dist",
                                max_dist = Inf,
                                na_value = NA,
                                maxneighbors = 1,
                                tie_method = c("first", "min.start", "min.end"),
                                overwrite = FALSE,
                                keep_order = TRUE,
                                intervals.set.out = NULL,
                                ...) {
    # Input validation
    if (is.null(intervals) || is.null(annotation_intervals)) {
        stop("Usage: gintervals.annotate(intervals, annotation_intervals, ...)", call. = FALSE)
    }

    tie_method <- match.arg(tie_method)

    # Normalize intervals.set.out to a name if provided
    intervals.set.out <- do.call(.gexpr2str, list(substitute(intervals.set.out)), envir = parent.frame())
    if (!is.null(intervals.set.out)) {
        fullpath <- .gintervals.check_new_set(intervals.set.out)
    }

    # Load intervals if needed
    intervals <- .gintervals.load_ext(intervals)
    annotation_intervals <- .gintervals.load_ext(annotation_intervals)

    # Ties are broken by start or end, which 2D intervals do not have
    if (tie_method != "first" && maxneighbors > 1 && !.gintervals.is1d(annotation_intervals)) {
        stop("tie_method other than \"first\" is supported for 1D intervals only", call. = FALSE)
    }

    # Original row order: used to group ties and, with keep_order, to restore the order
    intervals$.orig_order <- seq_len(nrow(intervals))

    # Determine annotation columns if not specified
    if (is.null(annotation_columns)) {
        basic_cols <- c(
            "chrom", "start", "end", "chrom1", "start1", "end1",
            "chrom2", "start2", "end2", "strand"
        )
        annotation_columns <- setdiff(colnames(annotation_intervals), basic_cols)
    }

    # Validate annotation columns exist
    missing_cols <- setdiff(annotation_columns, colnames(annotation_intervals))
    if (length(missing_cols) > 0) {
        stop(paste(
            "Annotation columns not found in annotation_intervals:",
            paste(missing_cols, collapse = ", ")
        ), call. = FALSE)
    }

    # Set up column names
    if (is.null(column_names)) {
        column_names <- annotation_columns
    } else if (length(column_names) != length(annotation_columns)) {
        stop("column_names must have same length as annotation_columns", call. = FALSE)
    }

    # Check for column name conflicts
    if (!overwrite) {
        existing_cols <- colnames(intervals)
        if (!is.null(dist_column) && dist_column %in% existing_cols) {
            stop(paste("Distance column", dist_column, "already exists. Use overwrite=TRUE or choose different name."), call. = FALSE)
        }

        conflicts <- intersect(column_names, existing_cols)
        if (length(conflicts) > 0) {
            stop(paste(
                "Annotation columns would overwrite existing columns:",
                paste(conflicts, collapse = ", "),
                ". Use overwrite=TRUE or provide different column_names."
            ), call. = FALSE)
        }
    }

    # Find neighbors using gintervals.neighbors
    neighbors_result <- gintervals.neighbors(intervals, annotation_intervals,
        maxneighbors = maxneighbors,
        na.if.notfound = TRUE,
        ...
    )

    # Handle empty result
    if (is.null(neighbors_result) || nrow(neighbors_result) == 0) {
        # Return original intervals with NA annotation columns, matching row count
        result <- intervals
        n <- nrow(result)
        for (i in seq_along(column_names)) {
            if (is.list(na_value) && column_names[i] %in% names(na_value)) {
                fill_val <- na_value[[column_names[i]]]
            } else {
                fill_val <- na_value
            }
            result[[column_names[i]]] <- rep(fill_val, n)
        }
        if (!is.null(dist_column)) {
            result[[dist_column]] <- rep(NA_real_, n)
        }
        result$.orig_order <- NULL
        return(result)
    }

    # The neighbors result is the query columns, then the annotation columns, then
    # the distance. Split it by position: shared column names get suffixes there,
    # and a name-based lookup would pick the query's column instead of the neighbor's.
    n1 <- ncol(intervals)
    n2 <- ncol(annotation_intervals)
    query_part <- neighbors_result[, seq_len(n1), drop = FALSE]
    colnames(query_part) <- colnames(intervals)
    ann_part <- neighbors_result[, n1 + seq_len(n2), drop = FALSE]
    colnames(ann_part) <- colnames(annotation_intervals)
    dist <- neighbors_result[[n1 + n2 + 1L]]

    # Apply tie-breaking if needed
    if (tie_method != "first" && maxneighbors > 1) {
        tie_col <- if (tie_method == "min.start") ann_part$start else ann_part$end
        # dist is signed when the annotations have a strand; ties are by absolute distance
        ord <- order(query_part$.orig_order, abs(dist), tie_col)
        query_part <- query_part[ord, , drop = FALSE]
        ann_part <- ann_part[ord, , drop = FALSE]
        dist <- dist[ord]
    }

    # Annotation values are replaced by na_value when no neighbor was found
    # (na.if.notfound gives an NA distance) or the neighbor is beyond max_dist
    to_fill <- is.na(dist)
    if (is.finite(max_dist)) {
        to_fill <- to_fill | abs(dist) > max_dist
    }

    annotations <- ann_part[, annotation_columns, drop = FALSE]
    colnames(annotations) <- column_names
    if (any(to_fill)) {
        for (i in seq_along(column_names)) {
            if (is.list(na_value) && column_names[i] %in% names(na_value)) {
                annotations[to_fill, i] <- na_value[[column_names[i]]]
            } else {
                annotations[to_fill, i] <- na_value
            }
        }
    }

    # Assemble: query columns (minus those overwritten), annotations, distance
    replaced <- c(column_names, dist_column)
    result <- cbind(query_part[, setdiff(colnames(query_part), replaced), drop = FALSE], annotations)
    if (!is.null(dist_column)) {
        result[[dist_column]] <- dist
    }
    rownames(result) <- NULL

    # Restore original order if requested
    if (keep_order) {
        result <- result[order(result$.orig_order), ]
    }
    result$.orig_order <- NULL

    # Handle intervals.set.out via helper
    return(.gintervals.save_set_or_return(result, intervals.set.out))
}
