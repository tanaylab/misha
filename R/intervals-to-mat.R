# gintervals.to_mat / gintervals.from_mat
# See dev/notes/features/2026-05-17_gintervals-to-mat-design.md

#' Convert intervals + values data.frame to an interval-indexed matrix
#'
#' Builds a numeric matrix of value columns whose rows are indexed by
#' intervals. The intervals are carried in \code{attr(mat, "intervals")} as the
#' authoritative identity; \code{rownames(mat)} are display-only (default
#' \code{"chrom:start-end"}) and are NEVER parsed back by
#' \code{\link{gintervals.from_mat}}. This avoids the round-trip corruption
#' that occurs when chrom names contain underscores or other separators.
#'
#' @param df data.frame with \code{chrom}, \code{start}, \code{end} and zero
#'   or more value columns. May contain an \code{intervalID} column (from
#'   \code{gextract}), which is kept in the attribute but excluded from the
#'   matrix.
#' @param id_col optional column name in \code{df} whose values become
#'   rownames. If \code{NULL} (default), rownames are
#'   \code{"chrom:start-end"}.
#' @param value_cols character vector of column names to use as matrix data.
#'   If \code{NULL} (default), auto-detect: all columns except \code{chrom},
#'   \code{start}, \code{end}, \code{intervalID}. Auto-detect errors if any
#'   selected column is non-numeric; pass \code{value_cols} explicitly to
#'   override.
#' @param labels if \code{TRUE} (default), set \code{rownames(mat)} to either
#'   \code{df[[id_col]]} (if \code{id_col} is supplied) or
#'   \code{"chrom:start-end"}. If \code{FALSE}, leave rownames \code{NULL} -
#'   useful in pipelines that don't need the display labels and would prefer to
#'   skip the construction cost on large inputs. When \code{FALSE}, the
#'   \code{id_col} argument is ignored.
#'
#' @return An \code{intervs_mat} object: a numeric matrix subclass with the
#'   intervals attached as \code{attr(., "intervals")}. Supports row/column
#'   subsetting (\code{[}) and \code{rbind()} while preserving the attribute.
#'
#' @seealso \code{\link{gintervals.from_mat}}
#'
#' @examples
#' df <- data.frame(
#'     chrom = c("chr1", "chr1", "chr2"),
#'     start = c(100L, 500L, 200L),
#'     end   = c(200L, 700L, 400L),
#'     t1    = c(1.5, 2.5, 3.5),
#'     t2    = c(10, 20, 30)
#' )
#' mat <- gintervals.to_mat(df)
#' rownames(mat)
#' # subset preserves intervals:
#' sub <- mat[c(1, 3), ]
#' attr(sub, "intervals")
#' # round-trip back to a data.frame:
#' gintervals.from_mat(sub)
#'
#' @export
gintervals.to_mat <- function(df, id_col = NULL, value_cols = NULL, labels = TRUE) {
    if (!is.data.frame(df)) {
        stop("`df` must be a data.frame", call. = FALSE)
    }
    required <- c("chrom", "start", "end")
    missing_cols <- setdiff(required, names(df))
    if (length(missing_cols) > 0) {
        stop(sprintf(
            "`df` is missing required interval column(s): %s",
            paste(missing_cols, collapse = ", ")
        ), call. = FALSE)
    }

    has_intervalID <- "intervalID" %in% names(df)
    identity_cols <- c(
        "chrom", "start", "end",
        if (has_intervalID) "intervalID"
    )

    if (is.null(value_cols)) {
        value_cols <- setdiff(names(df), identity_cols)
        non_numeric <- value_cols[!vapply(df[value_cols], is.numeric, logical(1))]
        if (length(non_numeric) > 0) {
            stop(sprintf(
                "Non-numeric value column(s): %s. Use `value_cols = ` to select numeric columns explicitly.",
                paste(non_numeric, collapse = ", ")
            ), call. = FALSE)
        }
    } else {
        missing_vals <- setdiff(value_cols, names(df))
        if (length(missing_vals) > 0) {
            stop(sprintf(
                "`value_cols` not found in `df`: %s",
                paste(missing_vals, collapse = ", ")
            ), call. = FALSE)
        }
    }

    mat <- as.matrix(df[, value_cols, drop = FALSE])

    if (labels) {
        if (is.null(id_col)) {
            rownames(mat) <- .Call(
                "C_intervals_coord_strings",
                df$chrom,
                as.integer(df$start),
                as.integer(df$end)
            )
        } else {
            if (!(id_col %in% names(df))) {
                stop(sprintf("`id_col` not found in `df`: %s", id_col), call. = FALSE)
            }
            rownames(mat) <- as.character(df[[id_col]])
        }
    }

    attr(mat, "intervals") <- df[, identity_cols, drop = FALSE]
    class(mat) <- c("intervs_mat", "matrix", "array")
    mat
}

#' Subset an intervs_mat preserving interval identity
#'
#' Row subset of an \code{intervs_mat} subsets \code{attr(., "intervals")} in
#' parallel. Column-only subset leaves the intervals attribute unchanged.
#' Single-row results with \code{drop = TRUE} return a plain numeric vector
#' (class dropped).
#'
#' @param x an \code{intervs_mat}
#' @param i row selector (logical, integer, or character matching rownames)
#' @param j column selector
#' @param ... unused
#' @param drop if \code{TRUE} and a single row is selected via \code{i},
#'   collapse to a named vector and drop the class. Defaults to \code{TRUE}
#'   to match base matrix \code{[} for row selection.
#'
#' @return An \code{intervs_mat} (still 2D), or a numeric vector (degenerate).
#'
#' @export
`[.intervs_mat` <- function(x, i, j, ..., drop = TRUE) {
    intervals <- attr(x, "intervals")
    bare <- unclass(x)
    attr(bare, "intervals") <- NULL

    i_missing <- missing(i)
    j_missing <- missing(j)

    # Always subset with drop=FALSE to preserve matrix dims.  We apply the
    # caller's `drop` ourselves at the end, but only along the row axis and
    # only when the caller actually supplied `i`.  Column-only subsets never
    # collapse rows regardless of `drop`.
    if (i_missing && j_missing) {
        out <- bare[, , drop = FALSE]
        row_selected <- FALSE
    } else if (i_missing) {
        out <- bare[, j, drop = FALSE]
        row_selected <- FALSE
    } else if (j_missing) {
        out <- bare[i, , drop = FALSE]
        intervals <- intervals[i, , drop = FALSE]
        rownames(intervals) <- NULL
        row_selected <- TRUE
    } else {
        out <- bare[i, j, drop = FALSE]
        intervals <- intervals[i, , drop = FALSE]
        rownames(intervals) <- NULL
        row_selected <- TRUE
    }

    # Collapse to a named vector only when the caller specified `i` and the
    # result has a single row and `drop = TRUE`.
    if (row_selected && drop && nrow(out) == 1L) {
        v <- out[1L, ]
        return(v)
    }

    attr(out, "intervals") <- intervals
    class(out) <- c("intervs_mat", "matrix", "array")
    out
}

#' Convert an interval-indexed matrix back to an intervals + values data.frame
#'
#' Inverse of \code{\link{gintervals.to_mat}}. Recovers intervals from
#' \code{attr(mat, "intervals")} if \code{mat} is an \code{intervs_mat}, or
#' from the explicit \code{intervals} argument if \code{mat} is a plain matrix.
#' Never parses rownames.
#'
#' @param mat an \code{intervs_mat} (produced by \code{gintervals.to_mat}) or a
#'   plain matrix (then \code{intervals} must be supplied).
#' @param intervals data.frame with \code{chrom}, \code{start}, \code{end}
#'   columns, required when \code{mat} is a plain matrix. Must satisfy
#'   \code{nrow(intervals) == nrow(mat)}; alignment is strictly positional.
#'   Must NOT be supplied when \code{mat} is already an \code{intervs_mat}.
#'
#' @return A data.frame with \code{chrom}, \code{start}, \code{end} (plus
#'   \code{intervalID} if it was present in the original input), followed by
#'   the value columns.
#'
#' @seealso \code{\link{gintervals.to_mat}}
#'
#' @examples
#' df <- data.frame(
#'     chrom = c("chr1", "chr2"),
#'     start = c(100L, 200L),
#'     end   = c(200L, 400L),
#'     t1    = c(1.0, 2.0)
#' )
#' mat <- gintervals.to_mat(df)
#' identical(gintervals.from_mat(mat)$t1, df$t1)
#'
#' # plain matrix path:
#' plain <- unclass(mat)
#' attr(plain, "intervals") <- NULL
#' gintervals.from_mat(plain, intervals = df[, c("chrom", "start", "end")])
#'
#' @export
gintervals.from_mat <- function(mat, intervals = NULL) {
    if (inherits(mat, "intervs_mat")) {
        if (!is.null(intervals)) {
            stop("`intervals` must not be supplied when `mat` is an intervs_mat (it already carries intervals).", call. = FALSE)
        }
        intervals <- attr(mat, "intervals")
    } else {
        if (is.null(intervals)) {
            stop("`mat` is a plain matrix; supply `intervals = ` (data.frame with chrom/start/end) to recover interval identity.", call. = FALSE)
        }
        if (!is.data.frame(intervals)) {
            stop("`intervals` must be a data.frame", call. = FALSE)
        }
        missing_cols <- setdiff(c("chrom", "start", "end"), names(intervals))
        if (length(missing_cols) > 0) {
            stop(sprintf(
                "`intervals` is missing required column(s): %s",
                paste(missing_cols, collapse = ", ")
            ), call. = FALSE)
        }
        if (nrow(intervals) != nrow(mat)) {
            stop(sprintf(
                "`intervals` has %d rows, `mat` has %d rows; they must match.",
                nrow(intervals), nrow(mat)
            ), call. = FALSE)
        }
    }

    bare <- unclass(mat)
    attr(bare, "intervals") <- NULL
    vals <- as.data.frame(bare, stringsAsFactors = FALSE)
    rownames(vals) <- NULL

    out <- cbind(intervals, vals, stringsAsFactors = FALSE)
    rownames(out) <- NULL
    out
}

#' Row-bind intervs_mat objects concatenating their intervals
#'
#' Concatenates the rows of the matrices AND the rows of their
#' \code{"intervals"} attributes. If any input is not an \code{intervs_mat},
#' falls back to plain matrix \code{rbind} (the result is a base matrix with
#' no intervals attribute).
#'
#' @param ... \code{intervs_mat} objects (and/or other matrix-like inputs).
#' @param deparse.level passed to base \code{rbind}.
#'
#' @return An \code{intervs_mat} if all inputs were \code{intervs_mat};
#'   otherwise a plain matrix.
#'
#' @export
rbind.intervs_mat <- function(..., deparse.level = 1) {
    args <- list(...)
    all_ours <- all(vapply(args, inherits, logical(1), what = "intervs_mat"))

    bare <- lapply(args, function(a) {
        b <- unclass(a)
        attr(b, "intervals") <- NULL
        b
    })
    out <- do.call(rbind, c(bare, list(deparse.level = deparse.level)))

    if (!all_ours) {
        return(out)
    }

    intervals <- do.call(rbind, lapply(args, attr, "intervals"))
    rownames(intervals) <- NULL
    attr(out, "intervals") <- intervals
    class(out) <- c("intervs_mat", "matrix", "array")
    out
}
