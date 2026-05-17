# gintervals.to_mat / gintervals.from_mat
# See dev/notes/features/2026-05-17_gintervals-to-mat-design.md

gintervals.to_mat <- function(df, id_col = NULL, value_cols = NULL) {
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

    if (is.null(id_col)) {
        rownames(mat) <- paste0(df$chrom, ":", df$start, "-", df$end)
    } else {
        if (!(id_col %in% names(df))) {
            stop(sprintf("`id_col` not found in `df`: %s", id_col), call. = FALSE)
        }
        rownames(mat) <- as.character(df[[id_col]])
    }

    attr(mat, "intervals") <- df[, identity_cols, drop = FALSE]
    class(mat) <- c("intervs_mat", "matrix", "array")
    mat
}

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
