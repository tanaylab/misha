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
