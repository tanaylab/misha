# Bin-based functions (gbins.quantiles, gbins.summary)

#' Calculates quantiles of a track expression for bins
#'
#' Calculates quantiles of a track expression for bins.
#'
#' This function is a binned version of 'gquantiles'. For each iterator
#' interval the value of 'bin_expr' is calculated and assigned to the
#' corresponding bin determined by 'breaks'. The quantiles of 'expr' are
#' calculated then separately for each bin.
#'
#' The bins can be multi-dimensional depending on the number of
#' 'bin_expr'-'breaks' pairs.
#'
#' The range of bins is determined by 'breaks' argument. For example:
#' 'breaks=c(x1, x2, x3, x4)' represents three different intervals (bins): (x1,
#' x2], (x2, x3], (x3, x4].
#'
#' If 'include.lowest' is 'TRUE' the the lowest value will be included in the
#' first interval, i.e. in [x1, x2].
#'
#' @param ... pairs of track expressions ('bin_expr') that determines the bins and breaks that define the bins. See \code{\link{gdist}}.
#' @param expr track expression for which quantiles are calculated
#' @param percentiles an array of percentiles of quantiles in [0, 1] range
#' @param intervals genomic scope for which the function is applied.
#' @param include.lowest if 'TRUE', the lowest value of the range determined by
#' breaks is included
#' @param iterator track expression iterator. If 'NULL' iterator is determined
#' implicitly based on track expressions.
#' @param band track expression band. If 'NULL' no band is used.
#' @return Multi-dimensional array representing quantiles for each percentile
#' and bin.
#' @seealso \code{\link{gquantiles}}, \code{\link{gintervals.quantiles}},
#' \code{\link{gdist}}
#' @keywords ~quantiles ~percentiles
#' @examples
#' \dontshow{
#' options(gmax.processes = 2)
#' }
#'
#' gdb.init_examples()
#' gbins.quantiles("dense_track", c(0, 0.2, 0.4, 2), "sparse_track",
#'     percentiles = c(0.2, 0.5),
#'     intervals = gintervals(1),
#'     iterator = "dense_track"
#' )
#'
#' @export gbins.quantiles
gbins.quantiles <- function(..., expr = NULL, percentiles = 0.5, intervals = get("ALLGENOME", envir = .misha), include.lowest = FALSE, iterator = NULL, band = NULL) {
    args <- as.list(substitute(list(...)))[-1L]

    if (length(args) >= 0 && length(args) %% 2 != 0) {
        expr <- args[[length(args)]]
    }

    if (length(args) < 2 || is.null(substitute(expr))) {
        stop("Usage: gbins.quantiles([bin_expr, breaks]+, expr, percentiles = 0.5, intervals = .misha$ALLGENOME, include.lowest = FALSE, iterator = NULL, band = NULL)", call. = FALSE)
    }
    .gcheckroot()

    intervals <- rescue_ALLGENOME(intervals, as.character(substitute(intervals)))

    exprs <- c()
    breaks <- list()

    exprs <- append(exprs, do.call(.gexpr2str, list(substitute(expr)), envir = parent.frame()))
    for (i in (0:((length(args) - 1) / 2 - 1))) {
        exprs <- append(exprs, do.call(.gexpr2str, list(args[[i * 2 + 1]]), envir = parent.frame()))
        breaks[[length(breaks) + 1]] <- eval.parent(args[[i * 2 + 2]])
    }

    .iterator <- do.call(.giterator, list(substitute(iterator)), envir = parent.frame())

    res <- .gcall("gbins_quantiles", exprs, breaks, include.lowest, percentiles, intervals, .iterator, band, .misha_env())
    attr(res, "breaks") <- breaks
    res
}


#' Calculates summary statistics of a track expression for bins
#'
#' Calculates summary statistics of a track expression for bins.
#'
#' This function is a binned version of 'gsummary'. For each iterator interval
#' the value of 'bin_expr' is calculated and assigned to the corresponding bin
#' determined by 'breaks'. The summary statistics of 'expr' are calculated then
#' separately for each bin.
#'
#' The bins can be multi-dimensional depending on the number of
#' 'bin_expr'-'breaks' pairs.
#'
#' The range of bins is determined by 'breaks' argument. For example:
#' 'breaks=c(x1, x2, x3, x4)' represents three different intervals (bins): (x1,
#' x2], (x2, x3], (x3, x4].
#'
#' If 'include.lowest' is 'TRUE' the the lowest value will be included in the
#' first interval, i.e. in [x1, x2].
#'
#' @param ... pairs of track expressions ('bin_expr') that determines the bins and breaks that define the bins. See \code{\link{gdist}}.
#' @param expr track expression for which summary statistics is calculated
#' @param intervals genomic scope for which the function is applied
#' @param include.lowest if 'TRUE', the lowest value of the range determined by
#' breaks is included
#' @param iterator track expression iterator. If 'NULL' iterator is determined
#' implicitly based on track expressions.
#' @param band track expression band. If 'NULL' no band is used.
#' @return Multi-dimensional array representing summary statistics for each
#' bin.
#' @seealso \code{\link{gsummary}}, \code{\link{gintervals.summary}},
#' \code{\link{gdist}}
#' @keywords ~summary
#' @examples
#' \dontshow{
#' options(gmax.processes = 2)
#' }
#'
#' gdb.init_examples()
#' gbins.summary("dense_track", c(0, 0.2, 0.4, 2), "sparse_track",
#'     intervals = gintervals(1), iterator = "dense_track"
#' )
#'
#' @export gbins.summary
gbins.summary <- function(..., expr = NULL, intervals = get("ALLGENOME", envir = .misha), include.lowest = FALSE, iterator = NULL, band = NULL) {
    args <- as.list(substitute(list(...)))[-1L]

    if (length(args) >= 0 && length(args) %% 2 != 0) {
        expr <- args[[length(args)]]
    }

    if (length(args) < 2 || is.null(substitute(expr))) {
        stop("Usage: gbins.summary([expr, breaks]+, expr, intervals = .misha$ALLGENOME, include.lowest = FALSE, iterator = NULL, band = NULL)", call. = FALSE)
    }
    .gcheckroot()

    intervals <- rescue_ALLGENOME(intervals, as.character(substitute(intervals)))

    exprs <- c()
    breaks <- list()

    exprs <- append(exprs, do.call(.gexpr2str, list(substitute(expr)), envir = parent.frame()))
    for (i in (0:((length(args) - 1) / 2 - 1))) {
        exprs <- append(exprs, do.call(.gexpr2str, list(args[[i * 2 + 1]]), envir = parent.frame()))
        breaks[[length(breaks) + 1]] <- eval.parent(args[[i * 2 + 2]])
    }

    .iterator <- do.call(.giterator, list(substitute(iterator)), envir = parent.frame())

    res <- .gcall("gbins_summary", exprs, breaks, include.lowest, intervals, .iterator, band, .misha_env())
    attr(res, "breaks") <- breaks
    res
}
