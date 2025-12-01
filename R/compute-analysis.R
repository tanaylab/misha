# Statistical analysis functions (gsegment, gwilcox)
#' Divides track expression into segments
#'
#' Divides the values of track expression into segments by using Wilcoxon test.
#'
#' This function divides the values of track expression into segments, where
#' each segment size is at least of 'minsegment' size and the P-value of
#' comparing the segment with the first 'minsegment' values from the next
#' segment is at most 'maxpval'. Comparison is done using Wilcoxon (also known
#' as Mann-Whitney) test.
#'
#' If 'intervals.set.out' is not 'NULL' the result is saved as an intervals
#' set. Use this parameter if the result size exceeds the limits of the
#' physical memory.
#'
#' @param expr track expression
#' @param minsegment minimal segment size
#' @param maxpval maximal P-value that separates two adjacent segments
#' @param onetailed if 'TRUE', Wilcoxon test is performed one tailed, otherwise
#' two tailed
#' @param intervals genomic scope for which the function is applied
#' @param iterator track expression iterator of "fixed bin" type. If 'NULL'
#' iterator is determined implicitly based on track expression.
#' @param intervals.set.out intervals set name where the function result is
#' optionally outputted
#' @return If 'intervals.set.out' is 'NULL' a set of intervals where each
#' interval represents a segment.
#' @seealso \code{\link{gscreen}}, \code{\link{gwilcox}}
#' @keywords ~segment ~wilcoxon ~Mann-Whitney
#' @examples
#' \dontshow{
#' options(gmax.processes = 2)
#' }
#'
#' gdb.init_examples()
#' gsegment("dense_track", 5000, 0.0001)
#'
#' @export gsegment
gsegment <- function(expr = NULL, minsegment = NULL, maxpval = 0.05, onetailed = TRUE, intervals = NULL, iterator = NULL, intervals.set.out = NULL) {
    if (is.null(substitute(expr)) || is.null(minsegment)) {
        stop("Usage: gsegment(expr, minsegment, maxpval = 0.05, onetailed = TRUE, intervals = .misha$ALLGENOME, iterator = NULL, intervals.set.out = NULL)", call. = FALSE)
    }
    .gcheckroot()

    intervals <- rescue_ALLGENOME(intervals, as.character(substitute(intervals)))

    if (is.null(intervals)) {
        intervals <- get("ALLGENOME", envir = .misha)
    }

    exprstr <- do.call(.gexpr2str, list(substitute(expr)), envir = parent.frame())
    .iterator <- do.call(.giterator, list(substitute(iterator)), envir = parent.frame())
    intervals.set.out <- do.call(.gexpr2str, list(substitute(intervals.set.out)), envir = parent.frame())
    if (!is.null(intervals.set.out)) {
        fullpath <- .gintervals.check_new_set(intervals.set.out)
    }

    if (!onetailed) {
        maxpval <- maxpval / 2
    }

    # intervals can be NULL if piped with gscreen and the latter returns NULL
    success <- FALSE
    res <- NULL
    tryCatch(
        {
            if (!is.null(intervals)) {
                res <- .gcall("C_gsegment", exprstr, intervals, minsegment, stats::qnorm(maxpval), onetailed, .iterator, intervals.set.out, .misha_env())
                if (!is.null(intervals.set.out) && .gintervals.is_bigset(intervals.set.out, FALSE) && !.gintervals.needs_bigset(intervals.set.out)) {
                    .gintervals.big2small(intervals.set.out)
                }
            }
            success <- TRUE
        },
        finally = {
            if (!success && !is.null(intervals.set.out)) {
                unlink(fullpath, recursive = TRUE)
            }
        }
    )

    # refresh the list of GINTERVS, etc.
    if (!is.null(intervals.set.out)) {
        .gdb.add_intervals.set(intervals.set.out)
        retv <- 0 # suppress return value
    } else {
        res
    }
}


#' Calculates Wilcoxon test on sliding windows over track expression
#'
#' Calculates Wilcoxon test on sliding windows over the values of track
#' expression.
#'
#' This function runs a Wilcoxon test (also known as a Mann-Whitney test) over
#' the values of track expression in the two sliding windows having an
#' identical center. The sizes of the windows are specified by 'winsize1' and
#' 'winsize2'. 'gwilcox' returns intervals where the smaller window tested
#' against a larger window gives a P-value below 'maxpval'. The test can be one
#' or two tailed.
#'
#' 'what2find' argument controls what should be searched: peaks, lows or both.
#'
#' If 'intervals.set.out' is not 'NULL' the result is saved as an intervals
#' set. Use this parameter if the result size exceeds the limits of the
#' physical memory.
#'
#' @param expr track expression
#' @param winsize1 number of values in the first sliding window
#' @param winsize2 number of values in the second sliding window
#' @param maxpval maximal P-value
#' @param onetailed if 'TRUE', Wilcoxon test is performed one tailed, otherwise
#' two tailed
#' @param what2find if '-1', lows are searched. If '1', peaks are searched. If
#' '0', both peaks and lows are searched
#' @param intervals genomic scope for which the function is applied
#' @param iterator track expression iterator of "fixed bin" type. If 'NULL'
#' iterator is determined implicitly based on track expression.
#' @param intervals.set.out intervals set name where the function result is
#' optionally outputted
#' @return If 'intervals.set.out' is 'NULL' a data frame representing the
#' intervals with an additional 'pval' column where P-value is below 'maxpval'.
#' @seealso \code{\link{gscreen}}, \code{\link{gsegment}}
#' @keywords ~wilcoxon ~Mann-Whitney
#' @examples
#' \dontshow{
#' options(gmax.processes = 2)
#' }
#'
#' gdb.init_examples()
#' gwilcox("dense_track", 100000, 1000,
#'     maxpval = 0.01,
#'     what2find = 1
#' )
#'
#' @export gwilcox
gwilcox <- function(expr = NULL, winsize1 = NULL, winsize2 = NULL, maxpval = 0.05, onetailed = TRUE, what2find = 1, intervals = NULL, iterator = NULL, intervals.set.out = NULL) {
    if (is.null(substitute(expr)) || is.null(winsize1) || is.null(winsize2)) {
        stop("Usage: gwilcox(expr, winsize1, winsize2, maxpval = 0.05, onetailed = TRUE, what2find = 1 (-1=lows, 0=lows/highs, 1=highs), intervals = .misha$ALLGENOME, iterator = NULL, intervals.set.out = NULL)", call. = FALSE)
    }
    .gcheckroot()

    intervals <- rescue_ALLGENOME(intervals, as.character(substitute(intervals)))

    if (is.null(intervals)) {
        intervals <- get("ALLGENOME", envir = .misha)
    }

    exprstr <- do.call(.gexpr2str, list(substitute(expr)), envir = parent.frame())
    .iterator <- do.call(.giterator, list(substitute(iterator)), envir = parent.frame())
    intervals.set.out <- do.call(.gexpr2str, list(substitute(intervals.set.out)), envir = parent.frame())
    if (!is.null(intervals.set.out)) {
        fullpath <- .gintervals.check_new_set(intervals.set.out)
    }

    if (!onetailed) {
        maxpval <- maxpval / 2
    }

    # intervals can be NULL if piped with gscreen and the latter returns NULL
    success <- FALSE
    res <- NULL
    tryCatch(
        {
            if (!is.null(intervals)) {
                res <- .gcall("C_gwilcox", exprstr, intervals, winsize1, winsize2, qnorm(maxpval), onetailed, as.integer(what2find), .iterator, intervals.set.out, .misha_env())
                if (!is.null(intervals.set.out) && .gintervals.is_bigset(intervals.set.out, FALSE) && !.gintervals.needs_bigset(intervals.set.out)) {
                    .gintervals.big2small(intervals.set.out)
                }
            }
            success <- TRUE
        },
        finally = {
            if (!success && !is.null(intervals.set.out)) {
                unlink(fullpath, recursive = TRUE)
            }
        }
    )

    # refresh the list of GINTERVS, etc.
    if (!is.null(intervals.set.out)) {
        .gdb.add_intervals.set(intervals.set.out)
        retv <- 0 # suppress return value
    } else {
        res
    }
}


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
