#' Partitions the values of track expression
#'
#' Converts the values of track expression to intervals that match
#' corresponding bin.
#'
#' This function converts first the values of track expression into 1-based
#' bin's index according 'breaks' argument. It returns then the intervals with
#' the corresponding bin's index.
#'
#' The range of bins is determined by 'breaks' argument. For example:
#' 'breaks=c(x1, x2, x3, x4)' represents three different intervals (bins): (x1,
#' x2], (x2, x3], (x3, x4].
#'
#' If 'include.lowest' is 'TRUE' the the lowest value will be included in the
#' first interval, i.e. in [x1, x2].
#'
#' If 'intervals.set.out' is not 'NULL' the result is saved as an intervals
#' set. Use this parameter if the result size exceeds the limits of the
#' physical memory.
#'
#' @param expr track expression
#' @param breaks breaks that determine the bin
#' @param intervals genomic scope for which the function is applied
#' @param include.lowest if 'TRUE', the lowest value of the range determined by
#' breaks is included
#' @param iterator track expression iterator. If 'NULL' iterator is determined
#' implicitly based on track expression.
#' @param band track expression band. If 'NULL' no band is used.
#' @param intervals.set.out intervals set name where the function result is
#' optionally outputted
#' @return If 'intervals.set.out' is 'NULL' a set of intervals with an
#' additional column that indicates the corresponding bin index.
#' @seealso \code{\link{gscreen}}, \code{\link{gextract}},
#' \code{\link{glookup}}, \code{\link{gdist}}
#' @keywords ~partition
#' @examples
#' \dontshow{
#' options(gmax.processes = 2)
#' }
#'
#' gdb.init_examples()
#' breaks <- seq(0, 0.2, by = 0.05)
#' gpartition("dense_track", breaks, gintervals(1, 0, 5000))
#'
#' @export gpartition
gpartition <- function(expr = NULL, breaks = NULL, intervals = NULL, include.lowest = FALSE, iterator = NULL, band = NULL, intervals.set.out = NULL) {
    if (is.null(substitute(expr)) || is.null(breaks) || is.null(intervals)) {
        stop("Usage: gpartition(expr, breaks, intervals, include.lowest = FALSE, iterator = NULL, band = NULL, intervals.set.out = NULL)", call. = FALSE)
    }
    .gcheckroot()

    intervals <- rescue_ALLGENOME(intervals, as.character(substitute(intervals)))

    exprstr <- do.call(.gexpr2str, list(substitute(expr)), envir = parent.frame())
    .iterator <- do.call(.giterator, list(substitute(iterator)), envir = parent.frame())
    intervals.set.out <- do.call(.gexpr2str, list(substitute(intervals.set.out)), envir = parent.frame())
    if (!is.null(intervals.set.out)) {
        fullpath <- .gintervals.check_new_set(intervals.set.out)
    }

    # intervals can be NULL if piped with gscreen and the latter returns NULL
    success <- FALSE
    res <- NULL
    tryCatch(
        {
            if (!is.null(intervals)) {
                res <- .gcall("C_gpartition", intervals, exprstr, breaks, include.lowest, .iterator, band, intervals.set.out, .misha_env())
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
