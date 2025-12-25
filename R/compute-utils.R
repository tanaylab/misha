# Utility functions (glookup, gsample, gscreen)
#' Returns values from a lookup table based on track expression
#'
#' Evaluates track expression and translates the values into bin indices that
#' are used in turn to retrieve and return values from a lookup table.
#'
#' This function evaluates the track expression for all iterator intervals and
#' translates this value into an index based on the breaks. This index is then
#' used to address the lookup table and return the according value. More than
#' one 'expr'-'breaks' pair can be used. In that case 'lookup_table' is
#' addressed in a multidimensional manner, i.e. 'lookup_table[i1, i2, ...]'.
#'
#' The range of bins is determined by 'breaks' argument. For example: 'breaks =
#' c(x1, x2, x3, x4)' represents three different intervals (bins): (x1, x2],
#' (x2, x3], (x3, x4].
#'
#' If 'include.lowest' is 'TRUE' then the lowest value is included in the first
#' interval, i.e. in [x1, x2].
#'
#' 'force.binning' parameter controls what should be done when the value of
#' 'expr' exceeds the range determined by 'breaks'. If 'force.binning' is
#' 'TRUE' then values smaller than the minimal break will be translated to
#' index 1, and the values exceeding the maximal break will be translated to
#' index 'M-1' where 'M' is the number of breaks. If 'force.binning' is 'FALSE'
#' the out-of-range values will produce 'NaN' values.
#'
#' Regardless of 'force.binning' value if the value of 'expr' is 'NaN' then
#' result is 'NaN' too.
#'
#' The order inside the result might not be the same as the order of intervals.
#' Use 'intervalID' column to refer to the index of the original interval from
#' the supplied 'intervals'.
#'
#' If 'intervals.set.out' is not 'NULL' the result (without 'columnID' column)
#' is saved as an intervals set. Use this parameter if the result size exceeds
#' the limits of the physical memory.
#'
#' @param lookup_table a multi-dimensional array containing the values that are
#' returned by the function
#' @param ... pairs of 'expr', 'breaks' where 'expr' is a track expression and the breaks determine the bin
#' @param intervals genomic scope for which the function is applied
#' @param include.lowest if 'TRUE', the lowest value of the range determined by
#' breaks is included
#' @param force.binning if 'TRUE', the values smaller than the minimal break
#' will be translated to index 1, and the values that exceed the maximal break
#' will be translated to index N-1 where N is the number of breaks. If 'FALSE'
#' the out-of-range values will produce NaN values.
#' @param iterator track expression iterator. If 'NULL' iterator is determined
#' implicitly based on track expressions.
#' @param band track expression band. If 'NULL' no band is used.
#' @param intervals.set.out intervals set name where the function result is
#' optionally outputted
#' @return If 'intervals.set.out' is 'NULL' a set of intervals with additional
#' 'value' and 'columnID' columns.
#' @seealso \code{\link{gtrack.lookup}}, \code{\link{gextract}},
#' \code{\link{gpartition}}, \code{\link{gdist}}
#' @keywords ~lookup ~extract
#' @examples
#' \dontshow{
#' options(gmax.processes = 2)
#' }
#'
#' gdb.init_examples()
#'
#' ## one-dimensional lookup table
#' breaks1 <- seq(0.1, 0.2, length.out = 6)
#' glookup(1:5, "dense_track", breaks1, gintervals(1, 0, 200))
#'
#' ## two-dimensional lookup table
#' t <- array(1:15, dim = c(5, 3))
#' breaks2 <- seq(0.31, 0.37, length.out = 4)
#' glookup(
#'     t, "dense_track", breaks1, "2 * dense_track", breaks2,
#'     gintervals(1, 0, 200)
#' )
#'
#' @export glookup
glookup <- function(lookup_table = NULL, ..., intervals = NULL, include.lowest = FALSE, force.binning = TRUE, iterator = NULL, band = NULL, intervals.set.out = NULL) {
    args <- as.list(substitute(list(...)))[-1L]
    if (is.null(lookup_table) || length(args) < 2 || (!is.null(intervals) && length(args) %% 2 != 0) || (is.null(intervals) && length(args) %% 2 == 0)) {
        stop("Usage: glookup(lookup_table, [expr, breaks]+, intervals, include.lowest = FALSE, force.binning = TRUE, iterator = NULL, band = NULL, intervals.set.out = NULL)", call. = FALSE)
    }
    .gcheckroot()

    if (length(args) %% 2 != 0) {
        intervals <- eval.parent(args[[length(args)]])
    }

    intervals <- rescue_ALLGENOME(intervals, as.character(substitute(intervals)))

    exprs <- c()
    breaks <- list()

    for (i in (0:(length(args) / 2 - 1))) {
        exprs <- append(exprs, do.call(.gexpr2str, list(args[[i * 2 + 1]]), envir = parent.frame()))
        breaks[[length(breaks) + 1]] <- eval.parent(args[[i * 2 + 2]])
    }

    .iterator <- do.call(.giterator, list(substitute(iterator)), envir = parent.frame())
    intervals.set.out <- do.call(.gexpr2str, list(substitute(intervals.set.out)), envir = parent.frame())
    if (!is.null(intervals.set.out)) {
        fullpath <- .gintervals.check_new_set(intervals.set.out)
    }

    # intervals can be NULL if gextract is piped with gscreen and the latter returns NULL
    success <- FALSE
    res <- NULL
    tryCatch(
        {
            if (!is.null(intervals)) {
                res <- .gcall("gbintransform", intervals, exprs, breaks, include.lowest, force.binning, lookup_table, .iterator, band, intervals.set.out, .misha_env())
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


#' Returns samples from the values of track expression
#'
#' Returns a sample of the specified size from the values of track expression.
#'
#' This function returns a sample of the specified size from the values of
#' track expression. If 'n' is less than the total number of values, the data
#' is randomly sampled. The seed of the pseudo-random generator can be
#' controlled through 'grnd.seed' option.
#'
#' If 'n' is higher than the total number of values, all values are returned
#' (yet reshuffled).
#'
#' @param expr track expression
#' @param n a number of items to choose
#' @param intervals genomic scope for which the function is applied
#' @param iterator track expression iterator. If 'NULL' iterator is determined
#' implicitly based on track expression.
#' @param band track expression band. If 'NULL' no band is used.
#' @return An array that represent quantiles.
#' @seealso \code{\link{gextract}}
#' @keywords ~sample
#' @examples
#' \dontshow{
#' options(gmax.processes = 2)
#' }
#'
#' gdb.init_examples()
#' gsample("sparse_track", 10)
#'
#' @export gsample
gsample <- function(expr = NULL, n = NULL, intervals = NULL, iterator = NULL, band = NULL) {
    if (is.null(substitute(expr))) {
        stop("Usage: gsample(expr, n, intervals = .misha$ALLGENOME, iterator = NULL, band = NULL)", call. = FALSE)
    }
    .gcheckroot()

    intervals <- rescue_ALLGENOME(intervals, as.character(substitute(intervals)))

    if (is.null(intervals)) {
        intervals <- get("ALLGENOME", envir = .misha)
    }

    exprstr <- do.call(.gexpr2str, list(substitute(expr)), envir = parent.frame())
    .iterator <- do.call(.giterator, list(substitute(iterator)), envir = parent.frame())

    .gcall("C_gsample", exprstr, n, intervals, .iterator, band, .misha_env())
}


#' Finds intervals that match track expression
#'
#' Finds all intervals where track expression is 'TRUE'.
#'
#' This function finds all intervals where track expression's value is 'TRUE'.
#'
#' If 'intervals.set.out' is not 'NULL' the result is saved as an intervals
#' set. Use this parameter if the result size exceeds the limits of the
#' physical memory.
#'
#' @param expr logical track expression
#' @param intervals genomic scope for which the function is applied
#' @param iterator track expression iterator. If 'NULL' iterator is determined
#' implicitly based on track expression.
#' @param band track expression band. If 'NULL' no band is used.
#' @param intervals.set.out intervals set name where the function result is
#' optionally outputted
#' @return If 'intervals.set.out' is 'NULL' a set of intervals that match track
#' expression.
#' @seealso \code{\link{gsegment}}, \code{\link{gextract}}
#' @keywords ~screen ~interval ~intervals
#' @examples
#' \dontshow{
#' options(gmax.processes = 2)
#' }
#'
#' gdb.init_examples()
#' gscreen("dense_track > 0.2 & sparse_track < 0.4",
#'     iterator = "dense_track"
#' )
#'
#' @export gscreen
gscreen <- function(expr = NULL, intervals = NULL, iterator = NULL, band = NULL, intervals.set.out = NULL) {
    if (is.null(substitute(expr))) {
        stop("Usage: gscreen(expr, intervals = .misha$ALLGENOME, iterator = NULL, band = NULL, intervals.set.out = NULL)", call. = FALSE)
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

    success <- FALSE
    res <- NULL
    tryCatch(
        {
            if (.ggetOption("gmultitasking")) {
                res <- .gcall("gscreen_multitask", exprstr, intervals, .iterator, band, intervals.set.out, .misha_env())
            } else {
                res <- .gcall("C_gscreen", exprstr, intervals, .iterator, band, intervals.set.out, .misha_env())
            }

            if (!is.null(intervals.set.out) && .gintervals.is_bigset(intervals.set.out, FALSE) && !.gintervals.needs_bigset(intervals.set.out)) {
                .gintervals.big2small(intervals.set.out)
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
