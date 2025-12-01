# Track expression operations (mapply, quantiles, summary)

#' Applies a function to values of track expressions
#'
#' Applies a function to values of track expressions for each interval.
#'
#' This function evaluates track expressions for each interval from
#' 'intervals'. The resulted vectors are passed then as arguments to 'FUN'.
#'
#' If the intervals are one-dimensional and have an additional column named
#' 'strand' whose value is '-1', the values of the track expression are placed
#' to the vector in reverse order.
#'
#' The current interval index (1-based) is stored in 'GAPPLY.INTERVID' variable
#' that is available during the execution of 'gintervals.mapply'. There is no
#' guarantee about the order in which the intervals are processed. Do not rely
#' on any specific order and use 'GITERATOR.INTERVID' variable to detect the
#' current interval id.
#'
#' If 'enable.gapply.intervals' is 'TRUE', an additional variable
#' 'GAPPLY.INTERVALS' is defined during the execution of 'gintervals.mapply'.
#' This variable stores the current iterator intervals prior to track
#' expression evaluation. Please note that setting 'enable.gapply.intervals' to
#' 'TRUE' might severely affect the run-time of the function.
#'
#' Note: all the changes made in R environment by 'FUN' will be void if
#' multitasking mode is switched on. One should also refrain from performing
#' any other operations in 'FUN' that might be not "thread-safe" such as
#' updating files, etc. Please switch off multitasking ('options(gmultitasking
#' = FALSE)') if you wish to perform such operations.
#'
#' If 'intervals.set.out' is not 'NULL' the result is saved as an intervals
#' set. Use this parameter if the result size exceeds the limits of the
#' physical memory.
#'
#' @param FUN function to apply, found via 'match.fun'
#' @param ... track expressions whose values are used as arguments for 'FUN'
#' @param intervals intervals for which track expressions are calculated
#' @param enable.gapply.intervals if 'TRUE', then a variable 'GAPPLY.INTERVALS'
#' is available
#' @param iterator track expression iterator. If 'NULL' iterator is determined
#' implicitly based on track expressions.
#' @param band track expression band. If 'NULL' no band is used.
#' @param intervals.set.out intervals set name where the function result is
#' optionally outputted
#' @param colnames name of the column that contains the return values of 'FUN'.
#' Default is "value".
#' @return If 'intervals.set.out' is 'NULL' a data frame representing intervals
#' with an additional column that contains the return values of 'FUN'. The name
#' of this additional column is specified by the 'colnames' parameter.
#' @seealso \code{\link{mapply}}
#' @keywords ~apply ~mapply
#' @examples
#' \dontshow{
#' options(gmultitasking = FALSE)
#' }
#'
#' gdb.init_examples()
#' gintervals.mapply(
#'     max, "dense_track",
#'     gintervals(c(1, 2), 0, 10000)
#' )
#' gintervals.mapply(
#'     function(x, y) {
#'         max(x + y)
#'     }, "dense_track",
#'     "sparse_track", gintervals(c(1, 2), 0, 10000),
#'     iterator = "sparse_track"
#' )
#' # Using custom column name
#' gintervals.mapply(
#'     max, "dense_track",
#'     gintervals(c(1, 2), 0, 10000),
#'     colnames = "max_value"
#' )
#'
#' @export gintervals.mapply
gintervals.mapply <- function(FUN = NULL, ..., intervals = NULL, enable.gapply.intervals = FALSE, iterator = NULL, band = NULL, intervals.set.out = NULL, colnames = "value") {
    assign("GINTERVID", -1, envir = .misha)
    args <- as.list(substitute(list(...)))[-1L]
    if (is.null(intervals) && length(args) < 2 || !is.null(intervals) && length(args) < 1) {
        stop("Usage: gintervals.mapply(FUN, [expr]+, intervals, enable.gapply.intervals = FALSE, iterator = NULL, intervals.set.out = NULL, colnames = \"value\")", call. = FALSE)
    }
    .gcheckroot()

    intervals <- rescue_ALLGENOME(intervals, as.character(substitute(intervals)))

    if (is.null(intervals)) {
        intervals <- eval.parent(args[[length(args)]])
        args <- args[1:(length(args) - 1)]
    }

    tracks <- c()
    for (track in args) {
        tracks <- c(tracks, do.call(.gexpr2str, list(track), envir = parent.frame()))
    }

    .iterator <- do.call(.giterator, list(substitute(iterator)), envir = parent.frame())

    if (exists("GAPPLY.INTERVALS", envir = .misha)) {
        remove(list = "GAPPLY.INTERVALS", envir = .misha)
    }

    intervals.set.out <- do.call(.gexpr2str, list(substitute(intervals.set.out)), envir = parent.frame())


    if (.gintervals.is_bigset(intervals) || !is.null(intervals.set.out)) {
        res <- NULL

        INTERVALS_FUN <- function(intervals, intervals.set.out, envir) {
            intervals <- intervals[[1]]
            chrom_res <- .gcall("gmapply", intervals, FUN, tracks, enable.gapply.intervals, .iterator, band, FALSE, colnames, .misha_env())
            if (!is.null(chrom_res) && nrow(chrom_res) > 0) {
                if (is.null(intervals.set.out)) {
                    assign("res", c(get("res", envir = envir), list(chrom_res)), envir = envir)
                    .gverify_max_data_size(sum(unlist(lapply(get("res", envir), nrow))), arguments = "intervals.set.out")
                }
            }
            chrom_res
        }

        .gintervals.apply(gintervals.chrom_sizes(intervals), intervals, intervals.set.out, INTERVALS_FUN, intervals.set.out, environment())

        if (!is.null(res)) {
            res <- do.call(.grbind, res)
        } # much faster than calling rbind incrementally in FUN

        if (is.null(intervals.set.out)) {
            if (!is.null(res) && nrow(res)) {
                res
            } else {
                NULL
            }
        } else {
            retv <- 0
        } # suppress return value
    } else {
        if (.ggetOption("gmultitasking")) {
            .gcall("gmapply_multitask", intervals, FUN, tracks, enable.gapply.intervals, .iterator, band, TRUE, colnames, .misha_env())
        } else {
            .gcall("gmapply", intervals, FUN, tracks, enable.gapply.intervals, .iterator, band, TRUE, colnames, .misha_env())
        }
    }
}


#' Calculates quantiles of a track expression for intervals
#'
#' Calculates quantiles of a track expression for intervals.
#'
#' This function calculates quantiles of 'expr' for each interval in
#' 'intervals'.
#'
#' If 'intervals.set.out' is not 'NULL' the result is saved as an intervals
#' set. Use this parameter if the result size exceeds the limits of the
#' physical memory.
#'
#' @param expr track expression for which quantiles are calculated
#' @param percentiles an array of percentiles of quantiles in [0, 1] range
#' @param intervals set of intervals
#' @param iterator track expression iterator. If 'NULL' iterator is determined
#' implicitly based on track expressions.
#' @param band track expression band. If 'NULL' no band is used.
#' @param intervals.set.out intervals set name where the function result is
#' optionally outputted
#' @return If 'intervals.set.out' is 'NULL' a set of intervals with additional
#' columns representing quantiles for each percentile.
#' @seealso \code{\link{gquantiles}}, \code{\link{gbins.quantiles}}
#' @keywords ~quantiles ~percentiles
#' @examples
#' \dontshow{
#' options(gmax.processes = 2)
#' }
#'
#' gdb.init_examples()
#' intervs <- gintervals(c(1, 2), 0, 5000)
#' gintervals.quantiles("dense_track",
#'     percentiles = c(0.5, 0.3, 0.9), intervs
#' )
#'
#' @export gintervals.quantiles
gintervals.quantiles <- function(expr = NULL, percentiles = 0.5, intervals = NULL, iterator = NULL, band = NULL, intervals.set.out = NULL) {
    if (is.null(substitute(expr)) || is.null(intervals)) {
        stop("Usage: gintervals.quantiles(expr, percentiles = 0.5, intervals, iterator = NULL, band = NULL, intervals.set.out = NULL)", call. = FALSE)
    }
    .gcheckroot()

    intervals <- rescue_ALLGENOME(intervals, as.character(substitute(intervals)))

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
                res <- .gcall("gintervals_quantiles_multitask", intervals, exprstr, percentiles, .iterator, band, intervals.set.out, .misha_env())
            } else {
                res <- .gcall("gintervals_quantiles", intervals, exprstr, percentiles, .iterator, band, intervals.set.out, .misha_env())
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
    if (is.null(intervals.set.out)) {
        res
    } else {
        .gdb.add_intervals.set(intervals.set.out)
        retv <- 0 # suppress return value
    }
}


#' Calculates summary statistics of track expression for intervals
#'
#' Calculates summary statistics of track expression for intervals.
#'
#' This function returns summary statistics of a track expression for each
#' interval 'intervals': total number of bins, total number of bins whose value
#' is NaN, min, max, sum, mean and standard deviation of the values.
#'
#' If 'intervals.set.out' is not 'NULL' the result is saved as an intervals
#' set. Use this parameter if the result size exceeds the limits of the
#' physical memory.
#'
#' @param expr track expression
#' @param intervals set of intervals
#' @param iterator track expression iterator. If 'NULL' iterator is determined
#' implicitly based on track expression.
#' @param band track expression band. If 'NULL' no band is used.
#' @param intervals.set.out intervals set name where the function result is
#' optionally outputted
#' @return If 'intervals.set.out' is 'NULL' a set of intervals with additional
#' columns representing summary statistics for each percentile and interval.
#' @seealso \code{\link{gsummary}}, \code{\link{gbins.summary}}
#' @keywords ~summary ~statistics
#' @examples
#' \dontshow{
#' options(gmax.processes = 2)
#' }
#'
#' gdb.init_examples()
#' intervs <- gintervals(c(1, 2), 0, 5000)
#' gintervals.summary("dense_track", intervs)
#'
#' @export gintervals.summary
gintervals.summary <- function(expr = NULL, intervals = NULL, iterator = NULL, band = NULL, intervals.set.out = NULL) {
    if (is.null(substitute(expr)) || is.null(intervals)) {
        stop("Usage: gintervals.summary(expr, intervals, iterator = NULL, band = NULL, intervals.set.out = NULL)", call. = FALSE)
    }
    .gcheckroot()

    intervals <- rescue_ALLGENOME(intervals, as.character(substitute(intervals)))

    exprstr <- do.call(.gexpr2str, list(substitute(expr)), envir = parent.frame())
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
                res <- .gcall("gintervals_summary", exprstr, intervals, .iterator, band, intervals.set.out, .misha_env())
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
