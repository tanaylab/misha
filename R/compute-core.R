#' Returns evaluated track expression
#'
#' Returns the result of track expressions evaluation for each of the iterator
#' intervals.
#'
#' This function returns the result of track expressions evaluation for each of
#' the iterator intervals. The returned value is a set of intervals with an
#' additional column for each of the track expressions. This value can be used
#' as an input for any other function that accepts intervals. If the intervals
#' inside 'intervals' argument overlap gextract returns the overlapped
#' coordinate more than once.
#'
#' The order inside the result might not be the same as the order of intervals.
#' An additional column 'intervalID' is added to the return value. Use this
#' column to refer to the index of the original interval from the supplied
#' 'intervals'.
#'
#' If 'file' parameter is not 'NULL' the result is outputted to a tab-delimited
#' text file (without 'intervalID' column) rather than returned to the user.
#' This can be especially useful when the result is too big to fit into the
#' physical memory.  The resulted file can be used as an input for
#' 'gtrack.import' or 'gtrack.array.import' functions.
#'
#' If 'intervals.set.out' is not 'NULL' the result is saved as an intervals
#' set. Similarly to 'file' parameter 'intervals.set.out' can be useful to
#' overcome the limits of the physical memory.
#'
#' 'colnames' parameter controls the names of the columns that contain the
#' evaluated expressions. By default the column names match the track
#' expressions.
#'
#' @param ... track expression
#' @param intervals genomic scope for which the function is applied
#' @param colnames sets the columns names in the returned value. If 'NULL'
#' names are set to track expression.
#' @param iterator track expression iterator. If 'NULL' iterator is determined
#' implicitly based on track expressions.
#' @param band track expression band. If 'NULL' no band is used.
#' @param file file name where the function result is optionally outputted in
#' tab-delimited format
#' @param intervals.set.out intervals set name where the function result is
#' optionally outputted
#' @return If 'file' and 'intervals.set.out' are 'NULL' a set of intervals with
#' an additional column for each of the track expressions and 'columnID'
#' column.
#' @seealso \code{\link{gtrack.array.extract}}, \code{\link{gsample}},
#' \code{\link{gtrack.import}}, \code{\link{gtrack.array.import}},
#' \code{\link{glookup}}, \code{\link{gpartition}}, \code{\link{gdist}}
#' @keywords ~extract
#' @examples
#' \dontshow{
#' options(gmax.processes = 2)
#' }
#'
#' gdb.init_examples()
#'
#' ## get values of 'dense_track' for [0, 400), chrom 1
#' gextract("dense_track", gintervals(1, 0, 400))
#'
#' ## get values of 'rects_track' (a 2D track) for a 2D interval
#' gextract(
#'     "rects_track",
#'     gintervals.2d("chr1", 0, 4000, "chr2", 2000, 5000)
#' )
#'
#' @export gextract
gextract <- function(..., intervals = NULL, colnames = NULL, iterator = NULL, band = NULL, file = NULL, intervals.set.out = NULL) {
    args <- as.list(substitute(list(...)))[-1L]
    if (is.null(intervals) && length(args) < 2 || !is.null(intervals) && length(args) < 1) {
        stop("Usage: gextract([expr]+, intervals, colnames = NULL, iterator = NULL, band = NULL, file = NULL, intervals.set.out = NULL)", call. = FALSE)
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
                if (.ggetOption("gmultitasking")) {
                    res <- .gcall("gextract_multitask", intervals, tracks, colnames, .iterator, band, file, intervals.set.out, .misha_env())
                } else {
                    res <- .gcall("C_gextract", intervals, tracks, colnames, .iterator, band, file, intervals.set.out, .misha_env())
                }

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

    normalize_output <- function(x) {
        if (is.data.frame(x)) {
            return(.gnormalize_chrom_names(x))
        }
        if (is.list(x)) {
            return(lapply(x, normalize_output))
        }
        x
    }

    # Output normalization (chrom aliases) is only needed for per-chromosome DBs.
    # Indexed multi-contig stores canonical chrom names already; skip to avoid extra O(N) passes.
    if (!isTRUE(.ggetOption("gmulticontig.indexed_format"))) {
        res <- normalize_output(res)
    }

    # refresh the list of GINTERVS, etc.
    if (!is.null(intervals.set.out)) {
        .gdb.add_intervals.set(intervals.set.out)
        retv <- 0 # suppress return value
    } else if (!is.null(file)) {
        retv <- 0
    } # suppress return value
    else {
        res
    }
}


#' Calculates quantiles of a track expression
#'
#' Calculates the quantiles of a track expression for the given percentiles.
#'
#' This function calculates the quantiles for the given percentiles.
#'
#' If data size exceeds the limit (see: 'getOption(gmax.data.size)'), the data
#' is randomly sampled to fit the limit. A warning message is generated. The
#' seed of the pseudo-random generator can be controlled through 'grnd.seed'
#' option.
#'
#' Note: this function is capable to run in multitasking mode. Sampling may
#' vary according to the extent of multitasking. Since multitasking depends on
#' the number of available CPU cores, running the function on two different
#' machines might give different results. Please switch off multitasking if you
#' want to achieve identical results on any machine. For more information
#' regarding multitasking please refer "User Manual".
#'
#' @param expr track expression
#' @param percentiles an array of percentiles of quantiles in [0, 1] range
#' @param intervals genomic scope for which the function is applied
#' @param iterator track expression iterator. If 'NULL' iterator is determined
#' implicitly based on track expression.
#' @param band track expression band. If 'NULL' no band is used.
#' @return An array that represent quantiles.
#' @seealso \code{\link{gbins.quantiles}}, \code{\link{gintervals.quantiles}},
#' \code{\link{gdist}}
#' @keywords ~quantiles ~percentiles
#' @examples
#' \dontshow{
#' options(gmax.processes = 2)
#' }
#'
#' gdb.init_examples()
#' gquantiles("dense_track", c(0.1, 0.6, 0.8), gintervals(c(1, 2)))
#'
#' @export gquantiles

gquantiles <- function(expr = NULL, percentiles = 0.5, intervals = get("ALLGENOME", envir = .misha), iterator = NULL, band = NULL) {
    if (is.null(substitute(expr))) {
        stop("Usage: gquantiles(expr, percentiles = 0.5, intervals = .misha$ALLGENOME, iterator = NULL, band = NULL)", call. = FALSE)
    }
    .gcheckroot()

    intervals <- rescue_ALLGENOME(intervals, as.character(substitute(intervals)))

    exprstr <- do.call(.gexpr2str, list(substitute(expr)), envir = parent.frame())
    .iterator <- do.call(.giterator, list(substitute(iterator)), envir = parent.frame())

    if (.ggetOption("gmultitasking")) {
        res <- .gcall("gquantiles_multitask", intervals, exprstr, percentiles, .iterator, band, .misha_env())
    } else {
        res <- .gcall("C_gquantiles", intervals, exprstr, percentiles, .iterator, band, .misha_env())
    }
    res
}


#' Calculates summary statistics of track expression
#'
#' Calculates summary statistics of track expression.
#'
#' This function returns summary statistics of a track expression: total number
#' of bins, total number of bins whose value is NaN, min, max, sum, mean and
#' standard deviation of the values.
#'
#' @param expr track expression
#' @param intervals genomic scope for which the function is applied
#' @param iterator track expression iterator. If 'NULL' iterator is determined
#' implicitly based on track expression.
#' @param band track expression band. If 'NULL' no band is used.
#' @return An array that represents summary statistics.
#' @seealso \code{\link{gintervals.summary}}, \code{\link{gbins.summary}}
#' @keywords ~summary ~statistics
#' @examples
#' \dontshow{
#' options(gmax.processes = 2)
#' }
#'
#' gdb.init_examples()
#' gsummary("rects_track")
#'
#' @export gsummary
gsummary <- function(expr = NULL, intervals = NULL, iterator = NULL, band = NULL) {
    if (is.null(substitute(expr))) {
        stop("Usage: gsummary(expr, intervals = .misha$ALLGENOME, iterator = NULL, band = NULL)", call. = FALSE)
    }
    .gcheckroot()

    intervals <- rescue_ALLGENOME(intervals, as.character(substitute(intervals)))

    if (is.null(intervals)) {
        intervals <- get("ALLGENOME", envir = .misha)
    }

    exprstr <- do.call(.gexpr2str, list(substitute(expr)), envir = parent.frame())
    .iterator <- do.call(.giterator, list(substitute(iterator)), envir = parent.frame())

    if (.ggetOption("gmultitasking")) {
        res <- .gcall("gtracksummary_multitask", exprstr, intervals, .iterator, band, .misha_env())
    } else {
        res <- .gcall("gtracksummary", exprstr, intervals, .iterator, band, .misha_env())
    }
    res
}
