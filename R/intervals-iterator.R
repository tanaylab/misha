# Iterator functions

#' Creates a cartesian-grid iterator
#'
#' Creates a cartesian grid two-dimensional iterator that can be used by any
#' function that accepts an iterator argument.
#'
#' This function creates and returns a cartesian grid two-dimensional iterator
#' that can be used by any function that accepts an iterator argument.
#'
#' Assume 'centers1' and 'centers2' to be the central points of each interval
#' from 'intervals1' and 'intervals2', and 'C1', 'C2' to be two points from
#' 'centers1', 'centers2' accordingly. Assume also that the values in
#' 'expansion1' and 'expansion2' are unique and sorted.
#'
#' 'giterator.cartesian_grid' creates a set of all possible unique and
#' non-overlapping two-dimensional intervals of form: '(chrom1, start1, end1,
#' chrom2, start2, end2)'. Each '(chrom1, start1, end1)' is created by taking a
#' point 'C1' - '(chrom1, coord1)' and converting it to 'start1' and 'end1'
#' such that 'start1 == coord1+E1[i]', 'end1 == coord1+E1[i+1]', where 'E1[i]'
#' is one of the sorted 'expansion1' values. Overlaps between rectangles or
#' expansion beyond the limits of chromosome are avoided.
#'
#' 'min.band.idx' and 'max.band.idx' parameters control whether a pair of 'C1'
#' and 'C2' is skipped or not. If both of these parameters are not 'NULL' AND
#' if both 'C1' and 'C2' share the same chromosome AND the delta of indices of
#' 'C1' and 'C2' ('C1 index - C2 index') lays within '[min.band.idx,
#' max.band.idx]' range - only then the pair will be used to create the
#' intervals. Otherwise 'C1-C2' pair is filtered out. Note: if 'min.band.idx'
#' and 'max.band.idx' are not 'NULL', i.e. band indices filtering is applied,
#' then 'intervals2' parameter must be set to 'NULL'.
#'
#' @param intervals1 one-dimensional intervals
#' @param expansion1 an array of integers that define expansion around
#' intervals1 centers
#' @param intervals2 one-dimensional intervals. If 'NULL' then 'intervals2' is
#' considered to be equal to 'intervals1'
#' @param expansion2 an array of integers that define expansion around
#' intervals2 centers. If 'NULL' then 'expansion2' is considered to be equal to
#' 'expansion1'
#' @param min.band.idx,max.band.idx integers that limit iterator intervals to
#' band
#' @return A list containing the definition of cartesian iterator.
#' @seealso \code{\link{giterator.intervals}}
#' @keywords ~iterator ~cartesian
#' @examples
#' \dontshow{
#' options(gmax.processes = 2)
#' }
#'
#' gdb.init_examples()
#'
#' intervs1 <- gintervals(
#'     c(1, 1, 2), c(100, 300, 200),
#'     c(300, 500, 300)
#' )
#' intervs2 <- gintervals(
#'     c(1, 2, 2), c(400, 1000, 3000),
#'     c(800, 2000, 4000)
#' )
#' itr <- giterator.cartesian_grid(
#'     intervs1, c(-20, 100), intervs2,
#'     c(-40, -10, 50)
#' )
#' giterator.intervals(iterator = itr)
#'
#' itr <- giterator.cartesian_grid(intervs1, c(-20, 50, 100))
#' giterator.intervals(iterator = itr)
#'
#' itr <- giterator.cartesian_grid(intervs1, c(-20, 50, 100),
#'     min.band.idx = -1,
#'     max.band.idx = 0
#' )
#' giterator.intervals(iterator = itr)
#'
#' @export giterator.cartesian_grid
giterator.cartesian_grid <- function(intervals1 = NULL, expansion1 = NULL, intervals2 = NULL, expansion2 = NULL, min.band.idx = NULL, max.band.idx = NULL) {
    if (is.null(intervals1) || is.null(expansion1)) {
        stop("Usage: giterator.cartesian_grid(intervals1, expansion1, intervals2 = NULL, expansion2 = NULL, min.band.idx = NULL, max.band.idx = NULL)", call. = FALSE)
    }

    use.band.idx.limit <- !is.null(min.band.idx) && !is.null(max.band.idx)
    if (use.band.idx.limit) {
        if (min.band.idx > max.band.idx) {
            stop("min.band.idx exceeds max.band.idx", call. = FALSE)
        }

        if (!is.null(intervals2)) {
            stop("band.idx limit can only be used when intervals2 is set to NULL", call. = FALSE)
        }
    } else {
        min.band.idx <- 0
        max.band.idx <- 0
    }

    r <- list(
        intervals1 = intervals1, intervals2 = intervals2, expansion1 = expansion1, expansion2 = expansion2,
        band.idx = c(min.band.idx, max.band.idx, use.band.idx.limit)
    )
    class(r) <- "cartesian.grid"
    .gcall("gcheck_iterator", r, .misha_env())
    r
}


#' Returns iterator intervals
#'
#' Returns iterator intervals given track expression, scope, iterator and band.
#'
#' This function returns a set of intervals used by the iterator intervals for
#' the given track expression, genomic scope, iterator and band. Some functions
#' accept an iterator without accepting a track expression (like
#' 'gtrack.create_pwm_energy'). These functions generate the values for each
#' iterator interval by themselves. Use set 'expr' to 'NULL' to simulate the
#' work of these functions.
#'
#' If 'intervals.set.out' is not 'NULL' the result is saved as an intervals
#' set. Use this parameter if the result size exceeds the limits of the
#' physical memory.
#'
#' When 'interval_relative' is TRUE, bins are aligned to each input interval's
#' start position rather than chromosome position 0. This mode requires a
#' numeric iterator (binsize) and returns an additional 'intervalID' column
#' indicating which input interval spawned each bin.
#'
#' @param expr track expression
#' @param intervals genomic scope
#' @param iterator track expression iterator. If 'NULL' iterator is determined
#' implicitly based on track expression.
#' @param band track expression band. If 'NULL' no band is used.
#' @param intervals.set.out intervals set name where the function result is
#' optionally outputted
#' @param interval_relative if TRUE, and iterator is numeric, bins start at each interval's start
#' position instead of chromosome position 0. Returns intervalID column. Default: FALSE.
#' @param partial_bins how to handle partial bins at interval boundaries when
#' interval_relative is TRUE. One of "clip" (default, truncate last bin to
#' interval boundary), "exact" or "drop" (only output full-size bins).
#' @return If 'intervals.set.out' is 'NULL' a data frame representing iterator intervals.
#' When 'interval_relative' is TRUE, includes an 'intervalID' column.
#' @seealso \code{\link{giterator.cartesian_grid}}
#' @keywords ~iterator ~intervals
#' @examples
#' \dontshow{
#' options(gmax.processes = 2)
#' }
#'
#' gdb.init_examples()
#'
#' ## iterator is set implicitly to bin size of 'dense' track
#' giterator.intervals("dense_track", gintervals(1, 0, 200))
#'
#' ## iterator = 30
#' giterator.intervals("dense_track", gintervals(1, 0, 200), 30)
#'
#' ## iterator is an intervals set named 'annotations'
#' giterator.intervals("dense_track", .misha$ALLGENOME, "annotations")
#'
#' ## iterator is set implicitly to intervals of 'array_track' track
#' giterator.intervals("array_track", gintervals(1, 0, 200))
#'
#' ## iterator is a rectangle 100000 by 50000
#' giterator.intervals(
#'     "rects_track",
#'     gintervals.2d(chroms1 = 1, chroms2 = "chrX"),
#'     c(100000, 50000)
#' )
#'
#' ## interval_relative mode: bins aligned to each interval's start
#' intervs <- gintervals(1, c(100, 500), c(300, 700))
#' giterator.intervals(NULL, intervs, iterator = 50, interval_relative = TRUE)
#'
#' @export giterator.intervals
giterator.intervals <- function(expr = NULL, intervals = .misha$ALLGENOME, iterator = NULL, band = NULL,
                                intervals.set.out = NULL, interval_relative = FALSE,
                                partial_bins = c("clip", "exact", "drop")) {
    if (is.null(substitute(expr)) && is.null(substitute(iterator))) {
        stop("Usage: giterator.intervals(expr = NULL, intervals = .misha$ALLGENOME, iterator = NULL, band = NULL, intervals.set.out = NULL, interval_relative = FALSE, partial_bins = 'clip')", call. = FALSE)
    }

    intervals <- rescue_ALLGENOME(intervals, as.character(substitute(intervals)))

    if (is.null(substitute(expr))) {
        exprstr <- "0"
    } else {
        exprstr <- do.call(.gexpr2str, list(substitute(expr)), envir = parent.frame())
    }

    .iterator <- do.call(.giterator, list(substitute(iterator)), envir = parent.frame())
    intervals.set.out <- do.call(.gexpr2str, list(substitute(intervals.set.out)), envir = parent.frame())

    # Validate and convert partial_bins to integer code
    partial_bins <- match.arg(partial_bins)
    partial_bins_code <- match(partial_bins, c("clip", "exact", "drop")) - 1L # 0=clip, 1=exact, 2=drop
    if (partial_bins_code == 2L) partial_bins_code <- 1L # "drop" is alias for "exact"

    # Validate interval_relative mode
    if (interval_relative) {
        if (!is.numeric(.iterator) || length(.iterator) != 1) {
            stop("interval_relative mode requires a numeric iterator (binsize)", call. = FALSE)
        }
    }

    if (!is.null(intervals.set.out)) {
        fullpath <- .gintervals.check_new_set(intervals.set.out)
    }

    # intervals can be NULL if gextract is piped with gscreen and the latter returns NULL
    success <- FALSE
    res <- NULL
    tryCatch(
        {
            if (!is.null(intervals)) {
                res <- .gcall(
                    "giterator_intervals", exprstr, intervals, .iterator, band, intervals.set.out,
                    interval_relative, as.integer(partial_bins_code), .misha_env()
                )

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
