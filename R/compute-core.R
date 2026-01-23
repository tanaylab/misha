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

#' Calculates correlation between track expressions
#'
#' Calculates correlation between track expressions over iterator bins
#' inside the supplied genomic scope. Expressions are processed in pairs:
#' (expr1, expr2), (expr3, expr4), etc. Only bins where both expressions are
#' not NaN are used.
#'
#' @param expr1 first track expression
#' @param expr2 second track expression
#' @param ... additional track expressions, supplied as pairs (expr3, expr4, ...)
#' @param intervals genomic scope for which the function is applied
#' @param iterator track expression iterator. If 'NULL' iterator is determined
#' implicitly based on track expression.
#' @param band track expression band. If 'NULL' no band is used.
#' @param method correlation method to use. One of 'pearson' (default),
#' 'spearman' (approximate, memory-efficient), or 'spearman.exact' (exact,
#' requires O(n) memory where n is number of non-NaN pairs).
#' @param details if 'TRUE' returns summary statistics for each pair, otherwise
#' returns correlations only. For Pearson, includes n, n.na, mean1, mean2, sd1,
#' sd2, cov, cor. For Spearman methods, includes n, n.na, cor.
#' @param names optional names for the pairs. If supplied, length must match the
#' number of pairs.
#' @return If 'details' is 'FALSE', a numeric vector of correlations. If
#' 'details' is 'TRUE', a data frame with summary statistics for each pair.
#' @seealso \code{\link{gextract}}, \code{\link{gscreen}}, \code{\link{gsummary}}
#' @examples
#' \dontshow{
#' options(gmax.processes = 2)
#' }
#'
#' gdb.init_examples()
#' gcor("dense_track", "sparse_track", intervals = gintervals(1, 0, 10000), iterator = 1000)
#'
#' # Spearman correlation (approximate, memory-efficient)
#' gcor("dense_track", "sparse_track",
#'     intervals = gintervals(1, 0, 10000),
#'     iterator = 1000, method = "spearman"
#' )
#'
#' # Exact Spearman correlation
#' gcor("dense_track", "sparse_track",
#'     intervals = gintervals(1, 0, 10000),
#'     iterator = 1000, method = "spearman.exact"
#' )
#'
#' @export gcor
gcor <- function(expr1 = NULL, expr2 = NULL, ..., intervals = NULL, iterator = NULL, band = NULL, method = c("pearson", "spearman", "spearman.exact"), details = FALSE, names = NULL) {
    if (is.null(substitute(expr1)) || is.null(substitute(expr2))) {
        stop("Usage: gcor(expr1, expr2, ..., intervals = .misha$ALLGENOME, iterator = NULL, band = NULL, method = 'pearson', details = FALSE, names = NULL)", call. = FALSE)
    }
    .gcheckroot()

    intervals <- rescue_ALLGENOME(intervals, as.character(substitute(intervals)))

    eval_env <- parent.frame()
    args <- c(list(substitute(expr1)), list(substitute(expr2)), as.list(substitute(list(...)))[-1L])

    is_intervals_candidate <- function(x) {
        if (is.data.frame(x)) {
            return(.gintervals.is1d(x) || .gintervals.is2d(x))
        }
        if (is.character(x) && length(x) == 1) {
            gintervs <- get("GINTERVS", envir = .misha)
            if (x %in% gintervs) {
                return(TRUE)
            }
            if (.gintervals.is_bigset(x, FALSE)) {
                return(TRUE)
            }
        }
        FALSE
    }

    # Check for positional intervals (last argument) before defaulting to ALLGENOME
    if (is.null(intervals) && length(args) %% 2 != 0) {
        intervals_candidate <- eval(args[[length(args)]], eval_env)
        if (is_intervals_candidate(intervals_candidate)) {
            intervals <- intervals_candidate
            args <- args[-length(args)]
        } else {
            stop("gcor expects an even number of track expressions (pairs).", call. = FALSE)
        }
    }

    if (length(args) %% 2 != 0) {
        stop("gcor expects an even number of track expressions (pairs).", call. = FALSE)
    }

    # Default to ALLGENOME after positional intervals check
    if (is.null(intervals)) {
        intervals <- get("ALLGENOME", envir = .misha)
    }

    exprs <- vapply(args, function(arg) {
        do.call(.gexpr2str, list(arg), envir = eval_env)
    }, character(1))
    .iterator <- do.call(.giterator, list(substitute(iterator)), envir = eval_env)
    num_pairs <- length(exprs) / 2

    if (!is.null(names) && length(names) != num_pairs) {
        stop("names length must match the number of expression pairs.", call. = FALSE)
    }

    method <- match.arg(method)

    # Select C++ function based on method
    if (method == "pearson") {
        if (.ggetOption("gmultitasking")) {
            res <- .gcall("gtrackcor_multitask", exprs, intervals, .iterator, band, .misha_env())
        } else {
            res <- .gcall("gtrackcor", exprs, intervals, .iterator, band, .misha_env())
        }
    } else if (method == "spearman") {
        if (.ggetOption("gmultitasking")) {
            res <- .gcall("gtrackcor_spearman_multitask", exprs, intervals, .iterator, band, .misha_env())
        } else {
            res <- .gcall("gtrackcor_spearman", exprs, intervals, .iterator, band, .misha_env())
        }
    } else if (method == "spearman.exact") {
        if (.ggetOption("gmultitasking")) {
            res <- .gcall("gtrackcor_spearman_exact_multitask", exprs, intervals, .iterator, band, .misha_env())
        } else {
            res <- .gcall("gtrackcor_spearman_exact", exprs, intervals, .iterator, band, .misha_env())
        }
    }

    if (is.null(dim(res))) {
        stats_matrix <- matrix(res, nrow = 1, dimnames = list(NULL, names(res)))
    } else {
        stats_matrix <- res
    }

    pair_names <- if (is.null(names)) {
        paste(exprs[seq(1, length(exprs), by = 2)], exprs[seq(2, length(exprs), by = 2)], sep = "~")
    } else {
        names
    }

    if (isTRUE(details)) {
        stats_df <- as.data.frame(stats_matrix, stringsAsFactors = FALSE)
        rownames(stats_df) <- pair_names
        return(stats_df)
    }

    cor_vals <- stats_matrix[, "cor"]
    names(cor_vals) <- pair_names
    cor_vals
}
