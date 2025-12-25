# Distribution functions (gcis_decay, gdist)
#' Calculates distribution of contact distances
#'
#' Calculates distribution of contact distances.
#'
#' A 2D iterator interval '(chrom1, start1, end1, chrom2, start2, end2)' is
#' said to represent a contact between two 1D intervals I1 and I2: '(chrom1,
#' start1, end1)' and '(chrom2, start2, end2)'.
#'
#' For contacts where 'chrom1' equals to 'chrom2' and I1 is within source
#' intervals the function calculates the distribution of distances between I1
#' and I2. The distribution is calculated separately for intra-domain and
#' inter-domain contacts.
#'
#' An interval is within source intervals if the unification of all source
#' intervals fully overlaps it. 'src' intervals are allowed to contain
#' overlapping intervals.
#'
#' Two intervals I1 and I2 are within the same domain (intra-domain contact) if
#' among the domain intervals exists an interval that fully overlaps both I1
#' and I2. Otherwise the contact is considered to be inter-domain. 'domain'
#' must contain only non-overlapping intervals.
#'
#' The distance between I1 and I2 is the absolute distance between the centers
#' of these intervals, i.e.: '|(start1 + end1 - start2 - end2) / 2|'.
#'
#' The range of distances for which the distribution is calculated is defined
#' by 'breaks' argument. For example: 'breaks=c(x1, x2, x3, x4)' represents
#' three different intervals (bins): (x1, x2], (x2, x3], (x3, x4].
#'
#' If 'include.lowest' is 'TRUE' the the lowest value will be included in the
#' first interval, i.e. in [x1, x2]
#'
#' @param expr track expression
#' @param breaks breaks that determine the bin
#' @param src source intervals
#' @param domain domain intervals
#' @param intervals genomic scope for which the function is applied
#' @param include.lowest if 'TRUE', the lowest value of the range determined by
#' breaks is included
#' @param iterator 2D track expression iterator. If 'NULL' iterator is
#' determined implicitly based on track expressions.
#' @param band track expression band. If 'NULL' no band is used.
#' @return 2-dimensional vector representing the distribution of contact
#' distances for inter and intra domains.
#' @seealso \code{\link{gdist}}, \code{\link{gtrack.2d.import_contacts}}
#' @keywords ~contacts
#' @examples
#' \dontshow{
#' options(gmax.processes = 2)
#' }
#'
#' gdb.init_examples()
#'
#' src <- rbind(
#'     gintervals(1, 10, 100),
#'     gintervals(1, 200, 300),
#'     gintervals(1, 400, 500),
#'     gintervals(1, 600, 700),
#'     gintervals(1, 7000, 9100),
#'     gintervals(1, 9000, 18000),
#'     gintervals(1, 30000, 31000),
#'     gintervals(2, 1130, 15000)
#' )
#'
#' domain <- rbind(
#'     gintervals(1, 0, 483000),
#'     gintervals(2, 0, 300000)
#' )
#'
#' gcis_decay("rects_track", 50000 * (1:10), src, domain)
#'
#' @export gcis_decay
gcis_decay <- function(expr = NULL, breaks = NULL, src = NULL, domain = NULL, intervals = NULL, include.lowest = FALSE, iterator = NULL, band = NULL) {
    if (is.null(substitute(expr)) || is.null(breaks) || is.null(src) || is.null(domain)) {
        stop("Usage: gcis_decay(expr, breaks, src, domain, intervals = .misha$ALLGENOME, include.lowest = FALSE, iterator = NULL, band = NULL)", call. = FALSE)
    }
    .gcheckroot()

    intervals <- rescue_ALLGENOME(intervals, as.character(substitute(intervals)))

    if (is.null(intervals)) {
        intervals <- get("ALLGENOME", envir = .misha)
    }

    exprstr <- do.call(.gexpr2str, list(substitute(expr)), envir = parent.frame())
    .iterator <- do.call(.giterator, list(substitute(iterator)), envir = parent.frame())

    src <- .gnormalize_chrom_names(src)
    domain <- .gnormalize_chrom_names(domain)
    intervals <- .gnormalize_chrom_names(intervals)

    res <- .gcall("C_gcis_decay", exprstr, breaks, src, domain, intervals, include.lowest, .iterator, band, .misha_env())
    attr(res, "breaks") <- breaks
    res
}


#' Calculates distribution of track expressions
#'
#' Calculates distribution of track expressions' values over the given set of
#' bins.
#'
#' This function calculates the distribution of values of the numeric track
#' expressions over the given set of bins.
#'
#' The range of bins is determined by 'breaks' argument. For example:
#' 'breaks=c(x1, x2, x3, x4)' represents three different intervals (bins): (x1,
#' x2], (x2, x3], (x3, x4].
#'
#' If 'include.lowest' is 'TRUE' the the lowest value will be included in the
#' first interval, i.e. in [x1, x2]
#'
#' 'gdist' can work with any number of dimensions. If more than one
#' 'expr'-'breaks' pair is passed, the result is a multidimensional vector, and
#' an individual value can be accessed by [i1,i2,...,iN] notation, where 'i1'
#' is the first track and 'iN' is the last track expression.
#'
#' @param ... pairs of 'expr', 'breaks' where 'expr' is a track expression and the breaks determine the bin
#' @param intervals genomic scope for which the function is applied
#' @param include.lowest if 'TRUE', the lowest value of the range determined by
#' breaks is included
#' @param iterator track expression iterator. If 'NULL' iterator is determined
#' implicitly based on track expressions.
#' @param band track expression band. If 'NULL' no band is used.
#' @param dataframe return a data frame instead of an N-dimensional vector.
#' @param names names for track expressions in the returned dataframe (only relevant when \code{dataframe == TRUE})
#' @return N-dimensional vector where N is the number of 'expr'-'breaks' pairs. If \code{dataframe == TRUE} - a data frame with a column for each track expression and an additional column 'n' with counts.
#' @seealso \code{\link{gextract}}
#' @keywords ~distribution
#' @examples
#' \dontshow{
#' options(gmax.processes = 2)
#' }
#'
#' gdb.init_examples()
#'
#' ## calculate the distribution of dense_track for bins:
#' ## (0, 0.2], (0.2, 0.5] and (0.5, 1]
#' gdist("dense_track", c(0, 0.2, 0.5, 1))
#'
#' ## calculate two-dimensional distribution:
#' ## dense_track vs. sparse_track
#' gdist("dense_track", seq(0, 1, by = 0.1), "sparse_track",
#'     seq(0, 2, by = 0.2),
#'     iterator = 100
#' )
#'
#' @export gdist
gdist <- function(..., intervals = NULL, include.lowest = FALSE, iterator = NULL, band = NULL, dataframe = FALSE, names = NULL) {
    args <- as.list(substitute(list(...)))[-1L]
    if (length(args) < 2 || (length(args) %% 2 != 0 && (length(args) - 1) %% 2 != 0)) {
        stop("Usage: gdist([expr, breaks]+, intervals = .misha$ALLGENOME, include.lowest = FALSE, iterator = NULL, band = NULL, dataframe = FALSE, names = NULL)", call. = FALSE)
    }
    .gcheckroot()

    intervals <- rescue_ALLGENOME(intervals, as.character(substitute(intervals)))

    if (length(args) %% 2 != 0) {
        intervals <- eval.parent(args[[length(args)]])
    } else if (is.null(intervals)) {
        intervals <- get("ALLGENOME", envir = .misha)
    }

    exprs <- c()
    breaks <- list()

    for (i in (0:(length(args) / 2 - 1))) {
        exprs <- append(exprs, do.call(.gexpr2str, list(args[[i * 2 + 1]]), envir = parent.frame()))
        breaks[[length(breaks) + 1]] <- eval.parent(args[[i * 2 + 2]])
    }

    .iterator <- do.call(.giterator, list(substitute(iterator)), envir = parent.frame())

    intervals <- .gnormalize_chrom_names(intervals)

    if (.ggetOption("gmultitasking")) {
        res <- .gcall("gtrackdist_multitask", intervals, exprs, breaks, include.lowest, .iterator, band, .misha_env())
    } else {
        res <- .gcall("gtrackdist", intervals, exprs, breaks, include.lowest, .iterator, band, .misha_env())
    }
    attr(res, "breaks") <- breaks

    if (dataframe) {
        res <- as.data.frame.table(res)
        if (is.null(names)) {
            names <- exprs
        }
        colnames(res) <- c(names, "n")
    }

    res
}
