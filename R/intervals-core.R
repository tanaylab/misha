# Core interval creation and validation functions

.gintervals <- function(chroms, starts, ends, strands) {
    if (is.null(strands)) {
        intervals <- data.frame(chrom = .gchroms(chroms), start = starts, end = ends)
    } else {
        intervals <- data.frame(chrom = .gchroms(chroms), start = starts, end = ends, strand = strands)
    }

    numintervals <- nrow(intervals)

    maxends <- get("ALLGENOME", envir = .misha)[[1]]$end[match(as.character(intervals$chrom), get("ALLGENOME", envir = .misha)[[1]]$chrom)]
    maxidx <- intervals$end == -1
    intervals$end[maxidx] <- maxends[maxidx]

    err.intervs <- intervals[intervals$start < 0, ]
    if (nrow(err.intervs) > 0) {
        stop(sprintf("Invalid interval (%s, %g, %g): start coordinate is out of range", err.intervs$chrom[1], err.intervs$start[1], err.intervs$end[1]), call. = FALSE)
    }

    err.intervs <- intervals[intervals$end > maxends, ]
    if (nrow(err.intervs) > 0) {
        stop(sprintf("Invalid interval (%s, %g, %g): end coordinate exceeds chromosome boundaries", err.intervs$chrom[1], err.intervs$start[1], err.intervs$end[1]), call. = FALSE)
    }

    err.intervs <- intervals[intervals$start >= intervals$end, ]
    if (nrow(err.intervs) > 0) {
        stop(sprintf("Invalid interval (%s, %g, %g): start coordinate exceeds or equals to end coordinate", err.intervs$chrom[1], err.intervs$start[1], err.intervs$end[1]), call. = FALSE)
    }

    if (!is.null(strands)) {
        if (!is.numeric(intervals$strand)) {
            stop("Invalid strand values", call. = FALSE)
        }

        err.intervs <- intervals[intervals$strand != as.integer(intervals$strand) | intervals$strand < -1 | intervals$strand > 1, ]
        if (nrow(err.intervs) > 0) {
            stop(sprintf("Invalid strand value %g of interval (%s, %g, %g)", err.intervs$strand[1], err.intervs$chrom[1], err.intervs$start[1], err.intervs$end[1]))
        }
    }

    intervals
}

.gintervals.is1d <- function(intervals) {
    if (is.character(intervals)) {
        if (.gintervals.is_bigset(intervals)) {
            return(.gintervals.big.is1d(intervals))
        }
        intervals <- .gintervals.load(intervals)
    }
    all(colnames(intervals)[1:3] == c("chrom", "start", "end"))
}

.gintervals.is2d <- function(intervals) {
    if (is.character(intervals)) {
        if (.gintervals.is_bigset(intervals)) {
            return(.gintervals.big.is2d(intervals))
        }
        intervals <- .gintervals.load(intervals)
    }
    all(colnames(intervals)[1:6] == c("chrom1", "start1", "end1", "chrom2", "start2", "end2"))
}


#' Creates a set of 1D intervals
#'
#' Creates a set of 1D intervals.
#'
#' This function returns a set of one-dimensional intervals. The returned value
#' can be used in all functions that accept 'intervals' argument.
#'
#' One-dimensional intervals is a data frame whose first three columns are
#' 'chrom', 'start' and 'end'. Each row of the data frame represents a genomic
#' interval of the specified chromosome in the range of [start, end).
#' Additional columns can be presented in 1D intervals object yet these columns
#' must be added after the three obligatory ones.
#'
#' If 'strands' argument is not 'NULL' an additional column "strand" is added
#' to the intervals. The possible values of a strand can be '1' (plus strand),
#' '-1' (minus strand) or '0' (unknown).
#'
#' @param chroms chromosomes - an array of strings with or without "chr"
#' prefixes or an array of integers (like: '1' for "chr1")
#' @param starts an array of start coordinates
#' @param ends an array of end coordinates. If '-1' chromosome size is assumed.
#' @param strands 'NULL' or an array consisting of '-1', '0' or '1' values
#' @return A data frame representing the intervals.
#' @seealso \code{\link{gintervals.2d}}, \code{\link{gintervals.force_range}}
#' @keywords ~intervals
#' @examples
#' \dontshow{
#' options(gmax.processes = 2)
#' }
#'
#' gdb.init_examples()
#'
#' ## the following 3 calls produce identical results
#' gintervals(1)
#' gintervals("1")
#' gintervals("chrX")
#'
#' gintervals(1, 1000)
#' gintervals(c("chr2", "chrX"), 10, c(3000, 5000))
#'
#' @export gintervals
gintervals <- function(chroms = NULL, starts = 0, ends = -1, strands = NULL) {
    if (is.null(chroms)) {
        stop("Usage: gintervals(chroms, starts = 0, ends = -1, strands = NULL)", call. = FALSE)
    }
    .gcheckroot()

    intervals <- .gintervals(chroms, starts, ends, strands)
    .gcall("gintervsort", intervals, .misha_env())
}


#' Creates a set of 2D intervals
#'
#' Creates a set of 2D intervals.
#'
#' This function returns a set of two-dimensional intervals. The returned value
#' can be used in all functions that accept 'intervals' argument.
#'
#' Two-dimensional intervals is a data frame whose first six columns are
#' 'chrom1', 'start1', 'end1', 'chrom2', 'start2' and 'end2'. Each row of the
#' data frame represents two genomic intervals from two chromosomes in the
#' range of [start, end). Additional columns can be presented in 2D intervals
#' object yet these columns must be added after the six obligatory ones.
#'
#' @param chroms1 chromosomes1 - an array of strings with or without "chr"
#' prefixes or an array of integers (like: '1' for "chr1")
#' @param starts1 an array of start1 coordinates
#' @param ends1 an array of end1 coordinates. If '-1' chromosome size is
#' assumed.
#' @param chroms2 chromosomes2 - an array of strings with or without "chr"
#' prefixes or an array of integers (like: '1' for "chr1"). If 'NULL',
#' 'chroms2' is assumed to be equal to 'chroms1'.
#' @param starts2 an array of start2 coordinates
#' @param ends2 an array of end2 coordinates. If '-1' chromosome size is
#' assumed.
#' @return A data frame representing the intervals.
#' @seealso \code{\link{gintervals}}, \code{\link{gintervals.force_range}}
#' @keywords ~intervals
#' @examples
#' \dontshow{
#' options(gmax.processes = 2)
#' }
#'
#' gdb.init_examples()
#'
#' ## the following 3 calls produce identical results
#' gintervals.2d(1)
#' gintervals.2d("1")
#' gintervals.2d("chrX")
#'
#' gintervals.2d(1, 1000, 2000, "chrX", 400, 800)
#' gintervals.2d(c("chr2", "chrX"), 10, c(3000, 5000), 1)
#'
#' @export gintervals.2d
gintervals.2d <- function(chroms1 = NULL, starts1 = 0, ends1 = -1, chroms2 = NULL, starts2 = 0, ends2 = -1) {
    if (is.null(chroms1)) {
        stop("Usage: gintervals.2d(chroms1, starts1 = 0, ends1 = -1, chroms2 = NULL, starts2 = 0, ends2 = -1)", call. = FALSE)
    }
    .gcheckroot()

    if (is.null(chroms2)) {
        chroms2 <- chroms1
    }

    intervals1 <- .gintervals(chroms1, starts1, ends1, NULL)
    intervals2 <- .gintervals(chroms2, starts2, ends2, NULL)

    intervals <- data.frame(
        chrom1 = intervals1$chrom, start1 = intervals1$start, end1 = intervals1$end,
        chrom2 = intervals2$chrom, start2 = intervals2$start, end2 = intervals2$end
    )

    .gcall("gintervsort", intervals, .misha_env())
}


#' Returns 2D intervals that cover the whole genome
#'
#' Returns 2D intervals that cover the whole genome.
#'
#' This function returns a set of two-dimensional intervals that cover the
#' whole genome as it is defined by 'chrom_sizes.txt' file.
#'
#' @return A data frame representing the intervals.
#' @seealso \code{\link{gintervals.2d}}
#' @keywords ~genome ~chromosome ~chromosomes ~ALLGENOME
#' @export gintervals.2d.all
gintervals.2d.all <- function() {
    .gcheckroot()
    intervals2d <- get("ALLGENOME", envir = .misha)[[2]]

    # Check if 2D is deferred and generate on demand
    if (.is_2d_deferred(intervals2d)) {
        mode <- getOption("gmulticontig.2d.mode", "diagonal")
        n_contigs <- attr(intervals2d, "n_contigs")

        if (mode == "full") {
            warning(sprintf(
                "Generating full 2D genome with %d contigs (%d pairs). This may take time and use significant memory.",
                n_contigs, n_contigs * n_contigs
            ))
        }

        intervals <- get("ALLGENOME", envir = .misha)[[1]]
        intervals2d <- .generate_2d_on_demand(intervals, mode)

        # Cache if requested
        if (getOption("gmulticontig.2d.cache", FALSE)) {
            allgenome <- get("ALLGENOME", envir = .misha)
            allgenome[[2]] <- intervals2d
            assign("ALLGENOME", allgenome, envir = .misha)
        }
    }

    intervals2d
}


#' Returns 1D intervals that cover the whole genome
#'
#' Returns 1D intervals that cover the whole genome.
#'
#' This function returns a set of one-dimensional intervals that cover the
#' whole genome as it is defined by 'chrom_sizes.txt' file.
#'
#' @return A data frame representing the intervals.
#' @seealso \code{\link{gintervals}}
#' @keywords ~genome ~chromosome ~chromosomes ~ALLGENOME
#' @export gintervals.all
gintervals.all <- function() {
    .gcheckroot()
    get("ALLGENOME", envir = .misha)[[1]]
}
