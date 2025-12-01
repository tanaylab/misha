# Utility functions (random, coverage, etc.)

#' Generate random genome intervals
#'
#' Generate random genome intervals with a specified number of regions of a specified size.
#' This function samples intervals uniformly across the genome, weighted by chromosome length.
#'
#' @param size The size of the intervals to generate (in base pairs)
#' @param n The number of intervals to generate
#' @param dist_from_edge The minimum distance from the edge of the chromosome for a region to start or end (default: 3e6)
#' @param chromosomes The chromosomes to sample from (default: all chromosomes). Can be a character vector of chromosome names.
#' @param filter A set of intervals to exclude from sampling (default: NULL). Generated intervals will not overlap with these regions.
#'
#' @return A data.frame with columns chrom, start, and end representing genomic intervals
#'
#' @details
#' The function samples intervals randomly across the genome, with chromosomes weighted by their length.
#' Each interval is guaranteed to:
#' \itemize{
#'   \item Be of the specified size
#'   \item Start and end at least \code{dist_from_edge} bases away from chromosome boundaries
#'   \item Fall entirely within a single chromosome
#'   \item Not overlap with any intervals in the \code{filter} (if provided)
#' }
#'
#' When a filter is provided, the function pre-computes valid genome segments (regions not in the filter)
#' and samples from these segments. Note that this can be slow
#' when the filter contains many intervals.
#'
#' The function uses R's random number generator, so \code{set.seed()} can be used for reproducibility.
#'
#' This function is implemented in C++ for high performance and can generate millions of intervals quickly.
#'
#' @examples
#' \dontrun{
#' gdb.init_examples()
#'
#' # Generate 1000 random intervals of 100bp
#' intervals <- gintervals.random(100, 1000)
#' head(intervals)
#'
#' # Generate intervals only on chr1 and chr2
#' intervals <- gintervals.random(100, 1000, chromosomes = c("chr1", "chr2"))
#'
#' # Generate intervals avoiding specific regions
#' filter_regions <- gintervals(c("chr1", "chr2"), c(1000, 5000), c(2000, 6000))
#' intervals <- gintervals.random(100, 1000, filter = filter_regions)
#'
#' # Verify no overlaps with filter
#' overlaps <- gintervals.intersect(intervals, filter_regions)
#' nrow(overlaps) # Should be 0
#'
#' # For reproducibility
#' set.seed(123)
#' intervals1 <- gintervals.random(100, 100)
#' set.seed(123)
#' intervals2 <- gintervals.random(100, 100)
#' identical(intervals1, intervals2) # TRUE
#' }
#'
#' @export
gintervals.random <- function(size, n, dist_from_edge = 3e6, chromosomes = NULL, filter = NULL) {
    # Check that database is initialized
    .gcheckroot()

    # Validate inputs
    if (!is.numeric(size) || length(size) != 1 || size <= 0) {
        stop("size must be a positive number", call. = FALSE)
    }
    if (!is.numeric(n) || length(n) != 1 || n <= 0) {
        stop("n must be a positive number", call. = FALSE)
    }
    if (!is.numeric(dist_from_edge) || length(dist_from_edge) != 1 || dist_from_edge < 0) {
        stop("dist_from_edge must be a non-negative number", call. = FALSE)
    }

    # Validate filter if provided
    if (!is.null(filter)) {
        if (!is.data.frame(filter)) {
            stop("filter must be a data frame", call. = FALSE)
        }
        if (!all(c("chrom", "start", "end") %in% names(filter))) {
            stop("filter must have columns: chrom, start, end", call. = FALSE)
        }
        if (nrow(filter) > 0) {
            # Validate filter intervals
            if (any(filter$start < 0)) {
                stop("filter intervals must have start >= 0", call. = FALSE)
            }
            if (any(filter$start >= filter$end)) {
                stop("filter intervals must have start < end", call. = FALSE)
            }
            # Sort and unify overlapping filter intervals for efficiency
            filter <- filter[order(filter$chrom, filter$start), ]
            filter <- gintervals.canonic(filter)
        } else {
            filter <- NULL # Empty filter same as no filter
        }
    }

    # Get all chromosomes
    all_genome <- gintervals.all()

    # Filter by chromosomes if specified
    if (!is.null(chromosomes)) {
        if (!is.character(chromosomes)) {
            stop("chromosomes must be a character vector", call. = FALSE)
        }
        all_genome <- all_genome[all_genome$chrom %in% chromosomes, , drop = FALSE]
        if (nrow(all_genome) == 0) {
            stop("No chromosomes named ", paste(chromosomes, collapse = ", "), " found in the genome", call. = FALSE)
        }
        # Also filter the filter intervals to only include selected chromosomes
        if (!is.null(filter)) {
            filter <- filter[filter$chrom %in% chromosomes, , drop = FALSE]
            if (nrow(filter) == 0) {
                filter <- NULL
            }
        }
    }

    # Pre-filter: remove chromosomes that are too short (only if no filter)
    # With filter, C++ will handle this more intelligently
    if (is.null(filter)) {
        chrom_lengths <- all_genome$end - all_genome$start
        min_required_length <- size + 2 * dist_from_edge
        valid_chroms <- chrom_lengths >= min_required_length

        if (!any(valid_chroms)) {
            stop("No chromosomes are long enough for intervals of size ", size,
                " with dist_from_edge ", dist_from_edge,
                " (minimum required chromosome length: ", min_required_length, ")",
                call. = FALSE
            )
        }

        all_genome <- all_genome[valid_chroms, , drop = FALSE]
    }

    # Call C++ function
    result <- .Call("C_grandom_genome",
        as.integer(size),
        as.integer(n),
        as.numeric(dist_from_edge),
        all_genome,
        filter,
        PACKAGE = "misha"
    )

    # Force range to ensure intervals are within chromosome boundaries
    # (should already be satisfied, but this is a safety check)
    result <- gintervals.force_range(result)

    return(result)
}


#' Calculate total base pairs covered by intervals
#'
#' Returns the total number of base pairs covered by a set of intervals.
#'
#' This function first canonicalizes the intervals to remove overlaps and
#' touching intervals, then sums up the lengths of all resulting intervals.
#' Overlapping intervals are counted only once.
#'
#' @param intervals set of one-dimensional intervals
#' @return A single numeric value representing the total number of base pairs
#' covered by the intervals.
#' @seealso \code{\link{gintervals}}, \code{\link{gintervals.canonic}},
#' \code{\link{gintervals.coverage_fraction}}
#' @keywords ~coverage ~genomics
#' @examples
#' \dontshow{
#' options(gmax.processes = 2)
#' }
#'
#' gdb.init_examples()
#'
#' # Create some intervals
#' intervs <- gintervals(
#'     c("chr1", "chr1", "chr2"),
#'     c(100, 150, 1000),
#'     c(200, 250, 2000)
#' )
#'
#' # Calculate total bp covered
#' # Note: intervals [100,200) and [150,250) overlap,
#' # so total is (200-100) + (250-150) + (2000-1000) = 100 + 100 + 1000 = 1200
#' # But after canonicalization: [100,250) + [1000,2000) = 150 + 1000 = 1150
#' gintervals.covered_bp(intervs)
#'
#' @export gintervals.covered_bp
gintervals.covered_bp <- function(intervals = NULL) {
    if (is.null(intervals)) {
        stop("Usage: gintervals.covered_bp(intervals)", call. = FALSE)
    }

    # Handle empty intervals before rescue_ALLGENOME
    if (is.data.frame(intervals) && nrow(intervals) == 0) {
        return(0)
    }

    intervals <- rescue_ALLGENOME(intervals, as.character(substitute(intervals)))

    # Canonicalize to remove overlaps
    canonical <- gintervals.canonic(intervals, unify_touching_intervals = TRUE)

    if (is.null(canonical) || nrow(canonical) == 0) {
        return(0)
    }

    # Sum up the lengths
    sum(canonical$end - canonical$start)
}


#' Calculate fraction of genomic space covered by intervals
#'
#' Returns the fraction of a genomic space that is covered by a set of intervals.
#'
#' This function calculates what fraction of 'intervals2' is covered by
#' 'intervals1'. If 'intervals2' is NULL, it calculates the fraction of the
#' entire genome that is covered by 'intervals1'. Overlapping intervals in
#' either set are automatically unified before calculation.
#'
#' @param intervals1 set of one-dimensional intervals (the covering set)
#' @param intervals2 set of one-dimensional intervals to be covered (default:
#' NULL, meaning the entire genome)
#' @return A single numeric value between 0 and 1 representing the fraction of
#' 'intervals2' (or the genome) covered by 'intervals1'.
#' @seealso \code{\link{gintervals}}, \code{\link{gintervals.intersect}},
#' \code{\link{gintervals.covered_bp}}, \code{\link{gintervals.all}}
#' @keywords ~coverage ~genomics
#' @examples
#' \dontshow{
#' options(gmax.processes = 2)
#' }
#'
#' gdb.init_examples()
#'
#' # Create some intervals
#' intervs1 <- gscreen("dense_track > 0.15")
#' intervs2 <- gintervals(c("chr1", "chr2"), 0, c(100000, 100000))
#'
#' # Calculate fraction of intervs2 covered by intervs1
#' gintervals.coverage_fraction(intervs1, intervs2)
#'
#' # Calculate fraction of entire genome covered by intervs1
#' gintervals.coverage_fraction(intervs1)
#'
#' @export gintervals.coverage_fraction
gintervals.coverage_fraction <- function(intervals1 = NULL, intervals2 = NULL) {
    if (is.null(intervals1)) {
        stop("Usage: gintervals.coverage_fraction(intervals1, intervals2 = NULL)", call. = FALSE)
    }

    # Handle empty intervals1 before rescue_ALLGENOME
    if (is.data.frame(intervals1) && nrow(intervals1) == 0) {
        return(0)
    }

    intervals1 <- rescue_ALLGENOME(intervals1, as.character(substitute(intervals1)))

    # If intervals2 is NULL, use entire genome
    if (is.null(intervals2)) {
        intervals2 <- gintervals.all()
    } else {
        # Handle empty intervals2 before rescue_ALLGENOME
        if (is.data.frame(intervals2) && nrow(intervals2) == 0) {
            return(0)
        }
        intervals2 <- rescue_ALLGENOME(intervals2, as.character(substitute(intervals2)))
    }

    # Calculate total bp in intervals2
    total_bp <- gintervals.covered_bp(intervals2)

    if (total_bp == 0) {
        return(0)
    }

    # Calculate intersection
    intersection <- gintervals.intersect(intervals1, intervals2)

    if (is.null(intersection) || nrow(intersection) == 0) {
        return(0)
    }

    # Calculate covered bp in intersection
    covered_bp <- gintervals.covered_bp(intersection)

    # Return fraction
    covered_bp / total_bp
}
