#' Finds neighbors between two sets of intervals
#'
#' For each interval in 'intervals1', finds the closest intervals from 'intervals2'.
#' Distance directionality can be determined by either the strand of the target
#' intervals (intervals2, default) or the query intervals (intervals1). When no strand column
#' is present, all intervals are treated as positive strand (strand = 1).
#'
#' This function finds for each interval in 'intervals1' the closest
#' 'maxneighbors' intervals from 'intervals2'.
#'
#' For 1D intervals the distance must fall in the range of ['mindist',
#' 'maxdist'].
#'
#' Distance is defined as the number of base pairs between the the last base pair of the query interval
#' and the first base pair of the target interval.
#'
#' **Strand handling:** By default, distance directionality is determined by the
#' 'strand' column in 'intervals2' (if present). If 'use_intervals1_strand' is TRUE,
#' distance directionality is instead determined by the 'strand' column in 'intervals1'.
#' This is particularly useful for TSS analysis where you want upstream/downstream
#' distances relative to gene direction.
#'
#' **Distance calculation modes:**
#' \itemize{
#'   \item **use_intervals1_strand = FALSE (default):** Uses intervals2 strand for directionality
#'   \item **use_intervals1_strand = TRUE:** Uses intervals1 strand for directionality
#' }
#'
#' **Important:** When 'use_intervals1_strand = TRUE', distance signs are interpreted as:
#' \itemize{
#'   \item **+ strand queries:** Negative distances = upstream, Positive distances = downstream
#'   \item **- strand queries:** Negative distances = downstream, Positive distances = upstream
#' }
#'
#' For 2D intervals two distances are calculated and returned for each axis.
#' The distances must fall in the range of ['mindist1', 'maxdist1'] for axis 1
#' and ['mindist2', 'maxdist2'] for axis 2. For selecting the closest
#' 'maxneighbors' intervals Manhattan distance is used (i.e. dist1+dist2).
#'
#' **Note:** 'use_intervals1_strand' is not yet supported for 2D intervals.
#'
#' The names of the returned columns are made unique using
#' \code{make.unique(colnames(df), sep = "")}, assuming 'df' is the result.
#'
#' If 'intervals.set.out' is not 'NULL' the result is saved as an intervals
#' set. Use this parameter if the result size exceeds the limits of the
#' physical memory.
#'
#' @param intervals1,intervals2 intervals
#' @param maxneighbors maximal number of neighbors
#' @param mindist,maxdist distance range for 1D intervals
#' @param mindist1,maxdist1,mindist2,maxdist2 distance range for 2D intervals
#' @param na.if.notfound if 'TRUE' return 'NA' interval if no matching
#' neighbors were found, otherwise omit the interval in the answer
#' @param use_intervals1_strand if 'TRUE' use intervals1 strand column for
#' distance directionality instead of intervals2 strand. If intervals1 has no
#' strand column, all intervals are treated as positive strand (strand = 1).
#' Invalid strand values (not -1 or 1) will cause an error.
#' @param warn.ignored.strand if 'TRUE' (default) show warning when 'intervals1'
#' contains a strand column that will be ignored for distance calculation
#' @param intervals.set.out intervals set name where the function result is
#' optionally outputted
#' @return If 'intervals.set.out' is 'NULL' a data frame containing the pairs
#' of intervals from 'intervals1', intervals from 'intervals2' and an
#' additional column named 'dist' ('dist1' and 'dist2' for 2D intervals)
#' representing the distance between the corresponding intervals. The intervals
#' from intervals2 would be changed to 'chrom1', 'start1', and 'end1' and for
#' 2D intervals chrom11, start11, end11 and chrom22, start22, end22. If
#' 'na.if.notfound' is 'TRUE', the data frame contains all the intervals from
#' 'intervals1' including those for which no matching neighbor was found. For
#' the latter intervals an 'NA' neighboring interval is stated. If
#' 'na.if.notfound' is 'FALSE', the data frame contains only intervals from
#' 'intervals1' for which matching neighbor(s) was found.
#' @seealso \code{\link{gintervals}}, \code{\link{gintervals.neighbors.upstream}},
#' \code{\link{gintervals.neighbors.downstream}}
#' @keywords ~intervals ~annotate ~nearest ~neighbor ~neighbors ~TSS ~strand
#' @examples
#' \dontshow{
#' options(gmax.processes = 2)
#' }
#'
#' gdb.init_examples()
#'
#' # Basic intervals
#' intervs1 <- giterator.intervals("dense_track",
#'     gintervals(1, 0, 4000),
#'     iterator = 233
#' )
#' intervs2 <- giterator.intervals(
#'     "sparse_track",
#'     gintervals(1, 0, 2000)
#' )
#'
#' # Original behavior - no strand considerations
#' gintervals.neighbors(intervs1, intervs2, 10,
#'     mindist = -300,
#'     maxdist = 500
#' )
#'
#' # Add strand to intervals2 - affects distance directionality (original behavior)
#' intervs2$strand <- c(1, 1, -1, 1)
#' gintervals.neighbors(intervs1, intervs2, 10,
#'     mindist = -300,
#'     maxdist = 500
#' )
#'
#' # TSS analysis example - use intervals1 (TSS) strand for directionality
#' tss <- data.frame(
#'     chrom = c("chr1", "chr1", "chr1"),
#'     start = c(1000, 2000, 3000),
#'     end = c(1001, 2001, 3001),
#'     strand = c(1, -1, 1), # +, -, +
#'     gene = c("GeneA", "GeneB", "GeneC")
#' )
#'
#' features <- data.frame(
#'     chrom = "chr1",
#'     start = c(500, 800, 1200, 1800, 2200, 2800, 3200),
#'     end = c(600, 900, 1300, 1900, 2300, 2900, 3300),
#'     feature_id = paste0("F", 1:7)
#' )
#'
#' # Use TSS strand for distance directionality
#' result <- gintervals.neighbors(tss, features,
#'     maxneighbors = 2,
#'     mindist = -1000, maxdist = 1000,
#'     use_intervals1_strand = TRUE
#' )
#'
#' # Convenience functions for common TSS analysis
#' # Find upstream neighbors (negative distances for + strand genes)
#' upstream <- gintervals.neighbors.upstream(tss, features,
#'     maxneighbors = 2, maxdist = 1000
#' )
#'
#' # Find downstream neighbors (positive distances for + strand genes)
#' downstream <- gintervals.neighbors.downstream(tss, features,
#'     maxneighbors = 2, maxdist = 1000
#' )
#'
#' # Find both directions
#' both_directions <- gintervals.neighbors.directional(tss, features,
#'     maxneighbors_upstream = 1,
#'     maxneighbors_downstream = 1,
#'     maxdist = 1000
#' )
#'
#' @export gintervals.neighbors
gintervals.neighbors <- function(intervals1 = NULL, intervals2 = NULL, maxneighbors = 1, mindist = -1e+09, maxdist = 1e+09,
                                 mindist1 = -1e+09, maxdist1 = 1e+09, mindist2 = -1e+09, maxdist2 = 1e+09,
                                 na.if.notfound = FALSE, use_intervals1_strand = FALSE, warn.ignored.strand = TRUE, intervals.set.out = NULL) {
    if (is.null(intervals1) || is.null(intervals2)) {
        stop(paste("Usage: gintervals.neighbors(intervals1, intervals2, maxneighbors = 1, mindist = -1e+09, maxdist = 1e+09, ",
            "mindist1 = -1e+09, maxdist1 = 1e+09, mindist2 = -1e+09, maxdist2 = 1e+09, na.if.notfound = FALSE, ",
            "use_intervals1_strand = FALSE, warn.ignored.strand = TRUE, intervals.set.out = NULL)",
            sep = ""
        ), call. = FALSE)
    }

    if (is.null(colnames)) {
        intervals1name <- deparse(substitute(intervals1), width.cutoff = 500)[1]
        intervals2name <- deparse(substitute(intervals2), width.cutoff = 500)[1]
    }

    intervals.set.out <- do.call(.gexpr2str, list(substitute(intervals.set.out)), envir = parent.frame())

    # Optional strand warning
    if (warn.ignored.strand && !use_intervals1_strand &&
        is.data.frame(intervals1) && "strand" %in% colnames(intervals1)) {
        warning("intervals1 contains a 'strand' column that will be ignored for distance calculation. ",
            "Set use_intervals1_strand = TRUE to use intervals1 strand for directionality, ",
            "or warn.ignored.strand = FALSE to suppress this warning.",
            call. = FALSE
        )
    }

    if (.gintervals.is_bigset(intervals1) || .gintervals.is_bigset(intervals2) || !is.null(intervals.set.out)) {
        res <- NULL

        FUN <- function(intervals, intervals.set.out, envir) {
            intervals1 <- intervals[[1]]
            intervals2 <- intervals[[2]]
            chrom_res <- .gcall("gfind_neighbors", intervals1, intervals2, maxneighbors, mindist, maxdist, mindist1, maxdist1, mindist2, maxdist2, na.if.notfound, FALSE, use_intervals1_strand, .misha_env())
            if (!is.null(chrom_res) && is.null(intervals.set.out)) {
                assign("res", c(get("res", envir = envir), list(chrom_res)), envir = envir)
                .gverify_max_data_size(sum(unlist(lapply(get("res", envir), nrow))), arguments = "intervals.set.out")
            }
            chrom_res
        }

        if (na.if.notfound) {
            .gintervals.apply(gintervals.chrom_sizes(intervals1), list(intervals1, intervals2), intervals.set.out, FUN, intervals.set.out, environment())
        } else {
            chroms1 <- gintervals.chrom_sizes(intervals1)
            chroms1$size <- NULL
            chroms2 <- gintervals.chrom_sizes(intervals2)
            chroms2$size <- NULL
            .gintervals.apply(merge(chroms1, chroms2), list(intervals1, intervals2), intervals.set.out, FUN, intervals.set.out, environment())
        }

        if (!is.null(res)) {
            res <- do.call(.grbind, res)
        } # much faster than calling rbind incrementally in FUN

        if (is.null(intervals.set.out)) {
            if (!is.null(res) && nrow(res)) {
                repair_names(res)
            } else {
                NULL
            }
        } else {
            retv <- 0
        } # suppress return value
    } else {
        intervals1 <- .gintervals.load_ext(intervals1)
        intervals2 <- .gintervals.load_ext(intervals2)
        res <- .gcall("gfind_neighbors", intervals1, intervals2, maxneighbors, mindist, maxdist, mindist1, maxdist1, mindist2, maxdist2, na.if.notfound, TRUE, use_intervals1_strand, .misha_env())
        repair_names(res)
    }
}

repair_names <- function(dataframe) {
    # if there are duplicate names - add a numeric suffix
    if (length(unique(names(dataframe))) < length(names(dataframe))) {
        names(dataframe) <- make.unique(names(dataframe), sep = "")
    }
    return(dataframe)
}


#' Directional neighbor finding functions
#'
#' These functions find neighbors using query strand directionality, where
#' upstream/downstream directionality is determined by the strand of the query
#' intervals rather than the target intervals. This is particularly useful for
#' TSS analysis where you want distances relative to gene direction.
#'
#' **Distance interpretation:**
#' \itemize{
#'   \item **Positive strand queries:** upstream distances < 0, downstream distances > 0
#'   \item **Negative strand queries:** upstream distances > 0, downstream distances < 0
#' }
#' If no strand column is present, all intervals are treated as positive strand.
#'
#' @param query_intervals intervals with strand information (query intervals)
#' @param target_intervals intervals to search for neighbors
#' @param maxneighbors maximum number of neighbors per query interval (default: 1)
#' @param maxdist maximum distance to search (default: 1e+09)
#' @param maxneighbors_upstream maximum upstream neighbors per query interval (default: 1)
#' @param maxneighbors_downstream maximum downstream neighbors per query interval (default: 1)
#' @param ... additional arguments passed to \code{gintervals.neighbors}
#'
#' @return
#' \describe{
#'   \item{gintervals.neighbors.upstream}{data frame of upstream neighbors}
#'   \item{gintervals.neighbors.downstream}{data frame of downstream neighbors}
#'   \item{gintervals.neighbors.directional}{list with 'upstream' and 'downstream' components}
#' }
#'
#' @examples
#' \dontshow{
#' options(gmax.processes = 2)
#' }
#'
#' gdb.init_examples()
#'
#' # Create TSS intervals with strand information
#' tss <- data.frame(
#'     chrom = c("chr1", "chr1", "chr1"),
#'     start = c(1000, 2000, 3000),
#'     end = c(1001, 2001, 3001),
#'     strand = c(1, -1, 1), # +, -, +
#'     gene = c("GeneA", "GeneB", "GeneC")
#' )
#'
#' # Create regulatory features
#' features <- data.frame(
#'     chrom = "chr1",
#'     start = c(500, 800, 1200, 1800, 2200, 2800, 3200),
#'     end = c(600, 900, 1300, 1900, 2300, 2900, 3300),
#'     feature_id = paste0("F", 1:7)
#' )
#'
#' # Find upstream neighbors (promoter analysis)
#' upstream <- gintervals.neighbors.upstream(tss, features,
#'     maxneighbors = 2, maxdist = 1000
#' )
#' print(upstream)
#'
#' # Find downstream neighbors (gene body analysis)
#' downstream <- gintervals.neighbors.downstream(tss, features,
#'     maxneighbors = 2, maxdist = 5000
#' )
#' print(downstream)
#'
#' # Find both directions in one call
#' both <- gintervals.neighbors.directional(tss, features,
#'     maxneighbors_upstream = 1,
#'     maxneighbors_downstream = 1,
#'     maxdist = 1000
#' )
#' print(both$upstream)
#' print(both$downstream)
#'
#' @seealso \code{\link{gintervals.neighbors}}
#' @rdname directional-neighbors
#' @export
gintervals.neighbors.upstream <- function(query_intervals, target_intervals,
                                          maxneighbors = 1, maxdist = 1e+09, ...) {
    gintervals.neighbors(query_intervals, target_intervals,
        maxneighbors = maxneighbors,
        mindist = -maxdist, maxdist = 0,
        use_intervals1_strand = TRUE,
        warn.ignored.strand = FALSE, ...
    )
}

#' @rdname directional-neighbors
#' @export
gintervals.neighbors.downstream <- function(query_intervals, target_intervals,
                                            maxneighbors = 1, maxdist = 1e+09, ...) {
    gintervals.neighbors(query_intervals, target_intervals,
        maxneighbors = maxneighbors,
        mindist = 0, maxdist = maxdist,
        use_intervals1_strand = TRUE,
        warn.ignored.strand = FALSE, ...
    )
}

#' @rdname directional-neighbors
#' @export
gintervals.neighbors.directional <- function(query_intervals, target_intervals,
                                             maxneighbors_upstream = 1,
                                             maxneighbors_downstream = 1,
                                             maxdist = 1e+09, ...) {
    upstream <- gintervals.neighbors.upstream(query_intervals, target_intervals,
        maxneighbors = maxneighbors_upstream,
        maxdist = maxdist, ...
    )

    downstream <- gintervals.neighbors.downstream(query_intervals, target_intervals,
        maxneighbors = maxneighbors_downstream,
        maxdist = maxdist, ...
    )

    list(upstream = upstream, downstream = downstream)
}
