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
