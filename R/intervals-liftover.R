# Liftover and chain functions

#' Converts intervals from another assembly
#'
#' Converts intervals from another assembly to the current one.
#'
#' This function converts 'intervals' from another assembly to the current one.
#' Chain file instructs how the conversion of coordinates should be done. It
#' can be either a name of a chain file or a data frame in the same format as
#' returned by 'gintervals.load_chain' function.
#'
#' The converted intervals are returned. An additional column named
#' 'intervalID' is added to the resulted data frame. For each interval in the
#' resulted intervals it indicates the index of the original interval.
#'
#' Note: When passing a pre-loaded chain (data frame), overlap policies cannot
#' be specified - they are taken from the chain's attributes that were set
#' during loading. When passing a chain file path, policies can be specified
#' and will be used for loading.
#'
#' @param intervals intervals from another assembly
#' @param chain name of chain file or data frame as returned by
#' 'gintervals.load_chain'
#' @param src_overlap_policy policy for handling source overlaps: "error" (default), "keep", or "discard". "keep" allows one source interval to map to multiple target intervals, "discard" discards all source intervals that have overlaps and "error" throws an error if source overlaps are detected.
#' @param tgt_overlap_policy policy for handling target overlaps. One of:
#' \tabular{ll}{
#'   Policy \tab Description \cr
#'   error \tab Throws an error if any target overlaps are detected. \cr
#'   auto \tab Default. Alias for "auto_score". \cr
#'   auto_score \tab Resolves overlaps by segmenting the target region and selecting the best chain for each segment based on alignment score (highest score wins). Tie-breakers: longest span, then lowest chain_id. \cr
#'   auto_longer \tab Resolves overlaps by segmenting and selecting the chain with the longest span for each segment. Tie-breakers: highest score, then lowest chain_id. \cr
#'   auto_first \tab Resolves overlaps by segmenting and selecting the chain with the lowest chain_id for each segment. \cr
#'   keep \tab Preserves all overlapping intervals. \cr
#'   discard \tab Discards any chain interval that has a target overlap with another chain interval. \cr
#'   agg \tab Segments overlaps into smaller disjoint regions where each region contains all contributing chains, allowing downstream aggregation to process multiple values per region. \cr
#'   best_source_cluster \tab Best source cluster strategy based on source overlap. When multiple chains map a source interval, clusters them by source overlap: if chain source intervals overlap (indicating true duplications), all mappings are retained; if chain source intervals are disjoint (indicating conflicting/alternative mappings), only the cluster with the largest total target length is kept. \cr
#' }
#' @param min_score optional minimum alignment score threshold. Chains with scores below this value are filtered out. Useful for excluding low-quality alignments.
#' @param include_metadata logical; if TRUE, adds 'score' column to the output indicating the alignment score of the chain used for each mapping. Only applicable with "auto_score" or "auto" policy.
#' @param canonic logical; if TRUE, merges adjacent target intervals that originated from the same source interval (same intervalID) and same chain (same chain_id). This is useful when a source interval maps to multiple adjacent target blocks due to chain gaps.
#' @param value_col optional character string specifying the name of a numeric column in the intervals data frame to track through the liftover. When specified, this column's values are preserved in the output with the same column name. Use with multi_target_agg to aggregate values when multiple source intervals map to overlapping target regions.
#' @param multi_target_agg aggregation method to use when value_col is specified. One of: "mean", "median", "sum", "min", "max", "count", "first", "last", "nth", "max.coverage_len", "min.coverage_len", "max.coverage_frac", "min.coverage_frac". Default: "mean". Ignored when value_col is NULL.
#' @param params additional parameters for specific aggregation methods. Currently only used for "nth" aggregation, where it specifies which element to select (e.g., params = 2 for second element, or params = list(n = 2)).
#' @param na.rm logical; if TRUE (default), NA values are removed before aggregation. If FALSE, any NA in the values will cause the result to be NA. Only used when value_col is specified.
#' @param min_n optional minimum number of non-NA observations required for aggregation. If fewer observations are available, the result is NA. NULL (default) means no minimum. Only used when value_col is specified.
#' @return A data frame representing the converted intervals. For 1D intervals, always includes 'intervalID' (index of original interval) and 'chain_id' (identifier of the chain that produced the mapping) columns. The chain_id column is essential for distinguishing results when a source interval maps to multiple target regions via different chains (duplications). When include_metadata=TRUE, also includes 'score' column. When value_col is specified, includes the value column with its original name.
#' @seealso \code{\link{gintervals.load_chain}}, \code{\link{gtrack.liftover}},
#' \code{\link{gintervals}}
#' @keywords ~intervals ~liftover ~chain
#' @examples
#' \dontshow{
#' options(gmax.processes = 2)
#' }
#'
#' gdb.init_examples()
#' chainfile <- paste(.misha$GROOT, "data/test.chain", sep = "/")
#' intervs <- data.frame(
#'     chrom = "chr25", start = c(0, 7000),
#'     end = c(6000, 20000)
#' )
#' # Liftover with default policies
#' gintervals.liftover(intervs, chainfile)
#'
#' # Liftover keeping source overlaps (one source interval may map to multiple targets)
#' # gintervals.liftover(intervs, chainfile, src_overlap_policy = "keep")
#'
#' @export gintervals.liftover
gintervals.liftover <- function(intervals = NULL,
                                chain = NULL,
                                src_overlap_policy = "error",
                                tgt_overlap_policy = "auto",
                                min_score = NULL,
                                include_metadata = FALSE,
                                canonic = FALSE,
                                value_col = NULL,
                                multi_target_agg = c(
                                    "mean", "median", "sum", "min", "max", "count",
                                    "first", "last", "nth",
                                    "max.coverage_len", "min.coverage_len",
                                    "max.coverage_frac", "min.coverage_frac"
                                ),
                                params = NULL,
                                na.rm = TRUE,
                                min_n = NULL) {
    if (is.null(intervals) || is.null(chain)) {
        stop("Usage: gintervals.liftover(intervals, chain, src_overlap_policy = \"error\", tgt_overlap_policy = \"auto\", min_score = NULL, include_metadata = FALSE, canonic = FALSE, value_col = NULL, multi_target_agg = \"mean\", params = NULL, na.rm = TRUE, min_n = NULL)", call. = FALSE)
    }
    .gcheckroot()

    if (!is.logical(include_metadata) || length(include_metadata) != 1) {
        stop("include_metadata must be a single logical value", call. = FALSE)
    }

    intervals <- rescue_ALLGENOME(intervals, as.character(substitute(intervals)))

    if (is.character(chain)) {
        # Chain file path provided - validate and use the policies
        if (!src_overlap_policy %in% c("error", "keep", "discard")) {
            stop("src_overlap_policy must be 'error', 'keep', or 'discard'", call. = FALSE)
        }

        if (!tgt_overlap_policy %in% c("error", "auto", "auto_first", "auto_longer", "auto_score", "discard", "keep", "agg", "best_source_cluster", "best_cluster_union", "best_cluster_sum", "best_cluster_max")) {
            stop("tgt_overlap_policy must be 'error', 'auto', 'auto_first', 'auto_longer', 'auto_score', 'keep', 'discard', 'agg', 'best_source_cluster', 'best_cluster_union', 'best_cluster_sum', or 'best_cluster_max'", call. = FALSE)
        }

        if (!is.null(min_score) && (!is.numeric(min_score) || length(min_score) != 1)) {
            stop("min_score must be a single numeric value", call. = FALSE)
        }

        # Convert "auto" to "auto_score" alias
        if (tgt_overlap_policy == "auto") {
            tgt_overlap_policy <- "auto_score"
        }

        chain.intervs <- gintervals.load_chain(chain, src_overlap_policy, tgt_overlap_policy, min_score = min_score)
    } else {
        # Pre-loaded chain provided
        chain.intervs <- chain

        # Check if chain has policy attributes
        chain_src_policy <- attr(chain.intervs, "src_overlap_policy")
        chain_tgt_policy <- attr(chain.intervs, "tgt_overlap_policy")

        if (!is.null(chain_src_policy) && !is.null(chain_tgt_policy)) {
            # Chain has attributes - use them, error if user tries to override
            policies_set <- !missing(src_overlap_policy) || !missing(tgt_overlap_policy) || !missing(min_score)
            if (policies_set) {
                stop("When using a pre-loaded chain, overlap policies cannot be specified. Set policies when loading the chain with gintervals.load_chain().", call. = FALSE)
            }
            src_overlap_policy <- chain_src_policy
            tgt_overlap_policy <- chain_tgt_policy
        } else {
            # Chain doesn't have attributes (e.g., manually created) - validate and use user-specified policies
            if (!src_overlap_policy %in% c("error", "keep", "discard")) {
                stop("src_overlap_policy must be 'error', 'keep', or 'discard'", call. = FALSE)
            }

            if (!tgt_overlap_policy %in% c("error", "auto", "auto_first", "auto_longer", "auto_score", "discard", "keep", "agg", "best_source_cluster", "best_cluster_union", "best_cluster_sum", "best_cluster_max")) {
                stop("tgt_overlap_policy must be 'error', 'auto', 'auto_first', 'auto_longer', 'auto_score', 'keep', 'discard', 'agg', 'best_source_cluster', 'best_cluster_union', 'best_cluster_sum', or 'best_cluster_max'", call. = FALSE)
            }

            # Convert "auto" to "auto_score" alias
            if (tgt_overlap_policy == "auto") {
                tgt_overlap_policy <- "auto_score"
            }
        }
    }

    if (!is.logical(canonic) || length(canonic) != 1) {
        stop("canonic must be a single logical value", call. = FALSE)
    }

    # Validate value_col and aggregation parameters
    use_aggregation <- !is.null(value_col)

    if (use_aggregation) {
        # Validate value_col exists in intervals
        if (!is.character(value_col) || length(value_col) != 1) {
            stop("value_col must be a single character string specifying the column name", call. = FALSE)
        }

        if (!value_col %in% names(intervals)) {
            stop(sprintf("value_col '%s' not found in intervals", value_col), call. = FALSE)
        }

        # Validate aggregation parameters
        multi_target_agg <- match.arg(multi_target_agg)

        if (!is.logical(na.rm) || length(na.rm) != 1 || is.na(na.rm)) {
            stop("na.rm must be a single non-NA logical value", call. = FALSE)
        }

        if (!is.null(min_n)) {
            if (!is.numeric(min_n) || length(min_n) != 1 ||
                is.na(min_n) || min_n < 0 || min_n != as.integer(min_n)) {
                stop("min_n must be NULL or a non-negative integer", call. = FALSE)
            }
            min_n <- as.integer(min_n)
        }

        nth_param <- NA_integer_
        if (identical(multi_target_agg, "nth")) {
            if (is.null(params)) {
                stop("params must be supplied for 'nth' aggregation (e.g. params = 2 or params = list(n = 2))", call. = FALSE)
            }

            extract_n <- function(obj) {
                if (is.list(obj)) {
                    if (length(obj) == 0L) {
                        stop("params list must contain an element 'n' for 'nth'", call. = FALSE)
                    }
                    if (!is.null(names(obj)) && "n" %in% names(obj)) {
                        return(obj[["n"]])
                    }
                    if (length(obj) == 1L) {
                        return(obj[[1]])
                    }
                    stop("params must contain a single numeric value (or named 'n') for 'nth'", call. = FALSE)
                }
                obj
            }

            n_value <- extract_n(params)
            if (length(n_value) != 1L || is.na(n_value) || !is.numeric(n_value)) {
                stop("params for 'nth' must be a single positive integer", call. = FALSE)
            }
            nth_param <- as.integer(n_value)
            if (nth_param <= 0L) {
                stop("params for 'nth' must be a positive integer", call. = FALSE)
            }
        } else if (!is.null(params)) {
            stop(sprintf("params is only supported for 'nth' aggregation, not '%s'", multi_target_agg), call. = FALSE)
        }
    } else {
        # No aggregation - use defaults
        multi_target_agg <- "mean"
        nth_param <- NA_integer_
        min_n <- -1L # -1 means disabled
    }

    # Convert NA/NULL min_n to -1 (disabled) for C++
    if (is.null(min_n) || (length(min_n) == 1 && is.na(min_n))) {
        min_n <- -1L
    }

    .gcall("gintervs_liftover", intervals, chain.intervs, src_overlap_policy, tgt_overlap_policy, canonic, include_metadata, value_col, multi_target_agg, nth_param, na.rm, min_n, .misha_env())
}


#' Loads a named intervals set
#'
#' Loads a named intervals set.
#'
#' This function loads and returns intervals stored in a named intervals set.
#'
#' If intervals set contains 1D intervals and 'chrom' is not 'NULL' only the
#' intervals of the given chromosome are returned.
#'
#' Likewise if intervals set contains 2D intervals and 'chrom1', 'chrom2' are
#' not 'NULL' only the intervals of the given pair of chromosomes are returned.
#'
#' For big intervals sets 'chrom' parameter (1D case) / 'chrom1', 'chrom2'
#' parameters (2D case) must be specified. In other words: big intervals sets
#' can be loaded only by chromosome or chromosome pair.
#'
#' @param intervals.set name of an intervals set
#' @param chrom chromosome for 1D intervals set
#' @param chrom1 first chromosome for 2D intervals set
#' @param chrom2 second chromosome for 2D intervals set
#' @return A data frame representing the intervals.
#' @seealso \code{\link{gintervals.save}}, \code{\link{gintervals.is.bigset}},
#' \code{\link{gintervals.exists}}, \code{\link{gintervals.ls}},
#' \code{\link{gintervals}}, \code{\link{gintervals.2d}}
#' @keywords ~intervals
#' @examples
#' \dontshow{
#' options(gmax.processes = 2)
#' }
#'
#' gdb.init_examples()
#' gintervals.load("annotations")
#'
#' @export gintervals.load
gintervals.load <- function(intervals.set = NULL, chrom = NULL, chrom1 = NULL, chrom2 = NULL) {
    .gintervals.load_ext(intervals.set, chrom, chrom1, chrom2, TRUE)
}

#' Validate source chromosomes against source genome
#'
#' @param chain Chain data frame with source chromosome information
#' @param src_groot Path to source genome database
#' @return NULL (stops with error if validation fails)
#' @keywords internal
#' @noRd
.validate_source_chromosomes <- function(chain, src_groot) {
    # Save current genome state
    old_groot <- .misha$GROOT
    old_allgenome <- .misha$ALLGENOME
    old_chrom_alias <- if (exists("CHROM_ALIAS", envir = .misha)) .misha$CHROM_ALIAS else NULL

    # Ensure genome is restored even if error occurs
    on.exit(
        {
            .misha$GROOT <- old_groot
            .misha$ALLGENOME <- old_allgenome
            if (!is.null(old_chrom_alias)) {
                .misha$CHROM_ALIAS <- old_chrom_alias
            } else if (exists("CHROM_ALIAS", envir = .misha)) {
                rm("CHROM_ALIAS", envir = .misha)
            }
        },
        add = TRUE
    )

    # Temporarily switch to source genome for validation
    gdb.init(src_groot)

    # Create source intervals data frame
    src_intervals <- data.frame(
        chrom = chain$chromsrc,
        start = chain$startsrc,
        end = chain$endsrc,
        stringsAsFactors = FALSE
    )

    validated <- .gintervals(chain$chromsrc, chain$startsrc, chain$endsrc, chain$strandsrc)

    # Genome will be restored by on.exit
}


#' Loads assembly conversion table from a chain file
#'
#' Loads assembly conversion table from a chain file.
#'
#' This function reads a file in 'chain' format and returns assembly conversion
#' table that can be used in 'gtrack.liftover' and 'gintervals.liftover'.
#'
#' Source overlaps occur when the same source genome position maps to multiple
#' target genome positions. Target overlaps occur when multiple source positions
#' map to overlapping regions in the target genome.
#'
#' The 'src_overlap_policy' controls how source overlaps are handled:
#' \itemize{
#'   \item "error" (default): Throw an error if source overlaps are detected
#'   \item "keep": Keep all mappings, allowing one source to map to multiple targets
#'   \item "discard": Remove all chain intervals involved in source overlaps
#' }
#'
#' The 'tgt_overlap_policy' controls how target overlaps are handled:
#' \itemize{
#'   \item "error": Throw an error if target overlaps are detected
#'   \item "auto" (default) / "auto_first": Keep the first overlapping chain (original file order) by trimming or discarding later overlaps while keeping source/target lengths consistent
#'   \item "auto_longer": Keep the longer overlapping chain and trim/drop the shorter ones
#'   \item "discard": Remove all chain intervals involved in target overlaps
#'   \item "keep": Allow target overlaps to remain untouched (liftover must be able to handle them)
#' }
#'
#' @param file name of chain file
#' @param src_groot optional path to source genome database for validating source chromosomes and coordinates. If provided, the function temporarily switches to this database to verify that all source chromosomes exist and coordinates are within bounds, then restores the original database.
#' @return A data frame representing assembly conversion table with columns: chrom, start, end, strand, chromsrc, startsrc, endsrc, strandsrc, chain_id, score.
#' @seealso \code{\link{gintervals.liftover}}, \code{\link{gtrack.liftover}}
#' @keywords ~intervals ~liftover ~chain
#' @examples
#' \dontshow{
#' options(gmax.processes = 2)
#' }
#'
#' gdb.init_examples()
#' chainfile <- paste(.misha$GROOT, "data/test.chain", sep = "/")
#' # Load chain file with default policies
#' gintervals.load_chain(chainfile)
#'
#' @inheritParams gintervals.liftover
#' @export gintervals.load_chain
gintervals.load_chain <- function(file = NULL, src_overlap_policy = "error", tgt_overlap_policy = "auto", src_groot = NULL, min_score = NULL) {
    if (is.null(file)) {
        stop("Usage: gintervals.load_chain(file, src_overlap_policy = \"error\", tgt_overlap_policy = \"auto\", src_groot = NULL, min_score = NULL)", call. = FALSE)
    }
    .gcheckroot()

    if (!src_overlap_policy %in% c("error", "keep", "discard")) {
        stop("src_overlap_policy must be 'error', 'keep', or 'discard'", call. = FALSE)
    }

    if (!tgt_overlap_policy %in% c("error", "auto", "auto_first", "auto_longer", "auto_score", "discard", "keep", "agg", "best_source_cluster", "best_cluster_union", "best_cluster_sum", "best_cluster_max")) {
        stop("tgt_overlap_policy must be 'error', 'auto', 'auto_first', 'auto_longer', 'auto_score', 'keep', 'discard', 'agg', 'best_source_cluster', 'best_cluster_union', 'best_cluster_sum', or 'best_cluster_max'", call. = FALSE)
    }

    if (!is.null(min_score) && (!is.numeric(min_score) || length(min_score) != 1)) {
        stop("min_score must be a single numeric value", call. = FALSE)
    }

    # Convert "auto" to "auto_score" alias
    if (tgt_overlap_policy == "auto") {
        tgt_overlap_policy <- "auto_score"
    }

    chain <- .gcall("gchain2interv", file, src_overlap_policy, tgt_overlap_policy, min_score, .misha_env())

    # Handle case where all chains were discarded
    if (is.null(chain) || nrow(chain) == 0) {
        chain <- data.frame(
            chrom = character(0),
            start = numeric(0),
            end = numeric(0),
            strand = numeric(0),
            chromsrc = character(0),
            startsrc = numeric(0),
            endsrc = numeric(0),
            strandsrc = numeric(0),
            chain_id = integer(0),
            score = numeric(0),
            stringsAsFactors = FALSE
        )
    }

    if (!is.null(src_groot)) {
        .validate_source_chromosomes(chain, src_groot)
    }

    # Store policies as attributes so they can be used during liftover
    attr(chain, "src_overlap_policy") <- src_overlap_policy
    attr(chain, "tgt_overlap_policy") <- tgt_overlap_policy
    if (!is.null(min_score)) {
        attr(chain, "min_score") <- min_score
    }

    chain
}

#' Transforms existing intervals to a chain format
#'
#' Transforms existing intervals to a chain format by validating required columns
#' and adding chain attributes.
#'
#' This function checks that the input intervals data frame has all the required
#' columns for a chain format and adds the necessary attributes. A chain format
#' requires both target coordinates (chrom, start, end, strand) and source
#' coordinates (chromsrc, startsrc, endsrc, strandsrc), as well as chain_id and
#' score columns.
#'
#' @param intervals a data frame with chain columns: chrom, start, end, strand,
#'   chromsrc, startsrc, endsrc, strandsrc, chain_id, score
#' @param src_overlap_policy source overlap policy: "error", "keep", or "discard"
#' @param tgt_overlap_policy target overlap policy: "error", "auto", "auto_first",
#'   "auto_longer", "auto_score", "discard", "keep", or "agg"
#' @param min_score optional minimum alignment score threshold
#' @return A data frame in chain format with chain attributes set
#' @seealso \code{\link{gintervals.load_chain}}, \code{\link{gintervals.liftover}}
#' @keywords ~intervals ~liftover ~chain
#' @examples
#' \dontshow{
#' options(gmax.processes = 2)
#' }
#'
#' gdb.init_examples()
#'
#' # Create a chain from existing intervals
#' chain_data <- data.frame(
#'     chrom = "chr1",
#'     start = 1000,
#'     end = 2000,
#'     strand = 0,
#'     chromsrc = "chr1",
#'     startsrc = 5000,
#'     endsrc = 6000,
#'     strandsrc = 0,
#'     chain_id = 1L,
#'     score = 1000.0
#' )
#' chain <- gintervals.as_chain(chain_data)
#'
#' @export gintervals.as_chain
gintervals.as_chain <- function(intervals = NULL, src_overlap_policy = "error", tgt_overlap_policy = "auto", min_score = NULL) {
    if (!is.data.frame(intervals)) {
        if (is.null(intervals)) {
            stop("Usage: gintervals.as_chain(intervals, src_overlap_policy = \"error\", tgt_overlap_policy = \"auto\", min_score = NULL)", call. = FALSE)
        }
        stop("intervals must be a data frame", call. = FALSE)
    }

    # Required columns for chain format
    required_cols <- c("chrom", "start", "end", "strand", "chromsrc", "startsrc", "endsrc", "strandsrc", "chain_id", "score")
    missing_cols <- setdiff(required_cols, colnames(intervals))
    if (length(missing_cols) > 0) {
        stop(sprintf("Missing required columns: %s", paste(missing_cols, collapse = ", ")), call. = FALSE)
    }

    # Validate policies
    if (!src_overlap_policy %in% c("error", "keep", "discard")) {
        stop("src_overlap_policy must be 'error', 'keep', or 'discard'", call. = FALSE)
    }

    if (!tgt_overlap_policy %in% c("error", "auto", "auto_first", "auto_longer", "auto_score", "discard", "keep", "agg")) {
        stop("tgt_overlap_policy must be 'error', 'auto', 'auto_first', 'auto_longer', 'auto_score', 'keep', 'discard', or 'agg'", call. = FALSE)
    }

    if (!is.null(min_score) && (!is.numeric(min_score) || length(min_score) != 1)) {
        stop("min_score must be a single numeric value", call. = FALSE)
    }

    # Convert "auto" to "auto_score" alias
    if (tgt_overlap_policy == "auto") {
        tgt_overlap_policy <- "auto_score"
    }

    # Validate data types
    if (!is.numeric(intervals$start) || !is.numeric(intervals$end)) {
        stop("start and end columns must be numeric", call. = FALSE)
    }

    if (!is.numeric(intervals$startsrc) || !is.numeric(intervals$endsrc)) {
        stop("startsrc and endsrc columns must be numeric", call. = FALSE)
    }

    if (!is.numeric(intervals$strand)) {
        stop("strand column must be numeric", call. = FALSE)
    }

    if (!is.numeric(intervals$strandsrc)) {
        stop("strandsrc column must be numeric", call. = FALSE)
    }

    if (!is.integer(intervals$chain_id) && !is.numeric(intervals$chain_id)) {
        stop("chain_id column must be integer or numeric", call. = FALSE)
    }

    if (!is.numeric(intervals$score)) {
        stop("score column must be numeric", call. = FALSE)
    }

    # Validate strand values
    invalid_strand <- intervals$strand != as.integer(intervals$strand) | intervals$strand < -1 | intervals$strand > 1
    if (any(invalid_strand, na.rm = TRUE)) {
        stop("strand values must be -1, 0, or 1", call. = FALSE)
    }

    invalid_strandsrc <- intervals$strandsrc != as.integer(intervals$strandsrc) | intervals$strandsrc < -1 | intervals$strandsrc > 1
    if (any(invalid_strandsrc, na.rm = TRUE)) {
        stop("strandsrc values must be -1, 0, or 1", call. = FALSE)
    }

    # Ensure chain_id is integer
    if (!is.integer(intervals$chain_id)) {
        intervals$chain_id <- as.integer(intervals$chain_id)
    }

    # Store policies as attributes so they can be used during liftover
    attr(intervals, "src_overlap_policy") <- src_overlap_policy
    attr(intervals, "tgt_overlap_policy") <- tgt_overlap_policy
    if (!is.null(min_score)) {
        attr(intervals, "min_score") <- min_score
    }

    intervals
}

.validate_source_chromosomes <- function(chain, src_groot) {
    # Save current genome state
    old_groot <- .misha$GROOT
    old_allgenome <- .misha$ALLGENOME
    old_chrom_alias <- if (exists("CHROM_ALIAS", envir = .misha)) .misha$CHROM_ALIAS else NULL

    # Ensure genome is restored even if error occurs
    on.exit(
        {
            .misha$GROOT <- old_groot
            .misha$ALLGENOME <- old_allgenome
            if (!is.null(old_chrom_alias)) {
                .misha$CHROM_ALIAS <- old_chrom_alias
            } else if (exists("CHROM_ALIAS", envir = .misha)) {
                rm("CHROM_ALIAS", envir = .misha)
            }
        },
        add = TRUE
    )

    # Temporarily switch to source genome for validation
    gdb.init(src_groot)

    # Create source intervals data frame
    src_intervals <- data.frame(
        chrom = chain$chromsrc,
        start = chain$startsrc,
        end = chain$endsrc,
        stringsAsFactors = FALSE
    )

    validated <- .gintervals(chain$chromsrc, chain$startsrc, chain$endsrc, chain$strandsrc)

    # Genome will be restored by on.exit
}
