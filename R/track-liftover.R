# Track liftover functions

#' Imports a track from another assembly
#'
#' Imports a track from another assembly.
#'
#' This function imports a track located in 'src.track.dir' of another assembly
#' to the current database. Chain file instructs how the conversion of
#' coordinates should be done. It can be either a name of a chain file or a
#' data frame in the same format as returned by 'gintervals.load_chain'
#' function. The name of the newly created track is specified by 'track'
#' argument and 'description' is added as a track attribute.
#'
#' Note: When passing a pre-loaded chain (data frame), overlap policies cannot
#' be specified - they are taken from the chain's attributes that were set
#' during loading. When passing a chain file path, policies can be specified
#' and will be used for loading. Aggregation parameters (multi_target_agg,
#' params, na.rm, min_n) can always be specified regardless of chain type.
#'
#' @param track name of a created track
#' @param description a character string description
#' @param src.track.dir path to the directory of the source track
#' @param multi_target_agg aggregation/selection policy for contributors that land on the same target locus. When multiple source intervals map to overlapping regions in the target genome (after applying tgt_overlap_policy), their values must be combined into a single value.
#' @param params additional parameters for aggregation (e.g., for "nth" aggregation)
#' @param na.rm logical indicating whether NA values should be removed before aggregation (default: TRUE)
#' @param min_n minimum number of non-NA values required for aggregation. If fewer values are available, the result will be NA.
#' @return None.
#'
#' @note
#' Terminology note for UCSC chain format users: In the UCSC chain format specification,
#' the fields prefixed with 't' (tName, tStart, tEnd, etc.) are called "target" or "reference",
#' while fields prefixed with 'q' (qName, qStart, qEnd, etc.) are called "query". However,
#' misha uses reversed terminology: UCSC's "target/reference" corresponds to misha's "source"
#' (chromsrc, startsrc, endsrc), and UCSC's "query" corresponds to misha's "target"
#' (chrom, start, end).
#'
#' @seealso \code{\link{gintervals.load_chain}},
#' \code{\link{gintervals.liftover}}
#' @keywords ~track ~liftover ~chain
#' @inheritParams gintervals.liftover
#' @export gtrack.liftover
gtrack.liftover <- function(track = NULL,
                            description = NULL,
                            src.track.dir = NULL,
                            chain = NULL,
                            src_overlap_policy = "error",
                            tgt_overlap_policy = "auto",
                            multi_target_agg = c(
                                "mean", "median", "sum", "min", "max", "count",
                                "first", "last", "nth",
                                "max.coverage_len", "min.coverage_len",
                                "max.coverage_frac", "min.coverage_frac"
                            ),
                            params = NULL,
                            na.rm = TRUE,
                            min_n = NULL,
                            min_score = NULL) {
    if (is.null(substitute(track)) || is.null(description) || is.null(src.track.dir) || is.null(chain)) {
        stop("Usage: gtrack.liftover(track, description, src.track.dir, chain, src_overlap_policy = \"error\", tgt_overlap_policy = \"auto\", ...)", call. = FALSE)
    }
    .gcheckroot()

    # Validate aggregation parameters (these can always be set)
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

    trackstr <- do.call(.gexpr2str, list(substitute(track)), envir = parent.frame())

    normalize_policy <- function(policy) {
        if (is.null(policy)) {
            return(NULL)
        }
        if (policy %in% c("auto", "auto_first")) {
            return("auto")
        }
        policy
    }

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

    .gconfirmtrackcreate(trackstr)
    trackdir <- .track_dir(trackstr)
    direxisted <- file.exists(trackdir)
    success <- FALSE
    tryCatch(
        {
            .gcall(
                "gtrack_liftover",
                trackstr,
                src.track.dir,
                chain.intervs,
                src_overlap_policy,
                tgt_overlap_policy,
                multi_target_agg,
                nth_param,
                na.rm,
                if (is.null(min_n)) NA_integer_ else min_n,
                min_score,
                .misha_env(),
                silent = TRUE
            )
            .gdb.add_track(trackstr)
            if (is.character(chain)) {
                .gtrack.attr.set(trackstr, "created.by", sprintf("gtrack.liftover(%s, description, \"%s\", \"%s\")", trackstr, src.track.dir, chain), TRUE)
            } else {
                .gtrack.attr.set(trackstr, "created.by", sprintf("gtrack.liftover(%s, description, \"%s\", chain)", trackstr, src.track.dir), TRUE)
            }
            .gtrack.attr.set(trackstr, "created.date", date(), TRUE)
            .gtrack.attr.set(trackstr, "created.user", Sys.getenv("USER"), TRUE)
            .gtrack.attr.set(trackstr, "description", description, TRUE)
            success <- TRUE

            # If database is indexed, automatically convert the track to indexed format
            # For empty tracks (no chromosome files), create an empty indexed track
            if (.gdb.is_indexed()) {
                track_has_files <- length(list.files(trackdir, pattern = "^[^.]")) > 0
                if (track_has_files) {
                    gtrack.convert_to_indexed(trackstr)
                } else {
                    # Create empty indexed track for empty tracks
                    # This ensures gtrack.info and gextract work correctly
                    .gcall("gtrack_create_empty_indexed", trackstr, .misha_env())
                }
            }
        },
        finally = {
            if (!success && !direxisted) {
                unlink(trackdir, recursive = TRUE)
                .gdb.rm_track(trackstr)
            }
        }
    )
    retv <- 0 # suppress return value
}
