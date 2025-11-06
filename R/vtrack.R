.gvtrack <- function(vtrack) {
    vtrackstr <- do.call(.gexpr2str, list(substitute(vtrack)), envir = parent.frame())
    if (!is.character(vtrackstr) || length(vtrackstr) != 1) {
        stop(sprintf("Virtual track must be specified as a character string"), call. = FALSE)
    }

    if (is.na(match(vtrackstr, gvtrack.ls()))) {
        stop(sprintf("Virtual track %s does not exist", vtrackstr), call. = FALSE)
    }
    vtrackstr
}

.gvtrack.get <- function(vtrackstr) {
    if (is.na(match(vtrackstr, gvtrack.ls()))) {
        stop(sprintf("Virtual track %s does not exist", vtrackstr), call. = FALSE)
    }

    gwd <- get("GWD", envir = .misha)
    get("GVTRACKS", envir = .misha)[[gwd]][[vtrackstr]]
}

.gvtrack.set <- function(vtrackstr, var) {
    if (exists("GVTRACKS", envir = .misha)) {
        gvtracks <- get("GVTRACKS", envir = .misha)
    } else {
        gvtracks <- list()
    }

    gwds <- names(gvtracks)
    if (!is.list(gvtracks) || (length(gvtracks) && !is.character(gwds)) || length(gvtracks) != length(gwds)) {
        stop("Invalid format of GVTRACKS variable.\nTo continue working with virtual tracks please remove this variable from the environment.", call. = FALSE)
    }

    gwd <- get("GWD", envir = .misha)
    idx1 <- match(gwd, gwds)
    if (is.na(idx1)) {
        gwds <- c(gwds, gwd)
        idx1 <- length(gwds)
        gvtracks[[idx1]] <- list()
    }

    vtracks <- gvtracks[[idx1]]
    names(gvtracks) <- gwds

    vtracknames <- names(vtracks)
    if (!is.list(vtracks) || (length(vtracks) && !is.character(vtracknames)) || length(vtracks) != length(vtracknames)) {
        stop("Invalid format of GVTRACKS variable.\nTo continue working with virtual tracks please remove this variable from the environment.", call. = FALSE)
    }

    idx2 <- match(vtrackstr, vtracknames)
    if (is.na(idx2)) {
        if (!is.na(match(vtrackstr, get("GTRACKS", envir = .misha)))) {
            stop(sprintf("Track %s already exists", vtrackstr), call. = FALSE)
        }

        if (!is.na(match(vtrackstr, get("GINTERVS", envir = .misha)))) {
            stop(sprintf("Interval %s already exists", vtrackstr), call. = FALSE)
        }

        if (.ggetOption(".gautocompletion", FALSE) && exists(vtrackstr, envir = .misha)) {
            stop(sprintf("Variable \"%s\" shadows the name of identically named virtual track.\nPlease remove this variable from the environment or switch off autocompletion mode.", vtrackstr), call. = FALSE)
        }

        vtracknames <- c(vtracknames, vtrackstr)
        idx2 <- length(vtracknames)
    }

    gvtracks[[idx1]][[idx2]] <- var

    if (!is.null(var)) {
        names(gvtracks[[idx1]]) <- vtracknames

        envir <- .misha_env()
        assign("GVTRACKS", gvtracks, envir$.misha)
        .gcall("gcheck_vtrack", vtrackstr, envir)
    }

    success <- FALSE
    old.gvtracks <- NULL
    if (exists("GVTRACKS", envir = .misha)) {
        old.gvtracks <- get("GVTRACKS", envir = .misha)
    }

    success <- FALSE
    tryCatch({
        assign("GVTRACKS", gvtracks, envir = .misha)
        success <- TRUE
    })
}


#' Creates a new virtual track
#'
#' Creates a new virtual track.
#'
#' This function creates a new virtual track named 'vtrack' with the given
#' source, function and parameters. 'src' can be either a track or intervals
#' (1D or 2D). Use the following table for a reference of all valid source,
#' function and parameters combinations:
#'
#' \emph{src = [Track], func = "avg", params = NULL} \cr Average track value in
#' iterator interval.
#'
#' \emph{src = [Track], func = "max", params = NULL} \cr Maximal track value in
#' iterator interval.
#'
#' \emph{src = [Track], func = "min", params = NULL} \cr Minimal track value in
#' iterator interval.
#'
#' \emph{src = ['Dense' / 'Sparse' / 'Array' track], func = "nearest", params =
#' NULL} \cr Mean track value in iterator interval. If there are no track
#' values covered by an iterator interator (can occur only in 'Sparse' track),
#' the nearest track value is returned.
#'
#' \emph{src = ['Dense' / 'Sparse' / 'Array' track], func = "stddev", params =
#' NULL} \cr Unbiased standard deviation of track values in iterator interval.
#'
#' \emph{src = ['Dense' / 'Sparse' / 'Array' track], func = "sum", params =
#' NULL} \cr Sum of track values in iterator interval.
#'
#' \emph{src = ['Dense' / 'Sparse' / 'Array' track], func = "quantile", params
#' = [Percentile in the range of [0, 1]]} \cr Quantile of track values in
#' iterator interval.
#'
#' \emph{src = ['Dense' track], func = "global.percentile", params = NULL} \cr
#' Percentile of an average track value in iterator interval relatively to all
#' values of the track.
#'
#' \emph{src = ['Dense' track], func = "global.percentile.max", params = NULL}
#' \cr Percentile of a maximal track value in iterator interval relatively to
#' all values of the track.
#'
#' \emph{src = ['Dense' track], func = "global.percentile.min", params = NULL}
#' \cr Percentile of a minimal track value in iterator interval relatively to
#' all values of the track.
#'
#' \emph{src = [2D track], func = "area", params = NULL} \cr Area covered by
#' iterator interval.
#'
#' \emph{src = [2D track], func = "weighted.sum", params = NULL} \cr Weighted
#' sum of values where each weight equals to the intersection area between the
#' iterator interval and the rectangle containing the value.
#'
#' \emph{src = [1D intervals], func = "distance", params = [Minimal distance
#' from center (default: 0)]} \cr Given the center 'C' of the current iterator
#' interval returns 'DC * X/2', where 'DC' is the normalized distance to the
#' center of the interval that contains 'C', and 'X' is the value of the
#' parameter. If no interval contains 'C' the resulted value is 'D + XXX/2'
#' where 'D' is the distance between 'C' and the edge of the closest interval.
#' Distance can be positive or negative depending on the position of the
#' coordinate relative to the interval and the strand (-1 or 1) of the
#' interval. Distance is always positive if 'strand' is '0' or if 'strand'
#' column is missing. Distance is 'NA' if no intervals exist for the current
#' chromosome.
#'
#' \emph{src = [1D intervals], func = "distance.center", params = NULL} \cr
#' Given the center 'C' of the current iterator interval returns 'NaN' if 'C'
#' is outside of the intervals, otherwise returns the distance between 'C' and
#' the center of the closest interval. Distance can be positive or negative
#' depending on the position of the coordinate relative to the interval and the
#' strand (-1 or 1) of the interval. Distance is always positive if 'strand' is
#' '0' or if 'strand' column is missing.
#'
#' \emph{src = [1D intervals], func = "coverage", params = NULL} \cr
#' For each iterator interval, calculates the fraction of its length that is covered by the
#' source intervals. Returns a value between 0 and 1. For example, if an iterator interval is [100,200]
#' and the source intervals cover positions 120-140 and 160-170, the coverage would be 0.3
#' ((20 + 10) / 100 = 0.3). Overlapping source intervals are first unified.
#'
#' \emph{src = [1D intervals], func = "neighbor.count", params = [Max distance >= 0]} \cr
#' Returns, for each iterator interval, the number of source intervals whose edge-to-edge distance
#' to the iterator interval is <= params. Equivalent to counting overlaps with source intervals expanded
#' by params on both sides. Overlapping sources are NOT unified (multiplicity preserved).
#'
#' \emph{func = "pwm", params = list(pssm = matrix, bidirect = TRUE,
#' prior = 0.01, extend = TRUE, spat_factor = NULL, spat_bin = NULL,
#' spat_min = NULL, spat_max = NULL)} \cr
#' Calculates total log-likelihood score of DNA sequence against PSSM.
#' Uses log-sum-exp over all positions. For bidirect=TRUE, scans both strands and
#' combines them per genomic start position (a per-position union).
#' Prior adds pseudocounts. The extend=TRUE parameter (default) pads the fetched
#' sequence by extending the END coordinate by (motif_length-1) bases for all strand modes.
#' This allows motifs anchored inside the iterator to be evaluated even when they extend
#' beyond the iterator boundary. At genomic position i, we always need sequence [i, i+motif_len),
#' regardless of strand (reverse strand computes RC of the same genomic positions). Only motif
#' anchors whose 0-based start lies inside the iterator contribute; the extra sequence provides
#' context only. With extend=FALSE, only motifs fully contained within the interval are scored.
#' Optional spatial weighting allows position-dependent weights.
#'
#' \emph{func = "pwm.max", params = list(pssm = matrix, bidirect = TRUE,
#' prior = 0.01, extend = TRUE, spat_factor = NULL, spat_bin = NULL,
#' spat_min = NULL, spat_max = NULL)} \cr
#' Returns maximum log-likelihood score of best PSSM match.
#' When \code{bidirect=TRUE}, both strands are scanned and the reported value is the
#' \emph{single best-strand} score after spatial weighting—the two strand scores are
#' compared, not summed. (This is a \emph{per-position union}, not per-strand accumulation.)
#' The \code{strand} parameter is ignored when \code{bidirect=TRUE}.
#' Prior adds pseudocounts. The \code{extend=TRUE} parameter (default) allows scoring
#' motifs whose start position falls within the interval, even if the motif extends
#' beyond the interval boundary. Neutral characters (\code{N}, \code{n}, \code{*} by default)
#' are scored with the mean log-probability of each PSSM column on both strands, so the
#' same penalty applies regardless of orientation.
#' Optional spatial weighting allows position-dependent weights.
#'
#' \emph{func = "pwm.max.pos", params = list(pssm = matrix, bidirect = TRUE,
#' prior = 0.01, extend = TRUE, spat_factor = NULL, spat_bin = NULL,
#' spat_min = NULL, spat_max = NULL)} \cr
#' Returns 1-based position of best PSSM match.
#' This position is 1-based relative to the start of the final scan window, which is determined after applying any iterator shifts (e.g., \code{sshift} from \code{\link{gvtrack.iterator}}) or the \code{spat_min} parameter. In case of ties, the position of the first match (most 5' / smallest coordinate) encountered during the scan is returned. The forward strand is checked before the reverse strand at each position, so it will win in a direct tie at the same location.
#' If bidirect=TRUE, the position would be positive if the best hit was at the
#' forward strand, and negative if it was at the reverse strand. When strand is
#' -1 the position is still according to the forward strand, but the hit is at
#' the end of the match.
#' Prior adds pseudocounts, extend=TRUE allows boundary scoring.
#' Optional spatial weighting allows position-dependent weights.
#'
#' \emph{func = "pwm.count", params = list(pssm = matrix, score.thresh = 0,
#' bidirect = TRUE, prior = 0.01, extend = TRUE, strand = 1)} \cr
#' Counts motif hits with score >= threshold inside each interval.
#' For \code{bidirect=FALSE}, only the strand given by \code{strand} is checked.
#' For \code{bidirect=TRUE}, both strands are evaluated at each position and the
#' \emph{combined} score is tested against the threshold. The default combination is
#' log-sum-exp (LSE), consistent with \code{pwm}. Each position contributes at most 1
#' to the count (per-position union), not a per-strand sum.
#' Returns the total number of passing positions and only counts anchors whose genomic
#' start lies inside the iterator (same 0-based half-open interval semantics as misha tracks).
#' Prior adds pseudocounts; \code{extend=TRUE} pads the fetched sequence exactly as in \code{pwm},
#' allowing boundary-spanning motifs on any strand to be evaluated without inflating the iterator
#' length. If spatial weights are provided, the threshold is applied to the spatially weighted
#' combined score.
#'
#' For all PWM functions:
#' \itemize{
#'   \item pssm: Position-specific scoring matrix (matrix or data frame with columns A,C,G,T containing frequencies; additional columns are allowed and will be ignored)
#'   \item bidirect: If TRUE, scans both strands and combines them per position (union).
#'     When TRUE, \code{strand} is ignored. If FALSE, only the strand given by \code{strand}
#'     is scanned (default: TRUE).
#'   \item prior: Pseudocount added to frequencies (default: 0.01). Set to 0 for no pseudocounts.
#'   \item extend: If TRUE, extends the fetched sequence so boundary-anchored motifs still
#'     have enough context (default: TRUE). All strand modes extend the END coordinate by
#'     (motif_length-1) bases. This is because at genomic position i, we always need sequence
#'     [i, i+motif_len) regardless of strand. Only anchors whose genomic start lies inside the
#'     iterator are scored/counted, regardless of the extra sequence fetched.
#'   \item neutral characters: By default \code{N}, \code{n}, and \code{*} are treated as unknown bases
#' and contribute the average log-probability of the corresponding PSSM column on both strands.
#'   \item strand: If 1, scans forward strand; if -1, scans reverse strand (default: 1).
#'     Ignored when \code{bidirect=TRUE}. For \code{pwm.max.pos}, when \code{strand == 1},
#'     the position of the best match is at the beginning of the match; when \code{strand == -1},
#'     the position is at the end of the match.
#'   \item score.thresh: Score threshold for pwm.count (default: 0). Only positions with
#' log-likelihood >= score.thresh are counted.
#'   \item spat_factor: Optional numeric vector of positive spatial weights (one per bin).
#' Weights are applied in log-space: weighted_score = log_likelihood + log(spat_factor[bin]).
#' If NULL (default), no spatial weighting is applied.
#'   \item spat_bin: Integer bin size in base pairs (required if spat_factor provided).
#' Must be > 0. Bins are 0-indexed from scan start: bin = floor(position / spat_bin).
#'   \item spat_min: Optional 1-based position defining start of scanning window (default: 1).
#' Use with spat_max to restrict scanning to a specific region.
#'   \item spat_max: Optional 1-based position defining end of scanning window (default: interval length).
#' Actual last scanned position may be earlier to accommodate motif length.
#' }
#'
#' \strong{Spatial Weighting:}
#' Enables position-dependent weighting for modeling positional biases. Bins are 0-indexed from the
#' scan start. When using gvtrack.iterator() shifts (e.g., sshift=-50, eshift=50), bins index from
#' the expanded scan window start, not the original interval. Both strands use the same bin at each
#' genomic position. Positions beyond the last bin use the last bin's weight. If the window size is
#' not divisible by spat_bin, the last bin will be smaller (e.g., 500bp window with 40bp bins creates
#' bins 0-11 of 40bp each, plus bin 12 of 20bp). Use spat_min and spat_max to restrict scanning to a
#' range divisible by spat_bin if needed.
#'
#' PWM parameters are accepted as list or individual parameters (see examples).
#'
#' \emph{func = "kmer.count", params = list(kmer = "ACGT", extend = TRUE, strand = 0)} \cr
#' Counts occurrences of the specified kmer in each interval. The extend=TRUE
#' parameter (default) allows counting kmers whose start position falls within
#' the interval, even if the kmer extends beyond the interval boundary. With
#' extend=FALSE, only kmers fully contained within the interval are counted.
#' The strand parameter can be 1 (forward strand), -1 (reverse strand), or 0 (both strands).
#'
#' \emph{func = "kmer.frac", params = list(kmer = "ACGT", extend = TRUE, strand = 0)} \cr
#' Calculates the fraction of possible positions in each interval that contain
#' the specified kmer. The extend=TRUE parameter (default) allows counting kmers
#' whose start position falls within the interval, even if the kmer extends beyond
#' the interval boundary. With extend=FALSE, only kmers fully contained within the
#' interval are counted. The strand parameter can be 1 (forward strand), -1
#' (reverse strand), or 0 (both strands).
#'
#' For kmer functions:
#' \itemize{
#'   \item kmer: The DNA sequence to count (case-insensitive)
#'   \item extend: If TRUE (default), counts kmers whose start position is within
#'         the interval, even if they extend beyond. If FALSE, only counts kmers
#'         fully contained within the interval.
#'   \item strand: If 1, counts kmers on forward strand; if -1, counts kmers on reverse strand. If
#'  0, counts kmers on both strands. Default is 0.
#' }
#'
#' Kmer parameters are accepted as list or individual parameters (see examples).
#' Note that for palindromic kmers, setting strand to 1 or -1 is recommended to avoid double counting.
#'
#' Modify iterator behavior with 'gvtrack.iterator' or 'gvtrack.iterator.2d'.
#'
#' @param vtrack virtual track name
#' @param src source (track/intervals). NULL for PWM functions
#' @param func function name (see above)
#' @param params function parameters (see above)
#' @param ... additional PWM parameters
#' @return None.
#' @seealso \code{\link{gvtrack.info}}, \code{\link{gvtrack.iterator}},
#' \code{\link{gvtrack.iterator.2d}}, \code{\link{gvtrack.array.slice}},
#' \code{\link{gvtrack.ls}}, \code{\link{gvtrack.rm}}
#' @keywords ~virtual
#' @examples
#' \dontshow{
#' options(gmax.processes = 2)
#' }
#'
#' gdb.init_examples()
#'
#' gvtrack.create("vtrack1", "dense_track", "max")
#' gvtrack.create("vtrack2", "dense_track", "quantile", 0.5)
#' gextract("dense_track", "vtrack1", "vtrack2",
#'     gintervals(1, 0, 10000),
#'     iterator = 1000
#' )
#'
#' gvtrack.create("vtrack3", "dense_track", "global.percentile")
#' gvtrack.create("vtrack4", "annotations", "distance")
#' gdist(
#'     "vtrack3", seq(0, 1, l = 10), "vtrack4",
#'     seq(-500, 500, 200)
#' )
#'
#' gvtrack.create("cov", "annotations", "coverage")
#' gextract("cov", gintervals(1, 0, 1000), iterator = 100)
#'
#' pssm <- matrix(
#'     c(
#'         0.7, 0.1, 0.1, 0.1, # Example PSSM
#'         0.1, 0.7, 0.1, 0.1,
#'         0.1, 0.1, 0.7, 0.1,
#'         0.1, 0.1, 0.7, 0.1,
#'         0.1, 0.1, 0.7, 0.1,
#'         0.1, 0.1, 0.7, 0.1
#'     ),
#'     ncol = 4, byrow = TRUE
#' )
#' colnames(pssm) <- c("A", "C", "G", "T")
#' gvtrack.create(
#'     "motif_score", NULL, "pwm",
#'     list(pssm = pssm, bidirect = TRUE, prior = 0.01)
#' )
#' gvtrack.create("max_motif_score", NULL, "pwm.max",
#'     pssm = pssm, bidirect = TRUE, prior = 0.01
#' )
#' gvtrack.create("max_motif_pos", NULL, "pwm.max.pos",
#'     pssm = pssm
#' )
#' gextract(
#'     c(
#'         "dense_track", "motif_score", "max_motif_score",
#'         "max_motif_pos"
#'     ),
#'     gintervals(1, 0, 10000),
#'     iterator = 500
#' )
#'
#' # Kmer counting examples
#' gvtrack.create("cg_count", NULL, "kmer.count", kmer = "CG", strand = 1)
#' gvtrack.create("cg_frac", NULL, "kmer.frac", kmer = "CG", strand = 1)
#' gextract(c("cg_count", "cg_frac"), gintervals(1, 0, 10000), iterator = 1000)
#'
#' gvtrack.create("at_pos", NULL, "kmer.count", kmer = "AT", strand = 1)
#' gvtrack.create("at_neg", NULL, "kmer.count", kmer = "AT", strand = -1)
#' gvtrack.create("at_both", NULL, "kmer.count", kmer = "AT", strand = 0)
#' gextract(c("at_pos", "at_neg", "at_both"), gintervals(1, 0, 10000), iterator = 1000)
#'
#' # GC content
#' gvtrack.create("g_frac", NULL, "kmer.frac", kmer = "G")
#' gvtrack.create("c_frac", NULL, "kmer.frac", kmer = "C")
#' gextract("g_frac + c_frac", gintervals(1, 0, 10000),
#'     iterator = 1000,
#'     colnames = "gc_content"
#' )
#'
#' # Spatial PWM examples
#' # Create a PWM with higher weight in the center of intervals
#' pssm <- matrix(
#'     c(
#'         0.7, 0.1, 0.1, 0.1,
#'         0.1, 0.7, 0.1, 0.1,
#'         0.1, 0.1, 0.7, 0.1,
#'         0.1, 0.1, 0.1, 0.7
#'     ),
#'     ncol = 4, byrow = TRUE
#' )
#' colnames(pssm) <- c("A", "C", "G", "T")
#'
#' # Spatial factors: low weight at edges, high in center
#' # For 200bp intervals with 40bp bins: bins 0, 40, 80, 120, 160
#' spatial_weights <- c(0.5, 1.0, 2.0, 1.0, 0.5)
#'
#' gvtrack.create(
#'     "spatial_pwm", NULL, "pwm",
#'     list(
#'         pssm = pssm,
#'         bidirect = TRUE,
#'         spat_factor = spatial_weights,
#'         spat_bin = 40L
#'     )
#' )
#'
#' # Compare with non-spatial PWM
#' gvtrack.create(
#'     "regular_pwm", NULL, "pwm",
#'     list(pssm = pssm, bidirect = TRUE)
#' )
#'
#' gextract(c("spatial_pwm", "regular_pwm"),
#'     gintervals(1, 0, 10000),
#'     iterator = 200
#' )
#'
#' # Using spatial parameters with iterator shifts
#' gvtrack.create(
#'     "spatial_extended", NULL, "pwm.max",
#'     pssm = pssm,
#'     spat_factor = c(0.5, 1.0, 2.0, 2.5, 2.0, 1.0, 0.5),
#'     spat_bin = 40L
#' )
#' # Scan window will be 280bp (100bp + 2*90bp)
#' gvtrack.iterator("spatial_extended", sshift = -90, eshift = 90)
#' gextract("spatial_extended", gintervals(1, 0, 10000), iterator = 100)
#'
#' # Using spat_min/spat_max to restrict scanning to a window
#' # For 500bp intervals, scan only positions 30-470 (440bp window)
#' gvtrack.create(
#'     "window_pwm", NULL, "pwm",
#'     pssm = pssm,
#'     bidirect = TRUE,
#'     spat_min = 30, # 1-based position
#'     spat_max = 470 # 1-based position
#' )
#' gextract("window_pwm", gintervals(1, 0, 10000), iterator = 500)
#'
#' # Combining spatial weighting with window restriction
#' # Scan positions 50-450 with spatial weights favoring the center
#' gvtrack.create(
#'     "window_spatial_pwm", NULL, "pwm",
#'     pssm = pssm,
#'     bidirect = TRUE,
#'     spat_factor = c(0.5, 1.0, 2.0, 2.5, 2.0, 1.0, 0.5, 1.0, 0.5, 0.5),
#'     spat_bin = 40L,
#'     spat_min = 50,
#'     spat_max = 450
#' )
#' gextract("window_spatial_pwm", gintervals(1, 0, 10000), iterator = 500)
#' @export gvtrack.create
gvtrack.create <- function(vtrack = NULL, src = NULL, func = NULL, params = NULL, ...) {
    if (is.null(substitute(vtrack))) {
        stop("Usage: gvtrack.create(vtrack, src, func = NULL, params = NULL, ...)", call. = FALSE)
    }
    if (is.null(substitute(src)) && !(func %in% c("pwm", "pwm.max", "pwm.max.pos", "pwm.count", "kmer.count", "kmer.frac"))) {
        stop("Usage: gvtrack.create(vtrack, src, func = NULL, params = NULL, ...)", call. = FALSE)
    }

    .gcheckroot()

    if (!is.null(func) && func %in% c("pwm", "pwm.max", "pwm.max.pos", "pwm.count")) {
        dots <- list(...)

        if (!is.null(params)) {
            # params list
            if (!is.list(params) || !("pssm" %in% names(params))) {
                stop("pwm function requires a list with at least 'pssm' matrix parameter")
            }
            dots <- params
        }

        if (!("pssm" %in% names(dots))) {
            stop("pwm function requires a 'pssm' matrix parameter")
        }
        pssm <- dots$pssm
        bidirect <- if (!is.null(dots$bidirect)) dots$bidirect else TRUE
        prior <- if (!is.null(dots$prior)) dots$prior else 0.01
        extend <- if (!is.null(dots$extend)) dots$extend else TRUE
        strand <- if (!is.null(dots$strand)) dots$strand else 1

        # Optional spatial parameters
        spat_factor <- dots$spat_factor
        spat_bin <- dots$spat_bin
        spat_min <- dots$spat_min
        spat_max <- dots$spat_max

        # Optional score threshold for pwm.count
        score.thresh <- if (!is.null(dots$score.thresh)) dots$score.thresh else 0

        pssm <- .coerce_pssm_matrix(
            pssm,
            numeric_msg = "PSSM must be a numeric matrix or data frame with numeric columns",
            ncol_msg = "PSSM must have columns named A, C, G, T",
            colnames_msg = "PSSM must have columns named A, C, G, T"
        )

        if (!is.numeric(prior) || prior < 0 || prior > 1) {
            stop("prior must be a number between 0 and 1")
        }

        if (!is.logical(bidirect)) {
            stop("bidirect must be TRUE or FALSE")
        }

        if (!is.logical(extend)) {
            stop("extend must be TRUE or FALSE")
        }

        if (strand != 1 && strand != -1) {
            stop("strand must be 1 or -1")
        }

        # Validate spatial parameters if provided
        if (!is.null(spat_factor)) {
            if (!is.numeric(spat_factor) || any(spat_factor <= 0)) {
                stop("spat_factor must be a numeric vector with all positive values")
            }
            if (is.null(spat_bin)) {
                spat_bin <- 1L
            }
            if (!is.numeric(spat_bin) || spat_bin <= 0) {
                stop("spat_bin must be a positive integer")
            }
            spat_bin <- as.integer(spat_bin)
        }

        # Set params with processed values
        params <- list(
            pssm = pssm,
            bidirect = bidirect,
            prior = prior,
            extend = extend,
            strand = strand,
            score.thresh = score.thresh
        )

        # Handle spat_min/spat_max coordinate conversion (independent of spatial factors)
        # Prego always trims sequences when spat_min/spat_max are provided
        # Misha uses scanning ranges, so we need to convert coordinates
        if (!is.null(spat_min)) {
            # Convert from 1-based R indexing to 0-based C++ indexing
            params$spat_min <- as.integer(spat_min - 1)
        }
        if (!is.null(spat_max)) {
            # Adjust spat_max to account for motif length and convert to 0-based indexing
            # spat_max defines the last base of the scanning window (1-based), but we need
            # the last valid start position for the motif (0-based)
            motif_length <- nrow(pssm)
            # First convert to 0-based, then adjust for motif length
            adjusted_spat_max <- (spat_max - 1) - motif_length + 1
            params$spat_max <- as.integer(adjusted_spat_max)
        }

        # Add spatial weighting parameters if provided
        if (!is.null(spat_factor)) {
            params$spat_factor <- spat_factor
            params$spat_bin <- spat_bin
        }
    } else if (!is.null(func) && func %in% c("kmer.count", "kmer.frac")) {
        # Check for kmer parameter
        dots <- list(...)

        if (!is.null(params)) {
            # Handle as list or string parameter
            if (is.list(params)) {
                # params list
                if (!is.list(params) || !("kmer" %in% names(params))) {
                    stop("kmer function requires a list with at least 'kmer' parameter")
                }
                kmer_params <- params
            }
        } else if ("kmer" %in% names(dots)) {
            # Use named parameters
            kmer_params <- dots
        } else {
            stop("kmer functions require a 'kmer' parameter")
        }

        kmer_params$extend <- if (!is.null(kmer_params$extend)) kmer_params$extend else TRUE
        kmer_params$strand <- if (!is.null(kmer_params$strand)) kmer_params$strand else 0

        # Validate required kmer parameter
        if (!("kmer" %in% names(kmer_params)) || !is.character(kmer_params$kmer) || length(kmer_params$kmer) != 1) {
            stop("kmer parameter must be a single string")
        }

        if (nchar(kmer_params$kmer) == 0) {
            stop("kmer sequence cannot be empty")
        }

        if (kmer_params$kmer == grevcomp(kmer_params$kmer) && kmer_params$strand == 0) {
            warning(paste0("kmer sequence '", kmer_params$kmer, "' is palindromic, please set strand to 1 or -1 to avoid double counting"))
        }

        params <- kmer_params
    } else if (!is.null(func) && func == "neighbor.count") {
        if (is.null(params)) {
            params <- 0
        }

        if (!is.numeric(params) || length(params) != 1 || is.na(params)) {
            stop("neighbor.count requires 'params' to be a single numeric value")
        }

        params <- as.numeric(params)

        if (params < 0) {
            stop("neighbor.count requires 'params' to be non-negative")
        }
    }

    vtrackstr <- do.call(.gexpr2str, list(substitute(vtrack)), envir = parent.frame())
    srcstr <- do.call(.gexpr2str, list(substitute(src)), envir = parent.frame())

    if (!is.na(match(vtrackstr, get("GTRACKS", envir = .misha)))) {
        stop(sprintf("Cannot create virtual track: regular track named %s already exists", vtrackstr), call. = FALSE)
    }

    if (!is.na(match(vtrackstr, get("GINTERVS", envir = .misha)))) {
        stop(sprintf("Cannot create virtual track: intervals named %s already exists", vtrackstr), call. = FALSE)
    }

    if (vtrackstr != make.names(vtrackstr)) {
        stop(sprintf("\"%s\" is not a syntactically valid name for a variable", vtrackstr), call. = FALSE)
    }

    var <- list()
    if (is.character(srcstr) && !is.na(match(srcstr, get("GTRACKS", envir = .misha)))) {
        var$src <- srcstr
    } else {
        var$src <- src
    }
    var$func <- func
    var$params <- params

    .gvtrack.set(vtrackstr, var)

    invisible(vtrackstr)
}


#' Returns the definition of a virtual track
#'
#' Returns the definition of a virtual track.
#'
#' This function returns the internal representation of a virtual track.
#'
#' @param vtrack virtual track name
#' @return Internal representation of a virtual track.
#' @seealso \code{\link{gvtrack.create}}
#' @keywords ~virtual
#' @examples
#' \dontshow{
#' options(gmax.processes = 2)
#' }
#'
#' gdb.init_examples()
#' gvtrack.create("vtrack1", "dense_track", "max")
#' gvtrack.info("vtrack1")
#'
#' @export gvtrack.info
gvtrack.info <- function(vtrack = NULL) {
    if (is.null(substitute(vtrack))) {
        stop("Usage: gvtrack.info(vtrack)", call. = FALSE)
    }
    .gcheckroot()

    vtrackstr <- do.call(.gvtrack, list(substitute(vtrack)), envir = parent.frame())
    info <- .gvtrack.get(vtrackstr)

    # If filter is present, add filter statistics
    if (!is.null(info$filter)) {
        filter_info <- tryCatch(
            {
                .gcall("C_get_filter_info", info$filter, .misha_env())
            },
            error = function(e) {
                NULL
            }
        )

        if (!is.null(filter_info)) {
            info$filter_stats <- filter_info
        }
    }

    info
}


#' Defines modification rules for a one-dimensional iterator in a virtual track
#'
#' Defines modification rules for a one-dimensional iterator in a virtual
#' track.
#'
#' This function defines modification rules for one-dimensional iterator
#' intervals in a virtual track.
#'
#' 'dim' converts a 2D iterator interval (chrom1, start1, end1, chrom2, start2,
#' end2) to a 1D interval. If 'dim' is '1' the interval is converted to
#' (chrom1, start1, end1). If 'dim' is '2' the interval is converted to
#' (chrom2, start2, end2). If 1D iterator is used 'dim' must be set to 'NULL'
#' or '0' (meaning: no conversion is made).
#'
#' Iterator interval's 'start' coordinate is modified by adding 'sshift'.
#' Similarly 'end' coordinate is altered by adding 'eshift'.
#'
#' @param vtrack virtual track name
#' @param dim use 'NULL' or '0' for 1D iterators. '1' converts 2D iterator to
#' (chrom1, start1, end1) , '2' converts 2D iterator to (chrom2, start2, end2)
#' @param sshift shift of 'start' coordinate
#' @param eshift shift of 'end' coordinate
#' @return None.
#' @seealso \code{\link{gvtrack.create}}, \code{\link{gvtrack.iterator.2d}}
#' @keywords ~virtual
#' @examples
#' \dontshow{
#' options(gmax.processes = 2)
#' }
#'
#' gdb.init_examples()
#'
#' gvtrack.create("vtrack1", "dense_track")
#' gvtrack.iterator("vtrack1", sshift = 200, eshift = 200)
#' gextract("dense_track", "vtrack1", gintervals(1, 0, 500))
#'
#' gvtrack.create("vtrack2", "dense_track")
#' gvtrack.iterator("vtrack2", dim = 1)
#' gextract("vtrack2", gintervals.2d(1, 0, 1000, 1, 0, -1),
#'     iterator = "rects_track"
#' )
#'
#' @export gvtrack.iterator
gvtrack.iterator <- function(vtrack = NULL, dim = NULL, sshift = 0, eshift = 0) {
    if (is.null(substitute(vtrack))) {
        stop("Usage: gvtrack.iterator(vtrack, dim = NULL, sshift = 0, eshift = 0)", call. = FALSE)
    }
    .gcheckroot()

    vtrackstr <- do.call(.gvtrack, list(substitute(vtrack)), envir = parent.frame())

    var <- .gvtrack.get(vtrackstr)
    itr <- list()
    itr$type <- "1d"
    itr$dim <- dim
    itr$sshift <- sshift
    itr$eshift <- eshift
    var$itr <- itr
    .gvtrack.set(vtrackstr, var)
    retv <- NULL
}


#' Defines modification rules for a two-dimensional iterator in a virtual track
#'
#' Defines modification rules for a two-dimensional iterator in a virtual
#' track.
#'
#' This function defines modification rules for one-dimensional iterator
#' intervals in a virtual track.
#'
#' Iterator interval's 'start1' coordinate is modified by adding 'sshift1'.
#' Similarly 'end1', 'start2', 'end2' coordinates are altered by adding
#' 'eshift1', 'sshift2' and 'eshift2' accordingly.
#'
#' @param vtrack virtual track name
#' @param sshift1 shift of 'start1' coordinate
#' @param eshift1 shift of 'end1' coordinate
#' @param sshift2 shift of 'start2' coordinate
#' @param eshift2 shift of 'end2' coordinate
#' @return None.
#' @seealso \code{\link{gvtrack.create}}, \code{\link{gvtrack.iterator}}
#' @keywords ~virtual
#' @examples
#' \dontshow{
#' options(gmax.processes = 2)
#' }
#'
#' gdb.init_examples()
#' gvtrack.create("vtrack1", "rects_track")
#' gvtrack.iterator.2d("vtrack1", sshift1 = 1000, eshift1 = 2000)
#' gextract(
#'     "rects_track", "vtrack1",
#'     gintervals.2d(1, 0, 5000, 2, 0, 5000)
#' )
#'
#' @export gvtrack.iterator.2d
gvtrack.iterator.2d <- function(vtrack = NULL, sshift1 = 0, eshift1 = 0, sshift2 = 0, eshift2 = 0) {
    if (is.null(substitute(vtrack))) {
        stop("Usage: gvtrack.iterator.2d(vtrack, sshift1 = 0, eshift1 = 0, sshift2 = 0, eshift2 = 0)", call. = FALSE)
    }
    .gcheckroot()

    vtrackstr <- do.call(.gvtrack, list(substitute(vtrack)), envir = parent.frame())

    var <- .gvtrack.get(vtrackstr)
    itr <- list()
    itr$type <- "2d"
    itr$sshift1 <- sshift1
    itr$eshift1 <- eshift1
    itr$sshift2 <- sshift2
    itr$eshift2 <- eshift2
    var$itr <- itr
    .gvtrack.set(vtrackstr, var)
    retv <- NULL
}


#' Returns a list of virtual track names
#'
#' Returns a list of virtual track names.
#'
#' This function returns a list of virtual tracks that exist in current R
#' environment that match the pattern (see 'grep'). If called without any
#' arguments all virtual tracks are returned.
#'
#' @param pattern,ignore.case,perl,fixed,useBytes see 'grep'
#' @return An array that contains the names of virtual tracks.
#' @seealso \code{\link{grep}}, \code{\link{gvtrack.create}},
#' \code{\link{gvtrack.rm}}
#' @keywords ~virtual ~ls
#' @examples
#' \dontshow{
#' options(gmax.processes = 2)
#' }
#'
#' gdb.init_examples()
#' gvtrack.create("vtrack1", "dense_track", "max")
#' gvtrack.create("vtrack2", "dense_track", "quantile", 0.5)
#' gvtrack.ls()
#' gvtrack.ls(pattern = "*2")
#'
#' @export gvtrack.ls
gvtrack.ls <- function(pattern = "", ignore.case = FALSE, perl = FALSE, fixed = FALSE, useBytes = FALSE) {
    .gcheckroot()

    if (!exists("GVTRACKS", envir = .misha)) {
        return(NULL)
    }

    gvtracks <- get("GVTRACKS", envir = .misha)
    gwds <- names(gvtracks)
    if (!is.list(gvtracks) || (length(gvtracks) && !is.character(gwds)) || length(gvtracks) != length(gwds)) {
        stop("Invalid format of GVTRACKS variable.\nTo continue working with virtual tracks please remove this variable from the environment.", call. = FALSE)
    }

    gwd <- get("GWD", envir = .misha)
    idx <- match(gwd, gwds)
    if (is.na(idx)) {
        return(NULL)
    }

    vtracks <- gvtracks[[idx]]
    vtracknames <- names(vtracks)
    if (!is.list(vtracks) || (length(vtracks) && !is.character(vtracknames)) || length(vtracks) != length(vtracknames)) {
        stop("Invalid format of GVTRACKS variable.\nTo continue working with virtual tracks please remove this variable from the environment.", call. = FALSE)
    }

    if (!length(vtracks)) {
        return(NULL)
    }

    if (pattern != "") {
        grep(pattern, vtracknames, value = TRUE, ignore.case = ignore.case, perl = perl, fixed = fixed, useBytes = useBytes)
    } else {
        vtracknames
    }
}


#' Deletes a virtual track
#'
#' Deletes a virtual track.
#'
#' This function deletes a virtual track from current R environment.
#'
#' @param vtrack virtual track name
#' @return None.
#' @seealso \code{\link{gvtrack.create}}, \code{\link{gvtrack.ls}}
#' @keywords ~virtual
#' @examples
#' \dontshow{
#' options(gmax.processes = 2)
#' }
#'
#' gdb.init_examples()
#' gvtrack.create("vtrack1", "dense_track", "max")
#' gvtrack.create("vtrack2", "dense_track", "quantile", 0.5)
#' gvtrack.ls()
#' gvtrack.rm("vtrack1")
#' gvtrack.ls()
#'
#' @export gvtrack.rm
gvtrack.rm <- function(vtrack = NULL) {
    if (is.null(substitute(vtrack))) {
        stop("Usage: gvtrack.rm(vtrack)", call. = FALSE)
    }
    .gcheckroot()

    vtrackstr <- do.call(.gvtrack, list(substitute(vtrack)), envir = parent.frame())
    .gvtrack.set(vtrackstr, NULL)
    retv <- NULL
}


#' Defines rules for a single value calculation of a virtual 'Array' track
#'
#' Defines how a single value within an interval is achieved for a virtual
#' track based on 'Array' track.
#'
#' A track (regular or virtual) used in a track expression is expected to
#' return one value for each track interval. 'Array' tracks store multiple
#' values per interval (one for each 'column') and hence if used in a track
#' expression one must define the way of how a single value should be deduced
#' from several ones.
#'
#' By default if an 'Array' track is used in a track expressions, its interval
#' value would be the average of all column values that are not NaN.
#' 'gvtrack.array.slice' allows to select specific columns and to specify the
#' function applied to their values.
#'
#' 'slice' parameter allows to choose the columns. Columns can be indicated by
#' their names or their indices. If 'slice' is 'NULL' the non-NaN values of all
#' track columns are used.
#'
#' 'func' parameter determines the function applied to the columns' values. Use
#' the following table for a reference of all valid functions and parameters
#' combinations:
#'
#' \emph{func = "avg", params = NULL} \cr Average of columns' values.
#'
#' \emph{func = "max", params = NULL} \cr Maximum of columns' values.
#'
#' \emph{func = "min", params = NULL} \cr Minimum of columns' values.
#'
#' \emph{func = "stdev", params = NULL} \cr Unbiased standard deviation of
#' columns' values.
#'
#' \emph{func = "sum", params = NULL} \cr Sum of columns' values.
#'
#' \emph{func = "quantile", params = [Percentile in the range of [0, 1]]} \cr
#' Quantile of columns' values.
#'
#' @param vtrack virtual track name
#' @param slice a vector of column names or column indices or 'NULL'
#' @param func,params see below
#' @return None.
#' @seealso \code{\link{gvtrack.create}},
#' \code{\link{gtrack.array.get_colnames}}, \code{\link{gtrack.array.extract}}
#' @keywords ~virtual ~array
#' @examples
#' \dontshow{
#' options(gmax.processes = 2)
#' }
#'
#' gdb.init_examples()
#' gvtrack.create("vtrack1", "array_track")
#' gvtrack.array.slice("vtrack1", c("col2", "col4"), "max")
#' gextract("vtrack1", gintervals(1, 0, 1000))
#'
#' @export gvtrack.array.slice
gvtrack.array.slice <- function(vtrack = NULL, slice = NULL, func = "avg", params = NULL) {
    if (is.null(substitute(vtrack))) {
        stop("Usage: gvtrack.array.slice(vtrack, slice = NULL, func = \"avg\", params = NULL)", call. = FALSE)
    }
    .gcheckroot()

    vtrackstr <- do.call(.gvtrack, list(substitute(vtrack)), envir = parent.frame())

    var <- .gvtrack.get(vtrackstr)
    .slice <- list()
    .slice$slice <- .gslice(var$src, slice)$slice
    .slice$func <- func
    .slice$params <- params
    var$slice <- .slice
    .gvtrack.set(vtrackstr, var)

    retv <- NULL
}


#' Attach or clear a genomic mask filter on a virtual track
#'
#' Attaches or clears a genomic mask filter on a virtual track. When a filter is attached,
#' the virtual track function is evaluated only over the unmasked regions (i.e., regions
#' not covered by the filter intervals).
#'
#' @param vtrack virtual track name
#' @param filter genomic mask to apply. Can be:
#'   \itemize{
#'     \item A data.frame with columns 'chrom', 'start', 'end' (intervals to mask)
#'     \item A character string naming an intervals set
#'     \item A character string naming a track (must be intervals-type track)
#'     \item A list of any combination of the above (all will be unified)
#'     \item NULL to clear the filter
#'   }
#' @details
#' The filter defines regions to \emph{exclude} from virtual track evaluation.
#' The virtual track function will be evaluated only on the complement of the filter.
#' Once a filter is attached to a virtual track, it applies to \strong{all subsequent extractions}
#' of that virtual track until explicitly cleared with \code{filter = NULL}.
#'
#' \strong{Order of Operations:}
#'
#' Filters are applied \emph{after} iterator modifiers (sshift/eshift/dim). The order is:
#' \enumerate{
#'   \item Apply iterator modifiers (gvtrack.iterator with sshift/eshift)
#'   \item Subtract mask from the modified intervals
#'   \item Evaluate virtual track function over unmasked regions
#' }
#'
#' \strong{Semantics by function type:}
#' \itemize{
#'   \item \strong{Aggregations (avg/sum/min/max/stddev/quantile):} Length-weighted over unmasked regions
#'   \item \strong{coverage:} Returns (covered length in unmasked region) / (total unmasked length)
#'   \item \strong{distance/distance.center:} Unaffected by mask (pure geometry)
#'   \item \strong{PWM/kmer:} Masked bases act as hard boundaries; matches cannot span masked regions.
#'         \strong{Important:} When \code{extend=TRUE} (the default), motifs at the boundaries of unmasked
#'         segments can use bases from the adjacent masked regions to complete the motif scoring.
#'         For example, if a 4bp motif starts at position 1998 in an unmasked region that ends at 2000,
#'         and positions 2000-2002 are masked, the motif will still be scored using the masked bases.
#'         In other words, motif matches \emph{starting positions} must be in unmasked regions,
#'         but the motif sequence itself can extend into masked regions when \code{extend=TRUE}.
#'         Set \code{extend=FALSE} to prevent any use of masked bases in scoring.
#' }
#'
#' \strong{Completely Masked Intervals:}
#' If an entire iterator interval is masked, the function returns \code{NA} (not 0).
#'
#'
#' @return None (invisibly).
#' @seealso \code{\link{gvtrack.create}}, \code{\link{gvtrack.iterator}}, \code{\link{gvtrack.info}}
#' @keywords ~virtual ~filter
#' @examples
#' \dontshow{
#' options(gmax.processes = 2)
#' }
#'
#' gdb.init_examples()
#'
#' ## Basic usage: Excluding specific regions
#' gvtrack.create("vtrack1", "dense_track", func = "avg")
#'
#' # Create intervals to mask (e.g., repetitive regions)
#' repeats <- gintervals(c(1, 1), c(100, 500), c(200, 600))
#'
#' # Attach filter - track will be evaluated excluding these regions
#' gvtrack.filter("vtrack1", filter = repeats)
#'
#' # Extract values - masked regions are excluded from calculation
#' result_filtered <- gextract("vtrack1", gintervals(1, 0, 1000))
#'
#' # Check filter info
#' gvtrack.info("vtrack1")
#'
#' # Clear the filter and compare
#' gvtrack.filter("vtrack1", filter = NULL)
#' result_unfiltered <- gextract("vtrack1", gintervals(1, 0, 1000))
#'
#' ## Using multiple filter sources (combined automatically)
#' centromeres <- gintervals(1, 10000, 15000)
#' telomeres <- gintervals(1, 0, 1000)
#' combined_mask <- list(repeats, centromeres, telomeres)
#'
#' gvtrack.filter("vtrack1", filter = combined_mask)
#' result_multi_filter <- gextract("vtrack1", gintervals(1, 0, 20000))
#'
#' ## Filters work with iterator modifiers
#' gvtrack.create("vtrack2", "dense_track", func = "sum")
#' gvtrack.filter("vtrack2", filter = repeats)
#' gvtrack.iterator("vtrack2", sshift = -50, eshift = 50)
#'
#' # Iterator shifts applied first, then mask subtracted
#' result_shifted <- gextract("vtrack2", gintervals(1, 1000, 2000), iterator = 100)
#'
#' @export
gvtrack.filter <- function(vtrack = NULL, filter = NULL) {
    if (is.null(substitute(vtrack))) {
        stop("Usage: gvtrack.filter(vtrack, filter = NULL)", call. = FALSE)
    }
    .gcheckroot()

    vtrackstr <- do.call(.gvtrack, list(substitute(vtrack)), envir = parent.frame())
    var <- .gvtrack.get(vtrackstr)

    # Clear filter if NULL
    if (is.null(filter)) {
        var$filter <- NULL
        .gvtrack.set(vtrackstr, var)
        return(invisible(NULL))
    }

    # Resolve filter to canonical intervals
    filter_df <- .resolve_filter_sources(filter)

    # Canonicalize: sort, merge overlaps, validate
    filter_df <- gintervals.canonic(filter_df, unify_touching_intervals = TRUE)

    # Generate cache key
    key <- .make_filter_key(filter_df, substitute(filter))

    # Ensure filter is compiled (calls C to register if needed)
    .ensure_filter_compiled(filter_df, key)

    # Attach filter key to vtrack
    var$filter <- key
    .gvtrack.set(vtrackstr, var)

    invisible(NULL)
}


# Helper: resolve filter sources to a unified data.frame
.resolve_filter_sources <- function(filter) {
    if (is.null(filter)) {
        return(NULL)
    }

    # If it's a list, union all elements
    if (is.list(filter) && !is.data.frame(filter)) {
        parts <- lapply(filter, .resolve_filter_sources)
        if (length(parts) == 0) {
            return(data.frame(chrom = character(0), start = integer(0), end = integer(0)))
        }
        result <- parts[[1]]
        if (length(parts) > 1) {
            for (i in 2:length(parts)) {
                result <- gintervals.union(result, parts[[i]])
            }
        }
        return(result)
    }

    # If it's a data.frame, validate columns
    if (is.data.frame(filter)) {
        if (!all(c("chrom", "start", "end") %in% names(filter))) {
            stop("Filter data frame must have columns: chrom, start, end", call. = FALSE)
        }
        return(filter[, c("chrom", "start", "end"), drop = FALSE])
    }

    # If it's a string, check if it's an intervals set or track
    if (is.character(filter) && length(filter) == 1) {
        # Check if it's an intervals set
        if (filter %in% gintervals.ls()) {
            return(gintervals.load(filter))
        }

        # Check if it's a track
        if (filter %in% gtrack.ls()) {
            track_info <- gtrack.info(filter)
            if (track_info$type %in% c("sparse", "intervals")) {
                # Extract intervals from track
                return(gextract(filter, gintervals.all()))
            } else {
                stop(sprintf("Track '%s' is not an intervals-type track", filter), call. = FALSE)
            }
        }

        stop(sprintf("Filter '%s' is neither an intervals set nor a track", filter), call. = FALSE)
    }

    stop("Filter must be a data.frame, intervals set name, track name, or list of these", call. = FALSE)
}


# Helper: generate a unique cache key for a filter
.make_filter_key <- function(filter_df, filter_expr) {
    # Get DB path and genome ID for cache key
    gwd <- get("GWD", envir = .misha)

    hash <- digest::digest(
        filter_df[, c("chrom", "start", "end")],
        algo = "xxhash64",
        serialize = TRUE
    )

    key <- sprintf(
        "filter_%s_%d_%s",
        gsub("[^a-zA-Z0-9]", "_", gwd),
        nrow(filter_df),
        hash
    )

    return(key)
}


# Helper: ensure filter is compiled in C registry
.ensure_filter_compiled <- function(filter_df, key) {
    # Call C function to register filter (idempotent - returns early if already registered)
    .gcall("C_register_filter", filter_df, key, .misha_env())
    invisible(NULL)
}
