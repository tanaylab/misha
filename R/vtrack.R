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

.set_vtrack_iterator_1d <- function(vtrackstr, dim = NULL, sshift = 0, eshift = 0) {
    var <- .gvtrack.get(vtrackstr)
    itr <- list()
    itr$type <- "1d"
    itr$dim <- dim
    itr$sshift <- sshift
    itr$eshift <- eshift
    var$itr <- itr
    .gvtrack.set(vtrackstr, var)
}


#' Creates a new virtual track
#'
#' Creates a new virtual track.
#'
#' This function creates a new virtual track named 'vtrack' with the given
#' source, function and parameters. 'src' can be either a track, intervals
#' (1D or 2D), or a data frame with intervals and a numeric value column
#' (value-based track). The tables below summarize the supported combinations.
#'
#' \strong{Value-based tracks}
#' Value-based tracks are data frames containing genomic intervals with associated
#' numeric values. They function as in-memory sparse tracks without requiring
#' track creation in the database. To create a value-based track, provide a data
#' frame with columns \code{chrom}, \code{start}, \code{end}, and one numeric
#' value column (any name is acceptable). Value-based tracks support all track-based
#' summarizer functions (e.g., \code{avg}, \code{min}, \code{max}, \code{sum},
#' \code{stddev}, \code{quantile}, \code{nearest}, \code{exists}, \code{size},
#' \code{first}, \code{last}, \code{sample}, and position functions), but do not
#' support overlapping intervals. They behave like sparse tracks in aggregation:
#' values are aggregated using count-based averaging (each interval contributes equally
#' regardless of length), not coverage-based averaging.
#'
#' \strong{Track-based summarizers}
#' \tabular{llll}{
#'   Source \tab func \tab params \tab Description \cr
#'   Track \tab avg \tab NULL \tab Average track value in the iterator interval. \cr
#'   Track (1D) \tab exists \tab vals (optional) \tab Returns 1 if any value exists (or specific vals if provided), 0 otherwise. \cr
#'   Track (1D) \tab first \tab NULL \tab First value in the iterator interval. \cr
#'   Track (1D) \tab last \tab NULL \tab Last value in the iterator interval. \cr
#'   Track \tab max \tab NULL \tab Maximum track value in the iterator interval. \cr
#'   Track \tab min \tab NULL \tab Minimum track value in the iterator interval. \cr
#'   Dense / Sparse / Array track \tab nearest \tab NULL \tab Average value inside the iterator; for sparse tracks with no samples in the interval, falls back to the closest sample outside the interval (by genomic distance). \cr
#'   Track (1D) \tab sample \tab NULL \tab Uniformly sampled source value from the iterator interval. \cr
#'   Track (1D) \tab size \tab NULL \tab Number of non-NaN values in the iterator interval. \cr
#'   Dense / Sparse / Array track \tab stddev \tab NULL \tab Unbiased standard deviation of values in the iterator interval. \cr
#'   Dense / Sparse / Array track \tab sum \tab NULL \tab Sum of values in the iterator interval. \cr
#'   Dense / Sparse / Array track \tab quantile \tab Percentile in [0, 1] \tab Quantile of values in the iterator interval. \cr
#'   Dense track \tab global.percentile \tab NULL \tab Percentile of the interval average relative to the full-track distribution. \cr
#'   Dense track \tab global.percentile.max \tab NULL \tab Percentile of the interval maximum relative to the full-track distribution. \cr
#'   Dense track \tab global.percentile.min \tab NULL \tab Percentile of the interval minimum relative to the full-track distribution. \cr
#' }
#'
#' \strong{Track position summarizers}
#' \tabular{llll}{
#'   Source \tab func \tab params \tab Description \cr
#'   Track (1D) \tab first.pos.abs \tab NULL \tab Absolute genomic coordinate of the first value. \cr
#'   Track (1D) \tab first.pos.relative \tab NULL \tab Zero-based position (relative to interval start) of the first value. \cr
#'   Track (1D) \tab last.pos.abs \tab NULL \tab Absolute genomic coordinate of the last value. \cr
#'   Track (1D) \tab last.pos.relative \tab NULL \tab Zero-based position (relative to interval start) of the last value. \cr
#'   Track (1D) \tab max.pos.abs \tab NULL \tab Absolute genomic coordinate of the maximum value inside the iterator interval. \cr
#'   Track (1D) \tab max.pos.relative \tab NULL \tab Zero-based position (relative to interval start) of the maximum value. \cr
#'   Track (1D) \tab min.pos.abs \tab NULL \tab Absolute genomic coordinate of the minimum value inside the iterator interval. \cr
#'   Track (1D) \tab min.pos.relative \tab NULL \tab Zero-based position (relative to interval start) of the minimum value. \cr
#'   Track (1D) \tab sample.pos.abs \tab NULL \tab Absolute genomic coordinate of a uniformly sampled value. \cr
#'   Track (1D) \tab sample.pos.relative \tab NULL \tab Zero-based position (relative to interval start) of a uniformly sampled value. \cr
#' }
#'
#' For \code{max.pos.relative}, \code{min.pos.relative}, \code{first.pos.relative}, \code{last.pos.relative}, \code{sample.pos.relative},
#' iterator modifiers (including \code{sshift} /
#' \code{eshift} and 1D projections generated via \code{gvtrack.iterator}) are
#' applied before the position is reported. In other words, the returned
#' coordinate is always 0-based and measured from the start of the iterator
#' interval after all modifier adjustments.
#'
#' \strong{Interval-based summarizers}
#' \tabular{llll}{
#'   Source \tab func \tab params \tab Description \cr
#'   1D intervals \tab distance \tab Minimal distance from center (default 0) \tab Signed distance using normalized formula when inside intervals, distance to edge when outside; see notes below for exact formula. \cr
#'   1D intervals \tab distance.center \tab NULL \tab Distance from iterator center to the closest interval center, \code{NA} if outside all intervals. \cr
#'   1D intervals \tab distance.edge \tab NULL \tab Edge-to-edge distance from iterator interval to closest source interval (like \code{gintervals.neighbors}); see notes below for strand handling. \cr
#'   1D intervals \tab coverage \tab NULL \tab Fraction of iterator length covered by source intervals (after unifying overlaps). \cr
#'   1D intervals \tab neighbor.count \tab Max distance (>= 0) \tab Number of source intervals whose edge-to-edge distance from the iterator interval is within params (no unification). \cr
#' }
#'
#' \strong{2D track summarizers}
#' \tabular{llll}{
#'   Source \tab func \tab params \tab Description \cr
#'   2D track \tab area \tab NULL \tab Area covered by intersections of track rectangles with the iterator interval. \cr
#'   2D track \tab weighted.sum \tab NULL \tab Weighted sum of values where each weight equals the intersection area. \cr
#' }
#'
#' \strong{Motif (PWM) summarizers}
#' \tabular{llll}{
#'   Source \tab func \tab Key params \tab Description \cr
#'   NULL (sequence) \tab pwm \tab pssm, bidirect, prior, extend, spat_* \tab Log-sum-exp score of motif likelihoods across all anchors inside the iterator interval. \cr
#'   NULL (sequence) \tab pwm.max \tab pssm, bidirect, prior, extend, spat_* \tab Maximum log-likelihood score among all anchors (per-position union across strands). \cr
#'   NULL (sequence) \tab pwm.max.pos \tab pssm, bidirect, prior, extend, spat_* \tab 1-based position of the best-scoring anchor (signed by strand when \code{bidirect = TRUE}); coordinates are always relative to the iterator interval after any \code{gvtrack.iterator()} shifts/extensions. \cr
#'   NULL (sequence) \tab pwm.count \tab pssm, score.thresh, bidirect, prior, extend, strand, spat_* \tab Count of anchors whose score exceeds \code{score.thresh} (per-position union). \cr
#' }
#'
#' \strong{K-mer summarizers}
#' \tabular{llll}{
#'   Source \tab func \tab Key params \tab Description \cr
#'   NULL (sequence) \tab kmer.count \tab kmer, extend, strand \tab Number of k-mer occurrences whose anchor lies inside the iterator interval. \cr
#'   NULL (sequence) \tab kmer.frac \tab kmer, extend, strand \tab Fraction of possible anchors within the interval that match the k-mer. \cr
#' }
#'
#' \strong{Masked sequence summarizers}
#' \tabular{llll}{
#'   Source \tab func \tab Key params \tab Description \cr
#'   NULL (sequence) \tab masked.count \tab NULL \tab Number of masked (lowercase) base pairs in the iterator interval. \cr
#'   NULL (sequence) \tab masked.frac \tab NULL \tab Fraction of base pairs in the iterator interval that are masked (lowercase). \cr
#' }
#'
#' The sections below provide additional notes for motif, interval, k-mer, and masked sequence functions.
#'
#' \strong{Motif (PWM) notes}
#' \itemize{
#'   \item \code{pssm}: Position-specific scoring matrix (matrix or data frame) with columns \code{A}, \code{C}, \code{G}, \code{T}; extra columns are ignored.
#'   \item \code{bidirect}: When TRUE (default), both strands are scanned and combined per genomic start (per-position union). The \code{strand} argument is ignored. When FALSE, only the strand specified by \code{strand} is scanned.
#'   \item \code{prior}: Pseudocount added to frequencies (default 0.01). Set to 0 to disable.
#'   \item \code{extend}: Extends the fetched sequence so boundary-anchored motifs retain full context (default TRUE). The END coordinate is padded by motif_length - 1 for all strand modes; anchors must still start inside the iterator.
#'   \item Neutral characters (\code{N}, \code{n}, \code{*}) contribute the mean log-probability of the corresponding PSSM column on both strands.
#'   \item \code{strand}: Used only when \code{bidirect = FALSE}; 1 scans the forward strand, -1 scans the reverse strand. For \code{pwm.max.pos}, strand = -1 reports the hit position at the end of the match (still relative to the forward orientation).
#'   \item \code{score.thresh}: Threshold for \code{pwm.count}. Anchors with log-likelihood >= \code{score.thresh} are counted; only one count per genomic start.
#'   \item Spatial weighting (\code{spat_factor}, \code{spat_bin}, \code{spat_min}, \code{spat_max}): optional position-dependent weights applied in log-space. Provide a positive numeric vector \code{spat_factor}; \code{spat_bin} (integer > 0) defines bin width; \code{spat_min}/\code{spat_max} restrict the scanning window.
#'   \item \code{pwm.max.pos}: Positions are reported 1-based relative to the final scan window (after iterator shifts and spatial trimming). Ties resolve to the most 5' anchor; the forward strand wins ties at the same coordinate. Values are signed when \code{bidirect = TRUE} (positive for forward, negative for reverse).
#' }
#'
#' \strong{Spatial weighting}
#' enables position-dependent weighting for modeling positional biases. Bins are 0-indexed from the
#' scan start. When using \code{gvtrack.iterator()} shifts (e.g., \code{sshift = -50}, \code{eshift = 50}), bins index from
#' the expanded scan window start, not the original interval. Both strands use the same bin at each
#' genomic position. Positions beyond the last bin reuse the final bin's weight. If the window size is
#' not divisible by \code{spat_bin}, the last bin is shorter (e.g., scanning 500 bp with 40 bp bins yields
#' bins 0-11 of 40 bp plus bin 12 of 20 bp). Use \code{spat_min} and \code{spat_max} to restrict scanning to a
#' range divisible by \code{spat_bin} if needed.
#'
#' PWM parameters can be supplied either as a single list (\code{params}) or via named arguments (see examples).
#'
#' \strong{Interval distance notes}
#'
#' \code{distance}: Given the center 'C' of the current iterator interval, returns 'DC * X/2' where 'DC' is the normalized distance to the center of the interval that contains 'C', and 'X' is the value of the parameter (default: 0). If no interval contains 'C', the result is 'D + X/2' where 'D' is the distance between 'C' and the edge of the closest interval.
#'
#' \code{distance.center}: Given the center 'C' of the current iterator interval, returns \code{NaN} if 'C' is outside of all intervals, otherwise returns the distance between 'C' and the center of the closest interval.
#'
#' \code{distance.edge}: Computes edge-to-edge distance from the iterator interval to the closest source interval, using the same calculation as \code{gintervals.neighbors}. Returns 0 for overlapping intervals. Distance sign depends on the strand column of source intervals; returns unsigned (absolute) distance if no strand column exists. Returns \code{NA} if no source intervals exist on the current chromosome.
#'
#' For \code{distance} and \code{distance.center}, distance can be positive or negative depending on the position of the coordinate relative to the interval and the strand (-1 or 1) of the interval. Distance is always positive if \code{strand = 0} or if the strand column is missing. The result is \code{NA} if no intervals exist for the current chromosome.
#'
#' \strong{Difference between distance functions:} The \code{distance} function measures from the \emph{center} of the iterator interval (a single coordinate point) to the closest \emph{edge} of source intervals when outside, or returns a normalized distance within the interval when inside. The \code{distance.center} function measures from the center of the iterator interval to the \emph{center} of source intervals. The \code{distance.edge} function measures \emph{edge-to-edge} distance between intervals, exactly like \code{gintervals.neighbors}. Use \code{distance.edge} when you need the same distance computation as \code{gintervals.neighbors} within a virtual track context.
#'
#' \strong{K-mer notes}
#' \itemize{
#'   \item \code{kmer}: DNA sequence (case-insensitive) to count.
#'   \item \code{extend}: If TRUE (default), counts kmers whose anchor lies in the interval even if the kmer extends beyond it; when FALSE, only kmers fully contained in the interval are considered.
#'   \item \code{strand}: 1 counts forward-strand occurrences, -1 counts reverse-strand occurrences, 0 counts both strands (default). For palindromic kmers, consider using 1 or -1 to avoid double counting.
#' }
#'
#' K-mer parameters can be supplied as a list or via named arguments (see examples).
#'
#' Modify iterator behavior with 'gvtrack.iterator' or 'gvtrack.iterator.2d'.
#'
#' @param vtrack virtual track name
#' @param src source (track/intervals). NULL for PWM functions. For value-based
#' tracks, provide a data frame with columns \code{chrom}, \code{start}, \code{end},
#' and one numeric value column. The data frame functions as an in-memory sparse
#' track and supports all track-based summarizer functions. Intervals must not overlap.
#' @param func function name (see above)
#' @param params function parameters (see above)
#' @param ... additional PWM parameters
#' @inheritParams gvtrack.iterator
#' @inheritParams gvtrack.filter
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
#' # Masked base pair counting
#' gvtrack.create("masked_count", NULL, "masked.count")
#' gvtrack.create("masked_frac", NULL, "masked.frac")
#' gextract(c("masked_count", "masked_frac"), gintervals(1, 0, 10000), iterator = 1000)
#'
#' # Combined with GC content (unmasked regions only)
#' gvtrack.create("gc", NULL, "kmer.frac", kmer = "G")
#' gextract("gc * (1 - masked_frac)",
#'     gintervals(1, 0, 10000),
#'     iterator = 1000,
#'     colnames = "gc_unmasked"
#' )
#'
#' # Value-based track examples
#' # Create a data frame with intervals and numeric values
#' intervals_with_values <- data.frame(
#'     chrom = "chr1",
#'     start = c(100, 300, 500),
#'     end = c(200, 400, 600),
#'     score = c(10, 20, 30)
#' )
#' # Use as value-based sparse track (functions like sparse track)
#' gvtrack.create("value_track", intervals_with_values, "avg")
#' gvtrack.create("value_track_max", intervals_with_values, "max")
#' gextract(c("value_track", "value_track_max"),
#'     gintervals(1, 0, 10000),
#'     iterator = 1000
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
#' @seealso \code{\link{gvtrack.iterator}}, \code{\link{gvtrack.iterator.2d}}, \code{\link{gvtrack.filter}}
#' @export gvtrack.create
gvtrack.create <- function(vtrack = NULL, src = NULL, func = NULL, params = NULL, dim = NULL, sshift = NULL, eshift = NULL, filter = NULL, ...) {
    if (is.null(substitute(vtrack))) {
        stop("Usage: gvtrack.create(vtrack, src, func = NULL, params = NULL, dim = NULL, sshift = NULL, eshift = NULL, filter = NULL, ...)", call. = FALSE)
    }
    if (is.null(substitute(src)) && !(func %in% c("pwm", "pwm.max", "pwm.max.pos", "pwm.count", "kmer.count", "kmer.frac", "masked.count", "masked.frac"))) {
        stop("Usage: gvtrack.create(vtrack, src, func = NULL, params = NULL, dim = NULL, sshift = NULL, eshift = NULL, filter = NULL, ...)", call. = FALSE)
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
    } else if (!is.null(func) && func %in% c("masked.count", "masked.frac")) {
        # Masked counting has no parameters - just validate function name
        # Any additional parameters in ... will be ignored with a warning
        dots <- list(...)
        if (length(dots) > 0) {
            warning("masked.count and masked.frac functions do not accept parameters; ignoring extra arguments")
        }
        params <- list()
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

    # Apply iterator shifts if specified
    if (!is.null(dim) || !is.null(sshift) || !is.null(eshift)) {
        .set_vtrack_iterator_1d(
            vtrackstr,
            dim = dim,
            sshift = ifelse(is.null(sshift), 0, sshift),
            eshift = ifelse(is.null(eshift), 0, eshift)
        )
    }

    # Apply filter if specified
    if (!is.null(filter)) {
        gvtrack.filter(vtrackstr, filter)
    }

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

    .set_vtrack_iterator_1d(vtrackstr, dim = dim, sshift = sshift, eshift = eshift)
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
