#' Find PWM motif hits in DNA sequences
#'
#' Returns every window of each sequence whose PWM score reaches
#' \code{score.thresh}: the windows \code{gseq.pwm(mode = "count")} counts,
#' with their positions, strand and score.
#'
#' Windows are scored exactly as \code{\link{gseq.pwm}} scores them in
#' \code{"count"} mode with its default \code{neutral_chars} and no spatial
#' weights: a natural-log likelihood (including the prior), with a window
#' holding \code{N} or any other non-ACGT character left out. Gap characters
#' are not skipped - they are non-ACGT characters like any other - so on a
#' sequence without gaps the number of hits per sequence equals
#' \code{gseq.pwm(mode = "count")}'s count. Windows that cannot reach the
#' threshold are skipped without being scored, which makes a strict threshold
#' fast; the hits are the same as scoring every window.
#'
#' Strand: under \code{bidirect = TRUE} a window is one hit, whose
#' \code{score} is the log-sum-exp of its forward and reverse-complement
#' scores (the value thresholded, as in \code{gseq.pwm} and \code{pwm.count}),
#' and whose \code{strand} is the strand that scores higher on its own
#' (\code{1} on a tie), as \code{gseq.pwm(mode = "pos")} names it. Under
#' \code{bidirect = FALSE} only the strand named by \code{strand} is read and
#' \code{strand} is that strand. \code{start} and \code{end} are the window's
#' positions on the input sequence whatever the strand.
#'
#' @inheritParams gseq.pwm
#' @param score.thresh single number; windows scoring at or above it are hits.
#' @param extend logical or integer; how far past \code{end_pos} a window may
#'   run, as in \code{\link{gseq.pwm}}.
#'
#' @return A data frame with one row per hit and columns \code{seq_idx}
#'   (1-based index into \code{seqs}), \code{start} and \code{end} (1-based,
#'   inclusive window positions on the sequence), \code{strand} (\code{1} or
#'   \code{-1}) and \code{score}, ordered by \code{seq_idx} and \code{start}.
#'   It has no rows when nothing passes.
#'
#' @seealso \code{\link{gseq.pwm}}, \code{\link{gintervals.pwm_hits}} for
#'   hits over genomic intervals.
#'
#' @examples
#' pssm <- matrix(
#'     c(
#'         0.7, 0.1, 0.1, 0.1,
#'         0.1, 0.7, 0.1, 0.1,
#'         0.1, 0.1, 0.7, 0.1,
#'         0.1, 0.1, 0.1, 0.7
#'     ),
#'     ncol = 4, byrow = TRUE, dimnames = list(NULL, c("A", "C", "G", "T"))
#' )
#' seqs <- c("TTACGTAAACGTTT", "GGGGGGGG")
#' gseq.pwm_hits(seqs, pssm, score.thresh = -3)
#'
#' # the same windows gseq.pwm(mode = "count") counts
#' gseq.pwm(seqs, pssm, mode = "count", score.thresh = -3)
#'
#' @export
gseq.pwm_hits <- function(seqs,
                          pssm,
                          score.thresh,
                          bidirect = TRUE,
                          strand = 1L,
                          start_pos = NULL,
                          end_pos = NULL,
                          extend = FALSE,
                          prior = 0.01) {
    if (missing(score.thresh)) {
        stop("gseq.pwm_hits() requires a 'score.thresh' argument. PWM scores are log-likelihoods, so there is no default that suits every PSSM - pick a threshold from the score distribution of your own matrix, e.g. with gseq.pwm(mode = \"max\").", call. = FALSE)
    }
    score.thresh <- .coerce_score_thresh(score.thresh)
    pssm <- .coerce_pssm_matrix(
        pssm,
        numeric_msg = "pssm must be a numeric matrix or data frame with numeric columns",
        ncol_msg = "pssm must have columns named A, C, G, T",
        colnames_msg = "pssm must have columns named A, C, G, T"
    )
    if (!is.logical(bidirect) || length(bidirect) != 1 || is.na(bidirect)) {
        stop("bidirect must be TRUE or FALSE", call. = FALSE)
    }
    strand <- as.integer(strand)
    if (length(strand) != 1 || is.na(strand) || !strand %in% c(-1L, 1L)) {
        if (!bidirect) {
            stop("gseq.pwm_hits(bidirect = FALSE) needs an explicit strand: 1 for the forward strand, -1 for the reverse.", call. = FALSE)
        }
        strand <- 1L
    }
    if (!is.numeric(prior) || length(prior) != 1 || is.na(prior) || prior < 0 || prior > 1) {
        stop("prior must be a number between 0 and 1", call. = FALSE)
    }

    seqs <- as.character(seqs)
    b <- .normalize_bounds(seqs, start_pos, end_pos, extend)

    out <- .Call(
        "C_gseq_pwm_hits",
        seqs,
        pssm,
        as.numeric(score.thresh),
        as.logical(bidirect),
        strand,
        as.integer(b$roi_start),
        as.integer(b$roi_end),
        extend,
        as.numeric(prior)
    )
    as.data.frame(out, stringsAsFactors = FALSE)
}

#' Find PWM motif hits in genomic intervals
#'
#' Returns, for each interval, the motif windows that reach
#' \code{score.thresh}: exactly the windows a \code{pwm.count} virtual track
#' with the same parameters counts over that interval, as intervals with their
#' strand and score. It replaces running \code{\link{gscreen}} over a 1 bp
#' iterator.
#'
#' A window is anchored at its first (leftmost) base, which must lie in the
#' interval; with \code{extend = TRUE} it may run up to w - 1 bases past the
#' interval's end (w being the motif length), as the PWM virtual tracks do.
#' Scores are natural-log likelihoods including the prior, computed as the
#' \code{pwm} virtual tracks compute them (in single precision), and a window
#' holding \code{N} or any other non-ACGT base is not scored. Windows that
#' cannot reach the threshold are skipped without being scored, which makes a
#' strict threshold fast; the hits are the same as scoring every window.
#'
#' Strand: under \code{bidirect = TRUE} a window is one hit, whose
#' \code{score} is the log-sum-exp of its forward and reverse-complement
#' scores (the value \code{pwm.count} thresholds), and whose \code{strand} is
#' the strand that scores higher on its own (\code{1} on a tie), as
#' \code{pwm.max.pos} names it. Under \code{bidirect = FALSE} only the strand
#' named by \code{strand} is read and \code{strand} is that strand. A
#' \code{strand} column of \code{intervals} is ignored, as by the virtual
#' tracks. \code{start} and \code{end} are genomic coordinates of the window
#' whatever the strand.
#'
#' @param intervals genomic intervals to scan
#' @param pssm numeric matrix or data frame with columns A, C, G, T; each row
#'   holds the probabilities or counts of the four bases at one motif position.
#' @param score.thresh single number; windows scoring at or above it are hits.
#' @param bidirect logical; if TRUE (default) both strands are scanned.
#' @param strand \code{1} or \code{-1}: the strand scanned when
#'   \code{bidirect = FALSE}.
#' @param prior numeric; pseudocount added to the frequencies (default 0.01).
#' @param extend logical; if TRUE (default) a window anchored near an
#'   interval's end may read bases past it.
#'
#' @return An intervals data frame with columns \code{chrom}, \code{start},
#'   \code{end} (the window), \code{strand}, \code{score} and
#'   \code{intervalID} (the row of \code{intervals} the window is anchored
#'   in), sorted by \code{intervalID} and \code{start}; \code{NULL} when
#'   nothing passes. Overlapping input intervals report a window once per
#'   interval it is anchored in.
#'
#' @seealso \code{\link{gseq.pwm_hits}}, \code{\link{gvtrack.create}} for the
#'   \code{pwm.count} virtual track.
#'
#' @examples
#' \dontshow{
#' options(gmax.processes = 2)
#' }
#'
#' gdb.init_examples()
#' pssm <- matrix(
#'     c(
#'         0.7, 0.1, 0.1, 0.1,
#'         0.1, 0.7, 0.1, 0.1,
#'         0.1, 0.1, 0.7, 0.1,
#'         0.1, 0.1, 0.1, 0.7
#'     ),
#'     ncol = 4, byrow = TRUE, dimnames = list(NULL, c("A", "C", "G", "T"))
#' )
#' intervs <- gintervals(1, 0, 10000)
#' hits <- gintervals.pwm_hits(intervs, pssm, score.thresh = -3)
#' head(hits)
#'
#' # the windows pwm.count counts
#' gvtrack.create("motif_count", NULL, "pwm.count", pssm = pssm, score.thresh = -3)
#' gextract("motif_count", intervs, iterator = intervs)
#'
#' @export
gintervals.pwm_hits <- function(intervals = NULL,
                                pssm = NULL,
                                score.thresh = NULL,
                                bidirect = TRUE,
                                strand = 1L,
                                prior = 0.01,
                                extend = TRUE) {
    if (is.null(intervals) || is.null(pssm)) {
        stop("Usage: gintervals.pwm_hits(intervals, pssm, score.thresh, bidirect = TRUE, strand = 1, prior = 0.01, extend = TRUE)", call. = FALSE)
    }
    if (is.null(score.thresh)) {
        stop("gintervals.pwm_hits() requires a 'score.thresh' argument. PWM scores are log-likelihoods, so there is no default that suits every PSSM - pick a threshold from the score distribution of your own matrix, e.g. with a 'pwm.max' virtual track.", call. = FALSE)
    }
    if (!is.logical(extend) || length(extend) != 1 || is.na(extend)) {
        stop("extend must be TRUE or FALSE", call. = FALSE)
    }
    if (!is.logical(bidirect) || length(bidirect) != 1 || is.na(bidirect)) {
        stop("bidirect must be TRUE or FALSE", call. = FALSE)
    }
    .gcheckroot()

    # The parameters a pwm.count virtual track would get, validated the same way
    params <- .vtrack_params_pwm("pwm.count", NULL, list(
        pssm = pssm, bidirect = bidirect, prior = prior, extend = extend,
        strand = strand, score.thresh = score.thresh
    ))

    intervals <- rescue_ALLGENOME(intervals, substitute(intervals))
    res <- .gcall("C_gintervals_pwm_hits", intervals, params, .misha_env())
    if (is.null(res)) {
        return(NULL)
    }
    res <- res[order(res$intervalID, res$start), , drop = FALSE]
    rownames(res) <- NULL
    res
}
