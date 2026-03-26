#' Get reverse complement of DNA sequence
#'
#' Takes a DNA sequence string and returns its reverse complement.
#'
#' @param seq A character vector containing DNA sequences (using A,C,G,T). Ignores other characters and NA values.
#' @return A character vector of the same length as the input, containing the reverse
#'         complement sequences
#' @examples
#' grevcomp("ACTG") # Returns "CAGT"
#' grevcomp(c("ACTG", "GGCC")) # Returns c("CAGT", "GGCC")
#' grevcomp(c("ACTG", NA, "GGCC")) # Returns c("CAGT", NA, "GGCC")
#'
#' @export
grevcomp <- function(seq) {
    if (!is.character(seq)) {
        stop("Sequence must be a character string")
    }
    rev_s <- .Call("C_revcomp", seq)
    if (anyNA(seq)) {
        rev_s[is.na(seq)] <- NA_character_
    }

    if (!is.null(names(seq))) {
        names(rev_s) <- names(seq)
    }
    return(rev_s)
}

#' Get reverse complement of DNA sequence
#'
#' Alias for \code{\link{grevcomp}}. Takes a DNA sequence string and returns its reverse complement.
#'
#' @return A character vector of the same length as the input, containing the reverse
#'         complement sequences
#' @examples
#' gseq.revcomp("ACTG") # Returns "CAGT"
#' gseq.revcomp(c("ACTG", "GGCC")) # Returns c("CAGT", "GGCC")
#'
#' @inheritParams grevcomp
#' @export
#' @seealso \code{\link{grevcomp}}, \code{\link{gseq.rev}}, \code{\link{gseq.comp}}
gseq.revcomp <- grevcomp

#' Reverse DNA sequence
#'
#' Takes a DNA sequence string and returns its reverse (without complementing).
#'
#' @param seq A character vector containing DNA sequences. Preserves case and handles NA values.
#' @return A character vector of the same length as the input, containing the reversed sequences
#' @examples
#' gseq.rev("ACTG") # Returns "GTCA"
#' gseq.rev(c("ACTG", "GGCC")) # Returns c("GTCA", "CCGG")
#' gseq.rev(c("ACTG", NA, "GGCC")) # Returns c("GTCA", NA, "CCGG")
#'
#' @export
#' @seealso \code{\link{gseq.revcomp}}, \code{\link{gseq.comp}}
gseq.rev <- function(seq) {
    if (!is.character(seq)) {
        stop("Sequence must be a character string")
    }
    rev_s <- .Call("C_rev", seq)
    if (anyNA(seq)) {
        rev_s[is.na(seq)] <- NA_character_
    }

    if (!is.null(names(seq))) {
        names(rev_s) <- names(seq)
    }
    return(rev_s)
}

#' Complement DNA sequence
#'
#' Takes a DNA sequence string and returns its complement (without reversing).
#'
#' @param seq A character vector containing DNA sequences (using A,C,G,T). Preserves case and handles NA values.
#' @return A character vector of the same length as the input, containing the complemented sequences
#' @examples
#' gseq.comp("ACTG") # Returns "TGAC"
#' gseq.comp(c("ACTG", "GGCC")) # Returns c("TGAC", "CCGG")
#' gseq.comp(c("ACTG", NA, "GGCC")) # Returns c("TGAC", NA, "CCGG")
#'
#' @export
#' @seealso \code{\link{gseq.revcomp}}, \code{\link{gseq.rev}}
gseq.comp <- function(seq) {
    if (!is.character(seq)) {
        stop("Sequence must be a character string")
    }
    comp_s <- .Call("C_comp", seq)
    if (anyNA(seq)) {
        comp_s[is.na(seq)] <- NA_character_
    }

    if (!is.null(names(seq))) {
        names(comp_s) <- names(seq)
    }
    return(comp_s)
}

# Helper function to normalize ROI bounds for sequence scoring
.normalize_bounds <- function(seqs, start_pos, end_pos, w, extend) {
    n <- length(seqs)
    L <- nchar(seqs, type = "chars", allowNA = TRUE, keepNA = TRUE)

    # Handle start_pos and end_pos (can be NULL, scalar, or vector)
    roi_start <- if (is.null(start_pos)) rep.int(1L, n) else rep_len(as.integer(start_pos), n)
    roi_end <- if (is.null(end_pos)) L else rep_len(as.integer(end_pos), n)

    # Handle extend parameter
    E <- if (isFALSE(extend)) {
        0L
    } else if (isTRUE(extend)) {
        w - 1L
    } else {
        as.integer(extend)
    }

    if (any(E < 0L, na.rm = TRUE)) {
        stop("extend must be FALSE, TRUE, or a non-negative integer")
    }

    # Compute start_min and start_max (1-based, inclusive)
    start_min <- pmax.int(1L, roi_start - E)
    start_max <- pmin.int(pmax.int(0L, L - w + 1L), roi_end - w + 1L + E)

    list(
        L = L,
        roi_start = roi_start,
        roi_end = roi_end,
        E = E,
        start_min = start_min,
        start_max = start_max
    )
}

#' Score DNA sequences with a PWM over a region of interest
#'
#' Scores full DNA sequences using a Position Weight Matrix (PWM) over a specified
#' region of interest (ROI). The ROI is defined by \code{start_pos} and \code{end_pos}
#' (1-based, inclusive), with optional extension controlled by \code{extend}.
#' All reported positions are on the full input sequence.
#'
#' @param seqs character vector of DNA sequences (A/C/G/T/N; case-insensitive)
#' @param pssm numeric matrix or data frame with columns named A, C, G, T (additional columns are allowed and will be ignored)
#' @param mode character; one of "lse", "max", "pos", or "count"
#' @param bidirect logical; if TRUE, scans both strands (default: TRUE)
#' @param strand integer; 1=forward, -1=reverse, 0=both strands (default: 0)
#' @param score.thresh numeric; score threshold for \code{mode="count"} (default: 0)
#' @param start_pos integer or NULL; 1-based inclusive start of ROI (default: 1)
#' @param end_pos integer or NULL; 1-based inclusive end of ROI (default: sequence length)
#' @param extend logical or integer; extension of allowed window starts (default: FALSE)
#' @param spat.factor numeric vector; spatial weighting factors (optional)
#' @param spat.bin integer; bin size for spatial weighting
#' @param spat.min numeric; start of scanning window
#' @param spat.max numeric; end of scanning window
#' @param return_strand logical; if TRUE and \code{mode="pos"}, returns data.frame with
#'   \code{pos} and \code{strand} columns
#' @param skip_gaps logical; if TRUE, treat gap characters as holes and skip them while
#'   scanning. Windows are w consecutive non-gap bases (default: TRUE)
#' @param gap_chars character vector; which characters count as gaps (default: c("-", "."))
#' @param neutral_chars character vector; bases treated as unknown and scored with the average
#'   log probability per position (default: c("N", "n", "*"))
#' @param neutral_chars_policy character string; how to treat neutral characters. One of
#'   \code{"average"} (default; use the column's mean log-probability), \code{"log_quarter"}
#'   (always use \code{log(1/4)}), or \code{"na"} (return NA when a neutral character is
#'   encountered in the scanning window).
#' @param prior numeric; pseudocount added to frequencies (default: 0.01). Set to 0 for no pseudocounts.
#'
#' @return Numeric vector (for "lse"/"max"/"count" modes), integer vector (for "pos" mode),
#'   or data.frame with \code{pos} and \code{strand} columns (for "pos" mode with
#'   \code{return_strand=TRUE}). Returns NA when no valid windows exist.
#'
#' @details
#' This function scores DNA sequences directly without requiring a genomics database.
#' For detailed documentation on PWM scoring modes, parameters, and spatial weighting,
#' see \code{\link{gvtrack.create}} (functions "pwm", "pwm.max", "pwm.max.pos", "pwm.count").
#'
#' The ROI (region of interest) is defined by \code{start_pos} and \code{end_pos}.
#' The \code{extend} parameter controls whether motif matches can extend beyond the ROI boundaries.
#'
#' When \code{skip_gaps=TRUE}, characters specified in \code{gap_chars} are treated as gaps.
#' Windows are defined as w consecutive non-gap bases. All positions (\code{pos}) are reported
#' as 1-based indices on the original full sequence (including gaps). \code{start_pos} and
#' \code{end_pos} are interpreted as physical coordinates on the full sequence.
#'
#' Neutral characters (\code{neutral_chars}, default \code{c("N", "n", "*")}) are treated as
#' unknown bases in both orientations. Each neutral contributes the mean log-probability of the
#' corresponding PSSM column, yielding identical penalties on forward and reverse strands without
#' hard-coded background scores. In \code{mode = "max"} the reported value is the single best
#' strand score after applying any spatial weights; forward and reverse contributions are not
#' aggregated. This matches the default behavior of the PWM virtual tracks (\code{pwm.max},
#' \code{pwm.max.pos}, etc.).
#'
#' @seealso \code{\link{gvtrack.create}} for detailed PWM parameter documentation
#'
#' @examples
#' \dontrun{
#' # Create a PSSM (position-specific scoring matrix) with frequency values
#' pssm <- matrix(
#'     c(
#'         0.7, 0.1, 0.1, 0.1, # Position 1: mostly A
#'         0.1, 0.7, 0.1, 0.1, # Position 2: mostly C
#'         0.1, 0.1, 0.7, 0.1, # Position 3: mostly G
#'         0.1, 0.1, 0.1, 0.7 # Position 4: mostly T
#'     ),
#'     ncol = 4, byrow = TRUE
#' )
#' colnames(pssm) <- c("A", "C", "G", "T")
#'
#' # Example sequences
#' seqs <- c("ACGTACGTACGT", "GGGGACGTCCCC", "TTTTTTTTTTT")
#'
#' # Score sequences using log-sum-exp (default mode)
#' gseq.pwm(seqs, pssm, mode = "lse")
#'
#' # Get maximum score
#' gseq.pwm(seqs, pssm, mode = "max")
#'
#' # Find position of best match
#' gseq.pwm(seqs, pssm, mode = "pos")
#'
#' # Find position with strand information
#' gseq.pwm(seqs, pssm, mode = "pos", bidirect = TRUE, return_strand = TRUE)
#'
#' # Count matches above threshold
#' gseq.pwm(seqs, pssm, mode = "count", score.thresh = 0.5)
#'
#' # Score only a region of interest
#' gseq.pwm(seqs, pssm, mode = "max", start_pos = 3, end_pos = 10)
#'
#' # Allow matches to extend beyond ROI boundaries
#' gseq.pwm(seqs, pssm, mode = "count", start_pos = 5, end_pos = 8, extend = TRUE)
#'
#' # Spatial weighting example: higher weight in the center
#' spatial_weights <- c(0.5, 1.0, 2.0, 1.0, 0.5)
#' gseq.pwm(seqs, pssm,
#'     mode = "lse",
#'     spat.factor = spatial_weights,
#'     spat.bin = 2
#' )
#' }
#'
#' @export
gseq.pwm <- function(seqs,
                     pssm,
                     mode = c("lse", "max", "pos", "count"),
                     bidirect = TRUE,
                     strand = 0L,
                     score.thresh = 0,
                     start_pos = NULL,
                     end_pos = NULL,
                     extend = FALSE,
                     spat.factor = NULL,
                     spat.bin = 1L,
                     spat.min = NULL,
                     spat.max = NULL,
                     return_strand = FALSE,
                     skip_gaps = TRUE,
                     gap_chars = c("-", "."),
                     neutral_chars = c("N", "n", "*"),
                     neutral_chars_policy = c("average", "log_quarter", "na"),
                     prior = 0.01) {
    # Validate inputs
    mode <- match.arg(mode)

    pssm <- .coerce_pssm_matrix(
        pssm,
        numeric_msg = "pssm must be a numeric matrix or data frame with numeric columns",
        ncol_msg = "pssm must have columns named A, C, G, T",
        colnames_msg = "pssm must have columns named A, C, G, T"
    )

    seqs <- as.character(seqs)

    w <- nrow(pssm)

    # Validate strand
    strand <- as.integer(strand)
    if (!strand %in% c(-1L, 0L, 1L)) {
        stop("strand must be -1, 0, or 1")
    }

    # If bidirect is TRUE, override strand to 0
    if (bidirect) {
        strand <- 0L
    }

    # Validate gap parameters
    skip_gaps <- as.logical(skip_gaps)[1]

    neutral_chars <- as.character(neutral_chars)
    if (length(neutral_chars) > 0) {
        if (any(nchar(neutral_chars) != 1)) {
            stop("neutral_chars must contain only single characters")
        }
    }

    neutral_chars_policy <- match.arg(neutral_chars_policy)

    if (skip_gaps) {
        if (!is.character(gap_chars) || length(gap_chars) == 0) {
            stop("gap_chars must be a non-empty character vector")
        }
        if (any(nchar(gap_chars) != 1)) {
            stop("gap_chars must be single characters")
        }
        if (length(gap_chars) != length(unique(gap_chars))) {
            stop("gap_chars must be distinct")
        }
    }

    # Validate prior parameter
    if (!is.numeric(prior) || prior < 0 || prior > 1) {
        stop("prior must be a number between 0 and 1")
    }

    # Normalize bounds
    b <- .normalize_bounds(seqs, start_pos, end_pos, w, extend)

    # Prepare spatial parameters as a list
    spat_params <- list(
        spat.factor = spat.factor,
        spat.bin = as.integer(spat.bin),
        spat.min = spat.min,
        spat.max = spat.max
    )

    # Call C++ implementation (multitask or sequential)
    if (.ggetOption("gmultitasking")) {
        out <- .Call(
            "C_gseq_pwm_multitask",
            seqs,
            pssm,
            mode,
            as.logical(bidirect),
            strand,
            as.numeric(score.thresh),
            as.integer(b$roi_start),
            as.integer(b$roi_end),
            extend,
            spat_params,
            as.logical(return_strand),
            skip_gaps,
            as.character(gap_chars),
            as.numeric(prior),
            as.character(neutral_chars),
            neutral_chars_policy,
            .misha_env()
        )
    } else {
        out <- .Call(
            "C_gseq_pwm",
            seqs,
            pssm,
            mode,
            as.logical(bidirect),
            strand,
            as.numeric(score.thresh),
            as.integer(b$roi_start),
            as.integer(b$roi_end),
            extend,
            spat_params,
            as.logical(return_strand),
            skip_gaps,
            as.character(gap_chars),
            as.numeric(prior),
            as.character(neutral_chars),
            neutral_chars_policy
        )
    }

    # For mode="pos" with return_strand=TRUE, ensure it's a proper data.frame
    if (mode == "pos" && return_strand && is.list(out)) {
        out <- as.data.frame(out, stringsAsFactors = FALSE)
    }

    return(out)
}

#' Show optimal edits to reach a PWM score threshold
#'
#' For each input sequence (or genomic interval), finds the optimal motif window
#' and the specific base changes needed to reach \code{score.thresh}. Returns a
#' long-format data frame with one row per edit.
#'
#' This is an investigation tool: use it on a small set of positions (e.g., from
#' \code{gscreen}) to see what mutations would activate latent binding sites.
#'
#' @param seqs character vector of DNA sequences, OR a data frame of genomic
#'   intervals (with columns chrom, start, end). When intervals are provided,
#'   sequences are extracted automatically via \code{gseq.extract}.
#' @param pssm numeric matrix or data frame with columns A, C, G, T. Each row
#'   is a motif position.
#' @param score.thresh numeric; target PWM log-likelihood score to reach.
#' @param max_edits integer or NULL; maximum number of edits to search. NULL
#'   means no cap. Default NULL.
#' @param max_indels integer or NULL; maximum number of insertions and deletions
#'   allowed (default NULL, substitutions only). When > 0, a banded
#'   Needleman-Wunsch DP is used. Edits are reported with \code{edit_type}
#'   "sub", "ins", or "del".
#' @param bidirect logical; scan both strands? Default TRUE.
#' @param prior numeric; pseudocount for PSSM frequencies. Default 0.01.
#' @param score.min numeric or NULL; skip windows with PWM score below this.
#'   Default NULL (no filter).
#' @param score.max numeric or NULL; skip windows with PWM score above this.
#'   Default NULL (no filter).
#' @param extend logical or integer; extend sequence for boundary motifs.
#'   Default TRUE.
#' @param strand integer; which strand to scan when \code{bidirect=FALSE}.
#'   1=forward, -1=reverse. Default 1.
#' @return A data frame (long format) with one row per edit, containing columns:
#' \describe{
#'   \item{seq_idx}{Index into input sequences/intervals (1-based)}
#'   \item{strand}{+1 (forward) or -1 (reverse strand)}
#'   \item{window_start}{1-based position of optimal window within sequence}
#'   \item{score_before}{PWM score before edits}
#'   \item{score_after}{PWM score after all edits}
#'   \item{n_edits}{Total number of edits needed}
#'   \item{edit_num}{Which edit this row represents (1, 2, ...)}
#'   \item{motif_col}{1-based position within the motif where the edit occurs}
#'   \item{ref}{Current base at this position}
#'   \item{alt}{Suggested replacement base}
#'   \item{gain}{Score improvement from this individual edit}
#'   \item{window_seq}{Motif-length sequence at the optimal window (as seen by PSSM, reverse-complemented if on reverse strand)}
#'   \item{mutated_seq}{Same sequence with all edits applied}
#' }
#'
#' When intervals are provided, additional columns \code{chrom}, \code{start},
#' \code{end} are included.
#'
#' Sequences already above the threshold produce a single row with
#' \code{n_edits = 0}. Unreachable sequences are omitted from the output.
#'
#' @examples
#' \dontshow{
#' options(gmax.processes = 2)
#' }
#'
#' gdb.init_examples()
#'
#' # Simple PSSM
#' pssm <- matrix(c(1, 0, 0, 0, 0, 1, 0, 0),
#'     nrow = 2,
#'     dimnames = list(NULL, c("A", "C", "G", "T"))
#' )
#'
#' # What edits are needed?
#' gseq.pwm_edits("CCGTACGT", pssm, score.thresh = -0.5, prior = 0)
#'
#' @export
#' @seealso \code{\link{gseq.pwm}}, \code{\link{gvtrack.create}}
gseq.pwm_edits <- function(seqs,
                           pssm,
                           score.thresh,
                           max_edits = NULL,
                           max_indels = NULL,
                           bidirect = TRUE,
                           prior = 0.01,
                           score.min = NULL,
                           score.max = NULL,
                           extend = TRUE,
                           strand = 1L) {
    # Handle interval input: extract sequences
    is_intervals <- is.data.frame(seqs) &&
        all(c("chrom", "start", "end") %in% colnames(seqs))

    intervals_df <- NULL
    if (is_intervals) {
        intervals_df <- seqs
        # Extend END only — matching vtrack's calculate_expanded_interval.
        # The motif window starts within the interval and extends rightward,
        # so we need extra sequence after the interval end.
        # With this extension, the valid window starts are exactly
        # offsets 0..(interval_length-1) in the extracted sequence.
        # When indels are enabled, deletions produce windows of length
        # L + max_indels, so we need max_indels additional bases beyond
        # the standard extension (matching PWMEditDistanceScorer which
        # calls calculate_expanded_interval with motif_len + max_indels).
        w <- nrow(pssm)
        indel_extra <- if (!is.null(max_indels)) as.integer(max_indels) else 0L
        ext <- if (isTRUE(extend)) w - 1L + indel_extra else if (isFALSE(extend)) 0L else as.integer(extend) + indel_extra
        extended_intervals <- intervals_df
        extended_intervals$end <- extended_intervals$end + ext
        # Clamp to chromosome boundaries
        chrom_sizes <- gintervals.all()
        for (i in seq_len(nrow(extended_intervals))) {
            chr <- as.character(extended_intervals$chrom[i])
            chr_row <- chrom_sizes[chrom_sizes$chrom == chr, , drop = FALSE]
            if (nrow(chr_row) > 0) {
                extended_intervals$end[i] <- min(extended_intervals$end[i], chr_row$end[1])
            }
        }
        seqs <- gseq.extract(extended_intervals)
        # Without indels: the extended sequence has length
        # interval_length + w - 1, and C++ scans positions 0..(seqlen-w)
        # = 0..(interval_length-1), giving exactly interval_length windows.
        # No ROI constraint needed.
        #
        # With indels: we extend by an extra indel_extra bases so that
        # deletion windows (length up to w + max_indels) at the last
        # interval position have enough sequence. But without ROI, the
        # C++ would scan extra window starts beyond the interval.
        # We constrain via ROI so that valid window starts are exactly
        # 0..(interval_length-1), matching the vtrack behavior.
        # ROI is 1-based: roi_end = interval_length + w - 1 means
        # the last valid start is roi_end - w = interval_length - 1.
        if (indel_extra > 0L) {
            interval_lengths <- intervals_df$end - intervals_df$start
            roi_start <- rep(1L, length(seqs))
            roi_end <- as.integer(interval_lengths + w - 1L)
        } else {
            roi_start <- NULL
            roi_end <- NULL
        }
    } else {
        seqs <- as.character(seqs)
        roi_start <- NULL
        roi_end <- NULL
    }

    # Validate PSSM
    pssm <- .coerce_pssm_matrix(
        pssm,
        numeric_msg = "pssm must be a numeric matrix or data frame with numeric columns",
        ncol_msg = "pssm must have columns named A, C, G, T",
        colnames_msg = "pssm must have columns named A, C, G, T"
    )

    # Validate parameters
    if (!is.numeric(score.thresh) || length(score.thresh) != 1) {
        stop("score.thresh must be a single numeric value")
    }
    if (!is.null(max_edits)) {
        max_edits <- as.integer(max_edits)
        if (max_edits < 1) stop("max_edits must be NULL or a positive integer")
    }
    if (!is.null(max_indels)) {
        max_indels <- as.integer(max_indels)
        if (max_indels < 0) stop("max_indels must be NULL or a non-negative integer")
    }
    if (!is.logical(bidirect)) stop("bidirect must be TRUE or FALSE")
    if (!is.numeric(prior) || prior < 0 || prior > 1) stop("prior must be between 0 and 1")
    if (strand != 1 && strand != -1) stop("strand must be 1 or -1")

    strand_mode <- if (bidirect) 0L else as.integer(strand)

    # When intervals were used, extension was already done in R — pass FALSE to C++
    # to avoid double-extension of the ROI window search range.
    # For bare sequences, pass the user's extend value.
    c_extend <- if (is_intervals) FALSE else extend

    # Call C++
    result <- .Call(
        "C_gseq_pwm_edits",
        seqs,
        pssm,
        as.numeric(score.thresh),
        max_edits,
        max_indels,
        as.logical(bidirect),
        strand_mode,
        as.numeric(prior),
        if (!is.null(roi_start)) as.integer(roi_start) else NULL,
        if (!is.null(roi_end)) as.integer(roi_end) else NULL,
        if (is.logical(c_extend)) as.logical(c_extend) else as.integer(c_extend),
        if (!is.null(score.min)) as.numeric(score.min) else NULL,
        if (!is.null(score.max)) as.numeric(score.max) else NULL
    )

    # Add interval columns if input was intervals
    if (is_intervals && nrow(result) > 0) {
        result$chrom <- intervals_df$chrom[result$seq_idx]
        result$start <- intervals_df$start[result$seq_idx]
        result$end <- intervals_df$end[result$seq_idx]
        # Reorder columns
        result <- result[, c(
            "seq_idx", "chrom", "start", "end",
            "strand", "window_start",
            "score_before", "score_after", "n_edits",
            "edit_num", "motif_col", "ref", "alt", "gain",
            "edit_type", "window_seq", "mutated_seq"
        )]
    }

    return(result)
}

#' Score DNA sequences with a k-mer over a region of interest
#'
#' Counts exact matches of a k-mer in DNA sequences over a specified region of interest
#' (ROI). The ROI is defined by \code{start_pos} and \code{end_pos} (1-based, inclusive),
#' with optional extension controlled by \code{extend}.
#'
#' @param seqs character vector of DNA sequences (A/C/G/T/N; case-insensitive)
#' @param kmer single character string containing the k-mer to search for (A/C/G/T only)
#' @param mode character; one of "count" or "frac"
#' @param strand integer; 1=forward, -1=reverse, 0=both strands (default: 0)
#' @param start_pos integer or NULL; 1-based inclusive start of ROI (default: 1)
#' @param end_pos integer or NULL; 1-based inclusive end of ROI (default: sequence length)
#' @param extend logical or integer; extension of allowed window starts (default: FALSE)
#' @param skip_gaps logical; if TRUE, treat gap characters as holes and skip them while
#'   scanning. Windows are k consecutive non-gap bases (default: TRUE)
#' @param gap_chars character vector; which characters count as gaps (default: c("-", "."))
#'
#' @return Numeric vector with counts (for "count" mode) or fractions (for "frac" mode).
#'   Returns 0 when sequence is too short or ROI is invalid.
#'
#' @details
#' This function counts k-mer occurrences in DNA sequences directly without requiring
#' a genomics database. For detailed documentation on k-mer counting parameters, see
#' \code{\link{gvtrack.create}} (functions "kmer.count" and "kmer.frac").
#'
#' The ROI (region of interest) is defined by \code{start_pos} and \code{end_pos}.
#' The \code{extend} parameter controls whether k-mer matches can extend beyond the ROI boundaries.
#' For palindromic k-mers, use \code{strand=1} or \code{-1} to avoid double counting.
#'
#' When \code{skip_gaps=TRUE}, characters specified in \code{gap_chars} are treated as gaps.
#' Windows are defined as k consecutive non-gap bases. The \code{frac} denominator counts the
#' number of possible logical starts (non-gap windows) in the region. \code{start_pos} and
#' \code{end_pos} are interpreted as physical coordinates on the full sequence.
#'
#' @seealso \code{\link{gvtrack.create}} for detailed k-mer parameter documentation
#'
#' @examples
#' \dontrun{
#' # Example sequences
#' seqs <- c("CGCGCGCGCG", "ATATATATAT", "ACGTACGTACGT")
#'
#' # Count CG dinucleotides on both strands
#' gseq.kmer(seqs, "CG", mode = "count", strand = 0)
#'
#' # Count on forward strand only
#' gseq.kmer(seqs, "CG", mode = "count", strand = 1)
#'
#' # Get CG fraction
#' gseq.kmer(seqs, "CG", mode = "frac", strand = 0)
#'
#' # Count in a specific region
#' gseq.kmer(seqs, "CG", mode = "count", start_pos = 2, end_pos = 8)
#'
#' # Allow k-mer to extend beyond ROI boundaries
#' gseq.kmer(seqs, "CG", mode = "count", start_pos = 2, end_pos = 8, extend = TRUE)
#'
#' # Calculate GC content by summing G and C fractions
#' g_frac <- gseq.kmer(seqs, "G", mode = "frac", strand = 1)
#' c_frac <- gseq.kmer(seqs, "C", mode = "frac", strand = 1)
#' gc_content <- g_frac + c_frac
#' gc_content
#'
#' # Compare AT counts on different strands
#' at_forward <- gseq.kmer(seqs, "AT", mode = "count", strand = 1)
#' at_reverse <- gseq.kmer(seqs, "AT", mode = "count", strand = -1)
#' at_both <- gseq.kmer(seqs, "AT", mode = "count", strand = 0)
#' data.frame(forward = at_forward, reverse = at_reverse, both = at_both)
#' }
#'
#' @export
gseq.kmer <- function(seqs,
                      kmer,
                      mode = c("count", "frac"),
                      strand = 0L,
                      start_pos = NULL,
                      end_pos = NULL,
                      extend = FALSE,
                      skip_gaps = TRUE,
                      gap_chars = c("-", ".")) {
    # Validate inputs
    mode <- match.arg(mode)

    if (!is.character(kmer) || length(kmer) != 1) {
        stop("kmer must be a single character string")
    }

    kmer <- toupper(as.character(kmer)[1L])
    if (!grepl("^[ACGT]+$", kmer)) {
        stop("kmer must contain only A, C, G, T characters")
    }

    seqs <- as.character(seqs)

    w <- nchar(kmer, type = "chars")

    # Validate strand
    strand <- as.integer(strand)
    if (!strand %in% c(-1L, 0L, 1L)) {
        stop("strand must be -1, 0, or 1")
    }

    # Validate gap parameters
    skip_gaps <- as.logical(skip_gaps)[1]
    if (skip_gaps) {
        if (!is.character(gap_chars) || length(gap_chars) == 0) {
            stop("gap_chars must be a non-empty character vector")
        }
        if (any(nchar(gap_chars) != 1)) {
            stop("gap_chars must be single characters")
        }
        if (length(gap_chars) != length(unique(gap_chars))) {
            stop("gap_chars must be distinct")
        }
    }

    # Normalize bounds
    b <- .normalize_bounds(seqs, start_pos, end_pos, w, extend)

    # Call C++ implementation (returns counts or fractions depending on mode)
    result <- .Call(
        "C_gseq_kmer",
        seqs,
        kmer,
        mode,
        strand,
        as.integer(b$roi_start),
        as.integer(b$roi_end),
        extend,
        skip_gaps,
        as.character(gap_chars)
    )

    return(result)
}

#' Compute k-mer distribution in genomic intervals
#'
#' Counts the occurrence of all k-mers (of size k) within the specified
#' genomic intervals, optionally excluding masked regions.
#'
#' @param intervals Genomic intervals to analyze
#' @param k Integer k-mer size (1-10). Default is 6.
#' @param mask Optional intervals to exclude from counting. Positions within
#'        the mask will not contribute to k-mer counts.
#'
#' @return A data frame with columns:
#'   \describe{
#'     \item{kmer}{Character string representing the k-mer sequence}
#'     \item{count}{Number of occurrences of this k-mer}
#'   }
#'   Only k-mers with count > 0 are included. K-mers containing N bases
#'   are not counted.
#'
#' @examples
#' gdb.init_examples()
#'
#' # Count all 6-mers in first 10kb of chr1
#' intervals <- data.frame(chrom = "chr1", start = 0, end = 10000)
#' kmer_dist <- gseq.kmer.dist(intervals, k = 6)
#' head(kmer_dist)
#'
#' # Count dinucleotides
#' dinucs <- gseq.kmer.dist(intervals, k = 2)
#' dinucs
#'
#' # Count with mask
#' mask <- data.frame(chrom = "chr1", start = 5000, end = 6000)
#' kmer_dist_masked <- gseq.kmer.dist(intervals, k = 6, mask = mask)
#'
#' @seealso \code{\link{gseq.extract}}, \code{\link{gseq.kmer}}
#' @export
gseq.kmer.dist <- function(intervals, k = 6L, mask = NULL) {
    .gcheckroot()

    # Validate k
    k <- as.integer(k)
    if (length(k) != 1 || is.na(k) || k < 1 || k > 10) {
        stop("k must be an integer between 1 and 10", call. = FALSE)
    }

    # Call C++ implementation
    result <- .gcall(
        "C_gseq_kmer_dist",
        intervals,
        k,
        mask,
        .misha_env()
    )

    result
}
