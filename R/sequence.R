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
#' @param pssm numeric matrix with 4 columns named A, C, G, T
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
#'   scanning. Windows are w consecutive non-gap bases (default: FALSE)
#' @param gap_chars character vector; which characters count as gaps (default: c("-", "."))
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
                     skip_gaps = FALSE,
                     gap_chars = c("-", ".")) {
    # Validate inputs
    mode <- match.arg(mode)

    if (!is.matrix(pssm) || !is.numeric(pssm)) {
        stop("pssm must be a numeric matrix")
    }
    if (ncol(pssm) != 4L) {
        stop("pssm must have exactly 4 columns")
    }
    if (!setequal(colnames(pssm), c("A", "C", "G", "T"))) {
        stop("pssm columns must be named A, C, G, T")
    }

    # Reorder columns to A, C, G, T
    pssm <- pssm[, c("A", "C", "G", "T"), drop = FALSE]

    seqs <- as.character(seqs)
    seqs <- toupper(seqs)

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

    # Prepare spatial parameters as a list
    spat_params <- list(
        spat.factor = spat.factor,
        spat.bin = as.integer(spat.bin),
        spat.min = spat.min,
        spat.max = spat.max
    )

    # Call C++ implementation
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
        as.character(gap_chars)
    )

    # For mode="pos" with return_strand=TRUE, ensure it's a proper data.frame
    if (mode == "pos" && return_strand && is.list(out)) {
        out <- as.data.frame(out, stringsAsFactors = FALSE)
    }

    return(out)
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
#'   scanning. Windows are k consecutive non-gap bases (default: FALSE)
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
                      skip_gaps = FALSE,
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
    seqs <- toupper(seqs)

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
