#' Shuffle sequences keeping their k-mer counts
#'
#' Draws random sequences with exactly the same k-mer counts as the input
#' sequences - the dinucleotide counts with the default \code{k = 2} - for use
#' as a composition-matched background, e.g. in motif enrichment.
#'
#' Each shuffle is drawn uniformly at random from all the sequences that keep
#' the k-mer counts and the first k-1 bases. The construction is the
#' Eulerian-path one of Altschul & Erickson (1985) and Kandel et al. (1996), as
#' in uShuffle (Jiang et al. 2008): the (k-1)-mers are the vertices of a graph
#' whose edges are the k-mers of the sequence, a uniformly random spanning tree
#' of "last exits" is drawn with Wilson's algorithm, and the remaining edges out
#' of each vertex are put in a random order. \code{k = 1} is a plain random
#' permutation of the bases. Because the k-mers form a path through the
#' (k-1)-mers, the last k-1 bases are kept too.
#'
#' Only A, C, G and T (in either case) are shuffled. Any other character - N,
#' an IUPAC code, a gap - stays at its position and splits the sequence into
#' segments, each shuffled on its own; the k-mers that touch such a character
#' are therefore unchanged, and the k-mer counts of the whole sequence,
#' N-containing k-mers included, are kept. Shuffling the N's into the sequence
#' instead would scatter assembly gaps across it. Bases are compared ignoring
#' case and every position keeps its case, so a soft-masked stretch stays
#' lowercase at the same coordinates (with new bases in it); \code{toupper()}
#' the input first to ignore case altogether. A segment of k bases or fewer has
#' only one arrangement and is returned as is.
#'
#' The random numbers come from misha's own generator (xoshiro256**), not from
#' R's: with a \code{seed} the result does not depend on, and does not change,
#' R's random number state, and is the same on every platform and the same as
#' PyMisha's \code{gseq_shuffle} for the same seed. Without a \code{seed} one is
#' drawn from R's generator, so \code{set.seed()} makes the call reproducible
#' (but the sequences then differ from PyMisha's, whose generator differs).
#' Each input sequence gets its own random stream, derived from the seed and its
#' position in \code{seqs}, so its shuffles do not depend on the sequences that
#' come after it.
#'
#' @param seqs a character vector of DNA sequences (\code{NA} gives \code{NA}),
#'   or an intervals data frame, whose sequences are read with
#'   \code{\link{gseq.extract}} (strand-aware).
#' @param k the k-mer length whose counts are kept, 1 to 32. Default 2
#'   (dinucleotides).
#' @param n the number of shuffles of each sequence.
#' @param seed a whole number (|seed| <= 2^53) seeding the shuffle, or
#'   \code{NULL} to draw one from R's random number generator.
#'
#' @return With \code{n = 1}, a character vector like \code{seqs} (names kept).
#'   With \code{n > 1}, a character matrix with a row per sequence and a column
#'   per shuffle.
#'
#' @references
#' Altschul SF, Erickson BW (1985). Significance of nucleotide sequence
#' alignments: a method for random sequence permutation that preserves
#' dinucleotide and codon usage. \emph{Mol Biol Evol} 2(6):526-538.
#'
#' Kandel D, Matias Y, Unger R, Winkler P (1996). Shuffling biological
#' sequences. \emph{Discrete Applied Mathematics} 71:171-185.
#'
#' Jiang M, Anderson J, Gillespie J, Mayne M (2008). uShuffle: a useful tool
#' for shuffling biological sequences while preserving the k-let counts.
#' \emph{BMC Bioinformatics} 9:192.
#'
#' @seealso \code{\link{gseq.kmer}}, \code{\link{gseq.kmer.dist}},
#'   \code{\link{gsynth.sample}} for sampling from a Markov model, which keeps
#'   the expected composition rather than the exact counts.
#'
#' @examples
#' seqs <- c(a = "ACGTTGCAACGTAGCTAGCTTTACG", b = "ccgNNacgtacgtTTGACG")
#' gseq.shuffle(seqs, seed = 1)
#'
#' # three shuffles of each sequence, keeping trinucleotide counts
#' gseq.shuffle(seqs, k = 3, n = 3, seed = 1)
#'
#' # the dinucleotide counts are unchanged
#' sh <- gseq.shuffle(seqs[["a"]], seed = 7)
#' identical(
#'     table(substring(sh, 1:24, 2:25)),
#'     table(substring(seqs[["a"]], 1:24, 2:25))
#' )
#'
#' @export
gseq.shuffle <- function(seqs, k = 2L, n = 1L, seed = NULL) {
    if (is.data.frame(seqs)) {
        seqs <- gseq.extract(seqs)
    }
    if (is.factor(seqs)) {
        seqs <- as.character(seqs)
    }
    if (is.logical(seqs) && all(is.na(seqs))) {
        seqs <- as.character(seqs)
    }
    if (!is.character(seqs)) {
        stop("seqs must be a character vector of sequences or an intervals data frame", call. = FALSE)
    }

    .whole <- function(x) is.numeric(x) && length(x) == 1L && !is.na(x) && is.finite(x) && x == round(x)
    if (!.whole(k) || k < 1 || k > 32) {
        stop("k must be a whole number between 1 and 32", call. = FALSE)
    }
    if (!.whole(n) || n < 1 || n > .Machine$integer.max) {
        stop("n must be a positive whole number", call. = FALSE)
    }
    if (is.null(seed)) {
        # 52 random bits from R's generator, so set.seed() reproduces the call
        seed <- floor(stats::runif(1) * 2^26) * 2^26 + floor(stats::runif(1) * 2^26)
    } else if (!.whole(seed) || abs(seed) > 2^53) {
        stop("seed must be NULL or a whole number between -2^53 and 2^53", call. = FALSE)
    }

    nms <- names(seqs)
    n <- as.integer(n)
    res <- .gcall("C_gseq_shuffle", unname(seqs), as.integer(k), n, as.double(seed))
    if (n == 1L) {
        names(res) <- nms
        return(res)
    }
    matrix(res,
        nrow = length(seqs), ncol = n, byrow = TRUE,
        dimnames = if (!is.null(nms)) list(nms, NULL)
    )
}
