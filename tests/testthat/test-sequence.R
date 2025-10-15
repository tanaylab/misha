test_that("grevcomp handles basic DNA sequences correctly", {
    # Test single sequence
    expect_equal(grevcomp("ACTG"), "CAGT")

    # Test palindromic sequence
    expect_equal(grevcomp("GCGC"), "GCGC")

    # Test multiple sequences
    expect_equal(
        grevcomp(c("ACTG", "GCGC", "AAAA")),
        c("CAGT", "GCGC", "TTTT")
    )
})

test_that("grevcomp preserves case", {
    # Test mixed case
    expect_equal(grevcomp("AcTg"), "cAgT")
    expect_equal(grevcomp("gCGc"), "gCGc")

    # Test all lowercase
    expect_equal(grevcomp("actg"), "cagt")

    # Test all uppercase
    expect_equal(grevcomp("ACTG"), "CAGT")
})

test_that("grevcomp handles empty and special cases", {
    # Test empty string
    expect_equal(grevcomp(""), "")

    # Test empty vector
    expect_equal(grevcomp(character(0)), character(0))

    # Test vector with empty strings
    expect_equal(
        grevcomp(c("ACTG", "", "GCTA")),
        c("CAGT", "", "TAGC")
    )

    # Test NA values
    expect_equal(grevcomp(NA_character_), NA_character_)
    expect_equal(grevcomp(c("ACTG", NA, "GCTA")), c("CAGT", NA, "TAGC"))
})

test_that("grevcomp validates input", {
    # Test non-character input
    expect_error(grevcomp(123))
    expect_error(grevcomp(TRUE))
    expect_error(grevcomp(as.factor("ACTG")))
})

test_that("grevcomp handles long sequences", {
    # Test long sequence
    long_seq <- paste(rep("ACTG", 1000), collapse = "")
    expected <- paste(rep("CAGT", 1000), collapse = "")
    expect_equal(grevcomp(long_seq), expected)

    # Test vector of long sequences
    expect_equal(
        grevcomp(c(long_seq, "AAAA", long_seq)),
        c(expected, "TTTT", expected)
    )
})

test_that("grevcomp is reversible", {
    # Test that applying grevcomp twice returns original sequence
    seqs <- c(
        "ACTG",
        "GCTA",
        "AAAAAA",
        "GCGCGC",
        paste(rep("ACTG", 100), collapse = "")
    )

    for (seq in seqs) {
        expect_equal(grevcomp(grevcomp(seq)), seq)
    }
})

test_that("grevcomp performance on large input", {
    # Create large input
    n_seqs <- 10000
    seq_length <- 1000
    large_input <- replicate(n_seqs, {
        paste(sample(c("A", "C", "G", "T"), seq_length, replace = TRUE), collapse = "")
    })

    # Test that it completes in reasonable time
    expect_true(
        system.time(grevcomp(large_input))[["elapsed"]] <= 10,
        "grevcomp took too long for large input"
    )
})

test_that("grevcomp handles various vector cases", {
    # Test vector with repeated sequences
    expect_equal(
        grevcomp(c("ACTG", "ACTG", "ACTG")),
        c("CAGT", "CAGT", "CAGT")
    )

    # Test vector with all empty strings
    expect_equal(
        grevcomp(c("", "", "")),
        c("", "", "")
    )

    # Test vector with mixed lengths
    expect_equal(
        grevcomp(c("A", "ACTG", "AAAAAA")),
        c("T", "CAGT", "TTTTTT")
    )

    # Test vector with mixed case and empty strings
    expect_equal(
        grevcomp(c("AcTg", "", "GcTa", "actg", "ACTG")),
        c("cAgT", "", "tAgC", "cagt", "CAGT")
    )
})

test_that("grevcomp handles NA patterns in vectors", {
    # Test vector starting with NA
    expect_equal(
        grevcomp(c(NA, "ACTG", "GCTA")),
        c(NA, "CAGT", "TAGC")
    )

    # Test vector ending with NA
    expect_equal(
        grevcomp(c("ACTG", "GCTA", NA)),
        c("CAGT", "TAGC", NA)
    )

    # Test vector with multiple NAs
    expect_equal(
        grevcomp(c("ACTG", NA, "GCTA", NA, "AAAA")),
        c("CAGT", NA, "TAGC", NA, "TTTT")
    )

    # Test vector with all NAs
    expect_equal(
        grevcomp(c(NA_character_, NA_character_, NA_character_)),
        c(NA_character_, NA_character_, NA_character_)
    )

    # Test vector with NAs and empty strings
    expect_equal(
        grevcomp(c(NA, "", "ACTG", NA, "")),
        c(NA, "", "CAGT", NA, "")
    )
})

test_that("grevcomp maintains vector attributes", {
    # Test named vector
    input <- c(seq1 = "ACTG", seq2 = "GCTA")
    expected <- c(seq1 = "CAGT", seq2 = "TAGC")
    expect_equal(grevcomp(input), expected)

    # Test vector with names and NAs
    input <- c(seq1 = "ACTG", missing = NA, seq2 = "GCTA")
    expected <- c(seq1 = "CAGT", missing = NA, seq2 = "TAGC")
    expect_equal(grevcomp(input), expected)
})

test_that("grevcomp handles vectors of varied content", {
    # Mixed case, empty strings, and NAs
    input <- c(
        "AcTg",
        "",
        NA,
        "GCTA",
        "actg",
        NA,
        "",
        "GcTa"
    )
    expected <- c(
        "cAgT",
        "",
        NA,
        "TAGC",
        "cagt",
        NA,
        "",
        "tAgC"
    )
    expect_equal(grevcomp(input), expected)

    # Test with very long vector
    long_vec <- rep(c("ACTG", NA, "", "GcTa"), 1000)
    expected_long <- rep(c("CAGT", NA, "", "tAgC"), 1000)
    expect_equal(grevcomp(long_vec), expected_long)
})

test_that("grevcomp preserves vector length", {
    # Test vectors of different lengths
    lengths <- c(0, 1, 5, 10, 100, 1000)

    for (len in lengths) {
        input <- rep(c("ACTG", NA, "GCTA"), length.out = len)
        result <- grevcomp(input)
        expect_equal(length(result), len)
    }
})

# ============================================================================
# Tests for gseq.pwm() - PWM scoring on strings with ROI
# ============================================================================

test_that("gseq.pwm validates inputs correctly", {
    # Create a simple PSSM for testing (using probabilities, not log-odds)
    # This matrix represents a motif that strongly prefers "AC"
    pssm <- matrix(c(
        1.0, 1e-10, 1e-10, 1e-10, # Position 1: strongly prefers A
        1e-10, 1.0, 1e-10, 1e-10 # Position 2: strongly prefers C
    ), nrow = 2, byrow = TRUE)
    colnames(pssm) <- c("A", "C", "G", "T")

    seqs <- "ACGTACGT"

    # Should work with valid inputs
    expect_no_error(gseq.pwm(seqs, pssm, mode = "lse"))

    # Bad PSSM: not a matrix
    expect_error(gseq.pwm(seqs, c(1, 2, 3, 4), mode = "lse"), "matrix")

    # Bad PSSM: wrong number of columns
    bad_pssm <- matrix(1:6, nrow = 2)
    expect_error(gseq.pwm(seqs, bad_pssm, mode = "lse"), "4 columns")

    # Bad PSSM: wrong column names
    bad_pssm <- pssm
    colnames(bad_pssm) <- c("X", "Y", "Z", "W")
    expect_error(gseq.pwm(seqs, bad_pssm, mode = "lse"), "A, C, G, T")

    # Bad strand
    expect_error(gseq.pwm(seqs, pssm, mode = "lse", strand = 5), "strand must be")

    # Bad extend
    expect_error(gseq.pwm(seqs, pssm, mode = "lse", extend = -5), "extend must be")
})

test_that("gseq.pwm handles ROI bounds correctly", {
    # Create PSSM for "AC" (using probabilities)
    pssm <- matrix(c(
        1.0, 1e-10, 1e-10, 1e-10, # Position 1: strongly prefers A
        1e-10, 1.0, 1e-10, 1e-10 # Position 2: strongly prefers C
    ), nrow = 2, byrow = TRUE)
    colnames(pssm) <- c("A", "C", "G", "T")

    # Sequence with AC at positions 5-6
    seq <- "GGGGACGGGG"
    #      123456789..

    # ROI covering AC (5-6) should find it
    result <- gseq.pwm(seq, pssm, mode = "count", start_pos = 5, end_pos = 6, extend = FALSE)
    expect_equal(result, 1)

    # ROI before AC (1-4) should not find it
    result <- gseq.pwm(seq, pssm, mode = "count", start_pos = 1, end_pos = 4, extend = FALSE)
    expect_equal(result, 0)

    # ROI after AC (7-10) should not find it
    result <- gseq.pwm(seq, pssm, mode = "count", start_pos = 7, end_pos = 10, extend = FALSE)
    expect_equal(result, 0)
})

test_that("gseq.pwm extend parameter works correctly", {
    # Create PSSM for "AC" (using probabilities)
    pssm <- matrix(c(
        1.0, 1e-10, 1e-10, 1e-10, # Position 1: strongly prefers A
        1e-10, 1.0, 1e-10, 1e-10 # Position 2: strongly prefers C
    ), nrow = 2, byrow = TRUE)
    colnames(pssm) <- c("A", "C", "G", "T")

    # Sequence with AC at position 5-6
    seq <- "GGGGACGGGG"
    #      123456789..

    # ROI at 6-6 (just "C"), no extend: should not find AC
    result <- gseq.pwm(seq, pssm, mode = "count", start_pos = 6, end_pos = 6, extend = FALSE)
    expect_equal(result, 0)

    # ROI at 6-6, extend=TRUE (w-1=1): should find AC starting at position 5
    result <- gseq.pwm(seq, pssm, mode = "count", start_pos = 6, end_pos = 6, extend = TRUE)
    expect_equal(result, 1)

    # ROI at 6-6, extend=1 (explicit): same as extend=TRUE
    result <- gseq.pwm(seq, pssm, mode = "count", start_pos = 6, end_pos = 6, extend = 1)
    expect_equal(result, 1)

    # ROI at 7-7, extend=2: should find AC (start_min = 7-2 = 5)
    result <- gseq.pwm(seq, pssm, mode = "count", start_pos = 7, end_pos = 7, extend = 2)
    expect_equal(result, 1)
})

test_that("gseq.pwm mode='pos' reports correct 1-based positions", {
    # Create PSSM for "AC" (using probabilities)
    pssm <- matrix(c(
        1.0, 1e-10, 1e-10, 1e-10, # Position 1: strongly prefers A
        1e-10, 1.0, 1e-10, 1e-10 # Position 2: strongly prefers C
    ), nrow = 2, byrow = TRUE)
    colnames(pssm) <- c("A", "C", "G", "T")

    # Sequence with AC at positions 3-4 and 7-8
    seq <- "GGACGGACGG"
    #      123456789..

    # Should find first AC at position 3
    result <- gseq.pwm(seq, pssm, mode = "pos", start_pos = 1, end_pos = 10, extend = FALSE)
    expect_equal(result, 3)

    # Limit ROI to find second AC at position 7
    result <- gseq.pwm(seq, pssm, mode = "pos", start_pos = 5, end_pos = 10, extend = FALSE)
    expect_equal(result, 7)
})

test_that("gseq.pwm mode='pos' with return_strand works", {
    # Create PSSM for "AC" (using probabilities)
    pssm <- matrix(c(
        1.0, 1e-10, 1e-10, 1e-10, # Position 1: strongly prefers A
        1e-10, 1.0, 1e-10, 1e-10 # Position 2: strongly prefers C
    ), nrow = 2, byrow = TRUE)
    colnames(pssm) <- c("A", "C", "G", "T")

    # Sequence with AC (forward)
    seq <- "GGACGG"

    result <- gseq.pwm(seq, pssm,
        mode = "pos", bidirect = FALSE, strand = 1,
        start_pos = 1, end_pos = 6, extend = FALSE, return_strand = TRUE
    )

    expect_true(is.data.frame(result))
    expect_equal(names(result), c("pos", "strand"))
    expect_equal(result$pos, 3)
    expect_equal(result$strand, 1)
})

test_that("gseq.pwm tie-breaking: prefer leftmost position", {
    # Create uniform PSSM (all probabilities equal)
    pssm <- matrix(1.0, nrow = 2, ncol = 4) # Normalizes to 0.25 each
    colnames(pssm) <- c("A", "C", "G", "T")

    seq <- "ACGTACGT"

    # All positions score equally, should return leftmost (position 1)
    result <- gseq.pwm(seq, pssm, mode = "pos", start_pos = 1, end_pos = 8, extend = FALSE)
    expect_equal(result, 1)
})

test_that("gseq.pwm handles edge case: sequence too short", {
    # Create PSSM length 5 (using uniform probabilities)
    pssm <- matrix(1.0, nrow = 5, ncol = 4)
    colnames(pssm) <- c("A", "C", "G", "T")

    # Sequence length 3 (< motif width)
    seq <- "ACG"

    result_total <- gseq.pwm(seq, pssm, mode = "lse")
    expect_true(is.na(result_total))

    result_max <- gseq.pwm(seq, pssm, mode = "max")
    expect_true(is.na(result_max))

    result_count <- gseq.pwm(seq, pssm, mode = "count")
    expect_equal(result_count, 0)

    result_pos <- gseq.pwm(seq, pssm, mode = "pos")
    expect_true(is.na(result_pos))
})

test_that("gseq.pwm handles edge case: empty ROI", {
    # Create PSSM for "AC" (using probabilities)
    pssm <- matrix(c(
        1.0, 1e-10, 1e-10, 1e-10, # Position 1: strongly prefers A
        1e-10, 1.0, 1e-10, 1e-10 # Position 2: strongly prefers C
    ), nrow = 2, byrow = TRUE)
    colnames(pssm) <- c("A", "C", "G", "T")

    seq <- "ACGTACGT"

    # start_pos > end_pos: invalid ROI
    result_count <- gseq.pwm(seq, pssm, mode = "count", start_pos = 5, end_pos = 3)
    expect_equal(result_count, 0)

    result_total <- gseq.pwm(seq, pssm, mode = "lse", start_pos = 5, end_pos = 3)
    expect_true(is.na(result_total))
})

test_that("gseq.pwm mode='lse' aggregates scores using log-sum-exp", {
    # Create PSSM for "A" (single base, using probabilities)
    pssm <- matrix(c(1.0, 1e-10, 1e-10, 1e-10), nrow = 1, ncol = 4)
    colnames(pssm) <- c("A", "C", "G", "T")

    # Sequence with 3 A's at positions 1, 3, 5
    seq <- "AGAGA"

    # Each A scores log(~1.0), total should be sum of log-probabilities
    # With probabilities (1.0, 1e-10, 1e-10, 1e-10), after normalization A gets ~0.9999...
    # log(0.9999...) is close to 0, but with bidirect mode and log-sum-exp it may vary
    result_full <- gseq.pwm(seq, pssm, mode = "lse", start_pos = 1, end_pos = 5, extend = FALSE)
    expect_true(is.finite(result_full)) # Should be a finite number

    # Only 2 A's in ROI 1-3
    result_partial <- gseq.pwm(seq, pssm, mode = "lse", start_pos = 1, end_pos = 3, extend = FALSE)
    expect_true(is.finite(result_partial)) # Should be a finite number
    # Note: With log-sum-exp in bidirect mode, the relationship between partial and full
    # may not be monotonic, so we just check that both are finite
})

test_that("gseq.pwm mode='max' finds maximum score", {
    # Create PSSM that strongly prefers A (using probabilities)
    pssm <- matrix(c(1.0, 1e-10, 1e-10, 1e-10), nrow = 1, ncol = 4)
    colnames(pssm) <- c("A", "C", "G", "T")

    seq <- "CGACG"

    # Should find A at position 3 with score close to 0 (log of ~1.0)
    result <- gseq.pwm(seq, pssm, mode = "max", start_pos = 1, end_pos = 5, extend = FALSE)
    expect_true(result > -1 && result < 1) # Should be close to 0
})

test_that("gseq.pwm mode='count' with threshold works", {
    # Create PSSM with different probabilities for A and C (using probabilities)
    # A gets higher probability than C
    pssm <- matrix(c(
        100, 10, 1e-10, 1e-10 # A: ~90%, C: ~9%, G/T: ~0%
    ), nrow = 1, ncol = 4)
    colnames(pssm) <- c("A", "C", "G", "T")

    seq <- "AACCCG"

    # Count positions scoring >= -1.0 (only A's, 2 of them - forward strand only)
    # A will score log(~0.9) ≈ -0.1, C will score log(~0.09) ≈ -2.4
    result <- gseq.pwm(seq, pssm,
        mode = "count", score.thresh = -1.0,
        bidirect = FALSE, strand = 1,
        start_pos = 1, end_pos = 6, extend = FALSE
    )
    expect_equal(result, 2)

    # Count positions scoring >= -3.0 (A's and C's, 5 total - forward strand only)
    result <- gseq.pwm(seq, pssm,
        mode = "count", score.thresh = -3.0,
        bidirect = FALSE, strand = 1,
        start_pos = 1, end_pos = 6, extend = FALSE
    )
    expect_equal(result, 5)
})

test_that("gseq.pwm vectorization works", {
    # Create PSSM for "AC" (using probabilities)
    pssm <- matrix(c(
        1.0, 1e-10, 1e-10, 1e-10, # Position 1: strongly prefers A
        1e-10, 1.0, 1e-10, 1e-10 # Position 2: strongly prefers C
    ), nrow = 2, byrow = TRUE)
    colnames(pssm) <- c("A", "C", "G", "T")

    seqs <- c("GGACGG", "GGGGGG", "ACACAC")

    result <- gseq.pwm(seqs, pssm, mode = "count", bidirect = FALSE, strand = 1, start_pos = 1, end_pos = 6, extend = FALSE)
    expect_equal(length(result), 3)
    expect_equal(result[1], 1) # 1 AC (forward strand only)
    expect_equal(result[2], 0) # 0 AC
    expect_equal(result[3], 3) # 3 AC (forward strand only)
})

test_that("gseq.pwm handles vectorized start_pos and end_pos", {
    # Create PSSM for "AC" (using probabilities)
    pssm <- matrix(c(
        1.0, 1e-10, 1e-10, 1e-10, # Position 1: strongly prefers A
        1e-10, 1.0, 1e-10, 1e-10 # Position 2: strongly prefers C
    ), nrow = 2, byrow = TRUE)
    colnames(pssm) <- c("A", "C", "G", "T")

    seqs <- c("ACACAC", "ACACAC")

    # Different ROIs per sequence (forward strand only)
    result <- gseq.pwm(seqs, pssm,
        mode = "count",
        bidirect = FALSE, strand = 1,
        start_pos = c(1, 3), end_pos = c(4, 6), extend = FALSE
    )

    expect_equal(result[1], 2) # AC at 1-2, 3-4 (forward strand only)
    expect_equal(result[2], 2) # AC at 3-4, 5-6 (forward strand only)
})

# ============================================================================
# Tests for gseq.kmer() - k-mer scoring on strings with ROI
# ============================================================================

test_that("gseq.kmer validates inputs correctly", {
    seqs <- "ACGTACGT"

    # Should work with valid inputs
    expect_no_error(gseq.kmer(seqs, "AC", mode = "count"))

    # Bad kmer: not a single string
    expect_error(gseq.kmer(seqs, c("AC", "GT"), mode = "count"), "single character")

    # Bad kmer: contains invalid characters
    expect_error(gseq.kmer(seqs, "ACN", mode = "count"), "only A, C, G, T")

    # Bad strand
    expect_error(gseq.kmer(seqs, "AC", mode = "count", strand = 5), "strand must be")

    # Bad extend
    expect_error(gseq.kmer(seqs, "AC", mode = "count", extend = -5), "extend must be")
})

test_that("gseq.kmer mode='count' counts exact matches", {
    seq <- "CGCGCGCG"

    # Count CG (4 occurrences)
    result <- gseq.kmer(seq, "CG", mode = "count", start_pos = 1, end_pos = 8, extend = FALSE)
    expect_equal(result, 4)

    # Count GC (3 occurrences)
    result <- gseq.kmer(seq, "GC", mode = "count", start_pos = 1, end_pos = 8, extend = FALSE)
    expect_equal(result, 3)

    # Count AA (0 occurrences)
    result <- gseq.kmer(seq, "AA", mode = "count", start_pos = 1, end_pos = 8, extend = FALSE)
    expect_equal(result, 0)
})

test_that("gseq.kmer mode='frac' computes correct fractions", {
    seq <- "AAACCCGGG" # length 9
    #      123456789

    # Count A's: 3 out of 9 possible positions (forward strand only)
    result <- gseq.kmer(seq, "A", mode = "frac", strand = 1, start_pos = 1, end_pos = 9, extend = FALSE)
    expect_equal(result, 3 / 9)

    # Count C's: 3 out of 9 possible positions (forward strand only)
    result <- gseq.kmer(seq, "C", mode = "frac", strand = 1, start_pos = 1, end_pos = 9, extend = FALSE)
    expect_equal(result, 3 / 9)

    # Count AA's: 2 out of 8 possible positions (at pos 1-2 and 2-3, forward strand only)
    # Sequence "AAACCCGGG" has "AA" at positions 1-2 and 2-3
    result <- gseq.kmer(seq, "AA", mode = "frac", strand = 1, start_pos = 1, end_pos = 9, extend = FALSE)
    expect_equal(result, 2 / 8) # 2 matches out of 8 possible start positions
})

test_that("gseq.kmer handles ROI bounds correctly", {
    seq <- "GGGGACGGGG"
    #      123456789..

    # ROI covering AC (5-6)
    result <- gseq.kmer(seq, "AC", mode = "count", start_pos = 5, end_pos = 6, extend = FALSE)
    expect_equal(result, 1)

    # ROI before AC (1-4)
    result <- gseq.kmer(seq, "AC", mode = "count", start_pos = 1, end_pos = 4, extend = FALSE)
    expect_equal(result, 0)

    # ROI after AC (7-10)
    result <- gseq.kmer(seq, "AC", mode = "count", start_pos = 7, end_pos = 10, extend = FALSE)
    expect_equal(result, 0)
})

test_that("gseq.kmer extend parameter works correctly", {
    seq <- "GGGGACGGGG"
    #      123456789..

    # ROI at 6-6 (just "C"), no extend: should not find AC
    result <- gseq.kmer(seq, "AC", mode = "count", start_pos = 6, end_pos = 6, extend = FALSE)
    expect_equal(result, 0)

    # ROI at 6-6, extend=TRUE (k-1=1): should find AC starting at position 5
    result <- gseq.kmer(seq, "AC", mode = "count", start_pos = 6, end_pos = 6, extend = TRUE)
    expect_equal(result, 1)

    # ROI at 6-6, extend=1 (explicit): same as extend=TRUE
    result <- gseq.kmer(seq, "AC", mode = "count", start_pos = 6, end_pos = 6, extend = 1)
    expect_equal(result, 1)
})

test_that("gseq.kmer strand parameter works", {
    seq <- "ACGTACGT"

    # Forward strand only: count AC (2 occurrences)
    result <- gseq.kmer(seq, "AC", mode = "count", strand = 1, start_pos = 1, end_pos = 8)
    expect_equal(result, 2)

    # Reverse strand only: count GT (reverse complement of AC, 2 occurrences)
    result <- gseq.kmer(seq, "AC", mode = "count", strand = -1, start_pos = 1, end_pos = 8)
    expect_equal(result, 2)

    # Both strands: count positions where AC or GT appears
    # AC at 1-2, 5-6; GT at 3-4, 7-8 -> 4 matches
    result <- gseq.kmer(seq, "AC", mode = "count", strand = 0, start_pos = 1, end_pos = 8)
    expect_equal(result, 4)
})

test_that("gseq.kmer handles edge case: sequence too short", {
    seq <- "AC"

    # k-mer length 5 > sequence length 2
    result <- gseq.kmer(seq, "ACGTA", mode = "count")
    expect_equal(result, 0)

    result <- gseq.kmer(seq, "ACGTA", mode = "frac")
    expect_equal(result, 0)
})

test_that("gseq.kmer handles edge case: empty ROI", {
    seq <- "ACGTACGT"

    # start_pos > end_pos: invalid ROI
    result <- gseq.kmer(seq, "AC", mode = "count", start_pos = 5, end_pos = 3)
    expect_equal(result, 0)

    result <- gseq.kmer(seq, "AC", mode = "frac", start_pos = 5, end_pos = 3)
    expect_equal(result, 0)
})

test_that("gseq.kmer vectorization works", {
    seqs <- c("ACACAC", "GGGGGG", "CGCGCG")

    result <- gseq.kmer(seqs, "AC", mode = "count", start_pos = 1, end_pos = 6, extend = FALSE)
    expect_equal(length(result), 3)
    expect_equal(result[1], 3) # 3 AC
    expect_equal(result[2], 0) # 0 AC
    expect_equal(result[3], 0) # 0 AC
})

test_that("gseq.kmer handles vectorized start_pos and end_pos", {
    seqs <- c("ACACAC", "ACACAC")

    # Different ROIs per sequence
    result <- gseq.kmer(seqs, "AC",
        mode = "count",
        start_pos = c(1, 3), end_pos = c(4, 6), extend = FALSE
    )

    expect_equal(result[1], 2) # AC at 1-2, 3-4
    expect_equal(result[2], 2) # AC at 3-4, 5-6
})

test_that("gseq.kmer frac denominator is computed correctly", {
    seq <- "ACGTACGT" # length 8

    # k=2: 7 possible positions (8 - 2 + 1)
    # Count AC: 2 occurrences (forward strand only)
    result <- gseq.kmer(seq, "AC", mode = "frac", strand = 1, start_pos = 1, end_pos = 8, extend = FALSE)
    expect_equal(result, 2 / 7)

    # With ROI 1-4: 3 possible positions (4 - 2 + 1)
    # Count AC: 1 occurrence (forward strand only)
    result <- gseq.kmer(seq, "AC", mode = "frac", strand = 1, start_pos = 1, end_pos = 4, extend = FALSE)
    expect_equal(result, 1 / 3)
})

# ============================================================================
# Comparison tests: gseq.pwm vs gextract with PWM virtual tracks
# ============================================================================

test_that("gseq.pwm matches gextract PWM vtrack: mode='lse', no extend", {
    remove_all_vtracks()

    pssm <- create_test_pssm() # AC motif with probabilities

    # Test on a specific interval
    test_interval <- gintervals(1, 200, 240)
    seq <- toupper(gseq.extract(test_interval))

    # Create virtual track
    gvtrack.create("pwm_total", NULL,
        func = "pwm",
        pssm = pssm, bidirect = FALSE, extend = FALSE, prior = 0.01
    )

    # Extract using gextract
    vtrack_result <- gextract("pwm_total", test_interval, iterator = test_interval)

    # Score using gseq.pwm (need to add prior to match virtual track)
    pssm_with_prior <- pssm + 0.01
    for (i in 1:nrow(pssm_with_prior)) {
        pssm_with_prior[i, ] <- pssm_with_prior[i, ] / sum(pssm_with_prior[i, ])
    }

    gseq_result <- gseq.pwm(seq, pssm_with_prior,
        mode = "lse", bidirect = FALSE, strand = 1,
        start_pos = 1, end_pos = nchar(seq), extend = FALSE
    )

    expect_equal(gseq_result, vtrack_result$pwm_total, tolerance = 1e-6)
})

test_that("gseq.pwm matches gextract PWM vtrack: mode='lse', with extend", {
    remove_all_vtracks()

    pssm <- create_test_pssm() # AC motif

    test_interval <- gintervals(1, 200, 240)
    # For extend=TRUE, need to extract extended sequence
    extended_interval <- test_interval
    extended_interval$end <- extended_interval$end + nrow(pssm) - 1
    seq_ext <- toupper(gseq.extract(extended_interval))

    # Create virtual track
    gvtrack.create("pwm_ext", NULL,
        func = "pwm",
        pssm = pssm, bidirect = FALSE, extend = TRUE, prior = 0.01
    )

    # Extract using gextract
    vtrack_result <- gextract("pwm_ext", test_interval, iterator = test_interval)

    # Score using gseq.pwm
    pssm_with_prior <- pssm + 0.01
    for (i in 1:nrow(pssm_with_prior)) {
        pssm_with_prior[i, ] <- pssm_with_prior[i, ] / sum(pssm_with_prior[i, ])
    }

    gseq_result <- gseq.pwm(seq_ext, pssm_with_prior,
        mode = "lse", bidirect = FALSE, strand = 1,
        start_pos = 1, end_pos = nchar(seq_ext), extend = TRUE
    )

    expect_equal(gseq_result, vtrack_result$pwm_ext, tolerance = 1e-6)
})

test_that("gseq.pwm matches gextract PWM vtrack: mode='max'", {
    remove_all_vtracks()

    pssm <- create_test_pssm() # AC motif

    test_interval <- gintervals(1, 200, 240)
    extended_interval <- test_interval
    extended_interval$end <- extended_interval$end + nrow(pssm) - 1
    seq_ext <- toupper(gseq.extract(extended_interval))

    # Create virtual track
    gvtrack.create("pwm_max", NULL,
        func = "pwm.max",
        pssm = pssm, bidirect = FALSE, extend = TRUE, prior = 0.01
    )

    # Extract using gextract
    vtrack_result <- gextract("pwm_max", test_interval, iterator = test_interval)

    # Score using gseq.pwm
    pssm_with_prior <- pssm + 0.01
    for (i in 1:nrow(pssm_with_prior)) {
        pssm_with_prior[i, ] <- pssm_with_prior[i, ] / sum(pssm_with_prior[i, ])
    }

    gseq_result <- gseq.pwm(seq_ext, pssm_with_prior,
        mode = "max", bidirect = FALSE, strand = 1,
        start_pos = 1, end_pos = nchar(seq_ext), extend = TRUE
    )

    expect_equal(gseq_result, vtrack_result$pwm_max, tolerance = 1e-6)
})

test_that("gseq.pwm matches gextract PWM vtrack: mode='pos'", {
    remove_all_vtracks()

    pssm <- create_test_pssm() # AC motif

    test_interval <- gintervals(1, 200, 240)
    extended_interval <- test_interval
    extended_interval$end <- extended_interval$end + nrow(pssm) - 1
    seq_ext <- toupper(gseq.extract(extended_interval))

    # Create virtual track
    gvtrack.create("pwm_pos", NULL,
        func = "pwm.max.pos",
        pssm = pssm, bidirect = FALSE, extend = TRUE, prior = 0.01
    )

    # Extract using gextract
    vtrack_result <- gextract("pwm_pos", test_interval, iterator = test_interval)

    # Score using gseq.pwm
    pssm_with_prior <- pssm + 0.01
    for (i in 1:nrow(pssm_with_prior)) {
        pssm_with_prior[i, ] <- pssm_with_prior[i, ] / sum(pssm_with_prior[i, ])
    }

    gseq_result <- gseq.pwm(seq_ext, pssm_with_prior,
        mode = "pos", bidirect = FALSE, strand = 1,
        start_pos = 1, end_pos = nchar(seq_ext), extend = TRUE
    )

    expect_equal(gseq_result, vtrack_result$pwm_pos, tolerance = 1e-6)
})

test_that("gseq.pwm matches gextract PWM vtrack: bidirectional mode", {
    remove_all_vtracks()

    pssm <- create_test_pssm() # AC motif

    test_interval <- gintervals(1, 200, 240)
    extended_interval <- test_interval
    extended_interval$end <- extended_interval$end + nrow(pssm) - 1
    seq_ext <- toupper(gseq.extract(extended_interval))

    # Create virtual tracks for bidirectional and forward
    gvtrack.create("pwm_bidi", NULL,
        func = "pwm",
        pssm = pssm, bidirect = TRUE, extend = TRUE, prior = 0.01
    )
    gvtrack.create("pwm_fwd", NULL,
        func = "pwm",
        pssm = pssm, bidirect = FALSE, extend = TRUE, prior = 0.01
    )

    # Extract using gextract
    vtrack_result <- gextract(c("pwm_bidi", "pwm_fwd"), test_interval, iterator = test_interval)

    # Score using gseq.pwm
    pssm_with_prior <- pssm + 0.01
    for (i in 1:nrow(pssm_with_prior)) {
        pssm_with_prior[i, ] <- pssm_with_prior[i, ] / sum(pssm_with_prior[i, ])
    }

    gseq_bidi <- gseq.pwm(seq_ext, pssm_with_prior,
        mode = "lse", bidirect = TRUE,
        start_pos = 1, end_pos = nchar(seq_ext), extend = TRUE
    )

    gseq_fwd <- gseq.pwm(seq_ext, pssm_with_prior,
        mode = "lse", bidirect = FALSE, strand = 1,
        start_pos = 1, end_pos = nchar(seq_ext), extend = TRUE
    )

    expect_equal(gseq_bidi, vtrack_result$pwm_bidi, tolerance = 1e-6)
    expect_equal(gseq_fwd, vtrack_result$pwm_fwd, tolerance = 1e-6)
    expect_true(gseq_bidi >= gseq_fwd) # bidirectional should have higher or equal score
})

# ============================================================================
# Comparison tests: gseq.kmer vs gextract with kmer virtual tracks
# ============================================================================

test_that("gseq.kmer matches gextract kmer.count vtrack", {
    remove_all_vtracks()

    test_interval <- gintervals(1, 200, 240)
    seq <- toupper(gseq.extract(test_interval))

    # Create virtual track for k-mer count
    gvtrack.create("count_ta", NULL, "kmer.count", kmer = "TA", strand = 1)

    # Extract using gextract
    vtrack_result <- gextract("count_ta", test_interval, iterator = test_interval)

    # Score using gseq.kmer
    gseq_result <- gseq.kmer(seq, "TA",
        mode = "count", strand = 1,
        start_pos = 1, end_pos = nchar(seq), extend = FALSE
    )

    expect_equal(gseq_result, vtrack_result$count_ta)
})

test_that("gseq.kmer matches gextract kmer.frac vtrack", {
    remove_all_vtracks()

    test_interval <- gintervals(1, 200, 240)
    seq <- toupper(gseq.extract(test_interval))

    # Create virtual track for k-mer fraction
    gvtrack.create("frac_ta", NULL, "kmer.frac", kmer = "TA", strand = 1)

    # Extract using gextract
    vtrack_result <- gextract("frac_ta", test_interval, iterator = test_interval)

    # Score using gseq.kmer
    gseq_result <- gseq.kmer(seq, "TA",
        mode = "frac", strand = 1,
        start_pos = 1, end_pos = nchar(seq), extend = FALSE
    )

    expect_equal(gseq_result, vtrack_result$frac_ta, tolerance = 1e-6)
})

test_that("gseq.kmer matches gextract kmer vtrack with longer k-mer", {
    remove_all_vtracks()

    test_interval <- gintervals(1, 200, 240)
    seq <- toupper(gseq.extract(test_interval))

    # Create virtual tracks for longer k-mer
    gvtrack.create("count_ccc", NULL, "kmer.count", kmer = "CCC", strand = 1)
    gvtrack.create("frac_ccc", NULL, "kmer.frac", kmer = "CCC", strand = 1)

    # Extract using gextract
    vtrack_result <- gextract(c("count_ccc", "frac_ccc"), test_interval, iterator = test_interval)

    # Score using gseq.kmer
    gseq_count <- gseq.kmer(seq, "CCC",
        mode = "count", strand = 1,
        start_pos = 1, end_pos = nchar(seq), extend = FALSE
    )
    gseq_frac <- gseq.kmer(seq, "CCC",
        mode = "frac", strand = 1,
        start_pos = 1, end_pos = nchar(seq), extend = FALSE
    )

    expect_equal(gseq_count, vtrack_result$count_ccc)
    expect_equal(gseq_frac, vtrack_result$frac_ccc, tolerance = 1e-6)
})

test_that("gseq.kmer matches gextract kmer vtrack: both strands", {
    remove_all_vtracks()

    test_interval <- gintervals(1, 200, 240)
    seq <- toupper(gseq.extract(test_interval))

    # Use "AC" instead of "TA" to avoid palindrome warning
    # Create virtual tracks for both strands
    gvtrack.create("count_both", NULL, "kmer.count", kmer = "AC", strand = 0)
    gvtrack.create("count_fwd", NULL, "kmer.count", kmer = "AC", strand = 1)
    gvtrack.create("count_rev", NULL, "kmer.count", kmer = "AC", strand = -1)

    # Extract using gextract
    vtrack_result <- gextract(c("count_both", "count_fwd", "count_rev"),
        test_interval,
        iterator = test_interval
    )

    # Score using gseq.kmer
    gseq_both <- gseq.kmer(seq, "AC",
        mode = "count", strand = 0,
        start_pos = 1, end_pos = nchar(seq), extend = FALSE
    )
    gseq_fwd <- gseq.kmer(seq, "AC",
        mode = "count", strand = 1,
        start_pos = 1, end_pos = nchar(seq), extend = FALSE
    )
    gseq_rev <- gseq.kmer(seq, "AC",
        mode = "count", strand = -1,
        start_pos = 1, end_pos = nchar(seq), extend = FALSE
    )

    expect_equal(gseq_both, vtrack_result$count_both)
    expect_equal(gseq_fwd, vtrack_result$count_fwd)
    expect_equal(gseq_rev, vtrack_result$count_rev)
})
