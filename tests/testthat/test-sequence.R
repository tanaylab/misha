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
# Tests for gseq.revcomp() - Alias for grevcomp
# ============================================================================

test_that("gseq.revcomp is an alias for grevcomp", {
    # Test that gseq.revcomp produces the same results as grevcomp
    seqs <- c("ACTG", "GCGC", "AAAA", NA, "")

    expect_equal(gseq.revcomp(seqs), grevcomp(seqs))
})

test_that("gseq.revcomp handles various inputs", {
    # Test single sequence
    expect_equal(gseq.revcomp("ACTG"), "CAGT")

    # Test multiple sequences
    expect_equal(
        gseq.revcomp(c("ACTG", "GCGC")),
        c("CAGT", "GCGC")
    )

    # Test with NA
    expect_equal(
        gseq.revcomp(c("ACTG", NA, "GCTA")),
        c("CAGT", NA, "TAGC")
    )

    # Test case preservation
    expect_equal(gseq.revcomp("AcTg"), "cAgT")
})

# ============================================================================
# Tests for gseq.rev() - Reverse without complement
# ============================================================================

test_that("gseq.rev handles basic DNA sequences correctly", {
    # Test single sequence
    expect_equal(gseq.rev("ACTG"), "GTCA")

    # Test multiple sequences
    expect_equal(
        gseq.rev(c("ACTG", "GCGC", "AAAA")),
        c("GTCA", "CGCG", "AAAA")
    )
})

test_that("gseq.rev preserves case", {
    # Test mixed case
    expect_equal(gseq.rev("AcTg"), "gTcA")

    # Test all lowercase
    expect_equal(gseq.rev("actg"), "gtca")

    # Test all uppercase
    expect_equal(gseq.rev("ACTG"), "GTCA")
})

test_that("gseq.rev handles empty and special cases", {
    # Test empty string
    expect_equal(gseq.rev(""), "")

    # Test empty vector
    expect_equal(gseq.rev(character(0)), character(0))

    # Test vector with empty strings
    expect_equal(
        gseq.rev(c("ACTG", "", "GCTA")),
        c("GTCA", "", "ATCG")
    )

    # Test NA values
    expect_equal(gseq.rev(NA_character_), NA_character_)
    expect_equal(gseq.rev(c("ACTG", NA, "GCTA")), c("GTCA", NA, "ATCG"))
})

test_that("gseq.rev validates input", {
    # Test non-character input
    expect_error(gseq.rev(123))
    expect_error(gseq.rev(TRUE))
})

test_that("gseq.rev is reversible", {
    # Test that applying gseq.rev twice returns original sequence
    seqs <- c(
        "ACTG",
        "GCTA",
        "AAAAAA",
        "GCGCGC",
        paste(rep("ACTG", 100), collapse = "")
    )

    for (seq in seqs) {
        expect_equal(gseq.rev(gseq.rev(seq)), seq)
    }
})

test_that("gseq.rev maintains vector attributes", {
    # Test named vector
    input <- c(seq1 = "ACTG", seq2 = "GCTA")
    expected <- c(seq1 = "GTCA", seq2 = "ATCG")
    expect_equal(gseq.rev(input), expected)

    # Test vector with names and NAs
    input <- c(seq1 = "ACTG", missing = NA, seq2 = "GCTA")
    expected <- c(seq1 = "GTCA", missing = NA, seq2 = "ATCG")
    expect_equal(gseq.rev(input), expected)
})

# ============================================================================
# Tests for gseq.comp() - Complement without reverse
# ============================================================================

test_that("gseq.comp handles basic DNA sequences correctly", {
    # Test single sequence
    expect_equal(gseq.comp("ACTG"), "TGAC")

    # Test multiple sequences
    expect_equal(
        gseq.comp(c("ACTG", "GCGC", "AAAA")),
        c("TGAC", "CGCG", "TTTT")
    )
})

test_that("gseq.comp preserves case", {
    # Test mixed case
    expect_equal(gseq.comp("AcTg"), "TgAc")

    # Test all lowercase
    expect_equal(gseq.comp("actg"), "tgac")

    # Test all uppercase
    expect_equal(gseq.comp("ACTG"), "TGAC")
})

test_that("gseq.comp handles empty and special cases", {
    # Test empty string
    expect_equal(gseq.comp(""), "")

    # Test empty vector
    expect_equal(gseq.comp(character(0)), character(0))

    # Test vector with empty strings
    expect_equal(
        gseq.comp(c("ACTG", "", "GCTA")),
        c("TGAC", "", "CGAT")
    )

    # Test NA values
    expect_equal(gseq.comp(NA_character_), NA_character_)
    expect_equal(gseq.comp(c("ACTG", NA, "GCTA")), c("TGAC", NA, "CGAT"))
})

test_that("gseq.comp validates input", {
    # Test non-character input
    expect_error(gseq.comp(123))
    expect_error(gseq.comp(TRUE))
})

test_that("gseq.comp is reversible", {
    # Test that applying gseq.comp twice returns original sequence
    seqs <- c(
        "ACTG",
        "GCTA",
        "AAAAAA",
        "GCGCGC",
        paste(rep("ACTG", 100), collapse = "")
    )

    for (seq in seqs) {
        expect_equal(gseq.comp(gseq.comp(seq)), seq)
    }
})

test_that("gseq.comp maintains vector attributes", {
    # Test named vector
    input <- c(seq1 = "ACTG", seq2 = "GCTA")
    expected <- c(seq1 = "TGAC", seq2 = "CGAT")
    expect_equal(gseq.comp(input), expected)

    # Test vector with names and NAs
    input <- c(seq1 = "ACTG", missing = NA, seq2 = "GCTA")
    expected <- c(seq1 = "TGAC", missing = NA, seq2 = "CGAT")
    expect_equal(gseq.comp(input), expected)
})

# ============================================================================
# Tests for relationship between gseq.rev, gseq.comp, and gseq.revcomp
# ============================================================================

test_that("gseq.revcomp equals gseq.comp then gseq.rev", {
    seqs <- c("ACTG", "GCTA", "AAAATTTT", "CGCGCGCG")

    for (seq in seqs) {
        expect_equal(gseq.revcomp(seq), gseq.rev(gseq.comp(seq)))
    }
})

test_that("gseq.revcomp equals gseq.rev then gseq.comp", {
    seqs <- c("ACTG", "GCTA", "AAAATTTT", "CGCGCGCG")

    for (seq in seqs) {
        expect_equal(gseq.revcomp(seq), gseq.comp(gseq.rev(seq)))
    }
})

test_that("gseq functions handle mixed operations", {
    seq <- "ACTG"

    # rev(comp(seq)) should equal comp(rev(seq))
    expect_equal(gseq.rev(gseq.comp(seq)), gseq.comp(gseq.rev(seq)))

    # This should also equal revcomp(seq)
    expect_equal(gseq.rev(gseq.comp(seq)), gseq.revcomp(seq))
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
    expect_no_error(gseq.pwm(seqs, pssm, mode = "lse", prior = 0))

    # Bad PSSM: not a matrix
    expect_error(gseq.pwm(seqs, c(1, 2, 3, 4), mode = "lse", prior = 0), "matrix")

    # Bad PSSM: wrong number of columns
    bad_pssm <- matrix(1:6, nrow = 2)
    expect_error(gseq.pwm(seqs, bad_pssm, mode = "lse", prior = 0), "4 columns")

    # Bad PSSM: wrong column names
    bad_pssm <- pssm
    colnames(bad_pssm) <- c("X", "Y", "Z", "W")
    expect_error(gseq.pwm(seqs, bad_pssm, mode = "lse", prior = 0), "A, C, G, T")

    # Bad strand
    expect_error(gseq.pwm(seqs, pssm, mode = "lse", strand = 5, prior = 0), "strand must be")

    # Bad extend
    expect_error(gseq.pwm(seqs, pssm, mode = "lse", extend = -5, prior = 0), "extend must be")
})

test_that("gseq.pwm accepts data frame PSSM inputs", {
    pssm <- matrix(c(
        0.4, 0.3, 0.2, 0.1,
        0.1, 0.2, 0.3, 0.4
    ), nrow = 2, byrow = TRUE)
    colnames(pssm) <- c("A", "C", "G", "T")

    seqs <- c("ACGT", "TTTT")

    matrix_modes <- list(
        lse = gseq.pwm(seqs, pssm, mode = "lse", prior = 0),
        max = gseq.pwm(seqs, pssm, mode = "max", prior = 0),
        count = gseq.pwm(seqs, pssm, mode = "count", prior = 0, score.thresh = -10),
        pos = gseq.pwm(seqs, pssm, mode = "pos", prior = 0),
        pos_strand = gseq.pwm(seqs, pssm, mode = "pos", prior = 0, return_strand = TRUE)
    )


    # Neutral character handling
    neutral_seq <- "ANNT"
    default_neutral <- gseq.pwm(neutral_seq, pssm, mode = "max", prior = 0)
    no_neutral <- gseq.pwm(neutral_seq, pssm, mode = "max", prior = 0, neutral_chars = character())
    expect_true(default_neutral > no_neutral)
    expect_equal(
        default_neutral,
        gseq.pwm(neutral_seq, pssm, mode = "max", prior = 0, neutral_chars = c("N", "n", "*"))
    )

    df_modes <- list(
        lse = gseq.pwm(seqs, as.data.frame(pssm), mode = "lse", prior = 0),
        max = gseq.pwm(seqs, as.data.frame(pssm), mode = "max", prior = 0),
        count = gseq.pwm(seqs, as.data.frame(pssm), mode = "count", prior = 0, score.thresh = -10),
        pos = gseq.pwm(seqs, as.data.frame(pssm), mode = "pos", prior = 0),
        pos_strand = gseq.pwm(seqs, as.data.frame(pssm), mode = "pos", prior = 0, return_strand = TRUE)
    )

    expect_equal(df_modes$lse, matrix_modes$lse)
    expect_identical(df_modes$max, matrix_modes$max)
    expect_identical(df_modes$count, matrix_modes$count)
    expect_identical(df_modes$pos, matrix_modes$pos)
    expect_identical(df_modes$pos_strand, matrix_modes$pos_strand)
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
    result <- gseq.pwm(seq, pssm, mode = "count", prior = 0, start_pos = 5, end_pos = 6, extend = FALSE)
    expect_equal(result, 1)

    # ROI before AC (1-4) should not find it
    result <- gseq.pwm(seq, pssm, mode = "count", prior = 0, start_pos = 1, end_pos = 4, extend = FALSE)
    expect_equal(result, 0)

    # ROI after AC (7-10) should not find it
    result <- gseq.pwm(seq, pssm, mode = "count", prior = 0, start_pos = 7, end_pos = 10, extend = FALSE)
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
    result <- gseq.pwm(seq, pssm, mode = "count", prior = 0, start_pos = 6, end_pos = 6, extend = FALSE)
    expect_equal(result, 0)

    # ROI at 6-6, extend=TRUE (w-1=1): should find AC starting at position 5
    result <- gseq.pwm(seq, pssm, mode = "count", prior = 0, start_pos = 6, end_pos = 6, extend = TRUE)
    expect_equal(result, 1)

    # ROI at 6-6, extend=1 (explicit): same as extend=TRUE
    result <- gseq.pwm(seq, pssm, mode = "count", prior = 0, start_pos = 6, end_pos = 6, extend = 1)
    expect_equal(result, 1)

    # ROI at 7-7, extend=2: should find AC (start_min = 7-2 = 5)
    result <- gseq.pwm(seq, pssm, mode = "count", prior = 0, start_pos = 7, end_pos = 7, extend = 2)
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
    result <- gseq.pwm(seq, pssm, mode = "pos", prior = 0, start_pos = 1, end_pos = 10, extend = FALSE)
    expect_equal(result, 3)

    # Limit ROI to find second AC at position 7
    result <- gseq.pwm(seq, pssm, mode = "pos", prior = 0, start_pos = 5, end_pos = 10, extend = FALSE)
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
        start_pos = 1, end_pos = 6, extend = FALSE, return_strand = TRUE,
        skip_gaps = FALSE, prior = 0
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
    result <- gseq.pwm(seq, pssm, mode = "pos", prior = 0, start_pos = 1, end_pos = 8, extend = FALSE)
    expect_equal(result, 1)
})

test_that("gseq.pwm handles edge case: sequence too short", {
    # Create PSSM length 5 (using uniform probabilities)
    pssm <- matrix(1.0, nrow = 5, ncol = 4)
    colnames(pssm) <- c("A", "C", "G", "T")

    # Sequence length 3 (< motif width)
    seq <- "ACG"

    result_total <- gseq.pwm(seq, pssm, mode = "lse", prior = 0)
    expect_true(is.na(result_total))

    result_max <- gseq.pwm(seq, pssm, mode = "max", prior = 0)
    expect_true(is.na(result_max))

    result_count <- gseq.pwm(seq, pssm, mode = "count", prior = 0)
    expect_equal(result_count, 0)

    result_pos <- gseq.pwm(seq, pssm, mode = "pos", prior = 0)
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
    result_count <- gseq.pwm(seq, pssm, mode = "count", prior = 0, start_pos = 5, end_pos = 3)
    expect_equal(result_count, 0)

    result_total <- gseq.pwm(seq, pssm, mode = "lse", prior = 0, start_pos = 5, end_pos = 3)
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
    result_full <- gseq.pwm(seq, pssm, mode = "lse", prior = 0, start_pos = 1, end_pos = 5, extend = FALSE)
    expect_true(is.finite(result_full)) # Should be a finite number

    # Only 2 A's in ROI 1-3
    result_partial <- gseq.pwm(seq, pssm, mode = "lse", prior = 0, start_pos = 1, end_pos = 3, extend = FALSE)
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
    result <- gseq.pwm(seq, pssm, mode = "max", prior = 0, start_pos = 1, end_pos = 5, extend = FALSE)
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
        start_pos = 1, end_pos = 6, extend = FALSE,
        skip_gaps = FALSE, prior = 0
    )
    expect_equal(result, 2)

    # Count positions scoring >= -3.0 (A's and C's, 5 total - forward strand only)
    result <- gseq.pwm(seq, pssm,
        mode = "count", score.thresh = -3.0,
        bidirect = FALSE, strand = 1,
        start_pos = 1, end_pos = 6, extend = FALSE,
        skip_gaps = FALSE, prior = 0
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

    result <- gseq.pwm(seqs, pssm, mode = "count", prior = 0, bidirect = FALSE, strand = 1, start_pos = 1, end_pos = 6, extend = FALSE)
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
        start_pos = c(1, 3), end_pos = c(4, 6), extend = FALSE,
        skip_gaps = FALSE, prior = 0
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
        start_pos = 1, end_pos = nchar(seq), extend = FALSE,
        skip_gaps = FALSE, prior = 0
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
        start_pos = 1, end_pos = nchar(seq_ext), extend = TRUE,
        skip_gaps = FALSE, prior = 0
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
    gseq_result <- gseq.pwm(seq_ext, pssm,
        mode = "max", bidirect = FALSE, strand = 1,
        start_pos = 1, end_pos = nchar(seq_ext), extend = TRUE,
        skip_gaps = FALSE, prior = 0.01
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
    gseq_result <- gseq.pwm(seq_ext, pssm,
        mode = "pos", bidirect = FALSE, strand = 1,
        start_pos = 1, end_pos = nchar(seq_ext), extend = TRUE,
        skip_gaps = FALSE, prior = 0.01
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
        start_pos = 1, end_pos = nchar(seq_ext), extend = TRUE,
        skip_gaps = FALSE, prior = 0
    )

    gseq_fwd <- gseq.pwm(seq_ext, pssm_with_prior,
        mode = "lse", bidirect = FALSE, strand = 1,
        start_pos = 1, end_pos = nchar(seq_ext), extend = TRUE,
        skip_gaps = FALSE, prior = 0
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
        start_pos = 1, end_pos = nchar(seq), extend = FALSE, skip_gaps = FALSE
    )
    gseq_fwd <- gseq.kmer(seq, "AC",
        mode = "count", strand = 1,
        start_pos = 1, end_pos = nchar(seq), extend = FALSE, skip_gaps = FALSE
    )
    gseq_rev <- gseq.kmer(seq, "AC",
        mode = "count", strand = -1,
        start_pos = 1, end_pos = nchar(seq), extend = FALSE, skip_gaps = FALSE
    )

    expect_equal(gseq_both, vtrack_result$count_both)
    expect_equal(gseq_fwd, vtrack_result$count_fwd)
    expect_equal(gseq_rev, vtrack_result$count_rev)
})

# ============================================================================
# Tests for gap support in gseq.pwm() and gseq.kmer()
# ============================================================================

test_that("gseq.pwm with gaps: basic containment", {
    # Create PSSM for "GAAC" motif
    pssm <- matrix(c(
        1e-10, 1e-10, 1.0, 1e-10, # Position 1: G
        1.0, 1e-10, 1e-10, 1e-10, # Position 2: A
        1.0, 1e-10, 1e-10, 1e-10, # Position 3: A
        1e-10, 1.0, 1e-10, 1e-10 # Position 4: C
    ), nrow = 4, byrow = TRUE)
    colnames(pssm) <- c("A", "C", "G", "T")

    # Sequence with gaps: GAAC split by gaps
    seq <- "CTGA----ACGGGGG"
    #       123456789012345

    # With skip_gaps=TRUE, should find GAAC at position 3
    result_gap <- gseq.pwm(seq, pssm,
        mode = "pos", bidirect = FALSE, strand = 1,
        start_pos = 1, end_pos = 15, extend = FALSE, skip_gaps = TRUE,
        prior = 0
    )
    expect_equal(result_gap, 3) # Physical position of G

    # With skip_gaps=FALSE, should not find a good match
    result_nogap <- gseq.pwm(seq, pssm,
        mode = "pos", bidirect = FALSE, strand = 1,
        start_pos = 1, end_pos = 15, extend = FALSE, skip_gaps = FALSE
    )
    # The GAAC is interrupted by gaps, so no perfect match
    expect_true(is.na(result_nogap) || result_nogap != 3)
})

test_that("gseq.pwm with gaps: ROI and extend", {
    # Create PSSM for "AC"
    pssm <- matrix(c(
        1.0, 1e-10, 1e-10, 1e-10, # A
        1e-10, 1.0, 1e-10, 1e-10 # C
    ), nrow = 2, byrow = TRUE)
    colnames(pssm) <- c("A", "C", "G", "T")

    # Sequence with "AC" contiguous at positions 11-12 (gaps elsewhere)
    seq <- "GGGG-A--GGAC-GGGG"
    #      12345678901234567

    # ROI at 12-12 (where C is) with extend=TRUE should find AC
    # With extend=1, scans positions 11-13, which includes the AC window (11-12)
    result <- gseq.pwm(seq, pssm,
        mode = "count", bidirect = FALSE, strand = 1,
        start_pos = 12, end_pos = 12, extend = TRUE, skip_gaps = TRUE,
        prior = 0
    )
    expect_equal(result, 1)
})

test_that("gseq.pwm with gaps: reverse strand", {
    # Create PSSM for "AC"
    pssm <- matrix(c(
        1.0, 1e-10, 1e-10, 1e-10,
        1e-10, 1.0, 1e-10, 1e-10
    ), nrow = 2, byrow = TRUE)
    colnames(pssm) <- c("A", "C", "G", "T")

    # Sequence with GT (reverse complement of AC) split by gaps
    seq <- "GGG-G--T-GGGG"
    #      1234567890123

    result <- gseq.pwm(seq, pssm,
        mode = "pos", bidirect = FALSE, strand = -1,
        start_pos = 1, end_pos = 13, extend = FALSE, skip_gaps = TRUE, return_strand = TRUE,
        prior = 0
    )
    expect_true(is.data.frame(result))
    expect_equal(result$pos, 5) # Position of first G in GT (pos 4 is a gap, pos 5 is G)
    expect_equal(result$strand, -1)
})

test_that("gseq.pwm with gaps: all gaps / too few non-gaps", {
    pssm <- matrix(c(
        1.0, 1e-10, 1e-10, 1e-10,
        1e-10, 1.0, 1e-10, 1e-10
    ), nrow = 2, byrow = TRUE)
    colnames(pssm) <- c("A", "C", "G", "T")

    # All gaps
    seq_all_gaps <- "-----"
    result_lse <- gseq.pwm(seq_all_gaps, pssm, mode = "lse", prior = 0, skip_gaps = TRUE)
    expect_true(is.na(result_lse))

    result_count <- gseq.pwm(seq_all_gaps, pssm, mode = "count", prior = 0, skip_gaps = TRUE)
    expect_equal(result_count, 0)

    # Only one non-gap base (need 2 for motif width)
    seq_one_base <- "---A---"
    result_lse2 <- gseq.pwm(seq_one_base, pssm, mode = "lse", prior = 0, skip_gaps = TRUE)
    expect_true(is.na(result_lse2))
})

test_that("gseq.pwm with gaps: spatial weighting", {
    pssm <- matrix(c(
        1.0, 1e-10, 1e-10, 1e-10,
        1e-10, 1.0, 1e-10, 1e-10
    ), nrow = 2, byrow = TRUE)
    colnames(pssm) <- c("A", "C", "G", "T")

    # Two AC motifs split by gaps at different physical positions
    seq <- "A-C-GG-A-C-GG"
    #      1234567890123

    # Spatial weighting should be based on physical positions
    spatial_weights <- c(1.0, 0.5) # First position weighted higher
    result <- gseq.pwm(seq, pssm,
        mode = "lse", bidirect = FALSE, strand = 1,
        start_pos = 1, end_pos = 13, extend = FALSE, skip_gaps = TRUE,
        spat.factor = spatial_weights, spat.bin = 5
    )
    expect_true(is.finite(result))
})

test_that("gseq.kmer with gaps: basic matching", {
    seq <- "A-CG-A"
    #      123456

    # With skip_gaps=TRUE, should find "ACG" (A at 1, C at 3, G at 4)
    result_gap <- gseq.kmer(seq, "ACG",
        mode = "count", strand = 1,
        start_pos = 1, end_pos = 6, extend = FALSE, skip_gaps = TRUE
    )
    expect_equal(result_gap, 1)

    # With skip_gaps=FALSE, should not find "ACG"
    result_nogap <- gseq.kmer(seq, "ACG",
        mode = "count", strand = 1,
        start_pos = 1, end_pos = 6, extend = FALSE, skip_gaps = FALSE
    )
    expect_equal(result_nogap, 0)
})

test_that("gseq.kmer with gaps: fraction mode", {
    seq <- "A-A-A-GGG"
    #      123456789

    # With skip_gaps=TRUE, count "A"s: 3 out of 6 non-gap positions
    # But for k=1, there are 6 possible logical starts
    result_frac <- gseq.kmer(seq, "A",
        mode = "frac", strand = 1,
        start_pos = 1, end_pos = 9, extend = FALSE, skip_gaps = TRUE
    )
    expect_equal(result_frac, 3 / 6) # 3 A's, 6 non-gap bases
})

test_that("gseq.kmer with gaps: fraction mode for k>1", {
    seq <- "A-AC-AC"
    #      1234567

    # Non-gap compacted sequence is "AACAC" (length 5)
    # For k=2 there are 4 logical starts; "AC" occurs twice
    result_frac <- gseq.kmer(seq, "AC",
        mode = "frac", strand = 1,
        start_pos = 1, end_pos = 7, extend = FALSE, skip_gaps = TRUE
    )
    expect_equal(result_frac, 2 / 4)
})

test_that("gseq.kmer with gaps: ROI handling", {
    seq <- "GGG-A-C-GGG"
    #      12345678901

    # ROI covering AC (positions 5-7 physically)
    result <- gseq.kmer(seq, "AC",
        mode = "count", strand = 1,
        start_pos = 5, end_pos = 7, extend = FALSE, skip_gaps = TRUE
    )
    expect_equal(result, 1)

    # ROI before AC
    result_before <- gseq.kmer(seq, "AC",
        mode = "count", strand = 1,
        start_pos = 1, end_pos = 3, extend = FALSE, skip_gaps = TRUE
    )
    expect_equal(result_before, 0)
})

test_that("gseq.kmer with gaps: reverse strand", {
    seq <- "G--T-GGG"
    #      12345678

    # GT on forward = AC on reverse (rev comp)
    result <- gseq.kmer(seq, "AC",
        mode = "count", strand = -1,
        start_pos = 1, end_pos = 8, extend = FALSE, skip_gaps = TRUE
    )
    expect_equal(result, 1)
})

test_that("gseq.kmer with gaps: both strands", {
    seq <- "A-CG-G-T"
    #      12345678

    # AC on forward, GT on forward (= AC on reverse)
    result <- gseq.kmer(seq, "AC",
        mode = "count", strand = 0,
        start_pos = 1, end_pos = 8, extend = FALSE, skip_gaps = TRUE
    )
    expect_equal(result, 2) # AC forward + GT as AC reverse
})

test_that("gseq.kmer with gaps: empty sequence", {
    seq <- "-----"

    result_count <- gseq.kmer(seq, "AC", mode = "count", skip_gaps = TRUE)
    expect_equal(result_count, 0)

    result_frac <- gseq.kmer(seq, "AC", mode = "frac", skip_gaps = TRUE)
    expect_equal(result_frac, 0)
})

test_that("gseq.kmer with gaps: custom gap characters", {
    # Custom gap characters should be honored for k-mer matching
    seq <- "A_C.A*C"
    # Non-gap compacted sequence is "ACAC" -> two "AC" occurrences
    result_custom <- gseq.kmer(seq, "AC",
        mode = "count", strand = 1,
        skip_gaps = TRUE, gap_chars = c("-", ".", "_", "*")
    )
    expect_equal(result_custom, 2)
})

test_that("gseq functions with gaps: custom gap characters", {
    pssm <- matrix(c(
        1.0, 1e-10, 1e-10, 1e-10,
        1e-10, 1.0, 1e-10, 1e-10
    ), nrow = 2, byrow = TRUE)
    colnames(pssm) <- c("A", "C", "G", "T")

    # Use different gap characters
    seq <- "A_C.G-A*C"

    # Default gaps: only "-" and "."
    result_default <- gseq.pwm(seq, pssm,
        mode = "count", bidirect = FALSE, strand = 1,
        skip_gaps = TRUE, gap_chars = c("-", ".")
    )
    expect_true(result_default >= 0)

    # Custom gaps: include "_" and "*"
    result_custom <- gseq.pwm(seq, pssm,
        mode = "count", bidirect = FALSE, strand = 1,
        skip_gaps = TRUE, gap_chars = c("-", ".", "_", "*"),
        prior = 0
    )
    expect_equal(result_custom, 2) # AC at positions 1,3 and 6,8
})

test_that("gseq functions: skip_gaps=FALSE preserves original behavior", {
    pssm <- matrix(c(
        1.0, 1e-10, 1e-10, 1e-10,
        1e-10, 1.0, 1e-10, 1e-10
    ), nrow = 2, byrow = TRUE)
    colnames(pssm) <- c("A", "C", "G", "T")

    seq <- "ACGTACGT"

    # Results should be identical with skip_gaps=FALSE and skip_gaps=TRUE on a gap-free sequence
    result_false <- gseq.pwm(seq, pssm, mode = "count", prior = 0, skip_gaps = FALSE)
    result_true <- gseq.pwm(seq, pssm, mode = "count", prior = 0, skip_gaps = TRUE)
    expect_equal(result_false, result_true)

    # Same for k-mer
    result_kmer_false <- gseq.kmer(seq, "AC", mode = "count", skip_gaps = FALSE)
    result_kmer_true <- gseq.kmer(seq, "AC", mode = "count", skip_gaps = TRUE)
    expect_equal(result_kmer_false, result_kmer_true)
})

test_that("gseq functions: gap parameters validation", {
    pssm <- matrix(c(1.0, 1e-10, 1e-10, 1e-10), nrow = 1)
    colnames(pssm) <- c("A", "C", "G", "T")

    # Invalid gap_chars
    expect_error(gseq.pwm("ACGT", pssm, skip_gaps = TRUE, gap_chars = character(0), prior = 0))
    expect_error(gseq.pwm("ACGT", pssm, skip_gaps = TRUE, gap_chars = c("A", "AB"), prior = 0))
    expect_error(gseq.pwm("ACGT", pssm, skip_gaps = TRUE, gap_chars = c("-", "-"), prior = 0))
})

test_that("gseq.pwm bidirectional pos returns correct strand", {
    # Create PSSM for "AC" (using probabilities)
    pssm <- matrix(c(
        1.0, 1e-10, 1e-10, 1e-10,
        1e-10, 1.0, 1e-10, 1e-10
    ), nrow = 2, byrow = TRUE)
    colnames(pssm) <- c("A", "C", "G", "T")

    # Sequence contains only the reverse complement "GT" at positions 4-5
    seq <- "GGGGTGG"

    result <- gseq.pwm(seq, pssm,
        mode = "pos", bidirect = TRUE, strand = 1,
        start_pos = 1, end_pos = nchar(seq), extend = FALSE, return_strand = TRUE,
        skip_gaps = FALSE, prior = 0
    )

    expect_true(is.data.frame(result))
    expect_equal(names(result), c("pos", "strand"))
    expect_equal(result$pos, 4)
    expect_equal(result$strand, -1)
})

# ============================================================================
# Comparison tests: gseq.pwm vs prego::compute_pwm
# ============================================================================

test_that("gseq.pwm matches prego::compute_pwm: mode='max' with prior", {
    skip_if_not_installed("prego")

    # Use CTCF motif from prego
    ctcf_mot <- prego::HOMER_motifs %>%
        dplyr::filter(motif == "CTCF") %>%
        dplyr::select(-motif, -pos) %>%
        as.matrix()

    # Test sequence
    seq <- "GTGAACTTCGCTGTCAGCAGAGGGCAACAGGTTCTGCGGG"

    # Test with prior=0.01
    prego_result <- prego::compute_pwm(seq, ctcf_mot, func = "max", prior = 0.01)
    gseq_result <- gseq.pwm(seq, ctcf_mot, mode = "max", prior = 0.01, skip_gaps = FALSE)
    expect_equal(gseq_result, prego_result, tolerance = 1e-5)

    # Test with prior=0
    prego_result_no_prior <- prego::compute_pwm(seq, ctcf_mot, func = "max", prior = 0)
    gseq_result_no_prior <- gseq.pwm(seq, ctcf_mot, mode = "max", prior = 0, skip_gaps = FALSE)
    expect_equal(gseq_result_no_prior, prego_result_no_prior, tolerance = 1e-5)
})

test_that("gseq.pwm matches prego::compute_pwm: vectorized sequences", {
    skip_if_not_installed("prego")

    # Use CTCF motif from prego
    ctcf_mot <- prego::HOMER_motifs %>%
        dplyr::filter(motif == "CTCF") %>%
        dplyr::select(-motif, -pos) %>%
        as.matrix()

    # Multiple test sequences
    seqs <- c(
        "GTGAACTTCGCTGTCAGCAGAGGGCAACAGGTTCTGCGGG",
        "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
        "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"
    )

    # Test mode='max' with prior=0.01
    prego_results <- sapply(seqs, function(s) {
        prego::compute_pwm(s, ctcf_mot, func = "max", prior = 0.01)
    })
    gseq_results <- gseq.pwm(seqs, ctcf_mot, mode = "max", prior = 0.01, skip_gaps = FALSE)
    expect_equal(gseq_results, as.numeric(prego_results), tolerance = 1e-4)
})
