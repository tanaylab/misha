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
