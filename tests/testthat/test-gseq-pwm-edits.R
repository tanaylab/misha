test_that("gseq.pwm_edits returns correct structure for bare sequences", {
    pssm <- matrix(
        c(
            0.9, 0.05, 0.025, 0.025,
            0.05, 0.9, 0.025, 0.025
        ),
        nrow = 2, byrow = TRUE,
        dimnames = list(NULL, c("A", "C", "G", "T"))
    )

    result <- gseq.pwm_edits("CCGTACGT", pssm, score.thresh = -0.5, prior = 0)

    expect_true(is.data.frame(result))
    expect_true(all(c(
        "seq_idx", "strand", "window_start", "score_before",
        "score_after", "n_edits", "edit_num", "motif_col",
        "ref", "alt", "gain", "window_seq", "mutated_seq"
    ) %in% colnames(result)))
    expect_true(nrow(result) > 0)
})

test_that("gseq.pwm_edits returns 0 edits when already above threshold", {
    pssm <- matrix(
        c(
            0.9, 0.05, 0.025, 0.025,
            0.05, 0.9, 0.025, 0.025
        ),
        nrow = 2, byrow = TRUE,
        dimnames = list(NULL, c("A", "C", "G", "T"))
    )

    # AC is the perfect match for this PSSM
    result <- gseq.pwm_edits("ACGTACGT", pssm, score.thresh = -5.0, prior = 0)

    expect_equal(result$n_edits[1], 0L)
    expect_equal(result$edit_num[1], 0L)
    expect_true(is.na(result$motif_col[1]))
})

test_that("gseq.pwm_edits identifies correct edits for single mutation", {
    # PSSM strongly prefers A at pos 1, C at pos 2
    pssm <- matrix(
        c(
            1, 0, 0, 0,
            0, 1, 0, 0
        ),
        nrow = 2, byrow = TRUE,
        dimnames = list(NULL, c("A", "C", "G", "T"))
    )

    # "TC" needs 1 edit: T->A at position 1
    result <- gseq.pwm_edits("TC", pssm,
        score.thresh = -0.01,
        prior = 0, bidirect = FALSE
    )

    one_edit_rows <- result[result$n_edits == 1, ]
    expect_equal(nrow(one_edit_rows), 1)
    expect_equal(one_edit_rows$motif_col, 1L)
    expect_equal(one_edit_rows$ref, "T")
    expect_equal(one_edit_rows$alt, "A")
})

test_that("gseq.pwm_edits window_seq and mutated_seq are correct", {
    pssm <- matrix(
        c(
            1, 0, 0, 0,
            0, 1, 0, 0
        ),
        nrow = 2, byrow = TRUE,
        dimnames = list(NULL, c("A", "C", "G", "T"))
    )

    result <- gseq.pwm_edits("TC", pssm,
        score.thresh = -0.01,
        prior = 0, bidirect = FALSE
    )

    expect_equal(result$window_seq[1], "TC")
    expect_equal(result$mutated_seq[1], "AC")
})

test_that("gseq.pwm_edits works with multiple sequences", {
    pssm <- matrix(
        c(
            0.9, 0.05, 0.025, 0.025,
            0.05, 0.9, 0.025, 0.025
        ),
        nrow = 2, byrow = TRUE,
        dimnames = list(NULL, c("A", "C", "G", "T"))
    )

    result <- gseq.pwm_edits(c("AC", "TT", "GG"), pssm,
        score.thresh = -1.0, prior = 0, bidirect = FALSE
    )

    # AC should be 0 edits, TT and GG should need edits
    expect_true(any(result$seq_idx == 1 & result$n_edits == 0))
    expect_true(any(result$seq_idx == 2 & result$n_edits >= 1))
    expect_true(any(result$seq_idx == 3 & result$n_edits >= 1))
})

test_that("gseq.pwm_edits with intervals includes chrom/start/end", {
    gdb.init_examples()
    pssm <- matrix(
        c(
            0.9, 0.05, 0.025, 0.025,
            0.05, 0.9, 0.025, 0.025
        ),
        nrow = 2, byrow = TRUE,
        dimnames = list(NULL, c("A", "C", "G", "T"))
    )

    intervals <- gintervals(1, 200, 210)
    result <- gseq.pwm_edits(intervals, pssm, score.thresh = -3.0, prior = 0)

    expect_true("chrom" %in% colnames(result))
    expect_true("start" %in% colnames(result))
    expect_true("end" %in% colnames(result))
    expect_true(nrow(result) > 0)
})

test_that("gseq.pwm_edits respects max_edits cap", {
    # Use non-zero PSSM so edits are optional (not mandatory from zero-prob)
    pssm <- matrix(
        c(
            0.9, 0.05, 0.025, 0.025,
            0.05, 0.9, 0.025, 0.025,
            0.025, 0.025, 0.9, 0.05,
            0.025, 0.025, 0.05, 0.9
        ),
        nrow = 4, byrow = TRUE,
        dimnames = list(NULL, c("A", "C", "G", "T"))
    )

    # TTTT needs multiple edits, but max_edits=1
    r1 <- gseq.pwm_edits("TTTT", pssm,
        score.thresh = -0.5,
        max_edits = 1L, prior = 0, bidirect = FALSE
    )
    # With max_edits=4, should find a solution
    r4 <- gseq.pwm_edits("TTTT", pssm,
        score.thresh = -0.5,
        max_edits = 4L, prior = 0, bidirect = FALSE
    )

    # max_edits=1 should return fewer (or no) results vs max_edits=4
    expect_true(nrow(r1) <= nrow(r4))
})

test_that("gseq.pwm_edits works with bidirectional scanning", {
    pssm <- matrix(
        c(
            1, 0, 0, 0,
            0, 1, 0, 0
        ),
        nrow = 2, byrow = TRUE,
        dimnames = list(NULL, c("A", "C", "G", "T"))
    )

    # "GT" reverse complement is "AC" — perfect match on reverse strand
    result <- gseq.pwm_edits("GT", pssm,
        score.thresh = -0.01,
        prior = 0, bidirect = TRUE
    )

    expect_equal(result$n_edits[1], 0L)
    expect_equal(result$strand[1], -1L)
})

test_that("gseq.pwm_edits score.min and score.max filters work", {
    pssm <- matrix(
        c(
            0.9, 0.05, 0.025, 0.025,
            0.05, 0.9, 0.025, 0.025
        ),
        nrow = 2, byrow = TRUE,
        dimnames = list(NULL, c("A", "C", "G", "T"))
    )

    # With score.max very low — should filter out ALL windows and return empty
    r_max <- gseq.pwm_edits("CCGTACGT", pssm,
        score.thresh = -0.5,
        prior = 0, score.max = -50.0
    )
    expect_equal(nrow(r_max), 0)

    # With score.min very high — should also filter out all and return empty
    r_min <- gseq.pwm_edits("CCGTACGT", pssm,
        score.thresh = -0.5,
        prior = 0, score.min = 0.0
    )
    expect_equal(nrow(r_min), 0)
})

test_that("gseq.pwm_edits returns empty data frame for empty input", {
    pssm <- matrix(
        c(
            1, 0, 0, 0,
            0, 1, 0, 0
        ),
        nrow = 2, byrow = TRUE,
        dimnames = list(NULL, c("A", "C", "G", "T"))
    )

    result <- gseq.pwm_edits(character(0), pssm, score.thresh = -1.0)

    expect_true(is.data.frame(result))
    expect_equal(nrow(result), 0)
    expect_true("window_seq" %in% colnames(result))
    expect_true("mutated_seq" %in% colnames(result))
})

test_that("gseq.pwm_edits mutated_seq has correct bases at edit positions", {
    pssm <- matrix(
        c(
            1, 0, 0, 0,
            0, 1, 0, 0,
            0, 0, 1, 0
        ),
        nrow = 3, byrow = TRUE,
        dimnames = list(NULL, c("A", "C", "G", "T"))
    )

    # "TTT" -> need to change all 3 to "ACG"
    result <- gseq.pwm_edits("TTT", pssm,
        score.thresh = -0.01,
        prior = 0, bidirect = FALSE
    )

    # All rows should have the same mutated_seq
    expect_true(all(result$mutated_seq == "ACG"))
    expect_true(all(result$window_seq == "TTT"))

    # Verify each edit individually
    for (i in seq_len(nrow(result))) {
        mc <- result$motif_col[i]
        expect_equal(substr(result$window_seq[i], mc, mc), result$ref[i])
        expect_equal(substr(result$mutated_seq[i], mc, mc), result$alt[i])
    }
})

test_that("gseq.pwm_edits gain values are positive for real edits", {
    pssm <- matrix(
        c(
            0.9, 0.05, 0.025, 0.025,
            0.05, 0.9, 0.025, 0.025
        ),
        nrow = 2, byrow = TRUE,
        dimnames = list(NULL, c("A", "C", "G", "T"))
    )

    result <- gseq.pwm_edits("TT", pssm,
        score.thresh = -0.5,
        prior = 0, bidirect = FALSE
    )

    real_edits <- result[result$edit_num > 0, ]
    expect_true(all(real_edits$gain > 0))
})

test_that("gseq.pwm_edits score_after > score_before when edits needed", {
    pssm <- matrix(
        c(
            0.9, 0.05, 0.025, 0.025,
            0.05, 0.9, 0.025, 0.025
        ),
        nrow = 2, byrow = TRUE,
        dimnames = list(NULL, c("A", "C", "G", "T"))
    )

    result <- gseq.pwm_edits("TT", pssm,
        score.thresh = -0.5,
        prior = 0, bidirect = FALSE
    )

    rows_with_edits <- result[result$n_edits > 0, ]
    if (nrow(rows_with_edits) > 0) {
        expect_true(all(rows_with_edits$score_after > rows_with_edits$score_before))
    }
})

test_that("gseq.pwm_edits validates parameters", {
    pssm <- matrix(c(1, 0, 0, 0),
        nrow = 1,
        dimnames = list(NULL, c("A", "C", "G", "T"))
    )

    expect_error(gseq.pwm_edits("ACGT", pssm, score.thresh = "bad"))
    expect_error(gseq.pwm_edits("ACGT", pssm, score.thresh = -1, max_edits = 0L))
    expect_error(gseq.pwm_edits("ACGT", pssm, score.thresh = -1, bidirect = "yes"))
})

test_that("gseq.pwm_edits with multiple intervals", {
    gdb.init_examples()
    pssm <- matrix(
        c(
            0.9, 0.05, 0.025, 0.025,
            0.05, 0.9, 0.025, 0.025
        ),
        nrow = 2, byrow = TRUE,
        dimnames = list(NULL, c("A", "C", "G", "T"))
    )

    intervals <- rbind(gintervals(1, 200, 210), gintervals(1, 300, 310))
    result <- gseq.pwm_edits(intervals, pssm, score.thresh = -3.0, prior = 0)

    expect_true(nrow(result) > 0)
    # Both intervals should have results
    expect_true(1 %in% result$seq_idx)
    expect_true(2 %in% result$seq_idx)
})

test_that("gseq.pwm_edits with pssm data frame", {
    pssm_df <- data.frame(
        A = c(0.9, 0.05), C = c(0.05, 0.9),
        G = c(0.025, 0.025), T = c(0.025, 0.025)
    )

    result <- gseq.pwm_edits("TT", pssm_df,
        score.thresh = -0.5,
        prior = 0, bidirect = FALSE
    )
    expect_true(is.data.frame(result))
    expect_true(nrow(result) > 0)
})
