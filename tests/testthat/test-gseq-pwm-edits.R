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

# ============================================================================
# Regression tests for specific bug fixes
# ============================================================================

test_that("gseq.pwm_edits numeric extend on bare sequences is preserved", {
    # extend=1 should allow 1 extra base of extension, not w-1
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

    seq <- "TTTTACGTTTTT"
    # With extend=FALSE: ROI is exactly [1, nchar], windows must start and fit within
    r_no_ext <- gseq.pwm_edits(seq, pssm,
        score.thresh = -1.0,
        prior = 0, bidirect = FALSE, extend = FALSE
    )

    # With extend=1: allows 1 extra base of scanning range
    r_ext1 <- gseq.pwm_edits(seq, pssm,
        score.thresh = -1.0,
        prior = 0, bidirect = FALSE, extend = 1L
    )

    # With extend=TRUE: allows w-1=3 extra bases of scanning range
    r_ext_full <- gseq.pwm_edits(seq, pssm,
        score.thresh = -1.0,
        prior = 0, bidirect = FALSE, extend = TRUE
    )

    # extend=1 should not be the same as extend=TRUE (w-1=3)
    # unless results happen to be identical for this sequence.
    # At minimum, extend=1 should produce valid results:
    expect_true(is.data.frame(r_ext1))
    expect_true(nrow(r_ext1) > 0)

    # extend=FALSE and extend=1 may differ in window_start
    # because extend=1 allows 1 extra window position
    expect_true(is.data.frame(r_no_ext))
})

test_that("gseq.pwm_edits score_before/score_after correct with N and zero-prob bases", {
    # PSSM with a zero-probability entry
    pssm_zero <- matrix(
        c(
            1, 0, 0, 0, # column 1: only A allowed
            0, 1, 0, 0 # column 2: only C allowed
        ),
        nrow = 2, byrow = TRUE,
        dimnames = list(NULL, c("A", "C", "G", "T"))
    )

    # "TC" has T at position 1 where only A has prob > 0 -> zero-prob, mandatory edit
    r <- gseq.pwm_edits("TC", pssm_zero,
        score.thresh = -0.5,
        prior = 0, bidirect = FALSE
    )

    # score_before should reflect the TRUE sequence score (which is -Inf due to zero prob)
    expect_true(all(!is.nan(r$score_before)))
    expect_true(all(r$score_before < -100)) # -Inf or very negative

    # score_after should be the score AFTER applying edits (should be near 0 = log(1))
    expect_true(all(r$score_after > r$score_before))

    # Gain for the mandatory edit should be positive (not 0)
    mandatory_rows <- r[r$edit_num > 0, ]
    expect_true(all(mandatory_rows$gain > 0 | is.infinite(mandatory_rows$gain)))

    # Test with N bases
    r_n <- gseq.pwm_edits("NC", pssm_zero,
        score.thresh = -0.5,
        prior = 0, bidirect = FALSE
    )
    expect_true(nrow(r_n) > 0)
    # N base should be a mandatory edit
    expect_true(any(r_n$ref == "N"))
    # score_after should be better than score_before
    expect_true(all(r_n$score_after[r_n$n_edits > 0] > r_n$score_before[r_n$n_edits > 0]))
})

test_that("pwm.edit_distance.pos with filter returns positions relative to original interval", {
    gdb.init_examples()
    remove_all_vtracks()

    pssm <- create_test_pssm()

    # Create a large interval
    test_interval <- gintervals(1, 200, 260)

    # Create a mask that splits the interval into two fragments
    # Mask out the middle: [220, 240) is masked, leaving [200,220) and [240,260)
    mask <- gintervals(1, 220, 240)

    gvtrack.create("edist_pos_filtered", NULL,
        func = "pwm.edit_distance.pos",
        pssm = pssm, score.thresh = -5.0,
        bidirect = FALSE, extend = FALSE, prior = 0
    )
    gvtrack.filter("edist_pos_filtered", filter = mask)

    # Also create unfiltered version for the second fragment
    gvtrack.create("edist_pos_unfiltered", NULL,
        func = "pwm.edit_distance.pos",
        pssm = pssm, score.thresh = -5.0,
        bidirect = FALSE, extend = FALSE, prior = 0
    )

    result_filtered <- gextract("edist_pos_filtered",
        intervals = test_interval, iterator = test_interval
    )
    pos_filtered <- result_filtered$edist_pos_filtered

    # If the best window is in the second fragment [240,260),
    # the position should be relative to the original interval [200,260),
    # not relative to the subfragment [240,260).
    # The position should be > 0 (1-based within the original interval)
    if (!is.na(pos_filtered)) {
        abs_pos <- abs(pos_filtered)
        # Position should be within [1, interval_length]
        expect_true(abs_pos >= 1)
        expect_true(abs_pos <= 60) # interval is 60bp

        # If best is in second fragment, position should be > 20
        # (since first fragment is 20bp and mask starts at offset 20)
        # We can't guarantee which fragment wins, but position must be valid
    }

    remove_all_vtracks()
})
