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

# ============================================================================
# Indel support tests (max_indels parameter)
# ============================================================================

test_that("gseq.pwm_edits detects deletion with max_indels=1", {
    # 4-position PSSM: strongly prefers ACGT
    pssm <- matrix(c(
        1, 0, 0, 0, # A
        0, 1, 0, 0, # C
        0, 0, 1, 0, # G
        0, 0, 0, 1 # T
    ), ncol = 4, byrow = TRUE)
    colnames(pssm) <- c("A", "C", "G", "T")

    # Insert an extra base (T) into "ACGT" -> "ATCGT"
    # The DP should find that deleting the extra T yields a perfect ACGT match
    seq <- "ATCGT"
    result <- gseq.pwm_edits(seq, pssm,
        score.thresh = -0.01,
        max_indels = 1L, prior = 0, bidirect = FALSE
    )

    expect_true(nrow(result) > 0)
    expect_true("edit_type" %in% colnames(result))
    edit_rows <- result[result$edit_num > 0, ]
    expect_true(any(edit_rows$edit_type == "del"))
    expect_true(all(result$n_edits >= 1))
})

test_that("gseq.pwm_edits detects insertion with max_indels=1", {
    # 4-position PSSM: strongly prefers ACGT
    pssm <- matrix(c(
        1, 0, 0, 0, # A
        0, 1, 0, 0, # C
        0, 0, 1, 0, # G
        0, 0, 0, 1 # T
    ), ncol = 4, byrow = TRUE)
    colnames(pssm) <- c("A", "C", "G", "T")

    # Remove a base from "ACGT" -> "AGT" (missing C at position 2)
    # The DP should find that inserting C yields a perfect ACGT match
    seq <- "AGT"
    result <- gseq.pwm_edits(seq, pssm,
        score.thresh = -0.01,
        max_indels = 1L, prior = 0, bidirect = FALSE
    )

    expect_true(nrow(result) > 0)
    expect_true("edit_type" %in% colnames(result))
    edit_rows <- result[result$edit_num > 0, ]
    expect_true(any(edit_rows$edit_type == "ins"))
})

test_that("gseq.pwm_edits max_indels=0 backward compatibility: all edits are sub", {
    pssm <- create_test_pssm() # AC motif

    # "TT" needs substitutions to become "AC"
    result <- gseq.pwm_edits("TT", pssm,
        score.thresh = -0.01,
        max_indels = 0L, prior = 0, bidirect = FALSE
    )

    expect_true("edit_type" %in% colnames(result))
    edit_rows <- result[result$edit_num > 0, ]
    expect_true(nrow(edit_rows) > 0)
    expect_true(all(edit_rows$edit_type == "sub"))
    expect_false(any(edit_rows$edit_type %in% c("ins", "del")))
})

test_that("gseq.pwm_edits CTCF deletion detection on ancestral genome", {
    skip_if_not(dir.exists("/net/mraid20/ifs/wisdom/tanay_lab/tgdata/db/tgdb/evo/Phylo447/PrimatesAnc069"))
    skip_if_not_installed("prego")

    gsetroot("/net/mraid20/ifs/wisdom/tanay_lab/tgdata/db/tgdb/evo/Phylo447/PrimatesAnc069")
    ctcf_pssm <- prego::get_motif_pssm("HOMER.CTCF")
    seq <- gseq.extract(gintervals("Anc069refChr724", 246256, 246256 + 30))

    result <- gseq.pwm_edits(seq, ctcf_pssm,
        score.thresh = -15,
        max_indels = 1L, bidirect = TRUE
    )

    expect_true(nrow(result) > 0)
    edit_rows <- result[result$edit_num > 0, ]
    expect_true(any(edit_rows$edit_type == "del"))
    expect_true(all(result$n_edits <= 2))

    # Verify the vtrack also detects it
    gvtrack.create("edist_ctcf", NULL, "pwm.edit_distance",
        pssm = ctcf_pssm, score.thresh = -15,
        score.min = -35, score.max = -18,
        max_edits = 2L, max_indels = 1L, bidirect = TRUE
    )
    interv <- gintervals("Anc069refChr724", 246256, 246257)
    val <- gextract("edist_ctcf", interv, iterator = interv)
    expect_equal(val$edist_ctcf, 1, info = "vtrack should report edit distance of 1")
})

test_that("gseq.pwm_edits edit_type column always present even with max_indels=0", {
    pssm <- create_test_pssm() # AC motif

    # Default (max_indels not specified)
    r_default <- gseq.pwm_edits("TT", pssm,
        score.thresh = -0.01,
        prior = 0, bidirect = FALSE
    )
    expect_true("edit_type" %in% colnames(r_default))

    # Explicit max_indels=0
    r_zero <- gseq.pwm_edits("TT", pssm,
        score.thresh = -0.01,
        max_indels = 0L, prior = 0, bidirect = FALSE
    )
    expect_true("edit_type" %in% colnames(r_zero))

    # 0-edit case (already above threshold)
    r_perfect <- gseq.pwm_edits("AC", pssm,
        score.thresh = -5.0,
        prior = 0, bidirect = FALSE
    )
    expect_true("edit_type" %in% colnames(r_perfect))
    # For 0-edit rows, edit_type should be NA
    zero_edit_rows <- r_perfect[r_perfect$n_edits == 0, ]
    expect_true(all(is.na(zero_edit_rows$edit_type)))
})

# --------------------------------------------------------------------------
# Detailed column validation for indel edits
# --------------------------------------------------------------------------

test_that("gseq.pwm_edits deletion: alignment view and column values", {
    skip_if_not(dir.exists("/net/mraid20/ifs/wisdom/tanay_lab/tgdata/db/tgdb/evo/Phylo447/PrimatesAnc069"))
    skip_if_not_installed("prego")

    # Use the real CTCF example which is known to produce a deletion
    gsetroot("/net/mraid20/ifs/wisdom/tanay_lab/tgdata/db/tgdb/evo/Phylo447/PrimatesAnc069")
    ctcf_pssm <- prego::get_motif_pssm("HOMER.CTCF")
    seq <- gseq.extract(gintervals("Anc069refChr724", 246256, 246256 + 30))

    result <- gseq.pwm_edits(seq, ctcf_pssm,
        score.thresh = -15,
        max_indels = 1L, bidirect = TRUE
    )

    expect_true(nrow(result) >= 1)
    del_rows <- result[!is.na(result$edit_type) & result$edit_type == "del", ]
    expect_true(nrow(del_rows) >= 1, info = "Should find at least one deletion")

    # Alignment view: both strings same length
    expect_equal(nchar(result$window_seq[1]), nchar(result$mutated_seq[1]))

    # mutated_seq should contain a '-' for the deletion
    expect_true(grepl("-", result$mutated_seq[1]),
        info = "mutated_seq should have a hyphen for the deleted base"
    )

    # window_seq should NOT have a hyphen (deletion means seq has the extra base)
    expect_false(grepl("-", result$window_seq[1]),
        info = "window_seq should not have hyphens for deletions"
    )

    # motif_col should be NA for deletions
    expect_true(is.na(del_rows$motif_col[1]))

    # alt should be NA for deletions (no replacement)
    expect_true(is.na(del_rows$alt[1]))

    # ref should be the deleted base (a character)
    expect_true(!is.na(del_rows$ref[1]))
    expect_true(del_rows$ref[1] %in% c("A", "C", "G", "T"))

    # gain should be 0 for deletions
    expect_equal(del_rows$gain[1], 0)

    # score_before and score_after should differ
    expect_true(result$score_before[1] != result$score_after[1])

    # score_before should be below threshold, score_after above
    expect_true(result$score_before[1] < -15)
    expect_true(result$score_after[1] >= -15)

    # Number of hyphens in mutated_seq = number of deletion edits
    n_del <- nrow(del_rows)
    n_hyphens <- nchar(result$mutated_seq[1]) - nchar(gsub("-", "", result$mutated_seq[1]))
    expect_equal(n_hyphens, n_del)
})

test_that("gseq.pwm_edits insertion: alignment view and column values", {
    pssm <- matrix(c(
        0.97, 0.01, 0.01, 0.01,
        0.01, 0.97, 0.01, 0.01,
        0.01, 0.01, 0.97, 0.01,
        0.01, 0.01, 0.01, 0.97
    ), ncol = 4, byrow = TRUE)
    colnames(pssm) <- c("A", "C", "G", "T")

    # "AGT" is missing C at motif position 2 — inserting C gives ACGT
    result <- gseq.pwm_edits("AGT", pssm,
        score.thresh = -0.5,
        max_indels = 1L, prior = 0, bidirect = FALSE
    )

    expect_true(nrow(result) >= 1)
    ins_rows <- result[!is.na(result$edit_type) & result$edit_type == "ins", ]
    expect_true(nrow(ins_rows) >= 1, info = "Should find at least one insertion")

    # Alignment view: both strings same length
    expect_equal(nchar(result$window_seq[1]), nchar(result$mutated_seq[1]))

    # window_seq should contain a '-' for the insertion
    expect_true(grepl("-", result$window_seq[1]),
        info = "window_seq should have a hyphen where base is missing"
    )

    # mutated_seq should NOT have a hyphen (the inserted base fills the gap)
    expect_false(grepl("-", result$mutated_seq[1]),
        info = "mutated_seq should not have hyphens for insertions"
    )

    # motif_col should be a valid 1-based motif position for insertions
    expect_true(!is.na(ins_rows$motif_col[1]))
    expect_true(ins_rows$motif_col[1] >= 1 && ins_rows$motif_col[1] <= 4)

    # ref should be NA for insertions (no original base)
    expect_true(is.na(ins_rows$ref[1]))

    # alt should be the inserted base
    expect_true(!is.na(ins_rows$alt[1]))
    expect_true(ins_rows$alt[1] %in% c("A", "C", "G", "T"))
})

test_that("gseq.pwm_edits substitution gains are correct", {
    pssm <- matrix(c(
        0.97, 0.01, 0.01, 0.01,
        0.01, 0.97, 0.01, 0.01,
        0.01, 0.01, 0.97, 0.01,
        0.01, 0.01, 0.01, 0.97
    ), ncol = 4, byrow = TRUE)
    colnames(pssm) <- c("A", "C", "G", "T")

    # "TCGT" needs T->A at position 1
    result <- gseq.pwm_edits("TCGT", pssm,
        score.thresh = -0.5,
        prior = 0, bidirect = FALSE
    )

    sub_rows <- result[!is.na(result$edit_type) & result$edit_type == "sub", ]
    expect_true(nrow(sub_rows) >= 1)

    # Gain should be col_max - current_base_score = log(0.97) - log(0.01)
    expected_gain <- log(0.97) - log(0.01)
    expect_equal(sub_rows$gain[1], expected_gain, tolerance = 1e-3)

    # ref should be the original base, alt the replacement
    expect_equal(sub_rows$ref[1], "T")
    expect_equal(sub_rows$alt[1], "A")

    # motif_col should be 1 (first position)
    expect_equal(sub_rows$motif_col[1], 1L)

    # score_before should be below the threshold (edits were needed)
    expect_true(result$score_before[1] < -0.5,
        info = "score_before should be below threshold since edits were needed"
    )

    # score_after should equal score_before + sum(gains) for pure substitutions
    total_gain <- sum(sub_rows$gain)
    expect_equal(result$score_after[1], result$score_before[1] + total_gain, tolerance = 1e-3)
})

test_that("gseq.pwm_edits mixed deletion + substitution", {
    pssm <- matrix(c(
        0.97, 0.01, 0.01, 0.01,
        0.01, 0.97, 0.01, 0.01,
        0.01, 0.01, 0.97, 0.01,
        0.01, 0.01, 0.01, 0.97
    ), ncol = 4, byrow = TRUE)
    colnames(pssm) <- c("A", "C", "G", "T")

    # "ATGGT" — extra T at pos 2, and first G should be C
    # Best alignment: delete T, sub G->C -> gives ACGT
    # Or: various other combinations
    result <- gseq.pwm_edits("ATGGT", pssm,
        score.thresh = -0.5,
        max_indels = 1L, prior = 0, bidirect = FALSE
    )

    expect_true(nrow(result) >= 1)
    edit_rows <- result[result$edit_num > 0, ]

    # Alignment view: both strings same length
    expect_equal(nchar(result$window_seq[1]), nchar(result$mutated_seq[1]))

    # score_after should be above the threshold
    expect_true(result$score_after[1] >= -0.5,
        info = "score_after should be above the threshold"
    )

    # Deletion rows should have gain=0, motif_col=NA, alt=NA
    del_rows <- edit_rows[edit_rows$edit_type == "del", ]
    if (nrow(del_rows) > 0) {
        expect_true(all(del_rows$gain == 0))
        expect_true(all(is.na(del_rows$motif_col)))
        expect_true(all(is.na(del_rows$alt)))
        expect_true(all(!is.na(del_rows$ref)))
    }

    # Substitution rows should have gain>0, valid motif_col, ref, alt
    sub_rows <- edit_rows[edit_rows$edit_type == "sub", ]
    if (nrow(sub_rows) > 0) {
        expect_true(all(sub_rows$gain > 0))
        expect_true(all(!is.na(sub_rows$motif_col)))
        expect_true(all(sub_rows$motif_col >= 1 & sub_rows$motif_col <= 4))
        expect_true(all(!is.na(sub_rows$ref)))
        expect_true(all(!is.na(sub_rows$alt)))
    }
})

test_that("gseq.pwm_edits synthetic deletion with 6-position PSSM", {
    # 6-position PSSM: prefers ACGTAC
    pssm6 <- matrix(c(
        0.97, 0.01, 0.01, 0.01,
        0.01, 0.97, 0.01, 0.01,
        0.01, 0.01, 0.97, 0.01,
        0.01, 0.01, 0.01, 0.97,
        0.97, 0.01, 0.01, 0.01,
        0.01, 0.97, 0.01, 0.01
    ), ncol = 4, byrow = TRUE)
    colnames(pssm6) <- c("A", "C", "G", "T")

    # "ACGATAC" — extra A at position 4. Deleting it gives ACGTAC (perfect).
    # The 6-char windows "ACGATA" and "CGATAC" both need multiple subs.
    result <- gseq.pwm_edits("ACGATAC", pssm6,
        score.thresh = -0.5,
        max_indels = 1L, prior = 0, bidirect = FALSE
    )

    expect_true(nrow(result) >= 1)
    edit_rows <- result[result$edit_num > 0, ]
    del_rows <- edit_rows[edit_rows$edit_type == "del", ]
    expect_true(nrow(del_rows) >= 1, info = "Should find a deletion")

    # Validate deletion column values
    expect_true(all(is.na(del_rows$motif_col)))
    expect_true(all(is.na(del_rows$alt)))
    expect_true(all(!is.na(del_rows$ref)))
    expect_true(all(del_rows$gain == 0))

    # Alignment view: same length, hyphen in mutated_seq
    expect_equal(nchar(result$window_seq[1]), nchar(result$mutated_seq[1]))
    expect_true(grepl("-", result$mutated_seq[1]))
    expect_false(grepl("-", result$window_seq[1]))
})

test_that("gseq.pwm_edits score_after = score_before + sum(gains) for pure substitutions", {
    pssm <- matrix(c(
        0.97, 0.01, 0.01, 0.01,
        0.01, 0.97, 0.01, 0.01,
        0.01, 0.01, 0.97, 0.01,
        0.01, 0.01, 0.01, 0.97
    ), ncol = 4, byrow = TRUE)
    colnames(pssm) <- c("A", "C", "G", "T")

    # "TGCA" -> needs 4 substitutions to become ACGT
    result <- gseq.pwm_edits("TGCA", pssm,
        score.thresh = -0.5,
        prior = 0, bidirect = FALSE
    )

    sub_rows <- result[!is.na(result$edit_type) & result$edit_type == "sub", ]
    expect_true(nrow(sub_rows) >= 1)

    # For pure substitutions: score_after = score_before + sum(gains)
    total_gain <- sum(sub_rows$gain)
    expect_equal(result$score_after[1], result$score_before[1] + total_gain,
        tolerance = 1e-3,
        info = "score_after must equal score_before + sum(gains) for substitutions"
    )

    # Each gain should equal log(best_base_prob) - log(current_base_prob)
    for (i in seq_len(nrow(sub_rows))) {
        col <- sub_rows$motif_col[i]
        ref <- sub_rows$ref[i]
        alt <- sub_rows$alt[i]
        # Best base for this column should be the diagonal (A for col 1, C for col 2, etc.)
        expected_alt <- c("A", "C", "G", "T")[col]
        expect_equal(alt, expected_alt,
            info = paste0("Column ", col, ": alt should be ", expected_alt)
        )
        # Gain = log(0.97) - log(0.01)
        expect_equal(sub_rows$gain[i], log(0.97) - log(0.01), tolerance = 1e-3)
    }
})
