# ============================================================================
# Tests for gseq.pwm_edits() with direction="below"
# ============================================================================

test_that("gseq.pwm_edits direction=below returns correct structure", {
    pssm <- matrix(
        c(
            0.97, 0.01, 0.01, 0.01,
            0.01, 0.97, 0.01, 0.01,
            0.01, 0.01, 0.97, 0.01,
            0.01, 0.01, 0.01, 0.97
        ),
        nrow = 4, byrow = TRUE,
        dimnames = list(NULL, c("A", "C", "G", "T"))
    )

    # "ACGT" is a perfect match — its score is well above any negative threshold
    result <- gseq.pwm_edits("ACGT", pssm,
        score.thresh = -1.0,
        prior = 0, bidirect = FALSE, direction = "below"
    )

    expect_true(is.data.frame(result))
    expect_true(all(c(
        "seq_idx", "strand", "window_start", "score_before",
        "score_after", "n_edits", "edit_num", "motif_col",
        "ref", "alt", "gain", "window_seq", "mutated_seq"
    ) %in% colnames(result)))
    expect_true(nrow(result) > 0)
    # Edits should be needed since score_before >> threshold
    expect_true(any(result$n_edits > 0))
})

test_that("gseq.pwm_edits direction=below returns 0 edits when already below threshold", {
    pssm <- matrix(
        c(
            0.97, 0.01, 0.01, 0.01,
            0.01, 0.97, 0.01, 0.01
        ),
        nrow = 2, byrow = TRUE,
        dimnames = list(NULL, c("A", "C", "G", "T"))
    )

    # "TT" scores very low against this PSSM (prefers AC)
    # With a threshold high enough that score is already below
    result <- gseq.pwm_edits("TT", pssm,
        score.thresh = -0.5,
        prior = 0, bidirect = FALSE, direction = "below"
    )

    expect_equal(result$n_edits[1], 0L)
    expect_equal(result$edit_num[1], 0L)
    expect_true(is.na(result$motif_col[1]))
})

test_that("gseq.pwm_edits direction=below: score_before > threshold and score_after <= threshold", {
    pssm <- matrix(
        c(
            0.97, 0.01, 0.01, 0.01,
            0.01, 0.97, 0.01, 0.01,
            0.01, 0.01, 0.97, 0.01,
            0.01, 0.01, 0.01, 0.97
        ),
        nrow = 4, byrow = TRUE,
        dimnames = list(NULL, c("A", "C", "G", "T"))
    )

    # "ACGT" is a perfect match — score is near 4*log(0.97) ~ -0.12
    result <- gseq.pwm_edits("ACGT", pssm,
        score.thresh = -1.0,
        prior = 0, bidirect = FALSE, direction = "below"
    )

    rows_with_edits <- result[result$n_edits > 0, ]
    expect_true(nrow(rows_with_edits) > 0, info = "Should need edits to bring score below threshold")

    # score_before should be above the threshold (that is why edits are needed)
    expect_true(all(rows_with_edits$score_before > -1.0),
        info = "score_before should be above threshold for below direction"
    )

    # score_after should be at or below the threshold
    expect_true(all(rows_with_edits$score_after <= -1.0 + 1e-9),
        info = "score_after should be at or below threshold for below direction"
    )
})

test_that("gseq.pwm_edits direction=below: edits reference bases match window_seq", {
    pssm <- matrix(
        c(
            0.97, 0.01, 0.01, 0.01,
            0.01, 0.97, 0.01, 0.01,
            0.01, 0.01, 0.97, 0.01
        ),
        nrow = 3, byrow = TRUE,
        dimnames = list(NULL, c("A", "C", "G", "T"))
    )

    # "ACG" is a perfect match
    result <- gseq.pwm_edits("ACG", pssm,
        score.thresh = -1.0,
        prior = 0, bidirect = FALSE, direction = "below"
    )

    edit_rows <- result[result$edit_num > 0, ]
    expect_true(nrow(edit_rows) > 0)

    for (i in seq_len(nrow(edit_rows))) {
        if (!is.na(edit_rows$edit_type[i]) && edit_rows$edit_type[i] == "sub") {
            mc <- edit_rows$motif_col[i]
            # ref should match the base at motif_col in window_seq
            expect_equal(
                substr(edit_rows$window_seq[i], mc, mc),
                edit_rows$ref[i],
                info = paste0("Edit ", i, ": ref should match window_seq at motif_col ", mc)
            )
        }
    }

    # Applying all edits to window_seq should produce mutated_seq
    # For substitution-only edits, build mutated manually
    ws <- result$window_seq[1]
    ms <- result$mutated_seq[1]
    ws_chars <- strsplit(ws, "")[[1]]
    for (i in seq_len(nrow(edit_rows))) {
        if (!is.na(edit_rows$edit_type[i]) && edit_rows$edit_type[i] == "sub") {
            mc <- edit_rows$motif_col[i]
            ws_chars[mc] <- edit_rows$alt[i]
        }
    }
    expect_equal(paste0(ws_chars, collapse = ""), ms,
        info = "Applying all edits to window_seq should produce mutated_seq"
    )
})

test_that("gseq.pwm_edits direction=below: gain values are negative for real edits", {
    pssm <- matrix(
        c(
            0.97, 0.01, 0.01, 0.01,
            0.01, 0.97, 0.01, 0.01,
            0.01, 0.01, 0.97, 0.01,
            0.01, 0.01, 0.01, 0.97
        ),
        nrow = 4, byrow = TRUE,
        dimnames = list(NULL, c("A", "C", "G", "T"))
    )

    # "ACGT" is a perfect match — edits needed to destroy it
    result <- gseq.pwm_edits("ACGT", pssm,
        score.thresh = -1.0,
        prior = 0, bidirect = FALSE, direction = "below"
    )

    real_edits <- result[result$edit_num > 0, ]
    expect_true(nrow(real_edits) > 0)

    sub_edits <- real_edits[!is.na(real_edits$edit_type) & real_edits$edit_type == "sub", ]
    if (nrow(sub_edits) > 0) {
        # In "below" direction, gain should be negative (score is being reduced)
        expect_true(all(sub_edits$gain < 0),
            info = "gain should be negative for below direction (score is being reduced)"
        )
    }
})

test_that("gseq.pwm_edits direction=below: score_after = score_before + sum(gains)", {
    pssm <- matrix(
        c(
            0.97, 0.01, 0.01, 0.01,
            0.01, 0.97, 0.01, 0.01,
            0.01, 0.01, 0.97, 0.01,
            0.01, 0.01, 0.01, 0.97
        ),
        nrow = 4, byrow = TRUE,
        dimnames = list(NULL, c("A", "C", "G", "T"))
    )

    result <- gseq.pwm_edits("ACGT", pssm,
        score.thresh = -1.0,
        prior = 0, bidirect = FALSE, direction = "below"
    )

    sub_rows <- result[!is.na(result$edit_type) & result$edit_type == "sub", ]
    if (nrow(sub_rows) > 0) {
        total_gain <- sum(sub_rows$gain)
        expect_equal(result$score_after[1], result$score_before[1] + total_gain,
            tolerance = 1e-3,
            info = "score_after must equal score_before + sum(gains) for substitutions"
        )
    }
})

test_that("gseq.pwm_edits direction=below with multiple sequences", {
    pssm <- matrix(
        c(
            0.97, 0.01, 0.01, 0.01,
            0.01, 0.97, 0.01, 0.01
        ),
        nrow = 2, byrow = TRUE,
        dimnames = list(NULL, c("A", "C", "G", "T"))
    )

    # "AC" is perfect match (above threshold), "TT" is poor match (likely below)
    result <- gseq.pwm_edits(c("AC", "TT"), pssm,
        score.thresh = -1.0,
        prior = 0, bidirect = FALSE, direction = "below"
    )

    # "AC" should need edits (score is above threshold)
    ac_rows <- result[result$seq_idx == 1, ]
    expect_true(any(ac_rows$n_edits > 0),
        info = "AC should need edits to bring score below threshold"
    )

    # "TT" should need 0 edits (score is already below threshold)
    tt_rows <- result[result$seq_idx == 2, ]
    expect_true(any(tt_rows$n_edits == 0),
        info = "TT should already be below threshold"
    )
})

test_that("gseq.pwm_edits direction=below matches manual_pwm_edit_distance_below", {
    pssm <- matrix(
        c(
            0.97, 0.01, 0.01, 0.01,
            0.01, 0.97, 0.01, 0.01,
            0.01, 0.01, 0.97, 0.01,
            0.01, 0.01, 0.01, 0.97
        ),
        nrow = 4, byrow = TRUE,
        dimnames = list(NULL, c("A", "C", "G", "T"))
    )

    seqs <- c("ACGT", "TCGT", "TTTT", "ACGA")
    threshold <- -5.0

    for (i in seq_along(seqs)) {
        result <- gseq.pwm_edits(seqs[i], pssm,
            score.thresh = threshold,
            prior = 0, bidirect = FALSE, direction = "below"
        )

        expected_edits <- manual_pwm_edit_distance_below(seqs[i], pssm, threshold, scan_all = TRUE)

        if (nrow(result) > 0) {
            actual_edits <- result$n_edits[1]
            if (is.na(expected_edits)) {
                # Unreachable — should not appear in output
                # (or seq was already below threshold -> 0 edits)
            } else {
                expect_equal(actual_edits, expected_edits,
                    info = paste0("Seq '", seqs[i], "': n_edits should match manual calculation")
                )
            }
        }
    }
})

test_that("gseq.pwm_edits direction=below with bidirectional scanning", {
    pssm <- matrix(
        c(
            0.97, 0.01, 0.01, 0.01,
            0.01, 0.97, 0.01, 0.01
        ),
        nrow = 2, byrow = TRUE,
        dimnames = list(NULL, c("A", "C", "G", "T"))
    )

    # "GT" forward strand scores ~-9.21 (already below -1.0 -> 0 edits)
    # "GT" revcomp "AC" scores ~-0.06 (above -1.0 -> needs edits)
    # With bidirect=TRUE, the code picks the strand needing fewest edits,
    # which is the forward strand with 0 edits.
    result <- gseq.pwm_edits("GT", pssm,
        score.thresh = -1.0,
        prior = 0, bidirect = TRUE, direction = "below"
    )

    expect_true(nrow(result) > 0)
    # Forward strand is already below threshold, so 0 edits on forward strand
    expect_equal(result$n_edits[1], 0L)
    expect_equal(result$strand[1], 1L)

    # Use a sequence where BOTH strands score above threshold
    # "AC" forward is ~-0.06 (above threshold), revcomp "GT" is ~-9.21 (below)
    # So again the best is 0 edits on the reverse strand.
    # Use a higher threshold to force both strands above it.
    result2 <- gseq.pwm_edits("AC", pssm,
        score.thresh = -0.01,
        prior = 0, bidirect = TRUE, direction = "below"
    )
    expect_true(nrow(result2) > 0)
    # Forward strand "AC" scores ~-0.06, reverse "GT" scores ~-9.21
    # Forward needs edits, reverse is already below -0.01
    # Best is reverse strand with 0 edits
    expect_equal(result2$n_edits[1], 0L)
})

test_that("gseq.pwm_edits direction=below with max_indels=1", {
    pssm <- matrix(
        c(
            0.97, 0.01, 0.01, 0.01,
            0.01, 0.97, 0.01, 0.01,
            0.01, 0.01, 0.97, 0.01,
            0.01, 0.01, 0.01, 0.97
        ),
        nrow = 4, byrow = TRUE,
        dimnames = list(NULL, c("A", "C", "G", "T"))
    )

    # "ACGT" is a perfect match — try to destroy it with indels allowed
    result <- gseq.pwm_edits("ACGT", pssm,
        score.thresh = -1.0,
        max_indels = 1L, prior = 0, bidirect = FALSE, direction = "below"
    )

    expect_true(is.data.frame(result))
    expect_true(nrow(result) > 0)
    expect_true("edit_type" %in% colnames(result))

    # Should find edits needed to bring score below threshold
    rows_with_edits <- result[result$n_edits > 0, ]
    expect_true(nrow(rows_with_edits) > 0)

    # score_after should be at or below threshold
    expect_true(all(rows_with_edits$score_after <= -1.0 + 1e-9),
        info = "score_after should be at or below threshold"
    )
})

test_that("gseq.pwm_edits direction=below: comparison with vtrack", {
    gdb.init("/home/aviezerl/hg38")
    remove_all_vtracks()

    pssm <- matrix(
        c(
            0.97, 0.01, 0.01, 0.01,
            0.01, 0.97, 0.01, 0.01,
            0.01, 0.01, 0.97, 0.01,
            0.01, 0.01, 0.01, 0.97
        ),
        nrow = 4, byrow = TRUE,
        dimnames = list(NULL, c("A", "C", "G", "T"))
    )

    # Find a region with a good match first
    test_interval <- gintervals("chr1", 10000, 10010)
    seq <- gseq.extract(test_interval)

    thresh <- -5.0

    # Use gseq.pwm_edits with direction=below
    result <- gseq.pwm_edits(seq, pssm,
        score.thresh = thresh,
        prior = 0, bidirect = FALSE, direction = "below"
    )

    # Create vtrack with direction=below
    gvtrack.create("edist_below", NULL, "pwm.edit_distance",
        pssm = pssm, score.thresh = thresh,
        bidirect = FALSE, prior = 0, extend = FALSE,
        direction = "below"
    )
    vtrack_result <- gextract("edist_below",
        intervals = test_interval, iterator = test_interval
    )

    # The n_edits from gseq.pwm_edits should match the vtrack result
    if (nrow(result) > 0 && !is.na(vtrack_result$edist_below)) {
        expect_equal(result$n_edits[1], as.integer(vtrack_result$edist_below),
            info = "n_edits from gseq.pwm_edits should match vtrack"
        )
    }

    remove_all_vtracks()
})

test_that("gseq.pwm_edits direction=below vs direction=above are complementary", {
    pssm <- matrix(
        c(
            0.97, 0.01, 0.01, 0.01,
            0.01, 0.97, 0.01, 0.01
        ),
        nrow = 2, byrow = TRUE,
        dimnames = list(NULL, c("A", "C", "G", "T"))
    )

    threshold <- -1.0

    # "AC" is a perfect match
    result_above <- gseq.pwm_edits("AC", pssm,
        score.thresh = threshold,
        prior = 0, bidirect = FALSE, direction = "above"
    )
    result_below <- gseq.pwm_edits("AC", pssm,
        score.thresh = threshold,
        prior = 0, bidirect = FALSE, direction = "below"
    )

    # "AC" is already above threshold: above should return 0 edits, below should need edits
    expect_equal(result_above$n_edits[1], 0L,
        info = "AC is already above threshold, direction=above should need 0 edits"
    )
    expect_true(result_below$n_edits[1] > 0,
        info = "AC is above threshold, direction=below should need edits"
    )

    # "TT" is well below threshold
    result_above_tt <- gseq.pwm_edits("TT", pssm,
        score.thresh = threshold,
        prior = 0, bidirect = FALSE, direction = "above"
    )
    result_below_tt <- gseq.pwm_edits("TT", pssm,
        score.thresh = threshold,
        prior = 0, bidirect = FALSE, direction = "below"
    )

    # "TT" is below threshold: above should need edits, below should return 0 edits
    expect_true(result_above_tt$n_edits[1] > 0,
        info = "TT is below threshold, direction=above should need edits"
    )
    expect_equal(result_below_tt$n_edits[1], 0L,
        info = "TT is already below threshold, direction=below should need 0 edits"
    )
})

test_that("gseq.pwm_edits direction=below with longer sequence picks best window", {
    pssm <- matrix(
        c(
            0.97, 0.01, 0.01, 0.01,
            0.01, 0.97, 0.01, 0.01,
            0.01, 0.01, 0.97, 0.01,
            0.01, 0.01, 0.01, 0.97
        ),
        nrow = 4, byrow = TRUE,
        dimnames = list(NULL, c("A", "C", "G", "T"))
    )

    # Sequence with a perfect ACGT match embedded in a longer context
    seq <- "TTTTACGTTTTT"
    threshold <- -1.0

    result <- gseq.pwm_edits(seq, pssm,
        score.thresh = threshold,
        prior = 0, bidirect = FALSE, direction = "below"
    )

    expect_true(nrow(result) > 0)
    # The best (highest-scoring) window should be found — the one containing ACGT
    # For direction=below, we want to destroy the BEST match
    rows_with_edits <- result[result$n_edits > 0, ]
    if (nrow(rows_with_edits) > 0) {
        # The window_seq should be the best-matching window
        expect_true(rows_with_edits$score_before[1] > threshold,
            info = "score_before should be above threshold (best window)"
        )
    }
})

test_that("gseq.pwm_edits direction=below with highly informative PSSM", {
    # PSSM with near-zero but non-zero probabilities (to avoid -Inf issues)
    pssm <- matrix(
        c(
            0.99, 0.003, 0.003, 0.004,
            0.003, 0.99, 0.003, 0.004
        ),
        nrow = 2, byrow = TRUE,
        dimnames = list(NULL, c("A", "C", "G", "T"))
    )

    # "AC" is a near-perfect match — score is 2*log(0.99) ~ -0.02
    # For direction=below, we need to make it score below threshold
    result <- gseq.pwm_edits("AC", pssm,
        score.thresh = -0.5,
        prior = 0, bidirect = FALSE, direction = "below"
    )

    expect_true(nrow(result) > 0)
    # Switching one position to a non-preferred base drops score dramatically
    rows_with_edits <- result[result$n_edits > 0, ]
    expect_true(nrow(rows_with_edits) > 0)
    expect_equal(rows_with_edits$n_edits[1], 1L,
        info = "One edit should suffice to bring score below threshold"
    )
    # Verify the score after is actually below threshold
    expect_true(rows_with_edits$score_after[1] <= -0.5 + 1e-9)
})
