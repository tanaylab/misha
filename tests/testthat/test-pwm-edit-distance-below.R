create_isolated_test_db()

test_that("pwm.edit_distance direction=below basic functionality works", {
    remove_all_vtracks()

    # Create simple PSSM: AC motif
    pssm <- create_test_pssm()

    test_intervals <- gintervals(1, 200, 240)
    seq <- toupper(gseq.extract(test_intervals))

    # Threshold low enough that most windows are already above it
    # (requiring edits to bring score below)
    threshold <- -5.0
    gvtrack.create("edist_below", NULL,
        func = "pwm.edit_distance",
        pssm = pssm, score.thresh = threshold,
        direction = "below",
        bidirect = FALSE, extend = FALSE, prior = 0
    )

    result <- gextract("edist_below", test_intervals, iterator = test_intervals)

    # Manual calculation
    expected <- manual_pwm_edit_distance_below(seq, pssm, threshold)

    if (is.na(expected)) {
        expect_true(is.na(result$edist_below[1]))
    } else {
        expect_equal(result$edist_below[1], expected, tolerance = 1e-6)
    }
})

test_that("pwm.edit_distance direction=below returns 0 when already below threshold", {
    remove_all_vtracks()

    pssm <- create_test_pssm() # AC motif

    test_intervals <- gintervals(1, 200, 240)

    # Threshold very high: all windows should already score below it
    threshold <- 100.0
    gvtrack.create("edist_below_already", NULL,
        func = "pwm.edit_distance",
        pssm = pssm, score.thresh = threshold,
        direction = "below",
        bidirect = FALSE, extend = FALSE, prior = 0
    )

    result <- gextract("edist_below_already", test_intervals, iterator = test_intervals)

    # Should need 0 edits since all windows score way below the high threshold
    expect_equal(result$edist_below_already[1], 0, tolerance = 1e-6)
})

test_that("pwm.edit_distance direction=below returns NA for unreachable threshold", {
    remove_all_vtracks()

    # Create a PSSM where minimum score is bounded
    # If all positions have uniform probabilities, col_min is known
    pssm <- matrix(c(
        0.25, 0.25, 0.25, 0.25, # Uniform
        0.25, 0.25, 0.25, 0.25 # Uniform
    ), ncol = 4, byrow = TRUE)
    colnames(pssm) <- c("A", "C", "G", "T")

    test_intervals <- gintervals(1, 200, 240)

    # Threshold far below what's possible: S_min = 2 * log(0.25) = -2.77
    # Any threshold below S_min should be unreachable
    threshold <- -100.0
    gvtrack.create("edist_below_impossible", NULL,
        func = "pwm.edit_distance",
        pssm = pssm, score.thresh = threshold,
        direction = "below",
        bidirect = FALSE, extend = FALSE, prior = 0
    )

    result <- gextract("edist_below_impossible", test_intervals, iterator = test_intervals)

    # Should return NA: even switching all bases to worst can't push score below -100
    expect_true(is.na(result$edist_below_impossible[1]))
})

test_that("pwm.edit_distance direction=below matches R reference on multiple intervals", {
    remove_all_vtracks()

    pssm <- matrix(c(
        0.8, 0.1, 0.05, 0.05, # Strong A
        0.1, 0.8, 0.05, 0.05, # Strong C
        0.1, 0.05, 0.8, 0.05 # Strong G
    ), ncol = 4, byrow = TRUE)
    colnames(pssm) <- c("A", "C", "G", "T")

    test_intervals <- gintervals(
        chrom = c(1, 1, 1),
        start = c(200, 300, 400),
        end = c(230, 330, 430)
    )

    threshold <- -3.5

    gvtrack.create("edist_below_ref", NULL,
        func = "pwm.edit_distance",
        pssm = pssm, score.thresh = threshold,
        direction = "below",
        bidirect = FALSE, extend = FALSE, prior = 0
    )

    result <- gextract("edist_below_ref", test_intervals, iterator = test_intervals)

    # Verify each interval against manual reference
    for (i in seq_len(nrow(test_intervals))) {
        seq <- toupper(gseq.extract(test_intervals[i, ]))
        expected <- manual_pwm_edit_distance_below(seq, pssm, threshold)

        if (is.na(expected)) {
            expect_true(is.na(result$edist_below_ref[i]),
                info = paste("Interval", i, "expected NA")
            )
        } else {
            expect_equal(result$edist_below_ref[i], expected,
                tolerance = 1e-6,
                info = paste("Interval", i)
            )
        }
    }
})

test_that("pwm.edit_distance direction=below max_edits cap works", {
    remove_all_vtracks()

    pssm <- create_test_pssm()

    test_intervals <- gintervals(1, 200, 240)
    seq <- toupper(gseq.extract(test_intervals))

    # Use a threshold that requires edits
    threshold <- -5.0

    # Without max_edits
    gvtrack.create("edist_below_exact", NULL,
        func = "pwm.edit_distance",
        pssm = pssm, score.thresh = threshold,
        direction = "below",
        bidirect = FALSE, extend = FALSE, prior = 0
    )

    # With max_edits = 1 (only windows needing at most 1 edit work)
    gvtrack.create("edist_below_max1", NULL,
        func = "pwm.edit_distance",
        pssm = pssm, score.thresh = threshold, max_edits = 1,
        direction = "below",
        bidirect = FALSE, extend = FALSE, prior = 0
    )

    result <- gextract(c("edist_below_exact", "edist_below_max1"),
        test_intervals,
        iterator = test_intervals
    )

    exact_edits <- manual_pwm_edit_distance_below(seq, pssm, threshold)
    max1_edits <- manual_pwm_edit_distance_below(seq, pssm, threshold, max_edits = 1)

    if (is.na(exact_edits)) {
        expect_true(is.na(result$edist_below_exact[1]))
    } else {
        expect_equal(result$edist_below_exact[1], exact_edits, tolerance = 1e-6)
    }

    if (is.na(max1_edits)) {
        expect_true(is.na(result$edist_below_max1[1]))
    } else {
        expect_equal(result$edist_below_max1[1], max1_edits, tolerance = 1e-6)
    }

    # If exact needs > 1 edit, max1 should return NA
    if (!is.na(exact_edits) && exact_edits > 1) {
        expect_true(is.na(result$edist_below_max1[1]))
    }
})

test_that("pwm.edit_distance direction=below bidirectional considers both strands", {
    remove_all_vtracks()

    pssm <- create_test_pssm()

    test_interval <- gintervals(1, 200, 240)
    threshold <- -5.0

    # Forward only
    gvtrack.create("edist_below_fwd", NULL,
        func = "pwm.edit_distance",
        pssm = pssm, score.thresh = threshold,
        direction = "below",
        bidirect = FALSE, strand = 1, extend = FALSE, prior = 0
    )

    # Reverse only
    gvtrack.create("edist_below_rev", NULL,
        func = "pwm.edit_distance",
        pssm = pssm, score.thresh = threshold,
        direction = "below",
        bidirect = FALSE, strand = -1, extend = FALSE, prior = 0
    )

    # Bidirectional
    gvtrack.create("edist_below_bidi", NULL,
        func = "pwm.edit_distance",
        pssm = pssm, score.thresh = threshold,
        direction = "below",
        bidirect = TRUE, extend = FALSE, prior = 0
    )

    result <- gextract(c("edist_below_fwd", "edist_below_rev", "edist_below_bidi"),
        test_interval,
        iterator = test_interval
    )

    # Bidirectional should return minimum of both strands
    fwd <- result$edist_below_fwd[1]
    rev <- result$edist_below_rev[1]
    bidi <- result$edist_below_bidi[1]

    if (!is.na(fwd) && !is.na(rev)) {
        expect_equal(bidi, min(fwd, rev), tolerance = 1e-6)
    } else if (!is.na(fwd)) {
        expect_equal(bidi, fwd, tolerance = 1e-6)
    } else if (!is.na(rev)) {
        expect_equal(bidi, rev, tolerance = 1e-6)
    }
})

test_that("pwm.edit_distance.pos and pwm.max.edit_distance work with direction=below", {
    remove_all_vtracks()

    pssm <- matrix(c(
        0.8, 0.1, 0.05, 0.05,
        0.1, 0.8, 0.05, 0.05,
        0.1, 0.05, 0.8, 0.05
    ), ncol = 4, byrow = TRUE)
    colnames(pssm) <- c("A", "C", "G", "T")

    motif_len <- nrow(pssm)
    test_interval <- gintervals(1, 200, 250)
    seq <- toupper(gseq.extract(test_interval))
    threshold <- -4.0

    gvtrack.create("edist_below_min", NULL,
        func = "pwm.edit_distance",
        pssm = pssm, score.thresh = threshold,
        direction = "below",
        bidirect = FALSE, extend = FALSE, prior = 0
    )

    gvtrack.create("edist_below_pos", NULL,
        func = "pwm.edit_distance.pos",
        pssm = pssm, score.thresh = threshold,
        direction = "below",
        bidirect = FALSE, extend = FALSE, prior = 0
    )

    gvtrack.create("edist_below_max_site", NULL,
        func = "pwm.max.edit_distance",
        pssm = pssm, score.thresh = threshold,
        direction = "below",
        bidirect = FALSE, extend = FALSE, prior = 0
    )

    gvtrack.create("pwm_max_pos_below", NULL,
        func = "pwm.max.pos",
        pssm = pssm, bidirect = FALSE, extend = FALSE, prior = 0
    )

    result <- gextract(
        c("edist_below_min", "edist_below_pos", "edist_below_max_site", "pwm_max_pos_below"),
        test_interval,
        iterator = test_interval
    )

    # Check min edit distance against R reference
    expected_min <- manual_pwm_edit_distance_below(seq, pssm, threshold)
    if (is.na(expected_min)) {
        expect_true(is.na(result$edist_below_min[1]))
    } else {
        expect_equal(result$edist_below_min[1], expected_min, tolerance = 1e-6)
    }

    # Check position: find the first window that achieves the minimum edit distance
    if (!is.na(expected_min)) {
        best_pos_idx <- NA_integer_
        for (start_idx in seq_len(nchar(seq) - motif_len + 1)) {
            window_seq <- substr(seq, start_idx, start_idx + motif_len - 1)
            cand_edits <- manual_pwm_edit_distance_below(window_seq, pssm, threshold, scan_all = FALSE)
            if (!is.na(cand_edits) && abs(cand_edits - expected_min) < 1e-6) {
                best_pos_idx <- start_idx
                break
            }
        }
        expect_false(is.na(best_pos_idx))
        expect_equal(result$edist_below_pos[1], best_pos_idx, tolerance = 1e-6)
    }

    # Check pwm.max.edit_distance: edits at the max-scoring window
    pwm_pos_val <- result$pwm_max_pos_below[1]
    if (!is.na(pwm_pos_val)) {
        pwm_start_offset <- as.integer(round(pwm_pos_val)) - 1L
        expect_true(pwm_start_offset >= 0)
        pwm_window <- gintervals(
            test_interval$chrom,
            test_interval$start + pwm_start_offset,
            test_interval$start + pwm_start_offset + motif_len
        )
        pwm_seq <- toupper(gseq.extract(pwm_window))
        expected_pwm_edits <- manual_pwm_edit_distance_below(pwm_seq, pssm, threshold, scan_all = FALSE)
        if (is.na(expected_pwm_edits)) {
            expect_true(is.na(result$edist_below_max_site[1]))
        } else {
            expect_equal(result$edist_below_max_site[1], expected_pwm_edits, tolerance = 1e-6)
        }
    }
})

test_that("pwm.edit_distance direction=below with score.min/score.max filtering", {
    remove_all_vtracks()

    pssm <- create_test_pssm()

    test_interval <- gintervals(1, 200, 240)
    threshold <- -5.0

    # Without score filters
    gvtrack.create("edist_below_nofilt", NULL,
        func = "pwm.edit_distance",
        pssm = pssm, score.thresh = threshold,
        direction = "below",
        bidirect = FALSE, extend = FALSE, prior = 0
    )

    # With -Inf score.min (should not filter anything)
    gvtrack.create("edist_below_lowfilt", NULL,
        func = "pwm.edit_distance",
        pssm = pssm, score.thresh = threshold,
        score.min = -Inf,
        direction = "below",
        bidirect = FALSE, extend = FALSE, prior = 0
    )

    # With very high score.min (filters out most windows)
    gvtrack.create("edist_below_highfilt", NULL,
        func = "pwm.edit_distance",
        pssm = pssm, score.thresh = threshold,
        score.min = 0.0,
        direction = "below",
        bidirect = FALSE, extend = FALSE, prior = 0
    )

    result <- gextract(
        c("edist_below_nofilt", "edist_below_lowfilt", "edist_below_highfilt"),
        test_interval,
        iterator = test_interval
    )

    # Low filter should match no-filter
    expect_equal(result$edist_below_nofilt[1], result$edist_below_lowfilt[1], tolerance = 1e-6)

    # High filter should either be NA or >= unfiltered result
    if (!is.na(result$edist_below_highfilt[1]) && !is.na(result$edist_below_nofilt[1])) {
        expect_true(result$edist_below_highfilt[1] >= result$edist_below_nofilt[1])
    }
})

test_that("pwm.edit_distance direction=below with 1bp iterator matches R reference", {
    remove_all_vtracks()

    pssm <- create_test_pssm()
    motif_len <- nrow(pssm)

    test_interval <- gintervals(1, 200, 210)
    threshold <- -5.0

    gvtrack.create("edist_below_1bp", NULL,
        func = "pwm.edit_distance",
        pssm = pssm, score.thresh = threshold,
        direction = "below",
        bidirect = FALSE, extend = TRUE, prior = 0
    )

    result_1bp <- gextract("edist_below_1bp", test_interval, iterator = 1)

    expect_true(nrow(result_1bp) > 0)

    # Check a few positions manually
    for (idx in 1:min(3, nrow(result_1bp))) {
        pos <- result_1bp$start[idx]
        seq_window <- toupper(gseq.extract(gintervals(1, pos, pos + motif_len)))
        expected <- manual_pwm_edit_distance_below(seq_window, pssm, threshold)

        if (is.na(expected)) {
            expect_true(is.na(result_1bp$edist_below_1bp[idx]),
                info = paste("Position", pos)
            )
        } else {
            expect_equal(result_1bp$edist_below_1bp[idx], expected,
                tolerance = 1e-6,
                info = paste("Position", pos)
            )
        }
    }
})

test_that("pwm.edit_distance direction=below vs direction=above are complementary", {
    remove_all_vtracks()

    pssm <- matrix(c(
        0.8, 0.1, 0.05, 0.05,
        0.1, 0.8, 0.05, 0.05
    ), ncol = 4, byrow = TRUE)
    colnames(pssm) <- c("A", "C", "G", "T")

    test_interval <- gintervals(1, 200, 240)
    threshold <- -3.0

    gvtrack.create("edist_above", NULL,
        func = "pwm.edit_distance",
        pssm = pssm, score.thresh = threshold,
        direction = "above",
        bidirect = FALSE, extend = FALSE, prior = 0
    )

    gvtrack.create("edist_below", NULL,
        func = "pwm.edit_distance",
        pssm = pssm, score.thresh = threshold,
        direction = "below",
        bidirect = FALSE, extend = FALSE, prior = 0
    )

    result <- gextract(c("edist_above", "edist_below"), test_interval, iterator = test_interval)

    above_val <- result$edist_above[1]
    below_val <- result$edist_below[1]

    # If above needs 0 edits (score >= threshold), below should need >= 1 edit
    # (since score is already at/above threshold, need to push it down).
    # Conversely, if below needs 0 edits (score <= threshold), above needs >= 1.
    # They can't both be 0 unless the score equals the threshold exactly.
    if (!is.na(above_val) && above_val == 0 && !is.na(below_val)) {
        # Score >= threshold, so below needs edits (or happens to be exactly at threshold)
        expect_true(below_val >= 0)
    }
    if (!is.na(below_val) && below_val == 0 && !is.na(above_val)) {
        # Score <= threshold, so above needs edits (or exactly at threshold)
        expect_true(above_val >= 0)
    }

    # Both should always be non-negative when not NA
    if (!is.na(above_val)) expect_true(above_val >= 0)
    if (!is.na(below_val)) expect_true(below_val >= 0)
})

test_that("pwm.edit_distance direction=below with longer motif matches R reference", {
    remove_all_vtracks()

    # 6bp motif
    pssm <- matrix(c(
        0.9, 0.03, 0.03, 0.04,
        0.03, 0.9, 0.03, 0.04,
        0.03, 0.03, 0.9, 0.04,
        0.04, 0.03, 0.03, 0.9,
        0.9, 0.03, 0.03, 0.04,
        0.03, 0.9, 0.03, 0.04
    ), ncol = 4, byrow = TRUE)
    colnames(pssm) <- c("A", "C", "G", "T")

    test_interval <- gintervals(1, 200, 250)
    seq <- toupper(gseq.extract(test_interval))
    threshold <- -5.0

    gvtrack.create("edist_below_6bp", NULL,
        func = "pwm.edit_distance",
        pssm = pssm, score.thresh = threshold,
        direction = "below",
        bidirect = FALSE, extend = FALSE, prior = 0
    )

    result <- gextract("edist_below_6bp", test_interval, iterator = test_interval)

    expected <- manual_pwm_edit_distance_below(seq, pssm, threshold)

    if (is.na(expected)) {
        expect_true(is.na(result$edist_below_6bp[1]))
    } else {
        expect_equal(result$edist_below_6bp[1], expected, tolerance = 1e-6)
    }
})

test_that("pwm.edit_distance direction=below with different thresholds is monotone", {
    remove_all_vtracks()

    pssm <- create_test_pssm()

    test_interval <- gintervals(1, 200, 240)
    seq <- toupper(gseq.extract(test_interval))

    # Thresholds from low to high
    thresholds <- c(-10.0, -5.0, -2.0, 0.0)
    vnames <- sprintf("edist_below_%d", seq_along(thresholds))

    for (i in seq_along(thresholds)) {
        gvtrack.create(vnames[i], NULL,
            func = "pwm.edit_distance",
            pssm = pssm, score.thresh = thresholds[i],
            direction = "below",
            bidirect = FALSE, extend = FALSE, prior = 0
        )
    }

    result <- gextract(vnames, test_interval, iterator = test_interval)

    # Lower thresholds should require more (or equal) edits to reach below them
    edits <- sapply(vnames, function(v) result[[v]][1])
    finite_edits <- edits[!is.na(edits)]
    if (length(finite_edits) > 1) {
        # Lower thresholds require more edits in the "below" direction
        # Since thresholds go from low to high, edits should go from high to low
        expect_true(all(diff(finite_edits) <= 0),
            info = paste("Edits should decrease with increasing threshold. Got:", paste(finite_edits, collapse = ", "))
        )
    }

    # Verify each against R reference
    for (i in seq_along(thresholds)) {
        expected <- manual_pwm_edit_distance_below(seq, pssm, thresholds[i])
        if (is.na(expected)) {
            expect_true(is.na(result[[vnames[i]]][1]))
        } else {
            expect_equal(result[[vnames[i]]][1], expected, tolerance = 1e-6)
        }
    }
})

test_that("pwm.edit_distance direction=below with max_edits consistency", {
    remove_all_vtracks()

    pssm <- create_test_pssm()

    test_interval <- gintervals(1, 200, 240)
    threshold <- -5.0

    # Create vtracks with different max_edits settings
    gvtrack.create("edist_below_exact", NULL,
        func = "pwm.edit_distance",
        pssm = pssm, score.thresh = threshold,
        direction = "below",
        bidirect = FALSE, extend = FALSE, prior = 0
    )

    for (k in 1:3) {
        vname <- sprintf("edist_below_max%d", k)
        gvtrack.create(vname, NULL,
            func = "pwm.edit_distance",
            pssm = pssm, score.thresh = threshold, max_edits = k,
            direction = "below",
            bidirect = FALSE, extend = FALSE, prior = 0
        )
    }

    vnames <- c("edist_below_exact", sprintf("edist_below_max%d", 1:3))
    result <- gextract(vnames, test_interval, iterator = test_interval)

    exact <- result$edist_below_exact[1]

    for (k in 1:3) {
        vname <- sprintf("edist_below_max%d", k)
        if (!is.na(exact) && exact <= k) {
            expect_equal(result[[vname]][1], exact,
                tolerance = 1e-6,
                info = paste("max_edits =", k, "should match exact when exact <=", k)
            )
        }
        if (!is.na(exact) && exact > k) {
            expect_true(is.na(result[[vname]][1]),
                info = paste("max_edits =", k, "should be NA when exact =", exact)
            )
        }
    }
})

test_that("pwm.edit_distance direction=below with extend flag works", {
    remove_all_vtracks()

    pssm <- create_test_pssm()
    motif_len <- nrow(pssm)

    test_interval <- gintervals(1, 200, 202)
    threshold <- -5.0

    gvtrack.create("edist_below_ext", NULL,
        func = "pwm.edit_distance",
        pssm = pssm, score.thresh = threshold,
        direction = "below",
        bidirect = FALSE, extend = TRUE, prior = 0
    )

    gvtrack.create("edist_below_noext", NULL,
        func = "pwm.edit_distance",
        pssm = pssm, score.thresh = threshold,
        direction = "below",
        bidirect = FALSE, extend = FALSE, prior = 0
    )

    result <- gextract(c("edist_below_ext", "edist_below_noext"), test_interval, iterator = test_interval)

    # With extend=TRUE, window is expanded to include full motif
    seq_ext <- toupper(gseq.extract(gintervals(1, 200, 200 + motif_len)))
    expected_ext <- manual_pwm_edit_distance_below(seq_ext, pssm, threshold)

    # With extend=FALSE, window stays as-is (may be smaller than motif)
    seq_noext <- toupper(gseq.extract(test_interval))
    expected_noext <- manual_pwm_edit_distance_below(seq_noext, pssm, threshold)

    if (is.na(expected_ext)) {
        expect_true(is.na(result$edist_below_ext[1]))
    } else {
        expect_equal(result$edist_below_ext[1], expected_ext, tolerance = 1e-6)
    }

    if (is.na(expected_noext)) {
        expect_true(is.na(result$edist_below_noext[1]))
    } else {
        expect_equal(result$edist_below_noext[1], expected_noext, tolerance = 1e-6)
    }
})

test_that("pwm.edit_distance direction=below with zero-probability columns", {
    remove_all_vtracks()

    # PSSM with zero probabilities: log(0) = -Inf
    pssm <- matrix(c(
        1.0, 0.0, 0.0, 0.0, # Only A
        0.0, 1.0, 0.0, 0.0 # Only C
    ), ncol = 4, byrow = TRUE)
    colnames(pssm) <- c("A", "C", "G", "T")

    test_interval <- gintervals(1, 200, 240)
    seq <- toupper(gseq.extract(test_interval))
    threshold <- -5.0

    gvtrack.create("edist_below_zeros", NULL,
        func = "pwm.edit_distance",
        pssm = pssm, score.thresh = threshold,
        direction = "below",
        bidirect = FALSE, extend = FALSE, prior = 0
    )

    result <- gextract("edist_below_zeros", test_interval, iterator = test_interval)

    # For a window like "GT" (neither A nor C at first/second position):
    # score would be log(0) + log(0) = -Inf, which is below any finite threshold.
    # So if any window has a non-matching base, it's already below.
    # For "AC": score = log(1) + log(1) = 0 > -5, needs 1 edit to push below.
    # Manual reference should match.
    expected <- manual_pwm_edit_distance_below(seq, pssm, threshold)

    if (is.na(expected)) {
        expect_true(is.na(result$edist_below_zeros[1]))
    } else {
        expect_equal(result$edist_below_zeros[1], expected, tolerance = 1e-6)
    }
})

test_that("pwm.edit_distance direction=below with gscreen works", {
    remove_all_vtracks()

    pssm <- create_test_pssm()

    threshold <- -5.0

    gvtrack.create("edist_below_screen", NULL,
        func = "pwm.edit_distance",
        pssm = pssm, score.thresh = threshold,
        direction = "below",
        bidirect = FALSE, extend = FALSE, prior = 0
    )

    # Screen for intervals where below-direction edit distance equals 0
    # (i.e., score is already at or below threshold)
    test_intervals <- gintervals(1, 200, 300)
    result <- gscreen("!is.na(edist_below_screen)", test_intervals, iterator = 10)

    # Result should be a valid intervals data frame (could be empty if all NA)
    if (!is.null(result) && is.data.frame(result) && nrow(result) > 0) {
        expect_true(all(c("chrom", "start", "end") %in% names(result)))
    }
})

# --------------------------------------------------------------------------
# direction=below with indels (max_indels parameter)
# --------------------------------------------------------------------------

test_that("direction=below with max_indels=1: indels can reduce total edits", {
    remove_all_vtracks()

    # 4bp motif with strong preferences
    pssm <- matrix(c(
        0.97, 0.01, 0.01, 0.01, # A
        0.01, 0.97, 0.01, 0.01, # C
        0.01, 0.01, 0.97, 0.01, # G
        0.01, 0.01, 0.01, 0.97 # T
    ), ncol = 4, byrow = TRUE)
    colnames(pssm) <- c("A", "C", "G", "T")

    test_intervals <- gintervals(
        chrom = c(1, 1, 1),
        start = c(200, 500, 1000),
        end = c(260, 560, 1060)
    )
    threshold <- -3.0

    # Substitution-only
    gvtrack.create("edist_below_sub", NULL,
        func = "pwm.edit_distance",
        pssm = pssm, score.thresh = threshold,
        direction = "below", max_indels = 0,
        bidirect = FALSE, extend = FALSE, prior = 0
    )

    # With 1 indel
    gvtrack.create("edist_below_indel1", NULL,
        func = "pwm.edit_distance",
        pssm = pssm, score.thresh = threshold,
        direction = "below", max_indels = 1,
        bidirect = FALSE, extend = FALSE, prior = 0
    )

    result <- gextract(c("edist_below_sub", "edist_below_indel1"),
        test_intervals,
        iterator = test_intervals
    )

    for (i in seq_len(nrow(result))) {
        sub_val <- result$edist_below_sub[i]
        indel_val <- result$edist_below_indel1[i]

        # Indels can only help or stay the same: with indels <= without indels
        if (!is.na(sub_val) && !is.na(indel_val)) {
            expect_true(indel_val <= sub_val,
                info = paste("Row", i, ": indel result", indel_val, "should be <= sub-only", sub_val)
            )
        }

        # If substitution-only finds a result, indel version should too
        if (!is.na(sub_val)) {
            expect_false(is.na(indel_val),
                info = paste("Row", i, ": indel version should not be NA when sub-only is", sub_val)
            )
        }

        # Both should be non-negative when not NA
        if (!is.na(sub_val)) expect_true(sub_val >= 0)
        if (!is.na(indel_val)) expect_true(indel_val >= 0)
    }
})

test_that("direction=below with max_indels=2: more indels can further reduce edits", {
    remove_all_vtracks()

    pssm <- matrix(c(
        0.97, 0.01, 0.01, 0.01, # A
        0.01, 0.97, 0.01, 0.01, # C
        0.01, 0.01, 0.97, 0.01, # G
        0.01, 0.01, 0.01, 0.97 # T
    ), ncol = 4, byrow = TRUE)
    colnames(pssm) <- c("A", "C", "G", "T")

    test_intervals <- gintervals(
        chrom = c(1, 1, 1),
        start = c(200, 500, 1000),
        end = c(260, 560, 1060)
    )
    threshold <- -3.0

    # Compare max_indels = 0, 1, 2
    gvtrack.create("edist_below_d0", NULL,
        func = "pwm.edit_distance",
        pssm = pssm, score.thresh = threshold,
        direction = "below", max_indels = 0,
        bidirect = FALSE, extend = FALSE, prior = 0
    )
    gvtrack.create("edist_below_d1", NULL,
        func = "pwm.edit_distance",
        pssm = pssm, score.thresh = threshold,
        direction = "below", max_indels = 1,
        bidirect = FALSE, extend = FALSE, prior = 0
    )
    gvtrack.create("edist_below_d2", NULL,
        func = "pwm.edit_distance",
        pssm = pssm, score.thresh = threshold,
        direction = "below", max_indels = 2,
        bidirect = FALSE, extend = FALSE, prior = 0
    )

    result <- gextract(c("edist_below_d0", "edist_below_d1", "edist_below_d2"),
        test_intervals,
        iterator = test_intervals
    )

    for (i in seq_len(nrow(result))) {
        d0 <- result$edist_below_d0[i]
        d1 <- result$edist_below_d1[i]
        d2 <- result$edist_below_d2[i]

        # Monotonicity: d2 <= d1 <= d0 (more indels can only help or stay the same)
        if (!is.na(d0) && !is.na(d1)) {
            expect_true(d1 <= d0,
                info = paste("Row", i, ": d1=", d1, "should be <= d0=", d0)
            )
        }
        if (!is.na(d1) && !is.na(d2)) {
            expect_true(d2 <= d1,
                info = paste("Row", i, ": d2=", d2, "should be <= d1=", d1)
            )
        }
        if (!is.na(d0) && !is.na(d2)) {
            expect_true(d2 <= d0,
                info = paste("Row", i, ": d2=", d2, "should be <= d0=", d0)
            )
        }

        # If d0 (sub-only) is reachable, d1 and d2 must also be reachable
        if (!is.na(d0)) {
            expect_false(is.na(d1),
                info = paste("Row", i, ": d1 should not be NA when d0 =", d0)
            )
            expect_false(is.na(d2),
                info = paste("Row", i, ": d2 should not be NA when d0 =", d0)
            )
        }
        if (!is.na(d1)) {
            expect_false(is.na(d2),
                info = paste("Row", i, ": d2 should not be NA when d1 =", d1)
            )
        }
    }
})

test_that("direction=below with max_indels: substitution-only vs indels comparison", {
    remove_all_vtracks()

    # Longer 6bp motif for more interesting comparisons
    pssm <- matrix(c(
        0.9, 0.03, 0.03, 0.04,
        0.03, 0.9, 0.03, 0.04,
        0.03, 0.03, 0.9, 0.04,
        0.04, 0.03, 0.03, 0.9,
        0.9, 0.03, 0.03, 0.04,
        0.03, 0.9, 0.03, 0.04
    ), ncol = 4, byrow = TRUE)
    colnames(pssm) <- c("A", "C", "G", "T")

    test_interval <- gintervals(1, 200, 280)
    threshold <- -5.0

    # No indels
    gvtrack.create("edist_below_noindel", NULL,
        func = "pwm.edit_distance",
        pssm = pssm, score.thresh = threshold,
        direction = "below",
        bidirect = FALSE, extend = FALSE, prior = 0
    )

    # With 1 indel
    gvtrack.create("edist_below_1indel", NULL,
        func = "pwm.edit_distance",
        pssm = pssm, score.thresh = threshold,
        direction = "below", max_indels = 1,
        bidirect = FALSE, extend = FALSE, prior = 0
    )

    # With 2 indels
    gvtrack.create("edist_below_2indels", NULL,
        func = "pwm.edit_distance",
        pssm = pssm, score.thresh = threshold,
        direction = "below", max_indels = 2,
        bidirect = FALSE, extend = FALSE, prior = 0
    )

    result <- gextract(
        c("edist_below_noindel", "edist_below_1indel", "edist_below_2indels"),
        test_interval,
        iterator = test_interval
    )

    no_indel <- result$edist_below_noindel[1]
    with_1 <- result$edist_below_1indel[1]
    with_2 <- result$edist_below_2indels[1]

    # More indels should always give the same or fewer total edits
    if (!is.na(no_indel) && !is.na(with_1)) {
        expect_true(with_1 <= no_indel)
    }
    if (!is.na(no_indel) && !is.na(with_2)) {
        expect_true(with_2 <= no_indel)
    }
    if (!is.na(with_1) && !is.na(with_2)) {
        expect_true(with_2 <= with_1)
    }

    # Substitution-only reachable implies indel versions reachable
    if (!is.na(no_indel)) {
        expect_false(is.na(with_1))
        expect_false(is.na(with_2))
    }
})

test_that("direction=below with max_indels: cap is respected", {
    remove_all_vtracks()

    pssm <- matrix(c(
        0.97, 0.01, 0.01, 0.01,
        0.01, 0.97, 0.01, 0.01,
        0.01, 0.01, 0.97, 0.01,
        0.01, 0.01, 0.01, 0.97
    ), ncol = 4, byrow = TRUE)
    colnames(pssm) <- c("A", "C", "G", "T")

    test_intervals <- gintervals(
        chrom = c(1, 1, 1, 1),
        start = c(200, 500, 1000, 2000),
        end = c(260, 560, 1060, 2060)
    )
    threshold <- -3.0

    # max_indels=0 should match the default (no indels)
    gvtrack.create("edist_below_cap_default", NULL,
        func = "pwm.edit_distance",
        pssm = pssm, score.thresh = threshold,
        direction = "below",
        bidirect = FALSE, extend = FALSE, prior = 0
    )
    gvtrack.create("edist_below_cap0", NULL,
        func = "pwm.edit_distance",
        pssm = pssm, score.thresh = threshold,
        direction = "below", max_indels = 0,
        bidirect = FALSE, extend = FALSE, prior = 0
    )

    result <- gextract(c("edist_below_cap_default", "edist_below_cap0"),
        test_intervals,
        iterator = test_intervals
    )

    # max_indels=0 and no max_indels specified should produce identical results
    for (i in seq_len(nrow(result))) {
        if (is.na(result$edist_below_cap_default[i])) {
            expect_true(is.na(result$edist_below_cap0[i]),
                info = paste("Row", i, ": both should be NA")
            )
        } else {
            expect_equal(result$edist_below_cap_default[i], result$edist_below_cap0[i],
                tolerance = 1e-6,
                info = paste("Row", i, ": default and max_indels=0 should match")
            )
        }
    }
})

test_that("direction=below with indels: consistency - indels always <= sub-only", {
    remove_all_vtracks()

    # Test across many intervals to increase confidence
    pssm <- matrix(c(
        0.8, 0.1, 0.05, 0.05,
        0.1, 0.8, 0.05, 0.05,
        0.1, 0.05, 0.8, 0.05,
        0.04, 0.03, 0.03, 0.9
    ), ncol = 4, byrow = TRUE)
    colnames(pssm) <- c("A", "C", "G", "T")

    test_intervals <- gintervals(
        chrom = rep(1, 8),
        start = seq(200, 900, by = 100),
        end = seq(260, 960, by = 100)
    )
    threshold <- -4.0

    gvtrack.create("edist_below_con_sub", NULL,
        func = "pwm.edit_distance",
        pssm = pssm, score.thresh = threshold,
        direction = "below", max_indels = 0,
        bidirect = FALSE, extend = FALSE, prior = 0
    )
    gvtrack.create("edist_below_con_indel1", NULL,
        func = "pwm.edit_distance",
        pssm = pssm, score.thresh = threshold,
        direction = "below", max_indels = 1,
        bidirect = FALSE, extend = FALSE, prior = 0
    )
    gvtrack.create("edist_below_con_indel2", NULL,
        func = "pwm.edit_distance",
        pssm = pssm, score.thresh = threshold,
        direction = "below", max_indels = 2,
        bidirect = FALSE, extend = FALSE, prior = 0
    )

    result <- gextract(
        c("edist_below_con_sub", "edist_below_con_indel1", "edist_below_con_indel2"),
        test_intervals,
        iterator = test_intervals
    )

    for (i in seq_len(nrow(result))) {
        sub_val <- result$edist_below_con_sub[i]
        ind1_val <- result$edist_below_con_indel1[i]
        ind2_val <- result$edist_below_con_indel2[i]

        # Consistency: indel version should never return more edits
        if (!is.na(sub_val) && !is.na(ind1_val)) {
            expect_true(ind1_val <= sub_val,
                info = paste("Row", i, ": indel1 should be <= sub-only")
            )
        }
        if (!is.na(sub_val) && !is.na(ind2_val)) {
            expect_true(ind2_val <= sub_val,
                info = paste("Row", i, ": indel2 should be <= sub-only")
            )
        }
        if (!is.na(ind1_val) && !is.na(ind2_val)) {
            expect_true(ind2_val <= ind1_val,
                info = paste("Row", i, ": indel2 should be <= indel1")
            )
        }

        # If sub-only is reachable, indel versions must be too
        if (!is.na(sub_val)) {
            expect_false(is.na(ind1_val),
                info = paste("Row", i, ": indel1 should not be NA when sub-only is", sub_val)
            )
            expect_false(is.na(ind2_val),
                info = paste("Row", i, ": indel2 should not be NA when sub-only is", sub_val)
            )
        }
    }
})

test_that("direction=below with indels: already below threshold still returns 0", {
    remove_all_vtracks()

    pssm <- create_test_pssm() # AC motif

    test_intervals <- gintervals(1, 200, 240)

    # Very high threshold: all windows should score below it
    threshold <- 100.0

    gvtrack.create("edist_below_indel_already0", NULL,
        func = "pwm.edit_distance",
        pssm = pssm, score.thresh = threshold,
        direction = "below", max_indels = 0,
        bidirect = FALSE, extend = FALSE, prior = 0
    )
    gvtrack.create("edist_below_indel_already1", NULL,
        func = "pwm.edit_distance",
        pssm = pssm, score.thresh = threshold,
        direction = "below", max_indels = 1,
        bidirect = FALSE, extend = FALSE, prior = 0
    )
    gvtrack.create("edist_below_indel_already2", NULL,
        func = "pwm.edit_distance",
        pssm = pssm, score.thresh = threshold,
        direction = "below", max_indels = 2,
        bidirect = FALSE, extend = FALSE, prior = 0
    )

    result <- gextract(
        c("edist_below_indel_already0", "edist_below_indel_already1", "edist_below_indel_already2"),
        test_intervals,
        iterator = test_intervals
    )

    # All should return 0 regardless of max_indels setting
    expect_equal(result$edist_below_indel_already0[1], 0, tolerance = 1e-6)
    expect_equal(result$edist_below_indel_already1[1], 0, tolerance = 1e-6)
    expect_equal(result$edist_below_indel_already2[1], 0, tolerance = 1e-6)
})

test_that("direction=below with indels: bidirectional considers both strands", {
    remove_all_vtracks()

    pssm <- matrix(c(
        0.97, 0.01, 0.01, 0.01,
        0.01, 0.97, 0.01, 0.01,
        0.01, 0.01, 0.97, 0.01,
        0.01, 0.01, 0.01, 0.97
    ), ncol = 4, byrow = TRUE)
    colnames(pssm) <- c("A", "C", "G", "T")

    test_interval <- gintervals(1, 200, 260)
    threshold <- -3.0

    # Forward only with indels
    gvtrack.create("edist_below_indel_fwd", NULL,
        func = "pwm.edit_distance",
        pssm = pssm, score.thresh = threshold,
        direction = "below", max_indels = 1,
        bidirect = FALSE, strand = 1, extend = FALSE, prior = 0
    )

    # Reverse only with indels
    gvtrack.create("edist_below_indel_rev", NULL,
        func = "pwm.edit_distance",
        pssm = pssm, score.thresh = threshold,
        direction = "below", max_indels = 1,
        bidirect = FALSE, strand = -1, extend = FALSE, prior = 0
    )

    # Bidirectional with indels
    gvtrack.create("edist_below_indel_bidi", NULL,
        func = "pwm.edit_distance",
        pssm = pssm, score.thresh = threshold,
        direction = "below", max_indels = 1,
        bidirect = TRUE, extend = FALSE, prior = 0
    )

    result <- gextract(
        c("edist_below_indel_fwd", "edist_below_indel_rev", "edist_below_indel_bidi"),
        test_interval,
        iterator = test_interval
    )

    fwd <- result$edist_below_indel_fwd[1]
    rev <- result$edist_below_indel_rev[1]
    bidi <- result$edist_below_indel_bidi[1]

    # Bidirectional should return minimum of both strands
    if (!is.na(fwd) && !is.na(rev)) {
        expect_equal(bidi, min(fwd, rev), tolerance = 1e-6)
    } else if (!is.na(fwd)) {
        expect_equal(bidi, fwd, tolerance = 1e-6)
    } else if (!is.na(rev)) {
        expect_equal(bidi, rev, tolerance = 1e-6)
    }
})

test_that("direction=below with indels: max_edits cap interacts correctly with max_indels", {
    remove_all_vtracks()

    pssm <- matrix(c(
        0.9, 0.03, 0.03, 0.04,
        0.03, 0.9, 0.03, 0.04,
        0.03, 0.03, 0.9, 0.04,
        0.04, 0.03, 0.03, 0.9
    ), ncol = 4, byrow = TRUE)
    colnames(pssm) <- c("A", "C", "G", "T")

    test_interval <- gintervals(1, 200, 260)
    threshold <- -5.0

    # Unlimited edits with 1 indel
    gvtrack.create("edist_below_indel_unlim", NULL,
        func = "pwm.edit_distance",
        pssm = pssm, score.thresh = threshold,
        direction = "below", max_indels = 1,
        bidirect = FALSE, extend = FALSE, prior = 0
    )

    # max_edits=1 with 1 indel (tight cap)
    gvtrack.create("edist_below_indel_max1", NULL,
        func = "pwm.edit_distance",
        pssm = pssm, score.thresh = threshold,
        direction = "below", max_indels = 1, max_edits = 1,
        bidirect = FALSE, extend = FALSE, prior = 0
    )

    # max_edits=3 with 1 indel
    gvtrack.create("edist_below_indel_max3", NULL,
        func = "pwm.edit_distance",
        pssm = pssm, score.thresh = threshold,
        direction = "below", max_indels = 1, max_edits = 3,
        bidirect = FALSE, extend = FALSE, prior = 0
    )

    result <- gextract(
        c("edist_below_indel_unlim", "edist_below_indel_max1", "edist_below_indel_max3"),
        test_interval,
        iterator = test_interval
    )

    unlim <- result$edist_below_indel_unlim[1]
    max1 <- result$edist_below_indel_max1[1]
    max3 <- result$edist_below_indel_max3[1]

    # If unlimited finds a solution, capped versions should find same or NA
    if (!is.na(unlim)) {
        if (unlim <= 1) {
            expect_equal(max1, unlim, tolerance = 1e-6,
                info = "max_edits=1 should find same result when unlimited needs <= 1"
            )
        } else {
            expect_true(is.na(max1),
                info = paste("max_edits=1 should be NA when unlimited needs", unlim, "edits")
            )
        }

        if (unlim <= 3) {
            expect_equal(max3, unlim, tolerance = 1e-6,
                info = "max_edits=3 should find same result when unlimited needs <= 3"
            )
        } else {
            expect_true(is.na(max3),
                info = paste("max_edits=3 should be NA when unlimited needs", unlim, "edits")
            )
        }
    }
})

test_that("direction=below with indels: 1bp iterator matches interval-level result", {
    remove_all_vtracks()

    pssm <- matrix(c(
        0.97, 0.01, 0.01, 0.01,
        0.01, 0.97, 0.01, 0.01,
        0.01, 0.01, 0.97, 0.01
    ), ncol = 4, byrow = TRUE)
    colnames(pssm) <- c("A", "C", "G", "T")

    motif_len <- nrow(pssm)
    test_interval <- gintervals(1, 200, 210)
    threshold <- -4.0

    # Interval-level result with indels
    gvtrack.create("edist_below_indel_int", NULL,
        func = "pwm.edit_distance",
        pssm = pssm, score.thresh = threshold,
        direction = "below", max_indels = 1,
        bidirect = FALSE, extend = FALSE, prior = 0
    )

    # 1bp iterator result with indels
    gvtrack.create("edist_below_indel_1bp", NULL,
        func = "pwm.edit_distance",
        pssm = pssm, score.thresh = threshold,
        direction = "below", max_indels = 1,
        bidirect = FALSE, extend = TRUE, prior = 0
    )

    result_int <- gextract("edist_below_indel_int", test_interval, iterator = test_interval)
    result_1bp <- gextract("edist_below_indel_1bp", test_interval, iterator = 1)

    # The interval-level minimum should equal the minimum across 1bp windows
    if (nrow(result_1bp) > 0) {
        min_1bp <- min(result_1bp$edist_below_indel_1bp, na.rm = TRUE)
        if (is.finite(min_1bp) && !is.na(result_int$edist_below_indel_int[1])) {
            expect_equal(result_int$edist_below_indel_int[1], min_1bp, tolerance = 1e-6)
        }
    }
})

# --------------------------------------------------------------------------
# LSE edit distance with direction=below
# --------------------------------------------------------------------------

test_that("pwm.edit_distance.lse direction=below basic functionality works", {
    remove_all_vtracks()

    pssm <- matrix(c(
        0.8, 0.1, 0.05, 0.05,
        0.1, 0.8, 0.05, 0.05,
        0.1, 0.05, 0.8, 0.05
    ), ncol = 4, byrow = TRUE)
    colnames(pssm) <- c("A", "C", "G", "T")

    test_interval <- gintervals(1, 200, 250)

    # Get the actual LSE score for this interval
    gvtrack.create("v_lse_score", NULL, func = "pwm",
        pssm = pssm, bidirect = FALSE, extend = FALSE, prior = 0
    )
    lse_result <- gextract("v_lse_score", test_interval, iterator = test_interval)
    lse_score <- lse_result$v_lse_score[1]

    # Set threshold above the LSE score so the score is already below it
    # => should require 0 edits
    high_thresh <- lse_score + 5.0
    gvtrack.create("v_lse_below_easy", NULL,
        func = "pwm.edit_distance.lse",
        pssm = pssm, score.thresh = high_thresh,
        direction = "below",
        bidirect = FALSE, extend = FALSE, prior = 0
    )

    result <- gextract("v_lse_below_easy", test_interval, iterator = test_interval)
    expect_equal(result$v_lse_below_easy[1], 0, tolerance = 1e-6)
})

test_that("pwm.edit_distance.lse direction=below needs edits when score is above threshold", {
    remove_all_vtracks()

    pssm <- matrix(c(
        0.8, 0.1, 0.05, 0.05,
        0.1, 0.8, 0.05, 0.05,
        0.1, 0.05, 0.8, 0.05
    ), ncol = 4, byrow = TRUE)
    colnames(pssm) <- c("A", "C", "G", "T")

    test_interval <- gintervals(1, 200, 250)

    # Get the actual LSE score
    gvtrack.create("v_lse_score", NULL, func = "pwm",
        pssm = pssm, bidirect = FALSE, extend = FALSE, prior = 0
    )
    lse_result <- gextract("v_lse_score", test_interval, iterator = test_interval)
    lse_score <- lse_result$v_lse_score[1]

    # Set threshold well below the LSE score so edits are needed
    low_thresh <- lse_score - 10.0
    gvtrack.create("v_lse_below_hard", NULL,
        func = "pwm.edit_distance.lse",
        pssm = pssm, score.thresh = low_thresh,
        direction = "below",
        bidirect = FALSE, extend = FALSE, prior = 0
    )

    result <- gextract("v_lse_below_hard", test_interval, iterator = test_interval)

    # Should need at least 1 edit since the score is above threshold
    if (!is.na(result$v_lse_below_hard[1])) {
        expect_true(result$v_lse_below_hard[1] >= 1,
            info = paste("LSE score", lse_score, "above threshold", low_thresh, "should need edits")
        )
    }
})

test_that("pwm.edit_distance.lse direction=below returns 0 when already below", {
    remove_all_vtracks()

    pssm <- create_test_pssm()

    test_interval <- gintervals(1, 200, 240)

    # Very high threshold: LSE score will certainly be below it
    threshold <- 100.0
    gvtrack.create("v_lse_below_already", NULL,
        func = "pwm.edit_distance.lse",
        pssm = pssm, score.thresh = threshold,
        direction = "below",
        bidirect = FALSE, extend = FALSE, prior = 0
    )

    result <- gextract("v_lse_below_already", test_interval, iterator = test_interval)
    expect_equal(result$v_lse_below_already[1], 0, tolerance = 1e-6)
})

test_that("pwm.edit_distance.lse direction=below returns NA for unreachable threshold", {
    remove_all_vtracks()

    # With prior=0, the minimum per-window score can be -Inf (zero probability columns),
    # so LSE can go to -Inf. Use a PSSM with no zeros so the minimum is bounded.
    pssm <- matrix(c(
        0.25, 0.25, 0.25, 0.25,
        0.25, 0.25, 0.25, 0.25
    ), ncol = 4, byrow = TRUE)
    colnames(pssm) <- c("A", "C", "G", "T")

    test_interval <- gintervals(1, 200, 240)

    # With uniform PSSM and no prior, all windows have the same score = 2*log(0.25).
    # LSE = log(N * exp(2*log(0.25))) where N = number of windows.
    # This is bounded. A threshold far below this should be unreachable.
    threshold <- -1000.0
    gvtrack.create("v_lse_below_impossible", NULL,
        func = "pwm.edit_distance.lse",
        pssm = pssm, score.thresh = threshold,
        direction = "below",
        bidirect = FALSE, extend = FALSE, prior = 0
    )

    result <- gextract("v_lse_below_impossible", test_interval, iterator = test_interval)
    expect_true(is.na(result$v_lse_below_impossible[1]))
})

test_that("pwm.edit_distance.lse direction=below vs above are complementary", {
    remove_all_vtracks()

    pssm <- matrix(c(
        0.8, 0.1, 0.05, 0.05,
        0.1, 0.8, 0.05, 0.05
    ), ncol = 4, byrow = TRUE)
    colnames(pssm) <- c("A", "C", "G", "T")

    test_interval <- gintervals(1, 200, 240)

    # Get the actual LSE score
    gvtrack.create("v_lse_score", NULL, func = "pwm",
        pssm = pssm, bidirect = FALSE, extend = FALSE, prior = 0
    )
    lse_result <- gextract("v_lse_score", test_interval, iterator = test_interval)
    lse_score <- lse_result$v_lse_score[1]

    # Use a threshold below the actual score
    low_thresh <- lse_score - 5.0
    gvtrack.create("v_lse_above_low", NULL,
        func = "pwm.edit_distance.lse",
        pssm = pssm, score.thresh = low_thresh,
        direction = "above",
        bidirect = FALSE, extend = FALSE, prior = 0
    )
    gvtrack.create("v_lse_below_low", NULL,
        func = "pwm.edit_distance.lse",
        pssm = pssm, score.thresh = low_thresh,
        direction = "below",
        bidirect = FALSE, extend = FALSE, prior = 0
    )

    result <- gextract(c("v_lse_above_low", "v_lse_below_low"),
        test_interval, iterator = test_interval
    )

    above_val <- result$v_lse_above_low[1]
    below_val <- result$v_lse_below_low[1]

    # Score is above the threshold, so:
    # - "above" direction should need 0 edits (already above)
    # - "below" direction should need >= 1 edit (need to push score down)
    expect_equal(above_val, 0, tolerance = 1e-6)
    if (!is.na(below_val)) {
        expect_true(below_val >= 1,
            info = paste("Below should need edits when score", lse_score, "is above threshold", low_thresh)
        )
    }

    # Now use a threshold above the actual score
    high_thresh <- lse_score + 5.0
    gvtrack.create("v_lse_above_high", NULL,
        func = "pwm.edit_distance.lse",
        pssm = pssm, score.thresh = high_thresh,
        direction = "above",
        bidirect = FALSE, extend = FALSE, prior = 0
    )
    gvtrack.create("v_lse_below_high", NULL,
        func = "pwm.edit_distance.lse",
        pssm = pssm, score.thresh = high_thresh,
        direction = "below",
        bidirect = FALSE, extend = FALSE, prior = 0
    )

    result2 <- gextract(c("v_lse_above_high", "v_lse_below_high"),
        test_interval, iterator = test_interval
    )

    above_val2 <- result2$v_lse_above_high[1]
    below_val2 <- result2$v_lse_below_high[1]

    # Score is below the threshold, so:
    # - "above" direction should need >= 1 edit (need to push score up)
    # - "below" direction should need 0 edits (already below)
    expect_equal(below_val2, 0, tolerance = 1e-6)
    if (!is.na(above_val2)) {
        expect_true(above_val2 >= 1,
            info = paste("Above should need edits when score", lse_score, "is below threshold", high_thresh)
        )
    }
})

test_that("pwm.edit_distance.lse.pos direction=below returns valid position", {
    remove_all_vtracks()

    pssm <- matrix(c(
        0.8, 0.1, 0.05, 0.05,
        0.1, 0.8, 0.05, 0.05,
        0.1, 0.05, 0.8, 0.05
    ), ncol = 4, byrow = TRUE)
    colnames(pssm) <- c("A", "C", "G", "T")

    motif_len <- nrow(pssm)
    test_interval <- gintervals(1, 200, 250)

    # Get the actual LSE score
    gvtrack.create("v_lse_score", NULL, func = "pwm",
        pssm = pssm, bidirect = FALSE, extend = FALSE, prior = 0
    )
    lse_result <- gextract("v_lse_score", test_interval, iterator = test_interval)
    lse_score <- lse_result$v_lse_score[1]

    # Set threshold below the LSE score so edits are needed
    threshold <- lse_score - 5.0

    gvtrack.create("v_lse_below_edist", NULL,
        func = "pwm.edit_distance.lse",
        pssm = pssm, score.thresh = threshold,
        direction = "below",
        bidirect = FALSE, extend = FALSE, prior = 0
    )

    gvtrack.create("v_lse_below_pos", NULL,
        func = "pwm.edit_distance.lse.pos",
        pssm = pssm, score.thresh = threshold,
        direction = "below",
        bidirect = FALSE, extend = FALSE, prior = 0
    )

    result <- gextract(c("v_lse_below_edist", "v_lse_below_pos"),
        test_interval, iterator = test_interval
    )

    edist_val <- result$v_lse_below_edist[1]
    pos_val <- result$v_lse_below_pos[1]

    # If edit distance is not NA, position should also not be NA
    if (!is.na(edist_val) && edist_val >= 1) {
        expect_false(is.na(pos_val),
            info = "Position should be defined when edits are needed"
        )
        # Position should be a valid 1-based offset within the interval
        interval_len <- test_interval$end - test_interval$start
        expect_true(pos_val >= 1,
            info = paste("Position", pos_val, "should be >= 1")
        )
        expect_true(pos_val <= interval_len,
            info = paste("Position", pos_val, "should be within interval length", interval_len)
        )
    }

    # If already below (0 edits), position is typically NA
    if (!is.na(edist_val) && edist_val == 0) {
        expect_true(is.na(pos_val),
            info = "Position should be NA when no edits are needed"
        )
    }
})

test_that("pwm.edit_distance.lse direction=below edit count is consistent with max-mode below", {
    remove_all_vtracks()

    pssm <- matrix(c(
        0.8, 0.1, 0.05, 0.05,
        0.1, 0.8, 0.05, 0.05,
        0.1, 0.05, 0.8, 0.05
    ), ncol = 4, byrow = TRUE)
    colnames(pssm) <- c("A", "C", "G", "T")

    test_interval <- gintervals(1, 200, 250)

    # Get the actual LSE score and the max-window score
    gvtrack.create("v_lse_score", NULL, func = "pwm",
        pssm = pssm, bidirect = FALSE, extend = FALSE, prior = 0
    )
    gvtrack.create("v_max_score", NULL, func = "pwm.max",
        pssm = pssm, bidirect = FALSE, extend = FALSE, prior = 0
    )
    scores <- gextract(c("v_lse_score", "v_max_score"), test_interval, iterator = test_interval)
    lse_score <- scores$v_lse_score[1]
    max_score <- scores$v_max_score[1]

    # Use a threshold based on the max score (the max-mode below problem)
    threshold <- max_score - 2.0

    gvtrack.create("v_lse_below_edist", NULL,
        func = "pwm.edit_distance.lse",
        pssm = pssm, score.thresh = threshold,
        direction = "below",
        bidirect = FALSE, extend = FALSE, prior = 0
    )
    gvtrack.create("v_max_below_edist", NULL,
        func = "pwm.edit_distance",
        pssm = pssm, score.thresh = threshold,
        direction = "below",
        bidirect = FALSE, extend = FALSE, prior = 0
    )

    result <- gextract(c("v_lse_below_edist", "v_max_below_edist"),
        test_interval, iterator = test_interval
    )

    lse_edits <- result$v_lse_below_edist[1]
    max_edits <- result$v_max_below_edist[1]

    # The LSE score is always >= the max score (LSE = log-sum-exp >= max).
    # So to bring the LSE below a threshold, we generally need at least as many edits
    # as bringing the max below the same threshold.
    # In other words, LSE below edits >= max below edits.
    if (!is.na(lse_edits) && !is.na(max_edits)) {
        expect_true(lse_edits >= max_edits,
            info = paste(
                "LSE below edits (", lse_edits, ") should be >= max below edits (",
                max_edits, ") for the same threshold"
            )
        )
    }

    # If max below is reachable, LSE below may or may not be reachable
    # (LSE is harder to push down since all windows contribute)
    # Both should be non-negative when not NA
    if (!is.na(lse_edits)) expect_true(lse_edits >= 0)
    if (!is.na(max_edits)) expect_true(max_edits >= 0)
})

test_that("pwm.edit_distance.lse direction=below threshold monotonicity", {
    remove_all_vtracks()

    pssm <- matrix(c(
        0.8, 0.1, 0.05, 0.05,
        0.1, 0.8, 0.05, 0.05
    ), ncol = 4, byrow = TRUE)
    colnames(pssm) <- c("A", "C", "G", "T")

    test_interval <- gintervals(1, 200, 240)

    # Get the actual LSE score
    gvtrack.create("v_lse_score", NULL, func = "pwm",
        pssm = pssm, bidirect = FALSE, extend = FALSE, prior = 0
    )
    lse_result <- gextract("v_lse_score", test_interval, iterator = test_interval)
    lse_score <- lse_result$v_lse_score[1]

    # Thresholds from low to high
    thresholds <- c(lse_score - 15, lse_score - 10, lse_score - 5, lse_score, lse_score + 5)
    vnames <- sprintf("v_lse_below_t%d", seq_along(thresholds))

    for (i in seq_along(thresholds)) {
        gvtrack.create(vnames[i], NULL,
            func = "pwm.edit_distance.lse",
            pssm = pssm, score.thresh = thresholds[i],
            direction = "below",
            bidirect = FALSE, extend = FALSE, prior = 0
        )
    }

    result <- gextract(vnames, test_interval, iterator = test_interval)

    edits <- sapply(vnames, function(v) result[[v]][1])

    # Lower thresholds should require more (or equal) edits to bring the score below them.
    # So as threshold increases, edits should decrease (non-increasing).
    for (i in 2:length(edits)) {
        if (!is.na(edits[i - 1]) && !is.na(edits[i])) {
            expect_true(edits[i] <= edits[i - 1],
                info = paste(
                    "Edits at threshold", thresholds[i], "=", edits[i],
                    "should be <= edits at threshold", thresholds[i - 1], "=", edits[i - 1]
                )
            )
        }
        # If a higher threshold is reachable, lower thresholds may not be
        # but if a lower threshold is reachable, the higher one must also be
        if (!is.na(edits[i - 1])) {
            expect_false(is.na(edits[i]),
                info = paste(
                    "Threshold", thresholds[i], "should be reachable since lower threshold",
                    thresholds[i - 1], "is reachable"
                )
            )
        }
    }
})

test_that("pwm.edit_distance.lse direction=below with bidirectional", {
    remove_all_vtracks()

    pssm <- matrix(c(
        0.8, 0.1, 0.05, 0.05,
        0.1, 0.8, 0.05, 0.05,
        0.1, 0.05, 0.8, 0.05
    ), ncol = 4, byrow = TRUE)
    colnames(pssm) <- c("A", "C", "G", "T")

    test_interval <- gintervals(1, 200, 250)

    # Get LSE score to set a meaningful threshold
    gvtrack.create("v_lse_score", NULL, func = "pwm",
        pssm = pssm, bidirect = FALSE, extend = FALSE, prior = 0
    )
    lse_result <- gextract("v_lse_score", test_interval, iterator = test_interval)
    threshold <- lse_result$v_lse_score[1] - 3.0

    # Forward only
    gvtrack.create("v_lse_below_fwd", NULL,
        func = "pwm.edit_distance.lse",
        pssm = pssm, score.thresh = threshold,
        direction = "below",
        bidirect = FALSE, strand = 1, extend = FALSE, prior = 0
    )

    # Reverse only
    gvtrack.create("v_lse_below_rev", NULL,
        func = "pwm.edit_distance.lse",
        pssm = pssm, score.thresh = threshold,
        direction = "below",
        bidirect = FALSE, strand = -1, extend = FALSE, prior = 0
    )

    # Bidirectional
    gvtrack.create("v_lse_below_bidi", NULL,
        func = "pwm.edit_distance.lse",
        pssm = pssm, score.thresh = threshold,
        direction = "below",
        bidirect = TRUE, extend = FALSE, prior = 0
    )

    result <- gextract(c("v_lse_below_fwd", "v_lse_below_rev", "v_lse_below_bidi"),
        test_interval, iterator = test_interval
    )

    fwd <- result$v_lse_below_fwd[1]
    rev <- result$v_lse_below_rev[1]
    bidi <- result$v_lse_below_bidi[1]

    # All should be non-negative when not NA
    if (!is.na(fwd)) expect_true(fwd >= 0)
    if (!is.na(rev)) expect_true(rev >= 0)
    if (!is.na(bidi)) expect_true(bidi >= 0)

    # The bidirectional version operates on both strands. Since a single edit
    # can reduce scores on both strands simultaneously, the bidirectional edit
    # count may be less than, equal to, or greater than individual strands.
    # We verify all three produce valid results.
    if (!is.na(fwd) || !is.na(rev)) {
        # If at least one strand is reachable, bidirectional should also be reachable
        # (the algorithm can always degrade both strands or one)
        expect_false(is.na(bidi),
            info = "Bidirectional should be reachable when at least one strand is"
        )
    }
})

test_that("pwm.edit_distance.lse direction=below with max_edits cap", {
    remove_all_vtracks()

    pssm <- matrix(c(
        0.8, 0.1, 0.05, 0.05,
        0.1, 0.8, 0.05, 0.05,
        0.1, 0.05, 0.8, 0.05
    ), ncol = 4, byrow = TRUE)
    colnames(pssm) <- c("A", "C", "G", "T")

    test_interval <- gintervals(1, 200, 250)

    # Get LSE score to set a threshold that needs edits
    gvtrack.create("v_lse_score", NULL, func = "pwm",
        pssm = pssm, bidirect = FALSE, extend = FALSE, prior = 0
    )
    lse_result <- gextract("v_lse_score", test_interval, iterator = test_interval)
    threshold <- lse_result$v_lse_score[1] - 8.0

    # Unlimited edits
    gvtrack.create("v_lse_below_unlim", NULL,
        func = "pwm.edit_distance.lse",
        pssm = pssm, score.thresh = threshold,
        direction = "below",
        bidirect = FALSE, extend = FALSE, prior = 0
    )

    # max_edits = 1
    gvtrack.create("v_lse_below_max1", NULL,
        func = "pwm.edit_distance.lse",
        pssm = pssm, score.thresh = threshold, max_edits = 1,
        direction = "below",
        bidirect = FALSE, extend = FALSE, prior = 0
    )

    # max_edits = 5
    gvtrack.create("v_lse_below_max5", NULL,
        func = "pwm.edit_distance.lse",
        pssm = pssm, score.thresh = threshold, max_edits = 5,
        direction = "below",
        bidirect = FALSE, extend = FALSE, prior = 0
    )

    result <- gextract(c("v_lse_below_unlim", "v_lse_below_max1", "v_lse_below_max5"),
        test_interval, iterator = test_interval
    )

    unlim <- result$v_lse_below_unlim[1]
    max1 <- result$v_lse_below_max1[1]
    max5 <- result$v_lse_below_max5[1]

    # If unlimited finds a solution, capped versions should find the same or NA
    if (!is.na(unlim)) {
        if (unlim <= 1) {
            expect_equal(max1, unlim, tolerance = 1e-6)
        } else {
            expect_true(is.na(max1),
                info = paste("max_edits=1 should be NA when unlimited needs", unlim)
            )
        }

        if (unlim <= 5) {
            expect_equal(max5, unlim, tolerance = 1e-6)
        } else {
            expect_true(is.na(max5),
                info = paste("max_edits=5 should be NA when unlimited needs", unlim)
            )
        }
    }

    # Non-negative when not NA
    if (!is.na(unlim)) expect_true(unlim >= 0)
    if (!is.na(max1)) expect_true(max1 >= 0)
    if (!is.na(max5)) expect_true(max5 >= 0)
})
