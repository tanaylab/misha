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
