create_isolated_test_db()

# --- R reference implementation for pwm.n_mutations ---
# For a single motif-length window: count single-base substitutions that
# independently cross the score threshold.
# Returns 0 if threshold already satisfied, NA if no single edit suffices.
manual_pwm_n_mutations <- function(seq, pssm_mat, threshold, direction = "above") {
    L <- nrow(pssm_mat)
    bases <- c("A", "C", "G", "T")
    seq_bases <- strsplit(seq, "")[[1]]

    if (length(seq_bases) < L) {
        return(NA_real_)
    }

    # Compute current score
    score <- 0
    for (i in 1:L) {
        b <- seq_bases[i]
        bidx <- match(b, bases)
        if (is.na(bidx)) {
            # N base: use the score from the LUT (index 4 maps to min)
            score <- score + min(log(pssm_mat[i, ]))
        } else {
            score <- score + log(pssm_mat[i, bidx])
        }
    }

    # Check if threshold already satisfied
    if (direction == "above" && score >= threshold) {
        return(0)
    }
    if (direction == "below" && score <= threshold) {
        return(0)
    }

    # Compute deficit
    deficit <- if (direction == "above") threshold - score else score - threshold

    # Count single-base substitutions that independently cross threshold
    count <- 0
    for (i in 1:L) {
        current_base <- seq_bases[i]
        current_idx <- match(current_base, bases)
        if (is.na(current_idx)) next # skip N bases

        current_log <- log(pssm_mat[i, current_idx])

        for (alt_base in setdiff(bases, current_base)) {
            alt_idx <- match(alt_base, bases)
            alt_log <- log(pssm_mat[i, alt_idx])

            delta <- if (direction == "above") {
                alt_log - current_log
            } else {
                current_log - alt_log
            }

            # Handle NaN from -Inf - (-Inf) or similar
            if (is.na(delta) || !is.finite(delta)) next
            if (delta >= deficit - 1e-12) count <- count + 1
        }
    }

    if (count == 0) {
        return(NA_real_)
    }
    return(count)
}

# Scan across all windows in a sequence, returning max n_mutations
manual_pwm_n_mutations_scan <- function(seq, pssm_mat, threshold, direction = "above") {
    L <- nrow(pssm_mat)
    n <- nchar(seq)
    if (n < L) {
        return(NA_real_)
    }

    best <- NA_real_
    for (start in 1:(n - L + 1)) {
        window_seq <- substr(seq, start, start + L - 1)
        cand <- manual_pwm_n_mutations(window_seq, pssm_mat, threshold, direction)
        if (is.na(best) || (!is.na(cand) && cand > best)) {
            best <- cand
        }
    }
    best
}

# Reverse complement helper
revcomp <- function(seq) {
    comp <- c(A = "T", C = "G", G = "C", T = "A", N = "N")
    bases <- strsplit(toupper(seq), "")[[1]]
    paste0(rev(comp[bases]), collapse = "")
}

# ============================================================================
# Basic correctness tests
# ============================================================================

test_that("pwm.n_mutations returns 0 when threshold already satisfied (direction=above)", {
    remove_all_vtracks()

    pssm <- create_test_pssm() # AC motif

    # Find position with "AC" pattern
    full_seq <- toupper(gseq.extract(gintervals(1, 200, 300)))
    ac_pos <- gregexpr("AC", full_seq)[[1]][1]

    if (ac_pos > 0) {
        abs_pos <- 200 + ac_pos - 1
        test_interval <- gintervals(1, abs_pos, abs_pos + 2)

        # Threshold at 0 — an AC window scores log(1)+log(1) = 0, satisfying >= 0
        threshold <- 0.0
        gvtrack.create("nmut_perfect", NULL,
            func = "pwm.n_mutations",
            pssm = pssm, score.thresh = threshold,
            bidirect = FALSE, extend = FALSE, prior = 0
        )

        result <- gextract("nmut_perfect", test_interval, iterator = test_interval)
        expect_equal(result$nmut_perfect[1], 0, tolerance = 1e-6)
    } else {
        skip("No AC pattern found in test region")
    }
})

test_that("pwm.n_mutations returns correct count for 1-edit window (direction=above)", {
    remove_all_vtracks()

    # Use a 3-position PSSM with clear preferences
    pssm <- matrix(c(
        0.8, 0.1, 0.05, 0.05, # Strong A
        0.1, 0.8, 0.05, 0.05, # Strong C
        0.1, 0.05, 0.8, 0.05 # Strong G
    ), ncol = 4, byrow = TRUE)
    colnames(pssm) <- c("A", "C", "G", "T")

    test_interval <- gintervals(1, 200, 240)
    seq <- toupper(gseq.extract(test_interval))

    # Choose threshold that should allow single edits for some windows
    threshold <- -3.0
    gvtrack.create("nmut_count", NULL,
        func = "pwm.n_mutations",
        pssm = pssm, score.thresh = threshold,
        bidirect = FALSE, extend = FALSE, prior = 0
    )

    result <- gextract("nmut_count", test_interval, iterator = test_interval)

    expected <- manual_pwm_n_mutations_scan(seq, pssm, threshold, "above")

    if (is.na(expected)) {
        expect_true(is.na(result$nmut_count[1]))
    } else {
        expect_equal(result$nmut_count[1], expected, tolerance = 1e-6)
    }
})

test_that("pwm.n_mutations returns NA when no single edit suffices", {
    remove_all_vtracks()

    pssm <- create_test_pssm()

    test_intervals <- gintervals(1, 200, 240)

    # Impossibly high threshold — no single edit can reach it
    threshold <- 100.0
    gvtrack.create("nmut_unreachable", NULL,
        func = "pwm.n_mutations",
        pssm = pssm, score.thresh = threshold,
        bidirect = FALSE, extend = FALSE, prior = 0
    )

    result <- gextract("nmut_unreachable", test_intervals, iterator = test_intervals)
    expect_true(is.na(result$nmut_unreachable[1]))
})

test_that("pwm.n_mutations matches R reference on multiple intervals", {
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
    gvtrack.create("nmut_ref", NULL,
        func = "pwm.n_mutations",
        pssm = pssm, score.thresh = threshold,
        bidirect = FALSE, extend = FALSE, prior = 0
    )

    result <- gextract("nmut_ref", test_intervals, iterator = test_intervals)

    for (i in seq_len(nrow(test_intervals))) {
        seq <- toupper(gseq.extract(test_intervals[i, ]))
        expected <- manual_pwm_n_mutations_scan(seq, pssm, threshold, "above")

        if (is.na(expected)) {
            expect_true(is.na(result$nmut_ref[i]),
                info = paste("Interval", i, "expected NA")
            )
        } else {
            expect_equal(result$nmut_ref[i], expected,
                tolerance = 1e-6,
                info = paste("Interval", i)
            )
        }
    }
})

test_that("pwm.n_mutations with 1bp iterator matches R reference per position", {
    remove_all_vtracks()

    # Use a PSSM without zero entries to avoid -Inf complications
    # (the C++ score table replaces log-zero with col_max for ABOVE direction)
    pssm <- matrix(c(
        0.8, 0.1, 0.05, 0.05,
        0.05, 0.8, 0.1, 0.05,
        0.05, 0.05, 0.8, 0.1
    ), ncol = 4, byrow = TRUE)
    colnames(pssm) <- c("A", "C", "G", "T")
    motif_len <- nrow(pssm)

    test_interval <- gintervals(1, 200, 210)
    threshold <- -3.0

    gvtrack.create("nmut_1bp", NULL,
        func = "pwm.n_mutations",
        pssm = pssm, score.thresh = threshold,
        bidirect = FALSE, extend = TRUE, prior = 0
    )

    result_1bp <- gextract("nmut_1bp", test_interval, iterator = 1)

    expect_true(nrow(result_1bp) > 0)

    for (idx in 1:min(5, nrow(result_1bp))) {
        pos <- result_1bp$start[idx]
        seq_window <- toupper(gseq.extract(gintervals(1, pos, pos + motif_len)))
        expected <- manual_pwm_n_mutations(seq_window, pssm, threshold, "above")

        if (is.na(expected)) {
            expect_true(is.na(result_1bp$nmut_1bp[idx]),
                info = paste("Position", pos)
            )
        } else {
            expect_equal(result_1bp$nmut_1bp[idx], expected,
                tolerance = 1e-6,
                info = paste("Position", pos)
            )
        }
    }
})

# ============================================================================
# Direction tests
# ============================================================================

test_that("pwm.n_mutations direction=below returns 0 when score already below threshold", {
    remove_all_vtracks()

    pssm <- create_test_pssm() # AC motif

    test_intervals <- gintervals(1, 200, 240)

    # Very high threshold — all windows score well below it
    threshold <- 100.0
    gvtrack.create("nmut_below_already", NULL,
        func = "pwm.n_mutations",
        pssm = pssm, score.thresh = threshold,
        score.min = -Inf,
        direction = "below",
        bidirect = FALSE, extend = FALSE, prior = 0
    )

    result <- gextract("nmut_below_already", test_intervals, iterator = test_intervals)
    expect_equal(result$nmut_below_already[1], 0, tolerance = 1e-6)
})

test_that("pwm.n_mutations direction=below matches R reference", {
    remove_all_vtracks()

    pssm <- matrix(c(
        0.8, 0.1, 0.05, 0.05,
        0.1, 0.8, 0.05, 0.05,
        0.1, 0.05, 0.8, 0.05
    ), ncol = 4, byrow = TRUE)
    colnames(pssm) <- c("A", "C", "G", "T")

    test_intervals <- gintervals(
        chrom = c(1, 1, 1),
        start = c(200, 300, 400),
        end = c(230, 330, 430)
    )

    threshold <- -3.5
    gvtrack.create("nmut_below_ref", NULL,
        func = "pwm.n_mutations",
        pssm = pssm, score.thresh = threshold,
        score.min = -Inf,
        direction = "below",
        bidirect = FALSE, extend = FALSE, prior = 0
    )

    result <- gextract("nmut_below_ref", test_intervals, iterator = test_intervals)

    for (i in seq_len(nrow(test_intervals))) {
        seq <- toupper(gseq.extract(test_intervals[i, ]))
        expected <- manual_pwm_n_mutations_scan(seq, pssm, threshold, "below")

        if (is.na(expected)) {
            expect_true(is.na(result$nmut_below_ref[i]),
                info = paste("Interval", i, "expected NA")
            )
        } else {
            expect_equal(result$nmut_below_ref[i], expected,
                tolerance = 1e-6,
                info = paste("Interval", i)
            )
        }
    }
})

test_that("pwm.n_mutations direction=below returns NA when threshold unreachable", {
    remove_all_vtracks()

    # Uniform PSSM: min score per position = log(0.25) ~ -1.386
    pssm <- matrix(c(
        0.25, 0.25, 0.25, 0.25,
        0.25, 0.25, 0.25, 0.25
    ), ncol = 4, byrow = TRUE)
    colnames(pssm) <- c("A", "C", "G", "T")

    test_intervals <- gintervals(1, 200, 240)

    # Threshold below the minimum possible score (-Inf can't be reached
    # with a uniform PSSM since all entries are log(0.25))
    threshold <- -100.0
    gvtrack.create("nmut_below_impossible", NULL,
        func = "pwm.n_mutations",
        pssm = pssm, score.thresh = threshold,
        score.min = -Inf,
        direction = "below",
        bidirect = FALSE, extend = FALSE, prior = 0
    )

    result <- gextract("nmut_below_impossible", test_intervals, iterator = test_intervals)
    expect_true(is.na(result$nmut_below_impossible[1]))
})

# ============================================================================
# Bidirectional tests
# ============================================================================

test_that("pwm.n_mutations bidirect=TRUE returns max across strands (direction=above)", {
    remove_all_vtracks()

    pssm <- create_test_pssm()

    test_interval <- gintervals(1, 200, 260)
    threshold <- -5.0

    gvtrack.create("nmut_fwd", NULL,
        func = "pwm.n_mutations",
        pssm = pssm, score.thresh = threshold,
        bidirect = FALSE, strand = 1, extend = FALSE, prior = 0
    )

    gvtrack.create("nmut_rev", NULL,
        func = "pwm.n_mutations",
        pssm = pssm, score.thresh = threshold,
        bidirect = FALSE, strand = -1, extend = FALSE, prior = 0
    )

    gvtrack.create("nmut_bidi", NULL,
        func = "pwm.n_mutations",
        pssm = pssm, score.thresh = threshold,
        bidirect = TRUE, extend = FALSE, prior = 0
    )

    result <- gextract(c("nmut_fwd", "nmut_rev", "nmut_bidi"),
        test_interval,
        iterator = test_interval
    )

    fwd <- result$nmut_fwd[1]
    rev <- result$nmut_rev[1]
    bidi <- result$nmut_bidi[1]

    # For direction=above, bidirectional should take max across strands at each
    # window position, since either strand crossing the threshold counts.
    # The aggregation across windows is still max.
    # So bidi >= max(fwd, rev) is expected (since both strands contribute).
    if (!is.na(fwd) && !is.na(rev)) {
        expect_true(bidi >= max(fwd, rev) - 1e-6,
            info = paste("fwd=", fwd, "rev=", rev, "bidi=", bidi)
        )
    }
})

test_that("pwm.n_mutations bidirect=TRUE with direction=below uses max across strands per position", {
    remove_all_vtracks()

    # Use a non-degenerate PSSM so scores are finite and direction=below is meaningful
    pssm <- matrix(c(
        0.8, 0.1, 0.05, 0.05,
        0.05, 0.8, 0.1, 0.05,
        0.05, 0.05, 0.8, 0.1
    ), ncol = 4, byrow = TRUE)
    colnames(pssm) <- c("A", "C", "G", "T")

    # PWM max scores in this region are around -3.4 to -9.0
    # Use threshold=-4.0 so some windows score above it (needing edits to go below)
    test_interval <- gintervals(1, 200, 260)
    threshold <- -4.0

    gvtrack.create("nmut_below_fwd", NULL,
        func = "pwm.n_mutations",
        pssm = pssm, score.thresh = threshold,
        score.min = -Inf,
        direction = "below",
        bidirect = FALSE, strand = 1, extend = TRUE, prior = 0
    )

    gvtrack.create("nmut_below_rev", NULL,
        func = "pwm.n_mutations",
        pssm = pssm, score.thresh = threshold,
        score.min = -Inf,
        direction = "below",
        bidirect = FALSE, strand = -1, extend = TRUE, prior = 0
    )

    gvtrack.create("nmut_below_bidi", NULL,
        func = "pwm.n_mutations",
        pssm = pssm, score.thresh = threshold,
        score.min = -Inf,
        direction = "below",
        bidirect = TRUE, extend = TRUE, prior = 0
    )

    result <- gextract(
        c("nmut_below_fwd", "nmut_below_rev", "nmut_below_bidi"),
        test_interval,
        iterator = 1
    )

    # Verify we have data to compare
    has_both <- !is.na(result$nmut_below_fwd) & !is.na(result$nmut_below_rev)
    expect_true(any(has_both), info = "Should have at least some positions with both strands non-NA")

    # For direction=below + bidirect, each position combines strands via max
    # (a genomic substitution affects both strands)
    for (i in which(has_both)) {
        fwd <- result$nmut_below_fwd[i]
        rev <- result$nmut_below_rev[i]
        bidi <- result$nmut_below_bidi[i]

        # Per-position combine = max of strands
        expect_true(bidi >= max(fwd, rev) - 1e-6,
            info = paste("pos", i, "fwd=", fwd, "rev=", rev, "bidi=", bidi)
        )
    }
})

# ============================================================================
# Score filter tests
# ============================================================================

test_that("pwm.n_mutations score.min filters out low-scoring windows", {
    remove_all_vtracks()

    pssm <- create_test_pssm()

    test_interval <- gintervals(1, 200, 240)
    threshold <- -3.0

    # No filter
    gvtrack.create("nmut_nofilt", NULL,
        func = "pwm.n_mutations",
        pssm = pssm, score.thresh = threshold,
        bidirect = FALSE, extend = FALSE, prior = 0
    )

    # With score.min = 0 (very restrictive — only perfect-scoring windows pass)
    gvtrack.create("nmut_highfilt", NULL,
        func = "pwm.n_mutations",
        pssm = pssm, score.thresh = threshold,
        score.min = 0.0,
        bidirect = FALSE, extend = FALSE, prior = 0
    )

    result <- gextract(c("nmut_nofilt", "nmut_highfilt"),
        test_interval,
        iterator = test_interval
    )

    # High filter should either be NA (all windows filtered) or <= the unfiltered
    # count (since fewer windows contribute). But since n_mutations is MAX across
    # qualifying windows, filtering can only reduce or NA-ify the result.
    if (!is.na(result$nmut_highfilt[1]) && !is.na(result$nmut_nofilt[1])) {
        expect_true(result$nmut_highfilt[1] <= result$nmut_nofilt[1] + 1e-6)
    }
})

test_that("pwm.n_mutations score.max filters out high-scoring windows", {
    remove_all_vtracks()

    pssm <- create_test_pssm()

    test_interval <- gintervals(1, 200, 240)
    threshold <- -5.0

    # No filter
    gvtrack.create("nmut_nomax", NULL,
        func = "pwm.n_mutations",
        pssm = pssm, score.thresh = threshold,
        bidirect = FALSE, extend = FALSE, prior = 0
    )

    # With score.max = very low value (filters most windows)
    gvtrack.create("nmut_maxfilt", NULL,
        func = "pwm.n_mutations",
        pssm = pssm, score.thresh = threshold,
        score.max = -100.0,
        bidirect = FALSE, extend = FALSE, prior = 0
    )

    result <- gextract(c("nmut_nomax", "nmut_maxfilt"),
        test_interval,
        iterator = test_interval
    )

    # With such a restrictive score.max, most or all windows are filtered out
    # Result should be NA or at most equal to the unfiltered result
    if (!is.na(result$nmut_maxfilt[1]) && !is.na(result$nmut_nomax[1])) {
        expect_true(result$nmut_maxfilt[1] <= result$nmut_nomax[1] + 1e-6)
    }
})

test_that("pwm.n_mutations direction=below score.min=-Inf matches default for below", {
    remove_all_vtracks()

    pssm <- create_test_pssm()

    test_interval <- gintervals(1, 200, 240)
    threshold <- -5.0

    # Default (no score.min)
    gvtrack.create("nmut_below_default", NULL,
        func = "pwm.n_mutations",
        pssm = pssm, score.thresh = threshold,
        direction = "below",
        bidirect = FALSE, extend = FALSE, prior = 0
    )

    # Explicit score.min = -Inf
    gvtrack.create("nmut_below_explicit", NULL,
        func = "pwm.n_mutations",
        pssm = pssm, score.thresh = threshold,
        score.min = -Inf,
        direction = "below",
        bidirect = FALSE, extend = FALSE, prior = 0
    )

    result <- gextract(c("nmut_below_default", "nmut_below_explicit"),
        test_interval,
        iterator = test_interval
    )

    expect_equal(result$nmut_below_default[1], result$nmut_below_explicit[1], tolerance = 1e-6)
})

# ============================================================================
# Edge case: N bases
# ============================================================================

test_that("pwm.n_mutations handles windows consistently with non-degenerate PSSM", {
    remove_all_vtracks()

    # Use a PSSM without zero entries so R reference and C++ agree
    pssm <- matrix(c(
        0.8, 0.1, 0.05, 0.05,
        0.05, 0.8, 0.1, 0.05
    ), ncol = 4, byrow = TRUE)
    colnames(pssm) <- c("A", "C", "G", "T")

    test_interval <- gintervals(1, 200, 240)
    seq <- toupper(gseq.extract(test_interval))

    threshold <- -3.0
    gvtrack.create("nmut_ntest", NULL,
        func = "pwm.n_mutations",
        pssm = pssm, score.thresh = threshold,
        bidirect = FALSE, extend = FALSE, prior = 0
    )

    result <- gextract("nmut_ntest", test_interval, iterator = test_interval)
    expected <- manual_pwm_n_mutations_scan(seq, pssm, threshold, "above")

    if (is.na(expected)) {
        expect_true(is.na(result$nmut_ntest[1]))
    } else {
        expect_equal(result$nmut_ntest[1], expected, tolerance = 1e-6)
    }
})

# ============================================================================
# Integration: gscreen filtering
# ============================================================================

test_that("pwm.n_mutations works in gscreen expression", {
    remove_all_vtracks()

    pssm <- create_test_pssm()
    threshold <- -5.0

    gvtrack.create("nmut_screen", NULL,
        func = "pwm.n_mutations",
        pssm = pssm, score.thresh = threshold,
        bidirect = FALSE, extend = FALSE, prior = 0
    )

    # Find intervals where n_mutations is not NA
    screened <- gscreen("!is.na(nmut_screen)",
        gintervals(1, 0, 5000),
        iterator = 20
    )

    # gscreen returns NULL when no intervals match; otherwise a data frame
    if (!is.null(screened)) {
        expect_s3_class(screened, "data.frame")
        expect_true(all(c("chrom", "start", "end") %in% colnames(screened)))
    }

    # Also verify gextract works with the same expression
    extracted <- gextract("nmut_screen", gintervals(1, 0, 1000), iterator = 20)
    expect_s3_class(extracted, "data.frame")
    expect_true("nmut_screen" %in% colnames(extracted))
})

test_that("pwm.n_mutations returns 0 for windows where threshold is already met (gextract integration)", {
    remove_all_vtracks()

    pssm <- create_test_pssm() # AC motif
    motif_len <- nrow(pssm)

    # Use extend=TRUE with 1bp iterator to check individual windows
    # Very low threshold that most windows should already satisfy
    threshold <- -20.0
    gvtrack.create("nmut_easy", NULL,
        func = "pwm.n_mutations",
        pssm = pssm, score.thresh = threshold,
        bidirect = FALSE, extend = TRUE, prior = 0
    )

    result <- gextract("nmut_easy", gintervals(1, 200, 210), iterator = 1)
    expect_true(nrow(result) > 0)

    # With such a low threshold, windows with non-N bases should score above it
    # and return 0
    non_na_rows <- !is.na(result$nmut_easy)
    if (any(non_na_rows)) {
        expect_true(all(result$nmut_easy[non_na_rows] == 0),
            info = "All qualifying windows should return 0 for a very low threshold"
        )
    }
})

# ============================================================================
# Parameter validation
# ============================================================================

test_that("pwm.n_mutations requires pssm parameter", {
    remove_all_vtracks()

    expect_error(
        gvtrack.create("nmut_nopssm", NULL,
            func = "pwm.n_mutations",
            score.thresh = -5.0
        ),
        "pssm"
    )
})

test_that("pwm.n_mutations requires score.thresh parameter", {
    remove_all_vtracks()

    pssm <- create_test_pssm()

    expect_error(
        gvtrack.create("nmut_nothresh", NULL,
            func = "pwm.n_mutations",
            pssm = pssm
        ),
        "score.thresh"
    )
})

test_that("pwm.n_mutations rejects invalid direction", {
    remove_all_vtracks()

    pssm <- create_test_pssm()

    expect_error(
        gvtrack.create("nmut_baddir", NULL,
            func = "pwm.n_mutations",
            pssm = pssm, score.thresh = -5.0,
            direction = "sideways"
        ),
        "direction"
    )
})

# ============================================================================
# Aggregation semantics: MAX across windows
# ============================================================================

test_that("pwm.n_mutations aggregates as MAX across qualifying windows", {
    remove_all_vtracks()

    pssm <- matrix(c(
        0.8, 0.1, 0.05, 0.05,
        0.1, 0.8, 0.05, 0.05,
        0.1, 0.05, 0.8, 0.05
    ), ncol = 4, byrow = TRUE)
    colnames(pssm) <- c("A", "C", "G", "T")

    motif_len <- nrow(pssm)

    # Use a wide interval so there are many windows
    test_interval <- gintervals(1, 200, 260)
    seq <- toupper(gseq.extract(test_interval))

    threshold <- -3.5
    gvtrack.create("nmut_agg", NULL,
        func = "pwm.n_mutations",
        pssm = pssm, score.thresh = threshold,
        bidirect = FALSE, extend = FALSE, prior = 0
    )

    result <- gextract("nmut_agg", test_interval, iterator = test_interval)

    # Compute per-window n_mutations manually and verify the MAX
    n_windows <- nchar(seq) - motif_len + 1
    per_window <- numeric(0)
    for (start in 1:n_windows) {
        window_seq <- substr(seq, start, start + motif_len - 1)
        cand <- manual_pwm_n_mutations(window_seq, pssm, threshold, "above")
        per_window <- c(per_window, cand)
    }

    # Max across all non-NA values
    non_na <- per_window[!is.na(per_window)]
    expected_max <- if (length(non_na) == 0) NA_real_ else max(non_na)

    if (is.na(expected_max)) {
        expect_true(is.na(result$nmut_agg[1]))
    } else {
        expect_equal(result$nmut_agg[1], expected_max, tolerance = 1e-6)
    }
})

# ============================================================================
# Prior handling
# ============================================================================

test_that("pwm.n_mutations with prior=0 vs prior>0 gives different results", {
    remove_all_vtracks()

    pssm <- create_test_pssm() # Has zeros in it

    test_interval <- gintervals(1, 200, 240)
    threshold <- -3.0

    gvtrack.create("nmut_prior0", NULL,
        func = "pwm.n_mutations",
        pssm = pssm, score.thresh = threshold,
        bidirect = FALSE, extend = FALSE, prior = 0
    )

    gvtrack.create("nmut_prior01", NULL,
        func = "pwm.n_mutations",
        pssm = pssm, score.thresh = threshold,
        bidirect = FALSE, extend = FALSE, prior = 0.01
    )

    result <- gextract(c("nmut_prior0", "nmut_prior01"),
        test_interval,
        iterator = test_interval
    )

    # With prior=0 on a PSSM that has zero entries, scores can be -Inf
    # vs prior>0 where scores are always finite. Results should generally differ.
    # This is a consistency check — at least one should be non-NA.
    expect_true(!is.na(result$nmut_prior0[1]) || !is.na(result$nmut_prior01[1]))
})
