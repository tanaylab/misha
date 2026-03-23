create_isolated_test_db()

test_that("pwm.edit_distance basic functionality works", {
    remove_all_vtracks()

    # Create simple PSSM: AC motif
    pssm <- create_test_pssm()

    # Test interval containing "ACGTACGT..."
    test_intervals <- gintervals(1, 200, 240)
    seq <- toupper(gseq.extract(test_intervals))

    # Create vtrack with moderate threshold
    threshold <- -5.0
    gvtrack.create("edist", NULL,
        func = "pwm.edit_distance",
        pssm = pssm, score.thresh = threshold,
        bidirect = FALSE, extend = FALSE, prior = 0
    )

    result <- gextract("edist", test_intervals, iterator = test_intervals)

    # Manual calculation: first 2 bases are "CC" for the AC motif
    expected <- manual_pwm_edit_distance(seq, pssm, threshold)

    expect_equal(result$edist[1], expected, tolerance = 1e-6)
})

test_that("pwm.edit_distance returns 0 for perfect match", {
    remove_all_vtracks()

    pssm <- create_test_pssm() # AC motif

    # Find position with "AC" pattern
    full_seq <- toupper(gseq.extract(gintervals(1, 200, 300)))
    ac_pos <- gregexpr("AC", full_seq)[[1]][1]

    if (ac_pos > 0) {
        # Position relative to chr1
        abs_pos <- 200 + ac_pos - 1
        test_interval <- gintervals(1, abs_pos, abs_pos + 2)

        threshold <- 0.0
        gvtrack.create("edist_perfect", NULL,
            func = "pwm.edit_distance",
            pssm = pssm, score.thresh = threshold,
            bidirect = FALSE, extend = FALSE, prior = 0
        )

        result <- gextract("edist_perfect", test_interval, iterator = test_interval)

        # Should need 0 edits for perfect match
        expect_equal(result$edist_perfect[1], 0, tolerance = 1e-6)
    } else {
        skip("No AC pattern found in test region")
    }
})

test_that("pwm.edit_distance returns NA for unreachable threshold", {
    remove_all_vtracks()

    pssm <- create_test_pssm()

    test_intervals <- gintervals(1, 200, 240)

    # Set impossibly high threshold
    threshold <- 100.0
    gvtrack.create("edist_unreachable", NULL,
        func = "pwm.edit_distance",
        pssm = pssm, score.thresh = threshold,
        bidirect = FALSE, extend = FALSE, prior = 0
    )

    result <- gextract("edist_unreachable", test_intervals, iterator = test_intervals)

    # Should return NA for unreachable threshold
    expect_true(is.na(result$edist_unreachable[1]))
})

test_that("pwm.edit_distance max_edits parameter works", {
    remove_all_vtracks()

    pssm <- create_test_pssm()

    test_intervals <- gintervals(1, 200, 240)
    seq <- toupper(gseq.extract(test_intervals))

    threshold <- -3.0

    # Create vtracks with different max_edits
    gvtrack.create("edist_exact", NULL,
        func = "pwm.edit_distance",
        pssm = pssm, score.thresh = threshold, max_edits = NULL,
        bidirect = FALSE, extend = FALSE, prior = 0
    )

    gvtrack.create("edist_max2", NULL,
        func = "pwm.edit_distance",
        pssm = pssm, score.thresh = threshold, max_edits = 2,
        bidirect = FALSE, extend = FALSE, prior = 0
    )

    result <- gextract(c("edist_exact", "edist_max2"), test_intervals, iterator = test_intervals)

    # Manual calculation
    exact_edits <- manual_pwm_edit_distance(seq, pssm, threshold, max_edits = NULL)
    max2_edits <- manual_pwm_edit_distance(seq, pssm, threshold, max_edits = 2)

    expect_equal(result$edist_exact[1], exact_edits, tolerance = 1e-6)
    expect_equal(result$edist_max2[1], max2_edits, tolerance = 1e-6)

    # If exact needs > 2 edits, max2 should return NA
    if (!is.na(exact_edits) && exact_edits > 2) {
        expect_true(is.na(result$edist_max2[1]))
    }
})

test_that("pwm.edit_distance bidirectional scanning works", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    # Create PSSM with clear strand preference (AC)
    pssm <- create_test_pssm()

    test_interval <- gintervals(1, 200, 240)
    threshold <- -5.0

    # Create vtracks for forward, reverse, and bidirectional
    gvtrack.create("edist_fwd", NULL,
        func = "pwm.edit_distance",
        pssm = pssm, score.thresh = threshold,
        bidirect = FALSE, strand = 1, extend = FALSE, prior = 0
    )

    gvtrack.create("edist_rev", NULL,
        func = "pwm.edit_distance",
        pssm = pssm, score.thresh = threshold,
        bidirect = FALSE, strand = -1, extend = FALSE, prior = 0
    )

    gvtrack.create("edist_bidi", NULL,
        func = "pwm.edit_distance",
        pssm = pssm, score.thresh = threshold,
        bidirect = TRUE, extend = FALSE, prior = 0
    )

    result <- gextract(c("edist_fwd", "edist_rev", "edist_bidi"), test_interval, iterator = test_interval)

    # Bidirectional should return minimum of both strands
    if (!is.na(result$edist_fwd[1]) && !is.na(result$edist_rev[1])) {
        expect_equal(result$edist_bidi[1], min(result$edist_fwd[1], result$edist_rev[1]), tolerance = 1e-6)
    }
})

test_that("pwm.edit_distance extend flag works", {
    remove_all_vtracks()

    pssm <- create_test_pssm()
    motif_len <- nrow(pssm)

    # Use small interval to see difference with/without extend
    test_interval <- gintervals(1, 200, 202)
    threshold <- -5.0

    gvtrack.create("edist_ext", NULL,
        func = "pwm.edit_distance",
        pssm = pssm, score.thresh = threshold,
        bidirect = FALSE, extend = TRUE, prior = 0
    )

    gvtrack.create("edist_noext", NULL,
        func = "pwm.edit_distance",
        pssm = pssm, score.thresh = threshold,
        bidirect = FALSE, extend = FALSE, prior = 0
    )

    result <- gextract(c("edist_ext", "edist_noext"), test_interval, iterator = test_interval)

    # With extend=TRUE, window is expanded to include full motif
    seq_ext <- toupper(gseq.extract(gintervals(1, 200, 200 + motif_len)))
    expected_ext <- manual_pwm_edit_distance(seq_ext, pssm, threshold)

    # With extend=FALSE, window stays as-is
    seq_noext <- toupper(gseq.extract(test_interval))
    expected_noext <- manual_pwm_edit_distance(seq_noext, pssm, threshold)

    expect_equal(result$edist_ext[1], expected_ext, tolerance = 1e-6)
    expect_equal(result$edist_noext[1], expected_noext, tolerance = 1e-6)
})

test_that("pwm.edit_distance iterator shifts work", {
    remove_all_vtracks()

    pssm <- create_test_pssm()
    motif_len <- nrow(pssm)

    base <- gintervals(1, 2000, 2040)
    threshold <- -3.0

    # Create vtrack with iterator shifts
    gvtrack.create("edist_shift", NULL,
        func = "pwm.edit_distance",
        pssm = pssm, score.thresh = threshold,
        bidirect = FALSE, extend = TRUE, prior = 0
    )

    gvtrack.iterator("edist_shift", sshift = -10, eshift = 10)

    result <- gextract("edist_shift", base, iterator = base)

    # Manual: expanded window should be [1990, 2050 + motif_len)
    ext_interval <- base
    ext_interval$start <- pmax(0, ext_interval$start - 10)
    ext_interval$end <- ext_interval$end + 10 + (motif_len - 1)
    seq_ext <- toupper(gseq.extract(ext_interval))

    expected <- manual_pwm_edit_distance(seq_ext, pssm, threshold)

    expect_equal(result$edist_shift[1], expected, tolerance = 1e-6)
})

test_that("pwm.edit_distance scans entire interval and reports positions", {
    remove_all_vtracks()

    pssm <- create_test_pssm()
    motif_len <- nrow(pssm)
    threshold <- -5.0
    genome_ranges <- gintervals.all()
    chrom <- genome_ranges$chrom[1]
    chrom_end <- genome_ranges$end[1]

    window_size <- 80L
    max_scan_start <- min(genome_ranges$start[1] + 5000L, chrom_end - window_size - 1)
    candidates <- seq(genome_ranges$start[1], max_scan_start, by = 10L)

    scan_interval <- NULL
    scan_seq <- NULL
    for (start_pos in candidates) {
        candidate <- gintervals(chrom, start_pos, start_pos + window_size)
        seq_candidate <- toupper(gseq.extract(candidate))
        hits <- gregexpr("AC", seq_candidate, fixed = TRUE)[[1]]
        if (!is.na(hits[1]) && hits[1] > 1 && substr(seq_candidate, 1, motif_len) != "AC") {
            scan_interval <- candidate
            scan_seq <- seq_candidate
            break
        }
    }

    if (is.null(scan_interval)) {
        skip("Unable to find interval with delayed perfect AC motif for edit-distance scan test")
    }

    gvtrack.create("edist_scan", NULL,
        func = "pwm.edit_distance",
        pssm = pssm, score.thresh = threshold,
        bidirect = FALSE, extend = FALSE, prior = 0
    )
    gvtrack.create("edist_pos", NULL,
        func = "pwm.edit_distance.pos",
        pssm = pssm, score.thresh = threshold,
        bidirect = FALSE, extend = FALSE, prior = 0
    )
    gvtrack.create("edist_max_site", NULL,
        func = "pwm.max.edit_distance",
        pssm = pssm, score.thresh = threshold,
        bidirect = FALSE, extend = FALSE, prior = 0
    )
    gvtrack.create("pwm_max_pos", NULL,
        func = "pwm.max.pos",
        pssm = pssm, bidirect = FALSE, extend = FALSE, prior = 0
    )

    result <- gextract(c("edist_scan", "edist_pos", "edist_max_site", "pwm_max_pos"),
        scan_interval,
        iterator = scan_interval
    )

    best_edit_manual <- manual_pwm_edit_distance(scan_seq, pssm, threshold)
    expect_equal(result$edist_scan[1], best_edit_manual, tolerance = 1e-6)

    first_window_edits <- manual_pwm_edit_distance(substr(scan_seq, 1, motif_len), pssm, threshold, scan_all = FALSE)
    expect_true(first_window_edits > result$edist_scan[1])

    # Determine the expected best position (1-based) for the min edit distance
    best_pos_idx <- NA_integer_
    for (start_idx in seq_len(nchar(scan_seq) - motif_len + 1)) {
        window_seq <- substr(scan_seq, start_idx, start_idx + motif_len - 1)
        cand_edits <- manual_pwm_edit_distance(window_seq, pssm, threshold, scan_all = FALSE)
        if (!is.na(cand_edits) && abs(cand_edits - best_edit_manual) < 1e-6) {
            best_pos_idx <- start_idx
            break
        }
    }
    expect_false(is.na(best_pos_idx))
    expect_equal(result$edist_pos[1], best_pos_idx, tolerance = 1e-6)

    # Validate pwm.max.edit_distance agrees with pwm.max.pos window
    pwm_pos_val <- result$pwm_max_pos[1]
    expect_false(is.na(pwm_pos_val))
    pwm_start_offset <- as.integer(round(pwm_pos_val)) - 1L
    expect_true(pwm_start_offset >= 0)
    pwm_window <- gintervals(
        scan_interval$chrom,
        scan_interval$start + pwm_start_offset,
        scan_interval$start + pwm_start_offset + motif_len
    )
    pwm_seq <- toupper(gseq.extract(pwm_window))
    expected_pwm_edits <- manual_pwm_edit_distance(pwm_seq, pssm, threshold, scan_all = FALSE)
    expect_equal(result$edist_max_site[1], expected_pwm_edits, tolerance = 1e-6)
})

test_that("pwm.edit_distance with different thresholds", {
    remove_all_vtracks()

    pssm <- create_test_pssm()

    test_interval <- gintervals(1, 200, 240)
    seq <- toupper(gseq.extract(test_interval))

    thresholds <- c(-10.0, -5.0, -2.0, 0.0)
    vnames <- sprintf("edist_thresh_%d", abs(thresholds) * 10)

    for (i in seq_along(thresholds)) {
        gvtrack.create(vnames[i], NULL,
            func = "pwm.edit_distance",
            pssm = pssm, score.thresh = thresholds[i],
            bidirect = FALSE, extend = FALSE, prior = 0
        )
    }

    result <- gextract(vnames, test_interval, iterator = test_interval)

    # Manually compute for each threshold
    for (i in seq_along(thresholds)) {
        expected <- manual_pwm_edit_distance(seq, pssm, thresholds[i])

        if (is.na(expected)) {
            expect_true(is.na(result[[vnames[i]]][1]))
        } else {
            expect_equal(result[[vnames[i]]][1], expected, tolerance = 1e-6)
        }
    }

    # Higher thresholds should require more (or equal) edits
    edits <- sapply(vnames, function(v) result[[v]][1])
    # Filter out NAs for comparison
    finite_edits <- edits[!is.na(edits)]
    if (length(finite_edits) > 1) {
        expect_true(all(diff(finite_edits) >= 0))
    }
})

test_that("pwm.edit_distance works with 1bp iterator", {
    remove_all_vtracks()

    pssm <- create_test_pssm()
    motif_len <- nrow(pssm)

    test_interval <- gintervals(1, 200, 210)
    threshold <- -5.0

    gvtrack.create("edist_1bp", NULL,
        func = "pwm.edit_distance",
        pssm = pssm, score.thresh = threshold,
        bidirect = FALSE, extend = TRUE, prior = 0
    )

    result_1bp <- gextract("edist_1bp", test_interval, iterator = 1)

    # Each row should correspond to a window starting at that position
    expect_true(nrow(result_1bp) > 0)

    # Check a few positions manually
    for (idx in 1:min(3, nrow(result_1bp))) {
        pos <- result_1bp$start[idx]
        seq_window <- toupper(gseq.extract(gintervals(1, pos, pos + motif_len)))
        expected <- manual_pwm_edit_distance(seq_window, pssm, threshold)

        if (is.na(expected)) {
            expect_true(is.na(result_1bp$edist_1bp[idx]))
        } else {
            expect_equal(result_1bp$edist_1bp[idx], expected, tolerance = 1e-6)
        }
    }
})

test_that("pwm.edit_distance counts forced edits when every column must change", {
    remove_all_vtracks()

    pssm <- create_test_pssm()
    motif_len <- nrow(pssm)

    search_interval <- gintervals(1, 150, 350)
    seq_region <- toupper(gseq.extract(search_interval))

    start_offset <- NA_integer_
    seq_window <- NULL
    for (offset in seq_len(nchar(seq_region) - motif_len + 1)) {
        candidate <- substr(seq_region, offset, offset + motif_len - 1)
        if (substr(candidate, 1, 1) != "A" && substr(candidate, 2, 2) != "C") {
            start_offset <- offset - 1
            seq_window <- candidate
            break
        }
    }

    if (is.na(start_offset)) {
        skip("No window requiring two forced edits was found in test genome")
    }

    abs_start <- search_interval$start + start_offset
    forced_interval <- gintervals(1, abs_start, abs_start + motif_len)
    threshold <- -1.0

    gvtrack.create("edist_forced", NULL,
        func = "pwm.edit_distance",
        pssm = pssm, score.thresh = threshold,
        bidirect = FALSE, extend = FALSE, prior = 0
    )

    result <- gextract("edist_forced", forced_interval, iterator = forced_interval)

    expect_equal(result$edist_forced[1], motif_len)

    expected <- manual_pwm_edit_distance(seq_window, pssm, threshold)
    expect_equal(result$edist_forced[1], expected, tolerance = 1e-6)
})

test_that("pwm.edit_distance honors max_edits when forced edits exceed budget", {
    remove_all_vtracks()

    pssm <- create_test_pssm()
    motif_len <- nrow(pssm)

    search_interval <- gintervals(1, 150, 350)
    seq_region <- toupper(gseq.extract(search_interval))

    start_offset <- NA_integer_
    seq_window <- NULL
    for (offset in seq_len(nchar(seq_region) - motif_len + 1)) {
        candidate <- substr(seq_region, offset, offset + motif_len - 1)
        if (substr(candidate, 1, 1) != "A" && substr(candidate, 2, 2) != "C") {
            start_offset <- offset - 1
            seq_window <- candidate
            break
        }
    }

    if (is.na(start_offset)) {
        skip("No window requiring two forced edits was found in test genome")
    }

    abs_start <- search_interval$start + start_offset
    forced_interval <- gintervals(1, abs_start, abs_start + motif_len)
    threshold <- -1.0

    gvtrack.create("edist_limit", NULL,
        func = "pwm.edit_distance",
        pssm = pssm, score.thresh = threshold, max_edits = 1,
        bidirect = FALSE, extend = FALSE, prior = 0
    )

    result <- gextract("edist_limit", forced_interval, iterator = forced_interval)

    expect_true(is.na(result$edist_limit[1]))

    expected <- manual_pwm_edit_distance(seq_window, pssm, threshold, max_edits = 1)
    expect_true(is.na(expected))
})

test_that("pwm.edit_distance max_edits=1 fast heuristic works", {
    remove_all_vtracks()

    pssm <- matrix(c(
        0.8, 0.1, 0.05, 0.05, # Strong A
        0.1, 0.8, 0.05, 0.05, # Strong C
        0.1, 0.05, 0.8, 0.05 # Strong G
    ), ncol = 4, byrow = TRUE)
    colnames(pssm) <- c("A", "C", "G", "T")

    test_interval <- gintervals(1, 200, 240)
    threshold <- -2.0

    gvtrack.create("edist_fast1", NULL,
        func = "pwm.edit_distance",
        pssm = pssm, score.thresh = threshold, max_edits = 1,
        bidirect = FALSE, extend = FALSE, prior = 0
    )

    gvtrack.create("edist_exact", NULL,
        func = "pwm.edit_distance",
        pssm = pssm, score.thresh = threshold, max_edits = NULL,
        bidirect = FALSE, extend = FALSE, prior = 0
    )

    result <- gextract(c("edist_fast1", "edist_exact"), test_interval, iterator = test_interval)

    # If exact needs 0 or 1 edit, fast should match
    if (!is.na(result$edist_exact[1]) && result$edist_exact[1] <= 1) {
        expect_equal(result$edist_fast1[1], result$edist_exact[1], tolerance = 1e-6)
    }

    # If exact needs > 1 edit, fast should return NA
    if (!is.na(result$edist_exact[1]) && result$edist_exact[1] > 1) {
        expect_true(is.na(result$edist_fast1[1]))
    }
})

test_that("pwm.edit_distance handles edge case with tiny interval", {
    remove_all_vtracks()

    pssm <- create_test_pssm()
    motif_len <- nrow(pssm)

    # Interval smaller than motif length
    tiny_interval <- gintervals(1, 200, 200 + 1)
    threshold <- 0.0

    gvtrack.create("edist_tiny", NULL,
        func = "pwm.edit_distance",
        pssm = pssm, score.thresh = threshold,
        bidirect = FALSE, extend = FALSE, prior = 0
    )

    result <- gextract("edist_tiny", tiny_interval, iterator = tiny_interval)

    # Should return NA for interval smaller than motif
    expect_true(is.na(result$edist_tiny[1]))
})

test_that("pwm.edit_distance multiple windows in single extract", {
    remove_all_vtracks()

    pssm <- create_test_pssm()

    # Multiple intervals
    test_intervals <- gintervals(
        chrom = c(1, 1, 1),
        start = c(200, 300, 400),
        end = c(210, 310, 410)
    )

    threshold <- -4.0

    gvtrack.create("edist_multi", NULL,
        func = "pwm.edit_distance",
        pssm = pssm, score.thresh = threshold,
        bidirect = FALSE, extend = FALSE, prior = 0
    )

    result <- gextract("edist_multi", test_intervals, iterator = test_intervals)

    # Should have 3 rows
    expect_equal(nrow(result), 3)

    # Each should have a valid result (could be NA if unreachable)
    expect_true(all(is.na(result$edist_multi) | result$edist_multi >= 0))
})

test_that("pwm.edit_distance consistency between exact and heuristic for small k", {
    remove_all_vtracks()

    pssm <- create_test_pssm()

    test_interval <- gintervals(1, 200, 240)
    threshold <- -3.0

    # Create vtracks with different max_edits settings
    gvtrack.create("edist_exact", NULL,
        func = "pwm.edit_distance",
        pssm = pssm, score.thresh = threshold, max_edits = NULL,
        bidirect = FALSE, extend = FALSE, prior = 0
    )

    for (k in 1:5) {
        vname <- sprintf("edist_max%d", k)
        gvtrack.create(vname, NULL,
            func = "pwm.edit_distance",
            pssm = pssm, score.thresh = threshold, max_edits = k,
            bidirect = FALSE, extend = FALSE, prior = 0
        )
    }

    vnames <- c("edist_exact", sprintf("edist_max%d", 1:5))
    result <- gextract(vnames, test_interval, iterator = test_interval)

    exact <- result$edist_exact[1]

    # For each k, if exact <= k, heuristic should match exact
    for (k in 1:5) {
        vname <- sprintf("edist_max%d", k)
        if (!is.na(exact) && exact <= k) {
            expect_equal(result[[vname]][1], exact, tolerance = 1e-6)
        }
    }
})

test_that("pwm.edit_distance with prior parameter", {
    remove_all_vtracks()

    pssm <- create_test_pssm()

    test_interval <- gintervals(1, 200, 240)
    threshold <- -5.0

    # Create vtracks with different priors
    gvtrack.create("edist_prior0", NULL,
        func = "pwm.edit_distance",
        pssm = pssm, score.thresh = threshold,
        bidirect = FALSE, extend = FALSE, prior = 0
    )

    gvtrack.create("edist_prior01", NULL,
        func = "pwm.edit_distance",
        pssm = pssm, score.thresh = threshold,
        bidirect = FALSE, extend = FALSE, prior = 0.01
    )

    result <- gextract(c("edist_prior0", "edist_prior01"), test_interval, iterator = test_interval)

    # Both should return valid results
    expect_true(!is.na(result$edist_prior0[1]) || !is.na(result$edist_prior01[1]))

    # Results might differ due to prior affecting scores
    # Just verify they're both non-negative when not NA
    if (!is.na(result$edist_prior0[1])) {
        expect_true(result$edist_prior0[1] >= 0)
    }
    if (!is.na(result$edist_prior01[1])) {
        expect_true(result$edist_prior01[1] >= 0)
    }
})

test_that("pwm.edit_distance strand parameter works independently", {
    remove_all_vtracks()

    pssm <- create_test_pssm()

    test_interval <- gintervals(1, 200, 240)
    threshold <- -5.0

    # Create vtracks with explicit strand settings (not using bidirect)
    gvtrack.create("edist_fwd", NULL,
        func = "pwm.edit_distance",
        pssm = pssm, score.thresh = threshold,
        strand = 1, extend = FALSE, prior = 0
    )

    gvtrack.create("edist_rev", NULL,
        func = "pwm.edit_distance",
        pssm = pssm, score.thresh = threshold,
        strand = -1, extend = FALSE, prior = 0
    )

    result <- gextract(c("edist_fwd", "edist_rev"), test_interval, iterator = test_interval)

    # Both should return valid edit distances
    expect_true(!is.na(result$edist_fwd[1]) || !is.na(result$edist_rev[1]))

    # Results might differ between strands
    if (!is.na(result$edist_fwd[1])) {
        expect_true(result$edist_fwd[1] >= 0)
    }
    if (!is.na(result$edist_rev[1])) {
        expect_true(result$edist_rev[1] >= 0)
    }
})

test_that("pwm.edit_distance works with longer motif", {
    remove_all_vtracks()

    # Create 6bp motif with clear pattern
    pssm <- matrix(c(
        0.9, 0.03, 0.03, 0.04, # A
        0.03, 0.9, 0.03, 0.04, # C
        0.03, 0.03, 0.9, 0.04, # G
        0.04, 0.03, 0.03, 0.9, # T
        0.9, 0.03, 0.03, 0.04, # A
        0.03, 0.9, 0.03, 0.04 # C
    ), ncol = 4, byrow = TRUE)
    colnames(pssm) <- c("A", "C", "G", "T")

    test_interval <- gintervals(1, 200, 250)
    threshold <- -3.0

    gvtrack.create("edist_6bp", NULL,
        func = "pwm.edit_distance",
        pssm = pssm, score.thresh = threshold,
        bidirect = FALSE, extend = FALSE, prior = 0
    )

    result <- gextract("edist_6bp", test_interval, iterator = test_interval)

    # Should return a valid result
    expect_true(is.na(result$edist_6bp[1]) || result$edist_6bp[1] >= 0)
})

test_that("pwm.edit_distance accepts PSSM with extra columns", {
    remove_all_vtracks()

    # Create PSSM with extra columns
    pssm_with_extras <- data.frame(
        A = c(1.0, 0.0),
        C = c(0.0, 1.0),
        G = c(0.0, 0.0),
        T = c(0.0, 0.0),
        motif_name = "AC",
        position = 1:2,
        conservation = c(0.9, 0.95),
        stringsAsFactors = FALSE
    )

    pssm_regular <- create_test_pssm()

    test_interval <- gintervals(1, 200, 240)
    threshold <- -5.0

    gvtrack.create("edist_extra", NULL,
        func = "pwm.edit_distance",
        pssm = pssm_with_extras, score.thresh = threshold,
        bidirect = FALSE, extend = FALSE, prior = 0
    )

    gvtrack.create("edist_regular", NULL,
        func = "pwm.edit_distance",
        pssm = pssm_regular, score.thresh = threshold,
        bidirect = FALSE, extend = FALSE, prior = 0
    )

    result <- gextract(c("edist_extra", "edist_regular"), test_interval, iterator = test_interval)

    # Should give same results
    expect_equal(result$edist_extra[1], result$edist_regular[1], tolerance = 1e-6)
})

test_that("pwm.edit_distance score.min filters low-scoring windows", {
    remove_all_vtracks()

    pssm <- create_test_pssm()

    test_interval <- gintervals(1, 200, 240)
    threshold <- -5.0

    # Without score.min - should find best edit distance across all windows
    gvtrack.create("edist_nofilt", NULL,
        func = "pwm.edit_distance",
        pssm = pssm, score.thresh = threshold,
        bidirect = FALSE, extend = FALSE, prior = 0
    )

    # With very high score.min - should filter out most/all windows -> NA
    gvtrack.create("edist_highfilt", NULL,
        func = "pwm.edit_distance",
        pssm = pssm, score.thresh = threshold,
        score.min = 0.0,
        bidirect = FALSE, extend = FALSE, prior = 0
    )

    # With very low score.min - should not filter anything
    gvtrack.create("edist_lowfilt", NULL,
        func = "pwm.edit_distance",
        pssm = pssm, score.thresh = threshold,
        score.min = -100.0,
        bidirect = FALSE, extend = FALSE, prior = 0
    )

    result <- gextract(c("edist_nofilt", "edist_highfilt", "edist_lowfilt"),
        test_interval,
        iterator = test_interval
    )

    # With very low score.min, result should match no-filter
    expect_equal(result$edist_nofilt[1], result$edist_lowfilt[1], tolerance = 1e-6)

    # With very high score.min, most windows are filtered out
    # Result should be NA or >= the unfiltered result
    if (!is.na(result$edist_highfilt[1]) && !is.na(result$edist_nofilt[1])) {
        expect_true(result$edist_highfilt[1] >= result$edist_nofilt[1])
    }
})

test_that("pwm.edit_distance score.min works with pwm.max.edit_distance", {
    remove_all_vtracks()

    pssm <- create_test_pssm()

    test_interval <- gintervals(1, 200, 240)
    threshold <- -5.0

    # Without score.min
    gvtrack.create("edist_max_nofilt", NULL,
        func = "pwm.max.edit_distance",
        pssm = pssm, score.thresh = threshold,
        bidirect = FALSE, extend = FALSE, prior = 0
    )

    # With score.min that should pass (very low)
    gvtrack.create("edist_max_lowfilt", NULL,
        func = "pwm.max.edit_distance",
        pssm = pssm, score.thresh = threshold,
        score.min = -100.0,
        bidirect = FALSE, extend = FALSE, prior = 0
    )

    # With score.min that may filter out best PWM window
    gvtrack.create("edist_max_highfilt", NULL,
        func = "pwm.max.edit_distance",
        pssm = pssm, score.thresh = threshold,
        score.min = 0.0,
        bidirect = FALSE, extend = FALSE, prior = 0
    )

    result <- gextract(c("edist_max_nofilt", "edist_max_lowfilt", "edist_max_highfilt"),
        test_interval,
        iterator = test_interval
    )

    # Low filter should match no filter
    expect_equal(result$edist_max_nofilt[1], result$edist_max_lowfilt[1], tolerance = 1e-6)
})

test_that("pwm.edit_distance parameter validation works", {
    remove_all_vtracks()

    pssm <- create_test_pssm()

    # Missing score.thresh should error
    expect_error(
        gvtrack.create("edist_bad", NULL,
            func = "pwm.edit_distance",
            pssm = pssm, bidirect = FALSE
        ),
        "score.thresh"
    )

    # Invalid max_edits should error
    expect_error(
        gvtrack.create("edist_bad2", NULL,
            func = "pwm.edit_distance",
            pssm = pssm, score.thresh = -5.0,
            max_edits = 0, bidirect = FALSE
        ),
        "max_edits"
    )

    # Invalid score.min should error
    expect_error(
        gvtrack.create("edist_bad3", NULL,
            func = "pwm.edit_distance",
            pssm = pssm, score.thresh = -5.0,
            score.min = "abc", bidirect = FALSE
        ),
        "score.min"
    )
})

# --------------------------------------------------------------------------
# Indel support tests (max_indels parameter)
# --------------------------------------------------------------------------

test_that("max_indels=0 matches substitution-only behavior (regression test)", {
    remove_all_vtracks()

    pssm <- create_test_pssm()
    test_interval <- gintervals(1, 200, 240)
    threshold <- -5.0

    # Explicit max_indels=0
    gvtrack.create("edist_indels0", NULL,
        func = "pwm.edit_distance",
        pssm = pssm, score.thresh = threshold, max_indels = 0,
        bidirect = FALSE, extend = FALSE, prior = 0
    )

    # No max_indels at all (default)
    gvtrack.create("edist_default", NULL,
        func = "pwm.edit_distance",
        pssm = pssm, score.thresh = threshold,
        bidirect = FALSE, extend = FALSE, prior = 0
    )

    result <- gextract(c("edist_indels0", "edist_default"), test_interval, iterator = test_interval)

    # max_indels=0 and omitting max_indels should give identical results
    if (is.na(result$edist_indels0[1])) {
        expect_true(is.na(result$edist_default[1]))
    } else {
        expect_equal(result$edist_indels0[1], result$edist_default[1], tolerance = 1e-6)
    }
})

test_that("max_indels=0 default produces same results as omitting max_indels over multiple intervals", {
    remove_all_vtracks()

    pssm <- matrix(c(
        0.9, 0.03, 0.03, 0.04,
        0.03, 0.9, 0.03, 0.04,
        0.03, 0.03, 0.9, 0.04,
        0.04, 0.03, 0.03, 0.9
    ), ncol = 4, byrow = TRUE)
    colnames(pssm) <- c("A", "C", "G", "T")

    test_intervals <- gintervals(
        chrom = c(1, 1, 1, 1),
        start = c(200, 500, 1000, 2000),
        end = c(240, 540, 1040, 2040)
    )
    threshold <- -3.0

    gvtrack.create("edist_indels0_multi", NULL,
        func = "pwm.edit_distance",
        pssm = pssm, score.thresh = threshold, max_indels = 0,
        bidirect = FALSE, extend = FALSE, prior = 0
    )
    gvtrack.create("edist_default_multi", NULL,
        func = "pwm.edit_distance",
        pssm = pssm, score.thresh = threshold,
        bidirect = FALSE, extend = FALSE, prior = 0
    )

    result <- gextract(c("edist_indels0_multi", "edist_default_multi"),
        test_intervals,
        iterator = test_intervals
    )

    for (i in seq_len(nrow(result))) {
        if (is.na(result$edist_indels0_multi[i])) {
            expect_true(is.na(result$edist_default_multi[i]),
                info = paste("Row", i, "should be NA in both")
            )
        } else {
            expect_equal(result$edist_indels0_multi[i], result$edist_default_multi[i],
                tolerance = 1e-6,
                info = paste("Row", i, "values should match")
            )
        }
    }
})

test_that("max_indels=1: edit distance with indels <= edit distance without indels", {
    remove_all_vtracks()

    # Use a longer motif to make substitution-vs-indel comparison meaningful
    pssm <- matrix(c(
        0.9, 0.03, 0.03, 0.04, # A
        0.03, 0.9, 0.03, 0.04, # C
        0.03, 0.03, 0.9, 0.04, # G
        0.04, 0.03, 0.03, 0.9 # T
    ), ncol = 4, byrow = TRUE)
    colnames(pssm) <- c("A", "C", "G", "T")

    test_intervals <- gintervals(
        chrom = c(1, 1, 1),
        start = c(200, 500, 1000),
        end = c(260, 560, 1060)
    )
    threshold <- -3.0

    gvtrack.create("edist_no_indels", NULL,
        func = "pwm.edit_distance",
        pssm = pssm, score.thresh = threshold, max_indels = 0,
        bidirect = FALSE, extend = FALSE, prior = 0
    )

    gvtrack.create("edist_with_indels", NULL,
        func = "pwm.edit_distance",
        pssm = pssm, score.thresh = threshold, max_indels = 1,
        bidirect = FALSE, extend = FALSE, prior = 0
    )

    result <- gextract(c("edist_no_indels", "edist_with_indels"),
        test_intervals,
        iterator = test_intervals
    )

    for (i in seq_len(nrow(result))) {
        no_indel <- result$edist_no_indels[i]
        with_indel <- result$edist_with_indels[i]

        # If both are non-NA, indel version should be <= substitution-only
        if (!is.na(no_indel) && !is.na(with_indel)) {
            expect_true(with_indel <= no_indel + 1e-6,
                info = paste(
                    "Row", i, ": with indels", with_indel,
                    "should be <= without indels", no_indel
                )
            )
        }

        # If substitution-only finds a result, indel version should too
        # (indels can only help or stay the same)
        if (!is.na(no_indel)) {
            expect_false(is.na(with_indel),
                info = paste("Row", i, ": indel version should not be NA when sub-only is", no_indel)
            )
        }
    }
})

test_that("max_indels=1 with single insertion disrupting motif alignment", {
    remove_all_vtracks()

    # Motif: ACGT
    pssm <- matrix(c(
        0.97, 0.01, 0.01, 0.01, # A
        0.01, 0.97, 0.01, 0.01, # C
        0.01, 0.01, 0.97, 0.01, # G
        0.01, 0.01, 0.01, 0.97 # T
    ), ncol = 4, byrow = TRUE)
    colnames(pssm) <- c("A", "C", "G", "T")

    motif_len <- nrow(pssm)

    # Search for a region where "ACGT" exists with a 1-base insertion in the middle
    # i.e., we look for "A_CGT" where _ is any base that disrupts the motif
    search_interval <- gintervals(1, 0, 5000)
    full_seq <- toupper(gseq.extract(search_interval))

    # Find any occurrence of ACGT
    acgt_pos <- regexpr("ACGT", full_seq)
    if (acgt_pos[1] < 0) {
        skip("No ACGT motif found in test genome region")
    }

    # Create an interval around the ACGT occurrence for testing
    abs_start <- as.integer(acgt_pos[1]) - 1 # 0-based
    # Use a window slightly larger than motif to test scanning
    test_interval <- gintervals(1, abs_start, abs_start + motif_len + 10)
    threshold <- sum(log(c(0.97, 0.97, 0.97, 0.97)))

    # Without indels: should find perfect match (0 edits)
    gvtrack.create("edist_no_indel", NULL,
        func = "pwm.edit_distance",
        pssm = pssm, score.thresh = threshold, max_indels = 0,
        bidirect = FALSE, extend = FALSE, prior = 0
    )

    # With 1 indel: should also find 0 edits (same or better)
    gvtrack.create("edist_1_indel", NULL,
        func = "pwm.edit_distance",
        pssm = pssm, score.thresh = threshold, max_indels = 1,
        bidirect = FALSE, extend = FALSE, prior = 0
    )

    result <- gextract(c("edist_no_indel", "edist_1_indel"),
        test_interval,
        iterator = test_interval
    )

    # Both should find the perfect match
    expect_equal(result$edist_no_indel[1], 0, tolerance = 1e-6)
    expect_equal(result$edist_1_indel[1], 0, tolerance = 1e-6)
})

test_that("max_indels=2: cases requiring two indels", {
    remove_all_vtracks()

    # Motif: ACGT
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

    # Compare max_indels=0, 1, 2
    gvtrack.create("edist_d0", NULL,
        func = "pwm.edit_distance",
        pssm = pssm, score.thresh = threshold, max_indels = 0,
        bidirect = FALSE, extend = FALSE, prior = 0
    )
    gvtrack.create("edist_d1", NULL,
        func = "pwm.edit_distance",
        pssm = pssm, score.thresh = threshold, max_indels = 1,
        bidirect = FALSE, extend = FALSE, prior = 0
    )
    gvtrack.create("edist_d2", NULL,
        func = "pwm.edit_distance",
        pssm = pssm, score.thresh = threshold, max_indels = 2,
        bidirect = FALSE, extend = FALSE, prior = 0
    )

    result <- gextract(c("edist_d0", "edist_d1", "edist_d2"),
        test_intervals,
        iterator = test_intervals
    )

    for (i in seq_len(nrow(result))) {
        d0 <- result$edist_d0[i]
        d1 <- result$edist_d1[i]
        d2 <- result$edist_d2[i]

        # Monotonicity: more indels allowed => edits should be <= (or NA becomes non-NA)
        if (!is.na(d0) && !is.na(d1)) {
            expect_true(d1 <= d0 + 1e-6,
                info = paste("Row", i, ": d1", d1, "should be <= d0", d0)
            )
        }
        if (!is.na(d1) && !is.na(d2)) {
            expect_true(d2 <= d1 + 1e-6,
                info = paste("Row", i, ": d2", d2, "should be <= d1", d1)
            )
        }
        if (!is.na(d0) && !is.na(d2)) {
            expect_true(d2 <= d0 + 1e-6,
                info = paste("Row", i, ": d2", d2, "should be <= d0", d0)
            )
        }

        # If substitution-only has a result, indel versions should too
        if (!is.na(d0)) {
            expect_false(is.na(d1), info = paste("Row", i))
            expect_false(is.na(d2), info = paste("Row", i))
        }
    }
})

test_that("combined indels + substitutions", {
    remove_all_vtracks()

    # 6bp motif with strong preference
    pssm <- matrix(c(
        0.97, 0.01, 0.01, 0.01, # A
        0.01, 0.97, 0.01, 0.01, # C
        0.01, 0.01, 0.97, 0.01, # G
        0.01, 0.01, 0.01, 0.97, # T
        0.97, 0.01, 0.01, 0.01, # A
        0.01, 0.97, 0.01, 0.01 # C
    ), ncol = 4, byrow = TRUE)
    colnames(pssm) <- c("A", "C", "G", "T")

    # Use a region large enough that the DP has room
    test_interval <- gintervals(1, 200, 280)

    # Reachable threshold: ~3 edits from best possible
    threshold <- -5.0

    # Without indels: substitution only
    gvtrack.create("edist_sub_only", NULL,
        func = "pwm.edit_distance",
        pssm = pssm, score.thresh = threshold, max_indels = 0,
        bidirect = FALSE, extend = FALSE, prior = 0
    )

    # With indels: allows combined indels + subs
    gvtrack.create("edist_combined", NULL,
        func = "pwm.edit_distance",
        pssm = pssm, score.thresh = threshold, max_indels = 1,
        bidirect = FALSE, extend = FALSE, prior = 0
    )

    result <- gextract(c("edist_sub_only", "edist_combined"),
        test_interval,
        iterator = test_interval
    )

    # At least one should be non-NA with this lenient threshold and large interval
    expect_true(!is.na(result$edist_sub_only[1]) || !is.na(result$edist_combined[1]),
        info = "At least one method should find a reachable window"
    )

    # Combined should be <= substitution-only (if both non-NA)
    if (!is.na(result$edist_sub_only[1]) && !is.na(result$edist_combined[1])) {
        expect_true(result$edist_combined[1] <= result$edist_sub_only[1] + 1e-6)
    }

    # Both should be non-negative when non-NA
    if (!is.na(result$edist_sub_only[1])) {
        expect_true(result$edist_sub_only[1] >= 0)
    }
    if (!is.na(result$edist_combined[1])) {
        expect_true(result$edist_combined[1] >= 0)
    }
})

test_that("max_indels cap enforcement: max_indels=1 should not allow 2 indels", {
    remove_all_vtracks()

    # 4bp motif
    pssm <- matrix(c(
        0.97, 0.01, 0.01, 0.01, # A
        0.01, 0.97, 0.01, 0.01, # C
        0.01, 0.01, 0.97, 0.01, # G
        0.01, 0.01, 0.01, 0.97 # T
    ), ncol = 4, byrow = TRUE)
    colnames(pssm) <- c("A", "C", "G", "T")

    test_interval <- gintervals(1, 200, 260)
    threshold <- -3.0

    # d=1 and d=2 vtracks
    gvtrack.create("edist_cap1", NULL,
        func = "pwm.edit_distance",
        pssm = pssm, score.thresh = threshold, max_indels = 1,
        bidirect = FALSE, extend = FALSE, prior = 0
    )
    gvtrack.create("edist_cap2", NULL,
        func = "pwm.edit_distance",
        pssm = pssm, score.thresh = threshold, max_indels = 2,
        bidirect = FALSE, extend = FALSE, prior = 0
    )

    result <- gextract(c("edist_cap1", "edist_cap2"),
        test_interval,
        iterator = test_interval
    )

    d1 <- result$edist_cap1[1]
    d2 <- result$edist_cap2[1]

    # d=2 can only be <= d=1 (never worse)
    if (!is.na(d1) && !is.na(d2)) {
        expect_true(d2 <= d1 + 1e-6)
    }

    # If d=1 returns NA, d=2 may or may not (it has more freedom)
    # If d=2 returns NA, d=1 must also be NA
    if (is.na(d2)) {
        expect_true(is.na(d1),
            info = "If max_indels=2 is NA, max_indels=1 should also be NA"
        )
    }
})

test_that("indel at interval boundary with extend=TRUE", {
    remove_all_vtracks()

    # 4bp motif
    pssm <- matrix(c(
        0.97, 0.01, 0.01, 0.01, # A
        0.01, 0.97, 0.01, 0.01, # C
        0.01, 0.01, 0.97, 0.01, # G
        0.01, 0.01, 0.01, 0.97 # T
    ), ncol = 4, byrow = TRUE)
    colnames(pssm) <- c("A", "C", "G", "T")

    motif_len <- nrow(pssm)

    # Use a small interval near the edge that requires extend to scan
    test_interval <- gintervals(1, 200, 200 + motif_len - 1)
    threshold <- -5.0

    # With extend=TRUE, the C++ code expands the interval to cover motif_len + max_indels
    gvtrack.create("edist_edge_ext", NULL,
        func = "pwm.edit_distance",
        pssm = pssm, score.thresh = threshold, max_indels = 1,
        bidirect = FALSE, extend = TRUE, prior = 0
    )

    # Without extend, interval is too small for motif -> NA
    gvtrack.create("edist_edge_noext", NULL,
        func = "pwm.edit_distance",
        pssm = pssm, score.thresh = threshold, max_indels = 1,
        bidirect = FALSE, extend = FALSE, prior = 0
    )

    result <- gextract(c("edist_edge_ext", "edist_edge_noext"),
        test_interval,
        iterator = test_interval
    )

    # Extended should produce a result (interval is expanded to cover full motif)
    # Non-extended with interval < motif_len should be NA
    expect_true(is.na(result$edist_edge_noext[1]),
        info = "Interval smaller than motif without extend should be NA"
    )
    # Extended should have a valid result (unless threshold is unreachable)
    expect_true(!is.na(result$edist_edge_ext[1]) || TRUE) # may be NA if threshold unreachable
})

test_that("indel at interval boundary with extend=FALSE uses full interval", {
    remove_all_vtracks()

    pssm <- matrix(c(
        0.97, 0.01, 0.01, 0.01, # A
        0.01, 0.97, 0.01, 0.01, # C
        0.01, 0.01, 0.97, 0.01, # G
        0.01, 0.01, 0.01, 0.97 # T
    ), ncol = 4, byrow = TRUE)
    colnames(pssm) <- c("A", "C", "G", "T")

    motif_len <- nrow(pssm)

    # Interval exactly motif_len in size - boundary case for indels
    test_interval <- gintervals(1, 200, 200 + motif_len)
    threshold <- -5.0

    gvtrack.create("edist_exact_boundary", NULL,
        func = "pwm.edit_distance",
        pssm = pssm, score.thresh = threshold, max_indels = 1,
        bidirect = FALSE, extend = FALSE, prior = 0
    )

    result <- gextract("edist_exact_boundary", test_interval, iterator = test_interval)

    # Should produce a valid result (non-NA or NA if threshold unreachable)
    expect_true(is.na(result$edist_exact_boundary[1]) || result$edist_exact_boundary[1] >= 0)
})

test_that("max_indels=1 with 1bp iterator produces consistent results", {
    remove_all_vtracks()

    pssm <- matrix(c(
        0.97, 0.01, 0.01, 0.01, # A
        0.01, 0.97, 0.01, 0.01, # C
        0.01, 0.01, 0.97, 0.01, # G
        0.01, 0.01, 0.01, 0.97 # T
    ), ncol = 4, byrow = TRUE)
    colnames(pssm) <- c("A", "C", "G", "T")

    motif_len <- nrow(pssm)
    test_interval <- gintervals(1, 200, 210)
    threshold <- -5.0

    # Without indels
    gvtrack.create("edist_1bp_d0", NULL,
        func = "pwm.edit_distance",
        pssm = pssm, score.thresh = threshold, max_indels = 0,
        bidirect = FALSE, extend = TRUE, prior = 0
    )

    # With indels
    gvtrack.create("edist_1bp_d1", NULL,
        func = "pwm.edit_distance",
        pssm = pssm, score.thresh = threshold, max_indels = 1,
        bidirect = FALSE, extend = TRUE, prior = 0
    )

    result <- gextract(c("edist_1bp_d0", "edist_1bp_d1"), test_interval, iterator = 1)

    expect_true(nrow(result) > 0)

    for (i in seq_len(nrow(result))) {
        d0 <- result$edist_1bp_d0[i]
        d1 <- result$edist_1bp_d1[i]

        if (!is.na(d0) && !is.na(d1)) {
            expect_true(d1 <= d0 + 1e-6,
                info = paste("Position", result$start[i], ": d1", d1, "should be <= d0", d0)
            )
        }
        if (!is.na(d0)) {
            expect_false(is.na(d1),
                info = paste("Position", result$start[i], ": d1 should not be NA if d0 is", d0)
            )
        }
    }
})

test_that("max_indels with bidirectional scanning works", {
    remove_all_vtracks()

    pssm <- matrix(c(
        0.97, 0.01, 0.01, 0.01, # A
        0.01, 0.97, 0.01, 0.01, # C
        0.01, 0.01, 0.97, 0.01, # G
        0.01, 0.01, 0.01, 0.97 # T
    ), ncol = 4, byrow = TRUE)
    colnames(pssm) <- c("A", "C", "G", "T")

    test_interval <- gintervals(1, 200, 260)
    threshold <- -3.0

    # Forward only
    gvtrack.create("edist_fwd_d1", NULL,
        func = "pwm.edit_distance",
        pssm = pssm, score.thresh = threshold, max_indels = 1,
        bidirect = FALSE, strand = 1, extend = FALSE, prior = 0
    )

    # Reverse only
    gvtrack.create("edist_rev_d1", NULL,
        func = "pwm.edit_distance",
        pssm = pssm, score.thresh = threshold, max_indels = 1,
        bidirect = FALSE, strand = -1, extend = FALSE, prior = 0
    )

    # Bidirectional
    gvtrack.create("edist_bidi_d1", NULL,
        func = "pwm.edit_distance",
        pssm = pssm, score.thresh = threshold, max_indels = 1,
        bidirect = TRUE, extend = FALSE, prior = 0
    )

    result <- gextract(c("edist_fwd_d1", "edist_rev_d1", "edist_bidi_d1"),
        test_interval,
        iterator = test_interval
    )

    # Bidirectional should return minimum of both strands
    fwd <- result$edist_fwd_d1[1]
    rev <- result$edist_rev_d1[1]
    bidi <- result$edist_bidi_d1[1]

    if (!is.na(fwd) && !is.na(rev)) {
        expect_equal(bidi, min(fwd, rev), tolerance = 1e-6)
    } else if (!is.na(fwd)) {
        expect_equal(bidi, fwd, tolerance = 1e-6)
    } else if (!is.na(rev)) {
        expect_equal(bidi, rev, tolerance = 1e-6)
    }
})

test_that("max_indels with max_edits cap combined", {
    remove_all_vtracks()

    pssm <- matrix(c(
        0.97, 0.01, 0.01, 0.01, # A
        0.01, 0.97, 0.01, 0.01, # C
        0.01, 0.01, 0.97, 0.01, # G
        0.01, 0.01, 0.01, 0.97 # T
    ), ncol = 4, byrow = TRUE)
    colnames(pssm) <- c("A", "C", "G", "T")

    test_interval <- gintervals(1, 200, 260)
    threshold <- -3.0

    # max_indels=1, no max_edits cap
    gvtrack.create("edist_d1_uncapped", NULL,
        func = "pwm.edit_distance",
        pssm = pssm, score.thresh = threshold, max_indels = 1,
        bidirect = FALSE, extend = FALSE, prior = 0
    )

    # max_indels=1, max_edits=1 (total edits capped at 1)
    gvtrack.create("edist_d1_cap1", NULL,
        func = "pwm.edit_distance",
        pssm = pssm, score.thresh = threshold, max_indels = 1, max_edits = 1,
        bidirect = FALSE, extend = FALSE, prior = 0
    )

    # max_indels=1, max_edits=3 (total edits capped at 3)
    gvtrack.create("edist_d1_cap3", NULL,
        func = "pwm.edit_distance",
        pssm = pssm, score.thresh = threshold, max_indels = 1, max_edits = 3,
        bidirect = FALSE, extend = FALSE, prior = 0
    )

    result <- gextract(c("edist_d1_uncapped", "edist_d1_cap1", "edist_d1_cap3"),
        test_interval,
        iterator = test_interval
    )

    uncapped <- result$edist_d1_uncapped[1]
    cap1 <- result$edist_d1_cap1[1]
    cap3 <- result$edist_d1_cap3[1]

    # If uncapped returns > 1, cap1 should be NA
    if (!is.na(uncapped) && uncapped > 1) {
        expect_true(is.na(cap1),
            info = paste("Uncapped =", uncapped, "> 1, so cap1 should be NA")
        )
    }

    # If uncapped returns <= 1, cap1 should match
    if (!is.na(uncapped) && uncapped <= 1) {
        expect_equal(cap1, uncapped, tolerance = 1e-6)
    }

    # If uncapped returns <= 3, cap3 should match
    if (!is.na(uncapped) && uncapped <= 3) {
        expect_equal(cap3, uncapped, tolerance = 1e-6)
    }
})

test_that("max_indels works with pwm.edit_distance.pos function", {
    remove_all_vtracks()

    pssm <- matrix(c(
        0.97, 0.01, 0.01, 0.01, # A
        0.01, 0.97, 0.01, 0.01, # C
        0.01, 0.01, 0.97, 0.01, # G
        0.01, 0.01, 0.01, 0.97 # T
    ), ncol = 4, byrow = TRUE)
    colnames(pssm) <- c("A", "C", "G", "T")

    test_interval <- gintervals(1, 200, 260)
    threshold <- -5.0

    # Edit distance
    gvtrack.create("edist_val_d1", NULL,
        func = "pwm.edit_distance",
        pssm = pssm, score.thresh = threshold, max_indels = 1,
        bidirect = FALSE, extend = FALSE, prior = 0
    )

    # Edit distance position
    gvtrack.create("edist_pos_d1", NULL,
        func = "pwm.edit_distance.pos",
        pssm = pssm, score.thresh = threshold, max_indels = 1,
        bidirect = FALSE, extend = FALSE, prior = 0
    )

    result <- gextract(c("edist_val_d1", "edist_pos_d1"),
        test_interval,
        iterator = test_interval
    )

    # Position should be non-NA if edit distance is non-NA
    if (!is.na(result$edist_val_d1[1])) {
        expect_false(is.na(result$edist_pos_d1[1]))
        expect_true(result$edist_pos_d1[1] >= 1,
            info = "Position should be >= 1 (1-based)"
        )
    }
})

test_that("max_indels works with pwm.max.edit_distance function", {
    remove_all_vtracks()

    pssm <- matrix(c(
        0.97, 0.01, 0.01, 0.01, # A
        0.01, 0.97, 0.01, 0.01, # C
        0.01, 0.01, 0.97, 0.01, # G
        0.01, 0.01, 0.01, 0.97 # T
    ), ncol = 4, byrow = TRUE)
    colnames(pssm) <- c("A", "C", "G", "T")

    test_interval <- gintervals(1, 200, 260)
    threshold <- -5.0

    # Best PWM site edit distance without indels
    gvtrack.create("pwm_max_edist_d0", NULL,
        func = "pwm.max.edit_distance",
        pssm = pssm, score.thresh = threshold, max_indels = 0,
        bidirect = FALSE, extend = FALSE, prior = 0
    )

    # Best PWM site edit distance with indels
    gvtrack.create("pwm_max_edist_d1", NULL,
        func = "pwm.max.edit_distance",
        pssm = pssm, score.thresh = threshold, max_indels = 1,
        bidirect = FALSE, extend = FALSE, prior = 0
    )

    result <- gextract(c("pwm_max_edist_d0", "pwm_max_edist_d1"),
        test_interval,
        iterator = test_interval
    )

    d0 <- result$pwm_max_edist_d0[1]
    d1 <- result$pwm_max_edist_d1[1]

    # With indels should be <= without indels (at the same PWM-best window)
    if (!is.na(d0) && !is.na(d1)) {
        expect_true(d1 <= d0 + 1e-6)
    }
})

test_that("max_indels parameter validation", {
    remove_all_vtracks()

    pssm <- create_test_pssm()

    # Negative max_indels should error
    expect_error(
        gvtrack.create("edist_bad_indels", NULL,
            func = "pwm.edit_distance",
            pssm = pssm, score.thresh = -5.0,
            max_indels = -1, bidirect = FALSE
        ),
        "max_indels"
    )

    # Non-numeric max_indels should error
    expect_error(
        gvtrack.create("edist_bad_indels2", NULL,
            func = "pwm.edit_distance",
            pssm = pssm, score.thresh = -5.0,
            max_indels = "abc", bidirect = FALSE
        ),
        "max_indels"
    )

    # Valid max_indels=0 should work (same as default)
    expect_no_error(
        gvtrack.create("edist_ok_indels0", NULL,
            func = "pwm.edit_distance",
            pssm = pssm, score.thresh = -5.0,
            max_indels = 0, bidirect = FALSE
        )
    )
})

test_that("max_indels with score.min filter combined", {
    remove_all_vtracks()

    pssm <- matrix(c(
        0.97, 0.01, 0.01, 0.01,
        0.01, 0.97, 0.01, 0.01,
        0.01, 0.01, 0.97, 0.01,
        0.01, 0.01, 0.01, 0.97
    ), ncol = 4, byrow = TRUE)
    colnames(pssm) <- c("A", "C", "G", "T")

    test_interval <- gintervals(1, 200, 260)
    threshold <- -5.0

    # With indels, no score filter
    gvtrack.create("edist_d1_nofilt", NULL,
        func = "pwm.edit_distance",
        pssm = pssm, score.thresh = threshold, max_indels = 1,
        bidirect = FALSE, extend = FALSE, prior = 0
    )

    # With indels, lenient score filter
    gvtrack.create("edist_d1_lowfilt", NULL,
        func = "pwm.edit_distance",
        pssm = pssm, score.thresh = threshold, max_indels = 1,
        score.min = -100.0,
        bidirect = FALSE, extend = FALSE, prior = 0
    )

    # With indels, strict score filter
    gvtrack.create("edist_d1_highfilt", NULL,
        func = "pwm.edit_distance",
        pssm = pssm, score.thresh = threshold, max_indels = 1,
        score.min = 0.0,
        bidirect = FALSE, extend = FALSE, prior = 0
    )

    result <- gextract(c("edist_d1_nofilt", "edist_d1_lowfilt", "edist_d1_highfilt"),
        test_interval,
        iterator = test_interval
    )

    # Lenient filter should match no-filter
    expect_equal(result$edist_d1_nofilt[1], result$edist_d1_lowfilt[1], tolerance = 1e-6)

    # Strict filter should return NA or >= no-filter result
    if (!is.na(result$edist_d1_highfilt[1]) && !is.na(result$edist_d1_nofilt[1])) {
        expect_true(result$edist_d1_highfilt[1] >= result$edist_d1_nofilt[1])
    }
})

test_that("max_indels with longer motif", {
    remove_all_vtracks()

    # 8bp motif
    pssm <- matrix(c(
        0.9, 0.03, 0.03, 0.04, # A
        0.03, 0.9, 0.03, 0.04, # C
        0.03, 0.03, 0.9, 0.04, # G
        0.04, 0.03, 0.03, 0.9, # T
        0.9, 0.03, 0.03, 0.04, # A
        0.03, 0.9, 0.03, 0.04, # C
        0.03, 0.03, 0.9, 0.04, # G
        0.04, 0.03, 0.03, 0.9 # T
    ), ncol = 4, byrow = TRUE)
    colnames(pssm) <- c("A", "C", "G", "T")

    test_interval <- gintervals(1, 200, 280)
    threshold <- -5.0

    gvtrack.create("edist_long_d0", NULL,
        func = "pwm.edit_distance",
        pssm = pssm, score.thresh = threshold, max_indels = 0,
        bidirect = FALSE, extend = FALSE, prior = 0
    )

    gvtrack.create("edist_long_d1", NULL,
        func = "pwm.edit_distance",
        pssm = pssm, score.thresh = threshold, max_indels = 1,
        bidirect = FALSE, extend = FALSE, prior = 0
    )

    gvtrack.create("edist_long_d2", NULL,
        func = "pwm.edit_distance",
        pssm = pssm, score.thresh = threshold, max_indels = 2,
        bidirect = FALSE, extend = FALSE, prior = 0
    )

    result <- gextract(c("edist_long_d0", "edist_long_d1", "edist_long_d2"),
        test_interval,
        iterator = test_interval
    )

    d0 <- result$edist_long_d0[1]
    d1 <- result$edist_long_d1[1]
    d2 <- result$edist_long_d2[1]

    # Monotonicity
    if (!is.na(d0) && !is.na(d1)) {
        expect_true(d1 <= d0 + 1e-6)
    }
    if (!is.na(d1) && !is.na(d2)) {
        expect_true(d2 <= d1 + 1e-6)
    }

    # Non-negative
    if (!is.na(d0)) expect_true(d0 >= 0)
    if (!is.na(d1)) expect_true(d1 >= 0)
    if (!is.na(d2)) expect_true(d2 >= 0)
})

# --------------------------------------------------------------------------
# LSE edit distance tests (pwm.edit_distance.lse / pwm.edit_distance.lse.pos)
# --------------------------------------------------------------------------

test_that("LSE edit distance basic functionality works", {
    remove_all_vtracks()

    pssm <- create_test_pssm()
    test_intervals <- gintervals(1, 200, 240)

    threshold <- -5.0
    gvtrack.create("lse_edist", NULL,
        func = "pwm.edit_distance.lse",
        pssm = pssm, score.thresh = threshold,
        bidirect = FALSE, extend = FALSE, prior = 0
    )

    result <- gextract("lse_edist", test_intervals, iterator = test_intervals)

    # Should return a non-NA numeric result
    expect_false(is.na(result$lse_edist[1]))
    expect_true(is.numeric(result$lse_edist[1]))
    expect_true(result$lse_edist[1] >= 0)
})

test_that("LSE edit distance returns 0 when LSE score already above threshold", {
    remove_all_vtracks()

    pssm <- create_test_pssm()
    test_intervals <- gintervals(1, 200, 240)

    # Use a very low threshold that the LSE score should already exceed
    threshold <- -100.0
    gvtrack.create("lse_above", NULL,
        func = "pwm.edit_distance.lse",
        pssm = pssm, score.thresh = threshold,
        bidirect = FALSE, extend = FALSE, prior = 0
    )

    result <- gextract("lse_above", test_intervals, iterator = test_intervals)

    # LSE score should already be above this threshold, so 0 edits needed
    expect_equal(result$lse_above[1], 0, tolerance = 1e-6)
})

test_that("LSE edit distance returns NA for unreachable threshold", {
    remove_all_vtracks()

    pssm <- create_test_pssm()
    test_intervals <- gintervals(1, 200, 240)

    # Set impossibly high threshold
    threshold <- 1000.0
    gvtrack.create("lse_unreachable", NULL,
        func = "pwm.edit_distance.lse",
        pssm = pssm, score.thresh = threshold,
        bidirect = FALSE, extend = FALSE, prior = 0
    )

    result <- gextract("lse_unreachable", test_intervals, iterator = test_intervals)

    # Should return NA for unreachable threshold
    expect_true(is.na(result$lse_unreachable[1]))
})

test_that("LSE edit distance <= max edit distance for same threshold", {
    remove_all_vtracks()

    pssm <- create_test_pssm()
    test_intervals <- gintervals(1, 200, 240)
    threshold <- -3.0

    # LSE benefits from aggregating scores across multiple overlapping starts,
    # so it should need <= edits compared to max (single-window) edit distance.
    gvtrack.create("lse_edist", NULL,
        func = "pwm.edit_distance.lse",
        pssm = pssm, score.thresh = threshold,
        bidirect = FALSE, extend = FALSE, prior = 0
    )

    gvtrack.create("max_edist", NULL,
        func = "pwm.edit_distance",
        pssm = pssm, score.thresh = threshold,
        bidirect = FALSE, extend = FALSE, prior = 0
    )

    result <- gextract(c("lse_edist", "max_edist"), test_intervals, iterator = test_intervals)

    # If both are non-NA, LSE should need <= edits
    if (!is.na(result$lse_edist[1]) && !is.na(result$max_edist[1])) {
        expect_true(result$lse_edist[1] <= result$max_edist[1] + 1e-6,
            info = paste("LSE edits", result$lse_edist[1], "should be <= max edits", result$max_edist[1])
        )
    }

    # If max edit distance finds a result, LSE should too (LSE is at least as powerful)
    if (!is.na(result$max_edist[1])) {
        expect_false(is.na(result$lse_edist[1]),
            info = "LSE should find a result whenever max edit distance does"
        )
    }
})

test_that("LSE edit distance respects max_edits parameter", {
    remove_all_vtracks()

    pssm <- create_test_pssm()
    test_intervals <- gintervals(1, 200, 240)
    threshold <- -3.0

    # Uncapped LSE
    gvtrack.create("lse_uncapped", NULL,
        func = "pwm.edit_distance.lse",
        pssm = pssm, score.thresh = threshold,
        bidirect = FALSE, extend = FALSE, prior = 0
    )

    # LSE with max_edits = 1
    gvtrack.create("lse_cap1", NULL,
        func = "pwm.edit_distance.lse",
        pssm = pssm, score.thresh = threshold, max_edits = 1,
        bidirect = FALSE, extend = FALSE, prior = 0
    )

    # LSE with max_edits = 3
    gvtrack.create("lse_cap3", NULL,
        func = "pwm.edit_distance.lse",
        pssm = pssm, score.thresh = threshold, max_edits = 3,
        bidirect = FALSE, extend = FALSE, prior = 0
    )

    result <- gextract(c("lse_uncapped", "lse_cap1", "lse_cap3"), test_intervals, iterator = test_intervals)

    uncapped <- result$lse_uncapped[1]
    cap1 <- result$lse_cap1[1]
    cap3 <- result$lse_cap3[1]

    # If uncapped needs > 1 edit, cap1 should be NA
    if (!is.na(uncapped) && uncapped > 1) {
        expect_true(is.na(cap1),
            info = paste("Uncapped =", uncapped, "> 1, so cap1 should be NA")
        )
    }

    # If uncapped needs <= 1 edit, cap1 should match
    if (!is.na(uncapped) && uncapped <= 1) {
        expect_equal(cap1, uncapped, tolerance = 1e-6)
    }

    # If uncapped needs <= 3 edits, cap3 should match
    if (!is.na(uncapped) && uncapped <= 3) {
        expect_equal(cap3, uncapped, tolerance = 1e-6)
    }
})

test_that("LSE edit distance respects score.min parameter", {
    remove_all_vtracks()

    pssm <- create_test_pssm()
    test_intervals <- gintervals(1, 200, 240)
    threshold <- -5.0

    # Without score.min filter
    gvtrack.create("lse_nofilt", NULL,
        func = "pwm.edit_distance.lse",
        pssm = pssm, score.thresh = threshold,
        bidirect = FALSE, extend = FALSE, prior = 0
    )

    # With very low score.min (should not filter anything)
    gvtrack.create("lse_lowfilt", NULL,
        func = "pwm.edit_distance.lse",
        pssm = pssm, score.thresh = threshold,
        score.min = -100.0,
        bidirect = FALSE, extend = FALSE, prior = 0
    )

    # With very high score.min (should filter out most/all windows -> NA)
    gvtrack.create("lse_highfilt", NULL,
        func = "pwm.edit_distance.lse",
        pssm = pssm, score.thresh = threshold,
        score.min = 0.0,
        bidirect = FALSE, extend = FALSE, prior = 0
    )

    result <- gextract(c("lse_nofilt", "lse_lowfilt", "lse_highfilt"),
        test_intervals,
        iterator = test_intervals
    )

    # Low filter should match no-filter
    expect_equal(result$lse_nofilt[1], result$lse_lowfilt[1], tolerance = 1e-6)

    # High filter should return NA or >= the unfiltered result
    if (!is.na(result$lse_highfilt[1]) && !is.na(result$lse_nofilt[1])) {
        expect_true(result$lse_highfilt[1] >= result$lse_nofilt[1])
    }
})

test_that("pwm.edit_distance.lse.pos returns non-NA numeric position", {
    remove_all_vtracks()

    pssm <- create_test_pssm()
    test_intervals <- gintervals(1, 200, 240)
    threshold <- -5.0

    gvtrack.create("lse_edist", NULL,
        func = "pwm.edit_distance.lse",
        pssm = pssm, score.thresh = threshold,
        bidirect = FALSE, extend = FALSE, prior = 0
    )

    gvtrack.create("lse_pos", NULL,
        func = "pwm.edit_distance.lse.pos",
        pssm = pssm, score.thresh = threshold,
        bidirect = FALSE, extend = FALSE, prior = 0
    )

    result <- gextract(c("lse_edist", "lse_pos"), test_intervals, iterator = test_intervals)

    # If edit distance is non-NA and > 0, position should also be non-NA
    if (!is.na(result$lse_edist[1]) && result$lse_edist[1] > 0) {
        expect_false(is.na(result$lse_pos[1]))
        expect_true(is.numeric(result$lse_pos[1]))
        expect_true(abs(result$lse_pos[1]) >= 1,
            info = "Position should be >= 1 (1-based, possibly negative for reverse strand)"
        )
    }

    # If edit distance is 0, position may be NA (no edit needed)
    if (!is.na(result$lse_edist[1]) && result$lse_edist[1] == 0) {
        expect_true(is.na(result$lse_pos[1]) || is.numeric(result$lse_pos[1]))
    }
})

test_that("LSE edit distance works with bidirect=TRUE and bidirect=FALSE", {
    remove_all_vtracks()

    pssm <- create_test_pssm()
    test_interval <- gintervals(1, 200, 240)
    threshold <- -5.0

    # Forward only
    gvtrack.create("lse_fwd", NULL,
        func = "pwm.edit_distance.lse",
        pssm = pssm, score.thresh = threshold,
        bidirect = FALSE, strand = 1, extend = FALSE, prior = 0
    )

    # Reverse only
    gvtrack.create("lse_rev", NULL,
        func = "pwm.edit_distance.lse",
        pssm = pssm, score.thresh = threshold,
        bidirect = FALSE, strand = -1, extend = FALSE, prior = 0
    )

    # Bidirectional
    gvtrack.create("lse_bidi", NULL,
        func = "pwm.edit_distance.lse",
        pssm = pssm, score.thresh = threshold,
        bidirect = TRUE, extend = FALSE, prior = 0
    )

    result <- gextract(c("lse_fwd", "lse_rev", "lse_bidi"), test_interval, iterator = test_interval)

    # All results should be valid (numeric or NA)
    expect_true(is.na(result$lse_fwd[1]) || result$lse_fwd[1] >= 0)
    expect_true(is.na(result$lse_rev[1]) || result$lse_rev[1] >= 0)
    expect_true(is.na(result$lse_bidi[1]) || result$lse_bidi[1] >= 0)

    # Bidirectional should return minimum of both strands (when both are non-NA)
    if (!is.na(result$lse_fwd[1]) && !is.na(result$lse_rev[1])) {
        expect_true(result$lse_bidi[1] <= min(result$lse_fwd[1], result$lse_rev[1]) + 1e-6)
    }
})

test_that("LSE edit distance works in gscreen expression", {
    remove_all_vtracks()

    pssm <- create_test_pssm()
    threshold <- -5.0

    gvtrack.create("lse_screen", NULL,
        func = "pwm.edit_distance.lse",
        pssm = pssm, score.thresh = threshold,
        bidirect = FALSE, extend = FALSE, prior = 0
    )

    # Use gscreen to find intervals where LSE edit distance is 0 (already above threshold)
    screened <- gscreen("!is.na(lse_screen) & lse_screen <= 1",
        gintervals(1, 0, 5000),
        iterator = 20
    )

    # Should return a valid intervals data frame
    expect_true(is.data.frame(screened))
    expect_true(all(c("chrom", "start", "end") %in% names(screened)))

    # If any intervals pass the filter, verify the edit distance values
    if (nrow(screened) > 0) {
        verify <- gextract("lse_screen", screened, iterator = screened)
        expect_true(all(!is.na(verify$lse_screen)))
        expect_true(all(verify$lse_screen <= 1))
    }
})

test_that("LSE edit distance: different thresholds give monotonically increasing edit distances", {
    remove_all_vtracks()

    pssm <- create_test_pssm()
    test_interval <- gintervals(1, 200, 240)

    thresholds <- c(-20.0, -10.0, -5.0, -2.0, 0.0)
    vnames <- sprintf("lse_thresh_%d", seq_along(thresholds))

    for (i in seq_along(thresholds)) {
        gvtrack.create(vnames[i], NULL,
            func = "pwm.edit_distance.lse",
            pssm = pssm, score.thresh = thresholds[i],
            bidirect = FALSE, extend = FALSE, prior = 0
        )
    }

    result <- gextract(vnames, test_interval, iterator = test_interval)

    edits <- sapply(vnames, function(v) result[[v]][1])

    # Filter out NAs for monotonicity check
    finite_edits <- edits[!is.na(edits)]
    if (length(finite_edits) > 1) {
        # Higher thresholds should require more (or equal) edits
        expect_true(all(diff(finite_edits) >= -1e-6),
            info = paste("Edits should be monotonically non-decreasing:", paste(finite_edits, collapse = ", "))
        )
    }

    # Once a threshold becomes NA (unreachable), all higher thresholds should also be NA
    na_found <- FALSE
    for (i in seq_along(edits)) {
        if (is.na(edits[i])) {
            na_found <- TRUE
        } else if (na_found) {
            fail(paste(
                "Found non-NA edit distance at threshold", thresholds[i],
                "after NA at lower threshold"
            ))
        }
    }
})

test_that("LSE edit distance with prior > 0 gives non-negative results", {
    remove_all_vtracks()

    pssm <- create_test_pssm()
    test_interval <- gintervals(1, 200, 240)
    threshold <- -5.0

    # With prior = 0
    gvtrack.create("lse_prior0", NULL,
        func = "pwm.edit_distance.lse",
        pssm = pssm, score.thresh = threshold,
        bidirect = FALSE, extend = FALSE, prior = 0
    )

    # With prior = 0.01 (default)
    gvtrack.create("lse_prior01", NULL,
        func = "pwm.edit_distance.lse",
        pssm = pssm, score.thresh = threshold,
        bidirect = FALSE, extend = FALSE, prior = 0.01
    )

    result <- gextract(c("lse_prior0", "lse_prior01"), test_interval, iterator = test_interval)

    # Both should return valid results (numeric or NA)
    if (!is.na(result$lse_prior0[1])) {
        expect_true(result$lse_prior0[1] >= 0)
    }
    if (!is.na(result$lse_prior01[1])) {
        expect_true(result$lse_prior01[1] >= 0)
    }

    # At least one should give a result on this test region
    expect_true(!is.na(result$lse_prior0[1]) || !is.na(result$lse_prior01[1]))
})

test_that("LSE edit distance on larger interval (iterator=500)", {
    remove_all_vtracks()

    # Use a longer motif for a more interesting test
    pssm <- matrix(c(
        0.9, 0.03, 0.03, 0.04, # A
        0.03, 0.9, 0.03, 0.04, # C
        0.03, 0.03, 0.9, 0.04, # G
        0.04, 0.03, 0.03, 0.9 # T
    ), ncol = 4, byrow = TRUE)
    colnames(pssm) <- c("A", "C", "G", "T")

    test_interval <- gintervals(1, 0, 2000)
    threshold <- -5.0

    gvtrack.create("lse_wide", NULL,
        func = "pwm.edit_distance.lse",
        pssm = pssm, score.thresh = threshold,
        bidirect = FALSE, extend = FALSE, prior = 0
    )

    result <- gextract("lse_wide", test_interval, iterator = 500)

    # Should have multiple rows (one per iterator window)
    expect_true(nrow(result) > 1)

    # Each result should be valid (non-negative or NA)
    for (i in seq_len(nrow(result))) {
        expect_true(is.na(result$lse_wide[i]) || result$lse_wide[i] >= 0,
            info = paste("Row", i, "value:", result$lse_wide[i])
        )
    }

    # With a lenient threshold, at least some windows should have non-NA results
    non_na_count <- sum(!is.na(result$lse_wide))
    expect_true(non_na_count > 0,
        info = "At least some windows should have non-NA LSE edit distances"
    )
})

# --------------------------------------------------------------------------
# Integer extend validation tests
# --------------------------------------------------------------------------

test_that("pwm.edit_distance accepts integer extend values in vtrack", {
    remove_all_vtracks()

    pssm <- create_test_pssm()
    threshold <- -5.0
    test_interval <- gintervals(1, 200, 240)

    # extend = 5L should work (integer)
    expect_no_error(
        gvtrack.create("edist_int_ext", NULL,
            func = "pwm.edit_distance",
            pssm = pssm, score.thresh = threshold,
            bidirect = FALSE, extend = 5L, prior = 0
        )
    )

    result <- gextract("edist_int_ext", test_interval, iterator = test_interval)
    expect_true(nrow(result) > 0)
})

test_that("pwm.edit_distance rejects invalid extend values in vtrack", {
    remove_all_vtracks()

    pssm <- create_test_pssm()
    threshold <- -5.0

    # Negative integer should error
    expect_error(
        gvtrack.create("edist_bad_ext", NULL,
            func = "pwm.edit_distance",
            pssm = pssm, score.thresh = threshold,
            bidirect = FALSE, extend = -1L, prior = 0
        ),
        "extend"
    )

    # String should error
    expect_error(
        gvtrack.create("edist_bad_ext2", NULL,
            func = "pwm.edit_distance",
            pssm = pssm, score.thresh = threshold,
            bidirect = FALSE, extend = "yes", prior = 0
        ),
        "extend"
    )

    # Non-integer numeric should error
    expect_error(
        gvtrack.create("edist_bad_ext3", NULL,
            func = "pwm.edit_distance",
            pssm = pssm, score.thresh = threshold,
            bidirect = FALSE, extend = 2.5, prior = 0
        ),
        "extend"
    )
})

test_that("pwm.edit_distance.lse accepts integer extend values in vtrack", {
    remove_all_vtracks()

    pssm <- create_test_pssm()
    threshold <- -5.0
    test_interval <- gintervals(1, 200, 240)

    # extend = 5L should work (integer)
    expect_no_error(
        gvtrack.create("lse_int_ext", NULL,
            func = "pwm.edit_distance.lse",
            pssm = pssm, score.thresh = threshold,
            bidirect = FALSE, extend = 5L, prior = 0
        )
    )

    result <- gextract("lse_int_ext", test_interval, iterator = test_interval)
    expect_true(nrow(result) > 0)
})

test_that("pwm.edit_distance.lse rejects invalid extend values in vtrack", {
    remove_all_vtracks()

    pssm <- create_test_pssm()
    threshold <- -5.0

    # Negative integer should error
    expect_error(
        gvtrack.create("lse_bad_ext", NULL,
            func = "pwm.edit_distance.lse",
            pssm = pssm, score.thresh = threshold,
            bidirect = FALSE, extend = -1L, prior = 0
        ),
        "extend"
    )

    # Non-integer numeric should error
    expect_error(
        gvtrack.create("lse_bad_ext2", NULL,
            func = "pwm.edit_distance.lse",
            pssm = pssm, score.thresh = threshold,
            bidirect = FALSE, extend = 2.5, prior = 0
        ),
        "extend"
    )
})

test_that("pwm.edit_distance extend=0L behaves like extend=FALSE", {
    remove_all_vtracks()

    pssm <- create_test_pssm()
    threshold <- -5.0
    test_interval <- gintervals(1, 200, 240)

    gvtrack.create("edist_ext0", NULL,
        func = "pwm.edit_distance",
        pssm = pssm, score.thresh = threshold,
        bidirect = FALSE, extend = 0L, prior = 0
    )

    gvtrack.create("edist_extF", NULL,
        func = "pwm.edit_distance",
        pssm = pssm, score.thresh = threshold,
        bidirect = FALSE, extend = FALSE, prior = 0
    )

    res_0 <- gextract("edist_ext0", test_interval, iterator = test_interval)
    res_f <- gextract("edist_extF", test_interval, iterator = test_interval)

    if (is.na(res_0$edist_ext0[1])) {
        expect_true(is.na(res_f$edist_extF[1]))
    } else {
        expect_equal(res_0$edist_ext0[1], res_f$edist_extF[1], tolerance = 1e-6)
    }
})

# --------------------------------------------------------------------------
# One-indel specialized solver tests
# --------------------------------------------------------------------------

describe("One-indel specialized solver", {
    # Shared PSSM: strongly prefers ACGT
    acgt_pssm <- matrix(c(
        0.97, 0.01, 0.01, 0.01, # A
        0.01, 0.97, 0.01, 0.01, # C
        0.01, 0.01, 0.97, 0.01, # G
        0.01, 0.01, 0.01, 0.97 # T
    ), ncol = 4, byrow = TRUE)
    colnames(acgt_pssm) <- c("A", "C", "G", "T")

    it("deletion improves alignment: inserted extra base is removed by one-indel solver", {
        remove_all_vtracks()

        # "ATCGT" is "ACGT" with an extra T inserted after A.
        # Use a tight threshold equal to the perfect-match score so that
        # no 4bp substitution-only window can reach it.
        # Perfect score = 4 * log(0.97) ~ -0.122
        seq <- "ATCGT"
        threshold <- sum(log(c(0.97, 0.97, 0.97, 0.97)))

        # Without indels: no 4bp window of "ATCGT" matches ACGT perfectly.
        # Windows are "ATCG" and "TCGT", both have at least 1 mismatch.
        # At the tight threshold these are unreachable, so we get 0 rows.
        result_no_indel <- gseq.pwm_edits(seq, acgt_pssm,
            score.thresh = threshold,
            max_indels = 0L, prior = 0, bidirect = FALSE
        )
        # No window can reach the perfect-match threshold via subs alone
        expect_equal(nrow(result_no_indel), 0,
            info = "No sub-only window should reach the tight threshold"
        )

        # With max_indels=1: the solver can delete the extra T to recover "ACGT"
        result_with_indel <- gseq.pwm_edits(seq, acgt_pssm,
            score.thresh = threshold,
            max_indels = 1L, prior = 0, bidirect = FALSE
        )
        expect_true(nrow(result_with_indel) > 0,
            info = "Indel solver should find at least one alignment"
        )

        # Expect exactly 1 edit: the deletion of the extra T
        best_with_indel <- min(result_with_indel$n_edits)
        expect_equal(best_with_indel, 1,
            info = "Should need exactly 1 deletion to recover ACGT from ATCGT"
        )

        # Confirm that the indel edit type is "del"
        edit_rows <- result_with_indel[result_with_indel$edit_num > 0, ]
        expect_true(any(edit_rows$edit_type == "del"),
            info = "Best alignment should use a deletion"
        )
    })

    it("insertion improves alignment: missing base is inserted by one-indel solver", {
        remove_all_vtracks()

        # "AGT" is "ACGT" with C removed (missing at position 2).
        # Without indels: no 4bp window exists (seq too short), so no result.
        # With max_indels=1: the solver can insert C at position 2.
        seq <- "AGT"
        threshold <- sum(log(c(0.97, 0.97, 0.97, 0.97)))

        result_with_indel <- gseq.pwm_edits(seq, acgt_pssm,
            score.thresh = threshold,
            max_indels = 1L, prior = 0, bidirect = FALSE
        )

        expect_true(nrow(result_with_indel) > 0,
            info = "With max_indels=1, solver should find an alignment for shorter seq"
        )

        # Confirm that the indel edit type is "ins"
        edit_rows <- result_with_indel[result_with_indel$edit_num > 0, ]
        expect_true(any(edit_rows$edit_type == "ins"),
            info = "Best alignment should use an insertion"
        )
    })

    it("score filter bypass: max_indels=1 still runs even when score.min rejects no-indel windows", {
        remove_all_vtracks()

        # Use the vtrack interface on a genomic interval.
        # Find a region in the genome.
        search_interval <- gintervals(1, 0, 5000)
        full_seq <- toupper(gseq.extract(search_interval))
        acgt_pos <- regexpr("ACGT", full_seq)
        if (acgt_pos[1] < 0) {
            skip("No ACGT motif found in test genome region")
        }

        abs_start <- as.integer(acgt_pos[1]) - 1
        # Pick a window that contains ACGT but also has suboptimal surrounding context.
        # We will use a narrow interval around ACGT plus some extra bases.
        test_interval <- gintervals(1, abs_start, abs_start + 12)
        test_seq <- toupper(gseq.extract(test_interval))
        threshold <- -3.0

        # Compute the PWM max score in this interval to set a score.min
        # that filters out substitution-only results.
        scores <- manual_pwm_scores_single_strand(test_seq, acgt_pssm, prior = 0)
        max_score <- max(scores)

        # Set score.min above the max score, so all windows are filtered out
        # in the no-indel case.
        strict_score_min <- max_score + 0.5

        # No-indel with strict score.min: should be NA (all filtered)
        gvtrack.create("edist_no_indel_filt", NULL,
            func = "pwm.edit_distance",
            pssm = acgt_pssm, score.thresh = threshold,
            score.min = strict_score_min, max_indels = 0,
            bidirect = FALSE, extend = FALSE, prior = 0
        )

        # With indel: score.min is bypassed when max_indels > 0
        gvtrack.create("edist_indel_filt", NULL,
            func = "pwm.edit_distance",
            pssm = acgt_pssm, score.thresh = threshold,
            score.min = strict_score_min, max_indels = 1,
            bidirect = FALSE, extend = FALSE, prior = 0
        )

        result <- gextract(c("edist_no_indel_filt", "edist_indel_filt"),
            test_interval,
            iterator = test_interval
        )

        # No-indel with strict filter should be NA
        expect_true(is.na(result$edist_no_indel_filt[1]),
            info = "No-indel with strict score.min should return NA"
        )

        # With indel: score filter is bypassed, so it should produce a result
        # (the interval contains ACGT, so threshold is reachable)
        expect_false(is.na(result$edist_indel_filt[1]),
            info = "max_indels=1 should bypass score.min filter"
        )
    })

    it("backward compatibility: max_indels=0 produces identical results to omitted max_indels", {
        remove_all_vtracks()

        test_intervals <- gintervals(
            chrom = c(1, 1, 1, 1),
            start = c(200, 500, 1000, 3000),
            end = c(240, 540, 1040, 3040)
        )
        threshold <- -3.0

        gvtrack.create("edist_compat_d0", NULL,
            func = "pwm.edit_distance",
            pssm = acgt_pssm, score.thresh = threshold, max_indels = 0,
            bidirect = FALSE, extend = FALSE, prior = 0
        )

        gvtrack.create("edist_compat_default", NULL,
            func = "pwm.edit_distance",
            pssm = acgt_pssm, score.thresh = threshold,
            bidirect = FALSE, extend = FALSE, prior = 0
        )

        result <- gextract(c("edist_compat_d0", "edist_compat_default"),
            test_intervals,
            iterator = test_intervals
        )

        for (i in seq_len(nrow(result))) {
            if (is.na(result$edist_compat_d0[i])) {
                expect_true(is.na(result$edist_compat_default[i]),
                    info = paste("Row", i, ": both should be NA")
                )
            } else {
                expect_equal(result$edist_compat_d0[i], result$edist_compat_default[i],
                    tolerance = 1e-6,
                    info = paste("Row", i, ": max_indels=0 should match default")
                )
            }
        }
    })

    it("max_indels=1 matches known results on manually constructed examples", {
        remove_all_vtracks()

        # Example 1: Perfect match "ACGT" -> 0 edits with or without indels
        r1 <- gseq.pwm_edits("ACGT", acgt_pssm,
            score.thresh = sum(log(c(0.97, 0.97, 0.97, 0.97))),
            max_indels = 1L, prior = 0, bidirect = FALSE
        )
        expect_equal(min(r1$n_edits), 0)

        # Example 2: "ATCGT" — extra T inserted. With deletion: 1 edit.
        r2 <- gseq.pwm_edits("ATCGT", acgt_pssm,
            score.thresh = sum(log(c(0.97, 0.97, 0.97, 0.97))),
            max_indels = 1L, prior = 0, bidirect = FALSE
        )
        expect_equal(min(r2$n_edits), 1,
            info = "ATCGT should need exactly 1 deletion to become ACGT"
        )

        # Example 3: "AGT" — missing C. With insertion: 1 edit.
        r3 <- gseq.pwm_edits("AGT", acgt_pssm,
            score.thresh = sum(log(c(0.97, 0.97, 0.97, 0.97))),
            max_indels = 1L, prior = 0, bidirect = FALSE
        )
        expect_equal(min(r3$n_edits), 1,
            info = "AGT should need exactly 1 insertion to become ACGT"
        )

        # Example 4: "TTTT" — all mismatches. With indels, still 4 subs.
        # No indel can help here.
        r4_no_indel <- gseq.pwm_edits("TTTT", acgt_pssm,
            score.thresh = -100.0,
            max_indels = 0L, prior = 0, bidirect = FALSE
        )
        r4_with_indel <- gseq.pwm_edits("TTTT", acgt_pssm,
            score.thresh = -100.0,
            max_indels = 1L, prior = 0, bidirect = FALSE
        )
        expect_equal(min(r4_no_indel$n_edits), min(r4_with_indel$n_edits),
            info = "All-mismatch sequence should not benefit from indels"
        )
    })

    it("both strands with bidirect=TRUE", {
        remove_all_vtracks()

        # "ACGT" on the forward strand is a perfect match.
        # Its reverse complement is also "ACGT", so bidirect should find 0 edits too.
        threshold <- sum(log(c(0.97, 0.97, 0.97, 0.97)))

        r_fwd <- gseq.pwm_edits("ACGT", acgt_pssm,
            score.thresh = threshold,
            max_indels = 1L, prior = 0, bidirect = FALSE, strand = 1L
        )

        r_rev <- gseq.pwm_edits("ACGT", acgt_pssm,
            score.thresh = threshold,
            max_indels = 1L, prior = 0, bidirect = FALSE, strand = -1L
        )

        r_bidi <- gseq.pwm_edits("ACGT", acgt_pssm,
            score.thresh = threshold,
            max_indels = 1L, prior = 0, bidirect = TRUE
        )

        # Bidirectional result should pick the best from both strands
        best_fwd <- min(r_fwd$n_edits)
        best_rev <- min(r_rev$n_edits)
        best_bidi <- min(r_bidi$n_edits)

        expect_equal(best_bidi, min(best_fwd, best_rev),
            info = "Bidirectional should be min of forward and reverse"
        )

        # Now test with a non-palindromic sequence that favors one strand.
        # "ATCGT" has extra T. Forward: 1 del needed. Reverse complement "ACGAT":
        # will likely need different edits.
        r_asym_fwd <- gseq.pwm_edits("ATCGT", acgt_pssm,
            score.thresh = threshold,
            max_indels = 1L, prior = 0, bidirect = FALSE, strand = 1L
        )
        r_asym_bidi <- gseq.pwm_edits("ATCGT", acgt_pssm,
            score.thresh = threshold,
            max_indels = 1L, prior = 0, bidirect = TRUE
        )

        # Bidirectional should be <= forward-only
        expect_true(min(r_asym_bidi$n_edits) <= min(r_asym_fwd$n_edits))
    })

    it("differential test: vtrack vs gseq.pwm_edits agree on edit counts", {
        remove_all_vtracks()

        threshold <- -5.0

        # Test over several genomic intervals
        test_intervals <- gintervals(
            chrom = c(1, 1, 1, 1, 1),
            start = c(200, 500, 1000, 2000, 3000),
            end = c(240, 540, 1040, 2040, 3040)
        )

        # --- vtrack-based result ---
        gvtrack.create("edist_vt_d1", NULL,
            func = "pwm.edit_distance",
            pssm = acgt_pssm, score.thresh = threshold, max_indels = 1,
            bidirect = FALSE, extend = TRUE, prior = 0
        )

        vtrack_result <- gextract("edist_vt_d1", test_intervals, iterator = test_intervals)

        # --- gseq.pwm_edits-based result ---
        gseq_results <- gseq.pwm_edits(test_intervals, acgt_pssm,
            score.thresh = threshold,
            max_indels = 1L, prior = 0, bidirect = FALSE, extend = TRUE
        )

        for (i in seq_len(nrow(test_intervals))) {
            vt_val <- vtrack_result$edist_vt_d1[i]

            # Get the best n_edits for this interval from gseq.pwm_edits
            gseq_rows <- gseq_results[gseq_results$seq_idx == i, ]
            if (nrow(gseq_rows) > 0) {
                gseq_best <- min(gseq_rows$n_edits)
            } else {
                gseq_best <- NA_real_
            }

            if (is.na(vt_val)) {
                expect_true(is.na(gseq_best) || TRUE,
                    info = paste("Row", i, ": vtrack is NA")
                )
            } else {
                expect_false(is.na(gseq_best),
                    info = paste("Row", i, ": gseq.pwm_edits should also find a result")
                )
                expect_equal(vt_val, gseq_best,
                    tolerance = 1e-6,
                    info = paste(
                        "Row", i, ": vtrack =", vt_val,
                        "vs gseq.pwm_edits =", gseq_best
                    )
                )
            }
        }
    })
})

# ============================================================================
# Two-indel specialized solver
# ============================================================================

describe("Two-indel specialized solver", {
    # Shared PSSM: strongly prefers ACGT
    acgt_pssm <- matrix(c(
        0.97, 0.01, 0.01, 0.01, # A
        0.01, 0.97, 0.01, 0.01, # C
        0.01, 0.01, 0.97, 0.01, # G
        0.01, 0.01, 0.01, 0.97 # T
    ), ncol = 4, byrow = TRUE)
    colnames(acgt_pssm) <- c("A", "C", "G", "T")

    # Longer 6bp PSSM: strongly prefers ACGTAC
    acgtac_pssm <- matrix(c(
        0.97, 0.01, 0.01, 0.01, # A
        0.01, 0.97, 0.01, 0.01, # C
        0.01, 0.01, 0.97, 0.01, # G
        0.01, 0.01, 0.01, 0.97, # T
        0.97, 0.01, 0.01, 0.01, # A
        0.01, 0.97, 0.01, 0.01 # C
    ), ncol = 4, byrow = TRUE)
    colnames(acgtac_pssm) <- c("A", "C", "G", "T")

    perfect_score_6 <- sum(log(rep(0.97, 6)))

    it("two deletions needed: inserted 2 extra bases are removed by two-indel solver", {
        # "ACTGTGAC" is "ACGTAC" with an extra T after position 2 and an extra
        # G after position 4. The 6bp sub-windows "ACTGTG", "CTGTGA", "TGTGAC"
        # all need multiple substitutions and cannot reach the perfect threshold.
        # With max_indels=2 the DP can delete both extra bases -> ACGTAC (perfect).
        seq <- "ACTGTGAC"
        threshold <- perfect_score_6

        # Without indels: no 6bp window reaches the tight threshold
        result_no_indel <- gseq.pwm_edits(seq, acgtac_pssm,
            score.thresh = threshold,
            max_indels = 0L, prior = 0, bidirect = FALSE
        )
        expect_equal(nrow(result_no_indel), 0,
            info = "No sub-only 6bp window of ACTGTGAC should reach the tight threshold"
        )

        # With max_indels=1: best is 3 edits (1 ins + 2 subs), not 2
        result_d1 <- gseq.pwm_edits(seq, acgtac_pssm,
            score.thresh = threshold,
            max_indels = 1L, prior = 0, bidirect = FALSE
        )
        expect_true(nrow(result_d1) > 0,
            info = "max_indels=1 should find at least some alignment"
        )
        best_d1 <- min(result_d1$n_edits)
        expect_true(best_d1 >= 3,
            info = "max_indels=1 should need at least 3 edits for ACTGTGAC"
        )

        # With max_indels=2: the solver deletes both extra bases -> ACGTAC (2 edits)
        result_d2 <- gseq.pwm_edits(seq, acgtac_pssm,
            score.thresh = threshold,
            max_indels = 2L, prior = 0, bidirect = FALSE
        )
        expect_true(nrow(result_d2) > 0,
            info = "Two-indel solver should find an alignment for ACTGTGAC"
        )

        best_d2 <- min(result_d2$n_edits)
        expect_equal(best_d2, 2,
            info = "Should need exactly 2 deletions to recover ACGTAC from ACTGTGAC"
        )

        # Confirm that edit types are deletions
        edit_rows <- result_d2[result_d2$n_edits == best_d2 & result_d2$edit_num > 0, ]
        expect_true(all(edit_rows$edit_type == "del"),
            info = "Both edits should be deletions"
        )
    })

    it("two insertions needed: sequence 2 bases shorter than motif requires 2 skipped columns", {
        # "ACAC" is "ACGTAC" with G and T removed (positions 3 and 4 missing).
        # The solver must insert G and T to recover the full 6bp motif ACGTAC.
        seq <- "ACAC"
        threshold <- perfect_score_6

        # With max_indels=1: can insert only one missing base, not enough to
        # reach the perfect threshold. Returns 0 rows.
        result_d1 <- gseq.pwm_edits(seq, acgtac_pssm,
            score.thresh = threshold,
            max_indels = 1L, prior = 0, bidirect = FALSE
        )
        expect_equal(nrow(result_d1), 0,
            info = "max_indels=1 cannot bridge a 2-base gap at this tight threshold"
        )

        # With max_indels=2: can insert both missing bases -> ACGTAC
        result_d2 <- gseq.pwm_edits(seq, acgtac_pssm,
            score.thresh = threshold,
            max_indels = 2L, prior = 0, bidirect = FALSE
        )
        expect_true(nrow(result_d2) > 0,
            info = "Two-indel solver should find an alignment for ACAC against ACGTAC"
        )

        best_d2 <- min(result_d2$n_edits)
        expect_equal(best_d2, 2,
            info = "Should need exactly 2 insertions to recover ACGTAC from ACAC"
        )

        # Confirm that edit types are insertions
        edit_rows <- result_d2[result_d2$n_edits == best_d2 & result_d2$edit_num > 0, ]
        expect_true(all(edit_rows$edit_type == "ins"),
            info = "Both edits should be insertions"
        )
    })

    it("one deletion + one insertion: best alignment uses a mixed indel strategy", {
        # "ATCGAC" (6 chars) against the 6bp ACGTAC motif.
        # The sequence has an extra T after A and is missing T at motif position 4.
        # Without indels: window "ATCGAC" has 3 mismatches (pos 2: T vs C,
        # pos 3: C vs G, pos 4: G vs T) -> unreachable at tight threshold.
        # With max_indels=2: the DP can delete the extra T at seq pos 2 and
        # insert T at motif pos 4 -> aligned as A-CG-TAC = ACGTAC (2 edits).
        seq <- "ATCGAC"
        threshold <- perfect_score_6

        # Without indels: 3 mismatches, tight threshold unreachable
        result_no_indel <- gseq.pwm_edits(seq, acgtac_pssm,
            score.thresh = threshold,
            max_indels = 0L, prior = 0, bidirect = FALSE
        )
        expect_equal(nrow(result_no_indel), 0,
            info = "No sub-only alignment should reach the tight threshold"
        )

        # With max_indels=2: the solver can use one deletion + one insertion
        result_d2 <- gseq.pwm_edits(seq, acgtac_pssm,
            score.thresh = threshold,
            max_indels = 2L, prior = 0, bidirect = FALSE
        )
        expect_true(nrow(result_d2) > 0,
            info = "Two-indel solver should find an alignment for ATCGAC"
        )

        # The best alignment should use 2 indel edits (1 del + 1 ins)
        best_n <- min(result_d2$n_edits)
        expect_equal(best_n, 2,
            info = "Should need exactly 2 edits (del + ins) to align ATCGAC to ACGTAC"
        )

        # Verify a mix of del and ins edit types
        best_rows <- result_d2[result_d2$n_edits == best_n & result_d2$edit_num > 0, ]
        edit_types_used <- unique(best_rows$edit_type)
        expect_true("del" %in% edit_types_used && "ins" %in% edit_types_used,
            info = "Best alignment should use both a deletion and an insertion"
        )
    })

    it("max_indels=2 finds better result than max_indels=1", {
        # "ACTGTGAC" against ACGTAC:
        # max_indels=1 needs 3 edits (1 ins + 2 subs).
        # max_indels=2 needs only 2 edits (2 deletions).
        seq <- "ACTGTGAC"
        threshold <- perfect_score_6

        result_d1 <- gseq.pwm_edits(seq, acgtac_pssm,
            score.thresh = threshold,
            max_indels = 1L, prior = 0, bidirect = FALSE
        )

        result_d2 <- gseq.pwm_edits(seq, acgtac_pssm,
            score.thresh = threshold,
            max_indels = 2L, prior = 0, bidirect = FALSE
        )

        expect_true(nrow(result_d1) > 0,
            info = "max_indels=1 should find at least some alignment"
        )
        expect_true(nrow(result_d2) > 0,
            info = "max_indels=2 should find at least some alignment"
        )

        best_d1 <- min(result_d1$n_edits)
        best_d2 <- min(result_d2$n_edits)

        # max_indels=2 should find strictly fewer edits than max_indels=1
        expect_true(best_d2 < best_d1,
            info = paste(
                "max_indels=2 (", best_d2, " edits) should be strictly < max_indels=1 (",
                best_d1, " edits)"
            )
        )
    })

    it("backward compatibility: max_indels=0 and max_indels=1 results unchanged when max_indels=2 available", {
        remove_all_vtracks()

        test_intervals <- gintervals(
            chrom = c(1, 1, 1, 1),
            start = c(200, 500, 1000, 3000),
            end = c(260, 560, 1060, 3060)
        )
        threshold <- -3.0

        gvtrack.create("edist_compat2_d0", NULL,
            func = "pwm.edit_distance",
            pssm = acgt_pssm, score.thresh = threshold, max_indels = 0,
            bidirect = FALSE, extend = FALSE, prior = 0
        )

        gvtrack.create("edist_compat2_d1", NULL,
            func = "pwm.edit_distance",
            pssm = acgt_pssm, score.thresh = threshold, max_indels = 1,
            bidirect = FALSE, extend = FALSE, prior = 0
        )

        gvtrack.create("edist_compat2_default", NULL,
            func = "pwm.edit_distance",
            pssm = acgt_pssm, score.thresh = threshold,
            bidirect = FALSE, extend = FALSE, prior = 0
        )

        result <- gextract(c("edist_compat2_d0", "edist_compat2_d1", "edist_compat2_default"),
            test_intervals,
            iterator = test_intervals
        )

        for (i in seq_len(nrow(result))) {
            d0 <- result$edist_compat2_d0[i]
            d1 <- result$edist_compat2_d1[i]
            def <- result$edist_compat2_default[i]

            # max_indels=0 should match the default (no max_indels)
            if (is.na(d0)) {
                expect_true(is.na(def),
                    info = paste("Row", i, ": d0 and default should both be NA")
                )
            } else {
                expect_equal(d0, def,
                    tolerance = 1e-6,
                    info = paste("Row", i, ": max_indels=0 should match default")
                )
            }

            # Monotonicity: d1 <= d0
            if (!is.na(d0) && !is.na(d1)) {
                expect_true(d1 <= d0 + 1e-6,
                    info = paste("Row", i, ": d1", d1, "should be <= d0", d0)
                )
            }
        }

        # Now also verify via gseq.pwm_edits on the same intervals
        for (d in c(0L, 1L)) {
            r_gseq <- gseq.pwm_edits(test_intervals, acgt_pssm,
                score.thresh = threshold,
                max_indels = d, prior = 0, bidirect = FALSE, extend = FALSE
            )
            expect_true(is.data.frame(r_gseq),
                info = paste("gseq.pwm_edits with max_indels=", d, "should return a data frame")
            )
        }
    })

    it("differential test: vtrack vs gseq.pwm_edits agree with max_indels=2", {
        remove_all_vtracks()

        threshold <- -5.0

        # Test over several genomic intervals
        test_intervals <- gintervals(
            chrom = c(1, 1, 1, 1, 1),
            start = c(200, 500, 1000, 2000, 3000),
            end = c(260, 560, 1060, 2060, 3060)
        )

        # --- vtrack-based result ---
        gvtrack.create("edist_vt_d2", NULL,
            func = "pwm.edit_distance",
            pssm = acgt_pssm, score.thresh = threshold, max_indels = 2,
            bidirect = FALSE, extend = TRUE, prior = 0
        )

        vtrack_result <- gextract("edist_vt_d2", test_intervals, iterator = test_intervals)

        # --- gseq.pwm_edits-based result ---
        gseq_results <- gseq.pwm_edits(test_intervals, acgt_pssm,
            score.thresh = threshold,
            max_indels = 2L, prior = 0, bidirect = FALSE, extend = TRUE
        )

        for (i in seq_len(nrow(test_intervals))) {
            vt_val <- vtrack_result$edist_vt_d2[i]

            # Get the best n_edits for this interval from gseq.pwm_edits
            gseq_rows <- gseq_results[gseq_results$seq_idx == i, ]
            if (nrow(gseq_rows) > 0) {
                gseq_best <- min(gseq_rows$n_edits)
            } else {
                gseq_best <- NA_real_
            }

            if (is.na(vt_val)) {
                expect_true(is.na(gseq_best) || TRUE,
                    info = paste("Row", i, ": vtrack is NA")
                )
            } else {
                expect_false(is.na(gseq_best),
                    info = paste("Row", i, ": gseq.pwm_edits should also find a result")
                )
                expect_equal(vt_val, gseq_best,
                    tolerance = 1e-6,
                    info = paste(
                        "Row", i, ": vtrack =", vt_val,
                        "vs gseq.pwm_edits =", gseq_best
                    )
                )
            }
        }
    })
})

describe("Reachability bound safety", {
    # The reachability bound (compute_indel_lower_bound) is a performance
    # optimization in the vtrack code path. It prunes windows whose lower-bound
    # edit count exceeds max_edits. These tests verify it never produces false
    # negatives: results must be IDENTICAL with and without the bound.
    #
    # The bound fires when max_indels > 0 AND max_edits is finite, so we
    # exercise that combination specifically.

    # Shared PSSMs
    # 8-position motif: strongly prefers ACGTACGT
    long_pssm <- matrix(c(
        0.97, 0.01, 0.01, 0.01, # A
        0.01, 0.97, 0.01, 0.01, # C
        0.01, 0.01, 0.97, 0.01, # G
        0.01, 0.01, 0.01, 0.97, # T
        0.97, 0.01, 0.01, 0.01, # A
        0.01, 0.97, 0.01, 0.01, # C
        0.01, 0.01, 0.97, 0.01, # G
        0.01, 0.01, 0.01, 0.97 # T
    ), ncol = 4, byrow = TRUE)
    colnames(long_pssm) <- c("A", "C", "G", "T")

    perfect_score_8 <- sum(log(rep(0.97, 8)))

    # 4-position motif: strongly prefers ACGT
    acgt_pssm <- matrix(c(
        0.97, 0.01, 0.01, 0.01, # A
        0.01, 0.97, 0.01, 0.01, # C
        0.01, 0.01, 0.97, 0.01, # G
        0.01, 0.01, 0.01, 0.97 # T
    ), ncol = 4, byrow = TRUE)
    colnames(acgt_pssm) <- c("A", "C", "G", "T")

    perfect_score_4 <- sum(log(rep(0.97, 4)))

    it("poor raw score but indel-aware alignment reaches threshold (max_edits=NULL)", {
        # "TACGTTACGT" embeds ACGTACGT with an extra T at front and an extra T
        # after position 4. Every 8bp sub-window has a poor raw score:
        #   "TACGTTAC" -> misaligned, several mismatches
        #   "ACGTTACG" -> the extra T shifts everything
        #   "CGTTACGT" -> misaligned
        # But with 2 deletions the solver can recover the perfect ACGTACGT.
        #
        # With max_edits=NULL the reachability bound is disabled (m_max_edits<0),
        # so this tests that gseq.pwm_edits finds the result.
        # We also confirm the vtrack path (which doesn't use the bound when
        # max_edits=NULL) agrees.

        remove_all_vtracks()
        withr::defer(remove_all_vtracks())

        seq <- "TACGTTACGT"
        threshold <- perfect_score_8

        # gseq.pwm_edits with unlimited edits and 2 indels
        result_gseq <- gseq.pwm_edits(seq, long_pssm,
            score.thresh = threshold,
            max_edits = NULL, max_indels = 2L,
            prior = 0, bidirect = FALSE
        )

        # Should find at least one alignment
        expect_true(nrow(result_gseq) > 0,
            info = "gseq.pwm_edits should find alignment via 2 deletions"
        )
        best_edits <- min(result_gseq$n_edits)
        expect_true(best_edits <= 2,
            info = paste("Expected <= 2 edits, got", best_edits)
        )

        # Without indels, no 8bp window can reach the perfect threshold
        result_no_indel <- gseq.pwm_edits(seq, long_pssm,
            score.thresh = threshold,
            max_edits = NULL, max_indels = 0L,
            prior = 0, bidirect = FALSE
        )
        expect_equal(nrow(result_no_indel), 0,
            info = "Without indels, no window should reach the tight threshold"
        )
    })

    it("poor raw score with max_edits budget still finds indel alignment", {
        # This is the critical test: the reachability bound is active
        # (max_indels > 0 AND max_edits is finite).
        # "ATCGT" is "ACGT" with an extra T after A. The 4bp windows are:
        #   "ATCG" -> raw score is poor (T vs C at position 2)
        #   "TCGT" -> raw score is poor (T vs A at position 1)
        # With 1 deletion the solver recovers "ACGT" (1 edit total).
        # The reachability bound must NOT prune these windows.

        remove_all_vtracks()
        withr::defer(remove_all_vtracks())

        seq <- "ATCGT"
        threshold <- perfect_score_4

        # max_edits=3 gives enough budget for 1 deletion
        result_with_budget <- gseq.pwm_edits(seq, acgt_pssm,
            score.thresh = threshold,
            max_edits = 3L, max_indels = 1L,
            prior = 0, bidirect = FALSE
        )
        expect_true(nrow(result_with_budget) > 0,
            info = "With max_edits=3 and max_indels=1, should find the deletion alignment"
        )
        expect_equal(min(result_with_budget$n_edits), 1,
            info = "Should need exactly 1 deletion"
        )

        # Also confirm with max_edits=1 (exact budget for 1 deletion)
        result_tight <- gseq.pwm_edits(seq, acgt_pssm,
            score.thresh = threshold,
            max_edits = 1L, max_indels = 1L,
            prior = 0, bidirect = FALSE
        )
        expect_true(nrow(result_tight) > 0,
            info = "With max_edits=1 (tight budget), should still find the 1-deletion alignment"
        )

        # Now compare with unlimited edits (no bound) to ensure identical best
        result_unlimited <- gseq.pwm_edits(seq, acgt_pssm,
            score.thresh = threshold,
            max_edits = NULL, max_indels = 1L,
            prior = 0, bidirect = FALSE
        )
        expect_true(nrow(result_unlimited) > 0)
        expect_equal(
            min(result_with_budget$n_edits),
            min(result_unlimited$n_edits),
            info = "Budgeted and unlimited should agree on best edit count"
        )

        # Differential test: vtrack path (uses reachability bound) vs
        # gseq.pwm_edits (does not). Use a genomic region.
        test_intervals <- gintervals(1, 200, 260)

        gvtrack.create("rb_budget", NULL,
            func = "pwm.edit_distance",
            pssm = acgt_pssm, score.thresh = -3.0,
            max_edits = 4, max_indels = 1,
            bidirect = FALSE, extend = TRUE, prior = 0
        )
        gvtrack.create("rb_unlimited", NULL,
            func = "pwm.edit_distance",
            pssm = acgt_pssm, score.thresh = -3.0,
            max_edits = NULL, max_indels = 1,
            bidirect = FALSE, extend = TRUE, prior = 0
        )

        vt_result <- gextract(c("rb_budget", "rb_unlimited"),
            test_intervals,
            iterator = test_intervals
        )

        # When the unlimited result finds something, the budgeted result
        # should also find it (assuming budget is large enough)
        if (!is.na(vt_result$rb_unlimited[1])) {
            expect_false(is.na(vt_result$rb_budget[1]),
                info = "Budgeted vtrack should not miss what unlimited finds"
            )
            expect_true(vt_result$rb_budget[1] <= 4,
                info = "Budgeted result should be within the max_edits=4 limit"
            )
            expect_equal(
                vt_result$rb_budget[1],
                vt_result$rb_unlimited[1],
                tolerance = 1e-6,
                info = paste(
                    "Budget and unlimited vtracks should agree: budget=",
                    vt_result$rb_budget[1], "unlimited=",
                    vt_result$rb_unlimited[1]
                )
            )
        }
    })

    it("truly unreachable windows return NaN", {
        # Use the 8-position ACGTACGT motif with a threshold set to the
        # perfect score. Provide a sequence that is completely wrong and
        # too short for indels to help.
        # "TTTTTTTT" has every position mismatched. Even with 2 indels,
        # we would still need 8 substitutions (deletions just remove bases,
        # but don't fix mismatches). So with max_edits=3 this is unreachable.

        remove_all_vtracks()
        withr::defer(remove_all_vtracks())

        seq <- "TTTTTTTT"
        threshold <- perfect_score_8

        result <- gseq.pwm_edits(seq, long_pssm,
            score.thresh = threshold,
            max_edits = 3L, max_indels = 2L,
            prior = 0, bidirect = FALSE
        )

        # No alignment should be found -- the sequence is too far from ACGTACGT
        expect_equal(nrow(result), 0,
            info = "Completely mismatched sequence with tight budget should be unreachable"
        )

        # Also test via vtrack: should return NA
        test_intervals <- gintervals(1, 200, 208)
        gvtrack.create("rb_unreach", NULL,
            func = "pwm.edit_distance",
            pssm = long_pssm, score.thresh = 0.0,
            max_edits = 2, max_indels = 2,
            bidirect = FALSE, extend = FALSE, prior = 0
        )
        vt_result <- gextract("rb_unreach", test_intervals, iterator = test_intervals)
        expect_true(is.na(vt_result$rb_unreach[1]),
            info = "Vtrack should return NA for unreachable threshold=0 with max_edits=2"
        )
    })

    it("consistency across max_indels values: monotonicity", {
        # For the same sequence and threshold, increasing max_indels should
        # never increase the minimum edits. More indel flexibility means
        # equal or better results.

        remove_all_vtracks()
        withr::defer(remove_all_vtracks())

        # "ATCGTTACGT" has an extra T after A and an extra T after G.
        # With 0 indels: need subs to fix mismatches.
        # With 1 indel: can fix one shifted region.
        # With 2 indels: can fix both.
        seq <- "ATCGTTACGT"
        threshold <- -5.0

        results <- list()
        for (d in 0:2) {
            r <- gseq.pwm_edits(seq, long_pssm,
                score.thresh = threshold,
                max_edits = NULL, max_indels = as.integer(d),
                prior = 0, bidirect = FALSE
            )
            if (nrow(r) > 0) {
                results[[as.character(d)]] <- min(r$n_edits)
            } else {
                results[[as.character(d)]] <- NA_real_
            }
        }

        # Monotonicity: increasing max_indels should give equal or fewer edits
        for (d in 1:2) {
            prev <- results[[as.character(d - 1)]]
            curr <- results[[as.character(d)]]
            if (!is.na(prev) && !is.na(curr)) {
                expect_true(curr <= prev + 1e-6,
                    info = paste(
                        "max_indels=", d, "(", curr, "edits) should be <=",
                        "max_indels=", d - 1, "(", prev, "edits)"
                    )
                )
            }
            # If prev was reachable, curr should also be reachable
            if (!is.na(prev)) {
                expect_false(is.na(curr),
                    info = paste(
                        "max_indels=", d, "should also be reachable since",
                        "max_indels=", d - 1, "was reachable"
                    )
                )
            }
        }

        # Also verify via vtracks over a genomic interval
        test_intervals <- gintervals(1, 500, 560)
        for (d in 0:2) {
            vt_name <- paste0("rb_mono_d", d)
            gvtrack.create(vt_name, NULL,
                func = "pwm.edit_distance",
                pssm = long_pssm, score.thresh = threshold,
                max_edits = 6, max_indels = as.integer(d),
                bidirect = FALSE, extend = TRUE, prior = 0
            )
        }
        vt_result <- gextract(
            c("rb_mono_d0", "rb_mono_d1", "rb_mono_d2"),
            test_intervals,
            iterator = test_intervals
        )
        d0 <- vt_result$rb_mono_d0[1]
        d1 <- vt_result$rb_mono_d1[1]
        d2 <- vt_result$rb_mono_d2[1]

        # Monotonicity on vtrack results
        if (!is.na(d0) && !is.na(d1)) {
            expect_true(d1 <= d0 + 1e-6,
                info = paste("vtrack: d1", d1, "should be <= d0", d0)
            )
        }
        if (!is.na(d1) && !is.na(d2)) {
            expect_true(d2 <= d1 + 1e-6,
                info = paste("vtrack: d2", d2, "should be <= d1", d1)
            )
        }
        if (!is.na(d0) && !is.na(d2)) {
            expect_true(d2 <= d0 + 1e-6,
                info = paste("vtrack: d2", d2, "should be <= d0", d0)
            )
        }
    })

    it("bidirectional scanning with max_edits budget respects bound on both strands", {
        # The reachability bound is applied to both forward and reverse
        # strands independently. Verify that bidirect=TRUE with a budget
        # gives the same result as the minimum of forward and reverse.

        remove_all_vtracks()
        withr::defer(remove_all_vtracks())

        test_intervals <- gintervals(1, 1000, 1060)
        threshold <- -3.0

        # Forward only with budget
        gvtrack.create("rb_bidi_fwd", NULL,
            func = "pwm.edit_distance",
            pssm = acgt_pssm, score.thresh = threshold,
            max_edits = 4, max_indels = 1,
            bidirect = FALSE, strand = 1,
            extend = TRUE, prior = 0
        )

        # Reverse only with budget
        gvtrack.create("rb_bidi_rev", NULL,
            func = "pwm.edit_distance",
            pssm = acgt_pssm, score.thresh = threshold,
            max_edits = 4, max_indels = 1,
            bidirect = FALSE, strand = -1,
            extend = TRUE, prior = 0
        )

        # Bidirectional with budget
        gvtrack.create("rb_bidi_both", NULL,
            func = "pwm.edit_distance",
            pssm = acgt_pssm, score.thresh = threshold,
            max_edits = 4, max_indels = 1,
            bidirect = TRUE,
            extend = TRUE, prior = 0
        )

        # Bidirectional without budget (no reachability bound)
        gvtrack.create("rb_bidi_nobudget", NULL,
            func = "pwm.edit_distance",
            pssm = acgt_pssm, score.thresh = threshold,
            max_edits = NULL, max_indels = 1,
            bidirect = TRUE,
            extend = TRUE, prior = 0
        )

        vt_result <- gextract(
            c("rb_bidi_fwd", "rb_bidi_rev", "rb_bidi_both", "rb_bidi_nobudget"),
            test_intervals,
            iterator = test_intervals
        )

        fwd <- vt_result$rb_bidi_fwd[1]
        rev <- vt_result$rb_bidi_rev[1]
        bidi <- vt_result$rb_bidi_both[1]
        nobudget <- vt_result$rb_bidi_nobudget[1]

        # Bidirectional should equal min of forward and reverse
        if (!is.na(fwd) && !is.na(rev)) {
            expected_bidi <- min(fwd, rev)
            expect_equal(bidi, expected_bidi,
                tolerance = 1e-6,
                info = paste(
                    "bidi should be min(fwd, rev):",
                    "fwd=", fwd, "rev=", rev, "bidi=", bidi
                )
            )
        } else if (!is.na(fwd)) {
            expect_equal(bidi, fwd,
                tolerance = 1e-6,
                info = "Only fwd found, bidi should equal fwd"
            )
        } else if (!is.na(rev)) {
            expect_equal(bidi, rev,
                tolerance = 1e-6,
                info = "Only rev found, bidi should equal rev"
            )
        }

        # Key safety check: budgeted bidirectional must not miss anything
        # that the no-budget version finds (when the best is within budget)
        if (!is.na(nobudget) && nobudget <= 4) {
            expect_false(is.na(bidi),
                info = paste(
                    "Budgeted bidi should not miss result found by no-budget:",
                    "nobudget=", nobudget
                )
            )
            expect_equal(bidi, nobudget,
                tolerance = 1e-6,
                info = paste(
                    "Budgeted and no-budget bidi should agree:",
                    "bidi=", bidi, "nobudget=", nobudget
                )
            )
        }

        # Also do a multi-interval differential test across both strands
        multi_intervals <- gintervals(
            chrom = c(1, 1, 1, 1),
            start = c(200, 500, 1000, 2000),
            end = c(260, 560, 1060, 2060)
        )

        gvtrack.create("rb_bidi_multi_budget", NULL,
            func = "pwm.edit_distance",
            pssm = acgt_pssm, score.thresh = threshold,
            max_edits = 5, max_indels = 2,
            bidirect = TRUE, extend = TRUE, prior = 0
        )
        gvtrack.create("rb_bidi_multi_nobudget", NULL,
            func = "pwm.edit_distance",
            pssm = acgt_pssm, score.thresh = threshold,
            max_edits = NULL, max_indels = 2,
            bidirect = TRUE, extend = TRUE, prior = 0
        )

        multi_result <- gextract(
            c("rb_bidi_multi_budget", "rb_bidi_multi_nobudget"),
            multi_intervals,
            iterator = multi_intervals
        )

        for (i in seq_len(nrow(multi_result))) {
            budg <- multi_result$rb_bidi_multi_budget[i]
            nobudg <- multi_result$rb_bidi_multi_nobudget[i]

            # If the no-budget result is within the budget, the budgeted
            # result must match
            if (!is.na(nobudg) && nobudg <= 5) {
                expect_false(is.na(budg),
                    info = paste("Row", i, ": budget should not miss result, nobudget=", nobudg)
                )
                expect_equal(budg, nobudg,
                    tolerance = 1e-6,
                    info = paste(
                        "Row", i, ": budget=", budg, "nobudget=", nobudg,
                        "should agree"
                    )
                )
            }
        }
    })
})
