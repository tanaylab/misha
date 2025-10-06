test_that("pwm.count counts hits above threshold", {
    remove_all_vtracks()

    pssm <- create_test_pssm() # AC motif

    test_intervals <- gintervals(1, 200, 240)
    seq <- toupper(gseq.extract(test_intervals)) # CCCTAACCCTAACCCTAACCCTAACCCTAACCCCTAACCC

    # Count AC occurrences manually
    ac_positions <- stringr::str_locate_all(seq, "AC")[[1]][, 1]
    expected_count <- length(ac_positions)

    # Create pwm.count virtual track with threshold 0
    gvtrack.create("count_hits", NULL,
        func = "pwm.count",
        pssm = pssm, bidirect = FALSE, extend = FALSE,
        prior = 0, score.thresh = 0
    )

    result <- gextract("count_hits", test_intervals, iterator = test_intervals)

    # With prior=0 and score.thresh=0, perfect matches should score 0 (log(1) = 0)
    # All other kmers should score -Inf
    expect_equal(result$count_hits[1], expected_count, ignore_attr = TRUE)
})

test_that("pwm.count respects score threshold", {
    remove_all_vtracks()

    pssm <- create_test_pssm() # AC motif

    test_intervals <- gintervals(1, 200, 240)

    # Create two tracks with different thresholds
    gvtrack.create("count_all", NULL,
        func = "pwm.count",
        pssm = pssm, bidirect = FALSE, extend = FALSE,
        prior = 0.01, score.thresh = -10
    )

    gvtrack.create("count_strict", NULL,
        func = "pwm.count",
        pssm = pssm, bidirect = FALSE, extend = FALSE,
        prior = 0.01, score.thresh = -1
    )

    result <- gextract(c("count_all", "count_strict"), test_intervals, iterator = test_intervals)

    # With lower threshold, we should count more hits
    expect_true(result$count_all[1] >= result$count_strict[1])

    # Both should be non-negative integers
    expect_true(result$count_all[1] >= 0)
    expect_true(result$count_strict[1] >= 0)
})

test_that("pwm.count works with bidirect=TRUE", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    pssm <- create_test_pssm() # AC motif

    test_interval <- gintervals(1, 200, 240)

    # Create tracks for forward, bidirectional
    gvtrack.create("count_fwd", NULL, "pwm.count",
        pssm = pssm, bidirect = FALSE, extend = FALSE,
        prior = 0.01, score.thresh = -5
    )

    gvtrack.create("count_bidi", NULL, "pwm.count",
        pssm = pssm, bidirect = TRUE, extend = FALSE,
        prior = 0.01, score.thresh = -5
    )

    result <- gextract(c("count_fwd", "count_bidi"), test_interval, iterator = test_interval)

    # Bidirectional should count at least as many as forward-only
    expect_true(result$count_bidi[1] >= result$count_fwd[1])
})

test_that("pwm.count honors extend parameter", {
    remove_all_vtracks()

    pssm <- create_test_pssm() # AC motif (length 2)

    # Short interval where extension matters
    test_intervals <- gintervals(1, 200, 210)

    gvtrack.create("count_noext", NULL, "pwm.count",
        pssm = pssm, bidirect = FALSE, extend = FALSE,
        prior = 0.01, score.thresh = -5
    )

    gvtrack.create("count_ext", NULL, "pwm.count",
        pssm = pssm, bidirect = FALSE, extend = TRUE,
        prior = 0.01, score.thresh = -5
    )

    result <- gextract(c("count_noext", "count_ext"), test_intervals, iterator = test_intervals)

    # With extend=TRUE, we might count one more position at the boundary
    expect_true(result$count_ext[1] >= result$count_noext[1])
})

test_that("pwm.count works with iterator shifts", {
    remove_all_vtracks()

    pssm <- create_test_pssm() # AC motif

    base <- gintervals(1, 2000, 2040)

    # Create track with iterator shifts
    gvtrack.create("count_shift", NULL, "pwm.count",
        pssm = pssm, bidirect = FALSE, extend = TRUE,
        prior = 0.01, score.thresh = -5
    )

    gvtrack.iterator("count_shift", sshift = -10, eshift = 10)

    result <- gextract("count_shift", base, iterator = base)

    # Should return non-negative count
    expect_true(result$count_shift[1] >= 0)
})

test_that("pwm.count returns 0 for very high threshold", {
    remove_all_vtracks()

    pssm <- create_test_pssm()

    test_intervals <- gintervals(1, 200, 240)

    # With very high threshold, nothing should pass
    gvtrack.create("count_impossible", NULL, "pwm.count",
        pssm = pssm, bidirect = FALSE, extend = FALSE,
        prior = 0.01, score.thresh = 100
    )

    result <- gextract("count_impossible", test_intervals, iterator = test_intervals)

    expect_equal(result$count_impossible[1], 0, ignore_attr = TRUE)
})

test_that("pwm.count works with strand parameter", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    pssm <- create_test_pssm()

    test_interval <- gintervals(1, 200, 240)

    # Create tracks for each strand
    gvtrack.create("count_plus", NULL, "pwm.count",
        pssm = pssm, bidirect = FALSE, extend = TRUE,
        prior = 0.01, score.thresh = -5, strand = 1
    )

    gvtrack.create("count_minus", NULL, "pwm.count",
        pssm = pssm, bidirect = FALSE, extend = TRUE,
        prior = 0.01, score.thresh = -5, strand = -1
    )

    result <- gextract(c("count_plus", "count_minus"), test_interval, iterator = test_interval)

    # Both should return non-negative counts
    expect_true(result$count_plus[1] >= 0)
    expect_true(result$count_minus[1] >= 0)
})

test_that("pwm.count matches manual counting for perfect matches", {
    remove_all_vtracks()

    # Create a simple pssm where only one specific sequence has score 0
    pssm <- create_test_pssm()

    test_intervals <- gintervals(1, 200, 240)
    seq <- toupper(gseq.extract(test_intervals))

    # Count "AC" occurrences
    ac_count <- stringr::str_count(seq, "AC")

    gvtrack.create("count_exact", NULL, "pwm.count",
        pssm = pssm, bidirect = FALSE, extend = FALSE,
        prior = 0, score.thresh = 0
    )

    result <- gextract("count_exact", test_intervals, iterator = test_intervals)

    # Should match manual count
    expect_equal(result$count_exact[1], ac_count, ignore_attr = TRUE)
})


test_that("pwm.count bidirectional equals sum of strand-specific counts", {
    remove_all_vtracks()

    pssm <- create_test_pssm()
    test_interval <- gintervals(1, 200, 300)

    gvtrack.create("count_plus", NULL, "pwm.count",
        pssm = pssm, bidirect = FALSE, extend = TRUE,
        prior = 0.01, score.thresh = -10, strand = 1
    )

    gvtrack.create("count_minus", NULL, "pwm.count",
        pssm = pssm, bidirect = FALSE, extend = TRUE,
        prior = 0.01, score.thresh = -10, strand = -1
    )

    gvtrack.create("count_bidi", NULL, "pwm.count",
        pssm = pssm, bidirect = TRUE, extend = TRUE,
        prior = 0.01, score.thresh = -10
    )

    result <- gextract(c("count_plus", "count_minus", "count_bidi"), test_interval, iterator = test_interval)

    expect_equal(result$count_bidi[1], result$count_plus[1] + result$count_minus[1], tolerance = 1e-8)
})

test_that("pwm.count spatial weighting can increase counts over non-spatial at positive threshold", {
    remove_all_vtracks()

    pssm <- create_test_pssm() # AC motif; with prior=0, perfect match score is 0
    test_interval <- gintervals(1, 200, 300)

    # Non-spatial with positive threshold: with prior=0 and score.thresh=log(2), forward logp=0 won't pass
    gvtrack.create("count_nospatial", NULL, "pwm.count",
        pssm = pssm, bidirect = FALSE, extend = TRUE,
        prior = 0, score.thresh = log(2)
    )

    # Spatial: positions in high-weight bins should pass (0 + log(2) >= log(2))
    gvtrack.create("count_spatial", NULL, "pwm.count",
        pssm = pssm, bidirect = FALSE, extend = TRUE,
        prior = 0, score.thresh = log(2),
        spat_factor = c(0.5, 1.0, 2.0, 1.0, 0.5),
        spat_bin = 20L
    )

    res <- gextract(c("count_nospatial", "count_spatial"), test_interval, iterator = test_interval)
    expect_true(res$count_spatial[1] >= res$count_nospatial[1])
})

test_that("pwm.count rejects non-positive spatial weights", {
    remove_all_vtracks()

    pssm <- create_test_pssm()

    expect_error(
        gvtrack.create("count_bad_spat", NULL, "pwm.count",
            pssm = pssm, bidirect = TRUE, extend = TRUE,
            spat_factor = c(1.0, 0.0, 1.0), spat_bin = 10L
        ),
        regexp = "spat_factor must be a numeric vector with all positive values"
    )
})

test_that("pwm.count prior reduces counts when threshold is fixed", {
    remove_all_vtracks()

    pssm <- create_test_pssm()
    test_interval <- gintervals(1, 200, 300)

    # Prior=0: perfect matches have logp=0 and pass threshold 0
    gvtrack.create("count_prior0", NULL, "pwm.count",
        pssm = pssm, bidirect = FALSE, extend = FALSE,
        prior = 0, score.thresh = 0
    )

    # Prior>0: perfect matches have logp<0, won't pass threshold 0
    gvtrack.create("count_prior01", NULL, "pwm.count",
        pssm = pssm, bidirect = FALSE, extend = FALSE,
        prior = 0.1, score.thresh = 0
    )

    out <- gextract(c("count_prior0", "count_prior01"), test_interval, iterator = test_interval)
    expect_true(out$count_prior0[1] >= out$count_prior01[1])
})

test_that("pwm.count batch path (many vtracks) returns consistent values", {
    remove_all_vtracks()

    pssm <- create_test_pssm()
    test_interval <- gintervals(1, 200, 300)

    # Create 5 identical pwm.count vtracks to trigger batch processing path
    for (i in 1:5) {
        gvtrack.create(paste0("count_rep_", i), NULL, "pwm.count",
            pssm = pssm, bidirect = FALSE, extend = TRUE,
            prior = 0.01, score.thresh = -5
        )
    }

    vnames <- paste0("count_rep_", 1:5)
    res <- gextract(vnames, test_interval, iterator = test_interval)

    # All columns should be equal (same inputs) and non-negative
    base <- res[[vnames[1]]][1]
    expect_true(base >= 0)
    for (nm in vnames[-1]) {
        expect_equal(res[[nm]][1], base, ignore_attr = TRUE)
    }
})

test_that("pwm.count shift equivalence: shifted vtrack over base equals unshifted over expanded iterator", {
    remove_all_vtracks()

    pssm <- create_test_pssm()

    base60 <- gintervals(1, 2100, 2160)
    base80 <- gintervals(1, 2090, 2170)

    gvtrack.create("pwmcount_shifted", NULL,
        func = "pwm.count",
        pssm = pssm, bidirect = FALSE, extend = TRUE, prior = 0.01, score.thresh = -5
    )
    gvtrack.iterator("pwmcount_shifted", sshift = -10, eshift = 10)

    s_shift <- gextract("pwmcount_shifted", base60, iterator = base60)

    gvtrack.create("pwmcount_unshifted", NULL,
        func = "pwm.count",
        pssm = pssm, bidirect = FALSE, extend = TRUE, prior = 0.01, score.thresh = -5
    )
    s_unshift <- gextract("pwmcount_unshifted", base80, iterator = base80)

    expect_equal(s_shift$pwmcount_shifted[1], s_unshift$pwmcount_unshifted[1], ignore_attr = TRUE)
})

test_that("pwm.count spatial + iterator shifts equivalence holds", {
    remove_all_vtracks()

    pssm <- create_test_pssm()

    base60 <- gintervals(1, 2200, 2260)
    base80 <- gintervals(1, 2190, 2270)

    gvtrack.create("pwmcount_spat_shifted", NULL,
        func = "pwm.count",
        pssm = pssm, bidirect = FALSE, extend = TRUE, prior = 0.01, score.thresh = -5,
        spat_factor = c(0.5, 1.0, 2.0, 1.0, 0.5),
        spat_bin = 20L
    )
    gvtrack.iterator("pwmcount_spat_shifted", sshift = -10, eshift = 10)
    a_shift <- gextract("pwmcount_spat_shifted", base60, iterator = base60)

    gvtrack.create("pwmcount_spat_unshifted", NULL,
        func = "pwm.count",
        pssm = pssm, bidirect = FALSE, extend = TRUE, prior = 0.01, score.thresh = -5,
        spat_factor = c(0.5, 1.0, 2.0, 1.0, 0.5),
        spat_bin = 20L
    )
    a_unshift <- gextract("pwmcount_spat_unshifted", base80, iterator = base80)

    expect_equal(a_shift$pwmcount_spat_shifted[1], a_unshift$pwmcount_spat_unshifted[1], ignore_attr = TRUE)
})
