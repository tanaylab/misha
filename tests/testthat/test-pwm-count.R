create_isolated_test_db()

# The PSSM these tests use, create_test_pssm() with prior = 0.01, is the 2-row AC motif, so
# a position can only score -0.0585 (AC), -4.6737 (one of the two bases matches) or -9.2888
# (neither). A threshold at or below -9.2888 - which is what score.thresh = -10 was - passes
# every position, so the count saturates at the interval length and the comparisons below
# hold whatever the counting path does: over chr1:200-300 the forward, reverse and
# bidirectional counts were all 100 out of 100, and "bidirectional equals the union of the
# strands" was true because all three were the ceiling.
#
# -5 is the only cut that separates all three score levels, and it is the one that leaves the
# three counts distinct and off the ceiling: 70 forward, 14 reverse, 84 bidirectional.
pwm_count_thresh <- -5

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
        prior = 0.01, score.thresh = pwm_count_thresh
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


test_that("pwm.count bidirectional equals per-position union (LSE) of strands", {
    remove_all_vtracks()

    pssm <- create_test_pssm()
    test_interval <- gintervals(1, 200, 300)

    gvtrack.create("count_plus", NULL, "pwm.count",
        pssm = pssm, bidirect = FALSE, extend = TRUE,
        prior = 0.01, score.thresh = pwm_count_thresh, strand = 1
    )

    gvtrack.create("count_minus", NULL, "pwm.count",
        pssm = pssm, bidirect = FALSE, extend = TRUE,
        prior = 0.01, score.thresh = pwm_count_thresh, strand = -1
    )

    gvtrack.create("count_bidi", NULL, "pwm.count",
        pssm = pssm, bidirect = TRUE, extend = TRUE,
        prior = 0.01, score.thresh = pwm_count_thresh
    )

    gvtrack.create("pwm_plus", NULL, "pwm",
        pssm = pssm, bidirect = FALSE, extend = TRUE,
        prior = 0.01, score.thresh = pwm_count_thresh, strand = 1
    )
    gvtrack.create("pwm_minus", NULL, "pwm",
        pssm = pssm, bidirect = FALSE, extend = TRUE,
        prior = 0.01, score.thresh = pwm_count_thresh, strand = -1
    )
    gvtrack.create("pwm_bidi", NULL, "pwm",
        pssm = pssm, bidirect = TRUE, extend = TRUE,
        prior = 0.01, score.thresh = pwm_count_thresh
    )
    pwm_result <- gextract(c("pwm_plus", "pwm_minus", "pwm_bidi"), test_interval, iterator = 1)
    n_pwm_plus <- sum(pwm_result$pwm_plus > pwm_count_thresh)
    n_pwm_minus <- sum(pwm_result$pwm_minus > pwm_count_thresh)
    # pwm_bidi is already LSE-combined per position; threshold it:
    n_pwm_bidi <- sum(pwm_result$pwm_bidi > pwm_count_thresh)

    result <- gextract(c("count_plus", "count_minus", "count_bidi"), test_interval, iterator = test_interval)

    expect_equal(result$count_plus[1], n_pwm_plus)
    expect_equal(result$count_minus[1], n_pwm_minus)
    # New union semantics: bidi count equals number of positions passing LSE-combined pwm
    expect_equal(result$count_bidi[1], n_pwm_bidi)
})

test_that("pwm.count: bidi equals union (LSE) and matches pwm thresholding", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    pssm <- create_test_pssm()
    test_interval <- gintervals(1, 200, 300)

    # Per-position PWM scores on each strand and bidirectional
    gvtrack.create("pwm_plus", NULL, "pwm",
        pssm = pssm, bidirect = FALSE, extend = TRUE,
        prior = 0.01, score.thresh = pwm_count_thresh, strand = 1
    )
    gvtrack.create("pwm_minus", NULL, "pwm",
        pssm = pssm, bidirect = FALSE, extend = TRUE,
        prior = 0.01, score.thresh = pwm_count_thresh, strand = -1
    )
    gvtrack.create("pwm_bidi", NULL, "pwm",
        pssm = pssm, bidirect = TRUE, extend = TRUE,
        prior = 0.01, score.thresh = pwm_count_thresh
    )

    # Iterator=1 gives per-base positions to threshold against
    pwm_result <- gextract(c("pwm_plus", "pwm_minus", "pwm_bidi"),
        test_interval,
        iterator = 1
    )
    n_pwm_plus <- sum(pwm_result$pwm_plus > pwm_count_thresh, na.rm = TRUE)
    n_pwm_minus <- sum(pwm_result$pwm_minus > pwm_count_thresh, na.rm = TRUE)
    n_pwm_bidi <- sum(pwm_result$pwm_bidi > pwm_count_thresh, na.rm = TRUE)

    # Strand-specific and bidirectional counts
    gvtrack.create("count_plus", NULL, "pwm.count",
        pssm = pssm, bidirect = FALSE, extend = TRUE,
        prior = 0.01, score.thresh = pwm_count_thresh, strand = 1
    )
    gvtrack.create("count_minus", NULL, "pwm.count",
        pssm = pssm, bidirect = FALSE, extend = TRUE,
        prior = 0.01, score.thresh = pwm_count_thresh, strand = -1
    )
    gvtrack.create("count_bidi", NULL, "pwm.count",
        pssm = pssm, bidirect = TRUE, extend = TRUE,
        prior = 0.01, score.thresh = pwm_count_thresh
    )

    result <- gextract(c("count_plus", "count_minus", "count_bidi"),
        test_interval,
        iterator = test_interval
    )

    # Strand-specific counts still match strand-specific PWM thresholding
    expect_equal(result$count_plus[1], n_pwm_plus, ignore_attr = TRUE)
    expect_equal(result$count_minus[1], n_pwm_minus, ignore_attr = TRUE)
    # New union semantics: bidi equals LSE-combined pwm thresholding per position
    expect_equal(result$count_bidi[1], n_pwm_bidi, ignore_attr = TRUE)
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
            score.thresh = 0,
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

test_that("pwm.count: strand=-1 matches pwm_minus (non-spatial, sliding path)", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    pssm <- create_test_pssm()
    test_interval <- gintervals(1, 200, 300)

    # Per-position minus-strand scores (reference)
    gvtrack.create("pwm_minus", NULL, "pwm",
        pssm = pssm, bidirect = FALSE, extend = TRUE,
        prior = 0.01, score.thresh = pwm_count_thresh, strand = -1
    )
    pwm <- gextract("pwm_minus", test_interval, iterator = 1)
    n_minus <- sum(pwm$pwm_minus > pwm_count_thresh, na.rm = TRUE)

    # Minus-only count (uses sliding path when non-spatial)
    gvtrack.create("count_minus", NULL, "pwm.count",
        pssm = pssm, bidirect = FALSE, extend = TRUE,
        prior = 0.01, score.thresh = pwm_count_thresh, strand = -1
    )
    res <- gextract("count_minus", test_interval, iterator = test_interval)

    expect_equal(res$count_minus[1], n_minus, ignore_attr = TRUE)
})

test_that("pwm.count: strand=-1 matches pwm_minus (spatial, non-sliding path; weights=1)", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    pssm <- create_test_pssm()
    test_interval <- gintervals(1, 200, 300)

    # Per-position minus-strand scores (reference)
    gvtrack.create("pwm_minus", NULL, "pwm",
        pssm = pssm, bidirect = FALSE, extend = TRUE,
        prior = 0.01, score.thresh = pwm_count_thresh, strand = -1
    )
    pwm <- gextract("pwm_minus", test_interval, iterator = 1)
    n_minus <- sum(pwm$pwm_minus > pwm_count_thresh, na.rm = TRUE)

    # Spatial path forced by spat_factor; weights=1 -> no effect on threshold
    gvtrack.create("count_minus_spat", NULL, "pwm.count",
        pssm = pssm, bidirect = FALSE, extend = TRUE,
        prior = 0.01, score.thresh = pwm_count_thresh, strand = -1,
        spat_factor = rep(1.0, 5), spat_bin = 20L
    )
    res <- gextract("count_minus_spat", test_interval, iterator = test_interval)

    expect_equal(res$count_minus_spat[1], n_minus, ignore_attr = TRUE)
})

test_that("pwm.count: bidirect ignores strand parameter (union semantics)", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    pssm <- create_test_pssm()
    test_interval <- gintervals(1, 200, 300)

    gvtrack.create("count_bidi_s1", NULL, "pwm.count",
        pssm = pssm, bidirect = TRUE, extend = TRUE,
        prior = 0.01, score.thresh = pwm_count_thresh, strand = 1
    )
    gvtrack.create("count_bidi_sneg1", NULL, "pwm.count",
        pssm = pssm, bidirect = TRUE, extend = TRUE,
        prior = 0.01, score.thresh = pwm_count_thresh, strand = -1
    )

    out <- gextract(c("count_bidi_s1", "count_bidi_sneg1"),
        test_interval,
        iterator = test_interval
    )

    expect_equal(out$count_bidi_s1[1], out$count_bidi_sneg1[1], ignore_attr = TRUE)
})

# ---------------------------------------------------------------------------
# score.thresh is mandatory for pwm.count
#
# pwm.count used to default score.thresh to 0. A misha PWM score is a plain
# log-likelihood (a sum of log p over the PSSM rows), so for any PSSM that is
# not perfectly deterministic every window scores below 0 and the default
# counted nothing at all - silently, everywhere. There is no value that is
# right for every PSSM, so the parameter is required, exactly as the
# pwm.edit_distance family already requires it.
# ---------------------------------------------------------------------------

test_that("pwm.count requires score.thresh (named arguments)", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    pssm <- create_test_pssm()

    expect_error(
        gvtrack.create("count_no_thresh", NULL, "pwm.count", pssm = pssm),
        "score.thresh"
    )
    expect_false("count_no_thresh" %in% gvtrack.ls())
})

test_that("pwm.count requires score.thresh (params list)", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    pssm <- create_test_pssm()

    expect_error(
        gvtrack.create("count_no_thresh_params", NULL, "pwm.count",
            params = list(pssm = pssm, bidirect = FALSE)
        ),
        "score.thresh"
    )
    expect_false("count_no_thresh_params" %in% gvtrack.ls())
})

test_that("pwm.count coerces character/factor score.thresh, rejects vector and logical", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    pssm <- create_test_pssm()
    test_interval <- gintervals(1, 200, 300)

    # A threshold read out of a config file or a read.csv column arrives as a
    # character, or as a factor if stringsAsFactors was left on. Both are
    # coerced, and must give the same answer as the number.
    gvtrack.create("count_chr", NULL, "pwm.count",
        pssm = pssm, bidirect = FALSE, extend = TRUE, prior = 0.01, score.thresh = "-5"
    )
    gvtrack.create("count_fct", NULL, "pwm.count",
        pssm = pssm, bidirect = FALSE, extend = TRUE, prior = 0.01, score.thresh = factor("-5")
    )
    gvtrack.create("count_num", NULL, "pwm.count",
        pssm = pssm, bidirect = FALSE, extend = TRUE, prior = 0.01, score.thresh = -5
    )
    res <- gextract(c("count_chr", "count_fct", "count_num"), test_interval, iterator = test_interval)
    expect_gt(res$count_num[1], 0)
    expect_equal(res$count_chr[1], res$count_num[1], ignore_attr = TRUE)
    # as.numeric() on a factor returns the level index, which for factor("-5")
    # is 1 - a threshold this PSSM never reaches, i.e. a silent zero.
    expect_equal(res$count_fct[1], res$count_num[1], ignore_attr = TRUE)

    # A vector would be silently truncated to its first element by the C++
    # layer, which is the failure mode this whole rule exists to close.
    expect_error(
        gvtrack.create("count_bad_thresh", NULL, "pwm.count", pssm = pssm, score.thresh = c(-10, -5)),
        "score.thresh must be a single value"
    )

    # Something that is not a number at all is still an error.
    expect_error(
        gvtrack.create("count_bad_thresh", NULL, "pwm.count", pssm = pssm, score.thresh = "loose"),
        "score.thresh must be a single number"
    )

    # A logical coerces to 1 or 0 under as.numeric(). Both are unreachable for
    # this PSSM, so accepting one would silently count nothing forever.
    expect_error(
        gvtrack.create("count_bad_thresh", NULL, "pwm.count", pssm = pssm, score.thresh = TRUE),
        "score.thresh must be a single number"
    )
    expect_error(
        gvtrack.create("count_bad_thresh", NULL, "pwm.count", pssm = pssm, score.thresh = FALSE),
        "score.thresh must be a single number"
    )
})

test_that("the other pwm functions still do not need score.thresh", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    pssm <- create_test_pssm()

    expect_silent(gvtrack.create("v_pwm", NULL, "pwm", pssm = pssm))
    expect_silent(gvtrack.create("v_pwm_max", NULL, "pwm.max", pssm = pssm))
    expect_silent(gvtrack.create("v_pwm_max_pos", NULL, "pwm.max.pos", pssm = pssm))
})

test_that("pwm.count with a workable score.thresh is unchanged", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    pssm <- create_test_pssm()
    test_interval <- gintervals(1, 200, 300)

    # Independent expectation: count the windows whose pwm score clears the
    # threshold, computed from the per-base pwm vtrack rather than from
    # pwm.count itself.
    gvtrack.create("pwm_ref", NULL, "pwm",
        pssm = pssm, bidirect = FALSE, extend = TRUE,
        prior = 0.01, strand = 1
    )
    per_base <- gextract("pwm_ref", test_interval, iterator = 1)
    expected <- sum(per_base$pwm_ref >= -5, na.rm = TRUE)

    # Anti-vacuity: the threshold must discriminate. A count that is 0 or that
    # equals the window count would make the comparison below pass no matter
    # what the counting path does - which is how both a threshold of 0 (counts
    # nothing) and a threshold of -10 (counts everything, on this PSSM) slipped
    # through the suite before.
    expect_gt(expected, 0)
    expect_lt(expected, nrow(per_base))

    gvtrack.create("count_ref", NULL, "pwm.count",
        pssm = pssm, bidirect = FALSE, extend = TRUE,
        prior = 0.01, strand = 1, score.thresh = -5
    )
    got <- gextract("count_ref", test_interval, iterator = test_interval)

    expect_equal(got$count_ref[1], expected, ignore_attr = TRUE)
})
