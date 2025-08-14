test_that("pwm vtrack works", {
    remove_all_vtracks()

    pssm <- create_test_pssm() # AC

    test_intervals <- gintervals(1, 200, 240)
    seq <- toupper(gseq.extract(test_intervals)) # CCCTAACCCTAACCCTAACCCTAACCCTAACCCCTAACCC
    seq_ext <- toupper(gseq.extract(gintervals(1, 200, 241))) # CCCTAACCCTAACCCTAACCCTAACCCTAACCCCTAACCCT

    gvtrack.create("e", NULL, func = "pwm", pssm = pssm, bidirect = FALSE, extend = FALSE, prior = 0.01)
    gvtrack.create("e_max", NULL, func = "pwm.max", pssm = pssm, bidirect = FALSE, extend = FALSE, prior = 0.01)
    gvtrack.create("max_pos", NULL, func = "pwm.max.pos", pssm = pssm, bidirect = FALSE, extend = TRUE, prior = 0.01)
    gvtrack.create("e_ext", NULL, func = "pwm", pssm = pssm, bidirect = FALSE, extend = TRUE, prior = 0.01)
    gvtrack.create("e_ext_no_prior", NULL, func = "pwm", pssm = pssm, bidirect = FALSE, extend = TRUE, prior = 0)

    scores <- gextract(c("e", "e_ext", "e_max", "max_pos", "e_ext_no_prior"), test_intervals, iterator = test_intervals)

    scores_1bp <- gextract(c("e_ext", "e_ext_no_prior"), test_intervals, iterator = 1)
    scores_man <- manual_pwm_scores_single_strand(seq_ext, pssm, 0.01)
    scores_man_no_prior <- manual_pwm_scores_single_strand(seq_ext, pssm, 0)

    expect_equal(scores_1bp$e_ext, scores_man, tolerance = 1e-6, ignore_attr = TRUE)
    expect_equal(scores$e, log_sum_exp(scores_man[-length(scores_man)]), tolerance = 1e-6, ignore_attr = TRUE)
    expect_equal(scores$e_ext, log_sum_exp(scores_man), tolerance = 1e-6, ignore_attr = TRUE)
    expect_equal(scores$e_max, max(scores_man), tolerance = 1e-5, ignore_attr = TRUE)
    expect_equal(scores$max_pos, which.max(scores_man), ignore_attr = TRUE)

    expect_equal(scores_1bp$e_ext_no_prior, scores_man_no_prior, tolerance = 1e-6, ignore_attr = TRUE)

    # make sure there is a 0 whenever there is 'AC' in the sequence
    AC_positions <- stringr::str_locate_all(seq_ext, "AC")[[1]][, 1]
    expect_true(all(scores_1bp$e_ext_no_prior[AC_positions] == 0))
    expect_true(all(scores_1bp$e_ext_no_prior[-AC_positions] == -Inf)) # the rest should be -Inf
    expect_equal(scores$max_pos[1], AC_positions[1])
})


test_that("pwm vtrack works with 10bp extract", {
    remove_all_vtracks()

    pssm <- create_test_pssm() # AC

    test_intervals <- gintervals(1, 200, 240)
    seq <- toupper(gseq.extract(test_intervals)) # CCCTAACCCTAACCCTAACCCTAACCCTAACCCCTAACCC
    seq_ext <- toupper(gseq.extract(gintervals(1, 200, 241))) # CCCTAACCCTAACCCTAACCCTAACCCTAACCCCTAACCCT

    gvtrack.create("e", NULL, func = "pwm", pssm = pssm, bidirect = FALSE, extend = TRUE, prior = 0.01)

    scores_10bp <- gextract("e", test_intervals, iterator = 10)
    scores_100bp <- gextract("e", test_intervals, iterator = 100)

    expect_equal(log_sum_exp(scores_10bp$e), scores_100bp$e, tolerance = 1e-6, ignore_attr = TRUE)
})

test_that("pwm vtrack bidirect returns the sum of the two strands", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    pssm <- create_test_pssm() # AC
    motif_length <- nrow(pssm)

    test_interval <- gintervals(1, 200, 240)

    # Get extended sequence to match extend=TRUE behavior
    test_interval_ext <- test_interval
    test_interval_ext$end <- test_interval_ext$end + motif_length - 1

    seq <- toupper(gseq.extract(test_interval_ext))
    seq_rc <- grevcomp(seq)

    # Create virtual tracks with different directionality settings
    gvtrack.create(
        "pwm_bidi", NULL, "pwm",
        list(pssm = pssm, bidirect = TRUE, extend = TRUE, prior = 0.01)
    )
    gvtrack.create(
        "pwm_fwd", NULL, "pwm",
        list(pssm = pssm, bidirect = FALSE, extend = TRUE, prior = 0.01)
    )

    # Extract scores
    scores <- gextract(c("pwm_bidi", "pwm_fwd"), test_interval, iterator = test_interval)

    # Calculate manual scores using extended sequence
    fwd_scores <- manual_pwm_scores_single_strand(seq, pssm, prior = 0.01)
    rev_scores <- manual_pwm_scores_single_strand(seq_rc, pssm, prior = 0.01)

    # Calculate total log-likelihoods
    manual_fwd_total <- log_sum_exp(fwd_scores)
    manual_bidi_total <- log_sum_exp(c(fwd_scores, rev_scores))

    # Test with appropriate tolerance
    expect_equal(scores$pwm_fwd[1], manual_fwd_total, tolerance = 1e-6)
    expect_equal(scores$pwm_bidi[1], manual_bidi_total, tolerance = 1e-6)
    expect_true(scores$pwm_bidi[1] >= scores$pwm_fwd[1])
})

test_that("pwm honors gvtrack.iterator shifts", {
    remove_all_vtracks()

    pssm <- create_test_pssm() # AC motif (length 2)

    # Base interval used as both intervals and iterator (single-bin)
    base <- gintervals(1, 2000, 2040)

    # Two identical PWM vtracks differing only by per-vtrack iterator shifts
    gvtrack.create("pwm_small_win", NULL, func = "pwm", pssm = pssm, bidirect = FALSE, extend = TRUE, prior = 0.01)
    gvtrack.create("pwm_large_win", NULL, func = "pwm", pssm = pssm, bidirect = FALSE, extend = TRUE, prior = 0.01)

    # Apply per-vtrack iterator shifts
    gvtrack.iterator("pwm_small_win", sshift = -10, eshift = 10)
    gvtrack.iterator("pwm_large_win", sshift = -1000, eshift = 1000)

    scores <- gextract(c("pwm_small_win", "pwm_large_win"), base, iterator = base)

    # Manually compute expected totals using the exact windows the scorer will read
    motif_len <- nrow(pssm)

    # NOTE: Shifts modify start and end independently (not center).
    # For PWM with extend=TRUE and strand=1, only the end is extended by (motif_len - 1) beyond the shifted end.
    ext_small <- base
    ext_small$start <- pmax(0, ext_small$start - 10)
    ext_small$end <- ext_small$end + 10 + (motif_len - 1)
    seq_small <- toupper(gseq.extract(ext_small))
    manual_small <- manual_pwm_scores_single_strand(seq_small, pssm, prior = 0.01)
    expect_equal(scores$pwm_small_win[1], log_sum_exp(manual_small), tolerance = 1e-6)

    ext_large <- base
    ext_large$start <- pmax(0, ext_large$start - 1000)
    ext_large$end <- ext_large$end + 1000 + (motif_len - 1)
    seq_large <- toupper(gseq.extract(ext_large))
    manual_large <- manual_pwm_scores_single_strand(seq_large, pssm, prior = 0.01)
    # Use a slightly looser tolerance here due to accumulated floating-point error on large windows
    expect_equal(scores$pwm_large_win[1], log_sum_exp(manual_large), tolerance = 1e-4)

    # The two totals should generally differ due to different windows
    expect_true(scores$pwm_small_win[1] != scores$pwm_large_win[1])
})

test_that("pwm scoring works correctly for forward and reverse strands", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    # Create PSSM with clear strand preference
    pssm <- matrix(c(
        0.8, 0.1, 0.05, 0.05, # Strong A
        0.1, 0.8, 0.05, 0.05, # Strong C
        0.8, 0.1, 0.05, 0.05 # Strong A
    ), ncol = 4, byrow = TRUE)
    colnames(pssm) <- c("A", "C", "G", "T")
    motif_length <- nrow(pssm)

    # Get intervals and sequence
    test_intervals <- gintervals(1, 200, 240)
    test_intervals_ext <- test_intervals
    test_intervals_ext$end <- test_intervals_ext$end + motif_length - 1

    # Get forward sequence with extension
    seq_fwd <- toupper(gseq.extract(test_intervals_ext))

    # Get reverse sequence with extension and reversed strand
    test_intervals_rev <- test_intervals_ext
    test_intervals_rev$strand <- -1
    seq_rev <- toupper(gseq.extract(test_intervals_rev)) # This will give reverse complement

    test_intervals_rev <- test_intervals
    test_intervals_rev$start <- test_intervals_rev$start - motif_length + 1
    test_intervals_rev$strand <- -1
    seq_rev_ext <- toupper(gseq.extract(test_intervals_rev))

    # Create tracks for forward, reverse, and bidirectional scanning
    gvtrack.create(
        "pwm_fwd", NULL, "pwm",
        list(pssm = pssm, bidirect = FALSE, extend = TRUE, prior = 0.01)
    )
    gvtrack.create(
        "pwm_rev", NULL, "pwm",
        list(pssm = pssm, bidirect = FALSE, extend = TRUE, prior = 0.01, strand = -1)
    )
    gvtrack.create(
        "pwm_bidi", NULL, "pwm",
        list(pssm = pssm, bidirect = TRUE, extend = TRUE, prior = 0.01)
    )

    # Extract scores for both strands
    scores_plus <- gextract(c("pwm_fwd", "pwm_rev", "pwm_bidi"),
        test_intervals,
        iterator = test_intervals
    )

    # Calculate manual scores
    fwd_scores <- manual_pwm_scores_single_strand(seq_fwd, pssm, prior = 0.01)
    rev_scores <- manual_pwm_scores_single_strand(seq_rev, pssm, prior = 0.01)
    rev_scors_ext <- manual_pwm_scores_single_strand(seq_rev_ext, pssm, prior = 0.01)

    # Calculate total log-likelihoods
    manual_fwd_total <- log_sum_exp(fwd_scores)
    manual_rev_total <- log_sum_exp(rev_scores)
    manual_rev_ext_total <- log_sum_exp(rev_scors_ext)
    manual_bidi_total <- log_sum_exp(c(fwd_scores, rev_scores))

    # Test all scores
    expect_equal(scores_plus$pwm_fwd[1], manual_fwd_total, tolerance = 1e-6)
    expect_equal(scores_plus$pwm_rev[1], manual_rev_ext_total, tolerance = 1e-6)
    expect_equal(scores_plus$pwm_bidi[1], manual_bidi_total, tolerance = 1e-6)

    # Test logical properties
    expect_true(scores_plus$pwm_bidi[1] >= scores_plus$pwm_fwd[1])
    expect_true(scores_plus$pwm_bidi[1] >= scores_plus$pwm_rev[1])
})

test_that("pwm honors iterator shifts across multiple magnitudes (extend=TRUE, fwd strand)", {
    remove_all_vtracks()

    pssm <- create_test_pssm() # AC motif (length 2)
    motif_len <- nrow(pssm)

    # Base interval used as both intervals and iterator (single-bin)
    base <- gintervals(1, 2000, 2040)

    # Build several vtracks with increasing iterator shifts
    shift_values <- c(0, 1, 5, 10, 25, 50, 100, 250, 1000)
    vnames <- sprintf("pwm_shift_%d", shift_values)

    for (i in seq_along(shift_values)) {
        s <- shift_values[i]
        gvtrack.create(vnames[i], NULL, func = "pwm", pssm = pssm, bidirect = FALSE, extend = TRUE, prior = 0.01)
        gvtrack.iterator(vnames[i], sshift = -s, eshift = s)
    }

    scores <- gextract(vnames, base, iterator = base)

    # Manually compute expectations for each shift
    for (i in seq_along(shift_values)) {
        s <- shift_values[i]
        ext <- base
        ext$start <- pmax(0, ext$start - s)
        ext$end <- ext$end + s + (motif_len - 1)
        seq_ext <- toupper(gseq.extract(ext))
        manual <- manual_pwm_scores_single_strand(seq_ext, pssm, prior = 0.01)

        expected <- log_sum_exp(manual)
        actual <- scores[[vnames[i]]][1]

        tol <- if (s >= 250) 1e-4 else 1e-6
        expect_equal(actual, expected, tolerance = tol)
    }
})

test_that("pwm honors iterator shifts without extension (extend=FALSE, fwd strand)", {
    remove_all_vtracks()

    pssm <- create_test_pssm() # AC motif (length 2)
    motif_len <- nrow(pssm)

    base <- gintervals(1, 2000, 2040)

    # Use a variety of shifts
    shift_values <- c(0, 5, 10, 50, 200)
    vnames <- sprintf("pwm_noext_shift_%d", shift_values)

    for (i in seq_along(shift_values)) {
        s <- shift_values[i]
        gvtrack.create(vnames[i], NULL, func = "pwm", pssm = pssm, bidirect = FALSE, extend = FALSE, prior = 0.01)
        gvtrack.iterator(vnames[i], sshift = -s, eshift = s)
    }

    scores <- gextract(vnames, base, iterator = base)

    for (i in seq_along(shift_values)) {
        s <- shift_values[i]
        # For extend=FALSE, no motif-length extension is applied
        ext <- base
        ext$start <- pmax(0, ext$start - s)
        ext$end <- ext$end + s
        seq_noext <- toupper(gseq.extract(ext))

        manual <- manual_pwm_scores_single_strand(seq_noext, pssm, prior = 0.01)
        expected <- log_sum_exp(manual)
        actual <- scores[[vnames[i]]][1]

        tol <- if (s >= 50) 1e-5 else 1e-6
        expect_equal(actual, expected, tolerance = tol)
    }
})

test_that("pwm bidirectional honors iterator shifts (extend=TRUE)", {
    remove_all_vtracks()

    pssm <- create_test_pssm() # AC motif (length 2)
    motif_len <- nrow(pssm)

    base <- gintervals(1, 200, 240)

    # A couple of representative shifts
    shift_values <- c(10, 200)
    vnames <- sprintf("pwm_bidi_shift_%d", shift_values)

    for (i in seq_along(shift_values)) {
        s <- shift_values[i]
        gvtrack.create(vnames[i], NULL, func = "pwm", pssm = pssm, bidirect = TRUE, extend = TRUE, prior = 0.01)
        gvtrack.iterator(vnames[i], sshift = -s, eshift = s)
    }

    scores <- gextract(vnames, base, iterator = base)

    for (i in seq_along(shift_values)) {
        s <- shift_values[i]
        ext <- base
        ext$start <- pmax(0, ext$start - s)
        ext$end <- ext$end + s + (motif_len - 1)
        seq_fwd <- toupper(gseq.extract(ext))
        seq_rev <- grevcomp(seq_fwd)

        fwd_scores <- manual_pwm_scores_single_strand(seq_fwd, pssm, prior = 0.01)
        rev_scores <- manual_pwm_scores_single_strand(seq_rev, pssm, prior = 0.01)
        expected <- log_sum_exp(c(fwd_scores, rev_scores))
        actual <- scores[[vnames[i]]][1]

        tol <- if (s >= 200) 1e-4 else 1e-6
        expect_equal(actual, expected, tolerance = tol)
    }
})

test_that("pwm.max honors iterator shifts (extend=TRUE, fwd strand)", {
    remove_all_vtracks()

    pssm <- create_test_pssm()
    motif_len <- nrow(pssm)

    base <- gintervals(1, 2000, 2040)

    shift_values <- c(0, 1, 10, 200)
    vnames <- sprintf("pwm_max_shift_%d", shift_values)

    for (i in seq_along(shift_values)) {
        s <- shift_values[i]
        gvtrack.create(vnames[i], NULL, func = "pwm.max", pssm = pssm, bidirect = FALSE, extend = TRUE, prior = 0.01)
        gvtrack.iterator(vnames[i], sshift = -s, eshift = s)
    }

    scores <- gextract(vnames, base, iterator = base)

    for (i in seq_along(shift_values)) {
        s <- shift_values[i]
        ext <- base
        ext$start <- pmax(0, ext$start - s)
        ext$end <- ext$end + s + (motif_len - 1)
        seq_ext <- toupper(gseq.extract(ext))
        manual <- manual_pwm_scores_single_strand(seq_ext, pssm, prior = 0.01)
        expected <- max(manual)
        actual <- scores[[vnames[i]]][1]
        # Single-precision arithmetic in the C++ scorer can introduce tiny differences vs R double
        tol <- if (s >= 200) 5e-6 else 5e-6
        expect_lt(abs(actual - expected), tol)
    }
})

test_that("pwm.max honors iterator shifts without extension (extend=FALSE, fwd strand)", {
    remove_all_vtracks()

    pssm <- create_test_pssm()

    base <- gintervals(1, 2000, 2040)

    shift_values <- c(0, 10, 200)
    vnames <- sprintf("pwm_max_noext_shift_%d", shift_values)

    for (i in seq_along(shift_values)) {
        s <- shift_values[i]
        gvtrack.create(vnames[i], NULL, func = "pwm.max", pssm = pssm, bidirect = FALSE, extend = FALSE, prior = 0.01)
        gvtrack.iterator(vnames[i], sshift = -s, eshift = s)
    }

    scores <- gextract(vnames, base, iterator = base)

    for (i in seq_along(shift_values)) {
        s <- shift_values[i]
        ext <- base
        ext$start <- pmax(0, ext$start - s)
        ext$end <- ext$end + s
        seq_noext <- toupper(gseq.extract(ext))
        manual <- manual_pwm_scores_single_strand(seq_noext, pssm, prior = 0.01)
        expected <- max(manual)
        actual <- scores[[vnames[i]]][1]
        tol <- if (s >= 50) 5e-6 else 5e-6
        expect_lt(abs(actual - expected), tol)
    }
})

test_that("pwm.max.pos honors iterator shifts (extend=TRUE, fwd strand)", {
    remove_all_vtracks()

    pssm <- create_test_pssm()
    motif_len <- nrow(pssm)

    base <- gintervals(1, 200, 240)

    shift <- 100
    gvtrack.create("pwm_max_pos_shift", NULL, func = "pwm.max.pos", pssm = pssm, bidirect = FALSE, extend = TRUE, prior = 0.01)
    gvtrack.iterator("pwm_max_pos_shift", sshift = -shift, eshift = shift)

    scores <- gextract("pwm_max_pos_shift", base, iterator = base)

    ext <- base
    ext$start <- pmax(0, ext$start - shift)
    ext$end <- ext$end + shift + (motif_len - 1)
    seq_ext <- toupper(gseq.extract(ext))
    manual <- manual_pwm_scores_single_strand(seq_ext, pssm, prior = 0.01)
    expected_pos <- which.max(manual)

    expect_equal(scores$pwm_max_pos_shift[1], expected_pos)
})

test_that("pwm.max bidirectional equals max of forward and reverse (extend=TRUE)", {
    remove_all_vtracks()

    pssm <- create_test_pssm()
    motif_len <- nrow(pssm)

    base <- gintervals(1, 200, 240)

    shift <- 200
    gvtrack.create("pwm_max_bidi", NULL, func = "pwm.max", pssm = pssm, bidirect = TRUE, extend = TRUE, prior = 0.01)
    gvtrack.iterator("pwm_max_bidi", sshift = -shift, eshift = shift)

    scores <- gextract("pwm_max_bidi", base, iterator = base)

    ext <- base
    ext$start <- pmax(0, ext$start - shift)
    ext$end <- ext$end + shift + (motif_len - 1)
    seq_fwd <- toupper(gseq.extract(ext))
    seq_rev <- grevcomp(seq_fwd)
    fwd_scores <- manual_pwm_scores_single_strand(seq_fwd, pssm, prior = 0.01)
    rev_scores <- manual_pwm_scores_single_strand(seq_rev, pssm, prior = 0.01)

    expected <- max(c(fwd_scores, rev_scores))
    # Slightly looser tolerance due to bidirectional scanning and single-precision accumulation
    expect_lt(abs(scores$pwm_max_bidi[1] - expected), 1e-3)
})

test_that("pwm total likelihood is monotonic with larger windows", {
    remove_all_vtracks()

    pssm <- create_test_pssm()
    motif_len <- nrow(pssm)

    base <- gintervals(1, 2000, 2040)

    shifts <- c(0, 10, 50, 200)
    vnames <- sprintf("pwm_total_shift_%d", shifts)

    for (i in seq_along(shifts)) {
        s <- shifts[i]
        gvtrack.create(vnames[i], NULL, func = "pwm", pssm = pssm, bidirect = FALSE, extend = TRUE, prior = 0.01)
        gvtrack.iterator(vnames[i], sshift = -s, eshift = s)
    }

    scores <- gextract(vnames, base, iterator = base)
    vals <- as.numeric(scores[1, vnames])
    # Non-decreasing with larger windows (superset of positions)
    expect_true(all(diff(vals) >= -1e-7))
})

test_that("iterator shift equals explicit iterator expansion for pwm (single-bin, extend=TRUE)", {
    remove_all_vtracks()

    pssm <- create_test_pssm()
    motif_len <- nrow(pssm)

    # Base 60bp bin
    base60 <- gintervals(1, 2000, 2060)
    # Expanded iterator by 10bp on each side => 80bp
    base80 <- gintervals(1, 1990, 2070)

    # Shifted vtrack over base60
    gvtrack.create("pwm_shifted", NULL, func = "pwm", pssm = pssm, bidirect = FALSE, extend = TRUE, prior = 0.01)
    gvtrack.iterator("pwm_shifted", sshift = -10, eshift = 10)
    s_shift <- gextract("pwm_shifted", base60, iterator = base60)

    # Unshifted vtrack over expanded iterator
    gvtrack.create("pwm_unshifted", NULL, func = "pwm", pssm = pssm, bidirect = FALSE, extend = TRUE, prior = 0.01)
    s_unshift <- gextract("pwm_unshifted", base80, iterator = base80)

    expect_lt(abs(s_shift$pwm_shifted[1] - s_unshift$pwm_unshifted[1]), 1e-6)
})

test_that("iterator shift equals explicit iterator expansion for pwm.max (single-bin, extend=TRUE)", {
    remove_all_vtracks()

    pssm <- create_test_pssm()

    base60 <- gintervals(1, 2100, 2160)
    base80 <- gintervals(1, 2090, 2170)

    gvtrack.create("pwmmax_shifted", NULL, func = "pwm.max", pssm = pssm, bidirect = FALSE, extend = TRUE, prior = 0.01)
    gvtrack.iterator("pwmmax_shifted", sshift = -10, eshift = 10)
    a_shift <- gextract("pwmmax_shifted", base60, iterator = base60)

    gvtrack.create("pwmmax_unshifted", NULL, func = "pwm.max", pssm = pssm, bidirect = FALSE, extend = TRUE, prior = 0.01)
    a_unshift <- gextract("pwmmax_unshifted", base80, iterator = base80)

    expect_lt(abs(a_shift$pwmmax_shifted[1] - a_unshift$pwmmax_unshifted[1]), 5e-6)
})

test_that("iterator shift equals explicit iterator expansion for pwm.max.pos (single-bin, extend=TRUE)", {
    remove_all_vtracks()

    pssm <- create_test_pssm()

    base60 <- gintervals(1, 2200, 2260)
    base80 <- gintervals(1, 2190, 2270)

    gvtrack.create("pwmmaxpos_shifted", NULL, func = "pwm.max.pos", pssm = pssm, bidirect = FALSE, extend = TRUE, prior = 0.01)
    gvtrack.iterator("pwmmaxpos_shifted", sshift = -10, eshift = 10)
    p_shift <- gextract("pwmmaxpos_shifted", base60, iterator = base60)

    gvtrack.create("pwmmaxpos_unshifted", NULL, func = "pwm.max.pos", pssm = pssm, bidirect = FALSE, extend = TRUE, prior = 0.01)
    p_unshift <- gextract("pwmmaxpos_unshifted", base80, iterator = base80)

    expect_equal(p_shift$pwmmaxpos_shifted[1], p_unshift$pwmmaxpos_unshifted[1])
})
