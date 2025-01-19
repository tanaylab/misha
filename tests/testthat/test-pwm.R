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
