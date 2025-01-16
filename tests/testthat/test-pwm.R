test_that("PWM total score (logsumexp) computation is correct", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    pssm <- create_test_pssm()
    motif_length <- nrow(pssm)

    # Create interval large enough for multiple PWM positions
    test_interval <- gintervals(1, 0, 100)

    # Create tracks for total likelihood (uses log sum exp) and max
    gvtrack.create(
        "pwm_total", NULL, "pwm",
        list(pssm = pssm, bidirect = FALSE, extend = TRUE)
    )
    gvtrack.create(
        "pwm_max", NULL, "pwm.max",
        list(pssm = pssm, bidirect = FALSE, extend = TRUE)
    )

    # Extract using iterator >= motif length to get complete PWM scores
    res <- gextract(c("pwm_total", "pwm_max"), test_interval, iterator = motif_length)

    expect_true(all(res$pwm_total >= res$pwm_max))
})

test_that("PWM scanning with extend works correctly", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    # Create simple PSSM that's easy to verify
    simple_pssm <- matrix(c(
        1.0, 0.0, 0.0, 0.0, # Only A
        0.0, 1.0, 0.0, 0.0 # Only C
    ), ncol = 4, byrow = TRUE)
    colnames(simple_pssm) <- c("A", "C", "G", "T")
    motif_length <- nrow(simple_pssm)

    # Test scanning with different interval sizes
    gvtrack.create(
        "pwm_scan", NULL, "pwm.max",
        list(pssm = simple_pssm, bidirect = FALSE, extend = TRUE)
    )

    # Test intervals of different sizes
    res_short <- gextract("pwm_scan", gintervals(1, 0, motif_length - 1),
        iterator = motif_length
    )
    res_exact <- gextract("pwm_scan", gintervals(1, 0, motif_length),
        iterator = motif_length
    )
    res_long <- gextract("pwm_scan", gintervals(1, 0, motif_length + 5),
        iterator = motif_length
    )

    # With extend=TRUE, all should have valid scores
    expect_false(all(is.na(res_short$pwm_scan)))
    expect_false(all(is.na(res_exact$pwm_scan)))
    expect_false(all(is.na(res_long$pwm_scan)))
})

test_that("PWM scanning without extend requires sufficient interval length", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    pssm <- create_test_pssm()
    motif_length <- nrow(pssm)

    gvtrack.create(
        "pwm_scan", NULL, "pwm.max",
        list(pssm = pssm, bidirect = FALSE, extend = FALSE)
    )

    # Test intervals shorter than motif
    res_short <- gextract("pwm_scan", gintervals(1, 0, motif_length - 1),
        iterator = motif_length
    )
    expect_true(all(is.na(res_short$pwm_scan)))

    # Test intervals longer than motif
    res_long <- gextract("pwm_scan", gintervals(1, 0, motif_length + 5),
        iterator = motif_length
    )
    expect_false(all(is.na(res_long$pwm_scan)))
})

test_that("Bidirectional PWM scanning finds best matches on both strands", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    # Create asymmetric PSSM that will give different scores on different strands
    asym_pssm <- matrix(c(
        0.9, 0.1, 0.0, 0.0, # Strong A preference
        0.0, 0.9, 0.1, 0.0 # Strong C preference
    ), ncol = 4, byrow = TRUE)
    colnames(asym_pssm) <- c("A", "C", "G", "T")
    motif_length <- nrow(asym_pssm)

    # Create tracks with and without bidirectional search
    gvtrack.create(
        "pwm_unidirect", NULL, "pwm.max",
        list(pssm = asym_pssm, bidirect = FALSE, extend = TRUE)
    )
    gvtrack.create(
        "pwm_bidirect", NULL, "pwm.max",
        list(pssm = asym_pssm, bidirect = TRUE, extend = TRUE)
    )
    gvtrack.create(
        "pwm_pos", NULL, "pwm.max.pos",
        list(pssm = asym_pssm, bidirect = TRUE, extend = TRUE)
    )

    # Test on sufficiently large interval to find matches
    res <- gextract(c("pwm_unidirect", "pwm_bidirect", "pwm_pos"),
        gintervals(1, 0, 50),
        iterator = motif_length
    )

    # Bidirectional search should sometimes find better scores
    expect_true(any(res$pwm_bidirect > res$pwm_unidirect))

    # Should see both positive and negative positions
    expect_true(any(res$pwm_pos > 0) && any(res$pwm_pos < 0))
})

test_that("PWM scanner correctly handles priors", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    pssm <- create_test_pssm()
    motif_length <- nrow(pssm)

    # Create same PWM with different priors
    gvtrack.create(
        "pwm_noprior", NULL, "pwm.max",
        list(pssm = pssm, bidirect = FALSE, prior = 0)
    )
    gvtrack.create(
        "pwm_prior", NULL, "pwm.max",
        list(pssm = pssm, bidirect = FALSE, prior = 0.1)
    )

    # Test on interval larger than motif
    res <- gextract(c("pwm_noprior", "pwm_prior"),
        gintervals(1, 0, motif_length + 10),
        iterator = motif_length
    )

    # Scores should differ with prior
    expect_false(all(res$pwm_noprior == res$pwm_prior))
})

test_that("Manual energy calculation matches PWM scorer results", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    # Create a simple PSSM that's easy to verify manually
    simple_pssm <- matrix(c(
        0.8, 0.1, 0.05, 0.05, # Strong A preference
        0.1, 0.8, 0.05, 0.05 # Strong C preference
    ), ncol = 4, byrow = TRUE)
    colnames(simple_pssm) <- c("A", "C", "G", "T")
    motif_length <- nrow(simple_pssm)

    # Create interval and get sequence
    test_interval <- gintervals(1, 0, 10)
    test_interval_ext <- test_interval
    test_interval_ext$end <- test_interval_ext$end + motif_length - 1
    seq <- toupper(gseq.extract(test_interval_ext))
    seq_str <- seq[1]

    # Create PWM tracks for different scoring modes
    gvtrack.create(
        "pwm_total", NULL, "pwm",
        list(pssm = simple_pssm, bidirect = FALSE, extend = TRUE, prior = 0)
    )
    gvtrack.create(
        "pwm_max", NULL, "pwm.max",
        list(pssm = simple_pssm, bidirect = FALSE, extend = TRUE, prior = 0)
    )
    gvtrack.create(
        "pwm_pos", NULL, "pwm.max.pos",
        list(pssm = simple_pssm, bidirect = FALSE, extend = TRUE, prior = 0)
    )

    # Extract PWM scores
    res <- gextract(c("pwm_total", "pwm_max", "pwm_pos"), test_interval, iterator = 1)
    res_interval <- gextract(c("pwm_total", "pwm_max", "pwm_pos"), test_interval, iterator = test_interval)

    # Manual calculation for first interval
    manual_scores <- numeric()
    for (i in 1:(nchar(seq_str) - motif_length + 1)) {
        subseq <- substr(seq_str, i, i + motif_length - 1)
        logp <- 0
        for (j in 1:motif_length) {
            base <- substr(subseq, j, j)
            base_idx <- switch(base,
                "A" = 1,
                "C" = 2,
                "G" = 3,
                "T" = 4
            )
            p <- simple_pssm[j, base_idx]
            logp <- logp + log(p)
        }
        manual_scores <- c(manual_scores, logp)
    }

    # Calculate total score using log-sum-exp
    manual_total <- log_sum_log_vec(manual_scores)
    manual_max <- max(manual_scores)

    # Find position of max score (1-based)
    manual_max_pos <- which.max(manual_scores)

    expect_equal(res$pwm_total, manual_scores, tolerance = 1e-6, ignore_attr = TRUE)
    expect_equal(res$pwm_max, manual_scores, tolerance = 1e-6, ignore_attr = TRUE)
    expect_equal(abs(res$pwm_pos), rep(1, nrow(res)), tolerance = 1e-6)

    # Compare aggregated scores
    expect_equal(res_interval$pwm_total, manual_total, tolerance = 1e-6)
    expect_equal(res_interval$pwm_max, manual_max, tolerance = 1e-6)
    expect_equal(abs(res_interval$pwm_pos), manual_max_pos, tolerance = 1e-6, ignore_attr = TRUE)
})


