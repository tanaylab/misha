test_that("misha PWM matches prego results - basic case without spatial, extend=FALSE", {
    skip_if_not_installed("prego")

    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    # Create test intervals on real genome
    test_intervals <- gintervals(1, 10000, 10100)
    test_intervals <- rbind(test_intervals, gintervals(1, 20000, 20100))
    test_intervals <- rbind(test_intervals, gintervals(2, 15000, 15100))

    # Extract sequences from genome
    seqs <- gseq.extract(test_intervals)

    # Create test PSSM - simple 4bp motif
    test_pssm <- data.frame(
        pos = 1:4,
        A = c(0.7, 0.1, 0.1, 0.1),
        C = c(0.1, 0.7, 0.1, 0.1),
        G = c(0.1, 0.1, 0.7, 0.1),
        T = c(0.1, 0.1, 0.1, 0.7)
    )
    pssm_mat <- as.matrix(test_pssm[, c("A", "C", "G", "T")])

    # Compute with prego
    prego_scores <- prego::compute_pwm(
        sequences = seqs,
        pssm = test_pssm,
        spat = NULL,
        bidirect = TRUE,
        prior = 0.01,
        func = "logSumExp"
    )

    # Compute with misha (extend=FALSE to match exact interval)
    gvtrack.create(
        "pwm_test", NULL, "pwm",
        list(
            pssm = pssm_mat,
            bidirect = TRUE,
            extend = FALSE,
            prior = 0.01
        )
    )

    result <- gextract("pwm_test", test_intervals, iterator = test_intervals)

    # Compare - should be very close
    expect_equal(result$pwm_test, prego_scores, tolerance = 1e-6)
    expect_regression(result, "pwm_basic_no_extend")
})

test_that("misha PWM matches prego results - basic case without spatial, extend=TRUE", {
    skip_if_not_installed("prego")

    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    # Create test intervals on real genome
    test_intervals <- gintervals(1, 10000, 10100)
    test_intervals <- rbind(test_intervals, gintervals(1, 20000, 20100))
    test_intervals <- rbind(test_intervals, gintervals(2, 15000, 15100))

    # Create test PSSM - simple 4bp motif
    test_pssm <- data.frame(
        pos = 1:4,
        A = c(0.7, 0.1, 0.1, 0.1),
        C = c(0.1, 0.7, 0.1, 0.1),
        G = c(0.1, 0.1, 0.7, 0.1),
        T = c(0.1, 0.1, 0.1, 0.7)
    )
    pssm_mat <- as.matrix(test_pssm[, c("A", "C", "G", "T")])
    motif_len <- nrow(pssm_mat)

    # Compute with misha using extend=TRUE (no iterator shifts needed)
    gvtrack.create(
        "pwm_test_extend", NULL, "pwm",
        list(
            pssm = pssm_mat,
            bidirect = TRUE,
            extend = TRUE,
            prior = 0.01
        )
    )

    result <- gextract("pwm_test_extend", test_intervals, iterator = test_intervals)

    # For extend=TRUE, misha automatically extends the END by (motif_len - 1)
    # Create extended intervals for prego to match this behavior
    extended_intervals <- test_intervals
    extended_intervals$end <- extended_intervals$end + (motif_len - 1)

    # Extract sequences from extended intervals
    extended_seqs <- gseq.extract(extended_intervals)

    # Compute with prego on extended sequences
    prego_scores <- prego::compute_pwm(
        sequences = extended_seqs,
        pssm = test_pssm,
        spat = NULL,
        bidirect = TRUE,
        prior = 0.01,
        func = "logSumExp"
    )

    # Compare - should be very close
    expect_equal(result$pwm_test_extend, prego_scores, tolerance = 1e-6)
    expect_regression(result, "pwm_basic_with_extend")
})

test_that("Misha PWM matches prego with iterator shifts", {
    skip_if_not_installed("prego")

    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    test_pssm <- data.frame(
        pos = 1:4,
        A = c(0.7, 0.1, 0.1, 0.1),
        C = c(0.1, 0.7, 0.1, 0.1),
        G = c(0.1, 0.1, 0.7, 0.1),
        T = c(0.1, 0.1, 0.1, 0.7)
    )

    pssm_mat <- as.matrix(test_pssm[, c("A", "C", "G", "T")])

    # Load prego
    gvtrack.create(
        "pwm_shift", NULL, "pwm",
        list(
            pssm = pssm_mat,
            bidirect = TRUE,
            extend = FALSE,
            prior = 0.01
        )
    )
    gvtrack.iterator("pwm_shift", sshift = -100, eshift = 100)

    result <- gextract("pwm_shift", gintervals(1, 1000, 2000), iterator = 200)
    test_intervals <- result %>%
        mutate(start = start - 100, end = end + 100)
    result_prego <- prego::compute_pwm(
        sequences = gseq.extract(test_intervals),
        pssm = test_pssm,
        spat = NULL,
        bidirect = TRUE,
        prior = 0.01
    )
    expect_equal(result$pwm_shift, result_prego, tolerance = 1e-5)
    expect_regression(result, "pwm_prego_regression_test_1")
})

test_that("misha PWM matches prego results - with spatial weighting, extend=FALSE", {
    skip_if_not_installed("prego")

    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    # Create test intervals (280bp to match spatial model)
    test_intervals <- gintervals(1, 10000, 10280)
    test_intervals <- rbind(test_intervals, gintervals(1, 20000, 20280))
    test_intervals <- rbind(test_intervals, gintervals(2, 15000, 15280))

    # Extract sequences from genome
    seqs <- gseq.extract(test_intervals)

    # Create test PSSM
    test_pssm <- data.frame(
        pos = 1:4,
        A = c(0.7, 0.1, 0.1, 0.1),
        C = c(0.1, 0.7, 0.1, 0.1),
        G = c(0.1, 0.1, 0.7, 0.1),
        T = c(0.1, 0.1, 0.1, 0.7)
    )
    pssm_mat <- as.matrix(test_pssm[, c("A", "C", "G", "T")])

    # Create spatial model
    spat_df <- data.frame(
        bin = seq(0, 240, by = 40),
        spat_factor = c(0.5, 1.0, 2.0, 2.5, 2.0, 1.0, 0.5)
    )

    # Compute with prego
    prego_scores <- prego::compute_pwm(
        sequences = seqs,
        pssm = test_pssm,
        spat = spat_df,
        bidirect = TRUE,
        prior = 0.01,
        func = "logSumExp"
    )

    # Compute with misha
    gvtrack.create(
        "pwm_spatial_test", NULL, "pwm",
        list(
            pssm = pssm_mat,
            bidirect = TRUE,
            extend = FALSE,
            prior = 0.01,
            spat_factor = spat_df$spat_factor,
            spat_bin = 40L
        )
    )

    result <- gextract("pwm_spatial_test", test_intervals, iterator = test_intervals)

    # Compare - should be very close
    expect_equal(result$pwm_spatial_test, prego_scores, tolerance = 1e-6)
    expect_regression(result, "pwm_spatial_no_extend")
})

test_that("misha PWM matches prego results - with spatial weighting, extend=TRUE", {
    skip_if_not_installed("prego")

    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    # Create test intervals (280bp for spatial model)
    test_intervals <- gintervals(1, 10000, 10280)
    test_intervals <- rbind(test_intervals, gintervals(1, 20000, 20280))
    test_intervals <- rbind(test_intervals, gintervals(2, 15000, 15280))

    # Create test PSSM
    test_pssm <- data.frame(
        pos = 1:4,
        A = c(0.7, 0.1, 0.1, 0.1),
        C = c(0.1, 0.7, 0.1, 0.1),
        G = c(0.1, 0.1, 0.7, 0.1),
        T = c(0.1, 0.1, 0.1, 0.7)
    )
    pssm_mat <- as.matrix(test_pssm[, c("A", "C", "G", "T")])
    motif_len <- nrow(pssm_mat)

    # Create spatial model
    # Extended range: 280 + (motif_len-1) = 283bp
    # With 40bp bins: need 8 bins to cover this
    spat_df <- data.frame(
        bin = seq(0, 280, by = 40),
        spat_factor = c(0.5, 1.0, 2.0, 2.5, 2.0, 1.0, 0.5, 0.5)
    )

    # Compute with misha using extend=TRUE (no iterator shifts)
    gvtrack.create(
        "pwm_spatial_extend", NULL, "pwm",
        list(
            pssm = pssm_mat,
            bidirect = TRUE,
            extend = TRUE,
            prior = 0.01,
            spat_factor = spat_df$spat_factor,
            spat_bin = 40L
        )
    )

    result <- gextract("pwm_spatial_extend", test_intervals, iterator = test_intervals)

    # For extend=TRUE, misha automatically extends the END by (motif_len - 1)
    extended_intervals <- test_intervals
    extended_intervals$end <- extended_intervals$end + (motif_len - 1)

    # Extract sequences from extended intervals
    extended_seqs <- gseq.extract(extended_intervals)

    # Compute with prego on extended sequences
    prego_scores <- prego::compute_pwm(
        sequences = extended_seqs,
        pssm = test_pssm,
        spat = spat_df,
        bidirect = TRUE,
        prior = 0.01,
        func = "logSumExp"
    )

    # Compare - should be very close
    expect_equal(result$pwm_spatial_extend, prego_scores, tolerance = 1e-6)
    expect_regression(result, "pwm_spatial_with_extend")
})

test_that("misha spatial PWM numerical agreement with prego - genome sequences", {
    skip_if_not_installed("prego")

    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    # Create test intervals (280bp to match spatial model)
    test_intervals <- gintervals(c(1, 1, 2), c(10000, 20000, 15000), c(10280, 20280, 15280))

    # Extract sequences from genome
    seqs <- gseq.extract(test_intervals)

    # Create a simple PSSM (2bp motif for speed)
    test_pssm <- data.frame(
        pos = 1:2,
        A = c(0.9, 0.05),
        C = c(0.05, 0.9),
        G = c(0.025, 0.025),
        T = c(0.025, 0.025)
    )
    pssm_mat <- as.matrix(test_pssm[, c("A", "C", "G", "T")])

    # Create spatial model
    spat_df <- data.frame(
        bin = seq(0, 240, by = 40),
        spat_factor = c(0.5, 1.0, 2.0, 2.5, 2.0, 1.0, 0.5)
    )

    # Compute using prego - no spatial
    prego_scores_nospatial <- prego::compute_pwm(
        sequences = seqs,
        pssm = test_pssm,
        spat = NULL,
        bidirect = FALSE,
        prior = 0.01,
        func = "logSumExp"
    )

    # Compute using prego - with spatial
    prego_scores_spatial <- prego::compute_pwm(
        sequences = seqs,
        pssm = test_pssm,
        spat = spat_df,
        bidirect = FALSE,
        prior = 0.01,
        func = "logSumExp"
    )

    # Compute with misha - no spatial
    gvtrack.create(
        "pwm_nospatial", NULL, "pwm",
        list(pssm = pssm_mat, bidirect = FALSE, extend = FALSE, prior = 0.01)
    )
    result_nospatial <- gextract("pwm_nospatial", test_intervals, iterator = test_intervals)

    # Compute with misha - with spatial
    gvtrack.create(
        "pwm_spatial", NULL, "pwm",
        list(
            pssm = pssm_mat,
            bidirect = FALSE,
            extend = FALSE,
            prior = 0.01,
            spat_factor = spat_df$spat_factor,
            spat_bin = 40L
        )
    )
    result_spatial <- gextract("pwm_spatial", test_intervals, iterator = test_intervals)

    # Compare misha with prego
    expect_equal(result_nospatial$pwm_nospatial, prego_scores_nospatial, tolerance = 1e-6)
    expect_equal(result_spatial$pwm_spatial, prego_scores_spatial, tolerance = 1e-6)

    # Validate that spatial factors change the scores
    expect_true(any(abs(prego_scores_spatial - prego_scores_nospatial) > 0.01))
    expect_true(any(abs(result_spatial$pwm_spatial - result_nospatial$pwm_nospatial) > 0.01))

    # Regression tracking
    expect_regression(result_nospatial, "pwm_genome_nospatial")
    expect_regression(result_spatial, "pwm_genome_spatial")
})

test_that("misha max PWM mode matches prego max mode, extend=FALSE", {
    skip_if_not_installed("prego")

    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    # Create test intervals
    test_intervals <- gintervals(1, 10000, 10100)
    test_intervals <- rbind(test_intervals, gintervals(1, 20000, 20100))
    test_intervals <- rbind(test_intervals, gintervals(2, 15000, 15100))

    # Extract sequences
    seqs <- gseq.extract(test_intervals)

    # Create test PSSM
    test_pssm <- data.frame(
        pos = 1:4,
        A = c(0.7, 0.1, 0.1, 0.1),
        C = c(0.1, 0.7, 0.1, 0.1),
        G = c(0.1, 0.1, 0.7, 0.1),
        T = c(0.1, 0.1, 0.1, 0.7)
    )
    pssm_mat <- as.matrix(test_pssm[, c("A", "C", "G", "T")])

    # Compute with prego
    prego_scores <- prego::compute_pwm(
        sequences = seqs,
        pssm = test_pssm,
        spat = NULL,
        bidirect = TRUE,
        prior = 0.01,
        func = "max"
    )

    # Compute with misha
    gvtrack.create(
        "pwm_max_test", NULL, "pwm.max",
        list(
            pssm = pssm_mat,
            bidirect = TRUE,
            extend = FALSE,
            prior = 0.01
        )
    )

    result <- gextract("pwm_max_test", test_intervals, iterator = test_intervals)

    # Compare
    expect_equal(result$pwm_max_test, prego_scores, tolerance = 1e-6)
    expect_regression(result, "pwm_max_no_extend")
})

test_that("misha max PWM mode matches prego max mode, extend=TRUE", {
    skip_if_not_installed("prego")

    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    # Create test intervals
    test_intervals <- gintervals(1, 10000, 10100)
    test_intervals <- rbind(test_intervals, gintervals(1, 20000, 20100))
    test_intervals <- rbind(test_intervals, gintervals(2, 15000, 15100))

    # Create test PSSM
    test_pssm <- data.frame(
        pos = 1:4,
        A = c(0.7, 0.1, 0.1, 0.1),
        C = c(0.1, 0.7, 0.1, 0.1),
        G = c(0.1, 0.1, 0.7, 0.1),
        T = c(0.1, 0.1, 0.1, 0.7)
    )
    pssm_mat <- as.matrix(test_pssm[, c("A", "C", "G", "T")])
    motif_len <- nrow(pssm_mat)

    # Compute with misha using extend=TRUE (no iterator shifts)
    gvtrack.create(
        "pwm_max_extend", NULL, "pwm.max",
        list(
            pssm = pssm_mat,
            bidirect = TRUE,
            extend = TRUE,
            prior = 0.01
        )
    )

    result <- gextract("pwm_max_extend", test_intervals, iterator = test_intervals)

    # For extend=TRUE, misha automatically extends the END by (motif_len - 1)
    extended_intervals <- test_intervals
    extended_intervals$end <- extended_intervals$end + (motif_len - 1)

    # Extract sequences from extended intervals
    extended_seqs <- gseq.extract(extended_intervals)

    # Compute with prego on extended sequences
    prego_scores <- prego::compute_pwm(
        sequences = extended_seqs,
        pssm = test_pssm,
        spat = NULL,
        bidirect = TRUE,
        prior = 0.01,
        func = "max"
    )

    # Compare
    expect_equal(result$pwm_max_extend, prego_scores, tolerance = 1e-6)
    expect_regression(result, "pwm_max_with_extend")
})

test_that("spatial binning is 0-indexed from scan start like prego", {
    # This test validates the spatial binning logic
    # In prego: bins are 0-indexed from the start of the scan range
    # Same bin for both strands
    skip_if_not_installed("prego")

    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    # Create test intervals (300bp)
    test_intervals <- gintervals(c(1, 2), c(10000, 15000), c(10300, 15300))

    # Extract sequences
    seqs <- gseq.extract(test_intervals)

    # Create a PSSM
    test_pssm <- data.frame(
        pos = 1:2,
        A = c(0.7, 0.1),
        C = c(0.1, 0.7),
        G = c(0.1, 0.1),
        T = c(0.1, 0.1)
    )
    pssm_mat <- as.matrix(test_pssm[, c("A", "C", "G", "T")])

    # Create spatial factors - higher weight in first bin
    spat_factors <- c(10.0, 1.0, 1.0)
    spat_bin <- 100L
    spat_df <- data.frame(
        bin = c(0, 100, 200),
        spat_factor = spat_factors
    )

    # Compute with prego
    prego_scores <- prego::compute_pwm(
        sequences = seqs,
        pssm = test_pssm,
        spat = spat_df,
        bidirect = FALSE,
        prior = 0.01,
        func = "logSumExp"
    )

    # Create vtrack
    gvtrack.create(
        "pwm_bintest", NULL, "pwm",
        list(
            pssm = pssm_mat,
            bidirect = FALSE,
            extend = FALSE,
            prior = 0.01,
            spat_factor = spat_factors,
            spat_bin = spat_bin
        )
    )

    # Extract on test intervals
    result <- gextract("pwm_bintest", test_intervals, iterator = test_intervals)

    # Compare with prego
    expect_equal(result$pwm_bintest, prego_scores, tolerance = 1e-6)

    # Verify scores are valid
    expect_false(any(is.na(result$pwm_bintest)))
    expect_false(any(is.infinite(result$pwm_bintest)))

    # Regression tracking
    expect_regression(result, "pwm_spatial_binning")
})

test_that("spatial parameters validation matches prego behavior", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    pssm <- create_test_pssm()

    # Both prego and misha should reject non-positive spatial factors
    expect_error(
        gvtrack.create(
            "pwm_bad", NULL, "pwm",
            list(
                pssm = pssm,
                spat_factor = c(-1, 1, 1),
                spat_bin = 10L
            )
        ),
        "positive"
    )

    expect_error(
        gvtrack.create(
            "pwm_bad", NULL, "pwm",
            list(
                pssm = pssm,
                spat_factor = c(0, 1, 1),
                spat_bin = 10L
            )
        ),
        "positive"
    )

    expect_error(
        gvtrack.create(
            "pwm_bad", NULL, "pwm",
            list(
                pssm = pssm,
                spat_factor = c(1, 1, 1),
                spat_bin = 0L
            )
        ),
        "positive"
    )
})

test_that("uniform spatial factors give same result as no spatial", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    pssm <- create_test_pssm()
    test_interval <- gintervals(1, 200, 300)

    # Uniform spatial factors (all 1.0) should give same result as no spatial
    uniform_spat <- rep(1.0, 10)
    spat_bin <- 10L

    gvtrack.create(
        "pwm_nospatial", NULL, "pwm",
        list(pssm = pssm, bidirect = FALSE, extend = TRUE, prior = 0.01)
    )

    gvtrack.create(
        "pwm_uniform", NULL, "pwm",
        list(
            pssm = pssm,
            bidirect = FALSE,
            extend = TRUE,
            prior = 0.01,
            spat_factor = uniform_spat,
            spat_bin = spat_bin
        )
    )

    scores <- gextract(c("pwm_nospatial", "pwm_uniform"), test_interval, iterator = test_interval)

    # Should be very close (accounting for log(1.0) = 0)
    expect_equal(scores$pwm_nospatial[1], scores$pwm_uniform[1], tolerance = 1e-5)
})

test_that("misha PWM exactly matches prego on genome sequences - no spatial, no extend", {
    skip_if_not_installed("prego")

    remove_all_vtracks()
    withr::defer(remove_all_vtracks())


    # Create test intervals
    test_intervals <- gintervals(1, 10000, 10100)
    test_intervals <- rbind(test_intervals, gintervals(1, 20000, 20100))
    test_intervals <- rbind(test_intervals, gintervals(2, 15000, 15100))

    # Extract sequences from genome
    seqs <- gseq.extract(test_intervals)

    # Create test PSSM
    test_pssm <- data.frame(
        pos = 1:4,
        A = c(0.7, 0.1, 0.1, 0.1),
        C = c(0.1, 0.7, 0.1, 0.1),
        G = c(0.1, 0.1, 0.7, 0.1),
        T = c(0.1, 0.1, 0.1, 0.7)
    )
    pssm_mat <- as.matrix(test_pssm[, c("A", "C", "G", "T")])

    # Compute with prego
    prego_scores <- prego::compute_pwm(
        sequences = seqs,
        pssm = test_pssm,
        spat = NULL,
        bidirect = TRUE,
        prior = 0.01,
        func = "logSumExp"
    )

    # Compute with misha (extend=FALSE to match exact interval)
    gvtrack.create(
        "pwm_test", NULL, "pwm",
        list(
            pssm = pssm_mat,
            bidirect = TRUE,
            extend = FALSE,
            prior = 0.01
        )
    )

    misha_scores <- gextract("pwm_test", test_intervals, iterator = test_intervals)$pwm_test

    # Compare - should be very close
    expect_equal(misha_scores, prego_scores, tolerance = 1e-6)
})

test_that("misha PWM exactly matches prego on genome sequences - with spatial, no extend", {
    skip_if_not_installed("prego")

    remove_all_vtracks()
    withr::defer(remove_all_vtracks())


    # Create test intervals
    test_intervals <- gintervals(1, 10000, 10280)
    test_intervals <- rbind(test_intervals, gintervals(1, 20000, 20280))
    test_intervals <- rbind(test_intervals, gintervals(2, 15000, 15280))

    # Extract sequences from genome
    seqs <- gseq.extract(test_intervals)

    # Create test PSSM
    test_pssm <- data.frame(
        pos = 1:4,
        A = c(0.7, 0.1, 0.1, 0.1),
        C = c(0.1, 0.7, 0.1, 0.1),
        G = c(0.1, 0.1, 0.7, 0.1),
        T = c(0.1, 0.1, 0.1, 0.7)
    )
    pssm_mat <- as.matrix(test_pssm[, c("A", "C", "G", "T")])

    # Create spatial model
    spat_df <- data.frame(
        bin = seq(0, 240, by = 40),
        spat_factor = c(0.5, 1.0, 2.0, 2.5, 2.0, 1.0, 0.5)
    )

    # Compute with prego
    prego_scores <- prego::compute_pwm(
        sequences = seqs,
        pssm = test_pssm,
        spat = spat_df,
        bidirect = TRUE,
        prior = 0.01,
        func = "logSumExp"
    )

    # Compute with misha
    gvtrack.create(
        "pwm_spatial", NULL, "pwm",
        list(
            pssm = pssm_mat,
            bidirect = TRUE,
            extend = FALSE,
            prior = 0.01,
            spat_factor = spat_df$spat_factor,
            spat_bin = 40L
        )
    )

    misha_scores <- gextract("pwm_spatial", test_intervals, iterator = test_intervals)$pwm_spatial

    # Compare - should be very close
    expect_equal(misha_scores, prego_scores, tolerance = 1e-6)
})

test_that("misha PWM max mode exactly matches prego on genome sequences", {
    skip_if_not_installed("prego")

    remove_all_vtracks()
    withr::defer(remove_all_vtracks())


    # Create test intervals
    test_intervals <- gintervals(1, 10000, 10100)
    test_intervals <- rbind(test_intervals, gintervals(1, 20000, 20100))
    test_intervals <- rbind(test_intervals, gintervals(2, 15000, 15100))

    # Extract sequences
    seqs <- gseq.extract(test_intervals)

    # Create test PSSM
    test_pssm <- data.frame(
        pos = 1:4,
        A = c(0.7, 0.1, 0.1, 0.1),
        C = c(0.1, 0.7, 0.1, 0.1),
        G = c(0.1, 0.1, 0.7, 0.1),
        T = c(0.1, 0.1, 0.1, 0.7)
    )
    pssm_mat <- as.matrix(test_pssm[, c("A", "C", "G", "T")])

    # Compute with prego
    prego_scores <- prego::compute_pwm(
        sequences = seqs,
        pssm = test_pssm,
        spat = NULL,
        bidirect = TRUE,
        prior = 0.01,
        func = "max"
    )

    # Compute with misha
    gvtrack.create(
        "pwm_max", NULL, "pwm.max",
        list(
            pssm = pssm_mat,
            bidirect = TRUE,
            extend = FALSE,
            prior = 0.01
        )
    )

    misha_scores <- gextract("pwm_max", test_intervals, iterator = test_intervals)$pwm_max

    # Compare
    expect_equal(misha_scores, prego_scores, tolerance = 1e-6)
})

test_that("misha PWM max with spatial exactly matches prego on genome sequences", {
    skip_if_not_installed("prego")

    remove_all_vtracks()
    withr::defer(remove_all_vtracks())


    # Create test intervals
    test_intervals <- gintervals(1, 10000, 10280)
    test_intervals <- rbind(test_intervals, gintervals(1, 20000, 20280))
    test_intervals <- rbind(test_intervals, gintervals(2, 15000, 15280))

    # Extract sequences
    seqs <- gseq.extract(test_intervals)

    # Create test PSSM
    test_pssm <- data.frame(
        pos = 1:4,
        A = c(0.7, 0.1, 0.1, 0.1),
        C = c(0.1, 0.7, 0.1, 0.1),
        G = c(0.1, 0.1, 0.7, 0.1),
        T = c(0.1, 0.1, 0.1, 0.7)
    )
    pssm_mat <- as.matrix(test_pssm[, c("A", "C", "G", "T")])

    # Create spatial model
    spat_df <- data.frame(
        bin = seq(0, 240, by = 40),
        spat_factor = c(0.5, 1.0, 2.0, 2.5, 2.0, 1.0, 0.5)
    )

    # Compute with prego using max
    prego_scores <- prego::compute_pwm(
        sequences = seqs,
        pssm = test_pssm,
        spat = spat_df,
        bidirect = TRUE,
        prior = 0.01,
        func = "max"
    )

    # Compute with misha
    gvtrack.create(
        "pwm_max_spatial", NULL, "pwm.max",
        list(
            pssm = pssm_mat,
            bidirect = TRUE,
            extend = FALSE,
            prior = 0.01,
            spat_factor = spat_df$spat_factor,
            spat_bin = 40L
        )
    )

    misha_scores <- gextract("pwm_max_spatial", test_intervals, iterator = test_intervals)$pwm_max_spatial

    # Compare
    expect_equal(misha_scores, prego_scores, tolerance = 1e-6)
})

test_that("misha PWM with extend parameter works correctly", {
    skip_if_not_installed("prego")

    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    # Create short interval (shorter than motif)
    test_interval <- gintervals(1, 10000, 10002) # Only 2bp, motif is 4bp

    # Extract sequence
    seq <- gseq.extract(test_interval)

    # Create test PSSM
    test_pssm <- data.frame(
        pos = 1:4,
        A = c(0.7, 0.1, 0.1, 0.1),
        C = c(0.1, 0.7, 0.1, 0.1),
        G = c(0.1, 0.1, 0.7, 0.1),
        T = c(0.1, 0.1, 0.1, 0.7)
    )
    pssm_mat <- as.matrix(test_pssm[, c("A", "C", "G", "T")])

    # With extend=FALSE, should return NA (interval too short)
    gvtrack.create(
        "pwm_no_extend", NULL, "pwm",
        list(
            pssm = pssm_mat,
            bidirect = TRUE,
            extend = FALSE,
            prior = 0.01
        )
    )

    score_no_extend <- gextract("pwm_no_extend", test_interval, iterator = test_interval)$pwm_no_extend
    expect_true(is.na(score_no_extend))

    # With extend=TRUE, should extend and return valid score
    gvtrack.create(
        "pwm_extend", NULL, "pwm",
        list(
            pssm = pssm_mat,
            bidirect = TRUE,
            extend = TRUE,
            prior = 0.01
        )
    )

    # Use iterator with extension
    gvtrack.iterator("pwm_extend", sshift = -10, eshift = 10)
    score_extend <- gextract("pwm_extend", test_interval, iterator = test_interval)$pwm_extend
    expect_false(is.na(score_extend))
})

test_that("misha PWM with spat_min/spat_max exactly matches prego gextract_pwm_old", {
    skip_if_not_installed("prego")

    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    # Create test data matching realistic usage (like Snai motif)
    # Use a longer motif (15bp) and spatial parameters with spat_min/spat_max
    set.seed(42)
    test_pssm <- data.frame(
        pos = 1:15,
        A = runif(15), C = runif(15), G = runif(15), T = runif(15)
    )
    test_pssm[, c("A", "C", "G", "T")] <- test_pssm[, c("A", "C", "G", "T")] /
        rowSums(test_pssm[, c("A", "C", "G", "T")])
    pssm_mat <- as.matrix(test_pssm[, c("A", "C", "G", "T")])

    # Spatial parameters: 11 bins of size 40, covering a 440bp window
    # This matches the pattern from the original bug report
    spat_factors <- c(0.02, 0.03, 0.10, 0.10, 0.11, 0.11, 0.11, 0.10, 0.10, 0.03, 0.02)
    spat_bin <- 40L
    spat_min <- 30L # 1-based R indexing
    spat_max <- 470L # 1-based R indexing (441bp window: positions 30-470)

    # Create test intervals (500bp long, so spat_min/max define scanning window)
    test_intervals <- gintervals(1, c(10000, 20000, 30000), c(10500, 20500, 30500))

    # Create misha vtrack with spatial parameters and range
    gvtrack.create(
        "pwm_spatial_range", NULL, "pwm",
        pssm = pssm_mat,
        bidirect = TRUE,
        spat_factor = spat_factors,
        spat_bin = spat_bin,
        spat_min = spat_min,
        spat_max = spat_max,
        prior = 0.01,
        extend = FALSE
    )

    # Extract with misha
    misha_scores <- gextract("pwm_spatial_range", test_intervals, iterator = test_intervals)

    # Extract with prego reference implementation
    spat_df <- data.frame(
        bin = seq(0, by = spat_bin, length.out = length(spat_factors)),
        spat_factor = spat_factors
    )

    test_pssm$motif <- "test"

    prego_scores <- prego::gextract_pwm_old(
        test_intervals,
        dataset = test_pssm,
        spat = spat_df,
        spat_min = spat_min,
        spat_max = spat_max,
        bidirect = TRUE,
        prior = 0.01
    )

    # Results should match within numerical precision
    # This test verifies the fix for:
    # 1. Correct 1-based to 0-based coordinate conversion
    # 2. Proper adjustment of spat_max for motif length
    # 3. Correct spatial bin alignment
    expect_equal(
        misha_scores$pwm_spatial_range,
        prego_scores$test,
        tolerance = 0.01,
        info = "Misha PWM with spat_min/spat_max should match prego gextract_pwm_old"
    )

    # Also verify without spatial parameters for baseline correctness
    gvtrack.create(
        "pwm_range_nospatial", NULL, "pwm",
        pssm = pssm_mat,
        bidirect = TRUE,
        spat_min = spat_min,
        spat_max = spat_max,
        prior = 0.01,
        extend = FALSE
    )

    misha_nospatial <- gextract("pwm_range_nospatial", test_intervals, iterator = test_intervals)

    prego_nospatial <- prego::gextract_pwm_old(
        test_intervals,
        dataset = test_pssm,
        spat = NULL,
        spat_min = spat_min,
        spat_max = spat_max,
        bidirect = TRUE,
        prior = 0.01
    )

    # Without spatial, should match very precisely
    expect_equal(
        misha_nospatial$pwm_range_nospatial,
        prego_nospatial$test,
        tolerance = 1e-6,
        info = "Misha PWM with spat_min/spat_max (no spatial) should match prego"
    )
})

test_that("misha PWM with sliding optimization matches prego - no spatial", {
    skip_if_not_installed("prego")

    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    # Create test PSSM
    test_pssm <- data.frame(
        pos = 1:4,
        A = c(0.7, 0.1, 0.1, 0.1),
        C = c(0.1, 0.7, 0.1, 0.1),
        G = c(0.1, 0.1, 0.7, 0.1),
        T = c(0.1, 0.1, 0.1, 0.7)
    )
    pssm_mat <- as.matrix(test_pssm[, c("A", "C", "G", "T")])
    motif_len <- nrow(pssm_mat)

    # Use a larger region to trigger sliding optimization
    # Query 5kb region with iterator=1 (stride-1 overlapping windows)
    query_region <- gintervals(1, 10000, 15000)

    # Create misha vtrack
    gvtrack.create(
        "pwm_sliding", NULL, "pwm",
        list(
            pssm = pssm_mat,
            bidirect = TRUE,
            extend = TRUE,
            prior = 0.01
        )
    )
    gvtrack.iterator("pwm_sliding", sshift = 0, eshift = 0)

    # Extract with sliding (iterator=1 means stride-1)
    misha_result <- gextract("pwm_sliding", query_region, iterator = 1)

    # For prego comparison, compute scores for each position
    # Each position needs an extended sequence
    prego_scores <- numeric(nrow(misha_result))
    for (i in seq_len(nrow(misha_result))) {
        interval <- gintervals(
            misha_result$chrom[i],
            misha_result$start[i],
            misha_result$end[i] + (motif_len - 1) # extend=TRUE adds motif_len-1
        )
        seq <- gseq.extract(interval)
        prego_scores[i] <- prego::compute_pwm(
            sequences = seq,
            pssm = test_pssm,
            spat = NULL,
            bidirect = TRUE,
            prior = 0.01,
            func = "logSumExp"
        )
    }

    # Compare sliding results with prego
    expect_equal(misha_result$pwm_sliding, prego_scores, tolerance = 1e-6)
})

test_that("misha PWM with sliding optimization matches prego - with spatial", {
    skip_if_not_installed("prego")

    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    # Create test PSSM
    test_pssm <- data.frame(
        pos = 1:4,
        A = c(0.7, 0.1, 0.1, 0.1),
        C = c(0.1, 0.7, 0.1, 0.1),
        G = c(0.1, 0.1, 0.7, 0.1),
        T = c(0.1, 0.1, 0.1, 0.7)
    )
    pssm_mat <- as.matrix(test_pssm[, c("A", "C", "G", "T")])
    motif_len <- nrow(pssm_mat)

    # Create spatial model: 280bp coverage
    spat_df <- data.frame(
        bin = seq(0, 240, by = 40),
        spat_factor = c(0.5, 1.0, 2.0, 2.5, 2.0, 1.0, 0.5)
    )

    # Query region: 1kb with sliding window
    # Each position will be scored within a 280bp context
    query_region <- gintervals(1, 10000, 11000)

    # Create misha vtrack with spatial
    gvtrack.create(
        "pwm_sliding_spatial", NULL, "pwm",
        list(
            pssm = pssm_mat,
            bidirect = TRUE,
            extend = TRUE,
            prior = 0.01,
            spat_factor = spat_df$spat_factor,
            spat_bin = 40L
        )
    )
    gvtrack.iterator("pwm_sliding_spatial", sshift = 0, eshift = 279) # 280bp window

    # Extract with sliding (iterator=1)
    misha_result <- gextract("pwm_sliding_spatial", query_region, iterator = 1)

    # For prego comparison - compute scores for each position
    prego_scores <- numeric(nrow(misha_result))
    for (i in seq_len(nrow(misha_result))) {
        # Get the extended interval (with motif extension)
        interval <- gintervals(
            misha_result$chrom[i],
            misha_result$start[i],
            misha_result$end[i] + 279 + (motif_len - 1) # eshift + extend
        )
        seq <- gseq.extract(interval)
        prego_scores[i] <- prego::compute_pwm(
            sequences = seq,
            pssm = test_pssm,
            spat = spat_df,
            bidirect = TRUE,
            prior = 0.01,
            func = "logSumExp"
        )
    }

    # Compare sliding results with prego
    expect_equal(misha_result$pwm_sliding_spatial, prego_scores, tolerance = 1e-5)
})

test_that("misha PWM.max with sliding optimization matches prego", {
    skip_if_not_installed("prego")

    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    # Create test PSSM
    test_pssm <- data.frame(
        pos = 1:4,
        A = c(0.7, 0.1, 0.1, 0.1),
        C = c(0.1, 0.7, 0.1, 0.1),
        G = c(0.1, 0.1, 0.7, 0.1),
        T = c(0.1, 0.1, 0.1, 0.7)
    )
    pssm_mat <- as.matrix(test_pssm[, c("A", "C", "G", "T")])
    motif_len <- nrow(pssm_mat)

    # Use smaller region for max mode (slower to compute in prego)
    query_region <- gintervals(1, 10000, 11000)

    # Create misha vtrack
    gvtrack.create(
        "pwm_max_sliding", NULL, "pwm.max",
        list(
            pssm = pssm_mat,
            bidirect = TRUE,
            extend = TRUE,
            prior = 0.01
        )
    )
    gvtrack.iterator("pwm_max_sliding", sshift = 0, eshift = 0)

    # Extract with sliding
    misha_result <- gextract("pwm_max_sliding", query_region, iterator = 1)

    # For prego comparison - compute max scores for each position
    prego_scores <- numeric(nrow(misha_result))
    for (i in seq_len(nrow(misha_result))) {
        interval <- gintervals(
            misha_result$chrom[i],
            misha_result$start[i],
            misha_result$end[i] + (motif_len - 1)
        )
        seq <- gseq.extract(interval)
        prego_scores[i] <- prego::compute_pwm(
            sequences = seq,
            pssm = test_pssm,
            spat = NULL,
            bidirect = TRUE,
            prior = 0.01,
            func = "max"
        )
    }

    # Compare
    expect_equal(misha_result$pwm_max_sliding, prego_scores, tolerance = 1e-6)
})

test_that("misha PWM with sliding optimization matches prego - minus strand with spatial", {
    skip_if_not_installed("prego")

    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    # Create test PSSM
    test_pssm <- data.frame(
        pos = 1:4,
        A = c(0.7, 0.1, 0.1, 0.1),
        C = c(0.1, 0.7, 0.1, 0.1),
        G = c(0.1, 0.1, 0.7, 0.1),
        T = c(0.1, 0.1, 0.1, 0.7)
    )
    pssm_mat <- as.matrix(test_pssm[, c("A", "C", "G", "T")])
    motif_len <- nrow(pssm_mat)

    # Create spatial model
    spat_df <- data.frame(
        bin = seq(0, 240, by = 40),
        spat_factor = c(0.5, 1.0, 2.0, 2.5, 2.0, 1.0, 0.5)
    )

    # Query region with sliding window
    query_region <- gintervals(1, 10000, 11000)

    # Create misha vtrack for minus strand with spatial
    gvtrack.create(
        "pwm_minus_sliding_spatial", NULL, "pwm",
        list(
            pssm = pssm_mat,
            bidirect = FALSE,
            strand = -1,
            extend = TRUE,
            prior = 0.01,
            spat_factor = spat_df$spat_factor,
            spat_bin = 40L
        )
    )
    gvtrack.iterator("pwm_minus_sliding_spatial", sshift = 0, eshift = 279)

    # Extract with sliding
    misha_result <- gextract("pwm_minus_sliding_spatial", query_region, iterator = 1)

    # For prego comparison - need reverse complement
    prego_scores <- numeric(nrow(misha_result))
    for (i in seq_len(nrow(misha_result))) {
        # Get the extended interval (with motif extension)
        # Follow the same pattern as the working test
        interval <- gintervals(
            misha_result$chrom[i],
            max(0L, misha_result$start[i] - (motif_len - 1)), # extend applies to start on minus strand
            misha_result$end[i] + 279 # eshift applies to interval end
        )
        seq_rc <- grevcomp(gseq.extract(interval))
        prego_scores[i] <- prego::compute_pwm(
            sequences = seq_rc,
            pssm = test_pssm,
            spat = spat_df,
            bidirect = FALSE,
            prior = 0.01,
            func = "logSumExp"
        )
    }

    expect_equal(misha_result$pwm_minus_sliding_spatial, prego_scores, tolerance = 1e-5)
})
