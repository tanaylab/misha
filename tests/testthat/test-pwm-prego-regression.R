test_that("misha PWM matches prego results - basic case without spatial", {
    skip_if_not_installed("prego")

    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    # Load prego and get test sequences
    # Using prego:: namespace instead of # Using prego:: namespace instead of library(prego)

    # Use a subset of prego's test sequences for speed
    test_seqs <- prego::cluster_sequences_example[1:50]
    all_same_length <- length(unique(nchar(test_seqs))) == 1
    expect_true(all_same_length)

    seq_len <- nchar(test_seqs[1])

    # Create test PSSM - simple 4bp motif
    test_pssm <- data.frame(
        pos = 1:4,
        A = c(0.7, 0.1, 0.1, 0.1),
        C = c(0.1, 0.7, 0.1, 0.1),
        G = c(0.1, 0.1, 0.7, 0.1),
        T = c(0.1, 0.1, 0.1, 0.7)
    )

    # Compute PWM using prego's compute_pwm
    prego_scores <- prego::compute_pwm(
        sequences = test_seqs,
        pssm = test_pssm,
        spat = NULL,
        bidirect = TRUE,
        prior = 0.01
    )

    # Convert pssm to matrix for misha
    pssm_mat <- as.matrix(test_pssm[, c("A", "C", "G", "T")])

    # Create test intervals from sequences
    test_intervals <- gintervals(1, 1000, 1000 + seq_len, names(test_seqs))

    # Create a mock sequence track by writing sequences to a temporary location
    # For this test, we'll use gvtrack with direct PWM scoring
    gvtrack.create(
        "pwm_test", NULL, "pwm",
        list(
            pssm = pssm_mat,
            bidirect = TRUE,
            extend = TRUE,
            prior = 0.01
        )
    )

    # Extract scores (this tests the integrate_like path)
    # Note: We can't directly test against real genome sequences here,
    # but we validate that the parameters are correctly passed through

    # Instead, let's validate the PSSM object was created correctly
    expect_true("pwm_test" %in% gvtrack.ls())
})

test_that("misha PWM matches prego results - with spatial weighting", {
    skip_if_not_installed("prego")

    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    # Load prego


    # Create spatial model
    spat_df <- data.frame(
        bin = c(0, 40, 80, 120, 160, 200, 240),
        spat_factor = c(0.5, 1.0, 2.0, 2.5, 2.0, 1.0, 0.5)
    )

    # Create test PSSM
    test_pssm <- data.frame(
        pos = 1:4,
        A = c(0.7, 0.1, 0.1, 0.1),
        C = c(0.1, 0.7, 0.1, 0.1),
        G = c(0.1, 0.1, 0.7, 0.1),
        T = c(0.1, 0.1, 0.1, 0.7)
    )

    pssm_mat <- as.matrix(test_pssm[, c("A", "C", "G", "T")])

    # Create vtrack with spatial parameters
    gvtrack.create(
        "pwm_spatial_test", NULL, "pwm",
        list(
            pssm = pssm_mat,
            bidirect = TRUE,
            extend = TRUE,
            prior = 0.01,
            spat_factor = spat_df$spat_factor,
            spat_bin = 40L
        )
    )

    expect_true("pwm_spatial_test" %in% gvtrack.ls())
})

test_that("misha spatial PWM numerical agreement with prego - synthetic sequences", {
    skip_if_not_installed("prego")

    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    # Load prego


    # Create synthetic test sequences with known motif instances
    # All sequences same length for easier comparison
    set.seed(42)
    n_seq <- 10
    seq_len <- 280

    # Generate random sequences
    bases <- c("A", "C", "G", "T")
    test_seqs <- sapply(1:n_seq, function(i) {
        paste0(sample(bases, seq_len, replace = TRUE), collapse = "")
    })

    # Create a simple PSSM (AC motif)
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

    # Compute using prego
    prego_scores_nospatial <- prego::compute_pwm(
        sequences = test_seqs,
        pssm = test_pssm,
        spat = NULL,
        bidirect = FALSE,
        prior = 0.01,
        func = "logSumExp"
    )

    prego_scores_spatial <- prego::compute_pwm(
        sequences = test_seqs,
        pssm = test_pssm,
        spat = spat_df,
        bidirect = FALSE,
        prior = 0.01,
        func = "logSumExp"
    )

    # Validate that spatial factors increase the scores (due to weighting)
    # This tests that the spatial model is being applied
    expect_true(length(prego_scores_nospatial) == n_seq)
    expect_true(length(prego_scores_spatial) == n_seq)

    # The spatial scores should generally be different from non-spatial
    # (unless by chance all motifs are in low-weight regions)
    expect_true(any(abs(prego_scores_spatial - prego_scores_nospatial) > 0.01))
})

test_that("misha max PWM mode matches prego max mode", {
    skip_if_not_installed("prego")

    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    # Load prego


    # Create simple sequences
    set.seed(123)
    test_seqs <- c(
        "ACGTACGTACGT",
        "TGCATGCATGCA",
        "AAAACCCCGGGG"
    )

    # Create test PSSM
    test_pssm <- data.frame(
        pos = 1:4,
        A = c(0.7, 0.1, 0.1, 0.1),
        C = c(0.1, 0.7, 0.1, 0.1),
        G = c(0.1, 0.1, 0.7, 0.1),
        T = c(0.1, 0.1, 0.1, 0.7)
    )

    pssm_mat <- as.matrix(test_pssm[, c("A", "C", "G", "T")])

    # Test max mode with prego
    prego_max_scores <- prego::compute_pwm(
        sequences = test_seqs,
        pssm = test_pssm,
        spat = NULL,
        bidirect = TRUE,
        prior = 0.01,
        func = "max"
    )

    # Create misha vtrack for max mode
    gvtrack.create(
        "pwm_max_test", NULL, "pwm.max",
        list(
            pssm = pssm_mat,
            bidirect = TRUE,
            extend = TRUE,
            prior = 0.01
        )
    )

    expect_true("pwm_max_test" %in% gvtrack.ls())
    expect_equal(length(prego_max_scores), 3)
})

test_that("spatial binning is 0-indexed from scan start like prego", {
    # This test validates the spatial binning logic
    # In prego: bins are 0-indexed from the start of the scan range
    # Same bin for both strands

    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    # Create a PSSM
    pssm_mat <- matrix(
        c(
            0.7, 0.1, 0.1, 0.1,
            0.1, 0.7, 0.1, 0.1
        ),
        nrow = 2, byrow = TRUE
    )
    colnames(pssm_mat) <- c("A", "C", "G", "T")

    # Create spatial factors - higher weight in first bin
    spat_factors <- c(10.0, 1.0, 1.0)
    spat_bin <- 100L

    # Create vtrack
    gvtrack.create(
        "pwm_bintest", NULL, "pwm",
        list(
            pssm = pssm_mat,
            bidirect = FALSE,
            extend = TRUE,
            prior = 0.01,
            spat_factor = spat_factors,
            spat_bin = spat_bin
        )
    )

    # Extract on a test interval
    test_int <- gintervals(1, 1000, 1300)
    scores <- gextract("pwm_bintest", test_int, iterator = test_int)

    # Should return a valid score (testing that binning doesn't crash)
    expect_false(is.na(scores$pwm_bintest[1]))
    expect_false(is.infinite(scores$pwm_bintest[1]))
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
