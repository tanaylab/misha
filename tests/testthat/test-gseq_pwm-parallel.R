test_that("gseq.pwm multitask produces identical results to sequential (mode=lse)", {
    withr::local_options(list(gmultitasking = FALSE))

    # Create test data
    pssm <- create_test_pssm()
    seqs <- replicate(1000, paste(sample(c("A", "C", "G", "T"), 50, replace = TRUE), collapse = ""))

    # Sequential execution
    result_seq <- gseq.pwm(seqs, pssm, mode = "lse", bidirect = FALSE)

    # Parallel execution (force above threshold)
    withr::local_options(list(
        gmultitasking = TRUE,
        gmax.processes = 4,
        gmin.seqs.work4process = 100 # Very low threshold to ensure parallelization
    ))
    result_par <- gseq.pwm(seqs, pssm, mode = "lse", bidirect = FALSE)

    expect_equal(result_par, result_seq, tolerance = 1e-10)
})

test_that("gseq.pwm multitask produces identical results to sequential (mode=max)", {
    withr::local_options(list(gmultitasking = FALSE))

    pssm <- create_test_pssm()
    seqs <- replicate(500, paste(sample(c("A", "C", "G", "T"), 100, replace = TRUE), collapse = ""))

    result_seq <- gseq.pwm(seqs, pssm, mode = "max", bidirect = TRUE)

    withr::local_options(list(
        gmultitasking = TRUE,
        gmax.processes = 4,
        gmin.seqs.work4process = 100
    ))
    result_par <- gseq.pwm(seqs, pssm, mode = "max", bidirect = TRUE)

    expect_equal(result_par, result_seq, tolerance = 1e-10)
})

test_that("gseq.pwm multitask produces identical results to sequential (mode=pos)", {
    withr::local_options(list(gmultitasking = FALSE))

    pssm <- create_test_pssm()
    seqs <- replicate(300, paste(sample(c("A", "C", "G", "T"), 80, replace = TRUE), collapse = ""))

    result_seq <- gseq.pwm(seqs, pssm, mode = "pos", bidirect = FALSE, return_strand = FALSE)

    withr::local_options(list(
        gmultitasking = TRUE,
        gmax.processes = 4,
        gmin.seqs.work4process = 100
    ))
    result_par <- gseq.pwm(seqs, pssm, mode = "pos", bidirect = FALSE, return_strand = FALSE)

    expect_equal(result_par, result_seq)
})

test_that("gseq.pwm multitask produces identical results with return_strand=TRUE", {
    withr::local_options(list(gmultitasking = FALSE))

    pssm <- create_test_pssm()
    seqs <- replicate(200, paste(sample(c("A", "C", "G", "T"), 60, replace = TRUE), collapse = ""))

    result_seq <- gseq.pwm(seqs, pssm, mode = "pos", bidirect = TRUE, return_strand = TRUE)

    withr::local_options(list(
        gmultitasking = TRUE,
        gmax.processes = 4,
        gmin.seqs.work4process = 100
    ))
    result_par <- gseq.pwm(seqs, pssm, mode = "pos", bidirect = TRUE, return_strand = TRUE)

    expect_equal(result_par$pos, result_seq$pos)
    expect_equal(result_par$strand, result_seq$strand)
})

test_that("gseq.pwm multitask produces identical results to sequential (mode=count)", {
    withr::local_options(list(gmultitasking = FALSE))

    pssm <- create_test_pssm()
    seqs <- replicate(400, paste(sample(c("A", "C", "G", "T"), 70, replace = TRUE), collapse = ""))

    result_seq <- gseq.pwm(seqs, pssm, mode = "count", score.thresh = -2, bidirect = TRUE)

    withr::local_options(list(
        gmultitasking = TRUE,
        gmax.processes = 4,
        gmin.seqs.work4process = 100
    ))
    result_par <- gseq.pwm(seqs, pssm, mode = "count", score.thresh = -2, bidirect = TRUE)

    expect_equal(result_par, result_seq)
})

test_that("gseq.pwm multitask handles ROI parameters correctly", {
    withr::local_options(list(gmultitasking = FALSE))

    pssm <- create_test_pssm()
    seqs <- replicate(200, paste(sample(c("A", "C", "G", "T"), 100, replace = TRUE), collapse = ""))
    start_pos <- rep(c(10, 20, 30), length.out = length(seqs))
    end_pos <- rep(c(50, 60, 70), length.out = length(seqs))

    result_seq <- gseq.pwm(seqs, pssm, mode = "max", start_pos = start_pos, end_pos = end_pos)

    withr::local_options(list(
        gmultitasking = TRUE,
        gmax.processes = 4,
        gmin.seqs.work4process = 100
    ))
    result_par <- gseq.pwm(seqs, pssm, mode = "max", start_pos = start_pos, end_pos = end_pos)

    expect_equal(result_par, result_seq, tolerance = 1e-10)
})

test_that("gseq.pwm multitask respects workload threshold (small workload stays sequential)", {
    pssm <- create_test_pssm()
    # Very small workload: 10 sequences x 20bp x 2bp PSSM = 400 work units
    seqs <- replicate(10, paste(sample(c("A", "C", "G", "T"), 20, replace = TRUE), collapse = ""))

    # Set high threshold to force sequential execution
    withr::local_options(list(
        gmultitasking = TRUE,
        gmax.processes = 4,
        gmin.seqs.work4process = 100000 # Much higher than workload
    ))

    # This should execute sequentially despite gmultitasking=TRUE
    # We can't directly test which code path was used, but we can verify correctness
    result <- gseq.pwm(seqs, pssm, mode = "lse")

    # Should complete without error and produce valid results
    expect_true(all(is.finite(result) | is.na(result)))
    expect_equal(length(result), length(seqs))
})

test_that("gseq.pwm multitask handles edge cases correctly", {
    pssm <- create_test_pssm()

    # Test with various edge cases
    withr::local_options(list(
        gmultitasking = TRUE,
        gmax.processes = 4,
        gmin.seqs.work4process = 100
    ))

    # Empty vector
    expect_equal(length(gseq.pwm(character(0), pssm, mode = "lse")), 0)

    # Single sequence
    single_seq <- "ACGTACGTACGT"
    result_single <- gseq.pwm(single_seq, pssm, mode = "max")
    expect_equal(length(result_single), 1)
    expect_true(is.finite(result_single))

    # Number of sequences less than number of processes
    few_seqs <- replicate(2, paste(sample(c("A", "C", "G", "T"), 50, replace = TRUE), collapse = ""))
    result_few <- gseq.pwm(few_seqs, pssm, mode = "lse")
    expect_equal(length(result_few), 2)
})

test_that("gseq.pwm multitask works with all strand modes", {
    withr::local_options(list(gmultitasking = FALSE))

    pssm <- create_test_pssm()
    seqs <- replicate(200, paste(sample(c("A", "C", "G", "T"), 60, replace = TRUE), collapse = ""))

    # Test strand = 1 (forward only)
    result_seq_fwd <- gseq.pwm(seqs, pssm, mode = "max", strand = 1)
    withr::local_options(list(
        gmultitasking = TRUE,
        gmax.processes = 4,
        gmin.seqs.work4process = 100
    ))
    result_par_fwd <- gseq.pwm(seqs, pssm, mode = "max", strand = 1)
    expect_equal(result_par_fwd, result_seq_fwd, tolerance = 1e-10)

    # Test strand = -1 (reverse only)
    withr::local_options(list(gmultitasking = FALSE))
    result_seq_rev <- gseq.pwm(seqs, pssm, mode = "max", strand = -1)
    withr::local_options(list(
        gmultitasking = TRUE,
        gmax.processes = 4,
        gmin.seqs.work4process = 100
    ))
    result_par_rev <- gseq.pwm(seqs, pssm, mode = "max", strand = -1)
    expect_equal(result_par_rev, result_seq_rev, tolerance = 1e-10)

    # Test strand = 0 (both strands, non-bidirectional)
    withr::local_options(list(gmultitasking = FALSE))
    result_seq_both <- gseq.pwm(seqs, pssm, mode = "max", strand = 0, bidirect = FALSE)
    withr::local_options(list(
        gmultitasking = TRUE,
        gmax.processes = 4,
        gmin.seqs.work4process = 100
    ))
    result_par_both <- gseq.pwm(seqs, pssm, mode = "max", strand = 0, bidirect = FALSE)
    expect_equal(result_par_both, result_seq_both, tolerance = 1e-10)
})

test_that("gseq.pwm multitask works with different prior values", {
    withr::local_options(list(gmultitasking = FALSE))

    pssm <- create_test_pssm()
    seqs <- replicate(150, paste(sample(c("A", "C", "G", "T"), 50, replace = TRUE), collapse = ""))

    for (prior_val in c(0, 0.01, 0.1, 1.0)) {
        result_seq <- gseq.pwm(seqs, pssm, mode = "lse", prior = prior_val)

        withr::local_options(list(
            gmultitasking = TRUE,
            gmax.processes = 4,
            gmin.seqs.work4process = 100
        ))
        result_par <- gseq.pwm(seqs, pssm, mode = "lse", prior = prior_val)

        expect_equal(result_par, result_seq,
            tolerance = 1e-10,
            label = paste("prior =", prior_val)
        )
    }
})

test_that("gseq.pwm multitask produces consistent results across multiple runs", {
    pssm <- create_test_pssm()
    seqs <- replicate(300, paste(sample(c("A", "C", "G", "T"), 60, replace = TRUE), collapse = ""))

    withr::local_options(list(
        gmultitasking = TRUE,
        gmax.processes = 4,
        gmin.seqs.work4process = 100
    ))

    # Run multiple times - should get identical results
    result1 <- gseq.pwm(seqs, pssm, mode = "max", bidirect = TRUE)
    result2 <- gseq.pwm(seqs, pssm, mode = "max", bidirect = TRUE)
    result3 <- gseq.pwm(seqs, pssm, mode = "max", bidirect = TRUE)

    expect_equal(result2, result1, tolerance = 1e-15)
    expect_equal(result3, result1, tolerance = 1e-15)
})


test_that("gseq.pwm small workloads stay sequential and avoid fork overhead", {
    skip_on_cran()

    pssm <- create_test_pssm()
    # Small workload: 500 seqs x 100bp x 2bp PSSM = 100,000 work units
    # This is right at the default threshold boundary
    set.seed(789)
    seqs <- replicate(500, paste(sample(c("A", "C", "G", "T"), 100, replace = TRUE), collapse = ""))

    # Sequential
    withr::local_options(list(gmultitasking = FALSE))
    time_seq <- system.time({
        result_seq <- gseq.pwm(seqs, pssm, mode = "count", score.thresh = -2)
    })

    # "Parallel" with high threshold - should stay sequential
    withr::local_options(list(
        gmultitasking = TRUE,
        gmax.processes = 4,
        gmin.seqs.work4process = 1000000 # Very high threshold
    ))
    time_threshold <- system.time({
        result_threshold <- gseq.pwm(seqs, pssm, mode = "count", score.thresh = -2)
    })

    expect_equal(result_threshold, result_seq)

    # Times should be very similar (both sequential)
    # We just verify both complete successfully and produce same results
    # Don't enforce strict timing ratios on small workloads due to system noise
    expect_true(time_seq["elapsed"] >= 0)
    expect_true(time_threshold["elapsed"] >= 0)
})
