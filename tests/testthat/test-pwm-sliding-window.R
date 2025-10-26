# Test sliding window optimization for PWM virtual tracks
# These tests verify that the fast path produces identical results to the legacy path

test_that("PWM sliding window gives same results as legacy - TOTAL_LIKELIHOOD", {
    pwm <- rbind(
        c(0.25, 0.25, 0.25, 0.25),
        c(0.8, 0.1, 0.05, 0.05),
        c(0.1, 0.7, 0.1, 0.1),
        c(0.1, 0.1, 0.7, 0.1),
        c(0.05, 0.05, 0.1, 0.8),
        c(0.25, 0.25, 0.25, 0.25)
    )
    colnames(pwm) <- c("A", "C", "G", "T")

    gvtrack.create("pwm_test", "seq", func = "pwm", params = list(pssm = pwm))

    # Test with dense iterator (step=1)
    # This should trigger the sliding window optimization
    result1 <- gextract("pwm_test", gintervals(1, 10000, 11000), iterator = 1)
    expect_regression(result1, "pwm_sliding_window_test_1")

    # Verify we got results
    expect_true(nrow(result1) > 0)
    expect_false(any(is.na(result1$pwm_test)))

    # Test that consecutive intervals give reasonable results
    # (the optimization should maintain correctness)
    expect_true(all(is.finite(result1$pwm_test)))

    gvtrack.rm("pwm_test")
})

test_that("PWM sliding window works with MAX_LIKELIHOOD mode", {
    pwm <- rbind(
        c(0.8, 0.1, 0.05, 0.05),
        c(0.1, 0.7, 0.1, 0.1),
        c(0.1, 0.1, 0.7, 0.1),
        c(0.05, 0.05, 0.1, 0.8)
    )
    colnames(pwm) <- c("A", "C", "G", "T")

    gvtrack.create("pwm_max_test", "seq", func = "pwm.max", params = list(pssm = pwm))
    result <- gextract("pwm_max_test", gintervals(1, 10000, 11000), iterator = 1)
    expect_regression(result, "pwm_sliding_window_test_2")

    expect_true(nrow(result) > 0)
    expect_false(any(is.na(result$pwm_max_test)))
    expect_true(all(is.finite(result$pwm_max_test)))

    gvtrack.rm("pwm_max_test")
})

test_that("PWM sliding window works with MOTIF_COUNT mode", {
    pwm <- rbind(
        c(0.9, 0.03, 0.03, 0.04),
        c(0.03, 0.9, 0.03, 0.04),
        c(0.03, 0.03, 0.9, 0.04),
        c(0.04, 0.03, 0.03, 0.9)
    )
    colnames(pwm) <- c("A", "C", "G", "T")

    gvtrack.create("pwm_count_test", "seq",
        func = "pwm.count",
        params = list(pssm = pwm, threshold = -10)
    )
    result <- gextract("pwm_count_test", gintervals(1, 10000, 11000), iterator = 1)
    expect_regression(result, "pwm_sliding_window_test_3")

    expect_true(nrow(result) > 0)
    expect_false(any(is.na(result$pwm_count_test)))
    expect_true(all(result$pwm_count_test >= 0))

    gvtrack.rm("pwm_count_test")
})

test_that("PWM sliding window cache invalidates on chromosome change", {
    pwm <- rbind(
        c(0.7, 0.1, 0.1, 0.1),
        c(0.1, 0.7, 0.1, 0.1),
        c(0.1, 0.1, 0.7, 0.1),
        c(0.1, 0.1, 0.1, 0.7)
    )
    colnames(pwm) <- c("A", "C", "G", "T")

    gvtrack.create("pwm_chrom_test", "seq", func = "pwm", params = list(pssm = pwm))

    # Extract from multiple chromosomes
    intervals <- rbind(
        gintervals(1, 1000, 2000),
        gintervals(2, 1000, 2000),
        gintervals(1, 2000, 3000)
    )

    result <- gextract("pwm_chrom_test", intervals, iterator = 1)
    expect_regression(result, "pwm_sliding_window_test_4")
    # Should have results for all three intervals (1000 rows each)
    expect_true(nrow(result) > 0)
    expect_false(any(is.na(result$pwm_chrom_test)))

    gvtrack.rm("pwm_chrom_test")
})

test_that("PWM sliding window handles different iterator steps", {
    pwm <- rbind(
        c(0.6, 0.15, 0.15, 0.1),
        c(0.15, 0.6, 0.15, 0.1),
        c(0.15, 0.15, 0.6, 0.1)
    )
    colnames(pwm) <- c("A", "C", "G", "T")

    gvtrack.create("pwm_iter_test", "seq", func = "pwm", params = list(pssm = pwm))

    # Test with different iterator steps
    result_1 <- gextract("pwm_iter_test", gintervals(1, 10000, 11000), iterator = 1)
    result_10 <- gextract("pwm_iter_test", gintervals(1, 10000, 11000), iterator = 10)
    result_100 <- gextract("pwm_iter_test", gintervals(1, 10000, 11000), iterator = 100)
    expect_regression(result_1, "pwm_sliding_window_test_5")
    expect_regression(result_10, "pwm_sliding_window_test_6")
    expect_regression(result_100, "pwm_sliding_window_test_7")

    # Different step sizes should give different number of results
    expect_true(nrow(result_1) > nrow(result_10))
    expect_true(nrow(result_10) > nrow(result_100))

    # But all should be valid
    expect_false(any(is.na(result_1$pwm_iter_test)))
    expect_false(any(is.na(result_10$pwm_iter_test)))
    expect_false(any(is.na(result_100$pwm_iter_test)))

    gvtrack.rm("pwm_iter_test")
})

test_that("PWM sliding window works with bidirectional PSSMs", {
    pwm <- rbind(
        c(0.8, 0.1, 0.05, 0.05),
        c(0.1, 0.7, 0.1, 0.1),
        c(0.1, 0.1, 0.7, 0.1),
        c(0.05, 0.05, 0.1, 0.8)
    )
    colnames(pwm) <- c("A", "C", "G", "T")

    # Test both strands
    gvtrack.create("pwm_bidir_test", "seq",
        func = "pwm",
        params = list(pssm = pwm, bidirect = TRUE)
    )

    result <- gextract("pwm_bidir_test", gintervals(1, 10000, 11000), iterator = 1)
    expect_regression(result, "pwm_sliding_window_test_8")
    expect_true(nrow(result) > 0)
    expect_false(any(is.na(result$pwm_bidir_test)))

    gvtrack.rm("pwm_bidir_test")
})


test_that("PWM sliding window matches spatial weighting baseline", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    pssm <- rbind(
        c(0.7, 0.1, 0.1, 0.1),
        c(0.1, 0.7, 0.1, 0.1),
        c(0.1, 0.1, 0.7, 0.1),
        c(0.1, 0.1, 0.1, 0.7)
    )
    colnames(pssm) <- c("A", "C", "G", "T")

    gvtrack.create(
        "pwm_slide_spatial", "seq",
        func = "pwm",
        params = list(
            pssm = pssm,
            spat_factor = c(0.2, 1.0, 3.0, 1.0, 0.5),
            spat_bin = 10L
        )
    )

    intervals <- rbind(
        gintervals(1, 10000, 10100),
        gintervals(1, 10001, 10101),
        gintervals(1, 10002, 10102),
        gintervals(1, 10003, 10103)
    )

    per_interval <- unlist(lapply(seq_len(nrow(intervals)), function(idx) {
        res <- gextract("pwm_slide_spatial", intervals[idx, , drop = FALSE], iterator = 1)
        res$pwm_slide_spatial
    }), use.names = FALSE)

    sliding_results <- gextract("pwm_slide_spatial", intervals, iterator = 1)$pwm_slide_spatial

    expect_equal(sliding_results, per_interval, tolerance = 1e-6)
    expect_regression(sliding_results, "pwm_sliding_window_test_9")
})

test_that("PWM sliding window works with iterator=20 and shifts", {
    # Create a longer PSSM (20bp) to match benchmark
    set.seed(42)
    pssm <- matrix(runif(80), nrow = 20, ncol = 4)
    pssm <- pssm / rowSums(pssm)
    colnames(pssm) <- c("A", "C", "G", "T")

    gvtrack.create("pwm_iter20_shift", "seq",
        func = "pwm",
        params = list(pssm = pssm, bidirect = TRUE, prior = 0.01)
    )
    gvtrack.iterator("pwm_iter20_shift", sshift = -250, eshift = 250)

    # Test with iterator=20 on a moderately sized region
    result <- gextract("pwm_iter20_shift", gintervals(1, 1000000, 1100000), iterator = 20)
    expect_regression(result, "pwm_sliding_window_test_10")
    # Should get approximately 100000/20 = 5000 positions (plus shifts)
    expect_true(nrow(result) > 4000)
    expect_false(any(is.na(result$pwm_iter20_shift)))
    expect_true(all(is.finite(result$pwm_iter20_shift)))

    gvtrack.rm("pwm_iter20_shift")
})

test_that("PWM sliding window works with iterator=20 without shifts", {
    set.seed(42)
    pssm <- matrix(runif(80), nrow = 20, ncol = 4)
    pssm <- pssm / rowSums(pssm)
    colnames(pssm) <- c("A", "C", "G", "T")

    gvtrack.create("pwm_iter20_no_shift", "seq",
        func = "pwm",
        params = list(pssm = pssm, bidirect = TRUE, prior = 0.01)
    )
    # No gvtrack.iterator call - uses default (no shifts)

    result <- gextract("pwm_iter20_no_shift", gintervals(1, 1000000, 1100000), iterator = 20)
    expect_regression(result, "pwm_sliding_window_test_11")
    # Should get approximately 100000/20 = 5000 positions
    expect_true(nrow(result) >= 4500)
    expect_true(nrow(result) <= 5500)
    expect_false(any(is.na(result$pwm_iter20_no_shift)))
    expect_true(all(is.finite(result$pwm_iter20_no_shift)))

    gvtrack.rm("pwm_iter20_no_shift")
})

test_that("PWM sliding window works with and without shifts", {
    set.seed(123)
    pssm <- matrix(runif(40), nrow = 10, ncol = 4)
    pssm <- pssm / rowSums(pssm)
    colnames(pssm) <- c("A", "C", "G", "T")

    # Create track with shifts
    gvtrack.create("pwm_with_shift", "seq",
        func = "pwm",
        params = list(pssm = pssm, bidirect = TRUE)
    )
    gvtrack.iterator("pwm_with_shift", sshift = -100, eshift = 100)

    # Create track without shifts
    gvtrack.create("pwm_no_shift", "seq",
        func = "pwm",
        params = list(pssm = pssm, bidirect = TRUE)
    )

    interval <- gintervals(1, 1000000, 1010000)

    result_with_shift <- gextract("pwm_with_shift", interval, iterator = 50)
    result_no_shift <- gextract("pwm_no_shift", interval, iterator = 50)
    expect_regression(result_with_shift, "pwm_sliding_window_test_12")
    expect_regression(result_no_shift, "pwm_sliding_window_test_13")
    # Both should return the same number of rows (shifts don't change output coordinates)
    expect_equal(nrow(result_with_shift), nrow(result_no_shift))

    # Both should cover the same coordinate range
    expect_equal(range(result_with_shift$start), range(result_no_shift$start))
    expect_equal(range(result_with_shift$end), range(result_no_shift$end))

    # All results should be valid
    expect_false(any(is.na(result_with_shift$pwm_with_shift)))
    expect_false(any(is.na(result_no_shift$pwm_no_shift)))

    gvtrack.rm("pwm_with_shift")
    gvtrack.rm("pwm_no_shift")
})

test_that("PWM sliding window handles different shift values", {
    set.seed(456)
    pssm <- matrix(runif(32), nrow = 8, ncol = 4)
    pssm <- pssm / rowSums(pssm)
    colnames(pssm) <- c("A", "C", "G", "T")

    # Test with large shifts (like benchmark: -250, 250)
    gvtrack.create("pwm_large_shift", "seq",
        func = "pwm",
        params = list(pssm = pssm)
    )
    gvtrack.iterator("pwm_large_shift", sshift = -250, eshift = 250)

    result_large <- gextract("pwm_large_shift", gintervals(1, 1000000, 1001000), iterator = 20)
    expect_regression(result_large, "pwm_sliding_window_test_14")
    # Should have ~50 positions (1000bp / 20)
    expect_true(nrow(result_large) >= 40)
    expect_false(any(is.na(result_large$pwm_large_shift)))

    # Test with small shifts
    gvtrack.create("pwm_small_shift", "seq",
        func = "pwm",
        params = list(pssm = pssm)
    )
    gvtrack.iterator("pwm_small_shift", sshift = -50, eshift = 50)

    result_small <- gextract("pwm_small_shift", gintervals(1, 1000000, 1001000), iterator = 20)
    expect_regression(result_small, "pwm_sliding_window_test_15")
    # Should also have ~50 positions (shifts don't change output coordinates)
    expect_true(nrow(result_small) >= 40)
    expect_false(any(is.na(result_small$pwm_small_shift)))

    # Both should give same number of results (shifts provide context, not extra output)
    expect_equal(nrow(result_large), nrow(result_small))

    gvtrack.rm("pwm_large_shift")
    gvtrack.rm("pwm_small_shift")
})

test_that("PWM sliding window with iterator=20 works across multiple intervals", {
    set.seed(789)
    pssm <- matrix(runif(48), nrow = 12, ncol = 4)
    pssm <- pssm / rowSums(pssm)
    colnames(pssm) <- c("A", "C", "G", "T")

    gvtrack.create("pwm_multi_iter20", "seq",
        func = "pwm",
        params = list(pssm = pssm)
    )
    gvtrack.iterator("pwm_multi_iter20", sshift = -100, eshift = 100)

    # Test with multiple consecutive intervals (similar to benchmark Test 4)
    intervals <- do.call(rbind, lapply(0:9, function(i) {
        gintervals(1, 1000000 + i * 10000, 1000000 + (i + 1) * 10000)
    }))

    result <- gextract("pwm_multi_iter20", intervals, iterator = 20)
    expect_regression(result, "pwm_sliding_window_test_16")

    # Should have results from all intervals
    expect_true(nrow(result) > 0)
    expect_false(any(is.na(result$pwm_multi_iter20)))
    expect_true(all(is.finite(result$pwm_multi_iter20)))

    # Each ~10kb interval with iterator=20 should give ~500 positions
    # 10 intervals * 500 = ~5000 positions (plus shift extensions)
    expect_true(nrow(result) >= 4000)

    gvtrack.rm("pwm_multi_iter20")
})

test_that("PWM sliding window with various iterator values and shifts", {
    set.seed(321)
    pssm <- matrix(runif(28), nrow = 7, ncol = 4)
    pssm <- pssm / rowSums(pssm)
    colnames(pssm) <- c("A", "C", "G", "T")

    gvtrack.create("pwm_var_iter", "seq",
        func = "pwm",
        params = list(pssm = pssm, bidirect = TRUE)
    )
    gvtrack.iterator("pwm_var_iter", sshift = -150, eshift = 150)

    interval <- gintervals(1, 1000000, 1050000)

    # Test multiple iterator values
    result_5 <- gextract("pwm_var_iter", interval, iterator = 5)
    result_20 <- gextract("pwm_var_iter", interval, iterator = 20)
    result_50 <- gextract("pwm_var_iter", interval, iterator = 50)
    result_100 <- gextract("pwm_var_iter", interval, iterator = 100)
    expect_regression(result_5, "pwm_sliding_window_test_17")
    expect_regression(result_20, "pwm_sliding_window_test_18")
    expect_regression(result_50, "pwm_sliding_window_test_19")
    expect_regression(result_100, "pwm_sliding_window_test_20")

    # Smaller iterator values should give more results
    expect_true(nrow(result_5) > nrow(result_20))
    expect_true(nrow(result_20) > nrow(result_50))
    expect_true(nrow(result_50) > nrow(result_100))

    # All should be valid
    expect_false(any(is.na(result_5$pwm_var_iter)))
    expect_false(any(is.na(result_20$pwm_var_iter)))
    expect_false(any(is.na(result_50$pwm_var_iter)))
    expect_false(any(is.na(result_100$pwm_var_iter)))

    gvtrack.rm("pwm_var_iter")
})

test_that("PWM sliding window MAX_LIKELIHOOD mode works with iterator=20 and shifts", {
    set.seed(654)
    pssm <- matrix(runif(32), nrow = 8, ncol = 4)
    pssm <- pssm / rowSums(pssm)
    colnames(pssm) <- c("A", "C", "G", "T")

    gvtrack.create("pwm_max_iter20", "seq",
        func = "pwm.max",
        params = list(pssm = pssm, bidirect = TRUE)
    )
    gvtrack.iterator("pwm_max_iter20", sshift = -200, eshift = 200)

    result <- gextract("pwm_max_iter20", gintervals(1, 1000000, 1100000), iterator = 20)
    expect_regression(result, "pwm_sliding_window_test_21")

    expect_true(nrow(result) > 0)
    expect_false(any(is.na(result$pwm_max_iter20)))
    expect_true(all(is.finite(result$pwm_max_iter20)))

    gvtrack.rm("pwm_max_iter20")
})

test_that("PWM sliding window MOTIF_COUNT mode works with iterator=20 and shifts", {
    set.seed(987)
    pssm <- matrix(runif(24), nrow = 6, ncol = 4)
    pssm <- pssm / rowSums(pssm)
    colnames(pssm) <- c("A", "C", "G", "T")

    gvtrack.create("pwm_count_iter20", "seq",
        func = "pwm.count",
        params = list(pssm = pssm, threshold = -10)
    )
    gvtrack.iterator("pwm_count_iter20", sshift = -100, eshift = 100)

    result <- gextract("pwm_count_iter20", gintervals(1, 1000000, 1100000), iterator = 20)
    expect_regression(result, "pwm_sliding_window_test_22")
    expect_true(nrow(result) > 0)
    expect_false(any(is.na(result$pwm_count_iter20)))
    expect_true(all(result$pwm_count_iter20 >= 0))

    gvtrack.rm("pwm_count_iter20")
})

test_that("PWM sliding window MAX_POS mode works with iterator=1", {
    set.seed(555)
    pssm <- matrix(runif(32), nrow = 8, ncol = 4)
    pssm <- pssm / rowSums(pssm)
    colnames(pssm) <- c("A", "C", "G", "T")

    gvtrack.create("pwm_pos_test", "seq",
        func = "pwm.max.pos",
        params = list(pssm = pssm, bidirect = TRUE)
    )

    result <- gextract("pwm_pos_test", gintervals(1, 10000, 11000), iterator = 1)
    expect_regression(result, "pwm_sliding_window_test_34")

    expect_true(nrow(result) > 0)
    expect_false(any(is.na(result$pwm_pos_test)))
    # Positions should be integers (1-based, positive for forward, negative for reverse)
    expect_true(all(result$pwm_pos_test == round(result$pwm_pos_test)))

    gvtrack.rm("pwm_pos_test")
})

test_that("PWM sliding window MAX_POS mode works with iterator=20 and shifts", {
    set.seed(666)
    pssm <- matrix(runif(40), nrow = 10, ncol = 4)
    pssm <- pssm / rowSums(pssm)
    colnames(pssm) <- c("A", "C", "G", "T")

    gvtrack.create("pwm_pos_iter20", "seq",
        func = "pwm.max.pos",
        params = list(pssm = pssm, bidirect = TRUE)
    )
    gvtrack.iterator("pwm_pos_iter20", sshift = -150, eshift = 150)
    result <- gextract("pwm_pos_iter20", gintervals(1, 1000000, 1100000), iterator = 20)
    expect_regression(result, "pwm_sliding_window_test_35")

    expect_true(nrow(result) > 0)
    expect_false(any(is.na(result$pwm_pos_iter20)))
    # Positions should be integers
    expect_true(all(result$pwm_pos_iter20 == round(result$pwm_pos_iter20)))

    gvtrack.rm("pwm_pos_iter20")
})

test_that("PWM sliding window MAX_POS mode without shifts works with iterator=50", {
    set.seed(777)
    pssm <- matrix(runif(28), nrow = 7, ncol = 4)
    pssm <- pssm / rowSums(pssm)
    colnames(pssm) <- c("A", "C", "G", "T")

    gvtrack.create("pwm_pos_no_shift", "seq",
        func = "pwm.max.pos",
        params = list(pssm = pssm, bidirect = TRUE)
    )

    result <- gextract("pwm_pos_no_shift", gintervals(1, 1000000, 1050000), iterator = 50)

    expect_regression(result, "pwm_sliding_window_test_36")

    expect_true(nrow(result) > 0)
    expect_false(any(is.na(result$pwm_pos_no_shift)))
    # All positions should be integers
    expect_true(all(result$pwm_pos_no_shift == round(result$pwm_pos_no_shift)))

    gvtrack.rm("pwm_pos_no_shift")
})

test_that("PWM sliding window MAX_POS mode with various iterators and shifts", {
    set.seed(888)
    pssm <- matrix(runif(36), nrow = 9, ncol = 4)
    pssm <- pssm / rowSums(pssm)
    colnames(pssm) <- c("A", "C", "G", "T")

    gvtrack.create("pwm_pos_var", "seq",
        func = "pwm.max.pos",
        params = list(pssm = pssm)
    )
    gvtrack.iterator("pwm_pos_var", sshift = -100, eshift = 100)

    interval <- gintervals(1, 1000000, 1050000)

    result_10 <- gextract("pwm_pos_var", interval, iterator = 10)
    result_100 <- gextract("pwm_pos_var", interval, iterator = 100)
    expect_regression(result_10, "pwm_sliding_window_test_37")
    expect_regression(result_100, "pwm_sliding_window_test_38")

    # Larger iterator should give fewer results
    expect_true(nrow(result_10) > nrow(result_100))

    # All should be valid integers
    expect_false(any(is.na(result_10$pwm_pos_var)))
    expect_false(any(is.na(result_100$pwm_pos_var)))
    expect_true(all(result_10$pwm_pos_var == round(result_10$pwm_pos_var)))
    expect_true(all(result_100$pwm_pos_var == round(result_100$pwm_pos_var)))

    gvtrack.rm("pwm_pos_var")
})

test_that("PWM shifts produce correct values - verification with small region", {
    # Create a simple, reproducible PSSM
    set.seed(111)
    pssm <- matrix(runif(32), nrow = 8, ncol = 4)
    pssm <- pssm / rowSums(pssm)
    colnames(pssm) <- c("A", "C", "G", "T")

    # Create track with shifts
    gvtrack.create("pwm_verify_shift", "seq",
        func = "pwm",
        params = list(pssm = pssm, bidirect = TRUE)
    )
    gvtrack.iterator("pwm_verify_shift", sshift = -100, eshift = 100)

    # Create track without shifts
    gvtrack.create("pwm_verify_no_shift", "seq",
        func = "pwm",
        params = list(pssm = pssm, bidirect = TRUE)
    )

    # Extract a small region with large iterator (sum of shift magnitudes)
    interval <- gintervals(1, 1000000, 1001000)

    result_with_shift <- gextract("pwm_verify_shift", interval, iterator = 200)
    result_no_shift <- gextract("pwm_verify_no_shift", interval, iterator = 200)
    expect_regression(result_with_shift, "pwm_sliding_window_test_23")
    expect_regression(result_no_shift, "pwm_sliding_window_test_24")

    # Both should have same coordinates
    expect_equal(result_with_shift$start, result_no_shift$start)
    expect_equal(result_with_shift$end, result_no_shift$end)

    # Values should differ because shifts provide different context
    # At least some positions should have different scores
    expect_false(all(abs(result_with_shift$pwm_verify_shift - result_no_shift$pwm_verify_no_shift) < 1e-10))

    # But all values should be valid
    expect_false(any(is.na(result_with_shift$pwm_verify_shift)))
    expect_false(any(is.na(result_no_shift$pwm_verify_no_shift)))

    gvtrack.rm("pwm_verify_shift")
    gvtrack.rm("pwm_verify_no_shift")
})

test_that("PWM shifts correctness - values differ from no-shift baseline", {
    # Use a known PSSM to verify correctness
    set.seed(222)
    pssm <- matrix(runif(40), nrow = 10, ncol = 4)
    pssm <- pssm / rowSums(pssm)
    colnames(pssm) <- c("A", "C", "G", "T")

    # Create track with shifts
    gvtrack.create("pwm_dense_shift", "seq",
        func = "pwm",
        params = list(pssm = pssm)
    )
    gvtrack.iterator("pwm_dense_shift", sshift = -50, eshift = 50)

    # Create track without shifts
    gvtrack.create("pwm_dense_no_shift", "seq",
        func = "pwm",
        params = list(pssm = pssm)
    )

    # Extract on a small region with same iterator for both
    interval <- gintervals(1, 1000500, 1000600)

    result_shift <- gextract("pwm_dense_shift", interval, iterator = 1)
    result_no_shift <- gextract("pwm_dense_no_shift", interval, iterator = 1)
    expect_regression(result_shift, "pwm_sliding_window_test_25")
    expect_regression(result_no_shift, "pwm_sliding_window_test_26")
    # Both should have same number of positions
    expect_equal(nrow(result_shift), nrow(result_no_shift))

    # With shifts, values should differ from without shifts (shifts provide context)
    # At least some positions should differ
    expect_false(all(abs(result_shift$pwm_dense_shift - result_no_shift$pwm_dense_no_shift) < 1e-10))

    # All values should be valid
    expect_false(any(is.na(result_shift$pwm_dense_shift)))
    expect_false(any(is.na(result_no_shift$pwm_dense_no_shift)))

    gvtrack.rm("pwm_dense_shift")
    gvtrack.rm("pwm_dense_no_shift")
})

test_that("PWM with shifts works correctly with large iterator values", {
    # Create PSSM
    set.seed(333)
    pssm <- matrix(runif(48), nrow = 12, ncol = 4)
    pssm <- pssm / rowSums(pssm)
    colnames(pssm) <- c("A", "C", "G", "T")

    # Create track with large shifts
    gvtrack.create("pwm_shift_consistent", "seq",
        func = "pwm",
        params = list(pssm = pssm, bidirect = TRUE)
    )
    gvtrack.iterator("pwm_shift_consistent", sshift = -200, eshift = 200)

    interval <- gintervals(1, 1000000, 1010000)

    # Extract with different iterators - all should produce valid results
    result_iter1 <- gextract("pwm_shift_consistent", interval, iterator = 1)
    result_iter100 <- gextract("pwm_shift_consistent", interval, iterator = 100)
    result_iter400 <- gextract("pwm_shift_consistent", interval, iterator = 400)
    expect_regression(result_iter1, "pwm_sliding_window_test_27")
    expect_regression(result_iter100, "pwm_sliding_window_test_28")
    expect_regression(result_iter400, "pwm_sliding_window_test_29")
    # All should return valid, non-NA results
    expect_false(any(is.na(result_iter1$pwm_shift_consistent)))
    expect_false(any(is.na(result_iter100$pwm_shift_consistent)))
    expect_false(any(is.na(result_iter400$pwm_shift_consistent)))

    # Larger iterators should give fewer results
    expect_true(nrow(result_iter1) > nrow(result_iter100))
    expect_true(nrow(result_iter100) > nrow(result_iter400))

    # All results should be finite
    expect_true(all(is.finite(result_iter1$pwm_shift_consistent)))
    expect_true(all(is.finite(result_iter100$pwm_shift_consistent)))
    expect_true(all(is.finite(result_iter400$pwm_shift_consistent)))

    gvtrack.rm("pwm_shift_consistent")
})

test_that("PWM with shifts at boundary positions returns valid results", {
    # Create a small PSSM
    set.seed(444)
    pssm <- matrix(runif(20), nrow = 5, ncol = 4)
    pssm <- pssm / rowSums(pssm)
    colnames(pssm) <- c("A", "C", "G", "T")

    # Create track with shifts
    gvtrack.create("pwm_boundary_test", "seq",
        func = "pwm",
        params = list(pssm = pssm)
    )
    sshift <- -30
    eshift <- 30
    gvtrack.iterator("pwm_boundary_test", sshift = sshift, eshift = eshift)

    # Extract a small interval with iterator = sum of shift magnitudes
    interval <- gintervals(1, 1000000, 1000200)
    result <- gextract("pwm_boundary_test", interval, iterator = abs(sshift) + abs(eshift))
    expect_regression(result, "pwm_sliding_window_test_30")
    # Verify we got results
    expect_true(nrow(result) > 0)
    expect_false(any(is.na(result$pwm_boundary_test)))
    expect_true(all(is.finite(result$pwm_boundary_test)))

    # Expected number of positions: 200bp / 60 ≈ 3-4 positions
    expect_true(nrow(result) >= 2)
    expect_true(nrow(result) <= 5)

    gvtrack.rm("pwm_boundary_test")
})

# Regression tests - comparing against saved results from old version
test_that("PWM regression: iterator=1 with shifts", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    set.seed(100)
    pssm <- matrix(runif(80), nrow = 20, ncol = 4)
    pssm <- pssm / rowSums(pssm)
    colnames(pssm) <- c("A", "C", "G", "T")

    gvtrack.create("pwm_reg1", "seq",
        func = "pwm",
        params = list(pssm = pssm, bidirect = TRUE, prior = 0.01)
    )
    gvtrack.iterator("pwm_reg1", sshift = -250, eshift = 250)

    result <- gextract("pwm_reg1", gintervals(1, 1000000, 1100000), iterator = 1)
    expect_regression(result, "pwm_sliding_window_test_31")
    expect_regression(result, "pwm_sliding_window_iter1_with_shifts")
})

test_that("PWM regression: iterator=20 with shifts", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    set.seed(100)
    pssm <- matrix(runif(80), nrow = 20, ncol = 4)
    pssm <- pssm / rowSums(pssm)
    colnames(pssm) <- c("A", "C", "G", "T")

    gvtrack.create("pwm_reg20", "seq",
        func = "pwm",
        params = list(pssm = pssm, bidirect = TRUE, prior = 0.01)
    )
    gvtrack.iterator("pwm_reg20", sshift = -250, eshift = 250)

    result <- gextract("pwm_reg20", gintervals(1, 1000000, 1100000), iterator = 20)
    expect_regression(result, "pwm_sliding_window_test_32")
    expect_regression(result, "pwm_sliding_window_iter20_with_shifts")
})

test_that("PWM regression: iterator=20 without shifts", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    set.seed(100)
    pssm <- matrix(runif(80), nrow = 20, ncol = 4)
    pssm <- pssm / rowSums(pssm)
    colnames(pssm) <- c("A", "C", "G", "T")

    gvtrack.create("pwm_reg_no_shift", "seq",
        func = "pwm",
        params = list(pssm = pssm, bidirect = TRUE, prior = 0.01)
    )

    result <- gextract("pwm_reg_no_shift", gintervals(1, 1000000, 1100000), iterator = 20)
    expect_regression(result, "pwm_sliding_window_test_33")
    expect_regression(result, "pwm_sliding_window_iter20_no_shifts")
})

test_that("PWM regression: iterator=100 with shifts", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    set.seed(100)
    pssm <- matrix(runif(80), nrow = 20, ncol = 4)
    pssm <- pssm / rowSums(pssm)
    colnames(pssm) <- c("A", "C", "G", "T")

    gvtrack.create("pwm_reg100", "seq",
        func = "pwm",
        params = list(pssm = pssm, bidirect = TRUE, prior = 0.01)
    )
    gvtrack.iterator("pwm_reg100", sshift = -200, eshift = 200)

    result <- gextract("pwm_reg100", gintervals(1, 1000000, 1100000), iterator = 100)
    expect_regression(result, "pwm_sliding_window_iter100_with_shifts")
})

test_that("PWM regression: MAX_LIKELIHOOD mode with shifts", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    set.seed(200)
    pssm <- matrix(runif(32), nrow = 8, ncol = 4)
    pssm <- pssm / rowSums(pssm)
    colnames(pssm) <- c("A", "C", "G", "T")

    gvtrack.create("pwm_max_reg", "seq",
        func = "pwm.max",
        params = list(pssm = pssm, bidirect = TRUE)
    )
    gvtrack.iterator("pwm_max_reg", sshift = -100, eshift = 100)

    result <- gextract("pwm_max_reg", gintervals(1, 1000000, 1050000), iterator = 20)

    expect_regression(result, "pwm_sliding_window_max_iter20_shifts")
})

test_that("PWM regression: MOTIF_COUNT mode with shifts", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    set.seed(300)
    pssm <- matrix(runif(24), nrow = 6, ncol = 4)
    pssm <- pssm / rowSums(pssm)
    colnames(pssm) <- c("A", "C", "G", "T")

    gvtrack.create("pwm_count_reg", "seq",
        func = "pwm.count",
        params = list(pssm = pssm, threshold = -10)
    )
    gvtrack.iterator("pwm_count_reg", sshift = -100, eshift = 100)

    result <- gextract("pwm_count_reg", gintervals(1, 1000000, 1050000), iterator = 20)

    expect_regression(result, "pwm_sliding_window_count_iter20_shifts")
})

test_that("PWM regression: multiple intervals with shifts", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    set.seed(400)
    pssm <- matrix(runif(48), nrow = 12, ncol = 4)
    pssm <- pssm / rowSums(pssm)
    colnames(pssm) <- c("A", "C", "G", "T")

    gvtrack.create("pwm_multi_reg", "seq",
        func = "pwm",
        params = list(pssm = pssm)
    )
    gvtrack.iterator("pwm_multi_reg", sshift = -100, eshift = 100)

    intervals <- do.call(rbind, lapply(0:9, function(i) {
        gintervals(1, 1000000 + i * 10000, 1000000 + (i + 1) * 10000)
    }))

    result <- gextract("pwm_multi_reg", intervals, iterator = 20)

    expect_regression(result, "pwm_sliding_window_multi_intervals_iter20_shifts")
})

test_that("PWM regression: MAX_POS mode iterator=1 with shifts", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    set.seed(500)
    pssm <- matrix(runif(40), nrow = 10, ncol = 4)
    pssm <- pssm / rowSums(pssm)
    colnames(pssm) <- c("A", "C", "G", "T")

    gvtrack.create("pwm_pos_reg1", "seq",
        func = "pwm.max.pos",
        params = list(pssm = pssm, bidirect = TRUE)
    )
    gvtrack.iterator("pwm_pos_reg1", sshift = -150, eshift = 150)

    result <- gextract("pwm_pos_reg1", gintervals(1, 1000000, 1050000), iterator = 1)
    expect_regression(result, "pwm_sliding_window_pos_iter1_shifts")
})

test_that("PWM regression: MAX_POS mode iterator=20 with shifts", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    set.seed(500)
    pssm <- matrix(runif(40), nrow = 10, ncol = 4)
    pssm <- pssm / rowSums(pssm)
    colnames(pssm) <- c("A", "C", "G", "T")

    gvtrack.create("pwm_pos_reg20", "seq",
        func = "pwm.max.pos",
        params = list(pssm = pssm, bidirect = TRUE)
    )
    gvtrack.iterator("pwm_pos_reg20", sshift = -150, eshift = 150)

    result <- gextract("pwm_pos_reg20", gintervals(1, 1000000, 1100000), iterator = 20)
    expect_regression(result, "pwm_sliding_window_pos_iter20_shifts")
})

test_that("PWM regression: MAX_POS mode iterator=20 without shifts", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    set.seed(500)
    pssm <- matrix(runif(40), nrow = 10, ncol = 4)
    pssm <- pssm / rowSums(pssm)
    colnames(pssm) <- c("A", "C", "G", "T")

    gvtrack.create("pwm_pos_reg_no_shift", "seq",
        func = "pwm.max.pos",
        params = list(pssm = pssm, bidirect = TRUE)
    )

    result <- gextract("pwm_pos_reg_no_shift", gintervals(1, 1000000, 1100000), iterator = 20)
    expect_regression(result, "pwm_sliding_window_pos_iter20_no_shifts")
})

test_that("PWM regression: MAX_POS mode iterator=100 with shifts", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    set.seed(500)
    pssm <- matrix(runif(40), nrow = 10, ncol = 4)
    pssm <- pssm / rowSums(pssm)
    colnames(pssm) <- c("A", "C", "G", "T")

    gvtrack.create("pwm_pos_reg100", "seq",
        func = "pwm.max.pos",
        params = list(pssm = pssm, bidirect = TRUE)
    )
    gvtrack.iterator("pwm_pos_reg100", sshift = -200, eshift = 200)

    result <- gextract("pwm_pos_reg100", gintervals(1, 1000000, 1100000), iterator = 100)
    expect_regression(result, "pwm_sliding_window_pos_iter100_shifts")
})

test_that("pwm (TOTAL_LIKELIHOOD): plus-strand sliding equals spatial (no-sliding) baseline", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    pssm <- create_test_pssm()

    n <- 40
    starts <- 2000 + 0:(n - 1)
    ends <- starts + 60L
    ivs <- gintervals(rep(1L, n), starts, ends)

    gvtrack.create("pwm_plus_slide", NULL, "pwm",
        pssm = pssm, bidirect = FALSE, strand = 1,
        extend = TRUE, prior = 0.01, score.thresh = -10
    )

    gvtrack.create("pwm_plus_spatial_ref", NULL, "pwm",
        pssm = pssm, bidirect = FALSE, strand = 1,
        extend = TRUE, prior = 0.01, score.thresh = -10,
        spat_factor = rep(1.0, 5), spat_bin = 20L
    )

    res <- gextract(c("pwm_plus_slide", "pwm_plus_spatial_ref"),
        ivs,
        iterator = ivs
    )

    expect_equal(res$pwm_plus_slide, res$pwm_plus_spatial_ref, tolerance = 1e-6)
})

test_that("pwm (TOTAL_LIKELIHOOD): minus-strand sliding equals spatial (no-sliding) baseline", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    pssm <- create_test_pssm()

    # Build a series of overlapping intervals shifted by 1bp to trigger sliding
    n <- 40
    starts <- 2000 + 0:(n - 1)
    ends <- starts + 60L
    ivs <- gintervals(rep(1L, n), starts, ends)

    # Non-spatial -> sliding path is enabled
    gvtrack.create("pwm_minus_slide", NULL, "pwm",
        pssm = pssm, bidirect = FALSE, strand = -1,
        extend = TRUE, prior = 0.01, score.thresh = -10
    )

    # Spatial with weights=1 -> no sliding path, but numerically identical
    gvtrack.create("pwm_minus_spatial_ref", NULL, "pwm",
        pssm = pssm, bidirect = FALSE, strand = -1,
        extend = TRUE, prior = 0.01, score.thresh = -10,
        spat_factor = rep(1.0, 5), spat_bin = 20L
    )

    res <- gextract(c("pwm_minus_slide", "pwm_minus_spatial_ref"),
        ivs,
        iterator = ivs
    )

    expect_equal(res$pwm_minus_slide, res$pwm_minus_spatial_ref, tolerance = 1e-6)
})

test_that("pwm.max (MAX_LIKELIHOOD): plus-strand sliding equals spatial (no-sliding) baseline", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    pssm <- create_test_pssm()

    n <- 40
    starts <- 2000 + 0:(n - 1)
    ends <- starts + 60L
    ivs <- gintervals(rep(1L, n), starts, ends)

    gvtrack.create("pwmmax_plus_slide", NULL, "pwm.max",
        pssm = pssm, bidirect = FALSE, strand = 1,
        extend = TRUE, prior = 0.01, score.thresh = -10
    )

    gvtrack.create("pwmmax_plus_spatial_ref", NULL, "pwm.max",
        pssm = pssm, bidirect = FALSE, strand = 1,
        extend = TRUE, prior = 0.01, score.thresh = -10,
        spat_factor = rep(1.0, 5), spat_bin = 20L
    )

    res <- gextract(c("pwmmax_plus_slide", "pwmmax_plus_spatial_ref"),
        ivs,
        iterator = ivs
    )

    expect_equal(res$pwmmax_plus_slide, res$pwmmax_plus_spatial_ref, tolerance = 1e-6)
})

test_that("pwm.max (MAX_LIKELIHOOD): minus-strand sliding equals spatial (no-sliding) baseline", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    pssm <- create_test_pssm()

    n <- 40
    starts <- 2100 + 0:(n - 1)
    ends <- starts + 60L
    ivs <- gintervals(rep(1L, n), starts, ends)

    # Non-spatial -> sliding path
    gvtrack.create("pwmmax_minus_slide", NULL, "pwm.max",
        pssm = pssm, bidirect = FALSE, strand = -1,
        extend = TRUE, prior = 0.01, score.thresh = -10
    )

    # Spatial weights=1 -> no sliding; should match
    gvtrack.create("pwmmax_minus_spatial_ref", NULL, "pwm.max",
        pssm = pssm, bidirect = FALSE, strand = -1,
        extend = TRUE, prior = 0.01, score.thresh = -10,
        spat_factor = rep(1.0, 5), spat_bin = 20L
    )

    res <- gextract(c("pwmmax_minus_slide", "pwmmax_minus_spatial_ref"),
        ivs,
        iterator = ivs
    )

    expect_equal(res$pwmmax_minus_slide, res$pwmmax_minus_spatial_ref, tolerance = 1e-8)
})

test_that("pwm.count (MOTIF_COUNT): plus-strand sliding equals spatial (no-sliding) baseline", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    pssm <- create_test_pssm()

    n <- 40
    starts <- 2000 + 0:(n - 1)
    ends <- starts + 60L
    ivs <- gintervals(rep(1L, n), starts, ends)

    gvtrack.create("count_plus_slide", NULL, "pwm.count",
        pssm = pssm, bidirect = FALSE, strand = 1,
        extend = TRUE, prior = 0.01, score.thresh = -10
    )

    gvtrack.create("count_plus_spatial_ref", NULL, "pwm.count",
        pssm = pssm, bidirect = FALSE, strand = 1,
        extend = TRUE, prior = 0.01, score.thresh = -10,
        spat_factor = rep(1.0, 5), spat_bin = 20L
    )

    res <- gextract(c("count_plus_slide", "count_plus_spatial_ref"),
        ivs,
        iterator = ivs
    )

    expect_equal(res$count_plus_slide, res$count_plus_spatial_ref, tolerance = 1e-8)
})

test_that("pwm.count (MOTIF_COUNT): minus-strand sliding equals spatial (no-sliding) baseline", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    pssm <- create_test_pssm()

    n <- 40
    starts <- 2200 + 0:(n - 1)
    ends <- starts + 60L
    ivs <- gintervals(rep(1L, n), starts, ends)

    # Non-spatial -> sliding path
    gvtrack.create("count_minus_slide", NULL, "pwm.count",
        pssm = pssm, bidirect = FALSE, strand = -1,
        extend = TRUE, prior = 0.01, score.thresh = -10
    )

    # Spatial weights=1 -> no sliding; should match
    gvtrack.create("count_minus_spatial_ref", NULL, "pwm.count",
        pssm = pssm, bidirect = FALSE, strand = -1,
        extend = TRUE, prior = 0.01, score.thresh = -10,
        spat_factor = rep(1.0, 5), spat_bin = 20L
    )

    res <- gextract(c("count_minus_slide", "count_minus_spatial_ref"),
        ivs,
        iterator = ivs
    )

    expect_equal(res$count_minus_slide, res$count_minus_spatial_ref, tolerance = 1e-8)
})

test_that("bidirect ignores strand parameter under sliding (union semantics)", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    pssm <- create_test_pssm()

    n <- 30
    starts <- 2300 + 0:(n - 1)
    ends <- starts + 50L
    ivs <- gintervals(rep(1L, n), starts, ends)

    # Non-spatial bidi tracks with different 'strand' values — should be identical
    gvtrack.create("count_bidi_s1_slide", NULL, "pwm.count",
        pssm = pssm, bidirect = TRUE, strand = 1,
        extend = TRUE, prior = 0.01, score.thresh = -10
    )
    gvtrack.create("count_bidi_sneg1_slide", NULL, "pwm.count",
        pssm = pssm, bidirect = TRUE, strand = -1,
        extend = TRUE, prior = 0.01, score.thresh = -10
    )

    out <- gextract(c("count_bidi_s1_slide", "count_bidi_sneg1_slide"),
        ivs,
        iterator = ivs
    )

    expect_equal(out$count_bidi_s1_slide, out$count_bidi_sneg1_slide, tolerance = 1e-8)
})
