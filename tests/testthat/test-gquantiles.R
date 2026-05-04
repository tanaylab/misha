create_isolated_test_db()

test_that("gquantiles with test.fixedbin", {
    intervs <- gscreen("test.fixedbin > 0.2", gintervals(c(1, 2), 0, -1))
    result <- gquantiles("test.fixedbin+0.2", percentile = c(0.5, 0.3, 0.2, 0.9), intervs)
    expect_regression(result, "gquantiles_fixedbin_result")
})

test_that("gquantiles with test.rects", {
    result <- gquantiles("test.rects", percentile = c(0.5, 0.3, 0.2, 0.9, 0.999), gintervals.2d(chroms1 = c(2, 3), chroms2 = c(2, 4)))
    expect_regression(result, "gquantiles_rects_result")
})

test_that("gquantiles with test.computed2d", {
    result <- gquantiles("test.computed2d", percentile = c(0.5, 0.3, 0.2, 0.9, 0.999), gintervals.2d(chroms1 = c(6, 5), chroms2 = c(8, 9)))
    expect_regression(result, "gquantiles_computed2d_result")
})

test_that("gquantiles with test.fixedbin without intervals", {
    intervs <- gscreen("test.fixedbin > 0.2", gintervals(c(1, 2), 0, -1))
    result <- gquantiles("test.fixedbin+0.2", percentile = c(0.5, 0.999))
    expect_regression(result, "gquantiles_fixedbin_no_intervals_result")
})

# Regression: with many percentiles (e.g. .gtrack.prepare.pvals requests
# ~12k), the single-process nth_element-on-suffix path was O(k * N), which
# could take hours on a dense full-genome track. The fix flips to a single
# O(N log N) sort once `targets.size()` exceeds ~2*log2(N), and must
# remain bit-identical to the multitask k-way merge result.
test_that("gquantiles single-process matches multitask for many percentiles", {
    pcts <- sort(unique(c(
        seq(0, 1, length.out = 200),
        seq(0.001, 0.999, length.out = 800)
    )))

    multitask <- .ggetOption("gmultitasking")
    on.exit(options(gmultitasking = multitask), add = TRUE)

    options(gmultitasking = FALSE)
    serial <- gquantiles("test.fixedbin",
        percentile = pcts,
        intervals = gintervals(c(1, 2), 0, -1)
    )

    options(gmultitasking = TRUE)
    parallel <- gquantiles("test.fixedbin",
        percentile = pcts,
        intervals = gintervals(c(1, 2), 0, -1)
    )

    expect_equal(unname(serial), unname(parallel))
})
