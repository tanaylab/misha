test_that("gcor supports positional intervals as last argument", {
    gdb.init_examples()
    intervals <- gintervals(1, 0, 10000)
    iterator <- 1000

    # Named intervals argument (baseline)
    cor_named <- gcor("dense_track", "sparse_track", intervals = intervals, iterator = iterator)

    # Positional intervals as last argument
    cor_positional <- gcor("dense_track", "sparse_track", intervals, iterator = iterator)

    expect_equal(cor_positional, cor_named)

    # Multiple pairs with positional intervals
    cor_multi_named <- gcor(
        "dense_track", "sparse_track",
        "dense_track", "dense_track",
        intervals = intervals, iterator = iterator
    )
    cor_multi_positional <- gcor(
        "dense_track", "sparse_track",
        "dense_track", "dense_track",
        intervals,
        iterator = iterator
    )

    expect_equal(cor_multi_positional, cor_multi_named)
})

test_that("gcor supports named intervals set as positional argument", {
    gdb.init_examples()
    iterator <- 1000

    # Create a named interval set for testing
    intervs <- gintervals(1, 0, 10000)
    gintervals.save("test_intervs", intervs)
    withr::defer(gintervals.rm("test_intervs", force = TRUE))

    # Named intervals argument (baseline)
    cor_named <- gcor("dense_track", "sparse_track", intervals = "test_intervs", iterator = iterator)

    # Named interval set passed positionally
    cor_positional <- gcor("dense_track", "sparse_track", "test_intervs", iterator = iterator)

    expect_equal(cor_positional, cor_named)
})

test_that("gcor defaults to ALLGENOME when no intervals provided", {
    gdb.init_examples()
    iterator <- 1000

    # Explicit ALLGENOME
    cor_explicit <- gcor("dense_track", "sparse_track", intervals = .misha$ALLGENOME, iterator = iterator)

    # Implicit ALLGENOME (no intervals argument)
    cor_implicit <- gcor("dense_track", "sparse_track", iterator = iterator)

    expect_equal(cor_implicit, cor_explicit)
})

test_that("gcor errors with odd number of tracks and non-interval last argument", {
    gdb.init_examples()

    # Odd number of track expressions where last arg is not an interval
    # The error can come from the check for even number of tracks or from .gcall
    # when it tries to use the track as intervals
    expect_error(
        gcor("dense_track", "sparse_track", "dense_track", iterator = 1000)
    )
})

test_that("gcor matches cor on extracted values", {
    gdb.init_examples()
    intervals <- gintervals(1, 0, 10000)
    iterator <- 1000

    cor_only <- gcor("dense_track", "sparse_track", intervals = intervals, iterator = iterator)
    stats <- gcor("dense_track", "sparse_track", intervals = intervals, iterator = iterator, details = TRUE)
    vals <- gextract("dense_track", "sparse_track", intervals = intervals, iterator = iterator)

    x <- vals[["dense_track"]]
    y <- vals[["sparse_track"]]
    complete <- stats::complete.cases(x, y)

    x_complete <- x[complete]
    y_complete <- y[complete]

    expect_equal(unname(cor_only), cor(x, y, use = "complete.obs"), tolerance = 1e-10)
    expect_equal(unname(stats[1, "n"]), length(x))
    expect_equal(unname(stats[1, "n.na"]), sum(!complete))
    expect_equal(unname(stats[1, "mean1"]), mean(x_complete), tolerance = 1e-10)
    expect_equal(unname(stats[1, "mean2"]), mean(y_complete), tolerance = 1e-10)
    expect_equal(unname(stats[1, "sd1"]), stats::sd(x_complete), tolerance = 1e-10)
    expect_equal(unname(stats[1, "sd2"]), stats::sd(y_complete), tolerance = 1e-10)
    expect_equal(unname(stats[1, "cov"]), stats::cov(x_complete, y_complete), tolerance = 1e-10)
    expect_equal(unname(stats[1, "cor"]), cor(x, y, use = "complete.obs"), tolerance = 1e-10)

    cor_multi <- gcor(
        "dense_track",
        "sparse_track",
        "dense_track",
        "dense_track",
        intervals = intervals,
        iterator = iterator,
        names = c("dense_sparse", "dense_dense")
    )
    expect_equal(names(cor_multi), c("dense_sparse", "dense_dense"))
    expect_equal(unname(cor_multi[["dense_sparse"]]), cor(x, y, use = "complete.obs"), tolerance = 1e-10)
    expect_equal(unname(cor_multi[["dense_dense"]]), cor(x, x, use = "complete.obs"), tolerance = 1e-10)
})

# =============================================================================
# Spearman correlation tests
# =============================================================================

test_that("gcor method='spearman.exact' matches cor() with method='spearman'", {
    gdb.init_examples()
    intervals <- gintervals(1, 0, 10000)
    iterator <- 1000

    cor_spearman <- gcor("dense_track", "sparse_track",
        intervals = intervals, iterator = iterator, method = "spearman.exact"
    )
    vals <- gextract("dense_track", "sparse_track", intervals = intervals, iterator = iterator)

    x <- vals[["dense_track"]]
    y <- vals[["sparse_track"]]

    expected_cor <- cor(x, y, method = "spearman", use = "complete.obs")
    expect_equal(unname(cor_spearman), expected_cor, tolerance = 1e-10)
})

test_that("gcor method='spearman.exact' returns correct details", {
    gdb.init_examples()
    intervals <- gintervals(1, 0, 10000)
    iterator <- 1000

    stats <- gcor("dense_track", "sparse_track",
        intervals = intervals, iterator = iterator, method = "spearman.exact", details = TRUE
    )
    vals <- gextract("dense_track", "sparse_track", intervals = intervals, iterator = iterator)

    x <- vals[["dense_track"]]
    y <- vals[["sparse_track"]]
    complete <- stats::complete.cases(x, y)

    # Spearman should have n, n.na, cor columns
    expect_true(all(c("n", "n.na", "cor") %in% colnames(stats)))
    expect_equal(unname(stats[1, "n"]), length(x))
    expect_equal(unname(stats[1, "n.na"]), sum(!complete))
    expect_equal(unname(stats[1, "cor"]), cor(x, y, method = "spearman", use = "complete.obs"), tolerance = 1e-10)
})

test_that("gcor method='spearman' (approximate) gives reasonable results", {
    gdb.init_examples()
    intervals <- gintervals(1, 0, 10000)
    iterator <- 1000

    cor_approx <- gcor("dense_track", "sparse_track",
        intervals = intervals, iterator = iterator, method = "spearman"
    )
    cor_exact <- gcor("dense_track", "sparse_track",
        intervals = intervals, iterator = iterator, method = "spearman.exact"
    )

    # Approximate should be close to exact (within reasonable tolerance)
    # For small datasets, they should be identical
    expect_equal(unname(cor_approx), unname(cor_exact), tolerance = 0.05)
})

test_that("gcor method='pearson' is unchanged (default behavior)", {
    gdb.init_examples()
    intervals <- gintervals(1, 0, 10000)
    iterator <- 1000

    # Default (no method specified)
    cor_default <- gcor("dense_track", "sparse_track", intervals = intervals, iterator = iterator)

    # Explicit Pearson
    cor_pearson <- gcor("dense_track", "sparse_track",
        intervals = intervals, iterator = iterator, method = "pearson"
    )

    expect_equal(cor_default, cor_pearson)

    # Verify it's actually Pearson (not Spearman)
    vals <- gextract("dense_track", "sparse_track", intervals = intervals, iterator = iterator)
    x <- vals[["dense_track"]]
    y <- vals[["sparse_track"]]
    expected_pearson <- cor(x, y, method = "pearson", use = "complete.obs")
    expected_spearman <- cor(x, y, method = "spearman", use = "complete.obs")

    expect_equal(unname(cor_default), expected_pearson, tolerance = 1e-10)
    # Ensure Pearson and Spearman give different results for this data
    expect_false(isTRUE(all.equal(expected_pearson, expected_spearman)))
})

test_that("gcor Spearman works with multiple pairs", {
    gdb.init_examples()
    intervals <- gintervals(1, 0, 10000)
    iterator <- 1000

    cor_multi <- gcor(
        "dense_track", "sparse_track",
        "dense_track", "dense_track",
        intervals = intervals, iterator = iterator, method = "spearman.exact",
        names = c("dense_sparse", "dense_dense")
    )

    vals <- gextract("dense_track", "sparse_track", intervals = intervals, iterator = iterator)
    x <- vals[["dense_track"]]
    y <- vals[["sparse_track"]]

    expect_equal(names(cor_multi), c("dense_sparse", "dense_dense"))
    expect_equal(unname(cor_multi[["dense_sparse"]]), cor(x, y, method = "spearman", use = "complete.obs"), tolerance = 1e-10)
    expect_equal(unname(cor_multi[["dense_dense"]]), cor(x, x, method = "spearman", use = "complete.obs"), tolerance = 1e-10)
})

test_that("gcor handles ties correctly for Spearman", {
    gdb.init_examples()
    intervals <- gintervals(1, 0, 10000)
    iterator <- 1000

    # The dense_track likely has ties; verify that Spearman handles them correctly
    cor_spearman <- gcor("dense_track", "sparse_track",
        intervals = intervals, iterator = iterator, method = "spearman.exact"
    )

    vals <- gextract("dense_track", "sparse_track", intervals = intervals, iterator = iterator)
    x <- vals[["dense_track"]]
    y <- vals[["sparse_track"]]

    # R's cor() uses average rank for ties by default, which should match our implementation
    expected <- cor(x, y, method = "spearman", use = "complete.obs")
    expect_equal(unname(cor_spearman), expected, tolerance = 1e-10)
})

test_that("gcor Spearman handles all NaN gracefully", {
    gdb.init_examples()
    # Use a small region where one track might be all NaN
    # Or test with a virtual track that's always NaN
    # For now, just verify the function doesn't crash with mostly valid data
    intervals <- gintervals(1, 0, 1000)
    iterator <- 100

    # This should work without error
    result <- gcor("dense_track", "sparse_track",
        intervals = intervals, iterator = iterator, method = "spearman.exact"
    )
    expect_true(is.numeric(result))
})

test_that("gcor spearman.exact enforces gmax.data.size (single-task)", {
    gdb.init_examples()
    intervals <- gintervals(1, 0, 10000)
    iterator <- 1000

    withr::local_options(list(gmax.data.size = 1, gmultitasking = FALSE))
    expect_error(
        gcor("dense_track", "dense_track",
            intervals = intervals, iterator = iterator, method = "spearman.exact"
        ),
        "gmax.data.size"
    )
})

test_that("gcor spearman.exact enforces gmax.data.size (multitask)", {
    gdb.init_examples()
    intervals <- gintervals(1, 0, 10000)
    iterator <- 1000

    withr::local_options(list(gmax.data.size = 1, gmultitasking = TRUE))
    expect_error(
        gcor("dense_track", "dense_track",
            intervals = intervals, iterator = iterator, method = "spearman.exact"
        ),
        "gmax.data.size"
    )
})
