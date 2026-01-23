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
    withr::defer(gintervals.rm("test_intervs"))

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
