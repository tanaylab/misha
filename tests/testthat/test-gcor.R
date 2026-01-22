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
