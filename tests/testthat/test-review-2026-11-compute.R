# Regression tests for the compute and virtual-track findings of the 2026-11 review of 5.14.4.

create_isolated_test_db()

test_that("gcor is accurate at a large offset, single process and multitask", {
    iv <- gintervals(c(1, 2), 0, 2e5)
    x <- gextract("test.fixedbin", iv, iterator = "test.fixedbin")$test.fixedbin
    x <- x[!is.na(x)] + 1e8
    y <- 2 * x - 1e8
    for (mt in c(FALSE, TRUE)) {
        withr::local_options(list(gmultitasking = mt, gmin.scope4process = 1000, gmax.processes = 8))
        r <- gcor("test.fixedbin + 1e8", "test.fixedbin * 2 + 1e8", intervals = iv, iterator = "test.fixedbin", details = TRUE)
        expect_equal(r[1, "sd1"], sd(x), tolerance = 1e-6)
        expect_equal(r[1, "sd2"], sd(y), tolerance = 1e-6)
        expect_equal(r[1, "cov"], cov(x, y), tolerance = 1e-6)
        expect_equal(r[1, "cor"], 1, tolerance = 1e-9)
    }
})

test_that("summaries give a NaN standard deviation when a value is infinite, as sd() does", {
    expr <- "ifelse(test.fixedbin > 0.15, Inf, test.fixedbin)"
    iv <- gintervals(1, 0, 1000)
    x <- gextract(expr, iv, iterator = 100, colnames = "v")$v
    expect_true(any(is.infinite(x)))
    expect_true(is.nan(sd(x)))

    expect_true(is.nan(gsummary(expr, iv, iterator = 100)[["Std dev"]]))
    expect_true(is.nan(gintervals.summary(expr, iv, iterator = 100)[["Std dev"]]))
    expect_true(is.nan(gbins.summary("test.fixedbin", breaks = c(0, 0.5), expr = expr, intervals = iv, iterator = 100)[1, "Std dev"]))
    # a finite set still has a finite one
    expect_equal(gsummary("test.fixedbin", iv, iterator = 100)[["Std dev"]], sd(gextract("test.fixedbin", iv, iterator = 100)$test.fixedbin))
})

test_that("2D gintervals.annotate with tie_method and maxneighbors > 1 gives a clear error", {
    q <- gintervals.2d(1, 100, 200, 1, 300, 400)
    a <- gintervals.2d(1, c(500, 1000), c(600, 1100), 1, c(700, 1200), c(800, 1300))
    a$lab <- c("x", "y")
    expect_error(
        gintervals.annotate(q, a, annotation_columns = "lab", maxneighbors = 2, tie_method = "min.start"),
        "1D intervals only"
    )
    expect_equal(gintervals.annotate(q, a, annotation_columns = "lab", maxneighbors = 2)$lab, c("x", "y"))
})
