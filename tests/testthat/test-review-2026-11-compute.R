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
