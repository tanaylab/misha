create_isolated_test_db()

test_that("lookup and extract with default binning", {
    tmptrack <- paste0("test.tmptrack_", sample(1:1e9, 1))
    m1 <- matrix(1:15, nrow = 5, ncol = 3)
    gtrack.rm(tmptrack, force = TRUE)
    withr::defer(gtrack.rm(tmptrack, force = TRUE))
    gtrack.lookup(tmptrack, "", m1, "test.fixedbin", seq(0.1, 0.2, length.out = 6), "test.sparse", seq(0.25, 0.48, length.out = 4), iterator = "test.fixedbin")
    r <- gextract(tmptrack, gintervals(c(1, 3)), colnames = "test.tmptrack")
    expect_regression(r, "gtrack.lookup.default_binning")
})

test_that("lookup and extract without force binning", {
    tmptrack <- paste0("test.tmptrack_", sample(1:1e9, 1))
    m1 <- matrix(1:15, nrow = 5, ncol = 3)
    gtrack.rm(tmptrack, force = TRUE)
    withr::defer(gtrack.rm(tmptrack, force = TRUE))
    gtrack.lookup(tmptrack, "", m1, "test.fixedbin", seq(0.1, 0.2, length.out = 6), "test.sparse", seq(0.25, 0.48, length.out = 4), force.binning = FALSE, iterator = "test.fixedbin")
    r <- gextract(tmptrack, gintervals(c(1, 3)), colnames = "test.tmptrack")
    expect_regression(r, "gtrack.lookup.no_force_binning")
})

test_that("lookup with rects and extract 2D intervals", {
    tmptrack <- paste0("test.tmptrack_", sample(1:1e9, 1))
    m1 <- matrix(1:15, nrow = 5, ncol = 3)
    gtrack.rm(tmptrack, force = TRUE)
    withr::defer(gtrack.rm(tmptrack, force = TRUE))
    gtrack.lookup(tmptrack, "", m1, "test.rects", seq(50, 100, length.out = 6), "test.rects / 2", seq(0, 40, length.out = 4), force.binning = FALSE)
    r <- gextract(tmptrack, gintervals.2d(chroms1 = c(3, 5), chroms2 = c(4, 6)), colnames = "test.tmptrack")
    expect_regression(r, "gtrack.lookup.2D_intervals")
})
