create_isolated_test_db()

test_that("gtrack.smooth with test.fixedbin using LINEAR_RAMP", {
    tmptrack <- paste0("test.tmptrack_", sample(1:1e9, 1))
    gtrack.rm(tmptrack, force = TRUE)
    withr::defer(gtrack.rm(tmptrack, force = TRUE))
    gtrack.smooth(tmptrack, "", "test.fixedbin", 10000, alg = "LINEAR_RAMP")
    r <- gextract(tmptrack, gintervals(c(1, 2), 0, 1000000), colnames = "test.tmptrack")
    expect_regression(r, "gtrack.smooth.fixedbin_LINEAR_RAMP")
})

test_that("gtrack.smooth with test.fixedbin using MEAN", {
    tmptrack <- paste0("test.tmptrack_", sample(1:1e9, 1))
    gtrack.rm(tmptrack, force = TRUE)
    withr::defer(gtrack.rm(tmptrack, force = TRUE))
    gtrack.smooth(tmptrack, "", "test.fixedbin", 10000, alg = "MEAN")
    r <- gextract(tmptrack, gintervals(c(1, 2), 0, 1000000), colnames = "test.tmptrack")
    expect_regression(r, "gtrack.smooth.fixedbin_MEAN")
})

test_that("gtrack.smooth with test.sparse using LINEAR_RAMP", {
    tmptrack <- paste0("test.tmptrack_", sample(1:1e9, 1))
    gtrack.rm(tmptrack, force = TRUE)
    withr::defer(gtrack.rm(tmptrack, force = TRUE))
    gtrack.smooth(tmptrack, "", "test.sparse", 10000, alg = "LINEAR_RAMP", iterator = 1000)
    r <- gextract(tmptrack, gintervals(c(1, 2), 0, 1000000), colnames = "test.tmptrack")
    expect_regression(r, "gtrack.smooth.sparse_LINEAR_RAMP")
})

test_that("gtrack.smooth with test.array using LINEAR_RAMP", {
    tmptrack <- paste0("test.tmptrack_", sample(1:1e9, 1))
    gtrack.rm(tmptrack, force = TRUE)
    withr::defer(gtrack.rm(tmptrack, force = TRUE))
    gtrack.smooth(tmptrack, "", "test.array", 10000, alg = "LINEAR_RAMP", iterator = 1000)
    r <- gextract(tmptrack, gintervals(c(1, 2), 0, 1000000), colnames = "test.tmptrack")
    expect_regression(r, "gtrack.smooth.array_LINEAR_RAMP")
})

test_that("gtrack.smooth with test.rects using LINEAR_RAMP", {
    tmptrack <- paste0("test.tmptrack_", sample(1:1e9, 1))
    expect_error(gtrack.smooth(tmptrack, "", "test.rects", 10000, alg = "LINEAR_RAMP"))
})
