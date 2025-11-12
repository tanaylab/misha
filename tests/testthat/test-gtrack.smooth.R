load_test_db()
test_that("gtrack.smooth with test.fixedbin using LINEAR_RAMP", {
    track_name <- random_track_name()
    gtrack.rm(track_name, force = TRUE)
    withr::defer(gtrack.rm(track_name, force = TRUE))
    gtrack.smooth(track_name, "", "test.fixedbin", 10000, alg = "LINEAR_RAMP")
    r <- gextract(track_name, gintervals(c(1, 2), 0, 1000000), colnames = "test.tmptrack")
    expect_regression(r, "gtrack.smooth.fixedbin_LINEAR_RAMP")
})

test_that("gtrack.smooth with test.fixedbin using MEAN", {
    track_name <- random_track_name()
    gtrack.rm(track_name, force = TRUE)
    withr::defer(gtrack.rm(track_name, force = TRUE))
    gtrack.smooth(track_name, "", "test.fixedbin", 10000, alg = "MEAN")
    r <- gextract(track_name, gintervals(c(1, 2), 0, 1000000), colnames = "test.tmptrack")
    expect_regression(r, "gtrack.smooth.fixedbin_MEAN")
})

test_that("gtrack.smooth with test.sparse using LINEAR_RAMP", {
    track_name <- random_track_name()
    gtrack.rm(track_name, force = TRUE)
    withr::defer(gtrack.rm(track_name, force = TRUE))
    gtrack.smooth(track_name, "", "test.sparse", 10000, alg = "LINEAR_RAMP", iterator = 1000)
    r <- gextract(track_name, gintervals(c(1, 2), 0, 1000000), colnames = "test.tmptrack")
    expect_regression(r, "gtrack.smooth.sparse_LINEAR_RAMP")
})

test_that("gtrack.smooth with test.array using LINEAR_RAMP", {
    track_name <- random_track_name()
    gtrack.rm(track_name, force = TRUE)
    withr::defer(gtrack.rm(track_name, force = TRUE))
    gtrack.smooth(track_name, "", "test.array", 10000, alg = "LINEAR_RAMP", iterator = 1000)
    r <- gextract(track_name, gintervals(c(1, 2), 0, 1000000), colnames = "test.tmptrack")
    expect_regression(r, "gtrack.smooth.array_LINEAR_RAMP")
})

test_that("gtrack.smooth with test.rects using LINEAR_RAMP", {
    track_name <- random_track_name()
    expect_error(gtrack.smooth(track_name, "", "test.rects", 10000, alg = "LINEAR_RAMP"))
})
