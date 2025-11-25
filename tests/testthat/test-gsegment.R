create_isolated_test_db()

test_that("gsegment with test.fixedbin", {
    expect_regression(gsegment("test.fixedbin", 10000, maxpval = 0.000001), "gsegment_fixedbin")
})

test_that("gsegment with test.sparse", {
    expect_error(gsegment("test.sparse", 10000, maxpval = 0.000001))
})

test_that("gsegment with test.array", {
    expect_error(gsegment("test.array", 10000, maxpval = 0.000001))
})

test_that("gsegment with test.rects", {
    expect_error(gsegment("test.rects", 10000, maxpval = 0.000001))
})

test_that("gsegment with test.computed2d", {
    expect_error(gsegment("test.computed2d", 10000, maxpval = 0.000001))
})

test_that("gsegment with modified test.fixedbin", {
    intervs <- gscreen("test.fixedbin > 0.2", gintervals(c(1, 2), 0, -1))
    expect_regression(gsegment("test.fixedbin*2", 10000, maxpval = 0.000001, intervals = intervs), "gsegment_fixedbin_mod")
})

test_that("gsegment with data size option and sampling for test.sparse", {
    temp_track_name <- paste0("test.tmptrack_", sample(1:1e9, 1))
    gintervals.rm(temp_track_name, force = TRUE)
    withr::defer(gintervals.rm(temp_track_name, force = TRUE))

    intervs <- gscreen("test.fixedbin > 0.14", gintervals(c(1, 2), 0, -1))
    set.seed(60427)
    intervs <- intervs[sample(nrow(intervs)), ]

    withr::with_options(c(gmax.data.size = 3200), {
        gsegment("test.sparse", 10000, maxpval = 0.0001, iterator = 50, intervals.set.out = temp_track_name)
    })
    r <- gintervals.load(temp_track_name)
    expect_regression(r, "gsegment_sparse_data_size")
})
