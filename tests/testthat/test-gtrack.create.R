create_isolated_test_db()

test_that("test for aaaaaaaaaaaaa.bbbbbbbbbbb", {
    withr::defer(gtrack.rm("aaaaaaaaaaaaa.bbbbbbbbbbb", force = TRUE))
    expect_error(gtrack.create("aaaaaaaaaaaaa.bbbbbbbbbbb", "", "test.fixedbin"))
    gdir.create("aaaaaaaaaaaaa")
    withr::defer(gdir.rm("aaaaaaaaaaaaa", recursive = TRUE, force = TRUE))
    gtrack.create("aaaaaaaaaaaaa.bbbbbbbbbbb", "", "test.fixedbin")
    r <- gextract("aaaaaaaaaaaaa.bbbbbbbbbbb", gintervals(c(1, 2), 0, 1000000))
    expect_regression(r, "gtrack.create.1")
})

test_that("test for test.tmptrack fixedbin+1", {
    tmptrack <- paste0("test.tmptrack_", sample(1:1e9, 1))
    gtrack.rm(tmptrack, force = TRUE)
    withr::defer(gtrack.rm(tmptrack, force = TRUE))
    gtrack.create(tmptrack, "", "test.fixedbin+1")
    r <- gextract(tmptrack, gintervals(c(1, 2), 0, 1000000), colnames = "test.tmptrack")
    expect_regression(r, "gtrack.create.2")
})

test_that("test for test.tmptrack with sparse iterator", {
    tmptrack <- paste0("test.tmptrack_", sample(1:1e9, 1))
    gtrack.rm(tmptrack, force = TRUE)
    withr::defer(gtrack.rm(tmptrack, force = TRUE))
    gtrack.create(tmptrack, "", "test.fixedbin+1", iterator = "test.sparse")
    r <- gextract(tmptrack, gintervals(c(1, 2), 0, 1000000), colnames = "test.tmptrack")
    expect_regression(r, "gtrack.create.3")
})

test_that("test for test.tmptrack with array iterator", {
    tmptrack <- paste0("test.tmptrack_", sample(1:1e9, 1))
    gtrack.rm(tmptrack, force = TRUE)
    withr::defer(gtrack.rm(tmptrack, force = TRUE))
    gtrack.create(tmptrack, "", "test.fixedbin+1", iterator = "test.array")
    r <- gextract(tmptrack, gintervals(c(1, 2), 0, 1000000), colnames = "test.tmptrack")
    expect_regression(r, "gtrack.create.4")
})

test_that("test for test.tmptrack rects+10", {
    tmptrack <- paste0("test.tmptrack_", sample(1:1e9, 1))
    gtrack.rm(tmptrack, force = TRUE)
    withr::defer(gtrack.rm(tmptrack, force = TRUE))
    gtrack.create(tmptrack, "", "test.rects+10")
    r <- gextract(tmptrack, gintervals.2d(chroms1 = c(2, 3), chroms2 = c(2, 4)), colnames = "test.tmptrack")
    expect_regression(r, "gtrack.create.5")
})

test_that("test for test.tmptrack with sparse intervals", {
    tmptrack <- paste0("test.tmptrack_", sample(1:1e9, 1))
    gtrack.rm(tmptrack, force = TRUE)
    intervs <- giterator.intervals("test.sparse", gintervals(c(1, 3, 4)))
    withr::defer(gtrack.rm(tmptrack, force = TRUE))
    gtrack.create(tmptrack, "", "test.fixedbin+1", iterator = intervs)
    r <- gextract(tmptrack, gintervals(c(1, 2), 0, 1000000), colnames = "test.tmptrack")
    expect_regression(r, "gtrack.create.6")
})

test_that("test for test.tmptrack with array intervals", {
    tmptrack <- paste0("test.tmptrack_", sample(1:1e9, 1))
    gtrack.rm(tmptrack, force = TRUE)
    intervs <- giterator.intervals("test.array", gintervals(c(1, 3, 4)))
    withr::defer(gtrack.rm(tmptrack, force = TRUE))
    gtrack.create(tmptrack, "", "test.fixedbin+1", iterator = intervs)
    r <- gextract(tmptrack, gintervals(c(1, 2), 0, 1000000), colnames = "test.tmptrack")
    expect_regression(r, "gtrack.create.7")
})

test_that("test for test.tmptrack with rects intervals", {
    tmptrack <- paste0("test.tmptrack_", sample(1:1e9, 1))
    gtrack.rm(tmptrack, force = TRUE)
    intervs <- giterator.intervals("test.rects", gintervals.2d(chroms1 = c(2, 3, 5), chroms2 = c(2, 4, 7)))
    withr::defer(gtrack.rm(tmptrack, force = TRUE))
    gtrack.create(tmptrack, "", "test.rects+10", iterator = intervs)
    r <- gextract(tmptrack, gintervals.2d(chroms1 = c(2, 3, 3), chroms2 = c(2, 3, 4)), colnames = "test.tmptrack")
    expect_regression(r, "gtrack.create.8")
})

test_that("test for sparse test.tmptrack with fixedbin", {
    tmptrack <- paste0("test.tmptrack_", sample(1:1e9, 1))
    gtrack.rm(tmptrack, force = TRUE)
    intervs <- gscreen("test.fixedbin > 0.2", gintervals(c(1, 2)))
    withr::defer(gtrack.rm(tmptrack, force = TRUE))
    gtrack.create_sparse(tmptrack, "", intervs, 1:dim(intervs)[1])
    r <- gextract(tmptrack, gintervals(c(1, 2), 0, 1000000), colnames = "test.tmptrack")
    expect_regression(r, "gtrack.create.9")
})

test_that("test for sparse test.tmptrack 2d", {
    tmptrack <- paste0("test.tmptrack_", sample(1:1e9, 1))
    expect_error(gtrack.create_sparse(tmptrack, "", gintervals.2d(c(1, 2)), 1:2))
})

test_that("test for 2d test.tmptrackv", {
    tmptrack <- paste0("test.tmptrack_", sample(1:1e9, 1))
    gtrack.rm(tmptrack, force = TRUE)
    intervs <- gscreen("test.rects > 80", gintervals.2d(c(1, 2)))
    withr::defer(gtrack.rm(tmptrack, force = TRUE))
    gtrack.2d.create(tmptrack, "", intervs, 1:dim(intervs)[1])
    r <- gextract(tmptrack, gintervals.2d(c(1, 2, 3)), colnames = "test.tmptrack")
    expect_regression(r, "gtrack.create.11")
})

test_that("test for 2d test.tmptrack", {
    tmptrack <- paste0("test.tmptrack_", sample(1:1e9, 1))
    intervs <- gintervals(c(1, 2))
    expect_error(gtrack.2d.create(tmptrack, "", intervs, 1:dim(intervs)[1]))
})
