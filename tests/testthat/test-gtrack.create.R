load_test_db()
test_that("test for aaaaaaaaaaaaa.bbbbbbbbbbb", {
    withr::defer(gtrack.rm("aaaaaaaaaaaaa.bbbbbbbbbbb", force = TRUE))
    expect_error(gtrack.create("aaaaaaaaaaaaa.bbbbbbbbbbb", "", "test.fixedbin"))
    gdir.create("aaaaaaaaaaaaa")
    withr::defer(gdir.rm("aaaaaaaaaaaaa", recursive = TRUE, force = TRUE))
    gtrack.create("aaaaaaaaaaaaa.bbbbbbbbbbb", "", "test.fixedbin")
    r <- gextract("aaaaaaaaaaaaa.bbbbbbbbbbb", gintervals(c(1, 2), 0, 1000000))
    expect_regression(r, "gtrack.create.1")
})

test_that("test for temp track fixedbin+1", {
    track_name <- random_track_name()
    gtrack.rm(track_name, force = TRUE)
    withr::defer(gtrack.rm(track_name, force = TRUE))
    gtrack.create(track_name, "", "test.fixedbin+1")
    r <- gextract(track_name, gintervals(c(1, 2), 0, 1000000), colnames = "test.tmptrack")
    expect_regression(r, "gtrack.create.2")
})

test_that("test for temp track with sparse iterator", {
    track_name <- random_track_name()
    gtrack.rm(track_name, force = TRUE)
    withr::defer(gtrack.rm(track_name, force = TRUE))
    gtrack.create(track_name, "", "test.fixedbin+1", iterator = "test.sparse")
    r <- gextract(track_name, gintervals(c(1, 2), 0, 1000000), colnames = "test.tmptrack")
    expect_regression(r, "gtrack.create.3")
})

test_that("test for temp track with array iterator", {
    track_name <- random_track_name()
    gtrack.rm(track_name, force = TRUE)
    withr::defer(gtrack.rm(track_name, force = TRUE))
    gtrack.create(track_name, "", "test.fixedbin+1", iterator = "test.array")
    r <- gextract(track_name, gintervals(c(1, 2), 0, 1000000), colnames = "test.tmptrack")
    expect_regression(r, "gtrack.create.4")
})

test_that("test for temp track rects+10", {
    track_name <- random_track_name()
    gtrack.rm(track_name, force = TRUE)
    withr::defer(gtrack.rm(track_name, force = TRUE))
    gtrack.create(track_name, "", "test.rects+10")
    r <- gextract(track_name, gintervals.2d(chroms1 = c(2, 3), chroms2 = c(2, 4)), colnames = "test.tmptrack")
    expect_regression(r, "gtrack.create.5")
})

test_that("test for temp track with sparse intervals", {
    track_name <- random_track_name()
    gtrack.rm(track_name, force = TRUE)
    intervs <- giterator.intervals("test.sparse", gintervals(c(1, 3, 4)))
    withr::defer(gtrack.rm(track_name, force = TRUE))
    gtrack.create(track_name, "", "test.fixedbin+1", iterator = intervs)
    r <- gextract(track_name, gintervals(c(1, 2), 0, 1000000), colnames = "test.tmptrack")
    expect_regression(r, "gtrack.create.6")
})

test_that("test for temp track with array intervals", {
    track_name <- random_track_name()
    gtrack.rm(track_name, force = TRUE)
    intervs <- giterator.intervals("test.array", gintervals(c(1, 3, 4)))
    withr::defer(gtrack.rm(track_name, force = TRUE))
    gtrack.create(track_name, "", "test.fixedbin+1", iterator = intervs)
    r <- gextract(track_name, gintervals(c(1, 2), 0, 1000000), colnames = "test.tmptrack")
    expect_regression(r, "gtrack.create.7")
})

test_that("test for temp track with rects intervals", {
    track_name <- random_track_name()
    gtrack.rm(track_name, force = TRUE)
    intervs <- giterator.intervals("test.rects", gintervals.2d(chroms1 = c(2, 3, 5), chroms2 = c(2, 4, 7)))
    withr::defer(gtrack.rm(track_name, force = TRUE))
    gtrack.create(track_name, "", "test.rects+10", iterator = intervs)
    r <- gextract(track_name, gintervals.2d(chroms1 = c(2, 3, 3), chroms2 = c(2, 3, 4)), colnames = "test.tmptrack")
    expect_regression(r, "gtrack.create.8")
})

test_that("test for sparse temp track with fixedbin", {
    track_name <- random_track_name()
    gtrack.rm(track_name, force = TRUE)
    intervs <- gscreen("test.fixedbin > 0.2", gintervals(c(1, 2)))
    withr::defer(gtrack.rm(track_name, force = TRUE))
    gtrack.create_sparse(track_name, "", intervs, 1:dim(intervs)[1])
    r <- gextract(track_name, gintervals(c(1, 2), 0, 1000000), colnames = "test.tmptrack")
    expect_regression(r, "gtrack.create.9")
})

test_that("test for sparse temp track 2d", {
    track_name <- random_track_name()
    expect_error(gtrack.create_sparse(track_name, "", gintervals.2d(c(1, 2)), 1:2))
})

test_that("test for 2d temp trackv", {
    track_name <- random_track_name()
    gtrack.rm(track_name, force = TRUE)
    intervs <- gscreen("test.rects > 80", gintervals.2d(c(1, 2)))
    withr::defer(gtrack.rm(track_name, force = TRUE))
    gtrack.2d.create(track_name, "", intervs, 1:dim(intervs)[1])
    r <- gextract(track_name, gintervals.2d(c(1, 2, 3)), colnames = "test.tmptrack")
    expect_regression(r, "gtrack.create.11")
})

test_that("test for 2d temp track", {
    track_name <- random_track_name()
    intervs <- gintervals(c(1, 2))
    expect_error(gtrack.2d.create(track_name, "", intervs, 1:dim(intervs)[1]))
})
