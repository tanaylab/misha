gdb.init_examples()
gdir.create("test", showWarnings = FALSE)

test_that("gdir.cd works", {
    # Create a test directory and some tracks in it
    gdir.create("test", showWarnings = FALSE)


    # Create some tracks in the test directory
    # Create multiple intervals for sparse track
    intervs_sparse <- rbind(
        gintervals(1, 0, 5000),
        gintervals(1, 10000, 15000),
        gintervals(2, 0, 5000)
    )
    gtrack.create_sparse(
        "test.sparse", "Test sparse track",
        intervs_sparse, c(1, 2, 3)
    )
    withr::defer(gtrack.rm("test.sparse", force = TRUE))

    # Create multiple intervals for fixedbin track
    intervs_fixedbin <- rbind(
        gintervals(1, 0, 10000),
        gintervals(2, 0, 10000)
    )
    gtrack.create_dense(
        "test.fixedbin", "Test fixedbin track",
        intervs_fixedbin, c(0.5, 0.7), 1000, 0
    )
    withr::defer(gtrack.rm("test.fixedbin", force = TRUE))

    # Create 2D intervals for rects track
    intervs_rects <- gintervals.2d(1, 0, 5000, 2, 0, 5000)
    gtrack.2d.create(
        "test.rects", "Test rects track",
        intervs_rects, c(10)
    )
    withr::defer(gtrack.rm("test.rects", force = TRUE))

    gdir.cd("test")
    r <- gtrack.ls()
    gdir.cd("..")

    # Check that we see the tracks we created
    expected <- c("sparse", "fixedbin", "rects")
    r <- r[!grepl("tmptrack", r)]
    expect_true(all(expected %in% r))
})

test_that("gdir.cd works (2)", {
    gdir.create("test", showWarnings = FALSE)


    gdir.cd("test")
    r <- gdir.cwd()
    gdir.cd("..")
    expect_equal(r, file.path(.misha$GROOT, "tracks", "test"))
})

test_that("gdir creates and removes directory correctly", {
    gdir.create("testdir")
    r1 <- dir(gdir.cwd())
    gdir.rm("testdir", force = TRUE)
    r2 <- dir(gdir.cwd())
    expect_true("testdir" %in% r1)
    expect_false("testdir" %in% r2)
})


test_that("gtrack.exists works correctly", {
    gdb.init_examples()
    gdir.create("test", showWarnings = FALSE)
    expect_false(gtrack.exists("aaaaaaaaaaa.nnnnnnnnnnnnnn"))

    # Create a track and verify it exists
    intervs <- rbind(
        gintervals(1, 0, 1000),
        gintervals(2, 0, 1000)
    )
    gtrack.create_sparse(
        "test.rects", "Test rects track",
        intervs, c(1, 2)
    )
    withr::defer(gtrack.rm("test.rects", force = TRUE))
    expect_true(gtrack.exists("test.rects"))
})

test_that("gtrack.ls works", {
    # Create some tracks
    gdir.create("test", showWarnings = FALSE)


    intervs_sparse <- rbind(
        gintervals(1, 0, 10000),
        gintervals(2, 0, 10000)
    )
    gtrack.create_sparse(
        "test.sparse", "Test sparse track",
        intervs_sparse, c(1, 2)
    )
    withr::defer(gtrack.rm("test.sparse", force = TRUE))

    intervs_fixedbin <- rbind(
        gintervals(1, 0, 10000),
        gintervals(2, 0, 10000)
    )
    gtrack.create_dense(
        "test.fixedbin", "Test fixedbin track",
        intervs_fixedbin, c(0.5, 0.7), 1000, 0
    )
    withr::defer(gtrack.rm("test.fixedbin", force = TRUE))

    intervs_rects <- gintervals.2d(1, 0, 5000, 2, 0, 5000)
    gtrack.2d.create(
        "test.rects", "Test rects track",
        intervs_rects, c(10)
    )
    withr::defer(gtrack.rm("test.rects", force = TRUE))

    # Test listing all tracks
    observed <- gtrack.ls()
    observed <- observed[!grepl("tmptrack", observed)]
    expect_true("test.sparse" %in% observed)
    expect_true("test.fixedbin" %in% observed)
    expect_true("test.rects" %in% observed)

    # Test listing tracks in a directory
    observed <- gtrack.ls("test")
    observed <- observed[!grepl("tmptrack", observed)]
    expect_true("test.sparse" %in% observed)
    expect_true("test.fixedbin" %in% observed)
    expect_true("test.rects" %in% observed)

    # Test filtering by non-existent attribute
    expect_null(gtrack.ls(blalaattr = "bubu"))
    expect_null(gtrack.ls("wig", created.by = "import"))
})

test_that("gtrack modify and extract for fixedbin", {
    load_test_db()
    tmptrack <- paste0("test.tmptrack_", sample(1:1e9, 1))
    gtrack.rm(tmptrack, force = TRUE)
    withr::defer(gtrack.rm(tmptrack, force = TRUE))
    gtrack.create(tmptrack, "", "test.fixedbin")
    intervs <- gscreen("test.fixedbin > 0.17 | is.na(test.fixedbin)", gintervals(c(1, 7)))
    gtrack.modify(tmptrack, "test.fixedbin + test.fixedbin", intervs)
    r <- gextract(tmptrack, gintervals(c(1, 2)), colnames = "test.tmptrack")
    expect_regression(r, "gtrack.modify_and_extract_for_fixedbin")
})

test_that("gtrack modify for sparse", {
    gdb.init_examples()
    gdir.create("test", showWarnings = FALSE)
    # Create a sparse track
    intervs_sparse <- rbind(
        gintervals(1, 0, 10000),
        gintervals(2, 0, 10000)
    )
    gtrack.create_sparse(
        "test.sparse", "Test sparse track",
        intervs_sparse, c(1, 2)
    )
    withr::defer(gtrack.rm("test.sparse", force = TRUE))

    # Create a fixedbin track for the expression
    intervs_fixedbin <- rbind(
        gintervals(1, 0, 10000),
        gintervals(2, 0, 10000)
    )
    gtrack.create_dense(
        "test.fixedbin", "Test fixedbin track",
        intervs_fixedbin, c(0.5, 0.7), 1000, 0
    )
    withr::defer(gtrack.rm("test.fixedbin", force = TRUE))

    tmptrack <- paste0("test.tmptrack_", sample(1:1e9, 1))
    gtrack.rm(tmptrack, force = TRUE)
    withr::defer(gtrack.rm(tmptrack, force = TRUE))
    gtrack.create(tmptrack, "", "test.sparse")
    expect_error(gtrack.modify(tmptrack, "test.fixedbin + test.fixedbin", gintervals(1, 1000, 2000)))
})

test_that("gtrack modify for rects", {
    gdb.init_examples()
    gdir.create("test", showWarnings = FALSE)
    # Create a rects track
    intervs_rects <- gintervals.2d(1, 0, 10000, 2, 0, 10000)
    gtrack.2d.create(
        "test.rects", "Test rects track",
        intervs_rects, c(10)
    )
    withr::defer(gtrack.rm("test.rects", force = TRUE))

    # Create a fixedbin track for the expression
    intervs_fixedbin <- rbind(
        gintervals(1, 0, 10000),
        gintervals(2, 0, 10000)
    )
    gtrack.create_dense(
        "test.fixedbin", "Test fixedbin track",
        intervs_fixedbin, c(0.5, 0.7), 1000, 0
    )
    withr::defer(gtrack.rm("test.fixedbin", force = TRUE))

    tmptrack <- paste0("test.tmptrack_", sample(1:1e9, 1))
    gtrack.rm(tmptrack, force = TRUE)
    withr::defer(gtrack.rm(tmptrack, force = TRUE))
    gtrack.create(tmptrack, "", "test.rects")
    expect_error(gtrack.modify(tmptrack, "test.fixedbin + test.fixedbin", gintervals.2d(1, 1000, 2000, 2, 3000, 4000)))
})

test_that("gtrack creation and removal for sparse", {
    gdb.init_examples()
    gdir.create("test", showWarnings = FALSE)
    tmptrack <- paste0("test.tmptrack_", sample(1:1e9, 1))
    gtrack.rm(tmptrack, force = TRUE)
    withr::defer(gtrack.rm(tmptrack, force = TRUE))
    intervs <- rbind(
        gintervals(1, 100, 2000),
        gintervals(2, 100, 2000)
    )
    gtrack.create_sparse(tmptrack, "", intervs, c(100, 200))
    r1 <- gtrack.ls()
    gtrack.rm(tmptrack, force = TRUE)
    r2 <- gtrack.ls()
    expect_true(tmptrack %in% r1)
    expect_false(tmptrack %in% r2)
})
