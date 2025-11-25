create_isolated_test_db()
test_that("gintervals.is.bigset checks bigintervs1d", {
    withr::with_options(list(gmax.data.size = 100), {
        expect_true(gintervals.is.bigset("bigintervs1d"))
    })
})

test_that("gintervals.is.bigset checks bigintervs2d", {
    withr::with_options(list(gmax.data.size = 100), {
        expect_true(gintervals.is.bigset("bigintervs2d"))
    })
})

test_that("gintervals.is.bigset checks test.tss", {
    withr::with_options(list(gmax.data.size = 100), {
        expect_false(gintervals.is.bigset("test.tss"))
    })
})

test_that("gintervals.rbind with saved data works", {
    temp_track_name <- paste0("test.tmptrack_", sample(1:1e9, 1))
    gintervals.rm(temp_track_name, force = TRUE)
    withr::defer(gintervals.rm(temp_track_name, force = TRUE))
    intervs1 <- gextract("test.fixedbin", gintervals(c(1, 2), 1000, 4000))
    intervs2 <- gextract("test.fixedbin", gintervals(c(2, "X"), 2000, 5000))
    gintervals.save(temp_track_name, intervs2)
    r <- gintervals.rbind(intervs1, temp_track_name)
    expect_true(!is.null(r))
    expect_regression(r, "gintervals.rbind.1")
})

test_that("gintervals.ls updates after save and remove operations", {
    temp_track_name <- paste0("test.tmptrack_", sample(1:1e9, 1))
    gintervals.rm(temp_track_name, force = TRUE)
    withr::defer(gintervals.rm(temp_track_name, force = TRUE))
    gintervals.save(temp_track_name, gintervals(c(1, 2)))
    r1 <- as.character(gintervals.ls())
    r2 <- as.character(gintervals.ls())
    # Filter out temporary track names and legacy test interval sets for regression comparison
    r1_filtered <- r1[!grepl("^test\\.tmptrack_", r1) & r1 != "test.testintervs" & r1 != "testintervs"]
    r2_filtered <- r2[!grepl("^test\\.tmptrack_", r2) & r2 != "test.testintervs" & r2 != "testintervs"]
    # Sort to ensure consistent ordering for comparison
    r1_filtered <- sort(r1_filtered)
    r2_filtered <- sort(r2_filtered)
    expect_regression(list(r1_filtered, r2_filtered), "gintervals.rbind.2")
})

test_that("gextract with removed interval gives an error", {
    temp_track_name <- paste0("test.tmptrack_", sample(1:1e9, 1))
    gintervals.rm(temp_track_name, force = TRUE)
    withr::defer(gintervals.rm(temp_track_name, force = TRUE))
    gintervals.save(temp_track_name, gintervals(c(1, 2), 1000, 2000))
    gintervals.rm(temp_track_name, force = TRUE)
    expect_error(gextract("test.fixedbin", temp_track_name))
})

test_that("gintervals.rm handles non-existent data without error when using force", {
    expect_silent(gintervals.rm("test.aaaaaaaaaaaaaaaaaaa", force = TRUE))
    expect_error(gintervals.rm("test.aaaaaaaaaaaaaaaaaaa"))
})

test_that("gintervals.ls reflects changes after save", {
    temp_track_name <- paste0("test.tmptrack_", sample(1:1e9, 1))
    gintervals.rm(temp_track_name, force = TRUE)
    withr::defer(gintervals.rm(temp_track_name, force = TRUE))
    r1 <- as.character(gintervals.ls())
    gintervals.save(temp_track_name, gintervals(c(1, 2), 1000, 2000))
    r2 <- as.character(gintervals.ls())
    # Filter out temporary track names and legacy test interval sets for regression comparison
    r1_filtered <- r1[!grepl("^test\\.tmptrack_", r1) & r1 != "test.testintervs" & r1 != "testintervs"]
    r2_filtered <- r2[!grepl("^test\\.tmptrack_", r2) & r2 != "test.testintervs" & r2 != "testintervs"]
    # Sort to ensure consistent ordering for comparison
    r1_filtered <- sort(r1_filtered)
    r2_filtered <- sort(r2_filtered)
    expect_regression(list(r1_filtered, r2_filtered), "gintervals.ls.1")
})

test_that("gextract works with saved intervals", {
    temp_track_name <- paste0("test.tmptrack_", sample(1:1e9, 1))
    gintervals.rm(temp_track_name, force = TRUE)
    withr::defer(gintervals.rm(temp_track_name, force = TRUE))
    gintervals.save(temp_track_name, gintervals(c(1, 2), 1000, 2000))
    r <- gextract("test.fixedbin", temp_track_name)
    expect_true(!is.null(r))
    expect_regression(r, "gintervals.save.1")
})

test_that("gscreen and gintervals.union works correctly", {
    temp_track_name <- paste0("test.tmptrack_", sample(1:1e9, 1))
    gintervals.rm(temp_track_name, force = TRUE)
    withr::defer(gintervals.rm(temp_track_name, force = TRUE))

    intervs1 <- gscreen("test.fixedbin > 0.1 & test.fixedbin < 0.3", gintervals(c(1, 2, 4, 8, 9), 0, -1))
    intervs2 <- gscreen("test.fixedbin < 0.2", gintervals(c(1, 2, 4, 7, 9), 0, -1))

    withr::with_options(list(gmax.data.size = 1000000), {
        gintervals.union(intervs1, intervs2, intervals.set.out = temp_track_name)
    })

    r <- gintervals.load(temp_track_name)
    expect_true(!is.null(r))
    expect_regression(r, "gscreen_and_gintervals.union.1")
})

test_that("gintervals.update with saved data using chrom 1", {
    temp_track_name <- paste0("test.tmptrack_", sample(1:1e9, 1))
    gintervals.rm(temp_track_name, force = TRUE)
    withr::defer(gintervals.rm(temp_track_name, force = TRUE))

    chrom_size_limit <- sum(gintervals.chrom_sizes("bigintervs1d")$size) + 100

    withr::with_options(list(gmax.data.size = chrom_size_limit), {
        gintervals.save(temp_track_name, "bigintervs1d")
    })

    expect_error(gintervals.update(temp_track_name, gintervals(1), chrom = 1))
})

test_that("gintervals.update with loaded data using chrom 1", {
    temp_track_name <- paste0("test.tmptrack_", sample(1:1e9, 1))
    gintervals.rm(temp_track_name, force = TRUE)
    withr::defer(gintervals.rm(temp_track_name, force = TRUE))

    chrom_size_limit <- sum(gintervals.chrom_sizes("bigintervs1d")$size) + 100

    withr::with_options(list(gmax.data.size = chrom_size_limit), {
        gintervals.save(temp_track_name, "bigintervs1d")
    })

    r <- gintervals.load(temp_track_name, chrom = 2)
    expect_error(gintervals.update(temp_track_name, r[c(2, 3), ], chrom1 = 1))
})

test_that("gintervals.update with loaded data using chrom 2", {
    temp_track_name <- paste0("test.tmptrack_", sample(1:1e9, 1))
    gintervals.rm(temp_track_name, force = TRUE)
    withr::defer(gintervals.rm(temp_track_name, force = TRUE))

    chrom_size_limit <- sum(gintervals.chrom_sizes("bigintervs1d")$size) + 100

    withr::with_options(list(gmax.data.size = chrom_size_limit), {
        gintervals.save(temp_track_name, "bigintervs1d")
    })

    r <- gintervals.load(temp_track_name, chrom = 2)
    gintervals.update(temp_track_name, r[c(2, 3), ], chrom = 2)
    result <- list(gintervals.load(temp_track_name, chrom = 2), gintervals.chrom_sizes(temp_track_name))
    expect_regression(result, "gintervals.update.3")
})

test_that("gintervals.update removes chrom 2 from saved data", {
    temp_track_name <- paste0("test.tmptrack_", sample(1:1e9, 1))
    gintervals.rm(temp_track_name, force = TRUE)
    withr::defer(gintervals.rm(temp_track_name, force = TRUE))

    chrom_size_limit <- sum(gintervals.chrom_sizes("bigintervs1d")$size) + 100

    withr::with_options(list(gmax.data.size = chrom_size_limit), {
        gintervals.save(temp_track_name, "bigintervs1d")
    })

    gintervals.update(temp_track_name, NULL, chrom = 2)
    result <- list(gintervals.load(temp_track_name, chrom = 2), gintervals.chrom_sizes(temp_track_name))
    expect_regression(result, "gintervals.update.4")
})

test_that("gintervals.update with saved 2d data using chrom1 1", {
    temp_track_name <- paste0("test.tmptrack_", sample(1:1e9, 1))
    gintervals.rm(temp_track_name, force = TRUE)
    withr::defer(gintervals.rm(temp_track_name, force = TRUE))

    chrom_size_limit <- sum(gintervals.chrom_sizes("bigintervs2d")$size) + 100

    withr::with_options(list(gmax.data.size = chrom_size_limit), {
        gintervals.save(temp_track_name, "bigintervs2d")
    })

    expect_error(gintervals.update(temp_track_name, gintervals.2d(1), chrom1 = 1))
})

test_that("gintervals.update with loaded 2d data using chrom 1", {
    temp_track_name <- paste0("test.tmptrack_", sample(1:1e9, 1))
    gintervals.rm(temp_track_name, force = TRUE)
    withr::defer(gintervals.rm(temp_track_name, force = TRUE))

    chrom_size_limit <- sum(gintervals.chrom_sizes("bigintervs2d")$size) + 100

    withr::with_options(list(gmax.data.size = chrom_size_limit), {
        gintervals.save(temp_track_name, "bigintervs2d")
    })

    r <- gintervals.load(temp_track_name, chrom1 = 2, chrom2 = 2)
    expect_error(gintervals.update(temp_track_name, r[c(2, 3), ], chrom = 1))
})

test_that("gintervals.update with loaded 2d data using chrom1 2 and chrom2 2", {
    temp_track_name <- paste0("test.tmptrack_", sample(1:1e9, 1))
    gintervals.rm(temp_track_name, force = TRUE)
    withr::defer(gintervals.rm(temp_track_name, force = TRUE))

    chrom_size_limit <- sum(gintervals.chrom_sizes("bigintervs2d")$size) + 100

    withr::with_options(list(gmax.data.size = chrom_size_limit), {
        gintervals.save(temp_track_name, "bigintervs2d")
    })

    r <- gintervals.load(temp_track_name, chrom1 = 2, chrom2 = 2)
    gintervals.update(temp_track_name, r[c(2, 3), ], chrom1 = 2, chrom2 = 2)
    result <- list(gintervals.load(temp_track_name, chrom1 = 2, chrom2 = 2), gintervals.chrom_sizes(temp_track_name))
    expect_regression(result, "gintervals.update.2d.3")
})

test_that("gintervals.update removes chrom1 2 and chrom2 2 from saved 2d data", {
    temp_track_name <- paste0("test.tmptrack_", sample(1:1e9, 1))
    gintervals.rm(temp_track_name, force = TRUE)
    withr::defer(gintervals.rm(temp_track_name, force = TRUE))

    chrom_size_limit <- sum(gintervals.chrom_sizes("bigintervs2d")$size) + 100

    withr::with_options(list(gmax.data.size = chrom_size_limit), {
        gintervals.save(temp_track_name, "bigintervs2d")
    })

    gintervals.update(temp_track_name, NULL, chrom1 = 2, chrom2 = 2)
    result <- list(gintervals.load(temp_track_name, chrom1 = 2, chrom2 = 2), gintervals.chrom_sizes(temp_track_name))
    expect_regression(result, "gintervals.update.2d.4")
})
