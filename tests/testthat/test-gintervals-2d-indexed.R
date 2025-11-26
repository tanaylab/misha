test_that("gintervals.2d.convert_to_indexed round-trips data", {
    skip_if_not_installed("withr")
    withr::local_options(list(gmulticontig.indexed_format = FALSE))
    old_gwd <- get("GWD", envir = .misha)
    old_root <- dirname(old_gwd)
    tmp_root <- withr::local_tempdir()
    utils::write.table(
        data.frame(chrom = c("1", "2", "3"), size = c(1e6, 1e6, 1e6)),
        file.path(tmp_root, "chrom_sizes.txt"),
        quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE
    )
    dir.create(file.path(tmp_root, "tracks", "test"), recursive = TRUE, showWarnings = FALSE)
    gdb.init(tmp_root)
    withr::defer(
        {
            if (dir.exists(old_root)) {
                try(gsetroot(old_root), silent = TRUE)
                try(gdb.reload(), silent = TRUE)
            }
        },
        envir = parent.frame()
    )
    gsetroot(tmp_root)
    gdb.reload()
    withr::defer(gintervals.rm("test.converted_2d", force = TRUE), envir = parent.frame())

    intervs <- gintervals.2d(
        chroms1 = c(1, 1, 2),
        starts1 = c(0, 10000, 0),
        ends1 = c(5000, 15000, 5000),
        chroms2 = c(1, 2, 2),
        starts2 = c(0, 0, 10000),
        ends2 = c(5000, 5000, 15000)
    )
    withr::with_options(list(gmax.data.size = 1), {
        misha::gintervals.save("test.converted_2d", intervs)
    })

    expect_no_error(misha::gintervals.2d.convert_to_indexed("test.converted_2d", remove.old = FALSE, force = TRUE))

    idx_path <- file.path(tmp_root, "tracks", "test", "converted_2d.interv", "intervals2d.idx")
    dat_path <- file.path(tmp_root, "tracks", "test", "converted_2d.interv", "intervals2d.dat")
    expect_true(file.exists(idx_path))
    expect_true(file.exists(dat_path))

    withr::with_options(list(gmax.data.size = 1e9), {
        loaded <- misha::gintervals.load("test.converted_2d")
        expect_equal(nrow(loaded), nrow(intervs))
        expect_equal(loaded$chrom1, intervs$chrom1)
        expect_equal(loaded$start1, intervs$start1)
        expect_equal(loaded$end1, intervs$end1)
        expect_equal(loaded$chrom2, intervs$chrom2)
        expect_equal(loaded$start2, intervs$start2)
        expect_equal(loaded$end2, intervs$end2)
    })
})

test_that("gintervals.2d.convert_to_indexed handles sparse chromosome pairs", {
    skip_if_not_installed("withr")
    withr::local_options(list(gmulticontig.indexed_format = FALSE))
    old_gwd <- get("GWD", envir = .misha)
    old_root <- dirname(old_gwd)
    tmp_root <- withr::local_tempdir()
    utils::write.table(
        data.frame(chrom = c("1", "2", "3"), size = c(1e6, 1e6, 1e6)),
        file.path(tmp_root, "chrom_sizes.txt"),
        quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE
    )
    dir.create(file.path(tmp_root, "tracks", "test"), recursive = TRUE, showWarnings = FALSE)
    gdb.init(tmp_root)
    withr::defer(
        {
            if (dir.exists(old_root)) {
                try(gsetroot(old_root), silent = TRUE)
                try(gdb.reload(), silent = TRUE)
            }
        },
        envir = parent.frame()
    )
    gsetroot(tmp_root)
    gdb.reload()
    withr::defer(gintervals.rm("test.sparse_pairs_2d", force = TRUE), envir = parent.frame())

    intervs <- gintervals.2d(
        chroms1 = c(1, 3),
        starts1 = c(0, 0),
        ends1 = c(1000, 1000),
        chroms2 = c(2, 3),
        starts2 = c(0, 0),
        ends2 = c(1000, 1000)
    )
    withr::with_options(list(gmax.data.size = 1), {
        misha::gintervals.save("test.sparse_pairs_2d", intervs)
    })

    expect_no_error(misha::gintervals.2d.convert_to_indexed("test.sparse_pairs_2d", remove.old = FALSE, force = TRUE))

    idx_path <- file.path(tmp_root, "tracks", "test", "sparse_pairs_2d.interv", "intervals2d.idx")
    dat_path <- file.path(tmp_root, "tracks", "test", "sparse_pairs_2d.interv", "intervals2d.dat")
    expect_true(file.exists(idx_path))
    expect_true(file.exists(dat_path))

    withr::with_options(list(gmax.data.size = 1e9), {
        loaded <- misha::gintervals.load("test.sparse_pairs_2d")
        expect_equal(loaded$chrom1, intervs$chrom1)
        expect_equal(loaded$chrom2, intervs$chrom2)
    })
})

test_that("indexed 2d intervals can be used as iterator in gextract", {
    create_isolated_test_db()
    gdir.create("temp", showWarnings = FALSE)
    intervs <- gextract("test.rects", gintervals.2d(chroms1 = c(2, 3), chroms2 = c(2, 4)), iterator = "test.rects")
    gintervals.save(intervals = intervs, intervals.set.out = "temp.intervs_2d")
    res <- gextract("test.rects", gintervals.2d(chroms1 = c(2, 3), chroms2 = c(2, 4)), iterator = "temp.intervs_2d")
    expect_equal(res, intervs)
})

test_that("indexed 2d intervals can be used as intervals in gextract", {
    create_isolated_test_db()
    gdir.create("temp", showWarnings = FALSE)
    intervs <- gintervals.2d(chroms1 = c(2, 3), chroms2 = c(2, 4))
    gintervals.save(intervals = intervs, intervals.set.out = "temp.intervs_2d_1")
    res1 <- gextract("test.rects", intervals = "temp.intervs_2d_1", iterator = "test.rects")
    res2 <- gextract("test.rects", intervals = intervs, iterator = "test.rects")
    expect_equal(res1, res2)
})
