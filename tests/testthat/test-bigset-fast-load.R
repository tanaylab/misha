# Tests for fast-path loading of indexed bigsets

test_that(".gintervals.is_indexed_bigset detects indexed bigsets correctly", {
    skip_if_not_installed("withr")
    withr::local_options(list(gmulticontig.indexed_format = FALSE))
    old_gwd <- get("GWD", envir = .misha)
    old_root <- dirname(old_gwd)
    tmp_root <- withr::local_tempdir()

    # Create minimal database
    utils::write.table(
        data.frame(chrom = c("chr1", "chr2", "chr3"), size = c(1e6, 1e6, 1e6)),
        file.path(tmp_root, "chrom_sizes.txt"),
        quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE
    )
    dir.create(file.path(tmp_root, "tracks", "test"), recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(tmp_root, "seq"), recursive = TRUE, showWarnings = FALSE)
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

    withr::defer(gintervals.rm("test.indexed_bigset", force = TRUE), envir = parent.frame())

    # Create intervals spanning multiple chromosomes
    intervs <- gintervals(
        chroms = c("chr1", "chr2", "chr3"),
        starts = c(0, 0, 0),
        ends = c(1000, 2000, 3000)
    )

    # Save as bigset
    withr::with_options(list(gmax.data.size = 1), {
        misha::gintervals.save("test.indexed_bigset", intervs)
    })

    # At this point it's a per-chromosome bigset, not indexed
    expect_false(misha:::.gintervals.is_indexed_bigset("test.indexed_bigset"))

    # Convert to indexed format
    misha::gintervals.convert_to_indexed("test.indexed_bigset", remove.old = TRUE, force = TRUE)

    # Now it should be detected as indexed
    expect_true(misha:::.gintervals.is_indexed_bigset("test.indexed_bigset"))
})

test_that(".gintervals.is_indexed_bigset returns FALSE for non-bigsets", {
    skip_if_not_installed("withr")
    withr::local_options(list(gmulticontig.indexed_format = FALSE))
    old_gwd <- get("GWD", envir = .misha)
    old_root <- dirname(old_gwd)
    tmp_root <- withr::local_tempdir()

    utils::write.table(
        data.frame(chrom = c("chr1", "chr2"), size = c(1e6, 1e6)),
        file.path(tmp_root, "chrom_sizes.txt"),
        quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE
    )
    dir.create(file.path(tmp_root, "tracks", "test"), recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(tmp_root, "seq"), recursive = TRUE, showWarnings = FALSE)
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

    withr::defer(gintervals.rm("test.small_intervs", force = TRUE), envir = parent.frame())

    # Create small intervals (won't be bigset)
    intervs <- gintervals(
        chroms = c("chr1", "chr2"),
        starts = c(0, 0),
        ends = c(100, 100)
    )

    misha::gintervals.save("test.small_intervs", intervs)

    # Should return FALSE for non-bigset
    expect_false(misha:::.gintervals.is_indexed_bigset("test.small_intervs"))
})

test_that(".gintervals.is_indexed_bigset returns FALSE when update files exist", {
    skip_if_not_installed("withr")
    withr::local_options(list(gmulticontig.indexed_format = FALSE))
    old_gwd <- get("GWD", envir = .misha)
    old_root <- dirname(old_gwd)
    tmp_root <- withr::local_tempdir()

    utils::write.table(
        data.frame(chrom = c("chr1", "chr2", "chr3"), size = c(1e6, 1e6, 1e6)),
        file.path(tmp_root, "chrom_sizes.txt"),
        quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE
    )
    dir.create(file.path(tmp_root, "tracks", "test"), recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(tmp_root, "seq"), recursive = TRUE, showWarnings = FALSE)
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

    withr::defer(gintervals.rm("test.updated_bigset", force = TRUE), envir = parent.frame())

    # Create intervals spanning multiple chromosomes
    intervs <- gintervals(
        chroms = c("chr1", "chr2", "chr3"),
        starts = c(0, 0, 0),
        ends = c(1000, 2000, 3000)
    )

    # Save as bigset
    withr::with_options(list(gmax.data.size = 1), {
        misha::gintervals.save("test.updated_bigset", intervs)
    })

    # Convert to indexed format
    misha::gintervals.convert_to_indexed("test.updated_bigset", remove.old = TRUE, force = TRUE)

    # Verify it's indexed
    expect_true(misha:::.gintervals.is_indexed_bigset("test.updated_bigset"))

    # Update the intervals (this creates per-chromosome files)
    update_intervs <- gintervals(
        chroms = c("chr1"),
        starts = c(5000),
        ends = c(6000)
    )
    misha::gintervals.update("test.updated_bigset", update_intervs, chrom = "chr1")

    # Now should return FALSE because update files exist
    expect_false(misha:::.gintervals.is_indexed_bigset("test.updated_bigset"))
})

test_that("indexed bigset load produces identical results to per-chromosome load", {
    skip_if_not_installed("withr")
    withr::local_options(list(gmulticontig.indexed_format = FALSE))
    old_gwd <- get("GWD", envir = .misha)
    old_root <- dirname(old_gwd)
    tmp_root <- withr::local_tempdir()

    utils::write.table(
        data.frame(chrom = c("chr1", "chr2", "chr3"), size = c(1e6, 1e6, 1e6)),
        file.path(tmp_root, "chrom_sizes.txt"),
        quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE
    )
    dir.create(file.path(tmp_root, "tracks", "test"), recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(tmp_root, "seq"), recursive = TRUE, showWarnings = FALSE)
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

    withr::defer(gintervals.rm("test.parity_bigset", force = TRUE), envir = parent.frame())

    # Create intervals with multiple chromosomes and extra columns
    intervs <- data.frame(
        chrom = factor(c("chr1", "chr1", "chr2", "chr2", "chr3")),
        start = c(0, 5000, 0, 3000, 1000),
        end = c(1000, 6000, 2000, 4000, 2000),
        score = c(1.5, 2.5, 3.5, 4.5, 5.5),
        name = c("a", "b", "c", "d", "e")
    )

    # Save as bigset (per-chromosome format)
    withr::with_options(list(gmax.data.size = 1), {
        misha::gintervals.save("test.parity_bigset", intervs)
    })

    # Load via slow path (per-chromosome)
    loaded_slow <- withr::with_options(list(gmax.data.size = 1e9), {
        misha::gintervals.load("test.parity_bigset")
    })

    # Convert to indexed format
    misha::gintervals.convert_to_indexed("test.parity_bigset", remove.old = TRUE, force = TRUE)

    # Load via fast path (indexed)
    loaded_fast <- withr::with_options(list(gmax.data.size = 1e9), {
        misha::gintervals.load("test.parity_bigset")
    })

    # Compare results
    expect_equal(nrow(loaded_fast), nrow(loaded_slow))
    expect_equal(as.character(loaded_fast$chrom), as.character(loaded_slow$chrom))
    expect_equal(loaded_fast$start, loaded_slow$start)
    expect_equal(loaded_fast$end, loaded_slow$end)
    expect_equal(loaded_fast$score, loaded_slow$score)
    expect_equal(loaded_fast$name, loaded_slow$name)

    # Verify factor levels are correct
    expect_equal(length(unique(loaded_fast$chrom)), 3)
})

test_that("2D indexed bigset load produces identical results to per-chromosome load", {
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
    dir.create(file.path(tmp_root, "seq"), recursive = TRUE, showWarnings = FALSE)
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

    withr::defer(gintervals.rm("test.parity_bigset_2d", force = TRUE), envir = parent.frame())

    # Create 2D intervals
    intervs <- gintervals.2d(
        chroms1 = c(1, 1, 2, 3),
        starts1 = c(0, 10000, 0, 0),
        ends1 = c(5000, 15000, 5000, 5000),
        chroms2 = c(1, 2, 2, 3),
        starts2 = c(0, 0, 10000, 0),
        ends2 = c(5000, 5000, 15000, 5000)
    )

    # Save as bigset (per-chromosome-pair format)
    withr::with_options(list(gmax.data.size = 1), {
        misha::gintervals.save("test.parity_bigset_2d", intervs)
    })

    # Load via slow path
    loaded_slow <- withr::with_options(list(gmax.data.size = 1e9), {
        misha::gintervals.load("test.parity_bigset_2d")
    })

    # Convert to indexed format
    misha::gintervals.2d.convert_to_indexed("test.parity_bigset_2d", remove.old = TRUE, force = TRUE)

    # Load via fast path
    loaded_fast <- withr::with_options(list(gmax.data.size = 1e9), {
        misha::gintervals.load("test.parity_bigset_2d")
    })

    # Compare results
    expect_equal(nrow(loaded_fast), nrow(loaded_slow))
    expect_equal(as.character(loaded_fast$chrom1), as.character(loaded_slow$chrom1))
    expect_equal(as.character(loaded_fast$chrom2), as.character(loaded_slow$chrom2))
    expect_equal(loaded_fast$start1, loaded_slow$start1)
    expect_equal(loaded_fast$end1, loaded_slow$end1)
    expect_equal(loaded_fast$start2, loaded_slow$start2)
    expect_equal(loaded_fast$end2, loaded_slow$end2)
})

# Restore the test database after all tests
suppressMessages(gdb.init("/net/mraid20/export/tgdata/db/tgdb/misha_test_db/"))
