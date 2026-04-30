# tests/testthat/test-track-copy-crossdb.R
test_that(".gdb.is_indexed_at and .gdb.chrom_names_at probe a db without loading it", {
    withr::with_tempdir({
        create_test_db("perchrom_db")
        expect_false(misha:::.gdb.is_indexed_at(normalizePath("perchrom_db")))
        expect_equal(
            misha:::.gdb.chrom_names_at(normalizePath("perchrom_db")),
            c("chr1", "chr2")
        )
    })
})

test_that(".gdb.is_indexed_at returns TRUE for an indexed db", {
    withr::with_tempdir({
        create_test_db("idx_db")
        gdb.init("idx_db")
        gdb.convert_to_indexed(force = TRUE, verbose = FALSE)
        expect_true(misha:::.gdb.is_indexed_at(normalizePath("idx_db")))
    })
})

test_that("split_indexed_to_per_chrom recreates per-chrom files identical to pre-conversion", {
    withr::with_tempdir({
        create_test_db("db_a")
        gsetroot("db_a")
        gtrack.create_sparse("t1", "test", gintervals(1, 0, 1000), 7)
        # Snapshot per-chrom files before conversion
        track_dir <- file.path(normalizePath("db_a"), "tracks", "t1.track")
        before <- list.files(track_dir)
        sizes_before <- file.info(file.path(track_dir, before))$size
        names(sizes_before) <- before

        gtrack.convert_to_indexed("t1")
        expect_true(file.exists(file.path(track_dir, "track.idx")))

        # Split back
        chrom_names <- misha:::.gdb.chrom_names_at(normalizePath("db_a"))
        misha:::.gtrack.split_indexed_to_per_chrom(track_dir, chrom_names, remove_indexed = TRUE)

        expect_false(file.exists(file.path(track_dir, "track.idx")))
        expect_false(file.exists(file.path(track_dir, "track.dat")))
        # Every original per-chrom file is back, byte-for-byte (compare sizes here; full bytes covered in later test)
        after <- list.files(track_dir)
        expect_setequal(after, before)
        sizes_after <- file.info(file.path(track_dir, after))$size
        names(sizes_after) <- after
        expect_equal(sizes_after[before], sizes_before[before])
    })
})

test_that("split_indexed_to_per_chrom is byte-identical to pre-conversion", {
    withr::with_tempdir({
        create_test_db("db_b")
        gsetroot("db_b")
        intervs <- gintervals(1, 0, 5000)
        gtrack.create_sparse("t2", "test", intervs, 42)
        track_dir <- file.path(normalizePath("db_b"), "tracks", "t2.track")
        files <- list.files(track_dir, full.names = FALSE)
        before <- setNames(lapply(file.path(track_dir, files), readBin, what = "raw", n = 1e8), files)

        gtrack.convert_to_indexed("t2")
        chrom_names <- misha:::.gdb.chrom_names_at(normalizePath("db_b"))
        misha:::.gtrack.split_indexed_to_per_chrom(track_dir, chrom_names, remove_indexed = TRUE)

        after <- setNames(lapply(file.path(track_dir, files), readBin, what = "raw", n = 1e8), files)
        expect_equal(after, before)
    })
})
