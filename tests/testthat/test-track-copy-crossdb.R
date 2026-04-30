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

test_that("split_indexed_to_per_chrom handles multiple non-empty chroms", {
    withr::with_tempdir({
        create_test_db("db_multi")
        gsetroot("db_multi")
        # Write data on BOTH chr1 and chr2 so the splitter has to produce two output files.
        intervs <- rbind(
            gintervals(1, 0, 1000),
            gintervals(2, 0, 1000)
        )
        gtrack.create_sparse("tm", "test", intervs, c(5, 5))
        track_dir <- file.path(normalizePath("db_multi"), "tracks", "tm.track")
        files_before <- list.files(track_dir)
        bytes_before <- setNames(
            lapply(file.path(track_dir, files_before), readBin, what = "raw", n = 1e8),
            files_before
        )

        gtrack.convert_to_indexed("tm")
        chrom_names <- misha:::.gdb.chrom_names_at(normalizePath("db_multi"))
        misha:::.gtrack.split_indexed_to_per_chrom(track_dir, chrom_names, remove_indexed = TRUE)

        # Both chrom files reappeared, byte-identical
        files_after <- list.files(track_dir)
        expect_setequal(files_after, files_before)
        bytes_after <- setNames(
            lapply(file.path(track_dir, files_after), readBin, what = "raw", n = 1e8),
            files_after
        )
        expect_equal(bytes_after[files_before], bytes_before[files_before])
        # Sanity: track values still extract correctly
        expect_equal(gextract("tm", gintervals(1, 0, 500))$tm[1], 5)
        expect_equal(gextract("tm", gintervals(2, 0, 500))$tm[1], 5)
    })
})

test_that("split_indexed_to_per_chrom with remove_indexed=FALSE keeps track.dat/idx", {
    withr::with_tempdir({
        create_test_db("db_keep")
        gsetroot("db_keep")
        gtrack.create_sparse("tk", "test", gintervals(1, 0, 1000), 3)
        track_dir <- file.path(normalizePath("db_keep"), "tracks", "tk.track")
        gtrack.convert_to_indexed("tk")
        expect_true(file.exists(file.path(track_dir, "track.idx")))

        chrom_names <- misha:::.gdb.chrom_names_at(normalizePath("db_keep"))
        misha:::.gtrack.split_indexed_to_per_chrom(track_dir, chrom_names, remove_indexed = FALSE)

        expect_true(file.exists(file.path(track_dir, "track.idx")))
        expect_true(file.exists(file.path(track_dir, "track.dat")))
        # Per-chrom files also produced
        expect_true("chr1" %in% list.files(track_dir))
    })
})

test_that("split_indexed_to_per_chrom errors on chromid out of range and preserves indexed pair", {
    withr::with_tempdir({
        create_test_db("db_oor")
        gsetroot("db_oor")
        # Create a track with data on BOTH chr1 and chr2 so that the index references
        # chrom_id 0 AND chrom_id 1. Passing only one chrom name will trigger the guard.
        intervs <- rbind(gintervals(1, 0, 1000), gintervals(2, 0, 1000))
        gtrack.create_sparse("t_oor", "test", intervs, c(1, 1))
        gtrack.convert_to_indexed("t_oor")
        track_dir <- file.path(normalizePath("db_oor"), "tracks", "t_oor.track")

        # Sanity: indexed pair exists
        expect_true(file.exists(file.path(track_dir, "track.idx")))
        expect_true(file.exists(file.path(track_dir, "track.dat")))

        # Pass only one chrom name -- must fail with a clear message
        expect_error(
            misha:::.gtrack.split_indexed_to_per_chrom(track_dir, "only_one_chrom",
                remove_indexed = TRUE
            ),
            "chrom_id"
        )

        # Indexed pair must still be intact (the splitter must NOT delete on error)
        expect_true(file.exists(file.path(track_dir, "track.idx")))
        expect_true(file.exists(file.path(track_dir, "track.dat")))

        # No leftover .tmp files
        expect_length(list.files(track_dir, pattern = "\\.tmp$"), 0)
    })
})

test_that("gtrack.copy with db= lands the track in the named dataset", {
    withr::with_tempdir({
        create_test_db("workdb")
        create_test_db("otherdb")
        gsetroot("workdb")
        gdataset.load(normalizePath("otherdb"))
        gtrack.create_sparse("src_t", "src", gintervals(1, 0, 1000), 9)

        gtrack.copy("src_t", "copied_t", db = normalizePath("otherdb"))

        expect_true(gtrack.exists("copied_t"))
        expect_equal(gtrack.dataset("copied_t"), normalizePath("otherdb"))
        expect_equal(gtrack.dataset("src_t"), normalizePath("workdb"))
        expect_equal(gextract("copied_t", gintervals(1, 0, 500))$copied_t[1], 9)
    })
})

test_that("gtrack.copy: per-chrom src to indexed dest converts on the fly", {
    withr::with_tempdir({
        create_test_db("perchrom_src")
        create_test_db("indexed_dest")
        gdb.init("indexed_dest")
        gdb.convert_to_indexed(force = TRUE, verbose = FALSE)
        gsetroot("perchrom_src")
        gdataset.load(normalizePath("indexed_dest"))
        gtrack.create_sparse("t", "src", gintervals(1, 0, 1000), 11)

        gtrack.copy("t", "t_copy", db = normalizePath("indexed_dest"))

        # Destination should be in indexed format
        dest_dir <- file.path(normalizePath("indexed_dest"), "tracks", "t_copy.track")
        expect_true(file.exists(file.path(dest_dir, "track.idx")))
        expect_true(file.exists(file.path(dest_dir, "track.dat")))
        # Per-chrom files should be gone
        expect_length(list.files(dest_dir, pattern = "^chr"), 0)
        # Values intact
        expect_equal(gextract("t_copy", gintervals(1, 0, 500))$t_copy[1], 11)
    })
})

test_that("gtrack.copy: indexed src to per-chrom dest splits on the fly", {
    withr::with_tempdir({
        create_test_db("indexed_src")
        create_test_db("perchrom_dest")
        gdb.init("indexed_src")
        gdb.convert_to_indexed(force = TRUE, verbose = FALSE)
        gtrack.create_sparse("t", "src", gintervals(1, 0, 1000), 22)
        gsetroot("perchrom_dest")
        gdataset.load(normalizePath("indexed_src"))

        gtrack.copy("t", "t_copy", db = normalizePath("perchrom_dest"))

        dest_dir <- file.path(normalizePath("perchrom_dest"), "tracks", "t_copy.track")
        expect_false(file.exists(file.path(dest_dir, "track.idx")))
        # Per-chrom files present (chr1 should exist; chr2 may not since data is only on chr1)
        expect_true(any(c("chr1", "chr2") %in% list.files(dest_dir)))
        expect_equal(gextract("t_copy", gintervals(1, 0, 500))$t_copy[1], 22)
    })
})

test_that("gtrack.copy: per-chrom dense src to indexed dest converts on the fly", {
    withr::with_tempdir({
        create_test_db("perchrom_src_dense")
        create_test_db("indexed_dest_dense")
        gdb.init("indexed_dest_dense")
        gdb.convert_to_indexed(force = TRUE, verbose = FALSE)
        gsetroot("perchrom_src_dense")
        gdataset.load(normalizePath("indexed_dest_dense"))

        # Dense track via gtrack.create with a fixed-bin iterator
        intervs <- gintervals(1, 0, 1000)
        gtrack.create("d", "dense", "1", iterator = 100)

        gtrack.copy("d", "d_copy", db = normalizePath("indexed_dest_dense"))

        dest_dir <- file.path(normalizePath("indexed_dest_dense"), "tracks", "d_copy.track")
        expect_true(file.exists(file.path(dest_dir, "track.idx")))
        expect_true(file.exists(file.path(dest_dir, "track.dat")))
        # Values intact for dense track
        result <- gextract("d_copy", gintervals(1, 0, 500), iterator = 100)
        expect_true(all(result$d_copy == 1))
    })
})

test_that("gtrack.copy drops chromosomes not present in destination, with a warning", {
    withr::with_tempdir({
        # src has chr1, chr2, chr3; dest has only chr1, chr2
        create_test_db("src3", chrom_sizes = data.frame(
            chrom = c("chr1", "chr2", "chr3"), size = c(10000, 10000, 10000)
        ))
        create_test_db("dest2", chrom_sizes = data.frame(
            chrom = c("chr1", "chr2"), size = c(10000, 10000)
        ))
        gsetroot("src3")
        gtrack.create_sparse(
            "t", "src",
            rbind(gintervals(1, 0, 100), gintervals(3, 0, 100)),
            c(5, 5)
        )

        expect_warning(
            gtrack.copy("t", "t_copy", db = normalizePath("dest2")),
            "chr3"
        )
        gsetroot("dest2")
        # Track exists; values for chr1 are 5, chr3 is gone (doesn't exist in dest)
        expect_true(gtrack.exists("t_copy"))
        expect_equal(gextract("t_copy", gintervals(1, 0, 100))$t_copy[1], 5)
    })
})

test_that("gtrack.copy: indexed -> indexed with different chrom order remaps via two-stage pipeline", {
    withr::with_tempdir({
        # src order: chr1, chr2; dest order: chr2, chr1
        create_test_db("src_idx", chrom_sizes = data.frame(
            chrom = c("chr1", "chr2"), size = c(10000, 10000)
        ))
        create_test_db("dest_idx", chrom_sizes = data.frame(
            chrom = c("chr2", "chr1"), size = c(10000, 10000)
        ))
        gdb.init("src_idx")
        gdb.convert_to_indexed(force = TRUE, verbose = FALSE)
        gtrack.create_sparse(
            "t", "src",
            rbind(gintervals(1, 0, 100), gintervals(2, 0, 100)),
            c(13, 13)
        )
        gdb.init("dest_idx")
        gdb.convert_to_indexed(force = TRUE, verbose = FALSE)
        gsetroot("src_idx")

        gtrack.copy("t", "t_copy", db = normalizePath("dest_idx"))

        gsetroot("dest_idx")
        expect_true(gtrack.exists("t_copy"))
        expect_equal(gextract("t_copy", gintervals(1, 0, 100))$t_copy[1], 13)
        expect_equal(gextract("t_copy", gintervals(2, 0, 100))$t_copy[1], 13)
    })
})

test_that("gtrack.copy: src 'chr1' -> dest '1' handles prefix variant via rename", {
    withr::with_tempdir({
        create_test_db("src_chrprefix", chrom_sizes = data.frame(
            chrom = c("chr1", "chr2"), size = c(10000, 10000)
        ))
        create_test_db("dest_noprefix", chrom_sizes = data.frame(
            chrom = c("1", "2"), size = c(10000, 10000)
        ))
        gsetroot("src_chrprefix")
        gtrack.create_sparse("t", "src", gintervals(1, 0, 100), 7)

        gtrack.copy("t", "t_copy", db = normalizePath("dest_noprefix"))

        gsetroot("dest_noprefix")
        expect_true(gtrack.exists("t_copy"))
        # Chrom is "1" (no prefix) in dest
        expect_equal(gextract("t_copy", gintervals("1", 0, 100))$t_copy[1], 7)
    })
})
