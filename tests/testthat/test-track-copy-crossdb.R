# tests/testthat/test-track-copy-crossdb.R

# Every test here gsetroot()s into a tempdir database that withr removes on exit, leaving
# GROOT dangling for whichever file the parallel worker picks up next.
restore_groot_on_exit()

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

test_that(".gdb.chrom_names_at gives the chrom id order gsetroot gives a per-chromosome db", {
    local_db_state()
    withr::with_tempdir({
        db <- create_db_with_unsorted_chrom_sizes("unsorted_db")
        names_at <- misha:::.gdb.chrom_names_at(db)
        gsetroot(db)
        expect_equal(names_at, as.character(gintervals.all()$chrom))
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

test_that("gtrack.copy: vector src with prefix dest", {
    withr::with_tempdir({
        create_test_db("a")
        create_test_db("b")
        gsetroot("a")
        gtrack.create_sparse("x", "x", gintervals(1, 0, 100), 1)
        gtrack.create_sparse("y", "y", gintervals(1, 0, 100), 2)

        out <- gtrack.copy(c("x", "y"), dest = "ns", db = normalizePath("b"))

        expect_setequal(out, c("ns.x", "ns.y"))
        gsetroot("b")
        expect_true(gtrack.exists("ns.x"))
        expect_true(gtrack.exists("ns.y"))
        expect_equal(gextract("ns.x", gintervals(1, 0, 100))$ns.x[1], 1)
        expect_equal(gextract("ns.y", gintervals(1, 0, 100))$ns.y[1], 2)
    })
})

test_that("gtrack.copy: vector src with NULL dest keeps names", {
    withr::with_tempdir({
        create_test_db("a")
        create_test_db("b")
        gsetroot("a")
        gtrack.create_sparse("x", "x", gintervals(1, 0, 100), 1)
        gtrack.create_sparse("y", "y", gintervals(1, 0, 100), 2)

        gtrack.copy(c("x", "y"), db = normalizePath("b"))

        gsetroot("b")
        expect_true(gtrack.exists("x"))
        expect_true(gtrack.exists("y"))
        expect_equal(gextract("x", gintervals(1, 0, 100))$x[1], 1)
        expect_equal(gextract("y", gintervals(1, 0, 100))$y[1], 2)
    })
})

test_that("gtrack.copy: overwrite=FALSE errors on existing dest, overwrite=TRUE replaces", {
    withr::with_tempdir({
        create_test_db("a")
        create_test_db("b")
        gsetroot("a")
        gtrack.create_sparse("x", "x", gintervals(1, 0, 100), 1)

        # First copy succeeds
        gtrack.copy("x", "x_copy", db = normalizePath("b"))

        # Second copy without overwrite errors
        expect_error(
            gtrack.copy("x", "x_copy", db = normalizePath("b")),
            "already exists"
        )

        # Replace source with different value, then overwrite=TRUE
        gtrack.rm("x", force = TRUE)
        gtrack.create_sparse("x", "x", gintervals(1, 0, 100), 99)

        gtrack.copy("x", "x_copy", db = normalizePath("b"), overwrite = TRUE)
        gsetroot("b")
        expect_equal(gextract("x_copy", gintervals(1, 0, 100))$x_copy[1], 99)
    })
})

test_that("gtrack.copy: .attrs and .vars survive cross-db split+pack roundtrip", {
    withr::with_tempdir({
        create_test_db("perchrom")
        create_test_db("indexed")
        gdb.init("indexed")
        gdb.convert_to_indexed(force = TRUE, verbose = FALSE)
        gsetroot("perchrom")
        gtrack.create_sparse("t_attr", "src", gintervals(1, 0, 100), 5)
        gtrack.attr.set("t_attr", "experiment", "foo")
        gtrack.var.set("t_attr", "extra_metadata", list(x = 1, y = "hello"))

        gtrack.copy("t_attr", "t_attr_copy", db = normalizePath("indexed"))

        gsetroot("indexed")
        expect_true(gtrack.exists("t_attr_copy"))
        expect_equal(gtrack.attr.get("t_attr_copy", "experiment"), "foo")
        expect_equal(
            gtrack.var.get("t_attr_copy", "extra_metadata"),
            list(x = 1, y = "hello")
        )
    })
})

test_that("split followed by pack reproduces the original indexed pair byte-for-byte", {
    withr::with_tempdir({
        create_test_db("rt")
        gsetroot("rt")
        gtrack.create_sparse(
            "t_rt", "src",
            rbind(gintervals(1, 0, 100), gintervals(2, 0, 100)),
            c(7, 7)
        )
        track_dir <- file.path(normalizePath("rt"), "tracks", "t_rt.track")
        gtrack.convert_to_indexed("t_rt")

        # Snapshot original indexed bytes
        idx_before <- readBin(file.path(track_dir, "track.idx"), "raw", n = 1e8)
        dat_before <- readBin(file.path(track_dir, "track.dat"), "raw", n = 1e8)

        chrom_names <- misha:::.gdb.chrom_names_at(normalizePath("rt"))
        misha:::.gtrack.split_indexed_to_per_chrom(track_dir, chrom_names, remove_indexed = TRUE)
        # Determine track type from a per-chrom file (sparse here)
        misha:::.gtrack.pack_per_chrom_to_indexed(track_dir, chrom_names, "sparse")

        idx_after <- readBin(file.path(track_dir, "track.idx"), "raw", n = 1e8)
        dat_after <- readBin(file.path(track_dir, "track.dat"), "raw", n = 1e8)

        expect_equal(idx_after, idx_before)
        expect_equal(dat_after, dat_before)
    })
})

test_that("gtrack.copy splits an indexed track of a per-chromosome db by the chrom ids gsetroot gives it", {
    local_db_state()
    withr::with_tempdir({
        # chrom_sizes.txt is unsorted, so the chrom ids of src follow the sorted names
        src <- create_db_with_unsorted_chrom_sizes("src")
        dest_perchrom <- create_db_with_unsorted_chrom_sizes("dest_perchrom")
        dest_indexed <- normalizePath(create_test_db("dest_indexed", chrom_sizes = data.frame(
            chrom = c("chr2", "chr10", "chr1", "chrX", "chr1_KI270706v1_random"),
            size = c(2000, 1500, 1000, 1200, 500)
        )))
        gdb.init(dest_indexed)
        gdb.convert_to_indexed(force = TRUE, verbose = FALSE)

        gsetroot(src)
        intervs <- gintervals.all()
        intervs$end <- 100
        gtrack.create_sparse("sp", "x", intervs, seq_len(nrow(intervs)))
        gtrack.convert_to_indexed("sp")
        src_vals <- gextract("sp", gintervals.all())
        expected <- setNames(src_vals$sp, as.character(src_vals$chrom))

        for (db in c(dest_perchrom, dest_indexed)) {
            gsetroot(src)
            gtrack.copy("sp", "sp_copy", db = db)
            gsetroot(db)
            res <- gextract("sp_copy", gintervals.all())
            expect_equal(setNames(res$sp_copy, as.character(res$chrom))[names(expected)], expected, info = db)
        }
    })
})

for (.seq_state in c("missing", "empty")) {
    test_that(sprintf("a dataset of a per-chromosome db with its seq/ %s numbers chromosomes as the db does", .seq_state), {
        local_db_state()
        td <- tempfile("ds_order_")
        dir.create(td)
        withr::defer(unlink(td, recursive = TRUE))
        db <- create_db_with_unsorted_chrom_sizes(file.path(td, "db"))
        gsetroot(db)

        iv <- gintervals(c("chr1", "chr2", "chr10", "chrX"), 0, 100)
        gtrack.create_sparse("t1", "x", iv, c(1, 2, 10, 23))
        suppressMessages(gtrack.convert_to_indexed("t1"))
        gtrack.2d.create("r2", "x", data.frame(chrom1 = "chr1", start1 = 10, end1 = 15, chrom2 = "chr2", start2 = 30, end2 = 35), 7)
        expected <- gextract("t1", gintervals.all())
        expected2d <- gextract("r2", gintervals.2d.all())

        ds <- file.path(td, "ds")
        suppressMessages(gdataset.save(ds, "d", tracks = "t1"))
        # a dataset moved away from its db, copied without its seq/ link or assembled by hand
        unlink(file.path(ds, "seq"))
        if (.seq_state == "empty") {
            dir.create(file.path(ds, "seq"))
        }
        gtrack.rm("t1", force = TRUE)
        suppressMessages(gdataset.load(ds))
        expect_equal(.gdb.chrom_names_at(ds), c("chr1", "chr1_KI270706v1_random", "chr10", "chr2", "chrX"))

        # an indexed track copied out of the dataset keeps its values on their chromosomes
        gtrack.copy("t1", "t1c")
        expect_equal(gextract("t1c", gintervals.all())$t1c, expected$t1)

        # a 2D track copied into the dataset
        gtrack.copy("r2", "r2c", db = ds)
        expect_equal(gextract("r2c", gintervals.2d.all())$r2c, expected2d$r2)
    })
}

test_that("a dataset of an indexed db without its seq/ link is indexed as the db is", {
    local_db_state()
    td <- tempfile("ds_indexed_")
    dir.create(td)
    withr::defer(unlink(td, recursive = TRUE))
    db <- create_db_with_unsorted_chrom_sizes(file.path(td, "db"))
    gsetroot(db)
    # tracks in per-chromosome and per-pair files, made before the database was converted
    gtrack.create_sparse("v2", "x", gintervals(c("chr1", "chr2"), 0, 100), c(7, 8))
    gtrack.2d.create("r2", "x", data.frame(chrom1 = "chr1", start1 = 10, end1 = 15, chrom2 = "chr2", start2 = 30, end2 = 35), 7)
    suppressMessages(gdb.convert_to_indexed(groot = db, force = TRUE, validate = FALSE))
    gsetroot(db)
    expect_true(.gdb.is_indexed_at(db))

    iv <- gintervals(c("chr1", "chr2"), 0, 100)
    gtrack.create_sparse("t1", "x", iv, c(1, 2))
    ds <- file.path(td, "ds")
    suppressMessages(gdataset.save(ds, "d", tracks = "t1"))
    unlink(file.path(ds, "seq"))
    gtrack.rm("t1", force = TRUE)
    suppressMessages(gdataset.load(ds))
    expect_true(.gdb.is_indexed_at(ds))
    # nor does a genome.seq without genome.idx make a sequence of its own
    dir.create(file.path(ds, "seq"))
    file.create(file.path(ds, "seq", "genome.seq"))
    expect_true(.gdb.is_indexed_at(ds))
    unlink(file.path(ds, "seq"), recursive = TRUE)

    # an indexed track is copied in the database's format; other tracks go in as they are
    gtrack.create_sparse("v1", "x", iv, c(5, 6))
    suppressMessages(gtrack.convert_to_indexed("v1"))
    gtrack.copy(c("v1", "v2", "r2"), "c", db = ds)
    expect_equal(gtrack.info("c.v1")$format, "indexed")
    expect_false(file.exists(file.path(ds, "tracks", "c", "v2.track", "track.idx")))
    expect_false(file.exists(file.path(ds, "tracks", "c", "r2.track", "track.idx")))
    expect_equal(gextract("c.v1", iv)$c.v1, c(5, 6))
    expect_equal(gextract("c.v2", iv)$c.v2, c(7, 8))
    expect_equal(gextract("c.r2", gintervals.2d.all())$c.r2, 7)
})

test_that("an unloaded database whose seq/ holds no sequence gives no chrom order to guess from", {
    local_db_state()
    td <- tempfile("noseq_")
    dir.create(td)
    withr::defer(unlink(td, recursive = TRUE))
    db <- create_db_with_unsorted_chrom_sizes(file.path(td, "db"))
    gsetroot(db)
    iv <- gintervals(c("chr1", "chr2"), 0, 100)
    gtrack.create_sparse("sp", "x", iv, c(1, 2))
    gtrack.create_sparse("spi", "x", iv, c(3, 4))
    suppressMessages(gtrack.convert_to_indexed("spi"))

    # a directory of the same database, not loaded, with tracks/ and chrom_sizes.txt but no seq/
    other <- file.path(td, "other")
    dir.create(file.path(other, "tracks"), recursive = TRUE)
    expect_true(file.copy(file.path(db, "chrom_sizes.txt"), other))
    expect_error(.gdb.chrom_names_at(other), "cannot be told")

    # a liftover of an indexed source track out of it needs the order
    src <- file.path(other, "tracks", "lift.track")
    dir.create(src)
    expect_error(misha:::.gtrack.liftover.src_chroms(src), "cannot be told")
    unlink(src, recursive = TRUE)
    # a genome.idx without genome.seq is no sequence either
    dir.create(file.path(other, "seq"))
    writeBin(raw(16), file.path(other, "seq", "genome.idx"))
    expect_error(.gdb.chrom_names_at(other), "cannot be told")
    unlink(file.path(other, "seq"), recursive = TRUE)

    # copies into it go by chromosome name: the source's order is known, as it is loaded
    gtrack.copy("sp", "sp_copy", db = other)
    gtrack.copy("spi", "spi_copy", db = other)
    suppressMessages(gdataset.load(other))
    expect_equal(gextract("sp_copy", iv)$sp_copy, c(1, 2))
    expect_equal(gextract("spi_copy", iv)$spi_copy, c(3, 4))
})

test_that("gtrack.copy into a database whose order cannot be told copies per-pair 2D tracks by name, and checks before overwriting", {
    local_db_state()
    td <- withr::local_tempdir()
    db <- create_db_with_unsorted_chrom_sizes(file.path(td, "db"))
    gsetroot(db)
    pairs <- data.frame(chrom1 = c("chr1", "chr10"), start1 = 10, end1 = 15, chrom2 = c("chr2", "chrX"), start2 = 30, end2 = 35)
    gtrack.2d.create("r2", "x", pairs, c(12, 1023))
    gtrack.2d.create("r2i", "x", pairs, c(5, 6))
    suppressMessages(gtrack.2d.convert_to_indexed("r2i"))
    # a tracks-only copy of the database: no seq/, unprefixed names
    other <- file.path(td, "other")
    dir.create(file.path(other, "tracks"), recursive = TRUE)
    expect_true(file.copy(file.path(db, "chrom_sizes.txt"), other))

    gtrack.copy("r2", "r2c", db = other)
    old <- sort(list.files(file.path(other, "tracks", "r2c.track")))
    # an indexed 2D track needs the order: it stops before the existing track is touched
    expect_error(gtrack.copy("r2i", "r2c", db = other, overwrite = TRUE), "cannot be told")
    expect_equal(sort(list.files(file.path(other, "tracks", "r2c.track"))), old)
    suppressMessages(gdataset.load(other))
    expect_equal(gextract("r2c", gintervals.2d.all())$r2c, c(12, 1023))

    # a chr-prefixed first name and a .seq file without the prefix still show the order
    gdb.unload()
    third <- file.path(td, "third")
    dir.create(file.path(third, "seq"), recursive = TRUE)
    writeLines(c("chrA\t10", paste0(c("B", "C", "D", "E"), "\t10")), file.path(third, "chrom_sizes.txt"))
    writeBin(charToRaw(strrep("A", 10)), file.path(third, "seq", "A.seq"))
    expect_equal(.gdb.chrom_names_at(third)[1], "chrA")
})

test_that("gtrack.copy copies a per-pair 2D track between chr-prefixed and unprefixed names of the same chromosomes", {
    local_db_state()
    td <- withr::local_tempdir()
    unprefixed <- create_db_with_unsorted_chrom_sizes(file.path(td, "unprefixed"))
    prefixed <- create_db_with_unsorted_chrom_sizes(file.path(td, "prefixed"))
    writeLines(paste0("chr", readLines(file.path(prefixed, "chrom_sizes.txt"))), file.path(prefixed, "chrom_sizes.txt"))
    pairs <- data.frame(chrom1 = c("chr1", "chr10"), start1 = 10, end1 = 15, chrom2 = c("chr2", "chrX"), start2 = 30, end2 = 35)
    for (from in c(prefixed, unprefixed)) {
        to <- setdiff(c(prefixed, unprefixed), from)
        gsetroot(from)
        gtrack.2d.create("r2", "x", pairs, c(12, 1023))
        gtrack.copy("r2", "r2c", db = to)
        gtrack.rm("r2", force = TRUE)
        gsetroot(to)
        res <- gextract("r2c", gintervals.2d.all())
        expect_equal(setNames(res$r2c, paste(res$chrom1, res$chrom2))[c("chr1 chr2", "chr10 chrX")], c("chr1 chr2" = 12, "chr10 chrX" = 1023), info = basename(from))
        gtrack.rm("r2c", force = TRUE)
    }
})

test_that("gtrack.copy with overwrite = TRUE into an unloaded database leaves the session's tracks out of its cache", {
    local_db_state()
    td <- withr::local_tempdir()
    src <- create_db_with_unsorted_chrom_sizes(file.path(td, "src"))
    dest <- create_db_with_unsorted_chrom_sizes(file.path(td, "dest"))
    gsetroot(src)
    gtrack.create_sparse("sp", "x", gintervals(c("chr1", "chr2"), 0, 10), c(1, 2))
    gtrack.create_sparse("other", "x", gintervals("chr1", 0, 10), 3)
    gtrack.copy("sp", "sp_c", db = dest)
    gtrack.copy("other", "sp_c", db = dest, overwrite = TRUE)
    cache <- file.path(dest, ".db.cache")
    if (file.exists(cache)) {
        expect_false(any(c("sp", "other") %in% unlist(readRDS(cache))))
    }
    gsetroot(dest)
    expect_equal(gtrack.ls(), "sp_c")
    expect_equal(gextract("sp_c", gintervals("chr1", 0, 10))$sp_c, 3)
})

test_that("gtrack.copy with overwrite = TRUE stops for the same directory under two names", {
    local_db_state()
    td <- withr::local_tempdir()
    db <- create_db_with_unsorted_chrom_sizes(file.path(td, "db"))
    gsetroot(db)
    gdir.create("sub", showWarnings = FALSE)
    gtrack.create_sparse("sub.t", "x", gintervals("chr1", 0, 10), 1)
    # tracks/alias is a link to tracks/sub: alias.t is sub.t
    expect_true(file.symlink(file.path(db, "tracks", "sub"), file.path(db, "tracks", "alias")))
    gdb.reload()
    expect_true(all(c("alias.t", "sub.t") %in% gtrack.ls()))
    expect_error(gtrack.copy("sub.t", "alias.t", overwrite = TRUE), "Source and destination are the same track")
    expect_equal(gextract("sub.t", gintervals("chr1", 0, 10))$sub.t, 1)
})

test_that("gtrack.copy with overwrite = TRUE keeps the existing track when the copy stops", {
    local_db_state()
    td <- withr::local_tempdir()
    mk_db <- function(name, chroms, indexed = FALSE) {
        path <- file.path(td, name)
        dir.create(file.path(path, "tracks"), recursive = TRUE)
        dir.create(file.path(path, "seq"))
        if (indexed) {
            fa <- tempfile(fileext = ".fa", tmpdir = td)
            writeLines(unlist(lapply(chroms, function(ch) c(paste0(">", ch), strrep("A", 1000)))), fa)
            invisible(misha:::.gcall("gseq_multifasta_import", fa, file.path(path, "seq", "genome.seq"), file.path(path, "seq", "genome.idx"), FALSE, misha:::.misha_env()))
        } else {
            for (ch in chroms) writeBin(charToRaw(strrep("A", 1000)), file.path(path, "seq", paste0(ch, ".seq")))
        }
        writeLines(paste(chroms, 1000, sep = "\t"), file.path(path, "chrom_sizes.txt"))
        normalizePath(path)
    }
    src <- mk_db("src", c("chr1", "chr2"))
    gsetroot(src)
    gtrack.create_sparse("s1", "x", gintervals(c("chr1", "chr2"), 0, 10), c(1, 2))
    pair <- data.frame(chrom1 = "chr1", start1 = 0, end1 = 10, chrom2 = "chr2", start2 = 0, end2 = 10)
    gtrack.2d.create("r2", "x", pair, 5)
    gtrack.2d.create("r2i", "x", pair, 6)
    suppressMessages(gtrack.2d.convert_to_indexed("r2i"))

    cases <- list(
        list(track = "r2", db = mk_db("renamed", c("chrA", "chrB")), error = "requires the same chromosome names"),
        list(track = "r2i", db = mk_db("reversed", c("chr2", "chr1"), indexed = TRUE), error = "requires identical chromosome order"),
        list(track = "s1", db = mk_db("disjoint", "chrZ"), error = "no chromosomes from source database are present"),
        list(track = "r2i", db = mk_db("per_chrom", c("chr1", "chr2")), error = "indexed 2D track .* into a per-chromosome database"),
        list(track = "r2", db = mk_db("indexed", c("chr1", "chr2"), indexed = TRUE), error = "format conversion to a non-active dataset")
    )
    for (case in cases) {
        old <- file.path(case$db, "tracks", "t.track")
        dir.create(old)
        writeLines("keep me", file.path(old, "sentinel"))
        expect_error(gtrack.copy(case$track, "t", db = case$db, overwrite = TRUE), case$error)
        expect_equal(list.files(old), "sentinel", info = basename(case$db))
    }
})

test_that("an unloaded database whose linked seq/ was converted since does not give a stale order", {
    local_db_state()
    td <- tempfile("stale_")
    dir.create(td)
    withr::defer(unlink(td, recursive = TRUE))
    parent <- create_db_with_unsorted_chrom_sizes(file.path(td, "parent"))
    # a database made from parent with its own chrom_sizes.txt and a link to its seq/
    child <- file.path(td, "child")
    dir.create(file.path(child, "tracks", "t1.track"), recursive = TRUE)
    expect_true(file.copy(file.path(parent, "chrom_sizes.txt"), child))
    expect_true(file.symlink(file.path(parent, "seq"), file.path(child, "seq")))
    expect_equal(.gdb.chrom_names_at(child), c("chr1", "chr1_KI270706v1_random", "chr10", "chr2", "chrX"))

    suppressMessages(gdb.convert_to_indexed(groot = parent, force = TRUE, validate = FALSE))
    gsetroot(parent)
    expect_error(.gdb.chrom_names_at(child), "does not match chrom_sizes.txt")
    expect_error(misha:::.gtrack.liftover.src_chroms(file.path(child, "tracks", "t1.track")), "does not match chrom_sizes.txt")
})
