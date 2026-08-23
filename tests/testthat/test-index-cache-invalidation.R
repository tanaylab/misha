# Regression tests for the process-static index cache invalidation bug.
#
# Background: misha keeps process-static caches keyed by directory path
# for the on-disk per-dir indexes:
#   - GenomeTrack::s_index_cache           (1D track.idx)
#   - TrackIndex2D::s_index_cache          (2D track.idx)
#   - GIntervalsBigSet1D::s_index_cache    (intervals.idx)
#   - GIntervalsBigSet2D::s_index_cache    (intervals2d.idx)
#
# Without explicit invalidation, the sequence
#   gtrack.rm("X"); gtrack.create_*("X", ...)
# on an indexed DB leaves the cache pointing at the previous lifecycle's
# index. Downstream readers see "indexed = true" and try to open a
# non-existent track.dat (or the wrong file layout).
#
# These tests exercise the rm-then-recreate cycle and verify the next
# read succeeds.

restore_groot_on_exit()

# Build a tiny per-chromosome DB then convert it to indexed format.
build_indexed_test_db <- function() {
    test_db <- tempfile("misha_idxcache_")
    dir.create(test_db)
    dir.create(file.path(test_db, "seq"))
    dir.create(file.path(test_db, "tracks"))
    chrom_sizes <- data.frame(
        chrom = c("chr1", "chr2", "chr3"),
        size  = c(20000L, 15000L, 10000L)
    )
    write.table(chrom_sizes, file.path(test_db, "chrom_sizes.txt"),
        row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE
    )
    for (i in seq_len(nrow(chrom_sizes))) {
        writeBin(
            charToRaw(paste(rep("A", chrom_sizes$size[i]), collapse = "")),
            file.path(test_db, "seq", paste0(chrom_sizes$chrom[i], ".seq"))
        )
    }
    gdb.convert_to_indexed(groot = test_db, force = TRUE, validate = FALSE)
    test_db
}

# Per-chromosome (non-indexed) test DB.
build_perchrom_test_db <- function() {
    test_db <- tempfile("misha_idxcache_pc_")
    dir.create(test_db)
    dir.create(file.path(test_db, "seq"))
    dir.create(file.path(test_db, "tracks"))
    chrom_sizes <- data.frame(
        chrom = c("chr1", "chr2", "chr3"),
        size  = c(20000L, 15000L, 10000L)
    )
    write.table(chrom_sizes, file.path(test_db, "chrom_sizes.txt"),
        row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE
    )
    for (i in seq_len(nrow(chrom_sizes))) {
        writeBin(
            charToRaw(paste(rep("A", chrom_sizes$size[i]), collapse = "")),
            file.path(test_db, "seq", paste0(chrom_sizes$chrom[i], ".seq"))
        )
    }
    test_db
}

make_sparse_data <- function(n = 200L, chrom = "chr1", chrom_size = 20000L,
                             width = 50L) {
    # Sparse intervals must be non-overlapping. Pick spacing >= width so
    # adjacent intervals don't collide; clamp n to what fits.
    span <- chrom_size - width
    spacing <- max(width + 1L, span %/% n)
    max_n <- as.integer(span %/% spacing)
    if (max_n < n) n <- max_n
    starts <- as.integer(seq.int(0L, by = spacing, length.out = n))
    data.frame(
        chrom = chrom, start = starts, end = starts + width,
        stringsAsFactors = FALSE
    )
}

test_that("gtrack.rm + gtrack.create_dense cycle works on indexed DB", {
    local_db_state()
    test_db <- build_indexed_test_db()
    withr::defer(unlink(test_db, recursive = TRUE))
    gsetroot(test_db)

    ivs <- make_sparse_data(n = 200L)
    vals <- runif(nrow(ivs))

    # Three rm-then-recreate cycles. Before the fix, rep 2's
    # gtrack.create_dense errored with
    #   "Cannot open .../<track>.track/track.dat: No such file or directory"
    # because the stale s_index_cache entry from rep 1 routed reads through
    # the indexed-format path on a directory that had been deleted.
    for (i in 1:3) {
        if ("idxcache_dense" %in% gtrack.ls()) {
            suppressMessages(gtrack.rm("idxcache_dense", force = TRUE))
        }
        expect_silent(
            gtrack.create_dense("idxcache_dense", "test",
                ivs, vals,
                binsize = 200L, defval = NaN
            )
        )
        # Read back through the normal extract path - this is what fails
        # if the cache is stale on the second iteration.
        ext <- gextract("idxcache_dense",
            intervals = data.frame(
                chrom = "chr1",
                start = 0L, end = 5000L
            ),
            iterator = 200L
        )
        expect_true(nrow(ext) > 0L,
            label = sprintf("rep %d: gextract returned no rows", i)
        )
    }
})

test_that("gtrack.rm + gtrack.create_sparse cycle works on indexed DB", {
    local_db_state()
    test_db <- build_indexed_test_db()
    withr::defer(unlink(test_db, recursive = TRUE))
    gsetroot(test_db)

    ivs <- make_sparse_data(n = 200L)
    vals <- runif(nrow(ivs))

    for (i in 1:3) {
        if ("idxcache_sparse" %in% gtrack.ls()) {
            suppressMessages(gtrack.rm("idxcache_sparse", force = TRUE))
        }
        expect_silent(
            gtrack.create_sparse("idxcache_sparse", "test", ivs, vals)
        )
        ext <- gextract("idxcache_sparse",
            intervals = data.frame(
                chrom = "chr1",
                start = 0L, end = 5000L
            )
        )
        expect_true(nrow(ext) > 0L,
            label = sprintf("rep %d: gextract returned no rows", i)
        )
    }
})

test_that("convert_to_indexed in place invalidates the cache", {
    local_db_state()
    test_db <- build_perchrom_test_db()
    withr::defer(unlink(test_db, recursive = TRUE))
    gsetroot(test_db)

    ivs <- make_sparse_data(n = 200L)
    vals <- runif(nrow(ivs))

    # Per-chrom create.
    gtrack.create_dense("idxcache_conv", "test",
        ivs, vals,
        binsize = 200L, defval = NaN
    )

    # Read once so the cache picks up the per-chrom state (nullptr entry).
    invisible(gextract("idxcache_conv",
        intervals = data.frame(
            chrom = "chr1",
            start = 0L, end = 5000L
        ),
        iterator = 200L
    ))

    # Convert per-chrom -> indexed. The cache MUST drop the prior nullptr
    # entry; otherwise readers keep probing per-chrom files that now
    # live behind track.dat + track.idx.
    expect_silent(gtrack.convert_to_indexed("idxcache_conv"))

    ext <- gextract("idxcache_conv",
        intervals = data.frame(
            chrom = "chr1",
            start = 0L, end = 5000L
        ),
        iterator = 200L
    )
    expect_true(nrow(ext) > 0L)
})

test_that("gintervals.rm + gintervals.save cycle works on indexed DB", {
    local_db_state()
    test_db <- build_indexed_test_db()
    withr::defer(unlink(test_db, recursive = TRUE))
    gsetroot(test_db)

    ivs <- make_sparse_data(n = 500L)

    for (i in 1:3) {
        if ("idxcache_iv" %in% gintervals.ls()) {
            suppressMessages(gintervals.rm("idxcache_iv", force = TRUE))
        }
        expect_silent(gintervals.save("idxcache_iv", ivs))
        loaded <- gintervals.load("idxcache_iv")
        expect_true(nrow(loaded) > 0L,
            label = sprintf("rep %d: gintervals.load returned no rows", i)
        )
    }
})

test_that("explicit cache invalidation is exposed and idempotent", {
    local_db_state()
    test_db <- build_indexed_test_db()
    withr::defer(unlink(test_db, recursive = TRUE))
    gsetroot(test_db)

    # The R-level helper must accept missing/empty/nonexistent paths
    # without erroring (it is called from cleanup paths where we don't
    # want to mask real errors).
    expect_silent(misha:::.gdb.invalidate_dir_cache(character(0)))
    expect_silent(misha:::.gdb.invalidate_dir_cache(""))
    expect_silent(misha:::.gdb.invalidate_dir_cache("/no/such/dir"))
    expect_silent(misha:::.gdb.invalidate_dir_cache(c(tempdir(), tempdir())))
})

# --- Mutation points other than rm/create ------------------------------------
#
# Any R-level operation that replaces the contents of a track or interval-set
# directory has to drop the cached index keyed by that path. The tests below
# each reproduce a confirmed wrong answer that appears when the corresponding
# .gdb.invalidate_dir_cache() call is removed.

test_that("gtrack.mv invalidates the cache of both source and destination", {
    local_db_state()
    test_db <- build_indexed_test_db()
    withr::defer(unlink(test_db, recursive = TRUE))
    gsetroot(test_db)

    ivs <- make_sparse_data(n = 190L)
    vals <- seq_len(nrow(ivs)) * 1.0
    gtrack.create_dense("mv_dense", "d", ivs, vals, binsize = 200L, defval = NaN)
    gtrack.create_sparse("mv_sparse", "s", ivs, vals)

    expect_equal(gtrack.info("mv_dense")$type, "dense")
    expect_equal(gtrack.info("mv_sparse")$type, "sparse")

    # Move the sparse track onto the name (and directory) the dense track
    # occupied. Without invalidation the cached index from "mv_dense" is
    # reused for the new occupant of that path: gtrack.info() answers
    # "dense" and gextract() errors with a bin-count mismatch.
    gtrack.mv("mv_dense", "mv_moved")
    gtrack.mv("mv_sparse", "mv_dense")

    expect_equal(gtrack.info("mv_dense")$type, "sparse")
    expect_equal(gtrack.info("mv_moved")$type, "dense")
    ext <- gextract("mv_dense",
        intervals = data.frame(chrom = "chr1", start = 0L, end = 5000L),
        iterator = 200L
    )
    expect_true(nrow(ext) > 0L)
})

test_that("gintervals.update big->small->big invalidates the bigset index", {
    local_db_state()
    test_db <- build_indexed_test_db()
    withr::defer(unlink(test_db, recursive = TRUE))
    gsetroot(test_db)

    withr::local_options(
        gmulticontig.indexed_format = TRUE,
        gmultitasking = FALSE
    )

    mk <- function(chrom, n, off = 0L, w = 20L, step = 100L) {
        s <- as.integer(off + step * seq_len(n) - step)
        gintervals(chrom, s, s + w)
    }

    # A gradient track so a summary over the set depends on which intervals
    # the C++ iterator actually visits.
    grad_ivs <- rbind(
        mk("chr1", 190L, step = 100L, w = 90L),
        mk("chr2", 140L, step = 100L, w = 90L)
    )
    gtrack.create_dense("bigset_grad", "g", grad_ivs,
        seq_len(nrow(grad_ivs)) * 1.0,
        binsize = 100L, defval = NaN
    )

    # Force the set to be a big set that is too large to load into R, so
    # reads go through GIntervalsBigSet1D (the class that caches the index)
    # rather than the fresh-load path used for small sets.
    withr::local_options(gbig.intervals.size = 50, gmax.data.size = 20)

    gintervals.save("bigset_upd", rbind(mk("chr1", 60L), mk("chr2", 15L)))
    expect_true(misha:::.gintervals.is_indexed_bigset("bigset_upd"))
    first <- gsummary("bigset_grad", intervals = "bigset_upd", iterator = "bigset_upd")
    expect_equal(as.numeric(first[1]), 75)

    # Shrink below the big-set threshold: .gintervals.big2small replaces the
    # directory with a flat file.
    gintervals.update("bigset_upd", NULL, chrom = "chr1")
    expect_false(misha:::.gintervals.is_bigset("bigset_upd"))

    # Grow back above the threshold: .gintervals.small2big writes a *new*
    # indexed big set at the same path, with different offsets. Without
    # invalidation the first read below fails with "unknown input format"
    # because the stale index seeks to lifecycle-1 offsets in the new
    # intervals.dat.
    gintervals.update("bigset_upd", mk("chr1", 80L, off = 5000L, w = 33L), chrom = "chr1")
    expect_true(misha:::.gintervals.is_indexed_bigset("bigset_upd"))

    second <- gsummary("bigset_grad", intervals = "bigset_upd", iterator = "bigset_upd")
    expect_equal(as.numeric(second[1]), 95)
})

test_that("gtrack.copy invalidates the destination's cached index", {
    local_db_state()
    test_db <- build_indexed_test_db()
    withr::defer(unlink(test_db, recursive = TRUE))
    gsetroot(test_db)

    ivs <- make_sparse_data(n = 190L)
    vals <- seq_len(nrow(ivs)) * 1.0
    gtrack.create_dense("cp_dense", "d", ivs, vals, binsize = 200L, defval = NaN)
    gtrack.create_sparse("cp_sparse", "s", ivs, vals)

    # Read the dense track so its index is cached.
    expect_true(nrow(gextract("cp_dense",
        intervals = data.frame(chrom = "chr1", start = 0L, end = 5000L),
        iterator = 200L
    )) > 0L)

    # Another process (a sibling R session, a shell rm -rf) removes the track
    # directory. This session's cache still holds the entry for that path.
    dense_dir <- misha:::.track_dir("cp_dense")
    system(sprintf("rm -rf %s", shQuote(dense_dir)))
    expect_false(dir.exists(dense_dir))

    # Copying a differently-shaped track onto the same path must not be read
    # through the stale entry.
    gtrack.copy("cp_sparse", "cp_dense", overwrite = TRUE)
    expect_equal(gtrack.info("cp_dense")$type, "sparse")
    ext <- gextract("cp_dense",
        intervals = data.frame(chrom = "chr1", start = 0L, end = 5000L),
        iterator = 200L
    )
    expect_true(nrow(ext) > 0L)
})
