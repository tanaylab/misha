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
