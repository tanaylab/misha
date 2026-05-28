load_test_db()

# Regression for the gtrack.copy SIGFPE / gtrack.info bin_size=0 bug.
#
# When an indexed dense track is packed with a per-chrom file list that does
# not include the genome's first chromosome (e.g. gtrack.copy with a
# destination DB whose leading chrom is not in the source), the resulting
# track.idx has length=0 for that chrom. The old GenomeTrackFixedBin::init_read
# returned early on the length=0 entry without populating m_bin_size, leaving
# it at the constructor default 0. gtrack.info then reported bin_size=0 and
# subsequent reads tripped the unguarded `interval.start / m_bin_size`
# division at GenomeTrackFixedBin.cpp:236/249/402/450/652/663 -> SIGFPE.
test_that("indexed dense track with empty leading chrom reports correct bin_size", {
    root <- tempfile("emptyleading_")
    dir.create(root)
    dir.create(file.path(root, "tracks"))
    dir.create(file.path(root, "seq"))
    writeLines(c("chr1\t10000", "chr2\t8000"), file.path(root, "chrom_sizes.txt"))
    file.create(file.path(root, "seq", "chr1.seq"))
    file.create(file.path(root, "seq", "chr2.seq"))
    gsetroot(root)

    # Build a dense track with data on chr2 only, by creating it normally then
    # removing chr1's per-chrom file. The pack will then see no chr1 file and
    # write a length=0 entry for it.
    gtrack.create_dense(
        track = "t", description = "x",
        intervals = data.frame(chrom = c("chr1", "chr2"), start = c(0L, 0L), end = c(10000L, 8000L)),
        values = c(1, 1),
        binsize = 20, defval = 0, func = "coverage"
    )
    trackdir <- file.path(root, "tracks", "t.track")
    file.remove(file.path(trackdir, "chr1")) # leave only chr2's per-chrom file
    expect_true(file.exists(file.path(trackdir, "chr2")))

    # Pack with both chroms in the dest chrom order; chr1 will get length=0.
    misha:::.gcall(
        "gtrack_pack_per_chrom_to_indexed",
        trackdir, c("chr1", "chr2"), "dense", misha:::.misha_env()
    )
    expect_true(file.exists(file.path(trackdir, "track.dat")))
    expect_true(file.exists(file.path(trackdir, "track.idx")))

    # gtrack.info reads the FIRST chrom by convention. Before the fix, that
    # chrom's length=0 entry skipped the bin_size read and returned 0. The fix
    # falls back to the first non-empty entry in the index.
    info <- gtrack.info("t")
    expect_equal(info$bin.size, 20L)

    # And reading chr2 (the populated chrom) must not be poisoned by the
    # earlier empty-chrom open that left m_bin_size at 0.
    out <- gextract("t", intervals = data.frame(chrom = "chr2", start = 0, end = 100), iterator = 20)
    expect_equal(nrow(out), 5L)
    expect_true(all(out$t == 1))
})
