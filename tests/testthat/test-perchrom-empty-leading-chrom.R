load_test_db()

# Per-chrom analog of the indexed empty-leading-chrom bug (#133).
#
# A per-chromosome dense track may legitimately lack a per-chrom file for the
# genome's first chromosome (e.g. a scaffold with no signal, or a partial
# track). gtrack.info() probes the FIRST genome chrom by convention; the old
# code resolved that chrom's filename and called init_read() on it directly,
# which errors with "No such file or directory" when the file is absent. Since
# gtrack.copy() internally calls gtrack.info(srcname), copying such a track
# failed outright. gtrack.info() must instead fall back to the first chrom that
# actually has data (bin_size is invariant across chroms).
test_that("per-chrom dense track with missing leading chrom reports correct bin_size", {
    root <- tempfile("pc_emptyleading_")
    dir.create(root)
    dir.create(file.path(root, "tracks"))
    dir.create(file.path(root, "seq"))
    writeLines(c("chr1\t10000", "chr2\t8000"), file.path(root, "chrom_sizes.txt"))
    file.create(file.path(root, "seq", "chr1.seq"))
    file.create(file.path(root, "seq", "chr2.seq"))
    gsetroot(root)

    gtrack.create_dense(
        track = "t", description = "x",
        intervals = data.frame(chrom = c("chr1", "chr2"), start = c(0L, 0L), end = c(10000L, 8000L)),
        values = c(1, 1),
        binsize = 20, defval = 0, func = "coverage"
    )
    trackdir <- file.path(root, "tracks", "t.track")
    # Remove chr1's per-chrom file: the track now lacks the genome's first chrom.
    file.remove(file.path(trackdir, "chr1"))
    expect_true(file.exists(file.path(trackdir, "chr2")))
    expect_false(file.exists(file.path(trackdir, "track.idx"))) # still per-chrom

    info <- gtrack.info("t")
    expect_equal(info$bin.size, 20L)
    expect_equal(info$format, "per-chromosome")
})

test_that("gtrack.copy of a per-chrom track missing its leading chrom succeeds", {
    root <- tempfile("pc_copy_emptyleading_")
    dir.create(root)
    dir.create(file.path(root, "tracks"))
    dir.create(file.path(root, "seq"))
    writeLines(c("chr1\t10000", "chr2\t8000"), file.path(root, "chrom_sizes.txt"))
    file.create(file.path(root, "seq", "chr1.seq"))
    file.create(file.path(root, "seq", "chr2.seq"))
    gsetroot(root)

    gtrack.create_dense(
        track = "src", description = "x",
        intervals = data.frame(chrom = c("chr1", "chr2"), start = c(0L, 0L), end = c(10000L, 8000L)),
        values = c(1, 1),
        binsize = 20, defval = 0, func = "coverage"
    )
    file.remove(file.path(root, "tracks", "src.track", "chr1"))

    expect_silent(suppressWarnings(gtrack.copy("src", "dst")))
    expect_equal(gtrack.info("dst")$bin.size, 20L)
})
