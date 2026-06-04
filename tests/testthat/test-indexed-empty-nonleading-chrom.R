load_test_db()

# Regression for the stale-mmap-window bug in GenomeTrackFixedBin (introduced by
# the mmap read path in commit 1cbfa801, "C++ optimization audit").
#
# Indexed dense tracks keep ONE GenomeTrackFixedBin object alive across all
# chromosomes (the persistent indexed backend in TrackExpressionVars). init_read
# sets up m_mmap_data / m_mmap_num_bins for a non-empty chromosome, but on a
# length-0 contig it returned early WITHOUT resetting those fields (the reset sat
# below the early return). read_next_bin / goto_bin bound the read on the stale
# m_mmap_num_bins (the previous chrom's bin count), not m_num_samples, so reading
# an empty contig AFTER a populated one returned the previous chromosome's values
# instead of NA.
#
# read_bins_bulk clamps on m_num_samples and was always safe, so single-function
# avg/sum/lse vtracks (mode 1/3) did not expose it. The bug surfaced on the paths
# that use read_next_bin/goto_bin: the general path (e.g. a max.pos vtrack) and
# mode 2 (avg+nearest), on single-bin iterators.
#
# Note the empty contig must be read in the SAME process right after the
# populated one, so the test runs single-process (gmultitasking = FALSE);
# otherwise different chroms can land in fresh worker backends with no stale
# state.

make_indexed_track_with_empty_trailing_chrom <- function(root, binsize = 20L,
                                                         chr1_size = 200L,
                                                         chr2_size = 200L) {
    dir.create(root)
    dir.create(file.path(root, "tracks"))
    dir.create(file.path(root, "seq"))
    writeLines(
        c(sprintf("chr1\t%d", chr1_size), sprintf("chr2\t%d", chr2_size)),
        file.path(root, "chrom_sizes.txt")
    )
    file.create(file.path(root, "seq", "chr1.seq"))
    file.create(file.path(root, "seq", "chr2.seq"))
    gsetroot(root)

    # Distinct per-bin values on chr1 so a stale read is unmistakable; chr2 will
    # become the length-0 entry. Build per-chrom, then drop chr2's file so the
    # pack writes a length-0 entry for it (mirror of the empty-leading-chrom
    # fixture, but with the empty chrom trailing a populated one).
    gtrack.create_dense(
        track = "t", description = "x",
        intervals = data.frame(
            chrom = c("chr1", "chr2"),
            start = c(0L, 0L),
            end = c(chr1_size, chr2_size)
        ),
        values = c(7, 7),
        binsize = binsize, defval = 0, func = "coverage"
    )
    trackdir <- file.path(root, "tracks", "t.track")
    file.remove(file.path(trackdir, "chr2")) # leave only chr1's per-chrom file
    expect_true(file.exists(file.path(trackdir, "chr1")))

    misha:::.gcall(
        "gtrack_pack_per_chrom_to_indexed",
        trackdir, c("chr1", "chr2"), "dense", misha:::.misha_env()
    )
    expect_true(file.exists(file.path(trackdir, "track.dat")))
    expect_true(file.exists(file.path(trackdir, "track.idx")))
    invisible(trackdir)
}

test_that("indexed track: empty contig after a populated one reads as NA (general path)", {
    old <- options(gmultitasking = FALSE)
    on.exit(options(old), add = TRUE)

    root <- tempfile("emptytrailing_")
    binsize <- 20L
    make_indexed_track_with_empty_trailing_chrom(root, binsize = binsize)

    remove_all_vtracks()
    # max.pos.abs routes through the generic read path (read_next_bin/goto_bin),
    # the one the stale mmap window corrupts.
    gvtrack.create("v_maxpos", "t", func = "max.pos.abs")

    scope <- data.frame(
        chrom = c("chr1", "chr2"),
        start = c(0L, 0L),
        end = c(200L, 200L)
    )
    res <- gextract("v_maxpos", intervals = scope, iterator = binsize)

    chr1 <- res[res$chrom == "chr1", ]
    chr2 <- res[res$chrom == "chr2", ]

    # chr1 has data in every bin -> max.pos is non-NA everywhere.
    expect_true(nrow(chr1) > 1)
    expect_false(any(is.na(chr1$v_maxpos)))

    # chr2 is an empty (length-0) contig: every bin must be NA, NOT chr1's
    # values leaked through the stale mmap window. Pre-fix, chr2 bins >= 1
    # returned non-NA positions.
    expect_true(nrow(chr2) > 1)
    expect_true(all(is.na(chr2$v_maxpos)))
})

test_that("indexed track: empty contig after a populated one reads as NA (avg + multi-function)", {
    old <- options(gmultitasking = FALSE)
    on.exit(options(old), add = TRUE)

    root <- tempfile("emptytrailing2_")
    binsize <- 20L
    make_indexed_track_with_empty_trailing_chrom(root, binsize = binsize)

    remove_all_vtracks()

    scope <- data.frame(
        chrom = c("chr1", "chr2"),
        start = c(0L, 0L),
        end = c(200L, 200L)
    )

    # Plain avg (single-function fast path, mode 3) was already safe via
    # read_bins_bulk; assert it stays correct as a control.
    avg_res <- gextract("t", intervals = scope, iterator = binsize)
    expect_false(any(is.na(avg_res$t[avg_res$chrom == "chr1"])))
    expect_true(all(is.na(avg_res$t[avg_res$chrom == "chr2"])))

    # avg + nearest in one expression set forces the avg/nearest fast path
    # (mode 2), which also uses read_next_bin/goto_bin.
    gvtrack.create("v_avg", "t", func = "avg")
    gvtrack.create("v_near", "t", func = "nearest")
    mn <- gextract("v_avg", "v_near", intervals = scope, iterator = binsize)
    expect_true(all(is.na(mn$v_avg[mn$chrom == "chr2"])))
    expect_true(all(is.na(mn$v_near[mn$chrom == "chr2"])))
})
