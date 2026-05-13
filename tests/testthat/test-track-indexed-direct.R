# Tests for the streaming indexed-direct fixed-bin writer
# (Phase 4b of the large-contig fixes). The headline guarantee is
# bit-for-bit equivalence with the legacy per-chrom + pack flow.

# Build a small indexed DB from a 5-chromosome FASTA and create the
# named track via gtrack.create. Returns the path to the track dir.
.make_indexed_db_with_track <- function(track_name, bin_size = 4) {
    test_fasta <- tempfile(fileext = ".fasta")
    cat(">chr1\nACTGACTGACTGACTGACTGACTGACTGACTG\n",
        ">chr2\nGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCC\n",
        ">chr3\nTATATATATATATATATATATATATATATATA\n",
        ">chr4\nAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n",
        ">chr5\nCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n",
        file = test_fasta, sep = ""
    )

    test_db <- tempfile(pattern = "misha_idx_direct_")

    # Caller owns cleanup via withr::defer.
    withr::with_options(list(gmulticontig.indexed_format = TRUE), {
        gdb.create(groot = test_db, fasta = test_fasta, verbose = FALSE)
        gdb.init(test_db)
    })

    # A trivial fixed-bin expression: position-derived counts.
    gtrack.create(track_name, "indexed-direct test", "1", iterator = bin_size)

    list(db = test_db, fasta = test_fasta)
}

test_that("indexed-direct fixed-bin matches convert-after-create bit-for-bit", {
    handles_a <- .make_indexed_db_with_track("idx_direct", bin_size = 4)
    withr::defer({
        unlink(handles_a$db, recursive = TRUE)
        unlink(handles_a$fasta)
    })
    direct_dir <- .track_dir("idx_direct")
    direct_dat <- readBin(file.path(direct_dir, "track.dat"),
        what = "raw",
        n = file.info(file.path(direct_dir, "track.dat"))$size + 1
    )
    direct_idx <- readBin(file.path(direct_dir, "track.idx"),
        what = "raw",
        n = file.info(file.path(direct_dir, "track.idx"))$size + 1
    )

    # Comparison artifact: build a second indexed DB, hide genome.idx
    # before calling gtrack.create, then convert. This forces the per-chrom
    # write + pack flow that produced track.dat/track.idx prior to Phase 4b.
    test_fasta <- tempfile(fileext = ".fasta")
    cat(">chr1\nACTGACTGACTGACTGACTGACTGACTGACTG\n",
        ">chr2\nGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCC\n",
        ">chr3\nTATATATATATATATATATATATATATATATA\n",
        ">chr4\nAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n",
        ">chr5\nCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n",
        file = test_fasta, sep = ""
    )
    test_db <- tempfile(pattern = "misha_idx_convert_")
    withr::defer({
        unlink(test_db, recursive = TRUE)
        unlink(test_fasta)
    })

    withr::with_options(list(gmulticontig.indexed_format = TRUE), {
        gdb.create(groot = test_db, fasta = test_fasta, verbose = FALSE)
        gdb.init(test_db)
    })

    # Hide genome.idx so .gdb.is_indexed() returns FALSE and
    # gtrack.create takes the per-chrom branch.
    seq_dir <- file.path(test_db, "seq")
    idx_path <- file.path(seq_dir, "genome.idx")
    hidden_path <- file.path(seq_dir, "genome.idx.hidden")
    expect_true(file.rename(idx_path, hidden_path))
    withr::defer({
        if (file.exists(hidden_path)) file.rename(hidden_path, idx_path)
    })
    gtrack.create("idx_via_convert", "via convert", "1", iterator = 4)

    # Restore idx, then run convert (mirroring the legacy two-pass flow).
    file.rename(hidden_path, idx_path)
    gdb.reload()
    gtrack.convert_to_indexed("idx_via_convert")

    convert_dir <- .track_dir("idx_via_convert")
    convert_dat <- readBin(file.path(convert_dir, "track.dat"),
        what = "raw",
        n = file.info(file.path(convert_dir, "track.dat"))$size + 1
    )
    convert_idx <- readBin(file.path(convert_dir, "track.idx"),
        what = "raw",
        n = file.info(file.path(convert_dir, "track.idx"))$size + 1
    )

    expect_identical(direct_dat, convert_dat)
    expect_identical(direct_idx, convert_idx)
})

test_that("indexed-direct: no per-chrom files left over", {
    handles <- .make_indexed_db_with_track("idx_clean", bin_size = 4)
    withr::defer({
        unlink(handles$db, recursive = TRUE)
        unlink(handles$fasta)
    })

    dir_contents <- list.files(.track_dir("idx_clean"), all.files = TRUE, no.. = TRUE)
    # Only track.dat + track.idx + .attributes should remain.
    expect_setequal(dir_contents, c("track.dat", "track.idx", ".attributes"))
})

test_that("indexed-direct produces a working track (basic functional checks)", {
    handles <- .make_indexed_db_with_track("idx_func", bin_size = 4)
    withr::defer({
        unlink(handles$db, recursive = TRUE)
        unlink(handles$fasta)
    })

    info <- gtrack.info("idx_func")
    expect_equal(info$type, "dense")

    val <- gextract("idx_func", intervals = gintervals.all(), iterator = 4)
    expect_gt(nrow(val), 0)
    # Track expression "1" -> every bin should be 1.
    expect_true(all(val$idx_func == 1))

    # Coarser iterator should still produce data.
    val_coarse <- gextract("idx_func", intervals = gintervals.all(), iterator = 8)
    expect_gt(nrow(val_coarse), 0)
})

test_that("indexed-direct track.idx has one entry per chromosome", {
    # Sanity: track.idx must contain N entries (N = number of chroms in
    # the genome), matching what gtrack_pack_per_chrom_to_indexed
    # produces. Header size + entry table size are deterministic.
    handles <- .make_indexed_db_with_track("idx_full", bin_size = 4)
    withr::defer({
        unlink(handles$db, recursive = TRUE)
        unlink(handles$fasta)
    })

    idx_path <- file.path(.track_dir("idx_full"), "track.idx")
    idx_size <- file.info(idx_path)$size

    # Header is 36 bytes (see TrackIndex.h). Each entry is 24 bytes.
    n_chroms <- nrow(gintervals.all())
    expect_equal(idx_size, 36 + 24 * n_chroms)
})
