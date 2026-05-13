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

# ---------------------------------------------------------------------
# Sparse (INTERVALS1D iterator) - Phase 4c.
# Builds a small indexed DB, creates a sparse iterator source via
# gtrack.create_sparse, then exercises the iterator-based gtrack.create
# (INTERVALS1D path) on both the indexed-direct and convert-after-create
# flows; asserts bit-for-bit equivalence.

# Build a small indexed DB and create a sparse "source" track via
# gtrack.create_sparse (which writes the indexed format directly via
# its own C path, unrelated to GenomeTrackIndexedWriter). Returns
# handles + the source track name. Source intervals deliberately skip
# some chroms to exercise the missing-chrom offset logic.
.make_indexed_db_with_sparse_source <- function(src_name = "src_sparse") {
    test_fasta <- tempfile(fileext = ".fasta")
    cat(">chr1\nACTGACTGACTGACTGACTGACTGACTGACTG\n",
        ">chr2\nGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCC\n",
        ">chr3\nTATATATATATATATATATATATATATATATA\n",
        ">chr4\nAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n",
        ">chr5\nCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n",
        file = test_fasta, sep = ""
    )

    test_db <- tempfile(pattern = "misha_idx_sparse_")

    withr::with_options(list(gmulticontig.indexed_format = TRUE), {
        gdb.create(groot = test_db, fasta = test_fasta, verbose = FALSE)
        gdb.init(test_db)
    })

    # Sparse source intervals on chr1, chr3, chr5 only - chr2 and chr4
    # are deliberately absent to exercise the missing-chrom offset path
    # in finalize().
    src_intervs <- data.frame(
        chrom = c("chr1", "chr1", "chr3", "chr5", "chr5"),
        start = c(0L, 16L, 8L, 0L, 20L),
        end = c(8L, 24L, 24L, 12L, 28L),
        stringsAsFactors = FALSE
    )
    gtrack.create_sparse(src_name, "sparse source", src_intervs,
        values = seq_len(nrow(src_intervs))
    )

    list(db = test_db, fasta = test_fasta, src = src_name)
}

test_that("indexed-direct sparse matches convert-after-create bit-for-bit", {
    handles_a <- .make_indexed_db_with_sparse_source("src_sparse")
    withr::defer({
        unlink(handles_a$db, recursive = TRUE)
        unlink(handles_a$fasta)
    })

    # Direct path: gtrack.create on the indexed DB. Note: the source
    # track's data is on chr1, chr3, chr5; the iterator scope is
    # gintervals.all() so chr2 and chr4 must be backfilled with
    # header-only entries.
    gtrack.create("idx_sparse_direct", "via direct", "src_sparse",
        iterator = "src_sparse"
    )

    direct_dir <- .track_dir("idx_sparse_direct")
    direct_dat <- readBin(file.path(direct_dir, "track.dat"),
        what = "raw",
        n = file.info(file.path(direct_dir, "track.dat"))$size + 1
    )
    direct_idx <- readBin(file.path(direct_dir, "track.idx"),
        what = "raw",
        n = file.info(file.path(direct_dir, "track.idx"))$size + 1
    )

    # Comparison artifact in a SECOND indexed DB: re-create source +
    # consumer by hiding genome.idx before the consumer create so the
    # per-chrom + convert path runs.
    handles_b <- .make_indexed_db_with_sparse_source("src_sparse")
    withr::defer({
        unlink(handles_b$db, recursive = TRUE)
        unlink(handles_b$fasta)
    })

    seq_dir <- file.path(handles_b$db, "seq")
    idx_path <- file.path(seq_dir, "genome.idx")
    hidden_path <- file.path(seq_dir, "genome.idx.hidden")
    expect_true(file.rename(idx_path, hidden_path))
    withr::defer({
        if (file.exists(hidden_path)) file.rename(hidden_path, idx_path)
    })
    gtrack.create("idx_sparse_via_convert", "via convert", "src_sparse",
        iterator = "src_sparse"
    )

    # Restore idx, then convert (mirroring legacy two-pass flow).
    file.rename(hidden_path, idx_path)
    gdb.reload()
    gtrack.convert_to_indexed("idx_sparse_via_convert")

    convert_dir <- .track_dir("idx_sparse_via_convert")
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

    # Functional equivalence (open in handles_b, which has both).
    intv <- gintervals.all()
    r_direct <- gextract("idx_sparse_via_convert", intervals = intv)
    expect_true(nrow(r_direct) > 0)
})

test_that("indexed-direct sparse: no per-chrom files left over", {
    handles <- .make_indexed_db_with_sparse_source("src_sparse_clean")
    withr::defer({
        unlink(handles$db, recursive = TRUE)
        unlink(handles$fasta)
    })

    gtrack.create("idx_sparse_clean", "no leftover", "src_sparse_clean",
        iterator = "src_sparse_clean"
    )

    dir_contents <- list.files(.track_dir("idx_sparse_clean"),
        all.files = TRUE, no.. = TRUE
    )
    expect_setequal(dir_contents, c("track.dat", "track.idx", ".attributes"))
})

test_that("indexed-direct sparse produces a working track", {
    handles <- .make_indexed_db_with_sparse_source("src_sparse_func")
    withr::defer({
        unlink(handles$db, recursive = TRUE)
        unlink(handles$fasta)
    })

    gtrack.create("idx_sparse_func", "functional", "src_sparse_func",
        iterator = "src_sparse_func"
    )

    info <- gtrack.info("idx_sparse_func")
    expect_equal(info$type, "sparse")

    val <- gextract("idx_sparse_func", intervals = gintervals.all())
    expect_gt(nrow(val), 0)
    expect_true(all(!is.na(val$idx_sparse_func)))
})

# ---------------------------------------------------------------------
# Liftover (gtrack.liftover) - Phase 6a.
# Builds a source DB with a sparse (or dense) track, an indexed target DB,
# and a chain file mapping a region of source -> target. Lifts the source
# track into the target DB twice: once with the target idx visible (direct
# indexed write path) and once with target idx hidden so the per-chrom +
# convert path runs. Asserts bit-for-bit equivalence of track.dat and
# track.idx.

# Build a small source DB with a sparse source track. Caller is
# responsible for the resulting db/fasta cleanup via withr::defer.
.make_liftover_source_sparse <- function(src_name) {
    test_fasta <- tempfile(fileext = ".fasta")
    cat(">chrSrc1\n", paste(rep("A", 64), collapse = ""), "\n",
        ">chrSrc2\n", paste(rep("C", 64), collapse = ""), "\n",
        file = test_fasta, sep = ""
    )

    test_db <- tempfile(pattern = "misha_liftover_src_sparse_")
    withr::with_options(list(gmulticontig.indexed_format = TRUE), {
        gdb.create(groot = test_db, fasta = test_fasta, verbose = FALSE)
        gdb.init(test_db)
    })

    src_intervs <- data.frame(
        chrom = c("chrSrc1", "chrSrc1", "chrSrc1", "chrSrc2", "chrSrc2"),
        start = c(0L, 16L, 40L, 0L, 32L),
        end = c(8L, 24L, 56L, 16L, 48L),
        stringsAsFactors = FALSE
    )
    gtrack.create_sparse(src_name, "sparse source", src_intervs,
        values = c(1.5, 2.25, 3.0, 4.125, 5.5)
    )

    list(db = test_db, fasta = test_fasta, src = src_name)
}

# Build an empty indexed target DB on >2 chroms (one of which the chain
# will not cover, exercising the header-only empty-chrom path).
.make_liftover_target_db <- function() {
    test_fasta <- tempfile(fileext = ".fasta")
    cat(">chr1\n", paste(rep("T", 96), collapse = ""), "\n",
        ">chr2\n", paste(rep("G", 96), collapse = ""), "\n",
        ">chr3\n", paste(rep("A", 96), collapse = ""), "\n",
        file = test_fasta, sep = ""
    )
    test_db <- tempfile(pattern = "misha_liftover_tgt_")
    withr::with_options(list(gmulticontig.indexed_format = TRUE), {
        gdb.create(groot = test_db, fasta = test_fasta, verbose = FALSE)
        gdb.init(test_db)
    })
    list(db = test_db, fasta = test_fasta)
}

.write_liftover_chain <- function(chain_file) {
    # chain 1: chrSrc1[0..56) + -> chr1[0..56) + (length 56)
    cat("chain 1000 chrSrc1 64 + 0 56 chr1 96 + 0 56 1\n",
        "56\n\n",
        sep = "", file = chain_file
    )
    # chain 2: chrSrc2[0..48) + -> chr2[0..48) + (length 48)
    # Note: chr3 is intentionally NOT covered by any chain - exercises the
    # empty-chrom header-only path in the indexed writer.
    cat("chain 1000 chrSrc2 64 + 0 48 chr2 96 + 0 48 2\n",
        "48\n\n",
        sep = "", file = chain_file, append = TRUE
    )
}

test_that("indexed-direct sparse liftover matches convert-after-create bit-for-bit", {
    # --- Source DB (shared by both runs) ---
    src_handles <- .make_liftover_source_sparse("src_for_liftover")
    withr::defer({
        unlink(src_handles$db, recursive = TRUE)
        unlink(src_handles$fasta)
    })
    src_track_dir <- file.path(src_handles$db, "tracks", "src_for_liftover.track")

    # --- Chain file (shared) ---
    chain_file <- tempfile(fileext = ".chain")
    withr::defer(unlink(chain_file))
    .write_liftover_chain(chain_file)

    # --- Direct path: liftover into a fresh indexed target DB ---
    tgt_a <- .make_liftover_target_db()
    withr::defer({
        unlink(tgt_a$db, recursive = TRUE)
        unlink(tgt_a$fasta)
    })
    gtrack.liftover("lifted_direct", "via direct", src_track_dir, chain_file)
    direct_dir <- .track_dir("lifted_direct")
    expect_true(file.exists(file.path(direct_dir, "track.dat")))
    expect_true(file.exists(file.path(direct_dir, "track.idx")))
    direct_dat <- readBin(file.path(direct_dir, "track.dat"),
        what = "raw",
        n = file.info(file.path(direct_dir, "track.dat"))$size + 1
    )
    direct_idx <- readBin(file.path(direct_dir, "track.idx"),
        what = "raw",
        n = file.info(file.path(direct_dir, "track.idx"))$size + 1
    )

    # --- Convert path: build a second target DB, hide genome.idx, lift,
    #     restore idx, convert.
    tgt_b <- .make_liftover_target_db()
    withr::defer({
        unlink(tgt_b$db, recursive = TRUE)
        unlink(tgt_b$fasta)
    })

    seq_dir <- file.path(tgt_b$db, "seq")
    idx_path <- file.path(seq_dir, "genome.idx")
    hidden_path <- file.path(seq_dir, "genome.idx.hidden")
    expect_true(file.rename(idx_path, hidden_path))
    withr::defer({
        if (file.exists(hidden_path)) file.rename(hidden_path, idx_path)
    })

    gdb.reload()
    gtrack.liftover("lifted_via_convert", "via convert", src_track_dir, chain_file)

    # Restore idx + reload + convert.
    file.rename(hidden_path, idx_path)
    gdb.reload()
    gtrack.convert_to_indexed("lifted_via_convert")

    convert_dir <- .track_dir("lifted_via_convert")
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

test_that("indexed-direct liftover: no per-chrom files left over", {
    src_handles <- .make_liftover_source_sparse("src_clean_liftover")
    withr::defer({
        unlink(src_handles$db, recursive = TRUE)
        unlink(src_handles$fasta)
    })
    src_track_dir <- file.path(src_handles$db, "tracks", "src_clean_liftover.track")

    chain_file <- tempfile(fileext = ".chain")
    withr::defer(unlink(chain_file))
    .write_liftover_chain(chain_file)

    tgt <- .make_liftover_target_db()
    withr::defer({
        unlink(tgt$db, recursive = TRUE)
        unlink(tgt$fasta)
    })

    gtrack.liftover("lifted_clean", "no leftover", src_track_dir, chain_file)

    dir_contents <- list.files(.track_dir("lifted_clean"),
        all.files = TRUE, no.. = TRUE
    )
    expect_setequal(dir_contents, c("track.dat", "track.idx", ".attributes"))
})

test_that("indexed-direct liftover: gextract returns identical values both ways", {
    src_handles <- .make_liftover_source_sparse("src_func_liftover")
    withr::defer({
        unlink(src_handles$db, recursive = TRUE)
        unlink(src_handles$fasta)
    })
    src_track_dir <- file.path(src_handles$db, "tracks", "src_func_liftover.track")

    chain_file <- tempfile(fileext = ".chain")
    withr::defer(unlink(chain_file))
    .write_liftover_chain(chain_file)

    # Run 1: direct
    tgt_a <- .make_liftover_target_db()
    withr::defer({
        unlink(tgt_a$db, recursive = TRUE)
        unlink(tgt_a$fasta)
    })
    gtrack.liftover("lifted_func_direct", "direct", src_track_dir, chain_file)
    res_direct <- gextract("lifted_func_direct", gintervals.all())

    # Run 2: hidden idx + convert
    tgt_b <- .make_liftover_target_db()
    withr::defer({
        unlink(tgt_b$db, recursive = TRUE)
        unlink(tgt_b$fasta)
    })
    seq_dir <- file.path(tgt_b$db, "seq")
    idx_path <- file.path(seq_dir, "genome.idx")
    hidden_path <- file.path(seq_dir, "genome.idx.hidden")
    expect_true(file.rename(idx_path, hidden_path))
    withr::defer({
        if (file.exists(hidden_path)) file.rename(hidden_path, idx_path)
    })
    gdb.reload()
    gtrack.liftover("lifted_func_convert", "convert", src_track_dir, chain_file)
    file.rename(hidden_path, idx_path)
    gdb.reload()
    gtrack.convert_to_indexed("lifted_func_convert")

    res_convert <- gextract("lifted_func_convert", gintervals.all())

    # Strip the (different) track-name column from each frame, compare the rest.
    res_direct$lifted_func_direct <- NULL
    res_convert$lifted_func_convert <- NULL
    expect_equal(
        res_direct[order(res_direct$chrom, res_direct$start), ],
        res_convert[order(res_convert$chrom, res_convert$start), ],
        ignore_attr = TRUE
    )
})
