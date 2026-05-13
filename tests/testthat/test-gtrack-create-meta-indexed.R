# Regression tests for indexed-aware shortcut in GTrackIntervalsFetcher::create_track_meta.
#
# When a sparse / arrays track lives in a misha DB and gets used as an
# intervals set, the C++ side lazily builds a .meta file by scanning every
# chromid. On indexed-format DBs there are no per-chrom files at all, so the
# old per-chrom probe (find_existing_1d_filename + access) wastes ~5 syscalls
# per chromid. The fix uses track.idx directly to enumerate chromids that
# have data (length > 0).
#
# These tests verify:
#  1. .meta produced on an indexed sparse track is functionally equivalent
#     to .meta produced on the same data in per-chrom format.
#  2. Meta creation completes for an indexed track where the overwhelming
#     majority of chromids have no data (the large-contig case).

setup_db <- function(num_chroms, chrom_size = 1e6, indexed = FALSE) {
    tmp_root <- withr::local_tempdir(.local_envir = parent.frame())
    chrom_sizes <- data.frame(
        chrom = paste0("chr", seq_len(num_chroms)),
        size = rep(chrom_size, num_chroms)
    )
    utils::write.table(
        chrom_sizes, file.path(tmp_root, "chrom_sizes.txt"),
        sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE
    )
    dir.create(file.path(tmp_root, "tracks"), showWarnings = FALSE)
    dir.create(file.path(tmp_root, "seq"), showWarnings = FALSE)
    for (chr in chrom_sizes$chrom) {
        writeLines(
            paste0(rep("A", chrom_size), collapse = ""),
            file.path(tmp_root, "seq", paste0(chr, ".seq"))
        )
    }

    old_gwd <- get("GWD", envir = .misha)
    old_root <- dirname(old_gwd)
    withr::defer(
        {
            if (dir.exists(old_root)) {
                try(gsetroot(old_root), silent = TRUE)
                try(gdb.reload(), silent = TRUE)
            }
        },
        envir = parent.frame()
    )
    withr::local_options(
        list(gmulticontig.indexed_format = indexed),
        .local_envir = parent.frame()
    )

    gdb.init(tmp_root)
    gsetroot(tmp_root)
    gdb.reload()

    tmp_root
}

test_that("gtrack.create_meta on indexed sparse track matches per-chrom result", {
    skip_if_not_installed("withr")

    # --- Per-chrom DB ---
    setup_db(num_chroms = 3, indexed = FALSE)
    intervs <- gintervals(
        chroms = c("chr1", "chr3"),
        starts = c(100, 500),
        ends   = c(200, 600)
    )
    gtrack.create_sparse("sparse_perchrom", "per-chrom sparse", intervs, c(1.5, 2.5))
    withr::defer(gtrack.rm("sparse_perchrom", force = TRUE), envir = parent.frame())

    # Force meta creation by reading .meta via the helper. The helper
    # invokes gtrack_create_meta on first call; subsequent calls hit the
    # cached .meta on disk.
    meta_perchrom <- misha:::.gintervals.big.meta("sparse_perchrom")

    # --- Indexed DB (same data) ---
    # Re-setup as a fresh DB and convert the sparse track to indexed format.
    setup_db(num_chroms = 3, indexed = FALSE)
    intervs2 <- gintervals(
        chroms = c("chr1", "chr3"),
        starts = c(100, 500),
        ends   = c(200, 600)
    )
    gtrack.create_sparse("sparse_indexed", "indexed sparse", intervs2, c(1.5, 2.5))
    withr::defer(gtrack.rm("sparse_indexed", force = TRUE), envir = parent.frame())
    gtrack.convert_to_indexed("sparse_indexed")

    # Sanity: the indexed track really has track.idx.
    idx_path <- file.path(misha:::.track_dir("sparse_indexed"), "track.idx")
    expect_true(file.exists(idx_path))

    meta_indexed <- misha:::.gintervals.big.meta("sparse_indexed")

    # The per-chrom stat rows for chroms with data should match.
    stats_p <- meta_perchrom$stats
    stats_i <- meta_indexed$stats
    expect_equal(nrow(stats_i), nrow(stats_p))

    # The "size" column counts intervals per chrom; on chrs without data it
    # must be 0 in both representations, and on chrs with data they must
    # agree exactly.
    expect_equal(stats_i$size, stats_p$size)
    expect_equal(stats_i$range, stats_p$range)
    expect_equal(stats_i$contains_overlaps, stats_p$contains_overlaps)
})

test_that("gtrack.create_meta scales to indexed track with many empty chroms", {
    skip_if_not_installed("withr")

    # 100 chroms, only 2 with data. Per-chrom path would do ~500 syscalls;
    # indexed path does ~2 contig reads + 1 stat. This is the large-contig
    # case in miniature - we just check correctness here, not timing.
    setup_db(num_chroms = 100, chrom_size = 1e5, indexed = FALSE)

    intervs <- gintervals(
        chroms = c("chr7", "chr42"),
        starts = c(0, 100),
        ends   = c(50, 200)
    )
    gtrack.create_sparse("sparse_sparse_idx", "mostly empty", intervs, c(1, 2))
    withr::defer(gtrack.rm("sparse_sparse_idx", force = TRUE), envir = parent.frame())
    gtrack.convert_to_indexed("sparse_sparse_idx")

    idx_path <- file.path(misha:::.track_dir("sparse_sparse_idx"), "track.idx")
    expect_true(file.exists(idx_path))

    # Triggers meta creation in C++ via the indexed shortcut.
    meta <- misha:::.gintervals.big.meta("sparse_sparse_idx")
    stats <- meta$stats

    # Only chr7 and chr42 should have non-zero size; everyone else 0.
    nonzero_chroms <- as.character(stats$chrom[stats$size > 0])
    expect_setequal(nonzero_chroms, c("chr7", "chr42"))
    expect_equal(sum(stats$size), 2)
})
