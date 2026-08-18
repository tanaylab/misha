# Regression tests for the two defects caused by the deferred whole-genome 2D
# intervals set (.misha$ALLGENOME[[2]]) that gdb.init/gdb.create store above
# `gmulticontig.2d.threshold` contigs:
#
#   * the 2D multitasking partitioner derived the chrom-pair universe from that
#     set, so with the placeholder it planned zero kids and every 2D entry point
#     silently returned NULL - even for a scope the caller passed explicitly;
#   * gtrack.liftover sized its per-chrom-pair buffer vector from that set and
#     indexed it with chromid1 * num_chroms + chromid2, i.e. wrote out of bounds
#     of a zero-length vector.
#
# Both are driven here by lowering the threshold on a scratch database, so they
# need no lab-specific >100-contig database.

flat_seq <- function(len = 2000) paste(rep("A", len), collapse = "")

# A database whose contig count exceeds `gmulticontig.2d.threshold`, so that
# ALLGENOME[[2]] holds the deferred placeholder rather than the N^2 grid.
scratch_db <- function(chrom_names, len = 2000) {
    setup_db(lapply(chrom_names, function(nm) sprintf(">%s\n%s\n", nm, flat_seq(len))),
        return_db = TRUE
    )
}

test_that("2D queries with an explicit 2D scope agree with and without multitasking on a deferred-2D database", {
    local_db_state()
    withr::local_options(list(gmulticontig.2d.threshold = 2))

    chroms <- c("chrMT1", "chrMT2", "chrMT3")
    scratch_db(chroms)

    # The database is above the threshold: the whole-genome 2D set is deferred.
    expect_equal(NROW(.misha$ALLGENOME[[2]]), 0)
    expect_gt(nrow(.misha$ALLGENOME[[1]]), 2)

    pairs <- expand.grid(chrom1 = chroms, chrom2 = chroms, stringsAsFactors = FALSE)
    intervs <- data.frame(
        chrom1 = pairs$chrom1, start1 = 100, end1 = 900,
        chrom2 = pairs$chrom2, start2 = 100, end2 = 900,
        stringsAsFactors = FALSE
    )
    gtrack.2d.create("mt2d", "2D track over every chrom pair", intervs, seq_len(nrow(intervs)))
    withr::defer(if (gtrack.exists("mt2d")) gtrack.rm("mt2d", force = TRUE))

    single <- withr::with_options(list(gmultitasking = FALSE), gextract("mt2d", intervals = "mt2d"))
    multi <- withr::with_options(list(gmultitasking = TRUE), gextract("mt2d", intervals = "mt2d"))

    expect_equal(nrow(single), nrow(intervs))
    expect_false(is.null(multi)) # returned NULL before the fix
    expect_equal(nrow(multi), nrow(single))
    expect_equal(sort(multi$mt2d), sort(single$mt2d))

    summary_single <- withr::with_options(list(gmultitasking = FALSE), gsummary("mt2d", intervals = "mt2d"))
    summary_multi <- withr::with_options(list(gmultitasking = TRUE), gsummary("mt2d", intervals = "mt2d"))
    expect_false(is.null(summary_multi))
    expect_equal(as.numeric(summary_multi), as.numeric(summary_single))
})

test_that("a 2D scope covering only some chrom pairs is partitioned correctly on a deferred-2D database", {
    local_db_state()
    withr::local_options(list(gmulticontig.2d.threshold = 2))

    chroms <- c("chrMS1", "chrMS2", "chrMS3", "chrMS4")
    scratch_db(chroms)
    expect_equal(NROW(.misha$ALLGENOME[[2]]), 0)

    # Only 3 of the 16 chrom pairs are populated, including a trans pair.
    intervs <- data.frame(
        chrom1 = c("chrMS1", "chrMS2", "chrMS4"),
        start1 = c(100, 100, 100), end1 = c(900, 900, 900),
        chrom2 = c("chrMS1", "chrMS3", "chrMS4"),
        start2 = c(100, 100, 100), end2 = c(900, 900, 900),
        stringsAsFactors = FALSE
    )
    gtrack.2d.create("mt2d_sparse", "sparse 2D track", intervs, c(7, 8, 9))
    withr::defer(if (gtrack.exists("mt2d_sparse")) gtrack.rm("mt2d_sparse", force = TRUE))

    multi <- withr::with_options(list(gmultitasking = TRUE), gextract("mt2d_sparse", intervals = "mt2d_sparse"))
    expect_false(is.null(multi))
    expect_equal(nrow(multi), 3)
    expect_equal(sort(multi$mt2d_sparse), c(7, 8, 9))
    # the trans pair must survive the partitioning
    expect_true(any(as.character(multi$chrom1) == "chrMS2" & as.character(multi$chrom2) == "chrMS3"))
})

# ---------------------------------------------------------------------------
# gtrack.liftover of a 2D track into a target database above the threshold
# ---------------------------------------------------------------------------

# Lift a fixed 3-rectangle 2D track from a 2-contig source database into a
# 3-contig target database, and return the lifted rectangles. `threshold`
# decides whether the target's whole-genome 2D set is deferred or materialised.
lift_2d_track <- function(tag, threshold) {
    withr::local_options(list(gmulticontig.2d.threshold = threshold), .local_envir = parent.frame())

    src_chroms <- paste0("chrS", tag, 1:2)
    src_db <- setup_db(lapply(src_chroms, function(nm) sprintf(">%s\n%s\n", nm, flat_seq())), return_db = TRUE)

    src_intervs <- data.frame(
        chrom1 = c(src_chroms[1], src_chroms[1], src_chroms[2]),
        start1 = c(100, 600, 100), end1 = c(300, 800, 300),
        chrom2 = c(src_chroms[1], src_chroms[2], src_chroms[2]),
        start2 = c(100, 600, 100), end2 = c(300, 800, 300),
        stringsAsFactors = FALSE
    )
    src_track <- paste0("src2d_", tag)
    gtrack.2d.create(src_track, "source 2D track", src_intervs, c(10, 20, 30))
    src_track_dir <- file.path(src_db, "tracks", paste0(src_track, ".track"))

    tgt_chroms <- paste0("chrT", tag, 1:3)
    chain <- new_chain_file()
    write_chain_entry(chain, src_chroms[1], 2000, "+", 0, 2000, tgt_chroms[1], 2000, "+", 0, 2000, 1)
    write_chain_entry(chain, src_chroms[2], 2000, "+", 0, 2000, tgt_chroms[2], 2000, "+", 0, 2000, 2)

    scratch_db(tgt_chroms)

    lifted <- paste0("lifted2d_", tag)
    gtrack.liftover(lifted, "lifted 2D track", src_track_dir, chain)

    scope <- gintervals.2d(
        chroms1 = c(tgt_chroms[1], tgt_chroms[1], tgt_chroms[2]), starts1 = 0, ends1 = 2000,
        chroms2 = c(tgt_chroms[1], tgt_chroms[2], tgt_chroms[2]), starts2 = 0, ends2 = 2000
    )
    res <- gextract(lifted, scope)
    res <- res[order(as.character(res$chrom1), as.character(res$chrom2), res$start1, res$start2), ]
    list(res = res, values = res[[lifted]], chroms = tgt_chroms)
}

test_that("gtrack.liftover of a 2D track into a deferred-2D target database writes the right rectangles", {
    local_db_state()

    # Before the fix this wrote through buffered_intervs[chromid1 * num_chroms +
    # chromid2] on a zero-length vector (SIGSEGV in this configuration).
    out <- lift_2d_track("A", threshold = 2)
    expect_equal(NROW(.misha$ALLGENOME[[2]]), 0)

    expect_equal(nrow(out$res), 3)
    expect_equal(out$values, c(10, 20, 30))
    expect_equal(as.character(out$res$chrom1), out$chroms[c(1, 1, 2)])
    expect_equal(as.character(out$res$chrom2), out$chroms[c(1, 2, 2)])
    expect_equal(out$res$start1, c(100, 600, 100))
    expect_equal(out$res$end1, c(300, 800, 300))
    expect_equal(out$res$start2, c(100, 600, 100))
    expect_equal(out$res$end2, c(300, 800, 300))
})

test_that("gtrack.liftover of a 2D track is unchanged when the whole-genome 2D grid is materialised", {
    local_db_state()

    out <- lift_2d_track("B", threshold = 1000)
    expect_equal(NROW(.misha$ALLGENOME[[2]]), 9) # 3 contigs, materialised grid

    expect_equal(nrow(out$res), 3)
    expect_equal(out$values, c(10, 20, 30))
    expect_equal(as.character(out$res$chrom1), out$chroms[c(1, 1, 2)])
    expect_equal(as.character(out$res$chrom2), out$chroms[c(1, 2, 2)])
})
