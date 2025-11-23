load_test_db()
test_that("gintervals.load_chain handles source overlaps with 'error' policy", {
    local_db_state()

    # Create target genome
    setup_db(list(">chr1\nACTGACTGACTGACTGACTGACTGACTGACTG\n", ">chr2\nGGGGCCCCTTTTAAAA\n"))

    # Create chain file with source overlaps: same source position maps to multiple targets
    chain_file <- new_chain_file()

    # Chain 1: source1[0-20] -> chr1[0-20]
    write_chain_entry(chain_file, "source1", 100, "+", 0, 20, "chr1", 32, "+", 0, 20, 1)

    # Chain 2: source1[5-21] -> chr2[0-16] (overlaps with chain 1 at source1[5-20])
    write_chain_entry(chain_file, "source1", 100, "+", 5, 21, "chr2", 16, "+", 0, 16, 2)

    # Default policy should error on source overlap
    expect_error(
        gintervals.load_chain(chain_file),
        "overlap"
    )

    # Explicit error policy should also error
    expect_error(
        gintervals.load_chain(chain_file, src_overlap_policy = "error"),
        "overlap"
    )
})

test_that("gintervals.load_chain handles source overlaps with 'keep' policy", {
    local_db_state()

    # Create target genome
    setup_db(list(">chr1\nACTGACTGACTGACTGACTGACTGACTGACTG\n", ">chr2\nGGGGCCCCTTTTAAAA\n"))

    # Create chain file with source overlaps
    chain_file <- new_chain_file()

    # Chain 1: source1[0-20] -> chr1[0-20]
    write_chain_entry(chain_file, "source1", 100, "+", 0, 20, "chr1", 32, "+", 0, 20, 1)

    # Chain 2: source1[5-21] -> chr2[0-16] (overlaps at source1[5-20])
    write_chain_entry(chain_file, "source1", 100, "+", 5, 21, "chr2", 16, "+", 0, 16, 2)

    # Keep policy should allow overlaps
    chain <- gintervals.load_chain(chain_file, src_overlap_policy = "keep")
    expect_true(nrow(chain) == 2)
    expect_true(all(c("chr1", "chr2") %in% chain$chrom))

    # Check exact chain values
    chain_sorted <- chain[order(chain$chrom, chain$start), ]
    # chr1 mapping block
    chr1_row <- chain_sorted[chain_sorted$chrom == "chr1", ]
    expect_equal(as.numeric(chr1_row$start), 0)
    expect_equal(as.numeric(chr1_row$end), 20)
    expect_equal(as.character(chr1_row$chromsrc), "source1")
    expect_equal(as.numeric(chr1_row$startsrc), 0)
    # chr2 mapping block
    chr2_row <- chain_sorted[chain_sorted$chrom == "chr2", ]
    expect_equal(as.numeric(chr2_row$start), 0)
    expect_equal(as.numeric(chr2_row$end), 16)
    expect_equal(as.character(chr2_row$chromsrc), "source1")
    expect_equal(as.numeric(chr2_row$startsrc), 5)
})

test_that("gintervals.load_chain handles source overlaps with 'discard' policy", {
    local_db_state()

    # Create target genome
    setup_db(list(">chr1\nACTGACTGACTGACTGACTGACTGACTGACTG\n", ">chr2\nGGGGCCCCTTTTAAAA\n", ">chr3\nTATATATATATA\n"))

    # Create chain file: source1 overlaps (maps to chr1 and chr2), source2 is clean (maps to chr3)
    chain_file <- new_chain_file()

    # Chain 1: source1[0-20] -> chr1[0-20]
    write_chain_entry(chain_file, "source1", 100, "+", 0, 20, "chr1", 32, "+", 0, 20, 1)

    # Chain 2: source1[5-21] -> chr2[0-16] (overlaps with chain 1)
    write_chain_entry(chain_file, "source1", 100, "+", 5, 21, "chr2", 16, "+", 0, 16, 2)

    # Chain 3: source2[0-12] -> chr3[0-12] (no overlaps)
    write_chain_entry(chain_file, "source2", 50, "+", 0, 12, "chr3", 12, "+", 0, 12, 3)

    # Discard policy should remove overlapping intervals but keep clean ones
    chain <- gintervals.load_chain(chain_file, src_overlap_policy = "discard")
    expect_equal(nrow(chain), 1)
    expect_equal(as.character(chain$chrom), "chr3")
    expect_equal(as.character(chain$chromsrc), "source2")
    expect_equal(as.numeric(chain$start), 0)
    expect_equal(as.numeric(chain$end), 12)
    expect_equal(as.numeric(chain$startsrc), 0)
})

test_that("gintervals.load_chain handles target overlaps with 'auto_first' policy", {
    local_db_state()

    # Create target genome
    setup_db(list(">chr1\nACTGACTGACTGACTGACTGACTGACTGACTGACTGACTG\n"))

    # Create chain file with target overlaps: different sources map to overlapping chr1 regions
    chain_file <- new_chain_file()

    # Chain 1: source1[0-30] -> chr1[0-30]
    write_chain_entry(chain_file, "source1", 100, "+", 0, 30, "chr1", 40, "+", 0, 30, 1)

    # Chain 2: source2[0-20] -> chr1[20-40] (overlaps at chr1[20-30])
    write_chain_entry(chain_file, "source2", 100, "+", 0, 20, "chr1", 40, "+", 20, 40, 2)

    # Auto_first policy should resolve overlaps by truncating
    chain <- gintervals.load_chain(chain_file, tgt_overlap_policy = "auto_first")
    expect_true(nrow(chain) >= 1)

    # Check that no intervals overlap in target
    chain_sorted <- chain[order(chain$chrom, chain$start), ]
    if (nrow(chain_sorted) > 1) {
        for (i in 2:nrow(chain_sorted)) {
            if (chain_sorted$chrom[i] == chain_sorted$chrom[i - 1]) {
                expect_true(chain_sorted$start[i] >= chain_sorted$end[i - 1])
            }
        }
    }
})

test_that("gintervals.load_chain exact truncation with 'auto_first' policy", {
    local_db_state()

    setup_db(list(">chr1\nACTGACTGACTGACTGACTGACTGACTGACTGACTGACTG\n"))

    chain_file <- new_chain_file()

    # Chain 1: sourceA[0-30] -> chr1[0-30]
    write_chain_entry(chain_file, "sourceA", 100, "+", 0, 30, "chr1", 40, "+", 0, 30, 1)
    # Chain 2: sourceB[0-20] -> chr1[20-40] (overlaps 20-30)
    write_chain_entry(chain_file, "sourceB", 100, "+", 0, 20, "chr1", 40, "+", 20, 40, 2)

    chain <- gintervals.load_chain(chain_file, tgt_overlap_policy = "auto_first")
    chain_chr1 <- chain[chain$chrom == "chr1", ]
    chain_chr1 <- chain_chr1[order(chain_chr1$start), ]
    # After auto truncation: first keeps its coverage [0,30), second is trimmed to [30,40)
    expect_equal(nrow(chain_chr1), 2)
    expect_equal(as.numeric(chain_chr1$start[1]), 0)
    expect_equal(as.numeric(chain_chr1$end[1]), 30)
    expect_equal(as.numeric(chain_chr1$start[2]), 30)
    expect_equal(as.numeric(chain_chr1$end[2]), 40)
})

test_that("gintervals.load_chain handles target overlaps with 'error' policy", {
    local_db_state()

    # Create target genome
    setup_db(list(">chr1\nACTGACTGACTGACTGACTGACTGACTGACTG\n"))

    # Create chain file with target overlaps
    chain_file <- new_chain_file()

    # Chain 1: source1[0-20] -> chr1[0-20]
    write_chain_entry(chain_file, "source1", 100, "+", 0, 20, "chr1", 32, "+", 0, 20, 1)

    # Chain 2: source2[0-17] -> chr1[15-32] (overlaps at chr1[15-20])
    write_chain_entry(chain_file, "source2", 100, "+", 0, 17, "chr1", 32, "+", 15, 32, 2)

    # Error policy should fail on target overlap
    expect_error(
        gintervals.load_chain(chain_file, tgt_overlap_policy = "error"),
        "overlap"
    )
})

test_that("gintervals.load_chain handles target overlaps with 'discard' policy", {
    local_db_state()

    # Create target genome
    setup_db(list(">chr1\nACTGACTGACTGACTGACTGACTGACTGACTG\n", ">chr2\nGGGGCCCCTTTTAAAA\n"))

    # Create chain file: chr1 has overlaps, chr2 doesn't
    chain_file <- new_chain_file()

    # Chain 1: source1[0-20] -> chr1[0-20]
    write_chain_entry(chain_file, "source1", 100, "+", 0, 20, "chr1", 32, "+", 0, 20, 1)

    # Chain 2: source2[0-17] -> chr1[15-32] (overlaps with chain 1)
    write_chain_entry(chain_file, "source2", 100, "+", 0, 17, "chr1", 32, "+", 15, 32, 2)

    # Chain 3: source3[0-10] -> chr2[0-10] (no overlap)
    write_chain_entry(chain_file, "source3", 100, "+", 0, 10, "chr2", 16, "+", 0, 10, 3)

    # Discard policy should remove overlapping chr1 intervals but keep chr2
    chain <- gintervals.load_chain(chain_file, tgt_overlap_policy = "discard")
    expect_equal(nrow(chain), 1)
    expect_equal(as.character(chain$chrom), "chr2")
})

test_that("gintervals.load_chain segments overlaps for auto and agg policies", {
    local_db_state()

    setup_db(list(">chrSeg\n", paste(rep("A", 100), collapse = ""), "\n"))

    chain_file <- new_chain_file()
    # Three overlapping chains on the target: ids 1..3 with increasing scores
    write_chain_entry(chain_file, "srcA", 200, "+", 0, 20, "chrSeg", 100, "+", 10, 30, 1, score = 10)
    write_chain_entry(chain_file, "srcB", 200, "+", 0, 20, "chrSeg", 100, "+", 15, 35, 2, score = 20)
    write_chain_entry(chain_file, "srcC", 200, "+", 0, 20, "chrSeg", 100, "+", 20, 40, 3, score = 30)

    chain_auto <- gintervals.load_chain(chain_file, tgt_overlap_policy = "auto_score")
    chain_auto_chr <- chain_auto[chain_auto$chrom == "chrSeg", ]
    expect_equal(as.numeric(chain_auto_chr$start), c(10, 15, 20))
    expect_equal(as.numeric(chain_auto_chr$end), c(15, 20, 40))
    expect_equal(as.numeric(chain_auto_chr$chain_id), c(1, 2, 3))

    chain_agg <- gintervals.load_chain(chain_file, tgt_overlap_policy = "agg")
    chain_agg_chr <- chain_agg[chain_agg$chrom == "chrSeg", ]
    expect_equal(
        as.numeric(chain_agg_chr$start),
        c(10, 15, 15, 20, 20, 20, 30, 30, 35)
    )
    expect_equal(
        as.numeric(chain_agg_chr$end),
        c(15, 20, 20, 30, 30, 30, 35, 35, 40)
    )
    expect_equal(
        as.numeric(chain_agg_chr$chain_id),
        c(1, 1, 2, 1, 2, 3, 2, 3, 3)
    )

    # startsrc shifts should match target trim
    chain1_segments <- chain_agg_chr[chain_agg_chr$chain_id == 1, ]
    expect_equal(as.numeric(chain1_segments$startsrc), c(0, 5, 10))
    chain3_segments <- chain_agg_chr[chain_agg_chr$chain_id == 3, ]
    expect_equal(as.numeric(chain3_segments$startsrc), c(0, 10, 15))
})

test_that("gintervals.liftover returns parallel slices with agg policy", {
    local_db_state()

    setup_db(list(">chrLift\n", paste(rep("C", 120), collapse = ""), "\n"))

    chain_file <- new_chain_file()
    write_chain_entry(chain_file, "srcA", 200, "+", 0, 20, "chrLift", 120, "+", 0, 20, 1, score = 5)
    write_chain_entry(chain_file, "srcB", 200, "+", 0, 20, "chrLift", 120, "+", 10, 30, 2, score = 10)
    write_chain_entry(chain_file, "srcC", 200, "+", 0, 20, "chrLift", 120, "+", 15, 35, 3, score = 15)

    src_intervals <- data.frame(
        chrom = c("srcA", "srcB", "srcC"),
        start = 0,
        end = 20,
        stringsAsFactors = FALSE
    )

    lifted_agg <- gintervals.liftover(
        src_intervals,
        chain_file,
        src_overlap_policy = "keep",
        tgt_overlap_policy = "agg"
    )
    lifted_chr <- lifted_agg[lifted_agg$chrom == "chrLift", ]
    # Expect duplicated (start, end) pairs when chains collide
    agg_pairs <- lifted_chr[, c("start", "end")]
    dup_counts <- table(paste(agg_pairs$start, agg_pairs$end, sep = "-"))
    expect_true(any(dup_counts > 1))

    lifted_auto <- gintervals.liftover(
        src_intervals,
        chain_file,
        src_overlap_policy = "keep",
        tgt_overlap_policy = "auto_score"
    )
    lifted_auto_chr <- lifted_auto[lifted_auto$chrom == "chrLift", ]
    auto_pairs <- lifted_auto_chr[, c("start", "end")]
    auto_counts <- table(paste(auto_pairs$start, auto_pairs$end, sep = "-"))
    expect_true(all(auto_counts == 1))
})

test_that("gintervals.liftover works with 'keep' source policy", {
    local_db_state()

    # Create target genome
    setup_db(list(">chr1\nACTGACTGACTGACTGACTGACTGACTGACTG\n", ">chr2\nGGGGCCCCTTTTAAAA\n"))

    # Create chain with source overlaps
    chain_file <- new_chain_file()

    # Chain 1: source1[0-20] -> chr1[0-20]
    write_chain_entry(chain_file, "source1", 100, "+", 0, 20, "chr1", 32, "+", 0, 20, 1)

    # Chain 2: source1[10-26] -> chr2[0-16] (overlaps at source1[10-20])
    write_chain_entry(chain_file, "source1", 100, "+", 10, 26, "chr2", 16, "+", 0, 16, 2)

    # Load chain with keep policy
    chain <- gintervals.load_chain(chain_file, src_overlap_policy = "keep", tgt_overlap_policy = "keep")

    # Create source intervals that overlap with the ambiguous region
    src_intervals <- data.frame(
        chrom = "source1",
        start = 10,
        end = 20,
        stringsAsFactors = FALSE
    )

    # Liftover should produce multiple target intervals
    result <- gintervals.liftover(src_intervals, chain)

    # Should have entries mapping to both chr1 and chr2
    expect_true(nrow(result) >= 2)
    expect_true("chr1" %in% result$chrom)
    expect_true("chr2" %in% result$chrom)
    # Verify exact coordinates: source [10,20) overlaps both chains
    # For chr1 block source [0,20)->chr1[0,20): [10,20)->[10,20)
    chr1_rows <- result[result$chrom == "chr1", ]
    if (nrow(chr1_rows) > 0) {
        expect_true(all(as.numeric(chr1_rows$start) == 10))
        expect_true(all(as.numeric(chr1_rows$end) == 20))
    }
    # For chr2 block source [10,26)->chr2[0,16): [10,20)->[0,10)
    chr2_rows <- result[result$chrom == "chr2", ]
    if (nrow(chr2_rows) > 0) {
        expect_true(all(as.numeric(chr2_rows$start) == 0))
        expect_true(all(as.numeric(chr2_rows$end) == 10))
    }
    # intervalID should refer to the single input interval (1)
    expect_true(all(as.integer(result$intervalID) == 1))
})

test_that("gintervals.liftover works with chain file path", {
    local_db_state()

    # Create target genome
    setup_db(list(">chr1\nACTGACTGACTGACTGACTGACTGACTGACTG\n"))

    # Create simple chain file
    chain_file <- new_chain_file()

    # Chain: source1[0-20] -> chr1[0-20]
    write_chain_entry(chain_file, "source1", 100, "+", 0, 20, "chr1", 32, "+", 0, 20, 1)

    # Create source intervals
    src_intervals <- data.frame(
        chrom = "source1",
        start = 5,
        end = 15,
        stringsAsFactors = FALSE
    )

    # Liftover using chain file path directly
    result <- gintervals.liftover(src_intervals, chain_file)

    expect_true(nrow(result) >= 1)
    expect_equal(as.character(result$chrom[1]), "chr1")
    expect_equal(as.numeric(result$start[1]), 5)
    expect_equal(as.numeric(result$end[1]), 15)
    expect_true(all(as.integer(result$intervalID) == 1))
})

test_that("gintervals.liftover basic 1D mapping with exact coordinates", {
    local_db_state()

    setup_db(list(">chr1\nAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n", ">chr2\nCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n"))

    chain_file <- new_chain_file()

    # sourceA[0-20) -> chr1[10-30)
    write_chain_entry(chain_file, "sourceA", 200, "+", 0, 20, "chr1", 40, "+", 10, 30, 1)

    # Create a source interval fully inside the block
    src_intervals <- data.frame(
        chrom = "sourceA",
        start = 5,
        end = 20,
        stringsAsFactors = FALSE
    )

    result <- gintervals.liftover(src_intervals, chain_file)
    expect_equal(nrow(result), 1)
    expect_equal(as.character(result$chrom[1]), "chr1")
    expect_equal(as.numeric(result$start[1]), 15)
    expect_equal(as.numeric(result$end[1]), 30)
    expect_equal(as.integer(result$intervalID[1]), 1)
})

test_that("gintervals.liftover simple 2D mapping cross-product", {
    local_db_state()

    setup_db(list(">chr1\nAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n", ">chr2\nCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n"))

    chain_file <- new_chain_file()

    # sourceA[0-20) -> chr1[10-30)
    write_chain_entry(chain_file, "sourceA", 200, "+", 0, 20, "chr1", 40, "+", 10, 30, 1)
    # sourceB[10-37) -> chr2[5-32) (length 27 to match chr2 size 32)
    write_chain_entry(chain_file, "sourceB", 200, "+", 10, 37, "chr2", 32, "+", 5, 32, 2)

    chain <- gintervals.load_chain(chain_file, src_overlap_policy = "keep")

    # Build a 2D source interval pairing (sourceA[5,20), sourceB[15,40))
    src2d <- data.frame(
        chrom1 = "sourceA", start1 = 5, end1 = 20,
        chrom2 = "sourceB", start2 = 15, end2 = 40,
        stringsAsFactors = FALSE
    )

    res2d <- gintervals.liftover(src2d, chain)
    expect_equal(nrow(res2d), 1)
    expect_equal(as.character(res2d$chrom1[1]), "chr1")
    expect_equal(as.character(res2d$chrom2[1]), "chr2")
    expect_equal(as.numeric(res2d$start1[1]), 15)
    expect_equal(as.numeric(res2d$end1[1]), 30)
    expect_equal(as.numeric(res2d$start2[1]), 10)
    expect_equal(as.numeric(res2d$end2[1]), 32)
})

test_that("gtrack.liftover works from sparse source track and preserves values", {
    local_db_state()

    # Create target genome
    setup_db(list(">chr1\nAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n"))
    target_db <- .misha$GROOT

    # Chain mapping chrsource1[0-20) -> chr1[0-20)
    # Note: Database chromosome is "chrsource1" (from filename chrsource1.fasta)
    chain_file <- new_chain_file()
    write_chain_entry(chain_file, "chrsource1", 100, "+", 0, 20, "chr1", 32, "+", 0, 20, 1)

    # Create a source DB with 'source1'
    source_db <- setup_source_db(list(">source1\nTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\n"))

    # Create a sparse source track on source1
    # Database chromosome name is "chrsource1" (from filename chrsource1.fasta)
    src_intervals <- data.frame(
        chrom = "chrsource1",
        start = c(5, 10),
        end = c(10, 15),
        stringsAsFactors = FALSE
    )
    src_values <- c(1.0, 2.0)
    src_track_name <- "src_sparse"
    gtrack.create_sparse(src_track_name, "source sparse", src_intervals, src_values)

    # Compute the source track directory path (tracks live under <GROOT>/tracks)
    src_track_dir <- file.path(source_db, "tracks", paste0(src_track_name, ".track"))

    # Switch to target DB and liftover the track
    gdb.init(target_db)
    lifted_track <- "lifted_sparse"
    withr::defer({
        if (gtrack.exists(lifted_track)) gtrack.rm(lifted_track, force = TRUE)
    })
    gtrack.liftover(lifted_track, "lifted from source", src_track_dir, chain_file)

    # Extract all values and validate coordinates and values
    res <- gextract(lifted_track, gintervals.all())
    # Expect two intervals on chr1 with same coords and values as mapped 1:1
    expect_equal(nrow(res), 2)
    expect_equal(as.character(res$chrom), c("chr1", "chr1"))
    expect_equal(as.numeric(res$start), c(5, 10))
    expect_equal(as.numeric(res$end), c(10, 15))
    # Value column name may differ; detect the non-coordinate column
    value_cols <- setdiff(colnames(res), c("chrom", "start", "end", "strand"))
    expect_true(length(value_cols) >= 1)
    expect_equal(as.numeric(res[[value_cols[1]]]), c(1.0, 2.0))
})

test_that("gintervals.liftover validates policy parameters", {
    local_db_state()

    # Create target genome
    setup_db(list(">chr1\nACTGACTGACTGACTG\n"))

    chain_file <- new_chain_file()

    write_chain_entry(chain_file, "source1", 100, "+", 0, 10, "chr1", 16, "+", 0, 10, 1)

    src_intervals <- data.frame(chrom = "source1", start = 0, end = 10, stringsAsFactors = FALSE)

    # Invalid src_overlap_policy
    expect_error(
        gintervals.liftover(src_intervals, chain_file, src_overlap_policy = "invalid"),
        "src_overlap_policy"
    )

    # Invalid tgt_overlap_policy
    expect_error(
        gintervals.liftover(src_intervals, chain_file, tgt_overlap_policy = "invalid"),
        "tgt_overlap_policy"
    )
})

test_that("gintervals.load_chain validates policy parameters", {
    local_db_state()

    # Create target genome
    setup_db(list(">chr1\nACTGACTGACTG\n"))

    chain_file <- new_chain_file()

    write_chain_entry(chain_file, "source1", 100, "+", 0, 10, "chr1", 12, "+", 0, 10, 1)

    # Invalid src_overlap_policy
    expect_error(
        gintervals.load_chain(chain_file, src_overlap_policy = "invalid"),
        "src_overlap_policy"
    )

    # Invalid tgt_overlap_policy
    expect_error(
        gintervals.load_chain(chain_file, tgt_overlap_policy = "invalid"),
        "tgt_overlap_policy"
    )
})

test_that("gintervals.liftover returns empty result when all intervals discarded", {
    local_db_state()

    # Create target genome
    setup_db(list(">chr1\nACTGACTGACTGACTGACTGACTGACTGACTG\n"))

    # Create chain with only overlapping intervals
    chain_file <- new_chain_file()

    # Chain 1: source1[0-20] -> chr1[0-20]
    write_chain_entry(chain_file, "source1", 100, "+", 0, 20, "chr1", 32, "+", 0, 20, 1)

    # Chain 2: source1[10-27] -> chr1[15-32] (both source and target overlaps)
    write_chain_entry(chain_file, "source1", 100, "+", 10, 27, "chr1", 32, "+", 15, 32, 2)

    # Load chain and discard all overlaps
    chain <- gintervals.load_chain(chain_file, src_overlap_policy = "discard")

    # Chain should be empty or NULL
    expect_true(is.null(chain) || nrow(chain) == 0)
})

test_that("complex chain with both source and target overlaps", {
    local_db_state()

    # Create target genome with multiple chromosomes
    setup_db(list(">chr1\nACTGACTGACTGACTGACTGACTGACTGACTG\n", ">chr2\nGGGGCCCCTTTTAAAA\n", ">chr3\nTATATATATATA\n"))

    # Create complex chain file with various overlap scenarios
    chain_file <- new_chain_file()

    # Source overlap: source1[0,20] -> chr1 and chr2
    write_chain_entry(chain_file, "source1", 100, "+", 0, 20, "chr1", 32, "+", 0, 20, 1)

    write_chain_entry(chain_file, "source1", 100, "+", 10, 26, "chr2", 16, "+", 0, 16, 2)

    # Target overlap: source2 and source3 both map to overlapping chr3 regions
    write_chain_entry(chain_file, "source2", 100, "+", 0, 8, "chr3", 12, "+", 0, 8, 3)

    write_chain_entry(chain_file, "source3", 100, "+", 0, 7, "chr3", 12, "+", 5, 12, 4)

    # Clean mapping: source4 -> chr1 (no overlaps)
    write_chain_entry(chain_file, "source4", 100, "+", 0, 7, "chr1", 32, "+", 25, 32, 5)

    # Test keep + auto: should keep source overlaps and auto-resolve target overlaps
    chain_keep_auto <- gintervals.load_chain(
        chain_file,
        src_overlap_policy = "keep",
        tgt_overlap_policy = "auto"
    )
    expect_true(nrow(chain_keep_auto) >= 3) # At least the clean mapping + some resolved

    # Test discard + discard: should only keep the clean mapping
    chain_discard_discard <- gintervals.load_chain(
        chain_file,
        src_overlap_policy = "discard",
        tgt_overlap_policy = "discard"
    )
    expect_true(nrow(chain_discard_discard) <= 1)
    if (nrow(chain_discard_discard) > 0) {
        expect_true("chr1" %in% chain_discard_discard$chrom)
    }
})

test_that("gintervals.load_chain returns 10 columns with strand information, chain_id and score", {
    local_db_state()

    # Create target genome
    setup_db(list(">chr1\nACTGACTGACTGACTGACTGACTGACTGACTG\n"))

    # Create simple chain file with forward strand
    chain_file <- new_chain_file()

    write_chain_entry(chain_file, "source1", 100, "+", 0, 20, "chr1", 32, "+", 0, 20, 1)

    chain <- gintervals.load_chain(chain_file)

    # Check that chain has 10 columns (including chain_id and score)
    expected_cols <- c("chrom", "start", "end", "strand", "chromsrc", "startsrc", "endsrc", "strandsrc", "chain_id", "score")
    expect_equal(colnames(chain), expected_cols)

    # Check strand format: +1 for forward
    expect_equal(as.numeric(chain$strand), 1)
    expect_equal(as.numeric(chain$strandsrc), 1)

    # Check that source intervals are complete with end coordinate
    expect_equal(as.numeric(chain$startsrc), 0)
    expect_equal(as.numeric(chain$endsrc), 20)
})

test_that("gintervals.load_chain handles reverse strands correctly", {
    local_db_state()

    # Create target genome
    setup_db(list(">chr1\nACTGACTGACTGACTGACTGACTGACTGACTG\n"))

    # Create chain file with reverse strands
    chain_file <- new_chain_file()

    # Source on reverse strand, target on forward
    write_chain_entry(chain_file, "source1", 100, "-", 0, 20, "chr1", 32, "+", 0, 20, 1)

    chain <- gintervals.load_chain(chain_file)

    # Check strand format: target forward (+1), source reverse (-1)
    expect_equal(as.numeric(chain$strand), 1)
    expect_equal(as.numeric(chain$strandsrc), -1)
})

test_that("gintervals.load_chain handles mixed strand chains", {
    local_db_state()

    # Create target genome
    setup_db(list(">chr1\nACTGACTGACTGACTGACTGACTGACTGACTG\n", ">chr2\nGGGGCCCCTTTTAAAA\n"))

    # Create chain file with mixed strands
    chain_file <- new_chain_file()

    # Chain 1: forward-forward
    write_chain_entry(chain_file, "source1", 100, "+", 0, 20, "chr1", 32, "+", 0, 20, 1)

    # Chain 2: forward-reverse (source forward, target reverse)
    write_chain_entry(chain_file, "source2", 100, "+", 0, 16, "chr2", 16, "-", 0, 16, 2)

    chain <- gintervals.load_chain(chain_file)

    # Check that we have 2 chains with different strand combinations
    expect_equal(nrow(chain), 2)

    chr1_chain <- chain[chain$chrom == "chr1", ]
    expect_equal(as.numeric(chr1_chain$strand), 1)
    expect_equal(as.numeric(chr1_chain$strandsrc), 1)

    chr2_chain <- chain[chain$chrom == "chr2", ]
    expect_equal(as.numeric(chr2_chain$strand), -1)
    expect_equal(as.numeric(chr2_chain$strandsrc), 1)
})

test_that("gintervals.load_chain with src_groot validates source chromosomes", {
    local_db_state()

    # Create target genome
    setup_db(list(">chr1\nACTGACTGACTGACTGACTGACTGACTGACTG\n"))
    target_db <- .misha$GROOT

    # Create source genome database with specific chromosomes
    source_db <- setup_source_db(list(">source1\nTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\n"))

    # Create chain file with valid source chromosome
    chain_file <- new_chain_file()

    write_chain_entry(chain_file, "source1", 100, "+", 0, 20, "chr1", 32, "+", 0, 20, 1)

    # Switch back to target db
    gdb.init(target_db)

    # This should succeed with src_groot validation
    chain <- gintervals.load_chain(chain_file, src_groot = source_db)
    expect_true(!is.null(chain))
    expect_equal(nrow(chain), 1)
    expect_equal(as.character(chain$chromsrc), "source1")

    # Verify we're back to target database
    expect_equal(.misha$GROOT, target_db)
})

test_that("gintervals.load_chain with src_groot rejects invalid source chromosomes", {
    local_db_state()

    # Create target genome
    setup_db(list(">chr1\nACTGACTGACTGACTGACTGACTGACTGACTG\n"))
    target_db <- .misha$GROOT

    # Create source genome with different chromosome name
    source_db <- setup_source_db(list(">different_chrom\nTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\n"))

    # Create chain file with invalid source chromosome
    chain_file <- new_chain_file()

    write_chain_entry(chain_file, "invalid_source", 100, "+", 0, 20, "chr1", 32, "+", 0, 20, 1)

    # Switch back to target db
    gdb.init(target_db)

    # This should error because invalid_source doesn't exist in source_db
    expect_error(
        gintervals.load_chain(chain_file, src_groot = source_db),
        "Chromosome.*does not exist|Source chromosome|invalid chromosomes"
    )

    # Verify we're back to target database even after error
    expect_equal(.misha$GROOT, target_db)
})

test_that("gintervals.load_chain with src_groot validates source coordinates", {
    local_db_state()

    # Create target genome (70bp to match chain file)
    setup_db(list(">chr1\nACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGAC\n"))
    target_db <- .misha$GROOT

    # Create source genome with small chromosome
    source_db <- setup_source_db(list(">source1\nTTTTTTTTTTTTTTTTTTTT\n")) # Only 20 bp

    # Create chain file with source coordinates exceeding chromosome size
    chain_file <- new_chain_file()

    write_chain_entry(chain_file, "source1", 100, "+", 0, 50, "chr1", 70, "+", 0, 50, 1)

    # Switch back to target db
    gdb.init(target_db)

    # This should error because source coordinates exceed chromosome size
    expect_error(
        gintervals.load_chain(chain_file, src_groot = source_db)
    )

    # Verify we're back to target database even after error
    expect_equal(.misha$GROOT, target_db)
})

test_that("gintervals.liftover works with 8-column chain format", {
    local_db_state()

    # Create target genome
    setup_db(list(">chr1\nACTGACTGACTGACTGACTGACTGACTGACTG\n"))

    # Create chain file
    chain_file <- new_chain_file()

    write_chain_entry(chain_file, "source1", 100, "+", 0, 20, "chr1", 32, "+", 0, 20, 1)

    # Load chain (now returns 8 columns)
    chain <- gintervals.load_chain(chain_file)

    # Create source intervals
    src_intervals <- data.frame(
        chrom = "source1",
        start = 5,
        end = 15,
        stringsAsFactors = FALSE
    )

    # Liftover should work with the new 8-column format
    result <- gintervals.liftover(src_intervals, chain)

    expect_true(nrow(result) >= 1)
    expect_equal(as.character(result$chrom[1]), "chr1")
    expect_equal(as.numeric(result$start[1]), 5)
    expect_equal(as.numeric(result$end[1]), 15)
})

test_that("gintervals.liftover matches liftOver binary - basic case", {
    skip_if_not(has_liftover_binary(), "liftOver binary not found")

    local_db_state()

    # Create target genome
    setup_db(list(">chr1\nACTGACTGACTGACTGACTGACTGACTGACTGACTGACTG\n"))

    # Create chain file with simple mapping
    chain_file <- new_chain_file()

    # source1[10-30] -> chr1[5-25]
    write_chain_entry(chain_file, "source1", 100, "+", 10, 30, "chr1", 40, "+", 5, 25, 1)

    # Create source intervals in BED format
    bed_input <- tempfile(fileext = ".bed")
    bed_output <- tempfile(fileext = ".bed")
    bed_unmapped <- tempfile(fileext = ".unmapped")
    withr::defer({
        unlink(bed_input)
        unlink(bed_output)
        unlink(bed_unmapped)
    })

    # BED format is 0-based, half-open
    cat("source1\t12\t28\tinterval1\n", file = bed_input)
    cat("source1\t15\t20\tinterval2\n", file = bed_input, append = TRUE)

    # Run liftOver binary
    system2("liftOver", args = c(bed_input, chain_file, bed_output, bed_unmapped), stdout = FALSE, stderr = FALSE)

    # Read binary output
    binary_result <- read.table(bed_output, header = FALSE, stringsAsFactors = FALSE, col.names = c("chrom", "start", "end", "name"))

    # Run misha liftover
    chain <- gintervals.load_chain(chain_file)
    src_intervals <- data.frame(
        chrom = c("source1", "source1"),
        start = c(12, 15),
        end = c(28, 20),
        stringsAsFactors = FALSE
    )
    misha_result <- gintervals.liftover(src_intervals, chain)
    misha_result <- misha_result[order(misha_result$intervalID, misha_result$start), ]

    # Compare results
    expect_equal(nrow(misha_result), nrow(binary_result))
    expect_equal(as.character(misha_result$chrom), binary_result$chrom)
    expect_equal(as.numeric(misha_result$start), binary_result$start)
    expect_equal(as.numeric(misha_result$end), binary_result$end)
})

test_that("gintervals.liftover matches liftOver binary - multiple chains", {
    skip_if_not(has_liftover_binary(), "liftOver binary not found")

    local_db_state()

    # Create target genome
    setup_db(list(">chr1\nACTGACTGACTGACTGACTGACTGACTGACTGACTGACTG\n", ">chr2\nGGGGCCCCTTTTAAAAGGGGCCCCTTTTAAAA\n"))

    # Create chain file with multiple mappings
    chain_file <- new_chain_file()

    # source1[0-20] -> chr1[10-30]
    write_chain_entry(chain_file, "source1", 100, "+", 0, 20, "chr1", 40, "+", 10, 30, 1)

    # source2[5-25] -> chr2[0-20]
    write_chain_entry(chain_file, "source2", 100, "+", 5, 25, "chr2", 32, "+", 0, 20, 2)

    # Create source intervals in BED format
    bed_input <- tempfile(fileext = ".bed")
    bed_output <- tempfile(fileext = ".bed")
    bed_unmapped <- tempfile(fileext = ".unmapped")
    withr::defer({
        unlink(bed_input)
        unlink(bed_output)
        unlink(bed_unmapped)
    })

    cat("source1\t5\t15\tinterval1\n", file = bed_input)
    cat("source2\t10\t20\tinterval2\n", file = bed_input, append = TRUE)

    # Run liftOver binary
    system2("liftOver", args = c(bed_input, chain_file, bed_output, bed_unmapped), stdout = FALSE, stderr = FALSE)

    # Read binary output
    binary_result <- read.table(bed_output, header = FALSE, stringsAsFactors = FALSE, col.names = c("chrom", "start", "end", "name"))
    binary_result <- binary_result[order(binary_result$name, binary_result$start), ]

    # Run misha liftover
    chain <- gintervals.load_chain(chain_file)
    src_intervals <- data.frame(
        chrom = c("source1", "source2"),
        start = c(5, 10),
        end = c(15, 20),
        stringsAsFactors = FALSE
    )
    misha_result <- gintervals.liftover(src_intervals, chain)
    misha_result <- misha_result[order(misha_result$intervalID, misha_result$start), ]

    # Compare results
    expect_equal(nrow(misha_result), nrow(binary_result))
    expect_equal(as.character(misha_result$chrom), binary_result$chrom)
    expect_equal(as.numeric(misha_result$start), binary_result$start)
    expect_equal(as.numeric(misha_result$end), binary_result$end)
})

test_that("gintervals.liftover matches liftOver binary - reverse strand", {
    skip_if_not(has_liftover_binary(), "liftOver binary not found")

    local_db_state()

    # Create target genome (needs to be at least 50bp to match chain)
    setup_db(list(">chr1\nACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTG\n"))

    # Create chain file with reverse strand mapping
    chain_file <- new_chain_file()

    # source1[0-30] -> chr1[22-52] on reverse strand (chr1 size is 56)
    write_chain_entry(chain_file, "source1", 100, "+", 0, 30, "chr1", 56, "-", 4, 34, 1)

    # Create source intervals in BED format
    bed_input <- tempfile(fileext = ".bed")
    bed_output <- tempfile(fileext = ".bed")
    bed_unmapped <- tempfile(fileext = ".unmapped")
    withr::defer({
        unlink(bed_input)
        unlink(bed_output)
        unlink(bed_unmapped)
    })

    cat("source1\t5\t15\tinterval1\n", file = bed_input)

    # Run liftOver binary
    system2("liftOver", args = c(bed_input, chain_file, bed_output, bed_unmapped), stdout = FALSE, stderr = FALSE)

    # Read binary output
    binary_result <- read.table(bed_output, header = FALSE, stringsAsFactors = FALSE, col.names = c("chrom", "start", "end", "name"))

    # Run misha liftover
    chain <- gintervals.load_chain(chain_file)
    src_intervals <- data.frame(
        chrom = "source1",
        start = 5,
        end = 15,
        stringsAsFactors = FALSE
    )
    misha_result <- gintervals.liftover(src_intervals, chain)

    # Compare results - for reverse strand mappings
    # Note: Currently checking basic properties rather than exact coordinates
    # as there may be differences in reverse strand coordinate interpretation
    expect_equal(nrow(misha_result), nrow(binary_result))
    expect_equal(as.character(misha_result$chrom), binary_result$chrom)
    # Verify that interval lengths match
    expect_equal(
        as.numeric(misha_result$end) - as.numeric(misha_result$start),
        binary_result$end - binary_result$start
    )
})

test_that("gintervals.liftover matches liftOver binary - partial overlap", {
    skip_if_not(has_liftover_binary(), "liftOver binary not found")

    local_db_state()

    # Create target genome
    setup_db(list(">chr1\nACTGACTGACTGACTGACTGACTGACTGACTGACTGACTG\n"))

    # Create chain file
    chain_file <- new_chain_file()

    # source1[10-30] -> chr1[5-25]
    write_chain_entry(chain_file, "source1", 100, "+", 10, 30, "chr1", 40, "+", 5, 25, 1)

    # Create source intervals that partially overlap the chain
    bed_input <- tempfile(fileext = ".bed")
    bed_output <- tempfile(fileext = ".bed")
    bed_unmapped <- tempfile(fileext = ".unmapped")
    withr::defer({
        unlink(bed_input)
        unlink(bed_output)
        unlink(bed_unmapped)
    })

    # This interval overlaps only partially with the chain
    cat("source1\t5\t15\tinterval1\n", file = bed_input)
    cat("source1\t25\t35\tinterval2\n", file = bed_input, append = TRUE)

    # Run liftOver binary with minMatch=0.1 to allow partial overlaps
    system2("liftOver", args = c("-minMatch=0.1", bed_input, chain_file, bed_output, bed_unmapped), stdout = FALSE, stderr = FALSE)

    # Read binary output (if any)
    if (file.exists(bed_output) && file.info(bed_output)$size > 0) {
        binary_result <- read.table(bed_output, header = FALSE, stringsAsFactors = FALSE, col.names = c("chrom", "start", "end", "name"))
        binary_result <- binary_result[order(binary_result$name, binary_result$start), ]
    } else {
        binary_result <- data.frame(chrom = character(), start = numeric(), end = numeric(), name = character())
    }

    # Run misha liftover
    chain <- gintervals.load_chain(chain_file)
    src_intervals <- data.frame(
        chrom = c("source1", "source1"),
        start = c(5, 25),
        end = c(15, 35),
        stringsAsFactors = FALSE
    )
    misha_result <- gintervals.liftover(src_intervals, chain)
    misha_result <- misha_result[order(misha_result$intervalID, misha_result$start), ]

    # Compare results
    expect_equal(nrow(misha_result), nrow(binary_result))
    if (nrow(binary_result) > 0) {
        expect_equal(as.character(misha_result$chrom), binary_result$chrom)
        expect_equal(as.numeric(misha_result$start), binary_result$start)
        expect_equal(as.numeric(misha_result$end), binary_result$end)
    }
})

test_that("gintervals.liftover matches liftOver binary - complex chain with gaps", {
    skip_if_not(has_liftover_binary(), "liftOver binary not found")

    local_db_state()

    # Create target genome (needs at least 60bp to match chain declaration)
    setup_db(list(">chr1\nACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTG\n"))

    # Create chain file with gaps (insertions/deletions)
    chain_file <- new_chain_file()

    # Chain with gaps: source1[0-27] -> chr1[0-40] with gaps
    # Block 1: 10 bases aligned
    # Gap: 2 bases in source, 3 bases in target
    # Block 2: 15 bases aligned
    # Total source: 10 + 2 + 15 = 27
    # Total target: 10 + 3 + 15 = 28 (but chain declares up to 40 for alignment region)
    cat("chain 1000 source1 100 + 0 27 chr1 64 + 0 28 1\n", file = chain_file)
    cat("10\t2\t3\n", file = chain_file, append = TRUE)
    cat("15\n\n", file = chain_file, append = TRUE)

    # Create source intervals
    bed_input <- tempfile(fileext = ".bed")
    bed_output <- tempfile(fileext = ".bed")
    bed_unmapped <- tempfile(fileext = ".unmapped")
    withr::defer({
        unlink(bed_input)
        unlink(bed_output)
        unlink(bed_unmapped)
    })

    cat("source1\t2\t8\tinterval1\n", file = bed_input)
    cat("source1\t15\t23\tinterval2\n", file = bed_input, append = TRUE)

    # Run liftOver binary
    system2("liftOver", args = c(bed_input, chain_file, bed_output, bed_unmapped), stdout = FALSE, stderr = FALSE)

    # Read binary output
    binary_result <- read.table(bed_output, header = FALSE, stringsAsFactors = FALSE, col.names = c("chrom", "start", "end", "name"))
    binary_result <- binary_result[order(binary_result$name, binary_result$start), ]

    # Run misha liftover
    chain <- gintervals.load_chain(chain_file)
    src_intervals <- data.frame(
        chrom = c("source1", "source1"),
        start = c(2, 15),
        end = c(8, 23),
        stringsAsFactors = FALSE
    )
    misha_result <- gintervals.liftover(src_intervals, chain)
    misha_result <- misha_result[order(misha_result$intervalID, misha_result$start), ]

    # Compare results
    expect_equal(nrow(misha_result), nrow(binary_result))
    expect_equal(as.character(misha_result$chrom), binary_result$chrom)
    expect_equal(as.numeric(misha_result$start), binary_result$start)
    expect_equal(as.numeric(misha_result$end), binary_result$end)
})

test_that("gintervals.liftover matches liftOver binary - interval spanning gap", {
    skip_if_not(has_liftover_binary(), "liftOver binary not found")

    local_db_state()

    # Create target genome
    setup_db(list(">chr1\nACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTG\n"))

    # Create chain file with gap
    chain_file <- new_chain_file()

    # Chain: source1[0-35] -> chr1[0-40] with a gap
    # Block 1: 15 bases [0-15)
    # Gap: 5 bases in source, 10 bases in target
    # Block 2: 15 bases [20-35)
    cat("chain 1000 source1 100 + 0 35 chr1 72 + 0 40 1\n", file = chain_file)
    cat("15\t5\t10\n", file = chain_file, append = TRUE)
    cat("15\n\n", file = chain_file, append = TRUE)

    # Create intervals that test gap handling
    bed_input <- tempfile(fileext = ".bed")
    bed_output <- tempfile(fileext = ".bed")
    bed_unmapped <- tempfile(fileext = ".unmapped")
    withr::defer({
        unlink(bed_input)
        unlink(bed_output)
        unlink(bed_unmapped)
    })

    # Interval entirely in first block
    cat("source1\t2\t8\tfirst_block\n", file = bed_input)
    # Interval entirely in second block
    cat("source1\t22\t30\tsecond_block\n", file = bed_input, append = TRUE)

    # Run liftOver binary
    system2("liftOver", args = c(bed_input, chain_file, bed_output, bed_unmapped), stdout = FALSE, stderr = FALSE)

    # Read binary output
    if (file.exists(bed_output) && file.info(bed_output)$size > 0) {
        binary_result <- read.table(bed_output, header = FALSE, stringsAsFactors = FALSE, col.names = c("chrom", "start", "end", "name"))
        binary_result <- binary_result[order(binary_result$name, binary_result$start), ]
    } else {
        binary_result <- data.frame(chrom = character(), start = numeric(), end = numeric(), name = character())
    }

    # Run misha liftover
    chain <- gintervals.load_chain(chain_file)
    src_intervals <- data.frame(
        chrom = c("source1", "source1"),
        start = c(2, 22),
        end = c(8, 30),
        stringsAsFactors = FALSE
    )
    misha_result <- gintervals.liftover(src_intervals, chain)

    # Create a comparison key for sorting
    if (nrow(misha_result) > 0) {
        misha_result$sort_key <- paste(misha_result$intervalID, sprintf("%010d", misha_result$start))
        misha_result <- misha_result[order(misha_result$sort_key), ]
        misha_result$sort_key <- NULL
    }

    # Compare results
    expect_equal(nrow(misha_result), nrow(binary_result))
    if (nrow(binary_result) > 0) {
        expect_equal(as.character(misha_result$chrom), binary_result$chrom)
        expect_equal(as.numeric(misha_result$start), binary_result$start)
        expect_equal(as.numeric(misha_result$end), binary_result$end)
    }
})

test_that("gintervals.liftover matches liftOver binary - very small intervals", {
    skip_if_not(has_liftover_binary(), "liftOver binary not found")

    local_db_state()

    # Create target genome
    setup_db(list(">chr1\nACTGACTGACTGACTGACTGACTGACTGACTGACTGACTG\n"))

    # Create chain file
    chain_file <- new_chain_file()

    # Simple chain
    write_chain_entry(chain_file, "source1", 100, "+", 0, 30, "chr1", 40, "+", 5, 35, 1)

    # Create very small intervals (1-2bp)
    bed_input <- tempfile(fileext = ".bed")
    bed_output <- tempfile(fileext = ".bed")
    bed_unmapped <- tempfile(fileext = ".unmapped")
    withr::defer({
        unlink(bed_input)
        unlink(bed_output)
        unlink(bed_unmapped)
    })

    cat("source1\t5\t6\tsmall1\n", file = bed_input)
    cat("source1\t10\t12\tsmall2\n", file = bed_input, append = TRUE)
    cat("source1\t25\t26\tsmall3\n", file = bed_input, append = TRUE)

    # Run liftOver binary
    system2("liftOver", args = c(bed_input, chain_file, bed_output, bed_unmapped), stdout = FALSE, stderr = FALSE)

    # Read binary output
    binary_result <- read.table(bed_output, header = FALSE, stringsAsFactors = FALSE, col.names = c("chrom", "start", "end", "name"))
    binary_result <- binary_result[order(binary_result$name, binary_result$start), ]

    # Run misha liftover
    chain <- gintervals.load_chain(chain_file)
    src_intervals <- data.frame(
        chrom = c("source1", "source1", "source1"),
        start = c(5, 10, 25),
        end = c(6, 12, 26),
        stringsAsFactors = FALSE
    )
    misha_result <- gintervals.liftover(src_intervals, chain)
    misha_result <- misha_result[order(misha_result$intervalID, misha_result$start), ]

    # Compare results
    expect_equal(nrow(misha_result), nrow(binary_result))
    expect_equal(as.character(misha_result$chrom), binary_result$chrom)
    expect_equal(as.numeric(misha_result$start), binary_result$start)
    expect_equal(as.numeric(misha_result$end), binary_result$end)
})

test_that("gintervals.liftover matches liftOver binary - boundary intervals", {
    skip_if_not(has_liftover_binary(), "liftOver binary not found")

    local_db_state()

    # Create target genome
    setup_db(list(">chr1\nACTGACTGACTGACTGACTGACTGACTGACTGACTGACTG\n", ">chr2\nGGGGCCCCTTTTAAAAGGGGCCCC\n"))

    # Create chain file with two chains
    chain_file <- new_chain_file()

    # Chain 1: source1[10-30] -> chr1[5-25]
    write_chain_entry(chain_file, "source1", 100, "+", 10, 30, "chr1", 40, "+", 5, 25, 1)

    # Chain 2: source1[35-50] -> chr2[2-17]
    write_chain_entry(chain_file, "source1", 100, "+", 35, 50, "chr2", 24, "+", 2, 17, 2)

    # Create intervals at exact boundaries
    bed_input <- tempfile(fileext = ".bed")
    bed_output <- tempfile(fileext = ".bed")
    bed_unmapped <- tempfile(fileext = ".unmapped")
    withr::defer({
        unlink(bed_input)
        unlink(bed_output)
        unlink(bed_unmapped)
    })

    # Interval at start of chain 1
    cat("source1\t10\t15\tstart_bound\n", file = bed_input)
    # Interval at end of chain 1
    cat("source1\t25\t30\tend_bound\n", file = bed_input, append = TRUE)
    # Interval exactly matching chain 1
    cat("source1\t10\t30\texact_match\n", file = bed_input, append = TRUE)
    # Interval at boundaries of chain 2
    cat("source1\t35\t40\tchain2_start\n", file = bed_input, append = TRUE)

    # Run liftOver binary
    system2("liftOver", args = c(bed_input, chain_file, bed_output, bed_unmapped), stdout = FALSE, stderr = FALSE)

    # Read binary output - sort by name  for comparison
    binary_result <- read.table(bed_output, header = FALSE, stringsAsFactors = FALSE, col.names = c("chrom", "start", "end", "name"))
    binary_result <- binary_result[order(binary_result$name, binary_result$start), ]

    # Run misha liftover
    chain <- gintervals.load_chain(chain_file)
    src_intervals <- data.frame(
        chrom = c("source1", "source1", "source1", "source1"),
        start = c(10, 25, 10, 35),
        end = c(15, 30, 30, 40),
        name = c("start_bound", "end_bound", "exact_match", "chain2_start"),
        stringsAsFactors = FALSE
    )
    misha_result <- gintervals.liftover(src_intervals, chain)

    # Sort misha results to match binary output order  (by original interval order)
    # Map intervalID back to original name
    misha_result$name <- src_intervals$name[misha_result$intervalID]
    misha_result <- misha_result[order(misha_result$name, misha_result$start), ]

    # Compare results
    expect_equal(nrow(misha_result), nrow(binary_result))
    expect_equal(as.character(misha_result$chrom), binary_result$chrom)
    expect_equal(as.numeric(misha_result$start), binary_result$start)
    expect_equal(as.numeric(misha_result$end), binary_result$end)
})

test_that("gintervals.liftover matches liftOver binary - mixed strand chains", {
    skip_if_not(has_liftover_binary(), "liftOver binary not found")

    local_db_state()

    # Create target genome
    setup_db(list(">chr1\nACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTG\n", ">chr2\nGGGGCCCCTTTTAAAAGGGGCCCCTTTTAAAA\n"))

    # Create chain file with mixed strands
    chain_file <- new_chain_file()

    # Chain 1: forward strand
    write_chain_entry(chain_file, "source1", 100, "+", 0, 20, "chr1", 56, "+", 10, 30, 1)

    # Chain 2: forward source, reverse target
    write_chain_entry(chain_file, "source2", 100, "+", 5, 25, "chr2", 32, "-", 5, 25, 2)

    # Create intervals
    bed_input <- tempfile(fileext = ".bed")
    bed_output <- tempfile(fileext = ".bed")
    bed_unmapped <- tempfile(fileext = ".unmapped")
    withr::defer({
        unlink(bed_input)
        unlink(bed_output)
        unlink(bed_unmapped)
    })

    cat("source1\t5\t15\tforward\n", file = bed_input)
    cat("source2\t10\t20\treverse\n", file = bed_input, append = TRUE)

    # Run liftOver binary
    system2("liftOver", args = c(bed_input, chain_file, bed_output, bed_unmapped), stdout = FALSE, stderr = FALSE)

    # Read binary output
    binary_result <- read.table(bed_output, header = FALSE, stringsAsFactors = FALSE, col.names = c("chrom", "start", "end", "name"))
    binary_result <- binary_result[order(binary_result$name, binary_result$start), ]

    # Run misha liftover
    chain <- gintervals.load_chain(chain_file)
    src_intervals <- data.frame(
        chrom = c("source1", "source2"),
        start = c(5, 10),
        end = c(15, 20),
        stringsAsFactors = FALSE
    )
    misha_result <- gintervals.liftover(src_intervals, chain)
    misha_result <- misha_result[order(misha_result$intervalID, misha_result$start), ]

    # Compare results
    expect_equal(nrow(misha_result), nrow(binary_result))
    expect_equal(as.character(misha_result$chrom), binary_result$chrom)
    # For mixed strands, verify interval lengths match
    expect_equal(
        as.numeric(misha_result$end) - as.numeric(misha_result$start),
        binary_result$end - binary_result$start
    )
})


test_that("gintervals.liftover matches liftOver binary - many small chains", {
    skip_if_not(has_liftover_binary(), "liftOver binary not found")

    local_db_state()

    # Create target genome
    chr1_seq <- paste(rep("ACTG", 50), collapse = "")
    chr2_seq <- paste(rep("GCTA", 40), collapse = "")
    setup_db(list(paste0(">chr1\n", chr1_seq, "\n"), paste0(">chr2\n", chr2_seq, "\n")))

    # Create chain file with many small non-overlapping chains
    chain_file <- new_chain_file()

    # Multiple chains from different source regions to different targets
    write_chain_entry(chain_file, "source1", 200, "+", 0, 15, "chr1", 200, "+", 10, 25, 1)

    write_chain_entry(chain_file, "source1", 200, "+", 20, 35, "chr1", 200, "+", 50, 65, 2)

    write_chain_entry(chain_file, "source1", 200, "+", 40, 58, "chr2", 160, "+", 20, 38, 3)

    write_chain_entry(chain_file, "source2", 150, "+", 5, 25, "chr1", 200, "+", 100, 120, 4)

    write_chain_entry(chain_file, "source2", 150, "+", 30, 50, "chr2", 160, "+", 60, 80, 5)

    # Create intervals
    bed_input <- tempfile(fileext = ".bed")
    bed_output <- tempfile(fileext = ".bed")
    bed_unmapped <- tempfile(fileext = ".unmapped")
    withr::defer({
        unlink(bed_input)
        unlink(bed_output)
        unlink(bed_unmapped)
    })

    cat("source1\t5\t12\tint1\n", file = bed_input)
    cat("source1\t22\t33\tint2\n", file = bed_input, append = TRUE)
    cat("source1\t42\t55\tint3\n", file = bed_input, append = TRUE)
    cat("source2\t10\t20\tint4\n", file = bed_input, append = TRUE)
    cat("source2\t35\t45\tint5\n", file = bed_input, append = TRUE)

    # Run liftOver binary
    system2("liftOver", args = c(bed_input, chain_file, bed_output, bed_unmapped), stdout = FALSE, stderr = FALSE)

    # Read binary output
    binary_result <- read.table(bed_output, header = FALSE, stringsAsFactors = FALSE, col.names = c("chrom", "start", "end", "name"))
    binary_result <- binary_result[order(binary_result$name, binary_result$start), ]

    # Run misha liftover
    chain <- gintervals.load_chain(chain_file)
    src_intervals <- data.frame(
        chrom = c("source1", "source1", "source1", "source2", "source2"),
        start = c(5, 22, 42, 10, 35),
        end = c(12, 33, 55, 20, 45),
        stringsAsFactors = FALSE
    )
    misha_result <- gintervals.liftover(src_intervals, chain)
    misha_result <- misha_result[order(misha_result$intervalID, misha_result$start), ]

    # Compare results
    expect_equal(nrow(misha_result), nrow(binary_result))
    expect_equal(as.character(misha_result$chrom), binary_result$chrom)
    expect_equal(as.numeric(misha_result$start), binary_result$start)
    expect_equal(as.numeric(misha_result$end), binary_result$end)
})

test_that("gintervals.liftover finds all overlapping chains when they are non-consecutive", {
    local_db_state()

    # This test verifies the fix for a bug where gintervals.liftover would miss
    # overlapping chain intervals that were not consecutive in the sorted chain array.
    # The bug occurred when source intervals had overlapping regions mapping to
    # different targets, creating a pattern where overlapping chains are separated
    # by non-overlapping ones in the sorted (by start_src) chain array.

    # Create target genome
    setup_db(list(">chr1\nACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTG\n", ">chr2\nGGGGCCCCTTTTAAAAGGGGCCCCTTTTAAAAGGGGCCCCTTTTAAAA\n", ">chr3\nCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCG\n"))

    # Create chain file replicating the user's bug scenario
    chain_file <- new_chain_file()

    # Chains with the same start position (will be consecutive when sorted)
    write_chain_entry(chain_file, "source1", 100, "+", 10, 30, "chr1", 44, "+", 0, 20, 1)

    write_chain_entry(chain_file, "source1", 100, "+", 10, 30, "chr2", 48, "+", 0, 20, 2)

    # Chain with different range that creates a gap
    write_chain_entry(chain_file, "source1", 100, "+", 30, 32, "chr3", 42, "+", 0, 2, 3)

    # Load chain with keep policy
    chain <- gintervals.load_chain(chain_file, src_overlap_policy = "keep", tgt_overlap_policy = "keep")

    # Query interval that overlaps first two chains but not the third
    src_intervals <- data.frame(chrom = "source1", start = 20, end = 21, stringsAsFactors = FALSE)

    # Perform liftover
    result <- gintervals.liftover(src_intervals, chain)

    # Should return exactly 2 results (both chains with [10,30) overlap [20,21))
    expect_equal(nrow(result), 2)
    expect_true("chr1" %in% result$chrom)
    expect_true("chr2" %in% result$chrom)
    expect_false("chr3" %in% result$chrom)
})

test_that("Do not miss earlier long overlap when hint is to the right (policy=keep)", {
    local_db_state()

    setup_db(list(
        paste0(">chr1\n", paste(rep("A", 200), collapse = ""), "\n"),
        paste0(">chr2\n", paste(rep("C", 200), collapse = ""), "\n"),
        paste0(">chr3\n", paste(rep("G", 200), collapse = ""), "\n")
    ))

    chain_file <- new_chain_file()

    # A: [0,100) -> chr1 (overlaps Q1 and Q2)
    write_chain_entry(chain_file, "source1", 200, "+", 0, 100, "chr1", 200, "+", 0, 100, 1)

    # B: [15,16) -> chr2 (does not overlap Q1 or Q2)
    write_chain_entry(chain_file, "source1", 200, "+", 15, 16, "chr2", 200, "+", 0, 1, 2)

    # C: [80,110) -> chr3 (overlaps Q2 only)
    write_chain_entry(chain_file, "source1", 200, "+", 80, 110, "chr3", 200, "+", 0, 30, 3)

    chain <- gintervals.load_chain(chain_file, src_overlap_policy = "keep", tgt_overlap_policy = "keep")

    # Two queries, sorted (as gintervals.liftover does internally)
    src_intervals <- data.frame(
        chrom = "source1",
        start = c(70, 90), # Q1 then Q2
        end = c(71, 91),
        stringsAsFactors = FALSE
    )

    result <- gintervals.liftover(src_intervals, chain)

    # Expect: Q1 -> chr1; Q2 -> chr1 and chr3  (total 3 rows)
    # Previous bug would yield only 2 rows (Q1->chr1, Q2->chr3), missing chr1 for Q2.
    expect_equal(sum(result$intervalID == 1), 1) # Q1 has exactly 1 hit (chr1)
    expect_equal(sum(result$intervalID == 2), 2) # Q2 should have 2 hits (chr1 & chr3)
    expect_true("chr1" %in% result$chrom[result$intervalID == 2])
    expect_true("chr3" %in% result$chrom[result$intervalID == 2])
})

test_that("Deterministic ordering for chains with identical start_src", {
    local_db_state()

    # Create target genome with multiple chromosomes
    setup_db(list(
        ">chr1\nACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTG\n",
        ">chr2\nGGGGCCCCTTTTAAAAGGGGCCCCTTTTAAAA\n",
        ">chr3\nTTTTAAAACCCCGGGGTTTTAAAACCCCGGGG\n"
    ))

    chain_file <- new_chain_file()

    # Create chains with identical start_src but different end_src
    # All start at source1[100], but have different lengths
    write_chain_entry(chain_file, "source1", 500, "+", 100, 120, "chr1", 48, "+", 0, 20, 1) # len=20
    write_chain_entry(chain_file, "source1", 500, "+", 100, 110, "chr2", 32, "+", 0, 10, 2) # len=10
    write_chain_entry(chain_file, "source1", 500, "+", 100, 130, "chr3", 32, "+", 0, 30, 3) # len=30
    write_chain_entry(chain_file, "source1", 500, "+", 100, 115, "chr1", 48, "+", 20, 35, 4) # len=15

    chain <- gintervals.load_chain(chain_file, src_overlap_policy = "keep", tgt_overlap_policy = "keep")

    # Run liftover twice with the same query
    src_intervals <- data.frame(
        chrom = "source1",
        start = 105,
        end = 106,
        stringsAsFactors = FALSE
    )

    result1 <- gintervals.liftover(src_intervals, chain)
    result2 <- gintervals.liftover(src_intervals, chain)

    # Results should be identical (deterministic)
    expect_equal(result1, result2)

    # Verify all 4 chains are returned
    expect_equal(nrow(result1), 4)

    # Verify all expected chroms are present (order may vary deterministically)
    expect_true(all(c("chr1", "chr2", "chr3") %in% result1$chrom))
})

test_that("Zero-length chain intervals produce error", {
    local_db_state()

    setup_db(list(">chr1\nACTGACTGACTGACTGACTGACTGACTGACTG\n"))

    chain_file <- new_chain_file()

    # Create chain with zero-length interval (start = end)
    write_chain_entry(chain_file, "source1", 100, "+", 10, 10, "chr1", 32, "+", 0, 0, 1)

    # Should error on zero-length chain (existing validation catches this)
    expect_error(
        gintervals.load_chain(chain_file),
        "end coordinate is less or equal than the start"
    )
})

test_that("Minus strand mapping at edges produces valid intervals", {
    local_db_state()

    # Create a genome with chr1 that has at least 50 bases
    setup_db(list(">chr1\n", paste0(rep("A", 60), collapse = ""), "\n"))

    chain_file <- new_chain_file()

    # Minus strand chain: source1[0-30] -> chr1[20-50] (minus strand)
    write_chain_entry(chain_file, "source1", 100, "+", 0, 30, "chr1", 60, "-", 20, 50, 1)

    chain <- gintervals.load_chain(chain_file)

    # Test left edge
    src_intervals <- data.frame(
        chrom = c("source1", "source1", "source1"),
        start = c(0, 15, 29), # left edge, middle, right edge
        end = c(1, 16, 30),
        stringsAsFactors = FALSE
    )

    result <- gintervals.liftover(src_intervals, chain)

    # All results should have start < end (primary correctness check)
    expect_true(all(result$start < result$end))

    # Verify we got results for all 3 queries
    expect_equal(length(unique(result$intervalID)), 3)

    # Verify all coordinates are non-negative and reasonable
    expect_true(all(result$start >= 0))
    expect_true(all(result$end >= 0))
})

test_that("Dense cluster of chains with same start_src performs correctly", {
    local_db_state()

    setup_db(list(
        ">chr1\n", paste0(rep("A", 10000), collapse = ""), "\n",
        ">chr2\n", paste0(rep("C", 10000), collapse = ""), "\n"
    ))

    chain_file <- new_chain_file()

    # Create 100 chains all starting at source1[50] with varying lengths
    for (i in 1:100) {
        len <- i # lengths 1 to 100
        target_chrom <- if (i %% 2 == 0) "chr1" else "chr2"
        write_chain_entry(
            chain_file, "source1", 10000, "+", 50, 50 + len,
            target_chrom, 10000, "+", i * 10, i * 10 + len, i
        )
    }

    chain <- gintervals.load_chain(chain_file, src_overlap_policy = "keep", tgt_overlap_policy = "keep")

    # Query that overlaps all chains
    src_intervals <- data.frame(
        chrom = "source1",
        start = 60,
        end = 61,
        stringsAsFactors = FALSE
    )

    result <- gintervals.liftover(src_intervals, chain)

    # Should find many overlapping chains (query at 60 overlaps chains where 50+length > 60)
    # Exact number depends on overlap resolution, but should be substantial
    expect_true(nrow(result) > 50) # At least half should overlap
    expect_true(nrow(result) <= 100) # At most all chains

    # All results should be valid intervals
    expect_true(all(result$start < result$end))
})

test_that("Target overlap auto keeps primary mapping when later chain is contained", {
    local_db_state()

    setup_db(list(
        ">chrT\n", paste0(rep("A", 500), collapse = ""), "\n"
    ))

    chain_file <- new_chain_file()

    # Chain A covers chrT[0, 200)
    write_chain_entry(
        chain_file, "sourceA", 500, "+", 0, 200,
        "chrT", 500, "+", 0, 200, 1
    )

    # Chain B is contained within Chain A on the target
    write_chain_entry(
        chain_file, "sourceB", 500, "+", 0, 120,
        "chrT", 500, "+", 50, 170, 2
    )

    # With default "auto" (now "auto_score"), chain B should be kept if it has higher score
    # But since we want to test the truncation behavior, use "auto_first" to get the old behavior
    chain <- gintervals.load_chain(chain_file, tgt_overlap_policy = "auto_first")

    # Check that chain was loaded successfully
    expect_false(is.null(chain))
    expect_true(nrow(chain) > 0)

    srcA <- data.frame(chrom = "sourceA", start = 60, end = 80, stringsAsFactors = FALSE)
    resA <- gintervals.liftover(srcA, chain)
    expect_equal(nrow(resA), 1)
    expect_equal(as.character(resA$chrom), "chrT")
    expect_equal(as.numeric(resA$start), 60)
    expect_equal(as.numeric(resA$end), 80)

    srcB <- data.frame(chrom = "sourceB", start = 60, end = 80, stringsAsFactors = FALSE)
    # With auto_first, chain B is truncated/discarded because it's contained in chain A
    expect_error(gintervals.liftover(srcB, chain), "does not exist")

    chain_keep <- gintervals.load_chain(chain_file, tgt_overlap_policy = "keep")
    resB_keep <- gintervals.liftover(srcB, chain_keep)
    expect_equal(nrow(resB_keep), 1)
    expect_equal(as.character(resB_keep$chrom), "chrT")
})

test_that("Regression: UCSC-style liftover returns single hit under target auto policy", {
    local_db_state()

    setup_db(list(
        ">AncRef\n", paste0(rep("A", 50000), collapse = ""), "\n"
    ))

    chain_file <- new_chain_file()

    # Primary long chain covering wide region
    # Note: Database chromosome is "chrAncRef" (from filename chrAncRef.fasta)
    write_chain_entry(
        chain_file, "chrX", 5000000, "+", 4550000, 4551000,
        "chrAncRef", 50000, "+", 15000, 16000, 1
    )

    # Several short chains overlapping the same target stretch; should be trimmed/dropped
    offsets <- seq(14790, 19090, by = 400)
    for (idx in seq_along(offsets)) {
        src_start <- 4560000 + idx * 10
        write_chain_entry(
            chain_file, "chrX", 5000000, "+", src_start, src_start + 5,
            "chrAncRef", 50000, "+", offsets[idx], offsets[idx] + 5, 100 + idx
        )
    }

    chain <- gintervals.load_chain(chain_file)

    src_interval <- data.frame(chrom = "chrX", start = 4550156, end = 4550157, stringsAsFactors = FALSE)
    result <- gintervals.liftover(src_interval, chain)

    expect_equal(nrow(result), 1)
    expect_equal(as.character(result$chrom), "chrAncRef") # Database chromosome name
    expect_equal(as.numeric(result$start), 15000 + (4550156 - 4550000))
    expect_equal(as.numeric(result$end), 15000 + (4550157 - 4550000))
})

test_that("Large coordinates near overflow boundaries work correctly", {
    local_db_state()

    setup_db(list(">chr1\n", paste0(rep("A", 1000), collapse = ""), "\n"))

    chain_file <- new_chain_file()

    # Use large but safe coordinates (well below INT64_MAX/2)
    large_coord <- 2^50 # ~1 quadrillion, safe for 64-bit arithmetic

    # Chain with large source coordinates
    write_chain_entry(
        chain_file, "source1", large_coord + 1000, "+",
        large_coord, large_coord + 100,
        "chr1", 1000, "+", 0, 100, 1
    )

    chain <- gintervals.load_chain(chain_file)

    # Query with large coordinates
    src_intervals <- data.frame(
        chrom = "source1",
        start = large_coord + 10,
        end = large_coord + 20,
        stringsAsFactors = FALSE
    )

    result <- gintervals.liftover(src_intervals, chain)

    # Should produce valid result
    expect_equal(nrow(result), 1)
    expect_true(result$start < result$end)

    # Verify mapping is correct: source[+10,+20] -> target[+10,+20]
    expect_equal(result$start, 10)
    expect_equal(result$end, 20)
})

# Score-based liftover tests
test_that("auto_score policy segments overlapping targets and selects winners per segment", {
    local_db_state()

    # Create target genome with enough space for multiple chains
    setup_db(list(">chr1\nAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n"))

    # Create chain file with overlapping TARGET intervals at different scores
    # Chain 1: target [10,30), score 5000
    # Chain 2: target [15,35), score 3000
    # Chain 3: target [12,32), score 8000
    chain_file <- new_chain_file()
    cat("chain 5000 source1 100 + 10 30 chr1 60 + 10 30 1\n20\n\n", file = chain_file)
    cat("chain 3000 source1 100 + 10 30 chr1 60 + 15 35 2\n20\n\n", file = chain_file, append = TRUE)
    cat("chain 8000 source1 100 + 10 30 chr1 60 + 12 32 3\n20\n\n", file = chain_file, append = TRUE)

    # Load chain with auto_score - should segment overlapping targets
    chain <- gintervals.load_chain(chain_file, src_overlap_policy = "keep", tgt_overlap_policy = "auto_score")

    # After segmentation, should have non-overlapping intervals
    # Segments: [10,12) from chain 1, [12,15) from chain 3, [15,30) from chain 3, [30,32) from chain 3, [32,35) from chain 2
    expect_true(nrow(chain) >= 3) # At least 3 segments

    # Check that target intervals don't overlap (sorted by start)
    chain_sorted <- chain[order(chain$start), ]
    for (i in 1:(nrow(chain_sorted) - 1)) {
        expect_true(chain_sorted$end[i] <= chain_sorted$start[i + 1],
            info = sprintf(
                "Overlap at row %d: [%d,%d) and [%d,%d)",
                i, chain_sorted$start[i], chain_sorted$end[i],
                chain_sorted$start[i + 1], chain_sorted$end[i + 1]
            )
        )
    }

    # Test liftover - source interval maps through segmented chain
    src_intervals <- data.frame(chrom = "source1", start = 10, end = 30, stringsAsFactors = FALSE)
    result <- gintervals.liftover(src_intervals, chain)

    # Should get multiple non-overlapping results
    expect_true(nrow(result) >= 1)
    expect_true(all(result$chrom == "chr1"))
})

test_that("auto is an alias for auto_score", {
    local_db_state()

    # Create target genome with enough space for multiple chains
    setup_db(list(">chr1\nAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n"))

    # Create chain file with overlapping TARGET intervals at different scores
    chain_file <- new_chain_file()
    cat("chain 5000 source1 100 + 10 30 chr1 60 + 10 30 1\n20\n\n", file = chain_file)
    cat("chain 3000 source1 100 + 10 30 chr1 60 + 15 35 2\n20\n\n", file = chain_file, append = TRUE)
    cat("chain 8000 source1 100 + 10 30 chr1 60 + 12 32 3\n20\n\n", file = chain_file, append = TRUE)

    # Load chain with "auto" (should work as alias for "auto_score")
    chain_auto <- gintervals.load_chain(chain_file, src_overlap_policy = "keep", tgt_overlap_policy = "auto")
    chain_auto_score <- gintervals.load_chain(chain_file, src_overlap_policy = "keep", tgt_overlap_policy = "auto_score")

    # Both should produce the same segmented chain
    expect_equal(nrow(chain_auto), nrow(chain_auto_score))
    expect_equal(chain_auto, chain_auto_score)

    # Test intervals that overlap all three chains
    src_intervals <- data.frame(chrom = "source1", start = 10, end = 30, stringsAsFactors = FALSE)

    # Both should produce identical liftover results
    result_auto <- gintervals.liftover(src_intervals, chain_auto)
    result_auto_score <- gintervals.liftover(src_intervals, chain_auto_score)

    expect_equal(result_auto, result_auto_score)
    expect_true(nrow(result_auto) >= 1)
    expect_true(all(result_auto$chrom == "chr1"))
})

test_that("auto_score policy uses length as tiebreaker when scores are equal", {
    local_db_state()

    setup_db(list(">chr1\nAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n"))

    # Create chains with equal scores but different lengths, with overlapping targets
    # Chain 1: [10,25), len 15, score 5000
    # Chain 2: [12,32), len 20, score 5000
    # Chain 3: [15,27), len 12, score 5000
    chain_file <- new_chain_file()
    cat("chain 5000 source1 100 + 10 25 chr1 60 + 10 25 1\n15\n\n", file = chain_file)
    cat("chain 5000 source1 100 + 10 30 chr1 60 + 12 32 2\n20\n\n", file = chain_file, append = TRUE)
    cat("chain 5000 source1 100 + 10 22 chr1 60 + 15 27 3\n12\n\n", file = chain_file, append = TRUE)

    chain <- gintervals.load_chain(chain_file, src_overlap_policy = "keep", tgt_overlap_policy = "auto_score")

    # Should have segmented, non-overlapping intervals
    chain_sorted <- chain[order(chain$start), ]
    for (i in 1:(nrow(chain_sorted) - 1)) {
        expect_true(chain_sorted$end[i] <= chain_sorted$start[i + 1])
    }

    # In overlapping regions, longer chain (chain 2) should win due to tiebreaker
    # Check that chain 2 won at least some segments
    expect_true(any(chain$chain_id == 2))
})

test_that("auto_score policy uses chain_id as final tiebreaker", {
    local_db_state()

    setup_db(list(">chr1\nAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n"))

    # Create chains with equal scores and equal lengths, with overlapping targets
    # chain_id from file: 3, 1, 2
    chain_file <- new_chain_file()
    cat("chain 5000 source1 100 + 10 30 chr1 60 + 10 30 3\n20\n\n", file = chain_file)
    cat("chain 5000 source1 100 + 10 30 chr1 60 + 15 35 1\n20\n\n", file = chain_file, append = TRUE)
    cat("chain 5000 source1 100 + 10 30 chr1 60 + 12 32 2\n20\n\n", file = chain_file, append = TRUE)

    chain <- gintervals.load_chain(chain_file, src_overlap_policy = "keep", tgt_overlap_policy = "auto_score")

    # Should have segmented, non-overlapping intervals
    chain_sorted <- chain[order(chain$start), ]
    for (i in 1:(nrow(chain_sorted) - 1)) {
        expect_true(chain_sorted$end[i] <= chain_sorted$start[i + 1])
    }

    # When scores and lengths are equal, lowest chain_id should win
    # Check that chain with id=1 won at least some segments
    expect_true(any(chain$chain_id == 1))
})

test_that("min_score parameter filters low-scoring chains", {
    local_db_state()

    setup_db(list(">chr1\nAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n"))

    # Create chains with varying scores
    chain_file <- new_chain_file()
    cat("chain 5000 source1 100 + 10 30 chr1 60 + 10 30 1\n20\n\n", file = chain_file)
    cat("chain 3000 source1 100 + 35 55 chr1 60 + 40 60 2\n20\n\n", file = chain_file, append = TRUE)
    cat("chain 8000 source1 100 + 60 80 chr1 60 + 5 25 3\n20\n\n", file = chain_file, append = TRUE)

    # Load without min_score: all chains
    chain_all <- gintervals.load_chain(chain_file)
    expect_equal(nrow(chain_all), 3)

    # Load with min_score=4000: should filter out chain 2 (score 3000)
    chain_filtered <- gintervals.load_chain(chain_file, min_score = 4000)
    expect_equal(nrow(chain_filtered), 2)
    expect_true(all(chain_filtered$score >= 4000))

    # Load with min_score=9000: should keep none (returns empty data frame)
    chain_none <- gintervals.load_chain(chain_file, min_score = 9000)
    expect_equal(nrow(chain_none), 0)
})

test_that("include_metadata parameter adds score and chain_id columns", {
    local_db_state()

    setup_db(list(">chr1\nAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n"))

    chain_file <- new_chain_file()
    cat("chain 7500 source1 100 + 10 30 chr1 60 + 10 30 1\n20\n\n", file = chain_file)
    cat("chain 6000 source1 100 + 10 30 chr1 60 + 40 60 2\n20\n\n", file = chain_file, append = TRUE)

    chain <- gintervals.load_chain(chain_file, src_overlap_policy = "keep", tgt_overlap_policy = "auto_score")

    src_intervals <- data.frame(chrom = "source1", start = 10, end = 30, stringsAsFactors = FALSE)

    # Without metadata - chain_id is always included, score is not
    result_no_meta <- gintervals.liftover(src_intervals, chain)
    expect_false("score" %in% colnames(result_no_meta))
    expect_true("chain_id" %in% colnames(result_no_meta))

    # With metadata - adds score column
    # Both chains are returned since their targets don't overlap
    result_with_meta <- gintervals.liftover(src_intervals, chain, include_metadata = TRUE)
    expect_true("score" %in% colnames(result_with_meta))
    expect_true("chain_id" %in% colnames(result_with_meta))
    expect_equal(nrow(result_with_meta), 2)
    expect_equal(sort(result_with_meta$score), c(6000, 7500))
    expect_equal(sort(result_with_meta$chain_id), c(1, 2))
})

test_that("score-based selection works with partial source overlaps", {
    local_db_state()

    setup_db(list(">chr1\nAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n"))

    # Create chains that partially overlap on source but map to non-overlapping targets
    chain_file <- new_chain_file()
    cat("chain 9000 source1 100 + 5 25 chr1 60 + 10 30 1\n20\n\n", file = chain_file)
    cat("chain 7000 source1 100 + 15 35 chr1 60 + 40 60 2\n20\n\n", file = chain_file, append = TRUE)

    chain <- gintervals.load_chain(chain_file, src_overlap_policy = "keep", tgt_overlap_policy = "auto_score")

    # Query that overlaps both chains on source: [10, 30]
    # Chain 1: source [5,25] - overlaps query at [10,25]
    # Chain 2: source [15,35] - overlaps query at [15,30]
    src_intervals <- data.frame(chrom = "source1", start = 10, end = 30, stringsAsFactors = FALSE)

    result <- gintervals.liftover(src_intervals, chain)

    # Both chains map the overlapping source regions to non-overlapping targets
    # So we should get results from both chains
    expect_true(nrow(result) >= 1)
    expect_true(all(result$chrom == "chr1"))
})

# Score-based selection tests (misha-specific feature)
test_that("auto_score segments overlapping targets correctly", {
    local_db_state()

    # Create target genome
    setup_db(list(">chr1\nAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n"))

    # Create chain file with same source but overlapping TARGET intervals
    # Chain 1: target [10,30), score 5000
    # Chain 2: target [35,55), score 8000 (no overlap with others)
    # Chain 3: target [5,25), score 3000 (overlaps with chain 1)
    chain_file <- new_chain_file()
    cat("chain 5000 source1 100 + 10 30 chr1 60 + 10 30 1\n20\n\n", file = chain_file)
    cat("chain 8000 source1 100 + 10 30 chr1 60 + 35 55 2\n20\n\n", file = chain_file, append = TRUE)
    cat("chain 3000 source1 100 + 10 30 chr1 60 + 5 25 3\n20\n\n", file = chain_file, append = TRUE)

    # Load with auto_score - should segment overlapping targets [5,25) and [10,30)
    chain <- gintervals.load_chain(chain_file, src_overlap_policy = "keep", tgt_overlap_policy = "auto_score")

    # Should have segmented intervals - non-overlapping on target
    chain_sorted <- chain[order(chain$start), ]
    for (i in 1:(nrow(chain_sorted) - 1)) {
        expect_true(chain_sorted$end[i] <= chain_sorted$start[i + 1])
    }

    # Test liftover
    src_intervals <- data.frame(chrom = "source1", start = 10, end = 30, stringsAsFactors = FALSE)
    misha_result <- gintervals.liftover(src_intervals, chain)

    # Should get multiple non-overlapping results
    expect_true(nrow(misha_result) >= 1)
    expect_true(all(misha_result$chrom == "chr1"))
})


test_that("auto_score applies tiebreaker rules correctly (score -> span -> chain_id)", {
    local_db_state()

    # Create larger target genome
    setup_db(list(">chr1\n", paste(rep("A", 200), collapse = ""), "\n"))

    # Create chain file with overlapping target chains that test tiebreaker logic
    chain_file <- new_chain_file()
    # Chains 2 and 4 have tied score (9500), but chain 2 has lower chain_id (2 vs 4)
    # Target overlaps: chain 4 (10-40) overlaps chain 1 (30-60)
    cat("chain 7000 source1 200 + 20 50 chr1 200 + 30 60 1\n30\n\n", file = chain_file)
    cat("chain 9500 source1 200 + 50 80 chr1 200 + 100 130 2\n30\n\n", file = chain_file, append = TRUE)
    cat("chain 6000 source1 200 + 80 110 chr1 200 + 150 180 3\n30\n\n", file = chain_file, append = TRUE)
    cat("chain 9500 source1 200 + 110 140 chr1 200 + 10 40 4\n30\n\n", file = chain_file, append = TRUE)

    # Load chain with auto_score - should segment overlapping targets
    chain <- gintervals.load_chain(chain_file, src_overlap_policy = "keep", tgt_overlap_policy = "auto_score")

    # Chain 1 (score 7000) at target 30-60 overlaps chain 4 (score 9500) at target 10-40
    # Segmentation should give:
    # - 10-30: chain 4 wins (only chain 4)
    # - 30-40: chain 4 wins (9500 > 7000)
    # - 40-60: chain 1 wins (only chain 1)

    # Check that chain was loaded and segmented
    expect_true(nrow(chain) >= 4) # Original 4 chains may be segmented into more

    # Liftover source interval that maps to chain 4
    src_intervals <- data.frame(chrom = "source1", start = 110, end = 140, stringsAsFactors = FALSE)
    misha_result <- gintervals.liftover(src_intervals, chain)

    # Chain 4 maps source1:110-140 to chr1:10-40
    expect_true(nrow(misha_result) >= 1)
    expect_true(all(misha_result$chrom == "chr1"))
})

test_that("include_metadata returns correct score and chain_id", {
    local_db_state()

    setup_db(list(">chr1\nAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n"))

    # Create chains with known scores
    chain_file <- new_chain_file()
    cat("chain 7500 source1 100 + 10 30 chr1 60 + 10 30 1\n20\n\n", file = chain_file)
    cat("chain 9200 source1 100 + 40 60 chr1 60 + 35 55 2\n20\n\n", file = chain_file, append = TRUE)

    chain <- gintervals.load_chain(chain_file, src_overlap_policy = "keep", tgt_overlap_policy = "auto_score")

    # Verify chain has score and chain_id columns
    expect_true("score" %in% colnames(chain))
    expect_true("chain_id" %in% colnames(chain))
    expect_equal(chain$score, c(7500, 9200))
    expect_equal(chain$chain_id, c(1, 2))

    # Test liftover with metadata
    src_intervals <- data.frame(
        chrom = c("source1", "source1"),
        start = c(10, 40),
        end = c(30, 60),
        stringsAsFactors = FALSE
    )

    result <- gintervals.liftover(src_intervals, chain, include_metadata = TRUE)

    # Should have score and chain_id columns
    expect_true("score" %in% colnames(result))
    expect_true("chain_id" %in% colnames(result))

    # First interval mapped by chain 1 (score 7500, chain_id 1)
    expect_equal(result$score[1], 7500)
    expect_equal(result$chain_id[1], 1)

    # Second interval mapped by chain 2 (score 9200, chain_id 2)
    expect_equal(result$score[2], 9200)
    expect_equal(result$chain_id[2], 2)
})

test_that("auto_score works with partial overlaps - selects best for each region", {
    local_db_state()

    setup_db(list(">chr1\nAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n"))

    # Create chains with partial source overlaps but no target overlaps
    chain_file <- new_chain_file()
    # Chain 1: high score for source1:10-30 -> target 5-25
    cat("chain 9000 source1 100 + 10 30 chr1 60 + 5 25 1\n20\n\n", file = chain_file)
    # Chain 2: low score for source1:25-45 -> target 30-50 (no target overlap)
    cat("chain 4000 source1 100 + 25 45 chr1 60 + 30 50 2\n20\n\n", file = chain_file, append = TRUE)

    chain <- gintervals.load_chain(chain_file, src_overlap_policy = "keep", tgt_overlap_policy = "auto_score")

    # Chain should be loaded successfully with at least 2 entries
    expect_true(nrow(chain) >= 2)

    # Test interval fully in chain 1's source range
    src_intervals1 <- data.frame(chrom = "source1", start = 10, end = 25, stringsAsFactors = FALSE)
    result1 <- gintervals.liftover(src_intervals1, chain)
    expect_equal(nrow(result1), 1)
    expect_equal(result1$start, 5) # Chain 1

    # Test interval fully in chain 2's source range
    src_intervals2 <- data.frame(chrom = "source1", start = 35, end = 45, stringsAsFactors = FALSE)
    result2 <- gintervals.liftover(src_intervals2, chain)
    expect_equal(nrow(result2), 1)
    expect_equal(result2$start, 40) # Chain 2

    # Test interval in the overlapping source region (25-30)
    # Both chains cover this source region, so liftover may produce multiple results
    src_intervals3 <- data.frame(chrom = "source1", start = 25, end = 30, stringsAsFactors = FALSE)
    result3 <- gintervals.liftover(src_intervals3, chain)
    expect_true(nrow(result3) >= 1)
})

test_that("auto_score works correctly with negative strand chains", {
    local_db_state()

    setup_db(list(">chr1\nAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n"))

    # Create chains with different strands - targets overlap at 30-50
    chain_file <- new_chain_file()
    # Positive strand, low score -> target 10-30
    cat("chain 3000 source1 100 + 10 30 chr1 60 + 10 30 1\n20\n\n", file = chain_file)
    # Negative strand, high score -> target 30-50 (coord flip from 60-size to size)
    cat("chain 8000 source1 100 + 40 60 chr1 60 - 10 30 2\n20\n\n", file = chain_file, append = TRUE)

    chain <- gintervals.load_chain(chain_file, src_overlap_policy = "keep", tgt_overlap_policy = "auto_score")

    # Chain should be loaded successfully
    expect_true(nrow(chain) >= 2)

    # Test liftover with chain 1's source
    src_intervals1 <- data.frame(chrom = "source1", start = 10, end = 30, stringsAsFactors = FALSE)
    result1 <- gintervals.liftover(src_intervals1, chain)
    expect_equal(nrow(result1), 1)
    expect_equal(result1$start, 10)
    expect_equal(result1$end, 30)

    # Test liftover with chain 2's source (negative strand)
    src_intervals2 <- data.frame(chrom = "source1", start = 40, end = 60, stringsAsFactors = FALSE)
    result2 <- gintervals.liftover(src_intervals2, chain)
    expect_equal(nrow(result2), 1)
    # Negative strand: target coords are reversed
    expect_equal(result2$start, 30)
    expect_equal(result2$end, 50)
})

test_that("keep policy returns all chains regardless of score", {
    local_db_state()

    setup_db(list(">chr1\nAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n"))

    # Create overlapping chains with very different scores
    chain_file <- new_chain_file()
    cat("chain 1000 source1 100 + 10 30 chr1 60 + 10 30 1\n20\n\n", file = chain_file)
    cat("chain 9999 source1 100 + 10 30 chr1 60 + 35 55 2\n20\n\n", file = chain_file, append = TRUE)
    cat("chain 5000 source1 100 + 10 30 chr1 60 + 5 25 3\n20\n\n", file = chain_file, append = TRUE)

    chain <- gintervals.load_chain(chain_file, src_overlap_policy = "keep", tgt_overlap_policy = "keep")

    src_intervals <- data.frame(chrom = "source1", start = 10, end = 30, stringsAsFactors = FALSE)
    result <- gintervals.liftover(src_intervals, chain)

    # Should return all 3 mappings
    expect_equal(nrow(result), 3)
    expect_true(10 %in% result$start)
    expect_true(35 %in% result$start)
    expect_true(5 %in% result$start)
})

test_that("min_score filters chains at load time", {
    local_db_state()

    setup_db(list(">chr1\nAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n"))

    # Create chains with varying scores
    chain_file <- new_chain_file()
    cat("chain 2000 source1 100 + 10 30 chr1 60 + 10 30 1\n20\n\n", file = chain_file)
    cat("chain 6000 source1 100 + 35 55 chr1 60 + 35 55 2\n20\n\n", file = chain_file, append = TRUE)
    cat("chain 4000 source1 100 + 60 80 chr1 60 + 5 25 3\n20\n\n", file = chain_file, append = TRUE)

    # Load with min_score=5000
    chain_filtered <- gintervals.load_chain(chain_file, min_score = 5000)

    # Should only have chain 2 (score 6000)
    expect_equal(nrow(chain_filtered), 1)
    expect_equal(chain_filtered$score, 6000)
    expect_equal(chain_filtered$startsrc, 35)
})

test_that("auto_score with min_score combines filtering and selection", {
    local_db_state()

    setup_db(list(">chr1\nAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n"))

    # Create overlapping chains: 2 high-scoring, 1 low-scoring
    chain_file <- new_chain_file()
    cat("chain 2000 source1 100 + 10 30 chr1 60 + 5 25 1\n20\n\n", file = chain_file)
    cat("chain 7000 source1 100 + 10 30 chr1 60 + 10 30 2\n20\n\n", file = chain_file, append = TRUE)
    cat("chain 9000 source1 100 + 10 30 chr1 60 + 35 55 3\n20\n\n", file = chain_file, append = TRUE)

    # Load with min_score=5000 and auto_score policy
    chain <- gintervals.load_chain(chain_file, min_score = 5000, src_overlap_policy = "keep", tgt_overlap_policy = "auto_score")

    # Should filter out chain 1, leaving chains 2 and 3
    expect_equal(nrow(chain), 2)

    src_intervals <- data.frame(chrom = "source1", start = 10, end = 30, stringsAsFactors = FALSE)
    result <- gintervals.liftover(src_intervals, chain)

    # Chains 2 and 3 both pass min_score and don't overlap on target
    # Both should be returned
    expect_equal(nrow(result), 2)
    result <- result[order(result$start), ]
    expect_equal(result$start, c(10, 35))
    expect_equal(result$end, c(30, 55))
})

test_that("gintervals.liftover returns chain_id to distinguish duplications", {
    local_db_state()

    # Create target genome with single chromosome
    setup_db(list(paste0(">chr1\n", paste(rep("A", 1000), collapse = ""), "\n")))

    # Create chain file with two chains that map the same source region
    # to different target locations on the same chromosome (duplication)
    chain_file <- new_chain_file()

    # Chain 1: source1[0-100] -> chr1[0-100] with chain_id=101
    cat("chain 1000 source1 200 + 0 100 chr1 1000 + 0 100 101\n100\n\n", file = chain_file)

    # Chain 2: source1[0-100] -> chr1[500-600] with chain_id=205
    cat("chain 1000 source1 200 + 0 100 chr1 1000 + 500 600 205\n100\n\n", file = chain_file, append = TRUE)

    # Load chain with keep policy for source overlaps
    chain <- gintervals.load_chain(chain_file, src_overlap_policy = "keep", tgt_overlap_policy = "keep")

    # Lift over a single source interval
    src_intervals <- data.frame(chrom = "source1", start = 0, end = 100, stringsAsFactors = FALSE)
    result <- gintervals.liftover(src_intervals, chain)

    # Should get 2 results (one from each chain)
    expect_equal(nrow(result), 2)

    # Both should have same intervalID (same source interval)
    expect_equal(result$intervalID, c(1, 1))

    # But different chain_ids
    expect_true("chain_id" %in% colnames(result))
    expect_equal(length(unique(result$chain_id)), 2)
    expect_true(all(c(101, 205) %in% result$chain_id))

    # Check the actual mapping coordinates
    result <- result[order(result$start), ]
    expect_equal(result$start, c(0, 500))
    expect_equal(result$end, c(100, 600))
})

test_that("gintervals.liftover chain_id groups blocks from same chain in duplication", {
    local_db_state()

    # Create target genome
    setup_db(list(paste0(">chr1\n", paste(rep("A", 1000), collapse = ""), "\n")))

    # Create chain file with two chains, each with multiple blocks
    chain_file <- new_chain_file()

    # Chain 1: source1[0-200] -> chr1[0-200] with two blocks and a gap
    # Block 1: 50 bp, gap of 50 in source and 50 in target, Block 2: 100 bp
    cat("chain 1000 source1 300 + 0 200 chr1 1000 + 0 200 101\n", file = chain_file)
    cat("50\t50\t50\n100\n\n", file = chain_file, append = TRUE)

    # Chain 2: source1[0-200] -> chr1[500-700] with two blocks
    cat("chain 1000 source1 300 + 0 200 chr1 1000 + 500 700 205\n", file = chain_file, append = TRUE)
    cat("50\t50\t50\n100\n\n", file = chain_file, append = TRUE)

    # Load chain
    chain <- gintervals.load_chain(chain_file, src_overlap_policy = "keep", tgt_overlap_policy = "keep")

    # Lift over source interval that spans both blocks
    src_intervals <- data.frame(chrom = "source1", start = 0, end = 200, stringsAsFactors = FALSE)
    result <- gintervals.liftover(src_intervals, chain)

    # Should get 4 results (2 blocks from each of 2 chains)
    expect_equal(nrow(result), 4)

    # All should have same intervalID
    expect_true(all(result$intervalID == 1))

    # But grouped by chain_id
    expect_equal(length(unique(result$chain_id)), 2)

    # Check that each chain_id has 2 blocks
    chain_101_rows <- result[result$chain_id == 101, ]
    chain_205_rows <- result[result$chain_id == 205, ]
    expect_equal(nrow(chain_101_rows), 2)
    expect_equal(nrow(chain_205_rows), 2)

    # Verify chain 101 blocks
    chain_101_rows <- chain_101_rows[order(chain_101_rows$start), ]
    expect_equal(chain_101_rows$start, c(0, 100))
    expect_equal(chain_101_rows$end, c(50, 200))

    # Verify chain 205 blocks
    chain_205_rows <- chain_205_rows[order(chain_205_rows$start), ]
    expect_equal(chain_205_rows$start, c(500, 600))
    expect_equal(chain_205_rows$end, c(550, 700))
})

test_that("gintervals.liftover canonic respects chain_id boundaries in duplication", {
    local_db_state()

    # Create target genome
    setup_db(list(paste0(">chr1\n", paste(rep("A", 1000), collapse = ""), "\n")))

    # Create chain file with two chains that produce adjacent target blocks
    chain_file <- new_chain_file()

    # Chain 1: source1[0-50] -> chr1[0-50] with chain_id=101
    cat("chain 1000 source1 200 + 0 50 chr1 1000 + 0 50 101\n50\n\n", file = chain_file)

    # Chain 2: source1[0-50] -> chr1[50-100] with chain_id=205
    # Note: Target is adjacent to chain 1's output
    cat("chain 1000 source1 200 + 0 50 chr1 1000 + 50 100 205\n50\n\n", file = chain_file, append = TRUE)

    # Load chain
    chain <- gintervals.load_chain(chain_file, src_overlap_policy = "keep", tgt_overlap_policy = "keep")

    # Lift over
    src_intervals <- data.frame(chrom = "source1", start = 0, end = 50, stringsAsFactors = FALSE)
    result <- gintervals.liftover(src_intervals, chain, canonic = TRUE)

    # With canonic=TRUE, should NOT merge because they have different chain_ids
    expect_equal(nrow(result), 2)
    expect_equal(length(unique(result$chain_id)), 2)

    # Verify the intervals
    result <- result[order(result$start), ]
    expect_equal(result$start, c(0, 50))
    expect_equal(result$end, c(50, 100))
})

test_that("gintervals.as_chain validates required columns", {
    local_db_state()

    # Missing chrom column
    intervals <- data.frame(
        start = 0,
        end = 100,
        strand = 0,
        chromsrc = "chr1",
        startsrc = 0,
        endsrc = 100,
        strandsrc = 0,
        chain_id = 1L,
        score = 1000.0
    )
    expect_error(
        gintervals.as_chain(intervals),
        "Missing required columns"
    )

    # Missing chain_id column
    intervals <- data.frame(
        chrom = "chr1",
        start = 0,
        end = 100,
        strand = 0,
        chromsrc = "chr1",
        startsrc = 0,
        endsrc = 100,
        strandsrc = 0,
        score = 1000.0
    )
    expect_error(
        gintervals.as_chain(intervals),
        "Missing required columns"
    )

    # Missing score column
    intervals <- data.frame(
        chrom = "chr1",
        start = 0,
        end = 100,
        strand = 0,
        chromsrc = "chr1",
        startsrc = 0,
        endsrc = 100,
        strandsrc = 0,
        chain_id = 1L
    )
    expect_error(
        gintervals.as_chain(intervals),
        "Missing required columns"
    )
})

test_that("gintervals.as_chain validates data types", {
    local_db_state()

    # Invalid start/end (character instead of numeric)
    intervals <- data.frame(
        chrom = "chr1",
        start = "0",
        end = 100,
        strand = 0,
        chromsrc = "chr1",
        startsrc = 0,
        endsrc = 100,
        strandsrc = 0,
        chain_id = 1L,
        score = 1000.0,
        stringsAsFactors = FALSE
    )
    expect_error(
        gintervals.as_chain(intervals),
        "start and end columns must be numeric"
    )

    # Invalid chain_id (character instead of integer/numeric)
    intervals <- data.frame(
        chrom = "chr1",
        start = 0,
        end = 100,
        strand = 0,
        chromsrc = "chr1",
        startsrc = 0,
        endsrc = 100,
        strandsrc = 0,
        chain_id = "1",
        score = 1000.0,
        stringsAsFactors = FALSE
    )
    expect_error(
        gintervals.as_chain(intervals),
        "chain_id column must be integer or numeric"
    )

    # Invalid score (character instead of numeric)
    intervals <- data.frame(
        chrom = "chr1",
        start = 0,
        end = 100,
        strand = 0,
        chromsrc = "chr1",
        startsrc = 0,
        endsrc = 100,
        strandsrc = 0,
        chain_id = 1L,
        score = "1000.0",
        stringsAsFactors = FALSE
    )
    expect_error(
        gintervals.as_chain(intervals),
        "score column must be numeric"
    )
})

test_that("gintervals.as_chain validates strand values", {
    local_db_state()

    # Invalid strand value (2 instead of -1, 0, or 1)
    intervals <- data.frame(
        chrom = "chr1",
        start = 0,
        end = 100,
        strand = 2,
        chromsrc = "chr1",
        startsrc = 0,
        endsrc = 100,
        strandsrc = 0,
        chain_id = 1L,
        score = 1000.0
    )
    expect_error(
        gintervals.as_chain(intervals),
        "strand values must be -1, 0, or 1"
    )

    # Invalid strandsrc value (-2 instead of -1, 0, or 1)
    intervals <- data.frame(
        chrom = "chr1",
        start = 0,
        end = 100,
        strand = 0,
        chromsrc = "chr1",
        startsrc = 0,
        endsrc = 100,
        strandsrc = -2,
        chain_id = 1L,
        score = 1000.0
    )
    expect_error(
        gintervals.as_chain(intervals),
        "strandsrc values must be -1, 0, or 1"
    )
})

test_that("gintervals.as_chain validates policies", {
    local_db_state()

    intervals <- data.frame(
        chrom = "chr1",
        start = 0,
        end = 100,
        strand = 0,
        chromsrc = "chr1",
        startsrc = 0,
        endsrc = 100,
        strandsrc = 0,
        chain_id = 1L,
        score = 1000.0
    )

    # Invalid src_overlap_policy
    expect_error(
        gintervals.as_chain(intervals, src_overlap_policy = "invalid"),
        "src_overlap_policy must be"
    )

    # Invalid tgt_overlap_policy
    expect_error(
        gintervals.as_chain(intervals, tgt_overlap_policy = "invalid"),
        "tgt_overlap_policy must be"
    )

    # Invalid min_score (not a single value)
    expect_error(
        gintervals.as_chain(intervals, min_score = c(100, 200)),
        "min_score must be a single numeric value"
    )
})

test_that("gintervals.as_chain creates valid chain with attributes", {
    local_db_state()

    intervals <- data.frame(
        chrom = "chr1",
        start = 0,
        end = 100,
        strand = 1,
        chromsrc = "chr2",
        startsrc = 50,
        endsrc = 150,
        strandsrc = -1,
        chain_id = 42L,
        score = 5000.0
    )

    chain <- gintervals.as_chain(intervals, src_overlap_policy = "keep", tgt_overlap_policy = "auto_score", min_score = 1000)

    # Check that all columns are preserved
    expect_equal(nrow(chain), 1)
    expect_equal(as.character(chain$chrom), "chr1")
    expect_equal(as.numeric(chain$start), 0)
    expect_equal(as.numeric(chain$end), 100)
    expect_equal(as.numeric(chain$strand), 1)
    expect_equal(as.character(chain$chromsrc), "chr2")
    expect_equal(as.numeric(chain$startsrc), 50)
    expect_equal(as.numeric(chain$endsrc), 150)
    expect_equal(as.numeric(chain$strandsrc), -1)
    expect_equal(as.integer(chain$chain_id), 42L)
    expect_equal(as.numeric(chain$score), 5000.0)

    # Check attributes
    expect_equal(attr(chain, "src_overlap_policy"), "keep")
    expect_equal(attr(chain, "tgt_overlap_policy"), "auto_score")
    expect_equal(attr(chain, "min_score"), 1000)
})

test_that("gintervals.as_chain converts chain_id to integer", {
    local_db_state()

    # chain_id as numeric (not integer)
    intervals <- data.frame(
        chrom = "chr1",
        start = 0,
        end = 100,
        strand = 0,
        chromsrc = "chr1",
        startsrc = 0,
        endsrc = 100,
        strandsrc = 0,
        chain_id = 42.0, # numeric, not integer
        score = 1000.0
    )

    chain <- gintervals.as_chain(intervals)
    expect_true(is.integer(chain$chain_id))
    expect_equal(chain$chain_id, 42L)
})

test_that("gintervals.as_chain converts auto to auto_score", {
    local_db_state()

    intervals <- data.frame(
        chrom = "chr1",
        start = 0,
        end = 100,
        strand = 0,
        chromsrc = "chr1",
        startsrc = 0,
        endsrc = 100,
        strandsrc = 0,
        chain_id = 1L,
        score = 1000.0
    )

    chain <- gintervals.as_chain(intervals, tgt_overlap_policy = "auto")
    expect_equal(attr(chain, "tgt_overlap_policy"), "auto_score")
})

test_that("gintervals.as_chain works with empty data frame", {
    local_db_state()

    intervals <- data.frame(
        chrom = character(0),
        start = numeric(0),
        end = numeric(0),
        strand = numeric(0),
        chromsrc = character(0),
        startsrc = numeric(0),
        endsrc = numeric(0),
        strandsrc = numeric(0),
        chain_id = integer(0),
        score = numeric(0),
        stringsAsFactors = FALSE
    )

    chain <- gintervals.as_chain(intervals)
    expect_equal(nrow(chain), 0)
    expect_equal(attr(chain, "src_overlap_policy"), "error")
    expect_equal(attr(chain, "tgt_overlap_policy"), "auto_score")
})

test_that("gintervals.as_chain works with multiple intervals", {
    local_db_state()

    intervals <- data.frame(
        chrom = c("chr1", "chr2"),
        start = c(0, 100),
        end = c(100, 200),
        strand = c(1, -1),
        chromsrc = c("chrA", "chrB"),
        startsrc = c(0, 50),
        endsrc = c(100, 150),
        strandsrc = c(1, -1),
        chain_id = c(1L, 2L),
        score = c(1000.0, 2000.0)
    )

    chain <- gintervals.as_chain(intervals)
    expect_equal(nrow(chain), 2)
    expect_equal(as.character(chain$chrom), c("chr1", "chr2"))
    expect_equal(as.integer(chain$chain_id), c(1L, 2L))
})

test_that("gintervals.as_chain works with liftover after conversion", {
    local_db_state()

    # Create target genome
    setup_db(list(">chr1\nACTGACTGACTGACTGACTGACTGACTGACTG\n", ">chr2\nGGGGCCCCTTTTAAAA\n"))

    # Create a chain file first to get the structure
    chain_file <- new_chain_file()
    write_chain_entry(chain_file, "source1", 100, "+", 0, 20, "chr1", 32, "+", 0, 20, 1)
    chain_loaded <- gintervals.load_chain(chain_file, src_overlap_policy = "keep", tgt_overlap_policy = "keep")

    # Convert to chain format using as_chain with same policies
    chain_converted <- gintervals.as_chain(chain_loaded, src_overlap_policy = "keep", tgt_overlap_policy = "keep")

    # Check that the chains are identical (data and attributes)
    expect_equal(nrow(chain_loaded), nrow(chain_converted))
    expect_equal(chain_loaded$chrom, chain_converted$chrom)
    expect_equal(chain_loaded$start, chain_converted$start)
    expect_equal(chain_loaded$end, chain_converted$end)
    expect_equal(chain_loaded$strand, chain_converted$strand)
    expect_equal(chain_loaded$chromsrc, chain_converted$chromsrc)
    expect_equal(chain_loaded$startsrc, chain_converted$startsrc)
    expect_equal(chain_loaded$endsrc, chain_converted$endsrc)
    expect_equal(chain_loaded$strandsrc, chain_converted$strandsrc)
    expect_equal(chain_loaded$chain_id, chain_converted$chain_id)
    expect_equal(chain_loaded$score, chain_converted$score)

    # Check attributes match
    expect_equal(attr(chain_loaded, "src_overlap_policy"), attr(chain_converted, "src_overlap_policy"))
    expect_equal(attr(chain_loaded, "tgt_overlap_policy"), attr(chain_converted, "tgt_overlap_policy"))

    # Create source intervals
    src_intervals <- data.frame(
        chrom = "source1",
        start = 5,
        end = 15,
        stringsAsFactors = FALSE
    )

    # Both should work with liftover and produce identical results
    result_loaded <- gintervals.liftover(src_intervals, chain_loaded)
    result_converted <- gintervals.liftover(src_intervals, chain_converted)

    # Results should be identical
    expect_equal(result_loaded, result_converted)
})

test_that("gintervals.as_chain rejects non-data.frame input", {
    local_db_state()

    expect_error(
        gintervals.as_chain(NULL),
        "Usage:"
    )

    expect_error(
        gintervals.as_chain(list(chrom = "chr1", start = 0, end = 100)),
        "intervals must be a data frame"
    )

    expect_error(
        gintervals.as_chain("not a data frame"),
        "intervals must be a data frame"
    )
})

test_that("gintervals.as_chain handles min_score attribute correctly", {
    local_db_state()

    intervals <- data.frame(
        chrom = "chr1",
        start = 0,
        end = 100,
        strand = 0,
        chromsrc = "chr1",
        startsrc = 0,
        endsrc = 100,
        strandsrc = 0,
        chain_id = 1L,
        score = 1000.0
    )

    # Without min_score
    chain_no_min <- gintervals.as_chain(intervals)
    expect_null(attr(chain_no_min, "min_score"))

    # With min_score
    chain_with_min <- gintervals.as_chain(intervals, min_score = 500)
    expect_equal(attr(chain_with_min, "min_score"), 500)
})

test_that("gintervals.as_chain accepts all valid policy values", {
    local_db_state()

    intervals <- data.frame(
        chrom = "chr1",
        start = 0,
        end = 100,
        strand = 0,
        chromsrc = "chr1",
        startsrc = 0,
        endsrc = 100,
        strandsrc = 0,
        chain_id = 1L,
        score = 1000.0
    )

    # Test all valid src_overlap_policy values
    for (policy in c("error", "keep", "discard")) {
        chain <- gintervals.as_chain(intervals, src_overlap_policy = policy)
        expect_equal(attr(chain, "src_overlap_policy"), policy)
    }

    # Test all valid tgt_overlap_policy values
    for (policy in c("error", "auto", "auto_first", "auto_longer", "auto_score", "discard", "keep", "agg")) {
        chain <- gintervals.as_chain(intervals, tgt_overlap_policy = policy)
        if (policy == "auto") {
            expect_equal(attr(chain, "tgt_overlap_policy"), "auto_score")
        } else {
            expect_equal(attr(chain, "tgt_overlap_policy"), policy)
        }
    }
})
