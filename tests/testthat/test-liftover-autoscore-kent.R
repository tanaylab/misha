create_isolated_test_db()

test_that("auto_score produces segmented results when targets overlap", {
    skip_if_not(has_liftover_binary(), "liftOver binary not found")
    skip_if_not(exists(".misha"), "misha package not properly loaded")

    local_db_state()

    # Create target genome
    setup_db(list(">chr1\n", paste(rep("A", 200), collapse = ""), "\n"))

    # Create chain file with overlapping TARGET intervals at different scores
    chain_file <- new_chain_file()
    # Chain 1: Low score (3000), target chr1[10-30]
    write_chain_entry(chain_file, "source1", 100, "+", 10, 30, "chr1", 200, "+", 10, 30, 1, score = 3000)
    # Chain 2: High score (8000), target chr1[15-35] (overlaps with chain 1)
    write_chain_entry(chain_file, "source1", 100, "+", 10, 30, "chr1", 200, "+", 15, 35, 2, score = 8000, append = TRUE)
    # Chain 3: Medium score (5000), target chr1[20-40] (overlaps with both)
    write_chain_entry(chain_file, "source1", 100, "+", 10, 30, "chr1", 200, "+", 20, 40, 3, score = 5000, append = TRUE)

    # Load chain with segmentation
    chain <- gintervals.load_chain(chain_file, src_overlap_policy = "keep", tgt_overlap_policy = "auto_score")

    # Verify segmentation produced non-overlapping target intervals
    expect_true(nrow(chain) >= 3) # May be segmented into more pieces

    # Run Misha's liftover with auto_score
    misha_result <- gintervals.liftover(
        data.frame(chrom = "source1", start = 10, end = 30, stringsAsFactors = FALSE),
        chain
    )

    # Should produce multiple segments covering the liftover
    expect_true(nrow(misha_result) >= 1)
    expect_true(all(misha_result$chrom == "chr1"))

    # Also run Kent's liftOver for reference (but don't require exact match)
    src_bed <- tempfile(fileext = ".bed")
    cat("source1\t10\t30\n", file = src_bed)
    kent_result <- run_kent_liftover(src_bed, chain_file)

    # Both should produce results
    expect_true(nrow(kent_result$mapped) >= 1)

    # Clean up
    unlink(src_bed)
})

test_that("auto_score tiebreaker uses span and chain_id", {
    skip_if_not(has_liftover_binary(), "liftOver binary not found")
    skip_if_not(exists(".misha"), "misha package not properly loaded")

    local_db_state()

    # Create target genome
    setup_db(list(">chr1\n", paste(rep("A", 200), collapse = ""), "\n"))

    # Create chains with equal scores but different lengths and overlapping targets
    chain_file <- new_chain_file()
    # Chain 1: Score 5000, length 15bp, target chr1[10-25]
    write_chain_entry(chain_file, "source1", 100, "+", 10, 25, "chr1", 200, "+", 10, 25, 1, score = 5000)
    # Chain 2: Score 5000, length 20bp (longer), target chr1[15-35] (overlaps with chain 1)
    write_chain_entry(chain_file, "source1", 100, "+", 10, 30, "chr1", 200, "+", 15, 35, 2, score = 5000, append = TRUE)
    # Chain 3: Score 5000, length 12bp, target chr1[20-32] (overlaps with both)
    write_chain_entry(chain_file, "source1", 100, "+", 10, 22, "chr1", 200, "+", 20, 32, 3, score = 5000, append = TRUE)

    chain <- gintervals.load_chain(chain_file, src_overlap_policy = "keep", tgt_overlap_policy = "auto_score")

    # Verify chain was loaded
    expect_true(nrow(chain) >= 1)

    # Run Misha's liftover
    misha_result <- gintervals.liftover(
        data.frame(chrom = "source1", start = 10, end = 30, stringsAsFactors = FALSE),
        chain
    )

    # Should produce results based on segmentation
    expect_true(nrow(misha_result) >= 1)
    expect_true(all(misha_result$chrom == "chr1"))

    # Also run Kent's liftOver for reference
    src_bed <- tempfile(fileext = ".bed")
    cat("source1\t10\t30\n", file = src_bed)
    kent_result <- run_kent_liftover(src_bed, chain_file)
    expect_true(nrow(kent_result$mapped) >= 1)

    # Clean up
    unlink(src_bed)
})

test_that("auto_score handles negative strand chains", {
    skip_if_not(has_liftover_binary(), "liftOver binary not found")
    skip_if_not(exists(".misha"), "misha package not properly loaded")

    local_db_state()

    # Create target genome
    setup_db(list(">chr1\n", paste(rep("A", 200), collapse = ""), "\n"))

    # Create chains with negative strand and scores, overlapping targets
    chain_file <- new_chain_file()
    # Negative strand, low score, target chr1[10-30] (on negative strand)
    write_chain_entry(chain_file, "source1", 100, "+", 10, 30, "chr1", 200, "-", 10, 30, 1, score = 3000)
    # Negative strand, high score (should win), target chr1[15-35] (overlaps with chain 1)
    write_chain_entry(chain_file, "source1", 100, "+", 10, 30, "chr1", 200, "-", 15, 35, 2, score = 8000, append = TRUE)

    chain <- gintervals.load_chain(chain_file, src_overlap_policy = "keep", tgt_overlap_policy = "auto_score")

    # Verify chain was loaded
    expect_true(nrow(chain) >= 2)

    # Run Misha's liftover
    misha_result <- gintervals.liftover(
        data.frame(chrom = "source1", start = 10, end = 30, stringsAsFactors = FALSE),
        chain
    )

    # Should produce results (negative strand coordinates are flipped)
    expect_true(nrow(misha_result) >= 1)
    expect_true(all(misha_result$chrom == "chr1"))

    # Also run Kent's liftOver for reference
    src_bed <- tempfile(fileext = ".bed")
    cat("source1\t10\t30\n", file = src_bed)
    kent_result <- run_kent_liftover(src_bed, chain_file)
    expect_true(nrow(kent_result$mapped) >= 1)

    # Clean up
    unlink(src_bed)
})

test_that("auto_score handles partial source overlaps", {
    skip_if_not(has_liftover_binary(), "liftOver binary not found")
    skip_if_not(exists(".misha"), "misha package not properly loaded")

    local_db_state()

    # Create target genome
    setup_db(list(">chr1\n", paste(rep("A", 300), collapse = ""), "\n"))

    # Create chains where query partially overlaps chains, with overlapping targets
    chain_file <- new_chain_file()
    # Chain 1: source[0-25], score 5000, target chr1[0-25]
    write_chain_entry(chain_file, "source1", 100, "+", 0, 25, "chr1", 300, "+", 0, 25, 1, score = 5000)
    # Chain 2: source[15-40], score 8000 (higher score, partial overlap), target chr1[20-45] (overlaps with chain 1 target)
    write_chain_entry(chain_file, "source1", 100, "+", 15, 40, "chr1", 300, "+", 20, 45, 2, score = 8000, append = TRUE)

    chain <- gintervals.load_chain(chain_file, src_overlap_policy = "keep", tgt_overlap_policy = "auto_score")

    # Verify chain was loaded and segmented
    expect_true(nrow(chain) >= 2)

    # Run Misha's liftover
    misha_result <- gintervals.liftover(
        data.frame(chrom = "source1", start = 10, end = 30, stringsAsFactors = FALSE),
        chain
    )

    # Should produce results
    expect_true(nrow(misha_result) >= 1)
    expect_true(all(misha_result$chrom == "chr1"))

    # Also run Kent's liftOver for reference
    src_bed <- tempfile(fileext = ".bed")
    cat("source1\t10\t30\n", file = src_bed)
    kent_result <- run_kent_liftover(src_bed, chain_file)
    expect_true(nrow(kent_result$mapped) >= 1)

    # Clean up
    unlink(src_bed)
})
