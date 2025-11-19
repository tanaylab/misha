load_test_db()

test_that("gintervals.liftover with canonic merges adjacent blocks from same interval", {
    local_db_state()

    # Create target genome with chr1 of 300 bp
    setup_db(list(paste0(">chr1\n", paste(rep("A", 300), collapse = ""), "\n")))

    # Create chain with internal gap that produces ADJACENT target blocks:
    # Source [0-200] maps with:
    #   Block 1: 50 bp aligned (source [0-50] -> target [0-50])
    #   Gap: 100 in source (dq), 0 in target (dt)
    #   Block 2: 50 bp aligned (source [150-200] -> target [50-100])
    # So target blocks are [0-50] and [50-100] - adjacent!
    chain_file <- new_chain_file()

    # Chain format: srcName srcSize srcStrand srcStart srcEnd tgtName tgtSize tgtStrand tgtStart tgtEnd
    # srcName = source (what we're lifting from), tgtName = target (in database)
    cat("chain 1000 source1 300 + 0 200 chr1 300 + 0 100 1\n", file = chain_file)
    cat("50\t100\t0\n", file = chain_file, append = TRUE) # block 50, gap: 100 in source, 0 in target
    cat("50\n\n", file = chain_file, append = TRUE) # final block 50

    # Test: single interval spanning both blocks
    src_intervals <- data.frame(
        chrom = "source1",
        start = 0,
        end = 200,
        stringsAsFactors = FALSE
    )

    # Load chain and run liftover
    chain <- gintervals.load_chain(chain_file)

    # Without canonic: should return 2 adjacent intervals
    result_no_canonic <- gintervals.liftover(src_intervals, chain)
    expect_equal(nrow(result_no_canonic), 2)
    expect_true("intervalID" %in% colnames(result_no_canonic))

    # Both should have same intervalID
    expect_equal(unique(result_no_canonic$intervalID), 1)

    # Verify they are adjacent: [0-50] and [50-100]
    result_no_canonic <- result_no_canonic[order(result_no_canonic$start), ]
    expect_equal(result_no_canonic$start, c(0, 50))
    expect_equal(result_no_canonic$end, c(50, 100))

    # With canonic=TRUE: should merge into single interval [0-100]
    result_canonic <- gintervals.liftover(src_intervals, chain, canonic = TRUE)
    expect_equal(nrow(result_canonic), 1)
    expect_equal(result_canonic$start, 0)
    expect_equal(result_canonic$end, 100)

    # intervalID and chain_id columns should be present
    expect_true("intervalID" %in% colnames(result_canonic))
    expect_true("chain_id" %in% colnames(result_canonic))
})

test_that("gintervals.liftover canonic does not merge intervals from different source intervals", {
    local_db_state()

    # Create target genome with chr1 of 300 bp
    setup_db(list(paste0(">chr1\n", paste(rep("A", 300), collapse = ""), "\n")))

    # Simple 1:1 mapping
    chain_file <- new_chain_file()
    write_chain_entry(chain_file, "source1", 300, "+", 0, 100, "chr1", 300, "+", 0, 100, 1)

    # Two separate source intervals that map to adjacent targets
    src_intervals <- data.frame(
        chrom = c("source1", "source1"),
        start = c(0, 50),
        end = c(50, 100),
        stringsAsFactors = FALSE
    )

    chain <- gintervals.load_chain(chain_file)

    # With canonic=TRUE: should still be 2 intervals (different intervalIDs)
    result_canonic <- gintervals.liftover(src_intervals, chain, canonic = TRUE)
    expect_equal(nrow(result_canonic), 2)

    # They should remain separate even though targets are adjacent [0-50] and [50-100]
    result_canonic <- result_canonic[order(result_canonic$start), ]
    expect_equal(result_canonic$start, c(0, 50))
    expect_equal(result_canonic$end, c(50, 100))
})

test_that("gintervals.liftover canonic handles non-adjacent blocks correctly", {
    local_db_state()

    # Create target genome with chr1 of 300 bp
    setup_db(list(paste0(">chr1\n", paste(rep("A", 300), collapse = ""), "\n")))

    # Chain with gap in BOTH source and target (non-adjacent targets)
    # Source [0-200] maps with:
    #   Block 1: 50 bp (source [0-50] -> target [0-50])
    #   Gap: 50 in target, 100 in source
    #   Block 2: 50 bp (source [150-200] -> target [100-150])
    # Target blocks [0-50] and [100-150] are NOT adjacent!
    chain_file <- new_chain_file()

    # Chain format: srcName = source (source1), tgtName = target (chr1)
    cat("chain 1000 source1 300 + 0 200 chr1 300 + 0 150 1\n", file = chain_file)
    cat("50\t100\t50\n", file = chain_file, append = TRUE) # gap: 100 in source, 50 in target
    cat("50\n\n", file = chain_file, append = TRUE)

    src_intervals <- data.frame(
        chrom = "source1",
        start = 0,
        end = 200,
        stringsAsFactors = FALSE
    )

    chain <- gintervals.load_chain(chain_file)

    # With canonic=TRUE: should still be 2 intervals (not adjacent)
    result_canonic <- gintervals.liftover(src_intervals, chain, canonic = TRUE)
    expect_equal(nrow(result_canonic), 2)

    result_canonic <- result_canonic[order(result_canonic$start), ]
    expect_equal(result_canonic$start, c(0, 100))
    expect_equal(result_canonic$end, c(50, 150))
})

test_that("gintervals.liftover canonic produces expected merged result", {
    # Note: Kent comparison removed because misha uses different gap format
    # (src_gap, tgt_gap) vs UCSC (tgt_gap, src_gap) in chain block lines.
    # Other tests already verify Kent compatibility with simple single-block chains.

    local_db_state()

    # Create target genome with chr1 of 300 bp
    setup_db(list(paste0(">chr1\n", paste(rep("A", 300), collapse = ""), "\n")))

    # Chain with adjacent target blocks
    chain_file <- new_chain_file()

    # Chain format: srcName = source (source1), tgtName = target (chr1)
    cat("chain 1000 source1 300 + 0 200 chr1 300 + 0 100 1\n", file = chain_file)
    cat("50\t100\t0\n", file = chain_file, append = TRUE) # gap: 100 in source, 0 in target
    cat("50\n\n", file = chain_file, append = TRUE)

    src_intervals <- data.frame(
        chrom = "source1",
        start = 0,
        end = 200,
        stringsAsFactors = FALSE
    )

    # Load chain and run misha liftover with canonic
    chain <- gintervals.load_chain(chain_file)
    misha_result <- gintervals.liftover(src_intervals, chain, canonic = TRUE)

    # Test that misha with canonic produces single merged interval
    expect_equal(nrow(misha_result), 1)
    expect_equal(misha_result$start, 0)
    expect_equal(misha_result$end, 100)
    expect_true("intervalID" %in% colnames(misha_result))
    expect_true("chain_id" %in% colnames(misha_result))
})

test_that("gintervals.liftover canonic merges multiple adjacent blocks", {
    local_db_state()

    # Create target genome with chr1 of 500 bp
    setup_db(list(paste0(">chr1\n", paste(rep("A", 500), collapse = ""), "\n")))

    # Chain with THREE adjacent target blocks
    # Source [0-350] maps with:
    #   Block 1: 50 bp (source [0-50] -> target [0-50])
    #   Gap: 0 in target, 100 in source
    #   Block 2: 50 bp (source [150-200] -> target [50-100])
    #   Gap: 0 in target, 100 in source
    #   Block 3: 50 bp (source [300-350] -> target [100-150])
    chain_file <- new_chain_file()

    # Chain format: srcName = source (source1), tgtName = target (chr1)
    cat("chain 1000 source1 500 + 0 350 chr1 500 + 0 150 1\n", file = chain_file)
    cat("50\t100\t0\n", file = chain_file, append = TRUE) # gap: 100 in source, 0 in target
    cat("50\t100\t0\n", file = chain_file, append = TRUE) # gap: 100 in source, 0 in target
    cat("50\n\n", file = chain_file, append = TRUE)

    src_intervals <- data.frame(
        chrom = "source1",
        start = 0,
        end = 350,
        stringsAsFactors = FALSE
    )

    chain <- gintervals.load_chain(chain_file)

    # Without canonic: should return 3 adjacent intervals
    result_no_canonic <- gintervals.liftover(src_intervals, chain)
    expect_equal(nrow(result_no_canonic), 3)

    # With canonic: should merge all 3 into single [0-150]
    result_canonic <- gintervals.liftover(src_intervals, chain, canonic = TRUE)
    expect_equal(nrow(result_canonic), 1)
    expect_equal(result_canonic$start, 0)
    expect_equal(result_canonic$end, 150)
})

test_that("gintervals.liftover default behavior unchanged (canonic=FALSE)", {
    local_db_state()

    # Create target genome with chr1 of 300 bp
    setup_db(list(paste0(">chr1\n", paste(rep("A", 300), collapse = ""), "\n")))

    # Chain with adjacent target blocks
    chain_file <- new_chain_file()

    # Chain format: srcName = source (source1), tgtName = target (chr1)
    cat("chain 1000 source1 300 + 0 200 chr1 300 + 0 100 1\n", file = chain_file)
    cat("50\t100\t0\n", file = chain_file, append = TRUE) # gap: 100 in source, 0 in target
    cat("50\n\n", file = chain_file, append = TRUE)

    src_intervals <- data.frame(
        chrom = "source1",
        start = 0,
        end = 200,
        stringsAsFactors = FALSE
    )

    chain <- gintervals.load_chain(chain_file)

    # Default (no canonic parameter): should behave like canonic=FALSE
    result_default <- gintervals.liftover(src_intervals, chain)
    expect_equal(nrow(result_default), 2)
    expect_true("intervalID" %in% colnames(result_default))

    # Explicit canonic=FALSE: same as default
    result_explicit <- gintervals.liftover(src_intervals, chain, canonic = FALSE)
    expect_equal(nrow(result_explicit), 2)
    expect_true("intervalID" %in% colnames(result_explicit))
})
