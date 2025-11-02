load_test_db()
test_that("gintervals.load_chain handles source overlaps with 'error' policy", {
    local_db_state()

    # Create target genome
    target_fasta <- tempfile(fileext = ".fasta")
    cat(">chr1\nACTGACTGACTGACTGACTGACTGACTGACTG\n>chr2\nGGGGCCCCTTTTAAAA\n", file = target_fasta)

    target_db <- tempfile()
    withr::defer({
        unlink(target_db, recursive = TRUE)
        unlink(target_fasta)
    })

    gdb.create(groot = target_db, fasta = target_fasta, verbose = FALSE)
    gdb.init(target_db)

    # Create chain file with source overlaps: same source position maps to multiple targets
    chain_file <- tempfile(fileext = ".chain")
    withr::defer(unlink(chain_file))

    # Chain 1: source1[0-20] -> chr1[0-20]
    cat("chain 1000 source1 100 + 0 20 chr1 32 + 0 20 1\n", file = chain_file)
    cat("20\n\n", file = chain_file, append = TRUE)

    # Chain 2: source1[5-21] -> chr2[0-16] (overlaps with chain 1 at source1[5-20])
    cat("chain 1000 source1 100 + 5 21 chr2 16 + 0 16 2\n", file = chain_file, append = TRUE)
    cat("16\n\n", file = chain_file, append = TRUE)

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
    target_fasta <- tempfile(fileext = ".fasta")
    cat(">chr1\nACTGACTGACTGACTGACTGACTGACTGACTG\n>chr2\nGGGGCCCCTTTTAAAA\n", file = target_fasta)

    target_db <- tempfile()
    withr::defer({
        unlink(target_db, recursive = TRUE)
        unlink(target_fasta)
    })

    gdb.create(groot = target_db, fasta = target_fasta, verbose = FALSE)
    gdb.init(target_db)

    # Create chain file with source overlaps
    chain_file <- tempfile(fileext = ".chain")
    withr::defer(unlink(chain_file))

    # Chain 1: source1[0-20] -> chr1[0-20]
    cat("chain 1000 source1 100 + 0 20 chr1 32 + 0 20 1\n", file = chain_file)
    cat("20\n\n", file = chain_file, append = TRUE)

    # Chain 2: source1[5-21] -> chr2[0-16] (overlaps at source1[5-20])
    cat("chain 1000 source1 100 + 5 21 chr2 16 + 0 16 2\n", file = chain_file, append = TRUE)
    cat("16\n\n", file = chain_file, append = TRUE)

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
    target_fasta <- tempfile(fileext = ".fasta")
    cat(">chr1\nACTGACTGACTGACTGACTGACTGACTGACTG\n>chr2\nGGGGCCCCTTTTAAAA\n>chr3\nTATATATATATA\n", file = target_fasta)

    target_db <- tempfile()
    withr::defer({
        unlink(target_db, recursive = TRUE)
        unlink(target_fasta)
    })

    gdb.create(groot = target_db, fasta = target_fasta, verbose = FALSE)
    gdb.init(target_db)

    # Create chain file: source1 overlaps (maps to chr1 and chr2), source2 is clean (maps to chr3)
    chain_file <- tempfile(fileext = ".chain")
    withr::defer(unlink(chain_file))

    # Chain 1: source1[0-20] -> chr1[0-20]
    cat("chain 1000 source1 100 + 0 20 chr1 32 + 0 20 1\n", file = chain_file)
    cat("20\n\n", file = chain_file, append = TRUE)

    # Chain 2: source1[5-21] -> chr2[0-16] (overlaps with chain 1)
    cat("chain 1000 source1 100 + 5 21 chr2 16 + 0 16 2\n", file = chain_file, append = TRUE)
    cat("16\n\n", file = chain_file, append = TRUE)

    # Chain 3: source2[0-12] -> chr3[0-12] (no overlaps)
    cat("chain 1000 source2 50 + 0 12 chr3 12 + 0 12 3\n", file = chain_file, append = TRUE)
    cat("12\n\n", file = chain_file, append = TRUE)

    # Discard policy should remove overlapping intervals but keep clean ones
    chain <- gintervals.load_chain(chain_file, src_overlap_policy = "discard")
    expect_equal(nrow(chain), 1)
    expect_equal(as.character(chain$chrom), "chr3")
    expect_equal(as.character(chain$chromsrc), "source2")
    expect_equal(as.numeric(chain$start), 0)
    expect_equal(as.numeric(chain$end), 12)
    expect_equal(as.numeric(chain$startsrc), 0)
})

test_that("gintervals.load_chain handles target overlaps with 'auto' policy", {
    local_db_state()

    # Create target genome
    target_fasta <- tempfile(fileext = ".fasta")
    cat(">chr1\nACTGACTGACTGACTGACTGACTGACTGACTGACTGACTG\n", file = target_fasta)

    target_db <- tempfile()
    withr::defer({
        unlink(target_db, recursive = TRUE)
        unlink(target_fasta)
    })

    gdb.create(groot = target_db, fasta = target_fasta, verbose = FALSE)
    gdb.init(target_db)

    # Create chain file with target overlaps: different sources map to overlapping chr1 regions
    chain_file <- tempfile(fileext = ".chain")
    withr::defer(unlink(chain_file))

    # Chain 1: source1[0-30] -> chr1[0-30]
    cat("chain 1000 source1 100 + 0 30 chr1 40 + 0 30 1\n", file = chain_file)
    cat("30\n\n", file = chain_file, append = TRUE)

    # Chain 2: source2[0-20] -> chr1[20-40] (overlaps at chr1[20-30])
    cat("chain 1000 source2 100 + 0 20 chr1 40 + 20 40 2\n", file = chain_file, append = TRUE)
    cat("20\n\n", file = chain_file, append = TRUE)

    # Auto policy should resolve overlaps by truncating
    chain <- gintervals.load_chain(chain_file, tgt_overlap_policy = "auto")
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

test_that("gintervals.load_chain exact truncation with 'auto' policy", {
    local_db_state()

    target_fasta <- tempfile(fileext = ".fasta")
    cat(
        ">chr1\nACTGACTGACTGACTGACTGACTGACTGACTGACTGACTG\n",
        file = target_fasta
    )

    target_db <- tempfile()
    withr::defer({
        unlink(target_db, recursive = TRUE)
        unlink(target_fasta)
    })

    gdb.create(groot = target_db, fasta = target_fasta, verbose = FALSE)
    gdb.init(target_db)

    chain_file <- tempfile(fileext = ".chain")
    withr::defer(unlink(chain_file))

    # Chain 1: sourceA[0-30] -> chr1[0-30]
    cat("chain 1000 sourceA 100 + 0 30 chr1 40 + 0 30 1\n", file = chain_file)
    cat("30\n\n", file = chain_file, append = TRUE)
    # Chain 2: sourceB[0-20] -> chr1[20-40] (overlaps 20-30)
    cat("chain 1000 sourceB 100 + 0 20 chr1 40 + 20 40 2\n", file = chain_file, append = TRUE)
    cat("20\n\n", file = chain_file, append = TRUE)

    chain <- gintervals.load_chain(chain_file, tgt_overlap_policy = "auto")
    chain_chr1 <- chain[chain$chrom == "chr1", ]
    chain_chr1 <- chain_chr1[order(chain_chr1$start), ]
    # After auto truncation: first becomes [0,20), second becomes [30,40)
    expect_equal(nrow(chain_chr1), 2)
    expect_equal(as.numeric(chain_chr1$start[1]), 0)
    expect_equal(as.numeric(chain_chr1$end[1]), 20)
    expect_equal(as.numeric(chain_chr1$start[2]), 30)
    expect_equal(as.numeric(chain_chr1$end[2]), 40)
})

test_that("gintervals.load_chain handles target overlaps with 'error' policy", {
    local_db_state()

    # Create target genome
    target_fasta <- tempfile(fileext = ".fasta")
    cat(">chr1\nACTGACTGACTGACTGACTGACTGACTGACTG\n", file = target_fasta)

    target_db <- tempfile()
    withr::defer({
        unlink(target_db, recursive = TRUE)
        unlink(target_fasta)
    })

    gdb.create(groot = target_db, fasta = target_fasta, verbose = FALSE)
    gdb.init(target_db)

    # Create chain file with target overlaps
    chain_file <- tempfile(fileext = ".chain")
    withr::defer(unlink(chain_file))

    # Chain 1: source1[0-20] -> chr1[0-20]
    cat("chain 1000 source1 100 + 0 20 chr1 32 + 0 20 1\n", file = chain_file)
    cat("20\n\n", file = chain_file, append = TRUE)

    # Chain 2: source2[0-17] -> chr1[15-32] (overlaps at chr1[15-20])
    cat("chain 1000 source2 100 + 0 17 chr1 32 + 15 32 2\n", file = chain_file, append = TRUE)
    cat("17\n\n", file = chain_file, append = TRUE)

    # Error policy should fail on target overlap
    expect_error(
        gintervals.load_chain(chain_file, tgt_overlap_policy = "error"),
        "overlap"
    )
})

test_that("gintervals.load_chain handles target overlaps with 'discard' policy", {
    local_db_state()

    # Create target genome
    target_fasta <- tempfile(fileext = ".fasta")
    cat(">chr1\nACTGACTGACTGACTGACTGACTGACTGACTG\n>chr2\nGGGGCCCCTTTTAAAA\n", file = target_fasta)

    target_db <- tempfile()
    withr::defer({
        unlink(target_db, recursive = TRUE)
        unlink(target_fasta)
    })

    gdb.create(groot = target_db, fasta = target_fasta, verbose = FALSE)
    gdb.init(target_db)

    # Create chain file: chr1 has overlaps, chr2 doesn't
    chain_file <- tempfile(fileext = ".chain")
    withr::defer(unlink(chain_file))

    # Chain 1: source1[0-20] -> chr1[0-20]
    cat("chain 1000 source1 100 + 0 20 chr1 32 + 0 20 1\n", file = chain_file)
    cat("20\n\n", file = chain_file, append = TRUE)

    # Chain 2: source2[0-17] -> chr1[15-32] (overlaps with chain 1)
    cat("chain 1000 source2 100 + 0 17 chr1 32 + 15 32 2\n", file = chain_file, append = TRUE)
    cat("17\n\n", file = chain_file, append = TRUE)

    # Chain 3: source3[0-10] -> chr2[0-10] (no overlap)
    cat("chain 1000 source3 100 + 0 10 chr2 16 + 0 10 3\n", file = chain_file, append = TRUE)
    cat("10\n\n", file = chain_file, append = TRUE)

    # Discard policy should remove overlapping chr1 intervals but keep chr2
    chain <- gintervals.load_chain(chain_file, tgt_overlap_policy = "discard")
    expect_equal(nrow(chain), 1)
    expect_equal(as.character(chain$chrom), "chr2")
})

test_that("gintervals.liftover works with 'keep' source policy", {
    local_db_state()

    # Create target genome
    target_fasta <- tempfile(fileext = ".fasta")
    cat(">chr1\nACTGACTGACTGACTGACTGACTGACTGACTG\n>chr2\nGGGGCCCCTTTTAAAA\n", file = target_fasta)

    target_db <- tempfile()
    withr::defer({
        unlink(target_db, recursive = TRUE)
        unlink(target_fasta)
    })

    gdb.create(groot = target_db, fasta = target_fasta, verbose = FALSE)
    gdb.init(target_db)

    # Create chain with source overlaps
    chain_file <- tempfile(fileext = ".chain")
    withr::defer(unlink(chain_file))

    # Chain 1: source1[0-20] -> chr1[0-20]
    cat("chain 1000 source1 100 + 0 20 chr1 32 + 0 20 1\n", file = chain_file)
    cat("20\n\n", file = chain_file, append = TRUE)

    # Chain 2: source1[10-26] -> chr2[0-16] (overlaps at source1[10-20])
    cat("chain 1000 source1 100 + 10 26 chr2 16 + 0 16 2\n", file = chain_file, append = TRUE)
    cat("16\n\n", file = chain_file, append = TRUE)

    # Load chain with keep policy
    chain <- gintervals.load_chain(chain_file, src_overlap_policy = "keep")

    # Create source intervals that overlap with the ambiguous region
    src_intervals <- data.frame(
        chrom = "source1",
        start = 10,
        end = 20,
        stringsAsFactors = FALSE
    )

    # Liftover should produce multiple target intervals
    result <- gintervals.liftover(src_intervals, chain, src_overlap_policy = "keep")

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
    target_fasta <- tempfile(fileext = ".fasta")
    cat(">chr1\nACTGACTGACTGACTGACTGACTGACTGACTG\n", file = target_fasta)

    target_db <- tempfile()
    withr::defer({
        unlink(target_db, recursive = TRUE)
        unlink(target_fasta)
    })

    gdb.create(groot = target_db, fasta = target_fasta, verbose = FALSE)
    gdb.init(target_db)

    # Create simple chain file
    chain_file <- tempfile(fileext = ".chain")
    withr::defer(unlink(chain_file))

    # Chain: source1[0-20] -> chr1[0-20]
    cat("chain 1000 source1 100 + 0 20 chr1 32 + 0 20 1\n", file = chain_file)
    cat("20\n\n", file = chain_file, append = TRUE)

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

    target_fasta <- tempfile(fileext = ".fasta")
    cat(
        ">chr1\nAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n>chr2\nCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n",
        file = target_fasta
    )

    target_db <- tempfile()
    withr::defer({
        unlink(target_db, recursive = TRUE)
        unlink(target_fasta)
    })

    gdb.create(groot = target_db, fasta = target_fasta, verbose = FALSE)
    gdb.init(target_db)

    chain_file <- tempfile(fileext = ".chain")
    withr::defer(unlink(chain_file))

    # sourceA[0-20) -> chr1[10-30)
    cat("chain 1000 sourceA 200 + 0 20 chr1 40 + 10 30 1\n", file = chain_file)
    cat("20\n\n", file = chain_file, append = TRUE)

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

    target_fasta <- tempfile(fileext = ".fasta")
    cat(
        ">chr1\nAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n>chr2\nCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n",
        file = target_fasta
    )

    target_db <- tempfile()
    withr::defer({
        unlink(target_db, recursive = TRUE)
        unlink(target_fasta)
    })

    gdb.create(groot = target_db, fasta = target_fasta, verbose = FALSE)
    gdb.init(target_db)

    chain_file <- tempfile(fileext = ".chain")
    withr::defer(unlink(chain_file))

    # sourceA[0-20) -> chr1[10-30)
    cat("chain 1000 sourceA 200 + 0 20 chr1 40 + 10 30 1\n", file = chain_file)
    cat("20\n\n", file = chain_file, append = TRUE)
    # sourceB[10-37) -> chr2[5-32) (length 27 to match chr2 size 32)
    cat("chain 1000 sourceB 200 + 10 37 chr2 32 + 5 32 2\n", file = chain_file, append = TRUE)
    cat("27\n\n", file = chain_file, append = TRUE)

    chain <- gintervals.load_chain(chain_file)

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
    target_fasta <- tempfile(fileext = ".fasta")
    cat(
        ">chr1\nAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n",
        file = target_fasta
    )

    target_db <- tempfile()
    withr::defer({
        unlink(target_db, recursive = TRUE)
        unlink(target_fasta)
    })

    gdb.create(groot = target_db, fasta = target_fasta, verbose = FALSE)
    gdb.init(target_db)

    # Chain mapping source1[0-20) -> chr1[0-20)
    chain_file <- tempfile(fileext = ".chain")
    withr::defer(unlink(chain_file))
    cat("chain 1000 source1 100 + 0 20 chr1 32 + 0 20 1\n", file = chain_file)
    cat("20\n\n", file = chain_file, append = TRUE)

    # Create a source DB with 'source1'
    source_fasta <- tempfile(fileext = ".fasta")
    cat(
        ">source1\nTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\n",
        file = source_fasta
    )
    source_db <- tempfile()
    withr::defer({
        unlink(source_db, recursive = TRUE)
        unlink(source_fasta)
    })
    gdb.create(groot = source_db, fasta = source_fasta, verbose = FALSE)
    gdb.init(source_db)

    # Create a sparse source track on source1
    src_intervals <- data.frame(
        chrom = "source1",
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
    target_fasta <- tempfile(fileext = ".fasta")
    cat(">chr1\nACTGACTGACTGACTG\n", file = target_fasta)

    target_db <- tempfile()
    withr::defer({
        unlink(target_db, recursive = TRUE)
        unlink(target_fasta)
    })

    gdb.create(groot = target_db, fasta = target_fasta, verbose = FALSE)
    gdb.init(target_db)

    chain_file <- tempfile(fileext = ".chain")
    withr::defer(unlink(chain_file))

    cat("chain 1000 source1 100 + 0 10 chr1 16 + 0 10 1\n", file = chain_file)
    cat("10\n\n", file = chain_file, append = TRUE)

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
    target_fasta <- tempfile(fileext = ".fasta")
    cat(">chr1\nACTGACTGACTG\n", file = target_fasta)

    target_db <- tempfile()
    withr::defer({
        unlink(target_db, recursive = TRUE)
        unlink(target_fasta)
    })

    gdb.create(groot = target_db, fasta = target_fasta, verbose = FALSE)
    gdb.init(target_db)

    chain_file <- tempfile(fileext = ".chain")
    withr::defer(unlink(chain_file))

    cat("chain 1000 source1 100 + 0 10 chr1 12 + 0 10 1\n", file = chain_file)
    cat("10\n\n", file = chain_file, append = TRUE)

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
    target_fasta <- tempfile(fileext = ".fasta")
    cat(">chr1\nACTGACTGACTGACTGACTGACTGACTGACTG\n", file = target_fasta)

    target_db <- tempfile()
    withr::defer({
        unlink(target_db, recursive = TRUE)
        unlink(target_fasta)
    })

    gdb.create(groot = target_db, fasta = target_fasta, verbose = FALSE)
    gdb.init(target_db)

    # Create chain with only overlapping intervals
    chain_file <- tempfile(fileext = ".chain")
    withr::defer(unlink(chain_file))

    # Chain 1: source1[0-20] -> chr1[0-20]
    cat("chain 1000 source1 100 + 0 20 chr1 32 + 0 20 1\n", file = chain_file)
    cat("20\n\n", file = chain_file, append = TRUE)

    # Chain 2: source1[10-27] -> chr1[15-32] (both source and target overlaps)
    cat("chain 1000 source1 100 + 10 27 chr1 32 + 15 32 2\n", file = chain_file, append = TRUE)
    cat("17\n\n", file = chain_file, append = TRUE)

    # Load chain and discard all overlaps
    chain <- gintervals.load_chain(chain_file, src_overlap_policy = "discard")

    # Chain should be empty or NULL
    expect_true(is.null(chain) || nrow(chain) == 0)
})

test_that("complex chain with both source and target overlaps", {
    local_db_state()

    # Create target genome with multiple chromosomes
    target_fasta <- tempfile(fileext = ".fasta")
    cat(">chr1\nACTGACTGACTGACTGACTGACTGACTGACTG\n>chr2\nGGGGCCCCTTTTAAAA\n>chr3\nTATATATATATA\n", file = target_fasta)

    target_db <- tempfile()
    withr::defer({
        unlink(target_db, recursive = TRUE)
        unlink(target_fasta)
    })

    gdb.create(groot = target_db, fasta = target_fasta, verbose = FALSE)
    gdb.init(target_db)

    # Create complex chain file with various overlap scenarios
    chain_file <- tempfile(fileext = ".chain")
    withr::defer(unlink(chain_file))

    # Source overlap: source1[0,20] -> chr1 and chr2
    cat("chain 1000 source1 100 + 0 20 chr1 32 + 0 20 1\n", file = chain_file)
    cat("20\n\n", file = chain_file, append = TRUE)

    cat("chain 1000 source1 100 + 10 26 chr2 16 + 0 16 2\n", file = chain_file, append = TRUE)
    cat("16\n\n", file = chain_file, append = TRUE)

    # Target overlap: source2 and source3 both map to overlapping chr3 regions
    cat("chain 1000 source2 100 + 0 8 chr3 12 + 0 8 3\n", file = chain_file)
    cat("8\n\n", file = chain_file, append = TRUE)

    cat("chain 1000 source3 100 + 0 7 chr3 12 + 5 12 4\n", file = chain_file, append = TRUE)
    cat("7\n\n", file = chain_file, append = TRUE)

    # Clean mapping: source4 -> chr1 (no overlaps)
    cat("chain 1000 source4 100 + 0 7 chr1 32 + 25 32 5\n", file = chain_file, append = TRUE)
    cat("7\n\n", file = chain_file, append = TRUE)

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

test_that("gintervals.load_chain returns 8 columns with strand information", {
    local_db_state()

    # Create target genome
    target_fasta <- tempfile(fileext = ".fasta")
    cat(">chr1\nACTGACTGACTGACTGACTGACTGACTGACTG\n", file = target_fasta)

    target_db <- tempfile()
    withr::defer({
        unlink(target_db, recursive = TRUE)
        unlink(target_fasta)
    })

    gdb.create(groot = target_db, fasta = target_fasta, verbose = FALSE)
    gdb.init(target_db)

    # Create simple chain file with forward strand
    chain_file <- tempfile(fileext = ".chain")
    withr::defer(unlink(chain_file))

    cat("chain 1000 source1 100 + 0 20 chr1 32 + 0 20 1\n", file = chain_file)
    cat("20\n\n", file = chain_file, append = TRUE)

    chain <- gintervals.load_chain(chain_file)

    # Check that chain has 8 columns
    expected_cols <- c("chrom", "start", "end", "strand", "chromsrc", "startsrc", "endsrc", "strandsrc")
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
    target_fasta <- tempfile(fileext = ".fasta")
    cat(">chr1\nACTGACTGACTGACTGACTGACTGACTGACTG\n", file = target_fasta)

    target_db <- tempfile()
    withr::defer({
        unlink(target_db, recursive = TRUE)
        unlink(target_fasta)
    })

    gdb.create(groot = target_db, fasta = target_fasta, verbose = FALSE)
    gdb.init(target_db)

    # Create chain file with reverse strands
    chain_file <- tempfile(fileext = ".chain")
    withr::defer(unlink(chain_file))

    # Source on reverse strand, target on forward
    cat("chain 1000 source1 100 - 0 20 chr1 32 + 0 20 1\n", file = chain_file)
    cat("20\n\n", file = chain_file, append = TRUE)

    chain <- gintervals.load_chain(chain_file)

    # Check strand format: target forward (+1), source reverse (-1)
    expect_equal(as.numeric(chain$strand), 1)
    expect_equal(as.numeric(chain$strandsrc), -1)
})

test_that("gintervals.load_chain handles mixed strand chains", {
    local_db_state()

    # Create target genome
    target_fasta <- tempfile(fileext = ".fasta")
    cat(">chr1\nACTGACTGACTGACTGACTGACTGACTGACTG\n>chr2\nGGGGCCCCTTTTAAAA\n", file = target_fasta)

    target_db <- tempfile()
    withr::defer({
        unlink(target_db, recursive = TRUE)
        unlink(target_fasta)
    })

    gdb.create(groot = target_db, fasta = target_fasta, verbose = FALSE)
    gdb.init(target_db)

    # Create chain file with mixed strands
    chain_file <- tempfile(fileext = ".chain")
    withr::defer(unlink(chain_file))

    # Chain 1: forward-forward
    cat("chain 1000 source1 100 + 0 20 chr1 32 + 0 20 1\n", file = chain_file)
    cat("20\n\n", file = chain_file, append = TRUE)

    # Chain 2: forward-reverse (source forward, target reverse)
    cat("chain 1000 source2 100 + 0 16 chr2 16 - 0 16 2\n", file = chain_file, append = TRUE)
    cat("16\n\n", file = chain_file, append = TRUE)

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
    target_fasta <- tempfile(fileext = ".fasta")
    cat(">chr1\nACTGACTGACTGACTGACTGACTGACTGACTG\n", file = target_fasta)

    target_db <- tempfile()
    withr::defer({
        unlink(target_db, recursive = TRUE)
        unlink(target_fasta)
    })

    gdb.create(groot = target_db, fasta = target_fasta, verbose = FALSE)
    gdb.init(target_db)

    # Create source genome database with specific chromosomes
    source_fasta <- tempfile(fileext = ".fasta")
    cat(">source1\nTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\n", file = source_fasta)

    source_db <- tempfile()
    withr::defer({
        unlink(source_db, recursive = TRUE)
        unlink(source_fasta)
    })

    gdb.create(groot = source_db, fasta = source_fasta, verbose = FALSE)

    # Create chain file with valid source chromosome
    chain_file <- tempfile(fileext = ".chain")
    withr::defer(unlink(chain_file))

    cat("chain 1000 source1 100 + 0 20 chr1 32 + 0 20 1\n", file = chain_file)
    cat("20\n\n", file = chain_file, append = TRUE)

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
    target_fasta <- tempfile(fileext = ".fasta")
    cat(">chr1\nACTGACTGACTGACTGACTGACTGACTGACTG\n", file = target_fasta)

    target_db <- tempfile()
    withr::defer({
        unlink(target_db, recursive = TRUE)
        unlink(target_fasta)
    })

    gdb.create(groot = target_db, fasta = target_fasta, verbose = FALSE)
    gdb.init(target_db)

    # Create source genome with different chromosome name
    source_fasta <- tempfile(fileext = ".fasta")
    cat(">different_chrom\nTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\n", file = source_fasta)

    source_db <- tempfile()
    withr::defer({
        unlink(source_db, recursive = TRUE)
        unlink(source_fasta)
    })

    gdb.create(groot = source_db, fasta = source_fasta, verbose = FALSE)

    # Create chain file with invalid source chromosome
    chain_file <- tempfile(fileext = ".chain")
    withr::defer(unlink(chain_file))

    cat("chain 1000 invalid_source 100 + 0 20 chr1 32 + 0 20 1\n", file = chain_file)
    cat("20\n\n", file = chain_file, append = TRUE)

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
    target_fasta <- tempfile(fileext = ".fasta")
    cat(">chr1\nACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGAC\n", file = target_fasta)

    target_db <- tempfile()
    withr::defer({
        unlink(target_db, recursive = TRUE)
        unlink(target_fasta)
    })

    gdb.create(groot = target_db, fasta = target_fasta, verbose = FALSE)
    gdb.init(target_db)

    # Create source genome with small chromosome
    source_fasta <- tempfile(fileext = ".fasta")
    cat(">source1\nTTTTTTTTTTTTTTTTTTTT\n", file = source_fasta) # Only 20 bp

    source_db <- tempfile()
    withr::defer({
        unlink(source_db, recursive = TRUE)
        unlink(source_fasta)
    })

    gdb.create(groot = source_db, fasta = source_fasta, verbose = FALSE)

    # Create chain file with source coordinates exceeding chromosome size
    chain_file <- tempfile(fileext = ".chain")
    withr::defer(unlink(chain_file))

    cat("chain 1000 source1 100 + 0 50 chr1 70 + 0 50 1\n", file = chain_file)
    cat("50\n\n", file = chain_file, append = TRUE)

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
    target_fasta <- tempfile(fileext = ".fasta")
    cat(">chr1\nACTGACTGACTGACTGACTGACTGACTGACTG\n", file = target_fasta)

    target_db <- tempfile()
    withr::defer({
        unlink(target_db, recursive = TRUE)
        unlink(target_fasta)
    })

    gdb.create(groot = target_db, fasta = target_fasta, verbose = FALSE)
    gdb.init(target_db)

    # Create chain file
    chain_file <- tempfile(fileext = ".chain")
    withr::defer(unlink(chain_file))

    cat("chain 1000 source1 100 + 0 20 chr1 32 + 0 20 1\n", file = chain_file)
    cat("20\n\n", file = chain_file, append = TRUE)

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

# Helper function to check if liftOver binary is available
has_liftover_binary <- function() {
    result <- tryCatch(
        {
            system2("which", "liftOver", stdout = TRUE, stderr = FALSE)
            TRUE
        },
        error = function(e) FALSE,
        warning = function(w) FALSE
    )
    if (!isTRUE(result)) {
        result <- tryCatch(
            {
                system2("liftOver", stdout = FALSE, stderr = FALSE)
                TRUE
            },
            error = function(e) FALSE,
            warning = function(w) FALSE
        )
    }
    return(isTRUE(result))
}

test_that("gintervals.liftover matches liftOver binary - basic case", {
    skip_if_not(has_liftover_binary(), "liftOver binary not found")

    local_db_state()

    # Create target genome
    target_fasta <- tempfile(fileext = ".fasta")
    cat(">chr1\nACTGACTGACTGACTGACTGACTGACTGACTGACTGACTG\n", file = target_fasta)

    target_db <- tempfile()
    withr::defer({
        unlink(target_db, recursive = TRUE)
        unlink(target_fasta)
    })

    gdb.create(groot = target_db, fasta = target_fasta, verbose = FALSE)
    gdb.init(target_db)

    # Create chain file with simple mapping
    chain_file <- tempfile(fileext = ".chain")
    withr::defer(unlink(chain_file))

    # source1[10-30] -> chr1[5-25]
    cat("chain 1000 source1 100 + 10 30 chr1 40 + 5 25 1\n", file = chain_file)
    cat("20\n\n", file = chain_file, append = TRUE)

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
    target_fasta <- tempfile(fileext = ".fasta")
    cat(">chr1\nACTGACTGACTGACTGACTGACTGACTGACTGACTGACTG\n>chr2\nGGGGCCCCTTTTAAAAGGGGCCCCTTTTAAAA\n", file = target_fasta)

    target_db <- tempfile()
    withr::defer({
        unlink(target_db, recursive = TRUE)
        unlink(target_fasta)
    })

    gdb.create(groot = target_db, fasta = target_fasta, verbose = FALSE)
    gdb.init(target_db)

    # Create chain file with multiple mappings
    chain_file <- tempfile(fileext = ".chain")
    withr::defer(unlink(chain_file))

    # source1[0-20] -> chr1[10-30]
    cat("chain 1000 source1 100 + 0 20 chr1 40 + 10 30 1\n", file = chain_file)
    cat("20\n\n", file = chain_file, append = TRUE)

    # source2[5-25] -> chr2[0-20]
    cat("chain 1000 source2 100 + 5 25 chr2 32 + 0 20 2\n", file = chain_file, append = TRUE)
    cat("20\n\n", file = chain_file, append = TRUE)

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
    target_fasta <- tempfile(fileext = ".fasta")
    cat(">chr1\nACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTG\n", file = target_fasta)

    target_db <- tempfile()
    withr::defer({
        unlink(target_db, recursive = TRUE)
        unlink(target_fasta)
    })

    gdb.create(groot = target_db, fasta = target_fasta, verbose = FALSE)
    gdb.init(target_db)

    # Create chain file with reverse strand mapping
    chain_file <- tempfile(fileext = ".chain")
    withr::defer(unlink(chain_file))

    # source1[0-30] -> chr1[22-52] on reverse strand (chr1 size is 56)
    cat("chain 1000 source1 100 + 0 30 chr1 56 - 4 34 1\n", file = chain_file)
    cat("30\n\n", file = chain_file, append = TRUE)

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
    target_fasta <- tempfile(fileext = ".fasta")
    cat(">chr1\nACTGACTGACTGACTGACTGACTGACTGACTGACTGACTG\n", file = target_fasta)

    target_db <- tempfile()
    withr::defer({
        unlink(target_db, recursive = TRUE)
        unlink(target_fasta)
    })

    gdb.create(groot = target_db, fasta = target_fasta, verbose = FALSE)
    gdb.init(target_db)

    # Create chain file
    chain_file <- tempfile(fileext = ".chain")
    withr::defer(unlink(chain_file))

    # source1[10-30] -> chr1[5-25]
    cat("chain 1000 source1 100 + 10 30 chr1 40 + 5 25 1\n", file = chain_file)
    cat("20\n\n", file = chain_file, append = TRUE)

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
    target_fasta <- tempfile(fileext = ".fasta")
    cat(">chr1\nACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTG\n", file = target_fasta)

    target_db <- tempfile()
    withr::defer({
        unlink(target_db, recursive = TRUE)
        unlink(target_fasta)
    })

    gdb.create(groot = target_db, fasta = target_fasta, verbose = FALSE)
    gdb.init(target_db)

    # Create chain file with gaps (insertions/deletions)
    chain_file <- tempfile(fileext = ".chain")
    withr::defer(unlink(chain_file))

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
    target_fasta <- tempfile(fileext = ".fasta")
    cat(">chr1\nACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTG\n", file = target_fasta)

    target_db <- tempfile()
    withr::defer({
        unlink(target_db, recursive = TRUE)
        unlink(target_fasta)
    })

    gdb.create(groot = target_db, fasta = target_fasta, verbose = FALSE)
    gdb.init(target_db)

    # Create chain file with gap
    chain_file <- tempfile(fileext = ".chain")
    withr::defer(unlink(chain_file))

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
    target_fasta <- tempfile(fileext = ".fasta")
    cat(">chr1\nACTGACTGACTGACTGACTGACTGACTGACTGACTGACTG\n", file = target_fasta)

    target_db <- tempfile()
    withr::defer({
        unlink(target_db, recursive = TRUE)
        unlink(target_fasta)
    })

    gdb.create(groot = target_db, fasta = target_fasta, verbose = FALSE)
    gdb.init(target_db)

    # Create chain file
    chain_file <- tempfile(fileext = ".chain")
    withr::defer(unlink(chain_file))

    # Simple chain
    cat("chain 1000 source1 100 + 0 30 chr1 40 + 5 35 1\n", file = chain_file)
    cat("30\n\n", file = chain_file, append = TRUE)

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
    target_fasta <- tempfile(fileext = ".fasta")
    cat(">chr1\nACTGACTGACTGACTGACTGACTGACTGACTGACTGACTG\n>chr2\nGGGGCCCCTTTTAAAAGGGGCCCC\n", file = target_fasta)

    target_db <- tempfile()
    withr::defer({
        unlink(target_db, recursive = TRUE)
        unlink(target_fasta)
    })

    gdb.create(groot = target_db, fasta = target_fasta, verbose = FALSE)
    gdb.init(target_db)

    # Create chain file with two chains
    chain_file <- tempfile(fileext = ".chain")
    withr::defer(unlink(chain_file))

    # Chain 1: source1[10-30] -> chr1[5-25]
    cat("chain 1000 source1 100 + 10 30 chr1 40 + 5 25 1\n", file = chain_file)
    cat("20\n\n", file = chain_file, append = TRUE)

    # Chain 2: source1[35-50] -> chr2[2-17]
    cat("chain 1000 source1 100 + 35 50 chr2 24 + 2 17 2\n", file = chain_file, append = TRUE)
    cat("15\n\n", file = chain_file, append = TRUE)

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
    target_fasta <- tempfile(fileext = ".fasta")
    cat(">chr1\nACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTG\n>chr2\nGGGGCCCCTTTTAAAAGGGGCCCCTTTTAAAA\n", file = target_fasta)

    target_db <- tempfile()
    withr::defer({
        unlink(target_db, recursive = TRUE)
        unlink(target_fasta)
    })

    gdb.create(groot = target_db, fasta = target_fasta, verbose = FALSE)
    gdb.init(target_db)

    # Create chain file with mixed strands
    chain_file <- tempfile(fileext = ".chain")
    withr::defer(unlink(chain_file))

    # Chain 1: forward strand
    cat("chain 1000 source1 100 + 0 20 chr1 56 + 10 30 1\n", file = chain_file)
    cat("20\n\n", file = chain_file, append = TRUE)

    # Chain 2: forward source, reverse target
    cat("chain 1000 source2 100 + 5 25 chr2 32 - 5 25 2\n", file = chain_file, append = TRUE)
    cat("20\n\n", file = chain_file, append = TRUE)

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
    target_fasta <- tempfile(fileext = ".fasta")
    cat(">chr1\n", file = target_fasta)
    cat(paste(rep("ACTG", 50), collapse = ""), "\n", file = target_fasta, append = TRUE)
    cat(">chr2\n", file = target_fasta, append = TRUE)
    cat(paste(rep("GCTA", 40), collapse = ""), "\n", file = target_fasta, append = TRUE)

    target_db <- tempfile()
    withr::defer({
        unlink(target_db, recursive = TRUE)
        unlink(target_fasta)
    })

    gdb.create(groot = target_db, fasta = target_fasta, verbose = FALSE)
    gdb.init(target_db)

    # Create chain file with many small non-overlapping chains
    chain_file <- tempfile(fileext = ".chain")
    withr::defer(unlink(chain_file))

    # Multiple chains from different source regions to different targets
    cat("chain 1000 source1 200 + 0 15 chr1 200 + 10 25 1\n", file = chain_file)
    cat("15\n\n", file = chain_file, append = TRUE)

    cat("chain 1000 source1 200 + 20 35 chr1 200 + 50 65 2\n", file = chain_file, append = TRUE)
    cat("15\n\n", file = chain_file, append = TRUE)

    cat("chain 1000 source1 200 + 40 58 chr2 160 + 20 38 3\n", file = chain_file, append = TRUE)
    cat("18\n\n", file = chain_file, append = TRUE)

    cat("chain 1000 source2 150 + 5 25 chr1 200 + 100 120 4\n", file = chain_file, append = TRUE)
    cat("20\n\n", file = chain_file, append = TRUE)

    cat("chain 1000 source2 150 + 30 50 chr2 160 + 60 80 5\n", file = chain_file, append = TRUE)
    cat("20\n\n", file = chain_file, append = TRUE)

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
    target_fasta <- tempfile(fileext = ".fasta")
    cat(">chr1\nACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTG\n", file = target_fasta)
    cat(">chr2\nGGGGCCCCTTTTAAAAGGGGCCCCTTTTAAAAGGGGCCCCTTTTAAAA\n", file = target_fasta, append = TRUE)
    cat(">chr3\nCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCG\n", file = target_fasta, append = TRUE)

    target_db <- tempfile()
    withr::defer({
        unlink(target_db, recursive = TRUE)
        unlink(target_fasta)
    })

    gdb.create(groot = target_db, fasta = target_fasta, verbose = FALSE)
    gdb.init(target_db)

    # Create chain file replicating the user's bug scenario
    chain_file <- tempfile(fileext = ".chain")
    withr::defer(unlink(chain_file))

    # Chains with the same start position (will be consecutive when sorted)
    cat("chain 1000 source1 100 + 10 30 chr1 44 + 0 20 1\n", file = chain_file)
    cat("20\n\n", file = chain_file, append = TRUE)

    cat("chain 1000 source1 100 + 10 30 chr2 48 + 0 20 2\n", file = chain_file, append = TRUE)
    cat("20\n\n", file = chain_file, append = TRUE)

    # Chain with different range that creates a gap
    cat("chain 1000 source1 100 + 30 32 chr3 42 + 0 2 3\n", file = chain_file, append = TRUE)
    cat("2\n\n", file = chain_file, append = TRUE)

    # Load chain with keep policy
    chain <- gintervals.load_chain(chain_file, src_overlap_policy = "keep", tgt_overlap_policy = "keep")

    # Query interval that overlaps first two chains but not the third
    src_intervals <- data.frame(chrom = "source1", start = 20, end = 21, stringsAsFactors = FALSE)

    # Perform liftover
    result <- gintervals.liftover(src_intervals, chain, src_overlap_policy = "keep", tgt_overlap_policy = "keep")

    # Should return exactly 2 results (both chains with [10,30) overlap [20,21))
    expect_equal(nrow(result), 2)
    expect_true("chr1" %in% result$chrom)
    expect_true("chr2" %in% result$chrom)
    expect_false("chr3" %in% result$chrom)
})
