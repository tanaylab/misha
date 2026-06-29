# Regression: map_interval's back() upper-bound guard used the last-by-start
# chain's own reach as if it bounded all chains on the source chromosome. With
# overlapping source chains (src_overlap_policy="keep") a wider EARLIER chain can
# reach past the query while the last-by-start chain ends before it, so the guard
# dropped the mapping. Same family as the 5.11.4 / 5.11.7 fixes, on the upper
# bound. (audit 2026-06-28, H5)

test_that("gintervals.liftover keeps a mapping reached only by a wide earlier overlapping chain", {
    local_db_state()

    # Target genome: chr1 long enough for the query's target, chr2 for the narrow chain.
    setup_db(list(
        paste0(">chr1\n", strrep("ACGT", 25), "\n"), # 100 bp
        paste0(">chr2\n", strrep("ACGT", 16), "\n") # 64 bp
    ))

    chain_file <- new_chain_file()
    # Wide chain A (smaller start_src): source1[0,100) -> chr1[0,100)
    write_chain_entry(chain_file, "source1", 100, "+", 0, 100, "chr1", 100, "+", 0, 100, 1)
    # Narrow chain B (larger start_src => this is back()): source1[10,20) -> chr2[50,60)
    write_chain_entry(chain_file, "source1", 100, "+", 10, 20, "chr2", 64, "+", 50, 60, 2)

    chain <- gintervals.load_chain(chain_file, src_overlap_policy = "keep", tgt_overlap_policy = "keep")

    # Query lies inside A (reaches 100) but past B (ends at 20).
    src_intervals <- data.frame(chrom = "source1", start = 50, end = 51, stringsAsFactors = FALSE)
    result <- gintervals.liftover(src_intervals, chain)

    # Before the fix the back() guard returned "no overlap" and the mapping was dropped (NULL).
    expect_false(is.null(result))
    expect_true("chr1" %in% result$chrom)
    chr1_rows <- result[result$chrom == "chr1", ]
    expect_equal(as.numeric(chr1_rows$start), 50)
    expect_equal(as.numeric(chr1_rows$end), 51)
})

test_that("ADVERSARIAL H5: relaxing the back() guard did not start mapping out-of-range queries", {
    local_db_state()
    setup_db(list(paste0(">chr1\n", strrep("ACGT", 50), "\n"))) # 200 bp

    chain_file <- new_chain_file()
    # Two non-overlapping source chains, both onto chr1.
    write_chain_entry(chain_file, "source1", 1000, "+", 100, 200, "chr1", 200, "+", 0, 100, 1)
    write_chain_entry(chain_file, "source1", 1000, "+", 300, 400, "chr1", 200, "+", 100, 200, 2)
    chain <- gintervals.load_chain(chain_file, src_overlap_policy = "keep", tgt_overlap_policy = "keep")

    lift <- function(s, e) {
        gintervals.liftover(data.frame(chrom = "source1", start = s, end = e, stringsAsFactors = FALSE), chain)
    }

    # Inside a chain -> maps; gap / before / after / wrong-chrom -> must stay unmapped.
    expect_false(is.null(lift(120, 130))) # inside chain 1
    expect_false(is.null(lift(350, 360))) # inside chain 2
    expect_null(lift(0, 50)) # before all chains
    expect_null(lift(220, 260)) # in the gap between chains
    expect_null(lift(500, 600)) # past all chains (the back() guard's job)
    expect_null(lift(950, 1000)) # far past all chains
})

test_that("ADVERSARIAL H5: minus-strand chains still lift correctly (back() guard change)", {
    local_db_state()
    setup_db(list(paste0(">chr1\n", strrep("ACGT", 50), "\n"))) # 200 bp

    chain_file <- new_chain_file()
    # source1[0,100) -> chr1[0,100) on the MINUS target strand
    write_chain_entry(chain_file, "source1", 200, "+", 0, 100, "chr1", 200, "-", 0, 100, 1)
    chain <- gintervals.load_chain(chain_file)

    res <- gintervals.liftover(
        data.frame(chrom = "source1", start = 10, end = 20, stringsAsFactors = FALSE), chain
    )
    expect_false(is.null(res))
    expect_equal(as.character(res$chrom[1]), "chr1")
    expect_equal(res$end[1] - res$start[1], 10) # width preserved under reversal
    expect_true(res$start[1] >= 0 && res$end[1] <= 200) # within the target chromosome
})

test_that("ADVERSARIAL M1: liftover value_col + canonic + aggregation returns finite output (scores-OOB path)", {
    local_db_state()
    setup_source_db(list(paste0(">source1\n", strrep("A", 400), "\n")))
    setup_db(list(paste0(">chrA\n", strrep("T", 400), "\n")))

    cf <- new_chain_file()
    # three sources mapping to OVERLAPPING target regions -> aggregation segments
    write_chain_entry(cf, "chrsource1", 400, "+", 0, 10, "chrA", 400, "+", 0, 10, 1)
    write_chain_entry(cf, "chrsource1", 400, "+", 10, 20, "chrA", 400, "+", 3, 13, 2)
    write_chain_entry(cf, "chrsource1", 400, "+", 20, 30, "chrA", 400, "+", 7, 17, 3)
    src <- data.frame(
        chrom = "chrsource1", start = c(0, 10, 20), end = c(10, 20, 30),
        value = c(1, 2, 3), stringsAsFactors = FALSE
    )

    # value_col + agg + canonic=TRUE, include_metadata=FALSE (default): the path where
    # `scores` was read out of bounds. Output must be correct, with no crash.
    res <- gintervals.liftover(src, cf, value_col = "value", multi_target_agg = "sum", canonic = TRUE)
    expect_false(is.null(res))
    expect_true("value" %in% names(res))
    expect_true(all(is.finite(res$value)))

    # include_metadata=TRUE path (unchanged by the fix) still emits a finite score column
    res2 <- gintervals.liftover(src, cf, tgt_overlap_policy = "auto_score", include_metadata = TRUE, canonic = TRUE)
    expect_true("score" %in% names(res2))
    expect_true(all(is.finite(res2$score)))
})
