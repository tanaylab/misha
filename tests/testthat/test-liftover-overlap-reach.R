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
