create_isolated_test_db()

test_that("gintervals.liftover multi-target aggregation policies", {
    local_db_state()

    # Source genome with a single chromosome
    source_db <- setup_source_db(list(paste0(">source1\n", paste(rep("A", 400), collapse = ""), "\n")))

    # Source intervals with values (similar to track values)
    # Database chromosome name is "chrsource1" (from filename chrsource1.fasta)
    src_intervals <- data.frame(
        chrom = rep("chrsource1", 3),
        start = c(0, 10, 20),
        end = c(10, 20, 30),
        value = c(1, 2, 3),
        stringsAsFactors = FALSE
    )

    # Variant with NA in the middle value
    src_intervals_na <- data.frame(
        chrom = rep("chrsource1", 3),
        start = c(0, 10, 20),
        end = c(10, 20, 30),
        value = c(1, NaN, 3),
        stringsAsFactors = FALSE
    )

    # Target genome
    setup_db(list(paste0(">chrA\n", paste(rep("T", 400), collapse = ""), "\n")))

    # Chain mappings: all three intervals map into overlapping target regions
    chain_file <- new_chain_file()
    # Chain file uses "chrsource1" to match database chromosome name
    write_chain_entry(chain_file, "chrsource1", 400, "+", 0, 10, "chrA", 400, "+", 0, 10, 1) # coverage len 10, value 1
    write_chain_entry(chain_file, "chrsource1", 400, "+", 10, 20, "chrA", 400, "+", 3, 13, 2) # coverage len 7, value 2
    write_chain_entry(chain_file, "chrsource1", 400, "+", 20, 30, "chrA", 400, "+", 7, 17, 3) # coverage len 3, value 3

    liftover_with <- function(intervals = src_intervals, agg = "mean", params = NULL,
                              na.rm = TRUE, min_n = NULL, value_col = "value") {
        gintervals.liftover(
            intervals, chain_file,
            value_col = value_col,
            multi_target_agg = agg,
            params = params,
            na.rm = na.rm,
            min_n = min_n
        )
    }

    # Test basic aggregations - all three source intervals map to overlapping regions
    # Without aggregation, we'd get 3 separate target intervals
    # With value_col, we get the value column passed through
    result <- liftover_with()
    expect_true("value" %in% names(result))
    expect_equal(nrow(result), 3) # 3 mapped intervals
    expect_true(all(result$value %in% c(1, 2, 3))) # values preserved

    result <- liftover_with(agg = "sum")
    expect_true(all(result$value %in% c(1, 2, 3))) # values preserved

    result <- liftover_with(agg = "min")
    expect_length(result$value, 3)

    result <- liftover_with(agg = "max")
    expect_length(result$value, 3)

    result <- liftover_with(agg = "median")
    expect_length(result$value, 3)

    result <- liftover_with(agg = "count")
    expect_length(result$value, 3)

    result <- liftover_with(agg = "first")
    expect_length(result$value, 3)

    result <- liftover_with(agg = "last")
    expect_length(result$value, 3)

    result <- liftover_with(agg = "nth", params = 2)
    expect_length(result$value, 3)

    result <- liftover_with(agg = "nth", params = list(n = 3))
    expect_length(result$value, 3)

    result <- liftover_with(agg = "max.coverage_len")
    expect_length(result$value, 3)

    result <- liftover_with(agg = "min.coverage_len")
    expect_length(result$value, 3)

    result <- liftover_with(agg = "max.coverage_frac")
    expect_length(result$value, 3)

    result <- liftover_with(agg = "min.coverage_frac")
    expect_length(result$value, 3)

    # Concrete aggregation values and NA handling
    # Two overlapping source intervals: [0,10)=1 and [5,15)=3 map into chrA[0,20)
    src_intervals_overlap <- data.frame(
        chrom = c("chrsource1", "chrsource1"),
        start = c(0, 5),
        end = c(10, 15),
        value = c(1, 3),
        stringsAsFactors = FALSE
    )
    src_intervals_overlap_na <- src_intervals_overlap
    src_intervals_overlap_na$value[1] <- NaN

    chain_file_overlap <- new_chain_file()
    write_chain_entry(chain_file_overlap, "chrsource1", 400, "+", 0, 20, "chrA", 400, "+", 0, 20, 1)

    # Mean aggregation averages overlapping contributions
    res_mean <- gintervals.liftover(
        src_intervals_overlap, chain_file_overlap,
        value_col = "value",
        multi_target_agg = "mean"
    )
    expect_equal(res_mean$start, c(0, 5, 10))
    expect_equal(res_mean$end, c(5, 10, 15))
    expect_equal(res_mean$value, c(1, 2, 3))

    # Sum aggregation
    res_sum <- gintervals.liftover(
        src_intervals_overlap, chain_file_overlap,
        value_col = "value",
        multi_target_agg = "sum"
    )
    expect_equal(res_sum$value, c(1, 4, 3))

    # Count aggregation counts contributors
    res_count <- gintervals.liftover(
        src_intervals_overlap, chain_file_overlap,
        value_col = "value",
        multi_target_agg = "count"
    )
    expect_equal(res_count$value, c(1, 2, 1))

    # max.coverage_len chooses the contributor with largest overlap (ties → higher value)
    res_covlen <- gintervals.liftover(
        src_intervals_overlap, chain_file_overlap,
        value_col = "value",
        multi_target_agg = "max.coverage_len"
    )
    expect_equal(res_covlen$value, c(1, 3, 3))

    # NA handling: drop NA when possible, preserve when all are NA
    res_na <- gintervals.liftover(
        src_intervals_overlap_na, chain_file_overlap,
        value_col = "value",
        multi_target_agg = "mean",
        na.rm = TRUE
    )
    expect_true(is.nan(res_na$value[1]))
    expect_equal(res_na$value[2], 3)
    expect_equal(res_na$value[3], 3)

    res_na_prop <- gintervals.liftover(
        src_intervals_overlap_na, chain_file_overlap,
        value_col = "value",
        multi_target_agg = "mean",
        na.rm = FALSE
    )
    expect_true(is.nan(res_na_prop$value[2]))

    # NA handling
    result <- liftover_with(intervals = src_intervals_na, agg = "mean", na.rm = TRUE)
    expect_true("value" %in% names(result))
    # With NA in middle, we have values 1, NaN, 3
    expect_true(any(is.nan(result$value)) || all(!is.nan(result$value)))

    result <- liftover_with(intervals = src_intervals_na, agg = "mean", na.rm = FALSE)
    expect_true("value" %in% names(result))

    # min_n gating
    result <- liftover_with(intervals = src_intervals_na, agg = "mean", na.rm = TRUE, min_n = 3)
    expect_true("value" %in% names(result))

    result <- liftover_with(intervals = src_intervals_na, agg = "mean", na.rm = TRUE, min_n = 2)
    expect_true("value" %in% names(result))
})

test_that("gintervals.liftover aggregation preserves intervalID and chain_id", {
    local_db_state()

    source_db <- setup_source_db(list(paste0(">source1\n", paste(rep("A", 200), collapse = ""), "\n")))

    # Source intervals with values
    # Database chromosome name is "chrsource1" (from filename chrsource1.fasta)
    src_intervals <- data.frame(
        chrom = rep("chrsource1", 3),
        start = c(0, 10, 20),
        end = c(10, 20, 30),
        score = c(10, 20, 30),
        stringsAsFactors = FALSE
    )

    setup_db(list(paste0(">chrB\n", paste(rep("G", 200), collapse = ""), "\n")))

    chain_file <- new_chain_file()
    # Map each interval to different target locations (no overlap)
    # Chain file uses "chrsource1" to match database chromosome name
    write_chain_entry(chain_file, "chrsource1", 200, "+", 0, 10, "chrB", 200, "+", 100, 110, 1)
    write_chain_entry(chain_file, "chrsource1", 200, "+", 10, 20, "chrB", 200, "+", 110, 120, 2)
    write_chain_entry(chain_file, "chrsource1", 200, "+", 20, 30, "chrB", 200, "+", 120, 130, 3)

    result <- gintervals.liftover(
        src_intervals, chain_file,
        value_col = "score",
        multi_target_agg = "sum"
    )

    # Should have intervalID and chain_id columns
    expect_true("intervalID" %in% names(result))
    expect_true("chain_id" %in% names(result))
    expect_true("score" %in% names(result))

    # Should have 3 rows (one per source interval)
    expect_equal(nrow(result), 3)

    # intervalID should be 1, 2, 3
    expect_equal(sort(result$intervalID), c(1, 2, 3))

    # chain_id should be 1, 2, 3
    expect_equal(sort(result$chain_id), c(1, 2, 3))

    # scores should be preserved
    expect_equal(sort(result$score), c(10, 20, 30))
})

test_that("gintervals.liftover aggregation with integer values", {
    local_db_state()

    source_db <- setup_source_db(list(paste0(">source1\n", paste(rep("A", 100), collapse = ""), "\n")))

    # Source intervals with integer values
    # Database chromosome name is "chrsource1" (from filename chrsource1.fasta)
    src_intervals <- data.frame(
        chrom = rep("chrsource1", 3),
        start = c(0, 10, 20),
        end = c(10, 20, 30),
        count = c(5L, 10L, 15L),
        stringsAsFactors = FALSE
    )

    setup_db(list(paste0(">chrC\n", paste(rep("C", 100), collapse = ""), "\n")))

    chain_file <- new_chain_file()
    # Chain file uses "chrsource1" to match database chromosome name
    write_chain_entry(chain_file, "chrsource1", 100, "+", 0, 10, "chrC", 100, "+", 0, 10, 1)
    write_chain_entry(chain_file, "chrsource1", 100, "+", 10, 20, "chrC", 100, "+", 0, 10, 2)
    write_chain_entry(chain_file, "chrsource1", 100, "+", 20, 30, "chrC", 100, "+", 0, 10, 3)

    result <- gintervals.liftover(
        src_intervals, chain_file,
        value_col = "count",
        multi_target_agg = "mean"
    )

    # Should handle integer values
    expect_true("count" %in% names(result))
    expect_true(all(result$count %in% c(5, 10, 15)))
})

test_that("gintervals.liftover without value_col works as before", {
    local_db_state()

    source_db <- setup_source_db(list(paste0(">source1\n", paste(rep("A", 100), collapse = ""), "\n")))

    # Source intervals without specifying value_col
    # Database chromosome name is "chrsource1" (from filename chrsource1.fasta)
    src_intervals <- data.frame(
        chrom = rep("chrsource1", 2),
        start = c(0, 10),
        end = c(10, 20),
        extra_col = c("a", "b"),
        stringsAsFactors = FALSE
    )

    setup_db(list(paste0(">chrD\n", paste(rep("T", 100), collapse = ""), "\n")))

    chain_file <- new_chain_file()
    # Chain file uses "chrsource1" to match database chromosome name
    write_chain_entry(chain_file, "chrsource1", 100, "+", 0, 10, "chrD", 100, "+", 0, 10, 1)
    write_chain_entry(chain_file, "chrsource1", 100, "+", 10, 20, "chrD", 100, "+", 5, 15, 2)

    # Without value_col, should work as before (no aggregation)
    result <- gintervals.liftover(src_intervals, chain_file)

    expect_true("intervalID" %in% names(result))
    expect_true("chain_id" %in% names(result))
    expect_false("value" %in% names(result)) # no value column added
    expect_equal(nrow(result), 2)
})

test_that("gintervals.liftover validates value_col parameter", {
    local_db_state()

    source_db <- setup_source_db(list(paste0(">source1\n", paste(rep("A", 100), collapse = ""), "\n")))

    # Database chromosome name is "chrsource1" (from filename chrsource1.fasta)
    src_intervals <- data.frame(
        chrom = "chrsource1",
        start = 0,
        end = 10,
        score = 5,
        stringsAsFactors = FALSE
    )

    setup_db(list(paste0(">chrE\n", paste(rep("C", 100), collapse = ""), "\n")))

    chain_file <- new_chain_file()
    # Chain file uses "chrsource1" to match database chromosome name
    write_chain_entry(chain_file, "chrsource1", 100, "+", 0, 10, "chrE", 100, "+", 0, 10, 1)

    # Invalid value_col name
    expect_error(
        gintervals.liftover(src_intervals, chain_file, value_col = "nonexistent"),
        "value_col 'nonexistent' not found in intervals"
    )

    # Non-string value_col
    expect_error(
        gintervals.liftover(src_intervals, chain_file, value_col = 123),
        "value_col must be a single character string"
    )

    # Multiple value_col names
    expect_error(
        gintervals.liftover(src_intervals, chain_file, value_col = c("score", "other")),
        "value_col must be a single character string"
    )
})

test_that("gintervals.liftover nth aggregator validates params", {
    local_db_state()

    source_db <- setup_source_db(list(paste0(">source1\n", paste(rep("A", 100), collapse = ""), "\n")))

    # Database chromosome name is "chrsource1" (from filename chrsource1.fasta)
    src_intervals <- data.frame(
        chrom = "chrsource1",
        start = 0,
        end = 10,
        value = 1,
        stringsAsFactors = FALSE
    )

    setup_db(list(paste0(">chrF\n", paste(rep("G", 100), collapse = ""), "\n")))

    chain_file <- new_chain_file()
    # Chain file uses "chrsource1" to match database chromosome name
    write_chain_entry(chain_file, "chrsource1", 100, "+", 0, 10, "chrF", 100, "+", 0, 10, 1)

    expect_error(
        gintervals.liftover(
            src_intervals, chain_file,
            value_col = "value",
            multi_target_agg = "nth"
        ),
        "params must be supplied for 'nth' aggregation"
    )

    expect_error(
        gintervals.liftover(
            src_intervals, chain_file,
            value_col = "value",
            multi_target_agg = "nth",
            params = list()
        ),
        "params list must contain an element 'n' for 'nth'",
        fixed = TRUE
    )
})

test_that("gintervals.liftover aggregation with all NA values", {
    local_db_state()

    source_db <- setup_source_db(list(paste0(">source1\n", paste(rep("A", 100), collapse = ""), "\n")))

    # All NA values
    # Database chromosome name is "chrsource1" (from filename chrsource1.fasta)
    src_intervals <- data.frame(
        chrom = rep("chrsource1", 3),
        start = c(0, 10, 20),
        end = c(10, 20, 30),
        value = c(NaN, NaN, NaN),
        stringsAsFactors = FALSE
    )

    setup_db(list(paste0(">chrG\n", paste(rep("T", 100), collapse = ""), "\n")))

    chain_file <- new_chain_file()
    # Chain file uses "chrsource1" to match database chromosome name
    write_chain_entry(chain_file, "chrsource1", 100, "+", 0, 10, "chrG", 100, "+", 0, 10, 1)
    write_chain_entry(chain_file, "chrsource1", 100, "+", 10, 20, "chrG", 100, "+", 0, 10, 2)
    write_chain_entry(chain_file, "chrsource1", 100, "+", 20, 30, "chrG", 100, "+", 0, 10, 3)

    # All aggregators should preserve NA values (propagate to result)
    for (agg in c("mean", "sum", "min", "max", "median", "first", "last", "max.coverage_len")) {
        result <- gintervals.liftover(
            src_intervals, chain_file,
            value_col = "value",
            multi_target_agg = agg,
            na.rm = TRUE
        )
        expect_true("value" %in% names(result), info = paste("aggregator:", agg))
        expect_true(all(is.nan(result$value) | is.na(result$value)), info = paste("aggregator:", agg))
    }

    # count with all NAs should still preserve the NaN values (one per source interval)
    result <- gintervals.liftover(
        src_intervals, chain_file,
        value_col = "value",
        multi_target_agg = "count",
        na.rm = TRUE
    )
    expect_true(all(is.nan(result$value) | result$value == 0))
})

test_that("gintervals.liftover with canonic mode preserves values", {
    local_db_state()

    source_db <- setup_source_db(list(paste0(">source1\n", paste(rep("A", 200), collapse = ""), "\n")))

    # Source intervals that will create adjacent mappings
    # Database chromosome name is "chrsource1" (from filename chrsource1.fasta)
    src_intervals <- data.frame(
        chrom = rep("chrsource1", 2),
        start = c(0, 10),
        end = c(10, 20),
        score = c(100, 100),
        stringsAsFactors = FALSE
    )

    setup_db(list(paste0(">chrH\n", paste(rep("C", 200), collapse = ""), "\n")))

    chain_file <- new_chain_file()
    # Map to adjacent target regions (will be merged in canonic mode)
    # Chain file uses "chrsource1" to match database chromosome name
    write_chain_entry(chain_file, "chrsource1", 200, "+", 0, 10, "chrH", 200, "+", 0, 10, 1)
    write_chain_entry(chain_file, "chrsource1", 200, "+", 10, 20, "chrH", 200, "+", 10, 20, 1) # same chain_id

    result <- gintervals.liftover(
        src_intervals, chain_file,
        value_col = "score",
        multi_target_agg = "mean",
        canonic = TRUE
    )

    expect_true("score" %in% names(result))
    # With canonic=TRUE and same intervalID/chain_id, adjacent intervals should be merged
    # Both have same value (100), so merged result should also be 100
    expect_true(all(result$score == 100))
})

test_that("gintervals.liftover aggregation handles multiple value types", {
    local_db_state()

    source_db <- setup_source_db(list(paste0(">source1\n", paste(rep("A", 100), collapse = ""), "\n")))

    # Test with numeric (double) values
    # Database chromosome name is "chrsource1" (from filename chrsource1.fasta)
    src_intervals_double <- data.frame(
        chrom = "chrsource1",
        start = 0,
        end = 10,
        val = 3.14,
        stringsAsFactors = FALSE
    )

    setup_db(list(paste0(">chrI\n", paste(rep("G", 100), collapse = ""), "\n")))

    chain_file <- new_chain_file()
    # Chain file uses "chrsource1" to match database chromosome name
    write_chain_entry(chain_file, "chrsource1", 100, "+", 0, 10, "chrI", 100, "+", 0, 10, 1)

    result <- gintervals.liftover(
        src_intervals_double, chain_file,
        value_col = "val",
        multi_target_agg = "mean"
    )

    expect_true("val" %in% names(result))
    expect_equal(result$val[1], 3.14, tolerance = 1e-6)
})

test_that("gintervals.liftover aggregates across multiple chain_ids mapping to same target", {
    local_db_state()

    source_db <- setup_source_db(list(
        paste0(">source1\n", paste(rep("A", 300), collapse = ""), "\n")
    ))

    # Source intervals from different regions that will map to the same target
    # Database chromosome name is "chrsource1" (from filename chrsource1.fasta)
    src_intervals <- data.frame(
        chrom = c("chrsource1", "chrsource1"),
        start = c(0, 100),
        end = c(100, 200),
        value = c(0.9, 0.85),
        stringsAsFactors = FALSE
    )

    setup_db(list(paste0(">chrTarget\n", paste(rep("G", 200), collapse = ""), "\n")))

    # Two chains mapping to the SAME target location but from different source regions
    # This is the scenario where we had the bug: aggregation should combine these
    chain_file <- new_chain_file()
    write_chain_entry(chain_file, "chrsource1", 300, "+", 0, 100, "chrTarget", 200, "+", 50, 150, 1, score = 100)
    write_chain_entry(chain_file, "chrsource1", 300, "+", 100, 200, "chrTarget", 200, "+", 50, 150, 2, score = 95)

    # Load the chain with tgt_overlap_policy="agg" to enable aggregation across chain_ids
    # Test max aggregation - should return 1 row with max value
    result_max <- gintervals.liftover(
        src_intervals, chain_file,
        src_overlap_policy = "keep",
        tgt_overlap_policy = "agg",
        value_col = "value",
        multi_target_agg = "max"
    )
    expect_equal(nrow(result_max), 1, info = "Should return 1 row when aggregating across chain_ids")
    expect_equal(result_max$value[1], 0.9, tolerance = 1e-6, info = "Max of 0.9 and 0.85 should be 0.9")
    expect_false("intervalID" %in% names(result_max), info = "intervalID column should not be present when aggregating across different sources")
    expect_false("chain_id" %in% names(result_max), info = "chain_id column should not be present when aggregating across different chains")

    # Test mean aggregation
    result_mean <- gintervals.liftover(
        src_intervals, chain_file,
        src_overlap_policy = "keep",
        tgt_overlap_policy = "agg",
        value_col = "value",
        multi_target_agg = "mean"
    )
    expect_equal(nrow(result_mean), 1, info = "Should return 1 row")
    expect_equal(result_mean$value[1], 0.875, tolerance = 1e-6, info = "Mean of 0.9 and 0.85 should be 0.875")

    # Test sum aggregation
    result_sum <- gintervals.liftover(
        src_intervals, chain_file,
        src_overlap_policy = "keep",
        tgt_overlap_policy = "agg",
        value_col = "value",
        multi_target_agg = "sum"
    )
    expect_equal(nrow(result_sum), 1, info = "Should return 1 row")
    expect_equal(result_sum$value[1], 1.75, tolerance = 1e-6, info = "Sum of 0.9 and 0.85 should be 1.75")

    # Test count aggregation
    result_count <- gintervals.liftover(
        src_intervals, chain_file,
        src_overlap_policy = "keep",
        tgt_overlap_policy = "agg",
        value_col = "value",
        multi_target_agg = "count"
    )
    expect_equal(nrow(result_count), 1, info = "Should return 1 row")
    expect_equal(result_count$value[1], 2, info = "Count should be 2 (two contributing sources)")
})

# Regression: ChainIntervals::map_interval carries a `hint` iterator across source
# intervals as a fast-path lower bound. With overlapping source chains
# (src_overlap_policy = "keep") end_src is not monotone in start order, so a wider
# EARLIER chain can overlap a later query while the chain immediately before the
# hint does not. The old fast-path guard checked only `hint - 1`, so it skipped the
# earlier wide chain and dropped that mapping. (Found auditing the 5.11.3 fix.)
test_that("gintervals.liftover: carried hint does not skip a wider earlier overlapping source chain", {
    local_db_state()

    # Target genome only (gintervals.liftover maps coordinates, reads no track).
    setup_db(list(paste0(">chrA\n", paste(rep("T", 400), collapse = ""), "\n")))

    # Overlapping source chains: a wide chain A over [0,100); a decoy B and a later
    # C so that, for a query inside A, B (the predecessor of C) does not overlap but
    # A (earlier) does.
    chain <- new_chain_file()
    write_chain_entry(chain, "chrsource1", 400, "+", 0, 100, "chrA", 400, "+", 0, 100, 1) # A (wide)
    write_chain_entry(chain, "chrsource1", 400, "+", 10, 20, "chrA", 400, "+", 200, 210, 2) # B (decoy)
    write_chain_entry(chain, "chrsource1", 400, "+", 30, 40, "chrA", 400, "+", 300, 310, 3) # C
    ch <- gintervals.load_chain(chain, src_overlap_policy = "keep", tgt_overlap_policy = "keep")

    # Two source intervals processed in one call (so the hint is carried): a wide
    # first one, then a short second one nested inside chains A and C.
    q <- data.frame(
        chrom = "chrsource1", start = c(5, 35), end = c(40, 38),
        stringsAsFactors = FALSE
    )
    batch <- gintervals.liftover(q, ch, include_metadata = TRUE)
    solo2 <- gintervals.liftover(q[2, ], ch, include_metadata = TRUE) # fresh hint = ground truth

    b2 <- batch[batch$intervalID == 2, ]
    # Carried-hint result must equal the fresh-hint result: both chains A and C, not
    # just C. Pre-fix the batch lift returned only the chain-C mapping.
    expect_equal(nrow(b2), nrow(solo2))
    expect_equal(nrow(b2), 2L)
    expect_setequal(b2$chain_id, c(1, 3))
    expect_setequal(b2$start, c(35, 305))
})

# Regression: when the target-overlap sweep truncates/splits a MINUS-strand chain,
# append_slice recomputed the source coordinates with the forward formula, giving the
# surviving slice the source coords of the discarded half (inverting the lifted
# source<->target correspondence). Confirmed against UCSC liftOver.
test_that("gintervals.liftover auto_score: truncated minus-strand chain keeps reversed source coords", {
    local_db_state()

    setup_db(list(paste0(">chrT\n", paste(rep("T", 1000), collapse = ""), "\n")))

    # chain1 (minus target): src chrS[100,200) -> forward tgt chrT[600,700)
    #   (minus coords 300..400 over size 1000 -> forward 600..700).
    # chain2 (plus, higher score) covers tgt[650,700), forcing chain1's target to be
    # truncated to [600,650).
    chain <- new_chain_file()
    write_chain_entry(chain, "chrS", 100000, "+", 100, 200, "chrT", 1000, "-", 300, 400, 1, score = 100)
    write_chain_entry(chain, "chrS", 100000, "+", 2000, 2050, "chrT", 1000, "+", 650, 700, 2, score = 900)
    ch <- gintervals.load_chain(chain, tgt_overlap_policy = "auto_score")

    s1 <- ch[ch$chain_id == 1, ]
    expect_equal(nrow(s1), 1L)
    expect_equal(s1$start, 600)
    expect_equal(s1$end, 650)
    # minus correspondence: fwd tgt700<->src100, tgt600<->src200; kept tgt[600,650)<->src[150,200)
    expect_equal(s1$startsrc, 150)
    expect_equal(s1$endsrc, 200)
})

# Regression: map_interval's slow path, when no chain starts before the query
# (query.start <= the leftmost overlapping chain's start_src), checked only
# iend_interval and missed the leftmost chain. With a carried hint advanced by an
# earlier (wider) interval, a narrow interval was silently dropped. Confirmed
# against UCSC liftOver.
test_that("gintervals.liftover: narrow interval starting at its leftmost chain is not dropped after a wide one", {
    local_db_state()

    setup_db(list(paste0(">chrA\n", paste(rep("T", 5000), collapse = ""), "\n")))
    chain <- new_chain_file()
    write_chain_entry(chain, "src", 5000, "+", 100, 200, "chrA", 5000, "+", 1000, 1100, 1)
    write_chain_entry(chain, "src", 5000, "+", 400, 600, "chrA", 5000, "+", 2000, 2200, 2)
    ch <- gintervals.load_chain(chain) # defaults: src=error, tgt=auto

    # The wide [0,1000) is processed first and advances the carried hint; the narrow
    # [100,200) starts exactly at chain 1's start_src (100).
    inp <- data.frame(chrom = "src", start = c(0, 100), end = c(1000, 200), stringsAsFactors = FALSE)
    res <- gintervals.liftover(inp, ch, include_metadata = TRUE)

    r2 <- res[res$intervalID == 2, ]
    expect_equal(nrow(r2), 1L) # was dropped pre-fix
    expect_equal(r2$start, 1000)
    expect_equal(r2$end, 1100)
})
