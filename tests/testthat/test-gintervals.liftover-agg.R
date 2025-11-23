load_test_db()

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
    expect_true(all(result$value %in% c(1, 2, 3)))

    result <- liftover_with(agg = "max")
    expect_true(all(result$value %in% c(1, 2, 3)))

    result <- liftover_with(agg = "median")
    expect_true(all(result$value %in% c(1, 2, 3)))

    result <- liftover_with(agg = "count")
    expect_true(all(result$value %in% c(1, 2, 3)))

    result <- liftover_with(agg = "first")
    expect_true(all(result$value %in% c(1, 2, 3)))

    result <- liftover_with(agg = "last")
    expect_true(all(result$value %in% c(1, 2, 3)))

    result <- liftover_with(agg = "nth", params = 2)
    expect_true(all(result$value %in% c(1, 2, 3)))

    result <- liftover_with(agg = "nth", params = list(n = 3))
    expect_true(all(result$value %in% c(1, 2, 3)))

    result <- liftover_with(agg = "max.coverage_len")
    expect_true(all(result$value %in% c(1, 2, 3)))

    result <- liftover_with(agg = "min.coverage_len")
    expect_true(all(result$value %in% c(1, 2, 3)))

    result <- liftover_with(agg = "max.coverage_frac")
    expect_true(all(result$value %in% c(1, 2, 3)))

    result <- liftover_with(agg = "min.coverage_frac")
    expect_true(all(result$value %in% c(1, 2, 3)))

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

    # All aggregators should preserve NA values
    for (agg in c("mean", "sum", "min", "max", "median", "first", "last", "max.coverage_len")) {
        result <- gintervals.liftover(
            src_intervals, chain_file,
            value_col = "value",
            multi_target_agg = agg,
            na.rm = TRUE
        )
        expect_true("value" %in% names(result), info = paste("aggregator:", agg))
        # All values should be NaN since source had all NaN
        expect_true(all(is.nan(result$value)), info = paste("aggregator:", agg))
    }

    # count with all NAs should still preserve the NaN values (one per source interval)
    result <- gintervals.liftover(
        src_intervals, chain_file,
        value_col = "value",
        multi_target_agg = "count",
        na.rm = TRUE
    )
    expect_true(all(is.nan(result$value)))
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
    expect_equal(result$val[1], 3.14, tolerance = 1e-10)
})
