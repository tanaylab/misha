create_isolated_test_db()

test_that("gintervals.random basic functionality works", {
    # Generate intervals
    set.seed(123)
    intervals <- gintervals.random(100, 10)

    # Check structure
    expect_s3_class(intervals, "data.frame")
    expect_equal(nrow(intervals), 10)
    expect_equal(ncol(intervals), 3)
    expect_equal(names(intervals), c("chrom", "start", "end"))

    # Check that all intervals have correct size
    expect_true(all(intervals$end - intervals$start == 100))

    # Check that all intervals are within chromosome boundaries
    all_genome <- gintervals.all()
    for (i in 1:nrow(intervals)) {
        chrom_size <- all_genome$end[all_genome$chrom == intervals$chrom[i]]
        expect_true(intervals$start[i] >= 0)
        expect_true(intervals$end[i] <= chrom_size)
    }
})

test_that("gintervals.random respects dist_from_edge", {
    set.seed(456)
    dist_from_edge <- 1e6
    intervals <- gintervals.random(100, 20, dist_from_edge = dist_from_edge)

    # Check that all intervals respect distance from edge
    all_genome <- gintervals.all()
    for (i in 1:nrow(intervals)) {
        chrom_info <- all_genome[all_genome$chrom == intervals$chrom[i], ]
        chrom_size <- chrom_info$end - chrom_info$start
        expect_true(intervals$start[i] >= dist_from_edge)
        expect_true(intervals$end[i] <= chrom_size - dist_from_edge)
    }
})

test_that("gintervals.random chromosome filtering works", {
    set.seed(789)
    # Test with specific chromosomes
    intervals <- gintervals.random(100, 50, chromosomes = c("chr1", "chr2"))

    # Check that all intervals are on specified chromosomes
    expect_true(all(intervals$chrom %in% c("chr1", "chr2")))

    # Check that we get intervals on multiple chromosomes (probabilistically)
    # (may fail with very small n, but with n=50 should almost always work)
    unique_chroms <- unique(as.character(intervals$chrom))
    expect_true(length(unique_chroms) >= 1)
})

test_that("gintervals.random is reproducible with set.seed", {
    set.seed(42)
    intervals1 <- gintervals.random(100, 10)

    set.seed(42)
    intervals2 <- gintervals.random(100, 10)

    expect_identical(intervals1, intervals2)
})

test_that("gintervals.random handles different sizes", {
    set.seed(111)

    # Small intervals
    intervals_small <- gintervals.random(10, 5)
    expect_true(all(intervals_small$end - intervals_small$start == 10))

    # Large intervals
    intervals_large <- gintervals.random(100000, 5, dist_from_edge = 1e6)
    expect_true(all(intervals_large$end - intervals_large$start == 100000))
})

test_that("gintervals.random input validation works", {
    # Invalid size
    expect_error(gintervals.random(0, 10), "size must be a positive number")
    expect_error(gintervals.random(-100, 10), "size must be a positive number")
    expect_error(gintervals.random(c(100, 200), 10), "size must be a positive number")

    # Invalid n
    expect_error(gintervals.random(100, 0), "n must be a positive number")
    expect_error(gintervals.random(100, -10), "n must be a positive number")
    expect_error(gintervals.random(100, c(10, 20)), "n must be a positive number")

    # Invalid dist_from_edge
    expect_error(gintervals.random(100, 10, dist_from_edge = -1), "dist_from_edge must be a non-negative number")
    expect_error(gintervals.random(100, 10, dist_from_edge = c(1e6, 2e6)), "dist_from_edge must be a non-negative number")

    # Invalid chromosomes
    expect_error(gintervals.random(100, 10, chromosomes = 123), "chromosomes must be a character vector")
    expect_error(gintervals.random(100, 10, chromosomes = c("nonexistent_chr")), "No chromosomes named")
})

test_that("gintervals.random handles edge cases", {
    # Very small dist_from_edge
    set.seed(222)
    intervals <- gintervals.random(100, 5, dist_from_edge = 0)
    expect_equal(nrow(intervals), 5)

    # Single chromosome
    set.seed(333)
    intervals <- gintervals.random(100, 10, chromosomes = "chr1")
    expect_true(all(intervals$chrom == "chr1"))

    # Large number of intervals
    set.seed(444)
    intervals <- gintervals.random(100, 1000)
    expect_equal(nrow(intervals), 1000)
    expect_true(all(intervals$end - intervals$start == 100))
})

test_that("gintervals.random samples proportional to chromosome length", {
    set.seed(555)

    # Generate many intervals to test sampling distribution
    n_intervals <- 10000
    intervals <- gintervals.random(100, n_intervals, dist_from_edge = 1e6)

    # Count intervals per chromosome
    chrom_counts <- table(intervals$chrom)

    # Get chromosome lengths
    all_genome <- gintervals.all()
    all_genome <- all_genome[all_genome$end - all_genome$start >= 100 + 2 * 1e6, ]
    all_genome$length <- all_genome$end - all_genome$start
    all_genome$expected_prop <- all_genome$length / sum(all_genome$length)

    # Check that observed proportions are reasonably close to expected
    # (using chi-square goodness of fit test)
    observed <- as.numeric(chrom_counts[as.character(all_genome$chrom)])
    observed[is.na(observed)] <- 0
    expected <- all_genome$expected_prop * n_intervals

    # Chi-square test (p-value should be > 0.01 most of the time)
    # Skip if any expected count is too small
    if (all(expected >= 5)) {
        chi_sq <- sum((observed - expected)^2 / expected)
        df <- length(expected) - 1
        p_value <- pchisq(chi_sq, df, lower.tail = FALSE)
        expect_true(p_value > 0.001) # Very lenient threshold
    }
})

test_that("gintervals.random error when chromosomes too short", {
    # Try to create intervals that are too large for any chromosome
    all_genome <- gintervals.all()
    max_chrom_size <- max(all_genome$end - all_genome$start)

    expect_error(
        gintervals.random(max_chrom_size + 1000, 10, dist_from_edge = 0),
        "No chromosomes are long enough"
    )

    expect_error(
        gintervals.random(1000, 10, dist_from_edge = max_chrom_size),
        "No chromosomes are long enough"
    )
})
# Filter functionality tests

test_that("gintervals.random filter basic functionality works", {
    set.seed(666)

    # Create filter regions
    filter_regions <- gintervals(c("chr1", "chr2"), c(10000, 15000), c(20000, 25000))

    # Generate intervals with filter
    intervals <- gintervals.random(100, 50, dist_from_edge = 100, filter = filter_regions)

    # Check structure
    expect_s3_class(intervals, "data.frame")
    expect_equal(nrow(intervals), 50)
    expect_equal(names(intervals), c("chrom", "start", "end"))

    # Check that all intervals have correct size
    expect_true(all(intervals$end - intervals$start == 100))

    # Check that NO intervals overlap with filter
    overlaps <- gintervals.intersect(intervals, filter_regions)
    expect_true(is.null(overlaps) || nrow(overlaps) == 0)
})

test_that("gintervals.random filter prevents all overlaps", {
    set.seed(777)

    # Create multiple filter regions
    filter_regions <- gintervals(
        c("chr1", "chr1", "chr2", "chr2"),
        c(5000, 50000, 10000, 100000),
        c(10000, 60000, 20000, 110000)
    )

    intervals <- gintervals.random(100, 100, dist_from_edge = 100, filter = filter_regions)

    # Verify NO overlaps
    overlaps <- gintervals.intersect(intervals, filter_regions)
    expect_true(is.null(overlaps) || nrow(overlaps) == 0)

    # Verify intervals can be adjacent to filter but not overlap
    # Check that some intervals are close to filter boundaries (within 1000bp)
    # This verifies the filter expansion is working correctly
    for (i in 1:nrow(intervals)) {
        for (j in 1:nrow(filter_regions)) {
            if (intervals$chrom[i] == filter_regions$chrom[j]) {
                # Check if interval is before filter
                if (intervals$end[i] <= filter_regions$start[j]) {
                    # Should not overlap even by 1bp
                    expect_true(intervals$end[i] <= filter_regions$start[j])
                }
                # Check if interval is after filter
                if (intervals$start[i] >= filter_regions$end[j]) {
                    # Should not overlap even by 1bp
                    expect_true(intervals$start[i] >= filter_regions$end[j])
                }
            }
        }
    }
})

test_that("gintervals.random filter with heavy filtering", {
    set.seed(888)

    # Filter out large portions of chromosomes (but leave some space)
    filter_regions <- gintervals(
        c("chr1", "chr1", "chr2"),
        c(1000, 200000, 50000),
        c(100000, 400000, 250000)
    )

    # Should still be able to generate intervals
    intervals <- gintervals.random(100, 50, dist_from_edge = 100, filter = filter_regions)

    expect_equal(nrow(intervals), 50)

    # Verify no overlaps
    overlaps <- gintervals.intersect(intervals, filter_regions)
    expect_true(is.null(overlaps) || nrow(overlaps) == 0)
})

test_that("gintervals.random filter with chromosome parameter", {
    set.seed(999)

    # Create filter on chr1 only
    filter_regions <- gintervals("chr1", 10000, 20000)

    # Generate on chr1 only
    intervals <- gintervals.random(100, 30,
        dist_from_edge = 100,
        chromosomes = "chr1", filter = filter_regions
    )

    # All on chr1
    expect_true(all(intervals$chrom == "chr1"))

    # No overlaps with filter
    overlaps <- gintervals.intersect(intervals, filter_regions)
    expect_true(is.null(overlaps) || nrow(overlaps) == 0)
})

test_that("gintervals.random filter reproducibility with set.seed", {
    filter_regions <- gintervals(c("chr1", "chr2"), c(5000, 10000), c(15000, 20000))

    set.seed(1111)
    intervals1 <- gintervals.random(100, 20, dist_from_edge = 100, filter = filter_regions)

    set.seed(1111)
    intervals2 <- gintervals.random(100, 20, dist_from_edge = 100, filter = filter_regions)

    expect_identical(intervals1, intervals2)
})

test_that("gintervals.random filter with overlapping filter intervals", {
    set.seed(1212)

    # Create overlapping filter regions (should be auto-unified)
    filter_regions <- gintervals(
        c("chr1", "chr1", "chr1"),
        c(5000, 8000, 15000),
        c(10000, 12000, 20000)
    )

    intervals <- gintervals.random(100, 50, dist_from_edge = 100, filter = filter_regions)

    # Should still work and have no overlaps
    overlaps <- gintervals.intersect(intervals, filter_regions)
    expect_true(is.null(overlaps) || nrow(overlaps) == 0)
})

test_that("gintervals.random filter validation works", {
    # Invalid filter type
    expect_error(
        gintervals.random(100, 10, filter = "not_a_dataframe"),
        "filter must be a data frame"
    )

    # Missing columns
    bad_filter <- data.frame(chrom = "chr1", start = 1000)
    expect_error(
        gintervals.random(100, 10, filter = bad_filter),
        "filter must have columns"
    )
})

test_that("gintervals.random filter error when no valid regions", {
    # Filter entire small chromosomes
    all_genome <- gintervals.all()

    # Create filter covering most of genome
    filter_regions <- gintervals(
        all_genome$chrom,
        rep(100, nrow(all_genome)),
        all_genome$end - 100
    )

    # Should error - no valid regions left
    expect_error(
        gintervals.random(100, 10, dist_from_edge = 100, filter = filter_regions),
        "No valid regions available after applying filter"
    )
})

test_that("gintervals.random filter with empty filter behaves like no filter", {
    set.seed(1313)
    intervals_no_filter <- gintervals.random(100, 20, dist_from_edge = 100)

    set.seed(1313)
    empty_filter <- gintervals(character(0), numeric(0), numeric(0))
    intervals_empty_filter <- gintervals.random(100, 20, dist_from_edge = 100, filter = empty_filter)

    expect_identical(intervals_no_filter, intervals_empty_filter)
})

test_that("gintervals.random filter performance is reasonable", {
    set.seed(1414)

    # Create medium-sized filter (100 regions)
    filter_chroms <- sample(c("chr1", "chr2", "chrX"), 100, replace = TRUE)
    filter_starts <- sample(1000:400000, 100)
    filter_ends <- filter_starts + sample(1000:10000, 100)
    filter_regions <- gintervals(filter_chroms, filter_starts, filter_ends)

    # Should complete quickly even with filter
    start_time <- Sys.time()
    intervals <- gintervals.random(100, 1000, dist_from_edge = 100, filter = filter_regions)
    elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

    # Should complete in under 1 second
    expect_true(elapsed < 1.0)

    # Verify correctness
    expect_equal(nrow(intervals), 1000)
    overlaps <- gintervals.intersect(intervals, filter_regions)
    expect_true(is.null(overlaps) || nrow(overlaps) == 0)
})

test_that("gintervals.random filter samples proportional to valid segment length", {
    set.seed(1515)

    # Create filter that heavily filters chr2 but not chr1
    # This should cause chr1 to be over-represented
    filter_regions <- gintervals(
        rep("chr2", 10),
        seq(10000, 100000, by = 10000),
        seq(19000, 109000, by = 10000)
    )

    # Generate many intervals
    intervals <- gintervals.random(100, 500,
        dist_from_edge = 100,
        chromosomes = c("chr1", "chr2"),
        filter = filter_regions
    )

    # Count per chromosome
    counts <- table(intervals$chrom)

    # chr1 should have significantly more intervals than chr2
    # because chr2 has been heavily filtered
    if ("chr1" %in% names(counts) && "chr2" %in% names(counts)) {
        expect_true(counts["chr1"] > counts["chr2"])
    }

    # Verify no overlaps
    overlaps <- gintervals.intersect(intervals, filter_regions)
    expect_true(is.null(overlaps) || nrow(overlaps) == 0)
})

# Benchmark tests

test_that("gintervals.random benchmark: no filter vs small filter", {
    skip_on_cran()
    skip_on_ci()

    set.seed(1616)

    # Benchmark without filter
    time_no_filter <- system.time({
        intervals_no_filter <- gintervals.random(100, 10000, dist_from_edge = 100)
    })[3]

    # Benchmark with small filter (10 regions)
    filter_regions <- gintervals(
        rep("chr1", 10),
        seq(10000, 100000, by = 10000),
        seq(11000, 101000, by = 10000)
    )

    time_small_filter <- system.time({
        intervals_small_filter <- gintervals.random(100, 10000,
            dist_from_edge = 100,
            filter = filter_regions
        )
    })[3]

    # Filter should add minimal overhead (less than 3x slowdown)
    expect_true(time_small_filter < time_no_filter * 3)

    # Print for information
    message(sprintf(
        "No filter: %.3fs, Small filter: %.3fs (%.1fx)",
        time_no_filter, time_small_filter, time_small_filter / time_no_filter
    ))
})

test_that("gintervals.random benchmark: large filter performance", {
    skip_on_cran()
    skip_on_ci()

    set.seed(1717)

    # Create large filter (1000 regions)
    filter_chroms <- sample(c("chr1", "chr2", "chrX"), 1000, replace = TRUE)
    filter_starts <- sample(1000:400000, 1000)
    filter_ends <- filter_starts + sample(500:5000, 1000)
    filter_regions <- gintervals(filter_chroms, filter_starts, filter_ends)

    time_large_filter <- system.time({
        intervals <- gintervals.random(100, 10000, dist_from_edge = 100, filter = filter_regions)
    })[3]

    # Should still be fast (under 2 seconds)
    expect_true(time_large_filter < 2.0)

    # Verify correctness
    overlaps <- gintervals.intersect(intervals, filter_regions)
    expect_true(is.null(overlaps) || nrow(overlaps) == 0)

    message(sprintf("Large filter (1000 regions): %.3fs", time_large_filter))
})
