create_isolated_test_db()

test_that("interval_relative creates bins from interval start", {
    # Create test intervals: [100, 300] and [500, 700] on chr1
    intervs <- gintervals(1, c(100, 500), c(300, 700))

    result <- giterator.intervals(NULL, intervs, iterator = 50, interval_relative = TRUE)

    # Check that intervalID column exists
    expect_true("intervalID" %in% names(result))
    expect_equal(ncol(result), 4) # chrom, start, end, intervalID

    # First interval [100, 300] should produce 4 bins of size 50: 100-150, 150-200, 200-250, 250-300
    first_interval_bins <- result[result$intervalID == 1, ]
    expect_equal(nrow(first_interval_bins), 4)
    expect_equal(first_interval_bins$start, c(100, 150, 200, 250))
    expect_equal(first_interval_bins$end, c(150, 200, 250, 300))

    # Second interval [500, 700] should produce 4 bins of size 50: 500-550, 550-600, 600-650, 650-700
    second_interval_bins <- result[result$intervalID == 2, ]
    expect_equal(nrow(second_interval_bins), 4)
    expect_equal(second_interval_bins$start, c(500, 550, 600, 650))
    expect_equal(second_interval_bins$end, c(550, 600, 650, 700))
})

test_that("interval_relative handles partial bins with clip mode", {
    # Interval of size 125 with binsize 50 -> 2 full + 1 partial
    intervs <- gintervals(1, 0, 125)

    # clip mode (default): should get 3 bins (last one is 25bp)
    result <- giterator.intervals(NULL, intervs, iterator = 50, interval_relative = TRUE, partial_bins = "clip")

    expect_equal(nrow(result), 3)
    expect_equal(result$start, c(0, 50, 100))
    expect_equal(result$end, c(50, 100, 125)) # Last bin clipped to 125
})

test_that("interval_relative handles partial bins with exact mode", {
    # Interval of size 125 with binsize 50 -> only 2 full bins
    intervs <- gintervals(1, 0, 125)

    # exact mode: should get 2 bins (partial dropped)
    result <- giterator.intervals(NULL, intervs, iterator = 50, interval_relative = TRUE, partial_bins = "exact")

    expect_equal(nrow(result), 2)
    expect_equal(result$start, c(0, 50))
    expect_equal(result$end, c(50, 100))
})

test_that("interval_relative handles partial bins with drop mode (alias for exact)", {
    intervs <- gintervals(1, 0, 125)

    result <- giterator.intervals(NULL, intervs, iterator = 50, interval_relative = TRUE, partial_bins = "drop")

    expect_equal(nrow(result), 2) # Same as exact
})

test_that("interval_relative intervalID tracking is correct", {
    # Multiple intervals across different chromosomes
    intervs <- gintervals(c(1, 2), c(0, 100), c(100, 300))

    result <- giterator.intervals(NULL, intervs, iterator = 50, interval_relative = TRUE)

    # First interval on chr1 [0, 100] produces 2 bins
    expect_equal(sum(result$intervalID == 1), 2)

    # Second interval on chr2 [100, 300] produces 4 bins
    expect_equal(sum(result$intervalID == 2), 4)
})

test_that("interval_relative requires numeric iterator", {
    intervs <- gintervals(1, 0, 100)

    # Should error when iterator is not numeric
    expect_error(giterator.intervals(NULL, intervs, iterator = "test.sparse", interval_relative = TRUE))
})

test_that("interval_relative works with single interval", {
    intervs <- gintervals(1, 1000, 2000)

    result <- giterator.intervals(NULL, intervs, iterator = 200, interval_relative = TRUE)

    expect_equal(nrow(result), 5)
    expect_equal(result$start[1], 1000)
    expect_equal(result$end[5], 2000)
    expect_true(all(result$intervalID == 1))
})

test_that("interval_relative handles intervals smaller than binsize in clip mode", {
    # Interval smaller than binsize
    intervs <- gintervals(1, 0, 30)

    result <- giterator.intervals(NULL, intervs, iterator = 50, interval_relative = TRUE, partial_bins = "clip")

    expect_equal(nrow(result), 1)
    expect_equal(result$start, 0)
    expect_equal(result$end, 30) # Single clipped bin
})

test_that("interval_relative handles intervals smaller than binsize in exact mode", {
    # Interval smaller than binsize - should produce no bins
    intervs <- gintervals(1, 0, 30)

    result <- giterator.intervals(NULL, intervs, iterator = 50, interval_relative = TRUE, partial_bins = "exact")

    # Result should be NULL or empty for intervals smaller than binsize in exact mode
    expect_true(is.null(result) || nrow(result) == 0)
})

test_that("interval_relative default behavior unchanged without flag", {
    intervs <- gintervals(1, 100, 300)

    # Without interval_relative, should use chromosome-aligned bins
    result_default <- giterator.intervals(NULL, intervs, iterator = 50)

    # Should NOT have intervalID column
    expect_false("intervalID" %in% names(result_default))
})

# Additional comprehensive tests

test_that("interval_relative preserves chromosome assignment correctly", {
    # Intervals on different chromosomes
    intervs <- gintervals(c(1, 2, 3), c(0, 100, 200), c(100, 200, 300))

    result <- giterator.intervals(NULL, intervs, iterator = 50, interval_relative = TRUE)

    # Check chromosomes are preserved
    chr1_bins <- result[result$intervalID == 1, ]
    chr2_bins <- result[result$intervalID == 2, ]
    chr3_bins <- result[result$intervalID == 3, ]

    expect_true(all(chr1_bins$chrom == "chr1"))
    expect_true(all(chr2_bins$chrom == "chr2"))
    expect_true(all(chr3_bins$chrom == "chr3"))
})

test_that("interval_relative handles many intervals efficiently", {
    # Create 100 small intervals
    n <- 100
    intervs <- gintervals(1, seq(0, by = 200, length.out = n), seq(100, by = 200, length.out = n))

    result <- giterator.intervals(NULL, intervs, iterator = 25, interval_relative = TRUE)

    # Each 100bp interval with 25bp bins = 4 bins per interval
    expect_equal(nrow(result), n * 4)
    expect_equal(max(result$intervalID), n)
    expect_equal(min(result$intervalID), 1)
})

test_that("interval_relative differs from chromosome-aligned bins", {
    # Interval that doesn't start at a multiple of binsize
    intervs <- gintervals(1, 123, 243)

    # Interval-relative: bins start at 123
    result_relative <- giterator.intervals(NULL, intervs, iterator = 20, interval_relative = TRUE)

    # Chromosome-aligned: bins start at nearest multiple of 20 (120)
    result_aligned <- giterator.intervals(NULL, intervs, iterator = 20)

    # First bin starts should differ
    expect_equal(result_relative$start[1], 123) # Starts at interval start
    expect_equal(result_aligned$start[1], 123) # Clipped to interval, but...
    expect_equal(result_aligned$end[1], 140) # ...ends at chromosome-aligned boundary

    # In relative mode, bins are evenly spaced from interval start
    expect_equal(result_relative$start, c(123, 143, 163, 183, 203, 223))
    expect_equal(result_relative$end, c(143, 163, 183, 203, 223, 243))
})

test_that("interval_relative handles exact division (no partial bins)", {
    # Interval size exactly divisible by binsize
    intervs <- gintervals(1, 0, 200)

    result_clip <- giterator.intervals(NULL, intervs, iterator = 50, interval_relative = TRUE, partial_bins = "clip")
    result_exact <- giterator.intervals(NULL, intervs, iterator = 50, interval_relative = TRUE, partial_bins = "exact")

    # Both should give same result when there are no partial bins
    expect_equal(nrow(result_clip), 4)
    expect_equal(nrow(result_exact), 4)
    expect_equal(result_clip$start, result_exact$start)
    expect_equal(result_clip$end, result_exact$end)
})

test_that("interval_relative handles binsize of 1", {
    intervs <- gintervals(1, 100, 105)

    result <- giterator.intervals(NULL, intervs, iterator = 1, interval_relative = TRUE)

    expect_equal(nrow(result), 5)
    expect_equal(result$start, 100:104)
    expect_equal(result$end, 101:105)
})

test_that("interval_relative handles large binsize relative to interval", {
    # Multiple intervals, some smaller than binsize
    intervs <- gintervals(1, c(0, 100, 200), c(50, 180, 500))

    result <- giterator.intervals(NULL, intervs, iterator = 100, interval_relative = TRUE, partial_bins = "clip")

    # First interval [0, 50]: 1 partial bin (50bp)
    # Second interval [100, 180]: 1 partial bin (80bp)
    # Third interval [200, 500]: 3 bins (100, 100, 100bp)
    expect_equal(sum(result$intervalID == 1), 1)
    expect_equal(sum(result$intervalID == 2), 1)
    expect_equal(sum(result$intervalID == 3), 3)

    # Check the partial bins have correct sizes
    expect_equal(result$end[result$intervalID == 1] - result$start[result$intervalID == 1], 50)
    expect_equal(result$end[result$intervalID == 2] - result$start[result$intervalID == 2], 80)
})

test_that("interval_relative with exact mode skips small intervals", {
    # Mix of intervals: some can produce full bins, some cannot
    intervs <- gintervals(1, c(0, 100, 200), c(30, 250, 500))

    result <- giterator.intervals(NULL, intervs, iterator = 50, interval_relative = TRUE, partial_bins = "exact")

    # First interval [0, 30]: too small, no bins
    # Second interval [100, 250]: 3 full bins (150bp / 50 = 3)
    # Third interval [200, 500]: 6 full bins (300bp / 50 = 6)
    expect_equal(sum(result$intervalID == 1), 0)
    expect_equal(sum(result$intervalID == 2), 3)
    expect_equal(sum(result$intervalID == 3), 6)
})

test_that("interval_relative works with floating point iterator", {
    intervs <- gintervals(1, 0, 100)

    # Should work with float that converts to integer
    result <- giterator.intervals(NULL, intervs, iterator = 25.0, interval_relative = TRUE)

    expect_equal(nrow(result), 4)
})

test_that("interval_relative handles non-overlapping intervals on same chromosome", {
    # Multiple disjoint intervals on same chromosome
    intervs <- gintervals(1, c(0, 1000, 5000), c(100, 1200, 5500))

    result <- giterator.intervals(NULL, intervs, iterator = 50, interval_relative = TRUE)

    # Each interval is processed independently
    expect_equal(length(unique(result$intervalID)), 3)

    # Check that bins don't cross interval boundaries
    for (id in 1:3) {
        bins <- result[result$intervalID == id, ]
        orig <- intervs[id, ]
        expect_true(all(bins$start >= orig$start))
        expect_true(all(bins$end <= orig$end))
    }
})

test_that("interval_relative assigns intervalID based on sorted order", {
    # Input intervals in non-sorted order
    intervs <- gintervals(c(3, 1, 2), c(0, 100, 50), c(100, 200, 150))

    result <- giterator.intervals(NULL, intervs, iterator = 50, interval_relative = TRUE)

    # intervalID is assigned based on sorted order (chr1, chr2, chr3)
    chr1_bins <- result[result$chrom == "chr1", ]
    chr2_bins <- result[result$chrom == "chr2", ]
    chr3_bins <- result[result$chrom == "chr3", ]

    # chr1 comes first in sorted order -> intervalID = 1
    expect_true(all(chr1_bins$intervalID == 1))
    # chr2 comes second in sorted order -> intervalID = 2
    expect_true(all(chr2_bins$intervalID == 2))
    # chr3 comes third in sorted order -> intervalID = 3
    expect_true(all(chr3_bins$intervalID == 3))
})

test_that("interval_relative bin sizes are consistent within each interval", {
    intervs <- gintervals(1, c(0, 200), c(175, 400))

    result <- giterator.intervals(NULL, intervs, iterator = 50, interval_relative = TRUE, partial_bins = "clip")

    # Check bin sizes for first interval (all except last should be 50)
    first_bins <- result[result$intervalID == 1, ]
    bin_sizes <- first_bins$end - first_bins$start
    expect_true(all(bin_sizes[-length(bin_sizes)] == 50))
    expect_equal(bin_sizes[length(bin_sizes)], 25) # Last bin is partial

    # Check bin sizes for second interval (all should be 50, exact fit)
    second_bins <- result[result$intervalID == 2, ]
    bin_sizes2 <- second_bins$end - second_bins$start
    expect_true(all(bin_sizes2 == 50))
})

test_that("interval_relative handles intervals at chromosome start", {
    intervs <- gintervals(1, 0, 100)

    result <- giterator.intervals(NULL, intervs, iterator = 30, interval_relative = TRUE)

    # Should produce: [0,30], [30,60], [60,90], [90,100]
    expect_equal(nrow(result), 4)
    expect_equal(result$start[1], 0)
    expect_equal(result$end[4], 100)
})
