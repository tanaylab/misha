test_that("value-based vtrack basic functionality works", {
    # Create simple intervals with values
    intervals_df <- data.frame(
        chrom = "chr1",
        start = c(100, 300, 500),
        end = c(200, 400, 600),
        score = c(10, 20, 30)
    )

    # Create value-based vtrack
    gvtrack.create("test.value.vt", src = intervals_df, func = "avg")

    # Extract over exact interval - should get exact value
    iter_int <- gintervals("chr1", 100, 200)
    result <- gextract("test.value.vt", intervals = iter_int, iterator = iter_int, colnames = "value")
    expect_equal(result$value, 10)

    # Extract over interval with no coverage - should get NA
    iter_int <- gintervals("chr1", 0, 50)
    result <- gextract("test.value.vt", intervals = iter_int, iterator = iter_int, colnames = "value")
    expect_true(is.na(result$value))

    # Extract over multiple intervals - count-based average (matching sparse tracks)
    # Interval 150-350 covers [150-200] of first (val=10) and [300-350] of second (val=20)
    # Average = (10 + 20) / 2 = 15
    iter_int <- gintervals("chr1", 150, 350)
    result <- gextract("test.value.vt", intervals = iter_int, iterator = iter_int, colnames = "value")
    expect_equal(result$value, 15)

    # Clean up
    gvtrack.rm("test.value.vt")
})

test_that("value-based vtrack min/max functions work", {
    intervals_df <- data.frame(
        chrom = "chr1",
        start = c(100, 300, 500),
        end = c(200, 400, 600),
        score = c(10, 20, 30)
    )

    # Test MIN function
    gvtrack.create("test.min.vt", src = intervals_df, func = "min")
    iter_int <- gintervals("chr1", 150, 350)
    result <- gextract("test.min.vt", intervals = iter_int, iterator = iter_int, colnames = "value")
    expect_equal(result$value, 10) # Min of 10 and 20
    gvtrack.rm("test.min.vt")

    # Test MAX function
    gvtrack.create("test.max.vt", src = intervals_df, func = "max")
    iter_int <- gintervals("chr1", 150, 350)
    result <- gextract("test.max.vt", intervals = iter_int, iterator = iter_int, colnames = "value")
    expect_equal(result$value, 20) # Max of 10 and 20
    gvtrack.rm("test.max.vt")
})

test_that("value-based vtrack sum function works", {
    intervals_df <- data.frame(
        chrom = "chr1",
        start = c(100, 300),
        end = c(200, 400),
        score = c(10, 20)
    )

    gvtrack.create("test.sum.vt", src = intervals_df, func = "sum")

    # Sum over interval [150-350]: count-based (like sparse tracks) = 10 + 20 = 30
    iter_int <- gintervals("chr1", 150, 350)
    result <- gextract("test.sum.vt", intervals = iter_int, iterator = iter_int, colnames = "value")
    expect_equal(result$value, 30)

    gvtrack.rm("test.sum.vt")
})

test_that("value-based vtrack stddev function works", {
    intervals_df <- data.frame(
        chrom = "chr1",
        start = c(100, 200),
        end = c(150, 250),
        score = c(10, 20)
    )

    gvtrack.create("test.stddev.vt", src = intervals_df, func = "stddev")

    # Stddev over both intervals
    iter_int <- gintervals("chr1", 100, 250)
    result <- gextract("test.stddev.vt", intervals = iter_int, iterator = iter_int, colnames = "value")
    # Count-based unbiased stddev (like sparse tracks): mean = (10 + 20)/2 = 15
    # Using formula: sqrt(sum(x^2)/(N-1) - N*(mean^2)/(N-1))
    # = sqrt((100 + 400)/1 - 2*(225)/1) = sqrt(500 - 450) = sqrt(50) ≈ 7.07
    expect_equal(result$value, sqrt(50), tolerance = 1e-6)

    gvtrack.rm("test.stddev.vt")
})

test_that("value-based vtrack quantile function works", {
    intervals_df <- data.frame(
        chrom = "chr1",
        start = c(100, 200, 300),
        end = c(150, 250, 350),
        score = c(10, 20, 30)
    )

    # Test median (0.5 percentile)
    gvtrack.create("test.quantile.vt", src = intervals_df, func = "quantile", params = 0.5)

    iter_int <- gintervals("chr1", 100, 350)
    result <- gextract("test.quantile.vt", intervals = iter_int, iterator = iter_int, colnames = "value")
    # Median of 10, 20, 30 (with equal coverage) should be 20
    expect_equal(result$value, 20)

    gvtrack.rm("test.quantile.vt")
})

test_that("value-based vtrack detects overlapping intervals", {
    # Create overlapping intervals
    intervals_df <- data.frame(
        chrom = "chr1",
        start = c(100, 150), # These overlap
        end = c(200, 250),
        score = c(10, 20)
    )

    # Should error because intervals overlap
    expect_error(
        gvtrack.create("test.overlap.vt", src = intervals_df, func = "avg"),
        regexp = "overlapping intervals"
    )
})

test_that("value-based vtrack works with iterator", {
    intervals_df <- data.frame(
        chrom = "chr1",
        start = c(0, 100, 200, 300, 400),
        end = c(50, 150, 250, 350, 450),
        score = c(1, 2, 3, 4, 5)
    )

    gvtrack.create("test.iter.vt", src = intervals_df, func = "avg")

    # Extract with 100bp iterator
    result <- gextract("test.iter.vt", gintervals("chr1", 0, 500), iterator = 100, colnames = "value")

    # First window [0-100]: covers [0-50] with val=1, [100-100] with val=2
    # = (1*50) / 50 = 1 (only 50bp coverage in 100bp window)
    expect_equal(result$value[1], 1)

    # Second window [100-200]: covers [100-150] with val=2, [200-200] with val=3
    # = (2*50) / 50 = 2
    expect_equal(result$value[2], 2)

    gvtrack.rm("test.iter.vt")
})

test_that("value-based vtrack works across multiple chromosomes", {
    intervals_df <- data.frame(
        chrom = c("chr1", "chr1", "chr2", "chr2"),
        start = c(100, 300, 100, 300),
        end = c(200, 400, 200, 400),
        score = c(10, 20, 30, 40)
    )

    gvtrack.create("test.multichrom.vt", src = intervals_df, func = "avg")

    # Extract from chr1
    iter_int <- gintervals("chr1", 150, 350)
    result <- gextract("test.multichrom.vt", intervals = iter_int, iterator = iter_int, colnames = "value")
    expect_equal(result$value, 15) # (10*50 + 20*50) / 100

    # Extract from chr2
    iter_int <- gintervals("chr2", 150, 350)
    result <- gextract("test.multichrom.vt", intervals = iter_int, iterator = iter_int, colnames = "value")
    expect_equal(result$value, 35) # (30*50 + 40*50) / 100

    gvtrack.rm("test.multichrom.vt")
})

test_that("value-based vtrack handles NA values", {
    intervals_df <- data.frame(
        chrom = "chr1",
        start = c(100, 200, 300),
        end = c(150, 250, 350),
        score = c(10, NA, 30)
    )

    gvtrack.create("test.na.vt", src = intervals_df, func = "avg")

    # Extract over all intervals - NA should be ignored
    iter_int <- gintervals("chr1", 100, 350)
    result <- gextract("test.na.vt", intervals = iter_int, iterator = iter_int, colnames = "value")
    # Only 10 and 30 count: (10*50 + 30*50) / 100 = 20
    expect_equal(result$value, 20)

    gvtrack.rm("test.na.vt")
})

test_that("value-based vtrack works in track expressions", {
    intervals_df <- data.frame(
        chrom = "chr1",
        start = c(100, 300),
        end = c(200, 400),
        score = c(10, 20)
    )

    gvtrack.create("test.expr.vt", src = intervals_df, func = "avg")

    # Use in expression
    iter_int <- gintervals("chr1", 150, 350)
    result <- gextract("test.expr.vt * 2", intervals = iter_int, iterator = iter_int, colnames = "value")
    expect_equal(result$value, 30) # 15 * 2

    # Use in expression with another vtrack
    gvtrack.create("test.expr2.vt", src = intervals_df, func = "max")
    result <- gextract("test.expr.vt + test.expr2.vt", intervals = iter_int, iterator = iter_int, colnames = "value")
    expect_equal(result$value, 35) # 15 + 20

    gvtrack.rm("test.expr.vt")
    gvtrack.rm("test.expr2.vt")
})

test_that("value-based vtrack exists function works", {
    intervals_df <- data.frame(
        chrom = "chr1",
        start = c(100, 300),
        end = c(200, 400),
        score = c(10, 20)
    )

    gvtrack.create("test.exists.vt", src = intervals_df, func = "exists")

    # Interval with data
    iter_int <- gintervals("chr1", 150, 350)
    result <- gextract("test.exists.vt", intervals = iter_int, iterator = iter_int, colnames = "value")
    expect_equal(result$value, 1)

    # Interval without data
    iter_int <- gintervals("chr1", 500, 600)
    result <- gextract("test.exists.vt", intervals = iter_int, iterator = iter_int, colnames = "value")
    expect_equal(result$value, 0)

    gvtrack.rm("test.exists.vt")
})

test_that("value-based vtrack size function works", {
    intervals_df <- data.frame(
        chrom = "chr1",
        start = c(100, 300, 500),
        end = c(200, 400, 600),
        score = c(10, 20, 30)
    )

    gvtrack.create("test.size.vt", src = intervals_df, func = "size")

    # Interval covering 2 data intervals
    iter_int <- gintervals("chr1", 150, 350)
    result <- gextract("test.size.vt", intervals = iter_int, iterator = iter_int, colnames = "value")
    expect_equal(result$value, 2)

    # Interval covering all 3 data intervals
    iter_int <- gintervals("chr1", 0, 1000)
    result <- gextract("test.size.vt", intervals = iter_int, iterator = iter_int, colnames = "value")
    expect_equal(result$value, 3)

    gvtrack.rm("test.size.vt")
})

test_that("value-based vtrack first/last functions work", {
    intervals_df <- data.frame(
        chrom = "chr1",
        start = c(100, 300, 500),
        end = c(200, 400, 600),
        score = c(10, 20, 30)
    )

    gvtrack.create("test.first.vt", src = intervals_df, func = "first")
    gvtrack.create("test.last.vt", src = intervals_df, func = "last")

    # Interval covering all 3
    iter_int <- gintervals("chr1", 0, 1000)
    first_result <- gextract("test.first.vt", intervals = iter_int, iterator = iter_int, colnames = "value")
    last_result <- gextract("test.last.vt", intervals = iter_int, iterator = iter_int, colnames = "value")
    expect_equal(first_result$value, 10) # First value
    expect_equal(last_result$value, 30) # Last value

    gvtrack.rm("test.first.vt")
    gvtrack.rm("test.last.vt")
})

test_that("value-based vtrack nearest function works", {
    intervals_df <- data.frame(
        chrom = "chr1",
        start = c(100, 300),
        end = c(200, 400),
        score = c(10, 20)
    )

    gvtrack.create("test.nearest.vt", src = intervals_df, func = "nearest")

    # Interval with data - should return avg
    iter_int <- gintervals("chr1", 150, 350)
    result <- gextract("test.nearest.vt", intervals = iter_int, iterator = iter_int, colnames = "value")
    expect_equal(result$value, 15) # Average like sparse tracks

    gvtrack.rm("test.nearest.vt")
})

test_that("value-based vtrack position functions work", {
    intervals_df <- data.frame(
        chrom = "chr1",
        start = c(100, 300, 500),
        end = c(200, 400, 600),
        score = c(10, 30, 20) # Middle one is max
    )

    # Test absolute positions
    gvtrack.create("test.min.pos.abs.vt", src = intervals_df, func = "min.pos.abs")
    gvtrack.create("test.max.pos.abs.vt", src = intervals_df, func = "max.pos.abs")
    gvtrack.create("test.first.pos.abs.vt", src = intervals_df, func = "first.pos.abs")
    gvtrack.create("test.last.pos.abs.vt", src = intervals_df, func = "last.pos.abs")

    iter_int <- gintervals("chr1", 0, 1000)

    min_pos_abs <- gextract("test.min.pos.abs.vt", intervals = iter_int, iterator = iter_int, colnames = "value")
    max_pos_abs <- gextract("test.max.pos.abs.vt", intervals = iter_int, iterator = iter_int, colnames = "value")
    first_pos_abs <- gextract("test.first.pos.abs.vt", intervals = iter_int, iterator = iter_int, colnames = "value")
    last_pos_abs <- gextract("test.last.pos.abs.vt", intervals = iter_int, iterator = iter_int, colnames = "value")

    expect_equal(min_pos_abs$value, 100) # Position of min value (10)
    expect_equal(max_pos_abs$value, 300) # Position of max value (30)
    expect_equal(first_pos_abs$value, 100) # Position of first interval
    expect_equal(last_pos_abs$value, 500) # Position of last interval

    # Test relative positions
    gvtrack.create("test.min.pos.rel.vt", src = intervals_df, func = "min.pos.relative")
    gvtrack.create("test.max.pos.rel.vt", src = intervals_df, func = "max.pos.relative")

    iter_int <- gintervals("chr1", 50, 1000)

    min_pos_rel <- gextract("test.min.pos.rel.vt", intervals = iter_int, iterator = iter_int, colnames = "value")
    max_pos_rel <- gextract("test.max.pos.rel.vt", intervals = iter_int, iterator = iter_int, colnames = "value")

    expect_equal(min_pos_rel$value, 50) # 100 - 50
    expect_equal(max_pos_rel$value, 250) # 300 - 50

    # Clean up
    gvtrack.rm("test.min.pos.abs.vt")
    gvtrack.rm("test.max.pos.abs.vt")
    gvtrack.rm("test.first.pos.abs.vt")
    gvtrack.rm("test.last.pos.abs.vt")
    gvtrack.rm("test.min.pos.rel.vt")
    gvtrack.rm("test.max.pos.rel.vt")
})

test_that("intervals with value column work with interval-based summarizers even when overlapping", {
    # Create overlapping intervals with a value column
    # When using interval-based functions, the value column should be ignored
    # and intervals should be treated as regular intervals (overlaps allowed)
    intervals_df <- data.frame(
        chrom = "chr1",
        start = c(100, 150, 300), # First two overlap
        end = c(200, 250, 400),
        score = c(10, 20, 30) # Value column (should be ignored for interval functions)
    )

    # Test coverage function - should work with overlapping intervals
    gvtrack.create("test.coverage.vt", src = intervals_df, func = "coverage")
    iter_int <- gintervals("chr1", 100, 400)
    result <- gextract("test.coverage.vt", intervals = iter_int, iterator = iter_int, colnames = "value")
    # Coverage: [100-200] and [150-250] unified to [100-250] = 150bp, plus [300-400] = 100bp
    # Total = 250bp out of 300bp = 0.833...
    expect_equal(result$value, 250 / 300, tolerance = 1e-6)
    gvtrack.rm("test.coverage.vt")

    # Test neighbor.count function - should work with overlapping intervals
    gvtrack.create("test.neighbor.vt", src = intervals_df, func = "neighbor.count", params = 10)
    iter_int <- gintervals("chr1", 250, 260) # Between second and third interval
    result <- gextract("test.neighbor.vt", intervals = iter_int, iterator = iter_int, colnames = "value")
    # Should count both overlapping intervals if within distance
    expect_true(result$value >= 0) # Should not error
    gvtrack.rm("test.neighbor.vt")

    # Test distance function - should work with overlapping intervals
    gvtrack.create("test.distance.vt", src = intervals_df, func = "distance")
    iter_int <- gintervals("chr1", 170, 180) # Center at 175, between first two intervals
    result <- gextract("test.distance.vt", intervals = iter_int, iterator = iter_int, colnames = "value")
    # Should calculate distance to nearest interval center
    expect_true(!is.na(result$value))
    gvtrack.rm("test.distance.vt")

    # Test distance.center function - note: this function requires non-overlapping intervals
    # but the key point is that it's treated as an interval-based function, not value-based
    # So the overlap check in the value-based path is skipped (which is the fix)
    # However, distance.center itself will check for overlaps and error if they exist
    expect_error(
        gvtrack.create("test.distance.center.vt", src = intervals_df, func = "distance.center"),
        regexp = "overlapping"
    )
})

test_that("intervals with value column and interval-based functions ignore value column", {
    # Verify that when using interval-based functions, the value column is ignored
    # and the intervals are treated as regular intervals

    # Create intervals with overlapping regions and different values
    intervals_df <- data.frame(
        chrom = "chr1",
        start = c(100, 150),
        end = c(200, 250),
        score = c(999, 888) # Values that would be used if treated as value-based
    )

    # With coverage, both intervals should contribute regardless of values
    gvtrack.create("test.coverage.ignore.vt", src = intervals_df, func = "coverage")
    iter_int <- gintervals("chr1", 100, 250)
    result <- gextract("test.coverage.ignore.vt", intervals = iter_int, iterator = iter_int, colnames = "value")
    # Coverage: [100-200] and [150-250] unified to [100-250] = 150bp out of 150bp = 1.0
    expect_equal(result$value, 1.0, tolerance = 1e-6)
    gvtrack.rm("test.coverage.ignore.vt")
})
