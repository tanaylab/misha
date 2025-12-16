create_isolated_test_db()

test_that("gintervals.normalize works with basic intervals", {
    # Test basic normalization with even size
    intervs <- gintervals(1, c(1000, 5000), c(2000, 6000))
    result <- gintervals.normalize(intervs, 500)

    # Expected: center of first interval is 1500, so normalized interval should be [1250, 1750]
    # Expected: center of second interval is 5500, so normalized interval should be [5250, 5750]
    expect_equal(nrow(result), 2)
    expect_equal(result$chrom, factor(c("chr1", "chr1"), levels = gintervals.all()$chrom))
    expect_equal(result$start, c(1250, 5250))
    expect_equal(result$end, c(1750, 5750))
})

test_that("gintervals.normalize works with odd size", {
    # Test normalization with odd size - should give exact size
    intervs <- gintervals(1, c(1000, 5000), c(3000, 8000))
    result <- gintervals.normalize(intervs, 499)

    # Expected: center of first interval is 2000, size=499 → [1751, 2250] (exactly 499bp)
    # Expected: center of second interval is 6500, size=499 → [6251, 6750] (exactly 499bp)
    expect_equal(nrow(result), 2)
    expect_equal(result$chrom, factor(c("chr1", "chr1"), levels = gintervals.all()$chrom))
    expect_equal(result$start, c(1751, 6251))
    expect_equal(result$end, c(2250, 6750))
    # Verify exact size
    expect_equal(result$end - result$start, c(499, 499))
})

test_that("gintervals.normalize handles chromosome boundaries", {
    # Test intervals near chromosome boundaries
    intervs <- gintervals(1, c(0, 247249700), c(100, 247249719))
    result <- gintervals.normalize(intervs, 1000)

    expect_equal(nrow(result), 2)
    expect_equal(result$chrom, factor(c("chr1", "chr1"), levels = gintervals.all()$chrom))
    expect_equal(result$start, c(0, 247249209))
    expect_equal(result$end, c(550, 247249719))
})

test_that("gintervals.normalize respects chromosome limits", {
    # Test that intervals don't exceed chromosome boundaries
    intervs <- gintervals(1, c(0, 247249700), c(50, 247249719))
    result <- gintervals.normalize(intervs, 200)

    # Check that no intervals exceed chromosome boundaries
    chrom_size <- 247249719 # chr1 size in test database
    expect_true(all(result$start >= 0))
    expect_true(all(result$end <= chrom_size))

    # Expected: center of first interval is 25, expansion is 100, but start cannot be negative
    # Expected: center of second interval is 247249709, expansion is 100, but end cannot exceed chromosome size
    expect_equal(nrow(result), 2)
    expect_equal(result$chrom, factor(c("chr1", "chr1"), levels = gintervals.all()$chrom))
    expect_equal(result$start, c(0, 247249609)) # max(25-100, 0) and max(247249709-100, 0)
    expect_equal(result$end, c(125, 247249719)) # min(25+100, chrom_size) and min(247249709+100, chrom_size)
})

test_that("gintervals.normalize handles empty intervals", {
    # Test with empty intervals set
    intervs <- data.frame(chrom = character(0), start = numeric(0), end = numeric(0))
    result <- gintervals.normalize(intervs, 500)
    expect_null(result)
})

test_that("gintervals.normalize works with multiple chromosomes", {
    # Test with intervals on multiple chromosomes
    intervs <- gintervals(c(1, 2), c(1000, 2000), c(2000, 3000))
    result <- gintervals.normalize(intervs, 800)

    # Expected: center of first interval is 1500, expansion is 400, so normalized interval should be [1100, 1900]
    # Expected: center of second interval is 2500, expansion is 400, so normalized interval should be [2100, 2900]
    expect_equal(nrow(result), 2)
    expect_equal(result$chrom, factor(c("chr1", "chr2"), levels = gintervals.all()$chrom))
    expect_equal(result$start, c(1100, 2100))
    expect_equal(result$end, c(1900, 2900))
})

test_that("gintervals.normalize preserves interval names", {
    # Test that additional columns are preserved
    intervs <- gintervals(c(1, 1), c(1000, 5000), c(2000, 6000), c(1, -1))
    intervs$name <- c("interval1", "interval2")
    intervs$name2 <- c("savta", "saba")
    result <- gintervals.normalize(intervs, 600)
    expect_true("name" %in% colnames(result))
    expect_true("strand" %in% colnames(result))
    expect_equal(result$name, c("interval1", "interval2"))
    expect_equal(result$name2, c("savta", "saba"))
    expect_equal(result$strand, c(1, -1))
    expect_equal(colnames(result), colnames(intervs))

    # Also check that the normalized intervals are correct
    # Expected: center of first interval is 1500, expansion is 300, so normalized interval should be [1200, 1800]
    # Expected: center of second interval is 5500, expansion is 300, so normalized interval should be [5200, 5800]
    expect_equal(result$start, c(1200, 5200))
    expect_equal(result$end, c(1800, 5800))
})

test_that("gintervals.normalize rejects invalid size parameters", {
    intervs <- gintervals(1, 1000, 2000)

    # Test negative size
    expect_error(gintervals.normalize(intervs, -100))

    # Test zero size
    expect_error(gintervals.normalize(intervs, 0))

    # Test non-numeric size
    expect_error(gintervals.normalize(intervs, "500"))

    # Test mismatched multiple sizes (1 interval with 2 sizes is now allowed as one-to-many)
    # But 1 interval with 3 sizes that are mismatched should work (one-to-many)
    # So we test a different mismatch case: 2 intervals with 3 sizes
    intervs2 <- gintervals(1, c(1000, 5000), c(2000, 6000))
    expect_error(gintervals.normalize(intervs2, c(500, 600, 700)))
})

test_that("gintervals.normalize rejects 2D intervals", {
    # Test that 2D intervals are rejected
    intervs2d <- gintervals.2d(1, 1000, 2000, 1, 3000, 4000)
    expect_error(gintervals.normalize(intervs2d, 500))
})

test_that("gintervals.normalize with very small size", {
    # Test with very small normalization size
    intervs <- gintervals(1, c(1000, 5000), c(2000, 6000))
    result <- gintervals.normalize(intervs, 10)

    # Expected: center of first interval is 1500, expansion is 5, so normalized interval should be [1495, 1505]
    # Expected: center of second interval is 5500, expansion is 5, so normalized interval should be [5495, 5505]
    expect_equal(nrow(result), 2)
    expect_equal(result$chrom, factor(c("chr1", "chr1"), levels = gintervals.all()$chrom))
    expect_equal(result$start, c(1495, 5495))
    expect_equal(result$end, c(1505, 5505))
})

test_that("gintervals.normalize with large size", {
    # Test with large normalization size
    intervs <- gintervals(1, c(1000000, 5000000), c(1001000, 5001000))
    result <- gintervals.normalize(intervs, 10000)

    # Expected: center of first interval is 1000500, expansion is 5000, so normalized interval should be [995500, 1005500]
    # Expected: center of second interval is 5000500, expansion is 5000, so normalized interval should be [4995500, 5005500]
    expect_equal(nrow(result), 2)
    expect_equal(result$chrom, factor(c("chr1", "chr1"), levels = gintervals.all()$chrom))
    expect_equal(result$start, c(995500, 4995500))
    expect_equal(result$end, c(1005500, 5005500))
})

test_that("gintervals.normalize works with intervals.set.out", {
    # Test saving result to intervals set
    temp_track_name <- paste0("test.tmptrack_", sample(1:1e9, 1))
    gintervals.rm(temp_track_name, force = TRUE)
    withr::defer(gintervals.rm(temp_track_name, force = TRUE))

    intervs <- gintervals(c(1, 2), c(1000, 2000), c(2000, 3000))
    gintervals.normalize(intervs, 800, intervals.set.out = temp_track_name)

    result <- gintervals.load(temp_track_name)

    # Expected: center of first interval is 1500, expansion is 400, so normalized interval should be [1100, 1900]
    # Expected: center of second interval is 2500, expansion is 400, so normalized interval should be [2100, 2900]
    expect_equal(nrow(result), 2)
    expect_equal(result$chrom, factor(c("chr1", "chr2"), levels = gintervals.all()$chrom))
    expect_equal(result$start, c(1100, 2100))
    expect_equal(result$end, c(1900, 2900))
})

test_that("gintervals.normalize maintains interval count", {
    # Test that the number of output intervals equals input intervals
    intervs <- data.frame(chrom = factor(c("chr1", "chr2", "chr1"), levels = gintervals.all()$chrom), start = c(1000, 2000, 5000), end = c(2000, 3000, 6000))
    result <- gintervals.normalize(intervs, 500)
    expect_equal(nrow(result), nrow(intervs))

    # Also verify the normalized intervals are correct
    # Expected: centers are 1500, 2500, 5500; expansion is 250
    expect_equal(result$chrom, factor(c("chr1", "chr2", "chr1"), levels = gintervals.all()$chrom))
    expect_equal(result$start, c(1250, 2250, 5250))
    expect_equal(result$end, c(1750, 2750, 5750))
})

# ===== Tests for vector of sizes =====

test_that("gintervals.normalize works with vector of sizes", {
    # Use intervals all on same chromosome, already sorted by start
    # Note: gintervals() sorts intervals by start position
    intervs <- gintervals(1, c(1000, 3000, 5000), c(2000, 4000, 6000))
    sizes <- c(500, 750, 1000)
    result <- gintervals.normalize(intervs, sizes)

    # Centers: 1500, 3500, 5500
    # Expansions: 250, 375, 500
    # Note: All on chr1, order matches sorted intervals
    expect_equal(nrow(result), 3)
    expect_equal(result$start, c(1250, 3125, 5000))
    expect_equal(result$end, c(1750, 3875, 6000))
})

test_that("gintervals.normalize vector of identical sizes matches single size", {
    intervs <- gintervals(1, c(1000, 5000), c(2000, 6000))

    # Use single size
    result_single <- gintervals.normalize(intervs, 500)

    # Use vector of same sizes
    result_vector <- gintervals.normalize(intervs, c(500, 500))

    expect_equal(result_single, result_vector)
})

test_that("gintervals.normalize one-to-many: single interval with multiple sizes", {
    interv <- gintervals(1, 1000, 2000)
    sizes <- c(500, 1000, 1500)
    result <- gintervals.normalize(interv, sizes)

    # Center: 1500
    # All intervals centered at 1500 with different sizes
    expect_equal(nrow(result), 3)
    expect_equal(result$start, c(1250, 1000, 750))
    expect_equal(result$end, c(1750, 2000, 2250))
})

test_that("gintervals.normalize rejects mismatched sizes", {
    intervs <- gintervals(1, c(1000, 5000, 3000), c(2000, 6000, 4000))

    # 3 intervals with 5 sizes should error
    expect_error(
        gintervals.normalize(intervs, c(500, 600, 700, 800, 900)),
        "must either match"
    )

    # 3 intervals with 2 sizes should error (neither 1 nor 3)
    expect_error(
        gintervals.normalize(intervs, c(500, 600)),
        "must either match"
    )
})

test_that("gintervals.normalize validates all size values are positive", {
    intervs <- gintervals(1, c(1000, 5000), c(2000, 6000))

    # Vector with negative value
    expect_error(
        gintervals.normalize(intervs, c(500, -100)),
        "positive"
    )

    # Vector with zero
    expect_error(
        gintervals.normalize(intervs, c(500, 0)),
        "positive"
    )
})

test_that("gintervals.normalize handles size=1", {
    intervs <- gintervals(1, c(1000, 5000), c(2000, 6000))
    result <- gintervals.normalize(intervs, 1)

    # Centers: 1500, 5500
    # Size=1 should create [center, center+1]
    expect_equal(nrow(result), 2)
    expect_equal(result$start, c(1500, 5500))
    expect_equal(result$end, c(1501, 5501))
})

test_that("gintervals.normalize handles size=1 in vector mode", {
    intervs <- gintervals(1, c(1000, 5000), c(2000, 6000))
    result <- gintervals.normalize(intervs, c(1, 500))

    # First interval: center=1500, size=1 -> [1500, 1501]
    # Second interval: center=5500, size=500 -> [5250, 5750]
    expect_equal(nrow(result), 2)
    expect_equal(result$start, c(1500, 5250))
    expect_equal(result$end, c(1501, 5750))
})

test_that("gintervals.normalize vector mode works across multiple chromosomes", {
    # Create intervals using data.frame to control exact order
    # Note: gintervals() would sort these, so we use data.frame directly
    intervs <- data.frame(
        chrom = factor(c("chr1", "chr1", "chr2"), levels = gintervals.all()$chrom),
        start = c(1000, 5000, 2000),
        end = c(2000, 6000, 3000)
    )
    sizes <- c(600, 400, 800) # Matches the sorted order of intervals
    result <- gintervals.normalize(intervs, sizes)

    # Verify each interval gets correct size
    expect_equal(nrow(result), 3)
    # All intervals already sorted: chr1 (1000-2000, size=600), chr1 (5000-6000, size=400), chr2 (2000-3000, size=800)
    # Centers: 1500, 5500, 2500
    # Expansions: 300, 200, 400
    expect_equal(result$start, c(1200, 5300, 2100))
    expect_equal(result$end, c(1800, 5700, 2900))
    expect_equal(as.character(result$chrom), c("chr1", "chr1", "chr2"))
})

test_that("gintervals.normalize preserves metadata with vector sizes", {
    intervs <- gintervals(c(1, 1), c(1000, 5000), c(2000, 6000), c(1, -1))
    intervs$name <- c("interval1", "interval2")
    intervs$score <- c(10.5, 20.3)
    result <- gintervals.normalize(intervs, c(600, 800))

    expect_true("name" %in% colnames(result))
    expect_true("strand" %in% colnames(result))
    expect_true("score" %in% colnames(result))
    expect_equal(result$name, c("interval1", "interval2"))
    expect_equal(result$strand, c(1, -1))
    expect_equal(result$score, c(10.5, 20.3))
    expect_equal(colnames(result), colnames(intervs))
})

test_that("gintervals.normalize preserves column order with vector sizes", {
    intervs <- gintervals(c(1, 1), c(1000, 5000), c(2000, 6000))
    intervs$col1 <- c("a", "b")
    intervs$col2 <- c(1, 2)
    intervs$col3 <- c(TRUE, FALSE)

    original_cols <- colnames(intervs)
    result <- gintervals.normalize(intervs, c(500, 600))

    expect_equal(colnames(result), original_cols)
})

test_that("gintervals.normalize vector mode respects chromosome boundaries", {
    # Test intervals near chromosome boundaries with different sizes
    intervs <- gintervals(1, c(0, 247249700), c(100, 247249719))
    sizes <- c(1000, 200)
    result <- gintervals.normalize(intervs, sizes)

    chrom_size <- 247249719 # chr1 size in test database
    expect_equal(nrow(result), 2)
    expect_true(all(result$start >= 0))
    expect_true(all(result$end <= chrom_size))

    # First interval: center=50, size=1000, expansion=500
    # -> [0, 550] (clamped at start)
    # Second interval: center=247249709.5, size=200, expansion=100
    # -> [247249609, 247249719] (clamped at end)
    expect_equal(result$start, c(0, 247249609))
    expect_equal(result$end, c(550, 247249719))
})

test_that("gintervals.normalize one-to-many preserves metadata", {
    interv <- gintervals(1, 1000, 2000, 1)
    interv$name <- "test_interval"
    interv$score <- 42.5

    result <- gintervals.normalize(interv, c(500, 1000))

    expect_equal(nrow(result), 2)
    expect_equal(result$name, c("test_interval", "test_interval"))
    expect_equal(result$score, c(42.5, 42.5))
    expect_equal(result$strand, c(1, 1))
})
