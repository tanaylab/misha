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
    # Test normalization with odd size
    intervs <- gintervals(1, c(1000, 5000), c(3000, 8000))
    result <- gintervals.normalize(intervs, 499)

    # Expected: center of first interval is 2000, expansion is 249, so normalized interval should be [1751, 2249]
    # Expected: center of second interval is 6500, expansion is 249, so normalized interval should be [6251, 6749]
    expect_equal(nrow(result), 2)
    expect_equal(result$chrom, factor(c("chr1", "chr1"), levels = gintervals.all()$chrom))
    expect_equal(result$start, c(1751, 6251))
    expect_equal(result$end, c(2249, 6749))
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

    # Test multiple sizes
    expect_error(gintervals.normalize(intervs, c(500, 600)))
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
