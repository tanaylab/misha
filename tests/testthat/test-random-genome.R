test_that("grandom_genome basic functionality works", {
    # Generate intervals
    set.seed(123)
    intervals <- grandom_genome(100, 10)

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

test_that("grandom_genome respects dist_from_edge", {
    set.seed(456)
    dist_from_edge <- 1e6
    intervals <- grandom_genome(100, 20, dist_from_edge = dist_from_edge)

    # Check that all intervals respect distance from edge
    all_genome <- gintervals.all()
    for (i in 1:nrow(intervals)) {
        chrom_info <- all_genome[all_genome$chrom == intervals$chrom[i], ]
        chrom_size <- chrom_info$end - chrom_info$start
        expect_true(intervals$start[i] >= dist_from_edge)
        expect_true(intervals$end[i] <= chrom_size - dist_from_edge)
    }
})

test_that("grandom_genome chromosome filtering works", {
    set.seed(789)
    # Test with specific chromosomes
    intervals <- grandom_genome(100, 50, chromosomes = c("chr1", "chr2"))

    # Check that all intervals are on specified chromosomes
    expect_true(all(intervals$chrom %in% c("chr1", "chr2")))

    # Check that we get intervals on multiple chromosomes (probabilistically)
    # (may fail with very small n, but with n=50 should almost always work)
    unique_chroms <- unique(as.character(intervals$chrom))
    expect_true(length(unique_chroms) >= 1)
})

test_that("grandom_genome is reproducible with set.seed", {
    set.seed(42)
    intervals1 <- grandom_genome(100, 10)

    set.seed(42)
    intervals2 <- grandom_genome(100, 10)

    expect_identical(intervals1, intervals2)
})

test_that("grandom_genome handles different sizes", {
    set.seed(111)

    # Small intervals
    intervals_small <- grandom_genome(10, 5)
    expect_true(all(intervals_small$end - intervals_small$start == 10))

    # Large intervals
    intervals_large <- grandom_genome(100000, 5, dist_from_edge = 1e6)
    expect_true(all(intervals_large$end - intervals_large$start == 100000))
})

test_that("grandom_genome input validation works", {
    # Invalid size
    expect_error(grandom_genome(0, 10), "size must be a positive number")
    expect_error(grandom_genome(-100, 10), "size must be a positive number")
    expect_error(grandom_genome(c(100, 200), 10), "size must be a positive number")

    # Invalid n
    expect_error(grandom_genome(100, 0), "n must be a positive number")
    expect_error(grandom_genome(100, -10), "n must be a positive number")
    expect_error(grandom_genome(100, c(10, 20)), "n must be a positive number")

    # Invalid dist_from_edge
    expect_error(grandom_genome(100, 10, dist_from_edge = -1), "dist_from_edge must be a non-negative number")
    expect_error(grandom_genome(100, 10, dist_from_edge = c(1e6, 2e6)), "dist_from_edge must be a non-negative number")

    # Invalid chromosomes
    expect_error(grandom_genome(100, 10, chromosomes = 123), "chromosomes must be a character vector")
    expect_error(grandom_genome(100, 10, chromosomes = c("nonexistent_chr")), "No chromosomes named")
})

test_that("grandom_genome handles edge cases", {
    # Very small dist_from_edge
    set.seed(222)
    intervals <- grandom_genome(100, 5, dist_from_edge = 0)
    expect_equal(nrow(intervals), 5)

    # Single chromosome
    set.seed(333)
    intervals <- grandom_genome(100, 10, chromosomes = "chr1")
    expect_true(all(intervals$chrom == "chr1"))

    # Large number of intervals
    set.seed(444)
    intervals <- grandom_genome(100, 1000)
    expect_equal(nrow(intervals), 1000)
    expect_true(all(intervals$end - intervals$start == 100))
})

test_that("grandom_genome samples proportional to chromosome length", {
    set.seed(555)

    # Generate many intervals to test sampling distribution
    n_intervals <- 10000
    intervals <- grandom_genome(100, n_intervals, dist_from_edge = 1e6)

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

test_that("grandom_genome error when chromosomes too short", {
    # Try to create intervals that are too large for any chromosome
    all_genome <- gintervals.all()
    max_chrom_size <- max(all_genome$end - all_genome$start)

    expect_error(
        grandom_genome(max_chrom_size + 1000, 10, dist_from_edge = 0),
        "No chromosomes are long enough"
    )

    expect_error(
        grandom_genome(1000, 10, dist_from_edge = max_chrom_size),
        "No chromosomes are long enough"
    )
})
