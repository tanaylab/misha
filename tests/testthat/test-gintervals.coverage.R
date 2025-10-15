test_that("gintervals.covered_bp works with simple non-overlapping intervals", {
    intervs <- gintervals(c("chr1", "chr1", "chr2"), c(100, 500, 1000), c(200, 700, 2000))
    # Expected: (200-100) + (700-500) + (2000-1000) = 100 + 200 + 1000 = 1300
    expect_equal(gintervals.covered_bp(intervs), 1300)
})

test_that("gintervals.covered_bp works with overlapping intervals", {
    intervs <- gintervals(c("chr1", "chr1"), c(100, 150), c(200, 250))
    # Overlapping intervals: [100,200) and [150,250)
    # After canonicalization: [100,250) = 150 bp
    expect_equal(gintervals.covered_bp(intervs), 150)
})

test_that("gintervals.covered_bp works with touching intervals", {
    intervs <- gintervals(c("chr1", "chr1"), c(100, 200), c(200, 300))
    # Touching intervals [100,200) and [200,300)
    # After canonicalization with unify_touching: [100,300) = 200 bp
    expect_equal(gintervals.covered_bp(intervs), 200)
})

test_that("gintervals.covered_bp works with completely overlapping intervals", {
    intervs <- gintervals(c("chr1", "chr1", "chr1"), c(100, 150, 120), c(300, 200, 180))
    # [100,300), [150,200), [120,180) - first one covers all others
    # After canonicalization: [100,300) = 200 bp
    expect_equal(gintervals.covered_bp(intervs), 200)
})

test_that("gintervals.covered_bp returns 0 for NULL canonical result", {
    # Create intervals that will result in 0 coverage after canonicalization
    # This tests the internal logic rather than the edge case of empty input
    intervs <- data.frame(chrom = character(0), start = numeric(0), end = numeric(0))
    expect_equal(gintervals.covered_bp(intervs), 0)
})

test_that("gintervals.covered_bp works with single interval", {
    intervs <- gintervals("chr1", 1000, 2000)
    expect_equal(gintervals.covered_bp(intervs), 1000)
})

test_that("gintervals.covered_bp works with actual track data", {
    intervs <- gscreen("test.fixedbin > 0.2", gintervals(1, 0, -1))
    bp <- gintervals.covered_bp(intervs)
    expect_true(bp > 0)
    expect_true(is.numeric(bp))
})

test_that("gintervals.coverage_fraction works with simple intervals", {
    # Create intervals1 that covers half of intervals2
    intervs1 <- gintervals("chr1", 0, 500)
    intervs2 <- gintervals("chr1", 0, 1000)
    expect_equal(gintervals.coverage_fraction(intervs1, intervs2), 0.5)
})

test_that("gintervals.coverage_fraction works with complete coverage", {
    intervs1 <- gintervals(c("chr1", "chr1"), c(0, 500), c(500, 1000))
    intervs2 <- gintervals("chr1", 200, 800)
    # intervs1 covers [0,500) and [500,1000) which completely covers [200,800)
    expect_equal(gintervals.coverage_fraction(intervs1, intervs2), 1.0)
})

test_that("gintervals.coverage_fraction works with no coverage", {
    intervs1 <- gintervals("chr1", 0, 100)
    intervs2 <- gintervals("chr1", 200, 300)
    expect_equal(gintervals.coverage_fraction(intervs1, intervs2), 0)
})

test_that("gintervals.coverage_fraction works with partial overlap", {
    intervs1 <- gintervals("chr1", 50, 150)
    intervs2 <- gintervals("chr1", 100, 200)
    # Intersection is [100,150) = 50 bp out of 100 bp total
    expect_equal(gintervals.coverage_fraction(intervs1, intervs2), 0.5)
})

test_that("gintervals.coverage_fraction with NULL intervals2 uses whole genome", {
    intervs1 <- gscreen("test.fixedbin > 0.2", gintervals(c(1, 2), 0, -1))
    frac <- gintervals.coverage_fraction(intervs1)
    expect_true(frac >= 0 && frac <= 1)
    expect_true(is.numeric(frac))
})

test_that("gintervals.coverage_fraction works across multiple chromosomes", {
    intervs1 <- gintervals(c("chr1", "chr2"), c(0, 0), c(1000, 1000))
    intervs2 <- gintervals(c("chr1", "chr2"), c(0, 0), c(500, 2000))
    # chr1: [0,1000) covers [0,500) completely = 500 bp
    # chr2: [0,1000) covers [0,2000) partially = 1000 bp out of 2000 bp
    # Total: 1500 bp covered out of 2500 bp = 0.6
    expect_equal(gintervals.coverage_fraction(intervs1, intervs2), 0.6)
})

test_that("gintervals.coverage_fraction handles overlapping intervals1", {
    # intervals1 has overlaps that should be unified
    intervs1 <- gintervals(c("chr1", "chr1"), c(0, 50), c(100, 150))
    intervs2 <- gintervals("chr1", 0, 200)
    # After unification: [0,150) covers 150 bp out of 200 bp
    expect_equal(gintervals.coverage_fraction(intervs1, intervs2), 0.75)
})

test_that("gintervals.coverage_fraction handles overlapping intervals2", {
    intervs1 <- gintervals("chr1", 0, 100)
    # intervals2 has overlaps that should be unified
    intervs2 <- gintervals(c("chr1", "chr1"), c(0, 50), c(150, 200))
    # intervals2 after unification: [0,200) = 200 bp
    # intervs1 covers [0,100) = 100 bp
    expect_equal(gintervals.coverage_fraction(intervs1, intervs2), 0.5)
})

test_that("gintervals.coverage_fraction returns 0 when intervals2 is empty", {
    intervs1 <- gintervals("chr1", 0, 100)
    intervs2 <- data.frame(chrom = character(0), start = numeric(0), end = numeric(0))
    expect_equal(gintervals.coverage_fraction(intervs1, intervs2), 0)
})

test_that("gintervals.coverage_fraction returns 0 when intervals1 is empty", {
    intervs1 <- data.frame(chrom = character(0), start = numeric(0), end = numeric(0))
    intervs2 <- gintervals("chr1", 0, 100)
    expect_equal(gintervals.coverage_fraction(intervs1, intervs2), 0)
})

test_that("gintervals.covered_bp and coverage_fraction work together", {
    intervs1 <- gintervals(c("chr1", "chr1"), c(100, 150), c(300, 400))
    intervs2 <- gintervals("chr1", 0, 1000)

    bp1 <- gintervals.covered_bp(intervs1)
    bp2 <- gintervals.covered_bp(intervs2)
    frac <- gintervals.coverage_fraction(intervs1, intervs2)

    # Verify the relationship
    expect_equal(frac, bp1 / bp2)
})

test_that("gintervals.coverage_fraction works with real track data", {
    intervs1 <- gscreen("test.fixedbin > 0.3", gintervals(c(1, 2), 0, -1))
    intervs2 <- gscreen("test.fixedbin > 0.1", gintervals(c(1, 2), 0, -1))

    frac <- gintervals.coverage_fraction(intervs1, intervs2)

    expect_true(frac >= 0 && frac <= 1)
    expect_true(is.numeric(frac))
    # intervs1 is more restrictive, so it should cover less than the full intervs2
    expect_true(frac < 1)
})

test_that("gintervals.coverage_fraction is commutative for symmetric sets", {
    intervs <- gintervals(c("chr1", "chr2"), c(100, 200), c(300, 400))

    # Coverage of intervs by intervs should be 1
    expect_equal(gintervals.coverage_fraction(intervs, intervs), 1)
})

test_that("error handling: gintervals.covered_bp requires intervals", {
    expect_error(gintervals.covered_bp(), "Usage")
})

test_that("error handling: gintervals.coverage_fraction requires intervals1", {
    expect_error(gintervals.coverage_fraction(), "Usage")
})
