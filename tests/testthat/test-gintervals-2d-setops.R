test_that("gintervals.2d.intersect clips overlapping rectangles pairwise", {
    a <- gintervals.2d("chr1", 0, 100, "chr2", 0, 100)
    b <- gintervals.2d("chr1", 50, 150, "chr2", 50, 150)
    res <- gintervals.2d.intersect(a, b)
    expect_equal(nrow(res), 1)
    expect_equal(res$start1, 50)
    expect_equal(res$end1, 100)
    expect_equal(res$start2, 50)
    expect_equal(res$end2, 100)
    expect_equal(colnames(res)[1:6], c("chrom1", "start1", "end1", "chrom2", "start2", "end2"))
})

test_that("gintervals.2d.intersect emits nothing unless both dimensions overlap", {
    a <- gintervals.2d("chr1", 0, 100, "chr2", 0, 100)
    b <- gintervals.2d("chr1", 50, 150, "chr2", 200, 300) # dim2 disjoint
    res <- gintervals.2d.intersect(a, b)
    expect_true(is.null(res) || nrow(res) == 0)
})

test_that("gintervals.2d.intersect only pairs matching chrom-pairs", {
    a <- gintervals.2d("chr1", 0, 100, "chr2", 0, 100)
    b <- gintervals.2d("chr1", 0, 100, "chr3", 0, 100) # different chrom2
    res <- gintervals.2d.intersect(a, b)
    expect_true(is.null(res) || nrow(res) == 0)
})

test_that("gintervals.2d.intersect does a pairwise cartesian within a chrom-pair", {
    a <- gintervals.2d(
        c("chr1", "chr1"), c(0, 200), c(100, 300),
        c("chr2", "chr2"), c(0, 200), c(100, 300)
    )
    b <- gintervals.2d("chr1", 50, 250, "chr2", 50, 250)
    res <- gintervals.2d.intersect(a, b)
    expect_equal(nrow(res), 2)
    expect_setequal(res$start1, c(50, 200))
    expect_setequal(res$end1, c(100, 250))
})

test_that("gintervals.2d.union concatenates and sorts without merging overlaps", {
    a <- gintervals.2d("chr1", 0, 100, "chr2", 0, 100)
    b <- gintervals.2d("chr1", 50, 150, "chr2", 50, 150)
    res <- gintervals.2d.union(a, b)
    expect_equal(nrow(res), 2)
    expect_equal(colnames(res)[1:6], c("chrom1", "start1", "end1", "chrom2", "start2", "end2"))
})

test_that("gintervals.2d.union sorts by (chrom1, start1, chrom2, start2)", {
    a <- gintervals.2d("chr1", 200, 300, "chr2", 0, 100)
    b <- gintervals.2d("chr1", 0, 100, "chr2", 0, 100)
    res <- gintervals.2d.union(a, b)
    expect_equal(res$start1, c(0, 200))
})

test_that("2D set ops reject 1D intervals", {
    iv1 <- gintervals("chr1", 0, 100)
    expect_error(gintervals.2d.intersect(iv1, iv1), "2D")
    expect_error(gintervals.2d.union(iv1, iv1), "2D")
})
