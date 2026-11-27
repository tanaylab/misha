create_isolated_test_db()

test_that("gintervals.promoters takes the window from the 5' base of each strand", {
    tss <- data.frame(
        chrom = "chr1",
        start = c(10000, 20000, 30000),
        end = c(10001, 20001, 30001),
        strand = c(1, -1, 0),
        name = c("p", "m", "both")
    )
    res <- gintervals.promoters(tss, upstream = 500, downstream = 50)
    # +: [tss - 500, tss + 50); -: 5' base is 20000, so [20001 - 50, 20001 + 500)
    # 0: covers both windows
    expect_equal(res$start, c(9500, 19951, 29500))
    expect_equal(res$end, c(10050, 20501, 30501))
    expect_equal(res$strand, c(1, -1, 0))
    expect_equal(res$name, c("p", "m", "both"))
})

test_that("gintervals.promoters uses the 5' end of longer intervals", {
    tx <- data.frame(chrom = "chr1", start = c(1000, 1000), end = c(5000, 5000), strand = c(1, -1))
    res <- gintervals.promoters(tx, upstream = 100, downstream = 10)
    expect_equal(res$start, c(900, 4990))
    expect_equal(res$end, c(1010, 5100))
})

test_that("gintervals.promoters clips at chromosome ends", {
    chr1_end <- gintervals.all()$end[1]
    tss <- data.frame(chrom = "chr1", start = c(100, chr1_end - 1), end = c(101, chr1_end), strand = c(1, -1))
    res <- gintervals.promoters(tss, upstream = 500, downstream = 50)
    expect_equal(res$start, c(0, chr1_end - 50))
    expect_equal(res$end, c(150, chr1_end))
})

test_that("gintervals.promoters requires strand 1, -1 or 0", {
    expect_error(gintervals.promoters(gintervals(1, 1000, 1001)), "strand column")
    bad <- data.frame(chrom = "chr1", start = 1000, end = 1001, strand = 2)
    expect_error(gintervals.promoters(bad), "strand must be")
    bad$strand <- NA
    expect_error(gintervals.promoters(bad), "strand must be")
})

test_that("gintervals.expand adds bp on both sides and clips", {
    res <- gintervals.expand(gintervals(1, c(1000, 5000), c(2000, 6000)), 100)
    expect_equal(res$start, c(900, 4900))
    expect_equal(res$end, c(2100, 6100))

    res <- gintervals.expand(gintervals(1, c(1000, 5000), c(2000, 6000)), c(10, 2000))
    expect_equal(res$start, c(990, 3000))
    expect_equal(res$end, c(2010, 8000))

    chr1_end <- gintervals.all()$end[1]
    res <- gintervals.expand(gintervals(1, c(50, chr1_end - 50), c(100, chr1_end)), 100)
    expect_equal(res$start, c(0, chr1_end - 150))
    expect_equal(res$end, c(200, chr1_end))

    expect_error(gintervals.expand(gintervals(1, 1000, 2000)), "Usage")
})

test_that("gintervals.remove_overlaps keeps one interval per overlap group in row order", {
    intervs <- data.frame(
        chrom = "chr1",
        start = c(11000, 100, 10000, 10500),
        end = c(12000, 200, 13000, 10600),
        data = c(10, 20, 30, 40)
    )
    expect_equal(gintervals.remove_overlaps(intervs)$data, c(10, 20))
    expect_equal(gintervals.remove_overlaps(intervs, select = "last")$data, c(20, 40))

    set.seed(1)
    picked <- replicate(20, gintervals.remove_overlaps(intervs, select = "random")$data, simplify = FALSE)
    expect_true(all(vapply(picked, length, integer(1)) == 2))
    expect_true(all(vapply(picked, function(x) 20 %in% x, logical(1))))
    expect_setequal(unique(unlist(picked)), c(10, 20, 30, 40))
})

test_that("gintervals.remove_overlaps follows unify_touching_intervals", {
    touching <- data.frame(chrom = "chr1", start = c(0, 100), end = c(100, 200))
    expect_equal(nrow(gintervals.remove_overlaps(touching)), 1)
    expect_equal(nrow(gintervals.remove_overlaps(touching, unify_touching_intervals = FALSE)), 2)
})
