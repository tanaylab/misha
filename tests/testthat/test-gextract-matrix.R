create_isolated_test_db()

# gextract output reshaped to intervals x bins x expressions, ordered by
# intervalID and start
gxmat_reference <- function(exprs, intervals, iterator = 1) {
    r <- gextract(exprs, intervals = intervals, iterator = iterator, colnames = paste0("v", seq_along(exprs)))
    r <- r[order(r$intervalID, r$start), ]
    num_bins <- (intervals$end[1] - intervals$start[1]) / iterator
    expect_equal(nrow(r), nrow(intervals) * num_bins)
    res <- array(NA_real_, dim = c(nrow(intervals), num_bins, length(exprs)))
    for (i in seq_along(exprs)) {
        res[, , i] <- matrix(r[[paste0("v", i)]], nrow = nrow(intervals), byrow = TRUE)
    }
    res
}

# unsorted windows on several chromosomes, with duplicated and overlapping rows
gxmat_windows <- function(width, num = 40, align = 1) {
    set.seed(17)
    chroms <- sample(c("chr1", "chr2", "chrX"), num, replace = TRUE)
    starts <- sample(0:2000, num, replace = TRUE) * align
    intervs <- data.frame(chrom = chroms, start = starts, end = starts + width)
    dup <- intervs[c(3, 3, 7), ]
    ovl <- intervs[c(5, 9), ]
    ovl$start <- ovl$start + align
    ovl$end <- ovl$end + align
    intervs <- rbind(intervs, dup, ovl)
    intervs[sample(nrow(intervs)), ]
}

test_that("gextract.matrix matches gextract per base on unsorted, duplicated and overlapping rows", {
    intervs <- gxmat_windows(200)
    ref <- gxmat_reference("test.fixedbin", intervs)
    m <- gextract.matrix("test.fixedbin", intervals = intervs)
    expect_true(is.matrix(m))
    expect_equal(dim(m), c(nrow(intervs), 200))
    expect_identical(m, ref[, , 1])
    withr::with_options(list(gmultitasking = FALSE), {
        expect_identical(gextract.matrix("test.fixedbin", intervals = intervs), ref[, , 1])
    })
    # positional intervals, as in gextract
    expect_identical(gextract.matrix("test.fixedbin", intervs), m)
})

test_that("gextract.matrix matches gextract with iterator > 1", {
    intervs <- gxmat_windows(500, align = 50)
    ref <- gxmat_reference("test.fixedbin", intervs, iterator = 50)
    m <- gextract.matrix("test.fixedbin", intervals = intervs, iterator = 50)
    expect_equal(dim(m), c(nrow(intervs), 10))
    expect_identical(m, ref[, , 1])
})

test_that("gextract.matrix keeps the NaN values of a sparse track", {
    intervs <- gxmat_windows(300)
    m <- gextract.matrix("test.sparse", intervals = intervs)
    expect_true(any(is.nan(m)))
    expect_true(any(!is.na(m)))
    expect_identical(m, gxmat_reference("test.sparse", intervs)[, , 1])
})

test_that("gextract.matrix gives a coverage mask from an intervals-based virtual track", {
    gvtrack.create("gxmat_cov", gscreen("test.fixedbin > 0.2", gintervals(c(1, 2, "X"))), func = "coverage")
    withr::defer(gvtrack.rm("gxmat_cov"))
    intervs <- gxmat_windows(250)
    m <- gextract.matrix("gxmat_cov", intervals = intervs)
    expect_setequal(unique(as.vector(m)), c(0, 1))
    expect_identical(m, gxmat_reference("gxmat_cov", intervs)[, , 1])

    intervs10 <- gxmat_windows(250, align = 10)
    m10 <- gextract.matrix("gxmat_cov", intervals = intervs10, iterator = 10)
    expect_identical(m10, gxmat_reference("gxmat_cov", intervs10, iterator = 10)[, , 1])
})

test_that("gextract.matrix returns an array for several expressions", {
    exprs <- c("test.fixedbin", "test.sparse", "test.fixedbin * 2")
    intervs <- gxmat_windows(100)
    ref <- gxmat_reference(exprs, intervs)
    a <- gextract.matrix("test.fixedbin", "test.sparse", "test.fixedbin * 2", intervals = intervs)
    expect_equal(dim(a), c(nrow(intervs), 100, 3))
    expect_equal(dimnames(a), list(NULL, NULL, exprs))
    expect_identical(unname(a), ref)

    a2 <- gextract.matrix(exprs, intervals = intervs, colnames = c("a", "b", "c"))
    expect_equal(dimnames(a2)[[3]], c("a", "b", "c"))
    expect_identical(unname(a2), ref)
    withr::with_options(list(gmultitasking = FALSE), {
        expect_identical(unname(gextract.matrix(exprs, intervals = intervs)), ref)
    })
})

test_that("gextract.matrix gives the same result in many chunks", {
    exprs <- c("test.fixedbin", "test.sparse")
    intervs <- gxmat_windows(100)
    ref <- gxmat_reference(exprs, intervs)
    # 3 intervals per gextract call
    withr::with_options(list(gmax.data.size = 350), {
        expect_identical(unname(gextract.matrix(exprs, intervals = intervs)), ref)
    })
})

test_that("gextract.matrix ignores strand and accepts an intervals set name", {
    intervs <- gxmat_windows(100)
    m <- gextract.matrix("test.fixedbin", intervals = intervs)
    stranded <- intervs
    stranded$strand <- -1
    expect_identical(gextract.matrix("test.fixedbin", intervals = stranded), m)

    gintervals.save("test.gxmat_windows", intervs)
    withr::defer(gintervals.rm("test.gxmat_windows", force = TRUE))
    saved <- gintervals.load("test.gxmat_windows")
    expect_identical(
        gextract.matrix("test.fixedbin", intervals = "test.gxmat_windows"),
        gextract.matrix("test.fixedbin", intervals = saved)
    )
})

test_that("gextract.matrix rejects invalid input", {
    intervs <- gintervals(c(1, 2), c(100, 200), c(200, 300))
    expect_error(gextract.matrix("test.fixedbin", intervals = gintervals(c(1, 2), c(100, 200), c(200, 350))), "same positive width")
    expect_error(gextract.matrix("test.fixedbin", intervals = intervs, iterator = 30), "does not divide")
    expect_error(gextract.matrix("test.fixedbin", intervals = intervs, iterator = "test.fixedbin"), "positive integer")
    expect_error(gextract.matrix("test.fixedbin", intervals = intervs, iterator = 2.5), "positive integer")
    expect_error(gextract.matrix("test.fixedbin", intervals = intervs, iterator = intervs), "positive integer")
    expect_error(gextract.matrix("test.fixedbin", intervals = gintervals(1, 105, 205), iterator = 10), "multiple of iterator")
    expect_error(gextract.matrix("test.fixedbin", intervals = intervs, colnames = c("a", "b")), "column names")
    expect_error(gextract.matrix("test.fixedbin", intervals = intervs[0, ]), "empty")
    expect_error(gextract.matrix("test.fixedbin", intervals = gintervals.2d(1, 0, 100, 1, 0, 100)), "1D")
    chrom_end <- gintervals.all()$end[gintervals.all()$chrom == "chr1"]
    expect_error(gextract.matrix("test.fixedbin", intervals = gintervals(1, chrom_end - 50, chrom_end + 50)), "chromosome boundaries")
})
