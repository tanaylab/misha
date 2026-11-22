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

# values of every bin of every window, from one gextract call per window with
# the window's bins, built here, as the iterator (the bins of one window touch
# but do not overlap, so gextract reports each one as it is)
gxmat_bin_reference <- function(exprs, intervals, iterator) {
    num_bins <- (intervals$end[1] - intervals$start[1]) / iterator
    res <- array(NA_real_, dim = c(nrow(intervals), num_bins, length(exprs)))
    for (i in seq_len(nrow(intervals))) {
        bin_starts <- intervals$start[i] + (seq_len(num_bins) - 1) * iterator
        bins <- data.frame(chrom = intervals$chrom[i], start = bin_starts, end = bin_starts + iterator)
        r <- gextract(exprs, intervals = intervals[i, c("chrom", "start", "end")], iterator = bins, colnames = paste0("v", seq_along(exprs)))
        r <- r[order(r$start), ]
        expect_equal(r$start, bin_starts)
        expect_equal(r$end, bin_starts + iterator)
        for (j in seq_along(exprs)) {
            res[i, , j] <- r[[paste0("v", j)]]
        }
    }
    res
}

test_that("gextract.matrix counts bins from each interval's start", {
    gvtrack.create("gxmat_cov", gscreen("test.fixedbin > 0.2", gintervals(c(1, 2, "X"), 0, 3000)), func = "coverage")
    withr::defer(gvtrack.rm("gxmat_cov"))
    exprs <- c("test.fixedbin", "test.sparse", "gxmat_cov")
    # stored values of the dense track, one per 50 bp track bin
    dense <- lapply(c("chr1", "chr2", "chrX"), function(chrom) gextract("test.fixedbin", gintervals(chrom, 0, 3000))$test.fixedbin)
    names(dense) <- c("chr1", "chr2", "chrX")
    intervs <- gxmat_windows(350)
    for (k in c(7, 10, 50)) {
        expect_true(any(intervs$start %% k != 0))
        ref <- withr::with_options(list(gmultitasking = FALSE), gxmat_bin_reference(exprs, intervs, k))
        a <- gextract.matrix(exprs, intervals = intervs, iterator = k)
        expect_equal(dim(a), c(nrow(intervs), 350 / k, 3))
        expect_identical(unname(a), ref)
        expect_identical(gextract.matrix("test.sparse", intervals = intervs, iterator = k), ref[, , 2])
        expect_true(any(is.nan(a[, , 2])) && any(!is.na(a[, , 2])))
        expect_true(all(a[, , 3] >= 0 & a[, , 3] <= 1) && any(a[, , 3] > 0 & a[, , 3] < 1))
        # three intervals per gextract call
        withr::with_options(list(gmax.data.size = 3 * 350 / k), {
            expect_identical(unname(gextract.matrix(exprs, intervals = intervs, iterator = k)), ref)
        })
        withr::with_options(list(gmultitasking = FALSE), {
            expect_identical(unname(gextract.matrix(exprs, intervals = intervs, iterator = k)), ref)
        })
        # a bin of the dense track gets the plain mean of the track bins it overlaps
        bin_starts <- outer(intervs$start, (seq_len(350 / k) - 1) * k, "+")
        direct <- matrix(mapply(function(chrom, s) {
            mean(dense[[chrom]][(s %/% 50):((s + k - 1) %/% 50) + 1], na.rm = TRUE)
        }, rep(intervs$chrom, 350 / k), as.vector(bin_starts)), nrow = nrow(intervs))
        expect_equal(a[, , 1], direct, tolerance = 1e-6)
    }
})

test_that("gextract.matrix gives gextract's numeric-iterator bins when every start is on the grid", {
    exprs <- c("test.fixedbin", "test.sparse")
    for (k in c(10, 50)) {
        intervs <- gxmat_windows(500, align = k)
        ref <- gxmat_reference(exprs, intervs, iterator = k)
        a <- gextract.matrix(exprs, intervals = intervs, iterator = k)
        expect_equal(dim(a), c(nrow(intervs), 500 / k, 2))
        expect_identical(unname(a), ref)
        # one interval off the grid: all rows go through the interval-relative bins
        mixed <- rbind(intervs, data.frame(chrom = "chr1", start = 3, end = 503))
        expect_identical(unname(gextract.matrix(exprs, intervals = mixed, iterator = k))[seq_len(nrow(intervs)), , ], ref)
    }
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
    expect_error(gextract.matrix("test.fixedbin", intervals = intervs, colnames = c("a", "b")), "column names")
    expect_error(gextract.matrix("test.fixedbin", intervals = intervs[0, ]), "empty")
    expect_error(gextract.matrix("test.fixedbin", intervals = gintervals.2d(1, 0, 100, 1, 0, 100)), "1D")
    chrom_end <- gintervals.all()$end[gintervals.all()$chrom == "chr1"]
    expect_error(gextract.matrix("test.fixedbin", intervals = gintervals(1, chrom_end - 50, chrom_end + 50)), "chromosome boundaries")
    expect_error(gextract.matrix("test.fixedbin", intervals = gintervals(1, chrom_end - 53, chrom_end + 47), iterator = 10), "chromosome boundaries")
})

test_that("gextract.matrix looks up the R variables of an expression in the caller's frame", {
    intervs <- gintervals(1, c(1000, 3000), c(1200, 3200))
    expected <- gextract.matrix("test.fixedbin", intervals = intervs, iterator = 20)
    # a variable of the calling function was not found
    f <- function() {
        thr <- 2
        gextract.matrix("test.fixedbin * thr", intervals = intervs, iterator = 20)
    }
    expect_equal(f(), expected * 2)
    # a variable named as one of gextract.matrix's locals took the local's value
    width <- 0.5
    expect_equal(gextract.matrix("test.fixedbin * width", intervals = intervs, iterator = 20), expected * 0.5)
    iterator <- 1000
    expect_equal(gextract.matrix("test.fixedbin + iterator", intervals = intervs, iterator = 20), expected + 1000)
})

test_that("gextract.matrix rejects coordinates that are not whole numbers", {
    intervs <- data.frame(chrom = "chr1", start = 1000.5, end = 1100.5)
    expect_error(gextract.matrix("test.fixedbin", intervals = intervs, iterator = 10), "whole numbers")
})
