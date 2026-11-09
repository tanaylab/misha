create_isolated_test_db()

test_that("gsample with test.sparse", {
    set.seed(60427)
    expect_regression(gsample("test.sparse", 100, gintervals(c(1, 2))), "gsample.sparse")
})

test_that("gsample with test.rects", {
    set.seed(60427)
    expect_regression(gsample("test.rects", 100, gintervals.2d(c(1, 2))), "gsample.rects")
})

test_that("gsample draws from .Random.seed and advances it", {
    withr::local_preserve_seed()
    st <- seq(0, 99900, by = 100)
    gvtrack.create("gsample_pos", data.frame(chrom = "chr1", start = st, end = st + 100, value = st), "avg")
    withr::defer(gvtrack.rm("gsample_pos"))
    iv <- gintervals(1, 0, 100000)

    # R's next draw does not repeat the numbers gsample used
    set.seed(2)
    x1 <- runif(3)
    set.seed(2)
    invisible(gsample("gsample_pos", 50, iv, iterator = 100))
    x2 <- runif(3)
    expect_false(identical(x1, x2))

    # a restored .Random.seed gives the same sample again
    set.seed(3)
    saved <- .Random.seed
    s1 <- gsample("gsample_pos", 50, iv, iterator = 100)
    invisible(runif(1000))
    assign(".Random.seed", saved, envir = .GlobalEnv)
    s2 <- gsample("gsample_pos", 50, iv, iterator = 100)
    expect_identical(s1, s2)

    # set.seed() still reproduces a sample
    set.seed(4)
    s3 <- gsample("gsample_pos", 50, iv, iterator = 100)
    set.seed(4)
    expect_identical(gsample("gsample_pos", 50, iv, iterator = 100), s3)
})

test_that("forked kids draw different random numbers", {
    withr::local_preserve_seed()
    withr::local_options(gmultitasking = TRUE, gmax.processes = 4, gmin.scope4process = 1000)
    # ten points per 1 kb window, valued 0-9, on three chromosomes: a "sample" vtrack returns
    # the index of the point it drew
    chroms <- c("chr1", "chr2", "chr3")
    bins <- seq(0, 19000, by = 1000)
    src <- do.call(rbind, lapply(chroms, function(chrom) {
        data.frame(chrom = chrom, start = rep(bins, each = 10) + rep(0:9 * 100, length(bins)))
    }))
    src$end <- src$start + 1
    src$value <- rep(0:9, nrow(src) / 10)
    gvtrack.create("kid_sample", src, "sample")
    withr::defer(gvtrack.rm("kid_sample"))
    windows <- do.call(rbind, lapply(chroms, function(chrom) gintervals(chrom, bins, bins + 1000)))

    set.seed(5)
    r <- gextract("kid_sample", gintervals(chroms, 0, 20000), iterator = windows)
    per_chrom <- split(r$kid_sample, as.character(r$chrom))
    # each chromosome goes to a kid of its own; kids that drew the same numbers would sample
    # the same points on every chromosome
    expect_false(identical(per_chrom$chr1, per_chrom$chr2))
    expect_false(identical(per_chrom$chr2, per_chrom$chr3))

    set.seed(5)
    expect_identical(gextract("kid_sample", gintervals(chroms, 0, 20000), iterator = windows), r)
})
