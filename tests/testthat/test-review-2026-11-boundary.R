# Regression tests for the R/C++ boundary findings of the 2026-11 review of 5.14.4:
# argument conversion and validation, and process state left after an error.

create_isolated_test_db()

scope <- gintervals(1, 0, 20000)

# The bin glookup() puts each value in: (b1, b2], ...; force.binning puts values outside the
# breaks into the first or the last bin.
lookup_bin <- function(vals, breaks) {
    pmin(pmax(findInterval(vals, breaks, left.open = TRUE), 1), length(breaks) - 1)
}

test_that("an integer or logical NA is stored and looked up as NaN", {
    withr::defer({
        for (tr in c("temp.review_na_sparse", "temp.review_na_2d", "temp.review_na_lookup")) {
            if (gtrack.exists(tr)) gtrack.rm(tr, force = TRUE)
        }
    })

    iv <- gintervals(c(1, 1, 2), c(10, 500, 100), c(20, 600, 300))
    gtrack.create_sparse("temp.review_na_sparse", "", iv, c(1L, NA, 3L))
    expect_equal(gextract("temp.review_na_sparse", iv)$temp.review_na_sparse, c(1, NaN, 3))

    iv2 <- gintervals.2d(c(1, 2), 0, 1000, c(1, 2), 2000, 3000)
    gtrack.2d.create("temp.review_na_2d", "", iv2, c(5L, NA))
    expect_equal(gextract("temp.review_na_2d", iv2)$temp.review_na_2d, c(5, NaN))

    vals <- gextract("test.fixedbin", scope)$test.fixedbin
    breaks <- c(0, 0.1, 0.2, 1)
    bin <- lookup_bin(vals, breaks)
    expected <- c(1, NaN, 3)[bin]
    expect_equal(glookup(c(1L, NA, 3L), "test.fixedbin", breaks, intervals = scope)$value, expected)
    expect_equal(glookup(c(TRUE, NA, FALSE), "test.fixedbin", breaks, intervals = scope)$value, c(1, NaN, 0)[bin])

    gtrack.lookup("temp.review_na_lookup", "", c(1L, NA, 3L), "test.fixedbin", breaks)
    expect_equal(gextract("temp.review_na_lookup", scope)$temp.review_na_lookup, expected)
})

test_that("integer arguments are accepted wherever a number is expected", {
    vals <- gextract("test.fixedbin", scope)$test.fixedbin

    # breaks
    ivals <- round(vals * 10)
    expr <- "round(test.fixedbin * 10)"
    expect_equal(glookup(1:7, expr, 0:7, intervals = scope)$value, (1:7)[lookup_bin(ivals, 0:7)])
    # gpartition: runs of 50 bp bins that fall in the same bin of (0, 1], ..., (6, 7]
    bin <- findInterval(ivals, 0:7, left.open = TRUE)
    bin[bin < 1 | bin > 7] <- NA
    runs <- rle(bin)
    ends <- cumsum(runs$lengths) * 50
    keep <- !is.na(runs$values)
    part <- gpartition(expr, 0:7, scope)
    expect_equal(part$start, (ends - runs$lengths * 50)[keep])
    expect_equal(part$end, ends[keep])
    expect_equal(part$bin, runs$values[keep])

    # percentiles
    expect_equal(unname(gquantiles("test.fixedbin", 1L, scope)), max(vals, na.rm = TRUE), tolerance = 1e-6)
    expect_equal(unname(unlist(gintervals.quantiles("test.fixedbin", 1L, scope)[, "1"])), max(vals, na.rm = TRUE), tolerance = 1e-6)

    # distances
    a <- data.frame(chrom = "chr1", start = c(100, 5000), end = c(200, 5100))
    b <- data.frame(chrom = "chr1", start = c(1000, 5300), end = c(1100, 5400))
    expect_equal(gintervals.neighbors(a, b, 1L, -1000L, 1000L)$dist, c(800, 200))
    expect_equal(gintervals.neighbors(a, b, mindist = 0L, maxdist = 500L)$dist, 200)

    # a FUN that returns an integer, and an expression that does
    expect_equal(
        gintervals.mapply(function(x) sum(x > 0.1, na.rm = TRUE), "test.fixedbin", intervals = scope)$value,
        sum(vals > 0.1, na.rm = TRUE)
    )
    expect_equal(gextract("as.integer(test.fixedbin * 10)", scope)[[4]], as.numeric(as.integer(vals * 10)))
    expect_equal(gextract("ifelse(test.fixedbin > 0.1, 1L, NA_integer_)", scope)[[4]], ifelse(vals > 0.1, 1, NaN))

    # window sizes and thresholds: the same result as the double form
    big <- gintervals(1, 0, 1e6)
    expect_equal(gsegment("test.fixedbin", 10000L, 0.0001, intervals = big), gsegment("test.fixedbin", 10000, 0.0001, intervals = big))
    expect_equal(
        gwilcox("test.fixedbin", 100000L, 1000L, maxpval = 0.01, intervals = big),
        gwilcox("test.fixedbin", 100000, 1000, maxpval = 0.01, intervals = big)
    )
    domain <- gintervals(1, 800000 * (0:5), 800000 * (0:5) + 400000)
    src <- gextract("test.sparse", gintervals(1))
    expect_equal(gcis_decay("test.rects", (0:20) * 1000L, src, domain), gcis_decay("test.rects", (0:20) * 1000, src, domain))

    withr::defer({
        for (tr in c("temp.review_smooth_i", "temp.review_smooth_d")) {
            if (gtrack.exists(tr)) gtrack.rm(tr, force = TRUE)
        }
    })
    gtrack.smooth("temp.review_smooth_i", "", "test.fixedbin", 1000L, weight_thr = 0L, alg = "MEAN")
    gtrack.smooth("temp.review_smooth_d", "", "test.fixedbin", 1000, weight_thr = 0, alg = "MEAN")
    expect_equal(gextract("temp.review_smooth_i", scope)[[4]], gextract("temp.review_smooth_d", scope)[[4]])
})

test_that("options(gmax.data.size) of Inf means no limit, [0, 1) counts as 1, NA or negative is an error", {
    withr::local_options(gmax.data.size = Inf)
    expect_equal(unname(gquantiles("test.fixedbin", 1, scope)), max(gextract("test.fixedbin", scope)$test.fixedbin, na.rm = TRUE), tolerance = 1e-6)
    expect_equal(nrow(gextract("test.fixedbin", scope)), 400)
    withr::local_options(gmax.data.size = 2^64)
    expect_equal(nrow(gextract("test.fixedbin", scope)), 400)
    for (small in list(0, 0.5)) {
        withr::local_options(gmax.data.size = small)
        # a sampling buffer of one value, not a crash
        expect_length(gquantiles("test.fixedbin", c(0.1, 0.5), scope), 2)
        expect_error(gextract("test.fixedbin", scope), "maximum allowed \\(1 ")
    }
    for (bad in list(NA_real_, -1, -Inf)) {
        withr::local_options(gmax.data.size = bad)
        expect_error(gquantiles("test.fixedbin", c(0.1, 0.5), scope), "gmax.data.size")
    }
})

test_that("a numeric chrom column is read as chromosome names; other chrom types are rejected without leaving misha's state behind", {
    # 1 and 2 are "1" and "2", which the chromosome aliases resolve to chr1 and chr2, as in gintervals(1, ...)
    chr <- data.frame(chrom = c("chr1", "chr2"), start = 0, end = 1000)
    expected <- gextract("test.fixedbin", chr)
    for (num in list(1:2, c(1, 2))) {
        expect_equal(gextract("test.fixedbin", transform(chr, chrom = num)), expected)
        expect_identical(.glifetime_counters(), c(ref_count = 0L, protect_count = 0L))
    }
    chr2d <- gintervals.2d(c("chr1", "chr2"), 0, 1e7, c("chr1", "chr2"), 0, 1e7)
    expected2d <- gextract("test.rects", chr2d)
    expect_gt(nrow(expected2d), 0)
    num2d <- transform(chr2d, chrom1 = c(1L, 2L), chrom2 = c(1, 2))
    expect_equal(gextract("test.rects", num2d), expected2d)

    bad <- list(
        data.frame(chrom = I(list("chr1")), start = 0, end = 100),
        data.frame(chrom = NA_integer_, start = 0, end = 100),
        data.frame(chrom = 1.5, start = 0, end = 100),
        c("chr1", "chr2")
    )
    for (intervs in bad) {
        expect_error(gextract("test.fixedbin", intervs), "Invalid format of interval")
        expect_identical(.glifetime_counters(), c(ref_count = 0L, protect_count = 0L))
    }
    iv2 <- gintervals.2d(1, 0, 100, 1, 0, 100)
    iv2$chrom2 <- I(list("chr1"))
    expect_error(gextract("test.rects", iv2), "Invalid format of intervals")
    expect_identical(.glifetime_counters(), c(ref_count = 0L, protect_count = 0L))
})

test_that("a multitasked misha call inside a multitasked gintervals.mapply FUN runs serially", {
    withr::local_options(gmultitasking = TRUE, gmax.processes = 4, gmin.scope4process = 100)
    inner <- gintervals(c(1, 2), 0, 1e6)
    expected <- unname(quantile(gextract("test.fixedbin", inner, iterator = 50)$test.fixedbin, 0.5, na.rm = TRUE, type = 1))
    iv <- giterator.intervals(intervals = gintervals(1, 0, 360000), iterator = 10000)
    res <- gintervals.mapply(
        function(x) gquantiles("test.fixedbin", 0.5, inner, iterator = 50)[[1]],
        "test.fixedbin",
        intervals = iv, iterator = 50
    )
    expect_equal(nrow(res), 36)
    expect_equal(res$value, rep(expected, 36), tolerance = 1e-6)
})

test_that("intervals given as a data.table or a grouped tibble are accepted", {
    skip_if_not_installed("data.table")
    a <- data.frame(chrom = c("chr1", "chr1", "chr2"), start = c(100, 5000, 300), end = c(200, 5100, 400), id = c("a", "b", "c"))
    b <- data.frame(chrom = c("chr1", "chr2"), start = c(1000, 10), end = c(1100, 20), nm = c("x", "y"))
    convs <- list(
        data.table = data.table::as.data.table,
        grouped_df = function(x) dplyr::group_by(tibble::as_tibble(x), chrom)
    )
    for (conv in convs) {
        res <- gintervals.neighbors(conv(a), conv(b))
        expect_equal(res$dist, c(800, 3900, 280))
        expect_equal(res$nm, c("x", "x", "y"))
        ann <- gintervals.annotate(conv(a), conv(b))
        expect_equal(ann$nm, c("x", "x", "y"))
        rb <- gintervals.rbind(conv(a), conv(a))
        expect_identical(class(rb), "data.frame")
        expect_equal(rb$id, rep(c("a", "b", "c"), 2))
        expect_null(attr(rb, "groups"))
    }
})

test_that("invalid scalar arguments and options are rejected with a clear error", {
    a <- data.frame(chrom = "chr1", start = c(100, 5000), end = c(200, 5100))
    b <- data.frame(chrom = "chr1", start = c(1000, 5300), end = c(1100, 5400))

    # an integer NA strand is not strand 0
    expect_error(gseq.extract(data.frame(chrom = "chr1", start = 0, end = 8, strand = c(1L, NA, -1L))), "strand")

    # NA where a TRUE/FALSE is needed
    withr::with_options(list(gmultitasking = NA), expect_error(gsummary("test.fixedbin", scope), "gmultitasking"))
    expect_error(gpartition("test.fixedbin", c(0, 0.1, 1), scope, include.lowest = NA), "include.lowest")
    expect_error(glookup(1:2, "test.fixedbin", c(0, 0.1, 1), intervals = scope, include.lowest = NA), "include.lowest")
    expect_error(gdist("test.fixedbin", c(0, 0.1, 1), intervals = scope, include.lowest = NA), "include.lowest")

    # NA, NaN or a missing value where a number is needed
    expect_error(gquantiles("test.fixedbin", NA_real_, scope), "range")
    expect_error(gintervals.neighbors(a, b, maxdist = numeric(0)), "maxdist")
    expect_error(gintervals.neighbors(a, b, maxdist = NA_real_), "maxdist")
    expect_error(gextract("test.fixedbin", scope, iterator = numeric(0)), "iterator")
    gvtrack.create("review_vt_shift", "test.fixedbin", "avg")
    for (shift in list(NaN, Inf, NA_real_)) {
        expect_error(gvtrack.iterator("review_vt_shift", sshift = shift), "sshift")
    }
    gvtrack.rm("review_vt_shift")
    expect_error(
        {
            gvtrack.create("review_vt_q", "test.fixedbin", "quantile", params = NA_real_)
            gextract("review_vt_q", scope)
        },
        "percentile"
    )

    # an infinite distance limit is no limit
    expect_equal(gintervals.neighbors(a, b, maxdist = Inf)$dist, c(800, 200))
    expect_equal(gintervals.neighbors(a, b, mindist = -Inf, maxdist = Inf)$dist, c(800, 200))
})

test_that("gintervals.load_chain accepts an integer min_score", {
    local_db_state()
    setup_db(list(">chr1\nAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n"))
    chain_file <- new_chain_file()
    cat("chain 5000 source1 100 + 10 30 chr1 60 + 10 30 1\n20\n\n", file = chain_file)
    cat("chain 3000 source1 100 + 35 55 chr1 60 + 40 60 2\n20\n\n", file = chain_file, append = TRUE)
    chain <- gintervals.load_chain(chain_file, min_score = 4000L)
    expect_equal(chain$score, 5000)
})

test_that("logical flags are TRUE or FALSE, not NA", {
    a <- data.frame(chrom = "chr1", start = c(100, 5000), end = c(200, 5100))
    b <- data.frame(chrom = "chr1", start = c(1000, 5300), end = c(1100, 5400))
    # (gsegment and gwilcox already stop on onetailed = NA in R, before the C++ check)
    expect_error(glookup(1:2, "test.fixedbin", c(0, 0.1, 1), intervals = scope, force.binning = NA), "force.binning")
    expect_error(gintervals.neighbors(a, b, na.if.notfound = NA), "na.if.notfound")
    expect_error(gintervals.neighbors(a, b, use_intervals1_strand = NA), "use_intervals1_strand")
    withr::defer(if (gtrack.exists("temp.review_smooth_na")) gtrack.rm("temp.review_smooth_na", force = TRUE))
    expect_error(gtrack.smooth("temp.review_smooth_na", "", "test.fixedbin", 1000, smooth_nans = NA), "nans")
})

test_that("computing a track's percentiles for a global.percentile vtrack keeps options(gmultitasking)", {
    withr::defer(if (gtrack.exists("temp.review_pv")) gtrack.rm("temp.review_pv", force = TRUE))
    gtrack.create_dense("temp.review_pv", "", gintervals(1, 0, 1e5), 1, 10000, 0)
    gvtrack.create("review_pv", "temp.review_pv", "global.percentile")
    withr::defer(gvtrack.rm("review_pv"))
    withr::local_options(gmultitasking = TRUE)
    gextract("review_pv", gintervals(1, 0, 1e5))
    expect_true(getOption("gmultitasking"))
})

test_that("a chain data frame with a numeric chromsrc is read; character strands are rejected without leaving state behind", {
    chain <- data.frame(
        chrom = "chr1", start = 1000, end = 2000, strand = 0, chromsrc = "1", startsrc = 0, endsrc = 1000,
        strandsrc = 0, chain_id = 1L, score = 1
    )
    src <- data.frame(chrom = "1", start = 100, end = 200)
    expected <- gintervals.liftover(src, chain)
    expect_equal(expected$start, 1100)
    for (num in list(1L, 1)) {
        expect_equal(gintervals.liftover(src, transform(chain, chromsrc = num)), expected)
        expect_identical(.glifetime_counters(), c(ref_count = 0L, protect_count = 0L))
    }
    expect_error(gintervals.liftover(src, transform(chain, strand = "+", strandsrc = "+")), "strand")
    expect_identical(.glifetime_counters(), c(ref_count = 0L, protect_count = 0L))
    expect_error(gintervals.liftover(src, transform(chain, chromsrc = NA_integer_)), "Invalid format")
    expect_identical(.glifetime_counters(), c(ref_count = 0L, protect_count = 0L))
})

test_that("options(gmax.data.size = Inf) is no limit for gbins.quantiles, gtrack.2d.import and a cartesian grid iterator", {
    vals <- gextract("test.fixedbin", scope)$test.fixedbin
    breaks <- c(0, 0.1, 0.2, 1)
    expected <- tapply(vals, cut(vals, breaks), max)
    withr::local_options(gmax.data.size = Inf)
    res <- gbins.quantiles("test.fixedbin", breaks, expr = "test.fixedbin", percentiles = 1, intervals = scope)
    expect_equal(unname(res[, 1]), unname(as.numeric(expected)), tolerance = 1e-6)

    withr::defer(if (gtrack.exists("temp.review_2d_inf")) gtrack.rm("temp.review_2d_inf", force = TRUE))
    set.seed(1)
    n <- 5000
    pts <- data.frame(chrom1 = "chr1", start1 = sample(1e6, n) - 1, chrom2 = "chr2", start2 = sample(1e6, n) - 1, value = seq_len(n))
    pts$end1 <- pts$start1 + 1
    pts$end2 <- pts$start2 + 1
    f <- tempfile(fileext = ".tsv")
    write.table(pts[, c("chrom1", "start1", "end1", "chrom2", "start2", "end2", "value")], f, sep = "\t", quote = FALSE, row.names = FALSE)
    gtrack.2d.import("temp.review_2d_inf", "", f)
    got <- gextract("temp.review_2d_inf", gintervals.2d(1, 0, 1e6, 2, 0, 1e6))
    expect_equal(sort(got$temp.review_2d_inf), as.numeric(seq_len(n)))

    # a limit above 2^31 used to be truncated to int by the 2D big set and cartesian grid iterators
    withr::local_options(gmax.data.size = 2^32)
    intervs <- gintervals(1, c(1000, 50000), c(2000, 51000))
    itr <- giterator.cartesian_grid(intervs, c(-100, 100))
    expect_equal(nrow(giterator.intervals(intervals = gintervals.2d(1), iterator = itr)), 4)
})
