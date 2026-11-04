# Regression tests for virtual track defects found in the November 2026 review of 5.14.4.
# Expected values are computed by hand or with base R, not with another misha code path.

vt_fix_db <- function() {
    set.seed(1)
    s1 <- paste(sample(c("A", "C", "G", "T"), 5000, TRUE), collapse = "")
    s2 <- paste(sample(c("A", "C", "G", "T"), 300, TRUE), collapse = "")
    setup_db(list(sprintf(">chr1\n%s\n", s1), sprintf(">chr2\n%s\n", s2)), return_db = TRUE)
}

vt_fix_pssm <- function() {
    set.seed(2)
    pssm <- matrix(runif(24), 6, 4)
    pssm <- pssm / rowSums(pssm)
    colnames(pssm) <- c("A", "C", "G", "T")
    pssm
}

test_that("a NaN bin of a dense track gives no value to exists, size and the position functions", {
    local_db_state()
    vt_fix_db()
    gtrack.create_dense("dn", "x", gintervals("chr1", c(0, 50, 100), c(50, 100, 150)), c(1, NaN, 3), 50, NaN)
    withr::defer(gtrack.rm("dn", force = TRUE))
    gvtrack.create("ex", "dn", "exists")
    gvtrack.create("sz", "dn", "size")
    gvtrack.create("mp", "dn", "max.pos.abs")
    gvtrack.create("fp", "dn", "first.pos.relative")
    withr::defer(gvtrack.clear())

    # iterator equal to the bin size, all functions in one call
    r <- gextract(c("dn", "ex", "sz", "mp", "fp"), gintervals("chr1", 0, 150), iterator = 50)
    expect_equal(r$ex, c(1, 0, 1))
    expect_equal(r$sz, c(1, 0, 1))
    expect_equal(r$mp, c(0, NaN, 100))
    expect_equal(r$fp, c(0, NaN, 0))

    # a position function alone, and an iterator smaller than the bin
    r <- gextract("mp", gintervals("chr1", 50, 100), iterator = 50)
    expect_true(is.na(r$mp))
    r <- gextract(c("ex", "sz", "mp"), gintervals("chr1", 40, 110), iterator = 10)
    in_nan_bin <- r$start >= 50 & r$start < 100
    expect_true(all(r$ex[in_nan_bin] == 0))
    expect_true(all(r$sz[in_nan_bin] == 0))
    expect_true(all(is.na(r$mp[in_nan_bin])))
    expect_equal(r$mp[!in_nan_bin], r$start[!in_nan_bin])
})

test_that("value-based vtracks read only the query's chromosome", {
    local_db_state()
    vt_fix_db()
    vals <- data.frame(chrom = c("chr1", "chr2"), start = c(1000, 10), end = c(1100, 20), score = c(5, 7))
    gvtrack.create("v_avg", vals, "avg")
    gvtrack.create("v_near", vals, "nearest")
    gvtrack.create("v_ex", vals, "exists")
    withr::defer(gvtrack.clear())
    q <- gintervals(c("chr1", "chr1", "chr2"), c(1000, 3000, 10), c(1100, 3100, 20))
    r <- gextract(c("v_avg", "v_near", "v_ex"), q, iterator = q)
    r <- r[order(r$chrom, r$start), ]
    expect_equal(r$v_avg, c(5, NaN, 7))
    expect_equal(r$v_near, c(5, 5, 7))
    expect_equal(r$v_ex, c(1, 0, 1))
})

test_that("a filter that splits a bin or an interval does not count its value twice", {
    local_db_state()
    vt_fix_db()
    gtrack.create_dense("dn", "x", gintervals("chr1", c(0, 50, 100), c(50, 100, 150)), c(1, NaN, 3), 50, NaN)
    gtrack.create_sparse("sp", "x", gintervals("chr1", c(313, 400), c(332, 420)), c(-1.5, 2))
    withr::defer({
        gtrack.rm("dn", force = TRUE)
        gtrack.rm("sp", force = TRUE)
    })
    mask <- gintervals("chr1", c(20, 320), c(30, 322))
    for (f in c("sum", "size", "lse")) {
        gvtrack.create(paste0("d_", f), "dn", f, filter = mask)
        gvtrack.create(paste0("s_", f), "sp", f, filter = mask)
    }
    gvtrack.create("d_q", "dn", "quantile", 0.5, filter = mask)
    withr::defer(gvtrack.clear())

    # bin [0,50) is split into [0,20) and [30,50) by the mask; its value 1 counts once
    iv <- gintervals("chr1", 0, 150)
    r <- gextract(c("d_sum", "d_size", "d_lse", "d_q"), iv, iterator = iv)
    expect_equal(r$d_sum, 4)
    expect_equal(r$d_size, 2)
    expect_equal(r$d_lse, log(exp(1) + exp(3)), tolerance = 1e-6)
    expect_equal(r$d_q, median(c(1, 3)))

    # sparse interval [313,332) is split by [320,322); [400,420) is the second value
    iv <- gintervals("chr1", 300, 450)
    r <- gextract(c("s_sum", "s_size", "s_lse"), iv, iterator = iv)
    expect_equal(r$s_sum, 0.5)
    expect_equal(r$s_size, 2)
    expect_equal(r$s_lse, log(exp(-1.5) + exp(2)), tolerance = 1e-6)

    # the same for a value-based vtrack
    vv <- data.frame(chrom = "chr1", start = c(313, 400), end = c(332, 420), score = c(-1.5, 2))
    gvtrack.create("v_sum", vv, "sum", filter = mask)
    gvtrack.create("v_size", vv, "size", filter = mask)
    gvtrack.create("v_avg", vv, "avg", filter = mask)
    withr::defer(gvtrack.clear())
    r <- gextract(c("v_sum", "v_size", "v_avg"), iv, iterator = iv)
    expect_equal(r$v_sum, 0.5)
    expect_equal(r$v_size, 2)
    expect_equal(r$v_avg, 0.25)
})

test_that("a value-based quantile with a filter uses every unmasked part", {
    local_db_state()
    vt_fix_db()
    vv <- data.frame(chrom = "chr1", start = seq(100, 600, 100), end = seq(150, 650, 100), score = c(1, 9, 2, 8, 3, 7))
    gvtrack.create("v_q", vv, "quantile", 0.5, filter = gintervals("chr1", 360, 380))
    withr::defer(gvtrack.clear())
    iv <- gintervals("chr1", 0, 1000)
    expect_equal(gextract("v_q", iv, iterator = iv)$v_q, median(vv$score))
})

test_that("edit distance positions are offset when a filter leaves one part that does not start the interval", {
    local_db_state()
    vt_fix_db()
    pssm <- vt_fix_pssm()
    gvtrack.create("pm", NULL, "pwm.max", pssm = pssm)
    gvtrack.create("p", NULL, "pwm", pssm = pssm)
    withr::defer(gvtrack.clear())
    part <- gintervals("chr1", 110, 140)
    whole <- gintervals("chr1", 100, 140)
    scores <- gextract(c("pm", "p"), part, iterator = part)
    # thresholds half a unit above the current score, so that at least one edit is needed
    thresh <- c(pwm.edit_distance.pos = scores$pm + 0.5, pwm.edit_distance.lse.pos = scores$p + 0.5)

    for (f in names(thresh)) {
        gvtrack.create("ed", NULL, f, pssm = pssm, score.thresh = thresh[[f]])
        unfiltered <- gextract("ed", part, iterator = part)$ed
        expect_false(is.na(unfiltered), info = f)
        gvtrack.filter("ed", gintervals("chr1", 100, 110))
        filtered <- gextract("ed", whole, iterator = whole)$ed
        expect_equal(abs(filtered), abs(unfiltered) + 10, info = f)
        expect_equal(sign(filtered), sign(unfiltered), info = f)
        gvtrack.rm("ed")
    }
})

test_that("gtrack.array.extract reads columns given in any order", {
    asc <- gtrack.array.extract("test.array", c("col1", "col5"), gintervals(c(1, 2)))
    desc <- gtrack.array.extract("test.array", c("col5", "col1"), gintervals(c(1, 2)))
    expect_equal(desc$col1, asc$col1)
    expect_equal(desc$col5, asc$col5)
    expect_true(any(!is.na(desc$col1)))
})

test_that("min.pos on an array track reports the position of the minimum", {
    local_db_state()
    vt_fix_db()
    f <- withr::local_tempfile(fileext = ".tsv")
    writeLines(c(
        "chrom\tstart\tend\ta\tb",
        "chr1\t0\t10\t5\t5",
        "chr1\t20\t30\t1\t1",
        "chr1\t40\t50\t3\t3"
    ), f)
    suppressMessages(gtrack.array.import("arr", "x", f))
    withr::defer(gtrack.rm("arr", force = TRUE))
    gvtrack.create("a_minp", "arr", "min.pos.abs")
    gvtrack.create("a_minr", "arr", "min.pos.relative")
    gvtrack.create("a_maxp", "arr", "max.pos.abs")
    withr::defer(gvtrack.clear())
    iv <- gintervals("chr1", 0, 100)
    r <- gextract(c("a_minp", "a_minr", "a_maxp"), iv, iterator = iv)
    expect_equal(r$a_minp, 20)
    expect_equal(r$a_minr, 20)
    expect_equal(r$a_maxp, 0)
})

test_that("stddev with a filter is the sample sd of the unmasked values, each counted once", {
    local_db_state()
    vt_fix_db()
    starts <- seq(0, 950, 50)
    gtrack.create_dense("dn2", "x", gintervals("chr1", starts, starts + 50), as.numeric(1:20), 50, NaN)
    gtrack.create_sparse("sp2", "x", gintervals("chr1", c(10, 100, 200, 300), c(60, 120, 220, 320)), c(4, 1, 7, 2))
    withr::defer({
        gtrack.rm("dn2", force = TRUE)
        gtrack.rm("sp2", force = TRUE)
    })
    withr::defer(gvtrack.clear())
    # [20,30) splits bin 1 (and the first sparse interval); [500,600) removes bins 11 and 12; [190,230) removes 7
    mask <- gintervals("chr1", c(20, 190, 500), c(30, 230, 600))
    gvtrack.create("d_sd", "dn2", "stddev", filter = mask)
    gvtrack.create("d_avg", "dn2", "avg", filter = mask)
    gvtrack.create("s_sd", "sp2", "stddev", filter = mask)
    vv <- data.frame(chrom = "chr1", start = c(10, 100, 200, 300), end = c(60, 120, 220, 320), score = c(4, 1, 7, 2))
    gvtrack.create("v_sd", vv, "stddev", filter = mask)
    iv <- gintervals("chr1", 0, 1000)

    expected_dense <- sd(c(1:10, 13:20))
    expect_equal(gextract("d_sd", iv, iterator = iv)$d_sd, expected_dense, tolerance = 1e-5)
    expect_equal(gextract(c("d_sd", "d_avg"), iv, iterator = iv)$d_sd, expected_dense, tolerance = 1e-5)
    expect_equal(gextract("s_sd", iv, iterator = iv)$s_sd, sd(c(4, 1, 2)), tolerance = 1e-5)
    expect_equal(gextract("v_sd", iv, iterator = iv)$v_sd, sd(c(4, 1, 2)), tolerance = 1e-5)
})

test_that("positions of sparse, value-based and array values are clipped to the iterator start", {
    local_db_state()
    vt_fix_db()
    gtrack.create_sparse("sp3", "x", gintervals("chr1", c(313, 400), c(332, 420)), c(5, 1))
    f <- withr::local_tempfile(fileext = ".tsv")
    writeLines(c("chrom\tstart\tend\ta", "chr1\t313\t332\t5", "chr1\t400\t420\t1"), f)
    suppressMessages(gtrack.array.import("arr3", "x", f))
    withr::defer({
        gtrack.rm("sp3", force = TRUE)
        gtrack.rm("arr3", force = TRUE)
    })
    withr::defer(gvtrack.clear())
    vv <- data.frame(chrom = "chr1", start = c(313, 400), end = c(332, 420), score = c(5, 1))
    iv <- gintervals("chr1", 318, 450)
    one <- gintervals("chr1", 318, 330)

    for (src in c("sp3", "vv")) {
        source_obj <- if (src == "vv") vv else src
        for (fn in c("first", "last", "min", "max", "sample")) {
            gvtrack.create("pa", source_obj, paste0(fn, ".pos.abs"))
            gvtrack.create("pr", source_obj, paste0(fn, ".pos.relative"))
            r <- gextract(c("pa", "pr"), iv, iterator = iv)
            expected <- c(first = 318, last = 400, min = 400, max = 318)
            if (fn != "sample") {
                expect_equal(r$pa, expected[[fn]], info = paste(src, fn))
                expect_equal(r$pr, expected[[fn]] - 318, info = paste(src, fn))
            }
            r <- gextract(c("pa", "pr"), one, iterator = one)
            expect_equal(r$pa, 318, info = paste(src, fn, "one value"))
            expect_equal(r$pr, 0, info = paste(src, fn, "one value"))
        }
    }
    gvtrack.create("a_max", "arr3", "max.pos.relative")
    gvtrack.create("a_min", "arr3", "min.pos.relative")
    r <- gextract(c("a_max", "a_min"), iv, iterator = iv)
    expect_equal(r$a_max, 0)
    expect_equal(r$a_min, 400 - 318)
})
