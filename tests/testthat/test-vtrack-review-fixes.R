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
