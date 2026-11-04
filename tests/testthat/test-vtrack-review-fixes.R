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
