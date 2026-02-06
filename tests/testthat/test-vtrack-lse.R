create_isolated_test_db()

# Reference R implementation of log-sum-exp (numerically stable)
r_lse <- function(x) {
    x <- x[!is.na(x)]
    if (length(x) == 0) {
        return(NA_real_)
    }
    m <- max(x)
    m + log(sum(exp(x - m)))
}

# Helper: compute LSE manually from raw track values over an interval
manual_lse <- function(track_name, interval_row) {
    interval <- gintervals(interval_row$chrom, interval_row$start, interval_row$end)
    values <- gextract(track_name, interval, iterator = track_name)
    track_vals <- values[[track_name]]
    r_lse(track_vals)
}

manual_lse_vec <- function(track_name, intervals) {
    vapply(seq_len(nrow(intervals)), function(i) {
        manual_lse(track_name, intervals[i, ])
    }, numeric(1))
}

# ============================================================
# Basic correctness
# ============================================================

test_that("lse vtrack function works on sparse track", {
    gvtrack.create("test.lse.vt", src = "test.sparse", func = "lse")
    gvtrack.create("test.sum.vt", src = "test.sparse", func = "sum")
    on.exit(
        {
            gvtrack.rm("test.lse.vt")
            gvtrack.rm("test.sum.vt")
        },
        add = TRUE
    )

    intervals <- gintervals.all()
    result <- gextract("test.lse.vt", intervals = intervals, iterator = intervals, colnames = "lse_value")
    result_sum <- gextract("test.sum.vt", intervals = intervals, iterator = intervals, colnames = "sum_value")

    expect_equal(nrow(result), nrow(result_sum))
    # LSE is NA iff sum is NA (both require at least one value)
    expect_equal(is.na(result$lse_value), is.na(result_sum$sum_value))
})

test_that("lse vtrack function correctness with known values", {
    intervals_df <- data.frame(
        chrom = "chr1",
        start = c(100, 200, 300),
        end = c(200, 300, 400),
        score = c(1.0, 2.0, 3.0)
    )

    gvtrack.create("test.lse.known", src = intervals_df, func = "lse")
    on.exit(gvtrack.rm("test.lse.known"), add = TRUE)

    iter_int <- gintervals("chr1", 100, 400)
    result <- gextract("test.lse.known", intervals = iter_int, iterator = iter_int, colnames = "value")

    expected <- log(exp(1.0) + exp(2.0) + exp(3.0))
    expect_equal(result$value, expected, tolerance = 1e-5)
})

test_that("lse of single value equals the value itself", {
    intervals_df <- data.frame(
        chrom = "chr1",
        start = 100,
        end = 200,
        score = 5.0
    )

    gvtrack.create("test.lse.single", src = intervals_df, func = "lse")
    on.exit(gvtrack.rm("test.lse.single"), add = TRUE)

    iter_int <- gintervals("chr1", 100, 200)
    result <- gextract("test.lse.single", intervals = iter_int, iterator = iter_int, colnames = "value")

    expect_equal(result$value, 5.0, tolerance = 1e-6)
})

test_that("lse returns NA for intervals with no data", {
    intervals_df <- data.frame(
        chrom = "chr1",
        start = 100,
        end = 200,
        score = 1.0
    )

    gvtrack.create("test.lse.empty", src = intervals_df, func = "lse")
    on.exit(gvtrack.rm("test.lse.empty"), add = TRUE)

    iter_int <- gintervals("chr1", 500, 600)
    result <- gextract("test.lse.empty", intervals = iter_int, iterator = iter_int, colnames = "value")

    expect_true(is.na(result$value))
})

test_that("lse with negative values (log-space)", {
    intervals_df <- data.frame(
        chrom = "chr1",
        start = c(100, 200),
        end = c(200, 300),
        score = c(-1.0, -2.0)
    )

    gvtrack.create("test.lse.neg", src = intervals_df, func = "lse")
    on.exit(gvtrack.rm("test.lse.neg"), add = TRUE)

    iter_int <- gintervals("chr1", 100, 300)
    result <- gextract("test.lse.neg", intervals = iter_int, iterator = iter_int, colnames = "value")

    expected <- log(exp(-1.0) + exp(-2.0))
    expect_equal(result$value, expected, tolerance = 1e-5)
})

# ============================================================
# Dense track
# ============================================================

test_that("lse works on dense track", {
    gvtrack.create("test.lse.dense", src = "test.fixedbin", func = "lse")
    on.exit(gvtrack.rm("test.lse.dense"), add = TRUE)

    intervals <- gintervals("chr1", 0, 1000)
    result <- gextract("test.lse.dense", intervals = intervals, iterator = intervals, colnames = "value")

    raw <- gextract("test.fixedbin", intervals = intervals, iterator = "test.fixedbin")
    vals <- raw$test.fixedbin[!is.na(raw$test.fixedbin)]
    if (length(vals) > 0) {
        expected <- r_lse(vals)
        expect_equal(result$value, expected, tolerance = 1e-3)
    }
})

test_that("lse on dense track matches manual computation across multiple intervals", {
    gvtrack.create("vt_lse_dense_multi", src = "test.fixedbin", func = "lse")
    on.exit(gvtrack.rm("vt_lse_dense_multi"), add = TRUE)

    intervals <- rbind(
        gintervals(1, 0, 500),
        gintervals(1, 500, 1000),
        gintervals(1, 1000, 2000)
    )

    result <- gextract("vt_lse_dense_multi", intervals = intervals, iterator = intervals, colnames = "value")
    manual <- manual_lse_vec("test.fixedbin", intervals)

    expect_equal(result$value, manual, tolerance = 1e-3)
})

test_that("lse on dense track with numeric iterator", {
    gvtrack.create("vt_lse_dense_iter", src = "test.fixedbin", func = "lse")
    on.exit(gvtrack.rm("vt_lse_dense_iter"), add = TRUE)

    result <- gextract("vt_lse_dense_iter", gintervals(1, 0, 500), iterator = 100, colnames = "value")

    # Each 100bp window should have its own LSE
    expect_equal(nrow(result), 5)
    # Verify first window manually
    raw <- gextract("test.fixedbin", gintervals(1, 0, 100), iterator = "test.fixedbin")
    vals <- raw$test.fixedbin[!is.na(raw$test.fixedbin)]
    if (length(vals) > 0) {
        expect_equal(result$value[1], r_lse(vals), tolerance = 1e-3)
    }
})

# ============================================================
# Sparse track: manual verification
# ============================================================

test_that("lse on sparse track matches manual computation", {
    gvtrack.create("vt_lse_sparse_verify", src = "test.sparse", func = "lse")
    on.exit(gvtrack.rm("vt_lse_sparse_verify"), add = TRUE)

    intervals <- rbind(
        gintervals(1, 0, 300),
        gintervals(1, 600, 1000),
        gintervals(1, 1200, 1500)
    )

    result <- gextract("vt_lse_sparse_verify", intervals = intervals, iterator = intervals, colnames = "value")
    manual <- manual_lse_vec("test.sparse", intervals)

    # Both should have same NA pattern
    expect_equal(is.na(result$value), is.na(manual))
    # Where both are non-NA, values should match
    both_valid <- !is.na(result$value) & !is.na(manual)
    if (any(both_valid)) {
        expect_equal(result$value[both_valid], manual[both_valid], tolerance = 1e-3)
    }
})

test_that("lse on sparse track with large iterator", {
    gvtrack.create("vt_lse_sparse_large", src = "test.sparse", func = "lse")
    on.exit(gvtrack.rm("vt_lse_sparse_large"), add = TRUE)

    # Large iterator covering entire chromosome
    result <- gextract("vt_lse_sparse_large", gintervals(c(1, 2)), iterator = 10000, colnames = "value")

    expect_true(nrow(result) > 0)
    # At least some values should be non-NA
    expect_true(any(!is.na(result$value)))
})

# ============================================================
# Iterator shifts
# ============================================================

test_that("lse with iterator shifts works", {
    intervals_df <- data.frame(
        chrom = "chr1",
        start = c(100, 200, 300),
        end = c(200, 300, 400),
        score = c(1.0, 2.0, 3.0)
    )

    gvtrack.create("test.lse.shift", src = intervals_df, func = "lse")
    gvtrack.iterator("test.lse.shift", sshift = -100, eshift = 100)
    on.exit(gvtrack.rm("test.lse.shift"), add = TRUE)

    # Iterator interval [200, 300] with shifts becomes [100, 400], covering all 3 values
    iter_int <- gintervals("chr1", 200, 300)
    result <- gextract("test.lse.shift", intervals = iter_int, iterator = iter_int, colnames = "value")

    expected <- log(exp(1.0) + exp(2.0) + exp(3.0))
    expect_equal(result$value, expected, tolerance = 1e-5)
})

test_that("lse with sshift only narrows covered data", {
    intervals_df <- data.frame(
        chrom = "chr1",
        start = c(100, 200, 300),
        end = c(200, 300, 400),
        score = c(1.0, 2.0, 3.0)
    )

    gvtrack.create("vt_lse_sshift", src = intervals_df, func = "lse")
    gvtrack.iterator("vt_lse_sshift", sshift = 100)
    on.exit(gvtrack.rm("vt_lse_sshift"), add = TRUE)

    # Iterator [100, 400] with sshift=100 becomes [200, 400], covering values 2.0 and 3.0
    iter_int <- gintervals("chr1", 100, 400)
    result <- gextract("vt_lse_sshift", intervals = iter_int, iterator = iter_int, colnames = "value")

    expected <- log(exp(2.0) + exp(3.0))
    expect_equal(result$value, expected, tolerance = 1e-5)
})

test_that("lse with eshift only extends covered data", {
    intervals_df <- data.frame(
        chrom = "chr1",
        start = c(100, 200, 300),
        end = c(200, 300, 400),
        score = c(1.0, 2.0, 3.0)
    )

    gvtrack.create("vt_lse_eshift", src = intervals_df, func = "lse")
    gvtrack.iterator("vt_lse_eshift", eshift = 200)
    on.exit(gvtrack.rm("vt_lse_eshift"), add = TRUE)

    # Iterator [100, 200] with eshift=200 becomes [100, 400], covering all 3
    iter_int <- gintervals("chr1", 100, 200)
    result <- gextract("vt_lse_eshift", intervals = iter_int, iterator = iter_int, colnames = "value")

    expected <- log(exp(1.0) + exp(2.0) + exp(3.0))
    expect_equal(result$value, expected, tolerance = 1e-5)
})

test_that("lse with shifts on dense track", {
    gvtrack.create("vt_lse_dense_shift", src = "test.fixedbin", func = "lse")
    gvtrack.iterator("vt_lse_dense_shift", sshift = -50, eshift = 50)
    on.exit(gvtrack.rm("vt_lse_dense_shift"), add = TRUE)

    iter_int <- gintervals(1, 200, 300)
    result <- gextract("vt_lse_dense_shift", intervals = iter_int, iterator = iter_int, colnames = "value")

    # Shifted interval is [150, 350]
    raw <- gextract("test.fixedbin", gintervals(1, 150, 350), iterator = "test.fixedbin")
    vals <- raw$test.fixedbin[!is.na(raw$test.fixedbin)]
    if (length(vals) > 0) {
        expect_equal(result$value, r_lse(vals), tolerance = 1e-3)
    }
})

# ============================================================
# Value-based tracks
# ============================================================

test_that("lse on value-based track works", {
    intervals_df <- data.frame(
        chrom = "chr1",
        start = c(0, 100, 200),
        end = c(100, 200, 300),
        score = c(0.5, 1.5, 2.5)
    )

    gvtrack.create("test.lse.val", src = intervals_df, func = "lse")
    on.exit(gvtrack.rm("test.lse.val"), add = TRUE)

    iter_int <- gintervals("chr1", 0, 300)
    result <- gextract("test.lse.val", intervals = iter_int, iterator = iter_int, colnames = "value")

    expected <- log(exp(0.5) + exp(1.5) + exp(2.5))
    expect_equal(result$value, expected, tolerance = 1e-5)
})

test_that("lse on value-based track across multiple chromosomes", {
    intervals_df <- data.frame(
        chrom = c("chr1", "chr1", "chr2", "chr2"),
        start = c(100, 300, 100, 300),
        end = c(200, 400, 200, 400),
        score = c(1.0, 2.0, 3.0, 4.0)
    )

    gvtrack.create("vt_lse_multichrom", src = intervals_df, func = "lse")
    on.exit(gvtrack.rm("vt_lse_multichrom"), add = TRUE)

    # chr1: lse(1, 2)
    iter_chr1 <- gintervals("chr1", 100, 400)
    res1 <- gextract("vt_lse_multichrom", intervals = iter_chr1, iterator = iter_chr1, colnames = "value")
    expect_equal(res1$value, log(exp(1) + exp(2)), tolerance = 1e-5)

    # chr2: lse(3, 4)
    iter_chr2 <- gintervals("chr2", 100, 400)
    res2 <- gextract("vt_lse_multichrom", intervals = iter_chr2, iterator = iter_chr2, colnames = "value")
    expect_equal(res2$value, log(exp(3) + exp(4)), tolerance = 1e-5)
})

test_that("lse on value-based track with partial overlap", {
    intervals_df <- data.frame(
        chrom = "chr1",
        start = c(100, 300, 500),
        end = c(200, 400, 600),
        score = c(1.0, 2.0, 3.0)
    )

    gvtrack.create("vt_lse_partial", src = intervals_df, func = "lse")
    on.exit(gvtrack.rm("vt_lse_partial"), add = TRUE)

    # Query only overlaps first two data intervals
    iter_int <- gintervals("chr1", 150, 350)
    result <- gextract("vt_lse_partial", intervals = iter_int, iterator = iter_int, colnames = "value")
    expect_equal(result$value, log(exp(1.0) + exp(2.0)), tolerance = 1e-5)
})

test_that("lse on value-based track handles NA values", {
    intervals_df <- data.frame(
        chrom = "chr1",
        start = c(100, 200, 300),
        end = c(200, 300, 400),
        score = c(1.0, NA, 3.0)
    )

    gvtrack.create("vt_lse_na", src = intervals_df, func = "lse")
    on.exit(gvtrack.rm("vt_lse_na"), add = TRUE)

    iter_int <- gintervals("chr1", 100, 400)
    result <- gextract("vt_lse_na", intervals = iter_int, iterator = iter_int, colnames = "value")

    # NA should be skipped; lse(1, 3)
    expected <- log(exp(1.0) + exp(3.0))
    expect_equal(result$value, expected, tolerance = 1e-5)
})

test_that("lse on value-based track with iterator windows", {
    intervals_df <- data.frame(
        chrom = "chr1",
        start = c(0, 100, 200, 300, 400),
        end = c(100, 200, 300, 400, 500),
        score = c(1.0, 2.0, 3.0, 4.0, 5.0)
    )

    gvtrack.create("vt_lse_iter_window", src = intervals_df, func = "lse")
    on.exit(gvtrack.rm("vt_lse_iter_window"), add = TRUE)

    # 200bp windows
    result <- gextract("vt_lse_iter_window", gintervals("chr1", 0, 400), iterator = 200, colnames = "value")

    # Window [0,200): lse(1, 2)
    expect_equal(result$value[1], log(exp(1) + exp(2)), tolerance = 1e-5)
    # Window [200,400): lse(3, 4)
    expect_equal(result$value[2], log(exp(3) + exp(4)), tolerance = 1e-5)
})

# ============================================================
# Numerical stability
# ============================================================

test_that("lse numerical stability with large values", {
    intervals_df <- data.frame(
        chrom = "chr1",
        start = c(100, 200),
        end = c(200, 300),
        score = c(100.0, 101.0)
    )

    gvtrack.create("test.lse.large", src = intervals_df, func = "lse")
    on.exit(gvtrack.rm("test.lse.large"), add = TRUE)

    iter_int <- gintervals("chr1", 100, 300)
    result <- gextract("test.lse.large", intervals = iter_int, iterator = iter_int, colnames = "value")

    # 101 + log(1 + exp(-1))
    expected <- 101.0 + log(1 + exp(-1.0))
    expect_equal(result$value, expected, tolerance = 1e-4)
})

test_that("lse numerical stability with very large values", {
    intervals_df <- data.frame(
        chrom = "chr1",
        start = c(100, 200, 300),
        end = c(200, 300, 400),
        score = c(80.0, 85.0, 80.0)
    )

    gvtrack.create("vt_lse_vlarge", src = intervals_df, func = "lse")
    on.exit(gvtrack.rm("vt_lse_vlarge"), add = TRUE)

    iter_int <- gintervals("chr1", 100, 400)
    result <- gextract("vt_lse_vlarge", intervals = iter_int, iterator = iter_int, colnames = "value")

    expected <- r_lse(c(80, 85, 80))
    expect_equal(result$value, expected, tolerance = 1e-3)
})

test_that("lse numerical stability with very negative values", {
    intervals_df <- data.frame(
        chrom = "chr1",
        start = c(100, 200, 300),
        end = c(200, 300, 400),
        score = c(-80.0, -85.0, -80.0)
    )

    gvtrack.create("vt_lse_vneg", src = intervals_df, func = "lse")
    on.exit(gvtrack.rm("vt_lse_vneg"), add = TRUE)

    iter_int <- gintervals("chr1", 100, 400)
    result <- gextract("vt_lse_vneg", intervals = iter_int, iterator = iter_int, colnames = "value")

    expected <- r_lse(c(-80, -85, -80))
    expect_equal(result$value, expected, tolerance = 1e-3)
})

test_that("lse with identical values", {
    intervals_df <- data.frame(
        chrom = "chr1",
        start = c(100, 200, 300),
        end = c(200, 300, 400),
        score = c(5.0, 5.0, 5.0)
    )

    gvtrack.create("vt_lse_identical", src = intervals_df, func = "lse")
    on.exit(gvtrack.rm("vt_lse_identical"), add = TRUE)

    iter_int <- gintervals("chr1", 100, 400)
    result <- gextract("vt_lse_identical", intervals = iter_int, iterator = iter_int, colnames = "value")

    # lse(x, x, x) = x + log(3)
    expected <- 5.0 + log(3)
    expect_equal(result$value, expected, tolerance = 1e-5)
})

test_that("lse with zero values", {
    intervals_df <- data.frame(
        chrom = "chr1",
        start = c(100, 200, 300),
        end = c(200, 300, 400),
        score = c(0.0, 0.0, 0.0)
    )

    gvtrack.create("vt_lse_zero", src = intervals_df, func = "lse")
    on.exit(gvtrack.rm("vt_lse_zero"), add = TRUE)

    iter_int <- gintervals("chr1", 100, 400)
    result <- gextract("vt_lse_zero", intervals = iter_int, iterator = iter_int, colnames = "value")

    # lse(0, 0, 0) = log(3 * exp(0)) = log(3)
    expect_equal(result$value, log(3), tolerance = 1e-5)
})

test_that("lse with mixed positive and negative values", {
    intervals_df <- data.frame(
        chrom = "chr1",
        start = c(100, 200, 300),
        end = c(200, 300, 400),
        score = c(-2.0, 0.0, 2.0)
    )

    gvtrack.create("vt_lse_mixed", src = intervals_df, func = "lse")
    on.exit(gvtrack.rm("vt_lse_mixed"), add = TRUE)

    iter_int <- gintervals("chr1", 100, 400)
    result <- gextract("vt_lse_mixed", intervals = iter_int, iterator = iter_int, colnames = "value")

    expected <- log(exp(-2) + exp(0) + exp(2))
    expect_equal(result$value, expected, tolerance = 1e-5)
})

test_that("lse with wide value range", {
    intervals_df <- data.frame(
        chrom = "chr1",
        start = c(100, 200),
        end = c(200, 300),
        score = c(-50.0, 50.0)
    )

    gvtrack.create("vt_lse_wide", src = intervals_df, func = "lse")
    on.exit(gvtrack.rm("vt_lse_wide"), add = TRUE)

    iter_int <- gintervals("chr1", 100, 300)
    result <- gextract("vt_lse_wide", intervals = iter_int, iterator = iter_int, colnames = "value")

    # When difference is huge, lse ~ max(values)
    # log(exp(-50) + exp(50)) = 50 + log(1 + exp(-100)) ~ 50
    expect_equal(result$value, 50.0, tolerance = 1e-4)
})

# ============================================================
# LSE property: always >= max(values)
# ============================================================

test_that("lse is always >= max of input values", {
    intervals_df <- data.frame(
        chrom = "chr1",
        start = c(100, 200, 300, 400),
        end = c(200, 300, 400, 500),
        score = c(-3.0, 1.0, -0.5, 2.0)
    )

    gvtrack.create("vt_lse_gemax", src = intervals_df, func = "lse")
    gvtrack.create("vt_max_gemax", src = intervals_df, func = "max")
    on.exit(
        {
            gvtrack.rm("vt_lse_gemax")
            gvtrack.rm("vt_max_gemax")
        },
        add = TRUE
    )

    iter_int <- gintervals("chr1", 100, 500)
    lse_result <- gextract("vt_lse_gemax", intervals = iter_int, iterator = iter_int, colnames = "value")
    max_result <- gextract("vt_max_gemax", intervals = iter_int, iterator = iter_int, colnames = "value")

    expect_gte(lse_result$value, max_result$value)
})

test_that("lse >= max holds over multiple intervals on dense track", {
    gvtrack.create("vt_lse_prop", src = "test.fixedbin", func = "lse")
    gvtrack.create("vt_max_prop", src = "test.fixedbin", func = "max")
    on.exit(
        {
            gvtrack.rm("vt_lse_prop")
            gvtrack.rm("vt_max_prop")
        },
        add = TRUE
    )

    intervals <- rbind(
        gintervals(1, 0, 500),
        gintervals(1, 500, 1000),
        gintervals(1, 1000, 2000)
    )

    lse_res <- gextract("vt_lse_prop", intervals = intervals, iterator = intervals, colnames = "value")
    max_res <- gextract("vt_max_prop", intervals = intervals, iterator = intervals, colnames = "value")

    both_valid <- !is.na(lse_res$value) & !is.na(max_res$value)
    if (any(both_valid)) {
        expect_true(all(lse_res$value[both_valid] >= max_res$value[both_valid]))
    }
})

# ============================================================
# Track expressions
# ============================================================

test_that("lse works in track expressions", {
    intervals_df <- data.frame(
        chrom = "chr1",
        start = c(100, 200),
        end = c(200, 300),
        score = c(1.0, 2.0)
    )

    gvtrack.create("vt_lse_expr", src = intervals_df, func = "lse")
    on.exit(gvtrack.rm("vt_lse_expr"), add = TRUE)

    iter_int <- gintervals("chr1", 100, 300)
    result <- gextract("vt_lse_expr * 2", intervals = iter_int, iterator = iter_int, colnames = "value")

    expected <- log(exp(1) + exp(2)) * 2
    expect_equal(result$value, expected, tolerance = 1e-5)
})

test_that("lse combined with other vtracks in expression", {
    intervals_df <- data.frame(
        chrom = "chr1",
        start = c(100, 200),
        end = c(200, 300),
        score = c(1.0, 2.0)
    )

    gvtrack.create("vt_lse_combo", src = intervals_df, func = "lse")
    gvtrack.create("vt_sum_combo", src = intervals_df, func = "sum")
    on.exit(
        {
            gvtrack.rm("vt_lse_combo")
            gvtrack.rm("vt_sum_combo")
        },
        add = TRUE
    )

    iter_int <- gintervals("chr1", 100, 300)
    result <- gextract("vt_lse_combo - vt_sum_combo", intervals = iter_int, iterator = iter_int, colnames = "value")

    lse_val <- log(exp(1) + exp(2))
    sum_val <- 3.0
    expect_equal(result$value, lse_val - sum_val, tolerance = 1e-5)
})

test_that("lse in conditional expression", {
    intervals_df <- data.frame(
        chrom = "chr1",
        start = c(100, 200),
        end = c(200, 300),
        score = c(1.0, 2.0)
    )

    gvtrack.create("vt_lse_cond", src = intervals_df, func = "lse")
    on.exit(gvtrack.rm("vt_lse_cond"), add = TRUE)

    iter_int <- gintervals("chr1", 100, 300)
    result <- gextract("ifelse(vt_lse_cond > 2, 1, 0)", intervals = iter_int, iterator = iter_int, colnames = "value")

    lse_val <- log(exp(1) + exp(2))
    expect_equal(result$value, ifelse(lse_val > 2, 1, 0))
})

# ============================================================
# gscreen / gsummary integration
# ============================================================

test_that("lse vtrack works with gscreen", {
    intervals_df <- data.frame(
        chrom = "chr1",
        start = c(100, 200, 300, 500, 600),
        end = c(200, 300, 400, 600, 700),
        score = c(1.0, 2.0, 3.0, 0.1, 0.2)
    )

    gvtrack.create("vt_lse_screen", src = intervals_df, func = "lse")
    on.exit(gvtrack.rm("vt_lse_screen"), add = TRUE)

    # Screen for intervals where lse > 2
    screened <- gscreen("vt_lse_screen > 2", gintervals("chr1", 0, 800), iterator = 200)
    expect_true(nrow(screened) > 0)

    # Verify: extract over screened intervals
    gextract_res <- gextract("vt_lse_screen", screened, iterator = screened, colnames = "value")
    expect_true(all(gextract_res$value > 2))
})

test_that("lse vtrack works with gsummary", {
    gvtrack.create("vt_lse_summary", src = "test.fixedbin", func = "lse")
    on.exit(gvtrack.rm("vt_lse_summary"), add = TRUE)

    summary_result <- gsummary("vt_lse_summary", gintervals(c(1, 2)), iterator = 5000)

    expect_true(!is.null(summary_result))
    expect_true("Mean" %in% names(summary_result) || "mean" %in% tolower(names(summary_result)))
})

# ============================================================
# NA pattern consistency with other functions
# ============================================================

test_that("lse NA pattern matches sum on dense track", {
    gvtrack.create("vt_lse_na_dense", src = "test.fixedbin", func = "lse")
    gvtrack.create("vt_sum_na_dense", src = "test.fixedbin", func = "sum")
    on.exit(
        {
            gvtrack.rm("vt_lse_na_dense")
            gvtrack.rm("vt_sum_na_dense")
        },
        add = TRUE
    )

    intervals <- gintervals(c(1, 2))
    lse_res <- gextract("vt_lse_na_dense", intervals, iterator = 2000, colnames = "value")
    sum_res <- gextract("vt_sum_na_dense", intervals, iterator = 2000, colnames = "value")

    expect_equal(is.na(lse_res$value), is.na(sum_res$value))
})

test_that("lse NA pattern matches sum on sparse track with intervals iterator", {
    gvtrack.create("vt_lse_na_sparse", src = "test.sparse", func = "lse")
    gvtrack.create("vt_sum_na_sparse", src = "test.sparse", func = "sum")
    on.exit(
        {
            gvtrack.rm("vt_lse_na_sparse")
            gvtrack.rm("vt_sum_na_sparse")
        },
        add = TRUE
    )

    intervals <- rbind(
        gintervals(1, 0, 500),
        gintervals(1, 500, 1000),
        gintervals(1, 5000, 6000),
        gintervals(1, 10000, 11000)
    )

    lse_res <- gextract("vt_lse_na_sparse", intervals, iterator = intervals, colnames = "value")
    sum_res <- gextract("vt_sum_na_sparse", intervals, iterator = intervals, colnames = "value")

    expect_equal(is.na(lse_res$value), is.na(sum_res$value))
})

# ============================================================
# Filters
# ============================================================

test_that("lse respects vtrack filter on sparse track", {
    gvtrack.create("vt_lse_filter_s", src = "test.sparse", func = "lse")
    gvtrack.create("vt_lse_nofilt_s", src = "test.sparse", func = "lse")

    # Restrictive filter on chr1
    gvtrack.filter("vt_lse_filter_s", filter = gintervals(1, 100, 200))
    on.exit(
        {
            gvtrack.rm("vt_lse_filter_s")
            gvtrack.rm("vt_lse_nofilt_s")
        },
        add = TRUE
    )

    iter_int <- gintervals(1, 0, 500)

    filt_res <- gextract("vt_lse_filter_s", iter_int, iterator = iter_int, colnames = "value")
    nofilt_res <- gextract("vt_lse_nofilt_s", iter_int, iterator = iter_int, colnames = "value")

    # Filtered has fewer (or equal) values, so LSE should be <= unfiltered
    if (!is.na(filt_res$value) && !is.na(nofilt_res$value)) {
        expect_lte(filt_res$value, nofilt_res$value)
    }
})

test_that("lse respects vtrack filter on dense track", {
    gvtrack.create("vt_lse_filter_dense", src = "test.fixedbin", func = "lse")
    gvtrack.create("vt_lse_nofilt_dense", src = "test.fixedbin", func = "lse")

    # Restrictive filter
    gvtrack.filter("vt_lse_filter_dense", filter = gintervals(1, 100, 200))
    on.exit(
        {
            gvtrack.rm("vt_lse_filter_dense")
            gvtrack.rm("vt_lse_nofilt_dense")
        },
        add = TRUE
    )

    iter_int <- gintervals(1, 0, 500)

    filt_res <- gextract("vt_lse_filter_dense", iter_int, iterator = iter_int, colnames = "value")
    nofilt_res <- gextract("vt_lse_nofilt_dense", iter_int, iterator = iter_int, colnames = "value")

    # Filtered should include fewer values, so LSE should be <= unfiltered
    if (!is.na(filt_res$value) && !is.na(nofilt_res$value)) {
        expect_lte(filt_res$value, nofilt_res$value)
    }
})

# ============================================================
# Edge cases
# ============================================================

test_that("lse with single bin on dense track", {
    gvtrack.create("vt_lse_single_bin", src = "test.fixedbin", func = "lse")
    on.exit(gvtrack.rm("vt_lse_single_bin"), add = TRUE)

    # Use the track's own bin size as iterator to get single-bin results
    result <- gextract("vt_lse_single_bin", gintervals(1, 0, 100), iterator = "test.fixedbin", colnames = "value")

    # For single bins, lse(x) = x
    raw <- gextract("test.fixedbin", gintervals(1, 0, 100), iterator = "test.fixedbin")
    both_valid <- !is.na(result$value) & !is.na(raw$test.fixedbin)
    if (any(both_valid)) {
        expect_equal(result$value[both_valid], raw$test.fixedbin[both_valid], tolerance = 1e-6)
    }
})

test_that("lse returns NA for empty value-based track query", {
    intervals_df <- data.frame(
        chrom = "chr1",
        start = 100,
        end = 200,
        score = 1.0
    )

    gvtrack.create("vt_lse_empty_query", src = intervals_df, func = "lse")
    on.exit(gvtrack.rm("vt_lse_empty_query"), add = TRUE)

    # Query a different chromosome
    iter_int <- gintervals("chr2", 100, 200)
    result <- gextract("vt_lse_empty_query", intervals = iter_int, iterator = iter_int, colnames = "value")

    expect_true(is.na(result$value))
})

test_that("lse on many consecutive NA intervals returns NA for each", {
    intervals_df <- data.frame(
        chrom = "chr1",
        start = 100,
        end = 200,
        score = 1.0
    )

    gvtrack.create("vt_lse_multi_na", src = intervals_df, func = "lse")
    on.exit(gvtrack.rm("vt_lse_multi_na"), add = TRUE)

    # Multiple query intervals, all outside data range
    query_ints <- rbind(
        gintervals("chr1", 500, 600),
        gintervals("chr1", 600, 700),
        gintervals("chr1", 700, 800)
    )

    result <- gextract("vt_lse_multi_na", intervals = query_ints, iterator = query_ints, colnames = "value")

    expect_equal(nrow(result), 3)
    expect_true(all(is.na(result$value)))
})

test_that("lse alternating data and no-data intervals", {
    intervals_df <- data.frame(
        chrom = "chr1",
        start = c(100, 500),
        end = c(200, 600),
        score = c(2.0, 3.0)
    )

    gvtrack.create("vt_lse_alt", src = intervals_df, func = "lse")
    on.exit(gvtrack.rm("vt_lse_alt"), add = TRUE)

    query_ints <- rbind(
        gintervals("chr1", 100, 200), # has data
        gintervals("chr1", 300, 400), # no data
        gintervals("chr1", 500, 600), # has data
        gintervals("chr1", 700, 800) # no data
    )

    result <- gextract("vt_lse_alt", intervals = query_ints, iterator = query_ints, colnames = "value")

    expect_equal(result$value[1], 2.0, tolerance = 1e-6) # lse(2) = 2
    expect_true(is.na(result$value[2]))
    expect_equal(result$value[3], 3.0, tolerance = 1e-6) # lse(3) = 3
    expect_true(is.na(result$value[4]))
})

test_that("lse with two values and known mathematical identity", {
    # lse(a, b) = max(a,b) + log(1 + exp(-|a-b|))
    a <- 3.0
    b <- 7.0
    intervals_df <- data.frame(
        chrom = "chr1",
        start = c(100, 200),
        end = c(200, 300),
        score = c(a, b)
    )

    gvtrack.create("vt_lse_identity", src = intervals_df, func = "lse")
    on.exit(gvtrack.rm("vt_lse_identity"), add = TRUE)

    iter_int <- gintervals("chr1", 100, 300)
    result <- gextract("vt_lse_identity", intervals = iter_int, iterator = iter_int, colnames = "value")

    expected <- max(a, b) + log(1 + exp(-abs(a - b)))
    expect_equal(result$value, expected, tolerance = 1e-5)
})

test_that("lse of many equal values equals value + log(n)", {
    n <- 10
    val <- 4.0
    intervals_df <- data.frame(
        chrom = "chr1",
        start = seq(0, by = 100, length.out = n),
        end = seq(100, by = 100, length.out = n),
        score = rep(val, n)
    )

    gvtrack.create("vt_lse_nequal", src = intervals_df, func = "lse")
    on.exit(gvtrack.rm("vt_lse_nequal"), add = TRUE)

    iter_int <- gintervals("chr1", 0, n * 100)
    result <- gextract("vt_lse_nequal", intervals = iter_int, iterator = iter_int, colnames = "value")

    expected <- val + log(n)
    expect_equal(result$value, expected, tolerance = 1e-4)
})

# ============================================================
# Monotonicity: adding more values can only increase LSE
# ============================================================

test_that("lse is monotonically non-decreasing as more values are added", {
    intervals_df <- data.frame(
        chrom = "chr1",
        start = c(100, 200, 300, 400),
        end = c(200, 300, 400, 500),
        score = c(1.0, 2.0, 3.0, 4.0)
    )

    gvtrack.create("vt_lse_mono", src = intervals_df, func = "lse")
    on.exit(gvtrack.rm("vt_lse_mono"), add = TRUE)

    # Increasing window sizes
    lse_values <- sapply(2:5, function(end_100) {
        iter_int <- gintervals("chr1", 100, end_100 * 100)
        gextract("vt_lse_mono", intervals = iter_int, iterator = iter_int, colnames = "value")$value
    })

    # Each value should be >= the previous one
    for (i in 2:length(lse_values)) {
        expect_gte(lse_values[i], lse_values[i - 1])
    }
})

# ============================================================
# Comparison: lse vs sum (lse >= sum when all values >= 0)
# ============================================================

test_that("lse equals max + log(1 + sum(exp(xi - max))) for known values", {
    vals <- c(1.0, 2.0, 3.0)
    intervals_df <- data.frame(
        chrom = "chr1",
        start = c(100, 200, 300),
        end = c(200, 300, 400),
        score = vals
    )

    gvtrack.create("vt_lse_formula", src = intervals_df, func = "lse")
    on.exit(gvtrack.rm("vt_lse_formula"), add = TRUE)

    iter_int <- gintervals("chr1", 100, 400)
    lse_res <- gextract("vt_lse_formula", intervals = iter_int, iterator = iter_int, colnames = "value")

    # Verify the log-sum-exp identity: lse(x) = max(x) + log(sum(exp(x - max(x))))
    m <- max(vals)
    expected <- m + log(sum(exp(vals - m)))
    expect_equal(lse_res$value, expected, tolerance = 1e-5)
})

# ============================================================
# Multiple vtracks simultaneously
# ============================================================

test_that("multiple lse vtracks can be used simultaneously", {
    df1 <- data.frame(
        chrom = "chr1",
        start = c(100, 200),
        end = c(200, 300),
        score = c(1.0, 2.0)
    )
    df2 <- data.frame(
        chrom = "chr1",
        start = c(100, 200),
        end = c(200, 300),
        score = c(3.0, 4.0)
    )

    gvtrack.create("vt_lse_multi1", src = df1, func = "lse")
    gvtrack.create("vt_lse_multi2", src = df2, func = "lse")
    on.exit(
        {
            gvtrack.rm("vt_lse_multi1")
            gvtrack.rm("vt_lse_multi2")
        },
        add = TRUE
    )

    iter_int <- gintervals("chr1", 100, 300)
    result <- gextract("vt_lse_multi1", "vt_lse_multi2", intervals = iter_int, iterator = iter_int)

    expect_equal(result$vt_lse_multi1, log(exp(1) + exp(2)), tolerance = 1e-5)
    expect_equal(result$vt_lse_multi2, log(exp(3) + exp(4)), tolerance = 1e-5)
})

test_that("lse and other functions work together", {
    df <- data.frame(
        chrom = "chr1",
        start = c(100, 200),
        end = c(200, 300),
        score = c(1.0, 2.0)
    )

    gvtrack.create("vt_lse_together", src = df, func = "lse")
    gvtrack.create("vt_sum_together", src = df, func = "sum")
    gvtrack.create("vt_avg_together", src = df, func = "avg")
    gvtrack.create("vt_max_together", src = df, func = "max")
    on.exit(
        {
            gvtrack.rm("vt_lse_together")
            gvtrack.rm("vt_sum_together")
            gvtrack.rm("vt_avg_together")
            gvtrack.rm("vt_max_together")
        },
        add = TRUE
    )

    iter_int <- gintervals("chr1", 100, 300)
    result <- gextract(
        "vt_lse_together", "vt_sum_together", "vt_avg_together", "vt_max_together",
        intervals = iter_int, iterator = iter_int
    )

    expect_equal(result$vt_lse_together, log(exp(1) + exp(2)), tolerance = 1e-5)
    expect_equal(result$vt_sum_together, 3.0)
    expect_equal(result$vt_avg_together, 1.5)
    expect_equal(result$vt_max_together, 2.0)
})

# ============================================================
# Dense track: multi-chromosome
# ============================================================

test_that("lse on dense track across multiple chromosomes", {
    gvtrack.create("vt_lse_dense_mc", src = "test.fixedbin", func = "lse")
    on.exit(gvtrack.rm("vt_lse_dense_mc"), add = TRUE)

    intervals <- rbind(
        gintervals(1, 0, 500),
        gintervals(2, 0, 500)
    )

    result <- gextract("vt_lse_dense_mc", intervals, iterator = intervals, colnames = "value")

    expect_equal(nrow(result), 2)
    # Verify per-chromosome against manual computation
    for (i in seq_len(nrow(result))) {
        if (!is.na(result$value[i])) {
            manual <- manual_lse("test.fixedbin", intervals[i, ])
            expect_equal(result$value[i], manual, tolerance = 1e-3)
        }
    }
})

# ============================================================
# Sparse track: full comparison against manual
# ============================================================

test_that("lse on sparse track full verification with manual lse", {
    gvtrack.create("vt_lse_full_sparse", src = "test.sparse", func = "lse")
    on.exit(gvtrack.rm("vt_lse_full_sparse"), add = TRUE)

    # Use sliding windows on chr1
    intervals <- rbind(
        gintervals(1, 0, 200),
        gintervals(1, 200, 400),
        gintervals(1, 400, 600),
        gintervals(1, 600, 800),
        gintervals(1, 800, 1000)
    )

    result <- gextract("vt_lse_full_sparse", intervals, iterator = intervals, colnames = "value")
    manual <- manual_lse_vec("test.sparse", intervals)

    expect_equal(is.na(result$value), is.na(manual))
    both_valid <- !is.na(result$value) & !is.na(manual)
    if (any(both_valid)) {
        expect_equal(result$value[both_valid], manual[both_valid], tolerance = 1e-3)
    }
})

# ============================================================
# PWM decomposition: LSE of fine-grained PWM = coarse PWM
# ============================================================

test_that("LSE of 20bp PWM values equals 200bp PWM directly", {
    # PWM vtrack internally computes LSE over per-position motif scores.
    # Therefore: PWM(200bp) = LSE(all positions in 200bp)
    #          = LSE(LSE(pos in 20bp_1), ..., LSE(pos in 20bp_10))
    # because exp(LSE(x)) = sum(exp(x_i)).
    #
    # We verify this by:
    # 1. Extracting PWM at 20bp resolution
    # 2. Storing those values in a value-based track
    # 3. Computing LSE over 200bp windows
    # 4. Comparing with direct PWM at 200bp resolution

    pssm <- create_test_pssm() # 2-position motif (AC)

    # Region on chr1 with coordinates divisible by both 20 and 200
    region <- gintervals(1, 10000, 11000)

    # Create PWM vtrack with extend=TRUE so boundary anchors are included
    gvtrack.create("vt_pwm_fine", NULL,
        func = "pwm",
        pssm = pssm, bidirect = FALSE, extend = TRUE, prior = 0.01
    )
    on.exit(gvtrack.rm("vt_pwm_fine"), add = TRUE)

    # Step 1: Extract PWM at 20bp resolution -> 50 windows
    fine <- gextract("vt_pwm_fine", region, iterator = 20)
    expect_equal(nrow(fine), 50)

    # Step 2: Extract PWM at 200bp resolution -> 5 windows
    coarse <- gextract("vt_pwm_fine", region, iterator = 200)
    expect_equal(nrow(coarse), 5)

    # Step 3: Build a value-based track from the 20bp PWM values
    fine_df <- data.frame(
        chrom = fine$chrom,
        start = fine$start,
        end = fine$end,
        score = fine$vt_pwm_fine
    )

    gvtrack.create("vt_lse_of_pwm", src = fine_df, func = "lse")
    on.exit(gvtrack.rm("vt_lse_of_pwm"), add = TRUE)

    # Step 4: Extract LSE at 200bp resolution
    lse_200 <- gextract("vt_lse_of_pwm", region, iterator = 200, colnames = "value")
    expect_equal(nrow(lse_200), 5)

    # Step 5: Compare - should be equal (up to float precision)
    expect_equal(lse_200$value, coarse$vt_pwm_fine, tolerance = 1e-4)
})
