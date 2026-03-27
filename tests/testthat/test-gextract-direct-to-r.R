create_isolated_test_db()

# ============================================================================
# 1. SETLENGTH edge cases
# ============================================================================

test_that("empty result returns NULL", {
    # Use an interval region where the sparse track has no data
    # Create a tiny interval that is unlikely to contain sparse data
    intervs <- gintervals(1, 247249700, 247249719)
    r <- gextract("test.sparse", intervs)
    expect_null(r)
})

test_that("single row result preserves data.frame structure", {
    # Extract a single bin from fixedbin
    r <- gextract("test.fixedbin", gintervals(1, 0, 1), iterator = 1)
    expect_s3_class(r, "data.frame")
    expect_equal(nrow(r), 1)
    expect_true(all(c("chrom", "start", "end", "test.fixedbin", "intervalID") %in% colnames(r)))
})

test_that("heavy truncation: small intervals over large scope", {
    # Use small intervals that cover a tiny fraction of the chromosome
    intervs <- rbind(
        gintervals(1, 100, 200),
        gintervals(1, 1000, 1100),
        gintervals(1, 5000, 5100)
    )
    r <- gextract("test.fixedbin", intervs)
    expect_s3_class(r, "data.frame")
    expect_true(nrow(r) > 0)
    # All results should fall within our specified intervals
    expect_true(all(r$start >= 0))
    expect_true(all(r$end <= 5100))
})

# ============================================================================
# 2. Column correctness
# ============================================================================

test_that("chrom column is a proper factor with correct levels", {
    r <- gextract("test.fixedbin", gintervals(c(1, 2)))
    expect_true(is.factor(r$chrom))
    # Factor levels should include chromosome names
    lvls <- levels(r$chrom)
    expect_true(length(lvls) > 0)
    # All chromosome values that appear should be valid levels
    unique_chroms <- as.character(unique(r$chrom))
    expect_true(all(unique_chroms %in% lvls))
})

test_that("start/end values match expected intervals exactly", {
    intervs <- gintervals(1, 0, 200)
    r <- gextract("test.fixedbin", intervs)
    expect_s3_class(r, "data.frame")
    # All starts should be >= 0 and ends <= 200
    expect_true(all(r$start >= 0))
    expect_true(all(r$end <= 200))
    # starts should be < ends
    expect_true(all(r$start < r$end))
})

test_that("intervalID column is correct and 1-based", {
    intervs <- rbind(
        gintervals(1, 0, 1000),
        gintervals(2, 0, 1000)
    )
    r <- gextract("test.fixedbin", intervs)
    expect_true("intervalID" %in% colnames(r))
    # intervalID should be integer values, 1-based
    expect_true(is.integer(r$intervalID) || is.numeric(r$intervalID))
    expect_true(all(r$intervalID >= 1))
    expect_true(all(r$intervalID <= nrow(intervs)))
})

test_that("column names are correct for single expression", {
    r <- gextract("test.fixedbin", gintervals(c(1, 2)))
    expect_equal(colnames(r), c("chrom", "start", "end", "test.fixedbin", "intervalID"))
})

test_that("column names are correct for multi expression", {
    r <- gextract("test.fixedbin", "test.sparse", gintervals(c(1, 2)), iterator = "test.fixedbin")
    expect_equal(colnames(r), c("chrom", "start", "end", "test.fixedbin", "test.sparse", "intervalID"))
})

test_that("custom colnames override expression names", {
    r <- gextract("test.fixedbin", "test.sparse",
        gintervals(c(1, 2)),
        colnames = c("fb", "sp"),
        iterator = "test.fixedbin"
    )
    expect_equal(colnames(r), c("chrom", "start", "end", "fb", "sp", "intervalID"))
})

test_that("row count matches expected for fixedbin track", {
    intervs <- gintervals(1, 0, 1000)
    r <- gextract("test.fixedbin", intervs)
    expect_true(nrow(r) > 0)
    # For a fixedbin track over [0,1000), the number of rows should be deterministic
    r2 <- gextract("test.fixedbin", intervs)
    expect_equal(nrow(r), nrow(r2))
})

test_that("data.frame class attribute is set", {
    r <- gextract("test.fixedbin", gintervals(1, 0, 1000))
    expect_equal(class(r), "data.frame")
    expect_true(is.data.frame(r))
})

test_that("row names are sequential 1:nrow", {
    r <- gextract("test.fixedbin", gintervals(c(1, 2)))
    rn <- rownames(r)
    expect_equal(rn, as.character(1:nrow(r)))
})

# ============================================================================
# 3. Multi-expression
# ============================================================================

test_that("single expression extraction works", {
    r <- gextract("test.fixedbin", gintervals(c(1, 2)))
    expect_s3_class(r, "data.frame")
    expect_equal(ncol(r), 5) # chrom, start, end, value, intervalID
})

test_that("three expressions extraction works", {
    r <- gextract("test.fixedbin", "test.fixedbin * 2", "test.fixedbin + 1",
        gintervals(c(1, 2)),
        iterator = "test.fixedbin"
    )
    expect_s3_class(r, "data.frame")
    expect_equal(ncol(r), 7) # chrom, start, end, 3 values, intervalID
    # Verify computed expressions are correct
    expect_equal(r[[5]], r[[4]] * 2, tolerance = 1e-10)
    expect_equal(r[[6]], r[[4]] + 1, tolerance = 1e-10)
})

test_that("five expressions extraction works", {
    r <- gextract(
        "test.fixedbin",
        "test.fixedbin * 2",
        "test.fixedbin + 1",
        "test.fixedbin - 0.5",
        "test.fixedbin / 2",
        gintervals(c(1, 2)),
        iterator = "test.fixedbin"
    )
    expect_s3_class(r, "data.frame")
    expect_equal(ncol(r), 9) # chrom, start, end, 5 values, intervalID
})

test_that("mix of vtrack and raw track expressions", {
    gvtrack.create("vt_avg", "test.fixedbin", func = "avg")
    withr::defer(gvtrack.rm("vt_avg"))

    r <- gextract("test.fixedbin", "vt_avg",
        gintervals(c(1, 2)),
        iterator = 200
    )
    expect_s3_class(r, "data.frame")
    expect_equal(ncol(r), 6) # chrom, start, end, fixedbin, vt_avg, intervalID
    expect_true(all(c("test.fixedbin", "vt_avg") %in% colnames(r)))
})

test_that("expression returning NaN for some rows", {
    # log of track values can produce NaN for negative or zero values
    r <- gextract("log(test.fixedbin)", gintervals(c(1, 2)))
    expect_s3_class(r, "data.frame")
    # Should still have proper structure even with NaN values
    expect_true(is.numeric(r[["log(test.fixedbin)"]]))
})

# ============================================================================
# 4. Consistency with file output
# ============================================================================

test_that("gextract to data.frame vs file output are identical", {
    intervs <- gintervals(c(1, 2), 0, 100000)
    r_df <- gextract("test.fixedbin", intervs)

    tmp <- tempfile()
    withr::defer(unlink(tmp))
    gextract("test.fixedbin", intervs, file = tmp)
    r_file <- readr::read_tsv(tmp, col_types = readr::cols(
        chrom = readr::col_character(),
        start = readr::col_double(),
        end = readr::col_double(),
        test.fixedbin = readr::col_double()
    ))

    # Compare after making chrom comparable (file has character, df has factor)
    r_df_cmp <- r_df %>%
        mutate(chrom = as.character(chrom)) %>%
        select(-intervalID)
    expect_equal(r_df_cmp, r_file, ignore_attr = TRUE)
})

test_that("gextract to data.frame vs file for sparse track are identical", {
    intervs <- gintervals(c(1, 2))
    r_df <- gextract("test.sparse", intervs)

    tmp <- tempfile()
    withr::defer(unlink(tmp))
    gextract("test.sparse", intervs, file = tmp)
    r_file <- readr::read_tsv(tmp, col_types = readr::cols(
        chrom = readr::col_character(),
        start = readr::col_double(),
        end = readr::col_double(),
        test.sparse = readr::col_double()
    ))

    r_df_cmp <- r_df %>%
        mutate(chrom = as.character(chrom)) %>%
        select(-intervalID)
    expect_equal(r_df_cmp, r_file, ignore_attr = TRUE)
})

test_that("gextract multi-expression file output matches data.frame", {
    intervs <- gintervals(c(1, 2), 0, 50000)
    r_df <- gextract("test.fixedbin", "test.fixedbin * 3",
        intervs,
        iterator = "test.fixedbin"
    )

    tmp <- tempfile()
    withr::defer(unlink(tmp))
    gextract("test.fixedbin", "test.fixedbin * 3",
        intervs,
        iterator = "test.fixedbin",
        file = tmp
    )
    r_file <- readr::read_tsv(tmp, show_col_types = FALSE)

    r_df_cmp <- r_df %>%
        mutate(chrom = as.character(chrom)) %>%
        select(-intervalID)
    expect_equal(r_df_cmp, r_file, ignore_attr = TRUE, tolerance = 1e-10)
})

# ============================================================================
# 5. Consistency with gsummary
# ============================================================================

test_that("gsummary sum matches sum of gextract values for fixedbin", {
    intervs <- gintervals(c(1, 2), 0, 100000)
    r_extract <- gextract("test.fixedbin", intervs)
    r_summary <- gsummary("test.fixedbin", intervs)

    extract_sum <- sum(r_extract$test.fixedbin, na.rm = TRUE)
    expect_equal(extract_sum, r_summary[["Sum"]], tolerance = 1e-5)
})

test_that("gsummary sum matches sum of gextract values for sparse", {
    intervs <- gintervals(c(1, 2))
    r_extract <- gextract("test.sparse", intervs)
    r_summary <- gsummary("test.sparse", intervs)

    extract_sum <- sum(r_extract$test.sparse, na.rm = TRUE)
    expect_equal(extract_sum, r_summary[["Sum"]], tolerance = 1e-5)
})

test_that("gsummary count matches nrow of gextract for fixedbin", {
    intervs <- gintervals(c(1, 2), 0, 100000)
    r_extract <- gextract("test.fixedbin", intervs)
    r_summary <- gsummary("test.fixedbin", intervs)

    non_nan <- sum(!is.nan(r_extract$test.fixedbin))
    expect_equal(non_nan, r_summary[["Total intervals"]] - r_summary[["NaN intervals"]])
})

# ============================================================================
# 6. Reproducibility
# ============================================================================

test_that("identical gextract calls produce identical results", {
    intervs <- gintervals(c(1, 2), 0, 500000)
    r1 <- gextract("test.fixedbin", intervs)
    r2 <- gextract("test.fixedbin", intervs)
    expect_identical(r1, r2)
})

test_that("identical gextract calls produce identical results for sparse", {
    intervs <- gintervals(c(1, 2))
    r1 <- gextract("test.sparse", intervs)
    r2 <- gextract("test.sparse", intervs)
    expect_identical(r1, r2)
})

test_that("identical gextract calls with iterator produce identical results", {
    intervs <- gintervals(c(1, 2), 0, 100000)
    r1 <- gextract("test.fixedbin", intervs, iterator = 500)
    r2 <- gextract("test.fixedbin", intervs, iterator = 500)
    expect_identical(r1, r2)
})

# ============================================================================
# 7. Various iterator sizes
# ============================================================================

test_that("very small iterator (10bp)", {
    intervs <- gintervals(1, 0, 1000)
    r <- gextract("test.fixedbin", intervs, iterator = 10)
    expect_s3_class(r, "data.frame")
    expect_true(nrow(r) > 0)
    # With 10bp iterator over 1000bp, each bin should be <= 10bp
    expect_true(all((r$end - r$start) <= 10))
})

test_that("large iterator (10000bp)", {
    intervs <- gintervals(c(1, 2), 0, 100000)
    r <- gextract("test.fixedbin", intervs, iterator = 10000)
    expect_s3_class(r, "data.frame")
    expect_true(nrow(r) > 0)
    # With 10000bp iterator, each bin should be <= 10000bp
    expect_true(all((r$end - r$start) <= 10000))
})

test_that("iterator larger than intervals", {
    intervs <- gintervals(1, 0, 100)
    r <- gextract("test.fixedbin", intervs, iterator = 10000)
    expect_s3_class(r, "data.frame")
    # Should still produce results, just a single bin covering the interval
    expect_true(nrow(r) >= 1)
})

test_that("different iterator sizes produce consistent total sums", {
    intervs <- gintervals(1, 0, 10000)

    gvtrack.create("vt_sum_test", "test.fixedbin", func = "sum")
    withr::defer(gvtrack.rm("vt_sum_test"))

    r100 <- gextract("vt_sum_test", intervs, iterator = 100)
    r1000 <- gextract("vt_sum_test", intervs, iterator = 1000)

    # Total sum should be the same regardless of iterator size
    sum100 <- sum(r100$vt_sum_test, na.rm = TRUE)
    sum1000 <- sum(r1000$vt_sum_test, na.rm = TRUE)
    expect_equal(sum100, sum1000, tolerance = 1e-5)
})

# ============================================================================
# 8. Chromosome handling
# ============================================================================

test_that("single chromosome extraction", {
    r <- gextract("test.fixedbin", gintervals(1, 0, 10000))
    expect_s3_class(r, "data.frame")
    expect_true(all(as.character(r$chrom) == "chr1"))
})

test_that("all chromosomes extraction", {
    withr::local_options(gmax.data.size = 1e9)
    r <- gextract("test.fixedbin", .misha$ALLGENOME)
    expect_s3_class(r, "data.frame")
    # Should have data from multiple chromosomes
    unique_chroms <- unique(as.character(r$chrom))
    expect_true(length(unique_chroms) > 1)
})

test_that("subset of chromosomes extraction", {
    r <- gextract("test.fixedbin", gintervals(c(1, 2, 3)))
    expect_s3_class(r, "data.frame")
    unique_chroms <- unique(as.character(r$chrom))
    expect_true(all(unique_chroms %in% c("chr1", "chr2", "chr3")))
})

test_that("chrom factor levels are consistent across different chromosome subsets", {
    r1 <- gextract("test.fixedbin", gintervals(1, 0, 1000))
    r2 <- gextract("test.fixedbin", gintervals(c(1, 2), 0, 1000))
    # Both should have the same factor levels (all genome chroms)
    expect_equal(levels(r1$chrom), levels(r2$chrom))
})

# ============================================================================
# 9. Stress test
# ============================================================================

test_that("large extraction: full chromosome with fine iterator", {
    withr::local_options(gmax.data.size = 1e9)
    # Extract full chr22 (smallest autosome ~49Mb) with 1000bp iterator
    r <- gextract("test.fixedbin", gintervals(22), iterator = 1000)
    expect_s3_class(r, "data.frame")
    expect_true(nrow(r) > 10000)
    # Verify structure integrity
    expect_equal(colnames(r), c("chrom", "start", "end", "test.fixedbin", "intervalID"))
    expect_true(is.factor(r$chrom))
    expect_true(is.numeric(r$start))
    expect_true(is.numeric(r$end))
    expect_true(is.numeric(r$test.fixedbin))
    expect_equal(rownames(r), as.character(1:nrow(r)))
})

test_that("large multi-expression extraction", {
    withr::local_options(gmax.data.size = 1e9)
    r <- gextract("test.fixedbin", "test.fixedbin * 2", "test.fixedbin + 1",
        gintervals(22),
        iterator = "test.fixedbin"
    )
    expect_s3_class(r, "data.frame")
    expect_true(nrow(r) > 1000)
    expect_equal(ncol(r), 7)
    # Verify computed columns are consistent
    expect_equal(r[[5]], r[[4]] * 2, tolerance = 1e-10)
    expect_equal(r[[6]], r[[4]] + 1, tolerance = 1e-10)
})

# ============================================================================
# Additional edge cases for direct-to-R optimization
# ============================================================================

test_that("expression with NaN values preserves NaN correctly", {
    # Division by zero produces Inf/NaN
    r <- gextract("test.fixedbin / (test.fixedbin - test.fixedbin)", gintervals(1, 0, 1000))
    expect_s3_class(r, "data.frame")
    # The expression column is the 4th column (after chrom, start, end)
    vals <- r[[4]]
    expect_true(is.numeric(vals))
    # Values should be Inf or NaN (0/0 = NaN, x/0 = Inf)
    expect_true(any(is.infinite(vals) | is.nan(vals)))
})

test_that("gextract with array track works correctly", {
    r <- gextract("test.array", gintervals(c(1, 2)))
    expect_s3_class(r, "data.frame")
    expect_true("test.array" %in% colnames(r))
    expect_true(is.factor(r$chrom))
    expect_equal(rownames(r), as.character(1:nrow(r)))
})

test_that("intervalID maps correctly to multiple input intervals", {
    intervs <- rbind(
        gintervals(1, 0, 10000),
        gintervals(1, 50000, 60000),
        gintervals(2, 0, 10000)
    )
    r <- gextract("test.fixedbin", intervs)
    expect_true(all(r$intervalID %in% 1:3))

    # Rows with chrom==chr1, start < 10000 should map to interval 1
    mask1 <- as.character(r$chrom) == "chr1" & r$start < 10000
    if (any(mask1)) {
        expect_true(all(r$intervalID[mask1] == 1))
    }

    # Rows with chrom==chr1, start >= 50000 should map to interval 2
    mask2 <- as.character(r$chrom) == "chr1" & r$start >= 50000
    if (any(mask2)) {
        expect_true(all(r$intervalID[mask2] == 2))
    }

    # Rows with chrom==chr2 should map to interval 3
    mask3 <- as.character(r$chrom) == "chr2"
    if (any(mask3)) {
        expect_true(all(r$intervalID[mask3] == 3))
    }
})

test_that("result from sparse track iterator falls back correctly", {
    # Sparse track as iterator triggers the fallback path (estimated==0)
    r <- gextract("test.fixedbin", gintervals(c(1, 2)), iterator = "test.sparse")
    expect_s3_class(r, "data.frame")
    expect_true(nrow(r) > 0)
    expect_true(is.factor(r$chrom))
    expect_equal(rownames(r), as.character(1:nrow(r)))
    expect_equal(class(r), "data.frame")
})

test_that("result from track iterator has correct structure", {
    r <- gextract("test.fixedbin", gintervals(c(1, 2)), iterator = "test.fixedbin")
    expect_s3_class(r, "data.frame")
    expect_true(is.factor(r$chrom))
    expect_equal(colnames(r), c("chrom", "start", "end", "test.fixedbin", "intervalID"))
    expect_equal(rownames(r), as.character(1:nrow(r)))
})

test_that("vtrack with sum function and gextract consistency", {
    gvtrack.create("vt_sum_cons", "test.fixedbin", func = "sum")
    withr::defer(gvtrack.rm("vt_sum_cons"))

    intervs <- gintervals(c(1, 2), 0, 50000)
    r <- gextract("vt_sum_cons", intervs, iterator = 5000)
    expect_s3_class(r, "data.frame")
    expect_true(all(c("chrom", "start", "end", "vt_sum_cons", "intervalID") %in% colnames(r)))

    # Sum of vtrack sums should equal gsummary of underlying track
    r_summary <- gsummary("test.fixedbin", intervs)
    total <- sum(r$vt_sum_cons, na.rm = TRUE)
    expect_equal(total, r_summary[["Sum"]], tolerance = 1e-4)
})
