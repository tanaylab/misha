test_that("round-trip preserves chrom/start/end and value columns", {
    df <- data.frame(
        chrom = c("chr1", "chr1", "chr2"),
        start = c(100L, 500L, 200L),
        end = c(200L, 700L, 400L),
        t1 = c(1.5, 2.5, 3.5),
        t2 = c(10.0, 20.0, 30.0),
        stringsAsFactors = FALSE
    )
    mat <- gintervals.to_mat(df)
    out <- gintervals.from_mat(mat)

    expect_equal(out$chrom, df$chrom)
    expect_equal(out$start, df$start)
    expect_equal(out$end, df$end)
    expect_equal(out$t1, df$t1)
    expect_equal(out$t2, df$t2)
    expect_type(out$start, "integer")
    expect_type(out$end, "integer")
})

test_that("id_col controls rownames; round-trip still works via attribute", {
    df <- data.frame(
        chrom = c("chr1", "chr2"),
        start = c(100L, 200L),
        end = c(200L, 400L),
        gene = c("FOO", "BAR"),
        t1 = c(1.0, 2.0),
        stringsAsFactors = FALSE
    )
    mat <- gintervals.to_mat(df, id_col = "gene", value_cols = "t1")
    expect_equal(rownames(mat), c("FOO", "BAR"))
    expect_equal(attr(mat, "intervals")$chrom, df$chrom)

    out <- gintervals.from_mat(mat)
    expect_equal(out$chrom, df$chrom)
    expect_equal(out$start, df$start)
    expect_equal(out$end, df$end)
    expect_equal(out$t1, df$t1)
})

test_that("default rownames are 'chrom:start-end'", {
    df <- data.frame(
        chrom = "chr1", start = 100L, end = 200L, t1 = 1.0,
        stringsAsFactors = FALSE
    )
    mat <- gintervals.to_mat(df)
    expect_equal(rownames(mat), "chr1:100-200")
})

test_that("missing id_col errors clearly", {
    df <- data.frame(
        chrom = "chr1", start = 100L, end = 200L, t1 = 1.0,
        stringsAsFactors = FALSE
    )
    expect_error(gintervals.to_mat(df, id_col = "nope"),
        "`id_col` not found",
        fixed = FALSE
    )
})

test_that("value_cols restricts which columns become matrix data", {
    df <- data.frame(
        chrom = "chr1", start = 100L, end = 200L,
        t1 = 1.0, t2 = 2.0, t3 = 3.0,
        stringsAsFactors = FALSE
    )
    mat <- gintervals.to_mat(df, value_cols = c("t1", "t3"))
    expect_equal(colnames(mat), c("t1", "t3"))
    expect_equal(as.numeric(mat[1, ]), c(1.0, 3.0))
})

test_that("auto-detect errors on non-numeric value column", {
    df <- data.frame(
        chrom = "chr1", start = 100L, end = 200L,
        gene = "FOO", t1 = 1.0,
        stringsAsFactors = FALSE
    )
    expect_error(gintervals.to_mat(df),
        "Non-numeric value column.*gene",
        fixed = FALSE
    )
})

test_that("value_cols missing in df errors clearly", {
    df <- data.frame(
        chrom = "chr1", start = 100L, end = 200L, t1 = 1.0,
        stringsAsFactors = FALSE
    )
    expect_error(gintervals.to_mat(df, value_cols = c("t1", "nope")),
        "`value_cols` not found",
        fixed = FALSE
    )
})

test_that("intervalID is excluded from value_cols by default but kept in attribute", {
    df <- data.frame(
        chrom = "chr1", start = 100L, end = 200L,
        intervalID = 1L, t1 = 1.0,
        stringsAsFactors = FALSE
    )
    mat <- gintervals.to_mat(df)
    expect_equal(colnames(mat), "t1")
    expect_true("intervalID" %in% names(attr(mat, "intervals")))

    out <- gintervals.from_mat(mat)
    expect_true("intervalID" %in% names(out))
    expect_equal(out$intervalID, df$intervalID)
})

test_that("from_mat with plain matrix + intervals arg works", {
    intervs <- data.frame(
        chrom = c("chr1", "chr1"),
        start = c(100L, 300L),
        end = c(200L, 400L),
        stringsAsFactors = FALSE
    )
    plain <- matrix(c(1.0, 2.0, 3.0, 4.0),
        nrow = 2,
        dimnames = list(NULL, c("a", "b"))
    )
    out <- gintervals.from_mat(plain, intervals = intervs)
    expect_equal(out$chrom, intervs$chrom)
    expect_equal(out$a, c(1.0, 2.0))
    expect_equal(out$b, c(3.0, 4.0))
})

test_that("from_mat without intervals on plain matrix errors clearly", {
    plain <- matrix(1:4, nrow = 2)
    expect_error(gintervals.from_mat(plain),
        "plain matrix",
        fixed = FALSE
    )
})

test_that("from_mat errors on row mismatch", {
    intervs <- data.frame(
        chrom = "chr1", start = 100L, end = 200L,
        stringsAsFactors = FALSE
    )
    plain <- matrix(1:4, nrow = 2)
    expect_error(gintervals.from_mat(plain, intervals = intervs),
        "must match",
        fixed = FALSE
    )
})

test_that("from_mat on intervs_mat + intervals arg errors (ambiguous)", {
    df <- data.frame(
        chrom = "chr1", start = 100L, end = 200L, t1 = 1.0,
        stringsAsFactors = FALSE
    )
    mat <- gintervals.to_mat(df)
    expect_error(gintervals.from_mat(mat, intervals = df),
        "intervs_mat",
        fixed = FALSE
    )
})

test_that("round-trip preserves chrom names containing underscores", {
    df <- data.frame(
        chrom = c("chr1", "chrUn_KI270442v1", "chrUn_GL000220v1"),
        start = c(100L, 5000L, 1L),
        end = c(200L, 6000L, 500L),
        t1 = c(1.0, 2.0, 3.0),
        stringsAsFactors = FALSE
    )
    mat <- gintervals.to_mat(df)
    out <- gintervals.from_mat(mat)
    expect_equal(out$chrom, df$chrom)
    expect_equal(out$start, df$start)
    expect_equal(out$end, df$end)
    expect_equal(out$t1, df$t1)
})
