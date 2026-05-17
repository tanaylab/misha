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

test_that("[ row subset keeps attr$intervals aligned", {
    df <- data.frame(
        chrom = c("chr1", "chr1", "chr2"),
        start = c(100L, 500L, 200L),
        end = c(200L, 700L, 400L),
        t1 = c(1, 2, 3), t2 = c(10, 20, 30),
        stringsAsFactors = FALSE
    )
    mat <- gintervals.to_mat(df)
    sub <- mat[c(1, 3), ]
    expect_s3_class(sub, "intervs_mat")
    expect_equal(nrow(sub), 2L)
    expect_equal(attr(sub, "intervals")$chrom, c("chr1", "chr2"))
    expect_equal(attr(sub, "intervals")$start, c(100L, 200L))
    expect_equal(rownames(sub), c("chr1:100-200", "chr2:200-400"))
})

test_that("[ column subset leaves attr$intervals unchanged", {
    df <- data.frame(
        chrom = "chr1", start = 100L, end = 200L,
        t1 = 1, t2 = 2, t3 = 3,
        stringsAsFactors = FALSE
    )
    mat <- gintervals.to_mat(df)
    sub <- mat[, c("t1", "t3")]
    expect_s3_class(sub, "intervs_mat")
    expect_equal(colnames(sub), c("t1", "t3"))
    expect_equal(attr(sub, "intervals"), attr(mat, "intervals"))
})

test_that("[ drop=TRUE on single row returns vector and drops class", {
    df <- data.frame(
        chrom = c("chr1", "chr2"),
        start = c(100L, 200L), end = c(200L, 400L),
        t1 = c(1, 2), t2 = c(10, 20),
        stringsAsFactors = FALSE
    )
    mat <- gintervals.to_mat(df)
    v <- mat[1, ]
    expect_false(inherits(v, "intervs_mat"))
    expect_null(attr(v, "intervals"))
    expect_equal(as.numeric(v), c(1, 10))
})

test_that("[ drop=FALSE on single row preserves class and attr", {
    df <- data.frame(
        chrom = c("chr1", "chr2"),
        start = c(100L, 200L), end = c(200L, 400L),
        t1 = c(1, 2),
        stringsAsFactors = FALSE
    )
    mat <- gintervals.to_mat(df)
    one <- mat[1, , drop = FALSE]
    expect_s3_class(one, "intervs_mat")
    expect_equal(nrow(one), 1L)
    expect_equal(attr(one, "intervals")$chrom, "chr1")
})

test_that("from_mat round-trips after row subset", {
    df <- data.frame(
        chrom = c("chr1", "chr2", "chr3"),
        start = c(100L, 200L, 300L),
        end = c(200L, 400L, 600L),
        t1 = c(1, 2, 3),
        stringsAsFactors = FALSE
    )
    mat <- gintervals.to_mat(df)
    out <- gintervals.from_mat(mat[c(3, 1), ])
    expect_equal(out$chrom, c("chr3", "chr1"))
    expect_equal(out$t1, c(3, 1))
})

test_that("rbind of two intervs_mat concatenates attr$intervals", {
    df1 <- data.frame(
        chrom = "chr1", start = 100L, end = 200L, t1 = 1,
        stringsAsFactors = FALSE
    )
    df2 <- data.frame(
        chrom = "chr2", start = 300L, end = 400L, t1 = 2,
        stringsAsFactors = FALSE
    )
    m1 <- gintervals.to_mat(df1)
    m2 <- gintervals.to_mat(df2)
    combined <- rbind(m1, m2)
    expect_s3_class(combined, "intervs_mat")
    expect_equal(nrow(combined), 2L)
    expect_equal(attr(combined, "intervals")$chrom, c("chr1", "chr2"))
    expect_equal(attr(combined, "intervals")$start, c(100L, 300L))

    out <- gintervals.from_mat(combined)
    expect_equal(out$chrom, c("chr1", "chr2"))
    expect_equal(out$t1, c(1, 2))
})

test_that("rbind of intervs_mat with plain matrix drops class (no false attr)", {
    df1 <- data.frame(
        chrom = "chr1", start = 100L, end = 200L, t1 = 1,
        stringsAsFactors = FALSE
    )
    m1 <- gintervals.to_mat(df1)
    plain <- matrix(2, nrow = 1, dimnames = list(NULL, "t1"))
    combined <- rbind(m1, plain)
    expect_false(inherits(combined, "intervs_mat"))
    expect_null(attr(combined, "intervals"))
})

test_that("labels = FALSE leaves rownames NULL", {
    df <- data.frame(
        chrom = c("chr1", "chr2"),
        start = c(100L, 200L),
        end = c(200L, 400L),
        t1 = c(1, 2),
        stringsAsFactors = FALSE
    )
    mat <- gintervals.to_mat(df, labels = FALSE)
    expect_null(rownames(mat))
    expect_equal(attr(mat, "intervals")$chrom, df$chrom)
    # round-trip still works (identity flows through attribute)
    out <- gintervals.from_mat(mat)
    expect_equal(out$chrom, df$chrom)
    expect_equal(out$t1, df$t1)
})

test_that("C++ helper output matches the R paste0 reference exactly", {
    df <- data.frame(
        chrom = c("chr1", "chrUn_KI270442v1", "chrX"),
        start = c(1L, 5000L, 100L),
        end = c(500L, 6000L, 200L),
        t1 = c(1, 2, 3),
        stringsAsFactors = FALSE
    )
    mat <- gintervals.to_mat(df) # default labels = TRUE, no id_col -> C path
    expect_equal(
        rownames(mat),
        paste0(df$chrom, ":", df$start, "-", df$end)
    )
})

test_that("C++ helper handles factor chrom column", {
    df <- data.frame(
        chrom = factor(c("chr1", "chr2", "chr1")),
        start = c(100L, 200L, 300L),
        end   = c(200L, 400L, 500L),
        t1    = c(1, 2, 3)
    )
    mat <- gintervals.to_mat(df)
    expect_equal(
        rownames(mat),
        c("chr1:100-200", "chr2:200-400", "chr1:300-500")
    )
})

test_that("labels = TRUE with id_col still uses the id_col (not the C path)", {
    df <- data.frame(
        chrom = "chr1", start = 100L, end = 200L,
        gene = "FOO", t1 = 1.0,
        stringsAsFactors = FALSE
    )
    mat <- gintervals.to_mat(df, id_col = "gene", value_cols = "t1")
    expect_equal(rownames(mat), "FOO")
})

test_that("explicit value_cols also errors on non-numeric column", {
    df <- data.frame(
        chrom = "chr1", start = 100L, end = 200L,
        gene = "FOO", t1 = 1.0,
        stringsAsFactors = FALSE
    )
    expect_error(gintervals.to_mat(df, value_cols = c("t1", "gene")),
        "Non-numeric value column.*gene",
        fixed = FALSE
    )
})

test_that("id_col is validated even when labels = FALSE", {
    df <- data.frame(
        chrom = "chr1", start = 100L, end = 200L, t1 = 1.0,
        stringsAsFactors = FALSE
    )
    expect_error(gintervals.to_mat(df, id_col = "nope", labels = FALSE),
        "`id_col` not found",
        fixed = FALSE
    )
})

test_that("head() and tail() dispatch through [ and preserve intervs_mat", {
    df <- data.frame(
        chrom = paste0("chr", 1:5),
        start = (1:5) * 100L,
        end = (1:5) * 100L + 50L,
        t1 = 1:5 * 1.0,
        stringsAsFactors = FALSE
    )
    mat <- gintervals.to_mat(df)
    h <- head(mat, 2)
    expect_s3_class(h, "intervs_mat")
    expect_equal(nrow(h), 2L)
    expect_equal(attr(h, "intervals")$chrom, c("chr1", "chr2"))

    t <- tail(mat, 2)
    expect_s3_class(t, "intervs_mat")
    expect_equal(nrow(t), 2L)
    expect_equal(attr(t, "intervals")$chrom, c("chr4", "chr5"))
})
