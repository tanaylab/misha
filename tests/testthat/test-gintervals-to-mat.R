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
