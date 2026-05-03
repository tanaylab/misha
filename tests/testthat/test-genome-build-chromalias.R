test_that(".parse_ucsc_chromalias reads header with leading # and returns one column per scheme", {
    f <- testthat::test_path("fixtures", "chrom-alias-mini.txt")
    df <- .parse_ucsc_chromalias(f)
    expect_named(df, c("ucsc", "assembly", "genbank", "refseq"))
    expect_equal(nrow(df), 4L)
    expect_equal(df$ucsc, c("chr1", "chr2", "chr3", "chrM"))
    expect_equal(df$refseq, c("NC_067374.1", "NC_067375.1", "NC_067376.1", "NC_067377.1"))
})

test_that(".parse_ucsc_chromalias handles gzipped input", {
    f <- testthat::test_path("fixtures", "chrom-alias-mini.txt")
    gz <- tempfile(fileext = ".txt.gz")
    on.exit(unlink(gz))
    src <- readBin(f, raw(), n = file.info(f)$size)
    con <- gzfile(gz, "wb")
    writeBin(src, con)
    close(con)
    df <- .parse_ucsc_chromalias(gz)
    expect_equal(nrow(df), 4L)
    expect_equal(df$ucsc[[1L]], "chr1")
})

test_that(".detect_alias_column returns the column with 100% coverage", {
    f <- testthat::test_path("fixtures", "chrom-alias-mini.txt")
    df <- .parse_ucsc_chromalias(f)
    expect_equal(.detect_alias_column(df, c("chr1", "chr2", "chr3", "chrM")), "ucsc")
    expect_equal(.detect_alias_column(df, c("1", "2", "3", "MT")), "assembly")
    expect_equal(.detect_alias_column(
        df,
        c("NC_067374.1", "NC_067375.1", "NC_067376.1", "NC_067377.1")
    ), "refseq")
})

test_that(".detect_alias_column returns NA when no column has 100% coverage", {
    f <- testthat::test_path("fixtures", "chrom-alias-mini.txt")
    df <- .parse_ucsc_chromalias(f)
    res <- .detect_alias_column(df, c("chr1", "chr2", "chr3", "chrM", "chrX")) # chrX absent
    expect_true(is.na(res))
    expect_true(!is.null(attr(res, "scores")))
    expect_equal(attr(res, "scores")[["ucsc"]], 4L) # 4 of 5 hits
})

test_that(".detect_alias_column ties broken by column order (first wins)", {
    df <- data.frame(a = c("x", "y"), b = c("x", "y"), stringsAsFactors = FALSE)
    expect_equal(.detect_alias_column(df, c("x", "y")), "a")
})

test_that(".detect_alias_column on empty target errors", {
    f <- testthat::test_path("fixtures", "chrom-alias-mini.txt")
    df <- .parse_ucsc_chromalias(f)
    expect_error(.detect_alias_column(df, character(0)), "empty target")
})

test_that(".translate_chroms is a no-op when source_col == groot_col", {
    f <- testthat::test_path("fixtures", "chrom-alias-mini.txt")
    alias_df <- .parse_ucsc_chromalias(f)
    rows <- data.frame(chrom = c("chr1", "chr2"), x = 1:2, stringsAsFactors = FALSE)
    out <- .translate_chroms(rows, "chrom", alias_df, "ucsc", "ucsc")
    expect_identical(out, rows)
})

test_that(".translate_chroms maps via alias table when columns differ", {
    f <- testthat::test_path("fixtures", "chrom-alias-mini.txt")
    alias_df <- .parse_ucsc_chromalias(f)
    rows <- data.frame(
        chrom = c("NC_067374.1", "NC_067375.1"), x = 1:2,
        stringsAsFactors = FALSE
    )
    out <- .translate_chroms(rows, "chrom", alias_df, "refseq", "ucsc")
    expect_equal(out$chrom, c("chr1", "chr2"))
    expect_equal(out$x, 1:2)
})
