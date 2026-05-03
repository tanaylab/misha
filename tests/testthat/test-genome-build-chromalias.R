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

test_that(".merge_chrom_aliases_tsv writes one row per (canonical, alias) when no existing file", {
    f <- testthat::test_path("fixtures", "chrom-alias-mini.txt")
    alias_df <- .parse_ucsc_chromalias(f)
    groot <- tempfile()
    dir.create(groot)
    on.exit(unlink(groot, recursive = TRUE))
    .merge_chrom_aliases_tsv(groot, alias_df, groot_col = "ucsc")
    out <- read.table(file.path(groot, "chrom_aliases.tsv"),
        sep = "\t",
        header = TRUE, stringsAsFactors = FALSE
    )
    # 4 contigs × 3 alias columns (assembly/genbank/refseq) = 12 rows.
    expect_equal(nrow(out), 12L)
    expect_setequal(out$canonical, c("chr1", "chr2", "chr3", "chrM"))
    expect_setequal(unique(out$source), c("assembly", "genbank", "refseq"))
})

test_that(".merge_chrom_aliases_tsv preserves existing rows and dedupes", {
    f <- testthat::test_path("fixtures", "chrom-alias-mini.txt")
    alias_df <- .parse_ucsc_chromalias(f)
    groot <- tempfile()
    dir.create(groot)
    on.exit(unlink(groot, recursive = TRUE))
    pre <- data.frame(
        canonical = "chr1", alias = "old_alias", source = "manual",
        stringsAsFactors = FALSE
    )
    write.table(pre, file.path(groot, "chrom_aliases.tsv"),
        sep = "\t", quote = FALSE, row.names = FALSE
    )
    .merge_chrom_aliases_tsv(groot, alias_df, groot_col = "ucsc")
    out <- read.table(file.path(groot, "chrom_aliases.tsv"),
        sep = "\t",
        header = TRUE, stringsAsFactors = FALSE
    )
    expect_true(any(out$alias == "old_alias" & out$canonical == "chr1"))
    expect_equal(nrow(out), 13L) # 12 new + 1 preserved
})

test_that(".merge_chrom_aliases_tsv warns on canonical conflict and keeps existing", {
    f <- testthat::test_path("fixtures", "chrom-alias-mini.txt")
    alias_df <- .parse_ucsc_chromalias(f)
    groot <- tempfile()
    dir.create(groot)
    on.exit(unlink(groot, recursive = TRUE))
    # Pre-existing maps NC_067374.1 -> chrFOO (conflicts with new chr1).
    pre <- data.frame(
        canonical = "chrFOO", alias = "NC_067374.1",
        source = "manual", stringsAsFactors = FALSE
    )
    write.table(pre, file.path(groot, "chrom_aliases.tsv"),
        sep = "\t", quote = FALSE, row.names = FALSE
    )
    expect_warning(
        .merge_chrom_aliases_tsv(groot, alias_df, groot_col = "ucsc"),
        "conflict"
    )
    out <- read.table(file.path(groot, "chrom_aliases.tsv"),
        sep = "\t",
        header = TRUE, stringsAsFactors = FALSE
    )
    # Existing chrFOO mapping preserved, new (chr1, NC_067374.1) skipped.
    expect_equal(out$canonical[out$alias == "NC_067374.1"], "chrFOO")
})
