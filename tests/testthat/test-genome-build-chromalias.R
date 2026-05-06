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
    expect_equal(.detect_alias_column(df, c("chr1", "chr2", "chr3", "chrM")), "ucsc",
        ignore_attr = TRUE
    )
    expect_equal(.detect_alias_column(df, c("1", "2", "3", "MT")), "assembly",
        ignore_attr = TRUE
    )
    expect_equal(
        .detect_alias_column(
            df,
            c("NC_067374.1", "NC_067375.1", "NC_067376.1", "NC_067377.1")
        ),
        "refseq",
        ignore_attr = TRUE
    )
})

test_that(".detect_alias_column with min_coverage < 1.0 picks highest partial match", {
    f <- testthat::test_path("fixtures", "chrom-alias-mini.txt")
    df <- .parse_ucsc_chromalias(f)
    # 4-of-5 in ucsc column (chrX absent); strict 100% gives NA, 0.5 threshold
    # picks ucsc.
    res <- .detect_alias_column(df,
        c("chr1", "chr2", "chr3", "chrM", "chrX"),
        min_coverage = 0.75
    )
    expect_equal(as.character(res), "ucsc")
    expect_equal(attr(res, "overlap"), 4 / 5)
    expect_equal(attr(res, "scores")[["ucsc"]], 4L)
})

test_that(".detect_alias_column with min_coverage = 1.0 stays strict", {
    f <- testthat::test_path("fixtures", "chrom-alias-mini.txt")
    df <- .parse_ucsc_chromalias(f)
    res <- .detect_alias_column(df,
        c("chr1", "chr2", "chr3", "chrM", "chrX"),
        min_coverage = 1.0
    )
    expect_true(is.na(res))
})

test_that(".detect_alias_column returns NA when no column has 100% coverage", {
    f <- testthat::test_path("fixtures", "chrom-alias-mini.txt")
    df <- .parse_ucsc_chromalias(f)
    res <- .detect_alias_column(df, c("chr1", "chr2", "chr3", "chrM", "chrX")) # chrX absent
    expect_true(is.na(res))
    expect_true(!is.null(attr(res, "scores")))
    expect_equal(attr(res, "scores")[["ucsc"]], 4L) # 4 of 5 hits
})

test_that(".detect_alias_column with chrom_lengths returns bp-weighted coverage", {
    # Bos-mutus-style scenario: groot has been renamed to genbank names where
    # available, but the 16 kb MT contig kept its refseq id because genbank
    # was empty for that row. So groot chroms span both name spaces. Each
    # column then covers exactly 1/2 by count (50%); by bp the genbank column
    # covers ~99.9984% (1 Gb of (1 Gb + 16 kb)) while refseq covers ~0.0016%.
    df <- data.frame(
        refseq = c("NC_006380.3", "NW_005392810.1"),
        genbank = c("", "JH880237.1"),
        stringsAsFactors = FALSE
    )
    target <- c("NC_006380.3", "JH880237.1")
    lengths <- c(16000, 1e9)

    # Count mode: both columns at 50%, both fail strict 0.99.
    res_count <- .detect_alias_column(df, target, min_coverage = 0.99)
    expect_true(is.na(res_count))

    # Bp mode: genbank at ~99.9984% passes, refseq at ~0.0016% fails.
    res_bp <- .detect_alias_column(df, target,
        min_coverage = 0.99, chrom_lengths = lengths
    )
    expect_equal(as.character(res_bp), "genbank")
    expect_equal(unname(attr(res_bp, "scores")[["genbank"]]), 1e9)
    expect_true(attr(res_bp, "overlap") > 0.99)
    expect_true(attr(res_bp, "bp_weighted"))
})

test_that(".detect_alias_column with chrom_lengths still rejects gross mismatches", {
    df <- data.frame(
        a = c("x", "y"),
        b = c("p", "q"),
        stringsAsFactors = FALSE
    )
    res <- .detect_alias_column(df, c("x", "y"),
        min_coverage = 0.99, chrom_lengths = c(1, 1)
    )
    expect_equal(as.character(res), "a")

    res2 <- .detect_alias_column(df, c("x", "z"),
        min_coverage = 0.99, chrom_lengths = c(1, 1)
    )
    # 'z' missing from both columns -> 50% bp coverage, fails 0.99.
    expect_true(is.na(res2))
})

test_that(".detect_alias_column errors when chrom_lengths length mismatches target", {
    df <- data.frame(a = c("x", "y"), stringsAsFactors = FALSE)
    expect_error(
        .detect_alias_column(df, c("x", "y"), chrom_lengths = c(1, 2, 3)),
        "chrom_lengths"
    )
})

test_that(".detect_alias_column ties broken by column order (first wins)", {
    df <- data.frame(a = c("x", "y"), b = c("x", "y"), stringsAsFactors = FALSE)
    expect_equal(.detect_alias_column(df, c("x", "y")), "a", ignore_attr = TRUE)
})

test_that(".detect_alias_column on empty target errors", {
    f <- testthat::test_path("fixtures", "chrom-alias-mini.txt")
    df <- .parse_ucsc_chromalias(f)
    expect_error(.detect_alias_column(df, character(0)), "empty target")
})

test_that(".length_match_fill maps empty alias rows to unique-length groot chroms", {
    # Two alias rows: row 1 has full canonical, row 2 has empty canonical and
    # needs length-fill. Lengths must be unique on both sides for a match.
    canonical <- c("JH880237.1", "")
    alias_lengths <- c(1e9, 16000)
    groot_chroms <- c("JH880237.1", "NC_006380.3")
    groot_lengths <- c(1e9, 16000)
    out <- .length_match_fill(canonical, alias_lengths, groot_chroms, groot_lengths)
    expect_equal(out, c("JH880237.1", "NC_006380.3"))
})

test_that(".length_match_fill leaves rows NA when length is ambiguous", {
    canonical <- c("", "", "")
    alias_lengths <- c(100, 100, 200)
    groot_chroms <- c("a", "b", "c")
    groot_lengths <- c(100, 100, 200)
    # Length 100 appears twice on both sides -> ambiguous; length 200 unique.
    out <- .length_match_fill(canonical, alias_lengths, groot_chroms, groot_lengths)
    expect_true(is.na(out[[1L]]) || !nzchar(out[[1L]]))
    expect_true(is.na(out[[2L]]) || !nzchar(out[[2L]]))
    expect_equal(out[[3L]], "c")
})

test_that(".length_match_fill is a no-op when no fills needed", {
    canonical <- c("X", "Y")
    out <- .length_match_fill(canonical, c(100, 200), c("X", "Y"), c(100, 200))
    expect_equal(out, c("X", "Y"))
})

test_that(".translate_chroms_per_row resolves any naming via cross-column lookup", {
    # alias_df has 3 cols + per-row canonical. A GFF mixing refseq (NW_X.1)
    # and a different namespace (chrM) should resolve to the right groot
    # canonical regardless.
    alias_df <- data.frame(
        refseq = c("NW_X.1", "NC_006380.3"),
        genbank = c("JH_X.1", ""),
        ucsc = c("NW_Xv1", "chrM"),
        canonical = c("JH_X.1", "NC_006380.3"),
        stringsAsFactors = FALSE
    )
    rows <- data.frame(
        chrom = c("NW_X.1", "chrM", "unknown"), x = 1:3,
        stringsAsFactors = FALSE
    )
    rev_idx <- .build_alias_rev_index(alias_df, "canonical")
    out <- .translate_chroms_per_row(rows, "chrom", rev_idx)
    expect_equal(out$chrom, c("JH_X.1", "NC_006380.3", NA_character_))
    expect_equal(out$x, 1:3)
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

test_that(".resolve_hub_target_col maps known friendly names", {
    cols <- c("ucsc", "assembly", "genbank", "refseq", "ncbi")
    expect_equal(.resolve_hub_target_col("ucsc", src_col = "refseq", alias_cols = cols), "ucsc")
    expect_equal(.resolve_hub_target_col("sequence_name", src_col = "refseq", alias_cols = cols), "assembly")
    # 'accession' = keep src_col (no rename).
    expect_equal(.resolve_hub_target_col("accession", src_col = "refseq", alias_cols = cols), "refseq")
})

test_that(".resolve_hub_target_col accepts arbitrary column names like 'genbank'", {
    cols <- c("ucsc", "assembly", "genbank", "refseq", "ncbi")
    expect_equal(.resolve_hub_target_col("genbank", src_col = "refseq", alias_cols = cols), "genbank")
    expect_equal(.resolve_hub_target_col("refseq", src_col = "ucsc", alias_cols = cols), "refseq")
})

test_that(".resolve_hub_target_col returns NA with a reason when column absent", {
    cols <- c("ucsc", "assembly", "refseq")
    res <- .resolve_hub_target_col("genbank", src_col = "refseq", alias_cols = cols)
    expect_true(is.na(res))
    expect_match(attr(res, "reason"), "no 'genbank' column")
    expect_match(attr(res, "reason"), "ucsc, assembly, refseq")
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
