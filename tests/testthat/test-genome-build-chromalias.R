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

test_that(".parse_ucsc_assembly_report parses 10-column NCBI format and normalizes 'na'", {
    f <- testthat::test_path("fixtures", "assembly-report-mini.txt")
    df <- .parse_ucsc_assembly_report(f)
    expect_equal(nrow(df), 5L)
    expect_true(all(c(
        "Sequence-Name", "GenBank-Accn", "RefSeq-Accn",
        "UCSC-style-name", "Sequence-Length"
    ) %in% names(df)))
    # 'na' placeholders collapse to "".
    expect_equal(df[["Assigned-Molecule"]][[5L]], "na") # not normalized: only GB/RS/UCSC are
    expect_equal(df[["GenBank-Accn"]], c(
        "CM00001.1", "CM00002.1", "CM00003.1",
        "AY999999.1", "JH99999.1"
    ))
    expect_equal(df[["UCSC-style-name"]][[5L]], "chrUn_JH99999v1")
})

test_that(".parse_ucsc_assembly_report strips CRLF carriage returns from fields", {
    # UCSC's mirrored copies use Windows line endings; readLines keeps the
    # trailing '\r' on the last column. Parser must strip it so downstream
    # %in% checks work.
    crlf_path <- tempfile(fileext = ".txt")
    on.exit(unlink(crlf_path))
    writeBin(charToRaw(paste0(
        "# Description: CRLF test\r\n",
        "# Sequence-Name\tGenBank-Accn\tUCSC-style-name\r\n",
        "1\tCM00001.1\tchr1\r\n",
        "MT\tAY99999.1\tchrM\r\n"
    )), crlf_path)
    df <- .parse_ucsc_assembly_report(crlf_path)
    expect_equal(df[["UCSC-style-name"]], c("chr1", "chrM"))
    # Verify no field carries trailing whitespace anywhere.
    expect_false(any(grepl("[[:space:]]$", unlist(df))))
})

test_that(".merge_assembly_report_into_alias keeps report columns independent (no suppression)", {
    # chromAlias has 4 rows; the report has 5 (one extra unplaced scaffold).
    # After merge, alias_df should grow to 5 rows and gain EVERY report
    # column under a normalized name -- including ucsc_style_name, even
    # though alias already has a 'ucsc' column. The two are different
    # naming conventions for unplaced scaffolds and must be scored
    # independently by .detect_alias_column.
    alias_df <- .parse_ucsc_chromalias(
        testthat::test_path("fixtures", "chrom-alias-mini.txt")
    )
    report_df <- .parse_ucsc_assembly_report(
        testthat::test_path("fixtures", "assembly-report-mini.txt")
    )
    merged <- .merge_assembly_report_into_alias(alias_df, report_df)
    expect_equal(nrow(merged), 5L)
    expect_true("ucsc_style_name" %in% names(merged))
    expect_true("sequence_role" %in% names(merged))
    expect_true("sequence_length" %in% names(merged))
    # ucsc_style_name carries the report's per-row UCSC names for ALL rows
    # (including the extra row appended from the report).
    expect_equal(
        merged[["ucsc_style_name"]],
        c("chr1", "chr2", "chr3", "chrM", "chrUn_JH99999v1")
    )
    # The original alias 'ucsc' column is untouched for matched rows and
    # empty for the appended row (no mixing of naming conventions).
    expect_equal(merged$ucsc, c("chr1", "chr2", "chr3", "chrM", ""))
    # Join key (refseq) is filled on the appended row for cross-ref.
    expect_equal(merged$refseq[[5L]], "NW_007777777.1")
})

test_that(".merge_assembly_report_into_alias appends a report-only row even when RefSeq-Accn is 'na'", {
    # Dog (GCF_000002285.5) case: the chromAlias has 147 rows and skips
    # the MT contig; the assembly_report carries MT with RefSeq-Accn='na'
    # and UCSC-style-name='na' but a real GenBank-Accn (CM023446.1) and
    # Sequence-Name='MT'. The merge must append this row using whatever
    # identifier the report supplies, populating the chromAlias-side
    # genbank/assembly columns from the report's equivalents.
    alias_df <- data.frame(
        refseq = c("NC_001.1", "NC_002.1"),
        genbank = c("CM01.1", "CM02.1"),
        ucsc = c("chr1", "chr2"),
        assembly = c("1", "2"),
        stringsAsFactors = FALSE
    )
    report_df <- data.frame(
        check.names = FALSE,
        "Sequence-Name" = c("1", "2", "MT"),
        "GenBank-Accn" = c("CM01.1", "CM02.1", "CM023446.1"),
        "RefSeq-Accn" = c("NC_001.1", "NC_002.1", ""),
        "UCSC-style-name" = c("chr1", "chr2", ""),
        stringsAsFactors = FALSE
    )
    merged <- .merge_assembly_report_into_alias(alias_df, report_df)
    expect_equal(nrow(merged), 3L)
    # Appended row has genbank/assembly filled from report equivalents but
    # ucsc stays empty (no UCSC-style-name to safely copy, and chromAlias's
    # ucsc convention may differ from report's anyway).
    expect_equal(merged$genbank[[3L]], "CM023446.1")
    expect_equal(merged$assembly[[3L]], "MT")
    expect_equal(merged$ucsc[[3L]], "")
    expect_equal(merged$refseq[[3L]], "")
})

test_that(".merge_assembly_report_into_alias preserves both columns when chromAlias and report disagree", {
    # Real-world rat case: chromAlias.ucsc uses RefSeq-derived
    # "chr1_NW_xxx_random"; report's UCSC-style-name uses GenBank-derived
    # "chr1_AABRxxx_random". Same row, same scaffold, different names. The
    # merge must keep BOTH so detect_alias_column can pick whichever
    # matches the groot.
    alias_df <- data.frame(
        refseq = c("NC_001.1", "NW_xxx.1"),
        ucsc = c("chr1", "chr1_NW_xxx_random"),
        genbank = c("CM01.1", "AABR_xxx.1"),
        stringsAsFactors = FALSE
    )
    report_df <- data.frame(
        check.names = FALSE,
        "Sequence-Name" = c("1", "1"),
        "RefSeq-Accn" = c("NC_001.1", "NW_xxx.1"),
        "GenBank-Accn" = c("CM01.1", "AABR_xxx.1"),
        "UCSC-style-name" = c("chr1", "chr1_AABRxxx_random"),
        stringsAsFactors = FALSE
    )
    merged <- .merge_assembly_report_into_alias(alias_df, report_df)
    # Both naming conventions preserved as separate columns.
    expect_equal(merged$ucsc[[2L]], "chr1_NW_xxx_random")
    expect_equal(merged$ucsc_style_name[[2L]], "chr1_AABRxxx_random")
})

test_that(".merge_assembly_report_into_alias preserves alias_df when no join key matches", {
    # alias_df has columns matching nothing in the report.
    alias_df <- data.frame(
        refseq = c("DIFFERENT1", "DIFFERENT2"),
        genbank = c("OTHER1", "OTHER2"),
        stringsAsFactors = FALSE
    )
    report_df <- .parse_ucsc_assembly_report(
        testthat::test_path("fixtures", "assembly-report-mini.txt")
    )
    out <- .merge_assembly_report_into_alias(alias_df, report_df)
    # Returns alias_df unchanged (no columns added, no rows appended).
    expect_equal(out, alias_df)
})

test_that(".merge_assembly_report_into_alias is a no-op on NULL or empty inputs", {
    alias_df <- data.frame(refseq = "NC_1", stringsAsFactors = FALSE)
    expect_null(.merge_assembly_report_into_alias(NULL, data.frame()))
    expect_equal(.merge_assembly_report_into_alias(alias_df, NULL), alias_df)
    expect_equal(.merge_assembly_report_into_alias(alias_df, data.frame()), alias_df)
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

test_that(".name_match_override picks a non-canonical column whose value matches groot", {
    # row 2: canonical="" (UCSC-style-name was 'na') but genbank=AABR07023006.1
    # is in groot under that bare accession. The override picks it up.
    canonical <- c("chr1", "")
    alias_df <- data.frame(
        ucsc_style_name = c("chr1", ""),
        genbank = c("CM01.1", "AABR07023006.1"),
        stringsAsFactors = FALSE
    )
    groot <- c("chr1", "AABR07023006.1")
    out <- .name_match_override(canonical, alias_df, "ucsc_style_name", groot)
    expect_equal(out, c("chr1", "AABR07023006.1"))
})

test_that(".name_match_override does not reuse a groot chrom already in canonical", {
    canonical <- c("chr1", "missing")
    alias_df <- data.frame(
        ucsc_style_name = c("chr1", "missing"),
        other = c("chr1", "chr1"), # would collide with row 1's canonical
        stringsAsFactors = FALSE
    )
    groot <- c("chr1")
    out <- .name_match_override(canonical, alias_df, "ucsc_style_name", groot)
    # row 2 stays "missing" because the only candidate (chr1) is already
    # canonical for row 1.
    expect_equal(out, c("chr1", "missing"))
})

test_that(".diagnose_unmapped_chroms categorizes each unmapped chrom", {
    alias_df <- data.frame(
        ucsc = c("chr1", "chr2"),
        genbank = c("CM01.1", "CM02.1"),
        stringsAsFactors = FALSE
    )
    # 3 unmapped groot chroms: one in another column, one length-ambiguous,
    # one absent.
    unmapped <- c("CM01.1", "amb_a", "ghost")
    groot_chroms <- c("chr1", "chr2", "CM01.1", "amb_a", "amb_b", "ghost")
    groot_lengths <- c(100, 200, 100, 50, 50, 999)
    alias_row_lengths <- c(100, 200) # only 2 alias rows
    out <- .diagnose_unmapped_chroms(
        unmapped, alias_df, "ucsc",
        groot_chroms, groot_lengths, alias_row_lengths, "/path/to/groot"
    )
    expect_match(out, "CM01.1.*present in alias column\\(s\\) genbank")
    expect_match(out, "ghost.*contig absent")
    expect_match(out, "/path/to/groot/chrom_aliases.tsv")
})

test_that(".length_match_override replaces canonical when canonical isn't in groot but length pairs uniquely", {
    # Rat MT scenario: canonical (ucsc_style_name) says "chrM" but the
    # groot has the MT contig under its GenBank accession "AY172581.1".
    # Length 16313 is unique on both sides -> safe to override.
    canonical <- c("chr1", "chrM", "chr2")
    alias_lengths <- c(100, 16313, 200)
    groot_chroms <- c("chr1", "AY172581.1", "chr2")
    groot_lengths <- c(100, 16313, 200)
    out <- .length_match_override(canonical, alias_lengths, groot_chroms, groot_lengths)
    expect_equal(out, c("chr1", "AY172581.1", "chr2"))
})

test_that(".length_match_override leaves canonical alone when length is ambiguous", {
    # Two alias rows at length 642 and four groot chroms at 642 -- can't
    # safely assign. Override must skip.
    canonical <- c("chr1", "chrM", "alias_x")
    alias_lengths <- c(100, 642, 642)
    groot_chroms <- c("chr1", "AABR_a.1", "AABR_b.1", "AABR_c.1", "AABR_d.1")
    groot_lengths <- c(100, 642, 642, 642, 642)
    out <- .length_match_override(canonical, alias_lengths, groot_chroms, groot_lengths)
    expect_equal(out, canonical)
})

test_that(".length_match_override does not reuse a groot chrom already in canonical", {
    # canonical[1] is "chr1" (in groot). canonical[2] is "chrM" (not in
    # groot). Length 100 is unique on both sides. But "chr1" already
    # claims the only groot chrom at length 100. Override must NOT replace
    # canonical[2] with "chr1" -- that would create a duplicate canonical.
    canonical <- c("chr1", "chrM")
    alias_lengths <- c(50, 100)
    groot_chroms <- c("chr1", "AY999.1")
    groot_lengths <- c(50, 100) # "chr1" already in canonical at length 50, AY999.1 at 100
    out <- .length_match_override(canonical, alias_lengths, groot_chroms, groot_lengths)
    # canonical[2] gets replaced with "AY999.1" (not yet in canonical).
    expect_equal(out, c("chr1", "AY999.1"))
})

test_that(".length_match_override is a no-op when canonical already matches groot", {
    canonical <- c("chr1", "chr2")
    alias_lengths <- c(100, 200)
    groot_chroms <- c("chr1", "chr2")
    groot_lengths <- c(100, 200)
    expect_equal(.length_match_override(
        canonical, alias_lengths,
        groot_chroms, groot_lengths
    ), canonical)
})

test_that(".length_match_fill is a no-op when no fills needed", {
    canonical <- c("X", "Y")
    out <- .length_match_fill(canonical, c(100, 200), c("X", "Y"), c(100, 200))
    expect_equal(out, c("X", "Y"))
})

test_that(".assign_target_chroms_per_row places every target via name match across columns", {
    # target_chroms appear scattered across columns (mixing namings is fine).
    df <- data.frame(
        ucsc = c("chr1", "chr2", "chr3", "chrM"),
        genbank = c("CM1", "CM2", "CM3", ""),
        refseq = c("NC1", "NC2", "NC3", "NCM"),
        stringsAsFactors = FALSE
    )
    target_chroms <- c("chr1", "CM2", "NC3", "chrM")
    target_lengths <- c(100, 200, 300, 16)
    alias_row_lengths <- c(100, 200, 300, 16)
    out <- .assign_target_chroms_per_row(
        df, target_chroms, target_lengths, alias_row_lengths
    )
    expect_equal(out, c("chr1", "CM2", "NC3", "chrM"))
})

test_that(".assign_target_chroms_per_row falls back to length pairing for unmatched targets", {
    # 'HALMT' isn't in any chromAlias column but its length pairs uniquely
    # with the empty row's alias length -> length fill places it.
    df <- data.frame(
        ucsc = c("chr1", "chr2", "chr3", ""),
        stringsAsFactors = FALSE
    )
    target_chroms <- c("chr1", "chr2", "chr3", "HALMT")
    target_lengths <- c(100, 200, 300, 16)
    alias_row_lengths <- c(100, 200, 300, 16)
    out <- .assign_target_chroms_per_row(
        df, target_chroms, target_lengths, alias_row_lengths
    )
    expect_equal(out, c("chr1", "chr2", "chr3", "HALMT"))
})

test_that(".assign_target_chroms_per_row errors when a target can't be placed", {
    # 'HALMT' isn't in chromAlias and its length is ambiguous (length 16
    # appears twice on alias side) -> unplaced -> stop.
    df <- data.frame(
        ucsc = c("chr1", "chr2", "", ""),
        stringsAsFactors = FALSE
    )
    target_chroms <- c("chr1", "chr2", "HALMT")
    target_lengths <- c(100, 200, 16)
    alias_row_lengths <- c(100, 200, 16, 16)
    expect_error(
        .assign_target_chroms_per_row(
            df, target_chroms, target_lengths, alias_row_lengths
        ),
        "HALMT"
    )
})

test_that(".assign_target_chroms_per_row leaves rows not in target_chroms empty", {
    # Extra alias row (smallContig with length 17) has no matching target;
    # it must come back empty (caller fills with src_col fallback if needed).
    df <- data.frame(
        ucsc = c("chr1", "chr2", "smallContig"),
        stringsAsFactors = FALSE
    )
    target_chroms <- c("chr1", "chr2")
    target_lengths <- c(100, 200)
    alias_row_lengths <- c(100, 200, 17)
    out <- .assign_target_chroms_per_row(
        df, target_chroms, target_lengths, alias_row_lengths
    )
    expect_equal(out, c("chr1", "chr2", ""))
})

test_that(".assign_target_chroms_per_row errors on input misalignment", {
    df <- data.frame(a = c("x", "y"), stringsAsFactors = FALSE)
    expect_error(
        .assign_target_chroms_per_row(
            df,
            target_chroms = c("X", "Y"),
            target_lengths = c(1, 2, 3),
            alias_row_lengths = c(1, 2)
        ),
        "target_lengths"
    )
    expect_error(
        .assign_target_chroms_per_row(
            df,
            target_chroms = c("X", "Y"),
            target_lengths = c(1, 2),
            alias_row_lengths = c(1, 2, 3)
        ),
        "alias_row_lengths"
    )
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

test_that(".merge_chrom_aliases_tsv tolerates pre-existing wide-format TSV (NCBI seq-only)", {
    # .build_seq_ncbi (pre-v5.6.30) wrote chrom_aliases.tsv in wide format:
    # (canonical, refseqAccession, genbankAccession, sequenceName, chrName, role, length).
    # When gdb.install_intervals runs against such a groot, the merge must
    # convert the wide rows to long and append the new alias rows.
    f <- testthat::test_path("fixtures", "chrom-alias-mini.txt")
    alias_df <- .parse_ucsc_chromalias(f)
    groot <- tempfile()
    dir.create(groot)
    on.exit(unlink(groot, recursive = TRUE))
    wide <- data.frame(
        canonical        = "chrZ",
        refseqAccession  = "NC_999999.1",
        genbankAccession = "CM999999.1",
        sequenceName     = "Z",
        chrName          = "Z",
        role             = "assembled-molecule",
        length           = 1000L,
        stringsAsFactors = FALSE
    )
    write.table(wide, file.path(groot, "chrom_aliases.tsv"),
        sep = "\t", quote = FALSE, row.names = FALSE
    )
    .merge_chrom_aliases_tsv(groot, alias_df, groot_col = "ucsc")
    out <- read.table(file.path(groot, "chrom_aliases.tsv"),
        sep = "\t", header = TRUE, stringsAsFactors = FALSE
    )
    # Existing wide row produces 4 long rows (chrZ -> NC_999999.1, CM999999.1, "Z", "Z").
    # The chrName/sequenceName rows dedupe to one alias ("Z") per (canonical, source).
    expect_true("NC_999999.1" %in% out$alias)
    expect_true("CM999999.1" %in% out$alias)
    expect_true(all(out$canonical[out$alias %in% c("NC_999999.1", "CM999999.1")] == "chrZ"))
    # New long-format rows from the merge are also present.
    expect_true(any(out$canonical == "chr1"))
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

test_that(".canonical_coverage reports bp-weighted coverage of canonical vs groot", {
    canonical <- c("chr1", "chr2", "chrX", "")
    groot_chroms <- c("chr1", "chr2", "chrX")
    groot_lengths <- c(1000, 2000, 500)
    expect_equal(
        .canonical_coverage(canonical, groot_chroms, groot_lengths),
        1.0
    )
})

test_that(".canonical_coverage handles partial coverage", {
    canonical <- c("chr1", "chrX") # chr2 missing entirely
    groot_chroms <- c("chr1", "chr2", "chrX")
    groot_lengths <- c(1000, 2000, 500)
    expect_equal(
        .canonical_coverage(canonical, groot_chroms, groot_lengths),
        1500 / 3500
    )
})

test_that(".canonical_coverage is 0 when nothing matches", {
    expect_equal(
        .canonical_coverage(c("foo", "bar"), c("chr1", "chr2"), c(100, 200)),
        0
    )
})

test_that(".coverage_gate returns winning column when threshold met", {
    df <- data.frame(
        ucsc = c("chr1", "chr2", "chrM"),
        genbank = c("CM00001.1", "CM00002.1", ""),
        stringsAsFactors = FALSE
    )
    target <- c("chr1", "chr2", "chrM")
    lengths <- c(1e8, 1e8, 1.6e4)
    col <- .coverage_gate(df, target, lengths,
        min_coverage = 0.99,
        label = "groot"
    )
    expect_equal(as.character(col), "ucsc")
})

test_that(".coverage_gate stops with diagnostic when no column meets threshold", {
    df <- data.frame(
        ucsc = c("chr1", "chr2", ""),
        genbank = c("CM00001.1", "", ""),
        stringsAsFactors = FALSE
    )
    target <- c("chr1", "chr2", "chrM")
    lengths <- c(1e8, 1e8, 1.6e4)
    expect_error(
        .coverage_gate(df, target, lengths, min_coverage = 1.0, label = "groot"),
        "no column with 100% bp coverage of groot"
    )
})

test_that(".coverage_gate is a no-op when alias_df is NULL", {
    expect_null(.coverage_gate(NULL, "chr1", 1, min_coverage = 1.0))
})
