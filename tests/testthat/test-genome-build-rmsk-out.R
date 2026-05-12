test_that(".parse_rm_out reads RepeatMasker .out skipping 3 header lines", {
    f <- testthat::test_path("fixtures", "rmsk-mini.out.gz")
    df <- .parse_rm_out(f, verbose = FALSE)
    expect_true(nrow(df) > 0L)
    expect_named(df, c("chrom", "start", "end", "strand", "name", "class", "family"))
    # Strand: + or C in source -> 1 or -1 in output.
    expect_true(all(df$strand %in% c(1L, -1L)))
    # 1-based RM coords -> 0-based misha coords: starts shifted down by 1.
    expect_true(all(df$end >= df$start))
    expect_true(all(df$start >= 0L))
})

test_that(".parse_rm_out splits class/family on /", {
    f <- testthat::test_path("fixtures", "rmsk-mini.out.gz")
    df <- .parse_rm_out(f, verbose = FALSE)
    # Some rows have class/family (e.g. SINE/Alu); some have just class.
    expect_true(any(!is.na(df$family)))
})

test_that(".parse_rm_out drops rows with < 11 fields silently with count", {
    # Build a synthetic .out with two header lines + one valid row + one short row.
    tmp <- tempfile(fileext = ".out")
    on.exit(unlink(tmp))
    writeLines(c(
        "   SW   perc perc perc  query  pos pos pos repeat class begin end left ID  ",
        "score   div del ins  sequence start end (left) repeat            position",
        "",
        "1234 0.0 0.0 0.0 chr1 100 200 (1000) + L1   LINE/L1  1 100 (0) 1",
        "short row"
    ), tmp)
    expect_warning(df <- .parse_rm_out(tmp, verbose = TRUE), "dropped 1 malformed")
    expect_equal(nrow(df), 1L)
})
