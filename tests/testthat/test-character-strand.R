create_isolated_test_db()

test_that("gintervals accepts +/- character strand", {
    num <- gintervals(c(1, 1, 1), c(0, 100, 200), c(50, 150, 250), c(1, -1, 0))
    chr <- gintervals(c(1, 1, 1), c(0, 100, 200), c(50, 150, 250), c("+", "-", "."))
    expect_equal(num, chr)
})

test_that("gintervals accepts factor strand", {
    num <- gintervals(c(1, 1, 1), c(0, 100, 200), c(50, 150, 250), c(1, -1, 0))
    fac <- gintervals(
        c(1, 1, 1), c(0, 100, 200), c(50, 150, 250),
        factor(c("+", "-", "."), levels = c("+", "-", ".", "*"))
    )
    expect_equal(num, fac)
})

test_that("gintervals accepts unknown-strand variants .,*,'' as 0", {
    res <- gintervals(c(1, 1, 1), c(0, 100, 200), c(50, 150, 250), c(".", "*", ""))
    expect_equal(res$strand, c(0, 0, 0))
})

test_that("gintervals rejects garbage strand strings", {
    expect_error(
        gintervals(c(1, 1), c(0, 100), c(50, 150), c("+", "foo")),
        "Invalid strand value"
    )
})

test_that("gintervals rejects strand integer outside {-1,0,1}", {
    expect_error(
        gintervals(c(1, 1), c(0, 100), c(50, 150), c(1, 2)),
        "Invalid strand value"
    )
})

test_that("data.frame with character strand flows through C++ converter", {
    df_num <- data.frame(
        chrom = "chr1", start = c(0L, 1000L, 2000L), end = c(500L, 1500L, 2500L),
        strand = c(1, -1, 0)
    )
    df_chr <- df_num
    df_chr$strand <- c("+", "-", ".")

    # gintervals.canonic flows through convert_rintervs in C++; both inputs must
    # produce identical output once the strand column has been normalized.
    out_num <- gintervals.canonic(df_num)
    out_chr <- gintervals.canonic(df_chr)
    expect_equal(out_num, out_chr)
})

test_that("gintervals.neighbors accepts character strand", {
    intervs1 <- data.frame(
        chrom = "chr1", start = c(1000L, 5000L), end = c(1100L, 5100L),
        strand = c("+", "-"), stringsAsFactors = FALSE
    )
    intervs2 <- data.frame(
        chrom = "chr1", start = c(2000L, 4000L), end = c(2100L, 4100L),
        strand = c("+", "-"), stringsAsFactors = FALSE
    )
    expect_silent(suppressWarnings(gintervals.neighbors(intervs1, intervs2, 1)))
})

test_that("data.frame with character strand reports invalid value via C++", {
    df <- data.frame(
        chrom = "chr1", start = c(0L, 1000L), end = c(500L, 1500L),
        strand = c("+", "wat"), stringsAsFactors = FALSE
    )
    expect_error(gintervals.canonic(df), "Bad strand")
})
