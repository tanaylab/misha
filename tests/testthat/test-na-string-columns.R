test_that("NA in a character column survives gintervals.rbind and gintervals.neighbors", {
    # DataFrameUtils used to copy character cells with
    # Rf_mkChar(CHAR(STRING_ELT(...))). CHAR(NA_STRING) is the C string "NA", so
    # the round trip turned a missing value into the literal string "NA" and
    # is.na() started reporting FALSE. Copy the CHARSXP directly instead.
    withr::local_options(list(gmultitasking = FALSE))

    a <- gintervals(1, 100, 200)
    a$lab <- NA_character_
    b <- gintervals(1, 300, 400)
    b$lab <- "real"

    r <- gintervals.rbind(a, b)
    expect_true(is.na(r$lab[1]))
    expect_identical(r$lab[2], "real")
    expect_false(identical(r$lab[1], "NA"))

    n <- gintervals.neighbors(gintervals(1, 150, 160), a)
    expect_true(is.na(n$lab[1]))
})
