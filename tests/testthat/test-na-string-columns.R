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

test_that("a track expression that parses to NULL errors instead of crashing a worker", {
    # m_eval_exprs uses R_NilValue as "this expression needs no R evaluation",
    # so an expression that PARSES to NULL looked like that sentinel, was never
    # evaluated, and the first read of its buffer segfaulted. Under multitasking
    # that surfaced as "Child process ended unexpectedly (signal 11)" - the
    # failure shape that can take the caller's tempdir and database with it.
    expect_error(
        gextract("NULL", gintervals(1, 0, 1000), iterator = 100),
        "evaluates to NULL"
    )

    withr::local_options(list(gmultitasking = TRUE, gmax.processes = 4, gmin.scope4process = 1))
    expect_error(
        gextract("NULL", gintervals(1, 0, 1000), iterator = 100),
        "evaluates to NULL"
    )
    # the session is still usable afterwards
    expect_equal(nrow(gextract("dense_track", gintervals(1, 0, 1000), iterator = 100)), 10)
})
