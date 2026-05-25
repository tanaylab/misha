test_that("gdb.unload clears the session state and a re-init restores it", {
    skip_if_not(exists("GROOT", envir = .misha, inherits = FALSE))
    root <- get("GROOT", envir = .misha)
    # Always restore the session so the rest of the suite keeps working.
    on.exit(
        {
            gsetroot(root)
            gdb.reload()
        },
        add = TRUE
    )

    gdb.unload()

    expect_false(exists("GROOT", envir = .misha, inherits = FALSE))
    expect_false(exists("ALLGENOME", envir = .misha, inherits = FALSE))
    # Package internals must survive an unload.
    expect_true(exists(".GLIBDIR", envir = .misha, inherits = FALSE))
    # Any genomic operation must now fail until gdb.init() is called.
    expect_error(gintervals.all(), "not set")
})

test_that("gdb.unload followed by gsetroot yields a working session", {
    skip_if_not(exists("GROOT", envir = .misha, inherits = FALSE))
    root <- get("GROOT", envir = .misha)
    on.exit(
        {
            gsetroot(root)
            gdb.reload()
        },
        add = TRUE
    )

    gdb.unload()
    gsetroot(root)
    expect_true(exists("GROOT", envir = .misha, inherits = FALSE))
    expect_gt(nrow(gintervals.all()), 0)
})
