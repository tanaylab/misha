create_isolated_test_db()

test_that("import with gmax data size option", {
    tmptrack <- paste0("test.tmptrack_", sample(1:1e9, 1))
    gtrack.rm(tmptrack, TRUE)
    withr::defer(gtrack.rm(tmptrack, force = TRUE))
    withr::with_options(list(gmax.data.size = 10000), {
        gtrack.2d.import(tmptrack, "aaa7", c("/net/mraid20/export/tgdata/db/tgdb/misha_snapshot/input_files/f4"))
    })
    r <- gextract(tmptrack, .misha$ALLGENOME, colnames = "test.tmptrack")
    expect_regression(r, "track.import_gmax_option")
})

test_that("import with attrs parameter - single attribute", {
    tmptrack <- paste0("test.tmptrack_", sample(1:1e9, 1))
    gtrack.rm(tmptrack, force = TRUE)
    withr::defer(gtrack.rm(tmptrack, force = TRUE))

    # Create a temporary file for testing
    temp_file <- tempfile(fileext = ".wig")
    writeLines(c(
        "track type=wiggle_0 name=\"test track\"",
        "fixedStep chrom=chr1 start=1 step=1",
        "1.0",
        "2.0",
        "3.0"
    ), temp_file)
    withr::defer(unlink(temp_file))

    # Import with single attribute
    attrs <- c("author" = "test_user")
    gtrack.import(tmptrack, "Test track", temp_file, binsize = 1, attrs = attrs)

    # Verify the attribute was set
    expect_equal(gtrack.attr.get(tmptrack, "author"), "test_user")
})
