create_isolated_test_db()
test_that("import with attrs parameter - error on unnamed attrs", {
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

    # Test error with unnamed vector
    attrs <- c("test_user", "1.0")
    expect_error(
        gtrack.import(tmptrack, "Test track", temp_file, binsize = 1, attrs = attrs),
        "attrs must be a named vector or list"
    )
})

test_that("import with attrs parameter - error on partially unnamed attrs", {
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

    # Test error with partially named vector
    attrs <- c("author" = "test_user", "1.0") # second element has no name
    expect_error(
        gtrack.import(tmptrack, "Test track", temp_file, binsize = 1, attrs = attrs),
        "attrs must be a named vector or list"
    )
})
