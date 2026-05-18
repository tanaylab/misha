create_isolated_test_db()
test_that("import with attrs parameter - multiple attributes", {
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

    # Import with multiple attributes
    attrs <- c("author" = "test_user", "version" = "1.0", "experiment" = "test_exp")
    gtrack.import(tmptrack, "Test track", temp_file, binsize = 1, attrs = attrs)

    # Verify all attributes were set
    expect_equal(gtrack.attr.get(tmptrack, "author"), "test_user")
    expect_equal(gtrack.attr.get(tmptrack, "version"), "1.0")
    expect_equal(gtrack.attr.get(tmptrack, "experiment"), "test_exp")
})

test_that("import with attrs parameter - list format", {
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

    # Import with attributes as named list
    attrs <- list("author" = "test_user", "version" = "2.0")
    gtrack.import(tmptrack, "Test track", temp_file, binsize = 1, attrs = attrs)

    # Verify attributes were set
    expect_equal(gtrack.attr.get(tmptrack, "author"), "test_user")
    expect_equal(gtrack.attr.get(tmptrack, "version"), "2.0")
})

test_that("malformed WIG error message includes underlying parser diagnostic", {
    tmptrack <- paste0("test.tmptrack_", sample(1:1e9, 1))
    gtrack.rm(tmptrack, force = TRUE)
    withr::defer(gtrack.rm(tmptrack, force = TRUE))

    # Simulate the concatenation-without-newline corruption: a value line glued
    # to a fixedStep header (e.g. "0.500fixedStep chrom=...").
    temp_file <- tempfile(fileext = ".wig")
    writeLines(c(
        "fixedStep chrom=chr1 start=1 step=1",
        "1.0",
        "2.0",
        "0.500fixedStep chrom=chr2 start=1 step=1",
        "3.0"
    ), temp_file)
    withr::defer(unlink(temp_file))

    # The error must surface the underlying WIG parser message (file + line)
    # rather than the generic "Unrecognized format" that hid the root cause.
    expect_error(
        gtrack.import(tmptrack, "", temp_file, binsize = 1),
        regexp = "WIG parser error.*Invalid format.*line"
    )
})
