test_that("test for dense track from intervals and values", {
    gtrack.rm("test.tmptrack", force = TRUE)
    withr::defer(gtrack.rm("test.tmptrack", force = TRUE))

    # Create intervals and values
    intervals <- gintervals(
        chrom = c(1, 1, 2),
        start = c(100, 500, 200),
        end = c(200, 600, 300)
    )
    values <- c(1.5, 2.7, 3.1)

    gtrack.create_dense("test.tmptrack", "Test dense track", intervals, values, 50, 0)

    # Extract data and check results
    r <- gextract("test.tmptrack", gintervals(c(1, 2), 0, 1000), colnames = "value")

    # Verify data is present
    expect_true(nrow(r) > 0)
    expect_true(any(r$value > 0))

    # Check track info - note that type is lowercase "dense" not "Dense"
    info <- gtrack.info("test.tmptrack")
    expect_equal(info$type, "dense")
    expect_equal(info$bin.size, 50)
})

test_that("test for overlapping intervals in dense track", {
    gtrack.rm("test.tmptrack", force = TRUE)
    withr::defer(gtrack.rm("test.tmptrack", force = TRUE))

    # Create intervals with overlapping regions
    intervals <- gintervals(
        chrom = c(1, 1, 1),
        start = c(100, 150, 200),
        end = c(250, 300, 350)
    )
    values <- c(1.0, 2.0, 3.0)

    gtrack.create_dense("test.tmptrack", "Test overlapping intervals", intervals, values, 100, 0)

    # Extract data and check results
    r <- gextract("test.tmptrack", gintervals(1, 0, 400), colnames = "value")

    # Verify overlapping intervals were processed
    expect_true(nrow(r) > 0)

    # Check if any value between overlapping intervals reflects weighted average
    # For overlapping regions, values should be between 1.0 and 3.0
    expect_true(any(r$value > 1.0 & r$value < 3.0))
})

test_that("test for values length not matching intervals", {
    # Create intervals with more intervals than values
    intervals <- gintervals(
        chrom = c(1, 1, 1),
        start = c(100, 200, 300),
        end = c(150, 250, 350)
    )
    values <- c(1.0, 2.0) # Only 2 values for 3 intervals

    expect_error(
        gtrack.create_dense("test.tmptrack", "Test mismatched lengths", intervals, values, 50, 0),
        "Length of values must match the number of intervals"
    )
})

test_that("test for default value in dense track", {
    gtrack.rm("test.tmptrack", force = TRUE)
    withr::defer(gtrack.rm("test.tmptrack", force = TRUE))

    # Create intervals with a gap
    intervals <- gintervals(
        chrom = c(1, 1),
        start = c(100, 300),
        end = c(150, 350)
    )
    values <- c(1.0, 2.0)

    # Create track with default value NaN
    defval <- NaN
    gtrack.create_dense("test.tmptrack", "Test default value", intervals, values, 50, defval)

    # Extract data and check results
    r <- gextract("test.tmptrack", gintervals(1, 0, 400), colnames = "value")

    # Check if there are regions with NaN (empty regions)
    expect_true(any(is.nan(r$value)))

    # Check that we have both data-containing regions and default regions
    expect_true(any(!is.nan(r$value)))
})

test_that("test for recreating test.fixedbin from extracted data", {
    # Extract data from test.fixedbin
    extracted_data <- gextract("test.fixedbin", gintervals.all(), colnames = "value")

    # Create intervals and values from the extracted data
    intervals <- gintervals(
        chrom = extracted_data$chrom,
        start = extracted_data$start,
        end = extracted_data$end
    )
    values <- extracted_data$value

    # Create a new track with the extracted data
    gtrack.rm("test.tmptrack", force = TRUE)
    withr::defer(gtrack.rm("test.tmptrack", force = TRUE))

    # Get track info to determine bin size
    info <- gtrack.info("test.fixedbin")
    bin_size <- info$bin.size

    # Create a new dense track with the extracted data
    gtrack.create_dense("test.tmptrack", "Recreated test.fixedbin", intervals, values, bin_size, 0)

    # Extract data from the recreated track
    recreated_data <- gextract("test.tmptrack", gintervals.all(), colnames = "value")

    # Compare the original and recreated data
    expect_equal(nrow(extracted_data), nrow(recreated_data))

    # Check track info
    new_info <- gtrack.info("test.tmptrack")
    expect_equal(new_info$type, "dense")
    expect_equal(new_info$bin.size, bin_size)

    expect_equal(gtrack.attr.get("test.tmptrack", "description"), "Recreated test.fixedbin")
    expect_equal(gtrack.attr.get("test.tmptrack", "binsize"), as.character(bin_size))
    expect_equal(gtrack.attr.get("test.tmptrack", "type"), "dense")
    expect_equal(gtrack.attr.get("test.tmptrack", "created.by"), "gtrack.create_dense(test.tmptrack, description, intervals, values, 50, 0)")
})
