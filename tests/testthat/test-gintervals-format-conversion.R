skip_if(getOption("gmulticontig.indexed_format", FALSE), "Indexed format not enabled, set gmulticontig.indexed_format = TRUE to run this test")

test_that("gintervals.convert_to_indexed fails for non-existent interval set", {
    expect_error(
        gintervals.convert_to_indexed("nonexistent_set"),
        "does not exist"
    )
})

test_that("gintervals.convert_to_indexed handles single-file format correctly", {
    # Clean up any existing test set
    test_set <- "test_1d_intervals"
    if (gintervals.exists(test_set)) gintervals.rm(test_set, force = TRUE)

    # Create a small test interval set (will be stored as single file)
    intervs <- gintervals(c(1, 1, 2), c(1000, 5000, 1000), c(2000, 6000, 3000))
    intervs$score <- rnorm(nrow(intervs))

    gintervals.save(test_set, intervs)

    # Load the original data for comparison
    original <- gintervals.load(test_set)

    # Try to convert - should skip since it's single-file format
    expect_message(
        gintervals.convert_to_indexed(test_set, remove.old = FALSE),
        "single-file format"
    )

    # Data should still be loadable and identical
    after_convert_attempt <- gintervals.load(test_set)
    expect_equal(original, after_convert_attempt)

    # Clean up
    gintervals.rm(test_set, force = TRUE)
})

test_that("gintervals.convert_to_indexed works with existing big set intervals", {
    # This test verifies convert function works with actual Big Set intervals
    # Note: Big Set intervals require large data (>1MB), so we skip if none exist

    # Look for existing Big Set intervals in the test database
    all_intervs <- gintervals.ls()

    # Find a big set interval (if any)
    big_set <- NULL
    for (intervset in all_intervs) {
        path <- gsub("\\.", "/", intervset)
        intervset_path <- paste0(get("GWD", envir = .misha), "/", path, ".interv")
        if (dir.exists(intervset_path)) {
            big_set <- intervset
            break
        }
    }

    skip_if(is.null(big_set), "No Big Set intervals available for testing")

    # Try converting it - should either convert or say already converted
    expect_no_error(gintervals.convert_to_indexed(big_set, remove.old = FALSE))

    # Data should still be loadable
    expect_no_error(gintervals.load(big_set))
})

test_that("gintervals.2d.convert_to_indexed fails for non-existent interval set", {
    expect_error(
        gintervals.2d.convert_to_indexed("nonexistent_2d_set"),
        "does not exist"
    )
})

test_that("gintervals.convert_to_indexed detects when intervals don't exist", {
    expect_error(
        gintervals.convert_to_indexed("definitely_nonexistent_set_12345"),
        "does not exist"
    )
})
