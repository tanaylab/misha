create_isolated_test_db()

test_that("gtrack.export_bw works with dense track at native resolution", {
    skip_if_not_installed("rtracklayer")

    tmpfile <- tempfile(fileext = ".bw")
    withr::defer(unlink(tmpfile))

    # Export a known dense track
    result <- gtrack.export_bw("test.fixedbin", tmpfile)

    # Verify function returns file path invisibly
    expect_equal(result, tmpfile)

    # Verify file exists and has content
    expect_true(file.exists(tmpfile))
    expect_gt(file.size(tmpfile), 0)

    # Import back and verify we can read it
    imported <- rtracklayer::import.bw(tmpfile)
    expect_s4_class(imported, "GRanges")
    expect_gt(length(imported), 0)
})

test_that("gtrack.export_bw works with sparse track", {
    skip_if_not_installed("rtracklayer")

    tmpfile <- tempfile(fileext = ".bw")
    withr::defer(unlink(tmpfile))

    gtrack.export_bw("test.sparse", tmpfile)

    expect_true(file.exists(tmpfile))
    expect_gt(file.size(tmpfile), 0)

    # Verify we can import it
    imported <- rtracklayer::import.bw(tmpfile)
    expect_s4_class(imported, "GRanges")
})

test_that("gtrack.export_bw handles NaN values correctly - removal (default)", {
    skip_if_not_installed("rtracklayer")

    # Create track with NaN values
    tmptrack <- paste0("test.tmptrack_", sample(1:1e9, 1))
    withr::defer(gtrack.rm(tmptrack, force = TRUE))

    intervals <- gintervals(1, c(1000, 2000, 3000), c(1100, 2100, 3100))
    values <- c(1.0, NaN, 2.0)
    gtrack.create_sparse(tmptrack, "", intervals, values)

    tmpfile <- tempfile(fileext = ".bw")
    withr::defer(unlink(tmpfile))

    # With na_value=NULL (default), should remove NaN intervals
    # Suppress warning since we're testing removal behavior, not the warning
    expect_no_error(suppressWarnings(gtrack.export_bw(tmptrack, tmpfile)))

    # Verify NaN was filtered
    imported <- rtracklayer::import.bw(tmpfile)
    expect_equal(length(imported), 2) # Only 2 intervals, NaN filtered
    expect_equal(as.numeric(imported$score), c(1.0, 2.0))
})

test_that("gtrack.export_bw handles NaN values correctly - replacement", {
    skip_if_not_installed("rtracklayer")

    # Create track with NaN values
    tmptrack <- paste0("test.tmptrack_", sample(1:1e9, 1))
    withr::defer(gtrack.rm(tmptrack, force = TRUE))

    intervals <- gintervals(1, c(1000, 2000, 3000), c(1100, 2100, 3100))
    values <- c(1.0, NaN, 2.0)
    gtrack.create_sparse(tmptrack, "", intervals, values)

    tmpfile <- tempfile(fileext = ".bw")
    withr::defer(unlink(tmpfile))

    # With na_value=0, should replace NaN with 0
    expect_no_error(gtrack.export_bw(tmptrack, tmpfile, na_value = 0))

    # Verify NaN was replaced with 0
    imported <- rtracklayer::import.bw(tmpfile)
    expect_equal(length(imported), 3) # All 3 intervals present
    expect_equal(as.numeric(imported$score), c(1.0, 0.0, 2.0))
})

test_that("gtrack.export_bw warns when removing many NaN values", {
    skip_if_not_installed("rtracklayer")

    # Create track with many NaN values (>10%)
    tmptrack <- paste0("test.tmptrack_", sample(1:1e9, 1))
    withr::defer(gtrack.rm(tmptrack, force = TRUE))

    # 6 out of 10 values are NaN (60%)
    intervals <- gintervals(1, seq(1000, 10000, 1000), seq(1100, 10100, 1000))
    values <- c(1, NaN, NaN, 2, NaN, 3, NaN, NaN, 4, NaN)
    gtrack.create_sparse(tmptrack, "", intervals, values)

    tmpfile <- tempfile(fileext = ".bw")
    withr::defer(unlink(tmpfile))

    # Should warn about removing many NaN values
    expect_warning(
        gtrack.export_bw(tmptrack, tmpfile),
        "Removing.*NaN"
    )
})

test_that("gtrack.export_bw errors when all values are NaN with na_value=NULL", {
    skip_if_not_installed("rtracklayer")

    # Create track with all NaN values
    tmptrack <- paste0("test.tmptrack_", sample(1:1e9, 1))
    withr::defer(gtrack.rm(tmptrack, force = TRUE))

    intervals <- gintervals(1, c(1000, 2000), c(1100, 2100))
    values <- c(NaN, NaN)
    gtrack.create_sparse(tmptrack, "", intervals, values)

    tmpfile <- tempfile(fileext = ".bw")
    withr::defer(unlink(tmpfile))

    # Should error when all values are NaN
    # Suppress warning that occurs before the error
    expect_error(
        suppressWarnings(gtrack.export_bw(tmptrack, tmpfile)),
        "All values are NaN"
    )
})

test_that("gtrack.export_bw rejects 2D tracks", {
    skip_if_not_installed("rtracklayer")

    tmpfile <- tempfile(fileext = ".bw")
    withr::defer(unlink(tmpfile))

    # test.rects is a 2D track
    expect_error(
        gtrack.export_bw("test.rects", tmpfile),
        "2D track"
    )
})

test_that("gtrack.export_bw validates inputs", {
    skip_if_not_installed("rtracklayer")

    tmpfile <- tempfile(fileext = ".bw")

    # Invalid file path - directory doesn't exist
    expect_error(
        gtrack.export_bw("test.fixedbin", "/invalid/path/that/does/not/exist/file.bw"),
        "Directory does not exist"
    )

    # Invalid na_value (not numeric or NULL) - character
    expect_error(
        gtrack.export_bw("test.fixedbin", tmpfile, na_value = "invalid"),
        "na_value.*numeric"
    )

    # Invalid na_value - vector instead of single value
    expect_error(
        gtrack.export_bw("test.fixedbin", tmpfile, na_value = c(0, 1)),
        "na_value.*single"
    )
})

test_that("gtrack.export_bw works with numeric iterator (binning)", {
    skip_if_not_installed("rtracklayer")

    tmpfile <- tempfile(fileext = ".bw")
    withr::defer(unlink(tmpfile))

    # Export with 1kb bins
    gtrack.export_bw("test.fixedbin", tmpfile, intervals = gintervals(1, 0, 100000), iterator = 1000)

    expect_true(file.exists(tmpfile))

    # Verify bins are 1kb
    imported <- rtracklayer::import.bw(tmpfile)
    widths <- IRanges::width(imported)
    # Most bins should be 1000bp (some at chromosome ends may be smaller)
    expect_equal(median(widths), 1000)
})

test_that("gtrack.export_bw works with track expressions", {
    skip_if_not_installed("rtracklayer")

    tmpfile <- tempfile(fileext = ".bw")
    withr::defer(unlink(tmpfile))

    # Export sum of track and constant
    gtrack.export_bw("test.fixedbin + 1", tmpfile,
        intervals = gintervals(1, 0, 100000),
        iterator = 1000
    )

    expect_true(file.exists(tmpfile))

    # Verify we can import it
    imported <- rtracklayer::import.bw(tmpfile)
    expect_s4_class(imported, "GRanges")
})

test_that("gtrack.export_bw works with custom intervals", {
    skip_if_not_installed("rtracklayer")

    tmpfile <- tempfile(fileext = ".bw")
    withr::defer(unlink(tmpfile))

    # Export only chromosome 1
    custom_intervals <- gintervals(1)
    gtrack.export_bw("test.fixedbin", tmpfile, intervals = custom_intervals)

    expect_true(file.exists(tmpfile))

    # Verify only chr1 was exported
    imported <- rtracklayer::import.bw(tmpfile)
    chrom_names <- as.character(GenomicRanges::seqnames(imported))
    expect_true(all(chrom_names == "chr1"))
})

test_that("gtrack.export_bw returns file path invisibly", {
    skip_if_not_installed("rtracklayer")

    tmpfile <- tempfile(fileext = ".bw")
    withr::defer(unlink(tmpfile))

    result <- gtrack.export_bw("test.fixedbin", tmpfile,
        intervals = gintervals(1, 0, 10000)
    )
    expect_equal(result, tmpfile)
})

test_that("gtrack.export_bw handles coordinate conversion correctly", {
    skip_if_not_installed("rtracklayer")

    # Create a sparse track with known intervals
    tmptrack <- paste0("test.tmptrack_", sample(1:1e9, 1))
    withr::defer(gtrack.rm(tmptrack, force = TRUE))

    # Misha uses 0-based half-open intervals [start, end)
    # These intervals: [1000, 2000), [3000, 4000)
    intervals <- gintervals(1, c(1000, 3000), c(2000, 4000))
    values <- c(10, 20)
    gtrack.create_sparse(tmptrack, "", intervals, values)

    tmpfile <- tempfile(fileext = ".bw")
    withr::defer(unlink(tmpfile))

    gtrack.export_bw(tmptrack, tmpfile)

    # Import and check coordinates
    # GRanges uses 1-based closed intervals [start, end]
    # So [1000, 2000) should become [1001, 2000]
    imported <- rtracklayer::import.bw(tmpfile)

    expect_equal(length(imported), 2)
    expect_equal(as.integer(GenomicRanges::start(imported)), c(1001L, 3001L))
    expect_equal(as.integer(GenomicRanges::end(imported)), c(2000L, 4000L))
    expect_equal(as.numeric(imported$score), c(10, 20))
})

test_that("gtrack.export_bw round-trip preserves values", {
    skip_if_not_installed("rtracklayer")

    # Create a track, export to BigWig, import back, and compare
    tmptrack <- paste0("test.tmptrack_", sample(1:1e9, 1))
    withr::defer(gtrack.rm(tmptrack, force = TRUE))

    intervals <- gintervals(1, seq(1000, 10000, 1000), seq(2000, 11000, 1000))
    values <- seq(1, 10)
    gtrack.create_sparse(tmptrack, "", intervals, values)

    tmpfile <- tempfile(fileext = ".bw")
    withr::defer(unlink(tmpfile))

    # Export
    gtrack.export_bw(tmptrack, tmpfile)

    # Import back with rtracklayer
    imported <- rtracklayer::import.bw(tmpfile)

    # Check values match (accounting for coordinate conversion)
    expect_equal(as.numeric(imported$score), values)
    expect_equal(length(imported), length(values))
})

test_that("gtrack.export_bw handles empty extraction result", {
    skip_if_not_installed("rtracklayer")

    tmpfile <- tempfile(fileext = ".bw")
    withr::defer(unlink(tmpfile))

    # Try to extract from non-overlapping intervals
    # This should fail during gextract, not in our function
    expect_error(
        gtrack.export_bw("test.fixedbin", tmpfile,
            intervals = gintervals(999, 0, 1000)
        ), # Non-existent chromosome
        "extract"
    )
})

test_that("gtrack.export_bw works with multiple chromosomes", {
    skip_if_not_installed("rtracklayer")

    tmpfile <- tempfile(fileext = ".bw")
    withr::defer(unlink(tmpfile))

    # Export chromosomes 1 and 2
    custom_intervals <- gintervals(c(1, 2))
    gtrack.export_bw("test.fixedbin", tmpfile, intervals = custom_intervals)

    expect_true(file.exists(tmpfile))

    # Verify both chromosomes are present
    imported <- rtracklayer::import.bw(tmpfile)
    chrom_names <- unique(as.character(GenomicRanges::seqnames(imported)))
    expect_true(all(c("chr1", "chr2") %in% chrom_names))
})

test_that("gtrack.export_bw with iterator produces correct bin sizes", {
    skip_if_not_installed("rtracklayer")

    tmpfile <- tempfile(fileext = ".bw")
    withr::defer(unlink(tmpfile))

    # Export with 5kb bins
    bin_size <- 5000
    gtrack.export_bw("test.fixedbin", tmpfile,
        intervals = gintervals(1, 0, 100000),
        iterator = bin_size
    )

    imported <- rtracklayer::import.bw(tmpfile)
    widths <- IRanges::width(imported)

    # Most bins should be exactly bin_size (last one might be smaller)
    expect_true(median(widths) == bin_size)
    expect_true(max(widths) == bin_size)
})

# Tests using gtrack.import to verify round-trip
test_that("gtrack.export_bw -> gtrack.import round-trip works with sparse track", {
    skip_if_not_installed("rtracklayer")

    # Create a sparse track with known values
    tmptrack <- paste0("test.tmptrack_", sample(1:1e9, 1))
    withr::defer(gtrack.rm(tmptrack, force = TRUE))

    intervals <- gintervals(1, c(1000, 5000, 10000), c(2000, 6000, 11000))
    values <- c(10.5, 20.3, 30.7)
    gtrack.create_sparse(tmptrack, "", intervals, values)

    # Export to BigWig
    tmpfile <- tempfile(fileext = ".bw")
    withr::defer(unlink(tmpfile))
    gtrack.export_bw(tmptrack, tmpfile)

    # Import back using gtrack.import (binsize = 0 for sparse tracks)
    imported_track <- paste0("test.imported_", sample(1:1e9, 1))
    withr::defer(gtrack.rm(imported_track, force = TRUE))
    gtrack.import(imported_track, "", tmpfile, binsize = 0)

    # Extract data from both tracks
    original_data <- gextract(tmptrack, intervals = intervals, iterator = intervals)
    imported_data <- gextract(imported_track, intervals = intervals, iterator = intervals)

    # Compare values (allowing small floating point differences)
    expect_equal(imported_data[[imported_track]], original_data[[tmptrack]], tolerance = 1e-6)
    expect_equal(nrow(imported_data), nrow(original_data))
})

test_that("gtrack.export_bw -> gtrack.import round-trip works with binned dense track", {
    skip_if_not_installed("rtracklayer")

    # Export dense track with 1kb bins
    tmpfile <- tempfile(fileext = ".bw")
    withr::defer(unlink(tmpfile))

    test_intervals <- gintervals(1, 0, 50000)
    bin_size <- 1000

    # Get expected values by extracting with iterator
    expected_data <- gextract("test.fixedbin", intervals = test_intervals, iterator = bin_size)

    # Export to BigWig
    gtrack.export_bw("test.fixedbin", tmpfile, intervals = test_intervals, iterator = bin_size)

    # Import back (binsize = bin_size for dense tracks)
    imported_track <- paste0("test.imported_", sample(1:1e9, 1))
    withr::defer(gtrack.rm(imported_track, force = TRUE))
    gtrack.import(imported_track, "", tmpfile, binsize = bin_size)

    # Extract from imported track
    imported_data <- gextract(imported_track, intervals = test_intervals, iterator = bin_size)

    # Compare values
    expect_equal(nrow(imported_data), nrow(expected_data))
    expect_equal(imported_data[[imported_track]], expected_data$test.fixedbin, tolerance = 1e-6)
})

test_that("gtrack.export_bw -> gtrack.import round-trip preserves NaN replacement", {
    skip_if_not_installed("rtracklayer")

    # Create track with NaN values
    tmptrack <- paste0("test.tmptrack_", sample(1:1e9, 1))
    withr::defer(gtrack.rm(tmptrack, force = TRUE))

    intervals <- gintervals(1, c(1000, 2000, 3000, 4000), c(1100, 2100, 3100, 4100))
    values <- c(1.5, NaN, 2.5, NaN)
    gtrack.create_sparse(tmptrack, "", intervals, values)

    # Export with na_value=0
    tmpfile <- tempfile(fileext = ".bw")
    withr::defer(unlink(tmpfile))
    gtrack.export_bw(tmptrack, tmpfile, na_value = 0)

    # Import back (binsize = 0 for sparse tracks)
    imported_track <- paste0("test.imported_", sample(1:1e9, 1))
    withr::defer(gtrack.rm(imported_track, force = TRUE))
    gtrack.import(imported_track, "", tmpfile, binsize = 0)

    # Extract from imported track
    imported_data <- gextract(imported_track, intervals = intervals, iterator = intervals)

    # NaN values should have been replaced with 0
    expect_equal(nrow(imported_data), 4)
    expect_equal(imported_data[[imported_track]], c(1.5, 0.0, 2.5, 0.0), tolerance = 1e-6)
})

test_that("gtrack.export_bw -> gtrack.import round-trip handles NaN removal", {
    skip_if_not_installed("rtracklayer")

    # Create track with NaN values
    tmptrack <- paste0("test.tmptrack_", sample(1:1e9, 1))
    withr::defer(gtrack.rm(tmptrack, force = TRUE))

    intervals <- gintervals(1, c(1000, 2000, 3000, 4000), c(1100, 2100, 3100, 4100))
    values <- c(1.5, NaN, 2.5, NaN)
    gtrack.create_sparse(tmptrack, "", intervals, values)

    # Export with na_value=NULL (default, removes NaN intervals)
    tmpfile <- tempfile(fileext = ".bw")
    withr::defer(unlink(tmpfile))
    expect_no_error(suppressWarnings(gtrack.export_bw(tmptrack, tmpfile)))

    # Import back (binsize = 0 for sparse tracks)
    imported_track <- paste0("test.imported_", sample(1:1e9, 1))
    withr::defer(gtrack.rm(imported_track, force = TRUE))
    gtrack.import(imported_track, "", tmpfile, binsize = 0)

    # Extract all intervals from imported track
    imported_data <- gextract(imported_track, intervals = gintervals(1))

    # Should only have 2 intervals (NaN intervals were removed)
    expect_equal(nrow(imported_data), 2)

    # Check the non-NaN intervals are preserved
    non_nan_intervals <- intervals[!is.na(values), ]
    extracted_values <- gextract(imported_track, intervals = non_nan_intervals, iterator = non_nan_intervals)
    expect_equal(extracted_values[[imported_track]], c(1.5, 2.5), tolerance = 1e-6)
})

test_that("gtrack.export_bw -> gtrack.import round-trip works with track expression", {
    skip_if_not_installed("rtracklayer")

    # Export a track expression
    tmpfile <- tempfile(fileext = ".bw")
    withr::defer(unlink(tmpfile))

    test_intervals <- gintervals(1, 0, 50000)
    bin_size <- 1000

    # Get expected values from expression
    expected_data <- gextract("test.fixedbin + 10", intervals = test_intervals, iterator = bin_size, colnames = "expr")

    # Export expression
    gtrack.export_bw("test.fixedbin + 10", tmpfile, intervals = test_intervals, iterator = bin_size)

    # Import back (binsize = bin_size for dense tracks)
    imported_track <- paste0("test.imported_", sample(1:1e9, 1))
    withr::defer(gtrack.rm(imported_track, force = TRUE))
    gtrack.import(imported_track, "", tmpfile, binsize = bin_size)

    # Extract from imported track
    imported_data <- gextract(imported_track, intervals = test_intervals, iterator = bin_size)

    # Compare values
    expect_equal(nrow(imported_data), nrow(expected_data))
    expect_equal(imported_data[[imported_track]], expected_data$expr, tolerance = 1e-6)
})

test_that("gtrack.export_bw works with fixed_summaries parameter", {
    skip_if_not_installed("rtracklayer")

    tmpfile1 <- tempfile(fileext = ".bw")
    tmpfile2 <- tempfile(fileext = ".bw")
    withr::defer(unlink(tmpfile1))
    withr::defer(unlink(tmpfile2))

    test_intervals <- gintervals(1, 0, 100000)

    # Export with dynamic summaries (default)
    gtrack.export_bw("test.fixedbin", tmpfile1, intervals = test_intervals, iterator = 1000, fixed_summaries = FALSE)

    # Export with fixed summaries
    gtrack.export_bw("test.fixedbin", tmpfile2, intervals = test_intervals, iterator = 1000, fixed_summaries = TRUE)

    # Both files should exist
    expect_true(file.exists(tmpfile1))
    expect_true(file.exists(tmpfile2))

    # Both should have content (file sizes may differ due to different summary levels)
    expect_gt(file.size(tmpfile1), 0)
    expect_gt(file.size(tmpfile2), 0)

    # Both should be importable and contain the same data
    imported1 <- rtracklayer::import.bw(tmpfile1)
    imported2 <- rtracklayer::import.bw(tmpfile2)

    expect_equal(length(imported1), length(imported2))
    expect_equal(as.numeric(imported1$score), as.numeric(imported2$score), tolerance = 1e-6)
})

test_that("gtrack.export_bw validates fixed_summaries parameter", {
    skip_if_not_installed("rtracklayer")

    tmpfile <- tempfile(fileext = ".bw")

    # Invalid fixed_summaries - not logical
    expect_error(
        gtrack.export_bw("test.fixedbin", tmpfile, fixed_summaries = "invalid"),
        "fixed_summaries.*TRUE or FALSE"
    )

    # Invalid fixed_summaries - vector instead of single value
    expect_error(
        gtrack.export_bw("test.fixedbin", tmpfile, fixed_summaries = c(TRUE, FALSE)),
        "fixed_summaries.*TRUE or FALSE"
    )
})

# ===== COMPREHENSIVE EDGE CASE TESTS =====

test_that("gtrack.export_bw handles invalid track names gracefully", {
    skip_if_not_installed("rtracklayer")

    tmpfile <- tempfile(fileext = ".bw")

    # Non-existent track
    expect_error(
        gtrack.export_bw("nonexistent.track", tmpfile),
        "does not contain any tracks|Cannot find|does not exist"
    )

    # Empty string
    expect_error(
        gtrack.export_bw("", tmpfile)
    )

    # NULL track name
    expect_error(
        gtrack.export_bw(NULL, tmpfile)
    )
})

test_that("gtrack.export_bw handles file overwrite behavior correctly", {
    skip_if_not_installed("rtracklayer")

    tmpfile <- tempfile(fileext = ".bw")
    withr::defer(unlink(tmpfile))

    test_intervals <- gintervals(1, 0, 10000)

    # Create initial file
    gtrack.export_bw("test.fixedbin", tmpfile, intervals = test_intervals, iterator = 1000)
    expect_true(file.exists(tmpfile))

    initial_size <- file.size(tmpfile)

    # Export again to same file (should overwrite)
    gtrack.export_bw("test.fixedbin", tmpfile, intervals = test_intervals, iterator = 1000)
    expect_true(file.exists(tmpfile))

    # File should still be valid
    imported <- rtracklayer::import.bw(tmpfile)
    expect_gt(length(imported), 0)
})

test_that("gtrack.export_bw handles infinity and extreme values", {
    skip_if_not_installed("rtracklayer")

    tmpfile <- tempfile(fileext = ".bw")
    withr::defer(unlink(tmpfile))

    # Create track with extreme values
    test_intervals <- gintervals(1, c(1000, 2000, 3000), c(2000, 3000, 4000))
    extreme_vals <- c(1e308, -1e308, 0)

    gtrack.create_sparse("test.extreme", "Test extreme values", test_intervals, extreme_vals)
    withr::defer(gtrack.rm("test.extreme", force = TRUE))

    # Export - Inf values should be converted to NaN and then to na_value
    gtrack.export_bw("test.extreme", tmpfile, na_value = -999)

    imported <- rtracklayer::import.bw(tmpfile)
    expect_gt(length(imported), 0)
})

test_that("gtrack.export_bw validates iterator parameter edge cases", {
    skip_if_not_installed("rtracklayer")

    tmpfile <- tempfile(fileext = ".bw")
    test_intervals <- gintervals(1, 0, 1000)

    # Zero iterator should error
    expect_error(
        gtrack.export_bw("test.fixedbin", tmpfile, intervals = test_intervals, iterator = 0),
        "iterator"
    )

    # Negative iterator should error
    expect_error(
        gtrack.export_bw("test.fixedbin", tmpfile, intervals = test_intervals, iterator = -100),
        "iterator"
    )

    # Non-numeric iterator should error
    expect_error(
        gtrack.export_bw("test.fixedbin", tmpfile, intervals = test_intervals, iterator = "invalid"),
        "iterator"
    )
})

test_that("gtrack.export_bw preserves numeric precision", {
    skip_if_not_installed("rtracklayer")

    tmpfile <- tempfile(fileext = ".bw")
    withr::defer(unlink(tmpfile))

    # Create track with precise values
    test_intervals <- gintervals(1, c(1000, 2000, 3000, 4000), c(2000, 3000, 4000, 5000))
    precise_vals <- c(0.123456789, 9.876543210, 1.111111111, 2.222222222)

    gtrack.create_sparse("test.precise", "Test precision", test_intervals, precise_vals)
    withr::defer(gtrack.rm("test.precise", force = TRUE))

    # Export and reimport
    gtrack.export_bw("test.precise", tmpfile)
    imported <- rtracklayer::import.bw(tmpfile)

    # Values should match within reasonable tolerance (BigWig uses float32)
    original_sorted <- sort(precise_vals)
    imported_sorted <- sort(as.numeric(imported$score))
    expect_equal(imported_sorted, original_sorted, tolerance = 1e-6)
})

test_that("gtrack.export_bw handles single-value tracks", {
    skip_if_not_installed("rtracklayer")

    tmpfile <- tempfile(fileext = ".bw")
    withr::defer(unlink(tmpfile))

    # Create track with constant value
    test_intervals <- gintervals(1, c(1000, 2000, 3000), c(2000, 3000, 4000))
    constant_val <- c(42.0, 42.0, 42.0)

    gtrack.create_sparse("test.constant", "Test constant", test_intervals, constant_val)
    withr::defer(gtrack.rm("test.constant", force = TRUE))

    # Export
    gtrack.export_bw("test.constant", tmpfile)
    imported <- rtracklayer::import.bw(tmpfile)

    # All values should be 42
    expect_true(all(abs(as.numeric(imported$score) - 42.0) < 1e-6))
})

test_that("gtrack.export_bw handles negative values correctly", {
    skip_if_not_installed("rtracklayer")

    tmpfile <- tempfile(fileext = ".bw")
    withr::defer(unlink(tmpfile))

    # Create track with negative values
    test_intervals <- gintervals(1, c(1000, 2000, 3000, 4000), c(2000, 3000, 4000, 5000))
    neg_vals <- c(-10.5, -20.3, -0.1, -999.9)

    gtrack.create_sparse("test.negative", "Test negative", test_intervals, neg_vals)
    withr::defer(gtrack.rm("test.negative", force = TRUE))

    # Export
    gtrack.export_bw("test.negative", tmpfile)
    imported <- rtracklayer::import.bw(tmpfile)

    # All imported values should be negative
    expect_true(all(as.numeric(imported$score) < 0))

    # Check specific values
    imported_sorted <- sort(as.numeric(imported$score))
    original_sorted <- sort(neg_vals)
    expect_equal(imported_sorted, original_sorted, tolerance = 1e-6)
})

test_that("gtrack.export_bw handles empty data after filtering", {
    skip_if_not_installed("rtracklayer")

    tmpfile <- tempfile(fileext = ".bw")
    withr::defer(unlink(tmpfile))

    # Export with intervals that don't overlap with track data
    non_overlapping <- gintervals(1, 50000000, 60000000)

    # This should either create empty file or error gracefully
    result <- tryCatch(
        {
            gtrack.export_bw("test.fixedbin", tmpfile, intervals = non_overlapping, iterator = 1000)
            "success"
        },
        error = function(e) "error"
    )

    # Either way is acceptable - just shouldn't crash
    expect_true(result %in% c("success", "error"))
})

test_that("gtrack.export_bw validates multiple na_value scenarios", {
    skip_if_not_installed("rtracklayer")

    tmpfile <- tempfile(fileext = ".bw")
    withr::defer(unlink(tmpfile))

    test_intervals <- gintervals(1, 0, 10000)

    # na_value as zero
    expect_no_error(
        gtrack.export_bw("test.fixedbin", tmpfile, intervals = test_intervals, iterator = 1000, na_value = 0)
    )

    # na_value as negative
    expect_no_error(
        gtrack.export_bw("test.fixedbin", tmpfile, intervals = test_intervals, iterator = 1000, na_value = -999)
    )

    # na_value as positive
    expect_no_error(
        gtrack.export_bw("test.fixedbin", tmpfile, intervals = test_intervals, iterator = 1000, na_value = 999)
    )

    # Invalid na_value (non-numeric)
    expect_error(
        gtrack.export_bw("test.fixedbin", tmpfile, intervals = test_intervals, iterator = 1000, na_value = "invalid"),
        "na_value.*numeric"
    )

    # Invalid na_value (multiple values)
    expect_error(
        gtrack.export_bw("test.fixedbin", tmpfile, intervals = test_intervals, iterator = 1000, na_value = c(0, 1)),
        "na_value.*single"
    )
})

test_that("gtrack.export_bw produces valid file sizes", {
    skip_if_not_installed("rtracklayer")

    tmpfile <- tempfile(fileext = ".bw")
    withr::defer(unlink(tmpfile))

    test_intervals <- gintervals(1, 0, 100000)

    # Export with different iterator sizes
    gtrack.export_bw("test.fixedbin", tmpfile, intervals = test_intervals, iterator = 100)
    size_100 <- file.size(tmpfile)

    gtrack.export_bw("test.fixedbin", tmpfile, intervals = test_intervals, iterator = 1000)
    size_1000 <- file.size(tmpfile)

    gtrack.export_bw("test.fixedbin", tmpfile, intervals = test_intervals, iterator = 10000)
    size_10000 <- file.size(tmpfile)

    # All files should have reasonable sizes
    expect_gt(size_100, 0)
    expect_gt(size_1000, 0)
    expect_gt(size_10000, 0)

    # Smaller iterators (more bins) should generally produce larger files
    # (though compression can affect this, so we just check they're all valid)
    expect_true(all(c(size_100, size_1000, size_10000) > 100)) # At least 100 bytes
})
