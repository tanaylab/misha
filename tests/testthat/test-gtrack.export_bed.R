load_test_db()

test_that("gtrack.export_bed exports sparse track to BED5 format", {
    tmpfile <- tempfile(fileext = ".bed")
    withr::defer(unlink(tmpfile))

    # Create sparse track
    test_intervals <- gintervals(c(1, 1, 2), c(1000, 5000, 2000), c(2000, 6000, 3000))
    test_values <- c(10.5, 20.3, 15.7)

    tmptrack <- paste0("test.bed_export_", sample(1:1e9, 1))
    gtrack.create_sparse(tmptrack, "Test sparse track", test_intervals, test_values)
    withr::defer(gtrack.rm(tmptrack, force = TRUE))

    # Export
    result <- gtrack.export_bed(tmptrack, tmpfile)

    # Should return file path invisibly
    expect_equal(result, tmpfile)
    expect_true(file.exists(tmpfile))

    # Read and verify BED5 format (chrom, start, end, name, score)
    bed_content <- read.table(tmpfile, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
    expect_equal(nrow(bed_content), 3)
    expect_equal(ncol(bed_content), 5) # BED5

    # Verify coordinates (no conversion - 0-based preserved)
    expect_equal(bed_content[, 2], test_intervals$start)
    expect_equal(bed_content[, 3], test_intervals$end)

    # Verify scores (use tolerance for float precision)
    expect_equal(bed_content[, 5], test_values, tolerance = 1e-6)
})

test_that("gtrack.export_bed exports dense track to BED5 format", {
    tmpfile <- tempfile(fileext = ".bed")
    withr::defer(unlink(tmpfile))

    # Export existing dense track
    test_intervals <- gintervals(1, 0, 10000)
    gtrack.export_bed("test.fixedbin", tmpfile, intervals = test_intervals)

    expect_true(file.exists(tmpfile))

    # Read and verify BED5 format
    bed_content <- read.table(tmpfile, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
    expect_gt(nrow(bed_content), 0)
    expect_equal(ncol(bed_content), 5) # BED5 (has track values)
})

test_that("gtrack.export_bed works with track expressions and iterator", {
    tmpfile <- tempfile(fileext = ".bed")
    withr::defer(unlink(tmpfile))

    # Export track expression with iterator
    test_intervals <- gintervals(1, 0, 10000)
    gtrack.export_bed("test.fixedbin * 2", tmpfile,
        intervals = test_intervals,
        iterator = 1000
    )

    expect_true(file.exists(tmpfile))

    bed_content <- read.table(tmpfile, sep = "\t", header = FALSE)
    expect_gt(nrow(bed_content), 0)
    expect_equal(ncol(bed_content), 5) # BED5
})

test_that("gtrack.export_bed handles NaN values with na_value parameter", {
    tmpfile1 <- tempfile(fileext = ".bed")
    tmpfile2 <- tempfile(fileext = ".bed")
    withr::defer(unlink(tmpfile1))
    withr::defer(unlink(tmpfile2))

    # Create track with NaN values
    test_intervals <- gintervals(1, c(1000, 2000, 3000), c(2000, 3000, 4000))
    test_values <- c(10.5, NaN, 15.7)

    tmptrack <- paste0("test.nan_", sample(1:1e9, 1))
    gtrack.create_sparse(tmptrack, "Test NaN", test_intervals, test_values)
    withr::defer(gtrack.rm(tmptrack, force = TRUE))

    # Export with default (NaN written as ".")
    gtrack.export_bed(tmptrack, tmpfile1)
    bed1 <- readLines(tmpfile1)
    expect_match(bed1[2], "\\.$") # Second line should end with "."

    # Export with na_value replacement
    gtrack.export_bed(tmptrack, tmpfile2, na_value = -999)
    bed2 <- read.table(tmpfile2, sep = "\t", header = FALSE)
    expect_equal(bed2[2, 5], -999)
})

test_that("gtrack.export_bed adds track header when requested", {
    tmpfile <- tempfile(fileext = ".bed")
    withr::defer(unlink(tmpfile))

    test_intervals <- gintervals(1, 0, 5000)

    # Export with track name only
    gtrack.export_bed("test.fixedbin", tmpfile,
        intervals = test_intervals,
        track_name = "MyTrack"
    )

    lines <- readLines(tmpfile)
    expect_match(lines[1], 'track name="MyTrack"')

    # Export with both track name and description
    tmpfile2 <- tempfile(fileext = ".bed")
    withr::defer(unlink(tmpfile2))

    gtrack.export_bed("test.fixedbin", tmpfile2,
        intervals = test_intervals,
        track_name = "MyTrack",
        description = "Test data"
    )

    lines2 <- readLines(tmpfile2)
    expect_match(lines2[1], 'track name="MyTrack" description="Test data"')
})

test_that("gtrack.export_bed round-trip preserves data", {
    tmpfile <- tempfile(fileext = ".bed")
    withr::defer(unlink(tmpfile))

    # Create test track
    test_intervals <- gintervals(c(1, 1, 2), c(1000, 5000, 2000), c(2000, 6000, 3000))
    test_values <- c(10.5, 20.3, 15.7)

    tmptrack1 <- paste0("test.roundtrip1_", sample(1:1e9, 1))
    gtrack.create_sparse(tmptrack1, "Test track", test_intervals, test_values)
    withr::defer(gtrack.rm(tmptrack1, force = TRUE))

    # Export to BED
    gtrack.export_bed(tmptrack1, tmpfile)

    # Import back
    tmptrack2 <- paste0("test.roundtrip2_", sample(1:1e9, 1))
    withr::defer(gtrack.rm(tmptrack2, force = TRUE))
    gtrack.import(tmptrack2, "", tmpfile, binsize = 0)

    # Compare original and reimported
    original <- gextract(tmptrack1, gintervals.all())
    reimported <- gextract(tmptrack2, gintervals.all())

    expect_equal(nrow(original), nrow(reimported))
    expect_equal(original$start, reimported$start)
    expect_equal(original$end, reimported$end)
    expect_equal(original[[5]], reimported[[5]], tolerance = 1e-6) # Values
})

test_that("gtrack.export_bed preserves coordinate system (no conversion)", {
    tmpfile <- tempfile(fileext = ".bed")
    withr::defer(unlink(tmpfile))

    # Create track with specific coordinates including edge cases
    test_intervals <- gintervals(c(1, 1, 2), c(0, 100, 1000), c(50, 200, 2000))
    test_values <- c(5.0, 10.0, 15.0)

    tmptrack <- paste0("test.coords_", sample(1:1e9, 1))
    gtrack.create_sparse(tmptrack, "Test coords", test_intervals, test_values)
    withr::defer(gtrack.rm(tmptrack, force = TRUE))

    # Export
    gtrack.export_bed(tmptrack, tmpfile)

    bed_content <- read.table(tmpfile, sep = "\t", header = FALSE)

    # Coordinates should be exactly as input (no +1 conversion like BigWig)
    expect_equal(bed_content[, 2], c(0, 100, 1000))
    expect_equal(bed_content[, 3], c(50, 200, 2000))
})

test_that("gtrack.export_bed handles multiple chromosomes", {
    tmpfile <- tempfile(fileext = ".bed")
    withr::defer(unlink(tmpfile))

    # Create track with multiple chromosomes
    test_intervals <- gintervals(
        c(1, 1, 2, 2, 3),
        c(1000, 5000, 2000, 6000, 1000),
        c(2000, 6000, 3000, 7000, 2000)
    )
    test_values <- c(10, 20, 30, 40, 50)

    tmptrack <- paste0("test.multichrom_", sample(1:1e9, 1))
    gtrack.create_sparse(tmptrack, "Test multi", test_intervals, test_values)
    withr::defer(gtrack.rm(tmptrack, force = TRUE))

    # Export
    gtrack.export_bed(tmptrack, tmpfile)

    bed_content <- read.table(tmpfile, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
    expect_equal(nrow(bed_content), 5)

    # Verify chromosomes preserved
    expect_equal(bed_content[, 1], c("chr1", "chr1", "chr2", "chr2", "chr3"))
})

test_that("gtrack.export_bed validates file parameter", {
    # Non-character file
    expect_error(
        gtrack.export_bed("test.fixedbin", 123),
        "must be a single character string"
    )

    # Multiple file paths
    expect_error(
        gtrack.export_bed("test.fixedbin", c("file1.bed", "file2.bed")),
        "must be a single character string"
    )

    # Non-existent directory
    expect_error(
        gtrack.export_bed("test.fixedbin", "/nonexistent_dir_12345/file.bed"),
        "Directory does not exist"
    )
})

test_that("gtrack.export_bed validates na_value parameter", {
    tmpfile <- tempfile(fileext = ".bed")

    # Non-numeric na_value
    expect_error(
        gtrack.export_bed("test.fixedbin", tmpfile, na_value = "invalid"),
        "must be a single numeric value"
    )

    # Multiple values
    expect_error(
        gtrack.export_bed("test.fixedbin", tmpfile, na_value = c(0, 1)),
        "must be a single numeric value"
    )

    # NA value
    expect_error(
        gtrack.export_bed("test.fixedbin", tmpfile, na_value = NA),
        "must be a single numeric value"
    )
})

test_that("gtrack.export_bed validates track_name and description", {
    tmpfile <- tempfile(fileext = ".bed")
    withr::defer(unlink(tmpfile))

    # Invalid track_name type
    expect_error(
        gtrack.export_bed("test.fixedbin", tmpfile, track_name = 123),
        "must be a single character string"
    )

    # Invalid description type
    expect_error(
        gtrack.export_bed("test.fixedbin", tmpfile,
            track_name = "test",
            description = 123
        ),
        "must be a single character string"
    )
})

test_that("gtrack.export_bed handles invalid track names", {
    tmpfile <- tempfile(fileext = ".bed")

    # Non-existent track
    expect_error(
        gtrack.export_bed("nonexistent.track.12345", tmpfile),
        "Failed to extract track data"
    )
})

test_that("gtrack.export_bed handles file overwrite", {
    tmpfile <- tempfile(fileext = ".bed")
    withr::defer(unlink(tmpfile))

    test_intervals <- gintervals(1, 0, 5000)

    # Write first file
    gtrack.export_bed("test.fixedbin", tmpfile, intervals = test_intervals, iterator = 1000)
    content1 <- readLines(tmpfile)
    size1 <- length(content1)

    # Overwrite with different iterator (different data)
    test_intervals2 <- gintervals(1, 0, 10000)
    gtrack.export_bed("test.fixedbin", tmpfile, intervals = test_intervals2, iterator = 500)
    content2 <- readLines(tmpfile)
    size2 <- length(content2)

    # Second export should have more lines (smaller iterator = more bins)
    expect_gt(size2, size1)
})

test_that("gtrack.export_bed works with custom intervals", {
    tmpfile <- tempfile(fileext = ".bed")
    withr::defer(unlink(tmpfile))

    # Export only specific region
    custom_intervals <- gintervals(1, 5000, 15000)
    gtrack.export_bed("test.fixedbin", tmpfile,
        intervals = custom_intervals,
        iterator = 1000
    )

    bed_content <- read.table(tmpfile, sep = "\t", header = FALSE)

    # All intervals should be within custom range
    expect_true(all(bed_content[, 2] >= 5000))
    expect_true(all(bed_content[, 3] <= 15000))
})

test_that("gtrack.export_bed works with virtual tracks", {
    tmpfile <- tempfile(fileext = ".bed")
    withr::defer(unlink(tmpfile))

    # Export virtual track with iterator (using masked.frac which exists)
    test_intervals <- gintervals(1, 1000, 5000)

    # Create virtual track first
    gvtrack.create("test_masked", NULL, "masked.frac")
    withr::defer(remove_all_vtracks())

    expect_no_error(
        gtrack.export_bed("test_masked", tmpfile,
            intervals = test_intervals,
            iterator = 500
        )
    )

    expect_true(file.exists(tmpfile))

    bed_content <- read.table(tmpfile, sep = "\t", header = FALSE)
    expect_gt(nrow(bed_content), 0)
    expect_equal(ncol(bed_content), 5) # BED5
})

test_that("gtrack.export_bed handles negative values", {
    tmpfile <- tempfile(fileext = ".bed")
    withr::defer(unlink(tmpfile))

    # Create track with negative values
    test_intervals <- gintervals(1, c(1000, 2000, 3000), c(2000, 3000, 4000))
    test_values <- c(-10.5, -20.3, -5.7)

    tmptrack <- paste0("test.negative_", sample(1:1e9, 1))
    gtrack.create_sparse(tmptrack, "Test negative", test_intervals, test_values)
    withr::defer(gtrack.rm(tmptrack, force = TRUE))

    # Export
    gtrack.export_bed(tmptrack, tmpfile)

    bed_content <- read.table(tmpfile, sep = "\t", header = FALSE)

    # Verify negative values preserved (use tolerance for float precision)
    expect_true(all(bed_content[, 5] < 0))
    expect_equal(bed_content[, 5], test_values, tolerance = 1e-6)
})

test_that("gtrack.export_bed produces reasonable file sizes", {
    tmpfile <- tempfile(fileext = ".bed")
    withr::defer(unlink(tmpfile))

    test_intervals <- gintervals(1, 0, 100000)

    # Export with different iterator sizes
    gtrack.export_bed("test.fixedbin", tmpfile, intervals = test_intervals, iterator = 1000)
    size1 <- file.size(tmpfile)

    gtrack.export_bed("test.fixedbin", tmpfile, intervals = test_intervals, iterator = 10000)
    size2 <- file.size(tmpfile)

    # All files should have reasonable sizes
    expect_gt(size1, 100) # Not empty
    expect_gt(size2, 100) # Not empty

    # Smaller iterator should produce larger file (more bins)
    expect_gt(size1, size2)
})

test_that("gtrack.export_bed handles regions with no track coverage", {
    tmpfile <- tempfile(fileext = ".bed")
    withr::defer(unlink(tmpfile))

    # Export with intervals that may have sparse/no track data
    # (gextract will still return bins, but values may be NaN)
    test_intervals <- gintervals(1, 50000000, 50010000)

    # Should complete without error (creates file with NaN values as ".")
    expect_no_error(
        gtrack.export_bed("test.fixedbin", tmpfile,
            intervals = test_intervals,
            iterator = 1000
        )
    )

    # File should be created
    expect_true(file.exists(tmpfile))
})

test_that("gtrack.export_bed exports numeric precision correctly", {
    tmpfile <- tempfile(fileext = ".bed")
    withr::defer(unlink(tmpfile))

    # Create track with precise decimal values
    test_intervals <- gintervals(1, c(1000, 2000, 3000), c(2000, 3000, 4000))
    test_values <- c(0.123456789, 9.876543210, 1.23456789)

    tmptrack <- paste0("test.precision_", sample(1:1e9, 1))
    gtrack.create_sparse(tmptrack, "Test precision", test_intervals, test_values)
    withr::defer(gtrack.rm(tmptrack, force = TRUE))

    # Export
    gtrack.export_bed(tmptrack, tmpfile)

    bed_content <- read.table(tmpfile, sep = "\t", header = FALSE)

    # Values should be preserved with reasonable precision
    expect_equal(bed_content[, 5], test_values, tolerance = 1e-6)
})

test_that("gtrack.export_bed auto-generates interval names", {
    tmpfile <- tempfile(fileext = ".bed")
    withr::defer(unlink(tmpfile))

    # Create small track
    test_intervals <- gintervals(1, c(1000, 2000, 3000), c(2000, 3000, 4000))
    test_values <- c(10, 20, 30)

    tmptrack <- paste0("test.names_", sample(1:1e9, 1))
    gtrack.create_sparse(tmptrack, "Test names", test_intervals, test_values)
    withr::defer(gtrack.rm(tmptrack, force = TRUE))

    # Export
    gtrack.export_bed(tmptrack, tmpfile)

    bed_content <- read.table(tmpfile, sep = "\t", header = FALSE, stringsAsFactors = FALSE)

    # Verify auto-generated names
    expect_equal(bed_content[, 4], c("interval_1", "interval_2", "interval_3"))
})
