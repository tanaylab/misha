load_test_db()

test_that("gintervals.export_bed exports basic BED3 format", {
    tmpfile <- tempfile(fileext = ".bed")
    withr::defer(unlink(tmpfile))

    # Create simple intervals
    test_intervals <- gintervals(c(1, 1, 2), c(1000, 5000, 2000), c(2000, 6000, 3000))

    # Export
    result <- gintervals.export_bed(test_intervals, tmpfile)

    # Should return file path invisibly
    expect_equal(result, tmpfile)
    expect_true(file.exists(tmpfile))

    # Read and verify content
    bed_content <- read.table(tmpfile, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
    expect_equal(nrow(bed_content), 3)
    expect_equal(ncol(bed_content), 3) # BED3

    # Verify coordinates (no conversion - 0-based preserved)
    expect_equal(bed_content[, 1], as.character(test_intervals$chrom))
    expect_equal(bed_content[, 2], test_intervals$start)
    expect_equal(bed_content[, 3], test_intervals$end)
})

test_that("gintervals.export_bed exports BED6 format with strand", {
    tmpfile <- tempfile(fileext = ".bed")
    withr::defer(unlink(tmpfile))

    # Create intervals with strand
    test_intervals <- gintervals(c(1, 1, 2), c(1000, 5000, 2000), c(2000, 6000, 3000))
    test_intervals$strand <- c(1, -1, 1)

    # Export
    gintervals.export_bed(test_intervals, tmpfile)

    # Read and verify content
    bed_content <- read.table(tmpfile, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
    expect_equal(nrow(bed_content), 3)
    expect_equal(ncol(bed_content), 6) # BED6

    # Verify strand conversion (1 -> "+", -1 -> "-")
    expect_equal(bed_content[, 6], c("+", "-", "+"))

    # Verify name and score placeholders
    expect_true(all(bed_content[, 4] == "."))
    expect_true(all(bed_content[, 5] == 0))
})

test_that("gintervals.export_bed adds track header when requested", {
    tmpfile <- tempfile(fileext = ".bed")
    withr::defer(unlink(tmpfile))

    test_intervals <- gintervals(1, c(1000, 2000), c(2000, 3000))

    # Export with track name only
    gintervals.export_bed(test_intervals, tmpfile, track_name = "MyRegions")

    lines <- readLines(tmpfile)
    expect_match(lines[1], 'track name="MyRegions"')
    expect_equal(length(lines), 3) # header + 2 data lines

    # Export with both track name and description
    tmpfile2 <- tempfile(fileext = ".bed")
    withr::defer(unlink(tmpfile2))

    gintervals.export_bed(test_intervals, tmpfile2,
        track_name = "MyRegions",
        description = "Test regions"
    )

    lines2 <- readLines(tmpfile2)
    expect_match(lines2[1], 'track name="MyRegions" description="Test regions"')
})

test_that("gintervals.export_bed works with saved intervals set", {
    tmpfile <- tempfile(fileext = ".bed")
    withr::defer(unlink(tmpfile))

    # Create and save intervals set
    test_intervals <- gintervals(c(1, 2), c(1000, 2000), c(2000, 3000))
    test_set_name <- paste0("test.export_bed_", sample(1:1e9, 1))

    gintervals.save(test_set_name, test_intervals)
    withr::defer(gintervals.rm(test_set_name, force = TRUE))

    # Export by name
    gintervals.export_bed(test_set_name, tmpfile)

    # Verify export
    bed_content <- read.table(tmpfile, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
    expect_equal(nrow(bed_content), 2)
    expect_equal(bed_content[, 2], test_intervals$start)
    expect_equal(bed_content[, 3], test_intervals$end)
})

test_that("gintervals.export_bed validates intervals parameter", {
    tmpfile <- tempfile(fileext = ".bed")

    # Missing required columns
    bad_intervals <- data.frame(chrom = 1, start = 100)
    expect_error(
        gintervals.export_bed(bad_intervals, tmpfile),
        "must have chrom, start, and end"
    )

    # Not a data.frame or valid name
    expect_error(
        gintervals.export_bed(list(chrom = 1, start = 100, end = 200), tmpfile),
        "must be a data.frame"
    )

    # Non-existent intervals set
    expect_error(
        gintervals.export_bed("nonexistent_set_12345", tmpfile),
        "does not exist"
    )
})

test_that("gintervals.export_bed validates file parameter", {
    test_intervals <- gintervals(1, 1000, 2000)

    # Non-character file
    expect_error(
        gintervals.export_bed(test_intervals, 123),
        "must be a single character string"
    )

    # Multiple file paths
    expect_error(
        gintervals.export_bed(test_intervals, c("file1.bed", "file2.bed")),
        "must be a single character string"
    )

    # Non-existent directory
    expect_error(
        gintervals.export_bed(test_intervals, "/nonexistent_dir_12345/file.bed"),
        "Directory does not exist"
    )
})

test_that("gintervals.export_bed validates track_name and description", {
    tmpfile <- tempfile(fileext = ".bed")
    withr::defer(unlink(tmpfile))

    test_intervals <- gintervals(1, 1000, 2000)

    # Invalid track_name type
    expect_error(
        gintervals.export_bed(test_intervals, tmpfile, track_name = 123),
        "must be a single character string"
    )

    # Invalid description type
    expect_error(
        gintervals.export_bed(test_intervals, tmpfile, track_name = "test", description = 123),
        "must be a single character string"
    )
})

test_that("gintervals.export_bed rejects empty intervals", {
    tmpfile <- tempfile(fileext = ".bed")

    # gintervals with empty vectors returns NULL, which triggers data.frame check
    empty_intervals <- gintervals(integer(0), integer(0), integer(0))

    expect_error(
        gintervals.export_bed(empty_intervals, tmpfile),
        "must be a data.frame"
    )

    # Also test with an actual empty data.frame
    empty_df <- data.frame(chrom = integer(0), start = integer(0), end = integer(0))
    expect_error(
        gintervals.export_bed(empty_df, tmpfile),
        "Cannot export empty intervals"
    )
})

test_that("gintervals.export_bed rejects 2D intervals", {
    tmpfile <- tempfile(fileext = ".bed")

    # Create 2D intervals
    intervals_2d <- gintervals.2d(1, 1000, 2000, 1, 3000, 4000)

    expect_error(
        gintervals.export_bed(intervals_2d, tmpfile),
        "Cannot export 2D intervals"
    )
})

test_that("gintervals.export_bed preserves coordinate system", {
    tmpfile <- tempfile(fileext = ".bed")
    withr::defer(unlink(tmpfile))

    # Create intervals with specific coordinates to verify no conversion
    test_intervals <- gintervals(
        c(1, 1, 2),
        c(0, 100, 1000), # Including 0 to test edge case
        c(50, 200, 2000)
    )

    gintervals.export_bed(test_intervals, tmpfile)

    bed_content <- read.table(tmpfile, sep = "\t", header = FALSE)

    # Coordinates should be exactly as input (no +1 conversion like BigWig)
    expect_equal(bed_content[, 2], c(0, 100, 1000))
    expect_equal(bed_content[, 3], c(50, 200, 2000))
})

test_that("gintervals.export_bed handles multiple chromosomes", {
    tmpfile <- tempfile(fileext = ".bed")
    withr::defer(unlink(tmpfile))

    # Create intervals across multiple chromosomes
    test_intervals <- gintervals(
        c(1, 1, 2, 2, 3),
        c(1000, 5000, 2000, 6000, 1000),
        c(2000, 6000, 3000, 7000, 2000)
    )

    gintervals.export_bed(test_intervals, tmpfile)

    bed_content <- read.table(tmpfile, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
    expect_equal(nrow(bed_content), 5)

    # Verify chromosomes are preserved
    expect_equal(bed_content[, 1], c("chr1", "chr1", "chr2", "chr2", "chr3"))
})

test_that("gintervals.export_bed handles file overwrite", {
    tmpfile <- tempfile(fileext = ".bed")
    withr::defer(unlink(tmpfile))

    test_intervals1 <- gintervals(1, 1000, 2000)
    test_intervals2 <- gintervals(c(1, 2), c(3000, 4000), c(4000, 5000))

    # Write first file
    gintervals.export_bed(test_intervals1, tmpfile)
    content1 <- readLines(tmpfile)
    expect_equal(length(content1), 1)

    # Overwrite with second file
    gintervals.export_bed(test_intervals2, tmpfile)
    content2 <- readLines(tmpfile)
    expect_equal(length(content2), 2)
})

test_that("gintervals.export_bed round-trip preserves data", {
    tmpfile <- tempfile(fileext = ".bed")
    withr::defer(unlink(tmpfile))

    # Create test intervals
    original_intervals <- gintervals(
        c(1, 1, 2),
        c(1000, 5000, 2000),
        c(2000, 6000, 3000)
    )

    # Export to BED
    gintervals.export_bed(original_intervals, tmpfile)

    # Import back using gtrack.import
    tmptrack <- paste0("test.reimport_", sample(1:1e9, 1))
    withr::defer(gtrack.rm(tmptrack, force = TRUE))

    gtrack.import(tmptrack, "", tmpfile, binsize = 0)

    # Extract reimported data
    reimported <- gextract(tmptrack, gintervals.all())

    # Compare coordinates (should match exactly - no conversion)
    expect_equal(nrow(reimported), nrow(original_intervals))
    expect_equal(reimported$start, original_intervals$start)
    expect_equal(reimported$end, original_intervals$end)
    expect_equal(as.character(reimported$chrom), as.character(original_intervals$chrom))
})

test_that("gintervals.export_bed handles large intervals sets", {
    tmpfile <- tempfile(fileext = ".bed")
    withr::defer(unlink(tmpfile))

    # Create large intervals set (1000 intervals)
    n <- 1000
    large_intervals <- gintervals(
        rep(1, n),
        seq(0, (n - 1) * 1000, 1000),
        seq(500, n * 1000 - 500, 1000)
    )

    # Should complete without error
    expect_no_error(gintervals.export_bed(large_intervals, tmpfile))

    # Verify file size
    expect_true(file.exists(tmpfile))
    expect_gt(file.size(tmpfile), 0)

    # Verify row count
    bed_content <- read.table(tmpfile, sep = "\t", header = FALSE)
    expect_equal(nrow(bed_content), n)
})

test_that("gintervals.export_bed handles special chromosome names", {
    tmpfile <- tempfile(fileext = ".bed")
    withr::defer(unlink(tmpfile))

    # Create intervals with various chromosome representations
    test_intervals <- data.frame(
        chrom = c("chr1", "chr2", "chrX"),
        start = c(1000, 2000, 3000),
        end = c(2000, 3000, 4000),
        stringsAsFactors = FALSE
    )

    gintervals.export_bed(test_intervals, tmpfile)

    bed_content <- read.table(tmpfile, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
    expect_equal(bed_content[, 1], c("chr1", "chr2", "chrX"))
})
