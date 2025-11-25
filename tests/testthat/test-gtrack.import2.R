create_isolated_test_db()
test_that("import with attrs parameter - NULL attrs works", {
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

    # Import with NULL attrs (should work normally)
    expect_no_error(
        gtrack.import(tmptrack, "Test track", temp_file, binsize = 1, attrs = NULL)
    )

    # Verify track was created successfully
    expect_true(gtrack.exists(tmptrack))
})

test_that("import with attrs parameter - attributes don't interfere with default attributes", {
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

    # Import with custom attributes
    attrs <- c("author" = "test_user", "custom_attr" = "custom_value")
    gtrack.import(tmptrack, "Test description", temp_file, binsize = 1, attrs = attrs)

    # Verify both custom and default attributes exist
    expect_equal(gtrack.attr.get(tmptrack, "author"), "test_user")
    expect_equal(gtrack.attr.get(tmptrack, "custom_attr"), "custom_value")
    expect_equal(gtrack.attr.get(tmptrack, "description"), "Test description")
    expect_true(nchar(gtrack.attr.get(tmptrack, "created.by")) > 0) # should contain creation info
    expect_true(nchar(gtrack.attr.get(tmptrack, "created.date")) > 0) # should contain creation date
})

test_that("import BED as sparse track", {
    gtrack.rm("test.bedtrack", force = TRUE)
    withr::defer(gtrack.rm("test.bedtrack", force = TRUE))

    # Create a temporary BED file
    bed_file <- tempfile(fileext = ".bed")
    writeLines(c(
        "track name=example",
        "chr1\t0\t2\tname1\t7",
        "chr1\t3\t5\tname2" # no score -> defaults to 1
    ), bed_file)
    withr::defer(unlink(bed_file))

    # Import BED; binsize is ignored for BED and track is created as sparse
    expect_no_error(
        gtrack.import("test.bedtrack", "BED import track", bed_file, binsize = 0)
    )

    # Basic sanity: track exists and extraction works on a small interval
    expect_true(gtrack.exists("test.bedtrack"))
    r <- gextract("test.bedtrack", gintervals(1, 0, 6))
    expect_true(nrow(r) > 0)
})

test_that("import BED as dense track with binsize", {
    gtrack.rm("test.bedtrack.dense", force = TRUE)
    withr::defer(gtrack.rm("test.bedtrack.dense", force = TRUE))

    bed_file <- tempfile(fileext = ".bed")
    writeLines(c(
        "chr1\t0\t10\tseg1\t2",
        "chr1\t10\t20\tseg2\t4"
    ), bed_file)
    withr::defer(unlink(bed_file))

    expect_no_error(
        gtrack.import("test.bedtrack.dense", "BED dense track", bed_file, binsize = 5, defval = 0)
    )
    info <- gtrack.info("test.bedtrack.dense")
    expect_equal(as.numeric(info$bin.size), 5)
})

test_that("import tab-delimited with header chrom/start/end/value (sparse)", {
    tmptrack <- paste0("test.tsvtrack_", sample(1:1e9, 1))
    gtrack.rm(tmptrack, force = TRUE)
    withr::defer(gtrack.rm(tmptrack, force = TRUE))

    # Create a temporary TSV file with required header and one value column
    tsv_file <- tempfile(fileext = ".tsv")
    writeLines(c(
        paste(c("chrom", "start", "end", "value"), collapse = "\t"),
        paste(c("chr1", 0, 3, 2.5), collapse = "\t"),
        paste(c("chr1", 4, 6, 1.0), collapse = "\t")
    ), tsv_file)
    withr::defer(unlink(tsv_file))

    # Import as sparse (binsize = 0)
    expect_no_error(
        gtrack.import(tmptrack, "TSV import track", tsv_file, binsize = 0)
    )
    expect_true(gtrack.exists(tmptrack))
})

test_that("import tab-delimited with multiple value columns fails", {
    tmptrack <- paste0("test.tsvbad_", sample(1:1e9, 1))
    gtrack.rm(tmptrack, force = TRUE)
    withr::defer(gtrack.rm(tmptrack, force = TRUE))

    # Create a TSV with two value columns -> should fail
    tsv_bad <- tempfile(fileext = ".tsv")
    writeLines(c(
        paste(c("chrom", "start", "end", "v1", "v2"), collapse = "\t"),
        paste(c("chr1", 0, 3, 2.5, 1.1), collapse = "\t")
    ), tsv_bad)
    withr::defer(unlink(tsv_bad))

    expect_error(
        gtrack.import(tmptrack, "bad tsv", tsv_bad, binsize = 0),
        "More than one value column appears"
    )
})
