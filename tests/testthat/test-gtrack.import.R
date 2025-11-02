load_test_db()
test_that("import and extract from s_7_export.txt", {
    track_name <- random_track_name()
    intervs <- gscreen("test.fixedbin > 0.1", gintervals(c(1, 2)))
    withr::defer(gtrack.rm(track_name, force = TRUE))
    gtrack.rm(track_name, force = TRUE)
    gtrack.import_mappedseq(track_name, "", "/net/mraid20/export/tgdata/db/tgdb/misha_snapshot/input_files/s_7_export.txt", remove.dups = FALSE)
    r <- gextract(track_name, intervs, colnames = "test.tmptrack")
    expect_regression(r, "track.import_mappedseq.s_7_export")
})

test_that("import and extract from sample-small.sam", {
    track_name <- random_track_name()
    intervs <- gscreen("test.fixedbin > 0.1", gintervals(c(1, 2)))
    gtrack.rm(track_name, force = TRUE)
    withr::defer(gtrack.rm(track_name, force = TRUE))
    gtrack.import_mappedseq(track_name, "", "/net/mraid20/export/tgdata/db/tgdb/misha_snapshot/input_files/sample-small.sam", cols.order = NULL, remove.dups = FALSE)
    r <- gextract(track_name, intervs, colnames = "test.tmptrack")
    expect_regression(r, "track.import_mappedseq.sample_small_sam")
})

test_that("import with pileup and binsize from s_7_export.txt", {
    track_name <- random_track_name()
    intervs <- gscreen("test.fixedbin > 0.1", gintervals(c(1, 2)))
    gtrack.rm(track_name, force = TRUE)
    withr::defer(gtrack.rm(track_name, force = TRUE))
    gtrack.import_mappedseq(track_name, "", "/net/mraid20/export/tgdata/db/tgdb/misha_snapshot/input_files/s_7_export.txt", remove.dups = FALSE, pileup = 180, binsize = 50)
    r <- gextract(track_name, intervs, colnames = "test.tmptrack")
    expect_regression(r, "track.import_pileup_binsize")
})

test_that("import with gmax data size option", {
    track_name <- random_track_name()
    gtrack.rm(track_name, TRUE)
    withr::defer(gtrack.rm(track_name, force = TRUE))
    withr::with_options(list(gmax.data.size = 10000), {
        gtrack.2d.import(track_name, "aaa7", c("/net/mraid20/export/tgdata/db/tgdb/misha_snapshot/input_files/f4"))
    })
    r <- gextract(track_name, .misha$ALLGENOME, colnames = "test.tmptrack")
    expect_regression(r, "track.import_gmax_option")
})

test_that("import with attrs parameter - single attribute", {
    track_name <- random_track_name()
    gtrack.rm(track_name, force = TRUE)
    withr::defer(gtrack.rm(track_name, force = TRUE))

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
    gtrack.import(track_name, "Test track", temp_file, binsize = 1, attrs = attrs)

    # Verify the attribute was set
    expect_equal(gtrack.attr.get(track_name, "author"), "test_user")
})

test_that("import with attrs parameter - multiple attributes", {
    track_name <- random_track_name()
    gtrack.rm(track_name, force = TRUE)
    withr::defer(gtrack.rm(track_name, force = TRUE))

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
    gtrack.import(track_name, "Test track", temp_file, binsize = 1, attrs = attrs)

    # Verify all attributes were set
    expect_equal(gtrack.attr.get(track_name, "author"), "test_user")
    expect_equal(gtrack.attr.get(track_name, "version"), "1.0")
    expect_equal(gtrack.attr.get(track_name, "experiment"), "test_exp")
})

test_that("import with attrs parameter - list format", {
    track_name <- random_track_name()
    gtrack.rm(track_name, force = TRUE)
    withr::defer(gtrack.rm(track_name, force = TRUE))

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
    gtrack.import(track_name, "Test track", temp_file, binsize = 1, attrs = attrs)

    # Verify attributes were set
    expect_equal(gtrack.attr.get(track_name, "author"), "test_user")
    expect_equal(gtrack.attr.get(track_name, "version"), "2.0")
})

test_that("import with attrs parameter - error on unnamed attrs", {
    track_name <- random_track_name()
    gtrack.rm(track_name, force = TRUE)
    withr::defer(gtrack.rm(track_name, force = TRUE))

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
        gtrack.import(track_name, "Test track", temp_file, binsize = 1, attrs = attrs),
        "attrs must be a named vector or list"
    )
})

test_that("import with attrs parameter - error on partially unnamed attrs", {
    track_name <- random_track_name()
    gtrack.rm(track_name, force = TRUE)
    withr::defer(gtrack.rm(track_name, force = TRUE))

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
        gtrack.import(track_name, "Test track", temp_file, binsize = 1, attrs = attrs),
        "attrs must be a named vector or list"
    )
})

test_that("import with attrs parameter - NULL attrs works", {
    track_name <- random_track_name()
    gtrack.rm(track_name, force = TRUE)
    withr::defer(gtrack.rm(track_name, force = TRUE))

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
        gtrack.import(track_name, "Test track", temp_file, binsize = 1, attrs = NULL)
    )

    # Verify track was created successfully
    expect_true(gtrack.exists(track_name))
})

test_that("import with attrs parameter - attributes don't interfere with default attributes", {
    track_name <- random_track_name()
    gtrack.rm(track_name, force = TRUE)
    withr::defer(gtrack.rm(track_name, force = TRUE))

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
    gtrack.import(track_name, "Test description", temp_file, binsize = 1, attrs = attrs)

    # Verify both custom and default attributes exist
    expect_equal(gtrack.attr.get(track_name, "author"), "test_user")
    expect_equal(gtrack.attr.get(track_name, "custom_attr"), "custom_value")
    expect_equal(gtrack.attr.get(track_name, "description"), "Test description")
    expect_true(nchar(gtrack.attr.get(track_name, "created.by")) > 0) # should contain creation info
    expect_true(nchar(gtrack.attr.get(track_name, "created.date")) > 0) # should contain creation date
})

test_that("import BED as sparse track", {
    gtrack.rm("temp.bedtrack", force = TRUE)
    withr::defer(gtrack.rm("temp.bedtrack", force = TRUE))

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
        gtrack.import("temp.bedtrack", "BED import track", bed_file, binsize = 0)
    )

    # Basic sanity: track exists and extraction works on a small interval
    expect_true(gtrack.exists("temp.bedtrack"))
    r <- gextract("temp.bedtrack", gintervals(1, 0, 6), colnames = "test.bedtrack")
    expect_true(nrow(r) > 0)
})

test_that("import BED as dense track with binsize", {
    gtrack.rm("temp.bedtrack.dense", force = TRUE)
    withr::defer(gtrack.rm("temp.bedtrack.dense", force = TRUE))

    bed_file <- tempfile(fileext = ".bed")
    writeLines(c(
        "chr1\t0\t10\tseg1\t2",
        "chr1\t10\t20\tseg2\t4"
    ), bed_file)
    withr::defer(unlink(bed_file))

    gdir.create("temp/bedtrack")
    expect_no_error(
        gtrack.import("temp.bedtrack.dense", "BED dense track", bed_file, binsize = 5, defval = 0)
    )
    info <- gtrack.info("temp.bedtrack.dense")
    expect_equal(as.numeric(info$bin.size), 5)
})

test_that("import tab-delimited with header chrom/start/end/value (sparse)", {
    gtrack.rm("temp.tsvtrack", force = TRUE)
    withr::defer(gtrack.rm("temp.tsvtrack", force = TRUE))

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
        gtrack.import("temp.tsvtrack", "TSV import track", tsv_file, binsize = 0)
    )
    expect_true(gtrack.exists("temp.tsvtrack"))
})

test_that("import tab-delimited with multiple value columns fails", {
    gtrack.rm("temp.tsvbad", force = TRUE)
    withr::defer(gtrack.rm("temp.tsvbad", force = TRUE))

    # Create a TSV with two value columns -> should fail
    tsv_bad <- tempfile(fileext = ".tsv")
    writeLines(c(
        paste(c("chrom", "start", "end", "v1", "v2"), collapse = "\t"),
        paste(c("chr1", 0, 3, 2.5, 1.1), collapse = "\t")
    ), tsv_bad)
    withr::defer(unlink(tsv_bad))

    expect_error(
        gtrack.import("temp.tsvbad", "bad tsv", tsv_bad, binsize = 0),
        "More than one value column appears"
    )
})
