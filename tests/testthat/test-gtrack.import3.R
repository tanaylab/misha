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

write_2d_intervals_file <- function(df, env = parent.frame()) {
    src <- tempfile(fileext = ".txt")
    withr::defer(unlink(src), envir = env)
    write.table(df, src, sep = "\t", row.names = FALSE, quote = FALSE)
    src
}

test_that("gtrack.2d.import rejects a negative start2 (rectangles)", {
    tmptrack <- paste0("test.tmptrack_", sample(1:1e9, 1))
    withr::defer(gtrack.rm(tmptrack, force = TRUE))

    src <- write_2d_intervals_file(data.frame(
        chrom1 = "chr1", start1 = 100, end1 = 200,
        chrom2 = "chr2", start2 = -5, end2 = 10, value = 1
    ))
    expect_error(gtrack.2d.import(tmptrack, "negative start2", src), "invalid format of start2 coordinate")
    expect_false(gtrack.exists(tmptrack))
})

test_that("gtrack.2d.import rejects a negative start2 (points)", {
    tmptrack <- paste0("test.tmptrack_", sample(1:1e9, 1))
    withr::defer(gtrack.rm(tmptrack, force = TRUE))

    src <- write_2d_intervals_file(data.frame(
        chrom1 = "chr1", start1 = 100, end1 = 101,
        chrom2 = "chr2", start2 = -1, end2 = 0, value = 1
    ))
    expect_error(gtrack.2d.import(tmptrack, "negative start2", src), "invalid format of start2 coordinate")
    expect_false(gtrack.exists(tmptrack))
})

test_that("gtrack.2d.import_contacts rejects a negative start2", {
    tmptrack <- paste0("test.tmptrack_", sample(1:1e9, 1))
    withr::defer(gtrack.rm(tmptrack, force = TRUE))

    src <- write_2d_intervals_file(data.frame(
        chrom1 = "chr1", start1 = 100, end1 = 101,
        chrom2 = "chr2", start2 = -5, end2 = 10, value = 1
    ))
    expect_error(gtrack.2d.import_contacts(tmptrack, "negative start2", src), "invalid format of start2 coordinate")
    expect_false(gtrack.exists(tmptrack))
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
