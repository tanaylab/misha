create_isolated_test_db()

# Reported by a user who built intervals with gintervals() (which sorts) while
# keeping values in the pre-sort order. create_sparse itself binds values by the
# original row index, so unsorted input is fine.
test_that("create_sparse binds values by original row order, not sorted order", {
    tmptrack <- paste0("test.tmptrack_", sample(1:1e9, 1))
    gtrack.rm(tmptrack, force = TRUE)
    withr::defer(gtrack.rm(tmptrack, force = TRUE))

    # deliberately out of canonical chrom order
    intervals <- data.frame(
        chrom = c("chr2", "chr1", "chr2", "chr1"),
        start = c(100, 200, 300, 400),
        end = c(101, 201, 301, 401),
        stringsAsFactors = FALSE
    )
    values <- c(2.1, 1.1, 2.3, 1.4)

    gtrack.create_sparse(tmptrack, "Test", intervals, values)

    r <- gextract(tmptrack, gintervals.all(), colnames = "value")
    expected <- intervals[order(intervals$chrom, intervals$start), ]
    expect_equal(as.character(r$chrom), expected$chrom)
    expect_equal(r$start, expected$start)
    expect_equal(r$value, values[order(intervals$chrom, intervals$start)], tolerance = 1e-6)
})

test_that("create_sparse takes values from the value column when values is omitted", {
    tmptrack <- paste0("test.tmptrack_", sample(1:1e9, 1))
    gtrack.rm(tmptrack, force = TRUE)
    withr::defer(gtrack.rm(tmptrack, force = TRUE))

    intervals <- data.frame(
        chrom = c("chr2", "chr1", "chr2"),
        start = c(100, 200, 300),
        end = c(101, 201, 301),
        value = c(2.1, 1.1, 2.3),
        stringsAsFactors = FALSE
    )

    gtrack.create_sparse(tmptrack, "Test", intervals)

    r <- gextract(tmptrack, gintervals.all(), colnames = "value")
    expect_equal(r$value, c(1.1, 2.1, 2.3), tolerance = 1e-6)
})

test_that("create_sparse errors when values is omitted and there is no value column", {
    intervals <- gintervals(1, 100, 200)
    expect_error(
        gtrack.create_sparse("test.no_value_col", "Test", intervals),
        "no \"value\" column"
    )
})
