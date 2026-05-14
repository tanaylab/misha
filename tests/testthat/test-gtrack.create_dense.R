create_isolated_test_db()

test_that("test for dense track from intervals and values", {
    tmptrack <- paste0("test.tmptrack_", sample(1:1e9, 1))
    gtrack.rm(tmptrack, force = TRUE)
    withr::defer(gtrack.rm(tmptrack, force = TRUE))

    # Create intervals and values
    intervals <- gintervals(
        chrom = c(1, 1, 2),
        start = c(100, 500, 200),
        end = c(200, 600, 300)
    )
    values <- c(1.5, 2.7, 3.1)

    gtrack.create_dense(tmptrack, "Test dense track", intervals, values, 50, 0)

    # Extract data and check results
    r <- gextract(tmptrack, gintervals(c(1, 2), 0, 1000), colnames = "value")

    # Verify data is present
    expect_true(nrow(r) > 0)
    expect_true(any(r$value > 0))

    # Check track info - note that type is lowercase "dense" not "Dense"
    info <- gtrack.info(tmptrack)
    expect_equal(info$type, "dense")
    expect_equal(info$bin.size, 50)
})

test_that("test for overlapping intervals in dense track", {
    tmptrack <- paste0("test.tmptrack_", sample(1:1e9, 1))
    gtrack.rm(tmptrack, force = TRUE)
    withr::defer(gtrack.rm(tmptrack, force = TRUE))

    # Create intervals with overlapping regions
    intervals <- gintervals(
        chrom = c(1, 1, 1),
        start = c(100, 150, 200),
        end = c(250, 300, 350)
    )
    values <- c(1.0, 2.0, 3.0)

    gtrack.create_dense(tmptrack, "Test overlapping intervals", intervals, values, 100, 0)

    # Extract data and check results
    r <- gextract(tmptrack, gintervals(1, 0, 400), colnames = "value")

    # Verify overlapping intervals were processed
    expect_true(nrow(r) > 0)

    # Check if any value between overlapping intervals reflects weighted average
    # For overlapping regions, values should be between 1.0 and 3.0
    expect_true(any(r$value > 1.0 & r$value < 3.0))
})

test_that("test for values length not matching intervals", {
    tmptrack <- paste0("test.tmptrack_", sample(1:1e9, 1))
    # Create intervals with more intervals than values
    intervals <- gintervals(
        chrom = c(1, 1, 1),
        start = c(100, 200, 300),
        end = c(150, 250, 350)
    )
    values <- c(1.0, 2.0) # Only 2 values for 3 intervals

    expect_error(
        gtrack.create_dense(tmptrack, "Test mismatched lengths", intervals, values, 50, 0),
        "Length of values must match the number of intervals"
    )
})

test_that("test for default value in dense track", {
    tmptrack <- paste0("test.tmptrack_", sample(1:1e9, 1))
    gtrack.rm(tmptrack, force = TRUE)
    withr::defer(gtrack.rm(tmptrack, force = TRUE))

    # Create intervals with a gap
    intervals <- gintervals(
        chrom = c(1, 1),
        start = c(100, 300),
        end = c(150, 350)
    )
    values <- c(1.0, 2.0)

    # Create track with default value NaN
    defval <- NaN
    gtrack.create_dense(tmptrack, "Test default value", intervals, values, 50, defval)

    # Extract data and check results
    r <- gextract(tmptrack, gintervals(1, 0, 400), colnames = "value")

    # Check if there are regions with NaN (empty regions)
    expect_true(any(is.nan(r$value)))

    # Check that we have both data-containing regions and default regions
    expect_true(any(!is.nan(r$value)))
})

test_that("test for recreating test.fixedbin from extracted data", {
    tmptrack <- paste0("test.tmptrack_", sample(1:1e9, 1))
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
    gtrack.rm(tmptrack, force = TRUE)
    withr::defer(gtrack.rm(tmptrack, force = TRUE))

    # Get track info to determine bin size
    info <- gtrack.info("test.fixedbin")
    bin_size <- info$bin.size

    # Create a new dense track with the extracted data
    gtrack.create_dense(tmptrack, "Recreated test.fixedbin", intervals, values, bin_size, 0)

    # Extract data from the recreated track
    recreated_data <- gextract(tmptrack, gintervals.all(), colnames = "value")

    # Compare the original and recreated data
    expect_equal(nrow(extracted_data), nrow(recreated_data))

    # Check track info
    new_info <- gtrack.info(tmptrack)
    expect_equal(new_info$type, "dense")
    expect_equal(new_info$bin.size, bin_size)

    expect_equal(gtrack.attr.get(tmptrack, "description"), "Recreated test.fixedbin")
    expect_equal(gtrack.attr.get(tmptrack, "binsize"), as.character(bin_size))
    expect_equal(gtrack.attr.get(tmptrack, "type"), "dense")
    expect_match(gtrack.attr.get(tmptrack, "created.by"), "gtrack.create_dense.*")
})

# ----------------------------------------------------------------------------
# func= parameter: per-bin aggregation
# ----------------------------------------------------------------------------

# Hand-computable fixture used across the func tests.
#   chr 1, binsize 20, three intervals on [0,100):
#     A: [ 0,  60)  value = 1
#     B: [40, 100)  value = 3
#     C: [50,  80)  value = 5
#
#   Bin coverage:
#     bin [ 0, 20):  A only           overlap A=20
#     bin [20, 40):  A only           overlap A=20
#     bin [40, 60):  A+B+C            overlap A=20, B=20, C=10
#     bin [60, 80):  B+C              overlap B=20, C=20
#     bin [80,100):  B only           overlap B=20
.func_fixture_intervals <- function() {
    gintervals(
        chrom = c(1, 1, 1),
        start = c(0, 40, 50),
        end = c(60, 100, 80)
    )
}
.func_fixture_values <- c(1, 3, 5)

# Helper: build a dense track with the fixture, extract on [0,100) as a vector
.func_make_and_extract <- function(tmptrack, func, defval = NaN) {
    gtrack.rm(tmptrack, force = TRUE)
    gtrack.create_dense(tmptrack, sprintf("func=%s", func),
        .func_fixture_intervals(), .func_fixture_values,
        binsize = 20, defval = defval, func = func
    )
    r <- gextract(tmptrack, gintervals(1, 0, 100), colnames = "value")
    r <- r[order(r$start), ]
    r$value
}

test_that("func='weighted.mean' is default and matches historical behavior", {
    tmptrack <- paste0("test.tmptrack_", sample(1:1e9, 1))
    withr::defer(gtrack.rm(tmptrack, force = TRUE))

    v_default <- {
        gtrack.rm(tmptrack, force = TRUE)
        gtrack.create_dense(tmptrack, "default",
            .func_fixture_intervals(), .func_fixture_values,
            binsize = 20, defval = NaN
        )
        r <- gextract(tmptrack, gintervals(1, 0, 100), colnames = "value")
        r <- r[order(r$start), ]
        r$value
    }
    v_explicit <- .func_make_and_extract(tmptrack, "weighted.mean")
    expect_equal(v_explicit, v_default)
    # Hand-computed per-bin: bin2 = (1*20 + 3*20 + 5*10) / (20+20+10) = 130/50 = 2.6
    #                        bin3 = (3*20 + 5*20) / 40 = 4
    expect_equal(v_default, c(1, 1, 2.6, 4, 3))
})

test_that("func='weighted.sum' returns sum(v*overlap) per bin", {
    tmptrack <- paste0("test.tmptrack_", sample(1:1e9, 1))
    withr::defer(gtrack.rm(tmptrack, force = TRUE))
    v <- .func_make_and_extract(tmptrack, "weighted.sum")
    # bin0:  1*20                 = 20
    # bin1:  1*20                 = 20
    # bin2:  1*20 + 3*20 + 5*10   = 130
    # bin3:  3*20 + 5*20          = 160
    # bin4:  3*20                 = 60
    expect_equal(v, c(20, 20, 130, 160, 60))
})

test_that("func='max' returns unweighted max of overlapping interval values", {
    tmptrack <- paste0("test.tmptrack_", sample(1:1e9, 1))
    withr::defer(gtrack.rm(tmptrack, force = TRUE))
    v <- .func_make_and_extract(tmptrack, "max")
    # bin2 sees A,B,C -> max(1,3,5) = 5
    expect_equal(v, c(1, 1, 5, 5, 3))
})

test_that("func='min' returns unweighted min of overlapping interval values", {
    tmptrack <- paste0("test.tmptrack_", sample(1:1e9, 1))
    withr::defer(gtrack.rm(tmptrack, force = TRUE))
    v <- .func_make_and_extract(tmptrack, "min")
    expect_equal(v, c(1, 1, 1, 3, 3))
})

test_that("func='count' returns number of overlapping intervals per bin", {
    tmptrack <- paste0("test.tmptrack_", sample(1:1e9, 1))
    withr::defer(gtrack.rm(tmptrack, force = TRUE))
    v <- .func_make_and_extract(tmptrack, "count")
    # bin2 sees A,B,C; bin3 sees B,C
    expect_equal(v, c(1, 1, 3, 2, 1))
})

test_that("func='count' empty bin returns 0 regardless of defval", {
    tmptrack <- paste0("test.tmptrack_", sample(1:1e9, 1))
    withr::defer(gtrack.rm(tmptrack, force = TRUE))
    # Place a single tiny interval on chr1 [0,5), binsize 20, extract [0,40)
    gtrack.rm(tmptrack, force = TRUE)
    intervs <- gintervals(chrom = 1, start = 0, end = 5)
    gtrack.create_dense(tmptrack, "count empty", intervs,
        values = 7, binsize = 20,
        defval = 99, func = "count"
    )
    r <- gextract(tmptrack, gintervals(1, 0, 40), colnames = "value")
    r <- r[order(r$start), ]
    expect_equal(r$value, c(1, 0))
})

test_that("func='median' returns overlap-weighted lower median", {
    tmptrack <- paste0("test.tmptrack_", sample(1:1e9, 1))
    withr::defer(gtrack.rm(tmptrack, force = TRUE))
    v <- .func_make_and_extract(tmptrack, "median")
    # bin0,1: only A (v=1)  -> 1
    # bin2:   A=20, B=20, C=10 -> sorted (1@20, 3@20, 5@10). total=50, half=25.
    #         acc 20 < 25, acc 40 >= 25 at value=3 -> 3
    # bin3:   B=20, C=20    -> sorted (3@20, 5@20). half=20. acc 20 >= 20 at 3 -> 3
    # bin4:   only B (v=3)  -> 3
    expect_equal(v, c(1, 1, 3, 3, 3))
})

test_that("defval matrix: weighted.sum with defval=0 vs NaN", {
    tmptrack <- paste0("test.tmptrack_", sample(1:1e9, 1))
    withr::defer(gtrack.rm(tmptrack, force = TRUE))

    # One interval [10, 30) value=5, binsize 20, bin [20,40) is partially covered.
    intervs <- gintervals(chrom = 1, start = 10, end = 30)

    gtrack.rm(tmptrack, force = TRUE)
    gtrack.create_dense(tmptrack, "ws nan", intervs, 5, binsize = 20, defval = NaN, func = "weighted.sum")
    r_nan <- gextract(tmptrack, gintervals(1, 0, 60), colnames = "value")
    r_nan <- r_nan[order(r_nan$start), ]
    # bin [ 0,20): overlap=10, sum=50
    # bin [20,40): overlap=10, sum=50
    # bin [40,60): empty, defval=NaN -> NaN
    expect_equal(r_nan$value[1:2], c(50, 50))
    expect_true(is.nan(r_nan$value[3]))

    gtrack.rm(tmptrack, force = TRUE)
    gtrack.create_dense(tmptrack, "ws 7", intervs, 5, binsize = 20, defval = 7, func = "weighted.sum")
    r_def <- gextract(tmptrack, gintervals(1, 0, 60), colnames = "value")
    r_def <- r_def[order(r_def$start), ]
    # bin [ 0,20): overlap=10, sum = 5*10 + 7*10 = 120
    # bin [20,40): overlap=10, sum = 5*10 + 7*10 = 120
    # bin [40,60): empty, defval=7 -> 7*20 = 140
    expect_equal(r_def$value, c(120, 120, 140))
})

test_that("defval matrix: max/min include defval synthetic on partial bins", {
    tmptrack <- paste0("test.tmptrack_", sample(1:1e9, 1))
    withr::defer(gtrack.rm(tmptrack, force = TRUE))

    intervs <- gintervals(chrom = 1, start = 10, end = 30)

    # max with defval=99 -> partial bins should reach 99
    gtrack.rm(tmptrack, force = TRUE)
    gtrack.create_dense(tmptrack, "max defval", intervs, 5, binsize = 20, defval = 99, func = "max")
    r <- gextract(tmptrack, gintervals(1, 0, 60), colnames = "value")
    r <- r[order(r$start), ]
    expect_equal(r$value, c(99, 99, 99)) # partial bins see 99; empty bin = 99

    # min with defval=-99 -> partial bins should reach -99
    gtrack.rm(tmptrack, force = TRUE)
    gtrack.create_dense(tmptrack, "min defval", intervs, 5, binsize = 20, defval = -99, func = "min")
    r <- gextract(tmptrack, gintervals(1, 0, 60), colnames = "value")
    r <- r[order(r$start), ]
    expect_equal(r$value, c(-99, -99, -99))
})

test_that("invalid func errors with allowed-values message", {
    tmptrack <- paste0("test.tmptrack_", sample(1:1e9, 1))
    intervs <- gintervals(chrom = 1, start = 0, end = 50)
    expect_error(
        gtrack.create_dense(tmptrack, "bad", intervs, 1, binsize = 20, func = "bogus"),
        "Invalid 'func'"
    )
})

test_that("NaN interval values are skipped regardless of func", {
    tmptrack <- paste0("test.tmptrack_", sample(1:1e9, 1))
    withr::defer(gtrack.rm(tmptrack, force = TRUE))

    intervs <- gintervals(chrom = c(1, 1), start = c(0, 20), end = c(40, 40))
    vals <- c(NaN, 7) # first interval has NaN value, should be skipped

    gtrack.rm(tmptrack, force = TRUE)
    gtrack.create_dense(tmptrack, "nan vals", intervs, vals, binsize = 20, defval = NaN, func = "count")
    r <- gextract(tmptrack, gintervals(1, 0, 40), colnames = "value")
    r <- r[order(r$start), ]
    expect_equal(r$value, c(0, 1)) # first bin has no real interval (NaN skipped)
})
