create_isolated_test_db()

# Tests to verify that value-based vtracks produce identical results to the original tracks
# they're derived from

test_that("value-based vtrack matches sparse track with avg function", {
    # Extract data from sparse track - get actual intervals with data, not whole chromosomes
    intervals_with_data <- gscreen("!is.na(test.sparse)", gintervals.all())
    extracted <- gextract("test.sparse", intervals = intervals_with_data, iterator = intervals_with_data)

    # Create value-based vtrack from extracted data
    gvtrack.create("test.sparse.vt", src = extracted, func = "avg")

    # Test with single interval
    test_int <- gintervals("chr1", 3000, 8000)
    track_result <- gextract("test.sparse", intervals = test_int, iterator = test_int, colnames = "value")
    vtrack_result <- gextract("test.sparse.vt", intervals = test_int, iterator = test_int, colnames = "value")

    expect_equal(track_result$value, vtrack_result$value, tolerance = 1e-6)

    # Test with iterator
    test_region <- gintervals("chr1", 0, 10000)
    track_iter <- gextract("test.sparse", intervals = test_region, iterator = 1000, colnames = "value")
    vtrack_iter <- gextract("test.sparse.vt", intervals = test_region, iterator = 1000, colnames = "value")

    expect_equal(track_iter$value, vtrack_iter$value, tolerance = 1e-6)

    # Test with multiple non-overlapping intervals
    test_ints <- gintervals(
        c("chr1", "chr1", "chr2"),
        c(1000, 5000, 2000),
        c(2000, 6000, 4000)
    )
    track_multi <- gextract("test.sparse", intervals = test_ints, iterator = test_ints, colnames = "value")
    vtrack_multi <- gextract("test.sparse.vt", intervals = test_ints, iterator = test_ints, colnames = "value")

    expect_equal(track_multi$value, vtrack_multi$value, tolerance = 1e-6)

    gvtrack.rm("test.sparse.vt")
})

test_that("value-based vtrack matches track-based vtrack for all aggregation functions", {
    # Extract data from sparse track
    intervals_with_data <- gscreen("!is.na(test.sparse)", gintervals.all())
    extracted <- gextract("test.sparse", intervals = intervals_with_data, iterator = intervals_with_data)

    # Test all aggregation functions by comparing track-based vtrack vs value-based vtrack
    # This is the correct comparison since track expressions like min(test.sparse) don't work in misha
    functions <- c("avg", "min", "max", "sum", "stddev", "quantile")

    test_int <- gintervals("chr1", 3000, 8000)

    for (func in functions) {
        # Create track-based vtrack (from original track)
        track_vtrack_name <- paste0("test.track.", func)
        gvtrack.create(track_vtrack_name, "test.sparse",
            func = func,
            params = if (func == "quantile") 0.5 else NULL
        )

        # Create value-based vtrack (from extracted data)
        value_vtrack_name <- paste0("test.value.", func)
        gvtrack.create(value_vtrack_name,
            src = extracted, func = func,
            params = if (func == "quantile") 0.5 else NULL
        )

        # Extract from both
        track_vtrack_result <- gextract(track_vtrack_name,
            intervals = test_int,
            iterator = test_int, colnames = "value"
        )
        value_vtrack_result <- gextract(value_vtrack_name,
            intervals = test_int,
            iterator = test_int, colnames = "value"
        )

        # They should match
        expect_equal(track_vtrack_result$value, value_vtrack_result$value,
            tolerance = 1e-6, info = paste("Function:", func)
        )

        # Clean up
        gvtrack.rm(track_vtrack_name)
        gvtrack.rm(value_vtrack_name)
    }
})

test_that("value-based vtrack works with gscreen like original track", {
    intervals_with_data <- gscreen("!is.na(test.sparse)", gintervals.all())
    extracted <- gextract("test.sparse", intervals = intervals_with_data, iterator = intervals_with_data)
    gvtrack.create("test.sparse.screen.vt", src = extracted, func = "avg")

    # Test gscreen with threshold
    test_region <- gintervals("chr1", 0, 10000)
    track_screen <- gscreen("test.sparse > 0.5", intervals = test_region, iterator = test_region)
    vtrack_screen <- gscreen("test.sparse.screen.vt > 0.5", intervals = test_region, iterator = test_region)

    # Both should return the same intervals (or both NULL)
    if (is.null(track_screen) && is.null(vtrack_screen)) {
        expect_true(TRUE)
    } else {
        expect_equal(nrow(track_screen), nrow(vtrack_screen))
        if (nrow(track_screen) > 0) {
            expect_equal(track_screen$start, vtrack_screen$start)
            expect_equal(track_screen$end, vtrack_screen$end)
        }
    }

    gvtrack.rm("test.sparse.screen.vt")
})

test_that("value-based vtrack works with track expressions", {
    intervals_with_data <- gscreen("!is.na(test.sparse)", gintervals.all())
    extracted <- gextract("test.sparse", intervals = intervals_with_data, iterator = intervals_with_data)
    gvtrack.create("test.expr.track.vt", src = extracted, func = "avg")

    test_int <- gintervals("chr1", 3000, 8000)

    # Test multiplication
    track_expr <- gextract("test.sparse * 2", intervals = test_int, iterator = test_int, colnames = "value")
    vtrack_expr <- gextract("test.expr.track.vt * 2", intervals = test_int, iterator = test_int, colnames = "value")
    expect_equal(track_expr$value, vtrack_expr$value, tolerance = 1e-6)

    # Test addition
    track_add <- gextract("test.sparse + 10", intervals = test_int, iterator = test_int, colnames = "value")
    vtrack_add <- gextract("test.expr.track.vt + 10", intervals = test_int, iterator = test_int, colnames = "value")
    expect_equal(track_add$value, vtrack_add$value, tolerance = 1e-6)

    gvtrack.rm("test.expr.track.vt")
})

test_that("value-based vtrack matches track-based vtrack for position and selector functions", {
    # Extract data from sparse track
    intervals_with_data <- gscreen("!is.na(test.sparse)", gintervals.all())
    extracted <- gextract("test.sparse", intervals = intervals_with_data, iterator = intervals_with_data)

    # Test position and selector functions
    functions <- c(
        "nearest", "exists", "size", "first", "last",
        "first.pos.abs", "last.pos.abs",
        "min.pos.abs", "max.pos.abs"
    )

    test_int <- gintervals("chr1", 3000, 8000)

    for (func in functions) {
        # Create track-based vtrack (from original track)
        track_vtrack_name <- paste0("test.track.", gsub("\\.", "_", func))
        gvtrack.create(track_vtrack_name, "test.sparse", func = func)

        # Create value-based vtrack (from extracted data)
        value_vtrack_name <- paste0("test.value.", gsub("\\.", "_", func))
        gvtrack.create(value_vtrack_name, src = extracted, func = func)

        # Extract from both
        track_vtrack_result <- gextract(track_vtrack_name,
            intervals = test_int,
            iterator = test_int, colnames = "value"
        )
        value_vtrack_result <- gextract(value_vtrack_name,
            intervals = test_int,
            iterator = test_int, colnames = "value"
        )

        # They should match
        expect_equal(track_vtrack_result$value, value_vtrack_result$value,
            tolerance = 1e-6, info = paste("Function:", func)
        )

        # Clean up
        gvtrack.rm(track_vtrack_name)
        gvtrack.rm(value_vtrack_name)
    }
})

test_that("value-based vtrack with iterator matches track-based vtrack", {
    # Extract data from sparse track
    intervals_with_data <- gscreen("!is.na(test.sparse)", gintervals.all())
    extracted <- gextract("test.sparse", intervals = intervals_with_data, iterator = intervals_with_data)

    # Test multiple functions with iterator
    functions <- c("avg", "min", "max", "sum", "stddev", "size", "first", "last")

    test_region <- gintervals("chr1", 0, 20000)

    for (func in functions) {
        # Create track-based vtrack
        track_vtrack_name <- paste0("test.iter.track.", gsub("\\.", "_", func))
        gvtrack.create(track_vtrack_name, "test.sparse", func = func)

        # Create value-based vtrack
        value_vtrack_name <- paste0("test.iter.value.", gsub("\\.", "_", func))
        gvtrack.create(value_vtrack_name, src = extracted, func = func)

        # Extract with 1000bp iterator
        track_result <- gextract(track_vtrack_name, intervals = test_region, iterator = 1000)
        value_result <- gextract(value_vtrack_name, intervals = test_region, iterator = 1000)

        # Compare all values
        expect_equal(nrow(track_result), nrow(value_result), info = paste("Function:", func, "- row count"))
        expect_equal(track_result[[track_vtrack_name]], value_result[[value_vtrack_name]],
            tolerance = 1e-6, info = paste("Function:", func, "- with iterator")
        )

        # Clean up
        gvtrack.rm(track_vtrack_name)
        gvtrack.rm(value_vtrack_name)
    }
})

test_that("value-based vtrack with multiple chromosomes matches track-based vtrack", {
    # Extract data from sparse track across multiple chromosomes
    intervals_with_data <- gscreen("!is.na(test.sparse)", gintervals.all())
    extracted <- gextract("test.sparse", intervals = intervals_with_data, iterator = intervals_with_data)

    # Test with intervals spanning multiple chromosomes
    test_ints <- gintervals(
        c("chr1", "chr1", "chr2", "chr2", "chr3"),
        c(0, 5000, 0, 10000, 0),
        c(5000, 10000, 10000, 20000, 5000)
    )

    functions <- c("avg", "min", "max", "sum", "first", "last", "size")

    for (func in functions) {
        # Create track-based vtrack
        track_vtrack_name <- paste0("test.multi.track.", gsub("\\.", "_", func))
        gvtrack.create(track_vtrack_name, "test.sparse", func = func)

        # Create value-based vtrack
        value_vtrack_name <- paste0("test.multi.value.", gsub("\\.", "_", func))
        gvtrack.create(value_vtrack_name, src = extracted, func = func)

        # Extract from multiple intervals
        track_result <- gextract(track_vtrack_name, intervals = test_ints, iterator = test_ints)
        value_result <- gextract(value_vtrack_name, intervals = test_ints, iterator = test_ints)

        # Compare all values
        expect_equal(nrow(track_result), nrow(value_result), info = paste("Function:", func))
        expect_equal(track_result[[track_vtrack_name]], value_result[[value_vtrack_name]],
            tolerance = 1e-6, info = paste("Function:", func, "- multiple chromosomes")
        )

        # Clean up
        gvtrack.rm(track_vtrack_name)
        gvtrack.rm(value_vtrack_name)
    }
})

test_that("value-based vtrack with quantile matches track-based vtrack for different percentiles", {
    # Extract data from sparse track
    intervals_with_data <- gscreen("!is.na(test.sparse)", gintervals.all())
    extracted <- gextract("test.sparse", intervals = intervals_with_data, iterator = intervals_with_data)

    test_int <- gintervals("chr1", 0, 20000)

    # Test multiple percentiles
    percentiles <- c(0.25, 0.5, 0.75, 0.9)

    for (pct in percentiles) {
        # Create track-based vtrack
        track_vtrack_name <- paste0("test.q.track.", gsub("\\.", "_", as.character(pct)))
        gvtrack.create(track_vtrack_name, "test.sparse", func = "quantile", params = pct)

        # Create value-based vtrack
        value_vtrack_name <- paste0("test.q.value.", gsub("\\.", "_", as.character(pct)))
        gvtrack.create(value_vtrack_name, src = extracted, func = "quantile", params = pct)

        # Extract
        track_result <- gextract(track_vtrack_name, intervals = test_int, iterator = test_int)
        value_result <- gextract(value_vtrack_name, intervals = test_int, iterator = test_int)

        # Compare
        expect_equal(track_result[[track_vtrack_name]], value_result[[value_vtrack_name]],
            tolerance = 1e-6, info = paste("Percentile:", pct)
        )

        # Clean up
        gvtrack.rm(track_vtrack_name)
        gvtrack.rm(value_vtrack_name)
    }
})

test_that("value-based vtrack edge cases match track-based vtrack", {
    # Extract data from sparse track
    intervals_with_data <- gscreen("!is.na(test.sparse)", gintervals.all())
    extracted <- gextract("test.sparse", intervals = intervals_with_data, iterator = intervals_with_data)

    # Create vtracks
    gvtrack.create("test.edge.track", "test.sparse", func = "avg")
    gvtrack.create("test.edge.value", src = extracted, func = "avg")

    # Test 1: Empty interval (no data)
    empty_int <- gintervals("chr1", 25000, 26000)
    track_empty <- gextract("test.edge.track", intervals = empty_int, iterator = empty_int)
    value_empty <- gextract("test.edge.value", intervals = empty_int, iterator = empty_int)
    expect_equal(is.nan(track_empty$test.edge.track), is.nan(value_empty$test.edge.value))

    # Test 2: Very small interval
    tiny_int <- gintervals("chr1", 100, 110)
    track_tiny <- gextract("test.edge.track", intervals = tiny_int, iterator = tiny_int)
    value_tiny <- gextract("test.edge.value", intervals = tiny_int, iterator = tiny_int)
    expect_equal(track_tiny$test.edge.track, value_tiny$test.edge.value, tolerance = 1e-6)

    # Test 3: Single data point coverage
    single_int <- gintervals("chr1", 0, 100)
    track_single <- gextract("test.edge.track", intervals = single_int, iterator = single_int)
    value_single <- gextract("test.edge.value", intervals = single_int, iterator = single_int)
    expect_equal(track_single$test.edge.track, value_single$test.edge.value, tolerance = 1e-6)

    # Clean up
    gvtrack.rm("test.edge.track")
    gvtrack.rm("test.edge.value")
})

test_that("value-based vtrack in complex expressions matches track-based vtrack", {
    # Extract data from sparse track
    intervals_with_data <- gscreen("!is.na(test.sparse)", gintervals.all())
    extracted <- gextract("test.sparse", intervals = intervals_with_data, iterator = intervals_with_data)

    # Create multiple vtracks with different functions
    gvtrack.create("test.expr.track.avg", "test.sparse", func = "avg")
    gvtrack.create("test.expr.track.min", "test.sparse", func = "min")
    gvtrack.create("test.expr.track.max", "test.sparse", func = "max")

    gvtrack.create("test.expr.value.avg", src = extracted, func = "avg")
    gvtrack.create("test.expr.value.min", src = extracted, func = "min")
    gvtrack.create("test.expr.value.max", src = extracted, func = "max")

    test_int <- gintervals("chr1", 3000, 8000)

    # Test complex expressions
    expressions <- c(
        "test.expr.track.avg * 2",
        "test.expr.track.min + test.expr.track.max",
        "(test.expr.track.max - test.expr.track.min) / 2",
        "test.expr.track.avg + test.expr.track.min * test.expr.track.max"
    )

    for (expr in expressions) {
        # Replace track names with value names
        value_expr <- gsub("track", "value", expr)

        # Extract with both expressions
        track_result <- gextract(expr, intervals = test_int, iterator = test_int, colnames = "value")
        value_result <- gextract(value_expr, intervals = test_int, iterator = test_int, colnames = "value")

        # Compare
        expect_equal(track_result$value, value_result$value,
            tolerance = 1e-6, info = paste("Expression:", expr)
        )
    }

    # Clean up
    gvtrack.rm("test.expr.track.avg")
    gvtrack.rm("test.expr.track.min")
    gvtrack.rm("test.expr.track.max")
    gvtrack.rm("test.expr.value.avg")
    gvtrack.rm("test.expr.value.min")
    gvtrack.rm("test.expr.value.max")
})
