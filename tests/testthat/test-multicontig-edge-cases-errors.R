load_test_db()
skip_if(!getOption("gmulticontig.indexed_format", FALSE))

# Edge cases and error handling tests for multi-contig implementation
# Tests for boundary conditions, error recovery, and robustness

# ============================================================================
# Track index file integrity tests
# ============================================================================

test_that("gtrack.info detects corrupted index file", {
    withr::defer(gtrack.rm("temp.corrupt_idx", force = TRUE))
    gtrack.rm("temp.corrupt_idx", force = TRUE)
    gtrack.create("temp.corrupt_idx", "", "test.fixedbin")
    gtrack.convert_to_indexed("temp.corrupt_idx")

    # Get track directory
    track_path <- file.path(.misha$GROOT, "tracks", "temp/corrupt_idx.track")
    idx_path <- file.path(track_path, "track.idx")
    # Corrupt the index by truncating it
    con <- file(idx_path, "r+b")
    seek(con, 0)
    writeBin(raw(10), con) # Write garbage
    close(con)
    # Should detect corruption
    expect_error(
        gtrack.info("temp.corrupt_idx"),
        "checksum|corrupt|invalid"
    )
})

test_that("missing track.dat file is detected", {
    withr::defer(gtrack.rm("temp.missing_dat", force = TRUE))

    gtrack.rm("temp.missing_dat", force = TRUE)
    gtrack.create("temp.missing_dat", "", "test.fixedbin")
    gtrack.convert_to_indexed("temp.missing_dat")
    # Remove track.dat
    track_path <- file.path(.misha$GROOT, "tracks", "temp/missing_dat.track")
    dat_path <- file.path(track_path, "track.dat")
    file.remove(dat_path)

    # Should fail
    expect_error(
        gextract("temp.missing_dat", gintervals(1, 0, 1000))
    )
})

# ============================================================================
# Boundary condition tests
# ============================================================================

test_that("extraction at chromosome boundaries works", {
    withr::defer(gtrack.rm("temp.chr_boundary", force = TRUE))

    gtrack.rm("temp.chr_boundary", force = TRUE)
    gtrack.create("temp.chr_boundary", "", "test.fixedbin")
    gtrack.convert_to_indexed("temp.chr_boundary")

    # Get chromosome sizes
    chr_sizes <- gintervals.all()

    # Extract from very end of chromosome 1
    chr1_size <- chr_sizes$end[chr_sizes$chrom == "chr1"]

    result1 <- gextract(
        "temp.chr_boundary",
        gintervals(1, chr1_size - 1000, chr1_size)
    )
    expect_true(nrow(result1) > 0)

    # Extract from very start of chromosome 2
    result2 <- gextract(
        "temp.chr_boundary",
        gintervals(2, 1, 1000)
    )
    expect_true(nrow(result2) > 0)
})

# ============================================================================
# Large file tests
# ============================================================================

test_that("very large bin_size works correctly", {
    skip_if_not(dir.exists(.misha$GROOT), "Genome root not available")

    withr::defer(gtrack.rm("temp.large_bins", force = TRUE))

    # Create track with very large bins
    gtrack.rm("temp.large_bins", force = TRUE)
    gtrack.create("temp.large_bins",
        description = "test",
        expr = "test.fixedbin",
        iterator = 1000000
    ) # 1MB bins

    gtrack.convert_to_indexed("temp.large_bins")

    info <- gtrack.info("temp.large_bins")
    expect_equal(info$bin.size, 1000000)

    result <- gextract("temp.large_bins", gintervals(1, 0, 10000000))
    expect_true(nrow(result) > 0)
})

test_that("handling tracks with thousands of intervals (sparse)", {
    withr::defer(gtrack.rm("temp.many_intervals", force = TRUE))

    gtrack.rm("temp.many_intervals", force = TRUE)
    # This test creates a sparse track with many intervals
    # First create intervals
    n_intervals <- 5000
    intervs <- gintervals(
        chrom = rep(1:3, length.out = n_intervals),
        start = seq(0, by = 1000, length.out = n_intervals),
        end = seq(100, by = 1000, length.out = n_intervals)
    )

    # Create sparse track on these intervals
    gtrack.create("temp.many_intervals",
        description = "test",
        expr = "runif(1)",
        iterator = intervs
    )

    gtrack.convert_to_indexed("temp.many_intervals")

    # Extract all
    result <- gextract("temp.many_intervals",
        gintervals(c(1, 2, 3)),
        iterator = "temp.many_intervals"
    )

    expect_true(nrow(result) > 1000)
})

# ============================================================================
# Concurrent access simulation
# ============================================================================

test_that("multiple simultaneous reads from same converted track", {
    withr::defer(gtrack.rm("temp.concurrent_read", force = TRUE))

    gtrack.rm("temp.concurrent_read", force = TRUE)
    gtrack.create("temp.concurrent_read", "", "test.fixedbin")
    gtrack.convert_to_indexed("temp.concurrent_read")

    # Simulate concurrent access by rapid sequential reads
    results <- list()
    for (i in 1:20) {
        results[[i]] <- gextract(
            "temp.concurrent_read",
            gintervals(sample(1:3, 1), 0, 10000)
        )
    }

    # All should succeed
    expect_equal(length(results), 20)
    for (r in results) {
        expect_true(nrow(r) > 0)
    }
})

# ============================================================================
# Numeric precision tests
# ============================================================================

test_that("very small values preserved in converted track", {
    withr::defer(gtrack.rm("temp.tiny_vals", force = TRUE))

    gtrack.rm("temp.tiny_vals", force = TRUE)
    gtrack.create("temp.tiny_vals", "", "test.fixedbin")

    # Set very small values
    tiny_val <- 1e-30
    gtrack.modify(
        "temp.tiny_vals",
        as.character(tiny_val),
        gintervals(1, 0, 1000)
    )

    vals_before <- gextract("temp.tiny_vals", gintervals(1, 0, 1000))

    gtrack.convert_to_indexed("temp.tiny_vals")

    vals_after <- gextract("temp.tiny_vals", gintervals(1, 0, 1000))

    # Should preserve small values (within float precision)
    expect_equal(vals_before$test.tiny_vals, vals_after$test.tiny_vals,
        tolerance = 1e-35
    )
})

test_that("negative and positive infinity handled", {
    withr::defer(gtrack.rm("temp.infinity", force = TRUE))

    gtrack.rm("temp.infinity", force = TRUE)
    gtrack.create("temp.infinity", "", "test.fixedbin")

    # Inf should be converted to NaN in misha
    gtrack.modify("temp.infinity", "Inf", gintervals(1, 0, 1000))
    gtrack.modify("temp.infinity", "-Inf", gintervals(1, 1000, 2000))

    vals_before <- gextract("temp.infinity", gintervals(1, 0, 2000))

    gtrack.convert_to_indexed("temp.infinity")

    vals_after <- gextract("temp.infinity", gintervals(1, 0, 2000))

    # Inf values should become NaN
    expect_true(all(is.na(vals_after$test.infinity)))
})

# ============================================================================
# Track attribute edge cases
# ============================================================================

test_that("tracks with no attributes convert correctly", {
    withr::defer(gtrack.rm("temp.no_attrs", force = TRUE))

    gtrack.rm("temp.no_attrs", force = TRUE)
    gtrack.create("temp.no_attrs", "", "test.fixedbin")
    # Remove all attributes if any
    attr_file_path <- file.path(.misha$GROOT, "tracks", "temp/no_attrs.track/.attributes")
    file.remove(attr_file_path)

    gtrack.convert_to_indexed("temp.no_attrs")

    # Should still work
    result <- gextract("temp.no_attrs", gintervals(1, 0, 1000))
    expect_true(nrow(result) > 0)
})

test_that("tracks with many attributes convert correctly", {
    withr::defer(gtrack.rm("temp.many_attrs", force = TRUE))

    gtrack.rm("temp.many_attrs", force = TRUE)
    gtrack.create("temp.many_attrs", "", "test.fixedbin")

    # Set many attributes
    for (i in 1:50) {
        gtrack.attr.set(
            "temp.many_attrs",
            paste0("key_", i),
            paste0("value_", i)
        )
    }

    attrs_before <- gtrack.attr.export("temp.many_attrs")

    gtrack.convert_to_indexed("temp.many_attrs")

    attrs_after <- gtrack.attr.export("temp.many_attrs")

    expect_equal(attrs_before, attrs_after)
})

# ============================================================================
# Interval edge cases
# ============================================================================

test_that("intervals with unusual metadata column names", {
    withr::defer(gintervals.rm("temp.weird_cols", force = TRUE))

    gintervals.rm("temp.weird_cols", force = TRUE)
    intervs <- gintervals(c(1, 2, 3), 0, 1000)
    intervs$`strange-name` <- 1:3
    intervs$`another weird name!` <- c("a", "b", "c")
    intervs$`123numeric` <- c(1.1, 2.2, 3.3)

    gintervals.save("temp.weird_cols", intervs)
    gintervals.convert_to_indexed("temp.weird_cols")

    loaded <- gintervals.load("temp.weird_cols")

    expect_true("strange-name" %in% names(loaded))
    expect_true("another weird name!" %in% names(loaded))
    expect_true("123numeric" %in% names(loaded))
})

test_that("intervals with NA values in metadata", {
    withr::defer(gintervals.rm("temp.na_meta", force = TRUE))

    gintervals.rm("temp.na_meta", force = TRUE)
    intervs <- gintervals(c(1, 2, 3), 0, 1000)
    intervs$score <- c(10, NA, 30)
    intervs$name <- c("a", NA, "c")

    gintervals.save("temp.na_meta", intervs)

    before <- gintervals.load("temp.na_meta")

    gintervals.convert_to_indexed("temp.na_meta")

    after <- gintervals.load("temp.na_meta")

    expect_equal(is.na(before$score), is.na(after$score))
    expect_equal(is.na(before$name), is.na(after$name))
})

# ============================================================================
# Error recovery tests
# ============================================================================

test_that("track usable after partial conversion failure", {
    withr::defer(gtrack.rm("temp.partial_fail", force = TRUE))

    gtrack.rm("temp.partial_fail", force = TRUE)
    gtrack.create("temp.partial_fail", "", "test.fixedbin")

    # Verify track works before
    before <- gextract("temp.partial_fail", gintervals(1, 0, 1000))
    expect_true(nrow(before) > 0)

    # Try to convert (should succeed normally)
    # In real scenario, could fail partway
    gtrack.convert_to_indexed("temp.partial_fail")

    # Track should still be usable
    after <- gextract("temp.partial_fail", gintervals(1, 0, 1000))
    expect_true(nrow(after) > 0)
})

# ============================================================================
# Chromosome name edge cases
# ============================================================================

test_that("tracks work with non-standard chromosome names", {
    # This test depends on genome having various chromosome names
    # Skip if not applicable
    skip_if_not(dir.exists(.misha$GROOT), "Genome root not available")

    gtrack.rm("temp.all_chroms", force = TRUE)
    chr_info <- gintervals.all()

    # Just verify we can work with all chromosomes
    withr::defer(gtrack.rm("temp.all_chroms", force = TRUE))

    gtrack.create("temp.all_chroms", "", "test.fixedbin")
    gtrack.convert_to_indexed("temp.all_chroms")

    # Try extracting from each chromosome
    for (i in 1:min(10, nrow(chr_info))) {
        chr <- chr_info$chrom[i]
        result <- gextract(
            "temp.all_chroms",
            gintervals(chr, 0, 1000)
        )
        expect_true(is.data.frame(result))
    }
})

# ============================================================================
# Memory and resource tests
# ============================================================================

test_that("converted track releases file handles properly", {
    withr::defer(gtrack.rm("temp.file_handles", force = TRUE))

    gtrack.rm("temp.file_handles", force = TRUE)
    gtrack.create("temp.file_handles", "", "test.fixedbin")
    gtrack.convert_to_indexed("temp.file_handles")

    # Perform many extractions
    for (i in 1:100) {
        result <- gextract(
            "temp.file_handles",
            gintervals(1, i * 100, i * 100 + 100)
        )
    }

    # Should not run out of file handles
    # If there's a leak, this would fail
    result <- gextract("temp.file_handles", gintervals(1, 0, 1000), colnames = "test.file_handles")
    expect_true(nrow(result) > 0)
})

# ============================================================================
# Special value handling
# ============================================================================

test_that("tracks with all NaN values work correctly", {
    withr::defer(gtrack.rm("temp.all_nan", force = TRUE))

    gtrack.rm("temp.all_nan", force = TRUE)
    gtrack.create("temp.all_nan", "", "test.fixedbin")

    # Set everything to NaN
    gtrack.modify("temp.all_nan", "NaN", gintervals(c(1, 2, 3)))

    vals_before <- gextract("temp.all_nan", gintervals(c(1, 2, 3)), colnames = "test.all_nan")

    gtrack.convert_to_indexed("temp.all_nan")

    vals_after <- gextract("temp.all_nan", gintervals(c(1, 2, 3)), colnames = "test.all_nan")

    expect_true(all(is.na(vals_before$test.all_nan)))
    expect_true(all(is.na(vals_after$test.all_nan)))
})

test_that("tracks with mixed finite and NaN values", {
    withr::defer(gtrack.rm("temp.mixed_nan", force = TRUE))

    gtrack.rm("temp.mixed_nan", force = TRUE)
    gtrack.create("temp.mixed_nan", "", "test.fixedbin")

    # Set some values to NaN, leave others
    gtrack.modify("temp.mixed_nan", "NaN", gintervals(1, 0, 5000))
    # Rest should have values from test.fixedbin

    vals_before <- gextract("temp.mixed_nan", gintervals(1, 0, 10000), colnames = "test.mixed_nan")

    gtrack.convert_to_indexed("temp.mixed_nan")

    vals_after <- gextract("temp.mixed_nan", gintervals(1, 0, 10000), colnames = "test.mixed_nan")

    expect_equal(
        is.na(vals_before$test.mixed_nan),
        is.na(vals_after$test.mixed_nan)
    )
})

# ============================================================================
# Index cache behavior
# ============================================================================

test_that("index cache works correctly across multiple tracks", {
    withr::defer(gtrack.rm("temp.cache1", force = TRUE))
    withr::defer(gtrack.rm("temp.cache2", force = TRUE))
    withr::defer(gtrack.rm("temp.cache3", force = TRUE))

    # Create and convert multiple tracks
    gtrack.create("temp.cache1", "", "test.fixedbin")
    gtrack.create("temp.cache2", "", "test.fixedbin")
    gtrack.create("temp.cache3", "", "test.fixedbin")

    gtrack.convert_to_indexed("temp.cache1")
    gtrack.convert_to_indexed("temp.cache2")
    gtrack.convert_to_indexed("temp.cache3")

    # Access all of them multiple times
    for (i in 1:5) {
        r1 <- gextract("temp.cache1", gintervals(1, 0, 1000))
        r2 <- gextract("temp.cache2", gintervals(1, 0, 1000))
        r3 <- gextract("temp.cache3", gintervals(1, 0, 1000))

        expect_true(nrow(r1) > 0)
        expect_true(nrow(r2) > 0)
        expect_true(nrow(r3) > 0)
    }
})
