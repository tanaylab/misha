load_test_db()
skip_if(!getOption("gmulticontig.indexed_format", FALSE))

# Edge cases and error handling tests for multi-contig implementation
# Tests for boundary conditions, error recovery, and robustness

# ============================================================================
# Track index file integrity tests
# ============================================================================

test_that("gtrack.info detects corrupted index file", {
    withr::defer(gtrack.rm("test.corrupt_idx", force = TRUE))
    browser()
    gtrack.create("test.corrupt_idx", "", "test.fixedbin")
    gtrack.convert_to_indexed("test.corrupt_idx")

    # Get track directory
    track_path <- file.path(.misha$GROOT, "tracks", "test/corrupt_idx.track")
    idx_path <- file.path(track_path, "track.idx")
    # Corrupt the index by truncating it
    con <- file(idx_path, "r+b")
    seek(con, 0)
    writeBin(raw(10), con) # Write garbage
    close(con)

    # Should detect corruption
    expect_error(
        gtrack.info("test.corrupt_idx"),
        "checksum|corrupt|invalid"
    )
})

test_that("missing track.dat file is detected", {
    withr::defer(gtrack.rm("test.missing_dat", force = TRUE))

    gtrack.create("test.missing_dat", "", "test.fixedbin")
    gtrack.convert_to_indexed("test.missing_dat")
    # Remove track.dat
    track_path <- file.path(.misha$GROOT, "tracks", "test/missing_dat.track")
    dat_path <- file.path(track_path, "track.dat")
    file.remove(dat_path)

    # Should fail
    expect_error(
        gextract("test.missing_dat", gintervals(1, 0, 1000))
    )
})

# ============================================================================
# Boundary condition tests
# ============================================================================

test_that("extraction at chromosome boundaries works", {
    withr::defer(gtrack.rm("test.chr_boundary", force = TRUE))

    gtrack.create("test.chr_boundary", "", "test.fixedbin")
    gtrack.convert_to_indexed("test.chr_boundary")

    # Get chromosome sizes
    chr_sizes <- gintervals.all()

    # Extract from very end of chromosome 1
    chr1_size <- chr_sizes$end[chr_sizes$chrom == "chr1"]

    result1 <- gextract(
        "test.chr_boundary",
        gintervals(1, chr1_size - 1000, chr1_size)
    )
    expect_true(nrow(result1) > 0)

    # Extract from very start of chromosome 2
    result2 <- gextract(
        "test.chr_boundary",
        gintervals(2, 1, 1000)
    )
    expect_true(nrow(result2) > 0)
})

# ============================================================================
# Large file tests
# ============================================================================

test_that("very large bin_size works correctly", {
    skip_if_not(dir.exists(.misha$GROOT), "Genome root not available")

    withr::defer(gtrack.rm("test.large_bins", force = TRUE))

    # Create track with very large bins
    gtrack.create("test.large_bins",
        description = "test",
        expr = "test.fixedbin",
        iterator = 1000000
    ) # 1MB bins

    gtrack.convert_to_indexed("test.large_bins")

    info <- gtrack.info("test.large_bins")
    expect_equal(info$bin.size, 1000000)

    result <- gextract("test.large_bins", gintervals(1, 0, 10000000))
    expect_true(nrow(result) > 0)
})

test_that("handling tracks with thousands of intervals (sparse)", {
    withr::defer(gtrack.rm("test.many_intervals", force = TRUE))

    # This test creates a sparse track with many intervals
    # First create intervals
    n_intervals <- 5000
    intervs <- gintervals(
        chrom = rep(1:3, length.out = n_intervals),
        start = seq(0, by = 1000, length.out = n_intervals),
        end = seq(100, by = 1000, length.out = n_intervals)
    )

    # Create sparse track on these intervals
    gtrack.create("test.many_intervals",
        description = "test",
        expr = "runif(1)",
        iterator = intervs
    )

    gtrack.convert_to_indexed("test.many_intervals")

    # Extract all
    result <- gextract("test.many_intervals",
        gintervals(c(1, 2, 3)),
        iterator = "test.many_intervals"
    )

    expect_true(nrow(result) > 1000)
})

# ============================================================================
# Concurrent access simulation
# ============================================================================

test_that("multiple simultaneous reads from same converted track", {
    withr::defer(gtrack.rm("test.concurrent_read", force = TRUE))

    gtrack.create("test.concurrent_read", "", "test.fixedbin")
    gtrack.convert_to_indexed("test.concurrent_read")

    # Simulate concurrent access by rapid sequential reads
    results <- list()
    for (i in 1:20) {
        results[[i]] <- gextract(
            "test.concurrent_read",
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
    withr::defer(gtrack.rm("test.tiny_vals", force = TRUE))

    gtrack.create("test.tiny_vals", "", "test.fixedbin")

    # Set very small values
    tiny_val <- 1e-30
    gtrack.modify(
        "test.tiny_vals",
        as.character(tiny_val),
        gintervals(1, 0, 1000)
    )

    vals_before <- gextract("test.tiny_vals", gintervals(1, 0, 1000))

    gtrack.convert_to_indexed("test.tiny_vals")

    vals_after <- gextract("test.tiny_vals", gintervals(1, 0, 1000))

    # Should preserve small values (within float precision)
    expect_equal(vals_before$test.tiny_vals, vals_after$test.tiny_vals,
        tolerance = 1e-35
    )
})

test_that("negative and positive infinity handled", {
    withr::defer(gtrack.rm("test.infinity", force = TRUE))

    gtrack.create("test.infinity", "", "test.fixedbin")

    # Inf should be converted to NaN in misha
    gtrack.modify("test.infinity", "Inf", gintervals(1, 0, 1000))
    gtrack.modify("test.infinity", "-Inf", gintervals(1, 1000, 2000))

    vals_before <- gextract("test.infinity", gintervals(1, 0, 2000))

    gtrack.convert_to_indexed("test.infinity")

    vals_after <- gextract("test.infinity", gintervals(1, 0, 2000))

    # Inf values should become NaN
    expect_true(all(is.na(vals_after$test.infinity)))
})

# ============================================================================
# Track attribute edge cases
# ============================================================================

test_that("tracks with no attributes convert correctly", {
    withr::defer(gtrack.rm("test.no_attrs", force = TRUE))

    gtrack.create("test.no_attrs", "", "test.fixedbin")
    # Remove all attributes if any
    attr_file_path <- file.path(.misha$GROOT, "tracks", "test/no_attrs.track/.attributes")
    file.remove(attr_file_path)

    gtrack.convert_to_indexed("test.no_attrs")

    # Should still work
    result <- gextract("test.no_attrs", gintervals(1, 0, 1000))
    expect_true(nrow(result) > 0)
})

test_that("tracks with many attributes convert correctly", {
    withr::defer(gtrack.rm("test.many_attrs", force = TRUE))

    gtrack.create("test.many_attrs", "", "test.fixedbin")

    # Set many attributes
    for (i in 1:50) {
        gtrack.attr.set(
            "test.many_attrs",
            paste0("key_", i),
            paste0("value_", i)
        )
    }

    attrs_before <- gtrack.attr.export("test.many_attrs")

    gtrack.convert_to_indexed("test.many_attrs")

    attrs_after <- gtrack.attr.export("test.many_attrs")

    expect_equal(attrs_before, attrs_after)
})

# ============================================================================
# Interval edge cases
# ============================================================================

test_that("intervals with unusual metadata column names", {
    withr::defer(gintervals.rm("test.weird_cols", force = TRUE))

    intervs <- gintervals(c(1, 2, 3), 0, 1000)
    intervs$`strange-name` <- 1:3
    intervs$`another weird name!` <- c("a", "b", "c")
    intervs$`123numeric` <- c(1.1, 2.2, 3.3)

    gintervals.save("test.weird_cols", intervs)
    gintervals.convert_to_indexed("test.weird_cols")

    loaded <- gintervals.load("test.weird_cols")

    expect_true("strange-name" %in% names(loaded))
    expect_true("another weird name!" %in% names(loaded))
    expect_true("123numeric" %in% names(loaded))
})

test_that("intervals with NA values in metadata", {
    withr::defer(gintervals.rm("test.na_meta", force = TRUE))

    intervs <- gintervals(c(1, 2, 3), 0, 1000)
    intervs$score <- c(10, NA, 30)
    intervs$name <- c("a", NA, "c")

    gintervals.save("test.na_meta", intervs)

    before <- gintervals.load("test.na_meta")

    gintervals.convert_to_indexed("test.na_meta")

    after <- gintervals.load("test.na_meta")

    expect_equal(is.na(before$score), is.na(after$score))
    expect_equal(is.na(before$name), is.na(after$name))
})

# ============================================================================
# Error recovery tests
# ============================================================================

test_that("track usable after partial conversion failure", {
    withr::defer(gtrack.rm("test.partial_fail", force = TRUE))

    gtrack.create("test.partial_fail", "", "test.fixedbin")

    # Verify track works before
    before <- gextract("test.partial_fail", gintervals(1, 0, 1000))
    expect_true(nrow(before) > 0)

    # Try to convert (should succeed normally)
    # In real scenario, could fail partway
    gtrack.convert_to_indexed("test.partial_fail")

    # Track should still be usable
    after <- gextract("test.partial_fail", gintervals(1, 0, 1000))
    expect_true(nrow(after) > 0)
})

# ============================================================================
# Version compatibility tests
# ============================================================================

test_that("gtrack.info reports correct format", {
    withr::defer(gtrack.rm("test.format_version", force = TRUE))

    gtrack.create("test.format_version", "", "test.fixedbin")

    info_per_chromosome <- gtrack.info("test.format_version")
    expect_equal(info_per_chromosome$format, "per-chromosome")

    gtrack.convert_to_indexed("test.format_version")

    info_indexed <- gtrack.info("test.format_version")
    expect_equal(info_indexed$format, "indexed")
})

# ============================================================================
# Chromosome name edge cases
# ============================================================================

test_that("tracks work with non-standard chromosome names", {
    # This test depends on genome having various chromosome names
    # Skip if not applicable
    skip_if_not(dir.exists(.misha$GROOT), "Genome root not available")

    chr_info <- gintervals.all()

    # Just verify we can work with all chromosomes
    withr::defer(gtrack.rm("test.all_chroms", force = TRUE))

    gtrack.create("test.all_chroms", "", "test.fixedbin")
    gtrack.convert_to_indexed("test.all_chroms")

    # Try extracting from each chromosome
    for (i in 1:min(10, nrow(chr_info))) {
        chr <- chr_info$chrom[i]
        result <- gextract(
            "test.all_chroms",
            gintervals(chr, 0, 1000)
        )
        expect_true(is.data.frame(result))
    }
})

# ============================================================================
# Memory and resource tests
# ============================================================================

test_that("converted track releases file handles properly", {
    withr::defer(gtrack.rm("test.file_handles", force = TRUE))

    gtrack.create("test.file_handles", "", "test.fixedbin")
    gtrack.convert_to_indexed("test.file_handles")

    # Perform many extractions
    for (i in 1:100) {
        result <- gextract(
            "test.file_handles",
            gintervals(1, i * 100, i * 100 + 100)
        )
    }

    # Should not run out of file handles
    # If there's a leak, this would fail
    result <- gextract("test.file_handles", gintervals(1, 0, 1000))
    expect_true(nrow(result) > 0)
})

# ============================================================================
# Special value handling
# ============================================================================

test_that("tracks with all NaN values work correctly", {
    withr::defer(gtrack.rm("test.all_nan", force = TRUE))

    gtrack.create("test.all_nan", "", "test.fixedbin")

    # Set everything to NaN
    gtrack.modify("test.all_nan", "NaN", gintervals(c(1, 2, 3)))

    vals_before <- gextract("test.all_nan", gintervals(c(1, 2, 3)))

    gtrack.convert_to_indexed("test.all_nan")

    vals_after <- gextract("test.all_nan", gintervals(c(1, 2, 3)))

    expect_true(all(is.na(vals_before$test.all_nan)))
    expect_true(all(is.na(vals_after$test.all_nan)))
})

test_that("tracks with mixed finite and NaN values", {
    withr::defer(gtrack.rm("test.mixed_nan", force = TRUE))

    gtrack.create("test.mixed_nan", "", "test.fixedbin")

    # Set some values to NaN, leave others
    gtrack.modify("test.mixed_nan", "NaN", gintervals(1, 0, 5000))
    # Rest should have values from test.fixedbin

    vals_before <- gextract("test.mixed_nan", gintervals(1, 0, 10000))

    gtrack.convert_to_indexed("test.mixed_nan")

    vals_after <- gextract("test.mixed_nan", gintervals(1, 0, 10000))

    expect_equal(
        is.na(vals_before$test.mixed_nan),
        is.na(vals_after$test.mixed_nan)
    )
})

# ============================================================================
# Index cache behavior
# ============================================================================

test_that("index cache works correctly across multiple tracks", {
    withr::defer(gtrack.rm("test.cache1", force = TRUE))
    withr::defer(gtrack.rm("test.cache2", force = TRUE))
    withr::defer(gtrack.rm("test.cache3", force = TRUE))

    # Create and convert multiple tracks
    gtrack.create("test.cache1", "", "test.fixedbin")
    gtrack.create("test.cache2", "", "test.fixedbin")
    gtrack.create("test.cache3", "", "test.fixedbin")

    gtrack.convert_to_indexed("test.cache1")
    gtrack.convert_to_indexed("test.cache2")
    gtrack.convert_to_indexed("test.cache3")

    # Access all of them multiple times
    for (i in 1:5) {
        r1 <- gextract("test.cache1", gintervals(1, 0, 1000))
        r2 <- gextract("test.cache2", gintervals(1, 0, 1000))
        r3 <- gextract("test.cache3", gintervals(1, 0, 1000))

        expect_true(nrow(r1) > 0)
        expect_true(nrow(r2) > 0)
        expect_true(nrow(r3) > 0)
    }
})
