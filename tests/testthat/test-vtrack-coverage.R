create_isolated_test_db()

test_that("Coverage virtual track works correctly", {
    # Basic coverage test - single overlap
    gvtrack.create("cov1", gintervals(1, 50, 100), "coverage")
    res <- gextract("cov1", gintervals(1, 0, 100), iterator = 50)
    expect_equal(res$cov1, c(0, 1)) # First 50bp: no overlap, Second 50bp: full overlap

    # Multiple source intervals
    gvtrack.create("cov2", rbind(
        gintervals(1, 150, 200), # 50bp interval
        gintervals(1, 320, 340) # 20bp interval
    ), "coverage")
    res <- gextract("cov2", rbind(
        gintervals(1, 0, 200), # Should have 50/200 coverage from first interval
        gintervals(1, 250, 500) # Should have 20/250 coverage from second interval
    ), iterator = 100)
    expect_equal(res$cov2, c(0, 0.5, 0, 0.2, 0))

    # Overlapping source intervals
    gvtrack.create("cov3", rbind(
        gintervals(1, 100, 200),
        gintervals(1, 150, 250) # Overlaps previous interval by 50bp
    ), "coverage")
    res <- gextract("cov3", gintervals(1, 0, 300), iterator = 100)
    expect_equal(res$cov3, c(0, 1, 0.5)) # Second bin fully covered despite overlap

    gvtrack.create("cov4", gintervals(1, 100, 250), "coverage")
    res <- gextract("cov4", gintervals(1, 0, 300), iterator = 100)
    expect_equal(res$cov4, c(0, 1, 0.5))
})

test_that("Virtual tracks work with test.sparse source", {
    src_intervs <- gextract("test.sparse", gintervals.all())

    gvtrack.create("test_cov", src_intervs, "coverage")
    res <- gextract("test_cov", gintervals(1, 0, 500), iterator = 100)
    expect_equal(nrow(res), 5)
})

test_that("Virtual tracks work with test.fixedbin source", {
    src_intervs <- gextract("test.fixedbin", gintervals.all())

    gvtrack.create("test_cov", src_intervs, "coverage")
    res <- gextract("test_cov", gintervals(1, 0, 500), iterator = 100)
    expect_equal(nrow(res), 5)
})

test_that("Coverage virtual track calculates basic coverage correctly", {
    # Test basic non-overlapping coverage
    source_intervals <- rbind(
        gintervals(1, 100, 200), # 100bp interval
        gintervals(1, 300, 400) # 100bp interval
    )
    gvtrack.create("basic_cov", source_intervals, "coverage")

    # Test single interval fully within source
    res <- gextract("basic_cov", gintervals(1, 120, 180), iterator = 60)
    expect_equal(res$basic_cov, 1.0) # Should be fully covered

    # Test interval partially overlapping source
    res <- gextract("basic_cov", gintervals.all(), iterator = gintervals(1, 250, 350))
    expect_equal(res$basic_cov, 0.5) # Half should be covered (50/100)

    # Test interval containing multiple source intervals
    res <- gextract("basic_cov", gintervals(1, 0, 500), iterator = 100)
    expected <- c(0, 1, 0, 1, 0) # Coverage per 100bp bin
    expect_equal(res$basic_cov, expected)
})

test_that("Coverage virtual track handles overlapping source intervals correctly", {
    # Create overlapping source intervals
    source_intervals <- rbind(
        gintervals(1, 100, 300), # 200bp interval
        gintervals(1, 200, 400) # 200bp interval overlapping by 100bp
    )
    gvtrack.create("overlap_cov", source_intervals, "coverage")

    # Test interval in region with overlap
    res <- gextract("overlap_cov", gintervals(1, 200, 300), iterator = 100)
    expect_equal(res$overlap_cov, 1.0) # Should be fully covered despite overlap

    # Test across entire region
    res <- gextract("overlap_cov", gintervals(1, 0, 500), iterator = 100)
    expected <- c(0, 1, 1, 1, 0) # Coverage for 100bp bins
    expect_equal(res$overlap_cov, expected)
})

test_that("Coverage virtual track handles edge cases correctly", {
    # Test single base intervals
    source_intervals <- rbind(
        gintervals(1, 100, 101),
        gintervals(1, 200, 201)
    )
    gvtrack.create("single_base_cov", source_intervals, "coverage")
    res <- gextract("single_base_cov", gintervals(1, 0, 300), iterator = 100)
    expected <- c(0, 0.01, 0.01) # One base out of 100 covered in bins 2 and 3
    expect_equal(res$single_base_cov, expected)

    # Test with multiple chromosomes
    source_intervals <- rbind(
        gintervals(1, 100, 200),
        gintervals(2, 100, 200)
    )
    gvtrack.create("multi_chrom_cov", source_intervals, "coverage")
    res <- gextract("multi_chrom_cov",
        rbind(
            gintervals(1, 0, 300),
            gintervals(2, 0, 300)
        ),
        iterator = 100
    )
    expected <- c(0, 1, 0, 0, 1, 0) # Coverage pattern repeats for each chromosome
    expect_equal(res$multi_chrom_cov, expected)
})

test_that("Coverage virtual track works with iterator modifiers", {
    source_intervals <- gintervals(1, 100, 200)
    gvtrack.create("mod_cov", source_intervals, "coverage")

    # Shift iterator bounds
    gvtrack.iterator("mod_cov", sshift = -50, eshift = 50)

    # Test that coverage is calculated on modified intervals
    res <- gextract("mod_cov", gintervals(1, 0, 300), iterator = 100)
    expected <- c(1 / 3, 0.5, 0.25)
    expect_equal(res$mod_cov, expected)
})

test_that("Coverage virtual track handles very large intervals efficiently", {
    # Create source intervals totaling 1Mb
    source_intervals <- do.call(rbind, lapply(1:100, function(i) {
        gintervals(1, (i - 1) * 10000, i * 10000)
    }))

    gvtrack.create("large_cov", source_intervals, "coverage")

    # Test with 100kb iterator intervals
    res <- gextract("large_cov", gintervals(1, 0, 1000000), iterator = 100000)
    expect_equal(length(res$large_cov), 10)
    expect_true(all(res$large_cov == 1.0)) # Should be fully covered
})

test_that("Coverage virtual track handles iterator intervals outside source range", {
    # Create source intervals in a specific range
    source_intervals <- gintervals(1, 1000, 2000)
    gvtrack.create("out_of_range_cov", source_intervals, "coverage")

    # Test completely outside range
    res <- gextract("out_of_range_cov", gintervals(1, 3000, 4000), iterator = 100)
    expect_equal(res$out_of_range_cov, rep(0, 10)) # Should be all zeros

    # Test partially outside range
    res <- gextract("out_of_range_cov", gintervals(1, 1500, 2500), iterator = 100)
    expected <- c(rep(1, 5), rep(0, 5)) # First half covered, second half not
    expect_equal(res$out_of_range_cov, expected)
})

test_that("Coverage virtual track handles many small, dispersed intervals", {
    # Create many small, dispersed intervals
    n_intervals <- 1000
    start_positions <- seq(1, 100000, length.out = n_intervals)

    # Each interval is 10bp long
    source_intervals <- do.call(rbind, lapply(start_positions, function(start) {
        gintervals(1, start, start + 10)
    }))

    gvtrack.create("dispersed_cov", source_intervals, "coverage")

    # Test with large iterator intervals
    res <- gextract("dispersed_cov", gintervals(1, 0, 100000), iterator = 10000)

    # Expect approximately 10% coverage in each bin (10bp * 100 intervals per bin / 10000bp)
    for (cov in res$dispersed_cov) {
        expect_true(abs(cov - 0.1) < 0.02) # Allow small deviation due to interval distribution
    }
})

test_that("Coverage virtual track calculates correct coverage with complex overlapping patterns", {
    # Create a complex pattern of overlapping intervals
    source_intervals <- rbind(
        gintervals(1, 100, 300), # 200bp base interval
        gintervals(1, 150, 250), # 100bp fully within first
        gintervals(1, 200, 350), # 150bp overlapping end of first
        gintervals(1, 275, 400) # 125bp overlapping third
    )

    gvtrack.create("complex_cov", source_intervals, "coverage")

    # Test across the entire region
    res <- gextract("complex_cov", gintervals(1, 0, 500), iterator = 100)

    # Expected coverage per 100bp:
    # 0-100: 0%
    # 100-200: 100% (covered by first interval, plus some additional overlap)
    # 200-300: 100% (covered by multiple overlapping intervals)
    # 300-400: 100% (covered by third and fourth intervals)
    # 400-500: 0%
    expected <- c(0, 1, 1, 1, 0)
    expect_equal(res$complex_cov, expected)
})

test_that("Coverage virtual track works with other virtual track functions", {
    # Create source intervals
    source_intervals <- gintervals(1, 100, 200)

    # Create two different virtual tracks on the same source
    gvtrack.create("distance_vtrack", source_intervals, "distance")
    gvtrack.create("coverage_vtrack", source_intervals, "coverage")

    # Extract both and compare
    res <- gextract(c("distance_vtrack", "coverage_vtrack"),
        gintervals(1, 0, 300),
        iterator = 100
    )

    # Coverage should be 0, 1, 0
    expect_equal(res$coverage_vtrack, c(0, 1, 0))

    # Distance should be positive for first bin, 0 for middle bin (as coordinate is inside interval),
    # and positive for last bin
    expect_true(res$distance_vtrack[1] > 0)
    expect_equal(res$distance_vtrack[2], 0)
    expect_true(res$distance_vtrack[3] > 0)
})

test_that("Coverage virtual track interacts correctly with gdist function", {
    # Create source intervals
    source_intervals <- rbind(
        gintervals(1, 100, 200),
        gintervals(1, 300, 400)
    )

    gvtrack.create("dist_cov", source_intervals, "coverage")

    # Use gdist to bin coverage values
    coverage_dist <- gdist("dist_cov",
        breaks = seq(0, 1, by = 0.2),
        intervals = gintervals(1, 0, 500),
        iterator = 50
    )

    # We expect:
    # 0.0-0.2: bins with little or no coverage
    # 0.2-0.4: bins with partial coverage
    # 0.4-0.6: bins with partial coverage
    # 0.6-0.8: bins with partial coverage
    # 0.8-1.0: bins with high or full coverage
    expect_true(all(coverage_dist >= 0))
    expect_equal(length(coverage_dist), 5) # 5 bins
})

test_that("Coverage virtual track performance is acceptable for large datasets", {
    skip_on_cran() # Skip on CRAN to avoid long test times

    # Create large test data: 10,000 intervals
    n_intervals <- 10000
    starts <- sample(1:1000000, n_intervals)
    ends <- starts + sample(100:1000, n_intervals, replace = TRUE)

    source_intervals <- do.call(rbind, lapply(1:n_intervals, function(i) {
        gintervals(1, starts[i], ends[i])
    }))

    # Time the creation of the virtual track
    start_time <- Sys.time()
    gvtrack.create("perf_cov", source_intervals, "coverage")
    creation_time <- difftime(Sys.time(), start_time, units = "secs")

    # Time a simple extraction
    start_time <- Sys.time()
    res <- gextract("perf_cov", gintervals(1, 0, 1000000), iterator = 100000)
    extraction_time <- difftime(Sys.time(), start_time, units = "secs")

    # Performance assertions - these thresholds should be adjusted based on real-world expectations
    expect_lt(as.numeric(creation_time), 10) # Creation should take less than 10 seconds
    expect_lt(as.numeric(extraction_time), 30) # Extraction should take less than 30 seconds
})

test_that("Coverage virtual track works with genomic scope", {
    # Create source intervals on two chromosomes
    source_intervals <- rbind(
        gintervals(1, 100, 200),
        gintervals(2, 100, 200)
    )

    gvtrack.create("scope_cov", source_intervals, "coverage")

    # Test with chromosome 1 only scope
    res_chr1 <- gextract("scope_cov", gintervals(1, 0, 300), iterator = 100)
    expect_equal(res_chr1$scope_cov, c(0, 1, 0))

    # Test with specific regions scope
    custom_scope <- rbind(
        gintervals(1, 150, 250),
        gintervals(2, 50, 150)
    )

    res_custom <- gextract("scope_cov", custom_scope, iterator = 50) %>% arrange(intervalID)
    expected <- c(1, 0, 0, 1)
    expect_equal(res_custom$scope_cov, expected)
})

test_that("Coverage virtual track works correctly with numeric iterators", {
    # Create source intervals
    source_intervals <- rbind(
        gintervals(1, 150, 250), # 100bp interval
        gintervals(1, 350, 450) # 100bp interval
    )

    gvtrack.create("numeric_cov", source_intervals, "coverage")

    # Test with a 100bp numeric iterator
    # Note: With a numeric iterator of 100, bins would be at positions:
    # [0-100, 100-200, 200-300, 300-400, 400-500, ...]
    res <- gextract("numeric_cov", gintervals(1, 0, 500), iterator = 100)

    # Expected coverage:
    # Bin 0-100: 0% coverage
    # Bin 100-200: 50% coverage (150-200 = 50bp covered)
    # Bin 200-300: 50% coverage (200-250 = 50bp covered)
    # Bin 300-400: 50% coverage (350-400 = 50bp covered)
    # Bin 400-500: 50% coverage (400-450 = 50bp covered)
    expected <- c(0, 0.5, 0.5, 0.5, 0.5)
    expect_equal(res$numeric_cov, expected)

    # Test with a different chromosome
    source_intervals_chr2 <- gintervals(2, 75, 225) # 150bp interval
    gvtrack.create("numeric_cov_chr2", source_intervals_chr2, "coverage")

    # Same 100bp numeric iterator but on chr2
    # Bins would be at positions: [0-100, 100-200, 200-300, ...]
    res_chr2 <- gextract("numeric_cov_chr2", gintervals(2, 0, 300), iterator = 100)

    # Expected coverage:
    # Bin 0-100: 25% coverage (75-100 = 25bp covered)
    # Bin 100-200: 100% coverage (100-200 = 100bp covered)
    # Bin 200-300: 25% coverage (200-225 = 25bp covered)
    expected_chr2 <- c(0.25, 1, 0.25)
    expect_equal(res_chr2$numeric_cov_chr2, expected_chr2)
})

test_that("Coverage virtual track handles non-aligned numeric iterators correctly", {
    # Create source interval at a specific position
    source_intervals <- gintervals(1, 125, 175) # 50bp interval
    gvtrack.create("nonaligned_cov", source_intervals, "coverage")

    # Test with a 100bp iterator
    # Bins would be at positions: [0-100, 100-200, ...]
    res <- gextract("nonaligned_cov", gintervals(1, 0, 200), iterator = 100)

    # Expected coverage:
    # Bin 0-100: 0% coverage
    # Bin 100-200: 50% coverage (125-175 = 50bp covered out of 100bp bin)
    expected <- c(0, 0.5)
    expect_equal(res$nonaligned_cov, expected)

    # Test with a different bin size (50bp)
    # Bins would be at positions: [0-50, 50-100, 100-150, 150-200, ...]
    res_50bp <- gextract("nonaligned_cov", gintervals(1, 0, 200), iterator = 50)

    # Expected coverage:
    # Bin 0-50: 0% coverage
    # Bin 50-100: 0% coverage
    # Bin 100-150: 50% coverage (125-150 = 25bp covered out of 50bp bin)
    # Bin 150-200: 50% coverage (150-175 = 25bp covered out of 50bp bin)
    expected_50bp <- c(0, 0, 0.5, 0.5)
    expect_equal(res_50bp$nonaligned_cov, expected_50bp)
})

test_that("Coverage virtual track works with non-divisible numeric iterators", {
    # Create source interval
    source_intervals <- gintervals(1, 100, 200) # 100bp interval
    gvtrack.create("nondiv_cov", source_intervals, "coverage")

    # Test with a 30bp iterator (not divisible into genome cleanly)
    # Bins would be at positions: [0-30, 30-60, 60-90, 90-120, 120-150, 150-180, 180-210, ...]
    res <- gextract("nondiv_cov", gintervals(1, 0, 210), iterator = 30)

    # Expected coverage:
    # Bin 0-30: 0% coverage
    # Bin 30-60: 0% coverage
    # Bin 60-90: 0% coverage
    # Bin 90-120: 33.3% coverage (100-120 = 20bp covered out of 30bp bin)
    # Bin 120-150: 100% coverage
    # Bin 150-180: 100% coverage
    # Bin 180-210: 66.7% coverage (180-200 = 20bp covered out of 30bp bin)
    expected <- c(0, 0, 0, 2 / 3, 1, 1, 2 / 3)
    # Using expect_equal with tolerance for floating point precision
    for (i in 1:length(expected)) {
        expect_equal(res$nondiv_cov[i], expected[i], tolerance = 0.01)
    }
})

test_that("Coverage virtual track handles crossing chromosome boundaries with numeric iterator", {
    # Create source intervals on two chromosomes
    source_intervals <- rbind(
        gintervals(1, 80, 120), # 40bp interval on chr1
        gintervals(2, 20, 60) # 40bp interval on chr2
    )
    gvtrack.create("cross_chr_cov", source_intervals, "coverage")

    # Test with 50bp iterator across two chromosomes
    # Bins would reset at chromosome boundaries
    # Chr1: [0-50, 50-100, 100-150, ...]
    # Chr2: [0-50, 50-100, ...]
    res <- gextract("cross_chr_cov",
        rbind(
            gintervals(1, 0, 150),
            gintervals(2, 0, 100)
        ),
        iterator = 50
    )

    # Expected coverage:
    # Chr1, Bin 0-50: 0% coverage
    # Chr1, Bin 50-100: 40% coverage (80-100 = 20bp covered out of 50bp bin)
    # Chr1, Bin 100-150: 40% coverage (100-120 = 20bp covered out of 50bp bin)
    # Chr2, Bin 0-50: 60% coverage (20-50 = 30bp covered out of 50bp bin)
    # Chr2, Bin 50-100: 20% coverage (50-60 = 10bp covered out of 50bp bin)
    expected <- c(0, 0.4, 0.4, 0.6, 0.2)
    expect_equal(res$cross_chr_cov, expected)
})

test_that("Coverage virtual track works with large numeric iterator", {
    # Create source intervals
    source_intervals <- rbind(
        gintervals(1, 10000, 20000), # 10,000bp interval
        gintervals(1, 50000, 60000) # 10,000bp interval
    )
    gvtrack.create("large_iter_cov", source_intervals, "coverage")

    # Test with a 50,000bp iterator
    # Bins would be at positions: [0-50000, 50000-100000, ...]
    res <- gextract("large_iter_cov", gintervals(1, 0, 100000), iterator = 50000)

    # Expected coverage:
    # Bin 0-50000: 20% coverage (10000-20000 = 10,000bp covered out of 50,000bp bin)
    # Bin 50000-100000: 20% coverage (50000-60000 = 10,000bp covered out of 50,000bp bin)
    expected <- c(0.2, 0.2)
    expect_equal(res$large_iter_cov, expected)
})

test_that("Coverage virtual track with numeric iterator respects interval boundaries", {
    # Create source interval
    source_intervals <- gintervals(1, 100, 200) # 100bp interval
    gvtrack.create("boundary_cov", source_intervals, "coverage")

    # Test with a query interval that doesn't align with iterator boundaries
    # Numeric iterator of 100bp would normally create bins at [0-100, 100-200, ...]
    # But our query is from 50-350, which doesn't align
    res <- gextract("boundary_cov", gintervals(1, 50, 350), iterator = 100)

    # The iterator should still calculate bins from the beginning of the chromosome:
    # [0-100, 100-200, 200-300, 300-400, ...]
    # But our query only includes portions of these bins: [50-100, 100-200, 200-300, 300-350]
    # Expected coverage for these bins:
    # Bin 50-100 (portion of 0-100): 0% coverage (no overlap with source)
    # Bin 100-200: 100% coverage (full overlap with source)
    # Bin 200-300: 0% coverage (no overlap with source)
    # Bin 300-350 (portion of 300-400): 0% coverage (no overlap with source)
    expected <- c(0, 1, 0, 0)
    expect_equal(res$boundary_cov, expected)
})

test_that("Coverage virtual track handles chromosome transitions correctly", {
    # Test case 1:
    # An interval on chr1 followed by interval on chr2, querying the chr2 interval
    interv1 <- gintervals(2, 10, 20, 1)
    interv2 <- rbind(gintervals(1, 50, 80, 1), interv1)

    gvtrack.create("cov_bug", src = interv2, func = "coverage")
    result <- gextract("cov_bug", intervals = interv1, iterator = interv1)

    # Should be 1.0 (complete coverage) not 0
    expect_equal(result$cov_bug, 1.0,
        info = "Coverage should be 1.0 when querying an interval that exactly matches a source interval on chr2"
    )

    # Test case 2: Multiple intervals across multiple chromosomes with numeric iterator
    # Testing behavior when crossing multiple chromosome boundaries
    multi_chr_source <- rbind(
        gintervals(1, 50, 150, 1),
        gintervals(2, 10, 20, 1),
        gintervals(3, 30, 40, 1)
    )

    query_intervals <- rbind(
        gintervals(1, 0, 200, 1),
        gintervals(2, 0, 30, 1),
        gintervals(3, 0, 50, 1)
    )

    gvtrack.create("multi_chr_cov", src = multi_chr_source, func = "coverage")
    result <- gextract("multi_chr_cov", intervals = query_intervals, iterator = 10)

    # Calculate expected results with numeric iterator starting at beginning of each chromosome:
    # Chr1: [0-10], [10-20], ...[190-200] = 20 bins, with bins [50-60]...[140-150] having coverage 1.0
    # Chr2: [0-10], [10-20], [20-30] = 3 bins, with bin [10-20] having coverage 1.0
    # Chr3: [0-10], [10-20], [20-30], [30-40], [40-50] = 5 bins, with bin [30-40] having coverage 1.0

    chr1_expected <- numeric(20)
    chr1_expected[6:15] <- 1.0 # Bins 50-150 have full coverage

    chr2_expected <- c(0, 1, 0) # Only bin [10-20] has coverage

    chr3_expected <- c(0, 0, 0, 1, 0) # Only bin [30-40] has coverage

    expected <- c(chr1_expected, chr2_expected, chr3_expected)

    expect_equal(result$multi_chr_cov, expected,
        info = "Coverage with numeric iterator should correctly handle transitions across multiple chromosomes"
    )

    # Test case 3: Reversed chromosome order in source
    # Testing behavior when source intervals are in reverse chromosome order
    reverse_chr_source <- rbind(
        gintervals(3, 30, 40, 1),
        gintervals(2, 10, 20, 1),
        gintervals(1, 50, 150, 1)
    )

    gvtrack.create("reverse_chr_cov", src = reverse_chr_source, func = "coverage")
    result <- gextract("reverse_chr_cov", intervals = query_intervals, iterator = 10)

    # Should produce identical results regardless of source interval order
    expect_equal(result$reverse_chr_cov, expected,
        info = "Coverage should work correctly regardless of chromosome ordering in source intervals"
    )

    # Test case 4: Back-and-forth chromosome transitions with interval iterator
    # Create series of iterator intervals that jump between chromosomes
    zigzag_intervals <- rbind(
        gintervals(1, 50, 60, 1),
        gintervals(2, 10, 20, 1),
        gintervals(1, 100, 110, 1),
        gintervals(2, 20, 30, 1)
    )

    gvtrack.create("zigzag_cov", src = multi_chr_source, func = "coverage")
    result <- gextract("zigzag_cov", intervals = zigzag_intervals, iterator = zigzag_intervals)

    # Expected results should be:
    # Chr1 [50-60] = 1.0
    # Chr2 [10-20] = 1.0
    # Chr1 [100-110] = 1.0
    # Chr2 [20-30] = 0.0
    expected_zigzag <- c(1.0, 1.0, 1.0, 0.0)

    expect_equal(result$zigzag_cov, expected_zigzag,
        info = "Coverage should handle back-and-forth chromosome transitions correctly"
    )

    # Test case 5: Compare with distance function (which doesn't have the bug)
    # This helps verify that both functions handle chromosome transitions similarly after the fix
    gvtrack.create("dist_test", src = interv2, func = "distance")
    gvtrack.create("cov_test", src = interv2, func = "coverage")

    # Get both distance and coverage for same interval
    result_dist <- gextract("dist_test", intervals = interv1, iterator = interv1)
    result_cov <- gextract("cov_test", intervals = interv1, iterator = interv1)

    # Distance should be 0 (interval is exactly matched)
    # Coverage should be 1.0 (complete coverage)
    expect_equal(result_dist$dist_test, 0,
        info = "Distance should be 0 for exact interval match"
    )
    expect_equal(result_cov$cov_test, 1.0,
        info = "Coverage should be 1.0 for exact interval match"
    )

    # Test case 6: Handling non-aligned numeric iterator
    # Test intervals that don't align with the numeric iterator boundaries
    nonaligned_source <- gintervals(1, 15, 35, 1)
    nonaligned_query <- gintervals(1, 0, 50, 1)

    gvtrack.create("nonaligned_cov", src = nonaligned_source, func = "coverage")
    result <- gextract("nonaligned_cov", intervals = nonaligned_query, iterator = 10)

    # Expected results:
    # [0-10]: 0.0
    # [10-20]: 0.5 (15-20 covered = 5/10)
    # [20-30]: 1.0
    # [30-40]: 0.5 (30-35 covered = 5/10)
    # [40-50]: 0.0
    expected_nonaligned <- c(0.0, 0.5, 1.0, 0.5, 0.0)

    expect_equal(result$nonaligned_cov, expected_nonaligned,
        info = "Coverage should handle non-aligned numeric iterator boundaries correctly"
    )

    # Test case 7: Empty chromosomes in between
    # Tests behavior when there are chromosomes with no intervals
    sparse_source <- rbind(
        gintervals(1, 10, 20, 1),
        gintervals(5, 10, 20, 1) # Skip chr 2,3,4
    )

    sparse_query <- rbind(
        gintervals(1, 0, 30, 1),
        gintervals(2, 0, 30, 1), # No intervals here
        gintervals(5, 0, 30, 1)
    )

    gvtrack.create("sparse_cov", src = sparse_source, func = "coverage")
    result <- gextract("sparse_cov", intervals = sparse_query, iterator = 10)

    # Expected:
    # Chr1 [0-10], [10-20], [20-30] = 0, 1, 0
    # Chr2 [0-10], [10-20], [20-30] = 0, 0, 0
    # Chr5 [0-10], [10-20], [20-30] = 0, 1, 0
    expected_sparse <- c(0, 1, 0, 0, 0, 0, 0, 1, 0)

    expect_equal(result$sparse_cov, expected_sparse,
        info = "Coverage should handle skipped/empty chromosomes correctly"
    )

    # Test case 8: Partial query windows with explicit iterator
    # Testing behavior with specific iterator windows
    partial_source <- rbind(
        gintervals(1, 50, 60, 1),
        gintervals(2, 10, 20, 1)
    )

    # Create explicit iterator intervals to get precisely the windows we want
    explicit_iterator <- rbind(
        gintervals(1, 45, 55, 1),
        gintervals(1, 55, 65, 1),
        gintervals(2, 15, 25, 1)
    )

    gvtrack.create("partial_cov", src = partial_source, func = "coverage")
    result <- gextract("partial_cov", intervals = explicit_iterator, iterator = explicit_iterator)

    # Expected:
    # Chr1 [45-55] = 0.5 (50-55 covered = 5/10)
    # Chr1 [55-65] = 0.5 (55-60 covered = 5/10)
    # Chr2 [15-25] = 0.5 (15-20 covered = 5/10)
    expected_partial <- c(0.5, 0.5, 0.5)

    expect_equal(result$partial_cov, expected_partial,
        info = "Coverage should handle partial query windows correctly"
    )
})
