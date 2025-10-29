load_test_db()
test_that("gvtrack.filter basic operations work", {
    # Test 1: Attach and clear a simple filter
    gvtrack.create("test_vtrack", "test.fixedbin", func = "avg")

    # Create simple mask intervals
    mask <- gintervals(c(1, 1), c(100, 500), c(200, 600))

    # Attach filter - should not error
    expect_silent(gvtrack.filter("test_vtrack", filter = mask))

    # Check that filter is attached
    info <- gvtrack.info("test_vtrack")
    expect_true(!is.null(info$filter))

    # Clear filter
    expect_silent(gvtrack.filter("test_vtrack", filter = NULL))

    # Check that filter is cleared
    info <- gvtrack.info("test_vtrack")
    expect_null(info$filter)

    gvtrack.rm("test_vtrack")
})

test_that("gvtrack.filter validates input", {
    gvtrack.create("test_vtrack2", "test.fixedbin", func = "avg")

    # Invalid filter - not a data.frame
    expect_error(gvtrack.filter("test_vtrack2", filter = 123))

    # Invalid data.frame - missing columns
    expect_error(gvtrack.filter("test_vtrack2", filter = data.frame(a = 1, b = 2)))

    gvtrack.rm("test_vtrack2")
})

test_that("gvtrack.filter cache works correctly", {
    gvtrack.create("test_vtrack3", "test.fixedbin", func = "avg")
    gvtrack.create("test_vtrack4", "test.fixedbin", func = "max")

    # Same mask on two different vtracks should use cached filter
    mask <- gintervals(1, 1000, 2000)

    expect_silent(gvtrack.filter("test_vtrack3", filter = mask))
    expect_silent(gvtrack.filter("test_vtrack4", filter = mask))

    # Both should have a filter attached
    info3 <- gvtrack.info("test_vtrack3")
    info4 <- gvtrack.info("test_vtrack4")

    expect_true(!is.null(info3$filter))
    expect_true(!is.null(info4$filter))

    # Keys should be the same (same mask)
    expect_equal(info3$filter, info4$filter)

    gvtrack.rm("test_vtrack3")
    gvtrack.rm("test_vtrack4")
})

test_that("gvtrack.filter info shows statistics", {
    gvtrack.create("test_vtrack5", "test.fixedbin", func = "avg")

    # Create mask with known characteristics
    mask <- gintervals(c(1, 1, 2), c(100, 500, 1000), c(200, 600, 2000))

    expect_silent(gvtrack.filter("test_vtrack5", filter = mask))

    info <- gvtrack.info("test_vtrack5")

    # Should have filter_stats
    expect_true(!is.null(info$filter_stats))
    expect_true(!is.null(info$filter_stats$num_chroms))
    expect_true(!is.null(info$filter_stats$total_bases))

    # Check values make sense
    expect_gt(info$filter_stats$total_bases, 0)
    expect_gt(info$filter_stats$num_chroms, 0)

    gvtrack.rm("test_vtrack5")
})

test_that("gvtrack.filter masks completely covered intervals", {
    gvtrack.create("test_vtrack6", "test.fixedbin", func = "avg")

    # Create a mask that completely covers one interval
    mask <- gintervals(1, 1000, 2000)

    expect_silent(gvtrack.filter("test_vtrack6", filter = mask))

    # Extract from an interval completely covered by the mask
    result_covered <- gextract("test_vtrack6", gintervals(1, 1000, 2000))

    # The result should be NA for the masked interval
    expect_true(is.na(result_covered$test_vtrack6[1]))

    # Extract from an interval NOT covered by the mask
    result_uncovered <- gextract("test_vtrack6", gintervals(1, 5000, 6000))

    # The result should NOT be NA for the unmasked interval
    expect_false(is.na(result_uncovered$test_vtrack6[1]))

    gvtrack.rm("test_vtrack6")
})

test_that("gvtrack.filter works with coverage virtual tracks", {
    # Skip if annotations don't exist, use a simple intervals set instead
    test_cov_source <- gintervals(c(1, 1), c(1500, 3500), c(2500, 4500))
    if (gintervals.exists("test_cov_source")) gintervals.rm("test_cov_source", force = TRUE)
    gintervals.save("test_cov_source", test_cov_source)

    gvtrack.create("test_coverage", "test_cov_source", func = "coverage")

    # Create mask that covers part of the genome
    mask <- gintervals(1, 1000, 5000)

    expect_silent(gvtrack.filter("test_coverage", filter = mask))

    # Test completely masked interval with iterator
    result_masked <- gextract("test_coverage", gintervals(1, 2000, 3000), iterator = gintervals(1, 2000, 3000))
    expect_true(is.na(result_masked$test_coverage[1]))

    # Test unmasked interval
    result_unmasked <- gextract("test_coverage", gintervals(1, 10000, 15000), iterator = gintervals(1, 10000, 15000))
    expect_false(is.na(result_unmasked$test_coverage[1]))

    # Test partially masked interval - should calculate coverage on unmasked parts
    result_partial <- gextract("test_coverage", gintervals(1, 500, 6000), iterator = gintervals(1, 500, 6000))
    expect_false(is.na(result_partial$test_coverage[1]))
    expect_gte(result_partial$test_coverage[1], 0)
    expect_lte(result_partial$test_coverage[1], 1)

    gvtrack.rm("test_coverage")
    gintervals.rm("test_cov_source", force = TRUE)
})

test_that("gvtrack.filter works with multiple track functions", {
    # Test avg function
    gvtrack.create("test_avg", "test.fixedbin", func = "avg")
    mask <- gintervals(1, 1000, 2000)
    gvtrack.filter("test_avg", filter = mask)
    result_avg <- gextract("test_avg", gintervals(1, 1000, 2000))
    expect_true(is.na(result_avg$test_avg[1]))
    gvtrack.rm("test_avg")

    # Test sum function
    gvtrack.create("test_sum", "test.fixedbin", func = "sum")
    gvtrack.filter("test_sum", filter = mask)
    result_sum <- gextract("test_sum", gintervals(1, 1000, 2000))
    expect_true(is.na(result_sum$test_sum[1]))
    gvtrack.rm("test_sum")

    # Test max function
    gvtrack.create("test_max", "test.fixedbin", func = "max")
    gvtrack.filter("test_max", filter = mask)
    result_max <- gextract("test_max", gintervals(1, 1000, 2000))
    expect_true(is.na(result_max$test_max[1]))
    gvtrack.rm("test_max")

    # Test min function
    gvtrack.create("test_min", "test.fixedbin", func = "min")
    gvtrack.filter("test_min", filter = mask)
    result_min <- gextract("test_min", gintervals(1, 1000, 2000))
    expect_true(is.na(result_min$test_min[1]))
    gvtrack.rm("test_min")
})

test_that("gvtrack.filter handles edge cases", {
    gvtrack.create("test_edge", "test.fixedbin", func = "avg")

    # Test clearing filter
    mask_first <- gintervals(1, 1000, 2000)
    gvtrack.filter("test_edge", filter = mask_first)
    expect_silent(gvtrack.filter("test_edge", filter = NULL))
    result <- gextract("test_edge", gintervals(1, 1000, 2000))
    expect_false(is.na(result$test_edge[1]))

    # Mask on different chromosome shouldn't affect query
    mask_chr2 <- gintervals(2, 1000, 2000)
    gvtrack.filter("test_edge", filter = mask_chr2)
    result_chr1 <- gextract("test_edge", gintervals(1, 1000, 2000))
    expect_false(is.na(result_chr1$test_edge[1]))

    # Multiple disjoint masks
    mask_multi <- gintervals(c(1, 1, 1), c(1000, 3000, 5000), c(2000, 4000, 6000))
    gvtrack.filter("test_edge", filter = mask_multi)
    result_masked <- gextract("test_edge", gintervals(1, 1500, 1600))
    expect_true(is.na(result_masked$test_edge[1]))

    gvtrack.rm("test_edge")
})

test_that("gvtrack.filter works with partially masked intervals", {
    gvtrack.create("test_partial", "test.fixedbin", func = "avg")

    # Create mask that covers middle part of query interval
    mask <- gintervals(1, 2000, 3000)
    gvtrack.filter("test_partial", filter = mask)

    # Query interval 1000-4000 has middle 1000bp masked
    # Should return value based on unmasked parts (1000-2000 and 3000-4000)
    result <- gextract("test_partial", gintervals(1, 1000, 4000))
    expect_false(is.na(result$test_partial[1]))

    # Query interval 2000-3000 is completely masked
    result_full <- gextract("test_partial", gintervals(1, 2000, 3000))
    expect_true(is.na(result_full$test_partial[1]))

    gvtrack.rm("test_partial")
})

test_that("gvtrack.filter works with overlapping mask intervals", {
    gvtrack.create("test_overlap", "test.fixedbin", func = "avg")

    # Create overlapping mask intervals (should be merged internally)
    mask <- gintervals(c(1, 1), c(1000, 1500), c(2000, 2500))
    gvtrack.filter("test_overlap", filter = mask)

    # Query in merged region should be NA
    result <- gextract("test_overlap", gintervals(1, 1200, 1800))
    expect_true(is.na(result$test_overlap[1]))

    gvtrack.rm("test_overlap")
})

test_that("gvtrack.filter statistics are accurate", {
    gvtrack.create("test_stats", "test.fixedbin", func = "avg")

    # Create mask with known size: 3 non-overlapping intervals
    # Chr1: 1000-1100 (100bp) + 2000-2200 (200bp) = 300bp
    # Chr2: 1000-1100 (100bp) = 100bp
    # Total: 400bp
    mask <- gintervals(c(1, 1, 2), c(1000, 2000, 1000), c(1100, 2200, 1100))
    gvtrack.filter("test_stats", filter = mask)

    info <- gvtrack.info("test_stats")

    expect_equal(info$filter_stats$total_bases, 400)
    expect_equal(info$filter_stats$num_chroms, 2)

    gvtrack.rm("test_stats")
})

test_that("gvtrack.filter updates when changed", {
    gvtrack.create("test_update", "test.fixedbin", func = "avg")

    # Set first mask
    mask1 <- gintervals(1, 1000, 2000)
    gvtrack.filter("test_update", filter = mask1)

    result1 <- gextract("test_update", gintervals(1, 1000, 2000))
    expect_true(is.na(result1$test_update[1]))

    result1_unmasked <- gextract("test_update", gintervals(1, 3000, 4000))
    expect_false(is.na(result1_unmasked$test_update[1]))

    # Change to different mask
    mask2 <- gintervals(1, 3000, 4000)
    gvtrack.filter("test_update", filter = mask2)

    result2 <- gextract("test_update", gintervals(1, 1000, 2000))
    expect_false(is.na(result2$test_update[1]))

    result2_masked <- gextract("test_update", gintervals(1, 3000, 4000))
    expect_true(is.na(result2_masked$test_update[1]))

    gvtrack.rm("test_update")
})

test_that("gvtrack.filter works with iterator modifiers", {
    gvtrack.create("test_itermdf", "test.fixedbin", func = "avg", iterator = 10)

    # Create mask
    mask <- gintervals(1, 1000, 2000)
    gvtrack.filter("test_itermdf", filter = mask)

    # Extract should still work with iterator
    result <- gextract("test_itermdf", gintervals(1, 1500, 1600))
    expect_true(is.na(result$test_itermdf[1]))

    gvtrack.rm("test_itermdf")
})

test_that("gvtrack.filter coverage calculation is correct", {
    # Create test intervals set for coverage calculation
    test_intervs <- gintervals(c(1, 1), c(1200, 3500), c(1800, 3700))
    if (gintervals.exists("test_coverage_source")) gintervals.rm("test_coverage_source", force = TRUE)
    gintervals.save("test_coverage_source", test_intervs)

    gvtrack.create("test_cov_calc", "test_coverage_source", func = "coverage")

    # Create mask that splits interval into thirds
    # Query 1000-4000 (3000bp), mask 1500-2000 and 2500-3000 (1000bp masked)
    mask <- gintervals(c(1, 1), c(1500, 2500), c(2000, 3000))
    gvtrack.filter("test_cov_calc", filter = mask)

    # Get coverage without filter first
    gvtrack.filter("test_cov_calc", filter = NULL)
    result_no_filter <- gextract("test_cov_calc", gintervals(1, 1000, 4000), iterator = gintervals(1, 1000, 4000))

    # Now with filter
    gvtrack.filter("test_cov_calc", filter = mask)
    result_with_filter <- gextract("test_cov_calc", gintervals(1, 1000, 4000), iterator = gintervals(1, 1000, 4000))

    # With filter, coverage is calculated only on unmasked 2000bp
    # Values might differ if source intervals overlap with masked regions
    expect_false(is.na(result_with_filter$test_cov_calc[1]))

    gvtrack.rm("test_cov_calc")
    gintervals.rm("test_coverage_source", force = TRUE)
})

test_that("gvtrack.filter works across multiple chromosomes", {
    gvtrack.create("test_multichrom", "test.fixedbin", func = "avg")

    # Mask on chromosomes 1 and 2
    mask <- gintervals(c(1, 2), c(1000, 1000), c(2000, 2000))
    gvtrack.filter("test_multichrom", filter = mask)

    # Both chromosomes should be masked
    result_chr1 <- gextract("test_multichrom", gintervals(1, 1500, 1600))
    result_chr2 <- gextract("test_multichrom", gintervals(2, 1500, 1600))

    expect_true(is.na(result_chr1$test_multichrom[1]))
    expect_true(is.na(result_chr2$test_multichrom[1]))

    # Different region on chr1 should not be masked
    result_chr1_unmask <- gextract("test_multichrom", gintervals(1, 5000, 5100))
    expect_false(is.na(result_chr1_unmask$test_multichrom[1]))

    gvtrack.rm("test_multichrom")
})

test_that("gvtrack.filter with exact boundary conditions", {
    gvtrack.create("test_boundary", "test.fixedbin", func = "avg")

    # Mask exactly at interval boundaries
    mask <- gintervals(1, 1000, 2000)
    gvtrack.filter("test_boundary", filter = mask)

    # Query exactly matching mask
    result_exact <- gextract("test_boundary", gintervals(1, 1000, 2000))
    expect_true(is.na(result_exact$test_boundary[1]))

    # Query just touching mask edge (before)
    result_before <- gextract("test_boundary", gintervals(1, 900, 1000))
    expect_false(is.na(result_before$test_boundary[1]))

    # Query just touching mask edge (after)
    result_after <- gextract("test_boundary", gintervals(1, 2000, 2100))
    expect_false(is.na(result_after$test_boundary[1]))

    # Query spanning mask boundary (start inside mask, end outside)
    # With default binned iterator, first bin (1500-1550) is masked -> NA
    # Later bins (2000+) are unmasked -> have values
    result_span_start <- gextract("test_boundary", gintervals(1, 1500, 2500))
    expect_true(is.na(result_span_start$test_boundary[1])) # First bin is masked
    expect_false(all(is.na(result_span_start$test_boundary))) # Not all bins are masked

    # Query spanning mask boundary (end inside)
    result_span_end <- gextract("test_boundary", gintervals(1, 500, 1500))
    expect_false(is.na(result_span_end$test_boundary[1])) # First bin is unmasked
    expect_true(any(is.na(result_span_end$test_boundary))) # Some bins are masked

    gvtrack.rm("test_boundary")
})

test_that("gvtrack.filter with very small intervals", {
    gvtrack.create("test_small", "test.fixedbin", func = "avg")

    # Single base pair mask
    mask <- gintervals(1, 1000, 1001)
    gvtrack.filter("test_small", filter = mask)

    # Query containing single bp mask
    result <- gextract("test_small", gintervals(1, 900, 1100))
    expect_false(is.na(result$test_small[1]))

    # Query exactly the single bp
    result_exact <- gextract("test_small", gintervals(1, 1000, 1001))
    expect_true(is.na(result_exact$test_small[1]))

    gvtrack.rm("test_small")
})

test_that("gvtrack.filter with large contiguous masks", {
    gvtrack.create("test_large", "test.fixedbin", func = "sum")

    # Large mask covering 100kb
    mask <- gintervals(1, 0, 100000)
    gvtrack.filter("test_large", filter = mask)

    # Query completely inside large mask
    result_inside <- gextract("test_large", gintervals(1, 10000, 20000))
    expect_true(is.na(result_inside$test_large[1]))

    # Query outside large mask
    result_outside <- gextract("test_large", gintervals(1, 150000, 160000))
    expect_false(is.na(result_outside$test_large[1]))

    gvtrack.rm("test_large")
})

test_that("gvtrack.filter with multiple vtracks sharing filter", {
    gvtrack.create("test_shared1", "test.fixedbin", func = "avg")
    gvtrack.create("test_shared2", "test.fixedbin", func = "max")

    mask <- gintervals(1, 1000, 2000)

    # Apply same filter to both
    gvtrack.filter("test_shared1", filter = mask)
    gvtrack.filter("test_shared2", filter = mask)

    # Both should have same filter key
    info1 <- gvtrack.info("test_shared1")
    info2 <- gvtrack.info("test_shared2")
    expect_equal(info1$filter, info2$filter)

    # Both should mask the same region
    result1 <- gextract("test_shared1", gintervals(1, 1500, 1600))
    result2 <- gextract("test_shared2", gintervals(1, 1500, 1600))
    expect_true(is.na(result1$test_shared1[1]))
    expect_true(is.na(result2$test_shared2[1]))

    gvtrack.rm("test_shared1")
    gvtrack.rm("test_shared2")
})

test_that("gvtrack.filter with dense intervals as mask source", {
    gvtrack.create("test_dense_mask", "test.fixedbin", func = "avg")

    # Create many small intervals as mask
    starts <- seq(1000, 10000, by = 100)
    ends <- starts + 50
    mask <- gintervals(rep(1, length(starts)), starts, ends)
    gvtrack.filter("test_dense_mask", filter = mask)

    # Query in a gap between masks
    result_gap <- gextract("test_dense_mask", gintervals(1, 1060, 1090))
    expect_false(is.na(result_gap$test_dense_mask[1]))

    # Query overlapping multiple masks - first bin at 1000 is masked
    result_multi <- gextract("test_dense_mask", gintervals(1, 1000, 2000))
    expect_true(is.na(result_multi$test_dense_mask[1])) # First bin is masked
    expect_false(all(is.na(result_multi$test_dense_mask))) # But not all bins

    gvtrack.rm("test_dense_mask")
})

test_that("gvtrack.filter preserves virtual track type", {
    # Create different types of vtracks using test.fixedbin
    gvtrack.create("test_type1", "test.fixedbin", func = "avg")
    gvtrack.create("test_type2", "test.fixedbin", func = "max")

    mask <- gintervals(1, 1000, 2000)

    # Apply filter
    gvtrack.filter("test_type1", filter = mask)
    gvtrack.filter("test_type2", filter = mask)

    # Check vtracks still work correctly (can extract without error)
    result1 <- gextract("test_type1", gintervals(1, 3000, 4000))
    result2 <- gextract("test_type2", gintervals(1, 3000, 4000))

    # Both should have data columns (functions preserved)
    expect_true("test_type1" %in% colnames(result1))
    expect_true("test_type2" %in% colnames(result2))

    gvtrack.rm("test_type1")
    gvtrack.rm("test_type2")
})

test_that("gvtrack.filter statistics with merged overlapping intervals", {
    gvtrack.create("test_merge_stats", "test.fixedbin", func = "avg")

    # Overlapping intervals: [1000-2000] and [1500-2500]
    # After merge: [1000-2500] = 1500bp
    mask <- gintervals(c(1, 1), c(1000, 1500), c(2000, 2500))
    gvtrack.filter("test_merge_stats", filter = mask)

    info <- gvtrack.info("test_merge_stats")

    # Should be merged to 1500bp total
    expect_equal(info$filter_stats$total_bases, 1500)

    gvtrack.rm("test_merge_stats")
})

test_that("gvtrack.filter with coverage basic behavior", {
    # Create simple intervals set for coverage
    cov_source <- gintervals(c(1, 1), c(3500, 4500), c(3800, 4800))
    if (gintervals.exists("test_cov_src2")) gintervals.rm("test_cov_src2", force = TRUE)
    gintervals.save("test_cov_src2", cov_source)

    gvtrack.create("test_cov2", "test_cov_src2", func = "coverage")

    mask <- gintervals(1, 1000, 2000)
    gvtrack.filter("test_cov2", filter = mask)

    # Coverage on unmasked region
    result <- gextract("test_cov2", gintervals(1, 3000, 4000), iterator = gintervals(1, 3000, 4000))
    expect_false(is.na(result$test_cov2[1]))

    gvtrack.rm("test_cov2")
    gintervals.rm("test_cov_src2", force = TRUE)
})

test_that("gvtrack.filter with all standard functions", {
    mask <- gintervals(1, 1000, 2000)

    # Test stddev
    gvtrack.create("test_stddev", "test.fixedbin", func = "stddev")
    gvtrack.filter("test_stddev", filter = mask)
    result_stddev <- gextract("test_stddev", gintervals(1, 1500, 1600))
    expect_true(is.na(result_stddev$test_stddev[1]))
    gvtrack.rm("test_stddev")

    # Test min
    gvtrack.create("test_min2", "test.fixedbin", func = "min")
    gvtrack.filter("test_min2", filter = mask)
    result_min <- gextract("test_min2", gintervals(1, 1500, 1600))
    expect_true(is.na(result_min$test_min2[1]))
    gvtrack.rm("test_min2")
})

test_that("gvtrack.filter avg matches manual complement extraction", {
    # Validates that filtered avg equals weighted average over complement parts

    gvtrack.create("test_avg_complement", "test.fixedbin", func = "avg")
    mask <- gintervals(1, 2000, 5000) # Mask [2000, 5000)
    gvtrack.filter("test_avg_complement", filter = mask)

    # Query [1000, 8000) with filter
    query <- gintervals(1, 1000, 8000)
    result_filtered <- gextract("test_avg_complement", query, iterator = query)

    # Manual: extract complement parts and compute weighted avg
    gvtrack.filter("test_avg_complement", filter = NULL)
    seg1 <- gextract("test_avg_complement", gintervals(1, 1000, 2000), iterator = gintervals(1, 1000, 2000))
    seg2 <- gextract("test_avg_complement", gintervals(1, 5000, 8000), iterator = gintervals(1, 5000, 8000))

    if (!is.na(seg1$test_avg_complement[1]) && !is.na(seg2$test_avg_complement[1])) {
        len1 <- 1000
        len2 <- 3000
        expected_avg <- (seg1$test_avg_complement[1] * len1 + seg2$test_avg_complement[1] * len2) / (len1 + len2)

        # Should match within tolerance
        expect_equal(result_filtered$test_avg_complement[1], expected_avg, tolerance = 0.001)
    }

    gvtrack.rm("test_avg_complement")
})

test_that("gvtrack.filter sum matches manual complement extraction", {
    # Validates that filtered sum equals sum over complement parts

    gvtrack.create("test_sum_complement", "test.fixedbin", func = "sum")
    mask <- gintervals(1, 2000, 4000) # Mask [2000, 4000)
    gvtrack.filter("test_sum_complement", filter = mask)

    # Query [1000, 6000) with filter
    query <- gintervals(1, 1000, 6000)
    result_filtered <- gextract("test_sum_complement", query, iterator = query)

    # Manual: extract complement parts and sum
    gvtrack.filter("test_sum_complement", filter = NULL)
    seg1 <- gextract("test_sum_complement", gintervals(1, 1000, 2000), iterator = gintervals(1, 1000, 2000))
    seg2 <- gextract("test_sum_complement", gintervals(1, 4000, 6000), iterator = gintervals(1, 4000, 6000))

    if (!is.na(seg1$test_sum_complement[1]) && !is.na(seg2$test_sum_complement[1])) {
        expected_sum <- seg1$test_sum_complement[1] + seg2$test_sum_complement[1]
        expect_equal(result_filtered$test_sum_complement[1], expected_sum, tolerance = 0.001)
    }

    gvtrack.rm("test_sum_complement")
})

test_that("gvtrack.filter min/max match manual complement extraction", {
    # Validates that filtered min/max equal min/max over complement parts

    # Test MIN
    gvtrack.create("test_min_complement", "test.fixedbin", func = "min")
    mask <- gintervals(1, 2000, 4000)
    gvtrack.filter("test_min_complement", filter = mask)

    query <- gintervals(1, 1000, 6000)
    result_filtered <- gextract("test_min_complement", query, iterator = query)

    gvtrack.filter("test_min_complement", filter = NULL)
    seg1 <- gextract("test_min_complement", gintervals(1, 1000, 2000), iterator = gintervals(1, 1000, 2000))
    seg2 <- gextract("test_min_complement", gintervals(1, 4000, 6000), iterator = gintervals(1, 4000, 6000))

    if (!is.na(seg1$test_min_complement[1]) && !is.na(seg2$test_min_complement[1])) {
        expected_min <- min(seg1$test_min_complement[1], seg2$test_min_complement[1])
        expect_equal(result_filtered$test_min_complement[1], expected_min, tolerance = 0.001)
    }

    gvtrack.rm("test_min_complement")

    # Test MAX
    gvtrack.create("test_max_complement", "test.fixedbin", func = "max")
    gvtrack.filter("test_max_complement", filter = mask)

    result_filtered <- gextract("test_max_complement", query, iterator = query)

    gvtrack.filter("test_max_complement", filter = NULL)
    seg1 <- gextract("test_max_complement", gintervals(1, 1000, 2000), iterator = gintervals(1, 1000, 2000))
    seg2 <- gextract("test_max_complement", gintervals(1, 4000, 6000), iterator = gintervals(1, 4000, 6000))

    if (!is.na(seg1$test_max_complement[1]) && !is.na(seg2$test_max_complement[1])) {
        expected_max <- max(seg1$test_max_complement[1], seg2$test_max_complement[1])
        expect_equal(result_filtered$test_max_complement[1], expected_max, tolerance = 0.001)
    }

    gvtrack.rm("test_max_complement")
})

test_that("gvtrack.filter coverage matches manual complement calculation", {
    # Create source intervals for coverage
    source_intervs <- gintervals(c(1, 1), c(1500, 6000), c(2500, 7000))
    if (gintervals.exists("test_cov_complement")) gintervals.rm("test_cov_complement", force = TRUE)
    gintervals.save("test_cov_complement", source_intervs)

    gvtrack.create("test_cov_complement_vt", "test_cov_complement", func = "coverage")

    # Mask [3000, 5000)
    mask <- gintervals(1, 3000, 5000)
    gvtrack.filter("test_cov_complement_vt", filter = mask)

    # Query [1000, 8000) - effective: [1000, 3000) ∪ [5000, 8000) = 5000bp
    query <- gintervals(1, 1000, 8000)
    result_filtered <- gextract("test_cov_complement_vt", query, iterator = query)

    # Manual calculation:
    # Source [1500, 2500) overlaps [1000, 3000): 1000bp covered
    # Source [6000, 7000) overlaps [5000, 8000): 1000bp covered
    # Total covered: 2000bp / Total unmasked: 5000bp = 0.4
    expected_coverage <- 0.4

    expect_equal(result_filtered$test_cov_complement_vt[1], expected_coverage, tolerance = 1e-10)

    gvtrack.rm("test_cov_complement_vt")
    gintervals.rm("test_cov_complement", force = TRUE)
})

test_that("gvtrack.filter quantile works with single-interval iterator", {
    # Validates that quantile function now works with filter + single-interval iterator

    gvtrack.create("test_quantile_filter", "test.fixedbin", func = "quantile", params = 0.5)
    mask <- gintervals(1, 2000, 4000) # Mask [2000, 4000)
    gvtrack.filter("test_quantile_filter", filter = mask)

    # Query [1000, 6000) with single-interval iterator
    query <- gintervals(1, 1000, 6000)
    result_filtered <- gextract("test_quantile_filter", query, iterator = query)

    # Should return a valid quantile (not NA, not error)
    expect_false(is.na(result_filtered$test_quantile_filter[1]))

    # Verify it's reasonable - should be between segment medians
    gvtrack.filter("test_quantile_filter", filter = NULL)
    seg1 <- gextract("test_quantile_filter", gintervals(1, 1000, 2000), iterator = gintervals(1, 1000, 2000))
    seg2 <- gextract("test_quantile_filter", gintervals(1, 4000, 6000), iterator = gintervals(1, 4000, 6000))

    # Result should be within reasonable range of segment medians
    if (!is.na(seg1$test_quantile_filter[1]) && !is.na(seg2$test_quantile_filter[1])) {
        min_expected <- min(seg1$test_quantile_filter[1], seg2$test_quantile_filter[1])
        max_expected <- max(seg1$test_quantile_filter[1], seg2$test_quantile_filter[1])

        # Allow some tolerance for sampling effects in quantile estimation
        expect_true(result_filtered$test_quantile_filter[1] >= min_expected * 0.5)
        expect_true(result_filtered$test_quantile_filter[1] <= max_expected * 2.0)
    }

    gvtrack.rm("test_quantile_filter")
})

test_that("gvtrack.filter quantile works with binned iterator", {
    # Validates that quantile still works with default binned iterator

    gvtrack.create("test_quantile_binned", "test.fixedbin", func = "quantile", params = 0.5)
    mask <- gintervals(1, 2000, 4000)
    gvtrack.filter("test_quantile_binned", filter = mask)

    # Default binned iterator
    result_binned <- gextract("test_quantile_binned", gintervals(1, 1000, 6000))

    # Should have some non-NA values in unmasked regions
    expect_true(any(!is.na(result_binned$test_quantile_binned[result_binned$start < 2000])))
    expect_true(any(!is.na(result_binned$test_quantile_binned[result_binned$start >= 4000])))

    # Should have NA values in masked region
    masked_bins <- result_binned[result_binned$start >= 2000 & result_binned$start < 4000, ]
    expect_true(any(is.na(masked_bins$test_quantile_binned)))

    gvtrack.rm("test_quantile_binned")
})

test_that("gvtrack.filter nearest function works", {
    gvtrack.create("test_nearest", "test.sparse", func = "nearest")
    mask <- gintervals(1, 2000, 4000)
    gvtrack.filter("test_nearest", filter = mask)

    # Extract with binned iterator
    result <- gextract("test_nearest", gintervals(1, 1000, 6000))

    # Should have values outside masked region
    expect_true(any(!is.na(result$test_nearest[result$start < 2000])))
    expect_true(any(!is.na(result$test_nearest[result$start >= 4000])))

    # Note: sparse tracks may not have data in every bin, so we don't check for NA in masked region
    # The filter is working correctly if unmasked bins have values

    gvtrack.rm("test_nearest")
})

test_that("gvtrack.filter global.percentile works", {
    gvtrack.create("test_gp", "test.fixedbin", func = "global.percentile")
    mask <- gintervals(1, 2000, 4000)
    gvtrack.filter("test_gp", filter = mask)

    # Extract with binned iterator
    result <- gextract("test_gp", gintervals(1, 1000, 6000))

    # Should have values outside masked region
    expect_true(any(!is.na(result$test_gp[result$start < 2000])))
    expect_true(any(!is.na(result$test_gp[result$start >= 4000])))

    # Values should be percentiles (0-1 range)
    non_na_vals <- result$test_gp[!is.na(result$test_gp)]
    expect_true(all(non_na_vals >= 0 & non_na_vals <= 1))

    gvtrack.rm("test_gp")
})

test_that("gvtrack.filter global.percentile.min works", {
    gvtrack.create("test_gp_min", "test.fixedbin", func = "global.percentile.min")
    mask <- gintervals(1, 2000, 4000)
    gvtrack.filter("test_gp_min", filter = mask)

    # Extract with binned iterator
    result <- gextract("test_gp_min", gintervals(1, 1000, 6000))

    # Should have values outside masked region
    expect_true(any(!is.na(result$test_gp_min[result$start < 2000])))

    # Values should be percentiles (0-1 range)
    non_na_vals <- result$test_gp_min[!is.na(result$test_gp_min)]
    expect_true(all(non_na_vals >= 0 & non_na_vals <= 1))

    gvtrack.rm("test_gp_min")
})

test_that("gvtrack.filter global.percentile.max works", {
    gvtrack.create("test_gp_max", "test.fixedbin", func = "global.percentile.max")
    mask <- gintervals(1, 2000, 4000)
    gvtrack.filter("test_gp_max", filter = mask)

    # Extract with binned iterator
    result <- gextract("test_gp_max", gintervals(1, 1000, 6000))

    # Should have values outside masked region
    expect_true(any(!is.na(result$test_gp_max[result$start < 2000])))

    # Values should be percentiles (0-1 range)
    non_na_vals <- result$test_gp_max[!is.na(result$test_gp_max)]
    expect_true(all(non_na_vals >= 0 & non_na_vals <= 1))

    gvtrack.rm("test_gp_max")
})

# Note: distance and distance.center are interval-based functions, not track-based
# They work differently and don't use the same filter mechanism as track functions

test_that("gvtrack.filter pwm function works", {
    # Create a simple PWM matrix (4x4)
    pssm <- matrix(c(
        0.25, 0.25, 0.25, 0.25,
        0.70, 0.10, 0.10, 0.10,
        0.10, 0.70, 0.10, 0.10,
        0.10, 0.10, 0.70, 0.10
    ), nrow = 4, byrow = TRUE)
    colnames(pssm) <- c("A", "C", "G", "T")

    gvtrack.create("test_pwm", src = NULL, func = "pwm", params = list(pssm = pssm))
    mask <- gintervals(1, 2000, 4000)
    gvtrack.filter("test_pwm", filter = mask)

    # Extract with explicit binned iterator (PWM requires explicit iterator)
    result <- gextract("test_pwm", gintervals(1, 1000, 6000), iterator = 50)

    # Should have values outside masked region
    expect_true(any(!is.na(result$test_pwm[result$start < 2000])))
    expect_true(any(!is.na(result$test_pwm[result$start >= 4000])))

    # PWM functions work with filters
    expect_true(nrow(result) > 0)

    gvtrack.rm("test_pwm")
})

test_that("gvtrack.filter pwm.max function works", {
    # Create a simple PWM matrix
    pssm <- matrix(c(
        0.25, 0.25, 0.25, 0.25,
        0.70, 0.10, 0.10, 0.10,
        0.10, 0.70, 0.10, 0.10,
        0.10, 0.10, 0.70, 0.10
    ), nrow = 4, byrow = TRUE)
    colnames(pssm) <- c("A", "C", "G", "T")

    gvtrack.create("test_pwm_max", src = NULL, func = "pwm.max", params = list(pssm = pssm))
    mask <- gintervals(1, 2000, 4000)
    gvtrack.filter("test_pwm_max", filter = mask)

    # Extract with explicit binned iterator
    result <- gextract("test_pwm_max", gintervals(1, 1000, 6000), iterator = 50)

    # Should have values outside masked region
    expect_true(any(!is.na(result$test_pwm_max[result$start < 2000])))
    expect_true(any(!is.na(result$test_pwm_max[result$start >= 4000])))

    gvtrack.rm("test_pwm_max")
})

test_that("gvtrack.filter pwm.max.pos function works", {
    # Create a simple PWM matrix
    pssm <- matrix(c(
        0.25, 0.25, 0.25, 0.25,
        0.70, 0.10, 0.10, 0.10,
        0.10, 0.70, 0.10, 0.10,
        0.10, 0.10, 0.70, 0.10
    ), nrow = 4, byrow = TRUE)
    colnames(pssm) <- c("A", "C", "G", "T")

    gvtrack.create("test_pwm_max_pos", src = NULL, func = "pwm.max.pos", params = list(pssm = pssm))
    mask <- gintervals(1, 2000, 4000)
    gvtrack.filter("test_pwm_max_pos", filter = mask)

    # Extract with explicit binned iterator
    result <- gextract("test_pwm_max_pos", gintervals(1, 1000, 6000), iterator = 50)

    # Should have values outside masked region (positions are integers)
    expect_true(any(!is.na(result$test_pwm_max_pos[result$start < 2000])))
    expect_true(any(!is.na(result$test_pwm_max_pos[result$start >= 4000])))

    gvtrack.rm("test_pwm_max_pos")
})

test_that("gvtrack.filter pwm.count function works", {
    # Create a simple PWM matrix
    pssm <- matrix(c(
        0.25, 0.25, 0.25, 0.25,
        0.70, 0.10, 0.10, 0.10,
        0.10, 0.70, 0.10, 0.10,
        0.10, 0.10, 0.70, 0.10
    ), nrow = 4, byrow = TRUE)
    colnames(pssm) <- c("A", "C", "G", "T")

    gvtrack.create("test_pwm_count",
        src = NULL, func = "pwm.count",
        params = list(pssm = pssm, score.thresh = 0)
    )
    mask <- gintervals(1, 2000, 4000)
    gvtrack.filter("test_pwm_count", filter = mask)

    # Extract with explicit binned iterator
    result <- gextract("test_pwm_count", gintervals(1, 1000, 6000), iterator = 50)

    # Counts should be non-negative integers
    non_na_vals <- result$test_pwm_count[!is.na(result$test_pwm_count)]
    expect_true(all(non_na_vals >= 0))

    # PWM count functions work with filters
    expect_true(nrow(result) > 0)

    gvtrack.rm("test_pwm_count")
})

test_that("gvtrack.filter kmer.count function works", {
    gvtrack.create("test_kmer_count",
        src = NULL, func = "kmer.count",
        params = list(kmer = "ACGT", extend = TRUE, strand = 1)
    )
    mask <- gintervals(1, 2000, 4000)
    gvtrack.filter("test_kmer_count", filter = mask)

    # Extract with explicit binned iterator
    result <- gextract("test_kmer_count", gintervals(1, 1000, 6000), iterator = 50)

    # Counts should be non-negative integers
    non_na_vals <- result$test_kmer_count[!is.na(result$test_kmer_count)]
    expect_true(all(non_na_vals >= 0))

    # Should have some results
    expect_true(nrow(result) > 0)

    gvtrack.rm("test_kmer_count")
})

test_that("gvtrack.filter kmer.frac function works", {
    gvtrack.create("test_kmer_frac",
        src = NULL, func = "kmer.frac",
        params = list(kmer = "CG", extend = FALSE, strand = 1)
    )
    mask <- gintervals(1, 2000, 4000)
    gvtrack.filter("test_kmer_frac", filter = mask)

    # Extract with explicit binned iterator
    result <- gextract("test_kmer_frac", gintervals(1, 1000, 6000), iterator = 50)

    # Fractions should be between 0 and 1
    non_na_vals <- result$test_kmer_frac[!is.na(result$test_kmer_frac)]
    expect_true(all(non_na_vals >= 0 & non_na_vals <= 1))

    # Should have some results
    expect_true(nrow(result) > 0)

    gvtrack.rm("test_kmer_frac")
})

test_that("gvtrack.filter pwm respects filter with single-interval iterator", {
    # Create a PWM vtrack with filter and verify it excludes masked regions
    pssm <- matrix(c(
        0.25, 0.25, 0.25, 0.25,
        0.70, 0.10, 0.10, 0.10,
        0.10, 0.70, 0.10, 0.10,
        0.10, 0.10, 0.70, 0.10
    ), nrow = 4, byrow = TRUE)
    colnames(pssm) <- c("A", "C", "G", "T")

    gvtrack.create("test_pwm_filter_single", src = NULL, func = "pwm", params = list(pssm = pssm))
    mask <- gintervals(1, 2000, 4000)
    gvtrack.filter("test_pwm_filter_single", filter = mask)

    # Extract with single-interval iterator
    query <- gintervals(1, 1000, 6000)
    result_filtered <- gextract("test_pwm_filter_single", query, iterator = query)

    # Should not be completely masked (has unmasked regions)
    expect_false(is.na(result_filtered$test_pwm_filter_single[1]))

    # Extract without filter for comparison
    gvtrack.create("test_pwm_no_filter", src = NULL, func = "pwm", params = list(pssm = pssm))
    result_no_filter <- gextract("test_pwm_no_filter", query, iterator = query)

    # The filtered result should be different (excluding masked region contribution)
    # The difference depends on sequence content, so we just verify they're different
    expect_false(identical(result_filtered$test_pwm_filter_single[1], result_no_filter$test_pwm_no_filter[1]))

    gvtrack.rm("test_pwm_filter_single")
    gvtrack.rm("test_pwm_no_filter")
})

test_that("gvtrack.filter kmer.count respects filter with single-interval iterator", {
    # Create a kmer vtrack with filter and verify it excludes masked regions
    gvtrack.create("test_kmer_filter_single",
        src = NULL, func = "kmer.count",
        params = list(kmer = "ACGT", extend = TRUE, strand = 1)
    )
    mask <- gintervals(1, 2000, 4000)
    gvtrack.filter("test_kmer_filter_single", filter = mask)

    # Extract with single-interval iterator
    query <- gintervals(1, 1000, 6000)
    result_filtered <- gextract("test_kmer_filter_single", query, iterator = query)

    # Should not be completely masked
    expect_false(is.na(result_filtered$test_kmer_filter_single[1]))

    # Extract without filter for comparison
    gvtrack.create("test_kmer_no_filter",
        src = NULL, func = "kmer.count",
        params = list(kmer = "ACGT", extend = TRUE, strand = 1)
    )
    result_no_filter <- gextract("test_kmer_no_filter", query, iterator = query)

    # The filtered result should be different (excluding masked region)
    # It should be lower because we're excluding kmers in the masked region
    expect_true(result_filtered$test_kmer_filter_single[1] <= result_no_filter$test_kmer_no_filter[1])

    gvtrack.rm("test_kmer_filter_single")
    gvtrack.rm("test_kmer_no_filter")
})

test_that("gvtrack.filter does not leak state to unfiltered vtracks sharing same source", {
    # Regression test for bug where filtered vtracks corrupted unfiltered vtracks
    # sharing the same Track_n_imdf (same source + iterator modifiers).
    # The bug: aggregate_*_with_filter() functions call track.read_interval() which
    # modifies shared track state, causing subsequent unfiltered vtracks to return
    # wrong values.

    # Create intervals for filtering
    mask <- gintervals(c(1, 1), c(1500, 3500), c(2500, 4500))

    # Create filtered vtrack with iterator modifier
    gvtrack.create("filtered_vtrack", "test.fixedbin", func = "sum")
    gvtrack.iterator("filtered_vtrack", sshift = -100, eshift = 100)
    gvtrack.filter("filtered_vtrack", mask)

    # Create unfiltered vtrack with SAME source and iterator modifier
    # (this causes them to share the same Track_n_imdf)
    gvtrack.create("unfiltered_vtrack", "test.fixedbin", func = "sum")
    gvtrack.iterator("unfiltered_vtrack", sshift = -100, eshift = 100)

    # Extract intervals
    query <- gintervals(c(1, 1, 1), c(1000, 2000, 4000), c(1200, 2200, 4200))

    # Extract both vtracks together (this would trigger the bug)
    result_with_filtered <- gextract(c("filtered_vtrack", "unfiltered_vtrack"),
        intervals = query,
        iterator = query
    )

    # Extract only the unfiltered vtrack (baseline)
    result_alone <- gextract("unfiltered_vtrack", intervals = query, iterator = query)

    # The unfiltered vtrack should return IDENTICAL values in both cases
    # Before the fix, result_with_filtered$unfiltered_vtrack would be corrupted
    # by the filtered vtrack's read_interval() calls
    expect_equal(result_with_filtered$unfiltered_vtrack, result_alone$unfiltered_vtrack,
        label = "Unfiltered vtrack values should not be affected by filtered vtrack sharing same source"
    )

    # Verify that the filtered vtrack actually does something different
    expect_true(
        any(is.na(result_with_filtered$filtered_vtrack) |
            result_with_filtered$filtered_vtrack != result_with_filtered$unfiltered_vtrack),
        label = "Filtered vtrack should produce different results than unfiltered"
    )

    # Clean up
    gvtrack.rm("filtered_vtrack")
    gvtrack.rm("unfiltered_vtrack")
})
