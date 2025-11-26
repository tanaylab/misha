create_isolated_test_db()

test_that("TSS directional neighbors work correctly", {
    # Create simple test data
    tss <- data.frame(
        chrom = c("chr1", "chr1", "chr1"),
        start = c(1000, 2000, 3000),
        end = c(1001, 2001, 3001),
        strand = c(1, -1, 1), # +, -, +
        gene = c("GeneA", "GeneB", "GeneC")
    )

    features <- data.frame(
        chrom = "chr1",
        start = c(500, 800, 1200, 1800, 2200, 2800, 3200),
        end = c(600, 900, 1300, 1900, 2300, 2900, 3300),
        feature_id = paste0("F", 1:7)
    )

    # Test basic use_intervals1_strand functionality
    result <- gintervals.neighbors(tss, features,
        maxneighbors = 3,
        mindist = -1000, maxdist = 1000,
        use_intervals1_strand = TRUE
    )

    # Should have results
    expect_true(nrow(result) > 0)

    # Check that distance signs are correctly calculated for different strands
    gene_a_results <- result[result$gene == "GeneA", ] # + strand gene at 1000
    gene_b_results <- result[result$gene == "GeneB", ] # - strand gene at 2000
    gene_c_results <- result[result$gene == "GeneC", ] # + strand gene at 3000

    # For + strand genes (A and C), upstream features should have negative distances
    # F1 (500-600) and F2 (800-900) are upstream of GeneA (1000)
    upstream_a <- gene_a_results[gene_a_results$feature_id %in% c("F1", "F2"), ]
    expect_true(all(upstream_a$dist < 0), "Upstream features should have negative distance for + strand genes")

    # F3 (1200-1300) is downstream of GeneA (1000)
    downstream_a <- gene_a_results[gene_a_results$feature_id == "F3", ]
    expect_true(all(downstream_a$dist > 0), "Downstream features should have positive distance for + strand genes")

    # For - strand gene B, upstream features should have positive distances
    # F4 (1800-1900) is upstream of GeneB (2000) on - strand (lower coordinates)
    upstream_b <- gene_b_results[gene_b_results$feature_id == "F4", ]
    if (nrow(upstream_b) > 0) {
        expect_true(all(upstream_b$dist > 0), "Upstream features should have positive distance for - strand genes")
    }

    # F5 (2200-2300) is downstream of GeneB (2000) on - strand (higher coordinates)
    downstream_b <- gene_b_results[gene_b_results$feature_id == "F5", ]
    if (nrow(downstream_b) > 0) {
        expect_true(all(downstream_b$dist < 0), "Downstream features should have negative distance for - strand genes")
    }
})

test_that("convenience functions work correctly", {
    # Create test data
    tss <- data.frame(
        chrom = c("chr1", "chr1"),
        start = c(1000, 2000),
        end = c(1001, 2001),
        strand = c(1, -1), # +, -
        gene = c("GeneA", "GeneB")
    )

    features <- data.frame(
        chrom = "chr1",
        start = c(500, 800, 1200, 1800, 2200, 2800),
        end = c(600, 900, 1300, 1900, 2300, 2900),
        feature_id = paste0("F", 1:6)
    )

    # Test upstream function
    upstream <- gintervals.neighbors.upstream(tss, features, maxneighbors = 2, maxdist = 1000)
    expect_true(nrow(upstream) > 0)
    expect_true(all(upstream$dist <= 0), "Upstream function should only return non-positive distances")

    # Test downstream function
    downstream <- gintervals.neighbors.downstream(tss, features, maxneighbors = 2, maxdist = 1000)
    expect_true(nrow(downstream) > 0)
    expect_true(all(downstream$dist >= 0), "Downstream function should only return non-negative distances")

    # Test directional function
    both <- gintervals.neighbors.directional(tss, features,
        maxneighbors_upstream = 1,
        maxneighbors_downstream = 1,
        maxdist = 1000
    )
    expect_true(is.list(both))
    expect_true("upstream" %in% names(both))
    expect_true("downstream" %in% names(both))
    expect_true(all(both$upstream$dist <= 0))
    expect_true(all(both$downstream$dist >= 0))
})

test_that("parameter validation works", {
    tss <- data.frame(
        chrom = "chr1", start = 1000, end = 1001, strand = 1
    )
    features <- data.frame(
        chrom = "chr1", start = 500, end = 600
    )

    # Test missing strand column - should work now (defaults to strand=1)
    tss_no_strand <- data.frame(
        chrom = "chr1", start = 1000, end = 1001
    )

    # This should now work without error, treating all intervals as strand=1
    expect_no_error(
        gintervals.neighbors(tss_no_strand, features, use_intervals1_strand = TRUE)
    )

    # Test convenience functions without strand column - should work and default to strand=1
    expect_no_error(
        gintervals.neighbors.upstream(tss_no_strand, features)
    )

    expect_no_error(
        gintervals.neighbors.downstream(tss_no_strand, features)
    )

    # Test invalid strand values (should now error)
    tss_invalid_strand <- data.frame(
        chrom = "chr1", start = 1000, end = 1001, strand = 2 # Invalid strand value
    )

    # Should error on invalid strand values
    expect_error(
        gintervals.neighbors(tss_invalid_strand, features, use_intervals1_strand = TRUE),
        "Invalid strand value.*Strand must be -1 or 1"
    )

    # Test invalid strand values in multiple rows
    tss_multiple_invalid <- data.frame(
        chrom = c("chr1", "chr1", "chr1"),
        start = c(1000, 2000, 3000),
        end = c(1001, 2001, 3001),
        strand = c(1, 0, -1) # 0 is invalid
    )

    # Should error on invalid strand values (0 at row 2)
    expect_error(
        gintervals.neighbors(tss_multiple_invalid, features, use_intervals1_strand = TRUE),
        "Invalid strand value.*Strand must be -1 or 1"
    )
})

test_that("backward compatibility maintained", {
    # Test that original behavior is unchanged when use_intervals1_strand = FALSE
    intervs1 <- data.frame(
        chrom = c("chr1", "chr1"),
        start = c(1000, 2000),
        end = c(1001, 2001),
        strand = c(1, -1) # This should be ignored with default parameters
    )

    intervs2 <- data.frame(
        chrom = "chr1",
        start = c(500, 1500),
        end = c(600, 1600),
        strand = c(1, -1) # This should be used for directionality
    )

    # Default behavior (should use intervals2 strand)
    result_default <- gintervals.neighbors(intervs1, intervs2,
        maxneighbors = 10,
        mindist = -1000, maxdist = 1000,
        warn.ignored.strand = FALSE
    )

    # Explicit FALSE (should be identical)
    result_false <- gintervals.neighbors(intervs1, intervs2,
        maxneighbors = 10,
        mindist = -1000, maxdist = 1000,
        use_intervals1_strand = FALSE,
        warn.ignored.strand = FALSE
    )

    # Results should be identical
    expect_equal(result_default, result_false)

    # New behavior should be different
    result_new <- gintervals.neighbors(intervs1, intervs2,
        maxneighbors = 10,
        mindist = -1000, maxdist = 1000,
        use_intervals1_strand = TRUE,
        warn.ignored.strand = FALSE
    )

    # Should have different distance values
    expect_false(all(result_default$dist == result_new$dist))
})

test_that("intervals without strand column work with use_intervals1_strand = TRUE", {
    # Test that intervals without a strand column work correctly
    # and are treated as if all have strand=1
    tss_no_strand <- data.frame(
        chrom = c("chr1", "chr1"),
        start = c(1000, 2000),
        end = c(1001, 2001)
        # No strand column - should default to strand=1
    )

    features <- data.frame(
        chrom = "chr1",
        start = c(500, 1500, 2500),
        end = c(600, 1600, 2600)
    )

    # This should work without error
    result <- gintervals.neighbors(tss_no_strand, features,
        maxneighbors = 1,
        mindist = -1000, maxdist = 1000,
        use_intervals1_strand = TRUE
    )

    # Should return results
    expect_true(nrow(result) > 0)

    # Since all intervals are treated as strand=1, distances should be consistent
    # with positive strand interpretation
    expect_true(all(!is.na(result$dist)))
})

test_that("strand column detection works correctly", {
    # This test verifies the fix for the hardcoded 'true' parameter bug
    # where get_strand_value was always called with has_strand_col=true

    # Test case 1: intervals with strand column
    tss_with_strand <- data.frame(
        chrom = c("chr1", "chr1"),
        start = c(1000, 2000),
        end = c(1001, 2001),
        strand = c(1, -1)
    )

    # Test case 2: intervals without strand column
    tss_no_strand <- data.frame(
        chrom = c("chr1", "chr1"),
        start = c(1000, 2000),
        end = c(1001, 2001)
        # No strand column
    )

    features <- data.frame(
        chrom = "chr1",
        start = c(500, 1500, 2500),
        end = c(600, 1600, 2600)
    )

    # Both should work without error
    result_with_strand <- gintervals.neighbors(tss_with_strand, features,
        maxneighbors = 1,
        mindist = -1000, maxdist = 1000,
        use_intervals1_strand = TRUE
    )

    result_no_strand <- gintervals.neighbors(tss_no_strand, features,
        maxneighbors = 1,
        mindist = -1000, maxdist = 1000,
        use_intervals1_strand = TRUE
    )

    # Both should return results
    expect_true(nrow(result_with_strand) > 0)
    expect_true(nrow(result_no_strand) > 0)

    # The distances should be different because one uses actual strand values
    # and the other defaults to all strand=1
    # For the first TSS at position 1000:
    # - With strand=1: upstream feature (500-600) should have negative distance
    # - With strand=-1: upstream feature (500-600) should have positive distance
    # - Without strand (defaults to 1): upstream feature should have negative distance

    # Find results for first TSS (position 1000)
    with_strand_tss1 <- result_with_strand[result_with_strand$start == 1000 &
        result_with_strand$start1 == 500, ]
    no_strand_tss1 <- result_no_strand[result_no_strand$start == 1000 &
        result_no_strand$start1 == 500, ]

    if (nrow(with_strand_tss1) > 0 && nrow(no_strand_tss1) > 0) {
        # Both should have negative distance (upstream for + strand)
        expect_true(with_strand_tss1$dist[1] < 0)
        expect_true(no_strand_tss1$dist[1] < 0)
        # Distances should be the same since both are + strand
        expect_equal(with_strand_tss1$dist[1], no_strand_tss1$dist[1])
    }

    # Find results for second TSS (position 2000)
    # With strand: strand=-1, so upstream (1500-1600) should have positive distance
    # Without strand: defaults to strand=1, so upstream should have negative distance
    with_strand_tss2 <- result_with_strand[result_with_strand$start == 2000 &
        result_with_strand$start1 == 1500, ]
    no_strand_tss2 <- result_no_strand[result_no_strand$start == 2000 &
        result_no_strand$start1 == 1500, ]

    if (nrow(with_strand_tss2) > 0 && nrow(no_strand_tss2) > 0) {
        # With strand=-1: upstream should have positive distance
        expect_true(with_strand_tss2$dist[1] > 0)
        # Without strand (default +1): upstream should have negative distance
        expect_true(no_strand_tss2$dist[1] < 0)
        # The signs should be different - this proves the fix works!
        expect_true(sign(with_strand_tss2$dist[1]) != sign(no_strand_tss2$dist[1]))
    }
})
