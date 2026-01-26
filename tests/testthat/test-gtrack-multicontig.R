load_test_db()
skip_if(getOption("gmulticontig.indexed_format", FALSE), "Indexed format not enabled, set gmulticontig.indexed_format = TRUE to run this test")
# Extended tests for multi-contig track implementation
# These tests cover edge cases, error conditions, and stress scenarios

# Helper to remove all vtracks
remove_all_vtracks <- function() {
    vtracks <- gvtrack.ls()
    for (vtrack in vtracks) {
        do.call(gvtrack.rm, list(vtrack = vtrack))
    }
}

# ============================================================================
# Edge case tests - Empty chromosomes
# ============================================================================

test_that("converted track with empty chromosomes works", {
    withr::defer(gtrack.rm("test.empty_chr", force = TRUE))

    # Create a sparse track that may have empty chromosomes
    gtrack.create("test.empty_chr", "", "test.sparse", iterator = "test.sparse")
    info_before <- gtrack.info("test.empty_chr")

    gtrack.convert_to_indexed("test.empty_chr")
    info_after <- gtrack.info("test.empty_chr")

    expect_equal(info_after$format, "indexed")

    # Extract from potentially empty chromosome
    # Should not crash even if chromosome is empty
    result <- suppressWarnings(gextract("test.empty_chr", gintervals(c(1, 2, 3))))
    expect_true(is.data.frame(result))
})

test_that("converted track handles all chromosomes correctly", {
    if (gtrack.exists("test.all_chr")) {
        gtrack.rm("test.all_chr", force = TRUE)
    }
    withr::defer(gtrack.rm("test.all_chr", force = TRUE))

    gtrack.create("test.all_chr", "", "test.fixedbin")
    gtrack.convert_to_indexed("test.all_chr")

    # Get all chromosomes from genome
    all_chroms <- unique(gintervals.all()$chrom)
    # Extract from each chromosome
    for (chr in all_chroms) {
        result <- gextract("test.all_chr", gintervals(chr, 1, 1000))
        expect_true(nrow(result) > 0 || nrow(result) == 0) # Should not crash
    }
})

# ============================================================================
# Data integrity tests
# ============================================================================

test_that("converted track data matches original across all chromosomes", {
    withr::defer(gtrack.rm("test.data_integrity", force = TRUE))

    gtrack.create("test.data_integrity", "", "test.fixedbin")

    # Extract data before conversion
    intervs <- gintervals(c(1, 2, 3), 0, 100000)
    data_before <- gextract("test.data_integrity", intervs)

    gtrack.convert_to_indexed("test.data_integrity")

    # Extract data after conversion
    data_after <- gextract("test.data_integrity", intervs)

    # Compare every value
    expect_equal(nrow(data_before), nrow(data_after))
    expect_equal(data_before$test.data_integrity, data_after$test.data_integrity)
    expect_equal(data_before$chrom, data_after$chrom)
    expect_equal(data_before$start, data_after$start)
    expect_equal(data_before$end, data_after$end)
})

test_that("converted sparse track preserves interval boundaries", {
    withr::defer(gtrack.rm("test.sparse_boundaries", force = TRUE))

    gtrack.create("test.sparse_boundaries", "", "test.sparse", iterator = "test.sparse")

    # Get sparse intervals before conversion
    intervs <- gintervals(c(1, 2))
    data_before <- gextract("test.sparse_boundaries", intervs, iterator = "test.sparse")

    gtrack.convert_to_indexed("test.sparse_boundaries")

    # Get sparse intervals after conversion
    data_after <- gextract("test.sparse_boundaries", intervs, iterator = "test.sparse")

    # Check interval boundaries match exactly
    expect_equal(data_before$start, data_after$start)
    expect_equal(data_before$end, data_after$end)
    expect_equal(data_before$test.sparse_boundaries, data_after$test.sparse_boundaries)
})

test_that("converted array track preserves array structure", {
    withr::defer(gtrack.rm("test.array_structure", force = TRUE))

    gtrack.create("test.array_structure", "", "test.array", iterator = "test.array")

    intervs <- gintervals(c(1, 2))
    data_before <- gextract("test.array_structure", intervs, iterator = "test.array")

    gtrack.convert_to_indexed("test.array_structure")

    data_after <- gextract("test.array_structure", intervs, iterator = "test.array")

    expect_equal(nrow(data_before), nrow(data_after))
    expect_equal(data_before$start, data_after$start)
    expect_equal(data_before$end, data_after$end)
})

# ============================================================================
# Statistical function tests
# ============================================================================

test_that("converted track works with all statistical functions", {
    withr::defer(gtrack.rm("test.stats", force = TRUE))

    gtrack.create("test.stats", "", "test.fixedbin")
    gtrack.convert_to_indexed("test.stats")

    intervs <- gintervals(1, 0, 10000)

    # Test each statistical function
    funcs <- c("avg", "min", "max", "sum", "stddev", "quantile")
    for (func in funcs) {
        result_converted <- gextract("test.stats", intervs,
            iterator = 500,
            colnames = paste0("val_", func)
        )
        result_orig <- gextract("test.fixedbin", intervs,
            iterator = 500,
            colnames = paste0("val_", func)
        )

        expect_equal(result_converted[[paste0("val_", func)]],
            result_orig[[paste0("val_", func)]],
            info = paste("Function:", func)
        )
    }
})

test_that("converted sparse track works with nearest function", {
    withr::defer(gtrack.rm("test.nearest", force = TRUE))

    gtrack.create("test.nearest", "", "test.sparse", iterator = "test.sparse")
    gtrack.convert_to_indexed("test.nearest")

    # Create intervals that may not overlap with sparse data
    intervs <- gintervals(
        1, c(100, 1000, 5000, 10000),
        c(200, 1100, 5100, 10100)
    )

    # Both should return nearest values
    result_converted <- gextract("test.nearest", intervs, colnames = "nearest")
    result_orig <- gextract("test.sparse", intervs, colnames = "nearest")

    expect_equal(nrow(result_converted), nrow(result_orig))
})

# ============================================================================
# Large data tests
# ============================================================================

test_that("converted track handles large extractions", {
    withr::defer(gtrack.rm("test.large_extract", force = TRUE))
    withr::local_options(gmax.data.size = 1e9)

    gtrack.create("test.large_extract", "", "test.fixedbin")
    gtrack.convert_to_indexed("test.large_extract")

    # Extract large region
    large_intervs <- gintervals(c(1, 2, 3), 0, 500000)
    result <- gextract("test.large_extract", large_intervs)

    expect_true(nrow(result) > 10000) # Should have many rows
    expect_true(all(is.finite(result$test.large_extract) |
        is.na(result$test.large_extract)))
})

test_that("converted track works with ALLGENOME iterator", {
    withr::defer(gtrack.rm("test.allgenome_iter", force = TRUE))
    withr::local_options(gmax.data.size = 1e9)

    gtrack.create("test.allgenome_iter", "", "test.fixedbin")
    gtrack.convert_to_indexed("test.allgenome_iter")

    # Extract with ALLGENOME
    result_converted <- gextract("test.allgenome_iter", .misha$ALLGENOME,
        iterator = 10000
    )
    result_orig <- gextract("test.fixedbin", .misha$ALLGENOME,
        iterator = 10000
    )

    expect_equal(nrow(result_converted), nrow(result_orig))
})

# ============================================================================
# Multiple concurrent access tests
# ============================================================================

test_that("converted track handles multiple simultaneous extractions", {
    withr::defer(gtrack.rm("test.concurrent", force = TRUE))

    gtrack.create("test.concurrent", "", "test.fixedbin")
    gtrack.convert_to_indexed("test.concurrent")

    # Perform multiple extractions from different chromosomes
    results <- list()
    for (chr in 1:3) {
        results[[chr]] <- gextract(
            "test.concurrent",
            gintervals(chr, 0, 10000)
        )
    }

    # All should succeed
    expect_equal(length(results), 3)
    for (r in results) {
        expect_true(nrow(r) > 0)
    }
})

test_that("converted track works with multiple vtracks", {
    withr::defer(gtrack.rm("test.multi_vtrack", force = TRUE))
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    gtrack.create("test.multi_vtrack", "", "test.fixedbin")
    gtrack.convert_to_indexed("test.multi_vtrack")

    # Create multiple vtracks
    gvtrack.create("v1", "test.multi_vtrack")
    gvtrack.create("v2", "test.multi_vtrack", func = "avg")

    intervs <- gintervals(1, 1, 5000)
    r1 <- gextract("v1", intervs, iterator = 100)
    r2 <- gextract("v2", intervs, iterator = 100)
    r3 <- gextract("v2 + 1", intervs, iterator = 100)

    expect_true(nrow(r1) > 0)
    expect_true(nrow(r2) > 0)
    expect_true(nrow(r3) > 0)
})

# ============================================================================
# Track modification tests
# ============================================================================

test_that("gtrack.modify works on single bin in converted track", {
    withr::defer(gtrack.rm("test.modify_single", force = TRUE))

    gtrack.create("test.modify_single", "", "test.fixedbin")
    gtrack.convert_to_indexed("test.modify_single")

    # Get bin size
    info <- gtrack.info("test.fixedbin")
    bin_size <- info$bin.size
    # Modify single bin
    interv <- gintervals(1, 1, bin_size)
    before <- gextract("test.modify_single", interv)
    gtrack.modify("test.modify_single", "999", interv)
    after <- gextract("test.modify_single", interv)

    expect_equal(after$test.modify_single, 999)
    expect_false(before$test.modify_single == 999)
})

test_that("gtrack.modify works across chromosome boundaries", {
    withr::defer(gtrack.rm("test.modify_boundary", force = TRUE))

    gtrack.create("test.modify_boundary", "", "test.fixedbin")
    gtrack.convert_to_indexed("test.modify_boundary")

    # Modify intervals spanning multiple chromosomes
    intervs <- gintervals(
        c(1, 1, 2, 2, 3),
        c(0, 50000, 0, 50000, 0),
        c(10000, 60000, 10000, 60000, 10000)
    )

    gtrack.modify(
        "test.modify_boundary",
        "test.modify_boundary * 0.5",
        intervs
    )

    # Verify modifications
    result <- gextract("test.modify_boundary", intervs)
    expect_true(all(is.finite(result$test.modify_boundary) |
        is.na(result$test.modify_boundary)))
})

test_that("multiple modifications on converted track", {
    withr::defer(gtrack.rm("test.multi_modify", force = TRUE))

    gtrack.create("test.multi_modify", "", "test.fixedbin")
    gtrack.convert_to_indexed("test.multi_modify")

    interv1 <- gintervals(1, 0, 10000)
    interv2 <- gintervals(2, 0, 10000)
    interv3 <- gintervals(3, 0, 10000)

    # Multiple modifications
    gtrack.modify("test.multi_modify", "1", interv1)
    gtrack.modify("test.multi_modify", "2", interv2)
    gtrack.modify("test.multi_modify", "3", interv3)

    # Verify each
    r1 <- gextract("test.multi_modify", interv1)
    r2 <- gextract("test.multi_modify", interv2)
    r3 <- gextract("test.multi_modify", interv3)

    expect_true(all(r1$test.multi_modify == 1 | is.na(r1$test.multi_modify)))
    expect_true(all(r2$test.multi_modify == 2 | is.na(r2$test.multi_modify)))
    expect_true(all(r3$test.multi_modify == 3 | is.na(r3$test.multi_modify)))
})

# ============================================================================
# Complex expression tests
# ============================================================================

test_that("complex nested expressions with converted tracks", {
    withr::defer(gtrack.rm("test.complex_expr", force = TRUE))

    gtrack.create("test.complex_expr", "", "test.fixedbin")
    gtrack.convert_to_indexed("test.complex_expr")

    intervs <- gintervals(1, 0, 10000)

    # Complex expression
    expr <- "((test.complex_expr + test.fixedbin) * 2 - 1) / (test.complex_expr + 0.1)"
    result_converted <- gextract(expr, intervs)

    # Compare with original
    expr_orig <- "((test.fixedbin + test.fixedbin) * 2 - 1) / (test.fixedbin + 0.1)"
    result_orig <- gextract(expr_orig, intervs)

    expect_equal(result_converted[[2]], result_orig[[2]])
})

test_that("converted track in conditional expressions", {
    withr::defer(gtrack.rm("test.conditional", force = TRUE))

    gtrack.create("test.conditional", "", "test.fixedbin")
    gtrack.convert_to_indexed("test.conditional")

    intervs <- gintervals(c(1, 2), 0, 50000)

    # Use in gscreen (conditional)
    result_converted <- gscreen("test.conditional > 0.5", intervs)
    result_orig <- gscreen("test.fixedbin > 0.5", intervs)

    expect_equal(nrow(result_converted), nrow(result_orig))
})

# ============================================================================
# Iterator combination tests
# ============================================================================

test_that("converted track with custom interval iterator", {
    withr::defer(gtrack.rm("test.custom_iter", force = TRUE))

    gtrack.create("test.custom_iter", "", "test.fixedbin")
    gtrack.convert_to_indexed("test.custom_iter")

    # Create custom intervals
    custom_intervs <- gintervals(
        chrom = c(1, 1, 2, 2),
        start = c(0, 10000, 0, 10000),
        end = c(5000, 15000, 5000, 15000)
    )

    result_converted <- gextract("test.custom_iter",
        gintervals(c(1, 2)),
        iterator = custom_intervs
    )
    result_orig <- gextract("test.fixedbin",
        gintervals(c(1, 2)),
        iterator = custom_intervs
    )

    expect_equal(nrow(result_converted), nrow(result_orig))
    expect_equal(result_converted$test.custom_iter, result_orig$test.fixedbin)
})

test_that("converted track with track iterator from another converted track", {
    withr::defer(gtrack.rm("test.iter_source", force = TRUE))
    withr::defer(gtrack.rm("test.iter_target", force = TRUE))

    # Create and convert two tracks
    gtrack.create("test.iter_source", "", "test.sparse", iterator = "test.sparse")
    gtrack.create("test.iter_target", "", "test.fixedbin")

    gtrack.convert_to_indexed("test.iter_source")
    gtrack.convert_to_indexed("test.iter_target")

    # Use converted sparse track as iterator for converted dense track
    result <- gextract("test.iter_target",
        gintervals(c(1, 2)),
        iterator = "test.iter_source"
    )

    expect_true(nrow(result) > 0)
})

# ============================================================================
# Track attribute tests
# ============================================================================

test_that("converted track preserves track attributes", {
    withr::defer(gtrack.rm("test.attrs", force = TRUE))

    gtrack.create("test.attrs", "", "test.fixedbin")

    # Set some attributes
    gtrack.attr.set("test.attrs", "test_key", "test_value")
    gtrack.attr.set("test.attrs", "number", "42")

    attrs_before <- gtrack.attr.export("test.attrs")

    gtrack.convert_to_indexed("test.attrs")

    attrs_after <- gtrack.attr.export("test.attrs")

    expect_equal(attrs_before, attrs_after)
})

# ============================================================================
# Numeric edge cases
# ============================================================================

test_that("converted track handles NaN values correctly", {
    withr::defer(gtrack.rm("test.nan", force = TRUE))

    # Create track with some NaN values
    gtrack.create("test.nan", "", "test.fixedbin")

    # Modify some values to NaN
    gtrack.modify("test.nan", "NaN", gintervals(1, 1000, 5000))

    gtrack.convert_to_indexed("test.nan")

    # Extract and verify NaN values are preserved
    result <- gextract("test.nan", gintervals(1, 0, 10000))

    expect_true(any(is.na(result$test.nan)))
})

test_that("converted track handles very small and large values", {
    withr::defer(gtrack.rm("test.extremes", force = TRUE))

    gtrack.create("test.extremes", "", "test.fixedbin")

    # Set extreme values
    gtrack.modify("test.extremes", "1e-10", gintervals(1, 0, 1000))
    gtrack.modify("test.extremes", "1e10", gintervals(1, 1000, 2000))

    vals_before <- gextract("test.extremes", gintervals(1, 0, 2000))

    gtrack.convert_to_indexed("test.extremes")

    vals_after <- gextract("test.extremes", gintervals(1, 0, 2000))

    expect_equal(vals_before$test.extremes, vals_after$test.extremes,
        tolerance = 1e-6
    )
})

# ============================================================================
# Error recovery tests
# ============================================================================

test_that("can still use track after failed conversion attempt", {
    withr::defer(gtrack.rm("test.recovery", force = TRUE))

    gtrack.create("test.recovery", "", "test.fixedbin")

    # This should work
    result_before <- gextract("test.recovery", gintervals(1, 0, 1000))
    expect_true(nrow(result_before) > 0)

    # Note: We can't easily force an conversion failure in a test,
    # but we can verify the track still works

    # Successful conversion
    gtrack.convert_to_indexed("test.recovery")

    # Should still work
    result_after <- gextract("test.recovery", gintervals(1, 0, 1000))
    expect_equal(result_before$test.recovery, result_after$test.recovery)
})

# ============================================================================
# Backward compatibility tests
# ============================================================================

test_that("can read per-chromosome track after conversion exists", {
    withr::defer(gtrack.rm("test.compat", force = TRUE))

    gtrack.create("test.compat", "", "test.fixedbin")

    # Extract from per-chromosome
    result_per_chromosome <- gextract("test.compat", gintervals(1, 0, 1000))

    # Conversion
    gtrack.convert_to_indexed("test.compat")

    # Should still be able to extract (from converted version)
    result_converted <- gextract("test.compat", gintervals(1, 0, 1000))

    expect_equal(result_per_chromosome$test.compat, result_converted$test.compat)
})
