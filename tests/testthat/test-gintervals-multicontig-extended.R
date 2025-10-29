load_test_db()
skip_if(getOption("gmulticontig.indexed_format", FALSE), "Indexed format not enabled, set gmulticontig.indexed_format = TRUE to run this test")
# # Extended tests for multi-contig intervals implementation
# # These tests cover edge cases, error conditions, and stress scenarios

# # ============================================================================
# # Basic 1D interval conversion tests
# # ============================================================================

# test_that("gintervals.convert_to_indexed works for 1D intervals", {
#     withr::defer(gintervals.rm("test.converted_1d", force = TRUE))

#     # Create test interval set
#     intervs <- gintervals(
#         chrom = c(1, 1, 2, 2, 3),
#         start = c(0, 10000, 0, 10000, 0),
#         end = c(5000, 15000, 5000, 15000, 5000)
#     )
#     gintervals.save("test.converted_1d", intervs)

#     # Check format before conversion
#     info_before <- gintervals.info("test.converted_1d")
#     expect_true("format" %in% names(info_before))
#     expect_equal(info_before$format, "per-chromosome")

#     # Conversion
#     gintervals.convert_to_indexed("test.converted_1d")

#     # Check format after conversion
#     info_after <- gintervals.info("test.converted_1d")
#     expect_equal(info_after$format, "indexed")

#     # Verify data integrity
#     loaded <- gintervals.load("test.converted_1d")
#     expect_equal(nrow(loaded), nrow(intervs))
#     expect_equal(loaded$chrom, intervs$chrom)
#     expect_equal(loaded$start, intervs$start)
#     expect_equal(loaded$end, intervs$end)
# })

# test_that("gintervals.convert_to_indexed handles empty chromosomes", {
#     withr::defer(gintervals.rm("test.empty_chr_1d", force = TRUE))

#     # Create intervals only on some chromosomes
#     intervs <- gintervals(
#         chrom = c(1, 1, 3), # Skip chromosome 2
#         start = c(0, 10000, 0),
#         end = c(5000, 15000, 5000)
#     )
#     gintervals.save("test.empty_chr_1d", intervs)

#     gintervals.convert_to_indexed("test.empty_chr_1d")

#     # Should still load correctly
#     loaded <- gintervals.load("test.empty_chr_1d")
#     expect_equal(nrow(loaded), nrow(intervs))
# })

# test_that("gintervals.convert_to_indexed preserves interval metadata", {
#     withr::defer(gintervals.rm("test.meta_1d", force = TRUE))

#     # Create intervals with metadata columns
#     intervs <- gintervals(
#         chrom = c(1, 2, 3),
#         start = c(0, 0, 0),
#         end = c(1000, 1000, 1000)
#     )
#     intervs$score <- c(10, 20, 30)
#     intervs$name <- c("a", "b", "c")

#     gintervals.save("test.meta_1d", intervs)

#     meta_before <- gintervals.load("test.meta_1d")

#     gintervals.convert_to_indexed("test.meta_1d")

#     meta_after <- gintervals.load("test.meta_1d")

#     expect_equal(meta_before$score, meta_after$score)
#     expect_equal(meta_before$name, meta_after$name)
# })

# # ============================================================================
# # 2D interval tests
# # ============================================================================

# test_that("gintervals.2d.convert_to_indexed works for 2D intervals", {
#     withr::defer(gintervals.rm("test.converted_2d", force = TRUE))

#     # Create test 2D interval set
#     intervs <- gintervals.2d(
#         chroms1 = c(1, 1, 2),
#         starts1 = c(0, 10000, 0),
#         ends1 = c(5000, 15000, 5000),
#         chroms2 = c(1, 2, 2),
#         starts2 = c(0, 0, 10000),
#         ends2 = c(5000, 5000, 15000)
#     )
#     gintervals.save("test.converted_2d", intervs)

#     # Check format before
#     info_before <- gintervals.info("test.converted_2d")
#     expect_equal(info_before$format, "per-pair")

#     # Convert
#     gintervals.2d.convert_to_indexed("test.converted_2d")

#     # Check format after
#     info_after <- gintervals.info("test.converted_2d")
#     expect_equal(info_after$format, "indexed")

#     # Verify data
#     loaded <- gintervals.load("test.converted_2d")
#     expect_equal(nrow(loaded), nrow(intervs))
#     expect_equal(loaded$chrom1, intervs$chrom1)
#     expect_equal(loaded$start1, intervs$start1)
#     expect_equal(loaded$end1, intervs$end1)
#     expect_equal(loaded$chrom2, intervs$chrom2)
#     expect_equal(loaded$start2, intervs$start2)
#     expect_equal(loaded$end2, intervs$end2)
# })

# test_that("gintervals.2d.convert_to_indexed handles sparse chromosome pairs", {
#     withr::defer(gintervals.rm("test.sparse_pairs", force = TRUE))

#     # Create 2D intervals with only some chromosome pairs
#     intervs <- gintervals.2d(
#         chroms1 = c(1, 3),
#         starts1 = c(0, 0),
#         ends1 = c(1000, 1000),
#         chroms2 = c(2, 3),
#         starts2 = c(0, 0),
#         ends2 = c(1000, 1000)
#     )
#     gintervals.save("test.sparse_pairs", intervs)

#     gintervals.2d.convert_to_indexed("test.sparse_pairs")

#     loaded <- gintervals.load("test.sparse_pairs")
#     expect_equal(nrow(loaded), nrow(intervs))
# })

# # ============================================================================
# # Data integrity tests
# # ============================================================================

# test_that("converted 1D intervals match original exactly", {
#     withr::defer(gintervals.rm("test.integrity_1d", force = TRUE))

#     # Create large interval set
#     set.seed(42)
#     n <- 1000
#     intervs <- gintervals(
#         chrom = sample(1:3, n, replace = TRUE),
#         start = sample(0:100000, n, replace = TRUE),
#         end = sample(0:100000, n, replace = TRUE)
#     )
#     intervs$end <- pmax(intervs$start + 100, intervs$end)
#     intervs <- intervs[order(intervs$chrom, intervs$start), ]

#     gintervals.save("test.integrity_1d", intervs)

#     before <- gintervals.load("test.integrity_1d")

#     gintervals.convert_to_indexed("test.integrity_1d")

#     after <- gintervals.load("test.integrity_1d")

#     expect_equal(before, after)
# })

# test_that("converted intervals work with gintervals.* functions", {
#     withr::defer(gintervals.rm("test.functions_1d", force = TRUE))

#     intervs <- gintervals(c(1, 2, 3), 0, 10000)
#     gintervals.save("test.functions_1d", intervs)
#     gintervals.convert_to_indexed("test.functions_1d")

#     # Test various functions
#     canonic <- gintervals.canonic("test.functions_1d")
#     expect_true(nrow(canonic) > 0)

#     # Use in operations
#     result <- gscreen("test.fixedbin > 0.5", intervals = "test.functions_1d")
#     expect_true(is.data.frame(result))
# })

# test_that("converted intervals work as iterator", {
#     withr::defer(gintervals.rm("test.iterator_1d", force = TRUE))

#     intervs <- gintervals(
#         chrom = c(1, 1, 2),
#         start = c(0, 10000, 0),
#         end = c(5000, 15000, 5000)
#     )
#     gintervals.save("test.iterator_1d", intervs)
#     gintervals.convert_to_indexed("test.iterator_1d")

#     # Use as iterator
#     result <- gextract("test.fixedbin",
#         gintervals(c(1, 2)),
#         iterator = "test.iterator_1d"
#     )

#     expect_equal(nrow(result), 3)
#     expect_equal(result$start, intervs$start)
#     expect_equal(result$end, intervs$end)
# })

# # ============================================================================
# # Large dataset tests
# # ============================================================================

# test_that("conversion handles large interval sets", {
#     withr::defer(gintervals.rm("test.large_1d", force = TRUE))
#     withr::local_options(gmax.data.size = 1e9)

#     # Create large interval set
#     set.seed(123)
#     n <- 10000
#     intervs <- gintervals(
#         chrom = sample(1:3, n, replace = TRUE),
#         start = sample(0:1000000, n, replace = TRUE),
#         end = sample(0:1000000, n, replace = TRUE)
#     )
#     intervs$end <- pmax(intervs$start + 100, intervs$end)

#     gintervals.save("test.large_1d", intervs)
#     gintervals.convert_to_indexed("test.large_1d")

#     loaded <- gintervals.load("test.large_1d")
#     expect_equal(nrow(loaded), n)
# })

# test_that("converted intervals work with large extraction", {
#     withr::defer(gintervals.rm("test.large_extract", force = TRUE))
#     withr::local_options(gmax.data.size = 1e9)

#     # Create many intervals
#     intervs <- gintervals(
#         chrom = rep(1:3, each = 100),
#         start = rep(seq(0, 990000, by = 10000), 3),
#         end = rep(seq(1000, 1000000, by = 10000), 3)
#     )

#     gintervals.save("test.large_extract", intervs)
#     gintervals.convert_to_indexed("test.large_extract")

#     # Extract large amount of data
#     result <- gextract("test.fixedbin",
#         iterator = "test.large_extract"
#     )

#     expect_equal(nrow(result), nrow(intervs))
# })

# # ============================================================================
# # Edge case tests
# # ============================================================================

# test_that("converted intervals handle overlapping intervals", {
#     withr::defer(gintervals.rm("test.overlapping", force = TRUE))

#     # Create overlapping intervals
#     intervs <- gintervals(
#         chrom = c(1, 1, 1),
#         start = c(0, 500, 1000),
#         end = c(1000, 1500, 2000)
#     )

#     gintervals.save("test.overlapping", intervs)
#     gintervals.convert_to_indexed("test.overlapping")

#     loaded <- gintervals.load("test.overlapping")
#     expect_equal(nrow(loaded), 3)
# })

# test_that("converted intervals handle adjacent intervals", {
#     withr::defer(gintervals.rm("test.adjacent", force = TRUE))

#     # Create adjacent intervals
#     intervs <- gintervals(
#         chrom = c(1, 1, 1),
#         start = c(0, 1000, 2000),
#         end = c(1000, 2000, 3000)
#     )

#     gintervals.save("test.adjacent", intervs)
#     gintervals.convert_to_indexed("test.adjacent")

#     loaded <- gintervals.load("test.adjacent")
#     expect_equal(nrow(loaded), 3)
#     expect_equal(loaded$start, intervs$start)
#     expect_equal(loaded$end, intervs$end)
# })

# test_that("converted intervals handle single-base intervals", {
#     withr::defer(gintervals.rm("test.single_base", force = TRUE))

#     # Create single-base intervals
#     intervs <- gintervals(
#         chrom = c(1, 2, 3),
#         start = c(100, 200, 300),
#         end = c(101, 201, 301)
#     )

#     gintervals.save("test.single_base", intervs)
#     gintervals.convert_to_indexed("test.single_base")

#     loaded <- gintervals.load("test.single_base")
#     expect_equal(loaded$start, intervs$start)
#     expect_equal(loaded$end, intervs$end)
# })

# # ============================================================================
# # Modification after conversion tests
# # ============================================================================

# test_that("can modify intervals after conversion", {
#     withr::defer(gintervals.rm("test.modify_after", force = TRUE))

#     intervs <- gintervals(c(1, 2), 0, 10000)
#     gintervals.save("test.modify_after", intervs)
#     gintervals.convert_to_indexed("test.modify_after")

#     # Add more intervals
#     new_intervs <- gintervals(3, 0, 10000)
#     combined <- rbind(gintervals.load("test.modify_after"), new_intervs)
#     gintervals.save("test.modify_after", combined, overwrite = TRUE)

#     loaded <- gintervals.load("test.modify_after")
#     expect_equal(nrow(loaded), 3)
# })

# # ============================================================================
# # Multiple interval set operations
# # ============================================================================

# test_that("converted intervals work with set operations", {
#     withr::defer(gintervals.rm("test.set_ops1", force = TRUE))
#     withr::defer(gintervals.rm("test.set_ops2", force = TRUE))

#     int1 <- gintervals(c(1, 2), 0, 10000)
#     int2 <- gintervals(c(2, 3), 0, 10000)

#     gintervals.save("test.set_ops1", int1)
#     gintervals.save("test.set_ops2", int2)

#     gintervals.convert_to_indexed("test.set_ops1")
#     gintervals.convert_to_indexed("test.set_ops2")

#     # Union
#     union_result <- gunion("test.set_ops1", "test.set_ops2")
#     expect_true(nrow(union_result) > 0)

#     # Intersection
#     intersect_result <- gintersect("test.set_ops1", "test.set_ops2")
#     expect_equal(nrow(intersect_result), 1) # Only chr2 overlaps

#     # Difference
#     diff_result <- gdiff("test.set_ops1", "test.set_ops2")
#     expect_true(nrow(diff_result) >= 0)
# })

# # ============================================================================
# # Attribute preservation tests
# # ============================================================================

# test_that("converted intervals preserve all column types", {
#     withr::defer(gintervals.rm("test.col_types", force = TRUE))

#     intervs <- gintervals(c(1, 2, 3), 0, 1000)
#     intervs$int_col <- c(1L, 2L, 3L)
#     intervs$num_col <- c(1.1, 2.2, 3.3)
#     intervs$char_col <- c("a", "b", "c")
#     intervs$logical_col <- c(TRUE, FALSE, TRUE)

#     gintervals.save("test.col_types", intervs)

#     before <- gintervals.load("test.col_types")

#     gintervals.convert_to_indexed("test.col_types")

#     after <- gintervals.load("test.col_types")

#     expect_equal(typeof(after$int_col), typeof(before$int_col))
#     expect_equal(typeof(after$num_col), typeof(before$num_col))
#     expect_equal(typeof(after$char_col), typeof(before$char_col))
#     expect_equal(typeof(after$logical_col), typeof(before$logical_col))

#     expect_equal(after$int_col, before$int_col)
#     expect_equal(after$num_col, before$num_col)
#     expect_equal(after$char_col, before$char_col)
#     expect_equal(after$logical_col, before$logical_col)
# })

# # ============================================================================
# # Complex workflow tests
# # ============================================================================

# test_that("converted intervals in complex gextract workflow", {
#     withr::defer(gintervals.rm("test.workflow", force = TRUE))

#     # Create intervals
#     intervs <- gintervals(c(1, 2, 3), 0, 50000)
#     gintervals.save("test.workflow", intervs)
#     gintervals.convert_to_indexed("test.workflow")

#     # Use in complex extraction
#     result <- gextract(
#         "test.fixedbin * 2 + test.sparse",
#         iterator = "test.workflow",
#         colnames = c("dense", "sparse")
#     )

#     expect_equal(nrow(result), 3)
#     expect_true("dense" %in% names(result))
#     expect_true("sparse" %in% names(result))
# })

# test_that("converted intervals work with gscreen conditions", {
#     withr::defer(gintervals.rm("test.screen_intervals", force = TRUE))

#     intervs <- gintervals(c(1, 2, 3), 0, 100000)
#     gintervals.save("test.screen_intervals", intervs)
#     gintervals.convert_to_indexed("test.screen_intervals")

#     # Screen within these intervals
#     result <- gscreen("test.fixedbin > 0.5",
#         intervals = "test.screen_intervals"
#     )

#     expect_true(is.data.frame(result))
#     # All results should be within our interval set
#     for (i in 1:nrow(result)) {
#         expect_true(result$chrom[i] %in% c(1, 2, 3))
#         expect_true(result$start[i] >= 0)
#         expect_true(result$end[i] <= 100000)
#     }
# })

# # ============================================================================
# # Performance and stress tests
# # ============================================================================

# test_that("converted intervals handle rapid sequential access", {
#     withr::defer(gintervals.rm("test.rapid_access", force = TRUE))

#     intervs <- gintervals(rep(1:3, each = 100), 0, 1000)
#     gintervals.save("test.rapid_access", intervs)
#     gintervals.convert_to_indexed("test.rapid_access")

#     # Rapid sequential loads
#     for (i in 1:10) {
#         loaded <- gintervals.load("test.rapid_access")
#         expect_equal(nrow(loaded), 300)
#     }
# })

# test_that("converted intervals work with nested operations", {
#     withr::defer(gintervals.rm("test.nested1", force = TRUE))
#     withr::defer(gintervals.rm("test.nested2", force = TRUE))

#     int1 <- gintervals(c(1, 2), 0, 10000)
#     int2 <- gintervals(c(2, 3), 5000, 15000)

#     gintervals.save("test.nested1", int1)
#     gintervals.save("test.nested2", int2)

#     gintervals.convert_to_indexed("test.nested1")
#     gintervals.convert_to_indexed("test.nested2")

#     # Nested set operations
#     union_set <- gunion("test.nested1", "test.nested2")
#     gintervals.save("test.nested_result", union_set)

#     # Use result in extraction
#     result <- gextract("test.fixedbin", iterator = union_set)
#     expect_true(nrow(result) > 0)

#     gintervals.rm("test.nested_result", force = TRUE)
# })

# # ============================================================================
# # Backward compatibility tests
# # ============================================================================

# test_that("can revert to per-chromosome format if needed", {
#     withr::defer(gintervals.rm("test.revert", force = TRUE))

#     intervs <- gintervals(c(1, 2, 3), 0, 10000)
#     gintervals.save("test.revert", intervs)

#     before_convert <- gintervals.load("test.revert")

#     gintervals.convert_to_indexed("test.revert", remove_old = FALSE)

#     after_convert <- gintervals.load("test.revert")

#     # Data should be identical
#     expect_equal(before_convert, after_convert)
# })

# # ============================================================================
# # Error handling tests
# # ============================================================================

# test_that("gintervals.convert_to_indexed fails gracefully for non-existent intervals", {
#     expect_error(
#         gintervals.convert_to_indexed("nonexistent_intervals_xyz"),
#         "does not exist|cannot find"
#     )
# })

# test_that("converted intervals handle missing chromosome gracefully", {
#     withr::defer(gintervals.rm("test.missing_chr", force = TRUE))

#     # Create intervals on subset of chromosomes
#     intervs <- gintervals(c(1, 3), 0, 10000) # Skip chr 2
#     gintervals.save("test.missing_chr", intervs)
#     gintervals.convert_to_indexed("test.missing_chr")

#     loaded <- gintervals.load("test.missing_chr")

#     # Should only have chr 1 and 3
#     expect_equal(sort(unique(loaded$chrom)), c(1, 3))
#     expect_false(2 %in% loaded$chrom)
# })
