create_isolated_test_db()
gdir.create("temp", showWarnings = FALSE)

skip_if(getOption("gmulticontig.indexed_format", FALSE) || gdb.info(.misha$GROOT)$format == "indexed", "Indexed format enbaled, set gmulticontig.indexed_format = FALSE to run this test")
# Tests for gtrack.convert_to_indexed() and indexed track format

# Helper to remove all vtracks
remove_all_vtracks <- function() {
    vtracks <- gvtrack.ls()
    for (vtrack in vtracks) {
        do.call(gvtrack.rm, list(vtrack = vtrack))
    }
}

# ============================================================================
# Basic conversion tests
# ============================================================================

test_that("gtrack.info reports format field correctly for per-chromosome tracks", {
    # Test with existing per-chromosome track
    info <- gtrack.info("test.fixedbin")
    expect_true("format" %in% names(info))
    expect_equal(info$format, "per-chromosome")

    # Test with sparse track
    info_sparse <- gtrack.info("test.sparse")
    expect_equal(info_sparse$format, "per-chromosome")

    # Test with array track
    info_array <- gtrack.info("test.array")
    expect_equal(info_array$format, "per-chromosome")
})

test_that("gtrack.convert_to_indexed fails for non-existent track", {
    expect_error(
        gtrack.convert_to_indexed("nonexistent_track_xyz"),
        "does not exist"
    )
})

test_that("gtrack.convert_to_indexed fails for 2D tracks", {
    expect_error(
        gtrack.convert_to_indexed("test.rects"),
        "only 1D tracks.*can be converted"
    )
})

test_that("gtrack.convert_to_indexed works for dense tracks", {
    withr::defer(gtrack.rm("temp.converted_dense", force = TRUE))
    gtrack.create("temp.converted_dense", "", "test.fixedbin")
    expect_equal(gtrack.info("temp.converted_dense")$format, "per-chromosome")

    gtrack.convert_to_indexed("temp.converted_dense")
    expect_equal(gtrack.info("temp.converted_dense")$format, "indexed")

    # Verify data is intact
    r <- gextract("temp.converted_dense", gintervals(1, 0, 1000))
    r_orig <- gextract("test.fixedbin", gintervals(1, 0, 1000))
    expect_equal(r$temp.converted_dense, r_orig$test.fixedbin)
})

test_that("gtrack.convert_to_indexed works for sparse tracks", {
    withr::defer(gtrack.rm("temp.converted_sparse", force = TRUE))
    gtrack.create("temp.converted_sparse", "", "test.sparse", iterator = "test.sparse")
    expect_equal(gtrack.info("temp.converted_sparse")$format, "per-chromosome")

    gtrack.convert_to_indexed("temp.converted_sparse")
    expect_equal(gtrack.info("temp.converted_sparse")$format, "indexed")

    # Verify data is intact
    r <- gextract("temp.converted_sparse", gintervals(c(1, 2)))
    r_orig <- gextract("test.sparse", gintervals(c(1, 2)))
    expect_equal(nrow(r), nrow(r_orig))
})

test_that("gtrack.convert_to_indexed works for array tracks", {
    withr::defer(gtrack.rm("temp.converted_array", force = TRUE))
    gtrack.create("temp.converted_array", "", "test.array", iterator = "test.array")
    expect_equal(gtrack.info("temp.converted_array")$format, "per-chromosome")

    gtrack.convert_to_indexed("temp.converted_array")
    expect_equal(gtrack.info("temp.converted_array")$format, "indexed")

    # Verify data is intact
    r <- gextract("temp.converted_array", gintervals(c(1, 2)))
    r_orig <- gextract("test.array", gintervals(c(1, 2)))
    expect_equal(nrow(r), nrow(r_orig))
})

# ============================================================================
# Track expression tests with converted tracks
# ============================================================================

test_that("track expressions work with converted dense tracks", {
    withr::defer(gtrack.rm("temp.converted_expr", force = TRUE))
    gtrack.create("temp.converted_expr", "", "test.fixedbin")
    gtrack.convert_to_indexed("temp.converted_expr")

    # Test simple arithmetic expressions
    r1 <- gextract("temp.converted_expr + 1", gintervals(1, 0, 5000))
    r_orig <- gextract("test.fixedbin + 1", gintervals(1, 0, 5000))
    expect_equal(r1[[2]], r_orig[[2]])

    r2 <- gextract("temp.converted_expr * 2", gintervals(1, 0, 5000))
    r_orig2 <- gextract("test.fixedbin * 2", gintervals(1, 0, 5000))
    expect_equal(r2[[2]], r_orig2[[2]])

    r3 <- gextract("2 * temp.converted_expr + 0.5", gintervals(c(1, 2), 0, 10000))
    r_orig3 <- gextract("2 * test.fixedbin + 0.5", gintervals(c(1, 2), 0, 10000))
    expect_equal(r3[[2]], r_orig3[[2]])
})

test_that("track expressions work with converted sparse tracks", {
    withr::defer(gtrack.rm("temp.converted_sparse_expr", force = TRUE))
    gtrack.create("temp.converted_sparse_expr", "", "test.sparse", iterator = "test.sparse")
    gtrack.convert_to_indexed("temp.converted_sparse_expr")

    # Extract with expression
    r1 <- gextract("temp.converted_sparse_expr * 3", gintervals(c(1, 2)), iterator = "test.sparse")
    r_orig <- gextract("test.sparse * 3", gintervals(c(1, 2)), iterator = "test.sparse")
    expect_equal(r1[[2]], r_orig[[2]])
})

test_that("mixed expressions with converted and non-converted tracks", {
    withr::defer(gtrack.rm("temp.converted_mixed", force = TRUE))
    gtrack.create("temp.converted_mixed", "", "test.fixedbin")
    gtrack.convert_to_indexed("temp.converted_mixed")

    # Expression mixing converted and non-converted tracks
    r1 <- gextract("temp.converted_mixed + test.fixedbin", gintervals(1, 0, 5000))
    r_orig <- gextract("test.fixedbin + test.fixedbin", gintervals(1, 0, 5000))
    expect_equal(r1[[2]], r_orig[[2]])

    r2 <- gextract("temp.converted_mixed * test.fixedbin + 0.1", gintervals(c(1, 2), 0, 10000))
    r_orig2 <- gextract("test.fixedbin * test.fixedbin + 0.1", gintervals(c(1, 2), 0, 10000))
    expect_equal(r2[[2]], r_orig2[[2]])
})

test_that("gtrack.create with expression using converted track", {
    withr::defer(gtrack.rm("temp.converted_source", force = TRUE))
    withr::defer(gtrack.rm("temp.created_from_converted", force = TRUE))

    gtrack.create("temp.converted_source", "", "test.fixedbin")
    gtrack.convert_to_indexed("temp.converted_source")

    # Create new track from converted track expression
    gtrack.create("temp.created_from_converted", "", "temp.converted_source * 2 + 5")
    r1 <- gextract("temp.created_from_converted", gintervals(1, 0, 5000))
    r_orig <- gextract("test.fixedbin * 2 + 5", gintervals(1, 0, 5000))
    expect_equal(r1[[2]], r_orig[[2]])
})

# ============================================================================
# Virtual track tests with converted tracks
# ============================================================================

test_that("vtrack based on converted dense track works", {
    withr::defer(gtrack.rm("temp.converted_vtrack", force = TRUE))
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    gtrack.create("temp.converted_vtrack", "", "test.fixedbin")
    gtrack.convert_to_indexed("temp.converted_vtrack")

    gvtrack.create("v_converted", "temp.converted_vtrack")
    r <- gextract("v_converted", gintervals(c(1, 2)), iterator = 100)
    r_orig <- gextract("test.fixedbin", gintervals(c(1, 2)), iterator = 100)
    expect_equal(r$v_converted, r_orig$test.fixedbin)
})

test_that("vtrack with avg function on converted track", {
    withr::defer(gtrack.rm("temp.converted_vtrack_avg", force = TRUE))
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    gtrack.create("temp.converted_vtrack_avg", "", "test.fixedbin")
    gtrack.convert_to_indexed("temp.converted_vtrack_avg")

    gvtrack.create("v_converted_avg", "temp.converted_vtrack_avg", func = "avg")
    r <- gextract("v_converted_avg", gintervals(c(1, 2)), iterator = 500)
    r_orig <- gextract("test.fixedbin", gintervals(c(1, 2)), iterator = 500)

    # Both should have same number of rows
    expect_equal(nrow(r), nrow(r_orig))
})

test_that("vtrack extracts expressions with converted track", {
    withr::defer(gtrack.rm("temp.converted_vtrack_expr", force = TRUE))
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    gtrack.create("temp.converted_vtrack_expr", "", "test.fixedbin")
    gtrack.convert_to_indexed("temp.converted_vtrack_expr")

    gvtrack.create("v_converted_expr", "temp.converted_vtrack_expr")
    # Use expression in gextract
    r <- gextract("v_converted_expr * 2 + 1", gintervals(c(1, 2)), iterator = 233)
    r_orig <- gextract("test.fixedbin * 2 + 1", gintervals(c(1, 2)), iterator = 233)
    expect_equal(r[[2]], r_orig[[2]])
})

test_that("vtrack expression mixing converted and non-converted tracks", {
    withr::defer(gtrack.rm("temp.converted_vtrack_mix", force = TRUE))
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    gtrack.create("temp.converted_vtrack_mix", "", "test.fixedbin")
    gtrack.convert_to_indexed("temp.converted_vtrack_mix")

    gvtrack.create("v_mixed", "temp.converted_vtrack_mix")
    # Use expression in gextract mixing vtrack with regular track
    r <- gextract("v_mixed + test.fixedbin", gintervals(c(1, 2)), iterator = 100)
    r_orig <- gextract("test.fixedbin + test.fixedbin", gintervals(c(1, 2)), iterator = 100)
    expect_equal(r[[2]], r_orig[[2]])
})

test_that("vtrack with sparse iterator on converted track", {
    withr::defer(gtrack.rm("temp.converted_vtrack_sparse_it", force = TRUE))
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    gtrack.create("temp.converted_vtrack_sparse_it", "", "test.fixedbin")
    gtrack.convert_to_indexed("temp.converted_vtrack_sparse_it")

    gvtrack.create("v_sparse_it", "temp.converted_vtrack_sparse_it")
    r <- gextract("v_sparse_it", gintervals(c(1, 2)), iterator = "test.sparse")
    r_orig <- gextract("test.fixedbin", gintervals(c(1, 2)), iterator = "test.sparse")
    expect_equal(nrow(r), nrow(r_orig))
})

test_that("vtrack with array iterator on converted track", {
    withr::defer(gtrack.rm("temp.converted_vtrack_array_it", force = TRUE))
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    gtrack.create("temp.converted_vtrack_array_it", "", "test.fixedbin")
    gtrack.convert_to_indexed("temp.converted_vtrack_array_it")

    gvtrack.create("v_array_it", "temp.converted_vtrack_array_it", func = "avg")
    r <- gextract("v_array_it", gintervals(c(1, 2)), iterator = "test.array")
    expect_true(nrow(r) > 0)
})

# ============================================================================
# gtrack.modify tests with converted tracks
# ============================================================================

test_that("gtrack.modify works on converted tracks", {
    withr::defer(gtrack.rm("temp.converted_modify", force = TRUE))

    gtrack.create("temp.converted_modify", "", "test.fixedbin")
    gtrack.convert_to_indexed("temp.converted_modify")

    intervs <- gintervals(1, 1000, 5000)
    before <- gextract("temp.converted_modify", intervs)

    gtrack.modify("temp.converted_modify", "temp.converted_modify * 2", intervs)

    after <- gextract("temp.converted_modify", intervs)
    expect_equal(after$temp.converted_modify, before$temp.converted_modify * 2)
})

test_that("gtrack.modify with expression on converted track", {
    withr::defer(gtrack.rm("temp.converted_modify_expr", force = TRUE))

    gtrack.create("temp.converted_modify_expr", "", "test.fixedbin")
    gtrack.convert_to_indexed("temp.converted_modify_expr")

    intervs <- gintervals(c(1, 2), 0, 10000)
    before <- gextract("temp.converted_modify_expr", intervs)

    gtrack.modify("temp.converted_modify_expr", "temp.converted_modify_expr + test.fixedbin", intervs)

    after <- gextract("temp.converted_modify_expr", intervs)
    expected <- gextract("test.fixedbin", intervs)
    expect_equal(after$temp.converted_modify_expr, before$temp.converted_modify_expr + expected$test.fixedbin)
})

test_that("gtrack.modify on multiple chromosomes of converted track", {
    withr::defer(gtrack.rm("temp.converted_modify_multi", force = TRUE))

    gtrack.create("temp.converted_modify_multi", "", "test.fixedbin")
    gtrack.convert_to_indexed("temp.converted_modify_multi")

    intervs <- gintervals(c(1, 2, 3), 0, 50000)
    before <- gextract("temp.converted_modify_multi", intervs)

    gtrack.modify("temp.converted_modify_multi", "temp.converted_modify_multi * 3 + 1", intervs)

    after <- gextract("temp.converted_modify_multi", intervs)
    expect_equal(after$temp.converted_modify_multi, before$temp.converted_modify_multi * 3 + 1, tolerance = 1e-6)
})

# ============================================================================
# Iterator tests with converted tracks
# ============================================================================

test_that("numeric iterator works with converted track", {
    withr::defer(gtrack.rm("temp.converted_iter", force = TRUE))

    gtrack.create("temp.converted_iter", "", "test.fixedbin")
    gtrack.convert_to_indexed("temp.converted_iter")

    r <- gextract("temp.converted_iter", gintervals(1, 0, 10000), iterator = 500)
    r_orig <- gextract("test.fixedbin", gintervals(1, 0, 10000), iterator = 500)
    expect_equal(nrow(r), nrow(r_orig))
    expect_equal(r[[2]], r_orig[[2]])
})

test_that("sparse track iterator works with converted track", {
    withr::defer({
        gtrack.rm("temp.converted_sparse_iter", force = TRUE)
        gdb.reload()
    })

    gtrack.create("temp.converted_sparse_iter", "", "test.fixedbin")
    gtrack.convert_to_indexed("temp.converted_sparse_iter")

    r <- gextract("temp.converted_sparse_iter", gintervals(c(1, 2)), iterator = "test.sparse")
    r_orig <- gextract("test.fixedbin", gintervals(c(1, 2)), iterator = "test.sparse")
    expect_equal(nrow(r), nrow(r_orig))
})

test_that("array track iterator works with converted track", {
    withr::defer(gtrack.rm("temp.converted_array_iter", force = TRUE))

    gtrack.create("temp.converted_array_iter", "", "test.fixedbin")
    gtrack.convert_to_indexed("temp.converted_array_iter")

    r <- gextract("temp.converted_array_iter", gintervals(c(1, 2)), iterator = "test.array")
    r_orig <- gextract("test.fixedbin", gintervals(c(1, 2)), iterator = "test.array")
    expect_equal(nrow(r), nrow(r_orig))
})

test_that("giterator.intervals works with converted track", {
    withr::defer(gtrack.rm("temp.converted_giter", force = TRUE))

    gtrack.create("temp.converted_giter", "", "test.fixedbin")
    gtrack.convert_to_indexed("temp.converted_giter")

    intervs <- giterator.intervals("test.sparse", gintervals(c(1, 2)))
    r <- gextract("temp.converted_giter", gintervals(c(1, 2)), iterator = intervs)
    expect_true(nrow(r) > 0)
})

# ============================================================================
# Complex scenario tests
# ============================================================================

test_that("gscreen works with converted tracks", {
    withr::defer(gtrack.rm("temp.converted_screen", force = TRUE))

    gtrack.create("temp.converted_screen", "", "test.fixedbin")
    gtrack.convert_to_indexed("temp.converted_screen")

    r <- gscreen("temp.converted_screen > 0.2", gintervals(c(1, 2)))
    r_orig <- gscreen("test.fixedbin > 0.2", gintervals(c(1, 2)))
    expect_equal(nrow(r), nrow(r_orig))
})

test_that("gtrack.create with sparse iterator from converted track", {
    withr::defer(gtrack.rm("temp.converted_source_sparse", force = TRUE))
    withr::defer(gtrack.rm("temp.created_sparse_from_converted", force = TRUE))

    gtrack.create("temp.converted_source_sparse", "", "test.fixedbin")
    gtrack.convert_to_indexed("temp.converted_source_sparse")

    gtrack.create("temp.created_sparse_from_converted", "", "temp.converted_source_sparse + 1",
        iterator = "test.sparse"
    )
    r <- gextract("temp.created_sparse_from_converted", gintervals(c(1, 2)))
    expect_true(nrow(r) > 0)
})

test_that("gtrack.create with array iterator from converted track", {
    withr::defer(gtrack.rm("temp.converted_source_array", force = TRUE))
    withr::defer(gtrack.rm("temp.created_array_from_converted", force = TRUE))

    gtrack.create("temp.converted_source_array", "", "test.fixedbin")
    gtrack.convert_to_indexed("temp.converted_source_array")

    gtrack.create("temp.created_array_from_converted", "", "temp.converted_source_array * 2",
        iterator = "test.array"
    )
    r <- gextract("temp.created_array_from_converted", gintervals(c(1, 2)))
    expect_true(nrow(r) > 0)
})

test_that("multiple operations on same converted track", {
    withr::defer(gtrack.rm("temp.converted_multi_ops", force = TRUE))
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    # Create and convert
    gtrack.create("temp.converted_multi_ops", "", "test.fixedbin")
    gtrack.convert_to_indexed("temp.converted_multi_ops")

    # Extract with expression
    r1 <- gextract("temp.converted_multi_ops * 2", gintervals(1, 0, 5000))
    expect_true(nrow(r1) > 0)

    # Create vtrack and use expression in gextract
    gvtrack.create("v_multi", "temp.converted_multi_ops")
    r2 <- gextract("v_multi + 0.5", gintervals(1, 0, 5000), iterator = 100)
    expect_true(nrow(r2) > 0)

    # Modify
    intervs <- gintervals(1, 1000, 3000)
    gtrack.modify("temp.converted_multi_ops", "temp.converted_multi_ops * 1.5", intervs)
    r3 <- gextract("temp.converted_multi_ops", intervs)
    expect_true(nrow(r3) > 0)

    # Extract again
    r4 <- gextract("temp.converted_multi_ops", gintervals(1, 0, 5000))
    expect_true(nrow(r4) > 0)
})

test_that("converted track with ALLGENOME", {
    withr::defer(gtrack.rm("temp.converted_allgenome", force = TRUE))
    withr::local_options(gmax.data.size = 1e9)

    gtrack.create("temp.converted_allgenome", "", "test.fixedbin")
    gtrack.convert_to_indexed("temp.converted_allgenome")

    r <- gextract("temp.converted_allgenome", .misha$ALLGENOME)
    r_orig <- gextract("test.fixedbin", .misha$ALLGENOME)
    expect_equal(nrow(r), nrow(r_orig))
})
