create_isolated_test_db()

# Tests for sparse track liftover overlapping interval merging
# These tests verify the fix for the bug where sparse track liftover
# would create overlapping intervals in the output track when different
# source intervals mapped to overlapping (but not identical) target regions.

test_that("gtrack.liftover merges overlapping target intervals in sparse tracks", {
    local_db_state()

    # Create source genome with one chromosome
    source_db <- setup_source_db(list(paste0(">source1\n", paste(rep("A", 300), collapse = ""), "\n")))

    # Create sparse track with a single non-overlapping source interval
    # This will map to OVERLAPPING target regions through two different chains
    src_intervals <- data.frame(
        chrom = "chrsource1",
        start = 50,
        end = 60,
        stringsAsFactors = FALSE
    )
    src_values <- 10
    gtrack.create_sparse("src_track", "Source track", src_intervals, src_values)

    src_track_dir <- file.path(source_db, "tracks", "src_track.track")

    # Create target genome
    setup_db(list(paste0(">chr1\n", paste(rep("T", 300), collapse = ""), "\n")))

    # Create two overlapping chains that map the same source interval to overlapping targets:
    # Chain 1: source[40,70) -> target[0,30), so source[50,60) maps to target[10,20)
    # Chain 2: source[45,75) -> target[12,42), so source[50,60) maps to target[17,27)
    # Target intervals [10,20) and [17,27) overlap at [17,20)!
    chain_file <- new_chain_file()

    write_chain_entry(chain_file, "chrsource1", 300, "+", 40, 70, "chr1", 300, "+", 0, 30, 1)
    write_chain_entry(chain_file, "chrsource1", 300, "+", 45, 75, "chr1", 300, "+", 12, 42, 2)

    lifted_track <- "lifted_track"
    withr::defer({
        if (gtrack.exists(lifted_track)) gtrack.rm(lifted_track, force = TRUE)
    })

    gtrack.liftover(lifted_track, "Lifted track", src_track_dir, chain_file,
        src_overlap_policy = "keep", tgt_overlap_policy = "keep",
        multi_target_agg = "max"
    )

    # Extract and validate - the key test is that gextract works (no overlapping interval error)
    result <- gextract(lifted_track, gintervals.all())

    # Should have merged the overlapping intervals
    expect_true(nrow(result) >= 1)

    # Verify no overlapping intervals in the output
    if (nrow(result) > 1) {
        result <- result[order(result$start), ]
        for (i in 1:(nrow(result) - 1)) {
            expect_true(result$end[i] <= result$start[i + 1],
                info = sprintf(
                    "Interval %d [%d,%d) overlaps with interval %d [%d,%d)",
                    i, result$start[i], result$end[i],
                    i + 1, result$start[i + 1], result$end[i + 1]
                )
            )
        }
    }
})

test_that("gtrack.liftover correctly aggregates values when merging overlapping intervals", {
    local_db_state()

    # Create source genome
    source_db <- setup_source_db(list(paste0(">source1\n", paste(rep("A", 200), collapse = ""), "\n")))

    # Create sparse track with intervals that will overlap in target
    src_intervals <- data.frame(
        chrom = c("chrsource1", "chrsource1", "chrsource1"),
        start = c(10, 30, 50),
        end = c(20, 40, 60),
        stringsAsFactors = FALSE
    )
    src_values <- c(100, 200, 300)
    gtrack.create_sparse("src_track", "Source track", src_intervals, src_values)

    src_track_dir <- file.path(source_db, "tracks", "src_track.track")

    # Create target genome
    setup_db(list(paste0(">chr1\n", paste(rep("T", 200), collapse = ""), "\n")))

    # Create a single chain that maps all source intervals to the SAME target region
    # This guarantees overlapping output intervals
    chain_file <- new_chain_file()

    # All three source intervals [10,20), [30,40), [50,60) map to chr1[0,10)
    # Through different chains with overlapping source regions
    write_chain_entry(chain_file, "chrsource1", 200, "+", 0, 30, "chr1", 200, "+", 0, 30, 1)
    write_chain_entry(chain_file, "chrsource1", 200, "+", 20, 50, "chr1", 200, "+", 10, 40, 2)
    write_chain_entry(chain_file, "chrsource1", 200, "+", 40, 70, "chr1", 200, "+", 20, 50, 3)

    # Test with max aggregation
    lifted_track_max <- "lifted_max"
    withr::defer({
        if (gtrack.exists(lifted_track_max)) gtrack.rm(lifted_track_max, force = TRUE)
    })

    gtrack.liftover(lifted_track_max, "Lifted max", src_track_dir, chain_file,
        src_overlap_policy = "keep", tgt_overlap_policy = "keep",
        multi_target_agg = "max"
    )

    result_max <- gextract(lifted_track_max, gintervals.all())
    expect_true(nrow(result_max) >= 1)

    # Verify no overlapping intervals
    if (nrow(result_max) > 1) {
        result_max <- result_max[order(result_max$start), ]
        for (i in 1:(nrow(result_max) - 1)) {
            expect_true(result_max$end[i] <= result_max$start[i + 1])
        }
    }

    # Test with mean aggregation
    lifted_track_mean <- "lifted_mean"
    withr::defer({
        if (gtrack.exists(lifted_track_mean)) gtrack.rm(lifted_track_mean, force = TRUE)
    })

    gtrack.liftover(lifted_track_mean, "Lifted mean", src_track_dir, chain_file,
        src_overlap_policy = "keep", tgt_overlap_policy = "keep",
        multi_target_agg = "mean"
    )

    result_mean <- gextract(lifted_track_mean, gintervals.all())
    expect_true(nrow(result_mean) >= 1)

    # Verify no overlapping intervals
    if (nrow(result_mean) > 1) {
        result_mean <- result_mean[order(result_mean$start), ]
        for (i in 1:(nrow(result_mean) - 1)) {
            expect_true(result_mean$end[i] <= result_mean$start[i + 1])
        }
    }
})

test_that("gtrack.liftover handles sparse track with many overlapping source-to-target mappings", {
    local_db_state()

    # Create source genome
    source_db <- setup_source_db(list(paste0(">source1\n", paste(rep("A", 1000), collapse = ""), "\n")))

    # Create sparse track with many small intervals
    n_intervals <- 20
    src_intervals <- data.frame(
        chrom = rep("chrsource1", n_intervals),
        start = seq(100, by = 10, length.out = n_intervals),
        end = seq(105, by = 10, length.out = n_intervals),
        stringsAsFactors = FALSE
    )
    src_values <- seq(1, n_intervals)
    gtrack.create_sparse("src_track", "Source track", src_intervals, src_values)

    src_track_dir <- file.path(source_db, "tracks", "src_track.track")

    # Create target genome
    setup_db(list(paste0(">chr1\n", paste(rep("T", 1000), collapse = ""), "\n")))

    # Create multiple overlapping chains that will cause many target overlaps
    chain_file <- new_chain_file()

    # Chain 1: wide mapping
    write_chain_entry(chain_file, "chrsource1", 1000, "+", 50, 400, "chr1", 1000, "+", 0, 350, 1)

    # Chain 2: overlapping mapping shifted by 50
    write_chain_entry(chain_file, "chrsource1", 1000, "+", 100, 450, "chr1", 1000, "+", 30, 380, 2)

    # Chain 3: another overlapping mapping
    write_chain_entry(chain_file, "chrsource1", 1000, "+", 150, 500, "chr1", 1000, "+", 60, 410, 3)

    lifted_track <- "lifted_many"
    withr::defer({
        if (gtrack.exists(lifted_track)) gtrack.rm(lifted_track, force = TRUE)
    })

    gtrack.liftover(lifted_track, "Lifted many", src_track_dir, chain_file,
        src_overlap_policy = "keep", tgt_overlap_policy = "keep",
        multi_target_agg = "max"
    )

    # The key test: gextract should work without "Invalid format" error
    result <- gextract(lifted_track, gintervals.all())

    expect_true(nrow(result) >= 1)

    # Verify no overlapping intervals in the output
    if (nrow(result) > 1) {
        result <- result[order(result$start), ]
        for (i in 1:(nrow(result) - 1)) {
            expect_true(result$end[i] <= result$start[i + 1],
                info = sprintf("Overlapping intervals found at positions %d and %d", i, i + 1)
            )
        }
    }
})

test_that("gtrack.liftover sparse track produces valid track that can be read without errors", {
    local_db_state()

    # This test specifically checks that the output track format is valid
    # (the original bug caused "Invalid format of a sparse track file" errors)

    source_db <- setup_source_db(list(paste0(">source1\n", paste(rep("A", 500), collapse = ""), "\n")))

    # Create sparse track
    src_intervals <- data.frame(
        chrom = c("chrsource1", "chrsource1", "chrsource1", "chrsource1"),
        start = c(50, 100, 150, 200),
        end = c(75, 125, 175, 225),
        stringsAsFactors = FALSE
    )
    src_values <- c(1.5, 2.5, 3.5, 4.5)
    gtrack.create_sparse("src_track", "Source track", src_intervals, src_values)

    src_track_dir <- file.path(source_db, "tracks", "src_track.track")

    setup_db(list(paste0(">chr1\n", paste(rep("T", 500), collapse = ""), "\n")))

    # Create chains that produce overlapping target intervals
    chain_file <- new_chain_file()
    write_chain_entry(chain_file, "chrsource1", 500, "+", 0, 150, "chr1", 500, "+", 0, 150, 1)
    write_chain_entry(chain_file, "chrsource1", 500, "+", 100, 250, "chr1", 500, "+", 80, 230, 2)

    lifted_track <- "lifted_valid"
    withr::defer({
        if (gtrack.exists(lifted_track)) gtrack.rm(lifted_track, force = TRUE)
    })

    gtrack.liftover(lifted_track, "Lifted valid", src_track_dir, chain_file,
        src_overlap_policy = "keep", tgt_overlap_policy = "keep",
        multi_target_agg = "mean"
    )

    # Verify track info works
    info <- gtrack.info(lifted_track)
    expect_equal(info$type, "sparse")

    # Verify gextract works on full chromosome
    result_full <- gextract(lifted_track, gintervals.all())
    expect_true(!is.null(result_full))

    # Verify gextract works on partial intervals
    result_partial <- gextract(lifted_track, gintervals("chr1", 50, 200))
    expect_true(!is.null(result_partial))

    # Verify gsummary works
    summary_result <- gsummary(lifted_track)
    expect_true(!is.null(summary_result))
})

test_that("gtrack.liftover merges adjacent intervals with same value", {
    local_db_state()

    # Test that adjacent intervals with the same value get merged

    source_db <- setup_source_db(list(paste0(">source1\n", paste(rep("A", 200), collapse = ""), "\n")))

    # Create sparse track with two adjacent intervals with SAME value
    src_intervals <- data.frame(
        chrom = c("chrsource1", "chrsource1"),
        start = c(10, 20),
        end = c(20, 30),
        stringsAsFactors = FALSE
    )
    src_values <- c(42, 42) # Same value
    gtrack.create_sparse("src_track", "Source track", src_intervals, src_values)

    src_track_dir <- file.path(source_db, "tracks", "src_track.track")

    setup_db(list(paste0(">chr1\n", paste(rep("T", 200), collapse = ""), "\n")))

    chain_file <- new_chain_file()
    write_chain_entry(chain_file, "chrsource1", 200, "+", 0, 50, "chr1", 200, "+", 0, 50, 1)

    lifted_track <- "lifted_adjacent"
    withr::defer({
        if (gtrack.exists(lifted_track)) gtrack.rm(lifted_track, force = TRUE)
    })

    gtrack.liftover(lifted_track, "Lifted adjacent", src_track_dir, chain_file)

    result <- gextract(lifted_track, gintervals.all())

    # Should have valid output
    expect_true(nrow(result) >= 1)

    # All values should be 42
    expect_true(all(result[[lifted_track]] == 42))
})

test_that("gtrack.liftover with sum aggregation correctly sums overlapping values", {
    local_db_state()

    source_db <- setup_source_db(list(paste0(">source1\n", paste(rep("A", 200), collapse = ""), "\n")))

    # Single source interval that will map to overlapping targets through multiple chains
    src_intervals <- data.frame(
        chrom = "chrsource1",
        start = 10,
        end = 20,
        stringsAsFactors = FALSE
    )
    src_values <- 100
    gtrack.create_sparse("src_track", "Source track", src_intervals, src_values)

    src_track_dir <- file.path(source_db, "tracks", "src_track.track")

    setup_db(list(paste0(">chr1\n", paste(rep("T", 200), collapse = ""), "\n")))

    # Two overlapping chains that map the same source interval to overlapping targets:
    # Chain 1: source[0,30) -> target[0,30), so source[10,20) maps to target[10,20)
    # Chain 2: source[5,35) -> target[8,38), so source[10,20) maps to target[13,23)
    # Target intervals [10,20) and [13,23) overlap at [13,20)!
    chain_file <- new_chain_file()
    write_chain_entry(chain_file, "chrsource1", 200, "+", 0, 30, "chr1", 200, "+", 0, 30, 1)
    write_chain_entry(chain_file, "chrsource1", 200, "+", 5, 35, "chr1", 200, "+", 8, 38, 2)

    lifted_track <- "lifted_sum"
    withr::defer({
        if (gtrack.exists(lifted_track)) gtrack.rm(lifted_track, force = TRUE)
    })

    gtrack.liftover(lifted_track, "Lifted sum", src_track_dir, chain_file,
        src_overlap_policy = "keep", tgt_overlap_policy = "keep",
        multi_target_agg = "sum"
    )

    result <- gextract(lifted_track, gintervals.all())

    expect_true(nrow(result) >= 1)

    # Verify no overlapping intervals
    if (nrow(result) > 1) {
        result <- result[order(result$start), ]
        for (i in 1:(nrow(result) - 1)) {
            expect_true(result$end[i] <= result$start[i + 1])
        }
    }
})
