interv2path <- function(intervname) {
    path <- gsub("\\.", "/", intervname)
    return(paste0(get("GWD", envir = .misha), "/", path, ".interv"))
}

test_that("2D intervals functionality is identical between Legacy and Indexed formats", {
    # 1. Setup a standard test environment
    gdb.init_examples()
    gdir.create("test", showWarnings = FALSE)

    # Define standard 2D intervals (need >10 intervals to trigger big set creation)
    intervs <- gintervals.2d(
        chroms1 = rep(c(1, 1, 2), each = 4),
        starts1 = rep(c(100, 200, 300, 400), 3),
        ends1 = rep(c(150, 250, 350, 450), 3),
        chroms2 = rep(c(1, 2, 2), each = 4),
        starts2 = rep(c(300, 100, 400, 500), 3),
        ends2 = rep(c(350, 150, 450, 550), 3)
    )

    # 2. Save both as Big Sets
    # We use withr::local_options to force 'gmax.data.size' to be small (10 intervals)
    # just for the scope of this test block. This ensures gintervals.save creates
    # Big Intervals Sets (directories) instead of small serialized files.

    # First, save the legacy format with indexed format temporarily disabled
    withr::with_options(list(gmax.data.size = 10, gmulticontig.indexed_format = FALSE), {
        gintervals.save("test.legacy_2d", intervs)
    })

    # Then save another copy for conversion to indexed
    withr::with_options(list(gmax.data.size = 10, gmulticontig.indexed_format = FALSE), {
        gintervals.save("test.indexed_2d", intervs)
    })

    # 3. Convert one copy to the new Indexed format
    misha::gintervals.2d.convert_to_indexed("test.indexed_2d", remove.old = TRUE, force = TRUE)

    # 4. Verify Physical Structure (Sanity Check)
    # We use misha's internal helper or construct the path manually as you defined
    legacy_dir <- interv2path("test.legacy_2d")

    expect_true(dir.exists(legacy_dir), info = "Legacy set should be a directory (Big Set)")
    expect_true(file.exists(file.path(legacy_dir, "chr1-chr1")))
    expect_false(file.exists(file.path(legacy_dir, "intervals2d.dat")))

    # Indexed should have the single DAT file
    indexed_dir <- interv2path("test.indexed_2d")
    expect_true(dir.exists(indexed_dir), info = "Indexed set should be a directory (Big Set)")
    expect_false(file.exists(file.path(indexed_dir, "chr1-chr1")))
    expect_true(file.exists(file.path(indexed_dir, "intervals2d.dat")))
    expect_true(file.exists(file.path(indexed_dir, "intervals2d.idx")))

    # 5. PARITY CHECK: Loading
    # Temporarily increase gmax.data.size to allow loading the full big set
    loaded_legacy <- withr::with_options(list(gmax.data.size = 1e9), {
        gintervals.load("test.legacy_2d")
    })
    loaded_indexed <- withr::with_options(list(gmax.data.size = 1e9), {
        gintervals.load("test.indexed_2d")
    })

    # Compare raw data frames (ignore row names)
    rownames(loaded_legacy) <- NULL
    rownames(loaded_indexed) <- NULL
    expect_equal(loaded_legacy, loaded_indexed, info = "Loading produces identical data frames")

    # 6. PARITY CHECK: Neighbors Analysis
    # This exercises the iterator logic
    nb_legacy <- gintervals.neighbors("test.legacy_2d", "test.legacy_2d", maxdist = 1000)
    nb_indexed <- gintervals.neighbors("test.indexed_2d", "test.indexed_2d", maxdist = 1000)

    # Normalize result (remove track names from columns for direct comparison)
    colnames(nb_legacy) <- gsub("test.legacy_2d", "track", colnames(nb_legacy))
    colnames(nb_indexed) <- gsub("test.indexed_2d", "track", colnames(nb_indexed))

    expect_equal(nb_legacy, nb_indexed, info = "gintervals.neighbors produces identical results")

    # Cleanup
    gintervals.rm("test.legacy_2d", force = TRUE)
    gintervals.rm("test.indexed_2d", force = TRUE)
})

test_that("2D intervals output regression", {
    gdb.init_examples()
    gdir.create("test", showWarnings = FALSE)

    # Create a deterministic 2D set
    intervs <- gintervals.2d(
        chroms1 = c(1, "X"), starts1 = c(1000, 2000), ends1 = c(2000, 3000),
        chroms2 = c(1, 2), starts2 = c(5000, 100), ends2 = c(6000, 500)
    )
    gintervals.save("test.reg_2d_parity", intervs)

    # 1. Check basic loading
    loaded <- gintervals.load("test.reg_2d_parity")
    expect_regression(loaded, "gintervals_load_2d_basic")

    # 2. Check gextract (2D track extraction)
    # We create a dummy 2D track for extraction
    gtrack.2d.create("test.track_2d", "Test Track", intervs, values = c(1.5, 2.5))

    extracted <- gextract("test.track_2d", intervs)
    expect_regression(extracted, "gextract_2d_basic")

    # Cleanup
    gtrack.rm("test.track_2d", force = TRUE)
    gintervals.rm("test.reg_2d_parity", force = TRUE)
})

test_that("2D intervals with extra columns maintain parity", {
    gdb.init_examples()
    gdir.create("test", showWarnings = FALSE)

    # Create intervals with additional columns
    intervs <- gintervals.2d(
        chroms1 = rep(c(1, 2), each = 6),
        starts1 = rep(c(100, 200, 300, 400, 500, 600), 2),
        ends1 = rep(c(150, 250, 350, 450, 550, 650), 2),
        chroms2 = rep(c(1, 2), each = 6),
        starts2 = rep(c(300, 400, 500, 600, 700, 800), 2),
        ends2 = rep(c(350, 450, 550, 650, 750, 850), 2)
    )

    # Add extra columns
    intervs$score <- runif(nrow(intervs))
    intervs$name <- paste0("interval_", seq_len(nrow(intervs)))
    intervs$category <- factor(sample(c("A", "B", "C"), nrow(intervs), replace = TRUE))

    # Save both formats
    withr::with_options(list(gmax.data.size = 10, gmulticontig.indexed_format = FALSE), {
        gintervals.save("test.legacy_2d_extra", intervs)
        gintervals.save("test.indexed_2d_extra", intervs)
    })

    gintervals.2d.convert_to_indexed("test.indexed_2d_extra", remove.old = TRUE, force = TRUE)

    # Load and compare
    loaded_legacy <- withr::with_options(list(gmax.data.size = 1e9), {
        gintervals.load("test.legacy_2d_extra")
    })
    loaded_indexed <- withr::with_options(list(gmax.data.size = 1e9), {
        gintervals.load("test.indexed_2d_extra")
    })

    rownames(loaded_legacy) <- NULL
    rownames(loaded_indexed) <- NULL

    expect_equal(loaded_legacy, loaded_indexed, info = "Extra columns preserved identically")
    expect_equal(colnames(loaded_legacy), colnames(loaded_indexed))
    expect_equal(class(loaded_legacy$category), class(loaded_indexed$category))

    # Cleanup
    gintervals.rm("test.legacy_2d_extra", force = TRUE)
    gintervals.rm("test.indexed_2d_extra", force = TRUE)
})

test_that("2D intervals with sparse chromosome pairs maintain parity", {
    gdb.init_examples()
    gdir.create("test", showWarnings = FALSE)

    # Create intervals with only some chromosome pairs populated
    # chr1-chr1: 4 intervals
    # chr1-chr2: empty
    # chr2-chr2: 8 intervals
    intervs <- gintervals.2d(
        chroms1 = c(rep(1, 4), rep(2, 8)),
        starts1 = c(seq(100, 400, 100), seq(100, 800, 100)),
        ends1 = c(seq(150, 450, 100), seq(150, 850, 100)),
        chroms2 = c(rep(1, 4), rep(2, 8)),
        starts2 = c(seq(500, 800, 100), seq(200, 900, 100)),
        ends2 = c(seq(550, 850, 100), seq(250, 950, 100))
    )

    # Save both formats
    withr::with_options(list(gmax.data.size = 10, gmulticontig.indexed_format = FALSE), {
        gintervals.save("test.legacy_2d_sparse", intervs)
        gintervals.save("test.indexed_2d_sparse", intervs)
    })

    gintervals.2d.convert_to_indexed("test.indexed_2d_sparse", remove.old = TRUE, force = TRUE)

    # Test loading specific chromosome pairs
    loaded_legacy_11 <- gintervals.load("test.legacy_2d_sparse", chrom1 = 1, chrom2 = 1)
    loaded_indexed_11 <- gintervals.load("test.indexed_2d_sparse", chrom1 = 1, chrom2 = 1)

    rownames(loaded_legacy_11) <- NULL
    rownames(loaded_indexed_11) <- NULL
    expect_equal(loaded_legacy_11, loaded_indexed_11, info = "chr1-chr1 identical")
    expect_equal(nrow(loaded_legacy_11), 4)

    # Test empty pair
    loaded_legacy_12 <- gintervals.load("test.legacy_2d_sparse", chrom1 = 1, chrom2 = 2)
    loaded_indexed_12 <- gintervals.load("test.indexed_2d_sparse", chrom1 = 1, chrom2 = 2)

    expect_equal(nrow(loaded_legacy_12), 0)
    expect_equal(nrow(loaded_indexed_12), 0)
    expect_equal(colnames(loaded_legacy_12), colnames(loaded_indexed_12))

    # Test chr2-chr2
    loaded_legacy_22 <- gintervals.load("test.legacy_2d_sparse", chrom1 = 2, chrom2 = 2)
    loaded_indexed_22 <- gintervals.load("test.indexed_2d_sparse", chrom1 = 2, chrom2 = 2)

    rownames(loaded_legacy_22) <- NULL
    rownames(loaded_indexed_22) <- NULL
    expect_equal(loaded_legacy_22, loaded_indexed_22, info = "chr2-chr2 identical")
    expect_equal(nrow(loaded_legacy_22), 8)

    # Cleanup
    gintervals.rm("test.legacy_2d_sparse", force = TRUE)
    gintervals.rm("test.indexed_2d_sparse", force = TRUE)
})

test_that("2D intervals gintervals.chrom_sizes produces identical results", {
    gdb.init_examples()
    gdir.create("test", showWarnings = FALSE)

    intervs <- gintervals.2d(
        chroms1 = rep(c(1, 1, 2), each = 4),
        starts1 = rep(c(100, 200, 300, 400), 3),
        ends1 = rep(c(150, 250, 350, 450), 3),
        chroms2 = rep(c(1, 2, 2), each = 4),
        starts2 = rep(c(300, 100, 400, 500), 3),
        ends2 = rep(c(350, 150, 450, 550), 3)
    )

    withr::with_options(list(gmax.data.size = 10, gmulticontig.indexed_format = FALSE), {
        gintervals.save("test.legacy_2d_sizes", intervs)
        gintervals.save("test.indexed_2d_sizes", intervs)
    })

    gintervals.2d.convert_to_indexed("test.indexed_2d_sizes", remove.old = TRUE, force = TRUE)

    sizes_legacy <- gintervals.chrom_sizes("test.legacy_2d_sizes")
    sizes_indexed <- gintervals.chrom_sizes("test.indexed_2d_sizes")

    expect_equal(sizes_legacy, sizes_indexed, info = "chrom_sizes identical")
    expect_equal(nrow(sizes_legacy), 3) # chr1-chr1, chr1-chr2, chr2-chr2
    expect_true(all(c("chrom1", "chrom2", "size") %in% colnames(sizes_legacy)))

    # Cleanup
    gintervals.rm("test.legacy_2d_sizes", force = TRUE)
    gintervals.rm("test.indexed_2d_sizes", force = TRUE)
})

test_that("2D intervals edge cases: single interval per pair", {
    gdb.init_examples()
    gdir.create("test", showWarnings = FALSE)

    # Single interval for each chromosome pair (minimal big set)
    intervs <- gintervals.2d(
        chroms1 = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
        starts1 = rep(100, 12),
        ends1 = rep(150, 12),
        chroms2 = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
        starts2 = seq(200, 1300, 100),
        ends2 = seq(250, 1350, 100)
    )

    withr::with_options(list(gmax.data.size = 10, gmulticontig.indexed_format = FALSE), {
        gintervals.save("test.legacy_2d_single", intervs)
        gintervals.save("test.indexed_2d_single", intervs)
    })

    gintervals.2d.convert_to_indexed("test.indexed_2d_single", remove.old = TRUE, force = TRUE)

    loaded_legacy <- withr::with_options(list(gmax.data.size = 1e9), {
        gintervals.load("test.legacy_2d_single")
    })
    loaded_indexed <- withr::with_options(list(gmax.data.size = 1e9), {
        gintervals.load("test.indexed_2d_single")
    })

    rownames(loaded_legacy) <- NULL
    rownames(loaded_indexed) <- NULL
    expect_equal(loaded_legacy, loaded_indexed)

    # Cleanup
    gintervals.rm("test.legacy_2d_single", force = TRUE)
    gintervals.rm("test.indexed_2d_single", force = TRUE)
})

test_that("2D intervals update operations maintain parity", {
    gdb.init_examples()
    gdir.create("test", showWarnings = FALSE)

    # Initial intervals
    intervs <- gintervals.2d(
        chroms1 = rep(c(1, 2), each = 6),
        starts1 = rep(seq(100, 600, 100), 2),
        ends1 = rep(seq(150, 650, 100), 2),
        chroms2 = rep(c(1, 2), each = 6),
        starts2 = rep(seq(200, 700, 100), 2),
        ends2 = rep(seq(250, 750, 100), 2)
    )

    withr::with_options(list(gmax.data.size = 10, gmulticontig.indexed_format = FALSE), {
        gintervals.save("test.legacy_2d_update", intervs)
        gintervals.save("test.indexed_2d_update", intervs)
    })

    gintervals.2d.convert_to_indexed("test.indexed_2d_update", remove.old = TRUE, force = TRUE)

    # Update chr1-chr1 with new intervals
    new_intervs <- gintervals.2d(
        chroms1 = rep(1, 3),
        starts1 = c(1000, 2000, 3000),
        ends1 = c(1100, 2100, 3100),
        chroms2 = rep(1, 3),
        starts2 = c(5000, 6000, 7000),
        ends2 = c(5100, 6100, 7100)
    )

    gintervals.update("test.legacy_2d_update", new_intervs, chrom1 = 1, chrom2 = 1)
    gintervals.update("test.indexed_2d_update", new_intervs, chrom1 = 1, chrom2 = 1)

    # Load and compare
    loaded_legacy <- withr::with_options(list(gmax.data.size = 1e9), {
        gintervals.load("test.legacy_2d_update")
    })
    loaded_indexed <- withr::with_options(list(gmax.data.size = 1e9), {
        gintervals.load("test.indexed_2d_update")
    })

    rownames(loaded_legacy) <- NULL
    rownames(loaded_indexed) <- NULL
    expect_equal(loaded_legacy, loaded_indexed, info = "After update, both formats identical")

    # Check that chr1-chr1 was updated
    loaded_11 <- gintervals.load("test.indexed_2d_update", chrom1 = 1, chrom2 = 1)
    expect_equal(nrow(loaded_11), 3)
    expect_equal(loaded_11$start1, c(1000, 2000, 3000))

    # Cleanup
    gintervals.rm("test.legacy_2d_update", force = TRUE)
    gintervals.rm("test.indexed_2d_update", force = TRUE)
})

test_that("2D intervals with all chromosome combinations", {
    gdb.init_examples()
    gdir.create("test", showWarnings = FALSE)

    # Create intervals for all possible chromosome pairs in test db (chr1, chr2, chrX)
    intervs <- gintervals.2d(
        chroms1 = c(rep(1, 5), rep(1, 5), rep(1, 5), rep(2, 5), rep(2, 5), rep("X", 5)),
        starts1 = rep(seq(100, 500, 100), 6),
        ends1 = rep(seq(150, 550, 100), 6),
        chroms2 = c(rep(1, 5), rep(2, 5), rep("X", 5), rep(2, 5), rep("X", 5), rep("X", 5)),
        starts2 = rep(seq(200, 600, 100), 6),
        ends2 = rep(seq(250, 650, 100), 6)
    )

    withr::with_options(list(gmax.data.size = 10, gmulticontig.indexed_format = FALSE), {
        gintervals.save("test.legacy_2d_allcomb", intervs)
        gintervals.save("test.indexed_2d_allcomb", intervs)
    })

    gintervals.2d.convert_to_indexed("test.indexed_2d_allcomb", remove.old = TRUE, force = TRUE)

    # Test all pairs
    pairs <- list(
        c(1, 1), c(1, 2), c(1, "X"),
        c(2, 2), c(2, "X"),
        c("X", "X")
    )

    for (pair in pairs) {
        loaded_legacy <- gintervals.load("test.legacy_2d_allcomb", chrom1 = pair[1], chrom2 = pair[2])
        loaded_indexed <- gintervals.load("test.indexed_2d_allcomb", chrom1 = pair[1], chrom2 = pair[2])

        rownames(loaded_legacy) <- NULL
        rownames(loaded_indexed) <- NULL
        expect_equal(loaded_legacy, loaded_indexed,
            info = sprintf("Pair %s-%s identical", pair[1], pair[2])
        )
        expect_equal(nrow(loaded_legacy), 5)
    }

    # Cleanup
    gintervals.rm("test.legacy_2d_allcomb", force = TRUE)
    gintervals.rm("test.indexed_2d_allcomb", force = TRUE)
})

test_that("2D intervals work as iterators in gextract", {
    gdb.init_examples()
    gdir.create("test", showWarnings = FALSE)

    # Create intervals
    intervs <- gintervals.2d(
        chroms1 = rep(c(1, 2), each = 6),
        starts1 = rep(seq(100, 600, 100), 2),
        ends1 = rep(seq(150, 650, 100), 2),
        chroms2 = rep(c(1, 2), each = 6),
        starts2 = rep(seq(200, 700, 100), 2),
        ends2 = rep(seq(250, 750, 100), 2)
    )

    # Create a 2D track
    values <- runif(nrow(intervs))
    gtrack.2d.create("test.track_2d_iter", "Test Track", intervs, values)

    # Save intervals in both formats
    withr::with_options(list(gmax.data.size = 10, gmulticontig.indexed_format = FALSE), {
        gintervals.save("test.legacy_2d_iter", intervs)
        gintervals.save("test.indexed_2d_iter", intervs)
    })

    gintervals.2d.convert_to_indexed("test.indexed_2d_iter", remove.old = TRUE, force = TRUE)

    # Use both as iterators in gextract
    result_legacy <- gextract("test.track_2d_iter",
        intervals = "test.legacy_2d_iter",
        iterator = "test.legacy_2d_iter"
    )
    result_indexed <- gextract("test.track_2d_iter",
        intervals = "test.indexed_2d_iter",
        iterator = "test.indexed_2d_iter"
    )

    # Compare results
    rownames(result_legacy) <- NULL
    rownames(result_indexed) <- NULL
    expect_equal(result_legacy, result_indexed, info = "gextract with iterator produces identical results")
    expect_equal(nrow(result_legacy), nrow(intervs))

    # Cleanup
    gtrack.rm("test.track_2d_iter", force = TRUE)
    gintervals.rm("test.legacy_2d_iter", force = TRUE)
    gintervals.rm("test.indexed_2d_iter", force = TRUE)
})

test_that("2D intervals work with gscreen", {
    gdb.init_examples()
    gdir.create("test", showWarnings = FALSE)

    # Create intervals and a track
    intervs <- gintervals.2d(
        chroms1 = rep(c(1, 2), each = 6),
        starts1 = rep(seq(100, 600, 100), 2),
        ends1 = rep(seq(150, 650, 100), 2),
        chroms2 = rep(c(1, 2), each = 6),
        starts2 = rep(seq(200, 700, 100), 2),
        ends2 = rep(seq(250, 750, 100), 2)
    )

    # Create track with values
    values <- seq(0.1, 1.2, length.out = nrow(intervs))
    gtrack.2d.create("test.track_2d_screen", "Test Track", intervs, values)

    # Save intervals in both formats
    withr::with_options(list(gmax.data.size = 10, gmulticontig.indexed_format = FALSE), {
        gintervals.save("test.legacy_2d_screen", intervs)
        gintervals.save("test.indexed_2d_screen", intervs)
    })

    gintervals.2d.convert_to_indexed("test.indexed_2d_screen", remove.old = TRUE, force = TRUE)

    # Use with gscreen
    screen_legacy <- gscreen("test.track_2d_screen > 0.5", intervals = "test.legacy_2d_screen")
    screen_indexed <- gscreen("test.track_2d_screen > 0.5", intervals = "test.indexed_2d_screen")

    # Compare results
    rownames(screen_legacy) <- NULL
    rownames(screen_indexed) <- NULL
    expect_equal(screen_legacy, screen_indexed, info = "gscreen produces identical results")
    expect_true(nrow(screen_legacy) > 0)
    expect_true(nrow(screen_legacy) < nrow(intervs))

    # Cleanup
    gtrack.rm("test.track_2d_screen", force = TRUE)
    gintervals.rm("test.legacy_2d_screen", force = TRUE)
    gintervals.rm("test.indexed_2d_screen", force = TRUE)
})

test_that("2D intervals work with virtual tracks", {
    gdb.init_examples()
    gdir.create("test", showWarnings = FALSE)

    # Create base intervals and tracks
    intervs <- gintervals.2d(
        chroms1 = rep(c(1, 2), each = 6),
        starts1 = rep(seq(100, 600, 100), 2),
        ends1 = rep(seq(150, 650, 100), 2),
        chroms2 = rep(c(1, 2), each = 6),
        starts2 = rep(seq(200, 700, 100), 2),
        ends2 = rep(seq(250, 750, 100), 2)
    )

    # Create two tracks for virtual track
    values1 <- seq(1, 12, 1)
    gtrack.2d.create("test.track_2d_vt1", "Track 1", intervs, values1)

    values2 <- seq(2, 24, 2)
    gtrack.2d.create("test.track_2d_vt2", "Track 2", intervs, values2)

    # Create virtual track from one track
    gvtrack.create("test.vtrack_2d", "test.track_2d_vt1")

    # Save intervals in both formats
    withr::with_options(list(gmax.data.size = 10, gmulticontig.indexed_format = FALSE), {
        gintervals.save("test.legacy_2d_vt", intervs)
        gintervals.save("test.indexed_2d_vt", intervs)
    })

    gintervals.2d.convert_to_indexed("test.indexed_2d_vt", remove.old = TRUE, force = TRUE)

    # Extract using virtual track combined with the other track in expression
    # Need explicit iterator when expression contains multiple 2D tracks
    result_legacy <- gextract("test.vtrack_2d + test.track_2d_vt2",
        intervals = "test.legacy_2d_vt",
        iterator = "test.legacy_2d_vt"
    )
    result_indexed <- gextract("test.vtrack_2d + test.track_2d_vt2",
        intervals = "test.indexed_2d_vt",
        iterator = "test.indexed_2d_vt"
    )

    # Compare results
    rownames(result_legacy) <- NULL
    rownames(result_indexed) <- NULL
    expect_equal(result_legacy, result_indexed, info = "Virtual track extraction identical")

    # Verify virtual track calculation (sum of the two tracks)
    # The column name will be the expression, not just the virtual track name
    expr_col <- colnames(result_legacy)[colnames(result_legacy) != "chrom1" & colnames(result_legacy) != "start1" & colnames(result_legacy) != "end1" & colnames(result_legacy) != "chrom2" & colnames(result_legacy) != "start2" & colnames(result_legacy) != "end2"]
    expect_equal(result_legacy[[expr_col[1]]], seq(3, 36, 3))

    # Cleanup
    gvtrack.rm("test.vtrack_2d")
    gtrack.rm("test.track_2d_vt1", force = TRUE)
    gtrack.rm("test.track_2d_vt2", force = TRUE)
    gintervals.rm("test.legacy_2d_vt", force = TRUE)
    gintervals.rm("test.indexed_2d_vt", force = TRUE)
})

test_that("2D intervals work with giterator.intervals", {
    gdb.init_examples()
    gdir.create("test", showWarnings = FALSE)

    # Create intervals
    intervs <- gintervals.2d(
        chroms1 = rep(c(1, 2), each = 6),
        starts1 = rep(seq(100, 600, 100), 2),
        ends1 = rep(seq(150, 650, 100), 2),
        chroms2 = rep(c(1, 2), each = 6),
        starts2 = rep(seq(200, 700, 100), 2),
        ends2 = rep(seq(250, 750, 100), 2)
    )

    # Save in both formats
    withr::with_options(list(gmax.data.size = 10, gmulticontig.indexed_format = FALSE), {
        gintervals.save("test.legacy_2d_giter", intervs)
        gintervals.save("test.indexed_2d_giter", intervs)
    })

    gintervals.2d.convert_to_indexed("test.indexed_2d_giter", remove.old = TRUE, force = TRUE)

    # Use with giterator.intervals
    iter_legacy <- giterator.intervals(iterator = "test.legacy_2d_giter")
    iter_indexed <- giterator.intervals(iterator = "test.indexed_2d_giter")

    # Compare
    rownames(iter_legacy) <- NULL
    rownames(iter_indexed) <- NULL
    expect_equal(iter_legacy, iter_indexed, info = "giterator.intervals produces identical results")
    expect_equal(nrow(iter_legacy), nrow(intervs))

    # Cleanup
    gintervals.rm("test.legacy_2d_giter", force = TRUE)
    gintervals.rm("test.indexed_2d_giter", force = TRUE)
})

test_that("2D intervals work with complex gextract expressions", {
    gdb.init_examples()
    gdir.create("test", showWarnings = FALSE)

    # Create intervals
    intervs <- gintervals.2d(
        chroms1 = rep(c(1, 2), each = 6),
        starts1 = rep(seq(100, 600, 100), 2),
        ends1 = rep(seq(150, 650, 100), 2),
        chroms2 = rep(c(1, 2), each = 6),
        starts2 = rep(seq(200, 700, 100), 2),
        ends2 = rep(seq(250, 750, 100), 2)
    )

    # Create tracks
    values <- runif(nrow(intervs))
    gtrack.2d.create("test.track_2d_expr", "Test Track", intervs, values)

    # Save intervals in both formats
    withr::with_options(list(gmax.data.size = 10, gmulticontig.indexed_format = FALSE), {
        gintervals.save("test.legacy_2d_expr", intervs)
        gintervals.save("test.indexed_2d_expr", intervs)
    })

    gintervals.2d.convert_to_indexed("test.indexed_2d_expr", remove.old = TRUE, force = TRUE)

    # Complex expression with both formats
    expr <- "test.track_2d_expr * 2 + 5"
    result_legacy <- gextract(expr, intervals = "test.legacy_2d_expr")
    result_indexed <- gextract(expr, intervals = "test.indexed_2d_expr")

    # Compare
    rownames(result_legacy) <- NULL
    rownames(result_indexed) <- NULL
    expect_equal(result_legacy, result_indexed, info = "Complex expression produces identical results")

    # Cleanup
    gtrack.rm("test.track_2d_expr", force = TRUE)
    gintervals.rm("test.legacy_2d_expr", force = TRUE)
    gintervals.rm("test.indexed_2d_expr", force = TRUE)
})

test_that("2D intervals work with colnames parameter in gextract", {
    gdb.init_examples()
    gdir.create("test", showWarnings = FALSE)

    # Create intervals
    intervs <- gintervals.2d(
        chroms1 = rep(c(1, 2), each = 6),
        starts1 = rep(seq(100, 600, 100), 2),
        ends1 = rep(seq(150, 650, 100), 2),
        chroms2 = rep(c(1, 2), each = 6),
        starts2 = rep(seq(200, 700, 100), 2),
        ends2 = rep(seq(250, 750, 100), 2)
    )

    # Create track
    values <- runif(nrow(intervs))
    gtrack.2d.create("test.track_2d_colnames", "Test Track", intervs, values)

    # Save in both formats
    withr::with_options(list(gmax.data.size = 10, gmulticontig.indexed_format = FALSE), {
        gintervals.save("test.legacy_2d_colnames", intervs)
        gintervals.save("test.indexed_2d_colnames", intervs)
    })

    gintervals.2d.convert_to_indexed("test.indexed_2d_colnames", remove.old = TRUE, force = TRUE)

    # Extract with custom column names
    result_legacy <- gextract("test.track_2d_colnames",
        intervals = "test.legacy_2d_colnames",
        colnames = "my_value"
    )
    result_indexed <- gextract("test.track_2d_colnames",
        intervals = "test.indexed_2d_colnames",
        colnames = "my_value"
    )

    # Compare
    rownames(result_legacy) <- NULL
    rownames(result_indexed) <- NULL
    expect_equal(result_legacy, result_indexed, info = "Custom colnames work identically")
    expect_true("my_value" %in% colnames(result_legacy))

    # Cleanup
    gtrack.rm("test.track_2d_colnames", force = TRUE)
    gintervals.rm("test.legacy_2d_colnames", force = TRUE)
    gintervals.rm("test.indexed_2d_colnames", force = TRUE)
})

test_that("2D intervals work with band parameter in gextract", {
    gdb.init_examples()
    gdir.create("test", showWarnings = FALSE)

    # Create intervals
    intervs <- gintervals.2d(
        chroms1 = rep(c(1, 2), each = 6),
        starts1 = rep(seq(100, 600, 100), 2),
        ends1 = rep(seq(150, 650, 100), 2),
        chroms2 = rep(c(1, 2), each = 6),
        starts2 = rep(seq(200, 700, 100), 2),
        ends2 = rep(seq(250, 750, 100), 2)
    )

    # Create track
    values <- runif(nrow(intervs))
    gtrack.2d.create("test.track_2d_band", "Test Track", intervs, values)

    # Save in both formats
    withr::with_options(list(gmax.data.size = 10, gmulticontig.indexed_format = FALSE), {
        gintervals.save("test.legacy_2d_band", intervs)
        gintervals.save("test.indexed_2d_band", intervs)
    })

    gintervals.2d.convert_to_indexed("test.indexed_2d_band", remove.old = TRUE, force = TRUE)

    # Extract with band parameter
    result_legacy <- gextract("test.track_2d_band",
        intervals = "test.legacy_2d_band",
        band = c(-50, 50)
    )
    result_indexed <- gextract("test.track_2d_band",
        intervals = "test.indexed_2d_band",
        band = c(-50, 50)
    )

    # Compare
    rownames(result_legacy) <- NULL
    rownames(result_indexed) <- NULL
    expect_equal(result_legacy, result_indexed, info = "Band parameter works identically")

    # Cleanup
    gtrack.rm("test.track_2d_band", force = TRUE)
    gintervals.rm("test.legacy_2d_band", force = TRUE)
    gintervals.rm("test.indexed_2d_band", force = TRUE)
})

test_that("2D intervals work with gpartition", {
    gdb.init_examples()
    gdir.create("test", showWarnings = FALSE)

    # Create intervals
    intervs <- gintervals.2d(
        chroms1 = rep(c(1, 2), each = 6),
        starts1 = rep(seq(100, 600, 100), 2),
        ends1 = rep(seq(150, 650, 100), 2),
        chroms2 = rep(c(1, 2), each = 6),
        starts2 = rep(seq(200, 700, 100), 2),
        ends2 = rep(seq(250, 750, 100), 2)
    )

    # Create track
    values <- seq(0.1, 1.2, length.out = nrow(intervs))
    gtrack.2d.create("test.track_2d_part", "Test Track", intervs, values)

    # Save in both formats
    withr::with_options(list(gmax.data.size = 10, gmulticontig.indexed_format = FALSE), {
        gintervals.save("test.legacy_2d_part", intervs)
        gintervals.save("test.indexed_2d_part", intervs)
    })

    gintervals.2d.convert_to_indexed("test.indexed_2d_part", remove.old = TRUE, force = TRUE)

    # Use with gpartition
    breaks <- c(0, 0.5, 1.0, 1.5)
    result_legacy <- gpartition("test.track_2d_part",
        breaks = breaks,
        intervals = "test.legacy_2d_part"
    )
    result_indexed <- gpartition("test.track_2d_part",
        breaks = breaks,
        intervals = "test.indexed_2d_part"
    )

    # Compare
    expect_equal(result_legacy, result_indexed, info = "gpartition produces identical results")

    # Cleanup
    gtrack.rm("test.track_2d_part", force = TRUE)
    gintervals.rm("test.legacy_2d_part", force = TRUE)
    gintervals.rm("test.indexed_2d_part", force = TRUE)
})

test_that("2D intervals work with multiple tracks in gextract", {
    gdb.init_examples()
    gdir.create("test", showWarnings = FALSE)

    # Create intervals
    intervs <- gintervals.2d(
        chroms1 = rep(c(1, 2), each = 6),
        starts1 = rep(seq(100, 600, 100), 2),
        ends1 = rep(seq(150, 650, 100), 2),
        chroms2 = rep(c(1, 2), each = 6),
        starts2 = rep(seq(200, 700, 100), 2),
        ends2 = rep(seq(250, 750, 100), 2)
    )
    values1 <- runif(nrow(intervs))
    values2 <- runif(nrow(intervs))
    values1 <- runif(nrow(intervs))
    values2 <- runif(nrow(intervs))
    values1 <- seq(1, 12, 1)
    gtrack.2d.create("test.track_2d_multi1", "Track 1", intervs, values1)

    values2 <- seq(2, 24, 2)
    gtrack.2d.create("test.track_2d_multi2", "Track 2", intervs, values2)

    # Save in both formats
    withr::with_options(list(gmax.data.size = 10, gmulticontig.indexed_format = FALSE), {
        gintervals.save("test.legacy_2d_multi", intervs)
        gintervals.save("test.indexed_2d_multi", intervs)
    })

    gintervals.2d.convert_to_indexed("test.indexed_2d_multi", remove.old = TRUE, force = TRUE)

    # Extract multiple tracks - need explicit iterator for multiple 2D tracks
    result_legacy <- gextract(c("test.track_2d_multi1", "test.track_2d_multi2"),
        intervals = "test.legacy_2d_multi",
        iterator = "test.legacy_2d_multi"
    )
    result_indexed <- gextract(c("test.track_2d_multi1", "test.track_2d_multi2"),
        intervals = "test.indexed_2d_multi",
        iterator = "test.indexed_2d_multi"
    )

    # Compare
    rownames(result_legacy) <- NULL
    rownames(result_indexed) <- NULL
    expect_equal(result_legacy, result_indexed, info = "Multiple track extraction identical")
    expect_true(all(c("test.track_2d_multi1", "test.track_2d_multi2") %in% colnames(result_legacy)))

    # Cleanup
    gtrack.rm("test.track_2d_multi1", force = TRUE)
    gtrack.rm("test.track_2d_multi2", force = TRUE)
    gintervals.rm("test.legacy_2d_multi", force = TRUE)
    gintervals.rm("test.indexed_2d_multi", force = TRUE)
})
# Shaman-style tests: operations used in real HiC analysis

test_that("2D intervals with gintervals.2d.band_intersect (shaman-style)", {
    gdb.init_examples()
    gdir.create("test", showWarnings = FALSE)

    # Create 2D intervals similar to HiC contact matrix
    intervs <- gintervals.2d(
        chroms1 = rep(c(1, 2), each = 6),
        starts1 = rep(seq(100, 600, 100), 2),
        ends1 = rep(seq(150, 650, 100), 2),
        chroms2 = rep(c(1, 2), each = 6),
        starts2 = rep(seq(200, 700, 100), 2),
        ends2 = rep(seq(250, 750, 100), 2)
    )

    # Save in both formats
    withr::with_options(list(gmax.data.size = 10, gmulticontig.indexed_format = FALSE), {
        gintervals.save("test.legacy_2d_band", intervs)
        gintervals.save("test.indexed_2d_band", intervs)
    })

    gintervals.2d.convert_to_indexed("test.indexed_2d_band", remove.old = TRUE, force = TRUE)

    # Apply band intersect (filter by distance range - key operation in shaman)
    # Band represents x-y distance. For these intervals, start2 > end1, so x-y is negative
    # The intervals have x-y ranging from approximately -150 to -50
    # Using a wider negative band to ensure all intervals are covered
    band <- c(-300, 0)
    result_legacy <- gintervals.2d.band_intersect("test.legacy_2d_band", band = band)
    result_indexed <- gintervals.2d.band_intersect("test.indexed_2d_band", band = band)

    # Compare
    rownames(result_legacy) <- NULL
    rownames(result_indexed) <- NULL
    expect_equal(result_legacy, result_indexed, info = "band_intersect produces identical results")
    expect_true(!is.null(result_legacy) && nrow(result_legacy) > 0, info = "band_intersect returns non-empty results")
    expect_true(nrow(result_legacy) <= nrow(intervs), info = "band_intersect filters or keeps all intervals")

    # Cleanup
    gintervals.rm("test.legacy_2d_band", force = TRUE)
    gintervals.rm("test.indexed_2d_band", force = TRUE)
})

test_that("2D intervals with negative band range (shaman cis-interactions style)", {
    gdb.init_examples()
    gdir.create("test", showWarnings = FALSE)

    # Create intervals
    intervs <- gintervals.2d(
        chroms1 = rep(c(1, 2), each = 6),
        starts1 = rep(seq(1000, 6000, 1000), 2),
        ends1 = rep(seq(1050, 6050, 1000), 2),
        chroms2 = rep(c(1, 2), each = 6),
        starts2 = rep(seq(2000, 7000, 1000), 2),
        ends2 = rep(seq(2050, 7050, 1000), 2)
    )

    # Save in both formats
    withr::with_options(list(gmax.data.size = 10, gmulticontig.indexed_format = FALSE), {
        gintervals.save("test.legacy_2d_negband", intervs)
        gintervals.save("test.indexed_2d_negband", intervs)
    })

    gintervals.2d.convert_to_indexed("test.indexed_2d_negband", remove.old = TRUE, force = TRUE)

    # Negative band (upstream contacts) - common in HiC analysis
    band <- c(-5000, -500)
    result_legacy <- gintervals.2d.band_intersect("test.legacy_2d_negband", band = band)
    result_indexed <- gintervals.2d.band_intersect("test.indexed_2d_negband", band = band)

    # Compare
    rownames(result_legacy) <- NULL
    rownames(result_indexed) <- NULL
    expect_equal(result_legacy, result_indexed, info = "Negative band filtering identical")

    # Cleanup
    gintervals.rm("test.legacy_2d_negband", force = TRUE)
    gintervals.rm("test.indexed_2d_negband", force = TRUE)
})

test_that("2D intervals with expand.grid pattern (shaman regional matrix style)", {
    gdb.init_examples()
    gdir.create("test", showWarnings = FALSE)

    # Simulate shaman's expand.grid pattern for creating all-vs-all 2D matrix
    bins_1d <- gintervals(
        c(1, 1, 1, 2, 2, 2),
        c(100, 200, 300, 100, 200, 300),
        c(150, 250, 350, 150, 250, 350)
    )

    # Create all pairwise combinations (like shaman does for regional analysis)
    g <- expand.grid(1:nrow(bins_1d), 1:nrow(bins_1d))
    g <- g[bins_1d$chrom[g$Var1] == bins_1d$chrom[g$Var2], ] # Same chromosome only

    intervs <- gintervals.2d(
        bins_1d$chrom[g$Var1], bins_1d$start[g$Var1], bins_1d$end[g$Var1],
        bins_1d$chrom[g$Var2], bins_1d$start[g$Var2], bins_1d$end[g$Var2]
    )

    # Save in both formats
    withr::with_options(list(gmax.data.size = 10, gmulticontig.indexed_format = FALSE), {
        gintervals.save("test.legacy_2d_expand", intervs)
        gintervals.save("test.indexed_2d_expand", intervs)
    })

    gintervals.2d.convert_to_indexed("test.indexed_2d_expand", remove.old = TRUE, force = TRUE)

    # Load and compare
    loaded_legacy <- withr::with_options(list(gmax.data.size = 1e9), {
        gintervals.load("test.legacy_2d_expand")
    })
    loaded_indexed <- withr::with_options(list(gmax.data.size = 1e9), {
        gintervals.load("test.indexed_2d_expand")
    })

    rownames(loaded_legacy) <- NULL
    rownames(loaded_indexed) <- NULL
    expect_equal(loaded_legacy, loaded_indexed, info = "Expand.grid pattern identical")

    # Test upper triangle filtering (another shaman pattern)
    band <- c(1, 1e9) # Upper triangle only
    upper_legacy <- gintervals.2d.band_intersect("test.legacy_2d_expand", band = band)
    upper_indexed <- gintervals.2d.band_intersect("test.indexed_2d_expand", band = band)

    rownames(upper_legacy) <- NULL
    rownames(upper_indexed) <- NULL
    expect_equal(upper_legacy, upper_indexed, info = "Upper triangle filtering identical")

    # Cleanup
    gintervals.rm("test.legacy_2d_expand", force = TRUE)
    gintervals.rm("test.indexed_2d_expand", force = TRUE)
})

test_that("2D intervals with per-chromosome loading (shaman per-chrom processing)", {
    gdb.init_examples()
    gdir.create("test", showWarnings = FALSE)

    # Create multi-chromosome 2D intervals
    intervs <- gintervals.2d(
        chroms1 = c(rep(1, 6), rep(2, 6)),
        starts1 = rep(seq(100, 600, 100), 2),
        ends1 = rep(seq(150, 650, 100), 2),
        chroms2 = c(rep(1, 6), rep(2, 6)),
        starts2 = rep(seq(200, 700, 100), 2),
        ends2 = rep(seq(250, 750, 100), 2)
    )

    # Save in both formats
    withr::with_options(list(gmax.data.size = 10, gmulticontig.indexed_format = FALSE), {
        gintervals.save("test.legacy_2d_perchrom", intervs)
        gintervals.save("test.indexed_2d_perchrom", intervs)
    })

    gintervals.2d.convert_to_indexed("test.indexed_2d_perchrom", remove.old = TRUE, force = TRUE)

    # Load each chromosome pair separately (shaman processes per chromosome)
    for (chrom in c(1, 2)) {
        loaded_legacy <- gintervals.load("test.legacy_2d_perchrom", chrom1 = chrom, chrom2 = chrom)
        loaded_indexed <- gintervals.load("test.indexed_2d_perchrom", chrom1 = chrom, chrom2 = chrom)

        rownames(loaded_legacy) <- NULL
        rownames(loaded_indexed) <- NULL
        expect_equal(loaded_legacy, loaded_indexed,
            info = sprintf("Per-chromosome loading identical for chr%d", chrom)
        )
        expect_equal(nrow(loaded_legacy), 6)
    }

    # Cleanup
    gintervals.rm("test.legacy_2d_perchrom", force = TRUE)
    gintervals.rm("test.indexed_2d_perchrom", force = TRUE)
})

test_that("2D intervals with distance calculations (shaman min_dist filtering)", {
    gdb.init_examples()
    gdir.create("test", showWarnings = FALSE)

    # Create intervals at various distances
    intervs <- gintervals.2d(
        chroms1 = rep(1, 12),
        starts1 = rep(seq(1000, 4000, 1000), 3),
        ends1 = rep(seq(1100, 4100, 1000), 3),
        chroms2 = rep(1, 12),
        starts2 = c(
            seq(1500, 4500, 1000), # Close contacts
            seq(3000, 6000, 1000), # Medium contacts
            seq(10000, 13000, 1000)
        ), # Far contacts
        ends2 = c(
            seq(1600, 4600, 1000),
            seq(3100, 6100, 1000),
            seq(10100, 13100, 1000)
        )
    )

    # Save in both formats
    withr::with_options(list(gmax.data.size = 10, gmulticontig.indexed_format = FALSE), {
        gintervals.save("test.legacy_2d_dist", intervs)
        gintervals.save("test.indexed_2d_dist", intervs)
    })

    gintervals.2d.convert_to_indexed("test.indexed_2d_dist", remove.old = TRUE, force = TRUE)

    # Filter by minimum distance (common shaman parameter)
    min_dist <- 1500
    max_dist <- 5000
    band <- c(min_dist, max_dist)

    result_legacy <- gintervals.2d.band_intersect("test.legacy_2d_dist", band = band)
    result_indexed <- gintervals.2d.band_intersect("test.indexed_2d_dist", band = band)

    rownames(result_legacy) <- NULL
    rownames(result_indexed) <- NULL
    expect_equal(result_legacy, result_indexed, info = "Distance filtering identical")

    # Verify filtering worked correctly
    expect_true(all(result_legacy$start2 - result_legacy$start1 >= min_dist))
    expect_true(all(result_legacy$start2 - result_legacy$start1 <= max_dist))

    # Cleanup
    gintervals.rm("test.legacy_2d_dist", force = TRUE)
    gintervals.rm("test.indexed_2d_dist", force = TRUE)
})

test_that("2D intervals large-scale operations (shaman-scale HiC data)", {
    gdb.init_examples()
    gdir.create("test", showWarnings = FALSE)

    # Create larger dataset similar to HiC bin size
    # Simulating 50 bins on each of 2 chromosomes = 2500 contacts per chr
    n_bins <- 50
    bins <- seq(1000, 50000, 1000)

    # Chr1 contacts
    g1 <- expand.grid(1:n_bins, 1:n_bins)
    intervs1 <- gintervals.2d(
        rep(1, nrow(g1)), bins[g1$Var1], bins[g1$Var1] + 100,
        rep(1, nrow(g1)), bins[g1$Var2], bins[g1$Var2] + 100
    )

    # Chr2 contacts
    intervs2 <- gintervals.2d(
        rep(2, nrow(g1)), bins[g1$Var1], bins[g1$Var1] + 100,
        rep(2, nrow(g1)), bins[g1$Var2], bins[g1$Var2] + 100
    )

    intervs <- rbind(intervs1, intervs2)

    # Save in both formats
    withr::with_options(list(gmax.data.size = 10, gmulticontig.indexed_format = FALSE), {
        gintervals.save("test.legacy_2d_large", intervs)
        gintervals.save("test.indexed_2d_large", intervs)
    })

    gintervals.2d.convert_to_indexed("test.indexed_2d_large", remove.old = TRUE, force = TRUE)

    # Test per-chromosome loading (common in shaman for large datasets)
    for (chrom in c(1, 2)) {
        loaded_legacy <- gintervals.load("test.legacy_2d_large", chrom1 = chrom, chrom2 = chrom)
        loaded_indexed <- gintervals.load("test.indexed_2d_large", chrom1 = chrom, chrom2 = chrom)

        rownames(loaded_legacy) <- NULL
        rownames(loaded_indexed) <- NULL
        expect_equal(nrow(loaded_legacy), n_bins * n_bins)
        expect_equal(loaded_legacy, loaded_indexed,
            info = sprintf("Large-scale chr%d identical", chrom)
        )
    }

    # Test band filtering on large dataset
    band <- c(1000, 10000)
    band_legacy <- gintervals.2d.band_intersect("test.legacy_2d_large", band = band)
    band_indexed <- gintervals.2d.band_intersect("test.indexed_2d_large", band = band)

    rownames(band_legacy) <- NULL
    rownames(band_indexed) <- NULL
    expect_equal(band_legacy, band_indexed, info = "Large-scale band filtering identical")

    # Cleanup
    gintervals.rm("test.legacy_2d_large", force = TRUE)
    gintervals.rm("test.indexed_2d_large", force = TRUE)
})
