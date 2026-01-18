# Tests for multi-database support using the dataset API
# These tests complement test-dataset.R by testing more advanced scenarios
# Note: create_test_db() helper is defined in helper-test_db.R

# ==============================================================================
# Migration tests: verify old multi-db patterns migrate to new API
# ==============================================================================

test_that("gsetroot() with multiple paths errors with helpful message", {
    withr::with_tempdir({
        create_test_db("db1")
        create_test_db("db2")

        expect_error(
            gsetroot(c("db1", "db2")),
            "gsetroot\\(\\) accepts a single database path.*Use gdataset\\.load\\(\\) to load additional datasets"
        )
    })
})

test_that("multi-db pattern: gsetroot(db1) + gdataset.load(db2)", {
    withr::with_tempdir({
        create_test_db("db1")
        gsetroot("db1")
        intervs <- gintervals(1, 0, 1000)
        gtrack.create_sparse("track1", "from db1", intervs, 1)

        create_test_db("db2")
        gsetroot("db2")
        gtrack.create_sparse("track2", "from db2", intervs, 2)

        # New pattern: working db + loaded dataset
        gsetroot("db1")
        gdataset.load("db2")

        expect_equal(length(gdataset.ls()), 2)
        expect_true("track1" %in% gtrack.ls())
        expect_true("track2" %in% gtrack.ls())
    })
})

# ==============================================================================
# Track operations across databases
# ==============================================================================

test_that("gtrack.info works for tracks in different databases", {
    withr::with_tempdir({
        create_test_db("working_db")
        gsetroot("working_db")
        intervs <- gintervals(1, 0, 1000)
        gtrack.create_sparse("working_track", "sparse track", intervs, 1)

        create_test_db("dataset1")
        gsetroot("dataset1")
        gtrack.create_sparse("ds_track", "sparse track 2", intervs, 2)

        gsetroot("working_db")
        gdataset.load("dataset1")

        # Info should work for both
        info1 <- gtrack.info("working_track")
        info2 <- gtrack.info("ds_track")

        expect_equal(info1$type, "sparse")
        expect_equal(info2$type, "sparse")
    })
})

test_that("gtrack.exists works for tracks in loaded datasets", {
    withr::with_tempdir({
        create_test_db("working_db")
        gsetroot("working_db")
        intervs <- gintervals(1, 0, 1000)
        gtrack.create_sparse("working_track", "working", intervs, 1)

        create_test_db("dataset1")
        gsetroot("dataset1")
        gtrack.create_sparse("ds_track", "dataset", intervs, 2)

        gsetroot("working_db")
        gdataset.load("dataset1")

        expect_true(gtrack.exists("working_track"))
        expect_true(gtrack.exists("ds_track"))
        expect_false(gtrack.exists("nonexistent"))
    })
})

test_that("gtrack.rm only affects correct database", {
    withr::with_tempdir({
        create_test_db("working_db")
        gsetroot("working_db")
        intervs <- gintervals(1, 0, 1000)
        gtrack.create_sparse("track_to_keep", "keep", intervs, 1)
        gtrack.create_sparse("track_to_delete", "delete", intervs, 2)

        # Delete track from working db
        gtrack.rm("track_to_delete", force = TRUE)

        # Verify it's gone
        expect_false(gtrack.exists("track_to_delete"))

        # Track to keep should still exist
        expect_true(gtrack.exists("track_to_keep"))

        # Verify persistence after reload
        gdb.reload()
        expect_false(gtrack.exists("track_to_delete"))
        expect_true(gtrack.exists("track_to_keep"))
    })
})

test_that("gtrack.mv works in single db context", {
    withr::with_tempdir({
        create_test_db("working_db")
        gsetroot("working_db")
        intervs <- gintervals(1, 0, 1000)
        gtrack.create_sparse("original", "original track", intervs, 42)

        # Rename track
        gtrack.mv("original", "renamed")

        expect_false(gtrack.exists("original"))
        expect_true(gtrack.exists("renamed"))
        expect_equal(gtrack.dataset("renamed"), normalizePath("working_db"))

        # Value should be preserved
        expect_equal(gextract("renamed", gintervals(1, 0, 500))$renamed[1], 42)
    })
})

test_that("gtrack.copy works in single db context", {
    withr::with_tempdir({
        create_test_db("working_db")
        gsetroot("working_db")
        intervs <- gintervals(1, 0, 1000)
        gtrack.create_sparse("source", "src", intervs, 123)

        gtrack.copy("source", "dest")

        expect_true(gtrack.exists("source"))
        expect_true(gtrack.exists("dest"))

        # Values should match
        expect_equal(
            gextract("source", gintervals(1, 0, 500))$source[1],
            gextract("dest", gintervals(1, 0, 500))$dest[1]
        )
    })
})

# ==============================================================================
# Extract and compute operations across databases
# ==============================================================================

test_that("gextract works with tracks from working db and loaded datasets", {
    withr::with_tempdir({
        create_test_db("working_db")
        gsetroot("working_db")
        intervs <- gintervals(1, 0, 5000)
        gtrack.create_sparse("working_track", "working", intervs, 10)

        create_test_db("dataset1")
        gsetroot("dataset1")
        gtrack.create_sparse("ds_track", "dataset", intervs, 20)

        gsetroot("working_db")
        gdataset.load("dataset1")

        # Extract from both tracks in single call
        result <- gextract(c("working_track", "ds_track"), gintervals(1, 0, 1000), iterator = 100)

        expect_true("working_track" %in% names(result))
        expect_true("ds_track" %in% names(result))
        expect_equal(result$working_track[1], 10)
        expect_equal(result$ds_track[1], 20)
    })
})

test_that("track expressions work across databases", {
    withr::with_tempdir({
        create_test_db("working_db")
        gsetroot("working_db")
        intervs <- gintervals(1, 0, 5000)
        gtrack.create_sparse("x", "working", intervs, 10)

        create_test_db("dataset1")
        gsetroot("dataset1")
        gtrack.create_sparse("y", "dataset", intervs, 3)

        gsetroot("working_db")
        gdataset.load("dataset1")

        # Expression using tracks from different databases
        result <- gextract("x + y", gintervals(1, 0, 1000), iterator = 100)
        expect_equal(result[["x + y"]][1], 13)
    })
})

test_that("gscreen works with tracks from different databases", {
    withr::with_tempdir({
        create_test_db("working_db")
        gsetroot("working_db")
        intervs <- gintervals(1, 0, 5000)
        gtrack.create_sparse("working_track", "working", intervs, 10)

        create_test_db("dataset1")
        gsetroot("dataset1")
        gtrack.create_sparse("ds_track", "dataset", intervs, 5)

        gsetroot("working_db")
        gdataset.load("dataset1")

        # Screen using track from working db
        result1 <- gscreen("working_track > 5", gintervals(1, 0, 1000))
        expect_true(nrow(result1) > 0)

        # Screen using track from dataset
        result2 <- gscreen("ds_track > 3", gintervals(1, 0, 1000))
        expect_true(nrow(result2) > 0)

        # Screen using expression combining both
        result3 <- gscreen("working_track + ds_track > 10", gintervals(1, 0, 1000), iterator = 100)
        expect_true(nrow(result3) > 0)
    })
})

test_that("gsummary works with tracks from different databases", {
    withr::with_tempdir({
        create_test_db("working_db")
        gsetroot("working_db")
        intervs <- gintervals(1, 0, 5000)
        gtrack.create_sparse("working_track", "working", intervs, 10)

        create_test_db("dataset1")
        gsetroot("dataset1")
        gtrack.create_sparse("ds_track", "dataset", intervs, 5)

        gsetroot("working_db")
        gdataset.load("dataset1")

        # Summary of track from working db
        sum1 <- gsummary("working_track", gintervals(1, 0, 1000))
        expect_equal(sum1[["Sum"]], 10)

        # Summary of track from dataset
        sum2 <- gsummary("ds_track", gintervals(1, 0, 1000))
        expect_equal(sum2[["Sum"]], 5)
    })
})

# ==============================================================================
# Virtual tracks across databases
# ==============================================================================

test_that("virtual tracks can reference tracks from any database", {
    withr::with_tempdir({
        create_test_db("working_db")
        gsetroot("working_db")
        intervs <- gintervals(1, 0, 5000)
        gtrack.create_sparse("working_track", "working", intervs, 10)

        create_test_db("dataset1")
        gsetroot("dataset1")
        gtrack.create_sparse("ds_track", "dataset", intervs, 20)

        gsetroot("working_db")
        gdataset.load("dataset1")

        gvtrack.create("vt_working", "working_track", "avg")
        gvtrack.create("vt_ds", "ds_track", "avg")
        withr::defer({
            gvtrack.rm("vt_working")
            gvtrack.rm("vt_ds")
        })

        result1 <- gextract("vt_working", gintervals(1, 0, 1000))
        result2 <- gextract("vt_ds", gintervals(1, 0, 1000))

        expect_equal(result1$vt_working[1], 10)
        expect_equal(result2$vt_ds[1], 20)
    })
})

test_that("virtual tracks remain global across database operations", {
    withr::with_tempdir({
        create_test_db("working_db")
        gsetroot("working_db")
        intervs <- gintervals(1, 0, 1000)
        gtrack.create_sparse("base_track", "base", intervs, 1)

        create_test_db("dataset1")

        # Create virtual track
        gvtrack.create("vtrack", "base_track", "avg")
        withr::defer(gvtrack.rm("vtrack"))

        # Load dataset
        gdataset.load("dataset1")

        # Virtual track should still be accessible
        expect_true("vtrack" %in% gvtrack.ls())
    })
})

# ==============================================================================
# Interval operations across databases
# ==============================================================================

test_that("gintervals operations work across databases", {
    withr::with_tempdir({
        create_test_db("working_db")
        gsetroot("working_db")
        int1 <- gintervals(1, 0, 1000)
        gintervals.save("working_intervals", int1)

        create_test_db("dataset1")
        gsetroot("dataset1")
        int2 <- gintervals(1, 500, 1500)
        gintervals.save("ds_intervals", int2)

        gsetroot("working_db")
        gdataset.load("dataset1")

        # Both interval sets should be visible
        all_intervals <- gintervals.ls()
        expect_true("working_intervals" %in% all_intervals)
        expect_true("ds_intervals" %in% all_intervals)

        # Verify the intervals exist and can be looked up
        expect_equal(gintervals.dataset("working_intervals"), normalizePath("working_db"))
        expect_equal(gintervals.dataset("ds_intervals"), normalizePath("dataset1"))
    })
})

# ==============================================================================
# Database reload behavior
# ==============================================================================

test_that("gdb.reload correctly refreshes state with loaded datasets", {
    withr::with_tempdir({
        create_test_db("working_db")
        gsetroot("working_db")
        intervs <- gintervals(1, 0, 1000)
        gtrack.create_sparse("track1", "t1", intervs, 1)

        create_test_db("dataset1")
        gsetroot("dataset1")
        gtrack.create_sparse("ds_track1", "ds1", intervs, 2)

        gsetroot("working_db")
        gdataset.load("dataset1")

        expect_equal(length(gtrack.ls()), 2)

        # Reload should preserve loaded datasets
        gdb.reload(rescan = TRUE)

        expect_equal(length(gtrack.ls()), 2)
        expect_true("track1" %in% gtrack.ls())
        expect_true("ds_track1" %in% gtrack.ls())
    })
})

# ==============================================================================
# Track attributes across databases
# ==============================================================================

test_that("track attributes work across databases", {
    withr::with_tempdir({
        create_test_db("working_db")
        gsetroot("working_db")
        intervs <- gintervals(1, 0, 1000)
        gtrack.create_sparse("working_track", "track in working db", intervs, 1)
        gtrack.attr.set("working_track", "custom_attr", "working_value")

        create_test_db("dataset1")
        gsetroot("dataset1")
        gtrack.create_sparse("ds_track", "track in dataset", intervs, 2)
        gtrack.attr.set("ds_track", "custom_attr", "ds_value")

        gsetroot("working_db")
        gdataset.load("dataset1")

        # Should read attributes from both databases
        expect_equal(gtrack.attr.get("working_track", "custom_attr"), "working_value")
        expect_equal(gtrack.attr.get("ds_track", "custom_attr"), "ds_value")
    })
})

# ==============================================================================
# Multiple datasets
# ==============================================================================

test_that("three or more datasets can be loaded", {
    withr::with_tempdir({
        create_test_db("working_db")
        gsetroot("working_db")
        intervs <- gintervals(1, 0, 1000)
        gtrack.create_sparse("track_a", "from working db", intervs, 1)

        create_test_db("dataset1")
        gsetroot("dataset1")
        gtrack.create_sparse("track_b", "from ds1", intervs, 2)

        create_test_db("dataset2")
        gsetroot("dataset2")
        gtrack.create_sparse("track_c", "from ds2", intervs, 3)

        create_test_db("dataset3")
        gsetroot("dataset3")
        gtrack.create_sparse("track_d", "from ds3", intervs, 4)

        # Connect working db + three datasets
        gsetroot("working_db")
        gdataset.load("dataset1")
        gdataset.load("dataset2")
        gdataset.load("dataset3")

        expect_equal(length(gdataset.ls()), 4)
        expect_equal(length(gtrack.ls()), 4)

        # Verify each track resolves to correct db
        expect_equal(gtrack.dataset("track_a"), normalizePath("working_db"))
        expect_equal(gtrack.dataset("track_b"), normalizePath("dataset1"))
        expect_equal(gtrack.dataset("track_c"), normalizePath("dataset2"))
        expect_equal(gtrack.dataset("track_d"), normalizePath("dataset3"))
    })
})

test_that("empty dataset in chain works", {
    withr::with_tempdir({
        create_test_db("working_db")
        gsetroot("working_db")
        intervs <- gintervals(1, 0, 1000)
        gtrack.create_sparse("track1", "working", intervs, 1)

        create_test_db("empty_dataset")
        # empty_dataset has no tracks

        create_test_db("dataset2")
        gsetroot("dataset2")
        gtrack.create_sparse("track2", "ds2", intervs, 2)

        gsetroot("working_db")
        gdataset.load("empty_dataset")
        gdataset.load("dataset2")

        expect_equal(length(gdataset.ls()), 3)
        tracks <- gtrack.ls()
        expect_equal(length(tracks), 2)
        expect_true("track1" %in% tracks)
        expect_true("track2" %in% tracks)
    })
})

# ==============================================================================
# Collision handling
# ==============================================================================

test_that("deletion of overriding track in working db still shows working db", {
    withr::with_tempdir({
        create_test_db("working_db")
        gsetroot("working_db")
        intervs <- gintervals(1, 0, 1000)
        gtrack.create_sparse("shared", "from working db", intervs, 100)

        create_test_db("dataset1")
        gsetroot("dataset1")
        gtrack.create_sparse("shared", "from dataset", intervs, 200)

        gsetroot("working_db")
        gdataset.load("dataset1", force = TRUE)

        # Working db wins, so value should be 100
        expect_equal(gextract("shared", gintervals(1, 0, 500))$shared[1], 100)
        expect_equal(gtrack.dataset("shared"), normalizePath("working_db"))

        # Delete from working db
        gtrack.rm("shared", force = TRUE)
        gdb.reload()

        # Now dataset1 should be visible
        expect_true(gtrack.exists("shared"))
        expect_equal(gtrack.dataset("shared"), normalizePath("dataset1"))
        expect_equal(gextract("shared", gintervals(1, 0, 500))$shared[1], 200)
    })
})

test_that("gdataset.ls(dataframe=TRUE) tracks counts are accurate with overlapping tracks", {
    withr::with_tempdir({
        create_test_db("working_db")
        gsetroot("working_db")
        intervs <- gintervals(1, 0, 1000)
        gtrack.create_sparse("unique_working", "working", intervs, 1)
        gtrack.create_sparse("shared", "shared from working", intervs, 10)

        create_test_db("dataset1")
        gsetroot("dataset1")
        gtrack.create_sparse("unique_ds", "ds", intervs, 2)
        gtrack.create_sparse("shared", "shared from ds", intervs, 20)

        gsetroot("working_db")
        gdataset.load("dataset1", force = TRUE)

        info <- gdataset.ls(dataframe = TRUE)

        # Working db should have 2 visible tracks (unique_working + shared)
        # because working db wins for collision
        expect_equal(info$tracks_visible[1], 2)

        # Dataset1 should have 1 visible track (unique_ds only, shared is shadowed)
        expect_equal(info$tracks_visible[2], 1)

        # Total visible tracks should be 3
        expect_equal(length(gtrack.ls()), 3)
    })
})

# ==============================================================================
# Path handling
# ==============================================================================

test_that("absolute and relative paths work correctly", {
    withr::with_tempdir({
        create_test_db("working_db")
        gsetroot("working_db")
        intervs <- gintervals(1, 0, 1000)
        gtrack.create_sparse("track1", "t1", intervs, 1)

        create_test_db("dataset1")
        gsetroot("dataset1")
        gtrack.create_sparse("track2", "t2", intervs, 2)

        gsetroot("working_db")

        # Mix absolute and relative paths
        abs_path <- normalizePath("dataset1")
        gdataset.load(abs_path)

        expect_equal(length(gtrack.ls()), 2)
        expect_true(gtrack.exists("track1"))
        expect_true(gtrack.exists("track2"))
    })
})

# ==============================================================================
# Error handling
# ==============================================================================

test_that("gsetroot gives clear error when directory doesn't exist", {
    expect_error(
        gsetroot("/this/path/does/not/exist"),
        "Database directory does not exist"
    )
})

test_that("gsetroot gives clear error when tracks/ subdirectory is missing", {
    withr::with_tempdir({
        dir.create("not_a_db/seq", recursive = TRUE)
        write.table(data.frame(chrom = "chr1", size = 1000), "not_a_db/chrom_sizes.txt",
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE
        )

        expect_error(
            gsetroot("not_a_db"),
            "does not contain a 'tracks' subdirectory"
        )
    })
})

test_that("gsetroot gives clear error when seq/ subdirectory is missing", {
    withr::with_tempdir({
        dir.create("not_a_db/tracks", recursive = TRUE)
        write.table(data.frame(chrom = "chr1", size = 1000), "not_a_db/chrom_sizes.txt",
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE
        )

        expect_error(
            gsetroot("not_a_db"),
            "does not contain a 'seq' subdirectory"
        )
    })
})

test_that("gsetroot gives clear error when chrom_sizes.txt is missing", {
    withr::with_tempdir({
        dir.create("not_a_db/tracks", recursive = TRUE)
        dir.create("not_a_db/seq", recursive = TRUE)

        expect_error(
            gsetroot("not_a_db"),
            "does not contain a chrom_sizes.txt file"
        )
    })
})

test_that("gsetroot detects when user passes database path as dir parameter", {
    withr::with_tempdir({
        create_test_db("db1")
        create_test_db("db2")

        # Common mistake: passing second database as 'dir' parameter
        expect_error(
            gsetroot("db1", "db2"),
            "looks like a misha database path.*To connect multiple databases, use.*gdataset\\.load"
        )
    })
})
