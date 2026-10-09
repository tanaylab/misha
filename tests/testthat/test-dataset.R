# This file re-roots into throwaway databases inside withr::with_tempdir(),
# which leaves .misha$GROOT pointing at a directory that no longer exists.
# Put a working root back for whatever runs next in this parallel worker -
# without it, test-gdb-unload.R (and anything else that reads the ambient
# root) fails on state this file left behind.
restore_groot_on_exit()

# Tests for dataset API
# Note: create_test_db() helper is defined in helper-test_db.R

# ==============================================================================
# Phase 1: gsetroot() - Single Database Only
# ==============================================================================

test_that("gsetroot() errors on vector input", {
    withr::with_tempdir({
        create_test_db("db1")
        create_test_db("db2")

        expect_error(
            gsetroot(c("db1", "db2")),
            "gsetroot\\(\\) accepts a single database path.*Use gdataset\\.load\\(\\) to load additional datasets"
        )
    })
})

test_that("gsetroot() works with single database (backward compatibility)", {
    withr::with_tempdir({
        create_test_db("single_db")

        gsetroot("single_db")

        # GROOT should be set
        expect_equal(.misha$GROOT, normalizePath("single_db"))

        # GDATASETS should be empty
        expect_equal(.misha$GDATASETS, character(0))
    })
})

# ==============================================================================
# Phase 2: Basic Dataset Loading
# ==============================================================================

test_that("gdataset.load() loads tracks from dataset", {
    withr::with_tempdir({
        # Create working db
        create_test_db("working_db")
        gsetroot("working_db")

        # Create dataset with a track
        create_test_db("dataset1")
        gsetroot("dataset1")
        intervs <- gintervals(1, 0, 1000)
        gtrack.create_sparse("dataset_track", "track from dataset", intervs, 5)

        # Switch back to working db and load dataset
        gsetroot("working_db")
        result <- gdataset.load("dataset1")

        # Track should be visible
        expect_true("dataset_track" %in% gtrack.ls())

        # Return value should show counts
        expect_equal(result$tracks, 1)
        expect_equal(result$intervals, 0)
        expect_equal(result$shadowed_tracks, 0)
        expect_equal(result$shadowed_intervals, 0)
    })
})

test_that("gdataset.load() loads intervals from dataset", {
    withr::with_tempdir({
        # Create working db
        create_test_db("working_db")
        gsetroot("working_db")

        # Create dataset with intervals
        create_test_db("dataset1")
        gsetroot("dataset1")
        intervs <- gintervals(1, 0, 1000)
        gintervals.save("dataset_intervals", intervs)

        # Switch back and load
        gsetroot("working_db")
        result <- gdataset.load("dataset1")

        # Intervals should be visible
        expect_true("dataset_intervals" %in% gintervals.ls())

        expect_equal(result$tracks, 0)
        expect_equal(result$intervals, 1)
    })
})

test_that("gdataset.load() loads tracks and intervals", {
    withr::with_tempdir({
        create_test_db("working_db")
        gsetroot("working_db")

        create_test_db("dataset1")
        gsetroot("dataset1")
        intervs <- gintervals(1, 0, 1000)
        gtrack.create_sparse("ds_track", "track", intervs, 1)
        gintervals.save("ds_intervals", intervs)

        gsetroot("working_db")
        result <- gdataset.load("dataset1")

        expect_true("ds_track" %in% gtrack.ls())
        expect_true("ds_intervals" %in% gintervals.ls())
        expect_equal(result$tracks, 1)
        expect_equal(result$intervals, 1)
    })
})

test_that("gdataset.load() verbose mode prints information", {
    withr::with_tempdir({
        create_test_db("working_db")
        gsetroot("working_db")

        create_test_db("dataset1")
        gsetroot("dataset1")
        intervs <- gintervals(1, 0, 1000)
        gtrack.create_sparse("track1", "t1", intervs, 1)

        gsetroot("working_db")

        expect_message(gdataset.load("dataset1", verbose = TRUE), "Loaded dataset")
    })
})

test_that("gdataset.load() reload (idempotency via unload+load)", {
    withr::with_tempdir({
        create_test_db("working_db")
        gsetroot("working_db")

        create_test_db("dataset1")
        gsetroot("dataset1")
        intervs <- gintervals(1, 0, 1000)
        gtrack.create_sparse("track1", "t1", intervs, 1)

        gsetroot("working_db")
        gdataset.load("dataset1")

        # Load again - should unload first then reload
        expect_no_error(gdataset.load("dataset1"))
        expect_true("track1" %in% gtrack.ls())
    })
})

test_that("gdataset.load() normalizes paths correctly", {
    withr::with_tempdir({
        create_test_db("working_db")
        gsetroot("working_db")

        create_test_db("dataset1")
        gsetroot("dataset1")
        intervs <- gintervals(1, 0, 1000)
        gtrack.create_sparse("track1", "t1", intervs, 1)

        gsetroot("working_db")

        # Load with relative path
        gdataset.load("dataset1")

        # Load with absolute path should be same dataset
        abs_path <- normalizePath("dataset1")
        expect_no_error(gdataset.load(abs_path))
    })
})

test_that("gdataset.load() loading working db path is no-op", {
    withr::with_tempdir({
        create_test_db("working_db")
        gsetroot("working_db")

        # Try to load working db as dataset
        result <- gdataset.load("working_db")

        # Should silently return zeros
        expect_equal(result$tracks, 0)
        expect_equal(result$intervals, 0)
        expect_equal(result$shadowed_tracks, 0)
        expect_equal(result$shadowed_intervals, 0)
    })
})

test_that("gdataset.load() errors without gsetroot() called", {
    withr::with_tempdir({
        create_test_db("dataset1")

        # Clear misha state
        rm(list = ls(envir = .misha), envir = .misha)

        expect_error(
            gdataset.load("dataset1"),
            "No working database.*Call gsetroot\\(\\) first"
        )
    })
})

test_that("gdataset.load() errors when path doesn't exist", {
    withr::with_tempdir({
        create_test_db("working_db")
        gsetroot("working_db")

        expect_error(
            gdataset.load("nonexistent_path"),
            "does not exist"
        )
    })
})

test_that("gdataset.load() errors when path has no tracks/ directory", {
    withr::with_tempdir({
        create_test_db("working_db")
        gsetroot("working_db")

        dir.create("not_a_dataset")

        expect_error(
            gdataset.load("not_a_dataset"),
            "tracks.*directory"
        )
    })
})

test_that("gdataset.load() errors when chrom_sizes.txt is missing", {
    withr::with_tempdir({
        create_test_db("working_db")
        gsetroot("working_db")

        dir.create("dataset1/tracks", recursive = TRUE)

        expect_error(
            gdataset.load("dataset1"),
            "chrom_sizes\\.txt"
        )
    })
})

test_that("gdataset.load() errors when chrom_sizes.txt doesn't match working db", {
    withr::with_tempdir({
        create_test_db("working_db", chrom_sizes = data.frame(chrom = "chr1", size = 10000))
        gsetroot("working_db")

        # Create dataset with different genome
        create_test_db("dataset1", chrom_sizes = data.frame(chrom = "chr1", size = 20000))

        expect_error(
            gdataset.load("dataset1"),
            "genome.*match"
        )
    })
})

# ==============================================================================
# Collision Handling
# ==============================================================================

test_that("gdataset.load() detects collision with working db tracks", {
    withr::with_tempdir({
        create_test_db("working_db")
        gsetroot("working_db")
        intervs <- gintervals(1, 0, 1000)
        gtrack.create_sparse("shared_track", "working", intervs, 1)

        create_test_db("dataset1")
        gsetroot("dataset1")
        gtrack.create_sparse("shared_track", "dataset", intervs, 2)

        gsetroot("working_db")

        expect_error(
            gdataset.load("dataset1"),
            "Cannot load dataset.*tracks 'shared_track'.*already exist in working database"
        )
    })
})

test_that("gdataset.load() detects collision with working db intervals", {
    withr::with_tempdir({
        create_test_db("working_db")
        gsetroot("working_db")
        intervs <- gintervals(1, 0, 1000)
        gintervals.save("shared_intervals", intervs)

        create_test_db("dataset1")
        gsetroot("dataset1")
        gintervals.save("shared_intervals", intervs)

        gsetroot("working_db")

        expect_error(
            gdataset.load("dataset1"),
            "Cannot load dataset.*interval sets 'shared_intervals'.*already exist in working database"
        )
    })
})

test_that("gdataset.load() force=TRUE allows working db to win for tracks", {
    withr::with_tempdir({
        create_test_db("working_db")
        gsetroot("working_db")
        intervs <- gintervals(1, 0, 1000)
        gtrack.create_sparse("shared_track", "working", intervs, 100)

        create_test_db("dataset1")
        gsetroot("dataset1")
        gtrack.create_sparse("shared_track", "dataset", intervs, 200)

        gsetroot("working_db")

        # Should succeed with force=TRUE
        result <- gdataset.load("dataset1", force = TRUE)

        # Working db should win
        expect_equal(gextract("shared_track", gintervals(1, 0, 500))$shared_track[1], 100)

        # Dataset track should be shadowed
        expect_equal(result$shadowed_tracks, 1)
        expect_equal(result$tracks, 0) # No visible tracks from dataset
    })
})

test_that("gdataset.load() force=TRUE allows working db to win for intervals", {
    withr::with_tempdir({
        create_test_db("working_db")
        gsetroot("working_db")
        intervs1 <- gintervals(1, 0, 1000)
        gintervals.save("shared_intervals", intervs1)

        create_test_db("dataset1")
        gsetroot("dataset1")
        intervs2 <- gintervals(1, 0, 2000)
        gintervals.save("shared_intervals", intervs2)

        gsetroot("working_db")

        result <- gdataset.load("dataset1", force = TRUE)

        # Working db should win
        loaded <- gintervals.load("shared_intervals")
        expect_equal(loaded$end[1], 1000)

        expect_equal(result$shadowed_intervals, 1)
    })
})

test_that("gdataset.load() detects dataset-to-dataset collision", {
    withr::with_tempdir({
        create_test_db("working_db")
        gsetroot("working_db")

        create_test_db("dataset1")
        gsetroot("dataset1")
        intervs <- gintervals(1, 0, 1000)
        gtrack.create_sparse("shared", "ds1", intervs, 1)

        create_test_db("dataset2")
        gsetroot("dataset2")
        gtrack.create_sparse("shared", "ds2", intervs, 2)

        gsetroot("working_db")
        gdataset.load("dataset1")

        expect_error(
            gdataset.load("dataset2"),
            "Cannot load dataset.*tracks 'shared'.*already exist in loaded dataset"
        )
    })
})

test_that("gdataset.load() force=TRUE allows later dataset to win", {
    withr::with_tempdir({
        create_test_db("working_db")
        gsetroot("working_db")

        create_test_db("dataset1")
        gsetroot("dataset1")
        intervs <- gintervals(1, 0, 1000)
        gtrack.create_sparse("shared", "ds1", intervs, 10)

        create_test_db("dataset2")
        gsetroot("dataset2")
        gtrack.create_sparse("shared", "ds2", intervs, 20)

        gsetroot("working_db")
        gdataset.load("dataset1")

        result <- gdataset.load("dataset2", force = TRUE)

        # Dataset2 should win
        expect_equal(gextract("shared", gintervals(1, 0, 500))$shared[1], 20)

        expect_equal(result$shadowed_tracks, 0) # Dataset2 doesn't shadow itself
        expect_equal(result$tracks, 1)
    })
})

# ==============================================================================
# gdataset.unload()
# ==============================================================================

test_that("gdataset.unload() removes dataset tracks", {
    withr::with_tempdir({
        create_test_db("working_db")
        gsetroot("working_db")

        create_test_db("dataset1")
        gsetroot("dataset1")
        intervs <- gintervals(1, 0, 1000)
        gtrack.create_sparse("ds_track", "track", intervs, 1)

        gsetroot("working_db")
        gdataset.load("dataset1")
        expect_true("ds_track" %in% gtrack.ls())

        gdataset.unload("dataset1")
        expect_false("ds_track" %in% gtrack.ls())
    })
})

test_that("gdataset.unload() validates path normalization", {
    withr::with_tempdir({
        create_test_db("working_db")
        gsetroot("working_db")

        create_test_db("dataset1")
        gsetroot("dataset1")
        intervs <- gintervals(1, 0, 1000)
        gtrack.create_sparse("track1", "t1", intervs, 1)

        gsetroot("working_db")
        gdataset.load("dataset1")

        # Unload with absolute path
        abs_path <- normalizePath("dataset1")
        expect_no_error(gdataset.unload(abs_path))
    })
})

test_that("gdataset.unload() with validate=FALSE silently ignores non-loaded paths", {
    withr::with_tempdir({
        create_test_db("working_db")
        gsetroot("working_db")

        # Should not error
        expect_no_error(gdataset.unload("/nonexistent/path", validate = FALSE))
    })
})

test_that("gdataset.unload() with validate=TRUE errors when path not loaded", {
    withr::with_tempdir({
        create_test_db("working_db")
        gsetroot("working_db")

        expect_error(
            gdataset.unload("/nonexistent/path", validate = TRUE),
            "not loaded"
        )
    })
})

test_that("gdataset.unload() restores shadowed tracks from working db", {
    withr::with_tempdir({
        create_test_db("working_db")
        gsetroot("working_db")
        intervs <- gintervals(1, 0, 1000)
        gtrack.create_sparse("shared", "working", intervs, 100)

        create_test_db("dataset1")
        gsetroot("dataset1")
        gtrack.create_sparse("shared", "dataset", intervs, 200)

        gsetroot("working_db")
        gdataset.load("dataset1", force = TRUE)

        # Working db wins, so value should be 100
        expect_equal(gextract("shared", gintervals(1, 0, 500))$shared[1], 100)

        # Unload - working db track still visible
        gdataset.unload("dataset1")
        expect_equal(gextract("shared", gintervals(1, 0, 500))$shared[1], 100)
    })
})

test_that("gdataset.unload() restores shadowed tracks from other datasets (working db priority)", {
    withr::with_tempdir({
        create_test_db("working_db")
        gsetroot("working_db")
        intervs <- gintervals(1, 0, 1000)
        gtrack.create_sparse("shared", "working", intervs, 100)

        create_test_db("dataset1")
        gsetroot("dataset1")
        gtrack.create_sparse("unique1", "ds1", intervs, 1)

        create_test_db("dataset2")
        gsetroot("dataset2")
        gtrack.create_sparse("shared", "ds2", intervs, 200)

        gsetroot("working_db")
        gdataset.load("dataset1")
        gdataset.load("dataset2", force = TRUE)

        # Working db should win
        expect_equal(gextract("shared", gintervals(1, 0, 500))$shared[1], 100)

        # Unload all datasets
        gdataset.unload("dataset1")
        gdataset.unload("dataset2")

        # Working db track still there
        expect_equal(gextract("shared", gintervals(1, 0, 500))$shared[1], 100)
    })
})

test_that("gdataset.unload() restores shadowed tracks in load order", {
    withr::with_tempdir({
        create_test_db("working_db")
        gsetroot("working_db")

        create_test_db("dataset1")
        gsetroot("dataset1")
        intervs <- gintervals(1, 0, 1000)
        gtrack.create_sparse("shared", "ds1", intervs, 10)

        create_test_db("dataset2")
        gsetroot("dataset2")
        gtrack.create_sparse("shared", "ds2", intervs, 20)

        gsetroot("working_db")
        gdataset.load("dataset1")
        gdataset.load("dataset2", force = TRUE)

        # Dataset2 wins
        expect_equal(gextract("shared", gintervals(1, 0, 500))$shared[1], 20)

        # Unload dataset2, dataset1 should become visible
        gdataset.unload("dataset2")
        expect_equal(gextract("shared", gintervals(1, 0, 500))$shared[1], 10)
    })
})

# ==============================================================================
# gdataset.save()
# ==============================================================================

test_that("gdataset.save() creates dataset with tracks", {
    withr::with_tempdir({
        create_test_db("working_db")
        gsetroot("working_db")
        intervs <- gintervals(1, 0, 1000)
        gtrack.create_sparse("track1", "t1", intervs, 1)
        gtrack.create_sparse("track2", "t2", intervs, 2)

        gdataset.save(
            path = "my_dataset",
            description = "Test dataset",
            tracks = c("track1", "track2")
        )

        # Check directory structure
        expect_true(dir.exists("my_dataset"))
        expect_true(dir.exists("my_dataset/tracks"))
        expect_true(file.exists("my_dataset/chrom_sizes.txt"))
        expect_true(file.exists("my_dataset/seq") || Sys.readlink("my_dataset/seq") != "")
        expect_true(file.exists("my_dataset/misha.yaml"))

        # Check tracks exist
        expect_true(dir.exists("my_dataset/tracks/track1.track"))
        expect_true(dir.exists("my_dataset/tracks/track2.track"))
    })
})

test_that("gdataset.save() creates dataset with intervals", {
    withr::with_tempdir({
        create_test_db("working_db")
        gsetroot("working_db")
        intervs1 <- gintervals(1, 0, 1000)
        intervs2 <- gintervals(2, 0, 2000)
        gintervals.save("intervals1", intervs1)
        gintervals.save("intervals2", intervs2)

        gdataset.save(
            path = "my_dataset",
            description = "Intervals dataset",
            intervals = c("intervals1", "intervals2")
        )

        # Intervals are stored under tracks/ directory
        expect_true(file.exists("my_dataset/tracks/intervals1.interv") ||
            dir.exists("my_dataset/tracks/intervals1.interv"))
    })
})

test_that("gdataset.save() creates dataset with tracks and intervals", {
    withr::with_tempdir({
        create_test_db("working_db")
        gsetroot("working_db")
        intervs <- gintervals(1, 0, 1000)
        gtrack.create_sparse("track1", "t1", intervs, 1)
        gintervals.save("intervals1", intervs)

        gdataset.save(
            path = "my_dataset",
            description = "Mixed dataset",
            tracks = "track1",
            intervals = "intervals1"
        )

        expect_true(dir.exists("my_dataset/tracks/track1.track"))
        expect_true(file.exists("my_dataset/tracks/intervals1.interv") ||
            dir.exists("my_dataset/tracks/intervals1.interv"))
    })
})

test_that("gdataset.save() with symlinks=TRUE creates symlinks", {
    skip_on_os("windows") # Symlinks may not work on Windows

    withr::with_tempdir({
        create_test_db("working_db")
        gsetroot("working_db")
        intervs <- gintervals(1, 0, 1000)
        gtrack.create_sparse("track1", "t1", intervs, 1)

        gdataset.save(
            path = "my_dataset",
            description = "Symlinked dataset",
            tracks = "track1",
            symlinks = TRUE
        )

        # Check track is symlink
        track_path <- "my_dataset/tracks/track1.track"
        expect_true(Sys.readlink(track_path) != "")
    })
})

test_that("gdataset.save() with copy_seq=TRUE copies seq directory", {
    withr::with_tempdir({
        create_test_db("working_db")
        gsetroot("working_db")
        intervs <- gintervals(1, 0, 1000)
        gtrack.create_sparse("track1", "t1", intervs, 1)

        gdataset.save(
            path = "my_dataset",
            description = "Dataset with copied seq",
            tracks = "track1",
            copy_seq = TRUE
        )

        # Check seq is a real directory (not symlink)
        expect_true(dir.exists("my_dataset/seq"))
        expect_equal(Sys.readlink("my_dataset/seq"), "")
    })
})

test_that("gdataset.save() creates valid misha.yaml", {
    withr::with_tempdir({
        create_test_db("working_db")
        gsetroot("working_db")
        intervs <- gintervals(1, 0, 1000)
        gtrack.create_sparse("track1", "t1", intervs, 1)

        gdataset.save(
            path = "my_dataset",
            description = "Test dataset",
            tracks = "track1"
        )

        yaml_data <- yaml::read_yaml("my_dataset/misha.yaml")

        expect_equal(yaml_data$description, "Test dataset")
        expect_true(!is.null(yaml_data$created))
        expect_true(!is.null(yaml_data$genome))
        expect_equal(yaml_data$track_count, 1)
    })
})

test_that("gdataset.save() errors when path already exists", {
    withr::with_tempdir({
        create_test_db("working_db")
        gsetroot("working_db")
        intervs <- gintervals(1, 0, 1000)
        gtrack.create_sparse("track1", "t1", intervs, 1)

        dir.create("existing_path")

        expect_error(
            gdataset.save(
                path = "existing_path",
                description = "Test",
                tracks = "track1"
            ),
            "already exists"
        )
    })
})

test_that("gdataset.save() errors when neither tracks nor intervals specified", {
    withr::with_tempdir({
        create_test_db("working_db")
        gsetroot("working_db")

        expect_error(
            gdataset.save(
                path = "my_dataset",
                description = "Empty dataset"
            ),
            "At least one of 'tracks' or 'intervals' must be specified"
        )
    })
})

test_that("gdataset.save() errors when track doesn't exist", {
    withr::with_tempdir({
        create_test_db("working_db")
        gsetroot("working_db")

        expect_error(
            gdataset.save(
                path = "my_dataset",
                description = "Test",
                tracks = "nonexistent_track"
            ),
            "does not exist"
        )
    })
})

test_that("gdataset.save() errors when interval set doesn't exist", {
    withr::with_tempdir({
        create_test_db("working_db")
        gsetroot("working_db")

        expect_error(
            gdataset.save(
                path = "my_dataset",
                description = "Test",
                intervals = "nonexistent_intervals"
            ),
            "does not exist"
        )
    })
})

# ==============================================================================
# gdataset.ls()
# ==============================================================================

test_that("gdataset.ls() returns working db and loaded datasets", {
    withr::with_tempdir({
        create_test_db("working_db")
        gsetroot("working_db")

        create_test_db("dataset1")
        create_test_db("dataset2")

        result <- gdataset.ls()

        # Should only have working db
        expect_equal(length(result), 1)
        expect_equal(result[1], normalizePath("working_db"))

        # Load datasets
        gsetroot("dataset1")
        intervs <- gintervals(1, 0, 1000)
        gtrack.create_sparse("track1", "t1", intervs, 1)

        gsetroot("working_db")
        gdataset.load("dataset1")

        result <- gdataset.ls()
        expect_equal(length(result), 2)
        expect_equal(result[1], normalizePath("working_db"))
        expect_true(normalizePath("dataset1") %in% result)
    })
})

test_that("gdataset.ls(dataframe=TRUE) returns detailed information", {
    withr::with_tempdir({
        create_test_db("working_db")
        gsetroot("working_db")
        intervs <- gintervals(1, 0, 1000)
        gtrack.create_sparse("working_track", "wt", intervs, 1)

        create_test_db("dataset1")
        gsetroot("dataset1")
        gtrack.create_sparse("ds_track", "dt", intervs, 2)

        gsetroot("working_db")
        gdataset.load("dataset1")

        result <- gdataset.ls(dataframe = TRUE)

        expect_true(is.data.frame(result))
        expect_true("path" %in% names(result))
        expect_true("tracks_total" %in% names(result))
        expect_true("tracks_visible" %in% names(result))
        expect_true("intervals_total" %in% names(result))
        expect_true("intervals_visible" %in% names(result))
        expect_true("has_metadata" %in% names(result))
        expect_true("writable" %in% names(result))

        # Working db should be writable
        expect_true(result$writable[1])

        # Dataset should not be writable
        expect_false(result$writable[2])
    })
})

# ==============================================================================
# gdataset.info()
# ==============================================================================

test_that("gdataset.info() returns metadata for dataset with misha.yaml", {
    withr::with_tempdir({
        create_test_db("working_db")
        gsetroot("working_db")
        intervs <- gintervals(1, 0, 1000)
        gtrack.create_sparse("track1", "t1", intervs, 1)

        gdataset.save(
            path = "my_dataset",
            description = "Test dataset",
            tracks = "track1"
        )

        info <- gdataset.info("my_dataset")

        expect_equal(info$description, "Test dataset")
        expect_equal(info$track_count, 1)
        expect_equal(info$interval_count, 0)
        expect_true(!is.null(info$genome))
        expect_false(info$is_loaded)
    })
})

test_that("gdataset.info() works for dataset without misha.yaml", {
    withr::with_tempdir({
        create_test_db("dataset1")
        gsetroot("dataset1")
        intervs <- gintervals(1, 0, 1000)
        gtrack.create_sparse("track1", "t1", intervs, 1)

        create_test_db("working_db")
        gsetroot("working_db")

        info <- gdataset.info("dataset1")

        expect_null(info$description)
        expect_null(info$author)
        expect_equal(info$track_count, 1)
        expect_false(info$is_loaded)
    })
})

test_that("gdataset.info() shows is_loaded=TRUE for loaded datasets", {
    withr::with_tempdir({
        create_test_db("working_db")
        gsetroot("working_db")

        create_test_db("dataset1")
        gsetroot("dataset1")
        intervs <- gintervals(1, 0, 1000)
        gtrack.create_sparse("track1", "t1", intervs, 1)

        gsetroot("working_db")
        gdataset.load("dataset1")

        info <- gdataset.info("dataset1")
        expect_true(info$is_loaded)
    })
})

# ==============================================================================
# gtrack.dataset() and gintervals.dataset()
# ==============================================================================

test_that("gtrack.dataset() returns working db path for working db tracks", {
    withr::with_tempdir({
        create_test_db("working_db")
        gsetroot("working_db")
        intervs <- gintervals(1, 0, 1000)
        gtrack.create_sparse("working_track", "wt", intervs, 1)

        result <- gtrack.dataset("working_track")
        expect_equal(result, normalizePath("working_db"))
    })
})

test_that("gtrack.dataset() returns dataset path for dataset tracks", {
    withr::with_tempdir({
        create_test_db("working_db")
        gsetroot("working_db")

        create_test_db("dataset1")
        gsetroot("dataset1")
        intervs <- gintervals(1, 0, 1000)
        gtrack.create_sparse("ds_track", "dt", intervs, 1)

        gsetroot("working_db")
        gdataset.load("dataset1")

        result <- gtrack.dataset("ds_track")
        expect_equal(result, normalizePath("dataset1"))
    })
})

test_that("gtrack.dataset() returns NA for non-existent tracks", {
    withr::with_tempdir({
        create_test_db("working_db")
        gsetroot("working_db")

        result <- gtrack.dataset("nonexistent")
        expect_true(is.na(result))
    })
})

test_that("gtrack.dataset() works on vectors", {
    withr::with_tempdir({
        create_test_db("working_db")
        gsetroot("working_db")
        intervs <- gintervals(1, 0, 1000)
        gtrack.create_sparse("track1", "t1", intervs, 1)

        create_test_db("dataset1")
        gsetroot("dataset1")
        gtrack.create_sparse("track2", "t2", intervs, 2)

        gsetroot("working_db")
        gdataset.load("dataset1")

        result <- gtrack.dataset(c("track1", "track2", "nonexistent"))

        expect_equal(result[1], normalizePath("working_db"))
        expect_equal(result[2], normalizePath("dataset1"))
        expect_true(is.na(result[3]))
    })
})

test_that("gintervals.dataset() returns correct paths", {
    withr::with_tempdir({
        create_test_db("working_db")
        gsetroot("working_db")
        intervs1 <- gintervals(1, 0, 1000)
        gintervals.save("working_intervals", intervs1)

        create_test_db("dataset1")
        gsetroot("dataset1")
        intervs2 <- gintervals(1, 0, 2000)
        gintervals.save("ds_intervals", intervs2)

        gsetroot("working_db")
        gdataset.load("dataset1")

        expect_equal(gintervals.dataset("working_intervals"), normalizePath("working_db"))
        expect_equal(gintervals.dataset("ds_intervals"), normalizePath("dataset1"))
        expect_true(is.na(gintervals.dataset("nonexistent")))
    })
})

# ==============================================================================
# gtrack.dbs() and gintervals.dbs()
# ==============================================================================

test_that("gtrack.dbs() shows all locations including shadowed", {
    withr::with_tempdir({
        create_test_db("working_db")
        gsetroot("working_db")
        intervs <- gintervals(1, 0, 1000)
        gtrack.create_sparse("shared", "working", intervs, 1)

        create_test_db("dataset1")
        gsetroot("dataset1")
        gtrack.create_sparse("shared", "dataset", intervs, 2)

        gsetroot("working_db")
        gdataset.load("dataset1", force = TRUE)

        dbs <- gtrack.dbs("shared")

        # Should show both locations
        expect_equal(length(dbs), 2)
        expect_true(normalizePath("working_db") %in% dbs)
        expect_true(normalizePath("dataset1") %in% dbs)
    })
})

test_that("gtrack.dbs(dataframe=TRUE) returns data frame", {
    withr::with_tempdir({
        create_test_db("working_db")
        gsetroot("working_db")
        intervs <- gintervals(1, 0, 1000)
        gtrack.create_sparse("track1", "t1", intervs, 1)

        df <- gtrack.dbs("track1", dataframe = TRUE)

        expect_true(is.data.frame(df))
        expect_true("track" %in% names(df))
        expect_true("db" %in% names(df))
    })
})

test_that("gintervals.dbs() shows all locations", {
    withr::with_tempdir({
        create_test_db("working_db")
        gsetroot("working_db")
        intervs <- gintervals(1, 0, 1000)
        gintervals.save("shared", intervs)

        create_test_db("dataset1")
        gsetroot("dataset1")
        gintervals.save("shared", intervs)

        gsetroot("working_db")
        gdataset.load("dataset1", force = TRUE)

        dbs <- gintervals.dbs("shared")

        expect_equal(length(dbs), 2)
        expect_true(normalizePath("working_db") %in% dbs)
        expect_true(normalizePath("dataset1") %in% dbs)
    })
})

# ==============================================================================
# gtrack.ls() and gintervals.ls() with db parameter
# ==============================================================================

test_that("gtrack.ls() filters by database with normalized paths", {
    withr::with_tempdir({
        create_test_db("working_db")
        gsetroot("working_db")
        intervs <- gintervals(1, 0, 1000)
        gtrack.create_sparse("working_track", "wt", intervs, 1)

        create_test_db("dataset1")
        gsetroot("dataset1")
        gtrack.create_sparse("ds_track", "dt", intervs, 2)

        gsetroot("working_db")
        gdataset.load("dataset1")

        # Filter by working db
        working_tracks <- gtrack.ls(db = "working_db")
        expect_true("working_track" %in% working_tracks)
        expect_false("ds_track" %in% working_tracks)

        # Filter by dataset with different path forms
        ds_tracks1 <- gtrack.ls(db = "dataset1")
        ds_tracks2 <- gtrack.ls(db = normalizePath("dataset1"))

        expect_equal(ds_tracks1, ds_tracks2)
        expect_true("ds_track" %in% ds_tracks1)
        expect_false("working_track" %in% ds_tracks1)
    })
})

test_that("gintervals.ls() filters by database", {
    withr::with_tempdir({
        create_test_db("working_db")
        gsetroot("working_db")
        intervs1 <- gintervals(1, 0, 1000)
        gintervals.save("working_intervals", intervs1)

        create_test_db("dataset1")
        gsetroot("dataset1")
        intervs2 <- gintervals(2, 0, 2000)
        gintervals.save("ds_intervals", intervs2)

        gsetroot("working_db")
        gdataset.load("dataset1")

        working_intervals <- gintervals.ls(db = "working_db")
        expect_true("working_intervals" %in% working_intervals)
        expect_false("ds_intervals" %in% working_intervals)

        ds_intervals <- gintervals.ls(db = "dataset1")
        expect_false("working_intervals" %in% ds_intervals)
        expect_true("ds_intervals" %in% ds_intervals)
    })
})

# ==============================================================================
# Complex Scenarios
# ==============================================================================

test_that("multiple datasets can be loaded simultaneously", {
    withr::with_tempdir({
        create_test_db("working_db")
        gsetroot("working_db")

        # Create multiple datasets
        for (i in 1:3) {
            create_test_db(paste0("dataset", i))
            gsetroot(paste0("dataset", i))
            intervs <- gintervals(1, 0, 1000)
            gtrack.create_sparse(paste0("track", i), paste0("t", i), intervs, i)
        }

        gsetroot("working_db")
        gdataset.load("dataset1")
        gdataset.load("dataset2")
        gdataset.load("dataset3")

        datasets <- gdataset.ls()
        expect_equal(length(datasets), 4) # working db + 3 datasets

        # All tracks should be visible
        tracks <- gtrack.ls()
        expect_true("track1" %in% tracks)
        expect_true("track2" %in% tracks)
        expect_true("track3" %in% tracks)
    })
})

test_that("dataset workflow: save, load, use", {
    withr::with_tempdir({
        # Create working db with tracks
        create_test_db("working_db")
        gsetroot("working_db")
        intervs <- gintervals(1, 0, 1000)
        gtrack.create_sparse("analysis_track", "analysis", intervs, 42)
        gintervals.save("analysis_intervals", intervs)

        # Save as dataset
        gdataset.save(
            path = "analysis_dataset",
            description = "My analysis results",
            tracks = "analysis_track",
            intervals = "analysis_intervals"
        )

        # Create new working db
        create_test_db("new_working_db")
        gsetroot("new_working_db")

        # Load the dataset
        gdataset.load("analysis_dataset")

        # Use the loaded data
        expect_true("analysis_track" %in% gtrack.ls())
        expect_true("analysis_intervals" %in% gintervals.ls())

        result <- gextract("analysis_track", gintervals(1, 0, 500))
        expect_equal(result$analysis_track[1], 42)
    })
})

# ==============================================================================
# Regression: gintervals.load() from outer database (GH bug report)
# Prior to fix, gintervals.load() / .gintervals.is_bigset() / .gintervals.big.meta()
# hardcoded GWD instead of resolving via GINTERVALS_DATASET, so intervals
# visible via gintervals.ls() could not be loaded.
# ==============================================================================

test_that("gintervals.load() works for small intervals from loaded dataset", {
    withr::with_tempdir({
        create_test_db("working_db")
        gsetroot("working_db")

        # Create dataset with a small (single-file) intervals set
        create_test_db("outer_db")
        gsetroot("outer_db")
        intervs <- gintervals(1, 0, 1000)
        gintervals.save("outer_intervals", intervs)

        # Switch to working db and load dataset
        gsetroot("working_db")
        gdataset.load("outer_db")

        # gintervals.ls() should find it
        expect_true("outer_intervals" %in% gintervals.ls())

        # gintervals.load() must also work (this was the reported bug)
        loaded <- gintervals.load("outer_intervals")
        expect_equal(nrow(loaded), 1)
        expect_equal(as.character(loaded$chrom[1]), "chr1")
        expect_equal(loaded$start[1], 0)
        expect_equal(loaded$end[1], 1000)
    })
})

test_that("gintervals.load() works for hierarchical intervals from loaded dataset", {
    withr::with_tempdir({
        create_test_db("working_db")
        gsetroot("working_db")

        # Create dataset with a hierarchically-named intervals set (dots → subdirs)
        create_test_db("outer_db")
        gsetroot("outer_db")
        # Create parent directory structure for hierarchical name
        dir.create(file.path("outer_db", "tracks", "intervs", "global"), recursive = TRUE)
        intervs <- gintervals(1, 100, 500)
        gintervals.save("intervs.global.test_set", intervs)

        gsetroot("working_db")
        gdataset.load("outer_db")

        expect_true("intervs.global.test_set" %in% gintervals.ls())

        loaded <- gintervals.load("intervs.global.test_set")
        expect_equal(nrow(loaded), 1)
        expect_equal(loaded$start[1], 100)
        expect_equal(loaded$end[1], 500)
    })
})

test_that("gintervals.load() works for big intervals from loaded dataset", {
    withr::with_tempdir({
        create_test_db("working_db")
        gsetroot("working_db")

        # Create dataset with a big (directory-based) intervals set
        create_test_db("outer_db")
        gsetroot("outer_db")
        # Create intervals spanning multiple chroms to trigger big-set storage
        intervs <- gintervals(c(1, 2), c(0, 0), c(1000, 2000))
        gintervals.save("big_outer", intervs)

        gsetroot("working_db")
        gdataset.load("outer_db")

        expect_true("big_outer" %in% gintervals.ls())

        loaded <- gintervals.load("big_outer")
        expect_equal(nrow(loaded), 2)
    })
})

test_that("gintervals.load() with chrom filter works for dataset intervals", {
    withr::with_tempdir({
        create_test_db("working_db")
        gsetroot("working_db")

        create_test_db("outer_db")
        gsetroot("outer_db")
        intervs <- gintervals(c(1, 2), c(0, 0), c(1000, 2000))
        gintervals.save("outer_multi", intervs)

        gsetroot("working_db")
        gdataset.load("outer_db")

        # Load only chr1
        loaded <- gintervals.load("outer_multi", chrom = "chr1")
        expect_equal(nrow(loaded), 1)
        expect_equal(as.character(loaded$chrom[1]), "chr1")
    })
})

# ==============================================================================
# gdataset.save() failure paths
#
# file.copy()/file.symlink() report failure by returning FALSE and at most a
# warning. Unchecked, gdataset.save() returned success after copying nothing
# and wrote a misha.yaml claiming a track_count the directory could not back
# up - and the half-built directory then blocked every retry, because the
# function refuses a path that already exists.
# ==============================================================================

test_that("gdataset.save() errors and cleans up when a track's data is missing", {
    withr::with_tempdir({
        create_test_db("working_db")
        gsetroot("working_db")
        intervs <- gintervals(1, 0, 1000)
        gtrack.create_sparse("track1", "t1", intervs, 1)

        # The track stays registered, but its data disappears (another
        # session removed it, a shell rm -rf, a broken NFS mount).
        track_dir <- misha:::.track_dir("track1")
        unlink(track_dir, recursive = TRUE)

        expect_error(
            gdataset.save(path = "ds", description = "Test", tracks = "track1"),
            "gdataset\\.save: track track1 is registered in .* but its data is missing"
        )
        # Nothing left behind: no directory, so in particular no misha.yaml
        # advertising a track that is not there.
        expect_false(dir.exists("ds"))
        expect_false(file.exists(file.path("ds", "misha.yaml")))
    })
})

test_that("gdataset.save() failure leaves the path free for a retry", {
    withr::with_tempdir({
        create_test_db("working_db")
        gsetroot("working_db")
        intervs <- gintervals(1, 0, 1000)
        gtrack.create_sparse("track1", "t1", intervs, 1)

        track_dir <- misha:::.track_dir("track1")
        stash <- file.path(tempdir(), "stashed.track")
        expect_true(file.rename(track_dir, stash))
        expect_error(gdataset.save(path = "ds", description = "Test", tracks = "track1"))

        # Put the data back and retry at the same path.
        expect_true(file.rename(stash, track_dir))
        expect_silent(gdataset.save(path = "ds", description = "Test", tracks = "track1"))
        expect_true(dir.exists("ds/tracks/track1.track"))
        expect_equal(yaml::read_yaml("ds/misha.yaml")$track_count, 1)
    })
})

test_that("gdataset.save() cleanup removes symlinks, not their targets", {
    withr::with_tempdir({
        create_test_db("working_db")
        gsetroot("working_db")
        intervs <- gintervals(1, 0, 1000)
        gtrack.create_sparse("track1", "t1", intervs, 1)
        gtrack.create_sparse("track2", "t2", intervs, 2)

        seq_dir <- file.path(normalizePath("working_db"), "seq")
        seq_before <- sort(list.files(seq_dir))
        track2_dir <- misha:::.track_dir("track2")
        track2_before <- sort(list.files(track2_dir))

        # track1's data is gone; the seq/ and track2 symlinks have already
        # been planted into the dataset directory when the failure hits.
        unlink(misha:::.track_dir("track1"), recursive = TRUE)
        expect_error(
            gdataset.save(
                path = "ds", description = "Test",
                tracks = c("track2", "track1"), symlinks = TRUE
            )
        )
        expect_false(dir.exists("ds"))
        # The working database must be untouched by the cleanup.
        expect_equal(sort(list.files(seq_dir)), seq_before)
        expect_equal(sort(list.files(track2_dir)), track2_before)
    })
})

test_that("gdataset.save() errors when an interval set's data is missing", {
    withr::with_tempdir({
        create_test_db("working_db")
        gsetroot("working_db")
        gintervals.save("ivs1", gintervals(1, 0, 1000))

        unlink(misha:::.intervals_dir("ivs1"), recursive = TRUE)

        expect_error(
            gdataset.save(path = "ds", description = "Test", intervals = "ivs1"),
            "gdataset\\.save: interval set ivs1 is registered in .* but its data is missing"
        )
        expect_false(dir.exists("ds"))
    })
})

test_that("gdataset.load refuses a dataset that numbers chromosomes differently and holds indexed tracks", {
    local_db_state()
    td <- withr::local_tempdir()
    # indexed, the same chrom_sizes.txt as `from`, chrom ids in chrom_sizes.txt order
    indexed_like <- function(from, path, contig = NULL) {
        dir.create(file.path(path, "tracks"), recursive = TRUE)
        dir.create(file.path(path, "seq"))
        cs <- utils::read.delim(file.path(from, "chrom_sizes.txt"), header = FALSE, colClasses = c("character", "numeric"))
        names <- if (is.null(contig)) cs$V1 else paste0(contig, seq_len(nrow(cs)))
        fa <- file.path(td, paste0(basename(path), ".fa"))
        writeLines(unlist(lapply(seq_len(nrow(cs)), function(k) c(paste0(">", names[k]), strrep("A", cs$V2[k])))), fa)
        invisible(misha:::.gcall("gseq_multifasta_import", fa, file.path(path, "seq", "genome.seq"), file.path(path, "seq", "genome.idx"), FALSE, misha:::.misha_env()))
        expect_true(file.copy(file.path(from, "chrom_sizes.txt"), path))
        normalizePath(path)
    }
    with_indexed_track <- function(db) {
        gsetroot(db)
        iv <- gintervals.all()
        iv$end <- 10
        gtrack.create_sparse("it", "x", iv, seq_len(nrow(iv)))
        suppressMessages(gtrack.convert_to_indexed("it"))
        expect_true(file.exists(file.path(db, "tracks", "it.track", "track.idx")))
        db
    }
    differently <- "numbers the chromosomes differently from the working database"
    # P: per-chromosome, chrom ids by the sorted names; I: indexed, in chrom_sizes.txt order
    p <- create_db_with_unsorted_chrom_sizes(file.path(td, "P"))
    i <- with_indexed_track(indexed_like(p, file.path(td, "I")))
    pi <- with_indexed_track(create_db_with_unsorted_chrom_sizes(file.path(td, "PI")))
    gsetroot(p)
    expect_error(gdataset.load(i), differently)
    expect_equal(get("GDATASETS", envir = misha:::.misha), character(0))
    gsetroot(i)
    expect_error(gdataset.load(pi), differently)

    # files keyed by chrom id of each kind, alone: an indexed track in a subdirectory, an indexed
    # interval set, an indexed 2D interval set
    for (kind in c("nested track", "interval set", "2D interval set")) {
        k <- create_db_with_unsorted_chrom_sizes(file.path(td, gsub(" ", "_", kind)))
        gsetroot(k)
        a <- gintervals.all()
        a$end <- 10
        if (kind == "nested track") {
            gdir.create("sub", showWarnings = FALSE)
            gtrack.create_sparse("sub.nt", "x", a, seq_len(nrow(a)))
            suppressMessages(gtrack.convert_to_indexed("sub.nt"))
            expect_true(file.exists(file.path(k, "tracks", "sub", "nt.track", "track.idx")))
        } else if (kind == "interval set") {
            withr::with_options(list(gbig.intervals.size = 2), gintervals.save("bv", a))
            suppressMessages(gintervals.convert_to_indexed("bv"))
            expect_true(file.exists(file.path(k, "tracks", "bv.interv", "intervals.idx")))
        } else {
            withr::with_options(list(gbig.intervals.size = 2), gintervals.save("bv2", gintervals.2d(a$chrom, 0, 10, a$chrom[c(2:nrow(a), 1)], 0, 10)))
            suppressMessages(gintervals.2d.convert_to_indexed("bv2"))
            expect_true(file.exists(file.path(k, "tracks", "bv2.interv", "intervals2d.idx")))
        }
        gsetroot(i)
        expect_error(gdataset.load(k), differently, info = kind)
    }

    # the same order with names that differ by "chr": an indexed database whose chrom_sizes.txt is
    # sorted, and a per-chromosome one with that file (chrom ids by the sorted names, with "chr")
    sorted_root <- file.path(td, "sorted_root")
    dir.create(file.path(sorted_root, "tracks"), recursive = TRUE)
    dir.create(file.path(sorted_root, "seq"))
    gsetroot(p)
    sorted <- gintervals.all()
    writeLines(paste(sub("^chr", "", sorted$chrom), sorted$end, sep = "\t"), file.path(sorted_root, "chrom_sizes.txt"))
    sorted_root <- indexed_like(sorted_root, file.path(td, "sorted_root_i"))
    sorted_ds <- file.path(td, "sorted_ds")
    dir.create(file.path(sorted_ds, "tracks"), recursive = TRUE)
    expect_true(file.copy(file.path(p, "seq"), sorted_ds, recursive = TRUE))
    expect_true(file.copy(file.path(sorted_root, "chrom_sizes.txt"), sorted_ds))
    with_indexed_track(normalizePath(sorted_ds))
    truth <- gextract("it", gintervals.all())
    gsetroot(sorted_root)
    expect_equal(as.character(gintervals.all()$chrom), sub("^chr", "", as.character(truth$chrom)))
    suppressMessages(gdataset.load(sorted_ds))
    expect_equal(gextract("it", gintervals.all())$it, truth$it)

    # an indexed track added to a dataset after its .db.cache was written is found
    late <- create_db_with_unsorted_chrom_sizes(file.path(td, "late"))
    gsetroot(i)
    suppressMessages(gdataset.load(late))
    expect_true(file.exists(file.path(late, ".db.cache")))
    expect_true(file.copy(file.path(pi, "tracks", "it.track"), file.path(late, "tracks"), recursive = TRUE))
    gsetroot(i)
    expect_error(gdataset.load(late), differently)

    # .seq files without the "chr" prefix give chrom_sizes.txt order
    u <- create_db_with_unsorted_chrom_sizes(file.path(td, "U"))
    for (f in list.files(file.path(u, "seq"), full.names = TRUE)) file.rename(f, file.path(dirname(f), sub("^chr", "", basename(f))))
    with_indexed_track(u)
    gsetroot(p)
    expect_error(gdataset.load(u), differently)

    # the working database's order is the session's, whatever its seq/ holds: sorted with the first
    # chromosome's .seq file missing, chrom_sizes.txt order with no sequence at all
    first_missing <- file.path(td, "first_missing")
    dir.create(file.path(first_missing, "tracks"), recursive = TRUE)
    dir.create(file.path(first_missing, "seq"))
    chroms6 <- c("2", "10", "X", "3", "4", "1")
    for (k in seq_along(chroms6)) writeBin(charToRaw(strrep("A", 1000 + k)), file.path(first_missing, "seq", paste0("chr", chroms6[k], ".seq")))
    writeLines(paste(chroms6, 1000 + seq_along(chroms6), sep = "\t"), file.path(first_missing, "chrom_sizes.txt"))
    i6 <- with_indexed_track(indexed_like(first_missing, file.path(td, "I6")))
    unlink(file.path(first_missing, "seq", "chr1.seq"))
    gsetroot(first_missing)
    expect_equal(as.character(gintervals.all()$chrom)[1], "chr1")
    expect_error(gdataset.load(i6), differently)
    for (kind in c("dangling", "empty")) {
        r <- file.path(td, paste0("R_", kind))
        dir.create(file.path(r, "tracks"), recursive = TRUE)
        dir.create(file.path(r, "seq"))
        expect_true(file.copy(file.path(p, "chrom_sizes.txt"), r))
        if (kind == "dangling") {
            for (f in list.files(file.path(p, "seq"))) file.symlink(file.path(td, "nowhere", f), file.path(r, "seq", f))
        }
        gsetroot(r)
        expect_equal(as.character(gintervals.all()$chrom)[1], "2", info = kind)
        expect_error(gdataset.load(pi), differently, info = kind)
        # one without files keyed by chrom id loads
        suppressMessages(gdataset.load(p))
        expect_equal(get("GDATASETS", envir = misha:::.misha), p, info = kind)
    }

    # both indexed: chrom_sizes.txt order on both sides, even with other contig names in the index
    d <- indexed_like(p, file.path(td, "D"), contig = "contig")
    gsetroot(i)
    expect_no_warning(suppressMessages(gdataset.load(d)))
    expect_equal(get("GDATASETS", envir = misha:::.misha), d)

    # the same order loads: a dataset on the working database's seq/, a copy, one without sequence
    gsetroot(p)
    gtrack.create_sparse("t", "x", gintervals("chr1", 0, 10), 1)
    linked <- file.path(td, "linked")
    suppressMessages(gdataset.save(linked, "d", tracks = "t"))
    copied <- file.path(td, "copied")
    suppressMessages(gdataset.save(copied, "d", tracks = "t", copy_seq = TRUE))
    bare <- file.path(td, "bare")
    suppressMessages(gdataset.save(bare, "d", tracks = "t"))
    unlink(file.path(bare, "seq"))
    gtrack.rm("t", force = TRUE)
    for (ds in c(linked, copied, bare)) {
        suppressMessages(gdataset.load(ds))
        expect_equal(gextract("t", gintervals("chr1", 0, 10))$t, 1, info = ds)
        gdataset.unload(ds)
    }
})

test_that("a dataset numbering chromosomes differently without indexed tracks loads, and keeps its own order", {
    local_db_state()
    td <- withr::local_tempdir()
    p <- create_db_with_unsorted_chrom_sizes(file.path(td, "P"))
    # I: indexed, the same chrom_sizes.txt, chrom ids in chrom_sizes.txt order, no tracks
    i <- file.path(td, "I")
    dir.create(file.path(i, "tracks"), recursive = TRUE)
    dir.create(file.path(i, "seq"))
    cs <- utils::read.delim(file.path(p, "chrom_sizes.txt"), header = FALSE, colClasses = c("character", "numeric"))
    fa <- file.path(td, "i.fa")
    writeLines(unlist(lapply(seq_len(nrow(cs)), function(k) c(paste0(">", cs$V1[k]), strrep("A", cs$V2[k])))), fa)
    invisible(misha:::.gcall("gseq_multifasta_import", fa, file.path(i, "seq", "genome.seq"), file.path(i, "seq", "genome.idx"), FALSE, misha:::.misha_env()))
    expect_true(file.copy(file.path(p, "chrom_sizes.txt"), i))
    i <- normalizePath(i)
    iv <- gintervals(c("1", "2", "10", "X"), 0, 10)

    # tracks and interval sets of P, read with I as the working database by chromosome name
    gsetroot(p)
    gtrack.create_sparse("pt", "x", iv, c(1, 2, 10, 99))
    gtrack.create_sparse("pi", "x", iv, c(1, 2, 10, 99))
    suppressMessages(gtrack.convert_to_indexed("pi"))
    gtrack.2d.create("p2", "x", data.frame(chrom1 = "chr1", start1 = 0, end1 = 10, chrom2 = "chr2", start2 = 0, end2 = 10), 5)
    gintervals.save("pv", iv)
    gintervals.save("pv2", gintervals.2d("chr1", 0, 10, "chr2", 0, 10))
    by_name <- function(track) {
        res <- gextract(track, iv)
        setNames(res[[track]], sub("^chr", "", as.character(res$chrom)))[c("1", "2", "10", "X")]
    }
    truth <- by_name("pt")
    ds <- file.path(td, "ds")
    suppressMessages(gdataset.save(ds, "d", tracks = c("pt", "p2"), intervals = c("pv", "pv2"), copy_seq = TRUE))
    gsetroot(i)
    suppressMessages(gdataset.load(ds))
    expect_equal(by_name("pt"), truth)
    expect_equal(sub("^chr", "", .gdb.chrom_names_at(ds))[1], "1")
    # nothing of it is written in the indexed format in this session
    expect_error(gtrack.convert_to_indexed("pt"), "in a dataset that numbers the chromosomes differently")
    expect_error(gtrack.2d.convert_to_indexed("p2"), "in a dataset that numbers the chromosomes differently")
    expect_error(gintervals.convert_to_indexed("pv"), "in a dataset that numbers the chromosomes differently")
    expect_error(gintervals.2d.convert_to_indexed("pv2"), "in a dataset that numbers the chromosomes differently")

    # the order kept is dropped with the dataset, by gdataset.unload(), gsetroot() and a gsetroot()
    # that fails
    expect_false(is.null(get("GDATASET_CHROMS", envir = misha:::.misha)[[ds]]))
    gdataset.unload(ds)
    expect_null(get("GDATASET_CHROMS", envir = misha:::.misha)[[ds]])
    suppressMessages(gdataset.load(ds))
    # loaded again while loaded: listed once
    suppressMessages(gdataset.load(ds))
    expect_equal(get("GDATASETS", envir = misha:::.misha), ds)
    gsetroot(i)
    expect_null(get("GDATASET_CHROMS", envir = misha:::.misha))
    suppressMessages(gdataset.load(ds))
    expect_error(gsetroot(i, dir = ""), "empty string")
    expect_null(get("GDATASET_CHROMS", envir = misha:::.misha))

    # copies into I loaded as a dataset of P, of an indexed and a per-chromosome track, by db = and
    # from inside its tracks/, are written in files named by chromosome: they read right in this
    # session and with I as the working database
    gsetroot(p)
    suppressMessages(gdataset.load(i))
    gtrack.copy("pi", "pic", db = i)
    gtrack.copy("pt", "ptc", db = i)
    gdir.cd(file.path(i, "tracks"))
    gtrack.copy("pi", "pic2")
    gdir.cd(file.path(p, "tracks"))
    for (t in c("pic", "ptc", "pic2")) {
        expect_false(file.exists(file.path(i, "tracks", paste0(t, ".track"), "track.idx")), info = t)
        expect_equal(by_name(t), truth, info = t)
    }
    gsetroot(i)
    expect_equal(by_name("pic"), truth)
    expect_equal(by_name("ptc"), truth)
    expect_equal(by_name("pic2"), truth)
    gsetroot(p)
    suppressMessages(gdataset.load(i))
    expect_equal(get("GDATASETS", envir = misha:::.misha), i)
})

test_that("gdataset.load refuses a dataset whose linked seq/ belongs to a database converted since", {
    local_db_state()
    td <- withr::local_tempdir()
    p <- create_db_with_unsorted_chrom_sizes(file.path(td, "P"))
    # I: indexed, the same chrom_sizes.txt, chrom ids in chrom_sizes.txt order
    i <- file.path(td, "I")
    dir.create(file.path(i, "tracks"), recursive = TRUE)
    dir.create(file.path(i, "seq"))
    cs <- utils::read.delim(file.path(p, "chrom_sizes.txt"), header = FALSE, colClasses = c("character", "numeric"))
    fa <- file.path(td, "i.fa")
    writeLines(unlist(lapply(seq_len(nrow(cs)), function(k) c(paste0(">", cs$V1[k]), strrep("A", cs$V2[k])))), fa)
    invisible(misha:::.gcall("gseq_multifasta_import", fa, file.path(i, "seq", "genome.seq"), file.path(i, "seq", "genome.idx"), FALSE, misha:::.misha_env()))
    expect_true(file.copy(file.path(p, "chrom_sizes.txt"), i))
    i <- normalizePath(i)
    # D: an indexed track of P (sorted chrom ids), saved with copy_seq = FALSE
    gsetroot(p)
    a <- gintervals.all()
    a$end <- 10
    gtrack.create_sparse("t", "x", a, seq_len(nrow(a)))
    suppressMessages(gtrack.convert_to_indexed("t"))
    d <- file.path(td, "D")
    suppressMessages(gdataset.save(d, "d", tracks = "t", copy_seq = FALSE))
    gtrack.rm("t", force = TRUE)
    gsetroot(i)
    expect_error(gdataset.load(d), "numbers the chromosomes differently")
    # P converted since: D's seq/ is indexed in P's order, which D's chrom_sizes.txt does not match
    suppressMessages(gdb.convert_to_indexed(groot = p, force = TRUE, validate = FALSE))
    expect_true(file.exists(file.path(d, "seq", "genome.idx")))
    gsetroot(i)
    expect_error(gdataset.load(d), "its seq/genome.idx does not match its chrom_sizes.txt")
})

test_that("gdataset.load gives no warning for a dataset's index names", {
    local_db_state()
    td <- withr::local_tempdir()
    p <- create_db_with_unsorted_chrom_sizes(file.path(td, "P"))
    # indexed, contigs named otherwise in genome.idx than in chrom_sizes.txt (as tgdb NZW_T2T)
    d <- file.path(td, "D")
    dir.create(file.path(d, "tracks"), recursive = TRUE)
    dir.create(file.path(d, "seq"))
    cs <- utils::read.delim(file.path(p, "chrom_sizes.txt"), header = FALSE, colClasses = c("character", "numeric"))
    fa <- file.path(td, "d.fa")
    writeLines(unlist(lapply(seq_len(nrow(cs)), function(k) c(paste0(">contig", k), strrep("A", cs$V2[k])))), fa)
    invisible(misha:::.gcall("gseq_multifasta_import", fa, file.path(d, "seq", "genome.seq"), file.path(d, "seq", "genome.idx"), FALSE, misha:::.misha_env()))
    expect_true(file.copy(file.path(p, "chrom_sizes.txt"), d))
    gsetroot(p)
    expect_no_warning(suppressMessages(gdataset.load(d)))
})
