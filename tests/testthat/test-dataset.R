# Tests for dataset API
# Following spec: dev/notes/2026-01-13-dataset-api-spec.md

# Helper to create a minimal database
create_test_db <- function(path, chrom_sizes = data.frame(chrom = c("chr1", "chr2"), size = c(10000, 10000))) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(path, "tracks"), showWarnings = FALSE)
    dir.create(file.path(path, "seq"), showWarnings = FALSE)

    write.table(chrom_sizes, file.path(path, "chrom_sizes.txt"),
        sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE
    )

    # Create dummy sequence files
    for (chr in chrom_sizes$chrom) {
        writeLines(
            paste0(rep("A", chrom_sizes$size[chrom_sizes$chrom == chr]), collapse = ""),
            file.path(path, "seq", paste0(chr, ".seq"))
        )
    }
}

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
