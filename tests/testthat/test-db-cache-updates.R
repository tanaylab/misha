copy_seed_db <- function(seed_root, dst) {
    dir.create(dst, recursive = TRUE, showWarnings = FALSE)
    entries <- list.files(seed_root,
        all.files = TRUE,
        full.names = TRUE,
        include.dirs = TRUE,
        no.. = TRUE
    )
    file.copy(entries, dst, recursive = TRUE, copy.mode = TRUE)
}

test_that("track cache stays updated across database switches without rescan", {
    gdb.init_examples()
    withr::defer(gdb.init_examples())

    seed_root <- get("GROOT", envir = .misha)
    tmp_base <- withr::local_tempdir()

    db1 <- file.path(tmp_base, "db1")
    db2 <- file.path(tmp_base, "db2")
    copy_seed_db(seed_root, db1)
    copy_seed_db(seed_root, db2)

    gsetroot(db1, rescan = TRUE)
    gtrack.create("tmp_track_db1", "tmp", "dense_track")
    withr::defer({
        if (gtrack.exists("tmp_track_db1")) {
            gtrack.rm("tmp_track_db1", force = TRUE)
        }
    }, envir = parent.frame())
    expect_true(gtrack.exists("tmp_track_db1"))

    gsetroot(db2, rescan = TRUE)
    expect_false("tmp_track_db1" %in% gtrack.ls())

    gsetroot(db1, rescan = FALSE)
    expect_true(gtrack.exists("tmp_track_db1"))

    gtrack.rm("tmp_track_db1", force = TRUE)
    gsetroot(db2, rescan = TRUE)
    gsetroot(db1, rescan = FALSE)
    expect_false(gtrack.exists("tmp_track_db1"))
})

test_that("interval cache stays updated across database switches without rescan", {
    gdb.init_examples()
    withr::defer(gdb.init_examples())

    seed_root <- get("GROOT", envir = .misha)
    tmp_base <- withr::local_tempdir()

    db1 <- file.path(tmp_base, "db1")
    db2 <- file.path(tmp_base, "db2")
    copy_seed_db(seed_root, db1)
    copy_seed_db(seed_root, db2)

    gsetroot(db1, rescan = TRUE)
    gintervals.save("tmp_interv_db1", gintervals(1, 0, 100))
    withr::defer({
        if ("tmp_interv_db1" %in% gintervals.ls()) {
            gintervals.rm("tmp_interv_db1", force = TRUE)
        }
    }, envir = parent.frame())
    expect_true("tmp_interv_db1" %in% gintervals.ls())

    gsetroot(db2, rescan = TRUE)
    expect_false("tmp_interv_db1" %in% gintervals.ls())

    gsetroot(db1, rescan = FALSE)
    expect_true("tmp_interv_db1" %in% gintervals.ls())

    gintervals.rm("tmp_interv_db1", force = TRUE)
    gsetroot(db2, rescan = TRUE)
    gsetroot(db1, rescan = FALSE)
    expect_false("tmp_interv_db1" %in% gintervals.ls())
})

test_that("dirty cache flag forces rescan on gsetroot", {
    gdb.init_examples()
    withr::defer(gdb.init_examples())

    seed_root <- get("GROOT", envir = .misha)
    tmp_base <- withr::local_tempdir()

    db1 <- file.path(tmp_base, "db1")
    db2 <- file.path(tmp_base, "db2")
    copy_seed_db(seed_root, db1)
    copy_seed_db(seed_root, db2)

    gsetroot(db1, rescan = TRUE)
    gtrack.create("tmp_track_dirty", "tmp", "dense_track")
    expect_true(gtrack.exists("tmp_track_dirty"))

    # Simulate an out-of-band change: remove track files and mark cache dirty
    unlink(file.path(db1, "tracks", "tmp_track_dirty.track"), recursive = TRUE)
    gdb.mark_cache_dirty()

    gsetroot(db2, rescan = TRUE)
    gsetroot(db1, rescan = FALSE)

    expect_false(gtrack.exists("tmp_track_dirty"))
    expect_false(file.exists(file.path(db1, ".db.cache.dirty")))
})

test_that("mark cache dirty helps detect manual interval deletions", {
    gdb.init_examples()
    withr::defer(gdb.init_examples())

    seed_root <- get("GROOT", envir = .misha)
    tmp_base <- withr::local_tempdir()

    db1 <- file.path(tmp_base, "db1")
    db2 <- file.path(tmp_base, "db2")
    copy_seed_db(seed_root, db1)
    copy_seed_db(seed_root, db2)

    gsetroot(db1, rescan = TRUE)
    gintervals.save("tmp_interv_dirty", gintervals(1, 0, 50))
    expect_true("tmp_interv_dirty" %in% gintervals.ls())

    # Remove interval set manually and mark cache dirty
    unlink(file.path(db1, "tracks", "tmp_interv_dirty.interv"), recursive = TRUE)
    gdb.mark_cache_dirty()

    gsetroot(db2, rescan = TRUE)
    gsetroot(db1, rescan = FALSE)

    expect_false("tmp_interv_dirty" %in% gintervals.ls())
    expect_false(file.exists(file.path(db1, ".db.cache.dirty")))
})

test_that("dirty flag contains timestamp for debugging", {
    gdb.init_examples()
    withr::defer(gdb.init_examples())

    seed_root <- get("GROOT", envir = .misha)
    tmp_base <- withr::local_tempdir()
    db1 <- file.path(tmp_base, "db1")
    copy_seed_db(seed_root, db1)

    gsetroot(db1, rescan = TRUE)

    # Mark cache dirty and verify timestamp exists
    gdb.mark_cache_dirty()
    dirty_path <- file.path(db1, ".db.cache.dirty")
    expect_true(file.exists(dirty_path))

    # Read the timestamp
    timestamp_line <- readLines(dirty_path, n = 1)
    expect_true(nchar(timestamp_line) > 0)

    # Verify it's a valid timestamp (can be parsed)
    timestamp <- tryCatch(
        as.POSIXct(timestamp_line),
        error = function(e) NULL
    )
    expect_false(is.null(timestamp))
})

test_that("cache write failures produce warnings", {
    # Test with invalid serialization (function objects can't be serialized)
    gdb.init_examples()
    withr::defer(gdb.init_examples())

    seed_root <- get("GROOT", envir = .misha)
    tmp_base <- withr::local_tempdir()
    db1 <- file.path(tmp_base, "db1")
    copy_seed_db(seed_root, db1)

    gsetroot(db1, rescan = TRUE)

    # Create a corrupted cache by making parent directory a file instead of directory
    cache_path <- file.path(db1, ".db.cache")
    parent_as_file <- file.path(dirname(cache_path), "blocking_file")

    # Create a file where a directory should be to force failure
    writeLines("block", parent_as_file)

    # Attempt to write cache to a path where parent is a file should produce a warning
    # This simulates a filesystem error
    expect_warning(
        result <- misha:::.gdb.cache_write_lists(c("track1"), c(), file.path(parent_as_file, "subdir")),
        "Failed to write database cache|Failed to mark cache dirty"
    )
    expect_false(result)
})

test_that("dirty flag operations handle failures gracefully", {
    gdb.init_examples()
    withr::defer(gdb.init_examples())

    seed_root <- get("GROOT", envir = .misha)
    tmp_base <- withr::local_tempdir()
    db1 <- file.path(tmp_base, "db1")
    copy_seed_db(seed_root, db1)

    gsetroot(db1, rescan = TRUE)

    # Mark cache dirty successfully
    result <- misha:::.gdb.cache_mark_dirty(db1)
    expect_true(result)
    dirty_path <- file.path(db1, ".db.cache.dirty")
    expect_true(file.exists(dirty_path))

    # Test that operations with invalid paths produce warnings
    invalid_path <- file.path(tmp_base, "nonexistent_parent", "db")

    # Create a file where directory should be
    blocking_file <- file.path(tmp_base, "blocking_file")
    writeLines("block", blocking_file)
    invalid_with_file_parent <- file.path(blocking_file, "db")

    # Should issue a warning and return FALSE when parent is a file
    expect_warning(
        result <- misha:::.gdb.cache_mark_dirty(invalid_with_file_parent),
        "Failed to mark cache dirty"
    )
    expect_false(result)
})

test_that("cache operations return correct success/failure status", {
    gdb.init_examples()
    withr::defer(gdb.init_examples())

    seed_root <- get("GROOT", envir = .misha)
    tmp_base <- withr::local_tempdir()
    db1 <- file.path(tmp_base, "db1")
    copy_seed_db(seed_root, db1)

    gsetroot(db1, rescan = TRUE)

    # Successful operations should return TRUE
    expect_true(misha:::.gdb.cache_mark_dirty(db1))
    expect_true(misha:::.gdb.cache_write_lists(c("track1"), c(), db1))
    expect_true(misha:::.gdb.cache_clear_dirty(db1))

    # Operations on NULL groot use current GROOT from environment (valid behavior)
    # So they should succeed when a database is set
    expect_true(misha:::.gdb.cache_mark_dirty(NULL))
    expect_true(misha:::.gdb.cache_write_lists(c("track1"), c(), NULL))
    expect_true(misha:::.gdb.cache_clear_dirty(NULL))

    # Operations on empty string groot should return FALSE
    suppressWarnings(expect_false(misha:::.gdb.cache_mark_dirty("")))
    expect_false(misha:::.gdb.cache_write_lists(c("track1"), c(), ""))
    expect_false(misha:::.gdb.cache_clear_dirty(""))

    # Test that clear_dirty returns TRUE when file doesn't exist
    # (already cleaned up above)
    expect_true(misha:::.gdb.cache_clear_dirty(db1))
})
