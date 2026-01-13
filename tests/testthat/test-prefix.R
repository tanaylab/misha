# Test cases for database prefix functionality

# Helper to create a test database with prefix
create_prefixed_test_db <- function(prefix = NULL, name = "testdb") {
    db_dir <- tempfile(pattern = paste0("misha_", name, "_"), tmpdir = tempdir())
    dir.create(db_dir, showWarnings = FALSE)
    dir.create(file.path(db_dir, "tracks"), showWarnings = FALSE)
    dir.create(file.path(db_dir, "seq"), showWarnings = FALSE)

    # Create minimal chrom_sizes.txt
    writeLines(c("chr1\t1000", "chr2\t2000"), file.path(db_dir, "chrom_sizes.txt"))

    # Create minimal seq files
    for (chrom in c("chr1", "chr2")) {
        chrom_dir <- file.path(db_dir, "seq", chrom)
        dir.create(chrom_dir, showWarnings = FALSE)
        # Create empty placeholder file
        writeLines("", file.path(chrom_dir, "seq.txt"))
    }

    # Create .misha config if prefix is specified
    if (!is.null(prefix)) {
        yaml::write_yaml(list(prefix = prefix), file.path(db_dir, ".misha"))
    }

    db_dir
}

# Helper to create a track directory structure in a database
create_test_track <- function(db_dir, trackname) {
    track_path <- file.path(db_dir, "tracks", paste0(gsub("\\.", "/", trackname), ".track"))
    dir.create(track_path, recursive = TRUE, showWarnings = FALSE)

    # Create minimal track files (sparse track structure)
    writeLines("sparse", file.path(track_path, "type"))
    invisible(track_path)
}

# Helper to create an intervals file in a database
create_test_intervals <- function(db_dir, intervname) {
    interv_path <- file.path(db_dir, "tracks", paste0(gsub("\\.", "/", intervname), ".interv"))
    # Create parent directory
    dir.create(dirname(interv_path), recursive = TRUE, showWarnings = FALSE)
    # Create minimal intervals file
    df <- data.frame(chrom = "chr1", start = 0, end = 100)
    saveRDS(df, interv_path)
    invisible(interv_path)
}


test_that("prefix parsing works correctly", {
    skip_if_not_installed("yaml")

    # Test parsing prefixed names
    parsed <- misha:::.gparse_prefixed_name("at@my.track")
    expect_equal(parsed$prefix, "at")
    expect_equal(parsed$name, "my.track")

    # Test parsing unprefixed names
    parsed <- misha:::.gparse_prefixed_name("my.track")
    expect_null(parsed$prefix)
    expect_equal(parsed$name, "my.track")

    # Test parsing names with multiple @
    parsed <- misha:::.gparse_prefixed_name("at@my@track")
    expect_equal(parsed$prefix, "at")
    expect_equal(parsed$name, "my@track")

    # Test NULL input
    parsed <- misha:::.gparse_prefixed_name(NULL)
    expect_null(parsed$prefix)
    expect_null(parsed$name)

    # Test empty prefix error
    expect_error(misha:::.gparse_prefixed_name("@track"), "non-empty")

    # Test empty name error
    expect_error(misha:::.gparse_prefixed_name("prefix@"), "non-empty")
})


test_that("ghas_prefix detects prefixes correctly", {
    skip_if_not_installed("yaml")

    expect_true(misha:::.ghas_prefix("at@my.track"))
    expect_false(misha:::.ghas_prefix("my.track"))
    expect_false(misha:::.ghas_prefix(NULL))
    expect_false(misha:::.ghas_prefix(""))
})


test_that("prefix validation works correctly", {
    skip_if_not_installed("yaml")

    # Valid prefixes
    expect_silent(misha:::.gdb.validate_prefix("at"))
    expect_silent(misha:::.gdb.validate_prefix("mydb"))
    expect_silent(misha:::.gdb.validate_prefix("db_v1"))
    expect_silent(misha:::.gdb.validate_prefix("A1"))

    # Invalid prefixes
    expect_error(misha:::.gdb.validate_prefix("1db"), "must start with a letter")
    expect_error(misha:::.gdb.validate_prefix("_db"), "must start with a letter")
    # @ is caught by the alphanumeric pattern check
    expect_error(misha:::.gdb.validate_prefix("at@db"), "alphanumeric")
    expect_error(misha:::.gdb.validate_prefix(""), "cannot be empty")
    expect_error(misha:::.gdb.validate_prefix(NULL), "must be a single string")

    # Warning for long prefix
    expect_warning(misha:::.gdb.validate_prefix("verylongprefix"), "longer than recommended")
})


test_that("single prefixed database works", {
    skip_if_not_installed("yaml")

    # Create a global database and a prefixed database
    global_db <- create_prefixed_test_db(prefix = NULL, name = "global")
    prefixed_db <- create_prefixed_test_db(prefix = "at", name = "prefixed")

    withr::defer({
        unlink(global_db, recursive = TRUE)
        unlink(prefixed_db, recursive = TRUE)
    })

    # Create tracks in each
    create_test_track(global_db, "global.track")
    create_test_track(prefixed_db, "prefixed.track")

    # Connect both databases
    gsetroot(c(global_db, prefixed_db))

    # Check gdb.prefixes returns correct mapping
    prefixes <- gdb.prefixes()
    expect_equal(length(prefixes), 1)
    expect_equal(names(prefixes), "at")
    expect_equal(unname(prefixes), prefixed_db)

    # Check gdb.config returns correct info
    config <- gdb.config("at")
    expect_equal(config$prefix, "at")

    # Check track listing
    tracks <- gtrack.ls()
    expect_true("global.track" %in% tracks)
    expect_true("at@prefixed.track" %in% tracks)

    # Check track resolution
    expect_equal(gtrack.db("global.track"), global_db)
    expect_equal(gtrack.db("at@prefixed.track"), prefixed_db)
})


test_that("duplicate prefixes are rejected", {
    skip_if_not_installed("yaml")

    # Create two databases with the same prefix
    db1 <- create_prefixed_test_db(prefix = "at", name = "db1")
    db2 <- create_prefixed_test_db(prefix = "at", name = "db2")

    withr::defer({
        unlink(db1, recursive = TRUE)
        unlink(db2, recursive = TRUE)
    })

    # Connecting should fail due to duplicate prefix
    expect_error(gsetroot(c(db1, db2)), "Duplicate prefix")
})


test_that("track creation with prefix works", {
    skip_if_not_installed("yaml")

    # Create a prefixed database
    prefixed_db <- create_prefixed_test_db(prefix = "at", name = "create_test")

    withr::defer({
        unlink(prefixed_db, recursive = TRUE)
    })

    # Create temp subdirectory for creating tracks
    dir.create(file.path(prefixed_db, "tracks", "temp"), showWarnings = FALSE)

    # Connect
    gsetroot(prefixed_db)

    # Create a sparse track with prefix
    intervs <- data.frame(chrom = "chr1", start = c(0, 100, 200), end = c(50, 150, 250))
    values <- c(1, 2, 3)
    gtrack.create_sparse("at@temp.new_track", "Test track", intervs, values)

    # Check it was created with the prefixed name
    tracks <- gtrack.ls()
    expect_true("at@temp.new_track" %in% tracks)

    # Check track exists and can be accessed
    expect_true(gtrack.exists("at@temp.new_track"))
    info <- gtrack.info("at@temp.new_track")
    expect_equal(info$type, "sparse")

    # Clean up
    gtrack.rm("at@temp.new_track", force = TRUE)
})


test_that("gtrack.ls filtering by prefix works", {
    skip_if_not_installed("yaml")

    global_db <- create_prefixed_test_db(prefix = NULL, name = "global2")
    prefixed_db <- create_prefixed_test_db(prefix = "at", name = "prefixed2")

    withr::defer({
        unlink(global_db, recursive = TRUE)
        unlink(prefixed_db, recursive = TRUE)
    })

    # Create tracks in each
    create_test_track(global_db, "global.track1")
    create_test_track(global_db, "global.track2")
    create_test_track(prefixed_db, "prefixed.track1")
    create_test_track(prefixed_db, "prefixed.track2")

    gsetroot(c(global_db, prefixed_db))

    # All tracks
    all_tracks <- gtrack.ls()
    expect_equal(length(all_tracks), 4)

    # Filter by prefix
    at_tracks <- gtrack.ls(db = "at")
    expect_equal(length(at_tracks), 2)
    expect_true(all(grepl("^at@", at_tracks)))

    # Filter by global database
    global_tracks <- gtrack.ls(db = global_db)
    expect_equal(length(global_tracks), 2)
    expect_false(any(grepl("@", global_tracks)))
})


test_that("cross-database track copy works", {
    skip_if_not_installed("yaml")

    # Use two prefixed databases for clear cross-db copy
    db1 <- create_prefixed_test_db(prefix = "src", name = "source_db")
    db2 <- create_prefixed_test_db(prefix = "dst", name = "dest_db")

    withr::defer({
        unlink(db1, recursive = TRUE)
        unlink(db2, recursive = TRUE)
    })

    # Create temp directories
    dir.create(file.path(db1, "tracks", "temp"), showWarnings = FALSE)
    dir.create(file.path(db2, "tracks", "temp"), showWarnings = FALSE)

    gsetroot(c(db1, db2))

    # Create a track in the source database using explicit prefix
    intervs <- data.frame(chrom = "chr1", start = c(0, 100), end = c(50, 150))
    values <- c(1, 2)
    gtrack.create_sparse("src@temp.source", "Source track", intervs, values)

    # Copy to destination database
    gtrack.copy("src@temp.source", "dst@temp.dest")

    # Verify the copy
    expect_true(gtrack.exists("dst@temp.dest"))
    expect_equal(gtrack.db("dst@temp.dest"), db2)
    expect_equal(gtrack.db("src@temp.source"), db1)

    # Clean up
    gtrack.rm("src@temp.source", force = TRUE)
    gtrack.rm("dst@temp.dest", force = TRUE)
})


test_that("gintervals.ls filtering by prefix works", {
    skip_if_not_installed("yaml")

    global_db <- create_prefixed_test_db(prefix = NULL, name = "global4")
    prefixed_db <- create_prefixed_test_db(prefix = "at", name = "prefixed4")

    withr::defer({
        unlink(global_db, recursive = TRUE)
        unlink(prefixed_db, recursive = TRUE)
    })

    # Create intervals in each
    create_test_intervals(global_db, "global.interv")
    create_test_intervals(prefixed_db, "prefixed.interv")

    gsetroot(c(global_db, prefixed_db))

    # All intervals
    all_intervs <- gintervals.ls()
    expect_true("global.interv" %in% all_intervs)
    expect_true("at@prefixed.interv" %in% all_intervs)

    # Filter by prefix
    at_intervs <- gintervals.ls(db = "at")
    expect_equal(length(at_intervs), 1)
    expect_equal(at_intervs, "at@prefixed.interv")
})


test_that("backward compatibility with no .misha file works", {
    skip_if_not_installed("yaml")

    # Create database without .misha file
    db_dir <- create_prefixed_test_db(prefix = NULL, name = "noconfig")

    withr::defer({
        unlink(db_dir, recursive = TRUE)
    })

    create_test_track(db_dir, "my.track")

    # Should work without .misha file
    gsetroot(db_dir)

    # Track names should not be prefixed
    tracks <- gtrack.ls()
    expect_true("my.track" %in% tracks)
    expect_false(any(grepl("@", tracks)))

    # gdb.prefixes should return empty
    prefixes <- gdb.prefixes()
    expect_equal(length(prefixes), 0)
})


test_that("prefix suggestions work in error messages", {
    skip_if_not_installed("yaml")

    global_db <- create_prefixed_test_db(prefix = NULL, name = "global5")
    prefixed_db <- create_prefixed_test_db(prefix = "at", name = "prefixed5")

    withr::defer({
        unlink(global_db, recursive = TRUE)
        unlink(prefixed_db, recursive = TRUE)
    })

    # Create track only in prefixed database
    create_test_track(prefixed_db, "unique.track")

    gsetroot(c(global_db, prefixed_db))

    # Looking for unprefixed track should suggest prefixed version
    error_msg <- tryCatch(
        gtrack.dbs("unique.track"),
        error = function(e) e$message
    )
    expect_match(error_msg, "Did you mean")
    expect_match(error_msg, "at@unique.track")
})


test_that("cross-database track move works", {
    skip_if_not_installed("yaml")

    db1 <- create_prefixed_test_db(prefix = "src", name = "move_source")
    db2 <- create_prefixed_test_db(prefix = "dst", name = "move_dest")

    withr::defer({
        unlink(db1, recursive = TRUE)
        unlink(db2, recursive = TRUE)
    })

    dir.create(file.path(db1, "tracks", "temp"), showWarnings = FALSE)
    dir.create(file.path(db2, "tracks", "temp"), showWarnings = FALSE)

    gsetroot(c(db1, db2))

    # Create a track in source
    intervs <- data.frame(chrom = "chr1", start = c(0, 100), end = c(50, 150))
    values <- c(1, 2)
    gtrack.create_sparse("src@temp.moveme", "Track to move", intervs, values)

    expect_true(gtrack.exists("src@temp.moveme"))

    # Move to destination
    gtrack.mv("src@temp.moveme", "dst@temp.moved")

    # Verify move
    expect_false(gtrack.exists("src@temp.moveme"))
    expect_true(gtrack.exists("dst@temp.moved"))
    expect_equal(gtrack.db("dst@temp.moved"), db2)

    # Clean up
    gtrack.rm("dst@temp.moved", force = TRUE)
})


test_that("multiple prefixed databases work together", {
    skip_if_not_installed("yaml")

    db1 <- create_prefixed_test_db(prefix = "alpha", name = "multi1")
    db2 <- create_prefixed_test_db(prefix = "beta", name = "multi2")
    db3 <- create_prefixed_test_db(prefix = "gamma", name = "multi3")

    withr::defer({
        unlink(db1, recursive = TRUE)
        unlink(db2, recursive = TRUE)
        unlink(db3, recursive = TRUE)
    })

    # Create tracks in each
    create_test_track(db1, "alpha.track")
    create_test_track(db2, "beta.track")
    create_test_track(db3, "gamma.track")

    gsetroot(c(db1, db2, db3))

    # Check all prefixes registered
    prefixes <- gdb.prefixes()
    expect_equal(length(prefixes), 3)
    expect_true("alpha" %in% names(prefixes))
    expect_true("beta" %in% names(prefixes))
    expect_true("gamma" %in% names(prefixes))

    # Check tracks
    tracks <- gtrack.ls()
    expect_true("alpha@alpha.track" %in% tracks)
    expect_true("beta@beta.track" %in% tracks)
    expect_true("gamma@gamma.track" %in% tracks)

    # Filter by each prefix
    expect_equal(gtrack.ls(db = "alpha"), "alpha@alpha.track")
    expect_equal(gtrack.ls(db = "beta"), "beta@beta.track")
    expect_equal(gtrack.ls(db = "gamma"), "gamma@gamma.track")
})


test_that("gdb.config returns database information", {
    skip_if_not_installed("yaml")

    db <- create_prefixed_test_db(prefix = "test", name = "config_test")

    withr::defer({
        unlink(db, recursive = TRUE)
    })

    # Add more config options
    yaml::write_yaml(
        list(
            prefix = "test",
            description = "Test database",
            author = "Test Author",
            version = "1.0"
        ),
        file.path(db, ".misha")
    )

    gsetroot(db)

    config <- gdb.config("test")
    expect_equal(config$prefix, "test")
    expect_equal(config$description, "Test database")
    expect_equal(config$author, "Test Author")
    expect_equal(config$version, "1.0")
})


test_that("gtrack.rm with prefix removes correct track", {
    skip_if_not_installed("yaml")

    db <- create_prefixed_test_db(prefix = "rm", name = "rm_test")

    withr::defer({
        unlink(db, recursive = TRUE)
    })

    dir.create(file.path(db, "tracks", "test"), showWarnings = FALSE)

    gsetroot(db)

    # Create multiple tracks
    intervs <- data.frame(chrom = "chr1", start = c(0), end = c(50))
    gtrack.create_sparse("rm@test.track1", "Track 1", intervs, 1)
    gtrack.create_sparse("rm@test.track2", "Track 2", intervs, 2)

    expect_true(gtrack.exists("rm@test.track1"))
    expect_true(gtrack.exists("rm@test.track2"))

    # Remove one track
    gtrack.rm("rm@test.track1", force = TRUE)

    expect_false(gtrack.exists("rm@test.track1"))
    expect_true(gtrack.exists("rm@test.track2"))

    # Clean up
    gtrack.rm("rm@test.track2", force = TRUE)
})


test_that("unknown prefix generates error for track info", {
    skip_if_not_installed("yaml")

    db <- create_prefixed_test_db(prefix = "known", name = "unknown_test")

    withr::defer({
        unlink(db, recursive = TRUE)
    })

    gsetroot(db)

    # gtrack.exists returns FALSE for unknown prefixes (doesn't error)
    expect_false(gtrack.exists("unknown@some.track"))

    # gtrack.info errors on unknown prefix
    expect_error(gtrack.info("unknown@some.track"), "Unknown database prefix")
})


test_that("prefix with metadata only (no prefix field) works as global", {
    skip_if_not_installed("yaml")

    db <- create_prefixed_test_db(prefix = NULL, name = "metadata_only")

    withr::defer({
        unlink(db, recursive = TRUE)
    })

    # Create .misha with metadata but no prefix
    yaml::write_yaml(
        list(
            description = "Database with metadata only",
            author = "Test Author"
        ),
        file.path(db, ".misha")
    )

    create_test_track(db, "my.track")

    gsetroot(db)

    # Should work as a global database
    prefixes <- gdb.prefixes()
    expect_equal(length(prefixes), 0)

    # Tracks should not be prefixed
    tracks <- gtrack.ls()
    expect_true("my.track" %in% tracks)
    expect_false(any(grepl("@", tracks)))
})


test_that("gintervals.save with prefix creates in correct database", {
    skip_if_not_installed("yaml")

    db1 <- create_prefixed_test_db(prefix = "int1", name = "interv_save1")
    db2 <- create_prefixed_test_db(prefix = "int2", name = "interv_save2")

    withr::defer({
        unlink(db1, recursive = TRUE)
        unlink(db2, recursive = TRUE)
    })

    dir.create(file.path(db1, "tracks", "test"), showWarnings = FALSE)
    dir.create(file.path(db2, "tracks", "test"), showWarnings = FALSE)

    gsetroot(c(db1, db2))

    # Create intervals in second database using prefix
    intervs <- data.frame(chrom = "chr1", start = c(0, 100), end = c(50, 150))
    gintervals.save("int2@test.intervals", intervs)

    # Verify
    all_intervs <- gintervals.ls()
    expect_true("int2@test.intervals" %in% all_intervs)

    # Filter by prefix
    int2_intervs <- gintervals.ls(db = "int2")
    expect_equal(length(int2_intervs), 1)
    expect_equal(int2_intervs, "int2@test.intervals")

    # Clean up
    gintervals.rm("int2@test.intervals", force = TRUE)
})


test_that("gtrack.exists handles prefixed and unprefixed names", {
    skip_if_not_installed("yaml")

    global_db <- create_prefixed_test_db(prefix = NULL, name = "exists_global")
    prefixed_db <- create_prefixed_test_db(prefix = "ex", name = "exists_prefixed")

    withr::defer({
        unlink(global_db, recursive = TRUE)
        unlink(prefixed_db, recursive = TRUE)
    })

    create_test_track(global_db, "global.track")
    create_test_track(prefixed_db, "prefixed.track")

    gsetroot(c(global_db, prefixed_db))

    # Global track
    expect_true(gtrack.exists("global.track"))
    expect_false(gtrack.exists("ex@global.track")) # Not in prefixed db

    # Prefixed track
    expect_true(gtrack.exists("ex@prefixed.track"))
    expect_false(gtrack.exists("prefixed.track")) # Not in global db
})


test_that("same track name in different databases works with prefixes", {
    skip_if_not_installed("yaml")

    db1 <- create_prefixed_test_db(prefix = "db1", name = "same_name1")
    db2 <- create_prefixed_test_db(prefix = "db2", name = "same_name2")

    withr::defer({
        unlink(db1, recursive = TRUE)
        unlink(db2, recursive = TRUE)
    })

    # Create track with same base name in both databases
    create_test_track(db1, "common.track")
    create_test_track(db2, "common.track")

    gsetroot(c(db1, db2))

    # Both should be accessible
    tracks <- gtrack.ls()
    expect_true("db1@common.track" %in% tracks)
    expect_true("db2@common.track" %in% tracks)

    # Each resolves to its own database
    expect_equal(gtrack.db("db1@common.track"), db1)
    expect_equal(gtrack.db("db2@common.track"), db2)
})


test_that("prefix name qualification and unqualification work", {
    skip_if_not_installed("yaml")

    db <- create_prefixed_test_db(prefix = "qual", name = "qualify_test")

    withr::defer({
        unlink(db, recursive = TRUE)
    })

    gsetroot(db)

    # Test qualification
    qualified <- misha:::.gqualify_name("my.track", db)
    expect_equal(qualified, "qual@my.track")

    # Test unqualification
    unqualified <- misha:::.gunqualify_name("qual@my.track")
    expect_equal(unqualified, "my.track")

    # Unqualify already unqualified name
    unqualified2 <- misha:::.gunqualify_name("my.track")
    expect_equal(unqualified2, "my.track")
})


test_that("gintervals.rm with prefix removes correct intervals", {
    skip_if_not_installed("yaml")

    db <- create_prefixed_test_db(prefix = "irm", name = "interv_rm_test")

    withr::defer({
        unlink(db, recursive = TRUE)
    })

    dir.create(file.path(db, "tracks", "test"), showWarnings = FALSE)

    gsetroot(db)

    # Create multiple interval sets
    intervs1 <- data.frame(chrom = "chr1", start = c(0), end = c(50))
    intervs2 <- data.frame(chrom = "chr1", start = c(100), end = c(150))
    gintervals.save("irm@test.intervs1", intervs1)
    gintervals.save("irm@test.intervs2", intervs2)

    expect_true(gintervals.exists("irm@test.intervs1"))
    expect_true(gintervals.exists("irm@test.intervs2"))

    # Remove one interval set
    gintervals.rm("irm@test.intervs1", force = TRUE)

    expect_false(gintervals.exists("irm@test.intervs1"))
    expect_true(gintervals.exists("irm@test.intervs2"))

    # Clean up
    gintervals.rm("irm@test.intervs2", force = TRUE)
})


test_that("gintervals.load with prefix loads correct data", {
    skip_if_not_installed("yaml")

    db1 <- create_prefixed_test_db(prefix = "load1", name = "interv_load1")
    db2 <- create_prefixed_test_db(prefix = "load2", name = "interv_load2")

    withr::defer({
        unlink(db1, recursive = TRUE)
        unlink(db2, recursive = TRUE)
    })

    dir.create(file.path(db1, "tracks", "test"), showWarnings = FALSE)
    dir.create(file.path(db2, "tracks", "test"), showWarnings = FALSE)

    gsetroot(c(db1, db2))

    # Create intervals with different data in each database
    intervs1 <- data.frame(chrom = "chr1", start = c(0, 10), end = c(5, 15))
    intervs2 <- data.frame(chrom = "chr1", start = c(100, 200), end = c(150, 250))
    gintervals.save("load1@test.data", intervs1)
    gintervals.save("load2@test.data", intervs2)

    # Verify both intervals were created
    all_intervs <- gintervals.ls()
    expect_true("load1@test.data" %in% all_intervs)
    expect_true("load2@test.data" %in% all_intervs)

    # Clean up
    gintervals.rm("load1@test.data", force = TRUE)
    gintervals.rm("load2@test.data", force = TRUE)
})


test_that("gtrack.db returns correct database for vectorized input", {
    skip_if_not_installed("yaml")

    db1 <- create_prefixed_test_db(prefix = "vec1", name = "vector_db1")
    db2 <- create_prefixed_test_db(prefix = "vec2", name = "vector_db2")

    withr::defer({
        unlink(db1, recursive = TRUE)
        unlink(db2, recursive = TRUE)
    })

    create_test_track(db1, "track1")
    create_test_track(db1, "track2")
    create_test_track(db2, "track3")

    gsetroot(c(db1, db2))

    # Test vectorized gtrack.db
    dbs <- gtrack.db(c("vec1@track1", "vec1@track2", "vec2@track3"))
    expect_equal(length(dbs), 3)
    expect_equal(dbs[1], db1)
    expect_equal(dbs[2], db1)
    expect_equal(dbs[3], db2)
})


test_that("gintervals.db returns correct database for prefixed intervals", {
    skip_if_not_installed("yaml")

    db1 <- create_prefixed_test_db(prefix = "idb1", name = "interv_db1")
    db2 <- create_prefixed_test_db(prefix = "idb2", name = "interv_db2")

    withr::defer({
        unlink(db1, recursive = TRUE)
        unlink(db2, recursive = TRUE)
    })

    create_test_intervals(db1, "intervs1")
    create_test_intervals(db2, "intervs2")

    gsetroot(c(db1, db2))

    # Test gintervals.db
    expect_equal(gintervals.db("idb1@intervs1"), db1)
    expect_equal(gintervals.db("idb2@intervs2"), db2)
})


test_that("prefix with special characters in track name works", {
    skip_if_not_installed("yaml")

    db <- create_prefixed_test_db(prefix = "sp", name = "special_chars")

    withr::defer({
        unlink(db, recursive = TRUE)
    })

    dir.create(file.path(db, "tracks", "sub"), showWarnings = FALSE)
    dir.create(file.path(db, "tracks", "sub", "nested"), showWarnings = FALSE)

    gsetroot(db)

    # Create track with nested namespace
    intervs <- data.frame(chrom = "chr1", start = c(0), end = c(50))
    gtrack.create_sparse("sp@sub.nested.track", "Nested track", intervs, 1)

    expect_true(gtrack.exists("sp@sub.nested.track"))
    expect_equal(gtrack.db("sp@sub.nested.track"), db)

    # Clean up
    gtrack.rm("sp@sub.nested.track", force = TRUE)
})


test_that("gdb.config with NULL returns current database config", {
    skip_if_not_installed("yaml")

    db <- create_prefixed_test_db(prefix = "curr", name = "current_config")

    withr::defer({
        unlink(db, recursive = TRUE)
    })

    yaml::write_yaml(
        list(prefix = "curr", description = "Current DB"),
        file.path(db, ".misha")
    )

    gsetroot(db)

    # gdb.config with NULL should return current database config
    config <- gdb.config(NULL)
    expect_equal(config$prefix, "curr")
    expect_equal(config$description, "Current DB")
})


test_that("gdb.config with path returns correct config", {
    skip_if_not_installed("yaml")

    db1 <- create_prefixed_test_db(prefix = "path1", name = "path_config1")
    db2 <- create_prefixed_test_db(prefix = "path2", name = "path_config2")

    withr::defer({
        unlink(db1, recursive = TRUE)
        unlink(db2, recursive = TRUE)
    })

    yaml::write_yaml(
        list(prefix = "path1", description = "First DB"),
        file.path(db1, ".misha")
    )
    yaml::write_yaml(
        list(prefix = "path2", description = "Second DB"),
        file.path(db2, ".misha")
    )

    gsetroot(c(db1, db2))

    # gdb.config with path
    config1 <- gdb.config(db1)
    config2 <- gdb.config(db2)

    expect_equal(config1$description, "First DB")
    expect_equal(config2$description, "Second DB")
})


test_that("prefix parsing handles edge cases", {
    skip_if_not_installed("yaml")

    # NA input
    parsed <- misha:::.gparse_prefixed_name(NA)
    expect_null(parsed$prefix)

    # Empty string
    parsed <- misha:::.gparse_prefixed_name("")
    expect_null(parsed$prefix)
    expect_equal(parsed$name, "")

    # Single character prefix
    parsed <- misha:::.gparse_prefixed_name("a@track")
    expect_equal(parsed$prefix, "a")
    expect_equal(parsed$name, "track")

    # Prefix with underscore
    parsed <- misha:::.gparse_prefixed_name("my_db@track.name")
    expect_equal(parsed$prefix, "my_db")
    expect_equal(parsed$name, "track.name")
})


test_that("gtrack.ls with db path filtering works", {
    skip_if_not_installed("yaml")

    # Use a single prefixed db with multiple tracks
    db <- create_prefixed_test_db(prefix = "filt", name = "filter_test")

    withr::defer({
        unlink(db, recursive = TRUE)
    })

    dir.create(file.path(db, "tracks", "ns"), showWarnings = FALSE)

    gsetroot(db)

    # Create actual tracks with data
    intervs <- data.frame(chrom = "chr1", start = c(0), end = c(50))
    gtrack.create_sparse("filt@ns.track1", "Track 1", intervs, 1)
    gtrack.create_sparse("filt@ns.track2", "Track 2", intervs, 2)

    # All tracks
    all_tracks <- gtrack.ls()
    expect_equal(length(all_tracks), 2)

    # Filter by prefix
    filt_tracks <- gtrack.ls(db = "filt")
    expect_equal(length(filt_tracks), 2)
    expect_true(all(grepl("^filt@", filt_tracks)))

    # Filter by database path
    path_tracks <- gtrack.ls(db = db)
    expect_equal(length(path_tracks), 2)

    # Clean up
    gtrack.rm("filt@ns.track1", force = TRUE)
    gtrack.rm("filt@ns.track2", force = TRUE)
})


test_that("track rename within prefixed database works", {
    skip_if_not_installed("yaml")

    db <- create_prefixed_test_db(prefix = "ren", name = "rename_test")

    withr::defer({
        unlink(db, recursive = TRUE)
    })

    dir.create(file.path(db, "tracks", "ns1"), showWarnings = FALSE)
    dir.create(file.path(db, "tracks", "ns2"), showWarnings = FALSE)

    gsetroot(db)

    # Create a track
    intervs <- data.frame(chrom = "chr1", start = c(0), end = c(50))
    gtrack.create_sparse("ren@ns1.original", "Original track", intervs, 1)

    expect_true(gtrack.exists("ren@ns1.original"))

    # Rename within same database (different namespace)
    gtrack.mv("ren@ns1.original", "ren@ns2.renamed")

    expect_false(gtrack.exists("ren@ns1.original"))
    expect_true(gtrack.exists("ren@ns2.renamed"))
    expect_equal(gtrack.db("ren@ns2.renamed"), db)

    # Clean up
    gtrack.rm("ren@ns2.renamed", force = TRUE)
})


test_that("gis_global_db correctly identifies database types", {
    skip_if_not_installed("yaml")

    global_db <- create_prefixed_test_db(prefix = NULL, name = "is_global_test1")
    prefixed_db <- create_prefixed_test_db(prefix = "isg", name = "is_global_test2")

    withr::defer({
        unlink(global_db, recursive = TRUE)
        unlink(prefixed_db, recursive = TRUE)
    })

    gsetroot(c(global_db, prefixed_db))

    expect_true(misha:::.gis_global_db(global_db))
    expect_false(misha:::.gis_global_db(prefixed_db))
})


test_that("prefix to db and db to prefix lookups work", {
    skip_if_not_installed("yaml")

    db <- create_prefixed_test_db(prefix = "lookup", name = "lookup_test")

    withr::defer({
        unlink(db, recursive = TRUE)
    })

    gsetroot(db)

    # Test prefix to db lookup
    result <- misha:::.gprefix_to_db("lookup")
    expect_equal(result, db)

    # Test unknown prefix
    result <- misha:::.gprefix_to_db("unknown")
    expect_null(result)

    # Test db to prefix lookup
    result <- misha:::.gdb_to_prefix(db)
    expect_equal(result, "lookup")
})


test_that("empty database with prefix works", {
    skip_if_not_installed("yaml")

    db <- create_prefixed_test_db(prefix = "empty", name = "empty_db")

    withr::defer({
        unlink(db, recursive = TRUE)
    })

    # Database with prefix but no tracks
    gsetroot(db)

    prefixes <- gdb.prefixes()
    expect_equal(names(prefixes), "empty")

    tracks <- gtrack.ls()
    expect_equal(length(tracks), 0)

    intervs <- gintervals.ls()
    expect_equal(length(intervs), 0)
})


test_that("gtrack.dbs returns all databases with track", {
    skip_if_not_installed("yaml")

    # Create two prefixed databases with same track name (different prefixes mean different qualified names)
    db1 <- create_prefixed_test_db(prefix = "dbs1", name = "dbs_test1")
    db2 <- create_prefixed_test_db(prefix = "dbs2", name = "dbs_test2")

    withr::defer({
        unlink(db1, recursive = TRUE)
        unlink(db2, recursive = TRUE)
    })

    create_test_track(db1, "shared.track")
    create_test_track(db2, "unique.track")

    gsetroot(c(db1, db2))

    # Each prefixed track exists only in one db
    dbs1 <- gtrack.dbs("dbs1@shared.track")
    expect_equal(length(dbs1), 1)
    expect_equal(unname(dbs1), db1)

    dbs2 <- gtrack.dbs("dbs2@unique.track")
    expect_equal(length(dbs2), 1)
    expect_equal(unname(dbs2), db2)
})


test_that("invalid .misha YAML file is handled gracefully", {
    skip_if_not_installed("yaml")

    db <- create_prefixed_test_db(prefix = NULL, name = "invalid_yaml")

    withr::defer({
        unlink(db, recursive = TRUE)
    })

    # Write invalid YAML
    writeLines("invalid: yaml: content: [", file.path(db, ".misha"))

    create_test_track(db, "my.track")

    # Should work with warning about invalid YAML
    expect_warning(gsetroot(db), "Failed to read")

    # Track should still be accessible as global
    tracks <- gtrack.ls()
    expect_true("my.track" %in% tracks)
})


test_that("prefix suggestions for intervals work", {
    skip_if_not_installed("yaml")

    global_db <- create_prefixed_test_db(prefix = NULL, name = "suggest_global")
    prefixed_db <- create_prefixed_test_db(prefix = "sug", name = "suggest_prefixed")

    withr::defer({
        unlink(global_db, recursive = TRUE)
        unlink(prefixed_db, recursive = TRUE)
    })

    # Create intervals only in prefixed database
    create_test_intervals(prefixed_db, "only.here")

    gsetroot(c(global_db, prefixed_db))

    # Test suggestion function directly
    suggestion <- misha:::.gsuggest_prefixed_intervals("only.here")
    expect_match(suggestion, "Did you mean")
    expect_match(suggestion, "sug@only.here")
})


test_that("global and prefixed databases order doesn't affect resolution", {
    skip_if_not_installed("yaml")

    global_db <- create_prefixed_test_db(prefix = NULL, name = "order_global")
    prefixed_db <- create_prefixed_test_db(prefix = "ord", name = "order_prefixed")

    withr::defer({
        unlink(global_db, recursive = TRUE)
        unlink(prefixed_db, recursive = TRUE)
    })

    create_test_track(global_db, "global.only")
    create_test_track(prefixed_db, "prefixed.only")

    # Test both orders
    gsetroot(c(global_db, prefixed_db))
    tracks1 <- sort(gtrack.ls())

    gsetroot(c(prefixed_db, global_db))
    tracks2 <- sort(gtrack.ls())

    # Same tracks should be available regardless of order
    expect_equal(tracks1, tracks2)
    expect_true("global.only" %in% tracks1)
    expect_true("ord@prefixed.only" %in% tracks1)
})


# Tests for dangerously_ignore_db_prefixes option

test_that("dangerously_ignore_db_prefixes=FALSE uses prefixes (default)", {
    skip_if_not_installed("yaml")

    prefixed_db <- create_prefixed_test_db(prefix = "ign", name = "ignore_test1")

    withr::defer({
        unlink(prefixed_db, recursive = TRUE)
    })

    create_test_track(prefixed_db, "my.track")
    gsetroot(prefixed_db)

    # Default: track requires prefix
    tracks <- gtrack.ls()
    expect_true("ign@my.track" %in% tracks)
    expect_false("my.track" %in% tracks)
})


test_that("dangerously_ignore_db_prefixes=TRUE ignores prefixes", {
    skip_if_not_installed("yaml")

    prefixed_db <- create_prefixed_test_db(prefix = "ign", name = "ignore_test2")

    withr::defer({
        unlink(prefixed_db, recursive = TRUE)
        options(misha.dangerously_ignore_db_prefixes = NULL)
    })

    create_test_track(prefixed_db, "my.track")

    # Set option and reload
    options(misha.dangerously_ignore_db_prefixes = TRUE)
    gsetroot(prefixed_db)

    # Track should be accessible without prefix
    tracks <- gtrack.ls()
    expect_true("my.track" %in% tracks)
    expect_false("ign@my.track" %in% tracks)
})


test_that("dangerously_ignore_db_prefixes=TRUE treats @ as literal character", {
    skip_if_not_installed("yaml")

    prefixed_db <- create_prefixed_test_db(prefix = "ign", name = "ignore_test3")

    withr::defer({
        unlink(prefixed_db, recursive = TRUE)
        options(misha.dangerously_ignore_db_prefixes = NULL)
    })

    dir.create(file.path(prefixed_db, "tracks", "test"), showWarnings = FALSE)

    options(misha.dangerously_ignore_db_prefixes = TRUE)
    gsetroot(prefixed_db)

    # Create track - prefix syntax should not work
    intervs <- data.frame(chrom = "chr1", start = c(0), end = c(50))
    gtrack.create_sparse("test.track", "Test", intervs, 1)

    # Track accessible by bare name
    expect_true(gtrack.exists("test.track"))
    expect_true("test.track" %in% gtrack.ls())

    # Clean up
    gtrack.rm("test.track", force = TRUE)
})


test_that("dangerously_ignore_db_prefixes with multiple databases uses last wins", {
    skip_if_not_installed("yaml")

    db1 <- create_prefixed_test_db(prefix = "first", name = "ignore_multi1")
    db2 <- create_prefixed_test_db(prefix = "second", name = "ignore_multi2")

    withr::defer({
        unlink(db1, recursive = TRUE)
        unlink(db2, recursive = TRUE)
        options(misha.dangerously_ignore_db_prefixes = NULL)
    })

    # Create track with same name in both databases
    create_test_track(db1, "shared.track")
    create_test_track(db2, "shared.track")

    options(misha.dangerously_ignore_db_prefixes = TRUE)
    gsetroot(c(db1, db2))

    # Should see unprefixed name only
    tracks <- gtrack.ls()
    expect_true("shared.track" %in% tracks)
    expect_false("first@shared.track" %in% tracks)
    expect_false("second@shared.track" %in% tracks)

    # Last connected database wins
    expect_equal(gtrack.db("shared.track"), db2)
})


test_that("dangerously_ignore_db_prefixes respects option scope", {
    skip_if_not_installed("yaml")

    prefixed_db <- create_prefixed_test_db(prefix = "scope", name = "ignore_scope")

    withr::defer({
        unlink(prefixed_db, recursive = TRUE)
    })

    create_test_track(prefixed_db, "my.track")

    # Default behavior
    gsetroot(prefixed_db)
    expect_true("scope@my.track" %in% gtrack.ls())

    # With local option
    withr::local_options(misha.dangerously_ignore_db_prefixes = TRUE)
    gsetroot(prefixed_db)
    expect_true("my.track" %in% gtrack.ls())

    # After local scope, back to default
    gsetroot(prefixed_db)
    expect_true("scope@my.track" %in% gtrack.ls())
})


# Tests for dangerously_allow_track_collision option

test_that("dangerously_allow_track_collision=FALSE rejects duplicate prefixes (default)", {
    skip_if_not_installed("yaml")

    db1 <- create_prefixed_test_db(prefix = "dup", name = "collision_test1")
    db2 <- create_prefixed_test_db(prefix = "dup", name = "collision_test2")

    withr::defer({
        unlink(db1, recursive = TRUE)
        unlink(db2, recursive = TRUE)
    })

    # Default: duplicate prefixes should error
    expect_error(gsetroot(c(db1, db2)), "Duplicate prefix")
})


test_that("dangerously_allow_track_collision=TRUE allows duplicate prefixes", {
    skip_if_not_installed("yaml")

    db1 <- create_prefixed_test_db(prefix = "dup", name = "collision_test3")
    db2 <- create_prefixed_test_db(prefix = "dup", name = "collision_test4")

    withr::defer({
        unlink(db1, recursive = TRUE)
        unlink(db2, recursive = TRUE)
        options(misha.dangerously_allow_track_collision = NULL)
    })

    # With option: duplicate prefixes allowed
    options(misha.dangerously_allow_track_collision = TRUE)
    expect_silent(gsetroot(c(db1, db2)))
})


test_that("dangerously_allow_track_collision uses last wins for duplicates", {
    skip_if_not_installed("yaml")

    db1 <- create_prefixed_test_db(prefix = "dup", name = "collision_test5")
    db2 <- create_prefixed_test_db(prefix = "dup", name = "collision_test6")

    withr::defer({
        unlink(db1, recursive = TRUE)
        unlink(db2, recursive = TRUE)
        options(misha.dangerously_allow_track_collision = NULL)
    })

    create_test_track(db1, "track1")
    create_test_track(db2, "track2")

    options(misha.dangerously_allow_track_collision = TRUE)
    gsetroot(c(db1, db2))

    # Both tracks accessible with same prefix
    expect_true("dup@track1" %in% gtrack.ls())
    expect_true("dup@track2" %in% gtrack.ls())

    # track2 from second db
    expect_equal(gtrack.db("dup@track2"), db2)
})


test_that("dangerously_allow_track_collision respects option scope", {
    skip_if_not_installed("yaml")

    db1 <- create_prefixed_test_db(prefix = "dup", name = "collision_scope1")
    db2 <- create_prefixed_test_db(prefix = "dup", name = "collision_scope2")

    withr::defer({
        unlink(db1, recursive = TRUE)
        unlink(db2, recursive = TRUE)
    })

    # Default: error
    expect_error(gsetroot(c(db1, db2)), "Duplicate prefix")

    # With local option: allowed
    withr::local_options(misha.dangerously_allow_track_collision = TRUE)
    expect_silent(gsetroot(c(db1, db2)))

    # After local scope, back to error
    expect_error(gsetroot(c(db1, db2)), "Duplicate prefix")
})


# Tests for both options combined

test_that("both dangerous options can be used together", {
    skip_if_not_installed("yaml")

    db1 <- create_prefixed_test_db(prefix = "a", name = "both_test1")
    db2 <- create_prefixed_test_db(prefix = "a", name = "both_test2")

    withr::defer({
        unlink(db1, recursive = TRUE)
        unlink(db2, recursive = TRUE)
        options(misha.dangerously_ignore_db_prefixes = NULL)
        options(misha.dangerously_allow_track_collision = NULL)
    })

    create_test_track(db1, "track.x")
    create_test_track(db2, "track.y")

    # Both options enabled
    options(misha.dangerously_ignore_db_prefixes = TRUE)
    options(misha.dangerously_allow_track_collision = TRUE)
    gsetroot(c(db1, db2))

    # Tracks accessible without prefix
    tracks <- gtrack.ls()
    expect_true("track.x" %in% tracks)
    expect_true("track.y" %in% tracks)
    expect_false(any(grepl("@", tracks)))
})


test_that("dangerous options work independently", {
    skip_if_not_installed("yaml")

    db1 <- create_prefixed_test_db(prefix = "ind1", name = "independent1")
    db2 <- create_prefixed_test_db(prefix = "ind2", name = "independent2")

    withr::defer({
        unlink(db1, recursive = TRUE)
        unlink(db2, recursive = TRUE)
        options(misha.dangerously_ignore_db_prefixes = NULL)
    })

    create_test_track(db1, "track1")
    create_test_track(db2, "track2")

    # Only ignore_db_prefixes enabled
    options(misha.dangerously_ignore_db_prefixes = TRUE)
    gsetroot(c(db1, db2))

    # Different prefixes still work, just without prefix in names
    tracks <- gtrack.ls()
    expect_true("track1" %in% tracks)
    expect_true("track2" %in% tracks)
})
