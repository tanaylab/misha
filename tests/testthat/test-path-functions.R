create_isolated_test_db()

test_that("gintervals.path returns correct path for single interval", {
    interv_name <- "annotations"
    path <- gintervals.path(interv_name)

    # Should return a character vector of length 1
    expect_type(path, "character")
    expect_length(path, 1)

    # Path should end with .interv
    expect_true(endsWith(path, ".interv"))

    # Path should contain the interval name
    expect_true(grepl(interv_name, path, fixed = TRUE))

    # Path should match expected format: GWD/interval_name.interv
    expected_path <- sprintf(
        "%s.interv",
        paste(get("GWD", envir = .misha), gsub("\\.", "/", interv_name), sep = "/")
    )
    expect_equal(path, expected_path)
})

test_that("gintervals.path works with vectorized input", {
    interv_names <- c("annotations", "coding")
    paths <- gintervals.path(interv_names)

    # Should return a character vector of correct length
    expect_type(paths, "character")
    expect_length(paths, length(interv_names))

    # Each path should end with .interv
    expect_true(all(endsWith(paths, ".interv")))

    # Verify each path matches expected format
    for (i in seq_along(interv_names)) {
        expected_path <- sprintf(
            "%s.interv",
            paste(get("GWD", envir = .misha), gsub("\\.", "/", interv_names[i]), sep = "/")
        )
        expect_equal(paths[i], expected_path)
    }
})

test_that("gintervals.path handles nested interval names correctly", {
    # Test with nested names (dots in the name)
    interv_names <- c("test.interv1", "test.interv2", "subdir.nested.interval")
    paths <- gintervals.path(interv_names)

    expect_type(paths, "character")
    expect_length(paths, length(interv_names))

    # Dots should be replaced with slashes in the path
    expect_true(grepl("test/interv1", paths[1], fixed = TRUE))
    expect_true(grepl("test/interv2", paths[2], fixed = TRUE))
    expect_true(grepl("subdir/nested/interval", paths[3], fixed = TRUE))

    # All should still end with .interv
    expect_true(all(endsWith(paths, ".interv")))
})

test_that("gintervals.path handles empty vector", {
    paths <- gintervals.path(character(0))
    expect_type(paths, "character")
    expect_length(paths, 0)
})

test_that("gintervals.path matches format used in .gintervals.load_file", {
    interv_name <- "annotations"
    path_from_function <- gintervals.path(interv_name)

    # The path should match what .gintervals.load_file would construct
    expected_intervfname <- sprintf(
        "%s.interv",
        paste(get("GWD", envir = .misha), gsub("\\.", "/", interv_name), sep = "/")
    )

    expect_equal(path_from_function, expected_intervfname)
})

test_that("gtrack.path returns correct path for single track", {
    tracks <- gtrack.ls()
    skip_if(length(tracks) == 0, "No tracks available for testing")

    track_name <- tracks[1]
    path <- gtrack.path(track_name)

    # Should return a character vector of length 1
    expect_type(path, "character")
    expect_length(path, 1)

    # Path should end with .track
    expect_true(endsWith(path, ".track"))

    # Path should contain the track name (dots replaced with slashes)
    track_path_part <- gsub("\\.", "/", track_name)
    expect_true(grepl(track_path_part, path, fixed = TRUE))

    # Path should match expected format: GWD/track_name.track
    expected_path <- sprintf(
        "%s.track",
        paste(get("GWD", envir = .misha), gsub("\\.", "/", track_name), sep = "/")
    )
    expect_equal(path, expected_path)
})

test_that("gtrack.path works with vectorized input", {
    tracks <- gtrack.ls()
    skip_if(length(tracks) < 2, "Need at least 2 tracks for vectorized testing")

    track_names <- tracks[1:min(3, length(tracks))]
    paths <- gtrack.path(track_names)

    # Should return a character vector of correct length
    expect_type(paths, "character")
    expect_length(paths, length(track_names))

    # Each path should end with .track
    expect_true(all(endsWith(paths, ".track")))

    # Verify each path matches expected format
    for (i in seq_along(track_names)) {
        expected_path <- sprintf(
            "%s.track",
            paste(get("GWD", envir = .misha), gsub("\\.", "/", track_names[i]), sep = "/")
        )
        expect_equal(paths[i], expected_path)
    }
})

test_that("gtrack.path handles nested track names correctly", {
    # Create a test track with nested name if possible
    tracks <- gtrack.ls()
    skip_if(length(tracks) == 0, "No tracks available for testing")

    # Test with existing track that might have dots
    track_with_dots <- tracks[grepl("\\.", tracks)][1]
    if (!is.na(track_with_dots)) {
        path <- gtrack.path(track_with_dots)

        # Dots should be replaced with slashes in the path
        track_path_part <- gsub("\\.", "/", track_with_dots)
        expect_true(grepl(track_path_part, path, fixed = TRUE))
        expect_true(endsWith(path, ".track"))
    }
})

test_that("gtrack.path handles empty vector", {
    paths <- gtrack.path(character(0))
    expect_type(paths, "character")
    expect_length(paths, 0)
})

test_that("gtrack.path matches format used in .track_dir", {
    tracks <- gtrack.ls()
    skip_if(length(tracks) == 0, "No tracks available for testing")

    track_name <- tracks[1]
    path_from_function <- gtrack.path(track_name)

    # The path should match what .track_dir would construct
    # (We can't call .track_dir directly as it's not exported, but we can construct it manually)
    expected_path <- sprintf(
        "%s.track",
        paste(get("GWD", envir = .misha), gsub("\\.", "/", track_name), sep = "/")
    )

    expect_equal(path_from_function, expected_path)
})

test_that("gtrack.path and gintervals.path return different extensions", {
    tracks <- gtrack.ls()
    intervals <- gintervals.ls()

    skip_if(
        length(tracks) == 0 || length(intervals) == 0,
        "Need both tracks and intervals for comparison"
    )

    track_path <- gtrack.path(tracks[1])
    interv_path <- gintervals.path(intervals[1])

    # Track paths should end with .track
    expect_true(endsWith(track_path, ".track"))

    # Interval paths should end with .interv
    expect_true(endsWith(interv_path, ".interv"))
})

test_that("both path functions handle large vectors efficiently", {
    tracks <- gtrack.ls()
    intervals <- gintervals.ls()

    # Test with all available tracks/intervals (could be many)
    if (length(tracks) > 0) {
        all_track_paths <- gtrack.path(tracks)
        expect_length(all_track_paths, length(tracks))
        expect_true(all(endsWith(all_track_paths, ".track")))
    }

    if (length(intervals) > 0) {
        all_interv_paths <- gintervals.path(intervals)
        expect_length(all_interv_paths, length(intervals))
        expect_true(all(endsWith(all_interv_paths, ".interv")))
    }
})

test_that("path functions work with big intervals sets", {
    # Test that path functions work with big intervals sets
    big_intervs <- gintervals.ls()
    big_intervs <- big_intervs[grep("^bigintervs", big_intervs)]

    skip_if(length(big_intervs) == 0, "No big intervals sets available for testing")

    paths <- gintervals.path(big_intervs[1])
    expect_type(paths, "character")
    expect_length(paths, 1)
    expect_true(endsWith(paths, ".interv"))

    # Should work with multiple big intervals
    if (length(big_intervs) > 1) {
        paths <- gintervals.path(big_intervs[1:min(2, length(big_intervs))])
        expect_length(paths, min(2, length(big_intervs)))
    }
})
