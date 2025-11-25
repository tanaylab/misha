create_isolated_test_db()

# Helper functions for manual computation
# These use the track's native iterator to match C++ behavior
manual_exists <- function(track_name, interval_row) {
    interval <- gintervals(interval_row$chrom, interval_row$start, interval_row$end)
    values <- gextract(track_name, interval, iterator = track_name)
    track_vals <- values[[track_name]]
    if (any(!is.na(track_vals))) 1.0 else 0.0
}

manual_exists_vec <- function(track_name, intervals) {
    vapply(seq_len(nrow(intervals)), function(i) {
        manual_exists(track_name, intervals[i, ])
    }, numeric(1))
}

manual_size <- function(track_name, interval_row) {
    interval <- gintervals(interval_row$chrom, interval_row$start, interval_row$end)
    values <- gextract(track_name, interval, iterator = track_name)
    track_vals <- values[[track_name]]
    sum(!is.na(track_vals))
}

manual_size_vec <- function(track_name, intervals) {
    vapply(seq_len(nrow(intervals)), function(i) {
        manual_size(track_name, intervals[i, ])
    }, numeric(1))
}

manual_first <- function(track_name, interval_row) {
    interval <- gintervals(interval_row$chrom, interval_row$start, interval_row$end)
    values <- gextract(track_name, interval, iterator = track_name)
    track_vals <- values[[track_name]]
    valid <- !is.na(track_vals)
    if (!any(valid)) {
        return(NA_real_)
    }
    track_vals[which(valid)[1]]
}

manual_first_vec <- function(track_name, intervals) {
    vapply(seq_len(nrow(intervals)), function(i) {
        manual_first(track_name, intervals[i, ])
    }, numeric(1))
}

manual_last <- function(track_name, interval_row) {
    interval <- gintervals(interval_row$chrom, interval_row$start, interval_row$end)
    values <- gextract(track_name, interval, iterator = track_name)
    track_vals <- values[[track_name]]
    valid <- !is.na(track_vals)
    if (!any(valid)) {
        return(NA_real_)
    }
    track_vals[rev(which(valid))[1]]
}

manual_last_vec <- function(track_name, intervals) {
    vapply(seq_len(nrow(intervals)), function(i) {
        manual_last(track_name, intervals[i, ])
    }, numeric(1))
}

manual_first_pos <- function(track_name, interval_row, relative = FALSE, relative_ref = NULL) {
    interval <- gintervals(interval_row$chrom, interval_row$start, interval_row$end)
    values <- gextract(track_name, interval, iterator = track_name)
    track_vals <- values[[track_name]]
    valid <- !is.na(track_vals)
    if (!any(valid)) {
        return(NA_real_)
    }
    first_pos <- values$start[which(valid)[1]]
    if (relative) {
        ref_start <- if (!is.null(relative_ref)) relative_ref else interval_row$start
        first_pos - ref_start
    } else {
        first_pos
    }
}

manual_first_pos_vec <- function(track_name, intervals, relative = FALSE, relative_ref = NULL) {
    vapply(seq_len(nrow(intervals)), function(i) {
        ref <- if (is.null(relative_ref)) NULL else relative_ref[i]
        manual_first_pos(track_name, intervals[i, ], relative = relative, relative_ref = ref)
    }, numeric(1))
}

manual_last_pos <- function(track_name, interval_row, relative = FALSE, relative_ref = NULL) {
    interval <- gintervals(interval_row$chrom, interval_row$start, interval_row$end)
    values <- gextract(track_name, interval, iterator = track_name)
    track_vals <- values[[track_name]]
    valid <- !is.na(track_vals)
    if (!any(valid)) {
        return(NA_real_)
    }
    last_pos <- values$start[rev(which(valid))[1]]
    if (relative) {
        ref_start <- if (!is.null(relative_ref)) relative_ref else interval_row$start
        last_pos - ref_start
    } else {
        last_pos
    }
}

manual_last_pos_vec <- function(track_name, intervals, relative = FALSE, relative_ref = NULL) {
    vapply(seq_len(nrow(intervals)), function(i) {
        ref <- if (is.null(relative_ref)) NULL else relative_ref[i]
        manual_last_pos(track_name, intervals[i, ], relative = relative, relative_ref = ref)
    }, numeric(1))
}

apply_interval_shift <- function(intervals, sshift = 0L, eshift = 0L) {
    shifted <- intervals
    shifted$start <- shifted$start + sshift
    shifted$end <- shifted$end + eshift
    shifted
}

# Tests for 'exists'
test_that("exists returns correct values for sparse track", {
    track_name <- "test.sparse"
    vtrack_name <- "vt_exists_sparse"
    intervals <- rbind(
        gintervals(1, 0, 300),
        gintervals(1, 600, 1000),
        gintervals(1, 1200, 1500)
    )

    gvtrack.create(vtrack_name, track_name, "exists")
    on.exit(gvtrack.rm(vtrack_name), add = TRUE)

    vt_res <- gextract(vtrack_name, intervals, iterator = intervals)
    manual <- manual_exists_vec(track_name, intervals)
    expect_equal(vt_res[[vtrack_name]], manual)
    # Check that at least some intervals have data
    expect_true(any(vt_res[[vtrack_name]] == 1))
})

test_that("exists returns 0 when no values exist", {
    track_name <- "test.sparse"
    vtrack_name <- "vt_exists_empty"
    # Create interval where track has no data
    intervals <- gintervals(1, 10000, 10100)

    gvtrack.create(vtrack_name, track_name, "exists")
    on.exit(gvtrack.rm(vtrack_name), add = TRUE)

    vt_res <- gextract(vtrack_name, intervals, iterator = intervals)
    expect_equal(vt_res[[vtrack_name]], 0)
})

# Tests for 'size'
test_that("size returns correct count for dense track", {
    track_name <- "test.fixedbin"
    vtrack_name <- "vt_size_dense"
    intervals <- rbind(
        gintervals(1, 0, 200),
        gintervals(1, 500, 900),
        gintervals(1, 1500, 2100)
    )

    gvtrack.create(vtrack_name, track_name, "size")
    on.exit(gvtrack.rm(vtrack_name), add = TRUE)

    vt_res <- gextract(vtrack_name, intervals, iterator = intervals)
    manual <- manual_size_vec(track_name, intervals)
    expect_equal(vt_res[[vtrack_name]], manual)
})

test_that("size returns 0 when no values exist", {
    track_name <- "test.sparse"
    vtrack_name <- "vt_size_empty"
    intervals <- gintervals(1, 10000, 10100)

    gvtrack.create(vtrack_name, track_name, "size")
    on.exit(gvtrack.rm(vtrack_name), add = TRUE)

    vt_res <- gextract(vtrack_name, intervals, iterator = intervals)
    expect_equal(vt_res[[vtrack_name]], 0)
})

# Tests for 'first'
test_that("first returns first value in sparse track", {
    track_name <- "test.sparse"
    vtrack_name <- "vt_first_sparse"
    intervals <- rbind(
        gintervals(1, 0, 300),
        gintervals(1, 600, 1000),
        gintervals(1, 1200, 1500)
    )

    gvtrack.create(vtrack_name, track_name, "first")
    on.exit(gvtrack.rm(vtrack_name), add = TRUE)

    vt_res <- gextract(vtrack_name, intervals, iterator = intervals)
    manual <- manual_first_vec(track_name, intervals)
    expect_equal(vt_res[[vtrack_name]], manual)
})

test_that("first returns NA when no values exist", {
    track_name <- "test.sparse"
    vtrack_name <- "vt_first_empty"
    intervals <- gintervals(1, 10000, 10100)

    gvtrack.create(vtrack_name, track_name, "first")
    on.exit(gvtrack.rm(vtrack_name), add = TRUE)

    vt_res <- gextract(vtrack_name, intervals, iterator = intervals)
    expect_true(is.na(vt_res[[vtrack_name]]))
})

# Tests for 'last'
test_that("last returns last value in sparse track", {
    track_name <- "test.sparse"
    vtrack_name <- "vt_last_sparse"
    intervals <- rbind(
        gintervals(1, 0, 300),
        gintervals(1, 600, 1000),
        gintervals(1, 1200, 1500)
    )

    gvtrack.create(vtrack_name, track_name, "last")
    on.exit(gvtrack.rm(vtrack_name), add = TRUE)

    vt_res <- gextract(vtrack_name, intervals, iterator = intervals)
    manual <- manual_last_vec(track_name, intervals)
    expect_equal(vt_res[[vtrack_name]], manual)
})

test_that("last returns NA when no values exist", {
    track_name <- "test.sparse"
    vtrack_name <- "vt_last_empty"
    intervals <- gintervals(1, 10000, 10100)

    gvtrack.create(vtrack_name, track_name, "last")
    on.exit(gvtrack.rm(vtrack_name), add = TRUE)

    vt_res <- gextract(vtrack_name, intervals, iterator = intervals)
    expect_true(is.na(vt_res[[vtrack_name]]))
})

# Tests for 'first.pos.abs'
test_that("first.pos.abs returns correct position for sparse track", {
    track_name <- "test.sparse"
    vtrack_name <- "vt_first_pos_abs_sparse"
    intervals <- rbind(
        gintervals(1, 0, 300),
        gintervals(1, 600, 1000),
        gintervals(1, 1200, 1500)
    )

    gvtrack.create(vtrack_name, track_name, "first.pos.abs")
    on.exit(gvtrack.rm(vtrack_name), add = TRUE)

    vt_res <- gextract(vtrack_name, intervals, iterator = intervals)
    manual <- manual_first_pos_vec(track_name, intervals, relative = FALSE)
    expect_equal(vt_res[[vtrack_name]], manual)
})

test_that("first.pos.abs honors iterator shifts", {
    track_name <- "test.fixedbin"
    vtrack_name <- "vt_first_pos_abs_shift"
    intervals <- rbind(
        gintervals(1, 100, 250),
        gintervals(1, 700, 900)
    )
    sshift <- -50L
    eshift <- 75L

    gvtrack.create(vtrack_name, track_name, "first.pos.abs")
    gvtrack.iterator(vtrack_name, sshift = sshift, eshift = eshift)
    on.exit(gvtrack.rm(vtrack_name), add = TRUE)

    vt_res <- gextract(vtrack_name, intervals, iterator = intervals)
    shifted <- apply_interval_shift(intervals, sshift = sshift, eshift = eshift)
    manual <- manual_first_pos_vec(track_name, shifted, relative = FALSE)
    expect_equal(vt_res[[vtrack_name]], manual)
})

# Tests for 'first.pos.relative'
test_that("first.pos.relative returns correct relative position", {
    track_name <- "test.sparse"
    vtrack_name <- "vt_first_pos_rel_sparse"
    intervals <- rbind(
        gintervals(1, 0, 300),
        gintervals(1, 600, 1000),
        gintervals(1, 1200, 1500)
    )

    gvtrack.create(vtrack_name, track_name, "first.pos.relative")
    on.exit(gvtrack.rm(vtrack_name), add = TRUE)

    vt_res <- gextract(vtrack_name, intervals, iterator = intervals)
    manual <- manual_first_pos_vec(track_name, intervals, relative = TRUE)
    expect_equal(vt_res[[vtrack_name]], manual)
})

test_that("first.pos.relative honors iterator shifts", {
    track_name <- "test.fixedbin"
    vtrack_name <- "vt_first_pos_rel_shift"
    intervals <- rbind(
        gintervals(1, 200, 320),
        gintervals(1, 600, 760)
    )
    sshift <- -30L
    eshift <- 60L

    gvtrack.create(vtrack_name, track_name, "first.pos.relative")
    gvtrack.iterator(vtrack_name, sshift = sshift, eshift = eshift)
    on.exit(gvtrack.rm(vtrack_name), add = TRUE)

    vt_res <- gextract(vtrack_name, intervals, iterator = intervals)
    shifted <- apply_interval_shift(intervals, sshift = sshift, eshift = eshift)
    # Relative positions are relative to the SHIFTED interval start (like max.pos.relative)
    manual <- manual_first_pos_vec(track_name, shifted, relative = TRUE, relative_ref = shifted$start)
    expect_equal(vt_res[[vtrack_name]], manual)
})

# Tests for 'last.pos.abs'
test_that("last.pos.abs returns correct position for sparse track", {
    track_name <- "test.sparse"
    vtrack_name <- "vt_last_pos_abs_sparse"
    intervals <- rbind(
        gintervals(1, 0, 300),
        gintervals(1, 600, 1000),
        gintervals(1, 1200, 1500)
    )

    gvtrack.create(vtrack_name, track_name, "last.pos.abs")
    on.exit(gvtrack.rm(vtrack_name), add = TRUE)

    vt_res <- gextract(vtrack_name, intervals, iterator = intervals)
    manual <- manual_last_pos_vec(track_name, intervals, relative = FALSE)
    expect_equal(vt_res[[vtrack_name]], manual)
})

test_that("last.pos.abs honors iterator shifts", {
    track_name <- "test.fixedbin"
    vtrack_name <- "vt_last_pos_abs_shift"
    intervals <- rbind(
        gintervals(1, 100, 250),
        gintervals(1, 700, 900)
    )
    sshift <- -50L
    eshift <- 75L

    gvtrack.create(vtrack_name, track_name, "last.pos.abs")
    gvtrack.iterator(vtrack_name, sshift = sshift, eshift = eshift)
    on.exit(gvtrack.rm(vtrack_name), add = TRUE)

    vt_res <- gextract(vtrack_name, intervals, iterator = intervals)
    shifted <- apply_interval_shift(intervals, sshift = sshift, eshift = eshift)
    manual <- manual_last_pos_vec(track_name, shifted, relative = FALSE)
    expect_equal(vt_res[[vtrack_name]], manual)
})

# Tests for 'last.pos.relative'
test_that("last.pos.relative returns correct relative position", {
    track_name <- "test.sparse"
    vtrack_name <- "vt_last_pos_rel_sparse"
    intervals <- rbind(
        gintervals(1, 0, 300),
        gintervals(1, 600, 1000),
        gintervals(1, 1200, 1500)
    )

    gvtrack.create(vtrack_name, track_name, "last.pos.relative")
    on.exit(gvtrack.rm(vtrack_name), add = TRUE)

    vt_res <- gextract(vtrack_name, intervals, iterator = intervals)
    manual <- manual_last_pos_vec(track_name, intervals, relative = TRUE)
    expect_equal(vt_res[[vtrack_name]], manual)
})

test_that("last.pos.relative honors iterator shifts", {
    track_name <- "test.fixedbin"
    vtrack_name <- "vt_last_pos_rel_shift"
    intervals <- rbind(
        gintervals(1, 200, 320),
        gintervals(1, 600, 760)
    )
    sshift <- -30L
    eshift <- 60L

    gvtrack.create(vtrack_name, track_name, "last.pos.relative")
    gvtrack.iterator(vtrack_name, sshift = sshift, eshift = eshift)
    on.exit(gvtrack.rm(vtrack_name), add = TRUE)

    vt_res <- gextract(vtrack_name, intervals, iterator = intervals)
    shifted <- apply_interval_shift(intervals, sshift = sshift, eshift = eshift)
    # Relative positions are relative to the SHIFTED interval start (like max.pos.relative)
    manual <- manual_last_pos_vec(track_name, shifted, relative = TRUE, relative_ref = shifted$start)
    expect_equal(vt_res[[vtrack_name]], manual)
})

# Tests for 'sample' - statistical properties
test_that("sample returns valid value from interval", {
    track_name <- "test.sparse"
    vtrack_name <- "vt_sample_sparse"
    intervals <- gintervals(1, 0, 1000)

    gvtrack.create(vtrack_name, track_name, "sample")
    on.exit(gvtrack.rm(vtrack_name), add = TRUE)

    # Get all values in the interval
    all_vals <- gextract(track_name, intervals, iterator = 1)
    valid_vals <- all_vals[[track_name]][!is.na(all_vals[[track_name]])]

    # Test that sampled values are valid
    set.seed(123)
    samples <- sapply(1:10, function(i) {
        gextract(vtrack_name, intervals, iterator = intervals)[[vtrack_name]]
    })

    # All samples should be in the valid set
    expect_true(all(samples %in% valid_vals))
    # At least some samples should be non-NA
    expect_true(any(!is.na(samples)))
})

test_that("sample is reproducible with set.seed", {
    track_name <- "test.sparse"
    vtrack_name <- "vt_sample_seed"
    intervals <- gintervals(1, 0, 1000)

    gvtrack.create(vtrack_name, track_name, "sample")
    on.exit(gvtrack.rm(vtrack_name), add = TRUE)

    set.seed(456)
    result1 <- gextract(vtrack_name, intervals, iterator = intervals)[[vtrack_name]]

    set.seed(456)
    result2 <- gextract(vtrack_name, intervals, iterator = intervals)[[vtrack_name]]

    expect_equal(result1, result2)
})

# Tests for 'sample.pos.abs' - statistical properties
test_that("sample.pos.abs returns valid position from interval", {
    track_name <- "test.sparse"
    vtrack_name <- "vt_sample_pos_abs"
    intervals <- gintervals(1, 0, 1000)

    gvtrack.create(vtrack_name, track_name, "sample.pos.abs")
    on.exit(gvtrack.rm(vtrack_name), add = TRUE)

    # Get all positions in the interval
    all_vals <- gextract(track_name, intervals, iterator = 1)
    valid_positions <- all_vals$start[!is.na(all_vals[[track_name]])]

    # Test that sampled positions are valid
    set.seed(789)
    positions <- sapply(1:10, function(i) {
        gextract(vtrack_name, intervals, iterator = intervals)[[vtrack_name]]
    })

    # All positions should be in the valid set
    expect_true(all(positions %in% valid_positions))
    # At least some positions should be non-NA
    expect_true(any(!is.na(positions)))
})

# Tests for 'sample.pos.relative'
test_that("sample.pos.relative returns relative position", {
    track_name <- "test.sparse"
    vtrack_name <- "vt_sample_pos_rel"
    intervals <- gintervals(1, 500, 1500)

    gvtrack.create(vtrack_name, track_name, "sample.pos.relative")
    on.exit(gvtrack.rm(vtrack_name), add = TRUE)

    set.seed(101)
    result <- gextract(vtrack_name, intervals, iterator = intervals)[[vtrack_name]]

    # Should be relative to start (>= 0)
    expect_gte(result, 0)
    # Should be within interval length
    expect_lte(result, 1000)
})

test_that("sample.pos.relative honors iterator shifts", {
    track_name <- "test.fixedbin"
    vtrack_name <- "vt_sample_pos_rel_shift"
    intervals <- gintervals(1, 200, 400)
    sshift <- -50L
    eshift <- 100L

    gvtrack.create(vtrack_name, track_name, "sample.pos.relative")
    gvtrack.iterator(vtrack_name, sshift = sshift, eshift = eshift)
    on.exit(gvtrack.rm(vtrack_name), add = TRUE)

    set.seed(102)
    result <- gextract(vtrack_name, intervals, iterator = intervals)[[vtrack_name]]

    # Should be relative to original interval start
    expect_gte(result, -50) # sshift
    expect_lte(result, 300) # interval length + eshift
})

# Edge case tests
test_that("all functions handle single-value intervals correctly", {
    track_name <- "test.fixedbin"
    intervals <- gintervals(1, 0, 10) # Small interval

    funcs <- c("exists", "size", "first", "last", "first.pos.abs", "last.pos.abs")

    for (func in funcs) {
        vtrack_name <- paste0("vt_single_", gsub("\\.", "_", func))
        gvtrack.create(vtrack_name, track_name, func)
        result <- gextract(vtrack_name, intervals, iterator = intervals)
        expect_true(nrow(result) == 1)
        gvtrack.rm(vtrack_name)
    }
})

test_that("position functions work with dense tracks", {
    track_name <- "test.fixedbin"
    vtrack_names <- c("first.pos.abs", "first.pos.relative", "last.pos.abs", "last.pos.relative")
    intervals <- rbind(
        gintervals(1, 0, 200),
        gintervals(1, 500, 700)
    )

    for (func in vtrack_names) {
        vtrack_name <- paste0("vt_dense_", gsub("\\.", "_", func))
        gvtrack.create(vtrack_name, track_name, func)

        result <- gextract(vtrack_name, intervals, iterator = intervals)
        # Should have results for both intervals
        expect_equal(nrow(result), 2)
        # Values should not all be NA
        expect_true(any(!is.na(result[[vtrack_name]])))

        gvtrack.rm(vtrack_name)
    }
})
