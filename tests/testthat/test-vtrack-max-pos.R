create_isolated_test_db()

manual_max_pos <- function(track_name, interval_row, relative = FALSE, relative_ref = NULL) {
    interval <- gintervals(interval_row$chrom, interval_row$start, interval_row$end)
    values <- gextract(track_name, interval, iterator = 1)
    track_vals <- values[[track_name]]
    valid <- !is.na(track_vals)
    if (!any(valid)) {
        return(NA_real_)
    }
    values <- values[valid, , drop = FALSE]
    track_vals <- values[[track_name]]
    max_val <- max(track_vals)
    best_pos <- values$start[which(track_vals == max_val)][1]
    if (relative) {
        ref_start <- if (!is.null(relative_ref)) relative_ref else interval_row$start
        best_pos - ref_start
    } else {
        best_pos
    }
}

manual_max_pos_vec <- function(track_name, intervals, relative = FALSE, relative_ref = NULL) {
    vapply(seq_len(nrow(intervals)), function(i) {
        ref <- if (is.null(relative_ref)) NULL else relative_ref[i]
        manual_max_pos(track_name, intervals[i, ], relative = relative, relative_ref = ref)
    }, numeric(1))
}

manual_min_pos <- function(track_name, interval_row, relative = FALSE, relative_ref = NULL) {
    interval <- gintervals(interval_row$chrom, interval_row$start, interval_row$end)
    values <- gextract(track_name, interval, iterator = 1)
    track_vals <- values[[track_name]]
    valid <- !is.na(track_vals)
    if (!any(valid)) {
        return(NA_real_)
    }
    values <- values[valid, , drop = FALSE]
    track_vals <- values[[track_name]]
    min_val <- min(track_vals)
    best_pos <- values$start[which(track_vals == min_val)][1]
    if (relative) {
        ref_start <- if (!is.null(relative_ref)) relative_ref else interval_row$start
        best_pos - ref_start
    } else {
        best_pos
    }
}

manual_min_pos_vec <- function(track_name, intervals, relative = FALSE, relative_ref = NULL) {
    vapply(seq_len(nrow(intervals)), function(i) {
        ref <- if (is.null(relative_ref)) NULL else relative_ref[i]
        manual_min_pos(track_name, intervals[i, ], relative = relative, relative_ref = ref)
    }, numeric(1))
}

apply_interval_shift <- function(intervals, sshift = 0L, eshift = 0L) {
    shifted <- intervals
    shifted$start <- shifted$start + sshift
    shifted$end <- shifted$end + eshift
    shifted
}

test_that("max.pos.abs matches manual argmax for dense tracks", {
    track_name <- "test.fixedbin"
    vtrack_name <- "vt_maxpos_abs_dense"
    intervals <- rbind(
        gintervals(1, 0, 200),
        gintervals(1, 500, 900),
        gintervals(1, 1500, 2100)
    )

    gvtrack.create(vtrack_name, track_name, "max.pos.abs")
    on.exit(gvtrack.rm(vtrack_name), add = TRUE)

    vt_res <- gextract(vtrack_name, intervals, iterator = intervals)
    manual <- manual_max_pos_vec(track_name, intervals, relative = FALSE)
    expect_equal(vt_res[[vtrack_name]], manual)
})

test_that("max.pos.abs honors iterator shifts", {
    track_name <- "test.fixedbin"
    vtrack_name <- "vt_maxpos_abs_dense_shift"
    intervals <- rbind(
        gintervals(1, 100, 250),
        gintervals(1, 700, 900)
    )
    sshift <- -50L
    eshift <- 75L

    gvtrack.create(vtrack_name, track_name, "max.pos.abs")
    gvtrack.iterator(vtrack_name, sshift = sshift, eshift = eshift)
    on.exit(gvtrack.rm(vtrack_name), add = TRUE)

    vt_res <- gextract(vtrack_name, intervals, iterator = intervals)
    shifted <- apply_interval_shift(intervals, sshift = sshift, eshift = eshift)
    manual <- manual_max_pos_vec(track_name, shifted, relative = FALSE)
    expect_equal(vt_res[[vtrack_name]], manual)
})

test_that("max.pos.relative matches manual argmax for sparse tracks", {
    track_name <- "test.sparse"
    vtrack_name <- "vt_maxpos_rel_sparse"
    intervals <- rbind(
        gintervals(1, 0, 300),
        gintervals(1, 600, 1000),
        gintervals(1, 1200, 1500)
    )

    gvtrack.create(vtrack_name, track_name, "max.pos.relative")
    on.exit(gvtrack.rm(vtrack_name), add = TRUE)

    vt_res <- gextract(vtrack_name, intervals, iterator = intervals)
    manual <- manual_max_pos_vec(track_name, intervals, relative = TRUE)
    expect_equal(vt_res[[vtrack_name]], manual)
})

test_that("max.pos.relative honors iterator shifts on dense tracks", {
    track_name <- "test.fixedbin"
    vtrack_name <- "vt_maxpos_rel_dense_shift"
    intervals <- rbind(
        gintervals(1, 200, 320),
        gintervals(1, 600, 760)
    )
    sshift <- -30L
    eshift <- 60L

    gvtrack.create(vtrack_name, track_name, "max.pos.relative")
    gvtrack.iterator(vtrack_name, sshift = sshift, eshift = eshift)
    on.exit(gvtrack.rm(vtrack_name), add = TRUE)

    vt_res <- gextract(vtrack_name, intervals, iterator = intervals)
    shifted <- apply_interval_shift(intervals, sshift = sshift, eshift = eshift)
    manual <- manual_max_pos_vec(track_name, shifted, relative = TRUE)
    expect_equal(vt_res[[vtrack_name]], manual)
})

test_that("min.pos.abs matches manual argmin for dense tracks", {
    track_name <- "test.fixedbin"
    vtrack_name <- "vt_minpos_abs_dense"
    intervals <- rbind(
        gintervals(1, 0, 200),
        gintervals(1, 500, 900),
        gintervals(1, 1500, 2100)
    )

    gvtrack.create(vtrack_name, track_name, "min.pos.abs")
    on.exit(gvtrack.rm(vtrack_name), add = TRUE)

    vt_res <- gextract(vtrack_name, intervals, iterator = intervals)
    manual <- manual_min_pos_vec(track_name, intervals, relative = FALSE)
    expect_equal(vt_res[[vtrack_name]], manual)
})

test_that("min.pos.abs honors iterator shifts", {
    track_name <- "test.fixedbin"
    vtrack_name <- "vt_minpos_abs_dense_shift"
    intervals <- rbind(
        gintervals(1, 100, 250),
        gintervals(1, 700, 900)
    )
    sshift <- -50L
    eshift <- 75L

    gvtrack.create(vtrack_name, track_name, "min.pos.abs")
    gvtrack.iterator(vtrack_name, sshift = sshift, eshift = eshift)
    on.exit(gvtrack.rm(vtrack_name), add = TRUE)

    vt_res <- gextract(vtrack_name, intervals, iterator = intervals)
    shifted <- apply_interval_shift(intervals, sshift = sshift, eshift = eshift)
    manual <- manual_min_pos_vec(track_name, shifted, relative = FALSE)
    expect_equal(vt_res[[vtrack_name]], manual)
})

test_that("min.pos.relative matches manual argmin for sparse tracks", {
    track_name <- "test.sparse"
    vtrack_name <- "vt_minpos_rel_sparse"
    intervals <- rbind(
        gintervals(1, 0, 300),
        gintervals(1, 600, 1000),
        gintervals(1, 1200, 1500)
    )

    gvtrack.create(vtrack_name, track_name, "min.pos.relative")
    on.exit(gvtrack.rm(vtrack_name), add = TRUE)

    vt_res <- gextract(vtrack_name, intervals, iterator = intervals)
    manual <- manual_min_pos_vec(track_name, intervals, relative = TRUE)
    expect_equal(vt_res[[vtrack_name]], manual)
})

test_that("min.pos.relative honors iterator shifts on dense tracks", {
    track_name <- "test.fixedbin"
    vtrack_name <- "vt_minpos_rel_dense_shift"
    intervals <- rbind(
        gintervals(1, 200, 320),
        gintervals(1, 600, 760)
    )
    sshift <- -30L
    eshift <- 60L

    gvtrack.create(vtrack_name, track_name, "min.pos.relative")
    gvtrack.iterator(vtrack_name, sshift = sshift, eshift = eshift)
    on.exit(gvtrack.rm(vtrack_name), add = TRUE)

    vt_res <- gextract(vtrack_name, intervals, iterator = intervals)
    shifted <- apply_interval_shift(intervals, sshift = sshift, eshift = eshift)
    manual <- manual_min_pos_vec(track_name, shifted, relative = TRUE)
    expect_equal(vt_res[[vtrack_name]], manual)
})
