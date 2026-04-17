# Phase 6 dispatch tests for gquantiles.

test_that("gquantiles single-expr default (fast=FALSE) unchanged", {
    gdb.init_examples()
    q <- gquantiles("dense_track", c(0.1, 0.5, 0.9), gintervals(c(1, 2)))
    expect_type(q, "double")
    expect_named(q)
    expect_equal(length(q), 3)
})

test_that("gquantiles single-expr with fast=TRUE returns named numeric", {
    gdb.init_examples()
    q <- gquantiles("dense_track", c(0.9, 0.99), fast = TRUE)
    expect_type(q, "double")
    expect_named(q)
    expect_equal(length(q), 2)
    expect_true(all(is.finite(q)))
    # Monotonicity: q[0.9] <= q[0.99]
    expect_true(q[1] <= q[2])
})

test_that("gquantiles multi-expr fast=TRUE returns data.frame", {
    gdb.init_examples()
    tracks <- c("dense_track", "subdir.dense_track2")
    df <- gquantiles(tracks, c(0.5, 0.9), fast = TRUE)
    expect_s3_class(df, "data.frame")
    expect_equal(colnames(df), c("track", "0.5", "0.9"))
    expect_equal(nrow(df), 2)
    expect_equal(df$track, tracks)
    expect_true(all(is.finite(df[["0.5"]])))
    expect_true(all(df[["0.5"]] <= df[["0.9"]]))
})

test_that("gquantiles multi-expr without fast=TRUE errors cleanly", {
    gdb.init_examples()
    expect_error(
        gquantiles(c("dense_track", "subdir.dense_track2"), 0.9),
        regexp = "multi-expression calls require fast=TRUE"
    )
})

test_that("gquantiles multi-expr with complex expression errors cleanly", {
    gdb.init_examples()
    # Complex expressions aren't fast-path eligible.
    expect_error(
        gquantiles(c("dense_track + 1", "subdir.dense_track2"), 0.9,
            fast = TRUE
        ),
        regexp = "fast-path not eligible"
    )
})

test_that("gquantiles single-expr fast=TRUE matches legacy on small examples", {
    gdb.init_examples()
    # On the tiny example db, StreamPercentiler doesn't sub-sample, so the
    # numbers should be close (not necessarily bit-equal because the window
    # vs bin-by-bin semantics differ for non-matching iterator/bin_size).
    # Here we use the track's natural bin size (50) as the iterator to
    # align with the legacy path's implicit per-bin iteration.
    q_fast <- gquantiles("dense_track", c(0.1, 0.5, 0.9),
        iterator = 50L, fast = TRUE
    )
    q_slow <- gquantiles("dense_track", c(0.1, 0.5, 0.9), iterator = 50L)
    expect_equal(unname(q_fast), unname(q_slow), tolerance = 0.01)
})

test_that("gquantiles vtrack vector works for fast=TRUE", {
    gdb.init_examples()
    gvtrack.create(vtrack = "vt_q1", src = "dense_track", func = "lse")
    gvtrack.create(
        vtrack = "vt_q2", src = "subdir.dense_track2",
        func = "lse"
    )
    gvtrack.iterator("vt_q1", sshift = -50, eshift = 50)
    gvtrack.iterator("vt_q2", sshift = -50, eshift = 50)
    df <- gquantiles(c("vt_q1", "vt_q2"), 0.9,
        iterator = 20L, fast = TRUE
    )
    gvtrack.rm("vt_q1")
    gvtrack.rm("vt_q2")
    expect_s3_class(df, "data.frame")
    expect_equal(nrow(df), 2)
    expect_true(all(is.finite(df[["0.9"]])))
})
