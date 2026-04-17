# Tests for gsummary fast path (Phase 3) on gdb.init_examples().

test_that("gsummary single-expr legacy path is unchanged when fast=FALSE", {
    gdb.init_examples()
    r <- gsummary("dense_track", fast = FALSE)
    expect_type(r, "double")
    expect_named(r)
    expect_true(all(c(
        "Total intervals", "NaN intervals",
        "Min", "Max", "Sum", "Mean", "Std dev"
    ) %in% names(r)))
})

test_that("gsummary single-expr with fast=TRUE still returns legacy-shaped row", {
    gdb.init_examples()
    # Single-expression calls always take the legacy path (even with
    # fast=TRUE), preserving bit-exact back-compat with pre-refactor
    # output. Multi-track is where the fast path activates.
    r <- gsummary("dense_track", fast = TRUE)
    expect_type(r, "double")
    expect_true(all(c(
        "Total intervals", "NaN intervals",
        "Min", "Max", "Sum", "Mean", "Std dev"
    ) %in% names(r)))
})

test_that("gsummary with vector of tracks returns a data.frame", {
    gdb.init_examples()
    tracks <- c("dense_track", "subdir.dense_track2")
    df <- gsummary(tracks, iterator = 20L, fast = TRUE)
    expect_s3_class(df, "data.frame")
    expect_equal(
        colnames(df),
        c("track", "n", "n_nan", "min", "max", "sum", "mean", "sd")
    )
    expect_equal(nrow(df), 2)
    expect_equal(df$track, tracks)
    expect_true(all(is.finite(df$mean)))
    expect_true(all(df$min <= df$mean))
    expect_true(all(df$mean <= df$max))
})

test_that("gsummary multi-track fast path matches per-track single calls", {
    gdb.init_examples()
    tracks <- c("dense_track", "subdir.dense_track2")
    # Omit iterator — fast path uses the tracks' native bin size (50bp
    # here), matching the slow path's implicit per-bin iteration.
    fast <- gsummary(tracks, fast = TRUE)
    # Compute reference per-track using the legacy slow path.
    ref <- do.call(rbind, lapply(tracks, function(t) {
        r <- gsummary(t, fast = FALSE)
        data.frame(
            track = t,
            n = r["Total intervals"],
            n_nan = r["NaN intervals"],
            min = r["Min"], max = r["Max"], sum = r["Sum"],
            mean = r["Mean"], sd = r["Std dev"],
            row.names = NULL, stringsAsFactors = FALSE
        )
    }))
    # Both paths iterate bin-by-bin for these dense tracks.
    # n = total scan positions (NaN + non-NaN); n_nan = NaN-aggregate
    # positions. Must match legacy exactly.
    expect_equal(fast$n, ref$n, tolerance = 0)
    expect_equal(fast$n_nan, ref$n_nan, tolerance = 0)
    expect_equal(fast$min, ref$min, tolerance = 1e-5)
    expect_equal(fast$max, ref$max, tolerance = 1e-5)
    expect_equal(fast$mean, ref$mean, tolerance = 1e-5)
    expect_equal(fast$sd, ref$sd, tolerance = 1e-5)
})

test_that("gsummary fast path accepts intervals restriction", {
    gdb.init_examples()
    allg <- get("ALLGENOME", envir = .misha)
    small <- data.frame(
        chrom = "1", start = 0, end = 2000,
        stringsAsFactors = FALSE
    )
    # Intervals as a data.frame must work.
    df <- gsummary(c("dense_track", "subdir.dense_track2"),
        iterator = 20L, intervals = small, fast = TRUE
    )
    expect_s3_class(df, "data.frame")
    expect_equal(nrow(df), 2)
    expect_true(all(is.finite(df$mean)))
})

test_that("gsummary falls through to slow path on complex expression", {
    gdb.init_examples()
    # A complex expression (arithmetic) — not fast-path eligible.
    # Suppress the one-time dispatch message.
    options(misha.quiet_dispatch = TRUE)
    r <- gsummary("dense_track + 1", fast = TRUE)
    expect_type(r, "double")
    expect_named(r)
    options(misha.quiet_dispatch = NULL)
})

test_that("gsummary multi-expr with fast=FALSE errors cleanly (Phase 5 scope)", {
    gdb.init_examples()
    expect_error(
        gsummary(c("dense_track", "subdir.dense_track2"), fast = FALSE),
        regexp = "multi-expression slow path not yet implemented"
    )
})

test_that("gsummary with vtrack (LSE window) — multi-expr fast path", {
    gdb.init_examples()
    gvtrack.create(vtrack = "vt_smry1", src = "dense_track", func = "lse")
    gvtrack.create(vtrack = "vt_smry2", src = "subdir.dense_track2", func = "lse")
    gvtrack.iterator("vt_smry1", sshift = -50, eshift = 50)
    gvtrack.iterator("vt_smry2", sshift = -50, eshift = 50)
    df <- gsummary(c("vt_smry1", "vt_smry2"), iterator = 20L, fast = TRUE)
    gvtrack.rm("vt_smry1")
    gvtrack.rm("vt_smry2")
    expect_s3_class(df, "data.frame")
    expect_equal(nrow(df), 2)
    expect_true(all(is.finite(df$mean)))
})
