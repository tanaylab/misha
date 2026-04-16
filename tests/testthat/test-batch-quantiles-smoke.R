# Synthetic smoke tests for the BatchTrackScan code path via glm_batch_quantiles.
# These run against gdb.init_examples() and don't require the mm10 motif db.
# They exercise the templated scan driver, per-(track, chrom) task queue,
# and TopKQuantile reducer in fallback mode (Phase 1) for both dense and
# sparse tracks.

test_that("glm_batch_quantiles works on a single dense example track", {
    gdb.init_examples()
    q <- glm_batch_quantiles(
        track_names = "dense_track",
        percentiles = 0.5,
        iterator = 50L,
        sshift = -50L,
        eshift = 50L,
        n_threads = 1L
    )
    expect_type(q, "double")
    expect_equal(length(q), 1)
    expect_true(is.finite(q))
    expect_named(q, "dense_track")
})

test_that("glm_batch_quantiles handles multiple percentiles (returns matrix)", {
    gdb.init_examples()
    q <- glm_batch_quantiles(
        track_names = "dense_track",
        percentiles = c(0.1, 0.5, 0.9),
        iterator = 50L,
        sshift = -50L,
        eshift = 50L,
        n_threads = 1L
    )
    expect_true(is.matrix(q))
    expect_equal(dim(q), c(1, 3))
    expect_equal(rownames(q), "dense_track")
    expect_true(all(is.finite(q)))
    # Quantiles must be monotone nondecreasing.
    expect_true(q[1, 1] <= q[1, 2])
    expect_true(q[1, 2] <= q[1, 3])
})

test_that("glm_batch_quantiles is deterministic across thread counts", {
    gdb.init_examples()
    tracks <- c("dense_track", "subdir.dense_track2")
    args <- list(
        track_names = tracks, percentiles = 0.75,
        iterator = 50L, sshift = -50L, eshift = 50L
    )
    q1 <- do.call(glm_batch_quantiles, c(args, list(n_threads = 1L)))
    q2 <- do.call(glm_batch_quantiles, c(args, list(n_threads = 2L)))
    expect_equal(q1, q2, tolerance = 0)
})

test_that("glm_batch_quantiles vs gquantiles(vtrack, LSE) — same ballpark", {
    # Note: on the tiny example db (~2K positions), StreamPercentiler vs
    # exact nth_element can differ substantially more than on full genomes
    # (the approximate tail estimator is noisy at small N). Use a loose
    # tolerance: this test confirms the batched path produces a value in
    # the same order of magnitude, not that it matches the legacy
    # approximate path. Tight correctness is covered by the motif-db test
    # (skipped without MISHA_MM10_MOTIF_DB) and by Phase 2 synthetic
    # fast-vs-slow parity tests.
    gdb.init_examples()
    track <- "dense_track"

    gvtrack.create("vt_ref", track, func = "lse")
    gvtrack.iterator("vt_ref", sshift = -100, eshift = 100)
    ref <- gquantiles("vt_ref", percentiles = 0.5, iterator = 50L)
    gvtrack.rm("vt_ref")

    fast <- glm_batch_quantiles(
        track_names = track, percentiles = 0.5,
        iterator = 50L, sshift = -100L, eshift = 100L, n_threads = 1L
    )

    expect_true(is.finite(ref) && is.finite(fast))
    expect_true(abs(unname(fast) - unname(ref)) < 1.0)
})

test_that("glm_batch_quantiles processes a sparse track", {
    gdb.init_examples()
    # sparse_track: quantile computation should succeed and return finite.
    q <- glm_batch_quantiles(
        track_names = "sparse_track",
        percentiles = 0.5,
        iterator = 100L,
        sshift = -50L,
        eshift = 50L,
        n_threads = 1L
    )
    expect_type(q, "double")
    expect_true(is.finite(q))
})
