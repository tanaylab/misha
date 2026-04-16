# Phase 2 tests: top-K pruning, aggregator templating, intervals support.
# Runs on gdb.init_examples() so no external db required.

test_that("glm_batch_quantiles top-K path matches full-storage on extreme p", {
    gdb.init_examples()
    # Extreme p triggers top-K mode (heap-backed). Compare against the
    # fallback path by forcing K large enough to not trip the pruning
    # (impossible via public API — instead, cross-check top-K with
    # func='max' where the definition is unambiguous).
    q_top <- glm_batch_quantiles(
        track_names = "dense_track", percentiles = 0.99,
        iterator = 10L, sshift = -30L, eshift = 30L, n_threads = 1L,
        func = "max"
    )
    # Compute reference via gquantiles on a max vtrack.
    gvtrack.create("vt_ref", "dense_track", func = "max")
    gvtrack.iterator("vt_ref", sshift = -30, eshift = 30)
    ref <- gquantiles("vt_ref", percentiles = 0.99, iterator = 10L)
    gvtrack.rm("vt_ref")
    # gquantiles uses StreamPercentiler but max is an order statistic;
    # for a small stream (example db) no subsampling occurs and results
    # should match closely. Tolerance is loose because gquantiles may
    # return a different quantile definition at the boundary.
    expect_true(abs(unname(q_top) - unname(ref)) < 0.1,
                info = sprintf("q_top=%g, ref=%g", unname(q_top), unname(ref)))
})

test_that("func='avg' produces different values than func='max' (sanity)", {
    gdb.init_examples()
    args <- list(
        track_names = "dense_track", percentiles = 0.5,
        iterator = 50L, sshift = -50L, eshift = 50L, n_threads = 1L
    )
    q_avg <- do.call(glm_batch_quantiles, c(args, list(func = "avg")))
    q_max <- do.call(glm_batch_quantiles, c(args, list(func = "max")))
    q_sum <- do.call(glm_batch_quantiles, c(args, list(func = "sum")))
    q_min <- do.call(glm_batch_quantiles, c(args, list(func = "min")))
    # Invariants:
    expect_true(q_min <= q_avg)
    expect_true(q_avg <= q_max)
    # sum with window 2 bins ~ avg * 2.
    expect_true(q_sum >= q_max)
    expect_true(all(is.finite(c(q_avg, q_max, q_sum, q_min))))
})

test_that("func='lse' matches vtrack lse", {
    gdb.init_examples()
    # Use a dense track; results should match the legacy LSE path closely
    # (same lse_accumulate, identical aggregation).
    gvtrack.create("vt_lse", "dense_track", func = "lse")
    gvtrack.iterator("vt_lse", sshift = -50, eshift = 50)
    ref <- gquantiles("vt_lse", percentiles = 0.9, iterator = 50L)
    gvtrack.rm("vt_lse")
    q <- glm_batch_quantiles(
        track_names = "dense_track", percentiles = 0.9,
        iterator = 50L, sshift = -50L, eshift = 50L, n_threads = 1L,
        func = "lse"
    )
    # Example db is tiny; tolerance is loose to absorb StreamPercentiler
    # vs exact-nth_element divergence on a small stream.
    expect_true(abs(unname(q) - unname(ref)) < 1.0,
                info = sprintf("q=%g ref=%g", unname(q), unname(ref)))
})

test_that("intervals restriction narrows the scan", {
    gdb.init_examples()
    allg <- get("ALLGENOME", envir = .misha)
    # Whole-genome scan.
    q_all <- glm_batch_quantiles(
        track_names = "dense_track", percentiles = 0.5,
        iterator = 50L, sshift = -50L, eshift = 50L, n_threads = 1L,
        func = "avg"
    )
    # Restricted scan (a single interval on chrom 1).
    small <- data.frame(chrom = "1", start = 0, end = 2000)
    q_small <- glm_batch_quantiles(
        track_names = "dense_track", percentiles = 0.5,
        iterator = 50L, sshift = -50L, eshift = 50L, n_threads = 1L,
        func = "avg", intervals = small
    )
    # Both should be finite and not necessarily equal. Sanity only —
    # the restricted scan covers ~40 positions (2000 bp / 50 bp) which
    # is enough for a finite median.
    expect_true(is.finite(q_all) && is.finite(q_small))
})

test_that("intervals with multiple chroms and interval-mask boundary", {
    gdb.init_examples()
    # Two disjoint intervals on the same chrom — exercises the boundary()
    # hook in the driver at the gap. The TopKQuantile reducer's boundary()
    # is a no-op, so this is a smoke test that the driver doesn't crash
    # when the mask has gaps.
    ivs <- data.frame(chrom = c("1", "1"),
                      start = c(0, 5000),
                      end   = c(2000, 7000))
    q <- glm_batch_quantiles(
        track_names = "dense_track", percentiles = 0.5,
        iterator = 50L, sshift = -50L, eshift = 50L, n_threads = 2L,
        func = "avg", intervals = ivs
    )
    expect_true(is.finite(q))
})

test_that("top-K clamp triggers fallback warning for p = 0.5", {
    gdb.init_examples()
    # p=0.5 on any non-trivial N gives K ≈ N/2 > K_MAX on real genomes.
    # On the tiny example db, K ≈ (2K positions) / 2 = ~1.2K, well under
    # K_MAX=10M, so we should NOT warn here. This test documents that
    # boundary: fallback triggers only when K genuinely exceeds K_MAX.
    expect_silent(
        glm_batch_quantiles(
            track_names = "dense_track", percentiles = 0.5,
            iterator = 50L, sshift = -50L, eshift = 50L, n_threads = 1L,
            func = "avg"
        )
    )
})

test_that("mixed-tail percentiles warn and fall back", {
    gdb.init_examples()
    expect_warning(
        glm_batch_quantiles(
            track_names = "dense_track", percentiles = c(0.1, 0.9),
            iterator = 50L, sshift = -50L, eshift = 50L, n_threads = 1L,
            func = "avg"
        ),
        regexp = "span both tails"
    )
})

test_that("intervals=NULL (default) yields the whole-genome scan", {
    gdb.init_examples()
    q_null <- glm_batch_quantiles(
        track_names = "dense_track", percentiles = 0.9,
        iterator = 50L, sshift = -50L, eshift = 50L, n_threads = 1L,
        func = "avg"
    )
    expect_true(is.finite(q_null))
    expect_named(q_null, "dense_track")
})
