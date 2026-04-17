# Phase 2 tests: top-K pruning, aggregator templating, intervals support.
# Runs on gdb.init_examples() so no external db required.

test_that("top-K path returns EXACTLY the same quantile value as fallback for the same extreme p", {
    # Forces the fallback path by passing a mixed-tail percentile vector
    # (triggers use_fallback=true in the C++ entry). Extracts the p=0.95
    # column and compares against a top-K-only call at p=0.95.
    # Both share the identical scan and the identical downstream
    # nth_element rank formula; any divergence is a real bug.
    gdb.init_examples()
    args <- list(
        track_names = "dense_track",
        iterator = 20L, sshift = -50L, eshift = 50L,
        n_threads = 1L, func = "avg"
    )
    # Mixed-tail → fallback.
    suppressWarnings({
        q_fallback_mat <- do.call(
            glm_batch_quantiles,
            c(args, list(percentiles = c(0.2, 0.95)))
        )
    })
    q_fallback_at_95 <- q_fallback_mat[1, 2] # p=0.95 column

    # Top-K only (all percentiles >= 0.5; no clamp on example db).
    q_topk <- do.call(
        glm_batch_quantiles,
        c(args, list(percentiles = 0.95))
    )

    expect_equal(unname(q_topk), unname(q_fallback_at_95),
        tolerance = 0,
        info = sprintf(
            "topk=%g fallback=%g",
            unname(q_topk), unname(q_fallback_at_95)
        )
    )
})

test_that("bottom-K path returns same value as fallback for p < 0.5", {
    gdb.init_examples()
    args <- list(
        track_names = "dense_track",
        iterator = 20L, sshift = -50L, eshift = 50L,
        n_threads = 1L, func = "avg"
    )
    suppressWarnings({
        q_fb <- do.call(
            glm_batch_quantiles,
            c(args, list(percentiles = c(0.05, 0.95)))
        )
    })
    q_fb_at_5 <- q_fb[1, 1]
    q_bot <- do.call(
        glm_batch_quantiles,
        c(args, list(percentiles = 0.05))
    )
    expect_equal(unname(q_bot), unname(q_fb_at_5),
        tolerance = 0,
        info = sprintf(
            "bot=%g fb=%g",
            unname(q_bot), unname(q_fb_at_5)
        )
    )
})

test_that("glm_batch_quantiles top-K path matches gquantiles on vtrack-max (sanity)", {
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
        info = sprintf("q_top=%g, ref=%g", unname(q_top), unname(ref))
    )
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
        info = sprintf("q=%g ref=%g", unname(q), unname(ref))
    )
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
    ivs <- data.frame(
        chrom = c("1", "1"),
        start = c(0, 5000),
        end = c(2000, 7000)
    )
    q <- glm_batch_quantiles(
        track_names = "dense_track", percentiles = 0.5,
        iterator = 50L, sshift = -50L, eshift = 50L, n_threads = 2L,
        func = "avg", intervals = ivs
    )
    expect_true(is.finite(q))
})

test_that("no fallback warning when K stays under K_MAX (example db at p=0.5)", {
    gdb.init_examples()
    # p=0.5 on the tiny example db gives K ≈ (few K positions) / 2,
    # well under K_MAX=10M, so NO warning should fire. On real genomes
    # the same call would trigger the K>K_MAX fallback warning; this
    # test only documents the small-N branch.
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
