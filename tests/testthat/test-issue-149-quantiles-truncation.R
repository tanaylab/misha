# Regression test for issue #149: non-multitask gintervals.quantiles() truncated a
# large scope to the first evaluation buffer (gbuf.size == 1000 bins). Root cause:
# estimate_num_bins() was called AFTER scanner.begin() in the non-multitask
# gintervals_quantiles; it iterates the scope object the scanner's iterator shares,
# leaving the scope cursor at end() so the scan stopped after the first batch.
# The bug only surfaces when the scope spans more than one evaluation buffer
# (> gbuf.size), so these tests deliberately use scopes well past 1000 intervals.

create_isolated_test_db()

# canonical order so the streaming (per-chrom) layout compares to the in-memory one
ord2d <- function(d) d[order(d$chrom1, d$start1, d$chrom2, d$start2), ]
ord1d <- function(d) d[order(d$chrom, d$start, d$end), ]

test_that("issue #149: non-multitask 2D quantiles do not truncate at the eval-buffer boundary", {
    full <- gscreen("test.rects > 40", gintervals.2d(c(1, 2), 0, -1))
    expect_gt(nrow(full), 1000) # must exceed gbuf.size to exercise the bug
    scope <- full[seq_len(5000), ] # >> 1000: spans several evaluation buffers

    # multitask in-memory is the known-good reference (it never shared the cursor)
    ref <- ord2d(gintervals.quantiles("test.rects", c(0.1, 0.5, 0.9), scope, iterator = c(1, 1)))
    expect_true(all(is.finite(ref[[7]]))) # sanity: reference is fully populated

    withr::local_options(list(gmultitasking = FALSE, gbig.intervals.size = 10))

    # non-multitask in-memory
    im <- ord2d(gintervals.quantiles("test.rects", c(0.1, 0.5, 0.9), scope, iterator = c(1, 1)))
    expect_equal(sum(is.finite(im[[7]])), nrow(scope))
    expect_equal(im, ref, ignore_attr = TRUE)

    # non-multitask streaming big-set (the path in the bug report)
    tmp <- paste0("test.q149_", sample(1:1e9, 1))
    withr::defer(gintervals.rm(tmp, force = TRUE))
    gintervals.quantiles("test.rects", c(0.1, 0.5, 0.9), scope, iterator = c(1, 1), intervals.set.out = tmp)
    big <- ord2d(gintervals.load(tmp))
    expect_equal(sum(is.finite(big[[7]])), nrow(scope))
    expect_equal(big, ref, ignore_attr = TRUE)
})

test_that("issue #149: non-multitask 1D quantiles do not truncate at the eval-buffer boundary", {
    scope <- gscreen("test.fixedbin > 0.2", gintervals(1, 0, -1))
    scope <- scope[seq_len(5000), ] # >> gbuf.size
    expect_gt(nrow(scope), 1000)

    ref <- ord1d(gintervals.quantiles("test.fixedbin", c(0.25, 0.5, 0.75), scope, iterator = 50))
    expect_true(all(is.finite(ref[[ncol(ref)]])))

    withr::local_options(list(gmultitasking = FALSE, gbig.intervals.size = 10))

    im <- ord1d(gintervals.quantiles("test.fixedbin", c(0.25, 0.5, 0.75), scope, iterator = 50))
    expect_equal(im, ref, ignore_attr = TRUE)

    tmp <- paste0("test.q149b_", sample(1:1e9, 1))
    withr::defer(gintervals.rm(tmp, force = TRUE))
    gintervals.quantiles("test.fixedbin", c(0.25, 0.5, 0.75), scope, iterator = 50, intervals.set.out = tmp)
    big <- ord1d(gintervals.load(tmp))
    expect_equal(big, ref, ignore_attr = TRUE)
})
