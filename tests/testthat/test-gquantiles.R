create_isolated_test_db()

test_that("gquantiles with test.fixedbin", {
    intervs <- gscreen("test.fixedbin > 0.2", gintervals(c(1, 2), 0, -1))
    result <- gquantiles("test.fixedbin+0.2", percentile = c(0.5, 0.3, 0.2, 0.9), intervs)
    expect_regression(result, "gquantiles_fixedbin_result")
})

test_that("gquantiles with test.rects", {
    result <- gquantiles("test.rects", percentile = c(0.5, 0.3, 0.2, 0.9, 0.999), gintervals.2d(chroms1 = c(2, 3), chroms2 = c(2, 4)))
    expect_regression(result, "gquantiles_rects_result")
})

test_that("gquantiles with test.computed2d", {
    result <- gquantiles("test.computed2d", percentile = c(0.5, 0.3, 0.2, 0.9, 0.999), gintervals.2d(chroms1 = c(6, 5), chroms2 = c(8, 9)))
    expect_regression(result, "gquantiles_computed2d_result")
})

test_that("gquantiles with test.fixedbin without intervals", {
    intervs <- gscreen("test.fixedbin > 0.2", gintervals(c(1, 2), 0, -1))
    result <- gquantiles("test.fixedbin+0.2", percentile = c(0.5, 0.999))
    expect_regression(result, "gquantiles_fixedbin_no_intervals_result")
})

# Regression: with many percentiles (e.g. .gtrack.prepare.pvals requests
# ~12k), the single-process nth_element-on-suffix path was O(k * N), which
# could take hours on a dense full-genome track. The fix flips to a single
# O(N log N) sort once `targets.size()` exceeds ~2*log2(N), and must
# remain bit-identical to the multitask k-way merge result.
test_that("gquantiles single-process matches multitask for many percentiles", {
    pcts <- sort(unique(c(
        seq(0, 1, length.out = 200),
        seq(0.001, 0.999, length.out = 800)
    )))

    multitask <- .ggetOption("gmultitasking")
    on.exit(options(gmultitasking = multitask), add = TRUE)

    options(gmultitasking = FALSE)
    serial <- gquantiles("test.fixedbin",
        percentile = pcts,
        intervals = gintervals(c(1, 2), 0, -1)
    )

    options(gmultitasking = TRUE)
    parallel <- gquantiles("test.fixedbin",
        percentile = pcts,
        intervals = gintervals(c(1, 2), 0, -1)
    )

    expect_equal(unname(serial), unname(parallel))
})

# Regression: in the multitasking merge the parent gathered the kids' tail buffers
# into lowest_vals/highest_vals and then partial_sort()ed and resize()d them to
# gquantile.edge.data.size without clamping to how many values were actually
# gathered. A kid's tail buffer holds at most its own share of the stream, so
# whenever the edge buffer was bigger than the stream could fill - which the
# default of 100000 is for any stream shorter than that - partial_sort() ran past
# the end of the vector and resize() padded the tail with zeroes. The quantiles
# came back wrong (zeroes, and values read out of uninitialised memory beyond the
# vector) when the overrun stayed inside the allocation, and the process aborted in
# the allocator when it did not.
#
# Every edge size below must give the exact answer: each kid hands over its whole
# share of the stream, and the parent's two tail buffers between them cover every
# rank once the edge buffer holds at least half the stream, so no sampling enters
# the result at all. A smaller edge buffer is legitimately approximate in the
# middle and is checked separately, for its exact tails and its shape.
test_that("multitasking gquantiles is exact whatever gquantile.edge.data.size is", {
    local_db_state()
    db <- withr::local_tempdir()
    create_test_db(db, chrom_sizes = data.frame(chrom = paste0("chr", 1:8), size = 20000L))
    suppressMessages(gsetroot(db))

    ivs <- do.call(rbind, lapply(paste0("chr", 1:8), function(chr) {
        gintervals(chr, seq(0, 19990, by = 10), seq(10, 20000, by = 10))
    }))
    set.seed(17)
    vals <- seq(0.01, 1, length.out = nrow(ivs))[sample(nrow(ivs))]
    expect_equal(length(vals), 16000)
    gtrack.create_sparse("q", "scratch quantile track", ivs, vals)

    pcts <- c(0, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 1)

    # Reference computed outside misha, and confirmed against the path that does no
    # sub-sampling at all (the whole stream fits in the reservoir).
    reference <- unname(quantile(vals, pcts, type = 7))
    withr::with_options(
        list(gmultitasking = FALSE, gmax.data.size = 1e8),
        # tolerance: track values are stored as 32-bit floats
        expect_equal(unname(gquantiles("q", pcts, gintervals.all())), reference, tolerance = 1e-6)
    )

    # gmax.data.size = 100 forces every kid to sub-sample, which is the only path
    # that reaches the merge. Two scopes: over all 8 chromosomes the gathered tail
    # buffer used to end up exactly at its allocated capacity, so any overrun left
    # the allocation and aborted; over 5 chromosomes it had spare capacity, so the
    # overrun stayed inside it and the wrong answer was returned silently.
    for (nchrom in c(8, 5)) {
        scope <- gintervals(paste0("chr", 1:nchrom), 0, 20000)
        expected <- withr::with_options(
            list(gmax.data.size = 1e8),
            unname(quantile(gextract("q", scope, colnames = "v")$v, pcts, type = 7))
        )

        for (edge in c(9000, 12000, 16000, 20000, 100000)) {
            withr::with_options(
                list(
                    gmultitasking = TRUE, gmin.scope4process = 1, gmax.processes = 8,
                    gmax.data.size = 100, gquantile.edge.data.size = edge
                ),
                {
                    res <- unname(suppressWarnings(gquantiles("q", pcts, scope)))
                    expect_equal(res, expected,
                        info = sprintf("%d chromosomes, gquantile.edge.data.size = %d", nchrom, edge)
                    )
                }
            )
        }

        # An edge buffer too small to cover the middle ranks still has to give the
        # exact extremes and a monotone, in-range curve in between.
        withr::with_options(
            list(
                gmultitasking = TRUE, gmin.scope4process = 1, gmax.processes = 8,
                gmax.data.size = 100, gquantile.edge.data.size = 10
            ),
            {
                res <- unname(suppressWarnings(gquantiles("q", pcts, scope)))
                expect_equal(res[c(1, length(res))], expected[c(1, length(expected))])
                expect_false(is.unsorted(res))
                expect_true(all(res >= expected[1] & res <= expected[length(expected)]))
            }
        )
    }
})
