create_isolated_test_db()

test_that("sliding window sum has no precision drift", {
    # Regression test for Kahan compensated summation fix in the sliding
    # window sum of GenomeTrackFixedBin.
    #
    # Before the fix, the incremental add/subtract sliding window could
    # accumulate floating-point drift (~1e-7) relative to a fresh sum.
    # With Kahan summation, the sliding result must be bit-identical to
    # a from-scratch computation for the same window.
    #
    # Test strategy:
    #   1. Create a vtrack with sum + wide shifts on a dense track.
    #   2. Extract with a contiguous iterator (triggers sliding window).
    #   3. Extract the same windows individually (each computed fresh).
    #   4. Assert the two are bit-identical (not merely close).

    withr::defer({
        for (vt in gvtrack.ls()) {
            gvtrack.rm(vt)
        }
    })

    # test.fixedbin is a dense track with 50bp bins.
    # Shifts of -2500/+2500 create a 5000bp window covering 100 bins,
    # so each sliding step adds 1 bin and removes 1 bin.
    gvtrack.create("vt_sliding", "test.fixedbin", func = "sum")
    gvtrack.iterator("vt_sliding", sshift = -2500, eshift = 2500)

    # --- Sliding extraction (contiguous iterator) ---
    query <- gintervals(1, 0, 200000)
    sliding <- gextract("vt_sliding", query, iterator = 50, colnames = "sliding_sum")

    # --- Fresh per-window extraction ---
    # Use a sample of positions spread across the range so any drift
    # that grows with distance from the start would be caught.
    sample_starts <- seq(0, 199950, by = 50)
    sample_ends <- sample_starts + 50
    sample_intervals <- gintervals(
        rep(1, length(sample_starts)), sample_starts, sample_ends
    )
    per_window <- gextract(
        "vt_sliding", sample_intervals,
        iterator = sample_intervals, colnames = "fresh_sum"
    )

    # Merge on genomic coordinates
    merged <- merge(sliding, per_window, by = c("chrom", "start", "end"))
    both_valid <- !is.na(merged$sliding_sum) & !is.na(merged$fresh_sum)

    expect_true(
        sum(both_valid) > 100,
        info = "Need enough non-NA windows to test precision"
    )

    # The sliding and fresh sums must be exactly equal -- no tolerance.
    # Any difference here means the sliding window accumulator has drifted.
    expect_identical(
        merged$sliding_sum[both_valid],
        merged$fresh_sum[both_valid]
    )
})

test_that("windows of two sizes slide and give what each window read afresh gives", {
    # A 60 bp iterator over test.fixedbin's 50 bp bins gives windows of two sizes, which slide
    # from one to the next (dropping and adding different numbers of bins) for sum, exists and
    # size, and for sum with avg, while the sum stays exact; lse slides between windows of one
    # size. Rows 600 bp apart are each read afresh. (lse with avg is left out: there a window read
    # afresh sums the LSE in float, a slid one in double.)
    withr::defer({
        for (vt in gvtrack.ls()) {
            gvtrack.rm(vt)
        }
    })
    sets <- list("sum", "lse", "exists", "size", c("sum", "lse", "exists", "size"), c("sum", "avg"), c("exists", "size"))
    end1 <- gintervals.all()$end[1]
    query <- rbind(gintervals(1, 0, 300000), gintervals(1, end1 - 30000, end1))
    for (shifts in list(c(-140, 173), c(-1000, 1000), c(7, 290))) {
        for (funcs in sets) {
            for (vt in gvtrack.ls()) {
                gvtrack.rm(vt)
            }
            vts <- paste0("vt_", funcs)
            for (i in seq_along(funcs)) {
                gvtrack.create(vts[i], "test.fixedbin", func = funcs[i])
                gvtrack.iterator(vts[i], sshift = shifts[1], eshift = shifts[2])
            }
            slid <- gextract(vts, query, iterator = 60)
            rows <- seq(1, nrow(slid), by = 10)
            far <- slid[rows, c("chrom", "start", "end")]
            fresh <- gextract(vts, far, iterator = far)
            label <- paste(paste(funcs, collapse = "+"), paste(shifts, collapse = ","))
            expect_gt(sum(!is.na(fresh[[vts[1]]])), 100)
            for (vt in vts) {
                expect_identical(slid[[vt]][rows], fresh[[vt]], label = paste(label, vt))
            }
        }
    }
})
