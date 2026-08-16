create_isolated_test_db()

# Three "peaks": the first two overlap (0-300 and 100-400), the third is
# separate. This is the shape of a unioned peak call or a promoter window
# around clustered TSSs.
overlapping_peaks <- function() gintervals(1, c(0, 100, 1000), c(300, 400, 1300))
disjoint_peaks <- function() gintervals(1, c(0, 400, 1000), c(300, 700, 1300))

test_that("an overlapping 1D iterator merges blocks and warns when the scope exposes it", {
    pk <- overlapping_peaks()
    expect_warning(
        res <- gextract("test.fixedbin", gintervals.all(), iterator = pk),
        "overlap"
    )
    # Behaviour is deliberately unchanged: the two overlapping peaks are
    # still reported as one merged block chr1:0-400.
    expect_equal(nrow(res), 2)
    expect_equal(res$start, c(0, 1000))
    expect_equal(res$end, c(400, 1300))
})

test_that("the warning names the merged block and the row-count collapse", {
    pk <- overlapping_peaks()
    w <- tryCatch(gextract("test.fixedbin", gintervals.all(), iterator = pk),
        warning = function(w) conditionMessage(w)
    )
    expect_match(w, "3 intervals were merged into 2", fixed = TRUE)
    expect_match(w, "chr1 0-300", fixed = TRUE)
    expect_match(w, "chr1 100-400", fixed = TRUE)
    expect_match(w, "chr1 0-400", fixed = TRUE)
    expect_match(w, "intervalID", fixed = TRUE)
})

test_that("passing the same intervals as the scope recovers one row per interval", {
    pk <- overlapping_peaks()
    expect_no_warning(res <- gextract("test.fixedbin", pk, iterator = pk))
    expect_equal(nrow(res), 3)
    # intervalID indexes the scope intervals, i.e. the peaks themselves
    expect_equal(res$intervalID, 1:3)
    expect_equal(res$start, pk$start)
    expect_equal(res$end, pk$end)
    # and each value is the mean over that peak alone, not over the merged block
    per_peak <- vapply(seq_len(nrow(pk)), function(i) {
        gsummary("test.fixedbin", pk[i, ])[["Mean"]]
    }, numeric(1))
    expect_equal(res$test.fixedbin, per_peak, tolerance = 1e-6)
})

test_that("a non-overlapping iterator is untouched: no warning, same rows", {
    nk <- disjoint_peaks()
    expect_no_warning(res <- gextract("test.fixedbin", gintervals.all(), iterator = nk))
    expect_equal(nrow(res), 3)
    expect_equal(res$start, nk$start)
    expect_equal(res$end, nk$end)
})

test_that("the merge warning is not limited to gextract", {
    pk <- overlapping_peaks()
    expect_warning(itr <- giterator.intervals("test.fixedbin", gintervals.all(), iterator = pk), "overlap")
    expect_equal(nrow(itr), 2)
    expect_warning(gsummary("test.fixedbin", gintervals.all(), iterator = pk), "overlap")
})

test_that("passing the intervals as the scope does not help where the scope is unified too", {
    pk <- overlapping_peaks()
    # Unlike gextract, giterator.intervals (and gscreen, gsummary, ...) unify their own
    # scope, so the merged block survives whatever scope is passed and the warning is
    # right to fire.
    expect_warning(itr <- giterator.intervals("test.fixedbin", pk, iterator = pk), "overlap")
    expect_equal(nrow(itr), 2)
    expect_equal(itr$end, c(400, 1300))
})

test_that("known gap: an interval nested in another one is dropped without a warning", {
    # The check is one-sided on purpose (see TrackExpressionScanner.cpp): it fires when a
    # row is wider than every input interval, which is what merging normally produces.
    # A nested interval disappears into its container without widening any row, so it is
    # still dropped silently. Pinned here so that closing the gap is a deliberate change.
    nested <- gintervals(1, c(0, 200), c(1000, 300))
    expect_no_warning(res <- gextract("test.fixedbin", gintervals.all(), iterator = nested))
    expect_equal(nrow(res), 1)
    expect_equal(res$end, 1000)
})

test_that("the warning is emitted once per call, not once per process", {
    pk <- overlapping_peaks()
    withr::local_options(gmultitasking = TRUE, gmin.scope4process = 1)
    ws <- character()
    withCallingHandlers(
        gextract("test.fixedbin", gintervals.all(), iterator = pk),
        warning = function(w) {
            ws <<- c(ws, conditionMessage(w))
            invokeRestart("muffleWarning")
        }
    )
    expect_equal(sum(grepl("were merged into", ws)), 1)
})
