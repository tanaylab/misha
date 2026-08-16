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

test_that("known gap: an interval nested after clipping is dropped without a warning", {
    # The check is one-sided on purpose (see TrackExpressionScanner.cpp): it fires when a
    # row is wider than every input interval, which is what merging normally produces.
    # An interval that ends up nested inside another one - in the input, or only after
    # the scope clips it - disappears without widening any row, so it is still dropped
    # silently. Pinned here so that closing the gap is a deliberate change.
    nested <- gintervals(1, c(0, 200), c(1000, 300))
    expect_no_warning(res <- gextract("test.fixedbin", gintervals.all(), iterator = nested))
    expect_equal(nrow(res), 1)
    expect_equal(res$end, 1000)

    # neither interval is nested in the input here; the scope is what makes them so
    clipped <- gintervals(1, c(0, 50), c(100, 150))
    expect_no_warning(res <- gextract("test.fixedbin", gintervals(1, 0, 100), iterator = clipped))
    expect_equal(nrow(res), 1)
    expect_equal(res$end, 100)
})

test_that("the warning survives gextract's track-parallel path", {
    withr::local_options(gmultitasking = TRUE, gmax.data.size = 1e9)
    st <- seq(0, by = 10000, length.out = 1200)
    scope <- gintervals(1, st, st + 10000)
    pk <- gintervals(1, c(0, 100, 20000), c(300, 400, 20300))
    exprs <- sprintf("test.fixedbin + %d", 0:8)

    # mclapply workers never deliver their warnings, so the parent has to re-raise
    # them. Assert the dispatch really goes track-parallel: on any other path this
    # test would pass without exercising anything.
    expect_equal(.gmultitasking_strategy(exprs, scope, pk, NULL, NULL, NULL), "tracks")

    ws <- character()
    res <- withCallingHandlers(
        do.call(gextract, c(as.list(exprs), list(intervals = scope, iterator = pk))),
        warning = function(w) {
            ws <<- c(ws, conditionMessage(w))
            invokeRestart("muffleWarning")
        }
    )
    expect_equal(sum(grepl("were merged into", ws)), 1)
    expect_equal(nrow(res), 2)
    expect_equal(ncol(res), length(exprs) + 4L)
})

test_that("a misha call inside gintervals.mapply's FUN does not swallow the warning", {
    withr::local_options(gmultitasking = FALSE)
    pk <- overlapping_peaks()
    expect_warning(
        gintervals.mapply(
            function(x) gsummary("test.fixedbin", gintervals(1, 0, 1000))[["Mean"]],
            "test.fixedbin", gintervals.all(),
            iterator = pk
        ),
        "overlap"
    )
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
