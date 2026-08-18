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

test_that("a nested interval is a merge too, and warns", {
    # A nested input widens no row - the block equals the interval it is nested in - so the
    # "is every row covered by an input interval" half of the test never sees it. It is the
    # other half that catches it: the nested peak is emitted nowhere on its own.
    # This is at least as common as partial overlap in the case the warning is for, a union
    # of peak calls from two conditions: a broad peak in A containing a narrow one in B.
    nested <- gintervals(1, c(0, 200), c(1000, 300))

    expect_warning(res <- gextract("test.fixedbin", gintervals.all(), iterator = nested), "overlap")
    expect_equal(nrow(res), 1)
    expect_equal(res$end, 1000)

    expect_warning(itr <- giterator.intervals("test.fixedbin", nested, iterator = nested), "overlap")
    expect_equal(nrow(itr), 1)

    expect_warning(scr <- gscreen("test.fixedbin > -1e9", nested, iterator = nested), "overlap")
    expect_equal(nrow(scr), 1)
})

test_that("the nested warning names the interval that disappeared", {
    nested <- gintervals(1, c(0, 200), c(1000, 300))
    w <- tryCatch(gextract("test.fixedbin", gintervals.all(), iterator = nested),
        warning = function(w) conditionMessage(w)
    )
    expect_match(w, "2 intervals were merged into 1", fixed = TRUE)
    expect_match(w, "chr1 200-300", fixed = TRUE)
    expect_match(w, "single row chr1 0-1000", fixed = TRUE)
})

test_that("an interval nested only after the scope clips it warns as well", {
    # Neither interval is nested in the input here; the scope is what makes them so, and the
    # first one is still swallowed by the merged block.
    clipped <- gintervals(1, c(0, 50), c(100, 150))
    expect_warning(res <- gextract("test.fixedbin", gintervals(1, 0, 100), iterator = clipped), "overlap")
    expect_equal(nrow(res), 1)
    expect_equal(res$end, 100)
})

test_that("remaining gap: an exact duplicate under a wider scope is still silent", {
    # Both copies are emitted verbatim, as the same single row, so no interval is missing
    # from the output even though two went in and one came out. Pinned so that closing this
    # last gap - which needs the opt-in strict mode, not a wider merge test - is deliberate.
    dup <- gintervals(1, c(0, 0), c(300, 300))
    expect_no_warning(res <- gextract("test.fixedbin", gintervals.all(), iterator = dup))
    expect_equal(nrow(res), 1)
    expect_equal(res$end, 300)
})

test_that("the recommended idiom stays silent for every overlap shape", {
    # gextract(t, pk, iterator = pk) is what the warning tells the user to do, so it has to
    # stay silent for the shapes a unioned peak call actually produces - and stay correct:
    # one row per peak, equal to that peak, valued over that peak alone.
    shapes <- list(
        partial = gintervals(1, c(0, 100), c(300, 400)),
        nested = gintervals(1, c(0, 100), c(1000, 300)),
        duplicate = gintervals(1, c(0, 0), c(300, 300)),
        chained = gintervals(1, c(0, 100, 200), c(200, 300, 400))
    )

    for (nm in names(shapes)) {
        pk <- shapes[[nm]]
        expect_no_warning(res <- gextract("test.fixedbin", pk, iterator = pk))
        expect_equal(nrow(res), nrow(pk), info = nm)
        expect_equal(res$start, pk$start, info = nm)
        expect_equal(res$end, pk$end, info = nm)
        expect_equal(res$intervalID, seq_len(nrow(pk)), info = nm)

        per_peak <- vapply(seq_len(nrow(pk)), function(i) {
            gsummary("test.fixedbin", pk[i, ])[["Mean"]]
        }, numeric(1))
        expect_equal(res$test.fixedbin, per_peak, tolerance = 1e-6, info = nm)
    }
})

test_that("the merge test survives a big-set scope, which it walks twice", {
    withr::defer(gintervals.rm("test.overlap_bigset_scope", force = TRUE))

    # A scope stored as a GIntervalsBigSet1D is read off disk, so the second walk the nested
    # case needs is a second full pass and a begin_iter() re-entry rather than a rewind of a
    # vector. Both halves of the test have to survive that.
    scope <- gintervals(1, c(0, 5000), c(2000, 6000))
    withr::with_options(list(gmax.data.size = 1), {
        gintervals.save("test.overlap_bigset_scope", scope)
    })
    expect_true(misha:::.gintervals.is_bigset("test.overlap_bigset_scope"))

    nested <- gintervals(1, c(0, 200), c(1000, 300))
    expect_warning(
        res <- gextract("test.fixedbin", "test.overlap_bigset_scope", iterator = nested),
        "overlap"
    )
    expect_equal(nrow(res), 1)
    expect_equal(res$end, 1000)

    disjoint <- gintervals(1, c(0, 5000), c(1000, 5500))
    expect_no_warning(gextract("test.fixedbin", "test.overlap_bigset_scope", iterator = disjoint))
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

test_that("a disjoint outer iterator does not silence an overlapping inner one", {
    withr::local_options(gmultitasking = FALSE)
    pk <- overlapping_peaks()
    nk <- disjoint_peaks()
    # The once-per-call key is consumed only by a call that has something to report.
    # The outer iterator here does not overlap, so it must leave the key for the inner
    # call, whose iterator does - otherwise the merge that actually happened is silent.
    #
    # The handler has to be installed inside FUN: gintervals.mapply evaluates FUN
    # through R_tryEval, which resets R's handler stack, so a handler around the
    # gintervals.mapply call never sees a condition raised inside it.
    seen <- character()
    gintervals.mapply(
        function(x) {
            withCallingHandlers(
                gsummary("test.fixedbin", gintervals(1, 0, 1000), iterator = pk)[["Mean"]],
                warning = function(w) {
                    seen <<- c(seen, conditionMessage(w))
                    invokeRestart("muffleWarning")
                }
            )
        },
        "test.fixedbin", gintervals(1, 0, 1300),
        iterator = nk
    )
    expect_equal(sum(grepl("were merged into", seen)), 1)
})

# A queued diagnostic describes the call that queued it. .gcall has to keep that true
# across an error or an interrupt, which is what the on.exit in it is for.

# gmax.data.size is verified after the expression iterator has been built, so the
# overlap warning is already queued by the time the call dies.
too_big_for_the_buffer <- function(expr) {
    withr::with_options(list(gmax.data.size = 1), expr)
}

test_that("a call that errors after queueing a diagnostic leaves nothing behind in .misha", {
    pk <- overlapping_peaks()
    expect_error(
        too_big_for_the_buffer(gextract("test.fixedbin", gintervals.all(), iterator = pk)),
        "exceeded the maximum"
    )
    expect_equal(grep("^\\.GPENDING", ls(.misha, all.names = TRUE), value = TRUE), character(0))
    expect_no_warning(gextract("test.fixedbin", gintervals(1, 0, 1000)))
})

test_that("the outer call's own diagnostic survives an error inside FUN", {
    withr::local_options(gmultitasking = FALSE)
    pk <- overlapping_peaks()
    # The inner call takes the outer queue aside while it runs; erroring must not be a
    # way to lose it. This is the shape a user writes without thinking about it:
    # tryCatch inside FUN so one bad interval does not abort the whole apply.
    expect_warning(
        gintervals.mapply(
            function(x) {
                tryCatch(
                    too_big_for_the_buffer(gextract("test.fixedbin", gintervals(1, 0, 1000))),
                    error = function(e) NA_real_
                )
                1
            },
            "test.fixedbin", gintervals(1, 0, 1300),
            iterator = pk
        ),
        "were merged into"
    )
})

test_that("a diagnostic queued by an inner call that errored is not raised by the outer one", {
    withr::local_options(gmultitasking = FALSE)
    pk <- overlapping_peaks()
    nk <- disjoint_peaks()
    # The outer iterator does not overlap, so the outer call has nothing to report. The
    # inner one does, but it errors before returning, and the user handled that error -
    # blaming the outer call for the inner call's merge names the wrong iterator.
    expect_no_warning(
        gintervals.mapply(
            function(x) {
                tryCatch(
                    too_big_for_the_buffer(
                        gextract("test.fixedbin", gintervals(1, 0, 1300), iterator = pk)
                    ),
                    error = function(e) NA_real_
                )
                1
            },
            "test.fixedbin", gintervals(1, 0, 1300),
            iterator = nk
        )
    )
    expect_equal(grep("^\\.GPENDING", ls(.misha, all.names = TRUE), value = TRUE), character(0))
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
