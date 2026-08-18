create_isolated_test_db()

# The C++ layer cannot raise a warning itself: an exiting handler - tryCatch(warning = ),
# or options(warn = 2) - longjmps out of the .Call, skipping ~RdbInitializer, so misha's
# reference count never returns to zero and its PROTECT counter is never restored. The
# text is queued instead (rdb::add_pending_diagnostic) and .gcall() raises it after the
# call has returned. These are the warnings that used to be raised from C++ directly.

# A reservoir of 100 values and 10-value tails overflow on a 2000-bin scope, so the
# quantile sampling warning fires in a fraction of a second instead of on a genome scan.
sampled <- function(expr) {
    withr::with_options(list(gmax.data.size = 100, gquantile.edge.data.size = 10), expr)
}

overlapping_peaks <- function() gintervals(1, c(0, 100, 1000), c(300, 400, 1300))

# How many times the once-per-call iterator-overlap warning fires in one call. It is the
# canary for a leaked reference count: the key is cleared by the outermost RdbInitializer,
# so a call that longjmped past its destructor leaves the key set for the rest of the
# session and every later call is silent.
overlap_warnings <- function() {
    ws <- character()
    withCallingHandlers(
        gextract("test.fixedbin", gintervals(1, 0, 2000), iterator = overlapping_peaks()),
        warning = function(w) {
            ws <<- c(ws, conditionMessage(w))
            invokeRestart("muffleWarning")
        }
    )
    sum(grepl("were merged into", ws))
}

pending_diagnostics <- function() {
    grep("^\\.GPENDING", ls(.misha, all.names = TRUE), value = TRUE)
}

test_that("the quantile sampling warnings still say what they said", {
    withr::local_options(gmultitasking = FALSE)
    expect_warning(
        sampled(gquantiles("test.fixedbin", 0.5, gintervals(1, 0, 100000))),
        "Data size \\(2000\\) exceeds the limit \\(100\\)"
    )
    expect_warning(
        sampled(gintervals.quantiles("test.fixedbin", 0.5, gintervals(1, c(0, 200000), c(100000, 300000)))),
        "Data size in one or more intervals exceeds the limit \\(100\\)"
    )
    expect_warning(
        sampled(gbins.quantiles("test.fixedbin", c(0, 1, 2),
            expr = "test.fixedbin", percentiles = 0.5, intervals = gintervals(1, 0, 100000)
        )),
        "Data size in one or more intervals exceeds the limit \\(100\\)"
    )
    expect_match(
        tryCatch(
            sampled(gquantiles("test.fixedbin", 0.5, gintervals(1, 0, 100000))),
            warning = function(w) conditionMessage(w)
        ),
        "The data was sampled to fit the limit and the resulted quantiles are hence approximate.\n(The limit can be controlled by gmax.data.size limit)",
        fixed = TRUE
    )
})

test_that("the multitasking quantile warnings say what they said, once per call", {
    withr::local_options(gmultitasking = TRUE, gmin.scope4process = 1)
    scope <- gintervals(c(1, 2, 3), 0, 100000)

    ws <- character()
    withCallingHandlers(sampled(gquantiles("test.fixedbin", 0.5, scope)),
        warning = function(w) {
            ws <<- c(ws, conditionMessage(w))
            invokeRestart("muffleWarning")
        }
    )
    expect_equal(sum(grepl("The data was sampled to fit the limit", ws)), 1)

    ws <- character()
    withCallingHandlers(sampled(gintervals.quantiles("test.fixedbin", 0.5, scope)),
        warning = function(w) {
            ws <<- c(ws, conditionMessage(w))
            invokeRestart("muffleWarning")
        }
    )
    expect_equal(sum(grepl("Data size in one or more intervals exceeds the limit", ws)), 1)
})

test_that("catching a quantile warning leaves the PROTECT stack and the call counter intact", {
    withr::local_options(gmultitasking = FALSE)
    # The warning is caught by an exiting handler, which is what used to longjmp past
    # ~RdbInitializer. Everything after it is the evidence that nothing was left behind:
    # the pending queue is empty, later calls still work, and - the visible symptom of a
    # leaked reference count - the once-per-call overlap warning still fires per call.
    expect_equal(overlap_warnings(), 1)

    caught <- tryCatch(
        sampled(gquantiles("test.fixedbin", 0.5, gintervals(1, 0, 100000))),
        warning = function(w) conditionMessage(w)
    )
    expect_match(caught, "Data size \\(2000\\) exceeds the limit \\(100\\)")

    expect_equal(pending_diagnostics(), character(0))
    expect_equal(nrow(gextract("test.fixedbin", gintervals(1, 0, 1000))), 20)
    expect_equal(overlap_warnings(), 1)
    expect_equal(overlap_warnings(), 1)
})

test_that("options(warn = 2) on a quantile warning leaves the PROTECT stack intact", {
    withr::local_options(gmultitasking = FALSE)
    expect_equal(overlap_warnings(), 1)

    expect_error(
        withr::with_options(
            list(warn = 2),
            sampled(gquantiles("test.fixedbin", 0.5, gintervals(1, 0, 100000)))
        ),
        "Data size \\(2000\\) exceeds the limit \\(100\\)"
    )

    expect_equal(pending_diagnostics(), character(0))
    expect_equal(nrow(gextract("test.fixedbin", gintervals(1, 0, 1000))), 20)
    expect_equal(overlap_warnings(), 1)
})

test_that("a multitasking quantile warning survives options(warn = 2)", {
    # Every forked worker reaches the same warning as the parent. A worker cannot deliver
    # one - its stderr goes to /dev/null - but under warn = 2 it used to turn the warning
    # into an error inside the fork and kill the call.
    withr::local_options(gmultitasking = TRUE, gmin.scope4process = 1)
    scope <- gintervals(c(1, 2, 3), 0, 100000)

    expect_error(
        withr::with_options(list(warn = 2), sampled(gintervals.quantiles("test.fixedbin", 0.5, scope))),
        "Data size in one or more intervals exceeds the limit"
    )

    expect_equal(pending_diagnostics(), character(0))
    expect_equal(nrow(gextract("test.fixedbin", gintervals(1, 0, 1000))), 20)
    expect_equal(overlap_warnings(), 1)
})

# ---------------------------------------------------------------------------------
# The zero-length interval warning, from the interval converter. It is the one site
# that fires more than once in a single call: one warning per named intervals set the
# call loads, which is what the pending queue's append path is for.

zero_length_set <- function(name, n_zero = 2) {
    intervs <- data.frame(
        chrom = factor(rep("chr1", n_zero + 1), levels = gintervals.all()$chrom),
        start = c(100, seq(500, by = 100, length.out = n_zero)),
        end = c(200, seq(500, by = 100, length.out = n_zero))
    )
    # gintervals.save() bumps zero-length intervals on the way out, so the only way to
    # get one on disk - which is how they reach users, from an old set or another tool -
    # is to serialize it directly.
    f <- file(file.path(get("GWD", envir = .misha), paste0(name, ".interv")), "wb")
    on.exit(close(f))
    serialize(intervs, f)
}

zero_length_set("zerolen")
zero_length_set("zerolen2", n_zero = 3)
gdb.reload(rescan = TRUE)

zero_length_warnings <- function(expr) {
    ws <- character()
    withCallingHandlers(force(expr),
        warning = function(w) {
            ws <<- c(ws, conditionMessage(w))
            invokeRestart("muffleWarning")
        }
    )
    ws[grepl("start == end", ws)]
}

test_that("the zero-length interval warning still says what it said", {
    withr::local_options(gmultitasking = FALSE)
    expect_equal(
        zero_length_warnings(gextract("test.fixedbin", "zerolen"))[1],
        "2 interval(s) had start == end and were extended by 1bp on load. misha does not support zero-length intervals; if this is unintended, fix the source data."
    )
})

test_that("one warning per intervals set the call loads, in the order they were loaded", {
    withr::local_options(gmultitasking = FALSE)
    # Two different sets in one call: the queue has to hold both, and keep their counts
    # apart - a single slot per severity would have dropped one of them.
    expect_equal(
        zero_length_warnings(gintervals.intersect("zerolen", "zerolen2")),
        c(
            "2 interval(s) had start == end and were extended by 1bp on load. misha does not support zero-length intervals; if this is unintended, fix the source data.",
            "3 interval(s) had start == end and were extended by 1bp on load. misha does not support zero-length intervals; if this is unintended, fix the source data."
        )
    )
    expect_equal(pending_diagnostics(), character(0))
})

test_that("catching the zero-length warning leaves the PROTECT stack intact", {
    withr::local_options(gmultitasking = FALSE)
    expect_equal(overlap_warnings(), 1)

    caught <- tryCatch(gextract("test.fixedbin", "zerolen"),
        warning = function(w) conditionMessage(w)
    )
    expect_match(caught, "start == end", fixed = TRUE)

    expect_equal(pending_diagnostics(), character(0))
    expect_equal(nrow(gextract("test.fixedbin", gintervals(1, 0, 1000))), 20)
    expect_equal(overlap_warnings(), 1)

    expect_error(
        withr::with_options(list(warn = 2), gextract("test.fixedbin", "zerolen")),
        "start == end"
    )
    expect_equal(pending_diagnostics(), character(0))
    expect_equal(overlap_warnings(), 1)
})

test_that("forked workers add nothing to what the zero-length warning already said", {
    # Every forked worker re-converts the iterator, and so hits this warning; the parent
    # converts it too, before it forks. Only the parent's conversions can reach the user,
    # so splitting the scope over more workers must not change what the caller sees.
    withr::local_options(gmultitasking = TRUE)
    scope <- gintervals(c(1, 2, 3), 0, 100000)

    many_kids <- withr::with_options(
        list(gmin.scope4process = 1),
        zero_length_warnings(gextract("test.fixedbin", scope, iterator = "zerolen"))
    )
    few_kids <- withr::with_options(
        list(gmin.scope4process = 1e9),
        zero_length_warnings(gextract("test.fixedbin", scope, iterator = "zerolen"))
    )
    expect_equal(many_kids, few_kids)
    expect_true(length(many_kids) >= 1)
    expect_true(all(many_kids == many_kids[1]))
    expect_equal(pending_diagnostics(), character(0))
})

test_that("a track-parallel worker's warning still reaches the caller", {
    # gextract's mclapply path converts the iterator inside each worker, which is not a
    # misha multitasking kid: the worker raises the queued text itself, its own
    # withCallingHandlers collects it, and the parent re-raises the distinct ones.
    withr::local_options(gmultitasking = TRUE, gmultitasking.strategy = "tracks")
    ws <- zero_length_warnings(
        gextract(c("test.fixedbin", "test.sparse"), gintervals(1, 0, 100000), iterator = "zerolen")
    )
    expect_equal(unique(ws), "2 interval(s) had start == end and were extended by 1bp on load. misha does not support zero-length intervals; if this is unintended, fix the source data.")
})

test_that("a warning queued by a call that then errors is not raised by a later one", {
    withr::local_options(gmultitasking = FALSE)
    # .gcall() drops the queue of a call that did not return. The zero-length warning is
    # queued while the scope is loaded, well before the result size is checked.
    expect_error(
        withr::with_options(
            list(gmax.data.size = 1),
            gextract("test.fixedbin", "zerolen")
        ),
        "exceeded the maximum"
    )
    expect_equal(pending_diagnostics(), character(0))
    expect_no_warning(gextract("test.fixedbin", gintervals(1, 0, 1000)))
})

# ---------------------------------------------------------------------------------
# Nested .gcall(): gintervals.mapply evaluates FUN inside its own .Call, so a misha
# call made from FUN nests inside a call that has already queued diagnostics of its
# own. This is what .gcall()'s take-the-outer-queue-aside-and-put-it-back exists for.
#
# The inner warnings have to be collected inside FUN: FUN is evaluated through
# R_tryEval, which runs it in a fresh top-level context, so a withCallingHandlers
# wrapped around the outer call does not see them.

zerolen_texts <- function() {
    list(
        outer = "2 interval(s) had start == end and were extended by 1bp on load. misha does not support zero-length intervals; if this is unintended, fix the source data.",
        inner = "3 interval(s) had start == end and were extended by 1bp on load. misha does not support zero-length intervals; if this is unintended, fix the source data."
    )
}

test_that("a nested .gcall raises only its own warnings, and the outer keeps its own", {
    withr::local_options(gmultitasking = FALSE)
    # Without the take-aside, the inner call raises the outer set's warning as well -
    # once per interval, since the outer queue would still be there each time - and the
    # outer call ends up raising whatever the last inner call left behind.
    txt <- zerolen_texts()
    n_intervals <- nrow(gintervals.load("zerolen"))
    inner_ws <- character()

    outer_ws <- zero_length_warnings(
        gintervals.mapply(
            function(x) {
                withCallingHandlers(
                    gsummary("test.fixedbin", "zerolen2"),
                    warning = function(w) {
                        inner_ws <<- c(inner_ws, conditionMessage(w))
                        invokeRestart("muffleWarning")
                    }
                )
                mean(x)
            },
            "test.fixedbin",
            intervals = "zerolen"
        )
    )

    # Each inner call raises its own set's warning and nothing else...
    expect_equal(inner_ws, rep(txt$inner, n_intervals))
    # ...and the outer call still raises its own, exactly once, after FUN has run.
    expect_equal(outer_ws, txt$outer)
    expect_equal(pending_diagnostics(), character(0))
})

test_that("an inner .gcall that errors drops its queue and leaves the outer one alone", {
    withr::local_options(gmultitasking = FALSE)
    # The inner call queues the zero-length warning while it loads the scope, then dies
    # on the result-size limit. Its queue goes with it; the outer call's does not.
    txt <- zerolen_texts()
    inner_ws <- character()

    outer_ws <- zero_length_warnings(
        gintervals.mapply(
            function(x) {
                withCallingHandlers(
                    try(
                        withr::with_options(
                            list(gmax.data.size = 1),
                            gextract("test.fixedbin", "zerolen2")
                        ),
                        silent = TRUE
                    ),
                    warning = function(w) {
                        inner_ws <<- c(inner_ws, conditionMessage(w))
                        invokeRestart("muffleWarning")
                    }
                )
                mean(x)
            },
            "test.fixedbin",
            intervals = "zerolen"
        )
    )

    expect_equal(inner_ws, character(0))
    expect_equal(outer_ws, txt$outer)
    expect_equal(pending_diagnostics(), character(0))
})
