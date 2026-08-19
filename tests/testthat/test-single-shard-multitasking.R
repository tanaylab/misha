create_isolated_test_db()

# A scope that lands entirely in one shard (a single chromosome, or a scope
# smaller than gmin.scope4process) makes prepare4multitasking() plan exactly one
# kid. Forking that kid buys nothing, so the multitask entry points fall back to
# running in-process. The failure this pins down is the tempting-but-wrong fix of
# having prepare4multitasking() return 0 for a one-kid plan: the callers read 0 as
# "empty scope" and every single-chromosome call would start returning NULL.

with_multitasking <- function(on, code) {
    old <- .ggetOption("gmultitasking")
    options(gmultitasking = on)
    on.exit(options(gmultitasking = old), add = TRUE)
    force(code)
}

single_chrom <- function() gintervals(1)
tiny_scope <- function() gintervals(1, 1000, 3000)
single_pair <- function() gintervals.2d(2, chroms2 = 2)

test_that("single-chromosome scopes return results, not NULL, under multitasking", {
    expect_false(is.null(with_multitasking(TRUE, gsummary("test.fixedbin", single_chrom()))))
    expect_false(is.null(with_multitasking(TRUE, gextract("test.fixedbin", single_chrom(), iterator = 1000))))
    expect_false(is.null(with_multitasking(TRUE, gquantiles("test.fixedbin", c(0.1, 0.5, 0.9), single_chrom()))))
    expect_false(is.null(with_multitasking(TRUE, gscreen("test.fixedbin > 0.2", single_chrom()))))
})

test_that("scopes below gmin.scope4process return results, not NULL, under multitasking", {
    expect_false(is.null(with_multitasking(TRUE, gsummary("test.fixedbin", tiny_scope()))))
    expect_false(is.null(with_multitasking(TRUE, gextract("test.fixedbin", tiny_scope(), iterator = 100))))
})

test_that("single chrom-pair 2D scopes return results, not NULL, under multitasking", {
    expect_false(is.null(with_multitasking(TRUE, gsummary("test.rects", single_pair()))))
    expect_false(is.null(with_multitasking(TRUE, gextract("test.rects", single_pair()))))
})

test_that("one-shard scopes give the same answer with and without multitasking", {
    cmp <- function(f) {
        expect_equal(with_multitasking(TRUE, f()), with_multitasking(FALSE, f()))
    }

    cmp(function() gsummary("test.fixedbin", single_chrom()))
    cmp(function() gsummary("test.fixedbin", single_chrom(), iterator = 100))
    cmp(function() gsummary("test.fixedbin", tiny_scope()))
    cmp(function() gsummary("test.rects", single_pair()))
    cmp(function() gextract("test.fixedbin", single_chrom(), iterator = 1000))
    cmp(function() gextract("test.fixedbin", tiny_scope(), iterator = 100))
    cmp(function() gextract("test.rects", single_pair()))
    cmp(function() gquantiles("test.fixedbin", c(0.1, 0.5, 0.9), single_chrom()))
    cmp(function() gintervals.quantiles("test.fixedbin", c(0.1, 0.5, 0.9), single_chrom()))
    cmp(function() gscreen("test.fixedbin > 0.2", single_chrom()))
    cmp(function() gdist("test.fixedbin", seq(0, 1, by = 0.1), intervals = single_chrom()))
    cmp(function() gcor("test.fixedbin", "test.fixedbin * 2", intervals = single_chrom()))
    cmp(function() gcor("test.fixedbin", "test.fixedbin * 2", intervals = single_chrom(), method = "spearman"))
    cmp(function() gcor("test.fixedbin", "test.fixedbin * 2", intervals = single_chrom(), method = "spearman.exact"))
    cmp(function() gintervals.summary("test.fixedbin", single_chrom()))
    cmp(function() gintervals.mapply(function(x) mean(x), "test.fixedbin", intervals = single_chrom(), iterator = 10000))
    cmp(function() glookup(c(10, 20, 30), "test.fixedbin", c(0, 0.3, 0.6, 1), intervals = single_chrom(), iterator = 100))
})

test_that("multi-chromosome scopes still agree with the serial path", {
    cmp <- function(f) {
        expect_equal(with_multitasking(TRUE, f()), with_multitasking(FALSE, f()))
    }

    all_chroms <- gintervals.all()
    cmp(function() gsummary("test.fixedbin", all_chroms))
    cmp(function() gextract("test.fixedbin", all_chroms, iterator = 10000))
    cmp(function() gquantiles("test.fixedbin", c(0.1, 0.5, 0.9), all_chroms))
    cmp(function() gscreen("test.fixedbin > 0.2", all_chroms))
})

test_that("an error inside the single-shard fallback leaves the session usable", {
    # The fallback runs the serial entry point, which reports errors to R through
    # Rf_error. Rf_error longjmps, so it must not be called while the multitasking
    # entry point still holds an RdbInitializer on the stack: its destructor would be
    # skipped and the library reference count would stay above zero, after which every
    # later multitasking call reuses stale shared memory and a stale kid index. This
    # showed up as "the whole session is broken after one failed call".
    all_chroms <- gintervals.all()
    reference <- gsummary("test.fixedbin", all_chroms)
    reference_rows <- nrow(gextract("test.fixedbin", all_chroms, iterator = 10000))

    for (i in 1:3) {
        withr::with_options(list(gmax.data.size = 100, gmultitasking = TRUE), {
            expect_error(gextract("test.fixedbin", gintervals(1), iterator = 10))
            expect_error(gscreen("test.fixedbin > 0", gintervals(1)))
        })

        # Multi-shard calls must still work, and still give the same answers.
        expect_equal(gsummary("test.fixedbin", all_chroms), reference)
        expect_equal(nrow(gextract("test.fixedbin", all_chroms, iterator = 10000)), reference_rows)
    }
})
