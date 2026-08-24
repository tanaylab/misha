load_test_db()

# Ctrl-C during a misha call, and what happens to the session afterwards.
#
# The mechanism these tests guard: RdbInitializer displaces R's SIGINT handler for the
# duration of a .Call, so a raw Rf_error inside that scope longjmps past the destructor
# that puts R's handler back - and the session's Ctrl-C is dead from then on, printing
# "CTL-C!" and doing nothing. Every test here therefore checks the session as well as
# the call.

# TRUE when R's own SIGINT handling is intact: a signal to ourselves raises an interrupt
# condition. If misha left its handler installed, the signal only sets misha's flag and
# Sys.sleep() runs to completion.
r_interrupt_works <- function() {
    tryCatch(
        {
            tools::pskill(Sys.getpid(), tools::SIGINT)
            Sys.sleep(1)
            FALSE
        },
        interrupt = function(i) TRUE
    )
}

test_that("an unusable gbuf.size is reported and leaves the session intact", {
    scope <- gintervals(1, 0, 100000)
    n <- nrow(gextract("test.fixedbin", scope))

    for (bad in list(2^32, 2^31, NA_real_, NA_integer_, Inf, 0, -5)) {
        withr::with_options(list(gbuf.size = bad), {
            expect_error(gextract("test.fixedbin", scope), "gbuf.size option must be a finite number")
        })
    }
    for (bad in list("abc", c(1000, 2000), TRUE)) {
        withr::with_options(list(gbuf.size = bad), {
            expect_error(gextract("test.fixedbin", scope), "gbuf.size option must be a single number")
        })
    }

    # The point of the fix: the rejected option used to reach R as a raw Rf_error from
    # VECTOR_ELT(NULL, ...), which never unwound misha's SIGINT handler.
    expect_equal(nrow(gextract("test.fixedbin", scope)), n)
    expect_true(r_interrupt_works())
})

test_that("gbuf.size accepts an integer and a real alike", {
    scope <- gintervals(1, 0, 100000)
    ref <- gextract("test.fixedbin", scope)
    for (good in list(1, 1L, 1000, 5000L, 1e5)) {
        withr::with_options(list(gbuf.size = good), {
            expect_equal(gextract("test.fixedbin", scope), ref)
        })
    }
})

test_that("an interrupt inside a track expression aborts the call, not the session", {
    # Deterministic: the expression signals us on its first evaluation, and the scanner's
    # next check_interrupt() throws. No timing involved.
    fired <- new.env(parent = emptyenv())
    fired$done <- FALSE
    e <- new.env(parent = globalenv())
    e$sigme <- function(x) {
        if (!fired$done) {
            fired$done <- TRUE
            tools::pskill(Sys.getpid(), tools::SIGINT)
        }
        x
    }
    # Single process on purpose: the expression has to run in this process for the signal
    # and the flag it sets to be ours.
    withr::with_options(list(gmultitasking = FALSE), {
        expect_error(
            eval(quote(gextract("sigme(test.fixedbin)", gintervals.all(), iterator = 10000)), e),
            "Command interrupted"
        )
    })
    expect_true(fired$done)

    expect_equal(nrow(gextract("test.fixedbin", gintervals(1, 0, 100000))), 2000)
    expect_true(r_interrupt_works())
})

test_that("Ctrl-C stops a running gseq.pwm scan", {
    skip_on_cran()
    # Timing-based, but only in one direction: the workload is far larger than the delay,
    # and if the machine still outran the signal the test skips rather than failing.
    delay <- 2
    n <- 40000
    seqs <- vapply(
        seq_len(n),
        function(i) paste(sample(c("A", "C", "G", "T"), 4000, TRUE), collapse = ""),
        character(1)
    )
    pssm <- matrix(rep(c(0.7, 0.1, 0.1, 0.1), 15), ncol = 4, byrow = TRUE)
    colnames(pssm) <- c("A", "C", "G", "T")

    sender <- sprintf(
        "sleep %g; kill -INT %d",
        delay, Sys.getpid()
    )
    t0 <- Sys.time()
    system2("/bin/sh", c("-c", shQuote(sender)), wait = FALSE)
    res <- withr::with_options(
        list(gmultitasking = FALSE),
        tryCatch(
            {
                gseq.pwm(seqs, pssm, mode = "max")
                "completed"
            },
            error = function(e) conditionMessage(e)
        )
    )
    elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    # Absorb a signal that landed after the call: outside a misha call it is R's, and an
    # uncaught interrupt would take the rest of the file with it.
    tryCatch(Sys.sleep(max(0, delay - elapsed) + 0.5), interrupt = function(i) NULL)

    if (identical(res, "completed") && elapsed < delay) {
        skip(sprintf("gseq.pwm finished in %.1fs, before the %gs signal", elapsed, delay))
    }
    # Completing despite a signal delivered mid-scan is the defect, not a reason to skip.
    expect_match(res, "Command interrupted")
    expect_lt(elapsed, delay + 5)

    expect_equal(nrow(gextract("test.fixedbin", gintervals(1, 0, 100000))), 2000)
    expect_true(r_interrupt_works())
})

test_that("gtrack.lookup reports a readable error when it fails before the first chromosome", {
    # The message used to be built with an uninitialised stack buffer, so it carried
    # arbitrary bytes and could not even be passed to nchar().
    gdir.create("temp", showWarnings = FALSE)
    # The expression goes in as a string: a bare symbol is resolved by .gexpr2str() in R
    # and never reaches the C++ scanner. This one parses, names no track, and so throws
    # while the iterator policy is being inferred - before any chromosome file is opened.
    msg <- tryCatch(
        {
            gtrack.lookup("temp.interrupt_lookup", "d", c(1, 2), "no_such_track_at_all", c(0, 1, 2))
            NA_character_
        },
        error = function(e) conditionMessage(e)
    )
    expect_false(is.na(msg))
    expect_true(validUTF8(msg))
    expect_silent(nchar(msg))
    expect_false(grepl("Error writing", msg, fixed = TRUE))
    expect_match(msg, "no_such_track_at_all")
})
