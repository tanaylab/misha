# misha runs some work on a few threads (run_in_threads, src/WorkerThreads.h). A thread that
# cannot be started must leave its share to the calling thread, not end R through
# std::terminate. A process limit forces it: with RLIMIT_NPROC at 1 the kernel refuses every
# new thread (and fork) of the user. Each case runs in a fresh R process, which sets the limit
# on itself once the session is ready.

skip_if_no_nproc_limit <- function() {
    skip_on_os(c("mac", "windows", "solaris"))
    skip_if(!nzchar(Sys.which("prlimit")), "prlimit is not available")
    # root is exempt from RLIMIT_NPROC
    skip_if(identical(Sys.info()[["effective_user"]], "root"), "root ignores the process limit")
}

run_limited_r <- function(groot, setup, limited) {
    script <- withr::local_tempfile(fileext = ".R")
    writeLines(c(
        sprintf("library(misha, lib.loc = %s)", deparse(dirname(find.package("misha")))),
        sprintf("gsetroot(%s)", deparse(groot)),
        setup,
        "invisible(system(sprintf('prlimit --pid %d --nproc=1', Sys.getpid())))",
        # the limit must hold: a fork now fails
        "if (suppressWarnings(system('true', ignore.stderr = TRUE)) == 0) { cat('NO_LIMIT\\n'); quit(save = 'no') }",
        limited
    ), script)
    out <- suppressWarnings(system2(file.path(R.home("bin"), "R"), c("--no-echo", "--no-save", "-f", shQuote(script)),
        stdout = TRUE, stderr = TRUE
    ))
    skip_if(any(grepl("^NO_LIMIT", out)), "the process limit did not take")
    out
}

test_that("intervals_join = 'intervals' works when no gather thread can start", {
    skip_if_no_nproc_limit()
    out <- run_limited_r(
        get("GROOT", envir = .misha),
        c(
            "options(gmultitasking = FALSE, gmax.processes = 8L)",
            "iv <- gintervals(c(1, 2), 0, 20e6)",
            "iv$label <- c('a', 'b')",
            "iv$k <- 1:2",
            # 400,000 rows and 5 columns of intervals: gathered on 5 threads
            "expected <- gextract('test.fixedbin', intervals = iv, iterator = 100, intervals_join = 'intervals')"
        ),
        c(
            "res <- gextract('test.fixedbin', intervals = iv, iterator = 100, intervals_join = 'intervals')",
            "cat('IDENTICAL', identical(res, expected), nrow(res), '\\n')"
        )
    )
    expect_true(any(grepl("^IDENTICAL TRUE 4e\\+05|^IDENTICAL TRUE 400000", out)), info = paste(out, collapse = "\n"))
})
