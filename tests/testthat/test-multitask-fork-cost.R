# Each multitasking process is a fork of the R session, and a fork copies the page tables
# of the session's memory. From a session holding many GB in 4 KB pages that takes
# seconds per call, so misha warns once per session (rdb::warn_if_forks_are_slow).
#
# Once per session means once per R process, so each case runs in a fresh one. The
# environment variable MISHA_TEST_FORK_4K_BYTES stands in for the bytes the session holds
# in 4 KB pages: forcing the real condition takes GB of memory and a host whose kernel
# gives R no huge pages.

run_in_fresh_r <- function(lines, env = character()) {
    script <- withr::local_tempfile(fileext = ".R")
    writeLines(c(
        sprintf("library(misha, lib.loc = %s)", deparse(dirname(find.package("misha")))),
        sprintf("gsetroot(%s)", deparse(get("GROOT", envir = .misha))),
        "options(gmax.processes = 4L, gmax.processes2core = 4L, gmin.scope4process = 1L, gmultitasking = TRUE)",
        "fork_warnings <- function(expr) {",
        "    n <- 0L",
        "    withCallingHandlers(expr, warning = function(w) {",
        "        if (grepl('4 KB pages', conditionMessage(w))) {",
        "            n <<- n + 1L",
        "            cat('TEXT', conditionMessage(w), '\\n')",
        "        }",
        "        invokeRestart('muffleWarning')",
        "    })",
        "    n",
        "}",
        lines
    ), script)
    system2(file.path(R.home("bin"), "R"), c("--no-echo", "--no-save", "-f", shQuote(script)),
        stdout = TRUE, stderr = TRUE, env = env
    )
}

counts <- function(out, key) as.integer(sub(paste0("^", key, " "), "", grep(paste0("^", key, " "), out, value = TRUE)))

test_that("a session holding much memory in 4 KB pages is warned once", {
    skip_if(parallel::detectCores() < 2L, "Multitasking test needs at least 2 cores")
    out <- run_in_fresh_r(c(
        "cat('FIRST', fork_warnings(gsummary('test.fixedbin', gintervals(c(1, 2)))), '\\n')",
        "cat('SECOND', fork_warnings(gsummary('test.fixedbin', gintervals(c(1, 2)))), '\\n')",
        "cat('EXTRACT', fork_warnings(gextract('test.fixedbin', gintervals(c(1, 2)), iterator = 1e5)), '\\n')"
    ), env = "MISHA_TEST_FORK_4K_BYTES=64000000000")
    expect_equal(counts(out, "FIRST"), 1L, info = paste(out, collapse = "\n"))
    expect_equal(counts(out, "SECOND"), 0L)
    expect_equal(counts(out, "EXTRACT"), 0L)
    text <- grep("^TEXT ", out, value = TRUE)
    expect_match(text, "This R session holds 64.0 GB of memory in 4 KB pages", fixed = TRUE)
    expect_match(text, "gmax.processes", fixed = TRUE)
    expect_match(text, "once per session", fixed = TRUE)
})

test_that("no fork-cost warning for a small session, few GB, or no multitasking", {
    skip_if(parallel::detectCores() < 2L, "Multitasking test needs at least 2 cores")
    call <- "fork_warnings(gsummary('test.fixedbin', gintervals(c(1, 2))))"
    # this session, as it is: small
    out <- run_in_fresh_r(sprintf("cat('N', %s, '\\n')", call))
    expect_equal(counts(out, "N"), 0L, info = paste(out, collapse = "\n"))
    # 8 GB in 4 KB pages with 4 processes: forks of about a second, below the threshold
    out <- run_in_fresh_r(sprintf("cat('N', %s, '\\n')", call), env = "MISHA_TEST_FORK_4K_BYTES=8000000000")
    expect_equal(counts(out, "N"), 0L)
    # nothing is forked
    out <- run_in_fresh_r(c("options(gmultitasking = FALSE)", sprintf("cat('N', %s, '\\n')", call)),
        env = "MISHA_TEST_FORK_4K_BYTES=64000000000"
    )
    expect_equal(counts(out, "N"), 0L)
})
