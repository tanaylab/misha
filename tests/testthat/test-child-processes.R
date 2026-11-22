load_test_db()

# misha collects only the child processes it started. It used to wait with waitpid(-1) at the
# end of every call, which also reaped any other finished child of the R process, and that
# child's owner then lost its exit status. A pipe() connection is such a child: close()
# returns its exit status.

# The status close() gives for a child that ended while `f` ran.
pipe_status_after <- function(f) {
    con <- pipe("exit 3", "r")
    readLines(con) # returns once the child closed its output, i.e. as it ends
    Sys.sleep(0.2)
    f()
    close(con)
}

test_that("a misha call leaves the exit status of other child processes to their owner", {
    expected <- pipe_status_after(function() NULL)
    expect_gt(expected, 0)

    expect_equal(pipe_status_after(function() gextract("test.fixedbin", gintervals(1, 0, 1000))), expected)
    withr::with_options(list(gmultitasking = TRUE), {
        expect_equal(pipe_status_after(function() gsummary("test.fixedbin")), expected)
    })
})
