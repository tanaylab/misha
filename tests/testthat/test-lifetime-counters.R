# C6, the last of the six invariants in the C++/R lifetime contract: RdbInitializer's two
# counters - s_ref_count, the number of entry-point scopes alive, and s_protect_counter, the
# depth of misha's own protection stack - are back at their entry values when a .Call returns,
# on every path including error and interrupt. ~RdbInitializer is what restores them, so
# anything that skips it breaks the invariant: an Rf_error raised inside the scope, or an entry
# point called inline from another one. The consequences are not subtle - a stuck s_ref_count
# means no later constructor re-initialises the library and no later destructor tears it down,
# so the next multitasking call reads a stale shared-memory arena with a stale kid index.
#
# Until .glifetime_counters() (R/utils.R) existed the invariant was only reachable through side
# effects of the teardown, mainly the process umask, which tests something adjacent to the
# claim rather than the claim, and became a fragile stand-in when 5.11.21 made misha's umask an
# option.

restore_groot_on_exit()

test_that("both counters are zero between calls", {
    expect_identical(.glifetime_counters(), c(ref_count = 0L, protect_count = 0L))
})

# Without this the file's other assertions would pass just as well against a stub that always
# returned zero. A track expression is evaluated by R from inside the entry point's scope, so
# the reading taken there must show the scope that is holding it.
test_that("a reading taken inside an entry point shows the live scope", {
    inner <- NULL
    probe <- function(x) {
        inner <<- .glifetime_counters()
        x
    }
    expect_s3_class(gextract("probe(test.fixedbin)", gintervals(1, 0, 1000)), "data.frame")
    expect_gte(inner[["ref_count"]], 1L)
    expect_identical(.glifetime_counters(), c(ref_count = 0L, protect_count = 0L))
})

test_that("a successful entry point leaves both counters where it found them", {
    before <- .glifetime_counters()
    scope <- gintervals(1, 0, 10000)

    invisible(gextract("test.fixedbin", scope))
    invisible(gscreen("test.fixedbin > 0.5", scope))
    invisible(gsummary("test.sparse", scope))
    invisible(gtrack.info("test.fixedbin"))

    expect_identical(.glifetime_counters(), before)
})

test_that("an entry point that fails leaves both counters where it found them", {
    before <- .glifetime_counters()
    scope <- gintervals(1, 0, 10000)

    expect_error(gextract("test.fixedbin + .no_such_object_", scope))
    expect_error(gscreen("test.fixedbin > 0.5 & .no_such_object_", scope))
    expect_error(gextract("no.such.track.at.all", scope))

    expect_identical(.glifetime_counters(), before)
})

test_that("a multitasking entry point leaves both counters where it found them", {
    withr::local_options(gmultitasking = TRUE, gmax.data.size = 1e9)
    before <- .glifetime_counters()
    scope <- gintervals(1:2)

    invisible(gextract("test.sparse", scope))
    expect_error(gextract("test.sparse + .no_such_object_", scope))

    expect_identical(.glifetime_counters(), before)
})
