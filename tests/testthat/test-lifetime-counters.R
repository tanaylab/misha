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

# The next misha call made outside any other one repairs a scope that a longjmp left behind
# (.gcall, RdbInitializer::repair_leaked_scope), so a reading taken after a sequence of calls
# shows only the last one. The tests below read the counters after each call.

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
    expect_identical(.glifetime_counters(), before)
    invisible(gscreen("test.fixedbin > 0.5", scope))
    expect_identical(.glifetime_counters(), before)
    invisible(gsummary("test.sparse", scope))
    expect_identical(.glifetime_counters(), before)
    invisible(gtrack.info("test.fixedbin"))
    expect_identical(.glifetime_counters(), before)
})

test_that("an entry point that fails leaves both counters where it found them", {
    before <- .glifetime_counters()
    scope <- gintervals(1, 0, 10000)

    expect_error(gextract("test.fixedbin + .no_such_object_", scope))
    expect_identical(.glifetime_counters(), before)
    expect_error(gscreen("test.fixedbin > 0.5 & .no_such_object_", scope))
    expect_identical(.glifetime_counters(), before)
    expect_error(gextract("no.such.track.at.all", scope))
    expect_identical(.glifetime_counters(), before)
})

test_that("a multitasking entry point leaves both counters where it found them", {
    withr::local_options(gmultitasking = TRUE, gmax.data.size = 1e9)
    before <- .glifetime_counters()
    scope <- gintervals(1:2)

    invisible(gextract("test.sparse", scope))
    expect_identical(.glifetime_counters(), before)
    expect_error(gextract("test.sparse + .no_such_object_", scope))
    expect_identical(.glifetime_counters(), before)
})

# The number of processes that evaluated the track expression of a gsummary() that may multitask
count_pids <- function() {
    pf <- tempfile()
    on.exit(unlink(pf))
    pidlog <- function(x) {
        cat(Sys.getpid(), "\n", file = pf, append = TRUE)
        x
    }
    withr::with_options(list(gmultitasking = TRUE, gmax.processes = 4, gmin.scope4process = 100), {
        gsummary("pidlog(test.fixedbin)", gintervals(c(1, 2), 0, 1e6))
    })
    length(unique(scan(pf, quiet = TRUE)))
}

# gintervals.rbind of a factor and a character chrom column raised an R error inside the .Call
# (SET_STRING_ELT on the factor's integer column). The longjmp skipped ~RdbInitializer, and from
# 5.16.0 every later call of the session then ran in one process, with the umask left at 007.
test_that("gintervals.rbind of a factor and a character chrom leaves the session as it was", {
    umask <- Sys.umask(NA)
    expect_gt(count_pids(), 1)
    r <- gintervals.rbind(gintervals(1, 0, 10), data.frame(chrom = "chr1", start = 5, end = 20))
    expect_equal(r$start, c(0, 5))
    expect_identical(.glifetime_counters(), c(ref_count = 0L, protect_count = 0L))
    expect_identical(Sys.umask(NA), umask)
    expect_gt(count_pids(), 1)
})

test_that("a misha call in an argument of another one is not nested in it", {
    # gseq.kmer.dist passes its intervals argument to .gcall unevaluated; the gscreen in it ran
    # at the raised depth, so in one process
    pf <- tempfile()
    on.exit(unlink(pf))
    pidlog <- function(x) {
        cat(Sys.getpid(), "\n", file = pf, append = TRUE)
        x
    }
    withr::with_options(list(gmultitasking = TRUE, gmax.processes = 4, gmin.scope4process = 100), {
        gseq.kmer.dist(gscreen("pidlog(test.fixedbin) > 0.2", gintervals(c(1, 2), 0, 1e6)), k = 2)
    })
    expect_gt(length(unique(scan(pf, quiet = TRUE))), 1)
})

test_that("a scope left behind by a longjmp is put back by the next misha call", {
    umask <- Sys.umask(NA)
    # C_test_leak_scope raises an R error inside its scope on purpose
    expect_error(.Call("C_test_leak_scope"), "R error inside an RdbInitializer scope")
    expect_identical(.glifetime_counters()[["ref_count"]], 1L)

    withr::with_options(list(gmultitasking = TRUE), expect_true(.ggetOption("gmultitasking")))
    expect_gt(count_pids(), 1)
    expect_identical(.glifetime_counters(), c(ref_count = 0L, protect_count = 0L))
    expect_identical(Sys.umask(NA), umask)
})

test_that("a direct grbind of columns of different types fails without leaving its scope", {
    a <- data.frame(chrom = factor("chr1"), start = 0, end = 10)
    b <- data.frame(chrom = "chr1", start = 5, end = 20)
    expect_error(.gcall("grbind", list(a, b), .misha_env()), "Cannot copy a data frame column of type character into a column of type factor")
    expect_identical(.glifetime_counters(), c(ref_count = 0L, protect_count = 0L))
    b$chrom <- factor("chr1", levels = c("chr2", "chr1"))
    expect_error(.gcall("grbind", list(a, b), .misha_env()), "differ in the levels of factor column 1")
    expect_identical(.glifetime_counters(), c(ref_count = 0L, protect_count = 0L))
})
