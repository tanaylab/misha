# When the shared-memory arena cannot be allocated, gextract_multitask() and
# gscreen_multitask() give up on multitasking and run their serial entry point instead.
# That serial run must happen *outside* the try block that owns the multitasking
# RdbInitializer (rdb::SingleShard, src/rdbutils.h): an error inside it reaches R through
# Rf_error, and the longjmp would skip the outer destructor, leaving misha's s_ref_count
# stuck above zero. From then on no constructor re-initialises the library and no
# destructor tears it down, so the second multitasking call afterwards reads a stale
# shared-memory arena with a stale kid index - wrong rows, then a dead child.
#
# The sibling path - the scope fitting in a single shard - is covered by
# test-single-shard-multitasking.R; this file covers the shared-memory failure.
#
# Forcing the allocation to fail, with no debug hook: over a sparse-track iterator neither
# entry point can estimate the record count, so distribute_task() sizes the arena from
# gmax.data.size alone. gmax.data.size = 1e13 asks mmap() for a few hundred terabytes, far
# past the 128 TB user address space, so the allocation fails outright - and the
# inflation-factor retry ladder, which only runs when there *is* an estimate, is skipped.
#
# The scope must span several chromosomes: with a single shard the entry points take the
# same tail for an unrelated reason and the fallback branch is never reached. gintervals(1:4)
# plans four children.

# Re-root at a database of this file's own. Under the parallel runner a file inherits
# whatever GROOT the previous file in the same worker left behind, and gintervals(1:4)
# then fails before the first expectation runs.
create_isolated_test_db()

SHM_FAIL_MAX_DATA_SIZE <- 1e13

test_that("gextract falls back to the serial path when shared memory cannot be allocated", {
    withr::local_options(gmultitasking = TRUE, gmax.data.size = 1e9)
    scope <- gintervals(1:4)
    expected <- gextract("test.sparse", scope)

    withr::local_options(gmax.data.size = SHM_FAIL_MAX_DATA_SIZE)
    expect_equal(gextract("test.sparse", scope), expected)
})

test_that("gscreen falls back to the serial path when shared memory cannot be allocated", {
    withr::local_options(gmultitasking = TRUE, gmax.data.size = 1e9)
    scope <- gintervals(1:4)
    expected <- gscreen("test.sparse > 0.5", scope)

    withr::local_options(gmax.data.size = SHM_FAIL_MAX_DATA_SIZE)
    expect_equal(gscreen("test.sparse > 0.5", scope), expected)
})

test_that("an error inside the serial fallback leaves the session usable", {
    withr::local_options(gmultitasking = TRUE, gmax.data.size = 1e9)
    scope <- gintervals(1:4)
    expected_extract <- gextract("test.sparse", scope)
    expected_screen <- gscreen("test.sparse > 0.5", scope)

    # The umask is the one piece of RdbInitializer's teardown that R can observe: the
    # outermost constructor sets umask(07), and only the matching destructor puts the old
    # value back. Pick a value misha's own cannot be mistaken for.
    old_umask <- Sys.umask("022")
    withr::defer(Sys.umask(old_umask))

    local({
        withr::local_options(gmax.data.size = SHM_FAIL_MAX_DATA_SIZE)
        expect_error(gextract("test.sparse + .no_such_object_", scope), "not found")
        expect_error(gscreen("test.sparse > 0.5 & .no_such_object_", scope), "not found")
    })

    expect_identical(Sys.umask(NA), as.octmode("022"))

    # The calls below are only meaningful once the destructor is known to have run - and
    # only safe: with s_ref_count stuck, the second multitasking call after the leak
    # indexes the previous call's arena out of bounds, and the child that dies leaves
    # through R's own shutdown, whose R_CleanTempDir() deletes the session tempdir and the
    # test database inside it. Guard them rather than take the rest of the file down.
    if (identical(Sys.umask(NA), as.octmode("022"))) {
        for (i in 1:2) {
            expect_equal(gextract("test.sparse", scope), expected_extract)
            expect_equal(gscreen("test.sparse > 0.5", scope), expected_screen)
        }
    }
})
