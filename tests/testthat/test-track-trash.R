test_that(".gdb.trash renames atomically and unlinks in background", {
    tmp <- withr::local_tempdir()
    target <- file.path(tmp, "mydir")
    dir.create(target)
    file.create(file.path(target, "a"))
    file.create(file.path(target, "b"))

    result <- .gdb.trash(target)

    expect_true(result)
    expect_false(file.exists(target))
    siblings <- list.files(tmp, all.files = TRUE, no.. = TRUE)
    trash <- siblings[grepl("^\\.trash\\.mydir\\.", siblings)]
    expect_length(trash, 1)

    # Background unlink eventually clears it. Wait up to 10s.
    deadline <- Sys.time() + 10
    while (Sys.time() < deadline && dir.exists(file.path(tmp, trash))) {
        Sys.sleep(0.1)
    }
    expect_false(dir.exists(file.path(tmp, trash)))
})

test_that(".gdb.trash on missing path is a no-op returning FALSE", {
    tmp <- withr::local_tempdir()
    expect_false(.gdb.trash(file.path(tmp, "nope")))
})

test_that(".gdb.trash sync mode unlinks before returning", {
    tmp <- withr::local_tempdir()
    target <- file.path(tmp, "syncdir")
    dir.create(target)
    file.create(file.path(target, "x"))
    expect_true(.gdb.trash(target, async = FALSE))
    expect_length(list.files(tmp, all.files = TRUE, no.. = TRUE), 0)
})

test_that(".gdb.trash_sweep_old removes only stale .trash.* siblings", {
    tmp <- withr::local_tempdir()
    fresh <- file.path(tmp, ".trash.foo.1234.aaa")
    stale <- file.path(tmp, ".trash.bar.5678.bbb")
    dir.create(fresh)
    dir.create(stale)
    Sys.setFileTime(stale, Sys.time() - as.difftime(48, units = "hours"))
    n <- .gdb.trash_sweep_old(tmp, max_age_hours = 24)
    expect_equal(n, 1L)
    expect_true(dir.exists(fresh))
    expect_false(dir.exists(stale))
})

test_that(".gdb.trash_sweep_old on missing parent is a no-op", {
    tmp <- withr::local_tempdir()
    expect_equal(.gdb.trash_sweep_old(file.path(tmp, "nope")), 0L)
})
