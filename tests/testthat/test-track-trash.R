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

test_that(".gdb.trash falls back to direct unlink when rename fails", {
    tmp <- withr::local_tempdir()
    target <- file.path(tmp, "renamefail")
    dir.create(target)
    file.create(file.path(target, "a"))

    local_mocked_bindings(
        file.rename = function(from, to) FALSE,
        .package = "base"
    )

    result <- .gdb.trash(target)
    expect_true(result)
    expect_false(file.exists(target))
    # No trash sibling because the fallback unlinks the original path directly
    expect_length(
        list.files(tmp, all.files = TRUE, no.. = TRUE, pattern = "^\\.trash\\."),
        0
    )
})

test_that(".gdb.trash called twice on same path returns FALSE second time", {
    tmp <- withr::local_tempdir()
    target <- file.path(tmp, "twice")
    dir.create(target)

    expect_true(.gdb.trash(target, async = FALSE))
    expect_false(file.exists(target))
    expect_false(.gdb.trash(target)) # gone, no-op
})

test_that(".gdb.trash returns FALSE when fallback unlink leaves residue", {
    tmp <- withr::local_tempdir()
    target <- file.path(tmp, "stubborn")
    dir.create(target)
    file.create(file.path(target, "child"))

    # Force the rename-fallback path and the unlink-keeps-file scenario.
    local_mocked_bindings(file.rename = function(from, to) FALSE, .package = "base")
    local_mocked_bindings(unlink = function(x, ...) 0L, .package = "base") # claim success but do nothing

    expect_false(.gdb.trash(target))
    expect_true(dir.exists(target)) # still there because unlink didn't actually remove it
})
