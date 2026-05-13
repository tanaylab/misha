create_isolated_test_db()

test_that(".gtrack.create_atomic writes through a tmp dir then renames", {
    suppressWarnings(gtrack.rm("atomic_ok", force = TRUE))
    withr::defer(suppressWarnings(gtrack.rm("atomic_ok", force = TRUE)))

    gtrack.create("atomic_ok", "test", "test.fixedbin")
    expect_true(gtrack.exists("atomic_ok"))

    parent <- dirname(.track_dir("atomic_ok"))
    # No leftover tmp siblings (trash is async; ignore .trash.* entries that
    # happen to embed a former .tmp. basename).
    all_entries <- list.files(parent, all.files = TRUE)
    tmp_leftovers <- all_entries[
        grepl("\\.tmp\\.", all_entries) & !grepl("^\\.trash\\.", all_entries)
    ]
    expect_length(tmp_leftovers, 0)
})

test_that("interrupted create leaves no visible track and no FS dir", {
    suppressWarnings(gtrack.rm("atomic_fail", force = TRUE))
    withr::defer(suppressWarnings(gtrack.rm("atomic_fail", force = TRUE)))

    # Reference a non-existent track to make C++ throw mid-creation
    expect_error(
        gtrack.create("atomic_fail", "test", "no_such_track_xyz"),
        regexp = "no_such_track_xyz|undefined|not found|does not exist",
        ignore.case = TRUE
    )

    expect_false(gtrack.exists("atomic_fail"))
    final_dir <- file.path(get("GWD", envir = .misha), "atomic_fail.track")
    expect_false(dir.exists(final_dir))

    gdb.reload(rescan = TRUE)
    expect_false("atomic_fail" %in% gtrack.ls())

    parent <- get("GWD", envir = .misha)
    # No lingering tmp siblings (a .trash.* rename is fine, async unlink
    # completes off-thread).
    all_entries <- list.files(parent, all.files = TRUE)
    tmp_leftovers <- all_entries[
        grepl("\\.atomic_fail\\.track\\.tmp\\.", all_entries) &
            !grepl("^\\.trash\\.", all_entries)
    ]
    expect_length(tmp_leftovers, 0)
})

test_that(".gconfirmtrackcreate still blocks re-create of an existing track", {
    suppressWarnings(gtrack.rm("atomic_dup", force = TRUE))
    withr::defer(suppressWarnings(gtrack.rm("atomic_dup", force = TRUE)))

    gtrack.create("atomic_dup", "test", "test.fixedbin")
    expect_error(
        gtrack.create("atomic_dup", "test", "test.fixedbin"),
        "already exists"
    )
})

test_that("tmp dir is invisible to gfind_tracks_n_intervals", {
    parent <- get("GWD", envir = .misha)
    fake <- file.path(parent, ".faketrack.tmp.1234.abc")
    suppressWarnings(unlink(fake, recursive = TRUE, force = TRUE))
    dir.create(fake)
    withr::defer(suppressWarnings(unlink(fake, recursive = TRUE, force = TRUE)))

    res <- .gcall("gfind_tracks_n_intervals", parent, .misha_env())
    expect_false(any(grepl("faketrack", res[[1]])))
})

test_that(".gtrack.create_atomic cleans up after C++ throws", {
    # Direct call to the helper with a create_fn that throws.
    suppressWarnings(gtrack.rm("atomic_throw", force = TRUE))
    withr::defer(suppressWarnings(gtrack.rm("atomic_throw", force = TRUE)))

    expect_error(
        .gtrack.create_atomic("atomic_throw", function() {
            stop("simulated mid-create failure")
        }),
        regexp = "simulated mid-create failure"
    )
    expect_false(dir.exists(file.path(get("GWD", envir = .misha), "atomic_throw.track")))

    parent <- get("GWD", envir = .misha)
    # tmp dir should be gone (.gdb.trash renamed it to .trash.* and unlinks
    # async). A trash entry is fine; a raw tmp is not.
    all_entries <- list.files(parent, all.files = TRUE)
    tmp_leftovers <- all_entries[
        grepl("\\.atomic_throw\\.track\\.tmp\\.", all_entries) &
            !grepl("^\\.trash\\.", all_entries)
    ]
    expect_length(tmp_leftovers, 0)
})
