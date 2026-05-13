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

test_that(".gtrack.create_atomic trashes a half-written tmp dir on failure", {
    suppressWarnings(gtrack.rm("atomic_partial", force = TRUE))
    withr::defer(suppressWarnings(gtrack.rm("atomic_partial", force = TRUE)))

    parent <- get("GWD", envir = .misha)

    expect_error(
        .gtrack.create_atomic("atomic_partial", function() {
            # Simulate what C++ would do: read+clear the override slot,
            # mkdir the tmp dir, write some partial files, then fail.
            tmp_dir <- get(".create_dir_override", envir = .misha)
            if (is.null(tmp_dir)) stop("override slot was not set")
            assign(".create_dir_override", NULL, envir = .misha) # mimic C++ clear-on-read
            dir.create(tmp_dir, recursive = TRUE)
            writeLines("partial1", file.path(tmp_dir, "chr1"))
            writeLines("partial2", file.path(tmp_dir, "chr2"))
            stop("simulated mid-create write failure")
        }),
        regexp = "simulated mid-create write failure"
    )

    # Final dir should not exist
    expect_false(dir.exists(file.path(parent, "atomic_partial.track")))

    # Tmp dir should be gone (either deleted or renamed-to-trash)
    tmp_left <- list.files(parent,
        all.files = TRUE,
        pattern = "^\\.atomic_partial\\.track\\.tmp\\."
    )
    expect_length(tmp_left, 0)

    # Wait for any in-flight trash unlink to settle (best-effort, async)
    Sys.sleep(0.5)
})

test_that(".gtrack.create_atomic loses gracefully when another session wins the rename race", {
    suppressWarnings(gtrack.rm("race_target", force = TRUE))
    withr::defer(suppressWarnings(gtrack.rm("race_target", force = TRUE)))

    parent <- get("GWD", envir = .misha)
    final_dir <- file.path(parent, "race_target.track")

    # Pre-create the final dir as if another session won the race after
    # .gconfirmtrackcreate but before our rename.
    suppressWarnings(unlink(final_dir, recursive = TRUE, force = TRUE))
    withr::defer(suppressWarnings(unlink(final_dir, recursive = TRUE, force = TRUE)))

    # Our create_fn writes a tmp dir; the helper should then refuse to
    # rename onto the existing final dir.
    create_fn <- function() {
        tmp_dir <- get(".create_dir_override", envir = .misha)
        if (is.null(tmp_dir)) stop("override slot was not set")
        assign(".create_dir_override", NULL, envir = .misha)
        dir.create(tmp_dir, recursive = TRUE)
        writeLines("data", file.path(tmp_dir, "x"))
        # NOW simulate the racing session creating the final dir.
        dir.create(final_dir)
        writeLines("winner", file.path(final_dir, "x"))
    }

    expect_error(
        .gtrack.create_atomic("race_target", create_fn),
        regexp = "Refusing to overwrite existing|already exists"
    )

    # Our tmp dir is gone (trashed); the winner's dir survives.
    tmp_left <- list.files(parent,
        all.files = TRUE,
        pattern = "^\\.race_target\\.track\\.tmp\\."
    )
    expect_length(tmp_left, 0)
    expect_true(dir.exists(final_dir))
    expect_equal(readLines(file.path(final_dir, "x")), "winner")
})
