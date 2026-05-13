test_that("gtrack.rm returns quickly and removes track", {
    skip_if_not(dir.exists("/net/mraid20/export/tgdata/db/tgdb/misha_test_db"))
    gsetroot("/net/mraid20/export/tgdata/db/tgdb/misha_test_db")
    on.exit(gdb.reload(), add = TRUE)

    # Create a track to remove
    suppressWarnings(gtrack.rm("test.rm_fast", force = TRUE))
    gtrack.create("test.rm_fast", "tmp", "test.fixedbin")
    on.exit(suppressWarnings(gtrack.rm("test.rm_fast", force = TRUE)), add = TRUE)

    trackdir <- .track_dir("test.rm_fast")
    expect_true(dir.exists(trackdir))

    t <- system.time(gtrack.rm("test.rm_fast", force = TRUE))
    expect_lt(t["elapsed"], 2)
    expect_false(dir.exists(trackdir))
    expect_false("test.rm_fast" %in% gtrack.ls())
})

test_that("gtrack.rm with force on non-existent track is silent", {
    skip_if_not(dir.exists("/net/mraid20/export/tgdata/db/tgdb/misha_test_db"))
    gsetroot("/net/mraid20/export/tgdata/db/tgdb/misha_test_db")
    expect_silent(gtrack.rm("nope_definitely_no_such_track", force = TRUE))
})

test_that("gtrack.rm uses .gdb.trash (creates trash sibling)", {
    skip_if_not(dir.exists("/net/mraid20/export/tgdata/db/tgdb/misha_test_db"))
    gsetroot("/net/mraid20/export/tgdata/db/tgdb/misha_test_db")
    on.exit(gdb.reload(), add = TRUE)

    suppressWarnings(gtrack.rm("test.rm_trash", force = TRUE))
    gtrack.create("test.rm_trash", "tmp", "test.fixedbin")
    on.exit(suppressWarnings(gtrack.rm("test.rm_trash", force = TRUE)), add = TRUE)

    trackdir <- .track_dir("test.rm_trash")
    parent <- dirname(trackdir)

    # Stub the async cleanup so the trash sibling is observable.
    local_mocked_bindings(
        system = function(command, ...) 0L,
        .package = "base"
    )

    gtrack.rm("test.rm_trash", force = TRUE)
    siblings <- list.files(parent,
        all.files = TRUE, no.. = TRUE,
        pattern = "^\\.trash\\.rm_trash\\.track\\."
    )
    expect_length(siblings, 1)
    # Clean up the trash sibling synchronously since we mocked the async path
    unlink(file.path(parent, siblings), recursive = TRUE, force = TRUE)
})

test_that("gtrack.rm force-branch reports failure when trash returns FALSE", {
    skip_if_not(dir.exists("/net/mraid20/export/tgdata/db/tgdb/misha_test_db"))
    gsetroot("/net/mraid20/export/tgdata/db/tgdb/misha_test_db")
    on.exit(gdb.reload(), add = TRUE)

    # Create a track and then drop it from GTRACKS so gtrack.rm enters the
    # force-branch (track not in GTRACKS but dir on disk).
    suppressWarnings(gtrack.rm("test_rm_failreport", force = TRUE))
    gtrack.create("test_rm_failreport", "tmp", "test.fixedbin")
    on.exit(suppressWarnings(gtrack.rm("test_rm_failreport", force = TRUE)),
        add = TRUE
    )

    trackdir <- .track_dir("test_rm_failreport")
    expect_true(dir.exists(trackdir))

    # Force the force-branch: remove from GTRACKS while the dir is still on
    # disk.
    gtracks <- get("GTRACKS", envir = .misha)
    assign("GTRACKS", setdiff(gtracks, "test_rm_failreport"), envir = .misha)

    local_mocked_bindings(.gdb.trash = function(...) FALSE)

    expect_message(
        gtrack.rm("test_rm_failreport", force = TRUE),
        regexp = "Failed to delete"
    )
})

test_that("gtrack.rm force-branch with no on-disk residue scrubs cache silently", {
    skip_if_not(dir.exists("/net/mraid20/export/tgdata/db/tgdb/misha_test_db"))
    gsetroot("/net/mraid20/export/tgdata/db/tgdb/misha_test_db")
    on.exit(gdb.reload(), add = TRUE)

    # Track that was never on disk but is also not in GTRACKS -> existed=FALSE
    # so .gdb.rm_track must still run (it's a cache scrub). The "Failed" message
    # must NOT appear.
    expect_silent(gtrack.rm("test_rm_neverexisted", force = TRUE))
})

test_that("gtrack.copy overwrite aborts when trash fails", {
    skip_if_not(dir.exists("/net/mraid20/export/tgdata/db/tgdb/misha_test_db"))
    gsetroot("/net/mraid20/export/tgdata/db/tgdb/misha_test_db")
    on.exit(gdb.reload(), add = TRUE)

    suppressWarnings(gtrack.rm("test_copy_src", force = TRUE))
    suppressWarnings(gtrack.rm("test_copy_dst", force = TRUE))
    gtrack.create("test_copy_src", "src", "test.fixedbin")
    gtrack.create("test_copy_dst", "dst", "test.fixedbin")
    on.exit(suppressWarnings({
        gtrack.rm("test_copy_src", force = TRUE)
        gtrack.rm("test_copy_dst", force = TRUE)
    }), add = TRUE)

    local_mocked_bindings(.gdb.trash = function(...) FALSE)

    expect_error(
        gtrack.copy("test_copy_src", "test_copy_dst", overwrite = TRUE),
        regexp = "Failed to remove existing destination|copy aborted"
    )
})
