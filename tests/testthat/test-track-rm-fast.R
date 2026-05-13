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
