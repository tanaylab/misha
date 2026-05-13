create_isolated_test_db()

test_that(".misha$.create_dir_override redirects track create dir", {
    # Clean slate
    suppressWarnings(gtrack.rm("test_override_target", force = TRUE))
    withr::defer(suppressWarnings(gtrack.rm("test_override_target", force = TRUE)))

    parent <- file.path(get("GWD", envir = .misha), ".test_override_tmp_dir")
    suppressWarnings(unlink(parent, recursive = TRUE, force = TRUE))
    withr::defer(suppressWarnings(unlink(parent, recursive = TRUE, force = TRUE)))

    # Sanity: default destination should not pre-exist.
    default_path <- file.path(get("GWD", envir = .misha), "test_override_target.track")
    expect_false(dir.exists(default_path))

    # Set the override; create_track_dir should mkdir THIS path instead
    # of the default test_override_target.track. The slot should be
    # cleared after consumption.
    assign(".create_dir_override", parent, envir = .misha)

    # The C++ writes into 'parent' (our override) but the trackname is
    # "test_override_target". R-level post-processing (.gdb.add_track,
    # attrs) likely fails because the expected track dir does not exist.
    # Catch any resulting error - we only care about the mkdir side
    # effect and slot consumption.
    suppressWarnings(try(
        gtrack.create("test_override_target", "test", "test.fixedbin"),
        silent = TRUE
    ))

    # Slot was consumed (cleared)
    override_after <- get(".create_dir_override", envir = .misha, inherits = FALSE)
    expect_true(is.null(override_after))

    # Override path was mkdir'd and has the per-chrom files
    expect_true(dir.exists(parent))
    expect_gt(length(list.files(parent)), 0)

    # Default path was NOT created
    expect_false(dir.exists(default_path))

    # Cache may still have a stale entry; reload.
    gdb.reload(rescan = TRUE)
})

test_that(".misha$.create_dir_override = NULL behaves as no override", {
    suppressWarnings(gtrack.rm("test_override_off", force = TRUE))
    withr::defer(suppressWarnings(gtrack.rm("test_override_off", force = TRUE)))

    assign(".create_dir_override", NULL, envir = .misha)

    gtrack.create("test_override_off", "test", "test.fixedbin")
    expect_true(gtrack.exists("test_override_off"))
})

test_that(".misha$.create_dir_override = empty string behaves as no override", {
    suppressWarnings(gtrack.rm("test_override_empty", force = TRUE))
    withr::defer(suppressWarnings(gtrack.rm("test_override_empty", force = TRUE)))

    # Length-0 character vector should be ignored (not a usable path)
    assign(".create_dir_override", character(0), envir = .misha)

    gtrack.create("test_override_empty", "test", "test.fixedbin")
    expect_true(gtrack.exists("test_override_empty"))
})

test_that(".create_dir_override is cleared even when create_track_dir errors", {
    # Pre-create the override path so create_track_dir's mkdir errors
    # with EEXIST. The C++ reads + clears the override BEFORE calling
    # mkdir, so the slot must be cleared even though the call errors.
    parent <- file.path(get("GWD", envir = .misha), ".override_mkdir_fail_target")
    suppressWarnings(unlink(parent, recursive = TRUE, force = TRUE))
    dir.create(parent)
    withr::defer(suppressWarnings(unlink(parent, recursive = TRUE, force = TRUE)))

    suppressWarnings(gtrack.rm("test_override_err_target", force = TRUE))
    withr::defer(suppressWarnings(gtrack.rm("test_override_err_target", force = TRUE)))

    assign(".create_dir_override", parent, envir = .misha)

    # mkdir errors because 'parent' already exists; the slot must
    # already be cleared by the time the error propagates.
    expect_error(
        gtrack.create("test_override_err_target", "test", "test.fixedbin"),
        regexp = "Cannot create track|exists"
    )

    override_after <- get0(".create_dir_override", envir = .misha, inherits = FALSE)
    expect_true(is.null(override_after))
})

test_that(".misha$.create_dir_override = NA_character_ behaves as no override", {
    # NA_character_ has length 1 but the C++ explicitly guards against
    # NA_STRING (see src/rdbutils.cpp::create_track_dir). The override
    # should be ignored and creation should succeed at the default path.
    # This is a regression guard: if someone removes the NA_STRING check
    # from the C++ the override would resolve to literal path "NA".
    suppressWarnings(gtrack.rm("test_override_na", force = TRUE))
    withr::defer(suppressWarnings(gtrack.rm("test_override_na", force = TRUE)))
    withr::defer(suppressWarnings(unlink("NA", recursive = TRUE, force = TRUE)))

    assign(".create_dir_override", NA_character_, envir = .misha)

    gtrack.create("test_override_na", "test", "test.fixedbin")
    expect_true(gtrack.exists("test_override_na"))

    # A literal "NA" directory must NOT have been created (in cwd or GWD).
    expect_false(dir.exists("NA"))
    expect_false(dir.exists(file.path(get("GWD", envir = .misha), "NA")))
})
