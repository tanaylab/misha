# Regressions found by turning on evaluation of the package vignettes.
# Each of these was a live defect that the vignettes documented or exercised
# but never actually ran.

restore_groot_on_exit()

test_that("gintervals.save() with the arguments swapped names the mistake", {
    gdb.init_examples()
    intervs <- gintervals(1, c(0, 250000), c(100000, 260000))
    expect_error(
        gintervals.save(intervs, "my_intervals_set"),
        "expected a single character string"
    )
    # And the right way round still works.
    gintervals.save("my_intervals_set", intervs)
    expect_true(gintervals.exists("my_intervals_set"))
    expect_equal(nrow(gintervals.load("my_intervals_set")), 2)
})

test_that("gtrack.dataset() keeps working after a track is created", {
    gdb.init_examples()
    # .gdb.add_track() used to rebuild GTRACK_DATASET from scratch, so the
    # first gtrack.create() of a session dropped every pre-existing track
    # from the map and gtrack.dataset() started returning NA for them.
    before <- gtrack.dataset("dense_track")
    expect_false(is.na(before))

    gtrack.create("vign_newtrack", "regression", "dense_track * 2")
    on.exit(gtrack.rm("vign_newtrack", force = TRUE), add = TRUE)

    expect_equal(gtrack.dataset("dense_track"), before)
    expect_equal(unname(gtrack.dataset("vign_newtrack")), before)
})

test_that("gintervals.dataset() keeps working after an intervals set is saved", {
    gdb.init_examples()
    before <- gintervals.dataset("annotations")
    expect_false(is.na(before))

    gintervals.save("vign_peaks", gintervals(1, 0, 1000))
    expect_equal(gintervals.dataset("annotations"), before)
    expect_equal(unname(gintervals.dataset("vign_peaks")), before)
})

test_that("gdataset.save() handles namespaced track and interval names", {
    gdb.init_examples()
    gdir.create("vignsub", showWarnings = FALSE)
    gintervals.save("vignsub.myintervs", gintervals(1, 0, 1000))

    ds <- file.path(tempdir(), "vign_dataset")
    unlink(ds, recursive = TRUE)
    on.exit(unlink(ds, recursive = TRUE), add = TRUE)

    # "subdir.dense_track2" lives at tracks/subdir/dense_track2.track. The
    # dotted name used to be pasted straight into the path, which produced an
    # empty tracks/subdir.dense_track2.track directory (plus a warning) and a
    # dataset whose namespaced track could not be read back.
    expect_silent(
        gdataset.save(
            path = ds, description = "namespaced",
            tracks = c("dense_track", "subdir.dense_track2"),
            intervals = c("annotations", "vignsub.myintervs")
        )
    )
    expect_true(dir.exists(file.path(ds, "tracks", "subdir", "dense_track2.track")))
    expect_true(file.exists(file.path(ds, "tracks", "vignsub", "myintervs.interv")))
    expect_false(dir.exists(file.path(ds, "tracks", "subdir.dense_track2.track")))

    # And the dataset is actually readable.
    linked <- file.path(tempdir(), "vign_linked_db")
    unlink(linked, recursive = TRUE)
    on.exit(unlink(linked, recursive = TRUE), add = TRUE)
    suppressMessages(gdb.create_linked(linked, parent = .misha$GROOT))
    gsetroot(linked)
    gdataset.load(ds)
    expect_true("subdir.dense_track2" %in% gtrack.ls())
    expect_equal(
        nrow(gextract("subdir.dense_track2", gintervals(1, 0, 300))),
        6
    )
    gdataset.unload(ds)
})
