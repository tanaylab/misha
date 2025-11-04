test_that("gtrack.create works with PWM virtual tracks on indexed databases", {
    # This is a regression test for a bug where gtrack.create would fail with
    # "Invalid format of intervals argument" when creating tracks from PWM virtual
    # tracks on indexed databases with many chromosomes.
    #
    # The bug occurred because empty intervals were passed during internal validation,
    # and the code didn't handle empty intervals gracefully.

    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    # Create an indexed database with multiple chromosomes
    test_fasta <- tempfile(fileext = ".fasta")
    cat(">chr1\nACTGACTGACTGACTG\n>chr2\nGGCCGGCCGGCCGGCC\n>chr3\nTATATATATATATATA\n",
        file = test_fasta
    )

    test_db <- tempfile()
    withr::defer({
        unlink(test_db, recursive = TRUE)
        unlink(test_fasta)
    })

    # Create indexed database
    withr::with_options(list(gmulticontig.indexed_format = TRUE), {
        gdb.create(groot = test_db, fasta = test_fasta, verbose = FALSE)
        gdb.init(test_db)

        # Create a simple PWM
        pssm <- create_test_pssm() # AC motif

        # Create PWM virtual track
        gvtrack.create("pwm_energy", NULL,
            func = "pwm",
            pssm = pssm,
            bidirect = TRUE,
            prior = 0.01
        )

        # Create a directory for the track
        gdir.create("test_tracks", showWarnings = FALSE)
        withr::defer(gdir.rm("test_tracks", force = TRUE, recursive = TRUE))

        # This should NOT fail with "Invalid format of intervals argument"
        # The bug occurred specifically when using gtrack.create with PWM vtracks
        # on indexed databases
        expect_no_error({
            gtrack.create(
                track = "test_tracks.pwm_result",
                description = "PWM track from vtrack",
                expr = "pwm_energy",
                iterator = gintervals("chr1", 0, 10)
            )
        })

        # Verify the track was created successfully
        expect_true(gtrack.exists("test_tracks.pwm_result"))

        # Verify the track info
        info <- gtrack.info("test_tracks.pwm_result")
        expect_equal(info$type, "sparse")
        expect_equal(info$dimensions, 1)
        expect_equal(info$format, "indexed")

        # Verify we can extract from the created track
        result <- gextract("test_tracks.pwm_result", gintervals("chr1", 0, 10))
        expect_true(nrow(result) > 0)
        expect_true("test_tracks.pwm_result" %in% names(result))
        expect_true(is.numeric(result$test_tracks.pwm_result))

        # Verify the values match the original virtual track
        vtrack_result <- gextract("pwm_energy", gintervals("chr1", 0, 10),
            iterator = gintervals("chr1", 0, 10)
        )

        # The track should have the same value as the virtual track for this interval
        expect_equal(result$test_tracks.pwm_result[1], vtrack_result$pwm_energy[1],
            tolerance = 1e-6
        )
    })
})

test_that("gtrack.create works with PWM virtual tracks on regular databases", {
    # Verify that the fix doesn't break regular (non-indexed) databases

    # Re-initialize test database after the indexed test
    gdb.init_examples()

    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    pssm <- create_test_pssm() # AC motif

    gvtrack.create("pwm_energy", NULL,
        func = "pwm",
        pssm = pssm,
        bidirect = TRUE,
        prior = 0.01
    )

    gdir.create("test_tracks", showWarnings = FALSE)
    withr::defer(gdir.rm("test_tracks", force = TRUE, recursive = TRUE))

    # This should work on regular databases too
    expect_no_error({
        gtrack.create(
            track = "test_tracks.pwm_regular",
            description = "PWM track from vtrack on regular DB",
            expr = "pwm_energy",
            iterator = gintervals(1, 200, 210)
        )
    })

    # Verify the track was created
    expect_true(gtrack.exists("test_tracks.pwm_regular"))

    # Verify we can extract from it
    result <- gextract("test_tracks.pwm_regular", gintervals(1, 200, 210))
    expect_true(nrow(result) > 0)
    expect_true(is.numeric(result$test_tracks.pwm_regular))
})
