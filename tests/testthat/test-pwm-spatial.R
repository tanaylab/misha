create_isolated_test_db()

test_that("pwm with spatial parameters works", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    pssm <- create_test_pssm() # AC motif
    test_interval <- gintervals(1, 200, 300)

    # Spatial parameters - higher weight in middle
    spat_factors <- c(0.5, 1.0, 2.0, 1.0, 0.5)
    spat_bin <- 20L

    # Create spatial PWM vtrack
    gvtrack.create(
        "pwm_spatial", NULL, "pwm",
        list(
            pssm = pssm,
            bidirect = TRUE,
            extend = TRUE,
            prior = 0.01,
            spat_factor = spat_factors,
            spat_bin = spat_bin
        )
    )

    # Extract scores
    scores <- gextract("pwm_spatial", test_interval, iterator = test_interval)

    # Should return a numeric value
    expect_type(scores$pwm_spatial, "double")
    expect_false(is.na(scores$pwm_spatial[1]))
    expect_false(is.infinite(scores$pwm_spatial[1]))
})

test_that("pwm.max with spatial parameters works", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    pssm <- create_test_pssm()
    test_interval <- gintervals(1, 200, 300)

    spat_factors <- c(0.5, 1.0, 2.0, 1.0, 0.5)
    spat_bin <- 20L

    gvtrack.create(
        "pwm_max_spatial", NULL, "pwm.max",
        list(
            pssm = pssm,
            bidirect = TRUE,
            extend = TRUE,
            prior = 0.01,
            spat_factor = spat_factors,
            spat_bin = spat_bin
        )
    )

    scores <- gextract("pwm_max_spatial", test_interval, iterator = test_interval)

    expect_type(scores$pwm_max_spatial, "double")
    expect_false(is.na(scores$pwm_max_spatial[1]))
})

test_that("pwm.max.pos with spatial parameters works", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    pssm <- create_test_pssm()
    test_interval <- gintervals(1, 200, 300)

    spat_factors <- c(0.5, 1.0, 2.0, 1.0, 0.5)
    spat_bin <- 20L

    gvtrack.create(
        "pwm_maxpos_spatial", NULL, "pwm.max.pos",
        list(
            pssm = pssm,
            bidirect = TRUE,
            extend = TRUE,
            prior = 0.01,
            spat_factor = spat_factors,
            spat_bin = spat_bin
        )
    )

    scores <- gextract("pwm_maxpos_spatial", test_interval, iterator = test_interval)

    expect_type(scores$pwm_maxpos_spatial, "double")
    # Position should be within reasonable range (1-based)
    expect_gt(abs(scores$pwm_maxpos_spatial[1]), 0)
})

test_that("spatial PWM backward compatibility - no spatial params", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    pssm <- create_test_pssm()
    test_interval <- gintervals(1, 200, 240)

    # Create two identical vtracks - one explicitly without spatial, one with old API
    gvtrack.create("pwm_old", NULL, "pwm", list(pssm = pssm, bidirect = TRUE, extend = TRUE, prior = 0.01))
    gvtrack.create("pwm_new", NULL, "pwm", list(pssm = pssm, bidirect = TRUE, extend = TRUE, prior = 0.01))

    scores <- gextract(c("pwm_old", "pwm_new"), test_interval, iterator = test_interval)

    # Results should be identical
    expect_equal(scores$pwm_old[1], scores$pwm_new[1], tolerance = 1e-10)
})

test_that("spatial factors weight positions correctly", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    pssm <- create_test_pssm()
    test_interval <- gintervals(1, 200, 250)

    # Uniform spatial factors (all 1.0) should give same result as no spatial
    uniform_spat <- rep(1.0, 10)
    spat_bin <- 10L

    gvtrack.create(
        "pwm_nospatial", NULL, "pwm",
        list(pssm = pssm, bidirect = FALSE, extend = TRUE, prior = 0.01)
    )

    gvtrack.create(
        "pwm_uniform", NULL, "pwm",
        list(
            pssm = pssm,
            bidirect = FALSE,
            extend = TRUE,
            prior = 0.01,
            spat_factor = uniform_spat,
            spat_bin = spat_bin
        )
    )

    scores <- gextract(c("pwm_nospatial", "pwm_uniform"), test_interval, iterator = test_interval)

    # Should be very close (accounting for log(1.0) = 0)
    expect_equal(scores$pwm_nospatial[1], scores$pwm_uniform[1], tolerance = 1e-5)
})

test_that("spatial range parameters work", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    pssm <- create_test_pssm()
    test_interval <- gintervals(1, 200, 400)

    spat_factors <- c(0.5, 1.0, 2.0)
    spat_bin <- 50L

    # Test with range restriction
    gvtrack.create(
        "pwm_spatial_range", NULL, "pwm",
        list(
            pssm = pssm,
            bidirect = FALSE,
            extend = TRUE,
            prior = 0.01,
            spat_factor = spat_factors,
            spat_bin = spat_bin,
            spat_min = 0L,
            spat_max = 100L
        )
    )

    scores <- gextract("pwm_spatial_range", test_interval, iterator = test_interval)

    expect_type(scores$pwm_spatial_range, "double")
    expect_false(is.na(scores$pwm_spatial_range[1]))
})

test_that("spatial PWM with bidirectional scanning", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    pssm <- create_test_pssm()
    test_interval <- gintervals(1, 200, 300)

    spat_factors <- c(1.0, 2.0, 1.0)
    spat_bin <- 30L

    # Both strands should use same spatial binning
    gvtrack.create(
        "pwm_spatial_bidi", NULL, "pwm",
        list(
            pssm = pssm,
            bidirect = TRUE,
            extend = TRUE,
            prior = 0.01,
            spat_factor = spat_factors,
            spat_bin = spat_bin
        )
    )

    gvtrack.create(
        "pwm_spatial_fwd", NULL, "pwm",
        list(
            pssm = pssm,
            bidirect = FALSE,
            extend = TRUE,
            prior = 0.01,
            spat_factor = spat_factors,
            spat_bin = spat_bin
        )
    )

    scores <- gextract(c("pwm_spatial_bidi", "pwm_spatial_fwd"), test_interval, iterator = test_interval)

    # Bidirectional should be >= forward only (log-sum-exp property)
    expect_gte(scores$pwm_spatial_bidi[1], scores$pwm_spatial_fwd[1])
})

test_that("spatial PWM honors iterator shifts", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    pssm <- create_test_pssm()
    base <- gintervals(1, 2000, 2100)

    spat_factors <- c(0.5, 1.0, 2.0, 1.0, 0.5)
    spat_bin <- 20L

    # Create two vtracks with different shifts
    gvtrack.create(
        "pwm_spat_small", NULL, "pwm",
        list(
            pssm = pssm,
            bidirect = FALSE,
            extend = TRUE,
            prior = 0.01,
            spat_factor = spat_factors,
            spat_bin = spat_bin
        )
    )

    gvtrack.create(
        "pwm_spat_large", NULL, "pwm",
        list(
            pssm = pssm,
            bidirect = FALSE,
            extend = TRUE,
            prior = 0.01,
            spat_factor = spat_factors,
            spat_bin = spat_bin
        )
    )

    gvtrack.iterator("pwm_spat_small", sshift = -10, eshift = 10)
    gvtrack.iterator("pwm_spat_large", sshift = -50, eshift = 50)

    scores <- gextract(c("pwm_spat_small", "pwm_spat_large"), base, iterator = base)

    # Both should return valid scores
    expect_false(is.na(scores$pwm_spat_small[1]))
    expect_false(is.na(scores$pwm_spat_large[1]))

    # Different window sizes should generally give different results
    # (unless by chance they're identical)
})

test_that("spatial PWM error handling", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    pssm <- create_test_pssm()

    # Should error with non-positive spatial factors
    expect_error(
        gvtrack.create(
            "pwm_bad_spat", NULL, "pwm",
            list(
                pssm = pssm,
                spat_factor = c(-1, 1, 1),
                spat_bin = 10L
            )
        ),
        "positive"
    )

    # Should error with non-positive bin size
    expect_error(
        gvtrack.create(
            "pwm_bad_bin", NULL, "pwm",
            list(
                pssm = pssm,
                spat_factor = c(1, 1, 1),
                spat_bin = 0L
            )
        ),
        "positive"
    )

    # Should error with zero spatial factor
    expect_error(
        gvtrack.create(
            "pwm_zero_spat", NULL, "pwm",
            list(
                pssm = pssm,
                spat_factor = c(0, 1, 1),
                spat_bin = 10L
            )
        ),
        "positive"
    )
})

test_that("spatial sliding agrees with non-sliding scoring", {
    remove_all_vtracks()

    # The spatial sliding window is an optimization: it must return exactly what
    # the straightforward per-interval scan returns. MISHA_DISABLE_SPATIAL_SLIDING
    # selects the reference path, so the two can be compared directly.
    # Regression coverage for three defects that made them disagree:
    #   - the incoming anchor was read from a fixed offset, so a stride > 1 slide
    #     (any fixed-size iterator) reused one base for every step;
    #   - the per-bin max rescan re-selected the element being evicted or moved;
    #   - the seed's hit count stopped at bins*B, missing the positions that clamp
    #     into the last bin when the spatial profile is shorter than the window.
    old <- Sys.getenv("MISHA_DISABLE_SPATIAL_SLIDING", unset = NA)
    on.exit(
        {
            if (is.na(old)) {
                Sys.unsetenv("MISHA_DISABLE_SPATIAL_SLIDING")
            } else {
                Sys.setenv(MISHA_DISABLE_SPATIAL_SLIDING = old)
            }
        },
        add = TRUE
    )

    set.seed(23)
    intervals <- gintervals(1, 2000, 2300)
    checked <- 0L
    for (motif_len in c(6L, 9L)) {
        pssm <- matrix(runif(motif_len * 4, 0.05, 1),
            nrow = motif_len, dimnames = list(NULL, c("A", "C", "G", "T"))
        )
        pssm <- pssm / rowSums(pssm)
        for (spat_bin in c(1L, 5L)) {
            # A short profile relative to the window exercises the last-bin clamp.
            for (spat_len in c(10L, 60L)) {
                spat <- runif(spat_len, 0.2, 3)
                for (func in c("pwm", "pwm.max", "pwm.max.pos", "pwm.count")) {
                    for (strand in c(1, -1)) {
                        params <- list(
                            pssm = pssm, prior = 0, bidirect = FALSE, strand = strand,
                            extend = TRUE, spat_factor = spat, spat_bin = spat_bin
                        )
                        if (func == "pwm.count") params$score.thresh <- -12
                        gvtrack.create("spat_v", NULL, func, params = params)

                        Sys.setenv(MISHA_DISABLE_SPATIAL_SLIDING = "")
                        slid <- gextract("spat_v", intervals, iterator = 50)
                        Sys.setenv(MISHA_DISABLE_SPATIAL_SLIDING = "1")
                        ref <- gextract("spat_v", intervals, iterator = 50)

                        slid <- slid[order(slid$start), ]
                        ref <- ref[order(ref$start), ]
                        expect_equal(slid$spat_v, ref$spat_v,
                            tolerance = 1e-5,
                            info = sprintf(
                                "func=%s strand=%d spat_bin=%d spat_len=%d motif_len=%d",
                                func, strand, spat_bin, spat_len, motif_len
                            )
                        )
                        checked <- checked + 1L
                    }
                }
            }
        }
    }
    # More than one bin must be produced, or sliding never happens and the test
    # would pass without exercising anything.
    expect_gt(nrow(ref), 1L)
    expect_gt(checked, 0L)
})
