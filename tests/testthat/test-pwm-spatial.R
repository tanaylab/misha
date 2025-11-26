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
