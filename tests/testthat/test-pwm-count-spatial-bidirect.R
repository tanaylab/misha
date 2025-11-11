test_that("pwm.count spatial sliding counts bidirectional hits once per position", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    pssm <- create_test_pssm()

    params <- list(
        pssm = pssm,
        bidirect = TRUE,
        strand = 1,
        extend = TRUE,
        prior = 0.01,
        score.thresh = -25,
        spat_factor = rep(1.0, 6),
        spat_bin = 15L
    )

    gvtrack.create("pwm_count_spat_slide", NULL, "pwm.count", params)

    withr::with_envvar(c(MISHA_DISABLE_SPATIAL_SLIDING = "1"), {
        gvtrack.create("pwm_count_spat_ref", NULL, "pwm.count", params)
    })

    n <- 30
    starts <- 2400 + 0:(n - 1)
    ends <- starts + 55L
    ivs <- gintervals(rep(1L, n), starts, ends)

    res <- gextract(c("pwm_count_spat_slide", "pwm_count_spat_ref"),
        ivs,
        iterator = ivs
    )

    expect_equal(res$pwm_count_spat_slide, res$pwm_count_spat_ref, tolerance = 1e-6)
})
