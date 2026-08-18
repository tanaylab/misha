create_isolated_test_db()

test_that("pwm.count spatial sliding counts bidirectional hits once per position", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    pssm <- create_test_pssm()

    # score.thresh = -25 was below the lowest score this PSSM can reach (-9.2888, see the
    # note at the top of test-pwm-count.R), so both paths counted every position and both
    # returned the interval length, 55, for all 30 intervals: they agreed because neither
    # could return anything else. At -5 they return 41-44, which varies across the intervals
    # and leaves room to differ in both directions.
    params <- list(
        pssm = pssm,
        bidirect = TRUE,
        strand = 1,
        extend = TRUE,
        prior = 0.01,
        score.thresh = -5,
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

    # Two paths agreeing on a constant they could not have missed proves nothing; assert the
    # counts actually vary before comparing them.
    expect_gt(length(unique(res$pwm_count_spat_slide)), 1)
    expect_lt(max(res$pwm_count_spat_slide), 55)
    expect_equal(res$pwm_count_spat_slide, res$pwm_count_spat_ref, tolerance = 1e-6)
})
