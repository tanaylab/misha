create_isolated_test_db()

# A gvtrack.filter on a pwm vtrack scores the filter's unmasked parts
# separately and reduces them. Every reference below is built from the same
# func and params with NO filter, evaluated over each unmasked part on its own,
# so the reduction is what is under test and not the per-part scoring.

filter_agg_pssm <- function() {
    pssm <- matrix(c(
        0.90, 0.04, 0.03, 0.03,
        0.03, 0.85, 0.06, 0.06,
        0.05, 0.05, 0.80, 0.10,
        0.08, 0.08, 0.09, 0.75,
        0.70, 0.10, 0.10, 0.10,
        0.12, 0.65, 0.12, 0.11
    ), ncol = 4, byrow = TRUE)
    colnames(pssm) <- c("A", "C", "G", "T")
    pssm
}

test_that("pwm vtracks aggregate across a filter's unmasked parts", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    pssm <- filter_agg_pssm()
    iv <- gintervals(1, 4000, 4300)
    mask <- gintervals(1, 4100, 4200) # punches the middle out, leaving two parts
    part_a <- gintervals(1, 4000, 4100)
    part_b <- gintervals(1, 4200, 4300)
    th <- -12

    for (fn in c("pwm", "pwm.max", "pwm.count")) {
        remove_all_vtracks()
        params <- list(
            pssm = pssm, bidirect = FALSE, strand = 1,
            extend = TRUE, prior = 0.01, score.thresh = th
        )
        gvtrack.create("t", NULL, fn, params = params)
        gvtrack.filter("t", filter = mask)
        gvtrack.create("t_part", NULL, fn, params = params)

        got <- gextract("t", iv, iterator = iv)$t
        parts <- c(
            gextract("t_part", part_a, iterator = part_a)$t_part,
            gextract("t_part", part_b, iterator = part_b)$t_part
        )

        expected <- switch(fn,
            # log-sum-exp, not a sum: summing log-likelihoods multiplies
            # probabilities, which is not "the score over this interval"
            pwm = log_sum_exp(parts),
            pwm.max = max(parts),
            pwm.count = sum(parts)
        )
        # the two parts must differ, or a max that picked the wrong one would
        # still agree, and a sum would be within tolerance of a log-sum-exp
        expect_gt(abs(diff(parts)), 1e-3)
        expect_equal(got, expected, tolerance = 1e-5, info = fn)
    }
})

test_that("pwm.max.pos aggregates across a filter's unmasked parts by score, not index", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    pssm <- filter_agg_pssm()
    iv <- gintervals(1, 4000, 4300)
    mask <- gintervals(1, 4100, 4200)
    parts <- list(gintervals(1, 4000, 4100), gintervals(1, 4200, 4300))
    params <- list(
        pssm = pssm, bidirect = FALSE, strand = 1,
        extend = TRUE, prior = 0.01
    )

    gvtrack.create("t_pos", NULL, "pwm.max.pos", params = params)
    gvtrack.filter("t_pos", filter = mask)
    gvtrack.create("p_pos", NULL, "pwm.max.pos", params = params)
    gvtrack.create("p_max", NULL, "pwm.max", params = params)

    got <- gextract("t_pos", iv, iterator = iv)$t_pos

    part_max <- vapply(parts, function(p) gextract("p_max", p, iterator = p)$p_max, numeric(1))
    part_pos <- vapply(parts, function(p) gextract("p_pos", p, iterator = p)$p_pos, numeric(1))

    # unambiguous winner, so the expected part is not a coin toss
    expect_gt(abs(diff(part_max)), 1e-3)
    win <- which.max(part_max)
    expected <- part_pos[win] + (parts[[win]]$start - iv$start)

    expect_equal(got, expected, info = "position, offset into the iterator interval")

    # The released code compared the parts by their part-RELATIVE position and
    # reported the larger one with no offset, so it named an anchor in whichever
    # part happened to carry the bigger local index, in that part's coordinate
    # frame rather than the interval's. On this geometry it returned 52 where
    # 216 is right, so the expect_equal above is the assertion that catches it.
    # The three below are range invariants that hold either way; they are here
    # to pin the output's meaning, not to detect this defect.
    expect_gt(got, 0)
    expect_lte(got, iv$end - iv$start)

    genomic_anchor_start <- iv$start + got - 1L
    expect_true(genomic_anchor_start < mask$start || genomic_anchor_start >= mask$end,
        info = paste("anchor at", genomic_anchor_start)
    )
})

test_that("a single unmasked pwm part is scored like the unfiltered survivor, offset when it is a position", {
    # A mask flush with one edge of the iterator interval leaves exactly one
    # unmasked part rather than two. pwm, pwm.max and pwm.count don't care: a
    # lone part is just scored directly. pwm.max.pos does: a mask that clips
    # the START moves the surviving part's own start away from the interval's
    # start, so the position the part reports relative to itself has to be
    # offset before it means anything relative to the original interval. A
    # mask clipping the END leaves offset zero, so that side alone could never
    # catch a missing offset.
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    pssm <- filter_agg_pssm()
    iv <- gintervals(1, 4000, 4300)
    th <- -12

    # The start clip is deliberately WIDER than the part it leaves (200bp of
    # mask, 100bp of part). A part-relative position is then always in 1..100,
    # so reporting one without its offset necessarily decodes to a coordinate
    # inside the mask - which makes the "outside the mask" assertion at the
    # bottom of the loop a second, independent detector of a missing offset,
    # not just an invariant that happens to hold.
    clips <- list(
        start = list(mask = gintervals(1, 4000, 4200), part = gintervals(1, 4200, 4300)),
        end = list(mask = gintervals(1, 4200, 4300), part = gintervals(1, 4000, 4200))
    )

    for (nm in names(clips)) {
        mask <- clips[[nm]]$mask
        part <- clips[[nm]]$part

        for (fn in c("pwm", "pwm.max", "pwm.count")) {
            remove_all_vtracks()
            params <- list(
                pssm = pssm, bidirect = FALSE, strand = 1,
                extend = TRUE, prior = 0.01, score.thresh = th
            )
            gvtrack.create("t", NULL, fn, params = params)
            gvtrack.filter("t", filter = mask)
            gvtrack.create("t_part", NULL, fn, params = params)

            got <- gextract("t", iv, iterator = iv)$t
            expected <- gextract("t_part", part, iterator = part)$t_part
            expect_equal(got, expected, tolerance = 1e-5, info = paste(nm, fn))
        }

        remove_all_vtracks()
        params <- list(pssm = pssm, bidirect = FALSE, strand = 1, extend = TRUE, prior = 0.01)
        gvtrack.create("t_pos", NULL, "pwm.max.pos", params = params)
        gvtrack.filter("t_pos", filter = mask)
        gvtrack.create("p_pos", NULL, "pwm.max.pos", params = params)

        got <- gextract("t_pos", iv, iterator = iv)$t_pos
        local_pos <- gextract("p_pos", part, iterator = part)$p_pos
        expected <- local_pos + (part$start - iv$start)
        expect_equal(got, expected, info = paste(nm, "pwm.max.pos value"))

        genomic_anchor_start <- iv$start + got - 1L
        expect_true(genomic_anchor_start < mask$start || genomic_anchor_start >= mask$end,
            info = paste(nm, "anchor at", genomic_anchor_start, "mask", mask$start, mask$end)
        )
    }
})

test_that("a fully masked pwm interval reports NA", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    pssm <- filter_agg_pssm()
    iv <- gintervals(1, 4000, 4300)

    for (fn in c("pwm", "pwm.max", "pwm.max.pos", "pwm.count")) {
        remove_all_vtracks()
        gvtrack.create("t", NULL, fn, params = list(
            pssm = pssm, bidirect = FALSE, strand = 1,
            extend = TRUE, prior = 0.01, score.thresh = -12
        ))
        gvtrack.filter("t", filter = iv)
        expect_true(is.na(gextract("t", iv, iterator = iv)$t), info = fn)
    }
})
