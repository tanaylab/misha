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

test_that("pwm.max.pos with spat_factor aggregates by score across a filter's parts", {
    # The spatial scorer answers pwm.max.pos from its own code path
    # (spat_answer_MAXPOS), separately from the non-spatial one, and each has
    # to publish the part's max score for the aggregation to compare parts by
    # score at all. Without spat_factor this file cannot reach that path, and
    # the geometry of the tests above answers its two parts through different
    # paths in an order that happens to give the right answer even when the
    # spatial one publishes nothing - so it takes a separate geometry, with a
    # winner the spatial path answers, to constrain that side.
    #
    # Verified to bite: with the max-score publication removed from
    # spat_answer_MAXPOS, this returns 236 - part B's local argmax of 36
    # offset by 200 - where 78, part A's, is correct.
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    pssm <- filter_agg_pssm()
    iv <- gintervals(1, 3000, 3300)
    mask <- gintervals(1, 3100, 3200)
    parts <- list(gintervals(1, 3000, 3100), gintervals(1, 3200, 3300))

    # unit weights, so the spatial factors leave the scores alone and the
    # reference can use the same params - only the code path differs
    params <- list(
        pssm = pssm, bidirect = FALSE, strand = 1,
        extend = TRUE, prior = 0.01,
        spat_factor = rep(1.0, 5), spat_bin = 20L
    )

    gvtrack.create("t_pos", NULL, "pwm.max.pos", params = params)
    gvtrack.filter("t_pos", filter = mask)
    gvtrack.create("p_pos", NULL, "pwm.max.pos", params = params)
    gvtrack.create("p_max", NULL, "pwm.max", params = params)

    got <- gextract("t_pos", iv, iterator = iv)$t_pos
    part_max <- vapply(parts, function(p) gextract("p_max", p, iterator = p)$p_max, numeric(1))
    part_pos <- vapply(parts, function(p) gextract("p_pos", p, iterator = p)$p_pos, numeric(1))

    # the winner must be part A, and by a clear margin: if part B won, or the
    # two tied, an aggregation that silently falls back to "first part wins"
    # would agree with the reference and this test would prove nothing
    expect_gt(abs(diff(part_max)), 1e-3)
    win <- which.max(part_max)
    expect_equal(win, 1L)

    expected <- part_pos[win] + (parts[[win]]$start - iv$start)
    expect_equal(got, expected, info = "spatial pwm.max.pos across two parts")

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

test_that("a pwm part the scorer cannot score is left out of the aggregation", {
    # Not an assembly gap, although that is the case the potts side cares
    # about: pwm charges an N the mean of its PSSM column, so an all-N part
    # scores -Inf (0 for pwm.count), and -Inf is the identity of both a
    # log-sum-exp and a maximum - it reduces correctly on its own. The part a
    # pwm scorer genuinely cannot score is one with no anchor to place at all:
    # narrower than the PSSM with extend = FALSE, where score_interval()
    # returns NaN.
    #
    # The mask leaves that part FIRST, which is the ordering that exposes it:
    # pwm.max seeded itself from the first part and could never recover once
    # that seed was NaN. pwm.count is order-independent - a NaN poisons its sum
    # from anywhere - and pwm.max.pos already skipped.
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    pssm <- filter_agg_pssm() # 6 rows
    iv <- gintervals(1, 4000, 4300)
    mask <- gintervals(1, 4003, 4200) # leaves a 3bp first part, < nrow(pssm)
    part_a <- gintervals(1, 4000, 4003)
    part_b <- gintervals(1, 4200, 4300)

    params <- list(
        pssm = pssm, bidirect = FALSE, strand = 1,
        extend = FALSE, prior = 0.01, score.thresh = -12
    )

    for (fn in c("pwm", "pwm.max", "pwm.count", "pwm.max.pos")) {
        remove_all_vtracks()
        gvtrack.create("t", NULL, fn, params = params)
        gvtrack.filter("t", filter = mask)
        gvtrack.create("t_part", NULL, fn, params = params)

        # the geometry has to actually produce one unscorable part and one
        # scorable one, or this test asserts nothing
        expect_true(is.na(gextract("t_part", part_a, iterator = part_a)$t_part), info = fn)
        b <- gextract("t_part", part_b, iterator = part_b)$t_part
        expect_false(is.na(b), info = fn)

        got <- gextract("t", iv, iterator = iv)$t
        expected <- if (fn == "pwm.max.pos") b + (part_b$start - iv$start) else b
        expect_equal(got, expected, tolerance = 1e-5, info = fn)
    }
})

test_that("a filtered pwm.count separates 'counted nothing' from 'nothing to count'", {
    # The same two answers score_potts_var() has to keep apart, and the same
    # rule:
    #
    #   every part narrower than the PSSM, so no anchor can be placed in any of
    #     them - nothing to count, NA. A 0 would be indistinguishable in the
    #     returned column from a bin that really was scanned and had no hits,
    #     which folds "excluded" into "no evidence" under any sum(), mean() or
    #     threshold downstream and hides the excluded bins from is.na().
    #   parts that hold anchors none of which is a match - an assembly gap -
    #     really did count, and found none, so 0.
    #
    # The second block is what stops this from being satisfied by "NA
    # everywhere".
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    pssm <- filter_agg_pssm() # 6 rows
    iv <- gintervals(1, 4000, 4300)
    narrow_mask <- gintervals(1, 4003, 4297) # leaves 3bp at each end

    gap <- gintervals(20, 0, 300)
    skip_if_not(
        grepl("^N+$", toupper(gseq.extract(gintervals(20, 0, 300L + nrow(pssm) - 1L)))),
        "no all-N interval at the start of chr20 in this fixture"
    )
    gap_mask <- gintervals(20, 100, 200)

    for (fn in c("pwm", "pwm.max", "pwm.max.pos", "pwm.count")) {
        remove_all_vtracks()
        gvtrack.create("t", NULL, fn, params = list(
            pssm = pssm, bidirect = FALSE, strand = 1,
            extend = FALSE, prior = 0.01, score.thresh = -12
        ))
        gvtrack.filter("t", filter = narrow_mask)
        expect_true(is.na(gextract("t", iv, iterator = iv)$t),
            info = paste("no anchor in any part", fn)
        )

        # An N part is scorable and merely hopeless - a PSSM charges an N the
        # mean of its column - so this is a real count of zero and must stay 0.
        remove_all_vtracks()
        gvtrack.create("t", NULL, fn, params = list(
            pssm = pssm, bidirect = FALSE, strand = 1,
            extend = TRUE, prior = 0.01, score.thresh = -12
        ))
        gvtrack.filter("t", filter = gap_mask)
        got <- gextract("t", gap, iterator = gap)$t
        if (fn == "pwm.count") {
            expect_equal(got, 0, info = paste("assembly gap", fn))
        } else {
            expect_false(is.na(got), info = paste("assembly gap", fn))
        }
    }
})
