create_isolated_test_db()

# gseq.pwm() and the pwm virtual tracks are two entry points into the same
# scorer, so under bidirect = TRUE they must give the same answer on the same
# anchors. These tests assert exactly that, mode by mode, rather than checking
# either side against a formula rewritten in R.
#
# A graded PSSM, not a consensus one: every column keeps probability on all
# four bases, so a window's forward and reverse-complement scores are both
# finite and different, and the two strand unions (log-sum-exp vs maximum) are
# genuinely distinguishable. A deterministic matrix would push one strand to
# the prior floor and hide the difference.
parity_pssm <- function() {
    m <- matrix(c(
        0.70, 0.10, 0.15, 0.05,
        0.05, 0.60, 0.25, 0.10,
        0.10, 0.15, 0.65, 0.10,
        0.20, 0.10, 0.10, 0.60,
        0.55, 0.20, 0.15, 0.10
    ), ncol = 4, byrow = TRUE)
    colnames(m) <- c("A", "C", "G", "T")
    m
}

parity_intervals <- function() {
    gintervals(
        c(1, 1, 1, 1, 2),
        c(200, 1000, 5000, 100000, 50000),
        c(260, 1080, 5040, 100200, 50200)
    )
}

# The vtrack scores in float, gseq.pwm in double, so the two agree to float
# precision and not further.
parity_tolerance <- 1e-5

test_that("gseq.pwm(mode = 'lse') equals the pwm vtrack under bidirect = TRUE", {
    remove_all_vtracks()
    pssm <- parity_pssm()
    ivs <- parity_intervals()

    gvtrack.create("p_lse", NULL,
        func = "pwm", pssm = pssm,
        bidirect = TRUE, extend = FALSE, prior = 0.01
    )
    v <- gextract("p_lse", ivs, iterator = ivs)
    v <- v[order(v$intervalID), ]

    g <- gseq.pwm(toupper(gseq.extract(ivs)), pssm,
        mode = "lse", bidirect = TRUE, extend = FALSE, prior = 0.01
    )

    expect_equal(g, v$p_lse, tolerance = parity_tolerance, ignore_attr = TRUE)
})

test_that("gseq.pwm(mode = 'max') equals the pwm.max vtrack under bidirect = TRUE", {
    remove_all_vtracks()
    pssm <- parity_pssm()
    ivs <- parity_intervals()

    gvtrack.create("p_max", NULL,
        func = "pwm.max", pssm = pssm,
        bidirect = TRUE, extend = FALSE, prior = 0.01
    )
    v <- gextract("p_max", ivs, iterator = ivs)
    v <- v[order(v$intervalID), ]

    g <- gseq.pwm(toupper(gseq.extract(ivs)), pssm,
        mode = "max", bidirect = TRUE, extend = FALSE, prior = 0.01
    )

    expect_equal(g, v$p_max, tolerance = parity_tolerance, ignore_attr = TRUE)
})

test_that("gseq.pwm equals the pwm vtracks with extend = TRUE", {
    remove_all_vtracks()
    pssm <- parity_pssm()
    w <- nrow(pssm)
    ivs <- parity_intervals()

    # extend = TRUE extends the END of the interval by w - 1 (see
    # GenomeSeqScorer::calculate_expanded_interval), so the equivalent
    # sequence-level call is the end-extended sequence scanned without
    # extension. That adds w - 1 = 4 anchors per interval, 20 in all.
    ext <- ivs
    ext$end <- ext$end + w - 1L
    seqs_ext <- toupper(gseq.extract(ext))

    gvtrack.create("p_max", NULL,
        func = "pwm.max", pssm = pssm,
        bidirect = TRUE, extend = TRUE, prior = 0.01
    )
    # mode = "max" cannot see the extension on these intervals: none of the 20
    # extra anchors beats the best one already inside, so the end-extended and
    # the plain sequence give the same five values and dropping the ext line
    # above would leave a max-only assertion passing. mode = "count" does see
    # it - 19 14 8 42 49 extended against 17 14 7 41 47 plain - so the count
    # line below is what holds the extend geometry, and the max line is here
    # for the strand union at extend = TRUE.
    #
    # The threshold is the one the extend = FALSE block justifies, and the 20
    # extra anchors do not disturb its gap: pooled over all 580 anchors the
    # nearest per-anchor union values are still -6.410838 and -6.350372 and the
    # nearest single-strand scores still -6.410507 and -6.361461.
    gvtrack.create("p_count", NULL,
        func = "pwm.count", pssm = pssm, score.thresh = -6.3776,
        bidirect = TRUE, extend = TRUE, prior = 0.01
    )
    v <- gextract(c("p_max", "p_count"), ivs, iterator = ivs)
    v <- v[order(v$intervalID), ]

    expect_equal(
        gseq.pwm(seqs_ext, pssm,
            mode = "max", bidirect = TRUE, extend = FALSE, prior = 0.01
        ),
        v$p_max,
        tolerance = parity_tolerance, ignore_attr = TRUE
    )
    expect_equal(
        gseq.pwm(seqs_ext, pssm,
            mode = "count", score.thresh = -6.3776,
            bidirect = TRUE, extend = FALSE, prior = 0.01
        ),
        v$p_count,
        tolerance = parity_tolerance, ignore_attr = TRUE
    )

    # And the extension is load-bearing: the plain interval's sequence is four
    # anchors short per interval and gives a different count.
    expect_false(isTRUE(all.equal(
        gseq.pwm(toupper(gseq.extract(ivs)), pssm,
            mode = "count", score.thresh = -6.3776,
            bidirect = TRUE, extend = FALSE, prior = 0.01
        ),
        v$p_count
    )))
})

test_that("gseq.pwm(mode = 'count') equals the pwm.count vtrack under bidirect = TRUE", {
    remove_all_vtracks()
    pssm <- parity_pssm()
    ivs <- parity_intervals()
    seqs <- toupper(gseq.extract(ivs))

    # -6.3776 sits in a gap of the score distribution of these intervals.
    # Pooled over all 560 anchors of parity_intervals(), the nearest per-anchor
    # log-sum-exp values are -6.410838 and -6.350372, and the nearest
    # single-strand scores -6.410507 and -6.361461, so neither side of the
    # comparison is deciding a near-tie. A threshold read off quantile() would
    # land on an exact tie instead - a window and its reverse complement
    # produce the same pair of numbers, so ties are common here - and then the
    # vtrack's float against gseq.pwm's double would decide the result.
    thresh <- -6.3776

    gvtrack.create("p_count", NULL,
        func = "pwm.count", pssm = pssm, score.thresh = thresh,
        bidirect = TRUE, extend = FALSE, prior = 0.01
    )
    v <- gextract("p_count", ivs, iterator = ivs)
    v <- v[order(v$intervalID), ]

    g <- gseq.pwm(seqs, pssm,
        mode = "count", score.thresh = thresh,
        bidirect = TRUE, extend = FALSE, prior = 0.01
    )

    expect_equal(g, v$p_count, tolerance = parity_tolerance, ignore_attr = TRUE)

    # Not saturated at either end, so the comparison above has something to say.
    n_anchors <- nchar(seqs) - nrow(pssm) + 1L
    expect_true(all(g > 0 & g < n_anchors))
})

test_that("a count is one hit per anchor, not one per strand", {
    remove_all_vtracks()
    pssm <- parity_pssm()
    ivs <- parity_intervals()
    seqs <- toupper(gseq.extract(ivs))
    n_anchors <- nchar(seqs) - nrow(pssm) + 1L

    # A threshold below every attainable score: with prior = 0.01 no score is
    # -Inf, so every anchor passes on both strands. One count per anchor is
    # therefore exactly the number of anchors; one count per strand would be
    # twice that.
    thresh <- -50

    gvtrack.create("p_count", NULL,
        func = "pwm.count", pssm = pssm, score.thresh = thresh,
        bidirect = TRUE, extend = FALSE, prior = 0.01
    )
    v <- gextract("p_count", ivs, iterator = ivs)
    v <- v[order(v$intervalID), ]

    g <- gseq.pwm(seqs, pssm,
        mode = "count", score.thresh = thresh,
        bidirect = TRUE, extend = FALSE, prior = 0.01
    )

    expect_equal(g, as.numeric(n_anchors), ignore_attr = TRUE)
    expect_equal(g, v$p_count, tolerance = parity_tolerance, ignore_attr = TRUE)
})

test_that("both scan paths take the same strand union as the vtracks", {
    remove_all_vtracks()
    pssm <- parity_pssm()
    ivs <- parity_intervals()
    seqs <- toupper(gseq.extract(ivs))
    thresh <- -6.3776

    # score_pwm_over_range() carries the mode branches twice: skip_gaps = TRUE
    # (the default, and what the tests above exercise) goes through the gap
    # projector, skip_gaps = FALSE through the contiguous scan. Both copies have
    # to agree with the vtrack, so both are pinned against it here.
    gvtrack.create("p_max", NULL,
        func = "pwm.max", pssm = pssm,
        bidirect = TRUE, extend = FALSE, prior = 0.01
    )
    gvtrack.create("p_count", NULL,
        func = "pwm.count", pssm = pssm, score.thresh = thresh,
        bidirect = TRUE, extend = FALSE, prior = 0.01
    )
    v <- gextract(c("p_max", "p_count"), ivs, iterator = ivs)
    v <- v[order(v$intervalID), ]

    for (sg in c(TRUE, FALSE)) {
        expect_equal(
            gseq.pwm(seqs, pssm,
                mode = "max", bidirect = TRUE, extend = FALSE,
                prior = 0.01, skip_gaps = sg
            ),
            v$p_max,
            tolerance = parity_tolerance, ignore_attr = TRUE
        )
        expect_equal(
            gseq.pwm(seqs, pssm,
                mode = "count", score.thresh = thresh, bidirect = TRUE,
                extend = FALSE, prior = 0.01, skip_gaps = sg
            ),
            v$p_count,
            tolerance = parity_tolerance, ignore_attr = TRUE
        )
    }

    # A gap the projector has to step over gives the ungapped answer back, so
    # the gap-path copy is reached with its projection actually doing something.
    with_gap <- paste0(substr(seqs, 1, 20), "---", substring(seqs, 21))
    expect_equal(
        gseq.pwm(with_gap, pssm,
            mode = "max", bidirect = TRUE, extend = FALSE, prior = 0.01
        ),
        v$p_max,
        tolerance = parity_tolerance, ignore_attr = TRUE
    )
    expect_equal(
        gseq.pwm(with_gap, pssm,
            mode = "count", score.thresh = thresh, bidirect = TRUE,
            extend = FALSE, prior = 0.01
        ),
        v$p_count,
        tolerance = parity_tolerance, ignore_attr = TRUE
    )
})

test_that("gseq.pwm(mode = 'pos') equals the pwm.max.pos vtrack under bidirect = TRUE", {
    remove_all_vtracks()
    pssm <- parity_pssm()
    ivs <- parity_intervals()

    gvtrack.create("p_pos", NULL,
        func = "pwm.max.pos", pssm = pssm,
        bidirect = TRUE, extend = FALSE, prior = 0.01
    )
    v <- gextract("p_pos", ivs, iterator = ivs)
    v <- v[order(v$intervalID), ]

    g <- gseq.pwm(toupper(gseq.extract(ivs)), pssm,
        mode = "pos", bidirect = TRUE, extend = FALSE, prior = 0.01,
        return_strand = TRUE
    )

    # The vtrack signs the position by strand; gseq.pwm returns the two apart.
    expect_equal(g$pos, abs(v$p_pos), ignore_attr = TRUE)
    expect_equal(g$strand, sign(v$p_pos), ignore_attr = TRUE)
})

test_that("a palindromic window is one hit scored on both strands", {
    # A two-base PSSM whose consensus AT is its own reverse complement, with
    # prior = 0: on the window "AT" the forward and reverse scores are both
    # log(1) = 0, so
    #   - the strand union is log(exp(0) + exp(0)) = log 2, not 0;
    #   - the anchor is one hit at score.thresh = 0, not two.
    pssm <- matrix(c(
        1, 0, 0, 0,
        0, 0, 0, 1
    ), ncol = 4, byrow = TRUE)
    colnames(pssm) <- c("A", "C", "G", "T")

    expect_equal(
        gseq.pwm("AT", pssm, mode = "max", bidirect = TRUE, prior = 0),
        log(2)
    )
    expect_equal(
        gseq.pwm("AT", pssm, mode = "count", score.thresh = 0, bidirect = TRUE, prior = 0),
        1
    )

    # Single-strand scoring is untouched: still the plain forward score.
    expect_equal(
        gseq.pwm("AT", pssm, mode = "max", bidirect = FALSE, strand = 1, prior = 0),
        0
    )
    expect_equal(
        gseq.pwm("AT", pssm, mode = "count", score.thresh = 0, bidirect = FALSE, strand = 1, prior = 0),
        1
    )
})
