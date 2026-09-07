create_isolated_test_db()

test_that("the pure-R oracle agrees with motifmodel::motif_score()", {
    skip_if_not(requireNamespace("motifmodel", quietly = TRUE) &&
        is.function(getExportedValue("motifmodel", "motif_score")))

    for (mode in c("full", "sparse", "none")) {
        m <- potts_ref_model(W = 8L, npair_mode = mode, seed = 7L)
        mm <- motifmodel::Potts(m$e, m$J, if (nrow(m$pairs)) m$pairs else NULL,
            intercept = m$intercept
        )
        set.seed(11L)
        seqs <- vapply(1:20, function(i) {
            paste(sample(POTTS_BASES, 8L, replace = TRUE), collapse = "")
        }, character(1))

        ref <- vapply(seqs, potts_ref_window, numeric(1), model = m)
        got <- motifmodel::motif_score(mm, seqs)

        # motifmodel's Potts() projects into the zero-sum gauge, which changes
        # the parameters but not the function, so the SCORES must match even
        # though the tables do not. That is the invariant worth pinning.
        expect_equal(as.numeric(got), as.numeric(ref),
            tolerance = 1e-8,
            info = paste("npair_mode =", mode)
        )
    }
})

test_that("the oracle's RC transform is self-consistent", {
    for (mode in c("full", "sparse", "none")) {
        m <- potts_ref_model(W = 6L, npair_mode = mode, seed = 3L)
        rc <- potts_ref_rc(m)
        set.seed(5L)
        seqs <- vapply(1:50, function(i) {
            paste(sample(POTTS_BASES, 6L, replace = TRUE), collapse = "")
        }, character(1))
        for (s in seqs) {
            expect_equal(potts_ref_window(s, rc),
                potts_ref_window(grevcomp(s), m),
                tolerance = 1e-12,
                ignore_attr = TRUE,
                info = paste(mode, s)
            )
        }
        # and it is an involution
        expect_equal(potts_ref_rc(rc)$e, m$e, tolerance = 1e-12)
        expect_equal(potts_ref_rc(rc)$intercept, m$intercept)
    }
})

test_that("gseq.potts scores a single window like the oracle", {
    for (W in c(4L, 8L, 20L)) {
        for (mode in c("full", "sparse", "none")) {
            m <- potts_ref_model(W = W, npair_mode = mode, seed = 13L)
            set.seed(17L)
            seqs <- vapply(1:30, function(i) {
                paste(sample(POTTS_BASES, W, replace = TRUE), collapse = "")
            }, character(1))
            ref <- vapply(seqs, potts_ref_window, numeric(1), model = m)

            # A single window: lse and max over one anchor are both that anchor.
            expect_equal(gseq.potts(seqs, m, mode = "max", bidirect = FALSE, strand = 1L),
                as.numeric(ref),
                tolerance = 1e-6,
                info = paste("W", W, mode)
            )
            expect_equal(gseq.potts(seqs, m, mode = "lse", bidirect = FALSE, strand = 1L),
                as.numeric(ref),
                tolerance = 1e-6,
                info = paste("W", W, mode)
            )
        }
    }
})

test_that("gseq.potts reverse strand equals the oracle's RC twin", {
    m <- potts_ref_model(W = 8L, npair_mode = "full", seed = 23L)
    rc <- potts_ref_rc(m)
    set.seed(29L)
    seqs <- vapply(1:40, function(i) {
        paste(sample(POTTS_BASES, 8L, replace = TRUE), collapse = "")
    }, character(1))

    # scoring the reverse strand of s == scoring the twin on s
    expect_equal(
        gseq.potts(seqs, m, mode = "max", bidirect = FALSE, strand = -1L),
        gseq.potts(seqs, rc, mode = "max", bidirect = FALSE, strand = 1L),
        tolerance = 1e-6
    )
    # == scoring the original model on revcomp(s)
    expect_equal(
        gseq.potts(seqs, m, mode = "max", bidirect = FALSE, strand = -1L),
        gseq.potts(grevcomp(seqs), m, mode = "max", bidirect = FALSE, strand = 1L),
        tolerance = 1e-6
    )
})

test_that("gseq.potts bidirect unions the two strands", {
    m <- potts_ref_model(W = 8L, npair_mode = "full", seed = 31L)
    set.seed(37L)
    seqs <- vapply(1:40, function(i) {
        paste(sample(POTTS_BASES, 8L, replace = TRUE), collapse = "")
    }, character(1))

    fwd <- gseq.potts(seqs, m, mode = "max", bidirect = FALSE, strand = 1L)
    rev <- gseq.potts(seqs, m, mode = "max", bidirect = FALSE, strand = -1L)

    # lse and max both lse-union the strands at one anchor (design doc's table)
    expect_equal(gseq.potts(seqs, m, mode = "max", bidirect = TRUE),
        vapply(seq_along(fwd), function(i) log_sum_exp(c(fwd[i], rev[i])), numeric(1)),
        tolerance = 1e-6
    )
    expect_equal(gseq.potts(seqs, m, mode = "lse", bidirect = TRUE),
        vapply(seq_along(fwd), function(i) log_sum_exp(c(fwd[i], rev[i])), numeric(1)),
        tolerance = 1e-6
    )
})

test_that("gseq.potts reduces over anchors in a longer sequence", {
    m <- potts_ref_model(W = 6L, npair_mode = "full", seed = 41L)
    set.seed(43L)
    seqs <- vapply(1:10, function(i) {
        paste(sample(POTTS_BASES, 40L, replace = TRUE), collapse = "")
    }, character(1))

    for (bi in c(TRUE, FALSE)) {
        ref_lse <- vapply(seqs, function(s) {
            log_sum_exp(potts_ref_anchors(s, m, bidirect = bi, strand = 1L, union = "lse"))
        }, numeric(1))
        ref_max <- vapply(seqs, function(s) {
            max(potts_ref_anchors(s, m, bidirect = bi, strand = 1L, union = "lse"))
        }, numeric(1))

        expect_equal(gseq.potts(seqs, m, mode = "lse", bidirect = bi, strand = 1L),
            as.numeric(ref_lse),
            tolerance = 1e-5, info = paste("bidirect", bi)
        )
        expect_equal(gseq.potts(seqs, m, mode = "max", bidirect = bi, strand = 1L),
            as.numeric(ref_max),
            tolerance = 1e-5, info = paste("bidirect", bi)
        )
    }
})

test_that("gseq.potts refuses to score a window with a non-ACGT base", {
    m <- potts_ref_model(W = 4L, npair_mode = "full", seed = 47L)

    # the only window is unscorable -> NA
    expect_true(is.na(gseq.potts("ACNT", m, mode = "max", bidirect = FALSE, strand = 1L)))
    expect_true(is.na(gseq.potts(NA_character_, m, mode = "max")))
    # too short for one window -> NA
    expect_true(is.na(gseq.potts("ACG", m, mode = "max")))

    # an N excludes its own anchors but not the clean ones
    s <- "ACGTACGTNACGTACGT"
    ref <- potts_ref_anchors(s, m, bidirect = FALSE, strand = 1L)
    expect_equal(gseq.potts(s, m, mode = "max", bidirect = FALSE, strand = 1L),
        max(ref, na.rm = TRUE),
        tolerance = 1e-6
    )
    expect_equal(gseq.potts(s, m, mode = "lse", bidirect = FALSE, strand = 1L),
        log_sum_exp(ref[!is.na(ref)]),
        tolerance = 1e-6
    )
})

test_that("gseq.potts accepts a whole motifmodel Potts as the model", {
    skip_if_not(requireNamespace("motifmodel", quietly = TRUE) &&
        is.function(getExportedValue("motifmodel", "motif_score")))

    m <- potts_ref_model(W = 8L, npair_mode = "full", seed = 53L)
    mm <- motifmodel::Potts(m$e, m$J, m$pairs, intercept = m$intercept)
    set.seed(59L)
    seqs <- vapply(1:20, function(i) {
        paste(sample(POTTS_BASES, 8L, replace = TRUE), collapse = "")
    }, character(1))

    # the object goes in verbatim - width, pair_strength, attr and link and all
    expect_equal(
        gseq.potts(seqs, mm, mode = "max", bidirect = FALSE, strand = 1L),
        as.numeric(motifmodel::motif_score(mm, seqs)),
        tolerance = 1e-6
    )
})

test_that("gseq.potts mode = pos reports a best-scoring window", {
    m <- potts_ref_model(W = 6L, npair_mode = "full", seed = 61L)
    set.seed(67L)
    seqs <- vapply(1:20, function(i) {
        paste(sample(POTTS_BASES, 30L, replace = TRUE), collapse = "")
    }, character(1))
    W <- nrow(m$e)

    # Assert the SCORE at the reported position, never the index. Two anchors
    # can tie exactly - under bidirect, a window and its reverse complement
    # always do - and then the index depends on tie-break order, which the
    # kernel and the oracle resolve differently for reasons of floating-point
    # summation rather than of logic.
    got_f <- gseq.potts(seqs, m, mode = "pos", bidirect = FALSE, strand = 1L)
    for (k in seq_along(seqs)) {
        a <- potts_ref_anchors(seqs[k], m, bidirect = FALSE, strand = 1L)
        i <- got_f[k]
        expect_true(i >= 1 && i <= length(a), info = seqs[k])
        expect_equal(as.numeric(a[i]), max(a),
            tolerance = 1e-6, ignore_attr = TRUE, info = seqs[k]
        )
    }

    # bidirect: the magnitude names a maximiser of the per-anchor strand
    # MAXIMUM (not the lse), and the sign names the winning strand there.
    got_b <- gseq.potts(seqs, m, mode = "pos", bidirect = TRUE)
    rc <- potts_ref_rc(m)
    for (k in seq_along(seqs)) {
        a <- potts_ref_anchors(seqs[k], m, bidirect = TRUE, strand = 1L, union = "max")
        i <- abs(got_b[k])
        expect_true(i >= 1 && i <= length(a), info = seqs[k])
        expect_equal(as.numeric(a[i]), max(a),
            tolerance = 1e-6, ignore_attr = TRUE, info = seqs[k]
        )
        win <- substr(seqs[k], i, i + W - 1L)
        f <- as.numeric(potts_ref_window(win, m))
        r <- as.numeric(potts_ref_window(win, rc))
        # Only assert the strand where the two strands actually differ; a
        # palindromic window ties and then either sign is correct.
        if (abs(f - r) > 1e-9) {
            expect_equal(sign(got_b[k]), if (r > f) -1 else 1, info = seqs[k])
        }
    }
})

test_that("gseq.potts mode = count needs a threshold and counts anchors", {
    m <- potts_ref_model(W = 6L, npair_mode = "full", seed = 71L)
    s <- paste(sample(POTTS_BASES, 200L, replace = TRUE), collapse = "")

    expect_error(gseq.potts(s, m, mode = "count"), "score.thresh")

    a <- potts_ref_anchors(s, m, bidirect = TRUE, strand = 1L, union = "lse")
    for (th in stats::quantile(a, c(0.1, 0.5, 0.9), names = FALSE)) {
        expect_equal(gseq.potts(s, m, mode = "count", score.thresh = th),
            sum(a >= th),
            info = paste("thresh", th)
        )
    }

    # a sequence too short for one window counts nothing rather than NA
    expect_equal(gseq.potts("ACG", m, mode = "count", score.thresh = 0), 0)
})

test_that("gseq.potts pos and max may disagree on which anchor wins", {
    # The strand union differs between the two - lse for max, maximum for pos -
    # so they can select different anchors. Inherited from the pwm family
    # (pwm.max passes combine_strands = true, pwm.max.pos does not) and pinned
    # here so neither half gets "corrected" alone.
    m <- potts_ref_model(W = 4L, npair_mode = "full", seed = 79L)
    set.seed(83L)
    disagreed <- FALSE
    for (i in 1:200) {
        s <- paste(sample(POTTS_BASES, 20L, replace = TRUE), collapse = "")
        lse_a <- potts_ref_anchors(s, m, bidirect = TRUE, union = "lse")
        max_a <- potts_ref_anchors(s, m, bidirect = TRUE, union = "max")

        # max mode reduces the LSE-union over anchors
        expect_equal(gseq.potts(s, m, mode = "max", bidirect = TRUE),
            max(lse_a),
            tolerance = 1e-6, ignore_attr = TRUE, info = s
        )

        # pos mode names a maximiser of the MAX-union
        p <- abs(gseq.potts(s, m, mode = "pos", bidirect = TRUE))
        expect_equal(as.numeric(max_a[p]), max(max_a),
            tolerance = 1e-6, ignore_attr = TRUE, info = s
        )

        # the wart itself: where the two unions rank anchors differently, the
        # anchor pos names scores strictly below the lse maximum, so max and
        # pos genuinely refer to different places.
        if (as.numeric(lse_a[p]) < max(lse_a) - 1e-6) {
            disagreed <- TRUE
        }
    }
    expect_true(disagreed,
        info = "no disagreement in 200 draws - the fixture no longer exercises the wart"
    )
})

test_that("gseq.potts mode = pos works on the reverse strand alone", {
    m <- potts_ref_model(W = 6L, npair_mode = "full", seed = 149L)
    set.seed(151L)
    seqs <- vapply(1:30, function(i) {
        paste(sample(POTTS_BASES, 30L, replace = TRUE), collapse = "")
    }, character(1))

    # bidirect = FALSE, strand = -1: the reported position is an unsigned
    # 1-based anchor index into the FORWARD sequence, and the anchor it names
    # must be a maximiser of the REVERSE-strand per-anchor scores. Asserting
    # the score rather than the index, for the tie reason established earlier.
    got <- gseq.potts(seqs, m, mode = "pos", bidirect = FALSE, strand = -1L)

    differed <- FALSE
    for (k in seq_along(seqs)) {
        rev_a <- potts_ref_anchors(seqs[k], m, bidirect = FALSE, strand = -1L)
        fwd_a <- potts_ref_anchors(seqs[k], m, bidirect = FALSE, strand = 1L)
        p <- got[k]
        expect_true(p >= 1 && p <= length(rev_a), info = seqs[k])
        expect_equal(as.numeric(rev_a[p]), max(rev_a),
            tolerance = 1e-6, ignore_attr = TRUE, info = seqs[k]
        )
        # The point of the test: a kernel that scored the forward strand here
        # would name a forward maximiser instead. Record whether the two
        # strands actually disagree about which anchor wins, so the assertion
        # above is known to be discriminating rather than accidentally
        # satisfied by a forward-scoring kernel.
        if (max(fwd_a) - as.numeric(fwd_a[p]) > 1e-6) {
            differed <- TRUE
        }
    }
    expect_true(differed,
        info = "the reverse-strand argmax never differed from the forward one - this fixture cannot distinguish the two strands"
    )
})
