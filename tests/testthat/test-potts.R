create_isolated_test_db()

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
    # NOT quantile(): at n = 195 a type-7 quantile lands exactly ON an anchor's
    # own score, so `a >= th` turns on rounding across the float/double
    # boundary - which is what test-potts-vtrack.R forbids by name. A midpoint
    # between two distinct sorted scores cannot tie with either.
    sa <- sort(unique(a))
    idx <- pmin(pmax(round(c(0.1, 0.5, 0.9) * (length(sa) - 1L)) + 1L, 1L), length(sa) - 1L)
    for (th in (sa[idx] + sa[idx + 1L]) / 2) {
        expect_equal(gseq.potts(s, m, mode = "count", score.thresh = th),
            sum(a >= th),
            info = paste("thresh", th)
        )
    }

    # "counted nothing" and "nothing to count" are different answers: a
    # sequence with windows but none of them scorable counts 0, a sequence too
    # short to hold one window - or NA - has nothing to count and is NA, for
    # count as for every other mode
    expect_equal(gseq.potts(strrep("N", 20L), m, mode = "count", score.thresh = 0), 0)
    expect_true(is.na(gseq.potts("ACG", m, mode = "count", score.thresh = 0)))
    expect_true(is.na(gseq.potts(NA_character_, m, mode = "count", score.thresh = 0)))
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

test_that("the blocked kernel agrees with the naive one on 1e6 random windows at W = 20, full pairwise", {
    # PottsModel.h's score_codes() only dispatches to score_codes_blocked()
    # when it has FEWER lookups than score_codes_naive() for this model - a
    # dense model like this one (W = 20, every pair present) is exactly the
    # case that wins, measured well past the point where it's worth it,
    # through gseq.potts() on an idle host. This is the equivalence check
    # that decision was conditioned on, kept as a permanent regression test -
    # C_potts_score_codes_cmp is a test-only entry point that scores an
    # integer code matrix with BOTH kernels directly, without a DNA string or
    # gseq.potts() in between.
    m <- potts_ref_model(W = 20L, npair_mode = "full", seed = 1L)
    params <- misha:::.potts_params(m,
        bidirect = TRUE, extend = FALSE, strand = 1L, score.thresh = 0,
        what = "kernel equivalence test"
    )

    set.seed(2L)
    W <- 20L
    n <- 1e6
    codes <- matrix(sample(0:3, W * n, replace = TRUE), nrow = W, ncol = n)
    storage.mode(codes) <- "integer"

    res <- .Call("C_potts_score_codes_cmp", params, codes)
    max_diff <- max(abs(res$naive - res$blocked))
    expect_true(max_diff < 1e-9, info = sprintf("max |naive - blocked| = %.3e over %d windows", max_diff, n))
})

# The trailing size-1 block build_blocked_tables() uses for an odd W has no
# coverage anywhere else in this file - every W above (4, 6, 8, 20) is even.
# Sweep W = 20/21 (even/odd) x every pairing density: score_codes_blocked()
# must agree with score_codes_naive() at EVERY density, even the ones where
# PottsModel.h's gate picks naive for production score_codes() - the
# equivalence check goes through C_potts_score_codes_cmp, which calls both
# kernels directly and does not go through that gate.
for (.W in c(20L, 21L)) {
    for (.npair_mode in c("full", "sparse", "none")) {
        test_that(sprintf("blocked kernel agrees with naive: W = %d, npair_mode = %s", .W, .npair_mode), {
            m <- potts_ref_model(W = .W, npair_mode = .npair_mode, seed = 3L)
            params <- misha:::.potts_params(m,
                bidirect = TRUE, extend = FALSE, strand = 1L, score.thresh = 0,
                what = "kernel equivalence test"
            )

            set.seed(4L)
            n <- 50000
            codes <- matrix(sample(0:3, .W * n, replace = TRUE), nrow = .W, ncol = n)
            storage.mode(codes) <- "integer"

            res <- .Call("C_potts_score_codes_cmp", params, codes)
            max_diff <- max(abs(res$naive - res$blocked))
            expect_true(max_diff < 1e-9,
                info = sprintf("W=%d %s: max |naive - blocked| = %.3e over %d windows", .W, .npair_mode, max_diff, n)
            )
        })
    }
}

test_that("a repeated pair contributes BOTH couplings, in both kernels", {
    # .coerce_potts_model() rejects a duplicate (p1, p2), so this cannot be
    # reached through gseq.potts() or a vtrack. PottsParams::parse() is a
    # .Call trust boundary all the same, and build_blocked_tables() SUMS a
    # repeat rather than letting the last one win - which is what
    # score_codes_naive()'s loop over k does naturally. Without this test the
    # sum can be replaced by an overwrite and the whole suite still passes.
    W <- 6L
    m <- potts_ref_model(W = W, npair_mode = "sparse", seed = 11L)
    single <- misha:::.potts_params(m,
        bidirect = TRUE, extend = FALSE, strand = 1L, score.thresh = 0,
        what = "duplicate pair test"
    )

    dup <- single
    dup$pairs <- rbind(dup$pairs, dup$pairs[1L, , drop = FALSE])
    dup$J <- rbind(dup$J, dup$J[1L, , drop = FALSE])

    set.seed(12L)
    n <- 2000L
    codes <- matrix(sample(0:3, W * n, replace = TRUE), nrow = W, ncol = n)
    storage.mode(codes) <- "integer"

    r1 <- .Call("C_potts_score_codes_cmp", single, codes)
    r2 <- .Call("C_potts_score_codes_cmp", dup, codes)

    # The duplicated block is counted twice, so the difference is exactly that
    # coupling's own contribution to each window - not zero, and not double the
    # whole score.
    J1 <- matrix(single$J[1L, ], 4L, 4L)
    extra <- J1[cbind(codes[single$pairs[1L, 1L], ] + 1L, codes[single$pairs[1L, 2L], ] + 1L)]
    expect_equal(r2$naive - r1$naive, extra)
    expect_equal(r2$blocked, r2$naive, tolerance = 1e-9)
    expect_false(isTRUE(all.equal(r2$naive, r1$naive)))
})

test_that("a labelled J block is reordered by name, like e's columns", {
    # .coerce_pssm_matrix() reorders `e` by column NAME, so reading a labelled
    # J block positionally would score a relabelled copy of the same model
    # differently, with no warning.
    m <- potts_ref_model(W = 4L, npair_mode = "sparse", seed = 21L)
    s <- c("ACGT", "TGCA", "GGTA")
    want <- gseq.potts(s, m, mode = "max", bidirect = FALSE, strand = 1L)

    o <- c(4L, 3L, 2L, 1L) # T, G, C, A
    relabelled <- m
    relabelled$e <- m$e[, o, drop = FALSE]
    relabelled$J <- lapply(m$J, function(Jk) {
        out <- Jk[o, o, drop = FALSE]
        dimnames(out) <- list(POTTS_BASES[o], POTTS_BASES[o])
        out
    })
    expect_equal(
        gseq.potts(s, relabelled, mode = "max", bidirect = FALSE, strand = 1L),
        want
    )

    # An unlabelled block stays positional, as before.
    bare <- m
    bare$J <- lapply(m$J, function(Jk) {
        dimnames(Jk) <- NULL
        Jk
    })
    expect_equal(gseq.potts(s, bare, mode = "max", bidirect = FALSE, strand = 1L), want)

    # Labelled with anything but the four bases is an error, not a guess.
    bad <- m
    bad$J <- lapply(m$J, function(Jk) {
        dimnames(Jk) <- list(c("A", "C", "G", "X"), POTTS_BASES)
        Jk
    })
    expect_error(
        gseq.potts(s, bad, mode = "max", bidirect = FALSE, strand = 1L),
        "named A, C, G, T"
    )

    # Half-labelled is too: it is as likely to be an accident as an intent.
    half <- m
    half$J <- lapply(m$J, function(Jk) {
        dimnames(Jk) <- list(POTTS_BASES, NULL)
        Jk
    })
    expect_error(
        gseq.potts(s, half, mode = "max", bidirect = FALSE, strand = 1L),
        "named A, C, G, T"
    )
})

test_that("a fractional pairs index is rejected, not truncated", {
    # storage.mode(pairs) <- "integer" truncates, and the range check that
    # follows would then validate the truncated value: 2.9 with W = 2 is both
    # fractional AND out of range, and used to be accepted as 2.
    m <- potts_ref_model(W = 4L, npair_mode = "sparse", seed = 31L)
    m$pairs <- matrix(c(1L, 2L), ncol = 2L)
    m$J <- m$J[1L]

    frac <- m
    frac$pairs <- matrix(c(1, 2.9), ncol = 2L)
    expect_error(
        gseq.potts("ACGT", frac, mode = "max", bidirect = FALSE, strand = 1L),
        "whole number"
    )

    # and the whole-number equivalent still passes
    ok <- m
    ok$pairs <- matrix(c(1, 2), ncol = 2L)
    expect_silent(gseq.potts("ACGT", ok, mode = "max", bidirect = FALSE, strand = 1L))
})

test_that("gseq.potts(bidirect = FALSE) requires an explicit strand", {
    # The signature default of 0 exists for parity with gseq.pwm(). Under
    # bidirect = FALSE it used to mean the forward strand silently, which hands
    # half an answer to anyone carrying gseq.pwm()'s "0 = both strands" over.
    m <- potts_ref_model(W = 4L, npair_mode = "sparse", seed = 41L)
    expect_error(gseq.potts("ACGT", m, mode = "max", bidirect = FALSE), "explicit strand")
    expect_silent(gseq.potts("ACGT", m, mode = "max", bidirect = FALSE, strand = 1L))
    expect_silent(gseq.potts("ACGT", m, mode = "max", bidirect = TRUE)) # 0 is fine here
})
