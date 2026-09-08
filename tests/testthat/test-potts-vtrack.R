create_isolated_test_db()

test_that("the potts params handler validates and normalises", {
    m <- potts_ref_model(W = 5L, npair_mode = "sparse", seed = 89L)

    p <- misha:::.vtrack_params_potts("potts", m, list())
    expect_equal(dim(p$e), c(5L, 4L))
    expect_equal(colnames(p$e), c("A", "C", "G", "T"))
    expect_equal(dim(p$J), c(nrow(m$pairs), 16L))
    expect_equal(dim(p$pairs), dim(m$pairs))
    expect_true(is.integer(p$pairs))
    expect_equal(p$intercept, m$intercept)
    expect_true(p$bidirect)
    expect_true(p$extend)

    # the flat J is column-major per block, so J[[k]][a, b] lands at (b-1)*4+a
    k <- 1L
    expect_equal(p$J[k, (3L - 1L) * 4L + 2L], m$J[[k]][2L, 3L])
})

test_that("the potts params handler tolerates a whole motifmodel Potts", {
    skip_if_not(requireNamespace("motifmodel", quietly = TRUE))
    m <- potts_ref_model(W = 5L, npair_mode = "full", seed = 97L)
    mm <- motifmodel::Potts(m$e, m$J, m$pairs, intercept = m$intercept)

    # width, pair_strength, attr and link ride along and are ignored
    expect_silent(p <- misha:::.vtrack_params_potts("potts", mm, list()))
    expect_equal(nrow(p$e), 5L)
})

test_that("the potts params handler rejects what it should", {
    m <- potts_ref_model(W = 5L, npair_mode = "full", seed = 101L)

    expect_error(misha:::.vtrack_params_potts("potts", list(J = m$J), list()), "'e'")

    bad <- m
    bad$width <- 7L
    expect_error(misha:::.vtrack_params_potts("potts", bad, list()), "width")

    bad <- m
    bad$pairs[1L, ] <- c(3L, 2L)
    expect_error(misha:::.vtrack_params_potts("potts", bad, list()), "lower position first")

    bad <- m
    bad$pairs[2L, ] <- bad$pairs[1L, ]
    expect_error(misha:::.vtrack_params_potts("potts", bad, list()), "twice")

    bad <- m
    bad$J <- bad$J[-1L]
    expect_error(misha:::.vtrack_params_potts("potts", bad, list()), "tables")

    # spatial parameters are not part of the potts family and must be named
    expect_error(
        misha:::.vtrack_params_potts("potts", NULL, c(m, list(spat_factor = c(1, 2)))),
        "spat_factor"
    )
    # and potts.count has no default threshold
    expect_error(misha:::.vtrack_params_potts("potts.count", m, list()), "score.thresh")
})

test_that("a potts model whose window score cannot fit a float is rejected", {
    # Every entry is finite, so the per-entry checks pass; what overflows is
    # the WINDOW. potts scores are reported in float, and the two scoring paths
    # then answer differently. Measured on the first model below before the
    # bound went in, over chr1:3000-3200 with a 20 bp iterator shift:
    # potts.max read Inf from a freshly seeded interval and NaN from a slid
    # one, and potts.max.pos 1 against NaN - RunningLogSumExp::value() and
    # RunningMaxDeque::value() report a non-finite running maximum as -inf, and
    # PottsScorer::slid_answer() maps that to NaN, which is exactly the
    # "nothing was scorable" answer it is supposed to mean.
    m <- potts_ref_model(W = 6L, npair_mode = "none", seed = 167L)

    over <- m
    over$e[] <- 1e38 # 6 * 1e38, past FLT_MAX
    expect_error(misha:::.vtrack_params_potts("potts", over, list()), "single-precision")
    # both entry points go through the one validator
    expect_error(gseq.potts("ACGTACGTAC", over, mode = "max"), "single-precision")

    # J and the intercept count towards the same bound
    over <- potts_ref_model(W = 4L, npair_mode = "full", seed = 167L)
    over$J <- lapply(over$J, function(Jk) {
        Jk[] <- 1e38
        Jk
    })
    expect_error(misha:::.vtrack_params_potts("potts", over, list()), "single-precision")

    over <- m
    over$intercept <- 1e39
    expect_error(misha:::.vtrack_params_potts("potts", over, list()), "single-precision")

    # A large model that does fit is left alone, so the bound is not merely
    # "no big numbers": 6 * 1e37 is inside FLT_MAX, and both scoring paths were
    # measured to agree on it over the same scan.
    ok <- m
    ok$e[] <- 1e37
    expect_silent(misha:::.vtrack_params_potts("potts", ok, list()))
})

test_that("potts and potts.max vtracks equal gseq.potts on the same sequence", {
    remove_all_vtracks()
    m <- potts_ref_model(W = 6L, npair_mode = "full", seed = 103L)

    iv <- gintervals(1, 200, 300)
    # extend = TRUE reads motif_len - 1 bases past the end, so the scored
    # sequence is longer than the interval. That is the sequence to compare on.
    seq_ext <- toupper(gseq.extract(gintervals(1, 200, 300 + 6L - 1L)))

    gvtrack.create("p_lse", NULL, "potts", params = c(m, list(bidirect = TRUE, extend = TRUE)))
    gvtrack.create("p_max", NULL, "potts.max", params = c(m, list(bidirect = TRUE, extend = TRUE)))

    got <- gextract(c("p_lse", "p_max"), iv, iterator = iv)
    expect_equal(got$p_lse, gseq.potts(seq_ext, m, mode = "lse", bidirect = TRUE),
        tolerance = 1e-5
    )
    expect_equal(got$p_max, gseq.potts(seq_ext, m, mode = "max", bidirect = TRUE),
        tolerance = 1e-5
    )
})

test_that("an order-1 potts vtrack reproduces the pwm family", {
    remove_all_vtracks()
    # A strictly positive PSSM: prior = 0 on a zero entry gives -Inf, and an
    # -Inf comparison proves nothing.
    set.seed(107L)
    pssm <- matrix(runif(8L * 4L, 0.05, 1),
        ncol = 4L,
        dimnames = list(NULL, c("A", "C", "G", "T"))
    )
    pssm <- pssm / rowSums(pssm)

    # An order-1 Potts whose e is log(pssm) is exactly a PWM with prior = 0.
    m <- list(
        e = log(pssm), J = list(),
        pairs = matrix(integer(0), 0L, 2L), intercept = 0
    )

    iv <- gintervals(1, 1000, 1400)
    for (bi in c(TRUE, FALSE)) {
        remove_all_vtracks()
        gvtrack.create("w_lse", NULL, "pwm",
            pssm = pssm, prior = 0, bidirect = bi, extend = TRUE, strand = 1
        )
        gvtrack.create("w_max", NULL, "pwm.max",
            pssm = pssm, prior = 0, bidirect = bi, extend = TRUE, strand = 1
        )
        gvtrack.create("t_lse", NULL, "potts",
            params = c(m, list(bidirect = bi, extend = TRUE, strand = 1))
        )
        gvtrack.create("t_max", NULL, "potts.max",
            params = c(m, list(bidirect = bi, extend = TRUE, strand = 1))
        )

        got <- gextract(c("w_lse", "w_max", "t_lse", "t_max"), iv, iterator = 10)
        expect_equal(got$t_lse, got$w_lse, tolerance = 1e-4, info = paste("bidirect", bi))
        expect_equal(got$t_max, got$w_max, tolerance = 1e-4, info = paste("bidirect", bi))
    }
})

test_that("potts vtracks honour extend", {
    remove_all_vtracks()
    m <- potts_ref_model(W = 20L, npair_mode = "full", seed = 109L)

    gvtrack.create("t_no", NULL, "potts.max", params = c(m, list(extend = FALSE)))
    gvtrack.create("t_yes", NULL, "potts.max", params = c(m, list(extend = TRUE)))

    short_iv <- gintervals(1, 500, 510) # 10 bp, narrower than W = 20
    got <- gextract(c("t_no", "t_yes"), short_iv, iterator = short_iv)
    expect_true(is.na(got$t_no))
    expect_false(is.na(got$t_yes))
})

test_that("potts.max.pos and potts.count work on a genome interval", {
    remove_all_vtracks()
    W <- 6L
    m <- potts_ref_model(W = W, npair_mode = "full", seed = 113L)

    iv <- gintervals(1, 2000, 2200)
    seq_ext <- toupper(gseq.extract(gintervals(1, 2000, 2200 + W - 1L)))
    a_max <- potts_ref_anchors(seq_ext, m, bidirect = TRUE, union = "max")
    gvtrack.create("t_pos", NULL, "potts.max.pos", params = c(m, list(extend = TRUE)))

    got <- gextract("t_pos", iv, iterator = iv)

    # Assert the SCORE at the reported anchor, not the argmax index: under
    # bidirect a window and its reverse complement give the same max-union, so
    # two anchors can tie exactly and either is a correct answer.
    p <- abs(got$t_pos)
    expect_true(p >= 1 && p <= length(a_max))
    expect_equal(as.numeric(a_max[p]), max(a_max), tolerance = 1e-5, ignore_attr = TRUE)

    # The sign names the winning strand, but only where the two strands differ:
    # a palindromic window ties and then either sign is right.
    win <- substr(seq_ext, p, p + W - 1L)
    f <- as.numeric(potts_ref_window(win, m))
    r <- as.numeric(potts_ref_window(win, potts_ref_rc(m)))
    if (abs(f - r) > 1e-9) {
        expect_equal(sign(got$t_pos), if (r > f) -1 else 1)
    }
})

test_that("potts.count counts the anchors at or above its threshold", {
    remove_all_vtracks()
    W <- 6L
    m <- potts_ref_model(W = W, npair_mode = "full", seed = 113L)

    iv <- gintervals(1, 2000, 2200)
    seq_ext <- toupper(gseq.extract(gintervals(1, 2000, 2200 + W - 1L)))
    a_lse <- potts_ref_anchors(seq_ext, m, bidirect = TRUE, union = "lse")
    srt <- sort(a_lse)
    n <- length(srt)

    # DO NOT replace this with quantile(a_lse, 0.9): a threshold at or near an
    # anchor's own score is not a testable quantity here. Under bidirect a
    # window and its reverse complement produce the same two per-strand numbers
    # swapped, so their unions are mathematically identical - and that holds for
    # the LOG-SUM-EXP union too, not only the maximum. a_lse therefore contains
    # exactly tied values, C++ and R round each tie independently, and an anchor
    # sitting ON the threshold lands either side of it: measured with a
    # quantile-derived threshold, R counted 21 of 200 anchors and C++ 19.
    # Threshold in the widest gap of the upper tail instead, assert the gap is
    # wide enough that no rounding can cross it, and pin both ends exactly
    # below.
    idx <- seq.int(floor(0.8 * n), n - 1L)
    j <- idx[which.max(srt[idx + 1L] - srt[idx])]
    th <- (srt[j] + srt[j + 1L]) / 2
    expect_gt(srt[j + 1L] - srt[j], 1e-3)

    # Below the minimum every scorable anchor counts, which pins the size of
    # the anchor set - the interval's length, not the extended fetch's.
    gvtrack.create("c_all", NULL, "potts.count",
        params = c(m, list(extend = TRUE, score.thresh = min(a_lse) - 1))
    )
    gvtrack.create("c_tail", NULL, "potts.count",
        params = c(m, list(extend = TRUE, score.thresh = th))
    )
    gvtrack.create("c_none", NULL, "potts.count",
        params = c(m, list(extend = TRUE, score.thresh = max(a_lse) + 1))
    )

    got <- gextract(c("c_all", "c_tail", "c_none"), iv, iterator = iv)
    expect_equal(got$c_all, n)
    expect_equal(got$c_tail, n - j)
    expect_equal(got$c_none, 0)
})

test_that("a potts vtrack over an all-N interval is NaN", {
    remove_all_vtracks()
    m <- potts_ref_model(W = 6L, npair_mode = "full", seed = 127L)

    # chr20 of the test fixture opens with an N run (assembly gap)
    iv <- gintervals(20, 0, 200)
    skip_if_not(
        grepl("^N+$", toupper(gseq.extract(iv))),
        "no all-N interval at the start of chr20 in this fixture"
    )

    gvtrack.create("t_lse", NULL, "potts", params = c(m, list(extend = FALSE)))
    gvtrack.create("t_cnt", NULL, "potts.count",
        params = c(m, list(extend = FALSE, score.thresh = 0))
    )
    got <- gextract(c("t_lse", "t_cnt"), iv, iterator = iv)
    expect_true(is.na(got$t_lse))
    expect_equal(got$t_cnt, 0)
})

test_that("a filtered potts vtrack aggregates instead of refusing", {
    remove_all_vtracks()
    m <- potts_ref_model(W = 6L, npair_mode = "sparse", seed = 131L)

    iv <- gintervals(1, 2000, 2200)
    mask <- gintervals(1, 2050, 2100)

    gvtrack.create("t_flt", NULL, "potts", params = c(m, list(extend = TRUE)))
    gvtrack.filter("t_flt", filter = mask)
    got <- gextract("t_flt", iv, iterator = iv)$t_flt

    # the two parts gvtrack.filter's complement leaves, scored on their own
    gvtrack.create("t_a", NULL, "potts", params = c(m, list(extend = TRUE)))
    gvtrack.create("t_b", NULL, "potts", params = c(m, list(extend = TRUE)))
    part_a <- gextract("t_a", gintervals(1, 2000, 2050), iterator = gintervals(1, 2000, 2050))$t_a
    part_b <- gextract("t_b", gintervals(1, 2100, 2200), iterator = gintervals(1, 2100, 2200))$t_b

    expect_equal(got, log_sum_exp(c(part_a, part_b)), tolerance = 1e-5)
})

test_that("potts.max.pos reports the argmax anchor when one strand is scored", {
    remove_all_vtracks()
    W <- 6L
    m <- potts_ref_model(W = W, npair_mode = "full", seed = 113L)

    iv <- gintervals(1, 2000, 2200)
    seq_ext <- toupper(gseq.extract(gintervals(1, 2000, 2200 + W - 1L)))
    a_fwd <- potts_ref_anchors(seq_ext, m, bidirect = FALSE, strand = 1L, union = "max")

    # bidirect = FALSE is what makes the argmax INDEX a testable quantity: with
    # only the forward strand scored there is no window/reverse-complement tie,
    # so the maximiser is unique and the position arithmetic is pinned exactly.
    expect_equal(sum(a_fwd == max(a_fwd)), 1L)

    gvtrack.create("t_pos", NULL, "potts.max.pos",
        params = c(m, list(bidirect = FALSE, strand = 1, extend = TRUE))
    )
    gvtrack.create("t_max", NULL, "potts.max",
        params = c(m, list(bidirect = FALSE, strand = 1, extend = TRUE))
    )

    got <- gextract(c("t_pos", "t_max"), iv, iterator = iv)
    # unsigned: the sign is only applied when bidirect
    expect_equal(got$t_pos, as.numeric(which.max(a_fwd)))
    expect_equal(got$t_max, max(a_fwd), tolerance = 1e-5, ignore_attr = TRUE)
})

test_that("the reverse orientation is scored and reported in forward coordinates", {
    remove_all_vtracks()
    W <- 6L
    m <- potts_ref_model(W = W, npair_mode = "full", seed = 113L)

    # The twin of the test above, and the ONLY suite coverage of PottsScorer's
    # m_strand == -1 half: the model/rc swap that reads the original reverse
    # strand off an already reverse-complemented target, and the reverse ->
    # forward coordinate mapping that turns a target index back into a
    # forward-strand position. bidirect = TRUE cannot reach either, because
    # .potts_params() clamps strand to 1 there (both strands are scored at
    # every anchor anyway), and C_gseq_potts cannot: it scores a string, so it
    # has no fetch, no reverse-complemented target and no position mapping.
    iv <- gintervals(1, 2000, 2200)
    seq_ext <- toupper(gseq.extract(gintervals(1, 2000, 2200 + W - 1L)))
    a_rev <- potts_ref_anchors(seq_ext, m, bidirect = FALSE, strand = -1L)

    # Same gapped threshold as the count test above, for the same reason: even
    # on one strand a repeated k-mer gives exactly tied anchor scores (8 of the
    # 200 here are duplicates), and a threshold landing on one is decided by
    # the last bit of a double in C++ and in R separately. Taken from the
    # middle of the distribution so the count cannot agree by saturating.
    srt <- sort(a_rev)
    n <- length(srt)
    idx <- seq.int(floor(0.2 * n), floor(0.8 * n))
    j <- idx[which.max(srt[idx + 1L] - srt[idx])]
    th <- (srt[j] + srt[j + 1L]) / 2
    expect_gt(srt[j + 1L] - srt[j], 1e-3)

    p <- c(m, list(bidirect = FALSE, strand = -1, extend = TRUE))
    gvtrack.create("r_lse", NULL, "potts", params = p)
    gvtrack.create("r_max", NULL, "potts.max", params = p)
    gvtrack.create("r_pos", NULL, "potts.max.pos", params = p)
    gvtrack.create("r_cnt", NULL, "potts.count", params = c(p, list(score.thresh = th)))
    gvtrack.create("r_all", NULL, "potts.count",
        params = c(p, list(score.thresh = min(a_rev) - 1))
    )
    gvtrack.create("r_none", NULL, "potts.count",
        params = c(p, list(score.thresh = max(a_rev) + 1))
    )

    got <- gextract(c("r_lse", "r_max", "r_pos", "r_cnt", "r_all", "r_none"),
        iv,
        iterator = iv
    )

    expect_equal(got$r_lse, log_sum_exp(a_rev), tolerance = 1e-5, ignore_attr = TRUE)
    expect_equal(got$r_max, max(a_rev), tolerance = 1e-5, ignore_attr = TRUE)
    expect_equal(got$r_cnt, n - j)
    expect_equal(got$r_cnt, sum(a_rev >= th))
    expect_equal(got$r_all, n)
    expect_equal(got$r_none, 0)

    # The position is reported in FORWARD-strand coordinates even though the
    # target was fetched reverse-complemented, so it indexes a_rev directly.
    # Asserted as the score at the reported anchor, not as an argmax index:
    # ties make the index itself undefined.
    pos <- got$r_pos
    expect_gt(pos, 0) # unsigned: the sign is only applied when bidirect
    expect_lte(pos, n)
    expect_equal(as.numeric(a_rev[pos]), max(a_rev), tolerance = 1e-5, ignore_attr = TRUE)
})

test_that("a potts vtrack honours gvtrack.iterator shifts", {
    remove_all_vtracks()
    m <- potts_ref_model(W = 6L, npair_mode = "full", seed = 113L)

    gvtrack.create("s_lse", NULL, "potts", params = c(m, list(extend = TRUE)))
    gvtrack.create("p_lse", NULL, "potts", params = c(m, list(extend = TRUE)))
    gvtrack.iterator("s_lse", sshift = -50, eshift = 50)

    narrow <- gintervals(1, 2050, 2150)
    wide <- gintervals(1, 2000, 2200)
    shifted <- gextract(c("s_lse", "p_lse"), narrow, iterator = narrow)
    plain <- gextract("p_lse", wide, iterator = wide)

    expect_equal(shifted$s_lse, plain$p_lse, tolerance = 1e-6)
    expect_false(isTRUE(all.equal(shifted$s_lse, shifted$p_lse)))
})

test_that("a bidirect potts vtrack ignores strand", {
    remove_all_vtracks()
    W <- 6L
    m <- potts_ref_model(W = W, npair_mode = "full", seed = 113L)

    all_chroms <- gintervals.all()
    chr21_end <- all_chroms$end[all_chroms$chrom == "chr21"]

    # Both intervals sit exactly where the two scan orders diverge: chr1 opens
    # with a TAACCC repeat, so many anchors tie to the last bit, and chr21's
    # tail has the end-only extension fully clipped by the chromosome end.
    #
    # `strand` cannot change WHICH anchors a bidirect scan covers - the range is
    # [0, tlen - W] in either orientation and maps onto the same genomic
    # windows - so potts, potts.max and potts.count agreed even before
    # .potts_params() clamped strand. potts.max.pos did not: strand = -1 fetches
    # the reverse-complemented target, walks the same anchors in the opposite
    # order, and first-wins tie-breaking then keeps a different, equally
    # maximal, one. Measured before the clamp: 2 vs 56 here, 9 vs 21 on the
    # chr21 tail, both pairs verified as exactly tied maximisers.
    ivs <- list(
        chr1_start = gintervals(1, 0, 60),
        chr21_tail = gintervals(21, chr21_end - 30, chr21_end)
    )
    funcs <- c("potts", "potts.max", "potts.max.pos", "potts.count")

    for (nm in names(ivs)) {
        iv <- ivs[[nm]]
        # A mid-range threshold, so potts.count cannot agree by saturating
        th <- stats::median(
            potts_ref_anchors(toupper(gseq.extract(iv)), m, bidirect = TRUE, union = "lse")
        )

        for (ext in c(TRUE, FALSE)) {
            remove_all_vtracks()
            base <- c(m, list(bidirect = TRUE, extend = ext, score.thresh = th))
            for (f in funcs) {
                gvtrack.create(paste0("p_", make.names(f)), NULL, f,
                    params = c(base, list(strand = 1))
                )
                gvtrack.create(paste0("m_", make.names(f)), NULL, f,
                    params = c(base, list(strand = -1))
                )
            }

            got <- gextract(
                c(paste0("p_", make.names(funcs)), paste0("m_", make.names(funcs))),
                iv,
                iterator = iv
            )
            for (f in funcs) {
                info <- paste(nm, "extend", ext, f)
                expect_identical(
                    got[[paste0("p_", make.names(f))]],
                    got[[paste0("m_", make.names(f))]],
                    info = info
                )
            }
            # not vacuous: a real count, and a real position
            expect_gt(got$p_potts.count, 0)
            expect_lt(got$p_potts.count, iv$end - iv$start)
            expect_gt(abs(got$p_potts.max.pos), 0)
        }
    }
})

# ---------------------------------------------------------------------------
# The sliding-window cache.
#
# Every test below compares the SLID path against the RE-SEEDED path: the same
# kernel with the same tie-breaking on both sides, so they must agree, and a
# disagreement is a cache bug and never a rounding artefact. The re-seeded side
# is one gextract() per position, each of which builds a fresh
# TrackExpressionVars and therefore a fresh, cold PottsScorer.
#
# Note what potts.max.pos returns: a 1-based index into the FETCHED target, not
# a genomic coordinate. A window whose argmax anchor does not move therefore
# reports a position one lower at every step, which is why RunningMaxDeque keys
# on absolute genomic positions and the cache has to map back through the
# current expanded.start.
# ---------------------------------------------------------------------------

# Compare a contiguous scan of `iv` against the same intervals fetched one at a
# time. `pad` is the gvtrack.iterator shift (NA for none) and `it` the iterator
# width, which together set the window size and the stride. `check_varies` asks
# that the column actually move, so a constant answer cannot satisfy the
# comparison by accident.
expect_potts_cache_agrees <- function(fn, params, pad, iv, info, it = 1,
                                      check_varies = TRUE) {
    remove_all_vtracks()
    gvtrack.create("t", NULL, fn, params = params)
    if (!is.na(pad)) {
        gvtrack.iterator("t", sshift = -pad, eshift = pad)
    }

    slid <- gextract("t", iv, iterator = it)
    # Fixed bins are aligned to absolute multiples of the bin size, so the span
    # can open and close with a partial bin. Asserted rather than assumed: a
    # scan that silently returned nothing would satisfy every comparison below.
    testthat::expect_equal(
        nrow(slid),
        floor((iv$end - 1) / it) - floor(iv$start / it) + 1,
        info = info
    )

    one_at_a_time <- vapply(rev(seq_len(nrow(slid))), function(k) {
        gextract("t", slid[k, 1:3], iterator = slid[k, 1:3])$t
    }, numeric(1))

    testthat::expect_equal(slid$t, rev(one_at_a_time), tolerance = 1e-6, info = info)
    if (check_varies) {
        testthat::expect_gt(length(unique(slid$t)), 1L)
    }
    slid$t
}

# A mid-range threshold, so potts.count cannot agree by saturating at 0 or at
# the window size.
potts_cache_thresh <- function(m, iv, pad, bidirect, strand) {
    W <- nrow(m$e)
    seq_ext <- toupper(gseq.extract(
        gintervals(iv$chrom, iv$start - pad, iv$end + pad + W - 1L)
    ))
    stats::median(potts_ref_anchors(seq_ext, m,
        bidirect = bidirect, strand = strand, union = "lse"
    ), na.rm = TRUE)
}

test_that("the sliding cache returns what the seeded path returns", {
    m <- potts_ref_model(W = 8L, npair_mode = "full", seed = 131L)
    iv <- gintervals(1, 5000, 5200)

    # Both orientations. The cache is keyed on strand_mode, and when the target
    # is reverse-complemented the anchor entering the window on a one-position
    # step arrives at the LOW target index, so the forward orientation alone
    # exercises only half of the slide arithmetic. bidirect = FALSE is the only
    # way in: .potts_params() clamps strand to 1 whenever bidirect is TRUE.
    orients <- list(
        forward = list(bidirect = TRUE, strand = 1),
        reverse = list(bidirect = FALSE, strand = -1)
    )

    # Four window/stride geometries:
    #   pad 250, iterator 1 - the regime the cache exists for, a 501-anchor
    #     window moving 1 bp. There the region's best anchor never leaves, so
    #     potts.max is constant and only the other three modes' comparisons say
    #     anything, hence the pass below.
    #   pad  40, iterator 1 - narrow enough that every mode's answer moves.
    #   pad  40, iterator 7 - stride 7 into an 81-anchor window, so the
    #     multi-anchor push loop runs with a stride well inside the window.
    #   pad  40, iterator 60 - stride 60 into a 141-anchor window, so most of
    #     the window is evicted and refilled on every step.
    #   no shift, iterator 10 - consecutive windows share no anchor at all, so
    #     the cache declines to keep one and every call answers directly. The
    #     comparison is trivially satisfied there; it is here to check that the
    #     un-populated path still returns the right number.
    geoms <- list(
        list(nm = "pad250/it1", pad = 250L, it = 1L, varies = FALSE),
        list(nm = "pad40/it1", pad = 40L, it = 1L, varies = TRUE),
        list(nm = "pad40/it7", pad = 40L, it = 7L, varies = TRUE),
        list(nm = "pad40/it60", pad = 40L, it = 60L, varies = TRUE),
        list(nm = "noshift/it10", pad = NA_integer_, it = 10L, varies = TRUE)
    )

    for (nm in names(orients)) {
        o <- orients[[nm]]
        for (g in geoms) {
            th <- potts_cache_thresh(
                m, iv, if (is.na(g$pad)) 0L else g$pad, o$bidirect, as.integer(o$strand)
            )
            for (fn in c("potts", "potts.max", "potts.max.pos", "potts.count")) {
                expect_potts_cache_agrees(
                    fn, c(m, o, list(extend = TRUE, score.thresh = th)),
                    g$pad, iv,
                    info = paste(nm, g$nm, fn),
                    it = g$it, check_varies = g$varies
                )
            }
        }
    }
})

test_that("the potts cache re-seeds where the chromosome end clips the fetch", {
    m <- potts_ref_model(W = 12L, npair_mode = "full", seed = 149L)
    pad <- 30L

    # The end-only extension gets progressively clipped by the chromosome end,
    # so the fetched target shortens and i_max moves while the interval keeps
    # its width. The slide guard has to notice and re-seed; if it slid anyway,
    # the reverse orientation would key its anchors off the wrong target length.
    chr21_end <- gintervals.all()$end[gintervals.all()$chrom == "chr21"]
    iv <- gintervals(21, chr21_end - 120, chr21_end)

    for (o in list(
        list(bidirect = TRUE, strand = 1),
        list(bidirect = FALSE, strand = -1)
    )) {
        for (fn in c("potts", "potts.max", "potts.max.pos", "potts.count")) {
            expect_potts_cache_agrees(
                fn, c(m, o, list(extend = TRUE, score.thresh = 0)),
                pad, iv,
                info = paste("chr21 tail", "bidirect", o$bidirect, fn),
                check_varies = FALSE
            )
        }
    }
})

test_that("the potts cache is invalidated on a chromosome change", {
    remove_all_vtracks()
    m <- potts_ref_model(W = 8L, npair_mode = "full", seed = 137L)
    gvtrack.create("t", NULL, "potts.max", params = c(m, list(extend = TRUE)))

    chroms <- gintervals.all()
    skip_if(nrow(chroms) < 2, "test database has one chromosome")

    iv <- rbind(
        gintervals(chroms$chrom[1], 3000, 3020),
        gintervals(chroms$chrom[2], 3000, 3020)
    )

    # Interleaved across the boundary, versus one chromosome at a time. A cache
    # that survives start_chrom() returns the first chromosome's window here.
    both <- gextract("t", iv, iterator = 1)
    sep <- rbind(
        gextract("t", iv[1, ], iterator = 1),
        gextract("t", iv[2, ], iterator = 1)
    )
    expect_equal(both$t, sep$t, tolerance = 1e-6)
})

# A model under which EVERY scorable anchor scores exactly the same: all-zero
# `e` and `J`, zero intercept. Nothing else makes a tie-break testable. With a
# random model an exact tie is a fixture accident, so a comparison that only
# ever sees distinct scores says nothing about which of two equal maximisers a
# path picks - and picking differently is precisely how the cached and the
# uncached paths can disagree while both look right.
potts_tied_model <- function(W) {
    full <- t(utils::combn(W, 2))
    storage.mode(full) <- "integer"
    dimnames(full) <- NULL
    list(
        e = matrix(0, nrow = W, ncol = 4L, dimnames = list(NULL, POTTS_BASES)),
        J = lapply(seq_len(nrow(full)), function(k) matrix(0, 4L, 4L)),
        pairs = full,
        intercept = 0
    )
}

test_that("the potts cache breaks a tie the same way whatever the scan's shape", {
    m <- potts_tied_model(8L)
    iv <- gintervals(1, 5000, 5200)

    # potts.max.pos is the only mode a tie-break can reach - the other three
    # reduce to a number that does not depend on WHICH maximiser is named - but
    # all four are compared, because a tied model is also the cheapest check
    # that nothing else in the cache depends on the scores being distinct.
    #
    # The geometries matter more than usual here. `noshift` is the case where
    # consecutive windows share no anchor, so the cache declines to keep one and
    # every call after the first answers from score_direct() - a SECOND code
    # path for the same question, and therefore a second chance to break a tie
    # differently from the deque. The one-interval-at-a-time reference always
    # takes the first call of a fresh scorer, so the two only agree if both
    # paths use one rule.
    geoms <- list(
        list(nm = "pad40/it1", pad = 40L, it = 1L),
        list(nm = "pad40/it7", pad = 40L, it = 7L),
        list(nm = "noshift/it10", pad = NA_integer_, it = 10L),
        list(nm = "noshift/it1", pad = NA_integer_, it = 1L)
    )
    orients <- list(
        forward = list(bidirect = TRUE, strand = 1),
        forward_uni = list(bidirect = FALSE, strand = 1),
        reverse = list(bidirect = FALSE, strand = -1)
    )

    for (nm in names(orients)) {
        for (g in geoms) {
            for (fn in c("potts", "potts.max", "potts.max.pos", "potts.count")) {
                expect_potts_cache_agrees(
                    fn, c(m, orients[[nm]], list(extend = TRUE, score.thresh = -1)),
                    g$pad, iv,
                    info = paste("tied", nm, g$nm, fn),
                    it = g$it, check_varies = FALSE
                )
            }
        }
    }
})

# The 0-based position where the leading N run of `chrom` ends, or NA.
potts_n_run_end <- function(chrom, limit = 3e5) {
    v <- strsplit(toupper(gseq.extract(gintervals(chrom, 0, limit))), "", fixed = TRUE)[[1L]]
    i <- which(v != "N")
    if (!length(i) || i[1L] == 1L) NA_integer_ else as.integer(i[1L] - 1L)
}

# The window shift that puts an N -> sequence boundary INSIDE one incoming
# batch, with an unscorable anchor pushed before a scorable one. NA if no shift
# in `pads` does.
#
# Geometry: the anchors of iterator bin [s, s + it) are
# [s - pad, s + it - 1 + pad]; misha aligns fixed bins to multiples of the bin
# size; an anchor is scorable only if it starts at or after `bnd`; and a slide
# pushes the window's top `it` anchors in ascending genomic order. So with T the
# top anchor of the arriving bin the case needs
#   T >= bnd            - the batch brings in a scorable anchor,
#   T - it < bnd        - the window held nothing scorable before the batch,
#   T - it + 1 < bnd    - and the batch pushes -Inf before that anchor,
# i.e. bnd <= T <= bnd + it - 2. That needs it >= 2 - at it == 2 the range is
# the single point T == bnd, where the batch is {bnd - 1, bnd}: an unscorable
# push followed by a scorable one, which is the whole requirement - and the
# right phase. An iterator = 1 scan cannot reach it however far it runs, because
# its batch is one anchor and so cannot push -Inf before a real value. The
# geometry below picks it = 7, which is stricter than necessary but leaves room
# for the phase search to find a shift.
potts_straddling_pad <- function(bnd, it, pads = 20:80) {
    for (pad in pads) {
        top <- seq(0L, 2L * bnd, by = it) - 1L + pad
        if (any(top >= bnd & top <= bnd + it - 2L)) {
            return(as.integer(pad))
        }
    }
    NA_integer_
}

test_that("a warm potts cache slides into and out of an all-N stretch", {
    m <- potts_ref_model(W = 8L, npair_mode = "full", seed = 139L)
    funcs <- c("potts", "potts.max", "potts.max.pos", "potts.count")

    # The one place the cached and the uncached paths can disagree SILENTLY. An
    # unscorable anchor enters the running structures as -Inf, the identity for
    # both a log-sum-exp and a maximum; when the window's last finite anchor is
    # evicted the aggregator goes to -Inf while the cache stays valid, and -Inf
    # is not the answer - NaN is, or 0 for potts.count. -Inf survives arithmetic
    # instead of poisoning it, so it is the harder failure to notice.
    #
    # Both directions, because they fail differently: sliding IN evicts the
    # window's last finite anchor, sliding OUT has to recover from a window that
    # is entirely -Inf.
    bnd20 <- potts_n_run_end(20)
    skip_if(is.na(bnd20), "chr20 does not open with an N run in this fixture")
    spans <- list(
        into_N = gintervals(1, 167240, 167320), # chr1 sequence ends at 167280
        out_of_N = gintervals(20, bnd20 - 40, bnd20 + 40)
    )
    for (nm in names(spans)) {
        iv <- spans[[nm]]
        sq <- toupper(gseq.extract(gintervals(
            iv$chrom, iv$start - 20L, iv$end + 20L + nrow(m$e) - 1L
        )))
        skip_if(!grepl("N", sq, fixed = TRUE), paste("no N run near", nm))
        skip_if(!grepl("[ACGT]", sq), paste("no clean sequence near", nm))

        for (fn in funcs) {
            v <- expect_potts_cache_agrees(
                fn, c(m, list(extend = TRUE, score.thresh = 0)),
                20L, iv,
                info = paste(nm, fn), it = 1L
            )
            # The span really does cross the boundary, and the answer where no
            # anchor is scorable is NaN - 0 for the count - never -Inf.
            info <- paste(nm, fn)
            if (fn == "potts.count") {
                expect_true(any(v == 0), info = info)
                expect_true(any(v > 0), info = info)
            } else {
                expect_true(any(is.na(v)), info = info)
                expect_true(any(is.finite(v)), info = info)
                expect_false(any(is.infinite(v)), info = info)
            }
        }
    }
})

test_that("a potts cache batch that straddles the end of an N run agrees too", {
    m <- potts_ref_model(W = 8L, npair_mode = "full", seed = 151L)
    it <- 7L

    # An iterator = 1 scan across the same boundary cannot reach this: the
    # arriving batch is one anchor, so it never pushes -Inf and then a real
    # value into the SAME batch while the window holds nothing scorable. That
    # ordering is what reaches RunningLogSumExp::push(-INFINITY) with its
    # running maximum already -inf, which computes exp(-inf + inf) = NaN into
    # the accumulator - and the NaN then survives into the first finite push,
    # so the interval reports NaN although it has scorable anchors. Measured at
    # 34 of 58 scanned positions before this was guarded.
    bnd <- potts_n_run_end(20)
    skip_if(is.na(bnd), "chr20 does not open with an N run in this fixture")
    pad <- potts_straddling_pad(bnd, it)
    skip_if(is.na(pad), "no window shift straddles this fixture's N boundary")

    iv <- gintervals(20, bnd - 4L * it - pad, bnd + 4L * it)
    for (fn in c("potts", "potts.max", "potts.max.pos", "potts.count")) {
        v <- expect_potts_cache_agrees(
            fn, c(m, list(extend = TRUE, score.thresh = 0)),
            pad, iv,
            info = paste("straddle", fn), it = it, check_varies = FALSE
        )
        info <- paste("straddle", fn)
        if (fn == "potts.count") {
            expect_true(any(v == 0), info = info)
            expect_true(any(v > 0), info = info)
        } else {
            expect_true(any(is.na(v)), info = info)
            expect_true(any(is.finite(v)), info = info)
        }
    }
})

test_that("potts vtracks aggregate across a filter's unmasked parts", {
    remove_all_vtracks()
    m <- potts_ref_model(W = 6L, npair_mode = "full", seed = 139L)

    iv <- gintervals(1, 4000, 4300)
    mask <- gintervals(1, 4100, 4200) # punches the middle out

    th <- 0
    for (fn in c("potts", "potts.max", "potts.count")) {
        remove_all_vtracks()
        gvtrack.create("t", NULL, fn,
            params = c(m, list(extend = TRUE, score.thresh = th))
        )
        gvtrack.filter("t", filter = mask)
        # scores each unmasked part on its own, with the same func and params
        # and no filter - the reference the filtered result is aggregated from
        gvtrack.create("t_part", NULL, fn,
            params = c(m, list(extend = TRUE, score.thresh = th))
        )
        got <- gextract("t", iv, iterator = iv)$t

        # the two unmasked parts, scored separately
        parts <- c(
            gextract("t_part", gintervals(1, 4000, 4100), iterator = gintervals(1, 4000, 4100))$t_part,
            gextract("t_part", gintervals(1, 4200, 4300), iterator = gintervals(1, 4200, 4300))$t_part
        )

        expected <- switch(fn,
            potts = log_sum_exp(parts),
            potts.max = max(parts),
            potts.count = sum(parts)
        )
        expect_equal(got, expected, tolerance = 1e-5, info = fn)
    }

    # a fully masked interval scores nothing - NaN for all four, potts.count
    # included. See the next test for why the count is not 0 here.
    for (fn in c("potts", "potts.max", "potts.max.pos", "potts.count")) {
        remove_all_vtracks()
        gvtrack.create("t", NULL, fn, params = c(m, list(extend = TRUE, score.thresh = th)))
        gvtrack.filter("t", filter = iv)
        got <- gextract("t", iv, iterator = iv)$t
        expect_true(is.na(got), info = fn)
    }
})

test_that("a filtered potts.count separates 'counted nothing' from 'nothing to count'", {
    # The count has two answers that look alike in a gextract() column and are
    # not the same claim:
    #
    #   the filter leaves parts, those parts hold anchors of the model, and
    #     none of the anchors is scorable - an assembly gap. The scan ran over
    #     a real set of candidates and none passed, so the count is 0.
    #   the filter leaves no part at all, or leaves only parts narrower than
    #     the model, so no window of it can be placed anywhere. There was no
    #     candidate to pass or fail. 0 would be indistinguishable from the row
    #     above, which folds "this bin was excluded" into "this bin had no
    #     hits" under any sum(), mean() or threshold downstream, and hides the
    #     excluded bins from is.na(). NA is the value R keeps for that, and it
    #     is what the other three funcs already report in the same situation.
    #
    # The third block is what stops this test from being satisfied by "NA
    # everywhere": it pins that a genuine count of zero is still 0.
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    W <- 6L
    m <- potts_ref_model(W = W, npair_mode = "full", seed = 173L)
    funcs <- c("potts", "potts.max", "potts.max.pos", "potts.count")

    iv <- gintervals(1, 4000, 4300)
    narrow_mask <- gintervals(1, 4003, 4297) # leaves 3 bp at each end, < W

    # chr20 opens with an assembly gap, so both surviving parts here hold
    # anchors and not one of them is scorable
    gap <- gintervals(20, 0, 300)
    skip_if_not(
        grepl("^N+$", toupper(gseq.extract(gintervals(20, 0, 300L + W - 1L)))),
        "no all-N interval at the start of chr20 in this fixture"
    )
    gap_mask <- gintervals(20, 100, 200)

    for (fn in funcs) {
        # no part at all
        remove_all_vtracks()
        gvtrack.create("t", NULL, fn, params = c(m, list(extend = TRUE, score.thresh = 0)))
        gvtrack.filter("t", filter = iv)
        expect_true(is.na(gextract("t", iv, iterator = iv)$t),
            info = paste("fully masked", fn)
        )

        # parts exist, none of them wide enough to hold a window. extend =
        # FALSE, or the fetch would be padded to the model's width and the
        # parts would score.
        remove_all_vtracks()
        gvtrack.create("t", NULL, fn, params = c(m, list(extend = FALSE, score.thresh = 0)))
        gvtrack.filter("t", filter = narrow_mask)
        expect_true(is.na(gextract("t", iv, iterator = iv)$t),
            info = paste("no anchor in any part", fn)
        )

        # parts exist and hold anchors, but the sequence is an assembly gap:
        # a real count of zero, which must stay 0
        remove_all_vtracks()
        gvtrack.create("t", NULL, fn, params = c(m, list(extend = TRUE, score.thresh = 0)))
        gvtrack.filter("t", filter = gap_mask)
        got <- gextract("t", gap, iterator = gap)$t
        if (fn == "potts.count") {
            expect_equal(got, 0, info = paste("assembly gap", fn))
        } else {
            expect_true(is.na(got), info = paste("assembly gap", fn))
        }
    }
})

test_that("potts.max.pos aggregates across a filter's unmasked parts by score, not index", {
    remove_all_vtracks()
    W <- 6L
    m <- potts_ref_model(W = W, npair_mode = "full", seed = 139L)

    iv <- gintervals(1, 4000, 4300)
    mask <- gintervals(1, 4100, 4200) # same punched-out middle as above

    seq_ext <- toupper(gseq.extract(gintervals(1, iv$start, iv$end + W - 1L)))
    a_max <- potts_ref_anchors(seq_ext, m, bidirect = TRUE, union = "max")

    # 1-based anchor indices (relative to iv$start) that the filter leaves
    # unmasked - the two parts gvtrack.filter's complement produces, expressed
    # as index ranges into a_max rather than as separate gextract calls, so
    # this does not depend on the C++ aggregation's part-vs-part tie-break.
    unmasked_idx <- c(
        seq_len(mask$start - iv$start),
        seq.int(mask$end - iv$start + 1L, length(a_max))
    )
    expected_max <- max(a_max[unmasked_idx])

    gvtrack.create("t_pos", NULL, "potts.max.pos", params = c(m, list(extend = TRUE)))
    gvtrack.filter("t_pos", filter = mask)
    got <- gextract("t_pos", iv, iterator = iv)$t_pos

    # Assert the SCORE at the reported anchor, and that the anchor sits in an
    # unmasked part - never the argmax index. Under bidirect a window and its
    # reverse complement give the same max-union, so two anchors (here,
    # potentially one per part) can tie exactly, and either is a correct
    # answer - see "potts.max.pos and potts.count work on a genome interval"
    # above, which ties this same window/revcomp degeneracy to the exact
    # values involved.
    p <- abs(got)
    expect_true(p %in% unmasked_idx, info = paste("p =", p))
    expect_equal(as.numeric(a_max[p]), expected_max, tolerance = 1e-5, ignore_attr = TRUE)

    win <- substr(seq_ext, p, p + W - 1L)
    f <- as.numeric(potts_ref_window(win, m))
    r <- as.numeric(potts_ref_window(win, potts_ref_rc(m)))
    if (abs(f - r) > 1e-9) {
        expect_equal(sign(got), if (r > f) -1 else 1)
    }
})

test_that("a single unmasked part is scored like the unfiltered survivor, offset when it is a position", {
    # A mask flush with one edge of the iterator interval leaves exactly one
    # unmasked part rather than two - the only geometry that exercises the
    # aggregation's one-part case, which the earlier tests in this file never
    # reached (they all mask out the MIDDLE, leaving two parts). potts,
    # potts.max and potts.count don't care: a lone part is just scored
    # directly. potts.max.pos does care: its answer is a position, and a mask
    # that clips the START of the interval moves the surviving part's own
    # start away from the interval's start, so the position the part reports
    # (relative to itself) has to be offset before it means anything relative
    # to the original interval - exactly the offset the two-or-more-parts case
    # already applies. A mask that clips the END leaves the surviving part
    # flush with the interval's start (offset zero), so that side alone could
    # never have caught a missing offset.
    remove_all_vtracks()
    m <- potts_ref_model(W = 6L, npair_mode = "full", seed = 139L)
    iv <- gintervals(1, 4000, 4300)
    th <- 0

    clips <- list(
        start = list(mask = gintervals(1, 4000, 4100), part = gintervals(1, 4100, 4300)),
        end = list(mask = gintervals(1, 4200, 4300), part = gintervals(1, 4000, 4200))
    )

    for (nm in names(clips)) {
        mask <- clips[[nm]]$mask
        part <- clips[[nm]]$part

        for (fn in c("potts", "potts.max", "potts.count")) {
            remove_all_vtracks()
            gvtrack.create("t", NULL, fn, params = c(m, list(extend = TRUE, score.thresh = th)))
            gvtrack.filter("t", filter = mask)
            gvtrack.create("t_part", NULL, fn, params = c(m, list(extend = TRUE, score.thresh = th)))

            got <- gextract("t", iv, iterator = iv)$t
            expected <- gextract("t_part", part, iterator = part)$t_part
            expect_equal(got, expected, tolerance = 1e-5, info = paste(nm, fn))
        }

        remove_all_vtracks()
        gvtrack.create("t_pos", NULL, "potts.max.pos", params = c(m, list(extend = TRUE)))
        gvtrack.filter("t_pos", filter = mask)
        gvtrack.create("p_pos", NULL, "potts.max.pos", params = c(m, list(extend = TRUE)))

        got <- gextract("t_pos", iv, iterator = iv)$t_pos
        local_pos <- gextract("p_pos", part, iterator = part)$p_pos
        offset <- part$start - iv$start
        expected <- if (local_pos > 0) {
            local_pos + offset
        } else if (local_pos < 0) {
            local_pos - offset
        } else {
            local_pos
        }
        expect_equal(got, expected, info = paste(nm, "potts.max.pos value"))

        # The assertion that actually would have caught the bug this test
        # guards against: a position that survives the filter must resolve to
        # a genomic coordinate OUTSIDE the mask. The old one-part fast path
        # returned the part-relative position unadjusted, which - for the
        # "start" clip - decodes to a coordinate INSIDE the masked region
        # (verified: iv = [4000,4300), mask = [4000,4100) reported -65,
        # decoding to iv$start + 65 - 1 = 4064, inside [4000,4100); the fixed
        # code reports -165, decoding to 4164, outside it).
        genomic_anchor_start <- iv$start + abs(got) - 1L
        outside_mask <- genomic_anchor_start < mask$start || genomic_anchor_start >= mask$end
        expect_true(outside_mask, info = paste(
            nm, "anchor at", genomic_anchor_start, "mask", mask$start, mask$end
        ))
    }
})

test_that("a potts interval with no anchor at all is NaN, potts.count included", {
    # Two situations that look alike and are not:
    #
    #   anchors exist, none of them scorable (an all-N assembly gap) - the scan
    #     ran and found nothing, so potts.count is 0 and the other three NaN;
    #   no anchor exists at all - nothing to count, so all four are NaN.
    #
    # The second row is reached at every chromosome end with the DEFAULT
    # extend = TRUE, not only through the narrow-interval-with-extend = FALSE
    # route the docs used to name: extend pads the fetch by W - 1, that padding
    # is clipped at the contig boundary, and the target then comes back shorter
    # than the model. potts.count answered 0 there, so every chromosome edge
    # read as a real count of zero.
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    W <- 6L
    m <- potts_ref_model(W = W, npair_mode = "full", seed = 163L)
    funcs <- c("potts", "potts.max", "potts.max.pos", "potts.count")

    chrom_end <- gintervals.all()
    chrom_end <- chrom_end$end[chrom_end$chrom == "chr20"]

    # 3 bp, so narrower than the model. At the chromosome end extend has
    # nowhere to pad and the fetch is short; 1 kb inside, the same 3 bp gets
    # its full W - 1 of padding and scores. The contrast is what shows the
    # answer comes from the clipped fetch and not merely from a narrow
    # interval.
    edge <- gintervals(20, chrom_end - 3L, chrom_end)
    inner <- gintervals(20, chrom_end - 1000L, chrom_end - 997L)
    skip_if_not(
        grepl("^[ACGT]+$", toupper(gseq.extract(gintervals(20, chrom_end - 1000L, chrom_end)))),
        "chr20 does not end in clean sequence in this fixture"
    )

    # all-N: chr20 opens with an assembly gap
    gap <- gintervals(20, 0, 200)
    skip_if_not(
        grepl("^N+$", toupper(gseq.extract(gap))),
        "no all-N interval at the start of chr20 in this fixture"
    )

    for (fn in funcs) {
        remove_all_vtracks()
        gvtrack.create("t", NULL, fn, params = c(m, list(extend = TRUE, score.thresh = 0)))

        # the 3 bp interval is scorable away from the boundary, so what makes
        # the edge unscorable is the clipped fetch
        expect_false(is.na(gextract("t", inner, iterator = inner)$t), info = paste("inner", fn))

        # no anchor at all -> NaN for all four
        expect_true(is.na(gextract("t", edge, iterator = edge)$t), info = paste("edge", fn))

        # anchors exist, none scorable -> 0 for the count, NaN for the rest
        got_gap <- gextract("t", gap, iterator = gap)$t
        if (fn == "potts.count") {
            expect_equal(got_gap, 0, info = paste("gap", fn))
        } else {
            expect_true(is.na(got_gap), info = paste("gap", fn))
        }
    }
})
