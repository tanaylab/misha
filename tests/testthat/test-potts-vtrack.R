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

test_that("a filtered potts vtrack refuses loudly", {
    remove_all_vtracks()
    m <- potts_ref_model(W = 6L, npair_mode = "sparse", seed = 131L)

    gvtrack.create("t_flt", NULL, "potts", params = c(m, list(extend = TRUE)))
    gvtrack.filter("t_flt", filter = gintervals(1, 2050, 2100))

    expect_error(
        gextract("t_flt", gintervals(1, 2000, 2200), iterator = 50),
        "not yet supported"
    )
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
