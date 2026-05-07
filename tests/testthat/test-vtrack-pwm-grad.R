test_that("oracle: flat PSSM gives g = 0 (lse, max, ism)", {
    pssm <- matrix(0.25, 3, 4, dimnames = list(NULL, c("A", "C", "G", "T")))
    expect_equal(manual_pwm_grad(pssm, "ACGTAA", "lse", bidirect = TRUE), 0)
    expect_equal(manual_pwm_grad(pssm, "ACGTAA", "max", bidirect = TRUE), 0)
    expect_equal(manual_pwm_grad_ism(pssm, "ACGTAA", "lse", bidirect = TRUE), 0)
})

test_that("oracle: L=1 PSSM, single strand, g_MAX = M[1,b_p] - min(M[1,])", {
    pssm <- matrix(c(0.7, 0.1, 0.1, 0.1), 1, 4,
        byrow = TRUE,
        dimnames = list(NULL, c("A", "C", "G", "T"))
    )
    M <- log((pssm + 0.01) / sum(pssm + 0.01))
    expected <- unname(M[1, 1] - min(M[1, ]))
    expect_equal(
        manual_pwm_grad(pssm, "ACGT", "max", bidirect = FALSE, strand = 1L),
        expected,
        tolerance = 1e-10
    )
})

test_that("oracle: linearized matches ISM when argmax is stable under flips (single strand)", {
    pssm <- matrix(
        c(
            0.97, 0.01, 0.01, 0.01,
            0.01, 0.97, 0.01, 0.01,
            0.01, 0.01, 0.97, 0.01
        ),
        3, 4,
        byrow = TRUE,
        dimnames = list(NULL, c("A", "C", "G", "T"))
    )
    seq <- "ACGAAAAA" # head anchor "ACG" perfect, others weak.
    g_lin <- manual_pwm_grad(pssm, seq, "max", bidirect = FALSE, strand = 1L)
    g_ism <- manual_pwm_grad_ism(pssm, seq, "max", bidirect = FALSE, strand = 1L)
    expect_equal(g_lin, g_ism, tolerance = 1e-8)
})

test_that("oracle: linearized lse with concentrated softmax matches max", {
    # Longer PSSM + uniform A-tail in seq pushes the non-head anchor scores
    # so far below the head anchor that the softmax concentrates at >= 1-1e-4
    # on the head, making g_LSE numerically indistinguishable from g_MAX.
    pssm <- matrix(
        c(
            0.97, 0.01, 0.01, 0.01,
            0.01, 0.97, 0.01, 0.01,
            0.01, 0.01, 0.97, 0.01,
            0.01, 0.01, 0.97, 0.01,
            0.01, 0.01, 0.97, 0.01
        ),
        5, 4,
        byrow = TRUE,
        dimnames = list(NULL, c("A", "C", "G", "T"))
    )
    seq <- "ACGGGAAAAAAAA" # head anchor "ACGGG" perfect, others weak.
    g_max <- manual_pwm_grad(pssm, seq, "max", bidirect = FALSE, strand = 1L)
    g_lse <- manual_pwm_grad(pssm, seq, "lse", bidirect = FALSE, strand = 1L)
    expect_equal(g_max, g_lse, tolerance = 1e-3)
})

test_that("oracle: bidirect on palindromic PSSM matches single-strand on palindrome", {
    pssm <- matrix(
        c(
            0.7, 0.1, 0.1, 0.1,
            0.1, 0.4, 0.4, 0.1,
            0.1, 0.1, 0.1, 0.7
        ),
        3, 4,
        byrow = TRUE,
        dimnames = list(NULL, c("A", "C", "G", "T"))
    )
    rc <- pssm[3:1, c(4, 3, 2, 1)]
    dimnames(rc) <- dimnames(pssm)
    expect_equal(rc, pssm)
    g <- manual_pwm_grad(pssm, "ATCAT", "lse", bidirect = TRUE)
    expect_true(is.finite(g) && g >= 0)
})

test_that("oracle: rc by RC-of-seq equals rc by RC-of-PSSM", {
    pssm <- matrix(
        c(
            0.7, 0.1, 0.1, 0.1,
            0.1, 0.6, 0.2, 0.1,
            0.1, 0.1, 0.6, 0.2
        ),
        3, 4,
        byrow = TRUE,
        dimnames = list(NULL, c("A", "C", "G", "T"))
    )
    seq <- "ACGCAT"
    g_strand_minus <- manual_pwm_grad(pssm, seq, "max", bidirect = FALSE, strand = -1L)
    expect_true(is.finite(g_strand_minus) && g_strand_minus >= 0)
})

test_that("oracle: rc-of-seq vs rc-of-PSSM give same per-anchor scores", {
    # Self-consistency: scoring fwd PSSM on RC(seq) at anchor a should equal
    # scoring rc PSSM on seq at the mirrored anchor (n - L - a + 1).
    pssm <- matrix(
        c(
            0.6, 0.2, 0.1, 0.1,
            0.1, 0.5, 0.3, 0.1,
            0.2, 0.1, 0.5, 0.2
        ),
        3, 4,
        byrow = TRUE,
        dimnames = list(NULL, c("A", "C", "G", "T"))
    )
    seq <- "ACGTACGTAC"
    L <- nrow(pssm)
    n <- nchar(seq)

    sc <- .oracle_anchor_scores(seq, pssm, prior = 0.01)

    # Independent path: rc score at anchor a == fwd score on RC(seq) at
    # position (n - L - (a - 1) + 1) - 1 + 1 = n - L - a + 2 (1-based).
    rc_seq <- .rc_seq(seq)
    fwd_on_rc <- manual_pwm_scores_single_strand(rc_seq, pssm, prior = 0.01)

    n_anchors <- n - L + 1L
    for (a in seq_len(n_anchors)) {
        mirror <- n_anchors - a + 1L
        expect_equal(sc$rc[a], fwd_on_rc[mirror], tolerance = 1e-12)
    }
})

test_that("oracle: ISM is always >= 0", {
    pssm <- matrix(
        c(
            0.7, 0.1, 0.1, 0.1,
            0.1, 0.6, 0.2, 0.1,
            0.1, 0.1, 0.6, 0.2
        ),
        3, 4,
        byrow = TRUE,
        dimnames = list(NULL, c("A", "C", "G", "T"))
    )
    for (seq in c("ACGCAT", "TTTTTT", "ACGACG", "GGGGGG")) {
        for (agg in c("lse", "max")) {
            g <- manual_pwm_grad_ism(pssm, seq, agg, bidirect = TRUE)
            expect_gte(g, -1e-10)
        }
    }
})

# ---- Engine tests for GRAD_MAX (Task 3) ----

create_isolated_test_db()

# Build a strong-consensus PSSM (0.97/0.01/0.01/0.01) whose column j has its
# peak at the j-th base of `bases`. `bases` must be a length-L character vector
# of A/C/G/T (no Ns).
.consensus_pssm <- function(bases) {
    L <- length(bases)
    stopifnot(all(bases %in% c("A", "C", "G", "T")))
    pssm <- matrix(0.01,
        nrow = L, ncol = 4,
        dimnames = list(NULL, c("A", "C", "G", "T"))
    )
    for (j in seq_len(L)) {
        idx <- match(bases[j], c("A", "C", "G", "T"))
        pssm[j, idx] <- 0.97
    }
    pssm
}

test_that("pwm.grad max: head=argmax => engine matches oracle col-0 diff", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    L <- 3L
    iter_W <- 20L # number of iterator-interval positions per pivot
    # Use a 1bp iterator interval; the iterator extension makes the vtrack
    # see a window of length iter_W (sshift=0, eshift=iter_W-1).
    pivot <- gintervals(1, 200, 201)
    # The engine fetches sequence over [pivot.start + sshift,
    #                                   pivot.end + eshift + (L-1))
    # because extend=TRUE adds (L-1) at the end. That's [200, 200 + iter_W + L - 1).
    seq_ext <- toupper(gseq.extract(gintervals(1, 200, 200 + iter_W + L - 1)))

    # Pick PSSM consensus = first L bases of seq_ext, so the head anchor wins.
    head_bases <- strsplit(substr(seq_ext, 1, L), "")[[1]]
    pssm <- .consensus_pssm(head_bases)

    gvtrack.create("g", NULL, "pwm.grad",
        pssm = pssm,
        aggregate = "max", bidirect = FALSE, strand = 1,
        extend = TRUE, prior = 0.01
    )
    gvtrack.iterator("g", sshift = 0, eshift = iter_W - 1)

    res <- gextract("g", iterator = 1, intervals = pivot)

    expected <- manual_pwm_grad(pssm, seq_ext, "max",
        bidirect = FALSE, strand = 1L, prior = 0.01
    )
    expect_equal(res$g[1], expected, tolerance = 1e-5, ignore_attr = TRUE)
    # Sanity: the head-wins case must produce a strictly positive gradient
    # (the head base is the consensus and column 0 is non-uniform).
    expect_gt(res$g[1], 0)
})

test_that("pwm.grad max: head not argmax => engine returns 0 (matches oracle)", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    L <- 3L
    iter_W <- 30L
    offset <- 10L # interior anchor index whose match wins
    pivot <- gintervals(1, 200, 201)
    seq_ext <- toupper(gseq.extract(gintervals(1, 200, 200 + iter_W + L - 1)))

    # Construct a PSSM whose consensus matches the interior anchor at `offset`,
    # NOT the head anchor at 0. Argmax should land at the interior, so head
    # gradient is 0 (provided the head bases are not coincidentally also a
    # consensus run, which we verify below).
    interior_bases <- strsplit(substr(seq_ext, offset + 1, offset + L), "")[[1]]
    pssm <- .consensus_pssm(interior_bases)

    # Guard the test setup: interior pattern must differ from head bases.
    head_bases <- strsplit(substr(seq_ext, 1, L), "")[[1]]
    skip_if(
        identical(interior_bases, head_bases),
        "test fixture: interior PSSM consensus collides with head bases; pick a different interval"
    )

    gvtrack.create("g", NULL, "pwm.grad",
        pssm = pssm,
        aggregate = "max", bidirect = FALSE, strand = 1,
        extend = TRUE, prior = 0.01
    )
    gvtrack.iterator("g", sshift = 0, eshift = iter_W - 1)

    res <- gextract("g", iterator = 1, intervals = pivot)

    expected <- manual_pwm_grad(pssm, seq_ext, "max",
        bidirect = FALSE, strand = 1L, prior = 0.01
    )
    expect_equal(res$g[1], expected, tolerance = 1e-5, ignore_attr = TRUE)
    # And it should be exactly 0 in this construction (argmax interior).
    expect_equal(expected, 0, tolerance = 1e-10)
})

test_that("pwm.grad max: bidirect=TRUE matches oracle on asymmetric PSSM", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    L <- 3L
    iter_W <- 20L
    pivot <- gintervals(1, 200, 201)
    seq_ext <- toupper(gseq.extract(gintervals(1, 200, 200 + iter_W + L - 1)))

    # Asymmetric PSSM (rc != self) so bidirect's strand-softmax is non-trivial.
    pssm <- matrix(
        c(
            0.7, 0.1, 0.1, 0.1,
            0.1, 0.6, 0.2, 0.1,
            0.1, 0.1, 0.5, 0.3
        ),
        L, 4,
        byrow = TRUE,
        dimnames = list(NULL, c("A", "C", "G", "T"))
    )

    gvtrack.create("g_bidir", NULL, "pwm.grad",
        pssm = pssm,
        aggregate = "max", bidirect = TRUE,
        extend = TRUE, prior = 0.01
    )
    gvtrack.iterator("g_bidir", sshift = 0, eshift = iter_W - 1)
    res <- gextract("g_bidir", iterator = 1, intervals = pivot)

    expected <- manual_pwm_grad(pssm, seq_ext, "max",
        bidirect = TRUE, prior = 0.01
    )
    expect_equal(res$g_bidir[1], expected, tolerance = 1e-5, ignore_attr = TRUE)
})

test_that("pwm.grad lse: bidirect=TRUE matches oracle on asymmetric PSSM", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    L <- 3L
    iter_W <- 20L
    pivot <- gintervals(1, 200, 201)
    seq_ext <- toupper(gseq.extract(gintervals(1, 200, 200 + iter_W + L - 1)))

    pssm <- matrix(
        c(
            0.7, 0.1, 0.1, 0.1,
            0.1, 0.6, 0.2, 0.1,
            0.1, 0.1, 0.5, 0.3
        ),
        L, 4,
        byrow = TRUE,
        dimnames = list(NULL, c("A", "C", "G", "T"))
    )

    gvtrack.create("g_bidir_lse", NULL, "pwm.grad",
        pssm = pssm,
        aggregate = "lse", bidirect = TRUE,
        extend = TRUE, prior = 0.01
    )
    gvtrack.iterator("g_bidir_lse", sshift = 0, eshift = iter_W - 1)
    res <- gextract("g_bidir_lse", iterator = 1, intervals = pivot)

    expected <- manual_pwm_grad(pssm, seq_ext, "lse",
        bidirect = TRUE, prior = 0.01
    )
    expect_equal(res$g_bidir_lse[1], expected, tolerance = 1e-5, ignore_attr = TRUE)
})

test_that("pwm.grad: bidirect on palindromic PSSM is consistent", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    # Palindromic: rc(pssm) == pssm. The bidirect comb score is fwd + log(2)
    # (offset, irrelevant for argmax). Strand softmax at any anchor is 0.5/0.5.
    pssm <- matrix(
        c(
            0.7, 0.1, 0.1, 0.1,
            0.1, 0.4, 0.4, 0.1,
            0.1, 0.1, 0.1, 0.7
        ),
        3, 4,
        byrow = TRUE,
        dimnames = list(NULL, c("A", "C", "G", "T"))
    )
    rc_test <- pssm[3:1, c(4, 3, 2, 1)]
    dimnames(rc_test) <- dimnames(pssm)
    expect_equal(rc_test, pssm) # confirm palindromic

    iter_W <- 20L
    pivot <- gintervals(1, 200, 201)
    seq_ext <- toupper(gseq.extract(gintervals(1, 200, 200 + iter_W + nrow(pssm) - 1)))

    gvtrack.create("g_bidir", NULL, "pwm.grad",
        pssm = pssm,
        aggregate = "lse", bidirect = TRUE,
        extend = TRUE, prior = 0.01
    )
    gvtrack.iterator("g_bidir", sshift = 0, eshift = iter_W - 1)
    res <- gextract("g_bidir", iterator = 1, intervals = pivot)

    expected <- manual_pwm_grad(pssm, seq_ext, "lse",
        bidirect = TRUE, prior = 0.01
    )
    expect_equal(res$g_bidir[1], expected, tolerance = 1e-5, ignore_attr = TRUE)
})

test_that("pwm.grad max: bidirect=FALSE strand=-1 matches oracle", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    L <- 3L
    iter_W <- 20L
    pivot <- gintervals(1, 200, 201)
    seq_ext <- toupper(gseq.extract(gintervals(1, 200, 200 + iter_W + L - 1)))

    # Build PSSM whose rc-consensus matches the first L bases of seq_ext, so
    # rc-PSSM scoring at the head anchor is strong.
    head_bases <- strsplit(substr(seq_ext, 1, L), "")[[1]]
    # rc-PSSM at column 0 sees complement(head[0]); we want a consensus PSSM
    # whose rc orientation has strong matches at our head bases. Simplest: use
    # the consensus PSSM for the actual head bases; the rc scoring will see
    # complement-of-head at L-1, L-2, ... so it will be a weak match. To make
    # rc scoring strong, use the rc of the consensus PSSM instead.
    fwd_pssm <- .consensus_pssm(head_bases)
    rc_target_pssm <- fwd_pssm[L:1, c(4, 3, 2, 1)]
    dimnames(rc_target_pssm) <- dimnames(fwd_pssm)

    gvtrack.create("g_minus", NULL, "pwm.grad",
        pssm = rc_target_pssm,
        aggregate = "max", bidirect = FALSE, strand = -1,
        extend = TRUE, prior = 0.01
    )
    gvtrack.iterator("g_minus", sshift = 0, eshift = iter_W - 1)
    res <- gextract("g_minus", iterator = 1, intervals = pivot)

    expected <- manual_pwm_grad(rc_target_pssm, seq_ext, "max",
        bidirect = FALSE, strand = -1L, prior = 0.01
    )
    expect_equal(res$g_minus[1], expected, tolerance = 1e-5, ignore_attr = TRUE)
})

test_that("pwm.grad lse: bidirect=FALSE strand=-1 matches oracle", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    L <- 3L
    iter_W <- 20L
    pivot <- gintervals(1, 200, 201)
    seq_ext <- toupper(gseq.extract(gintervals(1, 200, 200 + iter_W + L - 1)))

    head_bases <- strsplit(substr(seq_ext, 1, L), "")[[1]]
    fwd_pssm <- .consensus_pssm(head_bases)
    rc_target_pssm <- fwd_pssm[L:1, c(4, 3, 2, 1)]
    dimnames(rc_target_pssm) <- dimnames(fwd_pssm)

    gvtrack.create("g_minus_lse", NULL, "pwm.grad",
        pssm = rc_target_pssm,
        aggregate = "lse", bidirect = FALSE, strand = -1,
        extend = TRUE, prior = 0.01
    )
    gvtrack.iterator("g_minus_lse", sshift = 0, eshift = iter_W - 1)
    res <- gextract("g_minus_lse", iterator = 1, intervals = pivot)

    expected <- manual_pwm_grad(rc_target_pssm, seq_ext, "lse",
        bidirect = FALSE, strand = -1L, prior = 0.01
    )
    expect_equal(res$g_minus_lse[1], expected, tolerance = 1e-5, ignore_attr = TRUE)
})

# Task 4: GRAD_LSE single-strand fwd

test_that("pwm.grad lse: concentrated softmax (head=consensus) matches oracle", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    L <- 3L
    iter_W <- 20L
    pivot <- gintervals(1, 200, 201)
    seq_ext <- toupper(gseq.extract(gintervals(1, 200, 200 + iter_W + L - 1)))

    head_bases <- strsplit(substr(seq_ext, 1, L), "")[[1]]
    pssm <- .consensus_pssm(head_bases)

    gvtrack.create("g_lse", NULL, "pwm.grad",
        pssm = pssm,
        aggregate = "lse", bidirect = FALSE, strand = 1,
        extend = TRUE, prior = 0.01
    )
    gvtrack.iterator("g_lse", sshift = 0, eshift = iter_W - 1)
    res <- gextract("g_lse", iterator = 1, intervals = pivot)

    expected <- manual_pwm_grad(pssm, seq_ext, "lse",
        bidirect = FALSE, strand = 1L, prior = 0.01
    )
    expect_equal(res$g_lse[1], expected, tolerance = 1e-5, ignore_attr = TRUE)
    expect_gt(res$g_lse[1], 0)
})

test_that("pwm.grad lse: head not argmax => w_p small, g_lse < g_max", {
    # Same construction as the head-not-argmax test for max, but compare lse
    # vs max instead of just checking equality.
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    L <- 3L
    iter_W <- 30L
    offset <- 10L
    pivot <- gintervals(1, 200, 201)
    seq_ext <- toupper(gseq.extract(gintervals(1, 200, 200 + iter_W + L - 1)))

    interior_bases <- strsplit(substr(seq_ext, offset + 1, offset + L), "")[[1]]
    pssm <- .consensus_pssm(interior_bases)

    head_bases <- strsplit(substr(seq_ext, 1, L), "")[[1]]
    skip_if(
        identical(interior_bases, head_bases),
        "test fixture: interior PSSM consensus collides with head bases"
    )

    gvtrack.create("g_lse", NULL, "pwm.grad",
        pssm = pssm,
        aggregate = "lse", bidirect = FALSE, strand = 1,
        extend = TRUE, prior = 0.01
    )
    gvtrack.iterator("g_lse", sshift = 0, eshift = iter_W - 1)
    res <- gextract("g_lse", iterator = 1, intervals = pivot)

    expected <- manual_pwm_grad(pssm, seq_ext, "lse",
        bidirect = FALSE, strand = 1L, prior = 0.01
    )
    expect_equal(res$g_lse[1], expected, tolerance = 1e-5, ignore_attr = TRUE)
    # When the interior anchor is the LSE-dominant one, w_p (head) is small,
    # so g_lse stays close to 0 even though M[0, b_p] - worst is large.
    expect_lt(res$g_lse[1], 0.5)
})

test_that("pwm.grad lse: w_p computed correctly via softmax denominator", {
    # Hand-verifiable case: build a 2-anchor situation where head and one
    # interior anchor have known scores, then compare the engine's w_p * diff
    # against a direct computation.
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    L <- 3L
    iter_W <- 5L # exactly 5 anchors
    pivot <- gintervals(1, 200, 201)
    seq_ext <- toupper(gseq.extract(gintervals(1, 200, 200 + iter_W + L - 1)))

    head_bases <- strsplit(substr(seq_ext, 1, L), "")[[1]]
    pssm <- .consensus_pssm(head_bases)

    gvtrack.create("g_lse", NULL, "pwm.grad",
        pssm = pssm,
        aggregate = "lse", bidirect = FALSE, strand = 1,
        extend = TRUE, prior = 0.01
    )
    gvtrack.iterator("g_lse", sshift = 0, eshift = iter_W - 1)
    res <- gextract("g_lse", iterator = 1, intervals = pivot)

    expected <- manual_pwm_grad(pssm, seq_ext, "lse",
        bidirect = FALSE, strand = 1L, prior = 0.01
    )
    expect_equal(res$g_lse[1], expected, tolerance = 1e-5, ignore_attr = TRUE)

    # Cross-check: w_p must be in (0, 1] and g_lse <= g_max.
    g_max_oracle <- manual_pwm_grad(pssm, seq_ext, "max",
        bidirect = FALSE, strand = 1L, prior = 0.01
    )
    expect_lte(res$g_lse[1], g_max_oracle + 1e-8)
    expect_gt(res$g_lse[1], 0)
})

# Task 6: GRAD_MAX_ISM tests

test_that("pwm.grad.ism max: argmax-shift case differs from linearized", {
    # The test-DB seq at chr1:200-222 is periodic ("CCCTAACCC..."), so even
    # with a strong consensus PSSM matching the head, flipping the head base
    # doesn't drop the max because anchor 7 also matches the consensus.
    # Linearized would report ~3.9 (the col-0 diff); ISM correctly returns 0
    # because the argmax shifts under the flip. This is the canonical
    # ISM-vs-linearized divergence case.
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    L <- 3L
    iter_W <- 20L
    pivot <- gintervals(1, 200, 201)
    seq_ext <- toupper(gseq.extract(gintervals(1, 200, 200 + iter_W + L - 1)))

    head_bases <- strsplit(substr(seq_ext, 1, L), "")[[1]]
    pssm <- .consensus_pssm(head_bases)

    gvtrack.create("g_ism", NULL, "pwm.grad.ism",
        pssm = pssm, aggregate = "max",
        bidirect = FALSE, strand = 1, extend = TRUE, prior = 0.01
    )
    gvtrack.create("g_lin", NULL, "pwm.grad",
        pssm = pssm, aggregate = "max",
        bidirect = FALSE, strand = 1, extend = TRUE, prior = 0.01
    )
    gvtrack.iterator("g_ism", sshift = 0, eshift = iter_W - 1)
    gvtrack.iterator("g_lin", sshift = 0, eshift = iter_W - 1)
    res <- gextract(c("g_ism", "g_lin"), iterator = 1, intervals = pivot)

    expected_ism <- manual_pwm_grad_ism(pssm, seq_ext, "max",
        bidirect = FALSE, strand = 1L, prior = 0.01
    )
    expect_equal(res$g_ism[1], expected_ism, tolerance = 1e-5, ignore_attr = TRUE)
    # ISM = 0 (argmax shift); linearized > 0; they differ.
    expect_equal(expected_ism, 0, tolerance = 1e-10)
    expect_gt(res$g_lin[1], 1)
})

test_that("pwm.grad.ism max: matches oracle on weak-consensus PSSM", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    L <- 3L
    iter_W <- 15L
    pivot <- gintervals(1, 200, 201)
    seq_ext <- toupper(gseq.extract(gintervals(1, 200, 200 + iter_W + L - 1)))

    # Mid-strength PSSM so flipping can plausibly shift argmax.
    pssm <- matrix(
        c(
            0.5, 0.3, 0.1, 0.1,
            0.1, 0.5, 0.3, 0.1,
            0.1, 0.1, 0.5, 0.3
        ),
        L, 4,
        byrow = TRUE,
        dimnames = list(NULL, c("A", "C", "G", "T"))
    )

    gvtrack.create("g_ism", NULL, "pwm.grad.ism",
        pssm = pssm, aggregate = "max",
        bidirect = FALSE, strand = 1, extend = TRUE, prior = 0.01
    )
    gvtrack.iterator("g_ism", sshift = 0, eshift = iter_W - 1)
    res <- gextract("g_ism", iterator = 1, intervals = pivot)

    expected <- manual_pwm_grad_ism(pssm, seq_ext, "max",
        bidirect = FALSE, strand = 1L, prior = 0.01
    )
    expect_equal(res$g_ism[1], expected, tolerance = 1e-5, ignore_attr = TRUE)
    expect_gte(res$g_ism[1], -1e-7)
})

# Task 7: GRAD_LSE_ISM tests

test_that("pwm.grad.ism lse: matches oracle (concentrated softmax)", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    L <- 3L
    iter_W <- 20L
    pivot <- gintervals(1, 200, 201)
    seq_ext <- toupper(gseq.extract(gintervals(1, 200, 200 + iter_W + L - 1)))

    head_bases <- strsplit(substr(seq_ext, 1, L), "")[[1]]
    pssm <- .consensus_pssm(head_bases)

    gvtrack.create("g_ism_lse", NULL, "pwm.grad.ism",
        pssm = pssm, aggregate = "lse",
        bidirect = FALSE, strand = 1, extend = TRUE, prior = 0.01
    )
    gvtrack.iterator("g_ism_lse", sshift = 0, eshift = iter_W - 1)
    res <- gextract("g_ism_lse", iterator = 1, intervals = pivot)

    expected <- manual_pwm_grad_ism(pssm, seq_ext, "lse",
        bidirect = FALSE, strand = 1L, prior = 0.01
    )
    expect_equal(res$g_ism_lse[1], expected, tolerance = 1e-5, ignore_attr = TRUE)
})

test_that("pwm.grad.ism lse: matches oracle on mid-strength PSSM", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    L <- 3L
    iter_W <- 15L
    pivot <- gintervals(1, 200, 201)
    seq_ext <- toupper(gseq.extract(gintervals(1, 200, 200 + iter_W + L - 1)))

    pssm <- matrix(
        c(
            0.5, 0.3, 0.1, 0.1,
            0.1, 0.5, 0.3, 0.1,
            0.1, 0.1, 0.5, 0.3
        ),
        L, 4,
        byrow = TRUE,
        dimnames = list(NULL, c("A", "C", "G", "T"))
    )

    gvtrack.create("g_ism_lse", NULL, "pwm.grad.ism",
        pssm = pssm, aggregate = "lse",
        bidirect = FALSE, strand = 1, extend = TRUE, prior = 0.01
    )
    gvtrack.iterator("g_ism_lse", sshift = 0, eshift = iter_W - 1)
    res <- gextract("g_ism_lse", iterator = 1, intervals = pivot)

    expected <- manual_pwm_grad_ism(pssm, seq_ext, "lse",
        bidirect = FALSE, strand = 1L, prior = 0.01
    )
    expect_equal(res$g_ism_lse[1], expected, tolerance = 1e-5, ignore_attr = TRUE)
})

# Task 8: bidirect for ISM

test_that("pwm.grad.ism max: bidirect=TRUE matches oracle on asymmetric PSSM", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    L <- 3L
    iter_W <- 20L
    pivot <- gintervals(1, 200, 201)
    seq_ext <- toupper(gseq.extract(gintervals(1, 200, 200 + iter_W + L - 1)))

    pssm <- matrix(
        c(
            0.6, 0.2, 0.1, 0.1,
            0.1, 0.5, 0.3, 0.1,
            0.1, 0.1, 0.5, 0.3
        ),
        L, 4,
        byrow = TRUE,
        dimnames = list(NULL, c("A", "C", "G", "T"))
    )

    gvtrack.create("g_ism_bid", NULL, "pwm.grad.ism",
        pssm = pssm, aggregate = "max",
        bidirect = TRUE, extend = TRUE, prior = 0.01
    )
    gvtrack.iterator("g_ism_bid", sshift = 0, eshift = iter_W - 1)
    res <- gextract("g_ism_bid", iterator = 1, intervals = pivot)

    expected <- manual_pwm_grad_ism(pssm, seq_ext, "max",
        bidirect = TRUE, prior = 0.01
    )
    expect_equal(res$g_ism_bid[1], expected, tolerance = 1e-5, ignore_attr = TRUE)
})

test_that("pwm.grad.ism lse: bidirect=TRUE matches oracle on asymmetric PSSM", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    L <- 3L
    iter_W <- 20L
    pivot <- gintervals(1, 200, 201)
    seq_ext <- toupper(gseq.extract(gintervals(1, 200, 200 + iter_W + L - 1)))

    pssm <- matrix(
        c(
            0.6, 0.2, 0.1, 0.1,
            0.1, 0.5, 0.3, 0.1,
            0.1, 0.1, 0.5, 0.3
        ),
        L, 4,
        byrow = TRUE,
        dimnames = list(NULL, c("A", "C", "G", "T"))
    )

    gvtrack.create("g_ism_bid_lse", NULL, "pwm.grad.ism",
        pssm = pssm, aggregate = "lse",
        bidirect = TRUE, extend = TRUE, prior = 0.01
    )
    gvtrack.iterator("g_ism_bid_lse", sshift = 0, eshift = iter_W - 1)
    res <- gextract("g_ism_bid_lse", iterator = 1, intervals = pivot)

    expected <- manual_pwm_grad_ism(pssm, seq_ext, "lse",
        bidirect = TRUE, prior = 0.01
    )
    expect_equal(res$g_ism_bid_lse[1], expected, tolerance = 1e-5, ignore_attr = TRUE)
})

test_that("pwm.grad.ism max: strand=-1 matches oracle", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    L <- 3L
    iter_W <- 20L
    pivot <- gintervals(1, 200, 201)
    seq_ext <- toupper(gseq.extract(gintervals(1, 200, 200 + iter_W + L - 1)))

    head_bases <- strsplit(substr(seq_ext, 1, L), "")[[1]]
    fwd_pssm <- .consensus_pssm(head_bases)
    rc_target_pssm <- fwd_pssm[L:1, c(4, 3, 2, 1)]
    dimnames(rc_target_pssm) <- dimnames(fwd_pssm)

    gvtrack.create("g_ism_minus", NULL, "pwm.grad.ism",
        pssm = rc_target_pssm, aggregate = "max",
        bidirect = FALSE, strand = -1, extend = TRUE, prior = 0.01
    )
    gvtrack.iterator("g_ism_minus", sshift = 0, eshift = iter_W - 1)
    res <- gextract("g_ism_minus", iterator = 1, intervals = pivot)

    expected <- manual_pwm_grad_ism(rc_target_pssm, seq_ext, "max",
        bidirect = FALSE, strand = -1L, prior = 0.01
    )
    expect_equal(res$g_ism_minus[1], expected, tolerance = 1e-5, ignore_attr = TRUE)
})

# Task 9: spatial weighting

test_that("pwm.grad: spat_factor = all-1 matches non-spatial", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    L <- 3L
    iter_W <- 20L
    pivot <- gintervals(1, 200, 201)

    head_bases <- strsplit(
        toupper(gseq.extract(gintervals(1, 200, 200 + L))), ""
    )[[1]]
    pssm <- .consensus_pssm(head_bases)

    spat_one <- rep(1, iter_W)

    gvtrack.create("g_no_spat", NULL, "pwm.grad",
        pssm = pssm, aggregate = "lse",
        bidirect = FALSE, strand = 1, extend = TRUE, prior = 0.01
    )
    gvtrack.create("g_spat_one", NULL, "pwm.grad",
        pssm = pssm, spat_factor = spat_one, spat_bin = 1L,
        aggregate = "lse", bidirect = FALSE, strand = 1,
        extend = TRUE, prior = 0.01
    )
    gvtrack.iterator("g_no_spat", sshift = 0, eshift = iter_W - 1)
    gvtrack.iterator("g_spat_one", sshift = 0, eshift = iter_W - 1)

    res <- gextract(c("g_no_spat", "g_spat_one"),
        iterator = 1, intervals = pivot
    )
    expect_equal(res$g_no_spat[1], res$g_spat_one[1],
        tolerance = 1e-5, ignore_attr = TRUE
    )
})

test_that("pwm.grad lse: non-trivial spat_factor matches oracle", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    L <- 3L
    iter_W <- 10L
    pivot <- gintervals(1, 200, 201)
    seq_ext <- toupper(gseq.extract(gintervals(1, 200, 200 + iter_W + L - 1)))

    head_bases <- strsplit(substr(seq_ext, 1, L), "")[[1]]
    pssm <- .consensus_pssm(head_bases)

    # Down-weight late anchors so the LSE concentrates more on early ones,
    # which boosts the head's softmax weight.
    spat_factor <- c(rep(1, 5), rep(0.01, 5))

    gvtrack.create("g_spat", NULL, "pwm.grad",
        pssm = pssm, spat_factor = spat_factor, spat_bin = 1L,
        aggregate = "lse", bidirect = FALSE, strand = 1,
        extend = TRUE, prior = 0.01
    )
    gvtrack.iterator("g_spat", sshift = 0, eshift = iter_W - 1)
    res <- gextract("g_spat", iterator = 1, intervals = pivot)

    expected <- manual_pwm_grad(pssm, seq_ext, "lse",
        bidirect = FALSE, strand = 1L, prior = 0.01,
        spat_factor = spat_factor, spat_bin_size = 1L
    )
    expect_equal(res$g_spat[1], expected, tolerance = 1e-5, ignore_attr = TRUE)
})

test_that("pwm.grad.ism lse: spat_factor matches oracle", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    L <- 3L
    iter_W <- 10L
    pivot <- gintervals(1, 200, 201)
    seq_ext <- toupper(gseq.extract(gintervals(1, 200, 200 + iter_W + L - 1)))

    pssm <- matrix(
        c(
            0.5, 0.3, 0.1, 0.1,
            0.1, 0.5, 0.3, 0.1,
            0.1, 0.1, 0.5, 0.3
        ),
        L, 4,
        byrow = TRUE,
        dimnames = list(NULL, c("A", "C", "G", "T"))
    )

    spat_factor <- c(rep(1, 5), rep(0.05, 5))

    gvtrack.create("g_ism_spat", NULL, "pwm.grad.ism",
        pssm = pssm, spat_factor = spat_factor, spat_bin = 1L,
        aggregate = "lse", bidirect = FALSE, strand = 1,
        extend = TRUE, prior = 0.01
    )
    gvtrack.iterator("g_ism_spat", sshift = 0, eshift = iter_W - 1)
    res <- gextract("g_ism_spat", iterator = 1, intervals = pivot)

    expected <- manual_pwm_grad_ism(pssm, seq_ext, "lse",
        bidirect = FALSE, strand = 1L, prior = 0.01,
        spat_factor = spat_factor, spat_bin_size = 1L
    )
    expect_equal(res$g_ism_spat[1], expected, tolerance = 1e-5, ignore_attr = TRUE)
})

# Task 10: NA / composition / edge-case tests

test_that("pwm.grad: interval shorter than L returns NA", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    L <- 5L
    pivot <- gintervals(1, 200, 201)
    pssm <- .consensus_pssm(c("A", "C", "G", "T", "A"))

    # eshift only adds 2 bp -> total interval length = 3, < L = 5 -> NA.
    gvtrack.create("g_short", NULL, "pwm.grad",
        pssm = pssm, aggregate = "lse",
        bidirect = FALSE, strand = 1, extend = FALSE, prior = 0.01
    )
    gvtrack.iterator("g_short", sshift = 0, eshift = 2)
    res <- gextract("g_short", iterator = 1, intervals = pivot)
    expect_true(is.na(res$g_short[1]))
})

test_that("pwm.grad: composes with pwm in gextract", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    L <- 3L
    iter_W <- 20L
    pivot_block <- gintervals(1, 200, 210)
    seq_ext <- toupper(gseq.extract(gintervals(1, 200, 200 + iter_W + L - 1)))
    head_bases <- strsplit(substr(seq_ext, 1, L), "")[[1]]
    pssm <- .consensus_pssm(head_bases)

    gvtrack.create("score", NULL, "pwm",
        pssm = pssm, bidirect = FALSE, strand = 1,
        extend = TRUE, prior = 0.01
    )
    gvtrack.create("g", NULL, "pwm.grad",
        pssm = pssm, aggregate = "lse",
        bidirect = FALSE, strand = 1, extend = TRUE, prior = 0.01
    )
    gvtrack.iterator("score", sshift = 0, eshift = iter_W - 1)
    gvtrack.iterator("g", sshift = 0, eshift = iter_W - 1)

    res <- gextract(c("score", "g"), iterator = 1, intervals = pivot_block)
    expect_equal(nrow(res), 10L) # 10 1-bp positions
    expect_true(all(res$g >= -1e-7, na.rm = TRUE))
    # Score and gradient should both vary across positions.
    expect_gt(length(unique(res$score)), 1L)
})

test_that("pwm.grad: gscreen filter on gradient", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    L <- 3L
    iter_W <- 20L
    seq_ext <- toupper(gseq.extract(gintervals(1, 200, 200 + iter_W + L - 1)))
    head_bases <- strsplit(substr(seq_ext, 1, L), "")[[1]]
    pssm <- .consensus_pssm(head_bases)

    gvtrack.create("g", NULL, "pwm.grad",
        pssm = pssm, aggregate = "lse",
        bidirect = FALSE, strand = 1, extend = TRUE, prior = 0.01
    )
    gvtrack.iterator("g", sshift = 0, eshift = iter_W - 1)

    out <- gscreen("g > 0.1", iterator = 1, intervals = gintervals(1, 200, 220))
    expect_s3_class(out, "data.frame")
    # All returned intervals must satisfy the predicate.
    if (nrow(out) > 0) {
        vals <- gextract("g", iterator = 1, intervals = out)
        expect_true(all(vals$g > 0.1, na.rm = TRUE))
    }
})

test_that("pwm.grad: multi-chromosome intervals", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    L <- 3L
    iter_W <- 10L
    pssm <- .consensus_pssm(c("A", "C", "G"))

    gvtrack.create("g", NULL, "pwm.grad",
        pssm = pssm, aggregate = "lse",
        bidirect = FALSE, strand = 1, extend = TRUE, prior = 0.01
    )
    gvtrack.iterator("g", sshift = 0, eshift = iter_W - 1)

    intervals <- rbind(
        gintervals(1, 200, 210),
        gintervals(2, 200, 210)
    )
    res <- gextract("g", iterator = 1, intervals = intervals)
    expect_equal(nrow(res), 20L)
    expect_true(all(c(1, 2) %in% res$chrom1 |
        c("chr1", "chr2") %in% res$chrom))
})

test_that("pwm.grad: gsummary aggregation runs without error", {
    remove_all_vtracks()
    withr::defer(remove_all_vtracks())

    L <- 3L
    iter_W <- 10L
    pssm <- .consensus_pssm(c("A", "C", "G"))

    gvtrack.create("g", NULL, "pwm.grad",
        pssm = pssm, aggregate = "max",
        bidirect = FALSE, strand = 1, extend = TRUE, prior = 0.01
    )
    gvtrack.iterator("g", sshift = 0, eshift = iter_W - 1)

    s <- gsummary("g", iterator = 1, intervals = gintervals(1, 200, 250))
    expect_true(is.numeric(s))
    expect_true(length(s) >= 5L)
})
